using Pkg; Pkg.activate(".");
using DifferentialEquations
using DataFrames
using CSV
using SparseArrays
using LinearAlgebra
using JLD2
using Statistics
using Guadex

# =============================================================================
# --- Climate-driven simulations (all temperature projection scenarios) ---
#
# Runs the metacommunity model from the current year to 2045 with a year-by-year
# warming trajectory derived from the guadex_tw water-temperature projections,
# for every GCM x SSP combination (default 11 x 4 = 44 runs).  Each run writes
# the four reporting levels (sampling point, sub-basin, water body, whole basin)
# and the viz/ viewer files.
#
# Time step is daily (1 model time unit = 1 day); states are saved monthly so
# the outputs stay small while annual/seasonal metrics remain available.
#
# Settings come from the [run_climate_scenarios] section of parameters.toml.
# Individual settings can be overridden with the GUADEX_CLIMATE_* environment
# variables (see docs/climate_scenarios.md).
# =============================================================================

# --- Load run parameters (single source of truth: parameters.toml) ---
const _GUADEX_ROOT = isfile(joinpath(@__DIR__, "parameters.jl")) ? (@__DIR__) : dirname(@__DIR__)
include(joinpath(_GUADEX_ROOT, "parameters.jl"))
using .SimulationParameters
const GUADEX_PARAMS = SimulationParameters.load()
SimulationParameters.require_sections(GUADEX_PARAMS, "general", "inputs", "obstacles",
    "species", "subcatchments", "scenarios", "run_climate_scenarios")

const DAYS_PER_YEAR = Int(GUADEX_PARAMS["general"]["days_per_year"])

# Updated 2045 inputs, overridable per run.
const CEDEX_VAR_FILE = get(ENV, "GUADEX_CEDEX_VAR_FILE", GUADEX_PARAMS["inputs"]["cedex_var_file"])
const CEDEX_UTS_FILE = get(ENV, "GUADEX_CEDEX_UTS_FILE", GUADEX_PARAMS["inputs"]["cedex_esc_uts_file"])
const OBSTACLES_FILE = get(ENV, "GUADEX_OBSTACLES_FILE", GUADEX_PARAMS["inputs"]["obstacles_file"])
const OBSTACLE_MODE = Symbol(get(ENV, "GUADEX_OBSTACLE_MODE", GUADEX_PARAMS["obstacles"]["mode"]))
const OBSTACLE_MATCHING_TOLERANCE = parse(Float64, get(ENV, "GUADEX_OBSTACLE_TOLERANCE_M", string(GUADEX_PARAMS["obstacles"]["matching_tolerance_m"])))
const OBSTACLE_PASSABILITY = parse(Float64, get(ENV, "GUADEX_OBSTACLE_PASSABILITY", string(GUADEX_PARAMS["obstacles"]["upstream_passability"])))
const OBSTACLE_DOWNSTREAM_PASSABILITY = parse(Float64, get(ENV, "GUADEX_OBSTACLE_DOWNSTREAM_PASSABILITY", string(GUADEX_PARAMS["obstacles"]["downstream_passability"])))

# C4: dispersal-graph construction method from the [connectivity] section.
const CONNECTIVITY_METHOD = SimulationParameters.connectivity_method(GUADEX_PARAMS)

# E13-E16/E18: opt-in biological-assumption options (all legacy by default).
const BIOLOGICAL_OPTIONS = SimulationParameters.biological_options(GUADEX_PARAMS)

const NATIVE_SPECIES = String.(GUADEX_PARAMS["species"]["native"])
const INVASIVE_SPECIES = String.(GUADEX_PARAMS["species"]["invasive"])
const MIGRATORY_SPECIES = String.(get(GUADEX_PARAMS["species"], "migratory", String[]))
const UPSTREAM_COST = Float64(get(GUADEX_PARAMS["run_climate_scenarios"], "upstream_cost", 0.05))

# --- Climate run settings ---
function _climate_setting(key, default)
    return get(GUADEX_PARAMS["run_climate_scenarios"], key, default)
end

const PROJECTION_FILE = get(ENV, "GUADEX_CLIMATE_PROJECTIONS_FILE",
    _climate_setting("temperature_projections_file", "guadex_tw/outputs/tables/water_temp_future_2045.csv"))
const START_YEAR = parse(Int, get(ENV, "GUADEX_CLIMATE_START_YEAR", string(_climate_setting("start_year", 2026))))
const END_YEAR = parse(Int, get(ENV, "GUADEX_CLIMATE_END_YEAR", string(_climate_setting("end_year", 2045))))
const SAVE_INTERVAL_DAYS = parse(Float64, get(ENV, "GUADEX_CLIMATE_SAVE_DAYS", string(_climate_setting("save_interval_days", 365.0 / 12.0))))
const ELEVATION_SCALING = lowercase(get(ENV, "GUADEX_CLIMATE_ELEVATION_SCALING",
    string(_climate_setting("elevation_scaling", false)))) in ("1", "true", "yes")
const PRESENCE_THRESHOLD = Float64(_climate_setting("presence_threshold", 0.1))
const MAX_RUNS = parse(Int, get(ENV, "GUADEX_CLIMATE_MAX_RUNS", "0"))  # 0 = no cap (testing aid)
const FORCE_RERUN = lowercase(get(ENV, "GUADEX_CLIMATE_FORCE", "0")) in ("1", "true", "yes")
const MAKE_FIGURES = lowercase(get(ENV, "GUADEX_CLIMATE_PLOT",
    string(_climate_setting("make_figures", true)))) in ("1", "true", "yes")

# --- WP1-WP5 settings ---
const FORCING_MODE = lowercase(get(ENV, "GUADEX_CLIMATE_FORCING_MODE",
    string(_climate_setting("forcing_mode", "annual_mean"))))
const DAILY_FORCING_FILE = get(ENV, "GUADEX_CLIMATE_DAILY_FILE",
    string(_climate_setting("daily_forcing_file", "")))
const DAILY_FORCING_PER_GCM = lowercase(get(ENV, "GUADEX_CLIMATE_DAILY_PER_GCM",
    string(_climate_setting("daily_forcing_per_gcm", false)))) in ("1", "true", "yes")

# The per-GCM wide files written by guadex_tw/scripts/13_project_guadex_sites.py
# are named `<base>_<scenario>_<gcm>.csv`; they carry the GCM spread that the
# ensemble-median product deliberately removes.
function per_gcm_forcing_path(base::AbstractString, scenario::AbstractString, gcm::AbstractString)
    root, ext = splitext(base)
    return string(root, "_", scenario, "_", gcm, ext)
end
const BASELINE_PERIOD_START = parse(Int, get(ENV, "GUADEX_CLIMATE_BASELINE_START",
    string(_climate_setting("baseline_period_start", 1986))))
const BASELINE_PERIOD_END = parse(Int, get(ENV, "GUADEX_CLIMATE_BASELINE_END",
    string(_climate_setting("baseline_period_end", 2005))))
const SPIN_UP_ENABLED = lowercase(get(ENV, "GUADEX_CLIMATE_SPIN_UP",
    string(_climate_setting("spin_up", false)))) in ("1", "true", "yes")
const SPIN_UP_MAX_YEARS = parse(Int, get(ENV, "GUADEX_CLIMATE_SPIN_UP_YEARS",
    string(_climate_setting("spin_up_max_years", 50))))
const SPIN_UP_TOL = parse(Float64, get(ENV, "GUADEX_CLIMATE_SPIN_UP_TOL",
    string(_climate_setting("spin_up_tol", 1.0e-6))))
# E4: the composition (95th percentile of |Δ log N| over site × species cells)
# tolerance and the minimum number of burn-in years.  The default criterion is
# the stricter two-criteria rule; a run can opt into the legacy basin-only rule.
const SPIN_UP_COMPOSITION_TOL = parse(Float64, get(ENV, "GUADEX_CLIMATE_SPIN_UP_COMPOSITION_TOL",
    string(_climate_setting("spin_up_composition_tol", 1.0e-2))))
const SPIN_UP_MIN_YEARS = parse(Int, get(ENV, "GUADEX_CLIMATE_SPIN_UP_MIN_YEARS",
    string(_climate_setting("spin_up_min_years", 10))))
const SPIN_UP_CRITERION = Symbol(lowercase(get(ENV, "GUADEX_CLIMATE_SPIN_UP_CRITERION",
    string(_climate_setting("spin_up_criterion", "both")))))
const SPIN_UP_PROGRESS_EVERY = parse(Int, get(ENV, "GUADEX_CLIMATE_SPIN_UP_PROGRESS",
    string(_climate_setting("spin_up_progress_every", 50))))
# C3: projection route.  "full_burnin" (the default when the key is absent — so
# every other configuration keeps the E4 convergence machinery) integrates to the
# two-criteria stop rule above.  "interim_observed" is the team's interim route:
# start from the OBSERVED community and integrate only INTERIM_SPIN_UP_YEARS
# baseline years, reporting every scenario relative to the matched control.
const PROJECTION_ROUTE = projection_route(get(ENV, "GUADEX_CLIMATE_PROJECTION_ROUTE",
    string(_climate_setting("projection_route", ""))))
const INTERIM_SPIN_UP_YEARS = parse(Int, get(ENV, "GUADEX_CLIMATE_INTERIM_SPIN_UP_YEARS",
    string(_climate_setting("interim_spin_up_years", 3))))
const K_SCALING = parse(Float64, get(ENV, "GUADEX_CLIMATE_K_SCALING",
    string(_climate_setting("carrying_capacity_scaling", 1.0))))
# Multiplier linking observed total density to the base carrying capacity
# (`build_carrying_capacity`).  Legacy convention is 10x observed; set to 1.0 to
# treat the observed snapshot as the capacity level.  The effective multiplier is
# K_BASE_SCALING * K_SCALING.
const K_BASE_SCALING = parse(Float64, get(ENV, "GUADEX_CLIMATE_K_BASE",
    string(_climate_setting("carrying_capacity_base_scaling", 10.0))))
const OPTIMA_FRACTION = parse(Float64, get(ENV, "GUADEX_CLIMATE_OPTIMA_FRACTION",
    string(_climate_setting("thermal_optima_fraction", 0.5))))
const QUASI_EXTINCTION_Q = parse(Float64, get(ENV, "GUADEX_CLIMATE_QE_Q",
    string(_climate_setting("quasi_extinction_q", 0.1))))
const QUASI_EXTINCTION_PERSISTENCE = parse(Int, get(ENV, "GUADEX_CLIMATE_QE_PERSISTENCE",
    string(_climate_setting("quasi_extinction_persistence", 3))))

const STRESS_PARAMS = get(GUADEX_PARAMS, "temperature_stress", Dict{String,Any}())
const HEAT_STRESS_ENABLED = lowercase(get(ENV, "GUADEX_CLIMATE_HEAT_STRESS",
    string(get(STRESS_PARAMS, "enabled", false)))) in ("1", "true", "yes")
const HEAT_STRESS_K = parse(Float64, get(ENV, "GUADEX_CLIMATE_HEAT_STRESS_K",
    string(get(STRESS_PARAMS, "k", 0.0))))
const HEAT_STRESS_CALIBRATE = lowercase(get(ENV, "GUADEX_CLIMATE_HEAT_STRESS_CALIBRATE",
    string(get(STRESS_PARAMS, "calibrate", false)))) in ("1", "true", "yes")
const HEAT_STRESS_MAX_LOSS = parse(Float64, get(ENV, "GUADEX_CLIMATE_HEAT_STRESS_MAX_LOSS",
    string(get(STRESS_PARAMS, "max_annual_loss", 0.05))))

const CONTROL_RUN = lowercase(get(ENV, "GUADEX_CLIMATE_CONTROL",
    string(_climate_setting("control_run", false)))) in ("1", "true", "yes")
# C3 route validation: the interim route mandates a matched no-warming control
# (results are deltas) and a short fixed relaxation, so both are checked loudly
# rather than silently ignored.  The full-burn-in default is unconstrained here.
const INTERIM_ROUTE = PROJECTION_ROUTE === :interim_observed
if INTERIM_ROUTE
    (0 <= INTERIM_SPIN_UP_YEARS <= 5) ||
        error("interim_spin_up_years must be between 0 and 5 (got $INTERIM_SPIN_UP_YEARS)")
    CONTROL_RUN ||
        error("projection_route = \"interim_observed\" requires control_run = true " *
              "(scenario results are reported relative to the matched no-warming control)")
end
const CLIMATE_SCENARIOS = begin
    scenarios = String.(get(GUADEX_PARAMS["run_climate_scenarios"], "scenarios",
        ["ssp126", "ssp245", "ssp370", "ssp585"]))
    # A no-warming control (baseline climatology over the horizon) is reported
    # alongside the ensemble so scenario effects can be measured against it (WP4).
    CONTROL_RUN && !("control" in scenarios) ? vcat(scenarios, ["control"]) : scenarios
end
const DEFAULT_GCMS = ["ACCESS-CM2", "CMCC-CM2-SR5", "CNRM-ESM2-1", "EC-Earth3-Veg",
    "IITM-ESM", "KACE-1-0-G", "MIROC6", "MPI-ESM1-2-HR", "MRI-ESM2-0",
    "NorESM2-MM", "UKESM1-0-LL"]
const CLIMATE_GCMS = begin
    override = get(ENV, "GUADEX_CLIMATE_GCMS", "")
    isempty(override) ? String.(get(GUADEX_PARAMS["run_climate_scenarios"], "gcms", DEFAULT_GCMS)) :
        String.(split(override, ","))
end

const SIMULATION_YEARS = END_YEAR - START_YEAR + 1
const T_END = Float64(SIMULATION_YEARS * DAYS_PER_YEAR)
const YEAR_OFFSETS = collect(1:SIMULATION_YEARS)
const YEAR_LABELS = collect(START_YEAR:END_YEAR)

if get(ENV, "GUADEX_CONFIG_ONLY", "0") == "1"
    println("[parameters] climate-scenario configuration smoke test passed")
    println("  scenarios: $(join(CLIMATE_SCENARIOS, ", "))")
    println("  gcms: $(length(CLIMATE_GCMS))")
    println("  years: $START_YEAR-$END_YEAR (daily step, $(SAVE_INTERVAL_DAYS)-day saves)")
    println("  projections: $PROJECTION_FILE")
    println("  forcing mode: $FORCING_MODE")
    println("  heat stress: enabled=$HEAT_STRESS_ENABLED k=$HEAT_STRESS_K " *
            "calibrate=$HEAT_STRESS_CALIBRATE")
    println("  projection route: :$PROJECTION_ROUTE" *
            (INTERIM_ROUTE ? " (interim_spin_up_years=$INTERIM_SPIN_UP_YEARS; control_run=$CONTROL_RUN)" : ""))
    println("  spin-up: $SPIN_UP_ENABLED (max $(SPIN_UP_MAX_YEARS) yr, criterion=:$SPIN_UP_CRITERION, " *
            "composition_tol=$SPIN_UP_COMPOSITION_TOL, min_years=$SPIN_UP_MIN_YEARS)")
    println("  carrying-capacity base scaling: $(K_BASE_SCALING)x; WP4 scaling: $(K_SCALING)x; " *
            "effective: $(K_BASE_SCALING * K_SCALING)x")
    println("  optimum fraction: $OPTIMA_FRACTION")
    println("  biological options: absence=$(BIOLOGICAL_OPTIONS.absence_growth_fraction), " *
            "pool=$(BIOLOGICAL_OPTIONS.pool_capacity_mode), " *
            "dificil=$(BIOLOGICAL_OPTIONS.fishless_dificil_capacity), " *
            "nonrepro=$(BIOLOGICAL_OPTIONS.nonreproducing_local_growth), " *
            "fishfarm=$(BIOLOGICAL_OPTIONS.exclude_fishfarm_eel_records), " *
            "salinity=$(BIOLOGICAL_OPTIONS.salinity_envelope)")
    exit(0)
end

# =============================================================================
# --- Main ---
# =============================================================================

println("="^70)
println("Climate-driven simulations: $START_YEAR-$END_YEAR")
println("  scenarios: $(join(CLIMATE_SCENARIOS, ", "))")
println("  GCMs:      $(length(CLIMATE_GCMS))")
println("  save:      every $(SAVE_INTERVAL_DAYS) days (monthly by default)")
println("  elevation scaling: $ELEVATION_SCALING")
println("  forcing:   $FORCING_MODE")
println("  heat stress: $HEAT_STRESS_ENABLED (k=$(HEAT_STRESS_K))")
println("  spin-up:   $SPIN_UP_ENABLED (criterion=:$SPIN_UP_CRITERION)")
println("  route:     :$PROJECTION_ROUTE" *
        (INTERIM_ROUTE ? " ($INTERIM_SPIN_UP_YEARS baseline year(s); scenario-control deltas)" : ""))
println("  K multiplier: base $(K_BASE_SCALING)x * WP4 $(K_SCALING)x = " *
        "$(K_BASE_SCALING * K_SCALING)x observed density")
println("="^70)

projections = load_temperature_projections(PROJECTION_FILE)
available_scenarios = Set(String.(projections.scenario))
available_gcms = Set(String.(projections.gcm))

println("\nLoading base data (once)...")
data_base = prepare_ode_data(
    upstream_cost = UPSTREAM_COST,
    cedex_var_file = CEDEX_VAR_FILE,
    cedex_esc_uts_file = CEDEX_UTS_FILE,
    obstacles_file = OBSTACLES_FILE,
    obstacle_mode = OBSTACLE_MODE,
    obstacle_matching_tolerance = OBSTACLE_MATCHING_TOLERANCE,
    obstacle_passability = OBSTACLE_PASSABILITY,
    obstacle_downstream_passability = OBSTACLE_DOWNSTREAM_PASSABILITY,
    connectivity_method = CONNECTIVITY_METHOD,
    carrying_capacity_base_scaling = K_BASE_SCALING,
    absence_growth_fraction = BIOLOGICAL_OPTIONS.absence_growth_fraction,
    pool_capacity_mode = BIOLOGICAL_OPTIONS.pool_capacity_mode,
    fishless_dificil_capacity = BIOLOGICAL_OPTIONS.fishless_dificil_capacity,
    nonreproducing_local_growth = BIOLOGICAL_OPTIONS.nonreproducing_local_growth,
    exclude_fishfarm_eel_records = BIOLOGICAL_OPTIONS.exclude_fishfarm_eel_records,
    salinity_envelope = BIOLOGICAL_OPTIONS.salinity_envelope
)

n_sites = data_base.params.n_sites

# Initial conditions from the observed density matrix.  Align by site code
# (minor #3): `observed_density_matrix` joins on CODIGO and errors if a model
# site has no density row, instead of assuming row position matches site order.
u0 = observed_density_matrix(data_base.density_df, data_base.sites, data_base.species)
u0_flat = vec(u0)

# ---------------------------------------------------------------------------
# WP3/WP4: thermal optimum sweep, heat-stress slope, K scaling, daily forcing
# and the baseline spin-up.
# ---------------------------------------------------------------------------
params = data_base.params

if OPTIMA_FRACTION != 0.5
    optima = optimum_sweep_optima(data_base.species_chars_df, data_base.species;
        fraction=OPTIMA_FRACTION)
    params = set_thermal_optima(params, optima)
    println("WP3 optimum sweep: fraction=$OPTIMA_FRACTION -> $optima")
end

if K_SCALING != 1.0
    params = scale_carrying_capacity(params, K_SCALING)
    println("WP4 carrying-capacity scaling: $(K_SCALING)x")
end
println("Effective carrying-capacity multiplier vs observed density: " *
        "$(K_BASE_SCALING * K_SCALING)x")

# WP1: per-site daily forcing, loaded only in the daily mode.  A separate
# "historical" series provides the 1986-2005 baseline when present.
daily_forcing = nothing
baseline_temps = nothing
baseline_dates = nothing
forcing_matrix(df, sites; scenario) = is_wide_daily_forcing(df) ?
    wide_forcing_matrix(df, sites; scenario=scenario) :
    daily_forcing_matrix(df, sites; scenario=scenario)

if FORCING_MODE == "daily"
    isempty(DAILY_FORCING_FILE) && error("forcing_mode = daily requires daily_forcing_file")
    daily_forcing = load_daily_forcing_any(DAILY_FORCING_FILE)
    println("WP1 forcing layout: $(is_wide_daily_forcing(daily_forcing) ? "wide" : "long")")
    if "historical" in Set(String.(daily_forcing.scenario))
        baseline_dates, baseline_temps = forcing_matrix(daily_forcing, data_base.sites;
            scenario="historical")
        println("WP1 baseline series: $(length(baseline_dates)) days from 'historical'")
    else
        @warn "daily forcing has no 'historical' scenario; baseline will be taken from the " *
              "scenario series window $BASELINE_PERIOD_START-$BASELINE_PERIOD_END"
    end
end

# WP3: baseline-viability calibration of the shared heat-stress slope.
k_heat = HEAT_STRESS_ENABLED ? HEAT_STRESS_K : 0.0
if HEAT_STRESS_ENABLED && HEAT_STRESS_CALIBRATE
    if baseline_temps === nothing
        @warn "heat-stress calibration needs a daily baseline series; using k=$HEAT_STRESS_K"
    else
        energies = Float64[]
        for s in 1:params.n_species, i in 1:params.n_sites
            push!(energies, exceedance_energy(@view(baseline_temps[i, :]),
                params.thermal_upper_limits[s]))
        end
        k_heat = calibrate_heat_stress_rate(energies; max_annual_loss=HEAT_STRESS_MAX_LOSS)
        println("WP3 heat-stress calibration: k=$k_heat " *
                "(max annual loss $(HEAT_STRESS_MAX_LOSS))")
    end
end
params = set_heat_stress_rate(params, k_heat)
println("Heat-stress mortality: $(k_heat > 0 ? "k=$k_heat" : "disabled")")
data_base = merge(data_base, (params=params,))

# Baseline seasonal schedule shared by both spin-up routes (nothing => annual
# mean).  Spin up under the seasonal climatology, not the annual mean: the two
# equilibria differ and only the seasonal one matches the runs.
spin_schedule = nothing
if FORCING_MODE == "daily" && baseline_temps !== nothing
    spin_schedule = baseline_climatology_schedule(baseline_temps, baseline_dates, 1;
        days_per_year=Float64(DAYS_PER_YEAR))
end

# WP4/C3: initial state.  The projection route selects the strategy:
#   * "full_burnin" (default): integrate to the E4 two-criteria stop rule.
#   * "interim_observed": observed start + EXACTLY INTERIM_SPIN_UP_YEARS baseline
#     years, bypassing the E4 stop rule (this route only), with results reported
#     as scenario - control deltas.
baseline_species_density = nothing
spin_result = nothing
spin_performed = false
if INTERIM_ROUTE
    println("C3 interim projection route: observed start + $(INTERIM_SPIN_UP_YEARS) " *
            "baseline year(s) under $(spin_schedule === nothing ? "annual-mean" : "seasonal") " *
            "forcing (E4 stop rule bypassed; reporting scenario - control)")
    interim = interim_observed_spin_up(params; initial_state=u0_flat,
        schedule=spin_schedule, years=INTERIM_SPIN_UP_YEARS,
        days_per_year=Float64(DAYS_PER_YEAR), progress_every=SPIN_UP_PROGRESS_EVERY)
    u0_flat = interim.state
    spin_result = interim.spin
    spin_performed = interim.years > 0
    if spin_performed
        baseline_species_density = reshape(u0_flat, n_sites, params.n_species)
    end
    println("  interim spin-up: years=$(interim.years) (fixed, no convergence stop)")
elseif SPIN_UP_ENABLED
    println("WP4 spin-up under $(spin_schedule === nothing ? "annual-mean" : "seasonal") " *
            "baseline forcing (max $(SPIN_UP_MAX_YEARS) yr, tol=$SPIN_UP_TOL, " *
            "composition_tol=$SPIN_UP_COMPOSITION_TOL, min_years=$SPIN_UP_MIN_YEARS, " *
            "criterion=:$SPIN_UP_CRITERION)...")
    spin_result = spin_up(params; initial_state=u0_flat, schedule=spin_schedule,
        days_per_year=Float64(DAYS_PER_YEAR), max_years=SPIN_UP_MAX_YEARS, tol=SPIN_UP_TOL,
        composition_tol=SPIN_UP_COMPOSITION_TOL, min_years=SPIN_UP_MIN_YEARS,
        criterion=SPIN_UP_CRITERION, progress_every=SPIN_UP_PROGRESS_EVERY)
    u0_flat = spin_result.state
    baseline_species_density = reshape(u0_flat, n_sites, params.n_species)
    spin_performed = true
    println("  spin-up: years=$(spin_result.years), converged=$(spin_result.converged), " *
            "basin change=$(spin_result.last_basin_change), comp_q95=$(spin_result.last_composition_change) " *
            "(all_cells=$(spin_result.last_composition_change_all_cells), " *
            "active=$(spin_result.composition_active_cells)), " *
            "q95=$(spin_result.last_q95_change), max=$(spin_result.last_relative_change), " *
            "criteria=$(join(String.(spin_result.criteria), "+"))")
end

saveat = unique(vcat(collect(0.0:SAVE_INTERVAL_DAYS:T_END), T_END))

function positivity_condition(u, t, integrator)
    return any(x -> x < 0, u)
end
function positivity_affect!(integrator)
    integrator.u .= max.(integrator.u, 0.0)
end
# Clamp negatives but do NOT save a state on every trigger: under daily forcing
# and heat stress the callback fires many thousands of times, and
# `save_positions=(false, true)` would store a ~1.8 GB solution per run.  The
# monthly `saveat` grid already provides the reporting snapshots.
positivity_cb = DiscreteCallback(positivity_condition, positivity_affect!; save_positions=(false, false))

# Restrict path components to a safe character set so projection labels can
# never escape the results directory.
function safe_component(value)
    safe = replace(string(value), r"[^A-Za-z0-9._-]" => "_")
    safe = replace(safe, ".." => "__")
    return isempty(safe) ? "unnamed" : safe
end

function basin_row_from_csv(path, year)
    isfile(path) || return (NaN, NaN, NaN)
    df = CSV.read(path, DataFrame)
    ("year" in names(df)) || return (NaN, NaN, NaN)
    sub = df[df.year .== year, :]
    isempty(sub) && return (NaN, NaN, NaN)
    r = sub[1, :]
    return (r.mean_native_richness, r.mean_realised_richness_loss, r.mean_total_biomass)
end

# --- C7: persisted realised applied forcing ---------------------------------
# The annual-mean `warming` matrix is what the ODE is forced with, so its
# site-and-window means are the realised forcing (unlike `warming_end_degc`,
# which is a proxy read from the 6-gauge projection table).
function _json_number_field(text::AbstractString, key::AbstractString)
    m = match(Regex("\"" * key *
        "\"\\s*:\\s*(-?[0-9]+(?:\\.[0-9]+)?(?:[eE][-+]?[0-9]+)?|null)"), text)
    (m === nothing || m.captures[1] == "null") && return NaN
    return parse(Float64, m.captures[1])
end

function _json_int_field(text::AbstractString, key::AbstractString)
    value = _json_number_field(text, key)
    return isfinite(value) ? Int(round(value)) : 0
end

"""
    persisted_realised_forcing(run_dir)

Read the realised-forcing fields back from an existing run's
`export/run_metadata.json` (used on the digest-match skip path so the index can
be completed without re-running).  Returns `(early, late, n_sites)`, with `NaN`
when the fields are absent (runs written before C7).
"""
function persisted_realised_forcing(run_dir::AbstractString)
    path = joinpath(run_dir, "export", "run_metadata.json")
    isfile(path) || return (NaN, NaN, 0)
    text = read(path, String)
    return (_json_number_field(text, "realised_warming_2026_2045_mean_degc"),
        _json_number_field(text, "realised_warming_2036_2045_mean_degc"),
        _json_int_field(text, "realised_warming_n_sites"))
end

function ensure_index_columns!(df::DataFrame)
    for col in (:scenario, :gcm, :run_dir)
        hasproperty(df, col) || (df[!, col] = fill("", nrow(df)))
    end
    # Existing rows were produced from `basin_warming_curve`, so label the
    # retained proxy explicitly rather than leaving the source blank.
    hasproperty(df, :warming_end_source) ||
        (df[!, :warming_end_source] = fill("proxy_basin_warming_curve", nrow(df)))
    for col in (:start_year, :end_year, :realised_warming_n_sites)
        hasproperty(df, col) || (df[!, col] = fill(0, nrow(df)))
    end
    for col in (:warming_end_degc, :realised_warming_2026_2045_mean_degc,
                :realised_warming_2036_2045_mean_degc, :basin_native_richness,
                :basin_realised_richness_loss, :basin_total_biomass)
        hasproperty(df, col) || (df[!, col] = fill(NaN, nrow(df)))
    end
    return df
end

function index_row(scenario, gcm, warming_end, realised, richness, risk, biomass, run_dir)
    return (
        scenario = String(scenario), gcm = String(gcm),
        start_year = START_YEAR, end_year = END_YEAR,
        warming_end_degc = Float64(warming_end),
        warming_end_source = "proxy_basin_warming_curve",
        realised_warming_2026_2045_mean_degc = Float64(realised[1]),
        realised_warming_2036_2045_mean_degc = Float64(realised[2]),
        realised_warming_n_sites = Int(realised[3]),
        basin_native_richness = richness,
        basin_realised_richness_loss = risk,
        basin_total_biomass = biomass,
        run_dir = String(run_dir))
end

const CLIMATE_OUTPUT_DIR = get(ENV, "GUADEX_CLIMATE_OUTPUT_DIR",
    string(_climate_setting("output_dir", joinpath("results", "climate_scenarios"))))
base_output_dir = CLIMATE_OUTPUT_DIR
mkpath(base_output_dir)
index_path = joinpath(base_output_dir, "runs_index.csv")
index_df = isfile(index_path) ? CSV.read(index_path, DataFrame) : DataFrame(
    scenario=String[], gcm=String[], start_year=Int[], end_year=Int[],
    warming_end_degc=Float64[], warming_end_source=String[],
    realised_warming_2026_2045_mean_degc=Float64[],
    realised_warming_2036_2045_mean_degc=Float64[],
    realised_warming_n_sites=Int[],
    basin_native_richness=Float64[], basin_realised_richness_loss=Float64[],
    basin_total_biomass=Float64[], run_dir=String[])
# Old indexes written before C7 lack the realised-forcing columns; add them with
# NaN sentinels so appends keep working and the absence is explicit.
ensure_index_columns!(index_df)

run_number = 0
total_runs = length(CLIMATE_SCENARIOS) * length(CLIMATE_GCMS)

for scenario in CLIMATE_SCENARIOS
    is_control = scenario == "control"
    is_control || (scenario in available_scenarios ||
        (@warn "scenario $scenario absent from projections; skipping"; continue))
    scenario_gcms = is_control ? ["baseline"] : CLIMATE_GCMS
    for gcm in scenario_gcms
        is_control || (gcm in available_gcms ||
            (@warn "GCM $gcm absent from projections; skipping"; continue))
        global run_number += 1
        if MAX_RUNS > 0 && run_number > MAX_RUNS
            println("[cap] GUADEX_CLIMATE_MAX_RUNS=$MAX_RUNS reached; stopping after $MAX_RUNS run(s)")
            break
        end

        run_dir = joinpath(base_output_dir, safe_component(scenario), safe_component(gcm))
        curve_years, curve = is_control ?
            (collect(START_YEAR:END_YEAR), zeros(Float64, SIMULATION_YEARS)) :
            basin_warming_curve(projections, scenario, gcm, START_YEAR, END_YEAR)

        # E9/E10: identity of this run — resolved model state, resolved run
        # options and code version.  Used both to decide whether an existing
        # output can be trusted and to record provenance.
        per_gcm_digest = (FORCING_MODE == "daily" && DAILY_FORCING_PER_GCM) ?
            file_fingerprint(per_gcm_forcing_path(DAILY_FORCING_FILE, scenario, gcm)) : "na"
        resolved_config = Dict{String,Any}(
            "script" => "run_climate_scenarios.jl",
            "scenario" => scenario, "gcm" => gcm,
            "start_year" => START_YEAR, "end_year" => END_YEAR,
            "simulation_years" => SIMULATION_YEARS,
            "save_interval_days" => SAVE_INTERVAL_DAYS,
            "forcing_mode" => FORCING_MODE,
            "daily_forcing_file" => isempty(DAILY_FORCING_FILE) ? "" : basename(DAILY_FORCING_FILE),
            "daily_forcing_fingerprint" => isempty(DAILY_FORCING_FILE) ? "na" : file_fingerprint(DAILY_FORCING_FILE),
            "daily_forcing_per_gcm" => DAILY_FORCING_PER_GCM,
            "per_gcm_forcing_fingerprint" => per_gcm_digest,
            "baseline_period_start" => BASELINE_PERIOD_START,
            "baseline_period_end" => BASELINE_PERIOD_END,
            "spin_up" => SPIN_UP_ENABLED,
            "spin_up_max_years" => SPIN_UP_MAX_YEARS,
            "spin_up_tol" => SPIN_UP_TOL,
            "spin_up_composition_tol" => SPIN_UP_COMPOSITION_TOL,
            "spin_up_composition_active_floor" => Float64(spin_up_composition_active_floor()),
            "spin_up_min_years" => SPIN_UP_MIN_YEARS,
            "spin_up_criterion" => string(SPIN_UP_CRITERION),
            "projection_route" => string(PROJECTION_ROUTE),
            "interim_spin_up_years" => INTERIM_ROUTE ? INTERIM_SPIN_UP_YEARS : 0,
            "carrying_capacity_base_scaling" => K_BASE_SCALING,
            "carrying_capacity_scaling" => K_SCALING,
            "thermal_optima_fraction" => OPTIMA_FRACTION,
            "heat_stress_enabled" => HEAT_STRESS_ENABLED,
            "heat_stress_calibrate" => HEAT_STRESS_CALIBRATE,
            "heat_stress_max_loss" => HEAT_STRESS_MAX_LOSS,
            "upstream_cost" => UPSTREAM_COST,
            "obstacle_mode" => string(OBSTACLE_MODE),
            "connectivity_method" => string(CONNECTIVITY_METHOD),
            "obstacle_matching_tolerance_m" => OBSTACLE_MATCHING_TOLERANCE,
            "obstacle_passability" => OBSTACLE_PASSABILITY,
            "obstacle_downstream_passability" => OBSTACLE_DOWNSTREAM_PASSABILITY,
            "control_run" => CONTROL_RUN,
            "elevation_scaling" => ELEVATION_SCALING,
            "presence_threshold" => PRESENCE_THRESHOLD,
            "quasi_extinction_q" => QUASI_EXTINCTION_Q,
            "quasi_extinction_persistence" => QUASI_EXTINCTION_PERSISTENCE,
            "projections_file" => basename(PROJECTION_FILE),
            "projections_fingerprint" => file_fingerprint(PROJECTION_FILE),
            "warming_curve" => stable_digest(join(string.(vec(curve)), ",")),
            "code_version" => code_version(),
        )
        # E13-E16/E18: the resolved biological-assumption options are part of the
        # run identity, so changing one invalidates a cached run.
        merge!(resolved_config, SimulationParameters.biological_options_dict(GUADEX_PARAMS))
        run_digest = parameter_digest(data_base.params; options=resolved_config)

        # Resume support: reuse an existing output only when its stored digest
        # matches the current model state, options and code version; otherwise
        # re-run (a changed parameter, forcing or code version invalidates it).
        digest_path = joinpath(run_dir, "run_digest.txt")
        jld_path = joinpath(run_dir, "simulation_output.jld2")
        stored_digest = isfile(digest_path) ? strip(read(digest_path, String)) : ""
        if !FORCE_RERUN && isfile(jld_path) && stored_digest == run_digest
            println("\n[$run_number/$total_runs] $scenario / $gcm -> skipped (digest match)")
            if !any((index_df.scenario .== scenario) .& (index_df.gcm .== gcm))
                richness, risk, biomass = basin_row_from_csv(
                    joinpath(run_dir, "export", "levels", "level_basin.csv"), END_YEAR)
                # C7: recover the realised forcing already persisted in the run
                # metadata (NaN for runs written before C7; they need a re-run).
                realised = persisted_realised_forcing(run_dir)
                push!(index_df, index_row(scenario, gcm, curve[end], realised,
                    richness, risk, biomass, run_dir))
                CSV.write(index_path, index_df)
            end
            continue
        elseif !FORCE_RERUN && isfile(jld_path)
            println("\n[$run_number/$total_runs] $scenario / $gcm -> recomputing " *
                    "(no digest / digest mismatch)")
        end

        mkpath(run_dir)
        println("\n[$run_number/$total_runs] $scenario / $gcm -> $run_dir")

        warming = warming_matrix(n_sites, curve_years, curve;
            elevations=data_base.elevations, elevation_scaling=ELEVATION_SCALING)
        forcing_temps = nothing
        forcing_dates = nothing

        schedule = if FORCING_MODE == "daily" && !is_control
            src = daily_forcing
            if DAILY_FORCING_PER_GCM
                path = per_gcm_forcing_path(DAILY_FORCING_FILE, scenario, gcm)
                isfile(path) || error("per-GCM daily forcing not found: $path")
                src = load_daily_forcing_any(path)
            end
            forcing_dates, forcing_temps = forcing_matrix(src, data_base.sites;
                scenario=scenario)
            # C5/E1: reference the anomaly to the SAME per-site 1986-2005 mean
            # that defines the model level (`params.temperatures`, the corrected
            # `tw_baseline_mean`).  Then `level + anomaly` reconstructs the
            # corrected daily series exactly, with no offset (no per-GCM
            # historical-mean shift and no double correction).
            sched, _ = daily_temperature_schedule(forcing_temps;
                dates=forcing_dates,
                baseline_start=BASELINE_PERIOD_START, baseline_end=BASELINE_PERIOD_END,
                baseline_means=data_base.params.temperatures)
            years_out, warming_out = annual_mean_deltas_by_year(sched, forcing_dates)
            columns = [findfirst(==(y), years_out) for y in YEAR_LABELS]
            any(isnothing, columns) &&
                error("daily forcing ($scenario) does not cover every year in $START_YEAR-$END_YEAR")
            warming = warming_out[:, [something(c) for c in columns]]
            sched
        elseif is_control && FORCING_MODE == "daily" && baseline_temps !== nothing
            # No-warming control under the baseline daily climatology.
            forcing_dates, forcing_temps = baseline_dates, baseline_temps
            sched = baseline_climatology_schedule(baseline_temps, baseline_dates,
                SIMULATION_YEARS; days_per_year=Float64(DAYS_PER_YEAR))
            warming = annual_mean_deltas(sched; days_per_year=Float64(DAYS_PER_YEAR))
            sched
        else
            annual_mean_schedule(warming; days_per_year=Float64(DAYS_PER_YEAR))
        end
        # C7: realised applied warming anomaly over the 775 sites, from the
        # annual-mean `warming` matrix the ODE is forced with (2026-2045 and
        # 2036-2045).  This is the regressor for the dose-response analyses.
        realised = realised_warming_anomaly(warming, YEAR_LABELS;
            early_window=(START_YEAR, END_YEAR),
            late_window=(max(START_YEAR, END_YEAR - 9), END_YEAR))
        scheduled_params = ScheduledMetacommunityParams(data_base.params, schedule)

        prob = ODEProblem(metacommunity_ode_scheduled!, u0_flat, (0.0, T_END), scheduled_params)
        sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, saveat=saveat, callback=positivity_cb)
        println("  solved: $(length(sol.t)) snapshots, final basin delta = $(round(curve[end], digits=3)) degC")
        println("  realised applied warming: $(round(realised.mean_2026_2045, digits=3)) degC " *
                "($(START_YEAR)-$(END_YEAR)), $(round(realised.mean_2036_2045, digits=3)) degC " *
                "($(max(START_YEAR, END_YEAR - 9))-$(END_YEAR)) over $(realised.n_sites) sites")

        jldsave(joinpath(run_dir, "simulation_output.jld2");
            sol_t=sol.t, sol_u=sol.u,
            sites=data_base.sites, species=data_base.species,
            scenario=scenario, gcm=gcm,
            start_year=START_YEAR, end_year=END_YEAR,
            years=YEAR_LABELS, warming=warming,
            temperature_baseline=data_base.params.temperatures,
            habitat_suitability=data_base.params.habitat_suitability,
            save_interval_days=SAVE_INTERVAL_DAYS,
            elevation_scaling=ELEVATION_SCALING,
            forcing_mode=FORCING_MODE,
            heat_stress_k=k_heat,
            thermal_upper_limits=data_base.params.thermal_upper_limits,
            spin_up=spin_performed,
            projection_route=string(PROJECTION_ROUTE),
            interim_spin_up_years=INTERIM_ROUTE ? INTERIM_SPIN_UP_YEARS : 0,
            carrying_capacity_scaling=K_SCALING,
            carrying_capacity_base_scaling=K_BASE_SCALING,
            spin_up_criterion=string(SPIN_UP_CRITERION),
            spin_up_total_biomass_change=spin_performed ? spin_result.last_basin_change : nothing,
            spin_up_composition_q95_change=spin_performed ? spin_result.last_composition_change : nothing,
            spin_up_composition_q95_change_all_cells=spin_performed ? spin_result.last_composition_change_all_cells : nothing,
            spin_up_composition_active_cells=spin_performed ? spin_result.composition_active_cells : nothing,
            spin_up_criteria=spin_performed ? join(String.(spin_result.criteria), ";") : "",
            thermal_optima_fraction=OPTIMA_FRACTION,
            upstream_cost=UPSTREAM_COST,
            cedex_var_file=CEDEX_VAR_FILE,
            cedex_uts_file=CEDEX_UTS_FILE,
            obstacles_file=OBSTACLES_FILE,
            obstacle_mode=string(OBSTACLE_MODE),
            obstacle_matched_count=count(data_base.obstacle_mapping_diagnostics.matched),
            obstacle_total_count=nrow(data_base.obstacle_mapping_diagnostics),
            parameter_digest=run_digest)

        # E9/E10: full resolved configuration plus the parameter/code digest, so
        # the run can be audited against the committed config.
        run_metadata = Dict{String,Any}(
            "script" => "run_climate_scenarios.jl",
            "scenario" => scenario,
            "gcm" => gcm,
            "start_year" => START_YEAR,
            "end_year" => END_YEAR,
            "forcing_mode" => FORCING_MODE,
            "heat_stress_k" => k_heat,
            "spin_up" => spin_performed,
            # C3: which projection route produced the initial state, and (for the
            # interim route) how many baseline years were integrated.  Reported
            # so the corrected results are auditable as scenario - control deltas.
            "projection_route" => string(PROJECTION_ROUTE),
            "interim_spin_up_years" => INTERIM_ROUTE ? INTERIM_SPIN_UP_YEARS : 0,
            "carrying_capacity_scaling" => K_SCALING,
            "carrying_capacity_base_scaling" => K_BASE_SCALING,
            "spin_up_criterion" => string(SPIN_UP_CRITERION),
            "spin_up_converged" => spin_performed ? spin_result.converged : nothing,
            "spin_up_years" => spin_performed ? spin_result.years : 0,
            "spin_up_basin_change" => spin_performed ? spin_result.last_basin_change : nothing,
            "spin_up_total_biomass_change" => spin_performed ? spin_result.last_basin_change : nothing,
            "spin_up_composition_q95_change" => spin_performed ? spin_result.last_composition_change : nothing,
            "spin_up_composition_q95_change_all_cells" =>
                spin_performed ? spin_result.last_composition_change_all_cells : nothing,
            "spin_up_composition_active_cells" =>
                spin_performed ? spin_result.composition_active_cells : nothing,
            "spin_up_criteria" => spin_performed ? String.(spin_result.criteria) : String[],
            "spin_up_composition_tol" => SPIN_UP_COMPOSITION_TOL,
            "spin_up_composition_active_floor" => Float64(spin_up_composition_active_floor()),
            "spin_up_min_years" => SPIN_UP_MIN_YEARS,
            "thermal_optima_fraction" => OPTIMA_FRACTION,
            "baseline_period" => "$BASELINE_PERIOD_START-$BASELINE_PERIOD_END",
            "elevation_scaling" => ELEVATION_SCALING,
            "save_interval_days" => SAVE_INTERVAL_DAYS,
            "warming_end_degc" => curve[end],
            "warming_end_degc_source" => "proxy_basin_warming_curve",
            "warming_end_degc_note" =>
                "proxy curve from water_temp_future_2045.csv (6 gauges, window means); " *
                "NOT the applied per-site forcing. Use the realised_* fields instead.",
            "realised_warming_2026_2045_mean_degc" => realised.mean_2026_2045,
            "realised_warming_2036_2045_mean_degc" => realised.mean_2036_2045,
            "realised_warming_n_sites" => realised.n_sites,
            "realised_warming_source" =>
                FORCING_MODE == "daily" ? "daily_schedule" : "annual_mean_schedule",
            "projections_file" => PROJECTION_FILE,
            "obstacle_mode" => string(OBSTACLE_MODE),
            "connectivity_method" => string(CONNECTIVITY_METHOD),
            "parameter_digest" => run_digest,
            "code_version" => code_version(),
            # E13-E16/E18: resolved biological options (the team-chosen values
            # from the corrected config; parameters.toml keeps legacy defaults).
            "biological_options" => SimulationParameters.biological_options_dict(GUADEX_PARAMS),
            "resolved_config" => resolved_config)

        export_result = export_run_outputs(joinpath(run_dir, "export");
            sol_t=sol.t, sol_u=sol.u,
            sites=data_base.sites, species=data_base.species, site_df=data_base.site_df,
            crosswalk_path=joinpath(_GUADEX_ROOT, "data", "site_waterbody_crosswalk.csv"),
            native_species=NATIVE_SPECIES, invasive_species=INVASIVE_SPECIES,
            migratory_species=MIGRATORY_SPECIES,
            days_per_year=DAYS_PER_YEAR, threshold=PRESENCE_THRESHOLD,
            report_year_offsets=YEAR_OFFSETS, report_year_labels=YEAR_LABELS,
            temperature_baseline=data_base.params.temperatures,
            warming=warming,
            dams=data_base.dams, distance_matrix=data_base.distance_matrix,
            habitat_suitability=data_base.params.habitat_suitability,
            upstream_cost=UPSTREAM_COST,
            require_crosswalk=true,
            baseline_species_density=baseline_species_density,
            quasi_extinction_q=QUASI_EXTINCTION_Q,
            quasi_extinction_persistence=QUASI_EXTINCTION_PERSISTENCE,
            species_upper_limits=data_base.params.thermal_upper_limits,
            daily_forcing=forcing_temps === nothing ? nothing :
                (temps=forcing_temps, dates=forcing_dates),
            run_metadata=run_metadata)

        # Write the digest sidecar last, once the JLD2 and full export have
        # succeeded, so an interrupted run cannot be mistaken for a complete one.
        open(digest_path, "w") do io
            println(io, run_digest)
        end

        basin = export_result.level_tables[:basin]
        final_basin = basin[basin.year .== END_YEAR, :]
        duplicate = findall((index_df.scenario .== scenario) .& (index_df.gcm .== gcm))
        isempty(duplicate) || deleteat!(index_df, duplicate)
        push!(index_df, index_row(scenario, gcm, curve[end],
            (realised.mean_2026_2045, realised.mean_2036_2045, realised.n_sites),
            isempty(final_basin) ? NaN : final_basin[1, :mean_native_richness],
            isempty(final_basin) ? NaN : final_basin[1, :mean_realised_richness_loss],
            isempty(final_basin) ? NaN : final_basin[1, :mean_total_biomass],
            run_dir))
        CSV.write(index_path, index_df)
    end
    if MAX_RUNS > 0 && run_number >= MAX_RUNS
        break
    end
end

if nrow(index_df) > 0
    println("\n" * "="^70)
    println("Climate scenarios complete: $(nrow(index_df)) run(s) written under '$base_output_dir'")
    println("Run index: $index_path")
    println("="^70)
else
    println("\nNo climate runs were produced.")
end

if MAKE_FIGURES
    println("\nGenerating climate figures...")
    try
        plot_climate_scenario_figures(base_output_dir;
            figures_dir=joinpath(base_output_dir, "figures"))
    catch err
        @warn "climate figure generation failed" exception=err
    end
    println("\nGenerating climate diagnostic figures...")
    try
        plot_climate_diagnostics(base_output_dir;
            figures_dir=joinpath(base_output_dir, "figures"),
            species_chars_file=joinpath(_GUADEX_ROOT, "data", "ABIOTIC",
                "caracteristicas_peces_Guadalquivir_03-04-2018.csv"),
            native_codes=NATIVE_SPECIES)
    catch err
        @warn "climate diagnostic figure generation failed" exception=err
    end
else
    println("\nSkipping figures (GUADEX_CLIMATE_PLOT=0).")
end
