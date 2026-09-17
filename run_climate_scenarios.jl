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
const SPIN_UP_CRITERION = Symbol(lowercase(get(ENV, "GUADEX_CLIMATE_SPIN_UP_CRITERION",
    string(_climate_setting("spin_up_criterion", "basin")))))
const SPIN_UP_PROGRESS_EVERY = parse(Int, get(ENV, "GUADEX_CLIMATE_SPIN_UP_PROGRESS",
    string(_climate_setting("spin_up_progress_every", 50))))
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
const YEAR_OFFSETS = collect(0:(SIMULATION_YEARS - 1))
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
    println("  spin-up: $SPIN_UP_ENABLED (max $(SPIN_UP_MAX_YEARS) yr, criterion=:$SPIN_UP_CRITERION)")
    println("  carrying-capacity base scaling: $(K_BASE_SCALING)x; WP4 scaling: $(K_SCALING)x; " *
            "effective: $(K_BASE_SCALING * K_SCALING)x")
    println("  optimum fraction: $OPTIMA_FRACTION")
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
    carrying_capacity_base_scaling = K_BASE_SCALING
)

n_sites = data_base.params.n_sites

# Initial conditions from the observed density matrix.
density_cols = [Symbol("$(sp)_DEN") for sp in data_base.species]
density_df_filtered = filter(row -> row.CODIGO in data_base.sites, data_base.density_df)
u0 = Matrix(density_df_filtered[:, density_cols])
replace!(u0, NaN => 0.0)
u0 = max.(u0, 0.0)
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

# WP4: spin up to the baseline equilibrium once and reuse for every scenario.
baseline_species_density = nothing
if SPIN_UP_ENABLED
    spin_schedule = nothing
    if FORCING_MODE == "daily" && baseline_temps !== nothing
        # Spin up under the seasonal baseline climatology, not the annual mean:
        # the two equilibria differ and only the seasonal one matches the runs.
        spin_schedule = baseline_climatology_schedule(baseline_temps, baseline_dates, 1;
            days_per_year=Float64(DAYS_PER_YEAR))
    end
    println("WP4 spin-up under $(spin_schedule === nothing ? "annual-mean" : "seasonal") " *
            "baseline forcing (max $(SPIN_UP_MAX_YEARS) yr, tol=$SPIN_UP_TOL, " *
            "criterion=:$SPIN_UP_CRITERION)...")
    spin = spin_up(params; initial_state=u0_flat, schedule=spin_schedule,
        days_per_year=Float64(DAYS_PER_YEAR), max_years=SPIN_UP_MAX_YEARS, tol=SPIN_UP_TOL,
        criterion=SPIN_UP_CRITERION, progress_every=SPIN_UP_PROGRESS_EVERY)
    u0_flat = spin.state
    baseline_species_density = reshape(u0_flat, n_sites, params.n_species)
    println("  spin-up: years=$(spin.years), converged=$(spin.converged), " *
            "basin change=$(spin.last_basin_change), q95=$(spin.last_q95_change), " *
            "max=$(spin.last_relative_change)")
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
    return (r.mean_native_richness, r.mean_native_extinction_risk, r.mean_total_biomass)
end

const CLIMATE_OUTPUT_DIR = get(ENV, "GUADEX_CLIMATE_OUTPUT_DIR",
    string(_climate_setting("output_dir", joinpath("results", "climate_scenarios"))))
base_output_dir = CLIMATE_OUTPUT_DIR
mkpath(base_output_dir)
index_path = joinpath(base_output_dir, "runs_index.csv")
index_df = isfile(index_path) ? CSV.read(index_path, DataFrame) : DataFrame(
    scenario=String[], gcm=String[], start_year=Int[], end_year=Int[],
    warming_end_degc=Float64[], basin_native_richness=Float64[],
    basin_native_extinction_risk=Float64[], basin_total_biomass=Float64[],
    run_dir=String[])

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

        # Resume support: a completed run is skipped unless explicitly forced,
        # and its summary row is preserved in the incremental run index.
        if !FORCE_RERUN && isfile(joinpath(run_dir, "simulation_output.jld2"))
            println("\n[$run_number/$total_runs] $scenario / $gcm -> skipped (already exists)")
            if !any((index_df.scenario .== scenario) .& (index_df.gcm .== gcm))
                richness, risk, biomass = basin_row_from_csv(
                    joinpath(run_dir, "export", "levels", "level_basin.csv"), END_YEAR)
                push!(index_df, (scenario=scenario, gcm=gcm, start_year=START_YEAR,
                    end_year=END_YEAR, warming_end_degc=curve[end],
                    basin_native_richness=richness, basin_native_extinction_risk=risk,
                    basin_total_biomass=biomass, run_dir=run_dir))
                CSV.write(index_path, index_df)
            end
            continue
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
            # Prefer the baseline from the same source so per-GCM runs difference
            # against their own historical period (removing GCM bias).
            if "historical" in Set(String.(src.scenario))
                bdates, btemps = forcing_matrix(src, data_base.sites; scenario="historical")
            else
                bdates, btemps = baseline_dates, baseline_temps
            end
            sched, _ = daily_temperature_schedule(forcing_temps;
                dates=forcing_dates,
                baseline_temps=btemps, baseline_dates=bdates,
                baseline_start=BASELINE_PERIOD_START, baseline_end=BASELINE_PERIOD_END)
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
        scheduled_params = ScheduledMetacommunityParams(data_base.params, schedule)

        prob = ODEProblem(metacommunity_ode_scheduled!, u0_flat, (0.0, T_END), scheduled_params)
        sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, saveat=saveat, callback=positivity_cb)
        println("  solved: $(length(sol.t)) snapshots, final basin delta = $(round(curve[end], digits=3)) degC")

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
            spin_up=SPIN_UP_ENABLED,
            carrying_capacity_scaling=K_SCALING,
            carrying_capacity_base_scaling=K_BASE_SCALING,
            spin_up_criterion=string(SPIN_UP_CRITERION),
            thermal_optima_fraction=OPTIMA_FRACTION,
            upstream_cost=UPSTREAM_COST,
            cedex_var_file=CEDEX_VAR_FILE,
            cedex_uts_file=CEDEX_UTS_FILE,
            obstacles_file=OBSTACLES_FILE,
            obstacle_mode=string(OBSTACLE_MODE),
            obstacle_matched_count=count(data_base.obstacle_mapping_diagnostics.matched),
            obstacle_total_count=nrow(data_base.obstacle_mapping_diagnostics))

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
            run_metadata=Dict(
                "script" => "run_climate_scenarios.jl",
                "scenario" => scenario,
                "gcm" => gcm,
                "start_year" => START_YEAR,
                "end_year" => END_YEAR,
                "forcing_mode" => FORCING_MODE,
                "heat_stress_k" => k_heat,
                "spin_up" => SPIN_UP_ENABLED,
                "carrying_capacity_scaling" => K_SCALING,
                "carrying_capacity_base_scaling" => K_BASE_SCALING,
                "spin_up_criterion" => string(SPIN_UP_CRITERION),
                "spin_up_converged" => SPIN_UP_ENABLED ? spin.converged : nothing,
                "spin_up_years" => SPIN_UP_ENABLED ? spin.years : nothing,
                "spin_up_basin_change" => SPIN_UP_ENABLED ? spin.last_basin_change : nothing,
                "thermal_optima_fraction" => OPTIMA_FRACTION,
                "baseline_period" => "$BASELINE_PERIOD_START-$BASELINE_PERIOD_END",
                "elevation_scaling" => ELEVATION_SCALING,
                "save_interval_days" => SAVE_INTERVAL_DAYS,
                "warming_end_degc" => curve[end],
                "projections_file" => PROJECTION_FILE,
                "obstacle_mode" => string(OBSTACLE_MODE)))

        basin = export_result.level_tables[:basin]
        final_basin = basin[basin.year .== END_YEAR, :]
        duplicate = findall((index_df.scenario .== scenario) .& (index_df.gcm .== gcm))
        isempty(duplicate) || deleteat!(index_df, duplicate)
        push!(index_df, (
            scenario=scenario, gcm=gcm,
            start_year=START_YEAR, end_year=END_YEAR,
            warming_end_degc=curve[end],
            basin_native_richness=isempty(final_basin) ? NaN : final_basin[1, :mean_native_richness],
            basin_native_extinction_risk=isempty(final_basin) ? NaN : final_basin[1, :mean_native_extinction_risk],
            basin_total_biomass=isempty(final_basin) ? NaN : final_basin[1, :mean_total_biomass],
            run_dir=run_dir))
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
