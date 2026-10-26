# =============================================================================
# E7 — Global sensitivity runner (Latin hypercube + Sobol indices).
#
#   julia --project=. scripts/global_sensitivity_sobol.jl
#
# Reads the `[global_sensitivity]` block of parameters.toml (the agreed-range
# DEFAULT PROPOSAL; the team must confirm the ranges), builds the Saltelli
# A/B/AB/BA design, maps each sample onto the existing `MetacommunityParams`
# transforms (growth / K / interaction / thermal optimum / thermal breadth /
# dispersal / upstream cost / heat-stress k), integrates the ODE, and reports
# first- and total-order Sobol indices with bootstrap confidence intervals.
#
# Modes (env GUADEX_GS_MODE):
#   pilot      (default) reduced configuration: a subset of sites, one year and
#              a uniform warming stand-in; a few minutes.  CLEARLY a PILOT, not
#              the production index.
#   production full site set and the configured 2026-2045 daily warming
#              schedule/site set; N*(2D+2) model runs (expensive).
#
# Overrides: GUADEX_GS_N, GUADEX_GS_BOOTSTRAP, GUADEX_GS_MAX_SITES,
#            GUADEX_GS_SIM_YEARS, GUADEX_GS_SEED, GUADEX_GS_WARMING_C,
#            GUADEX_GS_HEAT_K, GUADEX_GS_OUTPUT_DIR, GUADEX_GS_FORCE_SPINUP.
#
# Nothing here changes the model equations or the run offsets: every knob is an
# existing model parameter, and the within-design ranking in
# scripts/plot_sensitivity_effects.jl is left untouched.
# =============================================================================

using Pkg; Pkg.activate(".")
using DifferentialEquations
using DataFrames
using CSV
using Statistics
using Printf
using Guadex

const _GUADEX_ROOT = isfile(joinpath(@__DIR__, "..", "parameters.jl")) ?
    dirname(@__DIR__) : @__DIR__
include(joinpath(_GUADEX_ROOT, "parameters.jl"))
include(joinpath(_GUADEX_ROOT, "sensitivity_core.jl"))
using .SimulationParameters

# --- Configuration ---------------------------------------------------------

const GUADEX_PARAMS = SimulationParameters.load()
SimulationParameters.require_sections(GUADEX_PARAMS, "general", "inputs", "obstacles",
    "species", "subcatchments", "obstacle_sensitivity", "temperature_stress",
    "global_sensitivity")

const GS = GUADEX_PARAMS["global_sensitivity"]
const SETTINGS = load_sensitivity_settings(GUADEX_PARAMS)

const MODE = lowercase(get(ENV, "GUADEX_GS_MODE", "pilot"))
MODE in ("pilot", "production") || error("GUADEX_GS_MODE must be 'pilot' or 'production'")

_env_int(key, default) = parse(Int, get(ENV, key, string(default)))
_env_float(key, default) = parse(Float64, get(ENV, key, string(default)))

const PILOT = MODE == "pilot"
const N_SAMPLES = _env_int("GUADEX_GS_N", PILOT ? GS["pilot_n_samples"] : GS["n_samples"])
const BOOTSTRAP = _env_int("GUADEX_GS_BOOTSTRAP",
    PILOT ? GS["pilot_bootstrap"] : GS["bootstrap"])
const SEED = _env_int("GUADEX_GS_SEED", GS["seed"])
const MAX_SITES = _env_int("GUADEX_GS_MAX_SITES", PILOT ? GS["pilot_max_sites"] : 0)
const SIM_YEARS = _env_int("GUADEX_GS_SIM_YEARS",
    PILOT ? GS["pilot_sim_years"] : SETTINGS.simulation_years)
const WARMING_C = _env_float("GUADEX_GS_WARMING_C", PILOT ? 2.0 : 0.0)
const OUTPUT_DIR = get(ENV, "GUADEX_GS_OUTPUT_DIR",
    joinpath("results", PILOT ? "global_sensitivity_pilot" : "global_sensitivity"))
const THRESHOLD = Float64(get(GS, "threshold", 0.1))
const METRIC_NAMES = String.(GS["metric_names"])

# --- Parameter block -------------------------------------------------------

function _gs_parameters(spec)
    names = String[]; lower = Float64[]; upper = Float64[]; defaults = Float64[]; rationale = String[]
    for entry in spec
        push!(names, String(entry["name"]))
        push!(lower, Float64(entry["lower"]))
        push!(upper, Float64(entry["upper"]))
        push!(defaults, Float64(entry["default"]))
        push!(rationale, String(get(entry, "rationale", "")))
        lower[end] <= upper[end] || error("global_sensitivity parameter $(names[end]): lower > upper")
    end
    length(names) >= 2 || error("global_sensitivity needs at least 2 parameters")
    return (; names, lower, upper, defaults, rationale)
end

const PARAMS = _gs_parameters(GS["parameters"])
const D = length(PARAMS.names)

# --- Base data and reduced configuration -----------------------------------

const CEDEX_VAR_FILE = get(ENV, "GUADEX_CEDEX_VAR_FILE", GUADEX_PARAMS["inputs"]["cedex_var_file"])
const CEDEX_UTS_FILE = get(ENV, "GUADEX_CEDEX_UTS_FILE", GUADEX_PARAMS["inputs"]["cedex_esc_uts_file"])
const OBSTACLES_FILE = get(ENV, "GUADEX_OBSTACLES_FILE", GUADEX_PARAMS["inputs"]["obstacles_file"])
const OBSTACLE_MODE = Symbol(get(ENV, "GUADEX_OBSTACLE_MODE", GUADEX_PARAMS["obstacles"]["mode"]))
const OBSTACLE_TOL = _env_float("GUADEX_OBSTACLE_TOLERANCE_M", GUADEX_PARAMS["obstacles"]["matching_tolerance_m"])
const OBSTACLE_PASS = _env_float("GUADEX_OBSTACLE_PASSABILITY", GUADEX_PARAMS["obstacles"]["upstream_passability"])
const OBSTACLE_DOWN = _env_float("GUADEX_OBSTACLE_DOWNSTREAM_PASSABILITY", GUADEX_PARAMS["obstacles"]["downstream_passability"])
const CONNECTIVITY_METHOD = SimulationParameters.connectivity_method(GUADEX_PARAMS)

println("="^72)
println("Global sensitivity (E7) — MODE=$(uppercase(MODE))" * (PILOT ? "  [PILOT, not production]" : ""))
println("  parameters:      $(join(PARAMS.names, ", "))")
println("  N (base):        $N_SAMPLES  (bootstrap $BOOTSTRAP, seed $SEED)")
println("  model runs:      $(N_SAMPLES * (2D + 2)) = N*(2D+2)")
println("  horizon:         $(SIM_YEARS) yr, uniform warming stand-in $(WARMING_C) degC")
println("  sites:           $(MAX_SITES == 0 ? "all" : "first $(MAX_SITES)")")
println("  output:          $OUTPUT_DIR")
println("="^72)

println("\nLoading base data (once)...")
data_base = prepare_ode_data(
    upstream_cost = SETTINGS.upstream_cost_default,
    cedex_var_file = CEDEX_VAR_FILE,
    cedex_esc_uts_file = CEDEX_UTS_FILE,
    obstacles_file = OBSTACLES_FILE,
    obstacle_mode = OBSTACLE_MODE,
    obstacle_matching_tolerance = OBSTACLE_TOL,
    obstacle_passability = OBSTACLE_PASS,
    obstacle_downstream_passability = OBSTACLE_DOWN,
    connectivity_method = CONNECTIVITY_METHOD,
    carrying_capacity_base_scaling = SETTINGS.carrying_capacity_base_scaling,
    carrying_capacity_scaling = SETTINGS.carrying_capacity_scaling,
    absence_growth_fraction = SETTINGS.absence_growth_fraction,
    pool_capacity_mode = SETTINGS.pool_capacity_mode,
    fishless_dificil_capacity = SETTINGS.fishless_dificil_capacity,
    nonreproducing_local_growth = SETTINGS.nonreproducing_local_growth,
    exclude_fishfarm_eel_records = SETTINGS.exclude_fishfarm_eel_records,
    salinity_envelope = SETTINGS.salinity_envelope)

NATIVE = String.(GUADEX_PARAMS["species"]["native"])
n_sites_full = data_base.params.n_sites

# Site subset (pilot): spread evenly over model order so the reduced network is
# not one subcatchment.  `MAX_SITES == 0` keeps every site.
site_idx = MAX_SITES == 0 || MAX_SITES >= n_sites_full ? collect(1:n_sites_full) :
    unique(round.(Int, range(1, n_sites_full, length=MAX_SITES)))

function _subset(data_base, idx)
    p = data_base.params
    params = MetacommunityParams(
        length(idx), p.n_species, p.interaction_matrix,
        p.dispersal_matrix[idx, idx], p.dispersal_scaling,
        p.intrinsic_growth_rates[idx, :], p.temperatures[idx],
        p.habitat_suitability[idx], p.thermal_optima, p.thermal_sigmas,
        p.carrying_capacity[idx], p.thermal_lower_limits, p.thermal_upper_limits,
        p.heat_stress_rate, p.interaction_inside_growth)
    return (
        params = params,
        sites = data_base.sites[idx],
        distance_matrix = Matrix(data_base.distance_matrix)[idx, idx],
        dams = Matrix(data_base.dams)[idx, idx],
        elevations = data_base.elevations[idx])
end

sub = _subset(data_base, site_idx)
n_sites = sub.params.n_sites

# Observed initial state restricted to the subset.  Reorder explicitly by site
# code and species code so the flattened state always matches `data_base.sites`.
density_filtered = filter(row -> row.CODIGO in data_base.sites, data_base.density_df)
site_pos = Dict(string(s) => i for (i, s) in enumerate(data_base.sites))
u0_full = zeros(n_sites_full, length(data_base.species))
for r in eachrow(density_filtered)
    pos = get(site_pos, string(r.CODIGO), 0)
    pos == 0 && continue
    for (j, sp) in enumerate(data_base.species)
        v = r[Symbol("$(sp)_DEN")]
        u0_full[pos, j] = (v === missing || isnan(Float64(v))) ? 0.0 : max(Float64(v), 0.0)
    end
end
u0_flat = vec(u0_full[site_idx, :])

# Base heat-stress slope: env override > opt-in calibration (WP3) > configured
# base.  Calibration is opt-in because it reads the (large) daily forcing file.
function _resolve_heat_k()
    haskey(ENV, "GUADEX_GS_HEAT_K") && return _env_float("GUADEX_GS_HEAT_K", 0.0)
    calibrate = lowercase(get(ENV, "GUADEX_GS_CALIBRATE_HEAT", "0")) in ("1", "true", "yes", "on")
    if calibrate && SETTINGS.heat_stress_enabled && SETTINGS.heat_stress_calibrate &&
       isfile(String(SETTINGS.daily_forcing_file))
        _, baseline_temps = load_baseline_forcing(SETTINGS.daily_forcing_file, data_base.sites)
        return resolve_heat_stress_rate(SETTINGS, data_base.params.thermal_upper_limits, baseline_temps)
    end
    return Float64(get(GS, "heat_stress_k_base", SETTINGS.heat_stress_k))
end
const HEAT_K_BASE = _resolve_heat_k()

base_params = sub.params
# Control = central parameter values under BASELINE (unwarmed) climate.  The
# `realised_richness_loss` response is measured against this matched control so
# the metric reflects the warming/parameter perturbation, not the observed
# snapshot's colonisation transient.
control_params = set_heat_stress_rate(base_params, HEAT_K_BASE)
base_params = WARMING_C == 0.0 ? control_params :
    with_temperature_baseline(control_params, control_params.temperatures .+ WARMING_C)

# Production: subset the configured warming schedule to the site set.  The pilot
# uses the static model with the uniform warming stand-in above.
schedule = nothing
if !PILOT
    model = parse_climate_models(get(ENV, "GUADEX_SENSITIVITY_CLIMATE_MODELS",
        GUADEX_PARAMS["obstacle_sensitivity"]["climate_models"]))[end]
    forcing = load_climate_model_forcing(SETTINGS.daily_forcing_file, data_base.sites,
        model.scenario, model.gcm;
        baseline_start=SETTINGS.baseline_period_start, baseline_end=SETTINGS.baseline_period_end,
        baseline_means=data_base.params.temperatures, days_per_year=SETTINGS.days_per_year,
        year_labels=SETTINGS.year_labels)
    schedule = TemperatureSchedule(forcing.schedule.deltas[site_idx, :], SETTINGS.days_per_year)
    println("  production warming: $(model.scenario)/$(model.gcm) over $(SETTINGS.start_year)-$(SETTINGS.end_year)")
end

println("  heat-stress k base: $(HEAT_K_BASE)")
println("  reduced network:   $(n_sites) sites × $(sub.params.n_species) species")

# --- Evaluator -------------------------------------------------------------

const T_END = Float64(SIM_YEARS * SETTINGS.days_per_year)
const POSITIVITY_CB = sensitivity_positivity_callback()

# Integrate one parameter set over the horizon and return the final state.
function _integrate(p, u0, sched)
    rhs_params = sched === nothing ? p : ScheduledMetacommunityParams(p, sched)
    rhs! = sched === nothing ? metacommunity_ode! : metacommunity_ode_scheduled!
    prob = ODEProblem(rhs!, u0, (0.0, T_END), rhs_params)
    sol = solve(prob, Tsit5(); reltol=1e-6, abstol=1e-6,
        save_everystep=false, save_start=true, save_end=true, callback=POSITIVITY_CB)
    return vec(Float64.(sol.u[end]))
end

# Matched control: central parameter values under baseline (unwarmed) climate.
const BASELINE_STATE = _integrate(control_params, u0_flat, nothing)

function evaluate(theta::AbstractVector)
    theta = Float64.(theta)
    p = base_params
    p = scale_intrinsic_growth_rate(p, theta[1])                 # growth multiplier
    p = scale_carrying_capacity(p, theta[2])                     # K scaling
    p = set_interaction_matrix(p, p.interaction_matrix .* theta[3])   # interaction strength
    p = set_thermal_optima(p, p.thermal_optima .+ theta[4])      # optimum shift (degC)
    p = set_thermal_sigma_multiplier(p, theta[5])                # thermal breadth
    p = set_heat_stress_rate(p, HEAT_K_BASE * theta[8])          # heat-stress k
    dispersal = precompute_dispersal_matrix(p.n_sites, sub.distance_matrix,
        sub.elevations, theta[7], sub.dams, data_base.species)   # upstream cost
    p = set_dispersal_matrix(p, dispersal)
    p = scale_dispersal(p, theta[6])                             # dispersal rate

    final_state = _integrate(p, u0_flat, schedule)
    metrics = basin_response_metrics(final_state, BASELINE_STATE, data_base.species;
        native=NATIVE, threshold=THRESHOLD)
    return [metrics.native_biomass, metrics.native_richness, metrics.realised_richness_loss]
end

# --- Design, run and report ------------------------------------------------

design = saltelli_design(N_SAMPLES, D; seed=SEED)
println("\nSaltelli design: n=$(N_SAMPLES), d=$D, model runs=$(n_model_runs(design))")
println("Evaluating $(n_model_runs(design)) model runs...")

started = time()
result = sobol_analyze(evaluate, design;
    lower=PARAMS.lower, upper=PARAMS.upper,
    bootstrap=BOOTSTRAP, seed=SEED,
    parameter_names=PARAMS.names, metric_names=METRIC_NAMES)
elapsed = round(time() - started, digits=1)

# --- Outputs ---------------------------------------------------------------

mkpath(OUTPUT_DIR)

range_df = DataFrame(
    name = PARAMS.names, lower = PARAMS.lower, upper = PARAMS.upper,
    default = PARAMS.defaults, rationale = PARAMS.rationale)
CSV.write(joinpath(OUTPUT_DIR, "parameter_ranges.csv"), range_df)

rows = NamedTuple[]
for (mi, metric) in enumerate(METRIC_NAMES)
    for (pi, name) in enumerate(PARAMS.names)
        push!(rows, (
            metric = metric, parameter = name,
            S1 = result.S1[pi, mi], S1_lo = result.S1_ci[pi, mi, 1], S1_hi = result.S1_ci[pi, mi, 2],
            ST = result.ST[pi, mi], ST_lo = result.ST_ci[pi, mi, 1], ST_hi = result.ST_ci[pi, mi, 2]))
    end
end
CSV.write(joinpath(OUTPUT_DIR, "sobol_indices.csv"), DataFrame(rows))

samples = saltelli_sample_matrix(design; lower=PARAMS.lower, upper=PARAMS.upper)
samples_df = DataFrame(samples, PARAMS.names)
samples_df.role = vcat(fill("A", N_SAMPLES), fill("B", N_SAMPLES),
    [fill("AB_$i", N_SAMPLES) for i in 1:D]..., [fill("BA_$i", N_SAMPLES) for i in 1:D]...)
CSV.write(joinpath(OUTPUT_DIR, "saltelli_samples.csv"), samples_df)

open(joinpath(OUTPUT_DIR, "run_summary.txt"), "w") do io
    println(io, "GuadeX E7 global sensitivity — mode=$(uppercase(MODE))")
    println(io, "pilot=$(PILOT) (pilot results are NOT the production index)")
    println(io, "N=$N_SAMPLES d=$D model_runs=$(n_model_runs(design)) bootstrap=$BOOTSTRAP seed=$SEED")
    println(io, "sites=$(n_sites) sim_years=$SIM_YEARS warming_c=$WARMING_C heat_stress_k_base=$HEAT_K_BASE")
    println(io, "elapsed_s=$elapsed")
    println(io, "parameters=$(join(PARAMS.names, ","))")
    println(io, "metrics=$(join(METRIC_NAMES, ","))")
end

println("\n" * "="^72)
println("Global sensitivity complete in $(elapsed)s")
println("  indices:  $(joinpath(OUTPUT_DIR, "sobol_indices.csv"))")
println("  ranges:   $(joinpath(OUTPUT_DIR, "parameter_ranges.csv"))")
println("  samples:  $(joinpath(OUTPUT_DIR, "saltelli_samples.csv"))")
println("="^72)

for (mi, metric) in enumerate(METRIC_NAMES)
    println("\n$(metric)")
    @printf("  %-30s %8s %8s\n", "parameter", "S1", "ST")
    for (pi, name) in enumerate(PARAMS.names)
        @printf("  %-30s %8.3f %8.3f\n", name, result.S1[pi, mi], result.ST[pi, mi])
    end
end
