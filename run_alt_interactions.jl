using Pkg; Pkg.activate(".")
using DifferentialEquations
using DataFrames
using CSV
using SparseArrays
using LinearAlgebra
using JLD2
using Dates
using Random
using Guadex

# =============================================================================
# --- Alternative-interaction obstacle sensitivity ---
#
# Extends `run_sensitivity_report.jl` with a robustness axis: the interaction
# matrix (empirical, random, invasive-favouring) and a thermal-sigma multiplier,
# on top of the same upstream-cost x passability sweep.  The climate background,
# daily forcing, seasonal burn-in, heat stress, K convention and species
# classification are shared with the climate scenarios.
#
# Settings come from `[obstacle_sensitivity]` and `[run_alt_interactions]` in
# parameters.toml.
# =============================================================================

const _GUADEX_ROOT = isfile(joinpath(@__DIR__, "parameters.jl")) ? (@__DIR__) : dirname(@__DIR__)
include(joinpath(_GUADEX_ROOT, "parameters.jl"))
include(joinpath(_GUADEX_ROOT, "sensitivity_core.jl"))
using .SimulationParameters
const GUADEX_PARAMS = SimulationParameters.load()
SimulationParameters.require_sections(GUADEX_PARAMS, "general", "inputs", "obstacles",
    "species", "subcatchments", "scenarios", "obstacle_sensitivity", "run_alt_interactions")

const SETTINGS = load_sensitivity_settings(GUADEX_PARAMS)
const THERMAL_SIGMA_MULTIPLIER = Float64(GUADEX_PARAMS["run_alt_interactions"]["thermal_sigma_multiplier"])

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

all_subcatchments = SimulationParameters.float_vector(GUADEX_PARAMS["subcatchments"]["all"])

const UPSTREAM_COSTS = SimulationParameters.float_vector(GUADEX_PARAMS["run_alt_interactions"]["upstream_costs"])
const PASSABILITY_NAMES = String.(GUADEX_PARAMS["run_alt_interactions"]["passability_scenarios"])
const PASSABILITY_SCENARIOS = SimulationParameters.scenario_library(
    all_subcatchments, GUADEX_PARAMS["scenarios"]["passability"], PASSABILITY_NAMES)
const CLIMATE_MODELS = parse_climate_models(get(ENV, "GUADEX_SENSITIVITY_CLIMATE_MODELS",
    GUADEX_PARAMS["obstacle_sensitivity"]["climate_models"]))
const MAX_RUNS = parse(Int, get(ENV, "GUADEX_SENSITIVITY_MAX_RUNS", "0"))
const FORCE_RERUN = lowercase(get(ENV, "GUADEX_SENSITIVITY_FORCE", "0")) in ("1", "true", "yes")

if get(ENV, "GUADEX_CONFIG_ONLY", "0") == "1"
    println("[parameters] configuration smoke test passed")
    println("  horizon:         $(SETTINGS.start_year)-$(SETTINGS.end_year) ($(SETTINGS.simulation_years) yr, daily)")
    println("  climate models:  $(join([climate_model_tag(m) for m in CLIMATE_MODELS], ", "))")
    println("  upstream costs:  $(join(UPSTREAM_COSTS, ", "))")
    println("  passability:     $(join(PASSABILITY_NAMES, ", "))")
    println("  sigma multiplier: $(THERMAL_SIGMA_MULTIPLIER)")
    println("  forcing file:    $(SETTINGS.daily_forcing_file)")
    println("  heat stress:     enabled=$(SETTINGS.heat_stress_enabled) calibrate=$(SETTINGS.heat_stress_calibrate)")
    println("  burn-in:         max $(SETTINGS.spin_up_max_years) yr, tol=$(SETTINGS.spin_up_tol), criterion=:$(SETTINGS.spin_up_criterion)")
    exit(0)
end

# =============================================================================
# --- Alternative interaction matrices ---
# =============================================================================

"""
    make_random_interaction_matrix(n_species, val_range)

Every off-diagonal cell gets a random negative value in `[-val_range, 0]`;
fixed seed for reproducibility.
"""
function make_random_interaction_matrix(n_species::Int, val_range::Float64=1.0)
    Random.seed!(42)
    mat = zeros(n_species, n_species)
    for i in 1:n_species, j in 1:n_species
        i == j && continue
        mat[i, j] = -rand() * val_range
    end
    return mat
end

"""
    make_invasive_favoring_matrix(n_species, native_idx, invasive_idx)

Invasives suppress natives (-0.8), natives have no effect on invasives (0.0),
invasives compete moderately (-0.3) and natives weakly (-0.1); unclassified
species are neutral.
"""
function make_invasive_favoring_matrix(n_species::Int, native_idx::Vector{Int}, invasive_idx::Vector{Int})
    mat = zeros(n_species, n_species)
    for inv in invasive_idx, nat in native_idx
        mat[nat, inv] = -0.8
    end
    for inv in invasive_idx, inv2 in invasive_idx
        inv == inv2 && continue
        mat[inv2, inv] = -0.3
    end
    for nat in native_idx, nat2 in native_idx
        nat == nat2 && continue
        mat[nat2, nat] = -0.1
    end
    return mat
end

println("="^70)
println("Alternative-interaction obstacle sensitivity (daily forcing + burn-in)")
println("  horizon:        $(SETTINGS.start_year)-$(SETTINGS.end_year)")
println("  climate models: $(join([climate_model_tag(m) for m in CLIMATE_MODELS], ", "))")
println("  sigma multiplier: $(THERMAL_SIGMA_MULTIPLIER)")
println("="^70)

# =============================================================================
# --- Base data, burn-in per interaction matrix ---
# =============================================================================

println("\nLoading base data (once)...")
data_base = prepare_ode_data(
    upstream_cost = SETTINGS.upstream_cost_default,
    cedex_var_file = CEDEX_VAR_FILE,
    cedex_esc_uts_file = CEDEX_UTS_FILE,
    obstacles_file = OBSTACLES_FILE,
    obstacle_mode = OBSTACLE_MODE,
    obstacle_matching_tolerance = OBSTACLE_MATCHING_TOLERANCE,
    obstacle_passability = OBSTACLE_PASSABILITY,
    obstacle_downstream_passability = OBSTACLE_DOWNSTREAM_PASSABILITY,
    carrying_capacity_base_scaling = SETTINGS.carrying_capacity_base_scaling,
    carrying_capacity_scaling = SETTINGS.carrying_capacity_scaling
)

n_sites = data_base.params.n_sites
n_species = data_base.params.n_species

base_params = data_base.params
if SETTINGS.thermal_optima_fraction != 0.5
    optima = optimum_sweep_optima(data_base.species_chars_df, data_base.species;
        fraction=SETTINGS.thermal_optima_fraction)
    base_params = set_thermal_optima(base_params, optima)
end

density_cols = [Symbol("$(sp)_DEN") for sp in data_base.species]
density_df_filtered = filter(row -> row.CODIGO in data_base.sites, data_base.density_df)
u0_obs = Matrix(density_df_filtered[:, density_cols])
replace!(u0_obs, NaN => 0.0)
u0_obs = max.(u0_obs, 0.0)

println("\nLoading baseline daily forcing ($(SETTINGS.daily_forcing_file))...")
baseline_dates, baseline_temps = load_baseline_forcing(SETTINGS.daily_forcing_file, data_base.sites)
println("  baseline series: $(length(baseline_dates)) days " *
        "($(SETTINGS.baseline_period_start)-$(SETTINGS.baseline_period_end))")

k_heat = resolve_heat_stress_rate(SETTINGS, base_params.thermal_upper_limits, baseline_temps)
println("Heat-stress mortality: $(k_heat > 0 ? "k=$k_heat" : "disabled")")

native_idx = classify_species_indices(data_base.species, NATIVE_SPECIES)
invasive_idx = classify_species_indices(data_base.species, INVASIVE_SPECIES)

original_interaction = base_params.interaction_matrix
matrices = [
    (name="original", matrix=original_interaction),
    (name="random", matrix=make_random_interaction_matrix(n_species, abs(minimum(original_interaction)))),
    (name="invasive_favoring", matrix=make_invasive_favoring_matrix(n_species, native_idx, invasive_idx)),
]

# One burn-in per interaction matrix + sigma: the equilibrium composition
# depends on the interaction structure, so every matrix is started from its own
# baseline equilibrium.
matrix_states = Dict{String,Any}()
for entry in matrices
    p = set_heat_stress_rate(set_interaction_matrix(base_params, entry.matrix), k_heat)
    if THERMAL_SIGMA_MULTIPLIER != 1.0
        p = set_thermal_sigma_multiplier(p, THERMAL_SIGMA_MULTIPLIER)
    end
    tag = "alt_$(entry.name)_sig$(THERMAL_SIGMA_MULTIPLIER)"
    spinup_extra = "mode=$(OBSTACLE_MODE)_tol=$(OBSTACLE_MATCHING_TOLERANCE)_" *
        "pass=$(OBSTACLE_PASSABILITY)_down=$(OBSTACLE_DOWNSTREAM_PASSABILITY)_" *
        "obs=$(basename(OBSTACLES_FILE))_matrix=$(entry.name)_sig=$(THERMAL_SIGMA_MULTIPLIER)"
    spin = if SETTINGS.spin_up
        get_or_compute_spinup(SETTINGS, p, vec(u0_obs), baseline_temps, baseline_dates;
            tag=tag, extra=spinup_extra)
    else
        println("Burn-in disabled ([obstacle_sensitivity].spin_up = false) for " *
                "matrix $(entry.name); starting from observed densities")
        observed_initial_state(u0_obs)
    end
    matrix_states[entry.name] = (params=p, spin=spin,
        baseline_species_density=reshape(spin.state, n_sites, n_species))
end

# The forcing depends only on (scenario, gcm), so load each per-GCM series once
# and reuse it across all interaction matrices.
forcing_by_model = Dict{String,Any}()
warming_end_by_model = Dict{String,Float64}()
for model in CLIMATE_MODELS
    model_tag = climate_model_tag(model)
    forcing = load_climate_model_forcing(SETTINGS.daily_forcing_file, data_base.sites,
        model.scenario, model.gcm;
        baseline_start=SETTINGS.baseline_period_start,
        baseline_end=SETTINGS.baseline_period_end,
        days_per_year=SETTINGS.days_per_year,
        year_labels=SETTINGS.year_labels)
    forcing_by_model[model_tag] = forcing
    warming_end_by_model[model_tag] = maximum(forcing.warming[:, end])
    println("  loaded daily forcing for $model_tag (" *
            "$(round(warming_end_by_model[model_tag], digits=3)) degC by $(SETTINGS.end_year))")
end

saveat = sensitivity_saveat(SETTINGS)
positivity_cb = sensitivity_positivity_callback()

# =============================================================================
# --- Main sweep ---
# =============================================================================

base_output_dir = get(ENV, "GUADEX_ALT_OUTPUT_DIR",
    joinpath(SETTINGS.output_dir, "alt_interactions"))
mkpath(base_output_dir)
index_path = joinpath(base_output_dir, "runs_index.csv")
index_df = read_sensitivity_index(index_path; with_matrix=true)

total_runs = length(matrices) * length(CLIMATE_MODELS) * length(UPSTREAM_COSTS) * length(PASSABILITY_NAMES)
current_run = 0
started = time()

for entry in matrices
    matrix_name = entry.name
    state = matrix_states[matrix_name]
    params = state.params
    spin = state.spin
    baseline_species_density = state.baseline_species_density

    println("\n" * "="^70)
    println("Interaction matrix: $matrix_name (sigma x $(THERMAL_SIGMA_MULTIPLIER))")
    println("="^70)

    for model in CLIMATE_MODELS
        model_tag = climate_model_tag(model)
        println("\n" * "-"^70)
        println("Climate model: $(model.scenario) / $(model.gcm)  [$matrix_name]")
        println("-"^70)

        forcing = forcing_by_model[model_tag]
        warming_end = warming_end_by_model[model_tag]
        println("  warming by $(SETTINGS.end_year): $(round(warming_end, digits=3)) degC")

        for uc in UPSTREAM_COSTS
            for pass_name in PASSABILITY_NAMES
                global current_run += 1
                if MAX_RUNS > 0 && current_run > MAX_RUNS
                    println("[cap] GUADEX_SENSITIVITY_MAX_RUNS=$MAX_RUNS reached; stopping")
                    current_run = MAX_RUNS
                    break
                end

                run_dir = joinpath(base_output_dir, matrix_name,
                    "sig_$(THERMAL_SIGMA_MULTIPLIER)", model_tag, "uc_$(uc)", pass_name)
                println("\n[$current_run/$total_runs] $matrix_name $model_tag uc=$(uc) pass=$(pass_name)")

                fingerprint = sensitivity_run_fingerprint(SETTINGS, model;
                    upstream_cost=uc, passability_scenario=pass_name,
                    interaction_matrix=matrix_name, sigma=THERMAL_SIGMA_MULTIPLIER)

                if !FORCE_RERUN && sensitivity_run_status(run_dir, fingerprint)
                    println("  skipped (complete and fingerprint matches; " *
                            "GUADEX_SENSITIVITY_FORCE=1 to rerun)")
                else
                    mkpath(run_dir)
                    passability_vector = build_dam_passability_vector(data_base.site_df, data_base.sites;
                        passability_per_subcatchment=PASSABILITY_SCENARIOS[pass_name])
                    sol, modified_dams, scenario_params = run_sensitivity_scenario(
                        data_base, params, forcing.schedule, spin.state, uc, passability_vector;
                        settings=SETTINGS, saveat=saveat, positivity_cb=positivity_cb)

                    jldsave(joinpath(run_dir, "simulation_output.jld2");
                        sol_t=sol.t, sol_u=sol.u,
                        sites=data_base.sites, species=data_base.species,
                        climate_scenario=model.scenario, gcm=model.gcm,
                        interaction_matrix_type=matrix_name,
                        thermal_sigma_multiplier=THERMAL_SIGMA_MULTIPLIER,
                        start_year=SETTINGS.start_year, end_year=SETTINGS.end_year,
                        upstream_cost=uc, passability_scenario=pass_name,
                        passability_vector=passability_vector,
                        simulation_years=SETTINGS.simulation_years,
                        warming_end_degc=warming_end,
                        temperature_baseline=params.temperatures, warming=forcing.warming,
                        heat_stress_k=k_heat,
                        spin_up=SETTINGS.spin_up, spin_up_years=spin.years, spin_up_converged=spin.converged,
                        carrying_capacity_base_scaling=SETTINGS.carrying_capacity_base_scaling,
                        carrying_capacity_scaling=SETTINGS.carrying_capacity_scaling,
                        thermal_optima_fraction=SETTINGS.thermal_optima_fraction,
                        cedex_var_file=CEDEX_VAR_FILE, cedex_uts_file=CEDEX_UTS_FILE,
                        obstacles_file=OBSTACLES_FILE, obstacle_mode=string(OBSTACLE_MODE),
                        obstacle_matching_tolerance_m=OBSTACLE_MATCHING_TOLERANCE,
                        obstacle_passability=OBSTACLE_PASSABILITY,
                        obstacle_downstream_passability=OBSTACLE_DOWNSTREAM_PASSABILITY,
                        obstacle_matched_count=count(data_base.obstacle_mapping_diagnostics.matched),
                        obstacle_total_count=nrow(data_base.obstacle_mapping_diagnostics))

                    export_sensitivity_run(run_dir;
                        data_base=data_base, sol=sol, params=scenario_params, forcing=forcing,
                        settings=SETTINGS, dams=modified_dams, upstream_cost=uc,
                        baseline_species_density=baseline_species_density,
                        native_species=NATIVE_SPECIES, invasive_species=INVASIVE_SPECIES,
                        migratory_species=MIGRATORY_SPECIES,
                        run_metadata=sensitivity_run_metadata(SETTINGS, model;
                            script="run_alt_interactions.jl",
                            upstream_cost=uc, passability_scenario=pass_name,
                            warming_end=warming_end, heat_stress_k=k_heat,
                            spin_years=spin.years, spin_converged=spin.converged,
                            obstacle_mode=OBSTACLE_MODE,
                            extras=Dict{String,Any}(
                                "interaction_matrix_type" => matrix_name,
                                "thermal_sigma_multiplier" => THERMAL_SIGMA_MULTIPLIER)))
                    mark_sensitivity_run_complete(run_dir, fingerprint)
                    println("  saved to: $run_dir")
                end

                row = index_row_values(run_dir, SETTINGS, model;
                    upstream_cost=uc, passability_scenario=pass_name, warming_end=warming_end,
                    matrix=matrix_name, sigma=THERMAL_SIGMA_MULTIPLIER)
                upsert_index_row!(index_df, row)
                CSV.write(index_path, index_df)
            end
            MAX_RUNS > 0 && current_run >= MAX_RUNS && break
        end
        MAX_RUNS > 0 && current_run >= MAX_RUNS && break
    end
    MAX_RUNS > 0 && current_run >= MAX_RUNS && break
end

elapsed = round(time() - started, digits=1)
println("\n" * "="^70)
println("Alternative-interaction sensitivity complete: $(nrow(index_df)) run(s) under '$base_output_dir' in $(elapsed)s")
println("Run index: $index_path")
println("="^70)
