using Pkg; Pkg.activate(".");
using DifferentialEquations
using DataFrames
using CSV
using SparseArrays
using LinearAlgebra
using JLD2
using Guadex

# =============================================================================
# --- Sensitivity Analysis: Temperature × Upstream Cost × Passability ---
# Runs simulations with different combinations of:
#   1. Temperature increase (0°C to 3°C gradient)
#   2. Upstream migration cost (upstream cost parameter)
#   3. Passability scenarios (per-subcatchment dam passability multipliers)
# No exploitation applied. 3-year simulations.
# Results are saved in organized folder structure.
# =============================================================================

# --- Load run parameters (single source of truth: parameters.toml) ---
const _GUADEX_ROOT = isfile(joinpath(@__DIR__, "parameters.jl")) ? (@__DIR__) : dirname(@__DIR__)
include(joinpath(_GUADEX_ROOT, "parameters.jl"))
using .SimulationParameters
const GUADEX_PARAMS = SimulationParameters.load()
SimulationParameters.require_sections(GUADEX_PARAMS, "general", "inputs", "obstacles", "species", "subcatchments", "scenarios", "run_sensitivity_report")

const DAYS_PER_YEAR = Int(GUADEX_PARAMS["general"]["days_per_year"])
const SIMULATION_YEARS = Int(GUADEX_PARAMS["run_sensitivity_report"]["simulation_years"])
const T_SPAN = (0.0, Float64(SIMULATION_YEARS * DAYS_PER_YEAR))

# Updated 2045 inputs.  CEDEX tables are loaded for provenance; the obstacle
# overlay is the configured network update used by this simulation.  Each value
# can be overridden per run via its GUADEX_* environment variable.
const CEDEX_VAR_FILE = get(ENV, "GUADEX_CEDEX_VAR_FILE", GUADEX_PARAMS["inputs"]["cedex_var_file"])
const CEDEX_UTS_FILE = get(ENV, "GUADEX_CEDEX_UTS_FILE", GUADEX_PARAMS["inputs"]["cedex_esc_uts_file"])
const OBSTACLES_FILE = get(ENV, "GUADEX_OBSTACLES_FILE", GUADEX_PARAMS["inputs"]["obstacles_file"])
const OBSTACLE_MODE = Symbol(get(ENV, "GUADEX_OBSTACLE_MODE", GUADEX_PARAMS["obstacles"]["mode"]))
const OBSTACLE_MATCHING_TOLERANCE = parse(Float64, get(ENV, "GUADEX_OBSTACLE_TOLERANCE_M", string(GUADEX_PARAMS["obstacles"]["matching_tolerance_m"])))
const OBSTACLE_PASSABILITY = parse(Float64, get(ENV, "GUADEX_OBSTACLE_PASSABILITY", string(GUADEX_PARAMS["obstacles"]["upstream_passability"])))
const OBSTACLE_DOWNSTREAM_PASSABILITY = parse(Float64, get(ENV, "GUADEX_OBSTACLE_DOWNSTREAM_PASSABILITY", string(GUADEX_PARAMS["obstacles"]["downstream_passability"])))

const NATIVE_SPECIES = String.(GUADEX_PARAMS["species"]["native"])
const INVASIVE_SPECIES = String.(GUADEX_PARAMS["species"]["invasive"])

all_subcatchments = SimulationParameters.float_vector(GUADEX_PARAMS["subcatchments"]["all"])

# =============================================================================
# --- Define Parameter Grid ---
# Grid values and scenario names are read from the [run_sensitivity_report]
# section of parameters.toml.
# =============================================================================

temperature_increases = SimulationParameters.float_vector(GUADEX_PARAMS["run_sensitivity_report"]["temperature_increases"])

upstream_costs = SimulationParameters.float_vector(GUADEX_PARAMS["run_sensitivity_report"]["upstream_costs"])

passability_scenarios = SimulationParameters.scenario_library(
    all_subcatchments, GUADEX_PARAMS["scenarios"]["passability"], GUADEX_PARAMS["run_sensitivity_report"]["passability_scenarios"])

report_years = Int.(SimulationParameters.float_vector(GUADEX_PARAMS["run_sensitivity_report"]["report_years"]))

# Stop before data loading when GUADEX_CONFIG_ONLY=1 (configuration smoke test).
if get(ENV, "GUADEX_CONFIG_ONLY", "0") == "1"
    println("[parameters] configuration smoke test passed")
    exit(0)
end

# =============================================================================
# --- Helper Functions ---
# =============================================================================

function build_dam_passability_vector(site_df, sites; passability_per_subcatchment::Dict{Float64, Float64}=Dict{Float64, Float64}())
    n_sites = length(sites)
    passability_vector = ones(n_sites)
    site_to_sc = Dict(row.CODIGO => row.CODIGO_S for row in eachrow(site_df))
    for (i, site) in enumerate(sites)
        if haskey(site_to_sc, site) && haskey(passability_per_subcatchment, site_to_sc[site])
            passability_vector[i] = passability_per_subcatchment[site_to_sc[site]]
        end
    end
    return passability_vector
end

function classify_species_indices(all_species, target_group)
    return [findfirst(==(sp), all_species) for sp in target_group if sp in all_species]
end

# =============================================================================
# --- Run Single Simulation ---
# =============================================================================

function run_single_simulation(data_base, temp_increase, upstream_cost, passability_dict;
    simulation_years=SIMULATION_YEARS)

    t_span = (0.0, Float64(simulation_years * DAYS_PER_YEAR))

    params_copy = MetacommunityParams(
        data_base.params.n_sites,
        data_base.params.n_species,
        copy(data_base.params.interaction_matrix),
        copy(data_base.params.dispersal_matrix),
        copy(data_base.params.dispersal_scaling),
        copy(data_base.params.intrinsic_growth_rates),
        copy(data_base.params.temperatures) .+ temp_increase,
        copy(data_base.params.habitat_suitability),
        copy(data_base.params.thermal_optima),
        copy(data_base.params.thermal_sigmas),
        copy(data_base.params.carrying_capacity)
    )

    passability_vector = build_dam_passability_vector(data_base.site_df, data_base.sites;
        passability_per_subcatchment=passability_dict)

    modified_dams = copy(data_base.dams)
    for j in 1:params_copy.n_sites
        for i in 1:params_copy.n_sites
            if i != j && modified_dams[i, j] < 1.0
                modified_dams[i, j] *= passability_vector[j]
                modified_dams[i, j] = min(1.0, modified_dams[i, j])
            end
        end
    end

    new_dispersal_matrix = precompute_dispersal_matrix(
        params_copy.n_sites,
        Matrix(data_base.distance_matrix),
        data_base.elevations,
        upstream_cost,
        modified_dams,
        data_base.species
    )

    params_final = MetacommunityParams(
        params_copy.n_sites,
        params_copy.n_species,
        params_copy.interaction_matrix,
        new_dispersal_matrix,
        params_copy.dispersal_scaling,
        params_copy.intrinsic_growth_rates,
        params_copy.temperatures,
        params_copy.habitat_suitability,
        params_copy.thermal_optima,
        params_copy.thermal_sigmas,
        params_copy.carrying_capacity
    )

    density_cols = [Symbol("$(sp)_DEN") for sp in data_base.species]
    density_df_filtered = filter(row -> row.CODIGO in data_base.sites, data_base.density_df)
    u0 = Matrix(density_df_filtered[:, density_cols])
    replace!(u0, NaN => 0.0)
    u0 = max.(u0, 0.0)
    u0_flat = vec(u0)

    prob = ODEProblem(metacommunity_ode!, u0_flat, t_span, params_final)

    function positivity_condition(u, t, integrator)
        any(x -> x < 0, u)
    end
    function positivity_affect!(integrator)
        integrator.u .= max.(integrator.u, 0.0)
    end
    positivity_cb = DiscreteCallback(positivity_condition, positivity_affect!; save_positions=(false, true))

    sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, saveat=0:1.0:t_span[2], callback=positivity_cb)

    return sol, passability_vector
end

# =============================================================================
# --- Main Sweep Loop ---
# =============================================================================

println("="^60)
println("Sensitivity Analysis: Temperature × Upstream Cost × Passability")
println("="^60)

base_output_dir = "results/sensitivity_temp_passability"
mkpath(base_output_dir)

println("\nLoading base data (once)...")
data_base = prepare_ode_data(
    upstream_cost = 0.05,
    cedex_var_file = CEDEX_VAR_FILE,
    cedex_esc_uts_file = CEDEX_UTS_FILE,
    obstacles_file = OBSTACLES_FILE,
    obstacle_mode = OBSTACLE_MODE,
    obstacle_matching_tolerance = OBSTACLE_MATCHING_TOLERANCE,
    obstacle_passability = OBSTACLE_PASSABILITY,
    obstacle_downstream_passability = OBSTACLE_DOWNSTREAM_PASSABILITY
)

native_idx = classify_species_indices(data_base.species, NATIVE_SPECIES)
invasive_idx = classify_species_indices(data_base.species, INVASIVE_SPECIES)

total_runs = length(temperature_increases) * length(upstream_costs) * length(passability_scenarios)
current_run = 0

for dt in temperature_increases
    for uc in upstream_costs
        for (pass_name, pass_dict) in passability_scenarios
            global current_run += 1
            run_label = "dT=$(dt)C_uc=$(uc)_pass=$(pass_name)"
            println("\n[$current_run/$total_runs] Running: $run_label")

            run_dir = joinpath(base_output_dir, "dT_$(dt)C", "uc_$(uc)", pass_name)
            mkpath(run_dir)

            sol, pass_vec = run_single_simulation(
                data_base, dt, uc, pass_dict
            )

            output_jld2 = joinpath(run_dir, "simulation_output.jld2")
            jldsave(output_jld2;
                sol_t = sol.t,
                sol_u = sol.u,
                sites = data_base.sites,
                species = data_base.species,
                temperature_increase = dt,
                upstream_cost = uc,
                simulation_years = SIMULATION_YEARS,
                passability_scenario = pass_name,
                passability_vector = pass_vec,
                cedex_var_file = CEDEX_VAR_FILE,
                cedex_uts_file = CEDEX_UTS_FILE,
                obstacles_file = OBSTACLES_FILE,
                obstacle_mode = string(OBSTACLE_MODE),
                obstacle_matching_tolerance_m = OBSTACLE_MATCHING_TOLERANCE,
                obstacle_passability = OBSTACLE_PASSABILITY,
                obstacle_downstream_passability = OBSTACLE_DOWNSTREAM_PASSABILITY,
                obstacle_matched_count = count(data_base.obstacle_mapping_diagnostics.matched),
                obstacle_total_count = nrow(data_base.obstacle_mapping_diagnostics)
            )

            fig_biomass = plot_avg_total_biomass(sol, data_base.sites, data_base.species)
            save_figure(fig_biomass, joinpath(run_dir, "avg_total_biomass.png"))

            fig_richness = plot_avg_species_richness(sol, data_base.sites, data_base.species)
            save_figure(fig_richness, joinpath(run_dir, "avg_species_richness.png"))

            fig_combined = plot_combined_analysis(sol, data_base.site_df, data_base.sites, data_base.species, data_base.distance_matrix)
            save_figure(fig_combined, joinpath(run_dir, "combined_analysis.png"))

            plot_richness_change_per_site(sol, data_base.species, data_base.sites, data_base.site_df,
                native_idx, invasive_idx, report_years,
                joinpath(run_dir, "richness_change_per_site.png"); days_per_year=DAYS_PER_YEAR)

            plot_richness_change_per_subcatchment(sol, data_base.species, data_base.sites, data_base.site_df,
                native_idx, invasive_idx, report_years,
                joinpath(run_dir, "richness_change_per_subcatchment.png"); days_per_year=DAYS_PER_YEAR)

            plot_richness_timeseries_grid(sol, data_base.species, data_base.sites, data_base.site_df,
                native_idx, invasive_idx, data_base.params.n_sites, data_base.params.n_species,
                joinpath(run_dir, "richness_timeseries_grid.png"); days_per_year=DAYS_PER_YEAR)

            println("  Saved to: $run_dir")
        end
    end
end

println("\n" * "="^60)
println("Sensitivity analysis complete. $total_runs runs saved to '$base_output_dir'")
println("="^60)
