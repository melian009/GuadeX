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

const DAYS_PER_YEAR = 365
const SIMULATION_YEARS = 3
const T_SPAN = (0.0, Float64(SIMULATION_YEARS * DAYS_PER_YEAR))

const NATIVE_SPECIES = ["AB", "AH", "SP", "PW", "LS", "SA", "IL", "CP", "IO"]
const INVASIVE_SPECIES = ["GH", "MS", "LG", "CC", "CG", "AM", "OM", "EL", "GL", "TT"]

all_subcatchments = [1.1, 1.2, 1.3, 1.4, 2.2, 3.0, 6.0, 7.0, 9.0, 10.0, 11.3, 11.5, 11.2, 11.1, 11.4, 12.0, 13.3, 13.2, 13.1, 14.0, 15.2, 15.1, 16.0, 17.0, 18.1, 18.3, 18.2, 19.0, 20.0, 21.0, 22.1, 22.4, 22.2, 23.0, 24.1, 24.2, 25.0, 26.1, 26.2, 26.3, 26.4, 26.5, 26.6, 27.0, 28.1, 28.2, 28.3, 28.4, 29.0, 30.1, 30.2, 30.3, 30.4, 30.5, 30.6, 30.7, 31.0, 32.4, 32.3, 32.1, 32.2, 33.0, 34.1, 34.2, 34.3, 34.4, 35.0, 36.1, 36.2, 36.3, 36.4, 37.0, 38.0, 39.0, 40.1, 40.2, 41.0]

# =============================================================================
# --- Define Parameter Grid ---
# =============================================================================

temperature_increases = [0.0, 1.0, 2.0, 3.0]

upstream_costs = [0.01, 0.05, 0.1, 0.5]

passability_scenarios = Dict{String, Dict{Float64, Float64}}(
    "baseline" => Dict{Float64, Float64}(a => 1.0 for a in all_subcatchments),
    "improved_passability" => Dict{Float64, Float64}(a => 1.5 for a in all_subcatchments),
    "reduced_passability" => Dict{Float64, Float64}(a => 0.5 for a in all_subcatchments),
    "blocked" => Dict{Float64, Float64}(a => 0.1 for a in all_subcatchments),
)

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
data_base = prepare_ode_data(upstream_cost = 0.05)

native_idx = classify_species_indices(data_base.species, NATIVE_SPECIES)
invasive_idx = classify_species_indices(data_base.species, INVASIVE_SPECIES)
report_years = [1, 2, 3]

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
                passability_vector = pass_vec
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
