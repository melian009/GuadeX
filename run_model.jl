using Pkg; Pkg.activate(".");
using DifferentialEquations
using DataFrames
using CSV
using SparseArrays
using LinearAlgebra
using JLD2
using Guadex


# --- Configuration ---
const UPSTREAM_COST = 0.05
const DISPERSAL_INTENSITY = 1.0

# Time configuration: 1 time unit = 1 day
const DAYS_PER_YEAR = 365
const SIMULATION_YEARS = 1
const T_SPAN = (0.0, Float64(SIMULATION_YEARS * DAYS_PER_YEAR)) # Time span in days

# =============================================================================
# --- Management Scenarios ---
# These vectors allow you to test the effect of management interventions
# on fish populations without modifying the ODE model structure.
# =============================================================================

"""
    build_exploitation_vector(site_df, sites; exploitation_per_subcatchment)

Create a vector of exploitation factors (0-1) for each site.
- Values < 1.0 reduce the intrinsic growth rate (simulating harvesting/mortality)
- Value of 0.5 means 50% exploitation (growth rate reduced by half)
- Value of 1.0 means no exploitation (baseline)

# Arguments
- `site_df`: DataFrame with site data including CODIGO and CODIGO_S columns
- `sites`: Vector of site codes in order
- `exploitation_per_subcatchment`: Dict mapping subcatchment ID to exploitation factor (0-1)
                                    Sites in unlisted subcatchments get factor 1.0 (no exploitation)
"""
function build_exploitation_vector(site_df, sites; exploitation_per_subcatchment::Dict{Float64, Float64}=Dict{Float64, Float64}())
    n_sites = length(sites)
    exploitation_vector = ones(n_sites)

    # Create site to subcatchment mapping
    site_to_sc = Dict(row.CODIGO => row.CODIGO_S for row in eachrow(site_df))

    for (i, site) in enumerate(sites)
        if haskey(site_to_sc, site) && haskey(exploitation_per_subcatchment, site_to_sc[site])
            exploitation_vector[i] = exploitation_per_subcatchment[site_to_sc[site]]
        end
    end

    return exploitation_vector
end

"""
    build_dam_passability_vector(site_df, sites; passability_per_subcatchment)

Create a vector of dam passability multipliers (0-1) for each site.
- Values < 1.0 reduce dam passability (simulating improved dam passage for fish)
- Value of 0.1 means very poor passability (dams block movement)
- Value of 1.0 means full passability (baseline, no management)
- Value of 2.0 would mean improved passability (fish ladder, etc.)

# Arguments
- `site_df`: DataFrame with site data including CODIGO and CODIGO_S columns
- `sites`: Vector of site codes in order
- `passability_per_subcatchment`: Dict mapping subcatchment ID to passability multiplier (0-1 or >1)
                                   Sites in unlisted subcatchments get factor 1.0 (baseline)
"""
function build_dam_passability_vector(site_df, sites; passability_per_subcatchment::Dict{Float64, Float64}=Dict{Float64, Float64}())
    n_sites = length(sites)
    passability_vector = ones(n_sites)

    # Create site to subcatchment mapping
    site_to_sc = Dict(row.CODIGO => row.CODIGO_S for row in eachrow(site_df))

    for (i, site) in enumerate(sites)
        if haskey(site_to_sc, site) && haskey(passability_per_subcatchment, site_to_sc[site])
            passability_vector[i] = passability_per_subcatchment[site_to_sc[site]]
        end
    end

    return passability_vector
end

"""
    apply_exploitation_to_growth_rates!(intrinsic_growth_rates, exploitation_vector)

Apply exploitation factors to intrinsic growth rates in-place.
Each column (species) is scaled by the site-specific exploitation factor.
"""
function apply_exploitation_to_growth_rates!(intrinsic_growth_rates, exploitation_vector)
    n_sites = length(exploitation_vector)
    @assert size(intrinsic_growth_rates, 1) == n_sites "Mismatch between exploitation vector and growth rate matrix"

    for i in 1:n_sites
        intrinsic_growth_rates[i, :] .*= exploitation_vector[i]
    end

    return intrinsic_growth_rates
end

"""
    apply_passability_to_dams!(dams, passability_vector)

Apply passability multipliers to dam passability matrix in-place.
The passability vector values replace the baseline values for each origin site.
"""
function apply_passability_to_dams!(dams, passability_vector)
    n_sites = length(passability_vector)
    @assert size(dams, 1) == n_sites "Mismatch between passability vector and dam matrix"

    for j in 1:n_sites  # origin sites (columns)
        for i in 1:n_sites  # destination sites (rows)
            if i != j
                # Scale the existing passability by the origin site's passability factor
                dams[i, j] *= passability_vector[j]
            end
        end
    end

    return dams
end

# =============================================================================
# --- Define Management Scenarios ---
# Modify these dictionaries to test different management interventions.
# Keys are subcatchment IDs (Float64), values are management factors.
# =============================================================================

all_subcatchments = [1.1, 1.2, 1.3, 1.4, 2.2, 3.0, 6.0, 7.0, 9.0, 10.0, 11.3, 11.5, 11.2, 11.1, 11.4, 12.0, 13.3, 13.2, 13.1, 14.0, 15.2, 15.1, 16.0, 17.0, 18.1, 18.3, 18.2, 19.0, 20.0, 21.0, 22.1, 22.4, 22.2, 23.0, 24.1, 24.2, 25.0, 26.1, 26.2, 26.3, 26.4, 26.5, 26.6, 27.0, 28.1, 28.2, 28.3, 28.4, 29.0, 30.1, 30.2, 30.3, 30.4, 30.5, 30.6, 30.7, 31.0, 32.4, 32.3, 32.1, 32.2, 33.0, 34.1, 34.2, 34.3, 34.4, 35.0, 36.1, 36.2, 36.3, 36.4, 37.0, 38.0, 39.0, 40.1, 40.2, 41.0]

# Example: Apply 30% exploitation (r reduced by 30%) in subcatchment 3.1
# and 50% exploitation in subcatchment 4.2
# Value of 0.7 means growth rate is multiplied by 0.7 (30% reduction)
# Value of 0.0 would mean complete closure (no reproduction)
example_exploitation_scenario = Dict{Float64, Float64}(
    # 3.1 => 0.7,   # 30% exploitation in subcatchment 3.1
    # 4.2 => 0.5,   # 50% exploitation in subcatchment 4.2
    a => 0.9 for a in all_subcatchments
)

# Example: Improve dam passability (multiplier) in subcatchments
# Value > 1.0 improves passability (e.g., 1.5 = 50% improvement via fish ladder)
# Value < 1.0 reduces passability (e.g., 0.5 = 50% reduction, blocked dam)
example_passability_scenario = Dict{Float64, Float64}(
    # 2.0 => 1.5,   # Improve passability by 50% in subcatchment 2.0
    # 3.1 => 2.0,   # Double passability (fish ladder) in subcatchment 3.1
    a => 1.0 for a in all_subcatchments
)

# --- Data Preparation ---
# This function loads all necessary data and returns a NamedTuple
# containing the MetacommunityParams struct and other data.
data = prepare_ode_data(
    upstream_cost = UPSTREAM_COST,
    dispersal_intensity = DISPERSAL_INTENSITY
)

# --- Apply Management Scenarios ---
println("\n--- Applying Management Scenarios ---")

# Build management vectors
exploitation_vector = build_exploitation_vector(data.site_df, data.sites;
    exploitation_per_subcatchment=example_exploitation_scenario)
passability_vector = build_dam_passability_vector(data.site_df, data.sites;
    passability_per_subcatchment=example_passability_scenario)

println("Exploitation vector (first 10 sites): $(exploitation_vector[1:min(10, end)])")
println("Passability vector (first 10 sites): $(passability_vector[1:min(10, end)])")

# Apply exploitation to growth rates (this modifies params in-place since intrinsic_growth_rates is a mutable matrix)
apply_exploitation_to_growth_rates!(data.params.intrinsic_growth_rates, exploitation_vector)
println("Exploitation factors applied to growth rates")

# For dam passability, we need to recompute the dispersal matrix since it was pre-computed
# First, apply the passability multiplier to the dams matrix
modified_dams = copy(data.dams)
for j in 1:data.params.n_sites  # origin sites (columns)
    for i in 1:data.params.n_sites  # destination sites (rows)
        if i != j
            # Scale the existing passability by the origin site's passability factor
            modified_dams[i, j] *= passability_vector[j]
        end
    end
end

# Recompute dispersal matrix with modified dam passability
println("Recomputing dispersal matrix with modified dam passability...")
new_dispersal_matrix = precompute_dispersal_matrix(
    data.params.n_sites,
    Matrix(data.distance_matrix),
    data.elevations,
    UPSTREAM_COST,
    DISPERSAL_INTENSITY,
    modified_dams
)

# Create new params with updated dispersal matrix
# Note: params is immutable, so we create a new instance
data_with_management = (
    params = MetacommunityParams(
        data.params.n_sites,
        data.params.n_species,
        data.params.interaction_matrix,
        new_dispersal_matrix,
        data.params.dispersal_scaling,
        data.params.intrinsic_growth_rates,  # Already modified in-place
        data.params.temperatures,
        data.params.habitat_suitability,
        data.params.thermal_optima,
        data.params.thermal_sigmas
    ),
    sites = data.sites,
    species = data.species,
    distance_matrix = data.distance_matrix,
    elevations = data.elevations,
    dams = modified_dams,
    site_df = data.site_df,
    species_chars_df = data.species_chars_df,
    density_df = data.density_df
)

println("Management scenarios applied successfully!")

# --- Initial Conditions ---
# The model expects u as a flattened vector where u[i, s] is the population of species s at site i.
# We extract the initial densities from the density_df.
# We need to ensure the order of species matches the order in the params.
# The params are built using the same species order as in the density data.
density_cols = [Symbol("$(sp)_DEN") for sp in data_with_management.species]
# Filter density_df to match sites
density_df_filtered = filter(row -> row.CODIGO in data_with_management.sites, data_with_management.density_df)
u0 = Matrix(density_df_filtered[:, density_cols])

# Flatten for the ODE solver
u0_flat = vec(u0)

# --- Solve ODE ---
# Define the ODE problem (using data_with_management which has modified params)
prob = ODEProblem(metacommunity_ode!, u0_flat, T_SPAN, data_with_management.params)

# Solve the ODE with adaptive time stepping (Tsit5)
# but save at regular daily intervals for consistent time series
# Using reltol=1e-6 and abstol=1e-6 for accurate solution
println("Starting simulation for $SIMULATION_YEARS years ($(SIMULATION_YEARS * DAYS_PER_YEAR) days)...")
sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, saveat=0:1.0:T_SPAN[2]) # Save at daily intervals
println("Simulation finished. Solution has $(length(sol.t)) time points.")

# The solution `sol` contains the time series of the population densities.
# You can access the results using `sol.u`.

# --- Visualization ---
println("\nGenerating visualizations...")

# Create output directory for figures
output_dir = "results/figures"
mkpath(output_dir)

# Plot ODE solution: population dynamics over time
println("  - Plotting ODE solution...")
fig_ode = plot_ode_solution(sol, data_with_management.sites, data_with_management.species,
    max_species_to_plot = 6,
    max_sites_to_plot = 4)
save_figure(fig_ode, joinpath(output_dir, "ode_solution.png"))

# Plot average total biomass over time
println("  - Plotting average total biomass...")
fig_biomass = plot_avg_total_biomass(sol, data_with_management.sites, data_with_management.species)
save_figure(fig_biomass, joinpath(output_dir, "avg_total_biomass.png"))

# Plot average species richness over time
println("  - Plotting average species richness...")
fig_richness = plot_avg_species_richness(sol, data_with_management.sites, data_with_management.species)
save_figure(fig_richness, joinpath(output_dir, "avg_species_richness.png"))

# Plot sites map with elevation coloring
println("  - Plotting sites map...")
fig_sites = plot_sites_map(data_with_management.site_df, color_by = :ALTITUD)
save_figure(fig_sites, joinpath(output_dir, "sites_map.png"))

# Plot site connectivity network
println("  - Plotting site connectivity network...")
fig_connectivity = plot_site_connectivity_map(data_with_management.site_df, data_with_management.distance_matrix, data_with_management.sites)
save_figure(fig_connectivity, joinpath(output_dir, "site_connectivity.png"))

# Plot subcatchment network structure
println("  - Plotting subcatchment network...")
fig_subcatchment = plot_subcatchment_network(data_with_management.site_df, data_with_management.sites)
save_figure(fig_subcatchment, joinpath(output_dir, "subcatchment_network.png"))

# Create combined analysis figure
println("  - Creating combined analysis...")
fig_combined = plot_combined_analysis(sol, data_with_management.site_df, data_with_management.sites, data_with_management.species, data_with_management.distance_matrix)
save_figure(fig_combined, joinpath(output_dir, "combined_analysis.png"))

println("\nAll figures saved to '$output_dir' directory")

# --- Save simulation output ---
output_jld2 = joinpath(output_dir, "simulation_output.jld2")
println("\nSaving simulation output to: $output_jld2")
jldsave(output_jld2;
    sol_t = sol.t,
    sol_u = sol.u,
    sites = data_with_management.sites,
    species = data_with_management.species,
    upstream_cost = UPSTREAM_COST,
    dispersal_intensity = DISPERSAL_INTENSITY,
    simulation_years = SIMULATION_YEARS,
    exploitation_vector = exploitation_vector,
    passability_vector = passability_vector
)
println("Simulation output saved.")
