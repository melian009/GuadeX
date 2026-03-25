using Pkg; Pkg.activate(".");
using DifferentialEquations
using DataFrames
using CSV
using SparseArrays
using LinearAlgebra
using Guadex


# --- Configuration ---
const UPSTREAM_COST = 0.01
const DISPERSAL_INTENSITY = 0.1
const T_SPAN = (0.0, 100.0) # Time span for simulation

# --- Data Preparation ---
# This function loads all necessary data and returns a NamedTuple
# containing the MetacommunityParams struct and other data.
data = prepare_ode_data(
    upstream_cost = UPSTREAM_COST,
    dispersal_intensity = DISPERSAL_INTENSITY
)

# --- Initial Conditions ---
# The model expects u as a flattened vector where u[i, s] is the population of species s at site i.
# We extract the initial densities from the density_df.
# We need to ensure the order of species matches the order in the params.
# The params are built using the same species order as in the density data.
density_cols = [Symbol("$(sp)_DEN") for sp in data.species]
# Filter density_df to match sites
density_df_filtered = filter(row -> row.CODIGO in data.sites, data.density_df)
u0 = Matrix(density_df_filtered[:, density_cols])

# Flatten for the ODE solver
u0_flat = vec(u0)

# --- Solve ODE ---
# Define the ODE problem
prob = ODEProblem(metacommunity_ode!, u0_flat, T_SPAN, data.params)

# Solve the ODE
println("Starting simulation...")
sol = solve(prob, Tsit5())
println("Simulation finished.")

# The solution `sol` contains the time series of the population densities.
# You can access the results using `sol.u`.

# --- Visualization ---
println("\nGenerating visualizations...")

# Create output directory for figures
output_dir = "results/figures"
mkpath(output_dir)

# Plot ODE solution: population dynamics over time
println("  - Plotting ODE solution...")
fig_ode = plot_ode_solution(sol, data.sites, data.species,
    max_species_to_plot = 6,
    max_sites_to_plot = 4)
save_figure(fig_ode, joinpath(output_dir, "ode_solution.png"))

# Plot total biomass over time
println("  - Plotting total biomass...")
fig_biomass = plot_total_biomass(sol, data.sites, data.species)
save_figure(fig_biomass, joinpath(output_dir, "total_biomass.png"))

# Plot species richness over time
println("  - Plotting species richness...")
fig_richness = plot_species_richness(sol, data.sites, data.species)
save_figure(fig_richness, joinpath(output_dir, "species_richness.png"))

# Plot sites map with elevation coloring
println("  - Plotting sites map...")
fig_sites = plot_sites_map(data.site_df, color_by = :ALTITUD)
save_figure(fig_sites, joinpath(output_dir, "sites_map.png"))

# Plot site connectivity network
println("  - Plotting site connectivity network...")
fig_connectivity = plot_site_connectivity_map(data.site_df, data.distance_matrix, data.sites)
save_figure(fig_connectivity, joinpath(output_dir, "site_connectivity.png"))

# Plot subcatchment network structure
println("  - Plotting subcatchment network...")
fig_subcatchment = plot_subcatchment_network(data.site_df, data.sites)
save_figure(fig_subcatchment, joinpath(output_dir, "subcatchment_network.png"))

# Create combined analysis figure
println("  - Creating combined analysis...")
fig_combined = plot_combined_analysis(sol, data.site_df, data.sites, data.species, data.distance_matrix)
save_figure(fig_combined, joinpath(output_dir, "combined_analysis.png"))

println("\nAll figures saved to '$output_dir' directory")
