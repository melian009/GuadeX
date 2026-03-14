using Pkg; Pkg.activate(".");
using DifferentialEquations
using DataFrames
using CSV
using SparseArrays
using LinearAlgebra
using Guadex


# --- Configuration ---
# Easily adjustable variables
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
