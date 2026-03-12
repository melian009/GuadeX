using Test
using Guadex
using SparseArrays

# Set up data paths
const data_dir = joinpath(@__DIR__, "../data")
const distance_file = joinpath(data_dir, "Matrix_distances_1037puntos_BRUTO_FINAL.csv")
const connectivity_file = joinpath(data_dir, "ConnectivityUTM.csv")

# Run test files
include("test_graph_construction.jl")
include("test_ode.jl")