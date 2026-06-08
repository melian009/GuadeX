using Test
using Guadex
using SparseArrays
using DataFrames
using CSV
using DifferentialEquations
using LinearAlgebra
using Statistics

# Shared constants for test files
const DATA_DIR = joinpath(@__DIR__, "../data")
const DISTANCE_FILE = joinpath(DATA_DIR, "Matrix_distances_1037puntos_BRUTO_FINAL.csv")
const CONNECTIVITY_FILE = joinpath(DATA_DIR, "ConnectivityUTM.csv")
const ENVIRONMENTAL_FILE = joinpath(DATA_DIR, "ABIOTIC", "Matriz_Ambiental_Data.csv")
const DENSITY_FILE = joinpath(DATA_DIR, "BIOTIC", "FishDensity_and_Juveniles_Matrix.csv")
const SPECIES_CHARS_FILE = joinpath(DATA_DIR, "ABIOTIC", "caracteristicas_peces_Guadalquivir_03-04-2018.csv")
const INTERACTION_FILE = joinpath(DATA_DIR, "BIOTIC", "Interacciones_peces_Guadalquivir_03-04-2018_ENG.csv")

include("test_graph_construction.jl")
include("test_ode.jl")
include("test_data_preparation.jl")
include("test_visualization.jl")
