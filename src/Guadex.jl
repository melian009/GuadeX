module Guadex

using Statistics
using DataFrames
using CSV
using Graphs
using GraphMakie
using CairoMakie
using NetworkLayout
using Colors
using DifferentialEquations
using LinearAlgebra
using SparseArrays
using Makie
using GeoMakie

# Export main functions
export build_stream_graph, get_graph_statistics, find_upstream_sites, find_downstream_sites, visualize_stream_graph, analyze_stream_connectivity, build_nearest_neighbor_connections!, build_threshold_connections!, build_mst_connections!, analyze_stream_connectivity
export plot_catchment_network
export metacommunity_ode!, MetacommunityParams, precompute_dispersal_matrix

# Export data preparation functions
export prepare_ode_data, save_ode_data
export load_species_characteristics, load_site_data, load_species_density_data
export load_interaction_matrix, build_distance_matrix
export build_elevation_vector, build_dam_passability_matrix
export extract_site_temperatures, extract_habitat_suitability
export build_intrinsic_growth_rates

# Export visualization functions
export plot_ode_solution, plot_total_biomass, plot_species_richness
export plot_sites_map, plot_site_connectivity_map, plot_subcatchment_network
export plot_combined_analysis, save_figure

include("graph_construction.jl")
include("visualize_graph.jl")
include("ode_model.jl")
include("data_preparation.jl")
include("visualization.jl")

end # module Guadex
