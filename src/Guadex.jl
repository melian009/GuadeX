module Guadex

using Statistics
using DataFrames
using CSV
using Graphs
using GraphMakie
using CairoMakie
using NetworkLayout
using Colors

# Export main functions
export build_stream_graph, add_flow_direction!, get_graph_statistics, find_upstream_sites, find_downstream_sites, visualize_stream_graph, analyze_stream_connectivity, build_nearest_neighbor_connections!, build_threshold_connections!, build_mst_connections!, analyze_stream_connectivity
export plot_catchment_network

include("graph_construction.jl")
include("visualize_graph.jl")


end # module Guadex
