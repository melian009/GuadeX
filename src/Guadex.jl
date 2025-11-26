module Guadex

using Statistics
using DataFrames
using CSV
using Graphs
using GraphPlot
using Colors
using ColorVectorSpace

# Export main functions
export build_stream_graph, add_flow_direction!, get_graph_statistics, find_upstream_sites, find_downstream_sites, visualize_stream_graph, analyze_stream_connectivity, build_nearest_neighbor_connections!, build_threshold_connections!, build_mst_connections!, analyze_stream_connectivity

include("graph_construction.jl")


end # module Guadex
