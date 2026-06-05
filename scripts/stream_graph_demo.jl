# Demonstration script for stream graph construction and visualization
using Pkg
Pkg.activate(".")

using Guadex
using DataFrames
using CSV

# Set up data paths
data_dir = joinpath(isdefined(Main, :__DIR__) ? Main.__DIR__ : pwd(), "data")
distance_file = joinpath(data_dir, "Matrix_distances_1037puntos_BRUTO_FINAL.csv")
connectivity_file = joinpath(data_dir, "ConnectivityUTM.csv")

println("=== Stream Graph Construction and Visualization Demo ===")
println("Distance file: $distance_file")
println("Connectivity file: $connectivity_file")

# Build the basic stream graph using nearest neighbors method
println("\n--- Building Stream Graph (Nearest Neighbors) ---")
graph, site_to_index, distance_data = build_stream_graph(distance_file, max_distance=Inf, connection_method=:nearest_neighbors, connectivity_file=connectivity_file)

# Get basic statistics
stats = get_graph_statistics(graph)
println("\nBasic Graph Statistics:")
println("  Number of nodes (sites): $(stats.num_nodes)")
println("  Number of edges: $(stats.num_edges)")
println("  Graph density: $(round(stats.density, digits=4))")
println("  Number of connected components: $(stats.num_components)")
println("  Largest component size: $(stats.largest_component_size)")


# # Compare different connection methods
# println("\n--- Comparing Connection Methods ---")

# # Method 1: Threshold distance
# println("\n1. Threshold Distance Method:")
# graph_threshold, _, _ = build_stream_graph(distance_file, max_distance=5000.0, connection_method=:threshold_distance)
# stats_threshold = get_graph_statistics(graph_threshold)
# println("   Nodes: $(stats_threshold.num_nodes), Edges: $(stats_threshold.num_edges)")

# # Method 2: Minimum spanning tree
# println("\n2. Minimum Spanning Tree Method:")
# graph_mst, _, _ = build_stream_graph(distance_file, connection_method=:minimum_spanning_tree)
# stats_mst = get_graph_statistics(graph_mst)
# println("   Nodes: $(stats_mst.num_nodes), Edges: $(stats_mst.num_edges)")

# # Method 3: All connections (original method)
# println("\n3. All Connections Method (Original):")
# graph_all, _, _ = build_stream_graph(distance_file, max_distance=10000.0, connection_method=:all_connections)
# stats_all = get_graph_statistics(graph_all)
# println("   Nodes: $(stats_all.num_nodes), Edges: $(stats_all.num_edges)")

# # Get updated statistics
# stats = get_graph_statistics(graph)
# println("\nUpdated Graph Statistics (with flow direction):")
# println("  Number of nodes (sites): $(stats.num_nodes)")
# println("  Number of edges: $(stats.num_edges)")
# println("  Graph density: $(round(stats.density, digits=4))")
# println("  Number of connected components: $(stats.num_components)")
# println("  Largest component size: $(stats.largest_component_size)")

# Analyze connectivity patterns
println("\n--- Analyzing Stream Connectivity ---")
analysis = analyze_stream_connectivity(graph, site_to_index)
println("  Critical nodes (high centrality): $(length(analysis.critical_nodes))")
println("  Source nodes (headwaters): $(length(analysis.source_nodes))")
println("  Sink nodes (mouths/confluences): $(length(analysis.sink_nodes))")
println("  Average indegree: $(round(analysis.avg_indegree, digits=3))")
println("  Average outdegree: $(round(analysis.avg_outdegree, digits=3))")

# Show some example nodes
println("\n--- Example Sites ---")
println("  Sample source nodes (headwaters): $(analysis.source_nodes[1:min(3, length(analysis.source_nodes))])")
println("  Sample sink nodes (mouths): $(analysis.sink_nodes[1:min(3, length(analysis.sink_nodes))])")

# Find upstream/downstream for a specific site
if length(site_to_index) > 0
    sample_site = collect(keys(site_to_index))[1]
    println("\n--- Upstream/Downstream Analysis for Site: $sample_site ---")
    upstream = find_upstream_sites(graph, site_to_index, sample_site)
    downstream = find_downstream_sites(graph, site_to_index, sample_site)
    println("  Upstream sites: $(length(upstream))")
    println("  Downstream sites: $(length(downstream))")

    if length(upstream) > 0
        println("  Sample upstream sites: $(upstream[1:min(3, length(upstream))])")
    end
    if length(downstream) > 0
        println("  Sample downstream sites: $(downstream[1:min(3, length(downstream))])")
    end
end

# Visualize the graph (optional - requires display backend)
println("\n--- Graph Visualization ---")
try
    # # Use plot_catchment_network for visualization
    labels = collect(keys(site_to_index))
    # # Extract coordinates for visualization
    # site_coords = Dict{String,Tuple{Float64,Float64}}()
    # for row in eachrow(CSV.read(connectivity_file, DataFrame))
    #     site_coords[row.CODIGO] = (row.UTMX, row.UTMY)
    # end
    # coordinates = [site_coords[label] for label in labels if haskey(site_coords, label)]
    coordinates = nothing  # For now, we skip fixed coordinates because subcatchments do not have coordinates, only sites have.
    labels = nothing # Skip labels for clarity in large graphs

    f = plot_catchment_network(graph; labels=labels, coordinates=coordinates)
    display(f)
    println("  Graph visualization created successfully!")
    println("  Note: To save the visualization, use CairoMakie.save or similar")
    # Save
    CairoMakie.save("stream_graph_visualization.png", f)
catch e
    println("  Visualization failed: $e")
    println("  This may be due to missing display backend or graphics libraries")
end

println("\n=== Demo Complete ===")
println("The stream graph has been successfully constructed with:")
println("  - $(stats.num_nodes) sampling sites as nodes")
println("  - $(stats.num_edges) directed edges representing stream flow")
println("  - Flow direction based on elevation data")
println("  - Connectivity analysis identifying critical network components")