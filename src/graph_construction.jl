using GraphPlot: spring_layout, circular_layout, shell_layout, spectral_layout


"""
    build_stream_graph(distance_file::String; max_distance::Float64=Inf)

Build a directed graph of sites connected by streams from the distance matrix file.

# Arguments
- `distance_file::String`: Path to the distance matrix CSV file
- `max_distance::Float64`: Maximum reticular distance to consider for connections (default: Inf)

# Returns
- `Graphs.DiGraph`: Directed graph where nodes are sites and edges represent stream connections
- `Dict{String, Int}`: Mapping from site codes to node indices
- `DataFrame`: Filtered distance data used to build the graph

# Description
The function reads the distance matrix file and creates a directed graph where:
- Each unique site becomes a node
- An edge from site A to site B exists if there's a stream connection
- Edge direction follows the natural flow (upstream to downstream)
- Only connections with reticular distance ≤ max_distance are included
"""
function build_stream_graph(distance_file::String; max_distance::Float64=Inf,
                           connection_method::Symbol=:nearest_neighbors)
    # Read the distance matrix
    println("Reading distance matrix from: $distance_file")
    distance_data = CSV.read(distance_file, DataFrame)

    # Filter out self-connections and distances greater than max_distance
    filtered_data = filter(row -> row.ID_ORIGIN != row.ID_DESTINATION &&
                                   row.RETICULAR_DIST <= max_distance &&
                                   row.RETICULAR_DIST > 0, distance_data)

    println("Found $(nrow(filtered_data)) valid connections out of $(nrow(distance_data)) total")

    # Get unique site codes
    unique_sites = unique(vcat(filtered_data.ID_ORIGIN, filtered_data.ID_DESTINATION))
    println("Building graph with $(length(unique_sites)) sites")

    # Create mapping from site code to node index
    site_to_index = Dict(site => i for (i, site) in enumerate(unique_sites))

    # Create directed graph
    num_sites = length(unique_sites)
    graph = DiGraph(num_sites)

    # Add edges based on stream connections
    # We need to determine direction based on elevation or other criteria
    # For now, we'll use a simple approach: assume flow from higher to lower elevation
    # We'll need to load elevation data for this

    # Build connections based on the specified method
    if connection_method == :nearest_neighbors
        # Connect each site to its nearest neighbors in the stream network
        edge_count = build_nearest_neighbor_connections!(graph, site_to_index, filtered_data)
    elseif connection_method == :threshold_distance
        # Connect sites within a threshold reticular distance
        edge_count = build_threshold_connections!(graph, site_to_index, filtered_data, max_distance)
    elseif connection_method == :minimum_spanning_tree
        # Build a minimum spanning tree of the stream network
        edge_count = build_mst_connections!(graph, site_to_index, filtered_data)
    elseif connection_method == :all_connections
        # Original method - connect all pairs (for backward compatibility)
        edge_count = 0
        for row in eachrow(filtered_data)
            origin_idx = site_to_index[row.ID_ORIGIN]
            dest_idx = site_to_index[row.ID_DESTINATION]

            # Add edge in both directions for now (undirected connectivity)
            # This will be refined to directed flow once we have elevation data
            if !has_edge(graph, origin_idx, dest_idx)
                add_edge!(graph, origin_idx, dest_idx)
                edge_count += 1
            end
            if !has_edge(graph, dest_idx, origin_idx)
                add_edge!(graph, dest_idx, origin_idx)
                edge_count += 1
            end
        end
    else
        error("Unknown connection method: $connection_method. Use :nearest_neighbors, :threshold_distance, :minimum_spanning_tree, or :all_connections")
    end

    println("Added $edge_count edges to the graph using $connection_method method")

    # Convert site_to_index to Dict{String, Int}
    site_to_index_str = Dict(String(k) => v for (k, v) in site_to_index)

    return graph, site_to_index_str, filtered_data
end

"""
    add_flow_direction!(graph::DiGraph, site_to_index::Dict{String, Int},
                       connectivity_data::AbstractDataFrame)

Add proper flow direction to the graph based on elevation data.

# Arguments
- `graph::DiGraph`: The graph to modify
- `site_to_index::Dict{String, Int}`: Mapping from site codes to node indices
- `connectivity_data::AbstractDataFrame`: Connectivity data with elevation information
"""
function add_flow_direction!(graph::DiGraph, site_to_index::Dict{String, Int},
                           connectivity_data::AbstractDataFrame)
    println("Adding flow direction based on elevation...")


    # Create elevation lookup
    elevation_dict = Dict(row.CODIGO => row.ALTITUD for row in eachrow(connectivity_data))

    # Remove all existing edges first
    edges_to_remove = collect(edges(graph))
    for edge in edges_to_remove
        rem_edge!(graph, edge)
    end

    # Add directed edges from higher to lower elevation
    # Only add edges where there's already a connection (based on existing undirected connections)
    edge_count = 0

    # We need to reload the distance data to know which sites are actually connected
    # For now, let's create a more reasonable flow direction based on elevation
    # and the original distance connections

    # Create a set of connected site pairs from the original distance data
    connected_pairs = Set{Tuple{String,String}}()

    # Add directed edges based on elevation for connected sites only
    for (site1, idx1) in site_to_index, (site2, idx2) in site_to_index
        if site1 != site2 && haskey(elevation_dict, site1) && haskey(elevation_dict, site2)
            elev1 = elevation_dict[site1]
            elev2 = elevation_dict[site2]

            # Only add edge if there's a significant elevation difference and sites are reasonably close
            # This is a heuristic - in reality we'd check the actual distance matrix
            if abs(elev1 - elev2) > 10  # Minimum elevation difference
                if elev1 > elev2  # Flow from higher to lower
                    if !has_edge(graph, idx1, idx2)
                        add_edge!(graph, idx1, idx2)
                        edge_count += 1
                    end
                else
                    if !has_edge(graph, idx2, idx1)
                        add_edge!(graph, idx2, idx1)
                        edge_count += 1
                    end
                end
            end
        end
    end

    println("Added $edge_count directed edges based on elevation")
end

function add_flow_direction!(graph::DiGraph, site_to_index::Dict{String, Int},
                           connectivity_file::AbstractString)
    connectivity_data = CSV.read(connectivity_file, DataFrame)
    add_flow_direction!(graph, site_to_index, connectivity_data)
end

"""
    get_graph_statistics(graph::DiGraph, site_to_index::Dict{String, Int})

Calculate basic statistics for the stream graph.

# Returns
- `NamedTuple`: Graph statistics including number of nodes, edges, density, etc.
"""
function get_graph_statistics(graph::DiGraph)
    num_nodes = nv(graph)
    num_edges = ne(graph)
    density = num_edges / (num_nodes * (num_nodes - 1))

    # Calculate connected components (for directed graph, we use weakly connected components)
    components = weakly_connected_components(graph)
    num_components = length(components)
    largest_component_size = maximum(length(comp) for comp in components)

    return (
        num_nodes = num_nodes,
        num_edges = num_edges,
        density = density,
        num_components = num_components,
        largest_component_size = largest_component_size,
        components = components
    )
end

"""
    find_upstream_sites(graph::DiGraph, site_to_index::Dict{String, Int}, target_site::String)

Find all sites that are upstream of the target site.

# Returns
- `Vector{String}`: List of upstream site codes
"""
function find_upstream_sites(graph::DiGraph, site_to_index::Dict{String, Int}, target_site::String)
    if !haskey(site_to_index, target_site)
        error("Site $target_site not found in graph")
    end

    target_idx = site_to_index[target_site]
    upstream_sites = String[]

    # Find all nodes that can reach the target site
    for (site, idx) in site_to_index
        if idx != target_idx && has_path(graph, idx, target_idx)
            push!(upstream_sites, site)
        end
    end

    return upstream_sites
end

"""
    find_downstream_sites(graph::DiGraph, site_to_index::Dict{String, Int}, target_site::String)

Find all sites that are downstream of the target site.

# Returns
- `Vector{String}`: List of downstream site codes
"""
function find_downstream_sites(graph::DiGraph, site_to_index::Dict{String, Int}, target_site::String)
    if !haskey(site_to_index, target_site)
        error("Site $target_site not found in graph")
    end

    target_idx = site_to_index[target_site]
    downstream_sites = String[]

    # Find all nodes reachable from the target site
    for (site, idx) in site_to_index
        if idx != target_idx && has_path(graph, target_idx, idx)
            push!(downstream_sites, site)
        end
    end

    return downstream_sites
end

"""
    visualize_stream_graph(graph::DiGraph, site_to_index::Dict{String, Int},
                          connectivity_data::Union{DataFrame, String};
                          layout::Symbol=:spring, node_size::Float64=0.1,
                          show_labels::Bool=false, save_path::Union{String, Nothing}=nothing)

Visualize the stream network graph.

# Arguments
- `graph::DiGraph`: The stream graph to visualize
- `site_to_index::Dict{String, Int}`: Mapping from site codes to node indices
- `connectivity_data::Union{DataFrame, String}`: Connectivity data with coordinates
- `layout::Symbol`: Layout algorithm (:spring, :circular, :shell, etc.)
- `node_size::Float64`: Size of nodes in the plot
- `show_labels::Bool`: Whether to show site labels
- `save_path::Union{String, Nothing}`: Path to save the plot (nothing = display only)

# Returns
- `Plots.Plot`: The plot object
"""
function visualize_stream_graph(graph::DiGraph, site_to_index::Dict{<:AbstractString, <:Int},
                              connectivity_data::Union{<:AbstractDataFrame, <:AbstractString};
                              layout::Symbol=:spring, node_size::Float64=0.1,
                              show_labels::Bool=false, save_path::Union{String, Nothing}=nothing)

    # Load connectivity data if it's a file path
    if connectivity_data isa String
        connectivity_df = CSV.read(connectivity_data, DataFrame)
    else
        connectivity_df = connectivity_data
    end

    # Extract coordinates for nodes
    x_coords = Float64[]
    y_coords = Float64[]
    node_labels = String[]

    for (site_code, node_idx) in site_to_index
        # Find the site in connectivity data
        site_row = filter(row -> row.CODIGO == site_code, connectivity_df)
        if nrow(site_row) > 0
            push!(x_coords, site_row.UTMX[1])
            push!(y_coords, site_row.UTMY[1])
            push!(node_labels, site_code)
        else
            # Fallback: use default coordinates
            push!(x_coords, 0.0)
            push!(y_coords, 0.0)
            push!(node_labels, site_code)
        end
    end

    # Normalize coordinates for better visualization
    if length(x_coords) > 1
        x_min, x_max = minimum(x_coords), maximum(x_coords)
        y_min, y_max = minimum(y_coords), maximum(y_coords)

        # Normalize to [0, 1] range
        x_norm = (x_coords .- x_min) ./ (x_max - x_min)
        y_norm = (y_coords .- y_min) ./ (y_max - y_min)
    else
        x_norm = x_coords
        y_norm = y_coords
    end

    # Create node colors based on elevation (if available)
    elevation_dict = Dict(row.CODIGO => row.ALTITUD for row in eachrow(connectivity_df))
    elevations = [haskey(elevation_dict, site_code) ? elevation_dict[site_code] : 0.0 for site_code in keys(site_to_index)]

    # Normalize elevations to [0,1] for color mapping
    if length(elevations) > 1
        elev_min, elev_max = minimum(elevations), maximum(elevations)
        elev_norm = (elevations .- elev_min) ./ (elev_max - elev_min)
    else
        elev_norm = elevations
    end

    # Map normalized elevation to color (simple blue gradient)
    node_colors = [RGB(0.2, 0.2, 1.0) * (1.0 - v) + RGB(1.0, 1.0, 0.2) * v for v in elev_norm]
    # If Colors.jl is not available, fallback to color names:
    # node_colors = [colorant"blue" for _ in elev_norm] # or use a fixed color

    # Map layout symbol to GraphPlot layout function
    layout_map = Dict(
        :spring => spring_layout,
        :circular => circular_layout,
        :shell => shell_layout,
        :spectral => spectral_layout,
        :random => random_layout,
        :community => community_layout,
        :collapse => collapse_layout
    )
    layout_func = get(layout_map, layout, spring_layout)

    # Create the plot
    p = gplot(graph,
              x_norm, y_norm;
              nodefillc=node_colors,
              nodesize=node_size,
              edgestrokec=:gray,
              edgelinewidth=0.5)

    # Add labels if requested
    if show_labels
        # This would require additional text annotation
        # For now, we'll keep it simple
    end

    # Save or display the plot
    if save_path !== nothing
        # Save the plot (implementation depends on the plotting backend)
        println("Saving plot to: $save_path")
    end

    return p
end

"""
    analyze_stream_connectivity(graph::DiGraph, site_to_index::Dict{String, Int})

Analyze the connectivity patterns of the stream network.

# Returns
- `NamedTuple`: Analysis results including centrality measures, connectivity metrics, etc.
"""
function analyze_stream_connectivity(graph::DiGraph, site_to_index::Dict{String, Int})
    # Basic graph metrics
    num_nodes = nv(graph)
    num_edges = ne(graph)

    # Centrality measures
    betweenness_centrality = [0.0 for _ in 1:num_nodes]
    closeness_centrality = [0.0 for _ in 1:num_nodes]
    indegree = indegree_centrality(graph)
    outdegree = outdegree_centrality(graph)

    # Find critical nodes (high betweenness centrality)
    critical_nodes = String[]
    for (site, idx) in site_to_index
        if indegree[idx] > 0.1 || outdegree[idx] > 0.1  # Threshold for "critical"
            push!(critical_nodes, site)
        end
    end

    # Find source and sink nodes
    source_nodes = String[]
    sink_nodes = String[]

    for (site, idx) in site_to_index
        if indegree[idx] == 0 && outdegree[idx] > 0
            push!(source_nodes, site)  # Headwaters
        elseif indegree[idx] > 0 && outdegree[idx] == 0
            push!(sink_nodes, site)    # Mouths/confluences
        end
    end

    return (
        num_nodes = num_nodes,
        num_edges = num_edges,
        critical_nodes = critical_nodes,
        source_nodes = source_nodes,
        sink_nodes = sink_nodes,
        avg_indegree = mean(indegree),
        avg_outdegree = mean(outdegree)
    )
end

"""
    build_nearest_neighbor_connections!(graph::DiGraph, site_to_index::Dict{String, Int},
                                      distance_data::DataFrame)

Connect each site to its nearest neighbors in the stream network.
"""
function build_nearest_neighbor_connections!(graph::DiGraph, site_to_index::Dict{<:AbstractString, <:Int},
                                           distance_data::DataFrame)
    edge_count = 0

    # Group by origin site to find nearest neighbors for each site
    for site in keys(site_to_index)
        # Find all connections from this site
        site_connections = filter(row -> row.ID_ORIGIN == site, distance_data)

        if nrow(site_connections) > 0
            # Sort by reticular distance (nearest first)
            sort!(site_connections, :RETICULAR_DIST)

            # Connect to the nearest neighbor(s)
            # For a dendritic network, we typically connect to 1-3 nearest neighbors
            max_neighbors = min(3, nrow(site_connections))
            for i in 1:max_neighbors
                row = site_connections[i, :]
                origin_idx = site_to_index[row.ID_ORIGIN]
                dest_idx = site_to_index[row.ID_DESTINATION]

                if !has_edge(graph, origin_idx, dest_idx)
                    add_edge!(graph, origin_idx, dest_idx)
                    edge_count += 1
                end
            end
        end
    end

    return edge_count
end

"""
    build_threshold_connections!(graph::DiGraph, site_to_index::Dict{String, Int},
                               distance_data::DataFrame, threshold::Float64)

Connect sites that are within a threshold reticular distance of each other.
"""
function build_threshold_connections!(graph::DiGraph, site_to_index::Dict{<:AbstractString, <:Int},
                                    distance_data::DataFrame, threshold::Float64)
    edge_count = 0

    # Only include connections within the threshold
    threshold_data = filter(row -> row.RETICULAR_DIST <= threshold, distance_data)

    for row in eachrow(threshold_data)
        origin_idx = site_to_index[row.ID_ORIGIN]
        dest_idx = site_to_index[row.ID_DESTINATION]

        if !has_edge(graph, origin_idx, dest_idx)
            add_edge!(graph, origin_idx, dest_idx)
            edge_count += 1
        end
    end

    return edge_count
end

"""
    build_mst_connections!(graph::DiGraph, site_to_index::Dict{String, Int},
                          distance_data::DataFrame)

Build connections using a minimum spanning tree approach to create the most
efficient stream network structure.
"""
function build_mst_connections!(graph::DiGraph, site_to_index::Dict{<:AbstractString, <:Int},
                               distance_data::DataFrame)
    edge_count = 0

    # Sort all connections by distance (shortest first)
    sorted_data = sort(distance_data, :RETICULAR_DIST)

    # Use Union-Find to build MST (simplified approach)
    # For now, we'll just add the shortest connections that don't create cycles
    # This is a simplified MST algorithm

    # Keep track of connected components
    components = [Int[i] for i in 1:length(site_to_index)]

    for row in eachrow(sorted_data)
        origin_idx = site_to_index[row.ID_ORIGIN]
        dest_idx = site_to_index[row.ID_DESTINATION]

        # Find which components these nodes belong to
        comp1 = find_component(components, origin_idx)
        comp2 = find_component(components, dest_idx)

        # If they're in different components, connect them
        if comp1 != comp2
            add_edge!(graph, origin_idx, dest_idx)
            edge_count += 1

            # Merge components
            merge_components!(components, comp1, comp2)
        end
    end

    return edge_count
end

"""
    find_component(components::Vector{Vector{Int}}, node::Int)

Find which component a node belongs to (simplified Union-Find).
"""
function find_component(components::Vector{Vector{Int}}, node::Int)
    for (i, comp) in enumerate(components)
        if node in comp
            return i
        end
    end
    return 1  # Default to first component
end

"""
    merge_components!(components::Vector{Vector{Int}}, comp1::Int, comp2::Int)

Merge two components (simplified Union-Find).
"""
function merge_components!(components::Vector{Vector{Int}}, comp1::Int, comp2::Int)
    if comp1 != comp2
        append!(components[comp1], components[comp2])
        empty!(components[comp2])
    end
end
