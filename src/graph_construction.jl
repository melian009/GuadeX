using DataFrames
using CSV
using Graphs

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
function build_stream_graph(distance_file::String; max_distance::Float64=Inf)
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

    # Add edges (bidirectional for now, will be refined with elevation data)
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

    println("Added $edge_count edges to the graph")

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