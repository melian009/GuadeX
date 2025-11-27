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
                           connection_method::Symbol=:nearest_neighbors,
                           connectivity_file::String="data/ConnectivityUTM.csv")
    # Read the distance matrix
    println("Reading distance matrix from: $distance_file")
    distance_data = CSV.read(distance_file, DataFrame)

    # Load subcatchment info from ConnectivityUTM.csv
    println("Reading subcatchment info from: $connectivity_file")
    connectivity_df = CSV.read(connectivity_file, DataFrame)
    site_to_subcatchment = Dict(row.CODIGO => row.CODIGO_S for row in eachrow(connectivity_df))

    # Filter site "1.30.20" from both datasets because it lacks reticular distance info
    distance_data = filter(row -> row.ID_ORIGIN != "1.30.20" && row.ID_DESTINATION != "1.30.20", distance_data)
    connectivity_df = filter(row -> row.CODIGO != "1.30.20", connectivity_df)

    # Filter out self-connections and distances greater than max_distance
    filtered_data = filter(row -> row.ID_ORIGIN != row.ID_DESTINATION &&
                                   0 < row.RETICULAR_DIST <= max_distance &&
                                   haskey(site_to_subcatchment, row.ID_ORIGIN),
                           distance_data)

    println("Found $(nrow(filtered_data)) valid connections out of $(nrow(distance_data)) total (with subcatchment info)")

    # Get unique site codes
    unique_sites = unique(vcat(filtered_data.ID_ORIGIN, filtered_data.ID_DESTINATION))
    println("Building graph with $(length(unique_sites)) sites")

    # Create mapping from site code to node index
    site_to_index = Dict(String(site) => i for (i, site) in enumerate(unique_sites))

    # Create directed graph
    num_sites = length(unique_sites)
    graph = DiGraph(num_sites)

    # Add edges based on stream connections
    # We need to determine direction based on elevation or other criteria
    # For now, we'll use a simple approach: assume flow from higher to lower elevation
    # We'll need to load elevation data for this

    # Build connections based on the specified method
    if connection_method == :nearest_neighbors
        edge_count = build_nearest_neighbor_connections!(graph, site_to_index, connectivity_df)
    elseif connection_method == :threshold_distance
        # edge_count = build_threshold_connections!(graph, site_to_index, filtered_data, max_distance, site_to_subcatchment)
        println("Method not yet tested: threshold_distance")
    elseif connection_method == :minimum_spanning_tree
        # edge_count = build_mst_connections!(graph, site_to_index, filtered_data, site_to_subcatchment)
        println("Method not yet tested: minimum_spanning_tree")
    elseif connection_method == :all_connections
        println("Method: all_connections not yet tested")
        # edge_count = 0
        # for row in eachrow(filtered_data)
        #     origin = row.ID_ORIGIN
        #     dest = row.ID_DESTINATION
        #     # Only connect if both sites are in the same subcatchment
        #     if site_to_subcatchment[origin] == site_to_subcatchment[dest]
        #         origin_idx = site_to_index[origin]
        #         dest_idx = site_to_index[dest]
        #         if !has_edge(graph, origin_idx, dest_idx)
        #             add_edge!(graph, origin_idx, dest_idx)
        #             edge_count += 1
        #         end
        #         if !has_edge(graph, dest_idx, origin_idx)
        #             add_edge!(graph, dest_idx, origin_idx)
        #             edge_count += 1
        #         end
        #     end
        # end
    else
        error("Unknown connection method: $connection_method. Use :nearest_neighbors, :threshold_distance, :minimum_spanning_tree, or :all_connections")
    end

    println("Added $edge_count edges to the graph using $connection_method method")

    return graph, site_to_index, filtered_data
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
function build_nearest_neighbor_connections!(graph::DiGraph, site_to_index::Dict{<:AbstractString, <:Int}, connectivity_df)
    edge_count = 0

    # Get all unique subcatchments
    subcatchments = unique(connectivity_df.CODIGO_S)
    # Add subcatchment nodes to graph if not present
    subcatchment_indices = Dict{Float64, Int}()
    next_idx = maximum(values(site_to_index)) + 1
    for sc in subcatchments
        if !haskey(site_to_index, string(sc))
            site_to_index[string(sc)] = next_idx
            subcatchment_indices[sc] = next_idx
            add_vertex!(graph)
            next_idx += 1
        else
            subcatchment_indices[sc] = site_to_index[string(sc)]
        end
    end

    # Connect sites within each subcatchment based on Dist.Guadalq.(m)
    sc_proxy_altitude = Dict{Float64, Float64}()  # To store proxy altitude for each subcatchment
    for sc in subcatchments
        # Get sites in this subcatchment
        sites_in_sc = filter(row -> row.CODIGO_S == sc, connectivity_df)
        # Sort sites by Dist.Guadalq.(m) (ascending: closest to river first)
        sorted_sites = sort(sites_in_sc, "Dist.Guadalq.(m)")
        # Connect sequentially (upstream to downstream)
        for i in 1:(nrow(sorted_sites)-1)
            site_up = sorted_sites.CODIGO[i+1]
            site_down = sorted_sites.CODIGO[i]
            idx_up = site_to_index[site_up]
            idx_down = site_to_index[site_down]
            if !has_edge(graph, idx_up, idx_down)
                add_edge!(graph, idx_up, idx_down)
                edge_count += 1
            end
        end
        # Connect the closest site to the river to the subcatchment node
        site_closest = sorted_sites.CODIGO[1]
        idx_closest = site_to_index[site_closest]
        idx_sc = subcatchment_indices[sc]
        if !has_edge(graph, idx_closest, idx_sc)
            add_edge!(graph, idx_closest, idx_sc)
            edge_count += 1
        end
        # Save proxy altitude for subcatchment
        sc_proxy_altitude[sc] = sorted_sites.ALTITUD[1]
    end

    # Connect subcatchment nodes sequentially by proxy altitude (highest to lowest)
    sorted_sc = sort(subcatchments, by=sc -> sc_proxy_altitude[sc], rev=true)
    for i in 1:(length(sorted_sc)-1)
        sc_up = sorted_sc[i]
        sc_down = sorted_sc[i+1]
        idx_up = subcatchment_indices[sc_up]
        idx_down = subcatchment_indices[sc_down]
        if !has_edge(graph, idx_up, idx_down)
            add_edge!(graph, idx_up, idx_down)
            edge_count += 1
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
                                    distance_data::DataFrame, threshold::Float64, site_to_subcatchment::Dict)
    edge_count = 0
    threshold_data = filter(row -> 0 < row.RETICULAR_DIST <= threshold &&
                                   haskey(site_to_subcatchment, row.ID_ORIGIN) &&
                                   haskey(site_to_subcatchment, row.ID_DESTINATION) &&
                                   site_to_subcatchment[row.ID_ORIGIN] == site_to_subcatchment[row.ID_DESTINATION],
                                   distance_data)
    for row in eachrow(threshold_data)
        origin = row.ID_ORIGIN
        dest = row.ID_DESTINATION
        # Only connect if both sites are in the same subcatchment
        origin_idx = site_to_index[origin]
        dest_idx = site_to_index[dest]
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
                               distance_data::DataFrame, site_to_subcatchment::Dict)
    edge_count = 0
    sorted_data = sort(distance_data, :RETICULAR_DIST)
    components = [Int[i] for i in 1:length(site_to_index)]
    for row in eachrow(sorted_data)
        origin = row.ID_ORIGIN
        dest = row.ID_DESTINATION
        # Only connect if both sites are in the same subcatchment
        if site_to_subcatchment[origin] == site_to_subcatchment[dest]
            origin_idx = site_to_index[origin]
            dest_idx = site_to_index[dest]
            comp1 = find_component(components, origin_idx)
            comp2 = find_component(components, dest_idx)
            if comp1 != comp2
                add_edge!(graph, origin_idx, dest_idx)
                edge_count += 1
                merge_components!(components, comp1, comp2)
            end
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
