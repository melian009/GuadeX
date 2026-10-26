@testset "Graph Construction" begin
    @testset "Basic Stream Graph" begin
        graph, site_to_index_str, distance_data = build_stream_graph(
            DISTANCE_FILE, max_distance=10000.0, connectivity_file=CONNECTIVITY_FILE
        )
        stats = Guadex.get_graph_statistics(graph)
        @test stats.num_nodes > 0
        @test stats.num_edges > 0
        # Density is the edge count over the n·(n−1) possible directed edges, so
        # it must be a genuine fraction, not merely non-negative.
        @test 0 < stats.density <= 1
        @test stats.num_components >= 1
        @test stats.largest_component_size > 0
        @test stats.largest_component_size <= stats.num_nodes

        # Topology: the builder never adds a self-loop and every edge joins two
        # distinct known nodes.
        @test all(e -> src(e) != dst(e), edges(graph))
        @test all(e -> 1 <= src(e) <= stats.num_nodes && 1 <= dst(e) <= stats.num_nodes,
            edges(graph))

        sample_sites = collect(keys(site_to_index_str))[1:min(5, length(site_to_index_str))]
        for site in sample_sites
            upstream = Guadex.find_upstream_sites(graph, site_to_index_str, site)
            downstream = Guadex.find_downstream_sites(graph, site_to_index_str, site)
            # A site is never its own up/down-stream neighbour and the two sets
            # can never contain every node (the target itself is excluded).
            @test !(site in upstream)
            @test !(site in downstream)
            @test length(upstream) <= stats.num_nodes - 1
            @test length(downstream) <= stats.num_nodes - 1
            # Every reported upstream site must actually reach the target, and
            # every downstream site must be reachable from it.
            t = site_to_index_str[site]
            @test all(u -> has_path(graph, site_to_index_str[u], t), upstream)
            @test all(d -> has_path(graph, t, site_to_index_str[d]), downstream)
        end
    end
end
