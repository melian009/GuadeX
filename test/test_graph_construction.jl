@testset "Graph Construction" begin
    @testset "Basic Stream Graph" begin
        graph, site_to_index_str, distance_data = build_stream_graph(
            DISTANCE_FILE, max_distance=10000.0, connectivity_file=CONNECTIVITY_FILE
        )
        stats = Guadex.get_graph_statistics(graph)
        @test stats.num_nodes > 0
        @test stats.num_edges > 0
        @test stats.density >= 0
        @test stats.num_components >= 1
        @test stats.largest_component_size > 0
        sample_sites = collect(keys(site_to_index_str))[1:min(5, length(site_to_index_str))]
        for site in sample_sites
            upstream = Guadex.find_upstream_sites(graph, site_to_index_str, site)
            downstream = Guadex.find_downstream_sites(graph, site_to_index_str, site)
            @test length(upstream) >= 0
            @test length(downstream) >= 0
        end
    end
end
