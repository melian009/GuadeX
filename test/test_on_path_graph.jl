# =============================================================================
# C4: dendritic dispersal graph reconstructed from the network-distance matrix.
#
# The preferred route (snap sites to SW_Line_4C_.shp and follow flow direction)
# is not feasible with the in-repo layer: it has no flow-direction attribute and
# no elevation, so the fallback on-path rule is used.  These tests pin the
# algorithm's topology, determinism, on-path consistency and reported lengths.
# =============================================================================

using TOML
using LinearAlgebra

# `SimulationParameters` is defined by test_parameters.jl (included earlier in
# runtests.jl); include it only if this file is run standalone.
if !isdefined(Main, :SimulationParameters)
    include(joinpath(@__DIR__, "..", "parameters.jl"))
end

"""
    _full_network_mst_km(distance_file, sites) -> Float64

Weight in km of the minimum spanning tree over every pair of `sites`, using the
reticular distances in `distance_file`.  Used only to check the reported MST
reference value (8,783 km); it is not a correctness criterion.
"""
function _full_network_mst_km(distance_file::String, sites::Vector{String})
    site_to_idx = Dict(s => i for (i, s) in enumerate(sites))
    n = length(sites)
    edges = Tuple{Float64, Int, Int}[]
    reader = CSV.File(distance_file; delim=';')
    for row in reader
        i = get(site_to_idx, row.ID_ORIGIN, 0)
        j = get(site_to_idx, row.ID_DESTINATION, 0)
        (i == 0 || j == 0 || i == j) && continue
        w = Float64(row.RETICULAR_DIST)
        (isfinite(w) && w > 0) || continue
        i < j && push!(edges, (w, i, j))
    end
    sort!(edges, by = e -> (e[1], e[2], e[3]))

    uf = collect(1:n)
    function find_root(x)
        while uf[x] != x
            uf[x] = uf[uf[x]]
            x = uf[x]
        end
        return x
    end

    total = 0.0
    for (w, a, b) in edges
        ra, rb = find_root(a), find_root(b)
        if ra != rb
            uf[ra] = rb
            total += w
        end
    end
    return total / 1000.0
end

@testset "C4 on-path river-network graph" begin

    site_df = Guadex.load_site_data(CONNECTIVITY_FILE, ENVIRONMENTAL_FILE)
    sites = String.(site_df.CODIGO)
    n = length(sites)
    d = Dict{String, Float64}(string(r.CODIGO) => Float64(coalesce(r."Dist.Guadalq.(m)", 0.0))
                              for r in eachrow(site_df))
    elevation = Dict{String, Float64}(string(r.CODIGO) => Float64(coalesce(r.ALTITUD, 500.0))
                                      for r in eachrow(site_df))

    result = Guadex.build_on_path_distance_matrix(DISTANCE_FILE, sites, d, elevation)
    diags = result.diagnostics

    @testset "topology: n-1 edges, connected, acyclic" begin
        @test diags.n_edges == n - 1
        @test diags.n_parent_edges + diags.n_root_join_edges == n - 1
        @test diags.n_roots == diags.n_root_join_edges + 1
        @test diags.nnz == 2 * (n - 1)          # both directions of every tree edge
        @test size(result.distance_matrix) == (n, n)
        @test all(iszero, LinearAlgebra.diag(result.distance_matrix))
        @test issymmetric(result.distance_matrix)

        # Union-find: every edge must join two different components (no cycle)
        # and all n sites must end up in a single component (spanning tree).
        uf = collect(1:n)
        find_root = function (x)
            while uf[x] != x
                uf[x] = uf[uf[x]]
                x = uf[x]
            end
            return x
        end
        cycles = 0
        for (u, v, _) in result.tree_edges
            ru, rv = find_root(u), find_root(v)
            if ru == rv
                cycles += 1
            else
                uf[ru] = rv
            end
        end
        @test cycles == 0
        @test length(Set(find_root(i) for i in 1:n)) == 1
    end

    @testset "on-path condition holds for every parent edge" begin
        seen_parent_edges = 0
        for (u, v, w) in result.tree_edges
            if result.parent[u] == v                 # parent (flow-path) edge
                seen_parent_edges += 1
                @test abs(w - (d[sites[u]] - d[sites[v]])) <= diags.rtol * w + diags.atol
            else                                      # root-join edge
                @test result.parent[u] == 0
            end
            # Flow direction: the source is never nearer to the mouth than the
            # sink (root joins with equal distance are allowed).
            @test d[sites[u]] + 1e-9 >= d[sites[v]]
        end
        @test seen_parent_edges == diags.n_parent_edges
    end

    @testset "elevation inversions are few and listed" begin
        inv = diags.parent_elevation_inversions
        @test length(inv) <= 10
        @test length(inv) == 4                     # pinned: 4 vs 238 in the old graph
        for (up_site, down_site, up_elev, down_elev) in inv
            @test up_elev < down_elev
            @test d[up_site] > d[down_site]
        end
        @test length(inv) < diags.n_parent_edges
    end

    @testset "deterministic across repeated runs" begin
        again = Guadex.build_on_path_distance_matrix(DISTANCE_FILE, sites, d, elevation)
        @test again.tree_edges == result.tree_edges
        @test again.parent == result.parent
        @test again.roots == result.roots
        @test again.diagnostics.total_length_km == diags.total_length_km
    end

    @testset "reported lengths (new tree vs old graph vs MST)" begin
        # Reviewer's prototype: 8,637 km parent edges + 5,235 km root MST.
        @test isapprox(diags.parent_length_km, 8636.94; atol=10.0)
        @test isapprox(diags.root_join_length_km, 5234.52; atol=10.0)
        @test isapprox(diags.total_length_km, 13871.46; atol=20.0)
        @test diags.total_length_km > diags.parent_length_km

        # Independent full-network minimum spanning tree on the same distances
        # (the review's 8,783 km).  Reported for comparison only; correctness is
        # judged by the topology/on-path tests above, not by proximity to the MST.
        mst_km = _full_network_mst_km(DISTANCE_FILE, sites)
        @test isapprox(mst_km, 8783.49; atol=10.0)
        @test diags.total_length_km > mst_km          # a spanning tree cannot beat the MST
    end

    @testset "legacy builder regression (unchanged default)" begin
        sc = Dict{String, String}(string(r.CODIGO) => string(r.CODIGO_S) for r in eachrow(site_df))
        legacy = Guadex.build_distance_matrix(DISTANCE_FILE, sites, sc, d, elevation)
        I, J, V = findnz(legacy)
        @test length(V) == 1548                       # both directions of 774 links
        pairs = Set{Tuple{Int, Int}}()
        for (i, j) in zip(I, J)
            push!(pairs, minmax(i, j))
        end
        @test length(pairs) == 774
        # 31,215 km in the review; the legacy outlet ties depend on Dict order, so
        # allow a small band while still separating it from the 8,783 km MST.
        @test isapprox(sum(V) / 2 / 1000, 31206.29; atol=100.0)
    end

    @testset "tie-breaking is deterministic (synthetic ties)" begin
        # Two downstream candidates at exactly the same network distance: the
        # smaller site code must win regardless of file/row order.
        mktempdir() do tmp
            path = joinpath(tmp, "dist.csv")
            rows = [
                "ID_ORIGIN;ID_DESTINATION;RETICULAR_DIST",
                "a;b;1000", "b;a;1000",
                "a;c;1000", "c;a;1000",
                "b;d;500", "d;b;500",
                "c;d;500", "d;c;500",
            ]
            write(path, join(rows, "\n") * "\n")
            dag = Dict("a" => 2000.0, "b" => 1000.0, "c" => 1000.0, "d" => 0.0)
            elev = Dict("a" => 300.0, "b" => 250.0, "c" => 200.0, "d" => 100.0)
            r1 = Guadex.build_on_path_distance_matrix(path, ["a", "b", "c", "d"], dag, elev)
            r2 = Guadex.build_on_path_distance_matrix(path, ["a", "b", "c", "d"], dag, elev)
            @test r1.tree_edges == r2.tree_edges
            # a -> b (1000) because b < c on ties; then the roots are MST-joined.
            @test (1, 2, 1000.0) in r1.tree_edges
            @test r1.diagnostics.n_edges == 3
        end
    end

    @testset "prepare_ode_data integration selects the method" begin
        data = Guadex.prepare_ode_data(connectivity_method=:on_path, upstream_cost=0.01)
        @test data.connectivity_method == :on_path
        @test data.connectivity_diagnostics !== nothing
        @test data.connectivity_diagnostics.n_edges == length(data.sites) - 1
        @test size(data.distance_matrix) == (length(data.sites), length(data.sites))
        # E12 interface: directed edge list + downstream-parent index.
        @test length(data.connectivity_tree_edges) == length(data.sites) - 1
        @test length(data.connectivity_parent) == length(data.sites)
        @test count(==(0), data.connectivity_parent) == data.connectivity_diagnostics.n_roots
    end

    @testset "corrected config selects the on-path method" begin
        corrected = TOML.parsefile(joinpath(@__DIR__, "..",
            "parameters_climate_scenarios_corrected.toml"))
        @test SimulationParameters.connectivity_method(corrected) == :on_path
        legacy_cfg = TOML.parsefile(joinpath(@__DIR__, "..", "parameters.toml"))
        @test SimulationParameters.connectivity_method(legacy_cfg) == :legacy
        k1x = TOML.parsefile(joinpath(@__DIR__, "..", "legacy",
            "parameters_climate_scenarios_k1x_burnin.toml"))
        @test SimulationParameters.connectivity_method(k1x) == :legacy
        # Unknown values fail loudly rather than silently falling back.
        @test_throws Exception SimulationParameters.connectivity_method(
            Dict{String, Any}("connectivity" => Dict{String, Any}("method" => "nonsense")))
    end

    @testset "E12 obstacles snap to the corrected tree edges" begin
        site_df = DataFrame(CODIGO=["a", "b", "c"], UTMX=[0.0, 100.0, 200.0],
            UTMY=[0.0, 0.0, 0.0])
        sites = ["a", "b", "c"]
        distances = sparse([2, 1, 3, 2], [1, 2, 2, 3],
            [100.0, 100.0, 100.0, 100.0], 3, 3)
        # Tree orientation and elevations deliberately disagree: the tree says
        # a is upstream of b, while b is the higher site.
        tree = [(1, 2, 100.0), (2, 3, 100.0)]
        elevations = [100.0, 200.0, 50.0]
        obs = DataFrame(ID_CLAVE_MGM=["o1", "o2"], coord_x_m=[50.0, 60.0],
            coord_y_m=[1.0, 1.0], ESTADO=["EX", "EX"])

        r = Guadex.build_obstacle_passability_matrix(obs, site_df, sites, distances,
            elevations; matching_tolerance=10.0, tree_edges=tree)
        @test r.metadata.use_tree_edges
        @test r.metadata.n_matched == 2
        # Barriers accumulate on the tree link a-b (0.1 * 0.1 / 0.5 * 0.5).
        @test r.passability[1, 2] ≈ 0.01     # a upstream -> b downstream
        @test r.passability[2, 1] ≈ 0.25
        # The corrected tree orientation, not the elevation rule, sets the
        # restricted direction (destination is the upstream site a).
        @test all(r.diagnostics.restricted_destination .== "a")
        @test all(r.diagnostics.restricted_origin .== "b")

        # Every matched link is a corrected tree edge.
        tree_pairs = Set((min(u, v), max(u, v)) for (u, v, _) in tree)
        for row in eachrow(r.diagnostics)
            row.matched || continue
            i = findfirst(==(row.edge_from), sites)
            j = findfirst(==(row.edge_to), sites)
            @test (min(i, j), max(i, j)) in tree_pairs
        end
    end

    @testset "E12 prepare_ode_data snaps obstacles on the on-path graph" begin
        data = Guadex.prepare_ode_data(connectivity_method=:on_path,
            obstacles_file=OBSTACLES_FILE, obstacle_mode=:overlay,
            obstacle_matching_tolerance=2000.0, upstream_cost=0.01)
        @test data.connectivity_method == :on_path
        @test data.obstacle_metadata !== nothing
        @test data.obstacle_metadata.use_tree_edges
        @test nrow(data.obstacle_mapping_diagnostics) == 1658
        # Matched links must be corrected tree edges, never artificial chords.
        tree_pairs = Set((min(u, v), max(u, v)) for (u, v, _) in data.connectivity_tree_edges)
        site_index = Dict(s => i for (i, s) in enumerate(data.sites))
        for row in eachrow(data.obstacle_mapping_diagnostics)
            row.matched || continue
            i = site_index[row.edge_from]
            j = site_index[row.edge_to]
            @test (min(i, j), max(i, j)) in tree_pairs
        end
        @test data.obstacle_metadata.n_applied > 0
        @test data.obstacle_metadata.n_status_nonoperational == 185
        @test data.obstacle_metadata.n_excluded_out_of_basin == 5
    end
end
