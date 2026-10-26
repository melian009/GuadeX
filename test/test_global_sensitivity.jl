using Test
using Guadex
using Statistics
using SparseArrays
using LinearAlgebra

# =============================================================================
# E7 — global sensitivity toolkit: Latin hypercube, Saltelli design, Sobol
# estimators validated against the analytic Ishigami function.
# =============================================================================

@testset "Global sensitivity toolkit (E7)" begin

    @testset "Latin hypercube: structure, marginals, determinism" begin
        n, d = 200, 3
        U = latin_hypercube(n, d; seed=42)
        @test size(U) == (n, d)
        @test all(0.0 .<= U .< 1.0)

        # One sample per equiprobable stratum in every dimension (exact LHS).
        for j in 1:d
            s = sort(U[:, j])
            @test all(floor(Int, s[k] * n) == k - 1 for k in 1:n)
        end

        # Determinism: same seed identical, different seed different.
        @test latin_hypercube(n, d; seed=42) == U
        @test latin_hypercube(n, d; seed=7) != U

        # Marginal range mapping onto the parameter box.
        lo = [-1.0, 0.0, 10.0]
        hi = [1.0, 5.0, 20.0]
        S = lhs_samples(n, lo, hi; seed=1)
        @test size(S) == (n, d)
        @test all(lo[j] <= S[i, j] <= hi[j] for i in 1:n, j in 1:d)

        # A bad box is rejected rather than silently sampled.
        @test_throws ErrorException latin_hypercube(0, 2; seed=1)
        @test_throws ErrorException rescale_columns(ones(3, 2), [0.0], [1.0, 2.0])
    end

    @testset "Saltelli A/B/AB/BA design" begin
        design = saltelli_design(16, 3; seed=1)
        @test size(design.A) == (16, 3)
        @test n_model_runs(design) == 16 * (2 * 3 + 2)

        # AB[i] differs from A only in column i (taken from B); BA symmetric.
        for i in 1:3
            @test design.AB[i][:, i] == design.B[:, i]
            @test design.BA[i][:, i] == design.A[:, i]
            for j in 1:3
                j == i && continue
                @test design.AB[i][:, j] == design.A[:, j]
                @test design.BA[i][:, j] == design.B[:, j]
            end
        end

        M = saltelli_sample_matrix(design)
        @test size(M) == (16 * (2 * 3 + 2), 3)
        lo = fill(-2.0, 3); hi = fill(2.0, 3)
        Ms = saltelli_sample_matrix(design; lower=lo, upper=hi)
        @test all(lo[j] <= Ms[i, j] <= hi[j] for i in 1:size(Ms, 1), j in 1:3)
    end

    @testset "Ishigami: S1/ST match the analytic values" begin
        ishigami(x) = sin(x[1]) + 7 * sin(x[2])^2 + 0.1 * x[3]^4 * sin(x[1])
        lower = fill(-pi, 3); upper = fill(pi, 3)

        design = saltelli_design(2^13, 3; seed=20261001)
        result = sobol_analyze(ishigami, design; lower=lower, upper=upper,
            bootstrap=300, seed=7)

        # Analytic variance and indices for a = 7, b = 0.1.
        V = 1 / 2 + 49 / 8 + 0.1 * pi^4 / 5 + 0.01 * pi^8 / 18
        S1_analytic = [(1 / 2 + 0.1 * pi^4 / 5 + 0.01 * pi^8 / 50) / V, (49 / 8) / V, 0.0]
        ST_analytic = [0.5576, 0.4424, 0.2437]

        @test result.n == 2^13
        @test result.n_model_runs == 2^13 * 8
        @test result.variance[1] ≈ V rtol = 0.05

        for i in 1:3
            @test result.S1[i, 1] ≈ S1_analytic[i] atol = 0.05
            @test result.ST[i, 1] ≈ ST_analytic[i] atol = 0.05
        end

        # Bootstrap intervals are ordered, finite and informative.
        for i in 1:3
            @test result.S1_ci[i, 1, 1] <= result.S1_ci[i, 1, 2]
            @test result.ST_ci[i, 1, 1] <= result.ST_ci[i, 1, 2]
            @test isfinite(result.S1_ci[i, 1, 1]) && isfinite(result.ST_ci[i, 1, 2])
            @test (result.ST_ci[i, 1, 2] - result.ST_ci[i, 1, 1]) < 0.5
        end
        # The dominant second parameter must bracket its estimate.
        @test result.S1_ci[2, 1, 1] <= result.S1[2, 1] <= result.S1_ci[2, 1, 2]
        # Non-degenerate widths: a permutation bootstrap would collapse these to
        # ~0 because resampling without replacement preserves every sample mean.
        @test (result.S1_ci[1, 1, 2] - result.S1_ci[1, 1, 1]) > 1e-4
        @test (result.ST_ci[3, 1, 2] - result.ST_ci[3, 1, 1]) > 1e-4

        # The bootstrap index generator must draw WITH replacement.
        boot = Guadex.bootstrap_indices(Guadex.GlobalSensitivityRNG(1), 100)
        @test length(boot) == 100
        @test all(1 .<= boot .<= 100)
        @test length(unique(boot)) < 100


        # Determinism of the full driver under a fixed seed.
        again = sobol_analyze(ishigami, saltelli_design(2^10, 3; seed=3);
            lower=lower, upper=upper, bootstrap=40, seed=11)
        repeat = sobol_analyze(ishigami, saltelli_design(2^10, 3; seed=3);
            lower=lower, upper=upper, bootstrap=40, seed=11)
        @test again.S1 == repeat.S1
        @test again.ST_ci == repeat.ST_ci
    end

    @testset "basin_response_metrics" begin
        species = ["A", "B", "C"]
        # Flattened states are `vec(sites × species)` (column-major), so build
        # the 2×3 matrices and vectorise them to keep site/species alignment.
        baseline = vec([1.0 1.0 1.0; 0.0 0.0 0.0])
        final = vec([1.0 0.0 0.0; 0.0 0.0 0.0])
        m = basin_response_metrics(final, baseline, species; native=["A", "B"], threshold=0.1)
        @test m.native_biomass ≈ 0.5
        @test m.native_richness ≈ 0.5
        @test m.realised_richness_loss ≈ 0.5
        # No native codes present -> explicit error, not a silent zero.
        @test_throws ErrorException basin_response_metrics(final, baseline, species;
            native=["ZZ"], threshold=0.1)
    end

    @testset "parameter transforms" begin
        p = MetacommunityParams(2, 2, [0.0 -0.3; -0.2 0.0],
            sparse([1, 2], [2, 1], [0.1, 0.1], 2, 2),
            [1.0, 1.0], [0.1 0.1; 0.1 0.1], [10.0, 10.0], [0.5, 0.5],
            [15.0, 15.0], [3.0, 3.0], [100.0, 100.0], [-Inf, -Inf], [Inf, Inf], 0.0)
        @test scale_intrinsic_growth_rate(p, 2.0).intrinsic_growth_rates ≈
              p.intrinsic_growth_rates .* 2.0
        @test scale_dispersal(p, 3.0).dispersal_matrix ≈ p.dispersal_matrix .* 3.0
        @test_throws ErrorException scale_intrinsic_growth_rate(p, 0.0)
        @test_throws ErrorException scale_dispersal(p, -1.0)
    end
end
