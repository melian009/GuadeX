using DifferentialEquations
using LinearAlgebra

# =============================================================================
# 1. Unit tests for gaussian_thermal_filter
# =============================================================================
@testset "gaussian_thermal_filter" begin
    # At optimum: filter = 1.0 regardless of sigma
    @test Guadex.gaussian_thermal_filter(20.0, 20.0, 3.0) ≈ 1.0
    @test Guadex.gaussian_thermal_filter(15.0, 15.0, 10.0) ≈ 1.0

    # At T = opt ± sigma: filter = exp(-1/2) ≈ 0.60653
    @test Guadex.gaussian_thermal_filter(25.0, 20.0, 5.0) ≈ exp(-0.5)  rtol=1e-6
    @test Guadex.gaussian_thermal_filter(15.0, 20.0, 5.0) ≈ exp(-0.5)  rtol=1e-6

    # At T = opt ± 2*sigma: filter = exp(-2) ≈ 0.13534
    @test Guadex.gaussian_thermal_filter(30.0, 20.0, 5.0) ≈ exp(-2.0) rtol=1e-6
    @test Guadex.gaussian_thermal_filter(10.0, 20.0, 5.0) ≈ exp(-2.0) rtol=1e-6

    # Symmetry
    @test Guadex.gaussian_thermal_filter(22.0, 20.0, 3.0) ≈
          Guadex.gaussian_thermal_filter(18.0, 20.0, 3.0)

    # Larger sigma → slower decay (same deviation gives higher filter)
    @test Guadex.gaussian_thermal_filter(23.0, 20.0, 6.0) >
          Guadex.gaussian_thermal_filter(23.0, 20.0, 3.0)
end

# =============================================================================
# 2. precompute_dispersal_matrix tests
# =============================================================================
@testset "precompute_dispersal_matrix" begin
    n_sites = 2
    distances = [0.0 1000.0; 1000.0 0.0]  # 1 km apart
    elevations = [10.0, 50.0]
    c = 0.01
    dams = ones(n_sites, n_sites)
    D_daily = 10.0 / 365.0  # median daily dispersal

    @testset "5-arg version (explicit D_daily)" begin
        mat = Guadex.precompute_dispersal_matrix(n_sites, distances, elevations, c, dams, D_daily)

        @test size(mat) == (2, 2)
        @test mat[1, 1] == 0.0  # diagonal = 0 (no self-dispersal)
        @test mat[2, 2] == 0.0

        # m_21: from site1 (e=10) to site2 (e=50), upstream
        # x = 1/(1 + c*(50-10)) = 1/1.4, d_km = 1.0
        expected_21 = D_daily * (1.0 / (1.0 + 0.01 * 40.0)) / 1.0
        @test mat[2, 1] ≈ expected_21

        # m_12: from site2 (e=50) to site1 (e=10), downstream
        # x = 1.0, d_km = 1.0
        expected_12 = D_daily * 1.0 / 1.0
        @test mat[1, 2] ≈ expected_12

        # Upstream is harder than downstream
        @test mat[2, 1] < mat[1, 2]
    end

    @testset "5-arg version (default D_daily)" begin
        # Test backwards-compatible 5-arg call uses default D_daily
        mat = Guadex.precompute_dispersal_matrix(n_sites, distances, elevations, c, dams)
        @test size(mat) == (2, 2)
        expected_12 = (10.0/365.0) * 1.0 / 1.0
        @test mat[1, 2] ≈ expected_12
    end

    @testset "6-arg version (with species_codes)" begin
        species_codes = ["AB", "AH", "MC"]
        mat = Guadex.precompute_dispersal_matrix(n_sites, distances, elevations, c, dams, species_codes)

        @test size(mat) == (2, 2)
        @test mat[1, 1] == 0.0

        # Should compute D_daily from species codes: rates = [0.3, 0.2, 75.0], median = 0.3
        # D_daily = 0.3/365 ≈ 0.000822
        D_from_species = 0.3 / 365.0
        expected_12 = D_from_species * 1.0 / 1.0
        @test mat[1, 2] ≈ expected_12
    end

    @testset "minimum distance clamping" begin
        # Very short distance → clamped to 0.001 km
        short_dist = [0.0 0.5; 0.5 0.0]  # 0.5 m
        mat = Guadex.precompute_dispersal_matrix(n_sites, short_dist, elevations, c, dams, D_daily)
        expected = D_daily * 1.0 / max(0.5/1000.0, 0.001)
        @test mat[1, 2] ≈ D_daily * 1.0 / 0.001
    end

    @testset "dam passability reduces dispersal" begin
        dams_blocked = ones(n_sites, n_sites)
        dams_blocked[1, 2] = 0.1  # dam from site2 to site1

        mat = Guadex.precompute_dispersal_matrix(n_sites, distances, elevations, c, dams_blocked, D_daily)
        expected_blocked = D_daily * 1.0 * 0.1 / 1.0
        @test mat[1, 2] ≈ expected_blocked
    end
end

# =============================================================================
# 3. MetacommunityParams construction
# =============================================================================
@testset "MetacommunityParams" begin
    n_sites, n_species = 3, 2
    distances = [0.0 1000.0 2000.0; 1000.0 0.0 3000.0; 2000.0 3000.0 0.0]
    elevations = [10.0, 50.0, 30.0]
    interaction_matrix = [-0.1 -0.3; -0.5 -0.1]
    intrinsic_growth = [0.01 0.02; 0.015 0.025; 0.012 0.022]
    D_daily = 10.0 / 365.0

    dispersal_mat = Guadex.precompute_dispersal_matrix(
        n_sites, distances, elevations, 0.01, ones(n_sites, n_sites), D_daily
    )

    p = Guadex.MetacommunityParams(
        n_sites, n_species,
        interaction_matrix,
        dispersal_mat,
        [1.0, 2.0],  # dispersal_scaling
        intrinsic_growth,
        [15.0, 18.0, 12.0],  # temperatures
        [0.8, 1.0, 0.6],     # habitat_suitability
        [19.0, 15.0],         # thermal_optima
        [3.67, 3.0],          # thermal_sigmas
        [100.0, 200.0, 150.0] # carrying_capacity
    )

    @test p.n_sites == 3
    @test p.n_species == 2
    @test p.interaction_matrix == interaction_matrix
    @test p.dispersal_scaling == [1.0, 2.0]
    @test p.thermal_optima == [19.0, 15.0]
    @test p.thermal_sigmas == [3.67, 3.0]
    @test p.carrying_capacity == [100.0, 200.0, 150.0]
end

# =============================================================================
# 4. metacommunity_ode! unit tests
# =============================================================================
@testset "metacommunity_ode! correctness" begin
    function setup_2x2_params()
        n_sites, n_species = 2, 2
        distances = [0.0 1000.0; 1000.0 0.0]
        elevations = [10.0, 10.0]  # same elevation → no upstream cost
        interaction = [-0.1 0.0; 0.0 -0.1]
        growth = [0.01 0.01; 0.01 0.01]
        D_daily = 10.0 / 365.0

        dispersal = Guadex.precompute_dispersal_matrix(
            n_sites, distances, elevations, 0.01, ones(n_sites, n_sites), D_daily
        )

        Guadex.MetacommunityParams(
            n_sites, n_species,
            interaction, dispersal,
            [1.0, 1.0],     # dispersal_scaling
            growth,
            [20.0, 20.0],   # temperatures
            [1.0, 1.0],     # habitat_suitability
            [20.0, 20.0],   # thermal_optima (at optimum)
            [5.0, 5.0],     # thermal_sigmas
            [30.0, 30.0]    # carrying_capacity
        )
    end

    @testset "growth increases population below capacity" begin
        p = setup_2x2_params()
        # Use zero interactions + symmetric populations to isolate logistic growth
        p_isolated = Guadex.MetacommunityParams(
            p.n_sites, p.n_species,
            [0.0 0.0; 0.0 0.0],  # no interactions
            p.dispersal_matrix, p.dispersal_scaling,
            p.intrinsic_growth_rates, p.temperatures,
            p.habitat_suitability, p.thermal_optima,
            p.thermal_sigmas, p.carrying_capacity
        )
        # Symmetric initial conditions → no net dispersal
        u0 = [5.0 3.0; 5.0 3.0]
        du = zeros(2, 2)
        Guadex.metacommunity_ode!(du, u0, p_isolated, 0.0)

        @test size(du) == (2, 2)
        @test all(isfinite.(du))
        # At optimum temp, env_filter = 1.0, r_eff = 0.01
        # total at each site = 8 < K=30 → logistic_term > 0
        # With symmetric pops, net dispersal ≈ 0
        @test du[1, 1] > 0.0
        @test du[1, 2] > 0.0
    end

    @testset "logistic term becomes negative above capacity" begin
        p = setup_2x2_params()
        p_no_interact = Guadex.MetacommunityParams(
            p.n_sites, p.n_species,
            [0.0 0.0; 0.0 0.0],
            p.dispersal_matrix, p.dispersal_scaling,
            p.intrinsic_growth_rates, p.temperatures,
            p.habitat_suitability, p.thermal_optima,
            p.thermal_sigmas, p.carrying_capacity
        )
        # Symmetric → no net dispersal. Total = 35 > K = 30 at each site.
        u0 = [20.0 15.0; 20.0 15.0]
        du = zeros(2, 2)
        Guadex.metacommunity_ode!(du, u0, p_no_interact, 0.0)

        # logistic_term = clamp(1 - 35/30, -1, 2) = -0.167
        # growth = N * r_eff * logistic_term < 0
        @test du[1, 1] < 0.0
        @test du[1, 2] < 0.0
        @test du[2, 1] < 0.0
    end

    @testset "thermal filter reduces growth away from optimum" begin
        p_base = setup_2x2_params()
        # Both params have zero interactions; only temperature differs
        make_params(T) = Guadex.MetacommunityParams(
            p_base.n_sites, p_base.n_species,
            [0.0 0.0; 0.0 0.0],  # no interactions in both
            p_base.dispersal_matrix, p_base.dispersal_scaling,
            p_base.intrinsic_growth_rates,
            [T, T], p_base.habitat_suitability,
            p_base.thermal_optima, p_base.thermal_sigmas,
            p_base.carrying_capacity
        )
        p_optimal = make_params(20.0)
        p_cold = make_params(10.0)

        # Symmetric initial conditions → no net dispersal
        u0 = [10.0 8.0; 10.0 8.0]
        du_opt = zeros(2, 2)
        du_cold = zeros(2, 2)

        Guadex.metacommunity_ode!(du_opt, u0, p_optimal, 0.0)
        Guadex.metacommunity_ode!(du_cold, u0, p_cold, 0.0)

        # Growth should be lower with suboptimal temperature
        @test du_cold[1, 1] < du_opt[1, 1]
        @test du_cold[1, 2] < du_opt[1, 2]
    end

    @testset "dispersal moves individuals between sites" begin
        n_sites, n_species = 2, 1
        distances = [0.0 1000.0; 1000.0 0.0]
        elevations = [10.0, 10.0]
        D_daily = 10.0 / 365.0
        dispersal = Guadex.precompute_dispersal_matrix(
            n_sites, distances, elevations, 0.01, ones(n_sites, n_sites), D_daily
        )

        p_disp = Guadex.MetacommunityParams(
            n_sites, n_species,
            [0.0;;],  # 1x1 interaction matrix (zero)
            dispersal,
            [1.0],    # dispersal_scaling
            [0.0; 0.0;;],  # zero growth (isolate dispersal)
            [20.0, 20.0],
            [1.0, 1.0],
            [20.0],
            [5.0],
            [100.0, 100.0]
        )

        # Unequal populations → net dispersal from high to low
        u0 = [100.0; 0.0]
        du = zeros(2, 1)
        Guadex.metacommunity_ode!(du, u0, p_disp, 0.0)

        # Site1 loses individuals, site2 gains
        @test du[1] < 0.0  # emigration from site1
        @test du[2] > 0.0  # immigration to site2
    end

    @testset "dispersal_scaling affects per-species rates" begin
        n_sites, n_species = 2, 2
        distances = [0.0 1000.0; 1000.0 0.0]
        elevations = [10.0, 10.0]
        D_daily = 10.0 / 365.0
        dispersal = Guadex.precompute_dispersal_matrix(
            n_sites, distances, elevations, 0.01, ones(n_sites, n_sites), D_daily
        )

        p = Guadex.MetacommunityParams(
            n_sites, n_species,
            [0.0 0.0; 0.0 0.0],
            dispersal,
            [0.1, 10.0],  # sp1 slow, sp2 fast
            [0.0 0.0; 0.0 0.0],  # zero growth
            [20.0, 20.0],
            [1.0, 1.0],
            [20.0, 20.0],
            [5.0, 5.0],
            [100.0, 100.0]
        )

        u0 = [100.0 100.0; 0.0 0.0]
        du = zeros(2, 2)
        Guadex.metacommunity_ode!(du, u0, p, 0.0)

        # Fast disperser (sp2) has larger magnitude derivative
        @test abs(du[2, 2]) > abs(du[2, 1])
    end

    @testset "zero population → zero derivative (from growth)" begin
        n_sites, n_species = 1, 1
        p = Guadex.MetacommunityParams(
            1, 1, [0.0;;], sparse([1], [1], [0.0], 1, 1),
            [1.0], [0.1;;], [20.0], [1.0], [20.0], [5.0], [100.0]
        )
        du = zeros(1, 1)
        Guadex.metacommunity_ode!(du, [0.0;;], p, 0.0)
        @test du[1] == 0.0
    end

    @testset "ANNUAL_DISPERSAL_RATES constant" begin
        rates = Guadex.ANNUAL_DISPERSAL_RATES
        @test "AB" in keys(rates)
        @test "GH" in keys(rates)
        @test "ST" in keys(rates)
        @test rates["AB"] == 0.3
        @test rates["MC"] == 75.0
        @test rates["GH"] == 25.0

        # Verify median is 10.0 (for documentation consistency)
        vals = collect(values(rates))
        @test Statistics.median(vals) ≈ 10.0
    end
end

# =============================================================================
# 5. Full ODE integration tests
# =============================================================================
@testset "Full ODE integration" begin
    @testset "Population grows from small initial condition" begin
        n_sites, n_species = 1, 1
        distances = [0.0;;]
        elevations = [10.0]
        D_daily = 10.0 / 365.0
        dispersal = Guadex.precompute_dispersal_matrix(
            n_sites, distances, elevations, 0.01, ones(1, 1), D_daily
        )

        p = Guadex.MetacommunityParams(
            1, 1, [0.0;;], dispersal,
            [1.0], [0.05;;],  # growth rate 0.05/day
            [20.0], [1.0], [20.0], [5.0], [100.0]
        )

        u0 = [1.0;;]
        tspan = (0.0, 50.0)
        prob = ODEProblem(Guadex.metacommunity_ode!, u0, tspan, p)
        sol = solve(prob, Tsit5(), saveat=10.0)

        @test length(sol.t) > 0
        # Population should grow toward carrying capacity
        @test sol.u[end][1] > sol.u[1][1]
        @test sol.u[end][1] <= 100.0  # capped by carrying capacity
    end

    @testset "Positivity preserved" begin
        n_sites, n_species = 2, 2
        distances = [0.0 1000.0; 1000.0 0.0]
        elevations = [10.0, 10.0]
        D_daily = 10.0 / 365.0
        dispersal = Guadex.precompute_dispersal_matrix(
            n_sites, distances, elevations, 0.01, ones(2, 2), D_daily
        )

        p = Guadex.MetacommunityParams(
            2, 2,
            [-0.1 0.0; 0.0 -0.1],  # weak competition
            dispersal,
            [1.0, 1.0],
            [0.01 0.01; 0.01 0.01],
            [20.0, 20.0], [1.0, 1.0],
            [20.0, 20.0], [5.0, 5.0],
            [30.0, 30.0]
        )

        u0 = [1.0 1.0; 1.0 1.0]
        tspan = (0.0, 100.0)
        prob = ODEProblem(Guadex.metacommunity_ode!, u0, tspan, p)

        positivity_cb = DiscreteCallback(
            (u, t, integrator) -> any(x -> x < 0, u),
            integrator -> (integrator.u .= max.(integrator.u, 0.0));
            save_positions=(false, true)
        )

        sol = solve(prob, Tsit5(), saveat=10.0, callback=positivity_cb)

        @test all(u -> u >= -1e-10, sol.u[end])  # allow tiny numerical error
    end

    @testset "invariant: all derivatives finite" begin
        n_sites, n_species = 2, 2
        distances = [0.0 1000.0; 1000.0 0.0]
        elevations = [10.0, 50.0]
        D_daily = 10.0 / 365.0
        dispersal = Guadex.precompute_dispersal_matrix(
            n_sites, distances, elevations, 0.01, ones(2, 2), D_daily
        )

        p = Guadex.MetacommunityParams(
            2, 2,
            [-0.5 -0.3; -0.8 -0.1],
            dispersal,
            [1.0, 1.0],
            [0.01 0.02; 0.03 0.04],
            [15.0, 20.0], [0.8, 1.0],
            [19.0, 15.0], [3.67, 3.0],
            [50.0, 100.0]
        )

        # Test at various state values
        test_states = [
            [10.0 5.0; 5.0 10.0],
            [50.0 50.0; 0.0 0.0],
            [0.0 100.0; 100.0 0.0],
            [0.0 0.0; 0.0 0.0],
        ]

        for u0 in test_states
            du = zeros(2, 2)
            Guadex.metacommunity_ode!(du, u0, p, 0.0)
            @test all(isfinite.(du))
        end
    end
end
