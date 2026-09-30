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

        # Interactions are zero here, so the corrected LV bracket reduces to
        # logistic_term = clamp(1 - 35/30, -1, 2) = -0.167 and growth is still
        # N * r_eff * logistic_term < 0 (issue C2 leaves this limit unchanged).
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
# 4b. Competitive Lotka-Volterra interaction form (issue C2)
# =============================================================================
# The corrected default places the interaction term INSIDE the intrinsic-growth
# bracket:  dN/dt = N * (r_eff * clamp(1 - total/K + sum_j alpha_sj N_j/K) - heat).
# The legacy form added it OUTSIDE growth.  These tests pin the corrected
# behaviour and guard the legacy switch.
@testset "interaction form: competitive Lotka-Volterra (C2)" begin
    # --- single species, one site, no dispersal, alpha = 0 ---------------
    # With interactions inside growth the single-species equilibrium is N* = K.
    @testset "single species equilibrium is K" begin
        p = Guadex.MetacommunityParams(
            1, 1, [0.0;;], sparse([1], [1], [0.0], 1, 1),
            [1.0], [0.05;;], [20.0], [1.0], [20.0], [5.0], [100.0])

        # The corrected form must be the constructor default.
        @test p.interaction_inside_growth
        @test p.interaction_inside_growth == true

        prob = ODEProblem(Guadex.metacommunity_ode!, [1.0;;], (0.0, 600.0), p)
        sol = solve(prob, Tsit5(), saveat=10.0)
        @test sol.u[end][1] ≈ 100.0 rtol=1e-3
    end

    # --- defining property: per-capita growth equals the LV bracket -------
    # The old form put the interaction term outside growth and fails this
    # identity, so it is the core regression guard for C2.
    @testset "per-capita growth matches r_eff*(1 - total/K + sum alpha N/K)" begin
        alpha = [-0.2 -0.3; -0.1 -0.2]
        growth = [0.1 0.2]  # 1 site x 2 species
        K = 100.0
        p = Guadex.MetacommunityParams(
            1, 2, alpha, sparse([1], [1], [0.0], 1, 1),
            [1.0, 1.0], growth, [20.0], [1.0], [20.0, 20.0], [5.0, 5.0], [K])

        u0 = [10.0 5.0]
        total = sum(u0)
        du = zeros(1, 2)
        Guadex.metacommunity_ode!(du, u0, p, 0.0)

        for s in 1:2
            interaction = sum(alpha[s, j] * u0[j] / K for j in 1:2)
            expected_percap = growth[1, s] * (1.0 - total / K + interaction)
            @test du[1, s] / u0[s] ≈ expected_percap rtol=1e-10
        end
    end

    # --- r_eff scaling ----------------------------------------------------
    # In the corrected form r_eff multiplies the whole bracket, so doubling
    # r_eff doubles the per-capita growth rate.  Under the legacy form the
    # interaction term was added outside r_eff, so this did NOT hold.
    @testset "per-capita growth doubles when r_eff doubles" begin
        alpha = [-0.2 -0.3; -0.1 -0.2]
        K = 100.0
        lv(g; form=true) = Guadex.MetacommunityParams(
            1, 2, alpha, sparse([1], [1], [0.0], 1, 1),
            [1.0, 1.0], g, [20.0], [1.0], [20.0, 20.0], [5.0, 5.0], [K],
            [-Inf, -Inf], [Inf, Inf], 0.0, form)

        u0 = [10.0 5.0]
        du1 = zeros(1, 2)
        du2 = zeros(1, 2)
        Guadex.metacommunity_ode!(du1, u0, lv([0.1 0.2]), 0.0)
        Guadex.metacommunity_ode!(du2, u0, lv([0.2 0.4]), 0.0)
        for s in 1:2
            @test du2[1, s] / du1[1, s] ≈ 2.0 rtol=1e-10
        end

        # Legacy form must NOT scale linearly: the additive interaction term
        # stays fixed when r_eff doubles.
        du1_leg = zeros(1, 2)
        du2_leg = zeros(1, 2)
        Guadex.metacommunity_ode!(du1_leg, u0, lv([0.1 0.2]; form=false), 0.0)
        Guadex.metacommunity_ode!(du2_leg, u0, lv([0.2 0.4]; form=false), 0.0)
        for s in 1:2
            ratio = du2_leg[1, s] / du1_leg[1, s]
            @test abs(ratio - 2.0) > 1e-6
        end
    end

    # --- legacy switch reproduces the old expression exactly --------------
    @testset "legacy switch reproduces interaction-outside-growth" begin
        alpha = [-0.2 -0.3; -0.1 -0.2]
        growth = [0.1 0.2]
        K = 100.0
        make(form) = Guadex.MetacommunityParams(
            1, 2, alpha, sparse([1], [1], [0.0], 1, 1),
            [1.0, 1.0], growth, [20.0], [1.0], [20.0, 20.0], [5.0, 5.0], [K],
            [-Inf, -Inf], [Inf, Inf], 0.0, form)

        u0 = [10.0 5.0]
        total = sum(u0)
        du_corrected = zeros(1, 2)
        du_legacy = zeros(1, 2)
        Guadex.metacommunity_ode!(du_corrected, u0, make(true), 0.0)
        Guadex.metacommunity_ode!(du_legacy, u0, make(false), 0.0)

        for s in 1:2
            interaction = sum(alpha[s, j] * u0[j] / K for j in 1:2)
            # Old form: r_eff * clamp(1 - total/K) + interaction, all outside.
            expected_legacy = growth[1, s] * clamp(1.0 - total / K, -1.0, 2.0) + interaction
            @test du_legacy[1, s] / u0[s] ≈ expected_legacy rtol=1e-10
            # And the two forms genuinely differ when alpha != 0.
            @test !isapprox(du_corrected[1, s], du_legacy[1, s])
        end
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

# =============================================================================
# 6. Time-varying temperature schedule
# =============================================================================
@testset "scheduled metacommunity ODE" begin
    n_sites, n_species = 2, 2
    distances = [0.0 1000.0; 1000.0 0.0]
    elevations = [10.0, 50.0]
    D_daily = 10.0 / 365.0
    dispersal = Guadex.precompute_dispersal_matrix(
        n_sites, distances, elevations, 0.01, ones(2, 2), D_daily
    )
    p = Guadex.MetacommunityParams(
        n_sites, n_species,
        [-0.1 0.0; 0.0 -0.1],
        dispersal,
        [1.0, 1.0],
        [0.01 0.01; 0.01 0.01],
        [20.0, 20.0], [1.0, 1.0],
        [20.0, 20.0], [5.0, 5.0],
        [50.0, 50.0]
    )

    u0 = [10.0 5.0; 5.0 10.0]

    # Zero anomalies must reproduce the static model exactly.
    zero_schedule = Guadex.TemperatureSchedule(zeros(2, 3), 365.0)
    scheduled = Guadex.ScheduledMetacommunityParams(p, zero_schedule)
    du_base = zeros(2, 2)
    du_zero = zeros(2, 2)
    Guadex.metacommunity_ode!(du_base, u0, p, 0.0)
    Guadex.metacommunity_ode_scheduled!(du_zero, u0, scheduled, 0.0)
    @test du_base ≈ du_zero

    # A warming anomaly away from the thermal optimum reduces growth.
    warm_schedule = Guadex.TemperatureSchedule([5.0 5.0 5.0; 5.0 5.0 5.0], 365.0)
    warm = Guadex.ScheduledMetacommunityParams(p, warm_schedule)
    du_warm = zeros(2, 2)
    Guadex.metacommunity_ode_scheduled!(du_warm, u0, warm, 0.0)
    @test du_warm[1, 1] < du_base[1, 1]
end

# =============================================================================
# 7. Heat-stress mortality (WP3)
# =============================================================================
@testset "heat-stress mortality" begin
    dispersal = sparse([1], [1], [0.0], 1, 1)
    make_params(T; k=0.0, upper=28.0) = Guadex.MetacommunityParams(
        1, 1, [0.0;;], dispersal, [1.0], [0.1;;],
        [T], [1.0], [20.0], [5.0], [100.0],
        [4.0], [upper], k)

    # The backwards-compatible constructor has no limits and no stress.
    default_params = Guadex.MetacommunityParams(
        1, 1, [0.0;;], dispersal, [1.0], [0.1;;],
        [30.0], [1.0], [20.0], [5.0], [100.0])
    @test default_params.heat_stress_rate == 0.0
    @test default_params.thermal_upper_limits == [Inf]

    u0 = [10.0;;]
    du_none = zeros(1, 1)
    du_heat = zeros(1, 1)

    # Below the empirical upper limit the term is exactly zero.
    Guadex.metacommunity_ode!(du_none, u0, make_params(25.0), 0.0)
    Guadex.metacommunity_ode!(du_heat, u0, make_params(25.0; k=0.01), 0.0)
    @test du_heat[1] ≈ du_none[1]

    # Above the limit the loss turns positive growth into a decline.
    Guadex.metacommunity_ode!(du_none, u0, make_params(30.0), 0.0)
    Guadex.metacommunity_ode!(du_heat, u0, make_params(30.0; k=0.01), 0.0)
    @test du_none[1] > 0.0
    @test du_heat[1] < 0.0

    # Monotone in the exceedance.
    du_hot = zeros(1, 1)
    Guadex.metacommunity_ode!(du_hot, u0, make_params(32.0; k=0.01), 0.0)
    @test du_hot[1] < du_heat[1]
end

# =============================================================================
# 8. E4: composition-aware burn-in convergence rule
# =============================================================================
@testset "spin_up convergence rule (E4)" begin
    # Synthetic 1-site × 2-species pair: total biomass is conserved (100 -> 100)
    # but the cells exchange mass, so composition is still drifting.
    prev = [50.0, 50.0]
    cur = [55.0, 45.0]
    basin_drift = abs(sum(cur) - sum(prev)) / sum(prev)
    comp_drift = Guadex.composition_q95_change(prev, cur, 1, 2)
    @test basin_drift < 1e-12          # total-biomass criterion is satisfied
    @test comp_drift > 1e-2            # composition criterion is not
    # The old basin-only rule would stop; the new both-criteria rule must not.
    @test Guadex.spin_up_converged(:basin, 20; min_years=10, tol=1e-6,
        composition_tol=1e-3, basin_change=basin_drift, composition_change=comp_drift,
        q95_change=comp_drift, max_change=comp_drift) == true
    @test Guadex.spin_up_converged(:both, 20; min_years=10, tol=1e-6,
        composition_tol=1e-3, basin_change=basin_drift, composition_change=comp_drift,
        q95_change=comp_drift, max_change=comp_drift) == false

    # A synthetic pair that has converged in both respects.
    prev_ok = [100.0, 100.0]
    cur_ok = [100.000001, 99.999999]
    basin_ok = abs(sum(cur_ok) - sum(prev_ok)) / sum(prev_ok)
    comp_ok = Guadex.composition_q95_change(prev_ok, cur_ok, 1, 2)
    @test basin_ok < 1e-6
    @test comp_ok < 1e-3
    @test Guadex.spin_up_converged(:both, 20; min_years=10, tol=1e-6,
        composition_tol=1e-3, basin_change=basin_ok, composition_change=comp_ok,
        q95_change=comp_ok, max_change=comp_ok) == true

    # Minimum years is respected even when both criteria are already met.
    @test Guadex.spin_up_converged(:both, 9; min_years=10, tol=1e-6,
        composition_tol=1e-3, basin_change=basin_ok, composition_change=comp_ok,
        q95_change=comp_ok, max_change=comp_ok) == false
    @test Guadex.spin_up_converged(:both, 10; min_years=10, tol=1e-6,
        composition_tol=1e-3, basin_change=basin_ok, composition_change=comp_ok,
        q95_change=comp_ok, max_change=comp_ok) == true

    # Definition/floor: a cell below the active floor in both years is excluded
    # from the robust statistic but still counted by the raw all-cells statistic.
    prev_flick = [10.0, 1e-6, 10.0]
    cur_flick = [10.0, 1e-5, 10.0]
    @test Guadex.composition_q95_change(prev_flick, cur_flick, 1, 3) == 0.0
    @test Guadex.composition_q95_change_all_cells(prev_flick, cur_flick, 1, 3) > 1.0
    diag = Guadex.composition_diagnostics(prev_flick, cur_flick, 1, 3)
    @test diag.active_cells == 2
    @test diag.q95 == 0.0
    @test diag.q95_all_cells > 1.0
    # A cell exactly at the floor is excluded (`>` not `>=`).
    @test Guadex.composition_q95_change([10.0, 0.05], [10.0, 0.1], 1, 2) == 0.0
    # With a floor below every cell the robust statistic equals the all-cells one.
    @test Guadex.composition_q95_change(prev, cur, 1, 2; active_floor=1e-12) ≈
        Guadex.composition_q95_change_all_cells(prev, cur, 1, 2)
    @test Guadex.spin_up_composition_active_floor() == 0.1

    # Robust converges while the raw all-cells statistic does not: one near-empty
    # flickering cell dominates the all-cells percentile but is excluded robustly.
    stable = fill(10.0, 20)
    prev_f = copy(stable); prev_f[1] = 1e-6
    cur_f = copy(stable); cur_f[1] = 1e-5
    robust_f = Guadex.composition_q95_change(prev_f, cur_f, 1, 20)
    all_f = Guadex.composition_q95_change_all_cells(prev_f, cur_f, 1, 20)
    @test robust_f < 1e-3
    @test all_f > 1e-3

    # Converse: the all-cells statistic converges while the robust one does not.
    # Only 1% of cells are active (above the floor) and genuinely drifting; the
    # other 99% are near-empty with negligible change, so the raw percentile is
    # tiny while the robust percentile over the active cells is large.
    prev_a = fill(1e-3, 1000); cur_a = fill(1e-3, 1000)
    prev_a[1:10] .= 10.0
    cur_a[1:10] .= 10.0 .* exp(1.0)
    robust_a = Guadex.composition_q95_change(prev_a, cur_a, 1, 1000)
    all_a = Guadex.composition_q95_change_all_cells(prev_a, cur_a, 1, 1000)
    @test all_a < 1e-3
    @test robust_a > 1e-3

    # End-to-end on a cheap 1-site × 1-species system: default is the stricter
    # two-criteria rule, min_years is enforced, and the legacy rule is available.
    dispersal = sparse([1], [1], [0.0], 1, 1)
    p = Guadex.MetacommunityParams(1, 1, [0.0;;], dispersal, [1.0], [0.05;;],
        [20.0], [1.0], [20.0], [5.0], [100.0])
    res = spin_up(p; initial_state = [1.0], max_years = 20, tol = 1e-5,
        composition_tol = 1e-5, min_years = 5)
    @test res.converged
    @test res.years >= 5
    @test res.criteria == [:basin, :composition]
    @test res.last_composition_change < 1e-5

    res_late = spin_up(p; initial_state = [1.0], max_years = 20, tol = 1e-5,
        composition_tol = 1e-5, min_years = 12)
    @test res_late.converged
    @test res_late.years >= 12

    res_legacy = spin_up(p; initial_state = [1.0], max_years = 20, tol = 1e-5,
        composition_tol = 1e-5, min_years = 2, criterion = :basin)
    @test res_legacy.converged
    @test res_legacy.criteria == [:basin]
end
