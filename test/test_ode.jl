@testset "ODE model test" begin
    n_sites = 2
    n_species = 2

    distances = [0.0 100.0; 100.0 0.0]
    elevations = [10.0, 50.0]
    interaction_matrix = [0.1 0.05; 0.05 0.1]
    intrinsic_growth_rates = [0.5 0.5; 0.4 0.4]
    upstream_cost = 0.01
    dam_passability = ones(n_sites, n_sites)

    dispersal_matrix = precompute_dispersal_matrix(
        n_sites,
        distances,
        elevations,
        upstream_cost,
        dam_passability
    )

    temperatures = [20.0, 20.0]
    habitat_suitability = [1.0, 1.0]
    thermal_optima = [20.0, 20.0]
    thermal_sigmas = [5.0, 5.0]
    carrying_capacity = [30.0, 30.0]

    p = MetacommunityParams(
        n_sites,
        n_species,
        interaction_matrix,
        dispersal_matrix,
        [1.0, 1.0],
        intrinsic_growth_rates,
        temperatures,
        habitat_suitability,
        thermal_optima,
        thermal_sigmas,
        carrying_capacity
    )

    u0 = [10.0 5.0; 5.0 10.0]
    du = zeros(size(u0))

    metacommunity_ode!(du, u0, p, 0.0)

    println("Testing ODE derivative calculation...")
    println("du: ", du)

    @test size(du) == (n_sites, n_species)
    @test all(isfinite.(du))

    m_12 = p.dispersal_matrix[2, 1]
    @test m_12 ≈ (1.0 / (1 + 0.01 * 40)) / 0.1

    m_21 = p.dispersal_matrix[1, 2]
    @test m_21 ≈ 1.0 / 0.1
end
