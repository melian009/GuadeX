include("../src/ode_model.jl")
using Test
using SparseArrays

function test_ode_structure()
    # Small test case: 2 sites, 2 species
    n_sites = 2
    n_species = 2

    # Distances: site 1 to 2 is 100m
    distances = [0.0 100.0; 100.0 0.0]

    # Elevations: site 2 is higher than site 1 (upstream)
    elevations = [10.0, 50.0]

    # Interaction matrix: species 1 and 2 compete
    interaction_matrix = [0.1 0.05; 0.05 0.1]

    # Growth rates
    intrinsic_growth_rates = [0.5 0.5; 0.4 0.4]

    dispersal_intensity = 0.1
    upstream_cost = 0.01

    # No dams
    dam_passability = ones(n_sites, n_sites)

    # Precompute dispersal matrix as required by the struct
    dispersal_matrix = precompute_dispersal_matrix(
        n_sites,
        distances,
        elevations,
        upstream_cost,
        dispersal_intensity,
        dam_passability
    )

    # Environmental parameters (missing in original test)
    temperatures = [20.0, 20.0]
    habitat_suitability = [1.0, 1.0]
    thermal_optima = [20.0, 20.0]
    thermal_sigmas = [5.0, 5.0]

    p = MetacommunityParams(
        n_sites,
        n_species,
        interaction_matrix,
        dispersal_matrix,
        intrinsic_growth_rates,
        temperatures,
        habitat_suitability,
        thermal_optima,
        thermal_sigmas
    )

    # Initial populations
    u0 = [10.0 5.0; 5.0 10.0] # (sites, species)
    du = zeros(size(u0))

    # Test the ODE function
    metacommunity_ode!(du, u0, p, 0.0)

    println("Testing ODE derivative calculation...")
    println("du: ", du)

    @test size(du) == (n_sites, n_species)
    @test all(isfinite.(du))

    # Check dispersal calculation manually for one case
    # m_12 (1 to 2, upstream): m * (1/d_12) * (1 / (1 + c * (e2 - e1)))
    # 0.1 * (1/100) * (1 / (1 + 0.01 * 40)) = 0.001 * (1 / 1.4) ≈ 0.0007142857
    m_12 = p.dispersal_matrix[2, 1] # Rate FROM 1 TO 2
    @test m_12 ≈ 0.1 * (1/100) * (1 / (1 + 0.01 * 40))

    # m_21 (2 to 1, downstream): m * (1/d_21) * 1
    # 0.1 * (1/100) * 1 = 0.001
    m_21 = p.dispersal_matrix[1, 2] # Rate FROM 2 TO 1
    @test m_21 ≈ 0.001

    println("Manual dispersal checks passed.")
end

test_ode_structure()
