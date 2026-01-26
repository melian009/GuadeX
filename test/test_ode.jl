include("ode_model.jl")
using Test

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

    p = MetacommunityParams(
        n_sites,
        n_species,
        distances,
        elevations,
        interaction_matrix,
        intrinsic_growth_rates,
        dispersal_intensity,
        upstream_cost,
        dam_passability
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
    # 0.1 * (1/100) * (1 / (1 + 0.01 * 40)) = 0.001 * (1 / 1.4) ≈ 0.000714
    m_12 = calculate_dispersal(1, 2, p)
    @test m_12 ≈ 0.1 * (1/100) * (1 / (1 + 0.01 * 40))

    # m_21 (2 to 1, downstream): m * (1/d_21) * 1
    # 0.1 * (1/100) * 1 = 0.001
    m_21 = calculate_dispersal(2, 1, p)
    @test m_21 ≈ 0.001

    println("Manual dispersal checks passed.")
end

test_ode_structure()
