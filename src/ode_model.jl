"""
    MetacommunityParams

Structure to hold parameters for the fish metacommunity ODE model.

# Fields
- `n_sites::Int`: Number of sites in the network.
- `n_species::Int`: Number of species in the metacommunity.
- `interaction_matrix::AbstractMatrix`: Asymmetric interaction matrix alpha_sj.
- `dispersal_matrix::AbstractMatrix`: Pre-calculated sparse matrix of dispersal rates m_ij.
- `intrinsic_growth_rates::AbstractMatrix`: Base intrinsic growth rates r_s(E_i) for each species at each site.
- `temperatures::AbstractVector`: Temperature T_i at each site.
- `habitat_suitability::AbstractVector`: Habitat suitability index h_i at each site.
- `thermal_optima::AbstractVector`: Optimal temperature for each species.
- `thermal_sigmas::AbstractVector`: Thermal tolerance (sigma) for each species.
"""
struct MetacommunityParams{T<:Real, M<:AbstractMatrix{T}, V<:AbstractVector{T}, S<:AbstractMatrix{T}}
    n_sites::Int
    n_species::Int
    interaction_matrix::M
    dispersal_matrix::S # Sparse matrix for performance
    intrinsic_growth_rates::M
    temperatures::V
    habitat_suitability::V
    thermal_optima::V
    thermal_sigmas::V
end

"""
    gaussian_thermal_filter(temp, opt, sigma)

Calculates the thermal suitability factor using a Gaussian window.
"""
function gaussian_thermal_filter(temp, opt, sigma)
    return exp(-(temp - opt)^2 / (2 * sigma^2))
end

"""
    metacommunity_ode!(du, u, p::MetacommunityParams, t)

Optimized ODE system for fish metacommunity dynamics.
u is a matrix or flattened vector where u[i, s] is the population of species s at site i.
"""
function metacommunity_ode!(du, u, p::MetacommunityParams, t)
    # Reshape u and du to (n_sites, n_species) for easier indexing
    # Note: In a production environment, consider using ComponentArrays to avoid reshape
    U = reshape(u, p.n_sites, p.n_species)
    dU = reshape(du, p.n_sites, p.n_species)

    # 1. Local Dynamics & Environmental Filtering
    for s in 1:p.n_species
        opt = p.thermal_optima[s]
        sigma = p.thermal_sigmas[s]

        for i in 1:p.n_sites
            # Environmental filtering: Gaussian thermal window * habitat suitability
            env_filter = gaussian_thermal_filter(p.temperatures[i], opt, sigma) * p.habitat_suitability[i]

            # Effective growth rate
            r_eff = p.intrinsic_growth_rates[i, s] * env_filter

            # Biotic interactions: sum_j alpha_sj * N_ij
            interaction_term = 0.0
            for j in 1:p.n_species
                interaction_term += p.interaction_matrix[s, j] * U[i, j]
            end

            # Local component: N_is * [r_eff - interaction_term]
            dU[i, s] = U[i, s] * (r_eff - interaction_term)
        end
    end

    # 2. Spatial Flux: dN_is/dt += sum_j [m_ji * N_js - m_ij * N_is]
    # This is equivalent to dU = dU + DispersalMatrix * U - diag(sum(DispersalMatrix, dims=1)) * U
    # But we iterate per species for clarity and to handle the sparse matrix efficiently.

    # Pre-calculate emigration sums if not already in params to avoid repeated work
    # For now, we use the sparse dispersal matrix directly.
    # dispersal_matrix[i, j] is the rate FROM j TO i.

    for s in 1:p.n_species
        species_pop = @view U[:, s]

        # Immigration: M * N_s (where M_ij is rate from j to i)
        # Emigration: N_is * sum_j(m_ji)

        # Using sparse matrix multiplication for immigration
        immigration = p.dispersal_matrix * species_pop

        for i in 1:p.n_sites
            # Emigration sum is the sum of the i-th column of the dispersal matrix
            # (all rates leaving site i)
            emigration_rate = 0.0
            # In a real implementation, this column sum should be pre-calculated in p
            # For this example, we assume the sparse matrix is structured for this.
            # Let's assume p.dispersal_matrix is CSC, so we can access columns.
            col_start = p.dispersal_matrix.colptr[i]
            col_end = p.dispersal_matrix.colptr[i+1] - 1
            for idx in col_start:col_end
                emigration_rate += p.dispersal_matrix.nzval[idx]
            end

            dU[i, s] += immigration[i] - emigration_rate * U[i, s]
        end
    end

    return nothing
end

"""
    precompute_dispersal_matrix(n_sites, distances, elevations, upstream_cost, dispersal_intensity, dam_passability)

Pre-calculates the sparse dispersal matrix M where M[i, j] is the rate FROM j TO i.
"""
function precompute_dispersal_matrix(n_sites, distances, elevations, c, m, dams)
    I = Int[]
    J = Int[]
    V = Float64[]

    for j in 1:n_sites
        for i in 1:n_sites
            d_ij = distances[i, j]
            if i == j || d_ij == 0 || isinf(d_ij)
                continue
            end

            # Rate FROM j TO i
            e_j = elevations[j]
            e_i = elevations[i]

            x = 1.0
            if e_i > e_j # Upstream
                x = 1.0 / (1.0 + c * (e_i - e_j))
            end

            rate = m * (1.0 / d_ij) * x * dams[j, i]

            push!(I, i)
            push!(J, j)
            push!(V, rate)
        end
    end

    return sparse(I, J, V, n_sites, n_sites)
end
