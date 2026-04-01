"""
    MetacommunityParams

Structure to hold parameters for the fish metacommunity ODE model.

# Fields
- `n_sites::Int`: Number of sites in the network.
- `n_species::Int`: Number of species in the metacommunity.
- `interaction_matrix::AbstractMatrix`: Asymmetric interaction matrix alpha_sj.
- `dispersal_matrix::AbstractMatrix`: Pre-calculated sparse matrix of base dispersal rates m_ij (in 1/day).
- `dispersal_scaling::AbstractVector`: Species-specific dispersal scaling factors to convert base rates to species-specific rates. Species with higher natural dispersal (e.g., migratory species) have higher scaling factors.
- `intrinsic_growth_rates::AbstractMatrix`: Base intrinsic growth rates r_s(E_i) for each species at each site (in 1/day).
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
    dispersal_scaling::V # Species-specific dispersal scaling
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
            # Calculates the total competitive/predatory pressure on species s at site i
            interaction_term = 0.0
            for j in 1:p.n_species
                interaction_term += p.interaction_matrix[s, j] * U[i, j]
            end

            # Local component: N_is * [r_eff + interaction_term]
            # This is logistic growth modified by interactions and environmental filtering
            dU[i, s] = U[i, s] * (r_eff + interaction_term)
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

        # Species-specific dispersal scaling factor
        # This scales the base dispersal matrix to species-specific rates
        # Higher values = more dispersive species (e.g., migratory fish)
        # Lower values = more sedentary species (e.g., small endemics)
        dispersal_scale = p.dispersal_scaling[s]

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

            # Apply species-specific dispersal scaling
            # Both immigration and emigration are scaled by the same factor
            dU[i, s] += dispersal_scale * (immigration[i] - emigration_rate * U[i, s])
        end
    end

    return nothing
end

"""
    precompute_dispersal_matrix(n_sites, distances, elevations, upstream_cost, dispersal_intensity, dam_passability)

Pre-calculates the sparse dispersal matrix M where M[i, j] is the rate FROM j TO i.

# Dispersal Rate Derivation
The dispersal rate (1/day) is calculated using a diffusion-inspired formula:
    rate = D * f_upstream * passability / distance

Where:
- D is the daily dispersal coefficient (km/day) for the species
- f_upstream is an elevation-based factor reducing upstream dispersal
- passability is a dam-based reduction factor (0-1)
- distance is the river distance between sites (km)

# Literature Sources for Daily Dispersal Coefficients
Dispersal rates from literature are typically reported as km/year and converted to daily:
    D_daily = D_annual / 365

## Native Species Dispersal Rates (km/year) and References:
- AB (Aphanius baeticus): < 0.5 km/year (fragmented habitats) (Ref 3)
- AH (Anaecypris hispanica): 0.1 - 0.3 km/year (highly sedentary) (Ref 11, 15, 16)
- SP (Squalius pyrenaicus): 1 - 5 km/year (Ref 11, 12)
- PW (Pseudochondrostoma willkommii): 15 - 40 km/year (potadromous, migratory) (Ref 9, 10)
- LS (Luciobarbus sclateri): 10 - 30 km/year (migratory) (Ref 6, 11)
- SA (Squalius alburnoides): 2 - 8 km/year (Ref 11)
- IL (Iberochondrostoma lemmingii): 0.5 - 2 km/year (Ref 11)
- CP (Cobitis paludica): 0.5 - 2 km/year (Ref 19)
- IO (Iberochondrostoma oretanum): < 1 km/year (fragmented) (Ref 5)

## Invasive Species Dispersal Rates (km/year):
- GH (Gambusia holbrooki): 8 - 42 km/year (high dispersal) (Ref 7, 10)
- MS (Micropterus salmoides): 10 - 25 km/year (Ref 2, 10)
- LG (Lepomis gibbosus): 8 - 42 km/year (Ref 5, 10)
- CC (Cyprinus carpio): 10 - 30 km/year (Ref 1)
- CG (Carassius gibelio): 5 - 15 km/year (Ref 18)
- AM (Ameiurus melas): 5 - 15 km/year (Ref 1)
- OM (Oncorhynchus mykiss): 5 - 20 km/year (Ref 1)
- EL (Esox lucius): 5 - 20 km/year (Ref 1)
- GL (Gobio lozanoi): 1 - 5 km/year (Ref 18)
- TT (Tinca tinca): 2 - 10 km/year (Ref 1)

## Other Species:
- AA (Anguilla anguilla): Limited by dams, < 10 km/year in Guadalquivir (Ref 9)
- MC (Mugil cephalus): 50 - 100 km/year (euryhaline, high mobility) (Ref 18)
- LR (Liza ramada): 50 - 100 km/year (euryhaline, high mobility) (Ref 18)
- ST (Salmo trutta): 5 - 20 km/year (typical for salmonids)

# References
See docs/Fish Growth and Dispersal Data Request.md for full citations.
"""
function precompute_dispersal_matrix(n_sites, distances, elevations, c, m, dams)
    I = Int[]
    J = Int[]
    V = Float64[]

    # Species-specific daily dispersal coefficients (km/day)
    # Converted from annual rates: D_daily = D_annual / 365
    # These are baseline diffusion coefficients for each species
    # Species with higher mobility (e.g., migratory) have higher values
    species_daily_dispersal = Dict{String, Float64}(
        # Native Endemics
        "AB" => 0.0005 / 365,   # < 0.5 km/year - highly fragmented (Ref 3)
        "AH" => 0.0008 / 365,  # ~0.3 km/year - highly sedentary (Ref 11, 15)
        "SP" => 0.008 / 365,    # ~3 km/year (Ref 11, 12)
        "PW" => 0.15 / 365,     # ~15-40 km/year - migratory (Ref 9, 10)
        "LS" => 0.055 / 365,    # ~10-30 km/year - migratory (Ref 6, 11)
        "SA" => 0.014 / 365,    # ~2-8 km/year (Ref 11)
        "IL" => 0.004 / 365,    # ~0.5-2 km/year (Ref 11)
        "CP" => 0.004 / 365,    # ~0.5-2 km/year (Ref 19)
        "IO" => 0.001 / 365,    # < 1 km/year - fragmented (Ref 5)

        # Invasive Species
        "GH" => 0.07 / 365,     # ~8-42 km/year - highly dispersive (Ref 7, 10)
        "MS" => 0.048 / 365,    # ~10-25 km/year (Ref 2, 10)
        "LG" => 0.07 / 365,     # ~8-42 km/year (Ref 5, 10)
        "CC" => 0.055 / 365,    # ~10-30 km/year (Ref 1)
        "CG" => 0.027 / 365,    # ~5-15 km/year (Ref 18)
        "AM" => 0.027 / 365,    # ~5-15 km/year (Ref 1)
        "OM" => 0.034 / 365,    # ~5-20 km/year (Ref 1)
        "EL" => 0.034 / 365,    # ~5-20 km/year (Ref 1)
        "GL" => 0.008 / 365,    # ~1-5 km/year (Ref 18)
        "TT" => 0.016 / 365,    # ~2-10 km/year (Ref 1)

        # Diadromous/Marine Species
        "AA" => 0.01 / 365,     # < 10 km/year - dam restricted (Ref 9)
        "MC" => 0.2 / 365,      # ~50-100 km/year - euryhaline (Ref 18)
        "LR" => 0.2 / 365,      # ~50-100 km/year - euryhaline (Ref 18)

        # Salmonids
        "ST" => 0.034 / 365,   # ~5-20 km/year (typical for salmonids)
    )

    for j in 1:n_sites
        for i in 1:n_sites
            d_ij = distances[i, j]
            if i == j || d_ij == 0 || isinf(d_ij)
                continue
            end

            # Rate FROM j TO i
            e_j = elevations[j]
            e_i = elevations[i]

            # Upstream dispersal penalty
            x = 1.0
            if e_i > e_j # Upstream - harder to disperse against flow
                x = 1.0 / (1.0 + c * (e_i - e_j))
            end

            # Convert distance from meters to kilometers
            d_km = d_ij / 1000.0

            # Base dispersal coefficient (km/day) - using default m as scaling factor
            # The m parameter acts as a system-wide modifier for dispersal intensity
            # Individual species rates are scaled by this factor
            D_base = m  # This is the dispersal_intensity parameter

            # Calculate rate with proper units: (km/day) / km = 1/day
            # This gives the dispersal rate in units of 1/day
            rate = D_base * x * dams[j, i] / max(d_km, 0.001)

            push!(I, i)
            push!(J, j)
            push!(V, rate)
        end
    end

    return sparse(I, J, V, n_sites, n_sites)
end
