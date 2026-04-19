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
- `carrying_capacity::AbstractVector`: Site-specific carrying capacity K_i for total biomass at each site.
"""
struct MetacommunityParams{T<:Real, M<:AbstractMatrix{T}, V<:AbstractVector{T}, S<:AbstractMatrix{T}}
    n_sites::Int
    n_species::Int
    interaction_matrix::M
    dispersal_matrix::S
    dispersal_scaling::V
    intrinsic_growth_rates::M
    temperatures::V
    habitat_suitability::V
    thermal_optima::V
    thermal_sigmas::V
    carrying_capacity::V
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

Optimized ODE system for fish metacommunity dynamics with logistic growth.
u is a matrix or flattened vector where u[i, s] is the population of species s at site i.

The model includes:
- Local logistic growth with carrying capacity K_i for total biomass at site i
- Environmental filtering (thermal niche and habitat suitability)
- Interspecific interactions (competition, predation)
- Species-specific dispersal between sites

The base dispersal matrix uses median species dispersal rates (5 km/year / 365).
Species-specific scaling factors (dispersal_scaling) adjust per-species to achieve
the correct relative dispersal rates based on literature values.
"""
function metacommunity_ode!(du, u, p::MetacommunityParams, t)
    U = reshape(u, p.n_sites, p.n_species)
    dU = reshape(du, p.n_sites, p.n_species)

    for s in 1:p.n_species
        opt = p.thermal_optima[s]
        sigma = p.thermal_sigmas[s]
        dispersal_scale = p.dispersal_scaling[s]

        for i in 1:p.n_sites
            total_biomass_i = 0.0
            for j in 1:p.n_species
                total_biomass_i += U[i, j]
            end

            env_filter = gaussian_thermal_filter(p.temperatures[i], opt, sigma) * p.habitat_suitability[i]
            r_eff = p.intrinsic_growth_rates[i, s] * env_filter

            interaction_term = 0.0
            for j in 1:p.n_species
                interaction_term += p.interaction_matrix[s, j] * U[i, j]
            end

            logistic_term = 1.0 - total_biomass_i / p.carrying_capacity[i]
            dU[i, s] = U[i, s] * (r_eff * logistic_term + interaction_term)
        end
    end

    for s in 1:p.n_species
        species_pop = @view U[:, s]
        dispersal_scale = p.dispersal_scaling[s]

        immigration = p.dispersal_matrix * species_pop

        for i in 1:p.n_sites
            emigration_rate = 0.0
            col_start = p.dispersal_matrix.colptr[i]
            col_end = p.dispersal_matrix.colptr[i+1] - 1
            for idx in col_start:col_end
                emigration_rate += p.dispersal_matrix.nzval[idx]
            end

            dU[i, s] += dispersal_scale * (immigration[i] - emigration_rate * U[i, s])
        end
    end

    return nothing
end

"""
    precompute_dispersal_matrix(n_sites, distances, elevations, c, dams, species_codes)

Pre-calculates the sparse dispersal matrix M where M[i, j] is the rate FROM j TO i.

# Dispersal Rate Derivation
The dispersal rate (1/day) is calculated using:
    rate = D_s * f_upstream * passability / distance

Where:
- D_s is the daily dispersal coefficient (km/day) for species s, derived from literature
- f_upstream is an elevation-based factor reducing upstream dispersal
- passability is a dam-based reduction factor (0-1)
- distance is the river distance between sites (km)

# Literature Sources for Daily Dispersal Coefficients
Dispersal rates from literature are reported as km/year and converted to daily:
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

Note: The base dispersal matrix uses the median species dispersal rate (5 km/year).
Species-specific dispersal is handled via the dispersal_scaling vector in MetacommunityParams,
which is applied per-species in the ODE function.
"""
function precompute_dispersal_matrix(n_sites, distances, elevations, c, dams)
    I = Int[]
    J = Int[]
    V = Float64[]

    for j in 1:n_sites
        for i in 1:n_sites
            d_ij = distances[i, j]
            if i == j || d_ij == 0 || isinf(d_ij)
                continue
            end

            e_j = elevations[j]
            e_i = elevations[i]

            x = 1.0
            if e_i > e_j
                x = 1.0 / (1.0 + c * (e_i - e_j))
            end

            d_km = d_ij / 1000.0

            rate = x * dams[j, i] / max(d_km, 0.001)

            push!(I, i)
            push!(J, j)
            push!(V, rate)
        end
    end

    return sparse(I, J, V, n_sites, n_sites)
end

function precompute_dispersal_matrix(n_sites, distances, elevations, c, dams, species_codes)
    return precompute_dispersal_matrix(n_sites, distances, elevations, c, dams)
end
