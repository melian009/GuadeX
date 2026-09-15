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
- `thermal_lower_limits::AbstractVector`: Empirical lower thermal limit per species
  (WP3); informational, no cold-stress term is applied.
- `thermal_upper_limits::AbstractVector`: Empirical upper thermal limit per species
  (WP3), above which the quadratic heat-stress loss is active.
- `heat_stress_rate::Real`: Shared heat-stress slope `k` (1/day/°C²). Zero
  disables the term entirely and reproduces the pre-WP3 model.

The 11-argument constructor (without the thermal limits and `k`) is retained for
backwards compatibility and fills `-Inf`/`+Inf` limits and `k = 0`, so existing
runs are unchanged.
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
    thermal_lower_limits::V
    thermal_upper_limits::V
    heat_stress_rate::T
end

"""
    MetacommunityParams(n_sites, n_species, interaction_matrix, dispersal_matrix,
                        dispersal_scaling, intrinsic_growth_rates, temperatures,
                        habitat_suitability, thermal_optima, thermal_sigmas,
                        carrying_capacity)

Backwards-compatible constructor: no empirical thermal limits and no heat-stress
mortality (`k = 0`), reproducing the model before WP3.
"""
function MetacommunityParams(n_sites::Int, n_species::Int, interaction_matrix,
        dispersal_matrix, dispersal_scaling, intrinsic_growth_rates, temperatures,
        habitat_suitability, thermal_optima, thermal_sigmas, carrying_capacity)
    T = float(promote_type(eltype(temperatures), eltype(thermal_optima),
        eltype(thermal_sigmas), eltype(dispersal_scaling), eltype(carrying_capacity)))
    vecT(x) = convert(Vector{T}, collect(T, x))
    return MetacommunityParams(n_sites, n_species, interaction_matrix, dispersal_matrix,
        vecT(dispersal_scaling), intrinsic_growth_rates, vecT(temperatures),
        vecT(habitat_suitability), vecT(thermal_optima), vecT(thermal_sigmas),
        vecT(carrying_capacity), fill(T(-Inf), n_species), fill(T(Inf), n_species), zero(T))
end

const ANNUAL_DISPERSAL_RATES = Dict{String, Float64}(
    "AB" => 0.3, "AH" => 0.2, "SP" => 3.0, "PW" => 27.5, "LS" => 20.0,
    "SA" => 5.0, "IL" => 1.25, "CP" => 1.25, "IO" => 0.5,
    "GH" => 25.0, "MS" => 17.5, "LG" => 25.0, "CC" => 20.0, "CG" => 10.0,
    "AM" => 10.0, "OM" => 12.5, "EL" => 12.5, "GL" => 3.0, "TT" => 6.0,
    "AA" => 5.0, "MC" => 75.0, "LR" => 75.0, "ST" => 12.5,
)

"""
    gaussian_thermal_filter(temp, opt, sigma)

Calculates the thermal suitability factor using a Gaussian window.
"""
function gaussian_thermal_filter(temp, opt, sigma)
    return exp(-(temp - opt)^2 / (2 * sigma^2))
end

"""
    _metacommunity_ode!(du, u, p::MetacommunityParams, delta_at, t)

Shared RHS for the metacommunity model.  `delta_at(i, t)` returns the
temperature anomaly (°C) applied to site `i` at time `t`; the static model
passes a constant-zero function and the climate model passes the interpolated
schedule.  Site temperatures and per-site total biomass are computed once per
site per RHS call rather than once per species.
"""
function _metacommunity_ode!(du, u, p::MetacommunityParams, delta_at, t)
    U = reshape(u, p.n_sites, p.n_species)
    dU = reshape(du, p.n_sites, p.n_species)

    for i in 1:p.n_sites
        temp_i = p.temperatures[i] + delta_at(i, t)
        total_biomass_i = 0.0
        for j in 1:p.n_species
            total_biomass_i += max(U[i, j], 0.0)
        end
        k_i = max(p.carrying_capacity[i], 1e-6)
        logistic_term = clamp(1.0 - total_biomass_i / k_i, -1.0, 2.0)

        for s in 1:p.n_species
            N_is = max(U[i, s], 0.0)
            env_filter = gaussian_thermal_filter(temp_i, p.thermal_optima[s], p.thermal_sigmas[s]) * p.habitat_suitability[i]
            r_eff = p.intrinsic_growth_rates[i, s] * env_filter

            interaction_term = 0.0
            for j in 1:p.n_species
                interaction_term += p.interaction_matrix[s, j] * max(U[i, j], 0.0)
            end
            interaction_term /= k_i

            # WP3: per-capita heat-stress loss, active only above the species'
            # empirical upper thermal limit.  Generic types are promoted so the
            # term is a no-op (0.0) for T<:Float32 etc.
            heat_stress = zero(temp_i)
            if p.heat_stress_rate > 0
                exceedance = temp_i - p.thermal_upper_limits[s]
                if exceedance > 0
                    heat_stress = p.heat_stress_rate * exceedance * exceedance
                end
            end

            dU[i, s] = N_is * (r_eff * logistic_term + interaction_term - heat_stress)
        end
    end

    # Dispersal.  `emigration` (column sums) is time-invariant and the
    # immigration buffer is reused across species to avoid per-species
    # allocations in the hot loop.
    immigration = zeros(Float64, p.n_sites)
    emigration = vec(sum(p.dispersal_matrix, dims=1))
    for s in 1:p.n_species
        species_pop = @view U[:, s]
        dispersal_scale = p.dispersal_scaling[s]

        mul!(immigration, p.dispersal_matrix, species_pop)

        for i in 1:p.n_sites
            dU[i, s] += dispersal_scale * (immigration[i] - emigration[i] * max(U[i, s], 0.0))
        end
    end

    return nothing
end

struct _ZeroDelta end
(::_ZeroDelta)(::Int, ::Real) = 0.0

"""
    metacommunity_ode!(du, u, p::MetacommunityParams, t)

Backwards-compatible static-temperature entry point.  The site temperatures in
`p.temperatures` are used unchanged.
"""
metacommunity_ode!(du, u, p::MetacommunityParams, t) =
    _metacommunity_ode!(du, u, p, _ZeroDelta(), t)

"""
    TemperatureSchedule

Year-by-year site temperature anomalies relative to the baseline
`MetacommunityParams.temperatures`.  `deltas[i, k]` is the anomaly (°C) for site
`i` at simulation-year node `k`; nodes are one simulation year (365 days) apart
and are linearly interpolated in time.
"""
struct TemperatureSchedule{T<:Real, M<:AbstractMatrix{T}}
    deltas::M
    days_per_year::T
end

TemperatureSchedule(deltas::AbstractMatrix, days_per_year::Real=365.0) =
    TemperatureSchedule(Float64.(deltas), Float64(days_per_year))

"""
    temperature_delta(schedule, i, t)

Interpolated temperature anomaly for site `i` at time `t` (days).
"""
@inline function temperature_delta(schedule::TemperatureSchedule, i::Int, t::Real)
    x = t / schedule.days_per_year
    n_nodes = size(schedule.deltas, 2)
    node = floor(Int, x) + 1
    node >= n_nodes && return @inbounds schedule.deltas[i, n_nodes]
    frac = x - (node - 1)
    return @inbounds schedule.deltas[i, node] * (1 - frac) + schedule.deltas[i, node + 1] * frac
end

"""
    ScheduledMetacommunityParams

Wraps static [`MetacommunityParams`](@ref) with a [`TemperatureSchedule`](@ref)
so the ODE can be solved under a year-by-year warming trajectory.
"""
struct ScheduledMetacommunityParams{P<:MetacommunityParams,S<:TemperatureSchedule}
    params::P
    schedule::S
end

struct _ScheduleDelta{S<:TemperatureSchedule}
    schedule::S
end

@inline (d::_ScheduleDelta)(i::Int, t::Real) = temperature_delta(d.schedule, i, t)

"""
    metacommunity_ode_scheduled!(du, u, p::ScheduledMetacommunityParams, t)

RHS of the metacommunity model with a time-varying (yearly interpolated)
temperature anomaly added to each site's baseline temperature.  Delegates to the
shared [`_metacommunity_ode!`](@ref) so the core model is defined once.
"""
function metacommunity_ode_scheduled!(du, u, p::ScheduledMetacommunityParams, t)
    return _metacommunity_ode!(du, u, p.params, _ScheduleDelta(p.schedule), t)
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

Note: The base dispersal matrix uses the median species dispersal rate (10 km/year).
Species-specific dispersal is handled via the dispersal_scaling vector in MetacommunityParams,
which is applied per-species in the ODE function.

The daily dispersal coefficient D_daily defaults to 10/365 km/day (the median of
literature-derived annual dispersal rates). Callers that know the exact species codes
should pass the dynamically computed median.
"""
function precompute_dispersal_matrix(n_sites, distances, elevations, c, dams, D_daily::Float64=10.0/365.0)
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

            rate = D_daily * x * dams[i, j] / max(d_km, 0.001)

            push!(I, i)
            push!(J, j)
            push!(V, rate)
        end
    end

    return sparse(I, J, V, n_sites, n_sites)
end

function precompute_dispersal_matrix(n_sites, distances, elevations, c, dams, species_codes)
    rates_for_median = Float64[get(ANNUAL_DISPERSAL_RATES, sp, 5.0) for sp in species_codes]
    D_daily = median(rates_for_median) / 365.0
    return precompute_dispersal_matrix(n_sites, distances, elevations, c, dams, D_daily)
end
