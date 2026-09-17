"""
    Equilibrium spin-up and carrying-capacity sensitivity (WP4).

The climate-scenario runs previously started from the observed 2026 density
snapshot, which is far below carrying capacity, so the first decades were
dominated by a colonisation transient rather than by the forcing.  WP4 replaces
that start with a state integrated to steady state under *baseline* forcing
(no trend, no seasonality).  The spin-up is computed once and reused by every
scenario because all scenarios share the same baseline.

These helpers are deliberately independent of the reporting pipeline so they can
be unit-tested and reused by both `run_climate_scenarios.jl` and the staged
experiment script.
"""

"""
    site_totals(state, n_sites, n_species)

Total (clamped, non-negative) biomass at each site for a flattened state.
"""
function site_totals(state::AbstractVector, n_sites::Int, n_species::Int)
    U = reshape(state, n_sites, n_species)
    return [sum(max(U[i, j], 0.0) for j in 1:n_species) for i in 1:n_sites]
end

"""
    spin_up(params; initial_state, days_per_year=365.0, max_years=50, tol=1e-6,
            solver=Tsit5(), reltol=1e-6, abstol=1e-6, schedule=nothing,
            criterion=:basin)

Integrate the static (baseline-forcing) metacommunity model in one-year blocks
until the chosen convergence measure falls below `tol`, or `max_years` blocks
have been integrated.

`criterion` selects the annual-change measure used to declare convergence:

- `:basin` (default) : relative change in **basin total** biomass,
  `|Σu(t+1) − Σu(t)| / |Σu(t)|`.  Robust and dominated by the aggregate, so it
  is the right stop rule for a burn-in whose purpose is to remove the
  community-filling transient; the former max-over-sites rule could not be met
  because a single near-empty site kept the relative change large.
- `:q95`             : 95th percentile of the per-site relative change, which
  additionally requires the bulk of sites (not just the basin mean) to settle.
- `:max`             : strictest legacy rule, the largest per-site relative
  change.  Kept for backwards compatibility and diagnostics; a single site with
  near-zero biomass can make it effectively unreachable.

Returns a named tuple:

- `state`                  : final flattened state (the equilibrium initial condition)
- `converged`              : whether `tol` was reached
- `years`                  : number of one-year blocks integrated
- `last_relative_change`   : final **max** per-site relative annual change (strict diagnostic)
- `last_basin_change`      : final basin-total relative annual change
- `last_q95_change`        : final 95th-percentile per-site relative annual change
- `site_biomass`           : final site total biomass vector
- `history`                : site total biomass at the end of each block
"""
function spin_up(params::MetacommunityParams;
        initial_state::AbstractVector,
        days_per_year::Real=365.0,
        max_years::Int=50,
        tol::Real=1e-6,
        solver=Tsit5(),
        reltol::Real=1e-6,
        abstol::Real=1e-6,
        schedule::Union{Nothing,TemperatureSchedule}=nothing,
        criterion::Symbol=:basin,
        progress_every::Int=0)
    max_years >= 1 || error("max_years must be >= 1")
    criterion in (:basin, :max, :q95) ||
        error("criterion must be :basin, :max or :q95 (got $criterion)")
    n_sites, n_species = params.n_sites, params.n_species
    length(initial_state) == n_sites * n_species ||
        error("initial_state has $(length(initial_state)) entries, expected $(n_sites * n_species)")

    u0 = vec(Float64.(initial_state))
    prev_total = site_totals(u0, n_sites, n_species)
    history = Vector{Vector{Float64}}()
    converged = false
    years = 0
    # Declared before the loop so the per-block assignments below update these
    # function-local variables instead of creating loop-local ones.
    change_max = Inf
    change_basin = Inf
    change_q95 = Inf

    # Keep the state non-negative and finite; long equilibrations can otherwise
    # drift into negative densities and, eventually, non-finite values.
    clamp_cb = DiscreteCallback(
        (u, t, integrator) -> any(x -> !isfinite(x) || x < 0, u),
        integrator -> begin
            @inbounds for i in eachindex(integrator.u)
                x = integrator.u[i]
                integrator.u[i] = (isfinite(x) && x > 0) ? x : 0.0
            end
        end;
        save_positions=(false, false))

    # When a one-year baseline schedule is supplied the equilibrium is computed
    # under seasonal forcing, matching the scenario runs (the static annual-mean
    # equilibrium is not the seasonal one; see the E1 seasonality gate).
    ode! = schedule === nothing ? metacommunity_ode! : metacommunity_ode_scheduled!
    rhs_params = schedule === nothing ? params :
        ScheduledMetacommunityParams(params, schedule)

    for _ in 1:max_years
        prob = ODEProblem(ode!, u0, (0.0, Float64(days_per_year)), rhs_params)
        # Only the end state is needed: saving every internal step for a
        # 775-site × 24-species state would exhaust memory.
        sol = solve(prob, solver; reltol=reltol, abstol=abstol,
            save_everystep=false, save_start=false, save_end=true, callback=clamp_cb)
        u0 = vec(Float64.(sol.u[end]))
        if !all(isfinite, u0)
            return (state=u0, converged=false, years=years + 1,
                last_relative_change=NaN,
                last_basin_change=NaN, last_q95_change=NaN,
                site_biomass=prev_total, history=history)
        end
        total = site_totals(u0, n_sites, n_species)
        rel_site = abs.(total .- prev_total) ./ max.(abs.(prev_total), 1e-12)
        change_max = maximum(rel_site)
        change_basin = abs(sum(total) - sum(prev_total)) /
            max(abs(sum(prev_total)), 1e-12)
        change_q95 = quantile(rel_site, 0.95)
        change = criterion === :basin ? change_basin :
                 criterion === :q95 ? change_q95 : change_max
        push!(history, total)
        years += 1
        if progress_every > 0 && years % progress_every == 0
            println("    burn-in yr $(years): basin=$(round(change_basin, sigdigits=4)), " *
                    "q95=$(round(change_q95, sigdigits=4)), max=$(round(change_max, sigdigits=4))")
        end
        if change < tol && years >= 2
            converged = true
            break
        end
        prev_total = total
    end

    return (
        state=u0,
        converged=converged,
        years=years,
        last_relative_change=change_max,
        last_basin_change=change_basin,
        last_q95_change=change_q95,
        site_biomass=prev_total,
        history=history,
    )
end

"""
    scale_carrying_capacity(params, factor)

Return a copy of `params` with every site's carrying capacity multiplied by
`factor`.  `K` enters the ODE both in the logistic term and as the denominator of
the interaction term, so scaling it scales absolute abundance without changing
the per-capita interaction strength (WP4 K sensitivity: 1×, 3×, 10×).
"""
function scale_carrying_capacity(params::MetacommunityParams, factor::Real)
    factor > 0 || error("carrying-capacity scaling factor must be > 0")
    return MetacommunityParams(
        params.n_sites, params.n_species, params.interaction_matrix,
        params.dispersal_matrix, params.dispersal_scaling, params.intrinsic_growth_rates,
        params.temperatures, params.habitat_suitability, params.thermal_optima,
        params.thermal_sigmas, params.carrying_capacity .* factor,
        params.thermal_lower_limits, params.thermal_upper_limits, params.heat_stress_rate)
end

"""
    set_thermal_optima(params, optima)

Return a copy of `params` with the per-species thermal optima replaced (WP3
optimum sweep).
"""
function set_thermal_optima(params::MetacommunityParams, optima::AbstractVector)
    length(optima) == params.n_species ||
        error("optima has $(length(optima)) entries, expected $(params.n_species)")
    return MetacommunityParams(
        params.n_sites, params.n_species, params.interaction_matrix,
        params.dispersal_matrix, params.dispersal_scaling, params.intrinsic_growth_rates,
        params.temperatures, params.habitat_suitability, Float64.(optima),
        params.thermal_sigmas, params.carrying_capacity,
        params.thermal_lower_limits, params.thermal_upper_limits, params.heat_stress_rate)
end

"""
    set_heat_stress_rate(params, k)

Return a copy of `params` with the shared heat-stress slope set to `k`
(`k = 0` disables the term).
"""
function set_heat_stress_rate(params::MetacommunityParams, k::Real)
    return MetacommunityParams(
        params.n_sites, params.n_species, params.interaction_matrix,
        params.dispersal_matrix, params.dispersal_scaling, params.intrinsic_growth_rates,
        params.temperatures, params.habitat_suitability, params.thermal_optima,
        params.thermal_sigmas, params.carrying_capacity,
        params.thermal_lower_limits, params.thermal_upper_limits, Float64(k))
end

"""
    with_temperature_baseline(params, temperatures)

Return a copy of `params` whose baseline site temperatures are replaced (used to
start from the `guadex_tw` baseline climatology instead of the legacy
`TEMP_MEDIA_SC` values).
"""
function with_temperature_baseline(params::MetacommunityParams, temperatures::AbstractVector)
    length(temperatures) == params.n_sites ||
        error("temperature baseline has $(length(temperatures)) entries, expected $(params.n_sites)")
    return MetacommunityParams(
        params.n_sites, params.n_species, params.interaction_matrix,
        params.dispersal_matrix, params.dispersal_scaling, params.intrinsic_growth_rates,
        Float64.(temperatures), params.habitat_suitability, params.thermal_optima,
        params.thermal_sigmas, params.carrying_capacity,
        params.thermal_lower_limits, params.thermal_upper_limits, params.heat_stress_rate)
end
