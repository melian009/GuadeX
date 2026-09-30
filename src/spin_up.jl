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

# Densities at or below this are treated as the same (extinct) cell when taking
# log-density differences, so an exactly-zero cell does not produce an infinite
# composition change.
const SPIN_UP_LOG_FLOOR = 1e-12

# A site × species cell is "active" (established enough to be part of the
# community) when its density exceeds this in at least one of the two compared
# years.  This is the project's presence/occupancy threshold — the same 0.1
# density units used for richness/occupancy elsewhere (`presence_threshold`) —
# so the robust composition criterion is measured over the community that is
# actually present, and near-extinct or flickering cells cannot dominate the
# 95th percentile.
const SPIN_UP_COMPOSITION_ACTIVE_FLOOR = 0.1

"""
    spin_up_composition_active_floor()

The density floor above which a site × species cell counts as active for the
robust composition criterion (see [`composition_q95_change`](@ref)).
"""
spin_up_composition_active_floor() = SPIN_UP_COMPOSITION_ACTIVE_FLOOR

function _check_spin_up_states(prev_state, state, n_sites, n_species)
    length(prev_state) == n_sites * n_species ||
        error("prev_state has $(length(prev_state)) entries, expected $(n_sites * n_species)")
    length(state) == n_sites * n_species ||
        error("state has $(length(state)) entries, expected $(n_sites * n_species)")
    return nothing
end

"""
    composition_q95_change(prev_state, state, n_sites, n_species;
                           active_floor=SPIN_UP_COMPOSITION_ACTIVE_FLOOR,
                           log_floor=SPIN_UP_LOG_FLOOR)

**Robust** composition change: the 95th percentile, over the *active* site ×
species cells, of the annual change in log density `|Δ log N|`.

A cell is active when its density exceeds `active_floor` in **at least one** of
the two compared years, so cells that are effectively extinct (including cells
that flicker around zero) are excluded and cannot dominate the percentile.  The
default floor is the project's presence threshold (0.1 density units).  Densities
are floored at `log_floor` before the log so an exactly-zero density in an
otherwise active cell does not produce an infinite change.  Returns `0.0` when
no cell is active.

The raw all-cells statistic is available separately as
[`composition_q95_change_all_cells`](@ref).
"""
function composition_q95_change(prev_state::AbstractVector, state::AbstractVector,
        n_sites::Int, n_species::Int;
        active_floor::Real=SPIN_UP_COMPOSITION_ACTIVE_FLOOR,
        log_floor::Real=SPIN_UP_LOG_FLOOR)
    _check_spin_up_states(prev_state, state, n_sites, n_species)
    deltas = Float64[]
    @inbounds for i in eachindex(prev_state)
        a_raw = Float64(prev_state[i]); b_raw = Float64(state[i])
        max(a_raw, b_raw) > active_floor || continue
        a = max(a_raw, log_floor)
        b = max(b_raw, log_floor)
        push!(deltas, abs(log(b) - log(a)))
    end
    isempty(deltas) && return 0.0
    return quantile(deltas, 0.95)
end

"""
    composition_q95_change_all_cells(prev_state, state, n_sites, n_species;
                                     log_floor=SPIN_UP_LOG_FLOOR)

Raw diagnostic companion to [`composition_q95_change`](@ref): the 95th
percentile of `|Δ log N|` over **all** site × species cells (no active floor).
Recorded for diagnostics but not used by the stop rule: on real data a handful of
near-empty, flickering cells keep it far above any practical tolerance.
"""
function composition_q95_change_all_cells(prev_state::AbstractVector, state::AbstractVector,
        n_sites::Int, n_species::Int; log_floor::Real=SPIN_UP_LOG_FLOOR)
    _check_spin_up_states(prev_state, state, n_sites, n_species)
    deltas = Vector{Float64}(undef, n_sites * n_species)
    @inbounds for i in eachindex(deltas)
        a = max(Float64(prev_state[i]), log_floor)
        b = max(Float64(state[i]), log_floor)
        deltas[i] = abs(log(b) - log(a))
    end
    return quantile(deltas, 0.95)
end

"""
    composition_diagnostics(prev_state, state, n_sites, n_species;
                            active_floor=SPIN_UP_COMPOSITION_ACTIVE_FLOOR,
                            log_floor=SPIN_UP_LOG_FLOOR)

One-pass burn-in diagnostics: `(q95, q95_all_cells, active_cells)`, where `q95`
is the robust active-cell statistic, `q95_all_cells` the raw all-cells statistic
and `active_cells` the number of cells above `active_floor` in at least one of
the two compared years.
"""
function composition_diagnostics(prev_state::AbstractVector, state::AbstractVector,
        n_sites::Int, n_species::Int;
        active_floor::Real=SPIN_UP_COMPOSITION_ACTIVE_FLOOR,
        log_floor::Real=SPIN_UP_LOG_FLOOR)
    _check_spin_up_states(prev_state, state, n_sites, n_species)
    active = Float64[]
    all_deltas = Float64[]
    @inbounds for i in eachindex(prev_state)
        a_raw = Float64(prev_state[i]); b_raw = Float64(state[i])
        a = max(a_raw, log_floor)
        b = max(b_raw, log_floor)
        delta = abs(log(b) - log(a))
        push!(all_deltas, delta)
        max(a_raw, b_raw) > active_floor && push!(active, delta)
    end
    return (
        q95 = isempty(active) ? 0.0 : quantile(active, 0.95),
        q95_all_cells = quantile(all_deltas, 0.95),
        active_cells = length(active),
    )
end

"""
    spin_up_criteria(criterion)

The individual convergence criteria that `criterion` requires, as a
`Vector{Symbol}`.  Recorded in run metadata so the stop rule is auditable.
"""
function spin_up_criteria(criterion::Symbol)
    criterion === :both && return [:basin, :composition]
    criterion === :basin && return [:basin]
    criterion === :q95 && return [:q95]
    criterion === :max && return [:max]
    error("criterion must be :both, :basin, :q95 or :max (got $criterion)")
end

"""
    spin_up_converged(criterion, years; min_years, tol, composition_tol,
                      basin_change, composition_change, q95_change, max_change)

Whether the burn-in stop rule is satisfied.  Always requires at least
`min_years` one-year blocks.  `:both` then requires **both** the basin-total
change below `tol` and the **robust** composition change (the active-cell
statistic from [`composition_q95_change`](@ref)) below `composition_tol`;
`:basin`, `:q95` and `:max` are the legacy single-criterion rules.
"""
function spin_up_converged(criterion::Symbol, years::Integer;
        min_years::Integer, tol::Real, composition_tol::Real,
        basin_change::Real, composition_change::Real,
        q95_change::Real, max_change::Real)
    years >= min_years || return false
    if criterion === :both
        return basin_change < tol && composition_change < composition_tol
    elseif criterion === :basin
        return basin_change < tol
    elseif criterion === :q95
        return q95_change < tol
    elseif criterion === :max
        return max_change < tol
    end
    error("criterion must be :both, :basin, :q95 or :max (got $criterion)")
end

"""
    spin_up(params; initial_state, days_per_year=365.0, max_years=50, tol=1e-6,
            composition_tol=1e-2, min_years=10, solver=Tsit5(), reltol=1e-6,
            abstol=1e-6, schedule=nothing, criterion=:both, progress_every=0)

Integrate the static (baseline-forcing) metacommunity model in one-year blocks
until the chosen convergence criteria are met, or `max_years` blocks have been
integrated.

The stop rule always requires at least `min_years` blocks.  `criterion` then
selects which annual-change measures must additionally fall below tolerance:

- `:both` (default) : **both** the basin-total relative biomass change
  `|Σu(t+1) − Σu(t)| / |Σu(t)|` below `tol` **and** the composition change below
  `composition_tol`.  The composition change is the **robust** statistic from
  [`composition_q95_change`](@ref): the 95th percentile of `|Δ log N|` over the
  site × species cells whose density exceeds the active floor (the project's
  0.1-unit presence threshold) in at least one of the two compared years.
  Restricting to active cells keeps near-extinct/flickering cells from dominating
  the percentile.  The raw all-cells statistic is still recorded as
  `last_composition_change_all_cells` for diagnostics.  Requiring both fixes the
  weak-criterion failure in which basin biomass has settled while the community
  composition keeps drifting.
- `:basin`          : basin-total relative change below `tol` only.  This is the
  legacy rule (pair it with `min_years=2` to reproduce the pre-E4 behaviour).
- `:q95`            : 95th percentile of the per-site relative **total** change
  below `tol` (legacy per-site diagnostic; not the composition criterion above).
- `:max`            : strictest legacy rule, the largest per-site relative total
  change below `tol`.  A single site with near-zero biomass can make it
  effectively unreachable.

Returns a named tuple:

- `state`                   : final flattened state (the equilibrium initial condition)
- `converged`               : whether every required criterion was met
- `years`                   : number of one-year blocks integrated
- `criteria`                : criteria required by the run (Vector{Symbol})
- `last_relative_change`    : final **max** per-site relative annual change (strict diagnostic)
- `last_basin_change`       : final basin-total relative annual change
- `last_composition_change` : final robust active-cell 95th-percentile `|Δ log N|`
- `last_composition_change_all_cells` : final raw all-cells 95th-percentile `|Δ log N|`
- `composition_active_cells` : number of active cells in the final composition statistic
- `last_q95_change`         : final 95th-percentile per-site relative total change
- `site_biomass`            : final site total biomass vector
- `history`                 : site total biomass at the end of each block
"""
function spin_up(params::MetacommunityParams;
        initial_state::AbstractVector,
        days_per_year::Real=365.0,
        max_years::Int=50,
        tol::Real=1e-6,
        composition_tol::Real=1e-2,
        min_years::Int=10,
        solver=Tsit5(),
        reltol::Real=1e-6,
        abstol::Real=1e-6,
        schedule::Union{Nothing,TemperatureSchedule}=nothing,
        criterion::Symbol=:both,
        progress_every::Int=0)
    max_years >= 1 || error("max_years must be >= 1")
    min_years >= 1 || error("min_years must be >= 1")
    criterion in (:both, :basin, :max, :q95) ||
        error("criterion must be :both, :basin, :max or :q95 (got $criterion)")
    n_sites, n_species = params.n_sites, params.n_species
    length(initial_state) == n_sites * n_species ||
        error("initial_state has $(length(initial_state)) entries, expected $(n_sites * n_species)")

    u0 = vec(Float64.(initial_state))
    prev_total = site_totals(u0, n_sites, n_species)
    # The full previous state is retained so the composition criterion can take
    # per-cell log-density differences, not just per-site totals.
    prev_state = copy(u0)
    history = Vector{Vector{Float64}}()
    converged = false
    years = 0
    criteria = spin_up_criteria(criterion)
    # Declared before the loop so the per-block assignments below update these
    # function-local variables instead of creating loop-local ones.
    change_max = Inf
    change_basin = Inf
    change_q95 = Inf
    change_composition = Inf
    change_composition_all = Inf
    composition_active = 0

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
                criteria=criteria,
                last_relative_change=NaN,
                last_basin_change=NaN, last_q95_change=NaN,
                last_composition_change=NaN,
                last_composition_change_all_cells=NaN,
                composition_active_cells=0,
                site_biomass=prev_total, history=history)
        end
        total = site_totals(u0, n_sites, n_species)
        rel_site = abs.(total .- prev_total) ./ max.(abs.(prev_total), 1e-12)
        change_max = maximum(rel_site)
        change_basin = abs(sum(total) - sum(prev_total)) /
            max(abs(sum(prev_total)), 1e-12)
        change_q95 = quantile(rel_site, 0.95)
        comp = composition_diagnostics(prev_state, u0, n_sites, n_species)
        change_composition = comp.q95
        change_composition_all = comp.q95_all_cells
        composition_active = comp.active_cells
        push!(history, total)
        years += 1
        if progress_every > 0 && years % progress_every == 0
            println("    burn-in yr $(years): basin=$(round(change_basin, sigdigits=4)), " *
                    "comp_q95=$(round(change_composition, sigdigits=4)) " *
                    "(all=$(round(change_composition_all, sigdigits=4)), " *
                    "active=$(composition_active)), " *
                    "q95=$(round(change_q95, sigdigits=4)), max=$(round(change_max, sigdigits=4))")
        end
        if spin_up_converged(criterion, years;
                min_years=min_years, tol=tol, composition_tol=composition_tol,
                basin_change=change_basin, composition_change=change_composition,
                q95_change=change_q95, max_change=change_max)
            converged = true
            break
        end
        prev_total = total
        prev_state = copy(u0)
    end

    return (
        state=u0,
        converged=converged,
        years=years,
        criteria=criteria,
        last_relative_change=change_max,
        last_basin_change=change_basin,
        last_q95_change=change_q95,
        last_composition_change=change_composition,
        last_composition_change_all_cells=change_composition_all,
        composition_active_cells=composition_active,
        site_biomass=prev_total,
        history=history,
    )
end

"""
    projection_route(value; default="full_burnin") -> Symbol

Resolve the C3 projection-route setting.  Accepted values are `"full_burnin"`
(the default; integrate to the E4 convergence criteria — see [`spin_up`](@ref))
and `"interim_observed"` (the team's interim route: start from the observed
community and integrate only a short fixed number of baseline years, then report
every scenario relative to the matched control).  An invalid value errors so a
typo cannot silently change the route; an absent/empty value falls back to
`default`, so existing configurations keep the full-burn-in behaviour.
"""
function projection_route(value; default::AbstractString="full_burnin")
    (value === nothing || string(value) == "") && return Symbol(default)
    route = Symbol(lowercase(string(value)))
    route in (:full_burnin, :interim_observed) ||
        error("projection_route must be 'full_burnin' or 'interim_observed' (got '$value')")
    return route
end

"""
    interim_observed_spin_up(params; initial_state, years, schedule=nothing,
                             days_per_year=365.0, ...) -> (state, years, spin)

C3 interim projection route.  Integrate the metacommunity model for **exactly**
`years` one-year baseline blocks from the observed `initial_state` (no
convergence stop: the E4 criterion is deliberately bypassed here and only here,
because the team chose a short fixed relaxation rather than a long burn-in).

`years = 0` performs no integration and returns the observed state unchanged (the
`spin` field is `nothing`).  For `years >= 1` the underlying [`spin_up`](@ref) is
called with `max_years = years`, `tol = 0` and `min_years = 1`, so it always runs
the full fixed number of blocks.  Returns `(state, years, spin)`; `spin` is the
full [`spin_up`](@ref) result for diagnostics/metadata (or `nothing`).
"""
function interim_observed_spin_up(params::MetacommunityParams;
        initial_state::AbstractVector,
        years::Integer,
        schedule::Union{Nothing,TemperatureSchedule}=nothing,
        days_per_year::Real=365.0,
        solver=Tsit5(),
        reltol::Real=1e-6,
        abstol::Real=1e-6,
        progress_every::Int=0)
    years >= 0 || error("interim spin-up years must be >= 0 (got $years)")
    n_expected = params.n_sites * params.n_species
    length(initial_state) == n_expected ||
        error("initial_state has $(length(initial_state)) entries, expected $n_expected")
    if years == 0
        return (state=vec(Float64.(initial_state)), years=0, spin=nothing)
    end
    spin = spin_up(params; initial_state=initial_state, schedule=schedule,
        days_per_year=days_per_year, max_years=Int(years), tol=0.0,
        composition_tol=0.0, min_years=1, criterion=:basin, progress_every=progress_every,
        solver=solver, reltol=reltol, abstol=abstol)
    return (state=spin.state, years=spin.years, spin=spin)
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
    set_dispersal_matrix(params, matrix)

Return a copy of `params` with the precomputed dispersal matrix replaced.  Used
by the obstacle/passability sweeps, where the dispersal network is rebuilt for
each scenario while every other parameter (including the empirical thermal
limits and heat-stress slope) is preserved.
"""
function set_dispersal_matrix(params::MetacommunityParams, matrix)
    return MetacommunityParams(
        params.n_sites, params.n_species, params.interaction_matrix,
        matrix, params.dispersal_scaling, params.intrinsic_growth_rates,
        params.temperatures, params.habitat_suitability, params.thermal_optima,
        params.thermal_sigmas, params.carrying_capacity,
        params.thermal_lower_limits, params.thermal_upper_limits, params.heat_stress_rate)
end

"""
    set_interaction_matrix(params, matrix)

Return a copy of `params` with the species interaction matrix replaced (used by
the alternative-interaction sensitivity).
"""
function set_interaction_matrix(params::MetacommunityParams, matrix)
    return MetacommunityParams(
        params.n_sites, params.n_species, matrix,
        params.dispersal_matrix, params.dispersal_scaling, params.intrinsic_growth_rates,
        params.temperatures, params.habitat_suitability, params.thermal_optima,
        params.thermal_sigmas, params.carrying_capacity,
        params.thermal_lower_limits, params.thermal_upper_limits, params.heat_stress_rate)
end

"""
    set_thermal_sigma_multiplier(params, factor)

Return a copy of `params` with every species' thermal tolerance scaled by
`factor` (narrowing sigma sharpens the thermal filter).
"""
function set_thermal_sigma_multiplier(params::MetacommunityParams, factor::Real)
    factor > 0 || error("thermal sigma multiplier must be > 0")
    return MetacommunityParams(
        params.n_sites, params.n_species, params.interaction_matrix,
        params.dispersal_matrix, params.dispersal_scaling, params.intrinsic_growth_rates,
        params.temperatures, params.habitat_suitability, params.thermal_optima,
        params.thermal_sigmas .* Float64(factor), params.carrying_capacity,
        params.thermal_lower_limits, params.thermal_upper_limits, params.heat_stress_rate)
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
