"""
    Temperature forcing, thermal-stress calibration and exposure diagnostics.

This file implements the infrastructure added by the GuadeX modelling
improvement plan:

* WP1 - loaders for the per-site daily water-temperature projections produced by
  the `guadex_tw` pipeline (one row per site/scenario/date, ensemble median).
* WP2 - builders that turn either an annual warming curve (the historical
  behaviour) or a daily series into a [`TemperatureSchedule`](@ref), plus a
  helper that recovers annual-mean anomalies from a daily schedule.
* WP3 - the empirical-optimum sweep and the baseline-viability calibration of the
  shared heat-stress slope `k`, and the exposure diagnostics that attribute a
  response to days above a species' empirical upper limit.

Nothing here changes the default behaviour: the annual-mean path is the default
and the heat-stress term is active only when `k > 0`.
"""

using Dates

# =============================================================================
# --- WP1: per-site daily forcing ---
# =============================================================================

const DAILY_FORCING_REQUIRED_COLUMNS = ["site_id", "scenario", "date"]

"""
    load_daily_temperature_forcing(path)

Load a per-site daily water-temperature forcing table.  The primary product of
the `guadex_tw` pipeline is a parquet/csv with columns

    site_id, scenario, date, tw_ensemble_median

and (optionally) `gcm` and `tw`.  Only the columns required to build a schedule
are validated; extra columns are preserved.
"""
function load_daily_temperature_forcing(path::AbstractString)
    isfile(path) || error("daily temperature forcing file not found: $path")
    df = CSV.read(path, DataFrame)
    missing_columns = setdiff(DAILY_FORCING_REQUIRED_COLUMNS, names(df))
    isempty(missing_columns) ||
        error("daily forcing file $path is missing columns: $(join(missing_columns, ", "))")
    value_columns = [c for c in ("tw_ensemble_median", "tw", "tw_mean") if c in names(df)]
    isempty(value_columns) &&
        error("daily forcing file $path has no temperature value column " *
              "(expected tw_ensemble_median, tw or tw_mean)")
    if !(eltype(df[!, :date]) <: Date)
        df[!, :date] = Date.(string.(df[!, :date]))
    end
    return df
end

"""
    daily_forcing_matrix(df, sites; scenario, value_col=:tw_ensemble_median, gcm=nothing)

Pivot one scenario (and optional GCM) of a daily forcing table into an
`n_sites × n_days` matrix aligned with `sites`, plus the sorted vector of dates.
Missing site/date combinations are `NaN`.
"""
function daily_forcing_matrix(df::DataFrame, sites::AbstractVector;
        scenario::AbstractString, value_col::Symbol=:tw_ensemble_median,
        gcm::Union{Nothing,AbstractString}=nothing)
    if !hasproperty(df, value_col)
        candidates = [c for c in ("tw_ensemble_median", "tw", "tw_mean") if hasproperty(df, c)]
        isempty(candidates) && error("daily forcing table has no recognised value column")
        value_col = Symbol(candidates[1])
    end
    sub = filter(row -> string(row.scenario) == scenario, df)
    if gcm !== nothing
        hasproperty(sub, :gcm) || error("gcm filter requested but table has no gcm column")
        sub = filter(row -> string(row.gcm) == gcm, sub)
    end
    isempty(sub) && error("no daily forcing rows for scenario=$scenario" *
                          (gcm === nothing ? "" : " gcm=$gcm"))

    dates = sort(unique(sub.date))
    date_index = Dict(d => k for (k, d) in enumerate(dates))
    site_index = Dict(string(s) => i for (i, s) in enumerate(sites))
    temps = fill(NaN, length(sites), length(dates))
    for row in eachrow(sub)
        i = get(site_index, string(row.site_id), nothing)
        i === nothing && continue
        temps[i, date_index[row.date]] = _to_float(row[value_col])
    end
    return dates, temps
end

"""
    is_wide_daily_forcing(df)

True when a daily forcing table is in the compact wide layout produced by
`guadex_tw/scripts/13_project_guadex_sites.py`: `date, scenario, <site columns>`
instead of the long `site_id, scenario, date, tw_ensemble_median` layout.
"""
is_wide_daily_forcing(df::DataFrame) =
    !hasproperty(df, :site_id) && hasproperty(df, :scenario) && hasproperty(df, :date)

"""
    load_daily_forcing_any(path)

Load a daily temperature forcing table in either the long layout
([`load_daily_temperature_forcing`](@ref)) or the wide layout
([`load_wide_daily_forcing`](@ref)).
"""
function load_daily_forcing_any(path::AbstractString)
    isfile(path) || error("daily temperature forcing file not found: $path")
    df = CSV.read(path, DataFrame)
    if hasproperty(df, :site_id)
        missing_columns = setdiff(DAILY_FORCING_REQUIRED_COLUMNS, names(df))
        isempty(missing_columns) ||
            error("daily forcing file $path is missing columns: $(join(missing_columns, ", "))")
        value_columns = [c for c in ("tw_ensemble_median", "tw", "tw_mean") if c in names(df)]
        isempty(value_columns) &&
            error("daily forcing file $path has no temperature value column")
    end
    if !(eltype(df[!, :date]) <: Date)
        df[!, :date] = Date.(string.(df[!, :date]))
    end
    return df
end

"""
    wide_forcing_matrix(df, sites; scenario)

Pivot one scenario of a wide daily forcing table into an `n_sites × n_days`
matrix aligned with `sites`, plus the sorted dates.  Site columns are matched by
their string label (e.g. `1.1.2`).
"""
function wide_forcing_matrix(df::DataFrame, sites::AbstractVector;
        scenario::AbstractString)
    sub = filter(row -> string(row.scenario) == scenario, df)
    isempty(sub) && error("wide daily forcing has no rows for scenario=$scenario")
    dates = Date.(sub[!, :date])
    temps = Matrix{Float64}(undef, length(sites), length(dates))
    for (i, site) in enumerate(sites)
        col = Symbol(string(site))
        hasproperty(sub, col) || error("wide daily forcing has no column for site '$(string(site))'")
        temps[i, :] = [_to_float(v) for v in sub[!, col]]
    end
    return dates, temps
end

# Missing/empty cells in a forcing CSV become `missing`; NaN is used internally.
_to_float(v) = v === missing ? NaN : Float64(v)

"""
    wide_baseline_means(df, sites; baseline_scenario="historical",
                        baseline_start, baseline_end)

Per-site baseline mean Tw from a wide forcing table (used when only the scalar
baseline, not the full historical daily series, is needed).
"""
function wide_baseline_means(df::DataFrame, sites::AbstractVector;
        baseline_scenario::AbstractString="historical")
    dates, temps = wide_forcing_matrix(df, sites; scenario=baseline_scenario)
    return [_finite_mean(@view temps[i, :]) for i in 1:length(sites)]
end

"""
    _finite_mean(values)

Mean of the finite entries of `values` (`NaN` when none are finite).  Used
instead of `skipmissing`, which does not skip `NaN` floats.
"""
function _finite_mean(values)
    total = 0.0
    n = 0
    for v in values
        if isfinite(v)
            total += v
            n += 1
        end
    end
    return n == 0 ? NaN : total / n
end

# =============================================================================
# --- WP2: schedule builders ---
# =============================================================================

"""
    daily_temperature_schedule(temps; dates, baseline_start, baseline_end)

Build a [`TemperatureSchedule`](@ref) from an `n_sites × n_days` matrix of
absolute daily water temperatures.  The anomaly stored in the schedule is
relative to each site's own mean over the baseline window `baseline_start` ..
`baseline_end` (inclusive, by calendar year), so already-realised warming is not
discarded and the site's observed temperature level is preserved when the
anomaly is added to `MetacommunityParams.temperatures`.

Returns `(schedule, baseline_means)`.
"""
function daily_temperature_schedule(temps::AbstractMatrix;
        dates::AbstractVector, baseline_start::Int, baseline_end::Int,
        baseline_temps::Union{Nothing,AbstractMatrix}=nothing,
        baseline_dates::Union{Nothing,AbstractVector}=nothing)
    size(temps, 2) == length(dates) ||
        error("temps has $(size(temps, 2)) day columns but $(length(dates)) dates")
    if baseline_temps === nothing
        mask = [baseline_start <= year(d) <= baseline_end for d in dates]
        any(mask) || error("no days in baseline window $baseline_start-$baseline_end; " *
                           "supply baseline_temps/baseline_dates for a separate baseline series")
        baseline_means = [_finite_mean(@view temps[i, mask]) for i in 1:size(temps, 1)]
    else
        size(baseline_temps, 1) == size(temps, 1) ||
            error("baseline_temps and temps must have the same number of sites")
        mask = [baseline_start <= year(d) <= baseline_end for d in baseline_dates]
        any(mask) || error("no days in baseline window $baseline_start-$baseline_end")
        baseline_means = [_finite_mean(@view baseline_temps[i, mask]) for i in 1:size(temps, 1)]
    end
    deltas = temps .- reshape(baseline_means, :, 1)
    schedule = TemperatureSchedule(deltas, 1.0)
    return schedule, baseline_means
end

"""
    annual_mean_schedule(deltas; days_per_year=365.0)

Build a [`TemperatureSchedule`](@ref) from an `n_sites × n_years` anomaly matrix,
reproducing the pre-WP2 annual-node behaviour.
"""
annual_mean_schedule(deltas::AbstractMatrix; days_per_year::Real=365.0) =
    TemperatureSchedule(deltas, Float64(days_per_year))

"""
    baseline_climatology_schedule(baseline_temps, baseline_dates, n_years;
                                  days_per_year=365.0)

Repeat the baseline-period daily climatology (per site, per day of year) for
`n_years` simulation years and return a daily [`TemperatureSchedule`](@ref)
whose anomalies are relative to each site's climatological mean.  This is the
no-trend baseline forcing used by the WP6 seasonality and calibration stages.
"""
function baseline_climatology_schedule(baseline_temps::AbstractMatrix,
        baseline_dates::AbstractVector, n_years::Int; days_per_year::Real=365.0)
    n_sites = size(baseline_temps, 1)
    n_days = Int(days_per_year)
    sums = zeros(n_sites, n_days)
    counts = zeros(Int, n_days)
    for (k, d) in enumerate(baseline_dates)
        doy = dayofyear(d)
        (1 <= doy <= n_days) || continue
        counts[doy] += 1
        for i in 1:n_sites
            v = baseline_temps[i, k]
            isfinite(v) && (sums[i, doy] += v)
        end
    end
    clim = fill(NaN, n_sites, n_days)
    for doy in 1:n_days
        counts[doy] > 0 && (clim[:, doy] .= sums[:, doy] ./ counts[doy])
    end
    tiled = repeat(clim, 1, n_years)
    baseline_means = [mean(skipmissing(@view clim[i, :])) for i in 1:n_sites]
    deltas = tiled .- reshape(baseline_means, :, 1)
    return TemperatureSchedule(deltas, 1.0)
end

"""
    annual_mean_deltas_by_year(schedule, dates)

Calendar-year mean anomaly matrix from a daily schedule, returning
`(years, matrix)` with `matrix` an `n_sites × length(years)` array.  Unlike
[`annual_mean_deltas`](@ref) this keeps the year labels, so callers can align the
reporting horizon exactly.
"""
function annual_mean_deltas_by_year(schedule::TemperatureSchedule, dates::AbstractVector)
    size(schedule.deltas, 2) == length(dates) ||
        error("schedule has $(size(schedule.deltas, 2)) nodes but $(length(dates)) dates")
    labels = sort(unique(year.(dates)))
    out = Matrix{Float64}(undef, size(schedule.deltas, 1), length(labels))
    for (k, y) in enumerate(labels)
        mask = [year(d) == y for d in dates]
        out[:, k] = vec(mean(schedule.deltas[:, mask], dims=2))
    end
    return labels, out
end

"""
    annual_mean_deltas(schedule; days_per_year=365.0)

Recover the `n_sites × n_years` annual-mean anomaly matrix from a daily (or
annual) schedule.  Used for reporting the projected annual temperature and for
exposure metrics.
"""
function annual_mean_deltas(schedule::TemperatureSchedule; days_per_year::Real=365.0)
    n_days = size(schedule.deltas, 2)
    if schedule.days_per_year == 1.0
        n_years = floor(Int, n_days / days_per_year)
        n_years <= 0 && return copy(schedule.deltas)
        out = Matrix{Float64}(undef, size(schedule.deltas, 1), n_years)
        for y in 1:n_years
            lo = round(Int, (y - 1) * days_per_year) + 1
            hi = min(round(Int, y * days_per_year), n_days)
            out[:, y] = vec(mean(schedule.deltas[:, lo:hi], dims=2))
        end
        return out
    end
    return copy(schedule.deltas)
end

# =============================================================================
# --- WP3: optimum sweep, heat-stress calibration, exposure ---
# =============================================================================

"""
    optimum_sweep_optima(species_chars_df, species_codes; margin=2.0, fraction=0.5)

Place every species' thermal optimum at `fraction` of the way across its
*empirical viable range*, inset by `margin` °C at each end:

    opt = (lower + margin) + fraction * ((upper - margin) - (lower + margin))

`fraction = 0` gives the cold edge, `0.5` the midpoint (the default) and `1` the
warm edge.  This is the documented sensitivity sweep required by WP3 (E5): the
midpoint is not calibrated, so the sign of the warming response must be reported
across the range.  Species without a finite empirical range keep 15 °C.
"""
function optimum_sweep_optima(species_chars_df, species_codes;
        margin::Real=2.0, fraction::Real=0.5)
    0 <= fraction <= 1 || error("fraction must be in [0, 1], got $fraction")
    lookup = Dict(lowercase(string(r.SP)) =>
        (hasproperty(r, :thermal_lower) ? r.thermal_lower : -Inf,
         hasproperty(r, :thermal_upper) ? r.thermal_upper : Inf)
        for r in eachrow(species_chars_df))
    optima = Float64[]
    for sp in species_codes
        info = get(lookup, lowercase(string(sp)), nothing)
        if info === nothing || !isfinite(info[1]) || !isfinite(info[2])
            push!(optima, 15.0)
            continue
        end
        lo = info[1] + margin
        hi = info[2] - margin
        if hi <= lo  # range narrower than 2*margin: fall back to the bounds
            lo, hi = info[1], info[2]
        end
        push!(optima, lo + fraction * (hi - lo))
    end
    return optima
end

"""
    exceedance_energy(temperatures, upper_limit; days_per_year=365.0)

Mean annual squared exceedance energy `Σ (T − T_upper)_+² / n_years` for one
site/species.  `temperatures` is a daily series covering a whole number of years
(with `NaN` allowed); the result is the sum divided by the number of years.
"""
function exceedance_energy(temperatures::AbstractVector, upper_limit::Real;
        days_per_year::Real=365.0)
    isfinite(upper_limit) || return 0.0
    n = length(temperatures)
    n == 0 && return 0.0
    energy = 0.0
    for t in temperatures
        (isfinite(t) && t > upper_limit) || continue
        energy += (t - upper_limit)^2
    end
    return energy / max(1.0, n / days_per_year)
end

"""
    calibrate_heat_stress_rate(exceedance_energies; max_annual_loss=0.05)

Baseline-viability calibration of the single shared heat-stress slope `k`
(WP3).  The per-capita annual survival under a quadratic heat-stress loss is
`exp(-k · E)`, with `E` the annual exceedance energy; requiring the worst
observed baseline species/site to keep at least `1 - max_annual_loss` survival
gives

    k = -log(1 - max_annual_loss) / max(E).

`exceedance_energies` is any iterable of annual exceedance energies (e.g. per
species/site under baseline forcing).  Returns `k` in 1/day/°C².
"""
function calibrate_heat_stress_rate(exceedance_energies; max_annual_loss::Real=0.05)
    0 < max_annual_loss < 1 || error("max_annual_loss must be in (0, 1)")
    worst = 0.0
    for e in exceedance_energies
        e = float(e)
        isfinite(e) && e > worst && (worst = e)
    end
    worst <= 0 && return 0.0
    return -log(1 - max_annual_loss) / worst
end

"""
    exposure_days(temperatures, upper_limit)

Number of days with `T > T_upper` for one series.
"""
function exposure_days(temperatures::AbstractVector, upper_limit::Real)
    isfinite(upper_limit) || return 0
    return count(t -> isfinite(t) && t > upper_limit, temperatures)
end

"""
    exposure_table(sites, species, temps_by_site, dates; upper_limits)

Per site × species × year exposure diagnostics: days above the empirical upper
limit and the annual squared exceedance energy.  `temps_by_site` is an
`n_sites × n_days` matrix aligned with `dates`.
"""
function exposure_table(sites::AbstractVector, species::AbstractVector,
        temps_by_site::AbstractMatrix, dates::AbstractVector;
        upper_limits::AbstractVector)
    size(temps_by_site, 1) == length(sites) ||
        error("temps_by_site has $(size(temps_by_site, 1)) rows but $(length(sites)) sites")
    size(temps_by_site, 2) == length(dates) ||
        error("temps_by_site has $(size(temps_by_site, 2)) columns but $(length(dates)) dates")
    labels = sort(unique(year.(dates)))
    buckets = [Int[] for _ in labels]
    for (k, d) in enumerate(dates)
        push!(buckets[searchsortedfirst(labels, year(d))], k)
    end

    rows = NamedTuple[]
    for (i, site) in enumerate(sites)
        for (s, sp) in enumerate(species)
            upper = Float64(upper_limits[s])
            for (k, y) in enumerate(labels)
                idx = buckets[k]
                n_days = count(t -> isfinite(t) && t > upper, @view temps_by_site[i, idx])
                energy = sum(max(t - upper, 0.0)^2 for t in @view temps_by_site[i, idx]
                             if isfinite(t))
                push!(rows, (year=Int(y), CODIGO=string(site), species=string(sp),
                    exposure_days=n_days, exceedance_energy=energy))
            end
        end
    end
    return DataFrame(rows)
end
