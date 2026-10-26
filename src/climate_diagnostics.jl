"""
    Climate-scenario diagnostics.

These figures explain *why* the climate ensemble behaves the way it does, as a
complement to the descriptive per-run / ensemble / across-scenario inventory in
`src/climate_figures.jl`.  Three diagnostics are produced from the exported
tables only (no model re-run):

* `forcing_and_response.png` - the warming forcing against the richness/biomass
  response, and the final-year warming-response relationship across all
  GCM x SSP runs.
* `thermal_niches.png` - species thermal response curves against the site
  temperature distribution and the end-of-horizon warming, which shows whether
  the forced warming moves sites towards or away from species optima.
* `community_filling.png` - the initial vs final per-site native richness
  distribution, which shows how much of the trend is colonisation of empty or
  species-poor sites rather than a climate response.
"""

# =============================================================================
# --- Basin series container ---
# =============================================================================

"""
    ClimateBasinSeries

Per-run basin-level series read directly from `export/levels/level_basin.csv`.
Unlike [`ClimateRun`](@ref) this keeps the climate-forcing column
(`mean_delta_temperature_c`) and the biomass columns, which the descriptive
figure metrics do not all need.
"""
struct ClimateBasinSeries
    scenario::String
    gcm::String
    run_dir::String
    years::Vector{Int}
    delta_temperature_c::Vector{Float64}
    native_richness::Vector{Float64}
    total_richness::Vector{Float64}
    invasive_richness::Vector{Float64}
    native_biomass::Vector{Float64}
    total_biomass::Vector{Float64}
    realised_richness_loss::Vector{Float64}
end

const _BASIN_SERIES_FIELDS = [
    (:delta_temperature_c, :mean_delta_temperature_c),
    (:native_richness, :mean_native_richness),
    (:total_richness, :mean_total_richness),
    (:invasive_richness, :mean_invasive_richness),
    (:native_biomass, :mean_native_biomass),
    (:total_biomass, :mean_total_biomass),
    (:realised_richness_loss, :mean_realised_richness_loss),
]

function _read_float_column(df, column::Symbol)
    hasproperty(df, column) || return fill(NaN, nrow(df))
    out = Vector{Float64}(undef, nrow(df))
    values = df[!, column]
    for i in 1:nrow(df)
        v = values[i]
        out[i] = (v isa Real && isfinite(float(v))) ? Float64(v) : NaN
    end
    return out
end

"""
    read_climate_basin_series(results_root)

Read the basin series of every discovered run, preferring `runs_index.csv` and
falling back to scanning the tree.  Runs without a basin table are skipped.
"""
function read_climate_basin_series(results_root::AbstractString)
    entries = _index_entries(results_root)
    isempty(entries) && (entries = _scan_entries(results_root))
    series = ClimateBasinSeries[]
    for entry in entries
        path = joinpath(entry.run_dir, "export", "levels", "level_basin.csv")
        isfile(path) || continue
        df = CSV.read(path, DataFrame)
        hasproperty(df, :year) || continue
        sort!(df, :year)
        years = Int[]
        for y in df.year
            ismissing(y) || push!(years, Int(y))
        end
        length(years) == nrow(df) || continue
        columns = Dict{Symbol,Vector{Float64}}()
        for (field, column) in _BASIN_SERIES_FIELDS
            columns[field] = _read_float_column(df, column)
        end
        push!(series, ClimateBasinSeries(String(entry.scenario), String(entry.gcm),
            String(entry.run_dir), years,
            columns[:delta_temperature_c], columns[:native_richness],
            columns[:total_richness], columns[:invasive_richness],
            columns[:native_biomass], columns[:total_biomass],
            columns[:realised_richness_loss]))
    end
    sort!(series; by=run -> (run.scenario, run.gcm))
    return series
end

# =============================================================================
# --- Helpers ---
# =============================================================================

const _DIAG_BAND = (0.10, 0.90)

function _series_envelope(runs::AbstractVector{<:ClimateBasinSeries}, field::Symbol,
        scenarios::AbstractVector{<:AbstractString})
    years = sort(unique(vcat([run.years for run in runs]...)))
    envelopes = Dict{String,NamedTuple}()
    for scenario in scenarios
        subset = [run for run in runs if run.scenario == scenario]
        pooled = Dict{Int,Vector{Float64}}()
        for run in subset
            values = getfield(run, field)
            for (year, value) in zip(run.years, values)
                isfinite(value) && push!(get!(pooled, year, Float64[]), value)
            end
        end
        band_lo, band_hi = _DIAG_BAND
        median_ = Float64[]
        lo = Float64[]
        hi = Float64[]
        for year in years
            values = get(pooled, year, Float64[])
            if isempty(values)
                push!(median_, NaN); push!(lo, NaN); push!(hi, NaN)
            else
                push!(median_, median(values))
                push!(lo, quantile(values, band_lo))
                push!(hi, quantile(values, band_hi))
            end
        end
        envelopes[scenario] = (years=years, median=median_, lo=lo, hi=hi, n=length(subset))
    end
    return envelopes
end

function _final_value(run::ClimateBasinSeries, field::Symbol)
    values = getfield(run, field)
    for i in length(values):-1:1
        isfinite(values[i]) && return values[i]
    end
    return NaN
end

function _pairwise_finite(xs, ys)
    x_out = Float64[]
    y_out = Float64[]
    for (x, y) in zip(xs, ys)
        (isfinite(x) && isfinite(y)) && (push!(x_out, Float64(x)); push!(y_out, Float64(y)))
    end
    return x_out, y_out
end

function _pearson(xs, ys)
    length(xs) < 2 && return NaN
    return cor(xs, ys)
end

# =============================================================================
# --- Scenario minus matched control (C3 interim route) ---
# =============================================================================

"""
    SCENARIO_CONTROL_METRICS

The response metrics differenced by [`scenario_minus_control`](@ref): native and
total richness, invasive richness, native and total biomass and the realised
native-richness-loss metric.
"""
const SCENARIO_CONTROL_METRICS = (
    :native_richness, :total_richness, :invasive_richness,
    :native_biomass, :total_biomass, :realised_richness_loss)

"""
    scenario_minus_control(scenario_run, control_run; metrics=SCENARIO_CONTROL_METRICS)

Per-metric annual `scenario − control` series for one scenario run against the
matched no-warming control run (`scenario == "control"`, `gcm == "baseline"`).

The two basin series are matched on their **shared years**; a year present in
only one run is dropped, and a metric is `NaN` wherever either side is missing or
non-finite.  Returns a `NamedTuple` with `years` plus one `<metric>_delta` vector
per requested metric, in the order given by `metrics`.

This is the reporting form used by the C3 interim projection route: those results
are provisional and are reported as deltas against the matched control, so any
residual spin-up/colonisation transient common to both runs cancels.
"""
function scenario_minus_control(scenario_run::ClimateBasinSeries,
        control_run::ClimateBasinSeries;
        metrics=SCENARIO_CONTROL_METRICS)
    control_index = Dict{Int,Int}()
    for (i, y) in enumerate(control_run.years)
        control_index[y] = i
    end
    years = Int[]
    scenario_positions = Int[]
    control_positions = Int[]
    for (i, y) in enumerate(scenario_run.years)
        j = get(control_index, y, 0)
        j == 0 && continue
        push!(years, y)
        push!(scenario_positions, i)
        push!(control_positions, j)
    end
    values = Vector{Float64}[]
    for metric in metrics
        scenario_values = getfield(scenario_run, metric)
        control_values = getfield(control_run, metric)
        delta = Vector{Float64}(undef, length(years))
        for (k, (i, j)) in enumerate(zip(scenario_positions, control_positions))
            a = Float64(scenario_values[i])
            b = Float64(control_values[j])
            delta[k] = (isfinite(a) && isfinite(b)) ? a - b : NaN
        end
        push!(values, delta)
    end
    names = Tuple(Symbol(metric, :_delta) for metric in metrics)
    return NamedTuple{(:years, names...)}((years, values...))
end

"""
    scenario_minus_control_table(runs; metrics=SCENARIO_CONTROL_METRICS,
                                 control_scenario="control") -> DataFrame

Tidy `scenario, gcm, year, metric, delta` table of every non-control run's
per-metric difference from the matched control run, for the reporting/dose-
response path.  Returns an empty `DataFrame` when the collection has no control
run, so callers can skip writing deltas without special-casing.
"""
function scenario_minus_control_table(runs::AbstractVector{<:ClimateBasinSeries};
        metrics=SCENARIO_CONTROL_METRICS,
        control_scenario::AbstractString="control")
    controls = [run for run in runs if run.scenario == control_scenario]
    isempty(controls) && return DataFrame()
    control = first(controls)
    rows = NamedTuple[]
    for run in runs
        run.scenario == control_scenario && continue
        deltas = scenario_minus_control(run, control; metrics=metrics)
        for (k, year) in enumerate(deltas.years)
            for metric in metrics
                push!(rows, (scenario=run.scenario, gcm=run.gcm, year=year,
                    metric=string(metric),
                    delta=getfield(deltas, Symbol(metric, :_delta))[k]))
            end
        end
    end
    isempty(rows) && return DataFrame()
    return DataFrame(rows)
end

# =============================================================================
# --- C7: dose-response with GCM ensemble structure, forcing axis, endpoints ---
# =============================================================================

_float_or_nan(v) = (v isa Real && isfinite(float(v))) ? Float64(v) : NaN

function _ols_slope_intercept(xs::AbstractVector, ys::AbstractVector)
    length(xs) == length(ys) || error("xs and ys must have the same length")
    n = length(xs)
    n == 0 && return (NaN, NaN)
    mx = mean(xs)
    my = mean(ys)
    sxx = sum(abs2, xs .- mx)
    sxy = sum((xs .- mx) .* (ys .- my))
    slope = sxx == 0 ? NaN : sxy / sxx
    return (slope, isfinite(slope) ? my - slope * mx : NaN)
end

function _ols_r2(xs::AbstractVector, ys::AbstractVector, slope::Real, intercept::Real)
    (isfinite(slope) && isfinite(intercept)) || return NaN
    yhat = intercept .+ slope .* xs
    ss_tot = sum(abs2, ys .- mean(ys))
    ss_res = sum(abs2, ys .- yhat)
    return ss_tot > 0 ? 1 - ss_res / ss_tot : NaN
end

"""
    dose_response(x, y, groups; min_group_n=2)

Fit the dose-response of a response `y` on the realised warming `x`, accounting
for the GCM ensemble structure **without adding a dependency**. `groups` labels
the GCM (or scenario × GCM) of each point.

The fit is deliberately a fixed-effect approximation to a per-GCM random
intercept: `y ~ x + GCM`, i.e. GCM enters as a categorical effect that absorbs
the between-GCM level shifts. `MixedModels` is *not* a dependency of this project
(see `Project.toml`), so it is not added merely for one diagnostic; the
fixed-effect slope and the per-GCM slope spread together carry the same
information (a random intercept and a fixed intercept coincide up to shrinkage).

Returns a named tuple with

* pooled OLS `slope`, `intercept`, `r2` (all points together);
* `fixed_effect_slope` / `fixed_effect_r2` from `y ~ x + GCM`;
* the per-group slopes (`group_labels`, `group_slopes`) and their spread
  (`slope_spread` = standard deviation, `slope_min`, `slope_max`);
* `between_group_sd`, the standard deviation of the per-GCM mean responses.

All arguments must have the same length; non-finite `(x, y)` pairs are dropped.
A group with fewer than `min_group_n` points or no spread in `x` gets a `NaN`
slope and is excluded from the spread.
"""
function dose_response(x::AbstractVector, y::AbstractVector, groups::AbstractVector;
        min_group_n::Int=2)
    length(x) == length(y) == length(groups) ||
        error("x, y and groups must have the same length")
    xs = Float64[]
    ys = Float64[]
    gs = String[]
    for (xi, yi, gi) in zip(x, y, groups)
        (xi isa Real && yi isa Real &&
         isfinite(float(xi)) && isfinite(float(yi))) || continue
        push!(xs, Float64(xi)); push!(ys, Float64(yi)); push!(gs, string(gi))
    end
    n = length(xs)
    empty_result = (
        slope=NaN, intercept=NaN, r2=NaN, n=0, n_groups=0,
        fixed_effect_slope=NaN, fixed_effect_r2=NaN,
        group_labels=String[], group_slopes=Float64[],
        slope_spread=NaN, slope_min=NaN, slope_max=NaN, between_group_sd=NaN)
    n < 2 && return empty_result

    pooled = _ols_slope_intercept(xs, ys)

    labels = sort(unique(gs))
    # GCM fixed effect: intercept + x + one dummy per non-reference GCM.
    Z = Matrix{Float64}(undef, n, 2 + max(0, length(labels) - 1))
    Z[:, 1] .= 1.0
    Z[:, 2] .= xs
    for (j, g) in enumerate(labels)
        j == 1 && continue
        Z[:, 1 + j] .= (gs .== g) .* 1.0
    end
    fe_slope = NaN
    fe_r2 = NaN
    if n > size(Z, 2)
        beta = Z \ ys
        fe_slope = beta[2]
        yhat = Z * beta
        ss_tot = sum(abs2, ys .- mean(ys))
        ss_res = sum(abs2, ys .- yhat)
        fe_r2 = ss_tot > 0 ? 1 - ss_res / ss_tot : NaN
    end

    group_slopes = Float64[]
    group_means = Float64[]
    for g in labels
        idx = findall(==(g), gs)
        push!(group_means, mean(ys[idx]))
        if length(idx) >= min_group_n && length(unique(xs[idx])) >= 2
            push!(group_slopes, _ols_slope_intercept(xs[idx], ys[idx])[1])
        else
            push!(group_slopes, NaN)
        end
    end
    finite_slopes = filter(isfinite, group_slopes)
    return (
        slope=pooled[1], intercept=pooled[2],
        r2=_ols_r2(xs, ys, pooled[1], pooled[2]),
        n=n, n_groups=length(labels),
        fixed_effect_slope=fe_slope, fixed_effect_r2=fe_r2,
        group_labels=labels, group_slopes=group_slopes,
        slope_spread=length(finite_slopes) >= 2 ? std(finite_slopes) : NaN,
        slope_min=isempty(finite_slopes) ? NaN : minimum(finite_slopes),
        slope_max=isempty(finite_slopes) ? NaN : maximum(finite_slopes),
        between_group_sd=length(group_means) >= 2 ? std(group_means) : NaN,
    )
end

"""
    realised_forcing_axis(index; window=:late, proxy_column=:warming_end_degc)

Return `(values, column, source, label)` for the dose-response regressor. Uses
the realised mean applied anomaly persisted by `run_climate_scenarios.jl`
(`realised_warming_2036_2045_mean_degc` for `window=:late`, the full-horizon
`realised_warming_2026_2045_mean_degc` for `window=:early`) when it is present
and finite. Falls back to the `warming_end_degc` proxy otherwise and labels both
the column and the axis so the fallback is visible in the figure.
"""
function realised_forcing_axis(index::DataFrame; window::Symbol=:late,
        proxy_column::Symbol=:warming_end_degc)
    column = window === :early ? :realised_warming_2026_2045_mean_degc :
        window === :late ? :realised_warming_2036_2045_mean_degc :
        throw(ArgumentError("window must be :early or :late"))
    window_label = window === :early ? "2026-2045" : "2036-2045"
    if hasproperty(index, column)
        values = [_float_or_nan(v) for v in index[!, column]]
        any(isfinite, values) && return (values=values, column=string(column),
            source="realised",
            label="Realised mean ΔT $window_label (°C)")
    end
    if hasproperty(index, proxy_column)
        values = [_float_or_nan(v) for v in index[!, proxy_column]]
        any(isfinite, values) && return (values=values, column=string(proxy_column),
            source="warming_end_degc_proxy",
            label="Warming by 2045 (°C) [proxy: basin_warming_curve, not applied forcing]")
    end
    throw(ArgumentError("index has neither a finite realised-forcing column " *
        "($column) nor a $proxy_column proxy; re-run run_climate_scenarios.jl"))
end

function _fallback_endpoints(cool::Tuple{String,String}, warm::Tuple{String,String},
        source::String)
    return (
        cool=(scenario=cool[1], gcm=cool[2], warming=NaN),
        warm=(scenario=warm[1], gcm=warm[2], warming=NaN),
        source=source, forcing_column="")
end

"""
    select_climate_endpoints(index; forcing_column, proxy_column, scenario_col,
                             gcm_col, fallback_cool, fallback_warm)

Re-select the cool (minimum realised forcing) and warm (maximum realised
forcing) runs from a runs index. Prefers the realised-forcing column written by
`run_climate_scenarios.jl` (default `:realised_warming_2036_2045_mean_degc`),
falls back to the `:warming_end_degc` proxy when the realised column is absent
or has no finite values, and finally to the configured pairs when neither is
available. `source` states which path was used.

Usage on the next re-run: after `run_climate_scenarios.jl` has populated the
realised columns, call

    endpoints = select_climate_endpoints(CSV.read(joinpath(output_dir, "runs_index.csv"), DataFrame))

and put `endpoints.cool.scenario`/`.gcm` and `endpoints.warm.scenario`/`.gcm`
into `[obstacle_sensitivity].climate_models` (or the
`GUADEX_SENSITIVITY_CLIMATE_MODELS` override). The current hard-coded pair is
only the fallback for runs written before the realised columns existed. For a
sensitivity index, pass `scenario_col=:climate_scenario`.
"""
function select_climate_endpoints(index::DataFrame;
        forcing_column::Symbol=:realised_warming_2036_2045_mean_degc,
        proxy_column::Symbol=:warming_end_degc,
        scenario_col::Symbol=:scenario, gcm_col::Symbol=:gcm,
        fallback_cool::Tuple{String,String}=("ssp126", "IITM-ESM"),
        fallback_warm::Tuple{String,String}=("ssp245", "UKESM1-0-LL"),
        exclude::Tuple{Vararg{String}}=("control",))
    hasproperty(index, scenario_col) ||
        throw(ArgumentError("index has no scenario column $scenario_col"))
    hasproperty(index, gcm_col) ||
        throw(ArgumentError("index has no gcm column $gcm_col"))
    nrow(index) == 0 && return _fallback_endpoints(fallback_cool, fallback_warm, "empty_index")
    keep = [string(r[scenario_col]) ∉ exclude for r in eachrow(index)]
    sub = index[keep, :]
    nrow(sub) == 0 && return _fallback_endpoints(fallback_cool, fallback_warm,
        "no_warming_runs")

    function _select(column::Symbol, source::String)
        pairs = [(_float_or_nan(r[column]), i) for (i, r) in enumerate(eachrow(sub))]
        finite = [p for p in pairs if isfinite(p[1])]
        length(finite) < 2 && return nothing
        sort!(finite; by=first)
        cool, warm = finite[1], finite[end]
        return (
            cool=(scenario=string(sub[cool[2], scenario_col]),
                gcm=string(sub[cool[2], gcm_col]), warming=cool[1]),
            warm=(scenario=string(sub[warm[2], scenario_col]),
                gcm=string(sub[warm[2], gcm_col]), warming=warm[1]),
            source=source, forcing_column=string(column))
    end

    result = hasproperty(sub, forcing_column) ? _select(forcing_column, "realised") : nothing
    if result === nothing && hasproperty(sub, proxy_column)
        result = _select(proxy_column, "warming_end_degc_proxy")
    end
    result === nothing &&
        (result = _fallback_endpoints(fallback_cool, fallback_warm, "configured_fallback"))
    return result
end

# =============================================================================
# --- Forcing vs response ---
# =============================================================================

"""
    plot_climate_forcing_response(runs, output_path)

Forcing-versus-response diagnostic in a 3 x 2 grid, with the SSP legend in its
own row below the panels.  `forcing` is the imposed basin-mean warming anomaly
(`delta_temperature_c`, relative to the 1986-2005 baseline): it is the driver
that the climate ensemble adds on top of the observed seasonal climatology, and
it is shown explicitly in panel (a) because the report otherwise plots only
absolute water temperature.  Panels (a)-(c) are annual trajectories (median
across the GCMs, 10-90% band); panels (d)-(f) are the final-year response
plotted against the *realised* warming of each run, which is the quantity used
to separate runs along the forcing axis.  Each point is one GCM x SSP run.
"""
function plot_climate_forcing_response(runs::AbstractVector{<:ClimateBasinSeries},
        output_path::AbstractString)
    report_theme!()
    isempty(runs) && return nothing
    scenarios = sort(unique(run.scenario for run in runs))
    warming = _series_envelope(runs, :delta_temperature_c, scenarios)
    richness = _series_envelope(runs, :native_richness, scenarios)
    biomass = _series_envelope(runs, :total_biomass, scenarios)

    fig = Figure(size=(1500, 1950))

    function _trajectory_axis(row, col, letter, ylabel, envelope)
        ax = Axis(fig[row, col]; width=520, height=430, xlabel="Year", ylabel=ylabel)
        panel_letter!(ax, letter)
        for scenario in scenarios
            env = envelope[scenario]
            (isempty(env.years) || all(isnan, env.median)) && continue
            color = scenario_color(scenario)
            band!(ax, env.years, env.lo, env.hi; color=(color, 0.12))
            lines!(ax, env.years, env.median; color=color, linewidth=2.4)
        end
        return ax
    end

    ax1 = _trajectory_axis(1, 1, "(a)", "Warming anomaly ΔT (°C)", warming)
    ax2 = _trajectory_axis(1, 2, "(b)", "Native richness (sp./site)", richness)
    ax3 = _trajectory_axis(2, 1, "(c)", "Total biomass (units/site)", biomass)
    # The four SSP trajectories are almost identical in the response panels; say
    # so explicitly, otherwise it looks like only one scenario was plotted.
    for ax in (ax2, ax3)
        report_annotation!(ax, "all SSPs overlap", position=:lb, color=(:black, 0.65))
    end

    # Final-year warming response across all runs.
    final_warming = [_final_value(run, :delta_temperature_c) for run in runs]
    final_richness = [_final_value(run, :native_richness) for run in runs]
    final_biomass = [_final_value(run, :total_biomass) for run in runs]
    final_risk = [_final_value(run, :realised_richness_loss) for run in runs]
    run_scenario = [run.scenario for run in runs]
    # The no-warming control is plotted (it is a useful reference point) but it
    # is excluded from the dose-response fits so they describe the warming runs,
    # matching the analysis elsewhere in the report.
    scenario_mask = [run.scenario != "control" for run in runs]
    run_gcm = [run.gcm for run in runs]

    function _response_axis(row, col, letter, ylabel, ys)
        ax = Axis(fig[row, col]; width=520, height=430,
            xlabel = "Realised ΔT by 2045 (°C)", ylabel = ylabel)
        panel_letter!(ax, letter)
        for scenario in scenarios
            idx = findall(==(scenario), run_scenario)
            x, y = _pairwise_finite(final_warming[idx], ys[idx])
            isempty(x) && continue
            scatter!(ax, x, y; color=(scenario_color(scenario), 0.85),
                markersize=11, strokewidth=0)
        end
        # Dose-response with GCM as a categorical fixed effect (no MixedModels
        # dependency); the per-GCM slope spread is shown next to the pooled slope.
        dose = dose_response(final_warming[scenario_mask], ys[scenario_mask],
            run_gcm[scenario_mask])
        if isfinite(dose.slope)
            text = "slope = $(round(dose.slope, digits=3))"
            isfinite(dose.slope_spread) &&
                (text *= " ± $(round(dose.slope_spread, digits=3)) (GCM SD)")
            text *= " per °C\nR² = $(round(dose.r2, digits=2)) " *
                    "(pooled, n=$(dose.n), GCMs=$(dose.n_groups))"
            report_annotation!(ax, text, position=:lt)
        end
        return ax
    end

    _response_axis(2, 2, "(d)", "Native richness (sp./site)", final_richness)
    _response_axis(3, 1, "(e)", "Total biomass (units/site)", final_biomass)
    _response_axis(3, 2, "(f)", "Richness loss (fraction)", final_risk)
    equal_panel_columns!(fig, 1, 1)

    Legend(fig[4, 1:2],
        [LineElement(color=scenario_color(s), linewidth=3.0) for s in scenarios],
        scenarios; orientation=:horizontal, framevisible=false,
        title="SSP pathway", titlefontsize=REPORT_LEGEND_FONTSIZE)
    return _save_climate_figure(fig, output_path)
end

# =============================================================================
# --- Thermal niches ---
# =============================================================================

function _read_species_thermal_table(species_chars_file::AbstractString)
    isfile(species_chars_file) ||
        error("species characteristics file not found: $species_chars_file")
    df = CSV.read(species_chars_file, DataFrame; delim=';')
    table = Dict{String,NamedTuple{(:optimum, :sigma, :label),Tuple{Float64,Float64,String}}}()
    for row in eachrow(df)
        code = lowercase(string(row.SP))
        params = parse_temperature_range_and_sigma(string(row.TEMPERATURE_C))
        table[code] = (optimum=params.optimum, sigma=params.sigma, label=string(row.SP))
    end
    return table
end

function _read_site_temperatures_baseline(site_table_path::AbstractString)
    isfile(site_table_path) || return Float64[]
    df = CSV.read(site_table_path, DataFrame)
    (hasproperty(df, :year) && hasproperty(df, :temperature_c)) || return Float64[]
    first_year = minimum(df.year)
    values = Float64[]
    for row in eachrow(df)
        Int(row.year) == first_year || continue
        v = row.temperature_c
        (v isa Real && isfinite(float(v))) && push!(values, Float64(v))
    end
    return values
end

"""
    plot_climate_thermal_niches(species_chars_file, site_table_path, output_path;
                               native_codes, extra_codes, warming_delta)

Species thermal response curves (`exp(-(T-opt)^2 / 2σ^2)`) for the richness
species plus optional context species, overlaid on the relative frequency of
baseline site temperatures.  Vertical lines mark the mean baseline temperature
and the mean end-of-horizon temperature after `warming_delta` °C.
"""
function plot_climate_thermal_niches(species_chars_file::AbstractString,
        site_table_path::AbstractString, output_path::AbstractString;
        native_codes::AbstractVector{<:AbstractString}=String[],
        extra_codes::AbstractVector{<:AbstractString}=String[],
        warming_delta::Real=1.0)
    report_theme!()
    thermal = _read_species_thermal_table(species_chars_file)
    site_temps = _read_site_temperatures_baseline(site_table_path)

    lo, hi = 8.0, 26.0

    # Match the report-plot rendering of Fig. 1 (fig01_thermal_niches): explicit
    # panel size, legend outside the axes, and no descriptive title (the LaTeX
    # caption carries the explanation).
    fig = Figure(size=(1700, 1050))
    ax = Axis(fig[1, 1]; width=760, height=650,
        xlabel="Water temperature (°C)",
        ylabel="Thermal suitability / rel. frequency")

    t_grid = range(lo - 0.5, hi + 1.0; length=300)

    # Site-temperature distribution as translucent bars behind the curves.
    if !isempty(site_temps)
        edges = collect(lo:1.0:hi)
        counts = zeros(Int, length(edges) - 1)
        for t in site_temps
            (lo <= t < hi) || continue
            counts[clamp(floor(Int, t - lo) + 1, 1, length(counts))] += 1
        end
        max_count = maximum(counts; init=0)
        max_count == 0 || barplot!(ax, edges[1:(end - 1)] .+ 0.5,
            counts ./ max_count; width=0.95, color=(:gray, 0.18),
            label="site temperatures (rel. freq.)")
    end

    # Many native species share the same thermal range (8-30 °C), so group
    # identical curves rather than drawing the same line several times.
    palette = [:dodgerblue, :seagreen, :darkorange, :firebrick, :mediumpurple,
        :teal, :goldenrod, :slateblue, :olive]
    groups = Dict{Tuple{Float64,Float64},Vector{String}}()
    for code in native_codes
        key = lowercase(string(code))
        haskey(thermal, key) || continue
        p = thermal[key]
        push!(get!(groups, (round(p.optimum, digits=4), round(p.sigma, digits=4)),
            String[]), uppercase(string(code)))
    end
    # Legend handles are built explicitly so the legend lives outside the axes,
    # exactly as in Fig. 1 rather than as an in-panel axis legend.
    handles = Any[PolyElement(color=(:gray, 0.35))]
    labels = String["Site-temperature distribution (relative frequency)"]
    for (i, (key, codes)) in enumerate(sort(collect(groups); by=first))
        optimum, sigma = key
        curve = [exp(-(t - optimum)^2 / (2 * sigma^2)) for t in t_grid]
        colour = palette[mod1(i, length(palette))]
        lines!(ax, t_grid, curve; color=colour, linewidth=2.6)
        push!(handles, LineElement(color=colour, linewidth=2.6))
        push!(labels, length(codes) == 1 ? "$(only(codes)): optimum $(round(optimum, digits=1)) °C" :
            "optimum $(round(optimum, digits=1)) °C ×$(length(codes)) ($(join(codes, ",")))")
    end

    # Context species that drive some local declines but are not in the metric.
    for code in extra_codes
        key = lowercase(string(code))
        haskey(thermal, key) || continue
        p = thermal[key]
        curve = [exp(-(t - p.optimum)^2 / (2 * p.sigma^2)) for t in t_grid]
        lines!(ax, t_grid, curve; color=(:black, 0.55), linewidth=2.0, linestyle=:dash)
        push!(handles, LineElement(color=(:black, 0.55), linewidth=2.0, linestyle=:dash))
        push!(labels, "$(uppercase(string(code))) (not in native metric)")
    end

    if !isempty(site_temps)
        baseline = mean(site_temps)
        vlines!(ax, [baseline]; color=:black, linewidth=2.0, linestyle=:dot)
        vlines!(ax, [baseline + warming_delta]; color=:crimson, linewidth=2.0, linestyle=:dot)
        text!(ax, baseline - 0.15, 1.06;
            text="Baseline mean $(round(baseline, digits=1)) °C",
            align=(:right, :bottom), fontsize=REPORT_ANNOTATION_FONTSIZE)
        text!(ax, baseline + warming_delta + 0.15, 1.06;
            text="+$(round(warming_delta, digits=2)) °C → $(round(baseline + warming_delta, digits=1)) °C",
            align=(:left, :bottom), fontsize=REPORT_ANNOTATION_FONTSIZE, color=:crimson)
    end
    ylims!(ax, 0.0, 1.12)

    Legend(fig[1, 2], handles, labels; framevisible=false,
        title="Native species (optimum)", titlefontsize=REPORT_LEGEND_FONTSIZE, nbanks=1)
    equal_panel_columns!(fig, 1, 0.55)
    resize_to_layout!(fig)
    return _save_climate_figure(fig, output_path)
end

# =============================================================================
# --- Community filling ---
# =============================================================================

function _site_table(run_dir::AbstractString)
    path = joinpath(run_dir, "export", "levels", "level_sampling_point.csv")
    isfile(path) || return nothing
    return CSV.read(path, DataFrame)
end

"""
    plot_climate_community_filling(run_dir, output_path)

Initial and final per-site native-richness distributions plus the site-by-site
change.  This isolates the colonisation transient from the warming response.
"""
function plot_climate_community_filling(run_dir::AbstractString,
        output_path::AbstractString)
    site = _site_table(run_dir)
    site === nothing && return nothing
    hasproperty(site, :native_richness) || return nothing

    years = sort(unique(Int.(site.year)))
    first_year, last_year = Base.first(years), Base.last(years)
    initial_rows = filter(row -> Int(row.year) == first_year, site)
    final_rows = filter(row -> Int(row.year) == last_year, site)

    joined = innerjoin(
        select(initial_rows, :CODIGO, :native_richness, :total_richness),
        select(final_rows, :CODIGO, :native_richness, :total_richness);
        on=:CODIGO, makeunique=true, order=:left)
    n_sites = nrow(joined)
    initial = Float64.(joined.native_richness)
    final = Float64.(joined.native_richness_1)
    change = final .- initial

    max_k = max(ceil(Int, maximum(vcat(initial, final))), 1)

    fig = Figure(size=(1250, 560))
    ax1 = Axis(fig[1, 1]; title="Per-site native richness distribution",
        xlabel="native species per site", ylabel="number of sites", titlesize=12)
    xs = 0:max_k
    counts_first = [count(==(k), round.(Int, initial)) for k in xs]
    counts_last = [count(==(k), round.(Int, final)) for k in xs]
    barplot!(ax1, xs .- 0.18, counts_first; width=0.34, color=(:steelblue, 0.8),
        label="$first_year")
    barplot!(ax1, xs .+ 0.18, counts_last; width=0.34, color=(:firebrick, 0.8),
        label="$last_year")
    axislegend(ax1; position=:rt, fontsize=10)

    ax2 = Axis(fig[1, 2]; title="Initial vs final richness per site",
        xlabel="native richness in $first_year", ylabel="native richness in $last_year",
        titlesize=12)
    jitter = 0.12 .* sin.(1:n_sites)
    scatter!(ax2, initial .+ jitter, final .- jitter;
        color=(:steelblue, 0.18), markersize=6)
    limit = max(max_k, 1) + 0.4
    lines!(ax2, [0, limit], [0, limit]; color=:black, linewidth=1.2, linestyle=:dash)
    up = count(>(0), change)
    down = count(<(0), change)
    flat = count(==(0), change)
    text!(ax2, 0.03, 0.97; space=:relative,
        text="$(up) sites gained, $(flat) unchanged, $(down) lost species",
        align=(:left, :top), fontsize=11)

    Label(fig[0, :],
        "Community filling: richness rises mainly where sites start empty or species-poor ($(run_dir))",
        fontsize=13, font=:bold)
    return _save_climate_figure(fig, output_path)
end

# =============================================================================
# --- Orchestration ---
# =============================================================================

"""
    plot_climate_diagnostics(results_root; figures_dir, species_chars_file,
        density_file, native_codes, extra_codes, representative_run)

Generate the three climate diagnostics under `<figures_dir>/diagnostics` and
return their paths.  `representative_run` defaults to the first complete run.
"""
function plot_climate_diagnostics(results_root::AbstractString=joinpath("results", "climate_scenarios");
        figures_dir::AbstractString=joinpath(results_root, "figures"),
        species_chars_file::AbstractString=joinpath("data", "ABIOTIC",
            "caracteristicas_peces_Guadalquivir_03-04-2018.csv"),
        density_file::AbstractString=joinpath("data", "BIOTIC",
            "FishDensity_and_Juveniles_Matrix.csv"),
        native_codes::AbstractVector{<:AbstractString}=["AB", "AH", "SP", "PW",
            "LS", "SA", "IL", "CP", "IO", "ST"],
        extra_codes::AbstractVector{<:AbstractString}=String[],
        representative_run::Union{Nothing,AbstractString}=nothing)
    runs = read_climate_basin_series(results_root)
    if isempty(runs)
        @warn "no climate runs found under $results_root; no diagnostics written"
        return String[]
    end
    diagnostic_dir = joinpath(figures_dir, "diagnostics")
    mkpath(diagnostic_dir)
    written = String[]

    # C3: scenario - matched-control deltas (written whenever a control run
    # exists, so the corrected interim-route results are reported as deltas).
    delta_table = scenario_minus_control_table(runs)
    if !isempty(delta_table)
        delta_path = joinpath(diagnostic_dir, "scenario_minus_control.csv")
        CSV.write(delta_path, delta_table)
        println("[climate_diagnostics] scenario - control deltas: $delta_path " *
                "($(nrow(delta_table)) rows)")
    end

    path = joinpath(diagnostic_dir, "forcing_and_response.png")
    saved = _safe_figure(() -> plot_climate_forcing_response(runs, path), path)
    saved === nothing || push!(written, saved)

    run_dir = representative_run
    if run_dir === nothing
        candidates = [run.run_dir for run in runs
                      if isfile(joinpath(run.run_dir, "export", "levels", "level_sampling_point.csv"))]
        run_dir = isempty(candidates) ? nothing : first(candidates)
    end

    if run_dir !== nothing
        site_table = joinpath(run_dir, "export", "levels", "level_sampling_point.csv")
        median_delta = median([_final_value(run, :delta_temperature_c) for run in runs])
        path = joinpath(diagnostic_dir, "thermal_niches.png")
        saved = _safe_figure(() -> plot_climate_thermal_niches(species_chars_file,
            site_table, path; native_codes=native_codes, extra_codes=extra_codes,
            warming_delta=median_delta), path)
        saved === nothing || push!(written, saved)

        path = joinpath(diagnostic_dir, "community_filling.png")
        saved = _safe_figure(() -> plot_climate_community_filling(run_dir, path), path)
        saved === nothing || push!(written, saved)
    end

    println("[climate_diagnostics] $(length(written)) diagnostic figure(s) written to $diagnostic_dir")
    return written
end
