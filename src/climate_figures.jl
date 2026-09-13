"""
    Climate-scenario figure module.

Reusable helpers that turn the four-level CSVs written by
`run_climate_scenarios.jl` into report figures.  Everything here reads the
exported tables (`export/levels/*.csv`) and, optionally, `run_metadata.json`;
the simulation state (`simulation_output.jld2`) is never touched, so figures can
be regenerated without re-running the model.

Reporting levels (see `docs/climate_scenarios.md`):

| level | GuadeX id | group column | metric column |
| :--- | :--- | :--- | :--- |
| `:site` | `CODIGO` | `CODIGO` (one row per site) | `<metric>` |
| `:subcatchment` | `CODIGO_S` | `subcatchment` | `mean_<metric>` |
| `:water_body` | `ID_masa` | `water_body` | `mean_<metric>` |
| `:basin` | `ES050` | `basin` | `mean_<metric>` |

Ensemble bands pool every finite value of a level/metric for a given year across
all selected runs.  At the basin level that is exactly the across-GCM spread; at
nested levels it also includes the spread across sub-basins / water bodies /
sites, so the sampling-point panels show the site-level 10-90% range.
"""

const CLIMATE_LEVELS = (:site, :subcatchment, :water_body, :basin)

const LEVEL_LABELS = Dict(
    :site => "Sampling point (CODIGO)",
    :subcatchment => "Sub-basin (CODIGO_S)",
    :water_body => "Water body (ID_masa)",
    :basin => "Whole basin (ES050)",
)

const LEVEL_FILES = Dict(
    :site => "level_sampling_point.csv",
    :subcatchment => "level_subcatchment.csv",
    :water_body => "level_water_body.csv",
    :basin => "level_basin.csv",
)

# Column that identifies the units inside each reporting level.
const LEVEL_GROUP_COLUMNS = Dict(
    :site => :CODIGO,
    :subcatchment => :subcatchment,
    :water_body => :water_body,
    :basin => :basin,
)

# The figure metrics are the subset of the canonical viewer metric definitions
# (`VIEWER_METRICS` in src/outputs.jl) rendered by the climate figures, so labels
# and units cannot drift between the viewer export and the figures.
const CLIMATE_METRIC_KEYS = (:native_richness, :invasive_richness,
    :native_extinction_risk, :total_biomass, :temperature_c)

const CLIMATE_METRICS = [definition for definition in VIEWER_METRICS
                         if definition[1] in CLIMATE_METRIC_KEYS]

const SCENARIO_COLORS = Dict(
    "ssp126" => :dodgerblue,
    "ssp245" => :seagreen,
    "ssp370" => :darkorange,
    "ssp585" => :firebrick,
)

scenario_color(scenario::AbstractString) = get(SCENARIO_COLORS, String(scenario), :gray30)

function _check_level(level::Symbol)
    level in CLIMATE_LEVELS || error("unknown climate level: $level (expected one of $CLIMATE_LEVELS)")
    return level
end

"""
    climate_metric_column(level, metric)

Column name that holds `metric` in the CSV of `level`.  Sampling-point (`:site`)
tables store the raw per-site metric; the aggregate levels prefix it with
`mean_`.
"""
function climate_metric_column(level::Symbol, metric::Symbol)
    _check_level(level)
    return level == :site ? metric : Symbol("mean_", metric)
end

# =============================================================================
# --- Run container ---
# =============================================================================

"""
    ClimateRun

One discovered run: identity, year range, completeness and the loaded level
tables (keyed by the level symbols in [`CLIMATE_LEVELS`](@ref)).
"""
mutable struct ClimateRun
    scenario::String
    gcm::String
    run_dir::String
    years::Vector{Int}
    start_year::Union{Int,Nothing}
    end_year::Union{Int,Nothing}
    complete::Bool
    tables::Dict{Symbol,Any}
end

# =============================================================================
# --- Table loading ---
# =============================================================================

function _needed_columns(level::Symbol)
    cols = Symbol[:year]
    level == :site || push!(cols, LEVEL_GROUP_COLUMNS[level])
    for (metric, _, _) in CLIMATE_METRICS
        push!(cols, climate_metric_column(level, metric))
    end
    return unique(cols)
end

"""
    load_level_table(path, level)

Read one level CSV, keeping only the columns used for figures.  Returns `nothing`
when the file does not exist.  Falls back to a full read if the requested subset
is unavailable (older or partial exports).
"""
function load_level_table(path::AbstractString, level::Symbol)
    isfile(path) || return nothing
    try
        return CSV.read(path, DataFrame; select=_needed_columns(level))
    catch
        return CSV.read(path, DataFrame)
    end
end

"""
    load_run_level_tables(run_dir)

Load the four level tables for one run directory.  Missing levels are `nothing`.
"""
function load_run_level_tables(run_dir::AbstractString)
    tables = Dict{Symbol,Any}()
    for level in CLIMATE_LEVELS
        path = joinpath(run_dir, "export", "levels", LEVEL_FILES[level])
        tables[level] = load_level_table(path, level)
    end
    return tables
end

function _table_years(df)
    (df === nothing || !hasproperty(df, :year) || nrow(df) == 0) && return Int[]
    years = Int[]
    for y in df.year
        ismissing(y) || push!(years, Int(y))
    end
    return sort(unique(years))
end

function _level_years(run::ClimateRun, level::Symbol)
    years = _table_years(get(run.tables, level, nothing))
    return isempty(years) ? run.years : years
end

# =============================================================================
# --- Discovery, metadata and completeness ---
# =============================================================================

function _int_or_nothing(value)
    (value === nothing || ismissing(value)) && return nothing
    try
        return Int(value)
    catch
        return nothing
    end
end

"""
    run_year_range(df)
    run_year_range(path)

`(minimum_year, maximum_year)` of a level table, or `nothing` when unavailable.
"""
function run_year_range(df::DataFrame)
    (hasproperty(df, :year) && nrow(df) > 0) || return nothing
    years = Int[]
    for y in df.year
        ismissing(y) || push!(years, Int(y))
    end
    isempty(years) && return nothing
    return (minimum(years), maximum(years))
end

function run_year_range(path::AbstractString)
    isfile(path) || return nothing
    return run_year_range(CSV.read(path, DataFrame))
end

"""
    read_run_metadata(run_dir)

Minimal, dependency-free reader for the four `run_metadata.json` fields used by
the figure pipeline (`scenario`, `gcm`, `start_year`, `end_year`).
"""
function read_run_metadata(run_dir::AbstractString)
    path = joinpath(run_dir, "export", "run_metadata.json")
    isfile(path) || return Dict{String,Any}()
    text = read(path, String)
    meta = Dict{String,Any}()
    for key in ("scenario", "gcm")
        m = match(Regex("\"" * key * "\"\\s*:\\s*\"([^\"]*)\""), text)
        m === nothing || (meta[key] = String(m.captures[1]))
    end
    for key in ("start_year", "end_year")
        m = match(Regex("\"" * key * "\"\\s*:\\s*(-?\\d+)"), text)
        m === nothing || (meta[key] = parse(Int, m.captures[1]))
    end
    return meta
end

function _subdirs(path::AbstractString)
    isdir(path) || return String[]
    return sort([joinpath(path, name) for name in readdir(path)
                 if isdir(joinpath(path, name))])
end

function _resolve_run_dir(results_root::AbstractString, scenario::AbstractString,
        gcm::AbstractString, stored::AbstractString)
    if !isempty(stored) && isdir(stored)
        return stored
    end
    if !isempty(stored)
        alt = replace(stored, '\\' => '/')
        isdir(alt) && return alt
    end
    return joinpath(results_root, scenario, gcm)
end

function _index_entries(results_root::AbstractString)
    path = joinpath(results_root, "runs_index.csv")
    isfile(path) || return NamedTuple[]
    df = CSV.read(path, DataFrame)
    hasproperty(df, :scenario) || return NamedTuple[]
    entries = NamedTuple[]
    for row in eachrow(df)
        scenario = string(row.scenario)
        gcm = hasproperty(df, :gcm) ? string(row.gcm) : ""
        run_dir = _resolve_run_dir(results_root, scenario, gcm,
            hasproperty(df, :run_dir) ? string(row.run_dir) : "")
        start_year = hasproperty(df, :start_year) ? _int_or_nothing(row.start_year) : nothing
        end_year = hasproperty(df, :end_year) ? _int_or_nothing(row.end_year) : nothing
        push!(entries, (scenario=scenario, gcm=gcm, run_dir=run_dir,
            start_year=start_year, end_year=end_year))
    end
    return entries
end

function _scan_entries(results_root::AbstractString)
    entries = NamedTuple[]
    for scenario_dir in _subdirs(results_root)
        scenario = basename(scenario_dir)
        scenario == "figures" && continue
        for gcm_dir in _subdirs(scenario_dir)
            basin_path = joinpath(gcm_dir, "export", "levels", "level_basin.csv")
            isfile(basin_path) || continue
            meta = read_run_metadata(gcm_dir)
            rng = run_year_range(basin_path)
            push!(entries, (scenario=get(meta, "scenario", scenario),
                gcm=get(meta, "gcm", basename(gcm_dir)), run_dir=gcm_dir,
                start_year=rng === nothing ? get(meta, "start_year", nothing) : rng[1],
                end_year=rng === nothing ? get(meta, "end_year", nothing) : rng[2]))
        end
    end
    return entries
end

function _is_complete(start_year::Union{Int,Nothing}, end_year::Union{Int,Nothing},
        reference_start::Union{Int,Nothing}, reference_end::Union{Int,Nothing})
    if reference_start !== nothing && start_year !== nothing && start_year > reference_start
        return false
    end
    if reference_end !== nothing && end_year !== nothing && end_year < reference_end
        return false
    end
    return true
end

"""
    discover_climate_runs(results_root; expected_end_year=nothing)

Discover runs under `results_root`, preferring `runs_index.csv` and falling back
to scanning the tree for `export/levels/level_basin.csv`.  Returns
`(runs, expected_end_year)`; each [`ClimateRun`](@ref) carries its loaded tables
and a `complete` flag (a run is incomplete when its year range does not cover the
reference horizon).
"""
function discover_climate_runs(results_root::AbstractString;
        expected_end_year::Union{Int,Nothing}=nothing)
    entries = _index_entries(results_root)
    isempty(entries) && (entries = _scan_entries(results_root))

    runs = ClimateRun[]
    starts = Int[]
    ends = Int[]
    for entry in entries
        tables = load_run_level_tables(entry.run_dir)
        basin_years = _table_years(get(tables, :basin, nothing))
        start_year = entry.start_year
        end_year = entry.end_year
        if !isempty(basin_years)
            start_year === nothing && (start_year = first(basin_years))
            end_year === nothing && (end_year = last(basin_years))
        end
        years = isempty(basin_years) ?
            (start_year === nothing || end_year === nothing ? Int[] : collect(start_year:end_year)) :
            basin_years
        start_year === nothing || push!(starts, start_year)
        end_year === nothing || push!(ends, end_year)
        push!(runs, ClimateRun(String(entry.scenario), String(entry.gcm),
            String(entry.run_dir), years, start_year, end_year, false, tables))
    end

    resolved_end = expected_end_year === nothing ?
        (isempty(ends) ? nothing : maximum(ends)) : expected_end_year
    reference_start = isempty(starts) ? nothing : minimum(starts)

    for run in runs
        run.complete = _is_complete(run.start_year, run.end_year,
            reference_start, resolved_end)
    end
    sort!(runs; by=run -> (run.scenario, run.gcm))
    return runs, resolved_end
end

"""
    load_climate_run(run_dir; scenario="", gcm="", expected_end_year=nothing)

Load a single run directory into a [`ClimateRun`](@ref) (used by the CLI and by
callers that already know the run folder).
"""
function load_climate_run(run_dir::AbstractString; scenario::AbstractString="",
        gcm::AbstractString="", expected_end_year::Union{Int,Nothing}=nothing)
    meta = read_run_metadata(run_dir)
    resolved_scenario = isempty(scenario) ? get(meta, "scenario", basename(dirname(run_dir))) : scenario
    resolved_gcm = isempty(gcm) ? get(meta, "gcm", basename(run_dir)) : gcm
    tables = load_run_level_tables(run_dir)
    basin_years = _table_years(get(tables, :basin, nothing))
    start_year = isempty(basin_years) ? get(meta, "start_year", nothing) : first(basin_years)
    end_year = isempty(basin_years) ? get(meta, "end_year", nothing) : last(basin_years)
    years = isempty(basin_years) ?
        (start_year === nothing || end_year === nothing ? Int[] : collect(start_year:end_year)) :
        basin_years
    complete = _is_complete(start_year, end_year, nothing, expected_end_year)
    return ClimateRun(String(resolved_scenario), String(resolved_gcm), String(run_dir),
        years, start_year, end_year, complete, tables)
end

"""
    incomplete_runs(runs)

Discovered runs whose year range does not cover the reference horizon.
"""
incomplete_runs(runs::AbstractVector{<:ClimateRun}) = [run for run in runs if !run.complete]

# =============================================================================
# --- Statistics (no rendering) ---
# =============================================================================

function _collect_year_values(df, col::Symbol)
    out = Dict{Int,Vector{Float64}}()
    (df === nothing || !hasproperty(df, :year) || !hasproperty(df, col)) && return out
    years = df.year
    values = df[!, col]
    for i in 1:nrow(df)
        y = years[i]
        ismissing(y) && continue
        v = values[i]
        (v isa Real && isfinite(float(v))) || continue
        push!(get!(out, Int(y), Float64[]), Float64(v))
    end
    return out
end

"""
    run_level_stats(run, level, metric)

Per-year statistics of `metric` within a single run and level.  `mean` averages
every finite unit value (sites, sub-basins, water bodies or the single basin);
`p10`/`p50`/`p90` describe the across-unit spread.  Missing metrics yield
`NaN` entries rather than an error.
"""
function run_level_stats(run::ClimateRun, level::Symbol, metric::Symbol)
    _check_level(level)
    years = _level_years(run, level)
    col = climate_metric_column(level, metric)
    year_values = _collect_year_values(get(run.tables, level, nothing), col)
    mean_ = Float64[]
    p10 = Float64[]
    p50 = Float64[]
    p90 = Float64[]
    for year in years
        values = get(year_values, year, Float64[])
        if isempty(values)
            push!(mean_, NaN); push!(p10, NaN); push!(p50, NaN); push!(p90, NaN)
        else
            push!(mean_, mean(values))
            push!(p10, quantile(values, 0.1))
            push!(p50, quantile(values, 0.5))
            push!(p90, quantile(values, 0.9))
        end
    end
    return (years=years, mean=mean_, p10=p10, p50=p50, p90=p90)
end

"""
    ensemble_level_stats(runs, level, metric)

Pool the finite values of every supplied run (and every unit for nested levels)
per year and return `median` plus `p10`/`p25`/`p75`/`p90` bands and the number of
contributing runs (`n_models`).
"""
function ensemble_level_stats(runs::AbstractVector, level::Symbol, metric::Symbol)
    _check_level(level)
    col = climate_metric_column(level, metric)
    pooled = Dict{Int,Vector{Float64}}()
    models = Dict{Int,Set{Int}}()
    yearset = Set{Int}()
    for (i, run) in enumerate(runs)
        for year in _level_years(run, level)
            push!(yearset, year)
        end
        for (year, values) in _collect_year_values(get(run.tables, level, nothing), col)
            isempty(values) && continue
            append!(get!(pooled, year, Float64[]), values)
            push!(get!(models, year, Set{Int}()), i)
        end
    end
    years = sort(collect(yearset))
    quantiles = Dict{Float64,Vector{Float64}}(p => Float64[] for p in (0.1, 0.25, 0.5, 0.75, 0.9))
    n_models = Int[]
    for year in years
        values = get(pooled, year, Float64[])
        if isempty(values)
            for p in keys(quantiles)
                push!(quantiles[p], NaN)
            end
            push!(n_models, 0)
        else
            for p in keys(quantiles)
                push!(quantiles[p], quantile(values, p))
            end
            push!(n_models, length(get(models, year, Set{Int}())))
        end
    end
    return (years=years, median=quantiles[0.5], p10=quantiles[0.1],
        p25=quantiles[0.25], p75=quantiles[0.75], p90=quantiles[0.9],
        n_models=n_models)
end

# =============================================================================
# --- Rendering ---
# =============================================================================

function _save_climate_figure(fig, path::AbstractString; px_per_unit::Real=2.0)
    mkpath(dirname(path))
    Makie.save(path, fig; px_per_unit=px_per_unit)
    println("Figure saved to: $path")
    return path
end

# Render one figure without letting a single failure abort the rest of the
# inventory; the failed path is reported and skipped.
function _safe_figure(render, path::AbstractString)
    try
        return render()
    catch err
        @warn "climate figure failed; continuing" path exception=(err, catch_backtrace())
        return nothing
    end
end

"""
    plot_climate_run_figure(run, output_path; incomplete=false)

One multi-panel overview per run: rows are the five reporting metrics, columns
are the four GuadeX levels.  Nested levels show the across-unit median with a
10-90% band (the site spread at the sampling-point level).
"""
function plot_climate_run_figure(run::ClimateRun, output_path::AbstractString;
        incomplete::Bool=false)
    fig = Figure(size=(1750, 2150))
    drew = false
    for (ri, (metric, label, unit)) in enumerate(CLIMATE_METRICS)
        for (ci, level) in enumerate(CLIMATE_LEVELS)
            ax = Axis(fig[ri, ci];
                title="$label — $(LEVEL_LABELS[level])", titlesize=10,
                xlabel=ri == length(CLIMATE_METRICS) ? "Year" : "",
                ylabel=ci == 1 ? "$label ($unit)" : "",
                xlabelsize=10, ylabelsize=9)
            stats = run_level_stats(run, level, metric)
            if !isempty(stats.years) && !all(isnan, stats.mean)
                drew = true
                if level != :basin
                    band!(ax, stats.years, stats.p10, stats.p90; color=(:steelblue, 0.20))
                    lines!(ax, stats.years, stats.p10; color=(:steelblue, 0.55), linewidth=0.8)
                    lines!(ax, stats.years, stats.p90; color=(:steelblue, 0.55), linewidth=0.8)
                end
                lines!(ax, stats.years, stats.mean; color=:black, linewidth=2)
            end
        end
    end
    drew || return nothing
    title = "GuadeX climate run — $(run.scenario) / $(run.gcm)"
    incomplete && (title *= "  [INCOMPLETE]")
    Label(fig[0, :], title, fontsize=18, font=:bold)
    return _save_climate_figure(fig, output_path)
end

"""
    plot_ensemble_level_figure(scenario, level, runs, output_path)

Ensemble figure for one scenario and level: five metric panels with the median
across the supplied GCMs and 25-75 / 10-90 bands.
"""
function plot_ensemble_level_figure(scenario::AbstractString, level::Symbol,
        runs::AbstractVector{<:ClimateRun}, output_path::AbstractString)
    _check_level(level)
    fig = Figure(size=(1950, 520))
    drew = false
    for (ci, (metric, label, unit)) in enumerate(CLIMATE_METRICS)
        ax = Axis(fig[1, ci]; title="$label ($unit)", titlesize=12,
            xlabel="Year", ylabel=ci == 1 ? LEVEL_LABELS[level] : "")
        stats = ensemble_level_stats(runs, level, metric)
        if !isempty(stats.years) && !all(isnan, stats.median)
            drew = true
            band!(ax, stats.years, stats.p10, stats.p90; color=(:steelblue, 0.18))
            band!(ax, stats.years, stats.p25, stats.p75; color=(:steelblue, 0.32))
            lines!(ax, stats.years, stats.median; color=:black, linewidth=2.2)
        end
    end
    drew || return nothing
    Label(fig[0, :],
        "Ensemble $scenario — $(LEVEL_LABELS[level]) ($(length(runs)) GCMs)\n" *
        "bands: 10-90% and 25-75% across GCMs (and sites/units at nested levels)",
        fontsize=14, font=:bold, justification=:center, halign=:center)
    return _save_climate_figure(fig, output_path)
end

"""
    plot_across_scenarios_figure(level, scenario_runs, output_path)

Across-scenario comparison for one level: five metric panels, one ensemble
median line per scenario with a light 10-90% band.
"""
function plot_across_scenarios_figure(level::Symbol,
        scenario_runs::AbstractDict, output_path::AbstractString)
    _check_level(level)
    scenarios = sort(collect(keys(scenario_runs)))
    fig = Figure(size=(2000, 520))
    drew = false
    legend_ax = nothing
    for (ci, (metric, label, unit)) in enumerate(CLIMATE_METRICS)
        ax = Axis(fig[1, ci]; title="$label ($unit)", titlesize=12,
            xlabel="Year", ylabel=ci == 1 ? LEVEL_LABELS[level] : "")
        for scenario in scenarios
            stats = ensemble_level_stats(scenario_runs[scenario], level, metric)
            (isempty(stats.years) || all(isnan, stats.median)) && continue
            drew = true
            color = scenario_color(scenario)
            band!(ax, stats.years, stats.p10, stats.p90; color=(color, 0.10))
            lines!(ax, stats.years, stats.median; color=color, linewidth=2.2, label=scenario)
        end
        ci == 1 && (legend_ax = ax)
    end
    drew || return nothing
    legend_ax === nothing || axislegend(legend_ax; position=:rb, nbanks=2, fontsize=10)
    Label(fig[0, :], "Across-scenario comparison — $(LEVEL_LABELS[level])",
        fontsize=15, font=:bold)
    return _save_climate_figure(fig, output_path)
end

"""
    plot_final_year_summary(scenario_runs, output_path; year)

Final-year summary: rows are levels, columns are metrics; every scenario is a
point (median) with 25-75 (thick) and 10-90 (thin) spreads across GCMs/units.
"""
function plot_final_year_summary(scenario_runs::AbstractDict,
        output_path::AbstractString; year::Int)
    scenarios = sort(collect(keys(scenario_runs)))
    isempty(scenarios) && return nothing
    fig = Figure(size=(440 * length(CLIMATE_METRICS) + 140,
        310 * length(CLIMATE_LEVELS) + 70))
    drew = false
    for (ri, level) in enumerate(CLIMATE_LEVELS)
        for (ci, (metric, label, unit)) in enumerate(CLIMATE_METRICS)
            ax = Axis(fig[ri, ci]; title="$label ($unit)", titlesize=10,
                xlabel=ri == length(CLIMATE_LEVELS) ? "Scenario" : "",
                ylabel=ci == 1 ? LEVEL_LABELS[level] : "",
                xticks=(1:length(scenarios), scenarios), xticklabelrotation=0.25)
            for (si, scenario) in enumerate(scenarios)
                stats = ensemble_level_stats(scenario_runs[scenario], level, metric)
                idx = findfirst(==(year), stats.years)
                (idx === nothing || isnan(stats.median[idx])) && continue
                drew = true
                color = scenario_color(scenario)
                lines!(ax, [si, si], [stats.p10[idx], stats.p90[idx]];
                    color=color, linewidth=1.2)
                lines!(ax, [si, si], [stats.p25[idx], stats.p75[idx]];
                    color=color, linewidth=4)
                scatter!(ax, [si], [stats.median[idx]]; color=color,
                    marker=:circle, markersize=8)
            end
        end
    end
    drew || return nothing
    Label(fig[0, :],
        "Year $year summary — median with 25-75% and 10-90% spreads " *
        "(GCMs, and sites/units at nested levels)",
        fontsize=14, font=:bold)
    return _save_climate_figure(fig, output_path)
end

# =============================================================================
# --- Orchestration ---
# =============================================================================

"""
    plot_climate_scenario_figures(results_root; figures_dir, expected_end_year)

Generate the full climate figure inventory from already-exported runs:

* `per_run/<scenario>__<gcm>.png` - five metrics x four levels overview.
* `ensemble/ensemble_<scenario>_<level>.png` - across-GCM ensemble per scenario.
* `comparison/across_scenarios_<level>.png` - across-scenario comparison.
* `summary_<year>.png` - final-year summary by level and metric.
* `figure_inventory.csv` - manifest of everything written.

Incomplete runs (for example a 2-year smoke run) are flagged in their per-run
title and excluded from every ensemble/summary.  Returns `figures_dir`.
"""
function plot_climate_scenario_figures(results_root::AbstractString=joinpath("results", "climate_scenarios");
        figures_dir::AbstractString=joinpath(results_root, "figures"),
        expected_end_year::Union{Int,Nothing}=nothing)
    runs, expected = discover_climate_runs(results_root;
        expected_end_year=expected_end_year)
    if isempty(runs)
        @warn "no climate runs found under $results_root; no figures written"
        return figures_dir
    end
    mkpath(figures_dir)

    complete = ClimateRun[run for run in runs if run.complete]
    for run in incomplete_runs(runs)
        @warn "incomplete climate run excluded from ensembles" scenario=run.scenario gcm=run.gcm years="$(run.start_year)-$(run.end_year)"
    end

    rows = NamedTuple[]

    per_run_dir = joinpath(figures_dir, "per_run")
    for run in runs
        path = joinpath(per_run_dir, "$(run.scenario)__$(run.gcm).png")
        saved = _safe_figure(() -> plot_climate_run_figure(run, path;
            incomplete=!run.complete), path)
        saved === nothing && continue
        push!(rows, (category="per_run", scenario=run.scenario, gcm=run.gcm,
            metric="all", level="all", path=saved,
            status=run.complete ? "complete" : "incomplete"))
    end

    by_scenario = Dict{String,Vector{ClimateRun}}()
    for run in complete
        push!(get!(by_scenario, run.scenario, ClimateRun[]), run)
    end

    ensemble_dir = joinpath(figures_dir, "ensemble")
    for scenario in sort(collect(keys(by_scenario)))
        for level in CLIMATE_LEVELS
            path = joinpath(ensemble_dir, "ensemble_$(scenario)_$(level).png")
            saved = _safe_figure(() -> plot_ensemble_level_figure(scenario, level,
                by_scenario[scenario], path), path)
            saved === nothing && continue
            push!(rows, (category="ensemble", scenario=scenario, gcm="",
                metric="all", level=string(level), path=saved, status="complete"))
        end
    end

    comparison_dir = joinpath(figures_dir, "comparison")
    for level in CLIMATE_LEVELS
        path = joinpath(comparison_dir, "across_scenarios_$(level).png")
        saved = _safe_figure(() -> plot_across_scenarios_figure(level, by_scenario, path), path)
        saved === nothing && continue
        push!(rows, (category="comparison", scenario="all", gcm="",
            metric="all", level=string(level), path=saved, status="complete"))
    end

    if expected !== nothing && !isempty(by_scenario)
        path = joinpath(figures_dir, "summary_$(expected).png")
        saved = _safe_figure(() -> plot_final_year_summary(by_scenario, path;
            year=Int(expected)), path)
        if saved !== nothing
            push!(rows, (category="summary", scenario="all", gcm="", metric="all",
                level="all", path=saved, status="complete"))
        end
    end

    inventory_path = joinpath(figures_dir, "figure_inventory.csv")
    if isempty(rows)
        CSV.write(inventory_path, DataFrame(category=String[], scenario=String[],
            gcm=String[], metric=String[], level=String[], path=String[], status=String[]))
    else
        CSV.write(inventory_path, DataFrame(rows))
    end

    println("\n[climate_figures] $(length(rows)) figure(s) written to $figures_dir")
    println("[climate_figures] inventory: $inventory_path")
    return figures_dir
end
