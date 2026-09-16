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
    native_extinction_risk::Vector{Float64}
end

const _BASIN_SERIES_FIELDS = [
    (:delta_temperature_c, :mean_delta_temperature_c),
    (:native_richness, :mean_native_richness),
    (:total_richness, :mean_total_richness),
    (:invasive_richness, :mean_invasive_richness),
    (:native_biomass, :mean_native_biomass),
    (:total_biomass, :mean_total_biomass),
    (:native_extinction_risk, :mean_native_extinction_risk),
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
            columns[:native_extinction_risk]))
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
# --- Forcing vs response ---
# =============================================================================

"""
    plot_climate_forcing_response(runs, output_path)

Six-panel diagnostic: warming trajectories, basin richness and biomass
trajectories, the final-year warming-response scatter and the (near-zero)
extinction-risk trajectory.  Richness scenarios are distinguished by the
canonical [`scenario_color`](@ref).
"""
function plot_climate_forcing_response(runs::AbstractVector{<:ClimateBasinSeries},
        output_path::AbstractString)
    isempty(runs) && return nothing
    scenarios = sort(unique(run.scenario for run in runs))
    warming = _series_envelope(runs, :delta_temperature_c, scenarios)
    richness = _series_envelope(runs, :native_richness, scenarios)
    biomass = _series_envelope(runs, :total_biomass, scenarios)
    risk = _series_envelope(runs, :native_extinction_risk, scenarios)

    fig = Figure(size=(1680, 1000))

    function _trajectory_axis(row, col, title, ylabel, envelope)
        ax = Axis(fig[row, col]; title=title, xlabel="Year", ylabel=ylabel,
            titlesize=12)
        for scenario in scenarios
            env = envelope[scenario]
            (isempty(env.years) || all(isnan, env.median)) && continue
            color = scenario_color(scenario)
            band!(ax, env.years, env.lo, env.hi; color=(color, 0.12))
            lines!(ax, env.years, env.median; color=color, linewidth=2.2, label=scenario)
        end
        return ax
    end

    ax1 = _trajectory_axis(1, 1, "Forcing: water-temperature change", "ΔT (°C)", warming)
    ax2 = _trajectory_axis(1, 2, "Response: mean native richness", "species / site", richness)
    ax3 = _trajectory_axis(1, 3, "Response: mean total biomass", "biomass / site", biomass)
    axislegend(ax1; position=:lt, nbanks=2, fontsize=9)
    # The four SSP trajectories are almost identical in the response panels; say
    # so explicitly, otherwise it looks like only one scenario was plotted.
    for ax in (ax2, ax3)
        text!(ax, 0.04, 0.06; space=:relative, text="all SSPs overlap",
            align=(:left, :bottom), fontsize=11, color=(:black, 0.65))
    end

    # Final-year warming response across all runs.
    final_warming = [_final_value(run, :delta_temperature_c) for run in runs]
    final_richness = [_final_value(run, :native_richness) for run in runs]
    final_biomass = [_final_value(run, :total_biomass) for run in runs]

    ax4 = Axis(fig[2, 1]; title="Final-year native richness vs warming",
        xlabel="ΔT by end of horizon (°C)", ylabel="mean native richness / site",
        titlesize=12)
    for scenario in scenarios
        idx = findall(==(scenario), [run.scenario for run in runs])
        xs = final_warming[idx]
        ys = final_richness[idx]
        x, y = _pairwise_finite(xs, ys)
        isempty(x) && continue
        scatter!(ax4, x, y; color=(scenario_color(scenario), 0.85), markersize=9,
            label=scenario)
    end
    xr, yr = _pairwise_finite(final_warming, final_richness)
    r_rich = _pearson(xr, yr)
    isfinite(r_rich) && text!(ax4, 0.04, 0.94; space=:relative,
        text="Pearson r = $(round(r_rich, digits=3))", align=(:left, :top), fontsize=11)

    ax5 = Axis(fig[2, 2]; title="Final-year total biomass vs warming",
        xlabel="ΔT by end of horizon (°C)", ylabel="mean total biomass / site",
        titlesize=12)
    for scenario in scenarios
        idx = findall(==(scenario), [run.scenario for run in runs])
        x, y = _pairwise_finite(final_warming[idx], final_biomass[idx])
        isempty(x) && continue
        scatter!(ax5, x, y; color=(scenario_color(scenario), 0.85), markersize=9)
    end
    xb, yb = _pairwise_finite(final_warming, final_biomass)
    r_bio = _pearson(xb, yb)
    isfinite(r_bio) && text!(ax5, 0.04, 0.94; space=:relative,
        text="Pearson r = $(round(r_bio, digits=3))", align=(:left, :top), fontsize=11)

    ax6 = Axis(fig[2, 3]; title="Extinction-risk metric (1 - richness / initial richness)",
        xlabel="Year", ylabel="mean relative richness loss", titlesize=12)
    for scenario in scenarios
        env = risk[scenario]
        (isempty(env.years) || all(isnan, env.median)) && continue
        lines!(ax6, env.years, env.median; color=scenario_color(scenario), linewidth=2.2)
    end
    risk_hi = 0.0
    for scenario in scenarios
        for value in risk[scenario].hi
            isfinite(value) && (risk_hi = max(risk_hi, value))
        end
    end
    ylims!(ax6, 0.0, max(0.05, risk_hi))

    Label(fig[0, :],
        "Climate diagnostics — forcing vs community response (median across " *
        "$(length(runs)) runs; bands are 10-90% across GCMs)",
        fontsize=15, font=:bold)
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
    thermal = _read_species_thermal_table(species_chars_file)
    site_temps = _read_site_temperatures_baseline(site_table_path)

    fig = Figure(size=(1150, 720))
    ax = Axis(fig[1, 1];
        title="Native thermal niches vs site temperature",
        xlabel="Water temperature (°C)", ylabel="Thermal suitability / relative frequency",
        titlesize=13)

    t_grid = range(minimum(site_temps; init=8.0) - 1, maximum(site_temps; init=20.0) + 3;
        length=300)

    # Site-temperature distribution as translucent bars behind the curves.
    if !isempty(site_temps)
        lo, hi = 8.0, 26.0
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
    for (i, (key, codes)) in enumerate(sort(collect(groups); by=first))
        optimum, sigma = key
        curve = [exp(-(t - optimum)^2 / (2 * sigma^2)) for t in t_grid]
        label = length(codes) == 1 ? "$(only(codes)): opt $(round(optimum, digits=1))" :
            "opt $(round(optimum, digits=1)) ×$(length(codes)) ($(join(codes, ",")))"
        lines!(ax, t_grid, curve; color=palette[mod1(i, length(palette))],
            linewidth=2.2, label=label)
    end

    # Context species that drive some local declines but are not in the metric.
    for (i, code) in enumerate(extra_codes)
        key = lowercase(string(code))
        haskey(thermal, key) || continue
        p = thermal[key]
        curve = [exp(-(t - p.optimum)^2 / (2 * p.sigma^2)) for t in t_grid]
        lines!(ax, t_grid, curve; color=(:black, 0.55), linewidth=2.0,
            linestyle=:dash, label="$(code) (not in native metric)")
    end

    if !isempty(site_temps)
        baseline = mean(site_temps)
        vlines!(ax, [baseline]; color=:black, linewidth=1.6, linestyle=:dot)
        vlines!(ax, [baseline + warming_delta]; color=:crimson, linewidth=1.6, linestyle=:dot)
        text!(ax, baseline - 0.1, 1.04; text="baseline $(round(baseline, digits=1))°C",
            align=(:right, :bottom), fontsize=10)
        text!(ax, baseline + warming_delta + 0.1, 1.04;
            text="+$(round(warming_delta, digits=2))°C → $(round(baseline + warming_delta, digits=1))°C",
            align=(:left, :bottom), fontsize=10, color=:crimson)
    end
    ylims!(ax, 0.0, 1.12)
    axislegend(ax; position=:rt, nbanks=2, fontsize=8, framevisible=true)

    Label(fig[0, :],
        "Thermal niches explain the warming response: most richness species optima lie above current water temperatures",
        fontsize=13, font=:bold)
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
