"""
    Output module: four-level reporting tables and 3-D viewer exports.

This module turns a solved population trajectory into the four reporting levels
requested by the GuadeX team (docs/Reporting/README.NewSimus_Sept26.md) and into
the JSON/CSV shapes consumed by the `viz/` 3-D explorer (viz/README.md).

Levels
------
1. sampling point : `CODIGO`
2. sub-basin      : `CODIGO_S` (the model subcatchment)
3. water body     : `ID_masa` (EUMASCod, resolved through
                    `data/site_waterbody_crosswalk.csv`)
4. whole basin    : `ES050` (all modelled sites)

Outputs written by [`export_run_outputs`](@ref):

    <run_dir>/levels/level_sampling_point.csv
    <run_dir>/levels/level_subcatchment.csv
    <run_dir>/levels/level_water_body.csv
    <run_dir>/levels/level_basin.csv
    <run_dir>/viewer/guadex_results_timeseries.csv
    <run_dir>/viewer/guadex_results_metrics.json
    <run_dir>/viewer/guadex_results_<metric>.json
    <run_dir>/viewer/level_subcatchment_timeseries.json
    <run_dir>/viewer/level_water_body_timeseries.json
    <run_dir>/viewer/level_basin_timeseries.json
    <run_dir>/run_metadata.json

The viewer files are keyed exactly as `viz/src/data/ResultsModel.js` expects:
`{name, unit, steps, data}` for time series and `{name, data}` for per-site
metrics, with site keys equal to `CODIGO`.
"""

const GUADEX_BASIN_ID = "ES050"
const GUADEX_UNASSIGNED_WATER_BODY = "UNASSIGNED"
const GUADEX_UNKNOWN_SUBCATCHMENT = "UNKNOWN"

# =============================================================================
# --- Minimal JSON serialisation (no extra dependency) ---
# =============================================================================

_json_escape(s::AbstractString) = replace(String(s), "\\" => "\\\\", "\"" => "\\\"",
    "\n" => "\\n", "\r" => "\\r", "\t" => "\\t")

function _json_number(x::Real)
    if isfinite(x)
        return isinteger(x) ? string(Int(x)) : string(Float64(x))
    end
    return "null"
end

_json_value(x::AbstractString) = "\"" * _json_escape(x) * "\""
_json_value(x::Real) = _json_number(x)
_json_value(x::Bool) = x ? "true" : "false"
_json_value(::Missing) = "null"
_json_value(::Nothing) = "null"

function _write_timeseries_json(path::AbstractString, name::AbstractString, unit::AbstractString,
                                steps::AbstractVector, data::AbstractDict{<:AbstractString,<:AbstractVector})
    io = IOBuffer()
    print(io, "{\n  \"name\": ", _json_value(name))
    print(io, ",\n  \"unit\": ", _json_value(unit))
    print(io, ",\n  \"steps\": [", join((_json_value(string(s)) for s in steps), ", "), "]")
    print(io, ",\n  \"data\": {")
    first_entry = true
    for key in sort(collect(keys(data)))
        first_entry || print(io, ",")
        first_entry = false
        print(io, "\n    ", _json_value(key), ": [",
            join((_json_value(v) for v in data[key]), ", "), "]")
    end
    print(io, "\n  }\n}\n")
    open(path, "w") do file
        write(file, take!(io))
    end
    return path
end

function _write_metrics_json(path::AbstractString, name::AbstractString,
                             data::AbstractDict{<:AbstractString,<:AbstractDict})
    io = IOBuffer()
    print(io, "{\n  \"name\": ", _json_value(name))
    print(io, ",\n  \"data\": {")
    first_entry = true
    for key in sort(collect(keys(data)))
        first_entry || print(io, ",")
        first_entry = false
        print(io, "\n    ", _json_value(key), ": {")
        first_field = true
        for field in sort(collect(keys(data[key])))
            first_field || print(io, ", ")
            first_field = false
            print(io, _json_value(string(field)), ": ", _json_value(data[key][field]))
        end
        print(io, "}")
    end
    print(io, "\n  }\n}\n")
    open(path, "w") do file
        write(file, take!(io))
    end
    return path
end

function _write_value_json(path::AbstractString, value)
    io = IOBuffer()
    _write_json_value(io, value, 0)
    print(io, "\n")
    open(path, "w") do file
        write(file, take!(io))
    end
    return path
end

function _write_json_value(io::IO, value, indent::Int)
    pad = "  "^indent
    if value isa AbstractDict
        print(io, "{")
        keys_sorted = sort(collect(keys(value)), by=string)
        for (i, key) in enumerate(keys_sorted)
            i == 1 ? print(io, "\n") : print(io, ",\n")
            print(io, pad, "  ", _json_value(string(key)), ": ")
            _write_json_value(io, value[key], indent + 1)
        end
        isempty(keys_sorted) ? print(io, "}") : print(io, "\n", pad, "}")
    elseif value isa AbstractVector
        print(io, "[")
        for (i, item) in enumerate(value)
            i == 1 ? print(io, "\n") : print(io, ",\n")
            print(io, pad, "  ")
            _write_json_value(io, item, indent + 1)
        end
        isempty(value) ? print(io, "]") : print(io, "\n", pad, "]")
    elseif value isa NamedTuple || value isa Tuple
        _write_json_value(io, Dict(string(k) => v for (k, v) in pairs(value)), indent)
    else
        print(io, _json_value(value))
    end
end

# =============================================================================
# --- Site -> level crosswalk ---
# =============================================================================

const CROSSWALK_COLUMNS = ["CODIGO", "CODIGO_S", "ID_masa", "water_body_name",
    "subzone", "zone"]

"""
    load_site_level_crosswalk(path)

Load the site -> (sub-basin, water body) crosswalk generated by
`viz/scripts/export-crosswalk.mjs`.  Missing files return an empty table so
callers can fall back to the model's own subcatchment field.
"""
function load_site_level_crosswalk(path::AbstractString)
    if !isfile(path)
        @warn "site crosswalk not found; water bodies will be UNASSIGNED" path
        return DataFrame(CODIGO=String[], CODIGO_S=String[], ID_masa=String[],
            water_body_name=String[], subzone=String[], zone=String[])
    end
    df = CSV.read(path, DataFrame; types=Dict(name => String for name in CROSSWALK_COLUMNS))
    for col in CROSSWALK_COLUMNS
        hasproperty(df, Symbol(col)) || (df[!, Symbol(col)] = fill("", nrow(df)))
        df[!, Symbol(col)] = [ismissing(v) ? "" : string(v) for v in df[!, Symbol(col)]]
    end
    return df
end

"""
    site_level_vectors(sites, site_df, crosswalk)

Return aligned vectors mapping every modelled site to its sub-basin, water body,
water-body name, zone, subzone and basin id.  The model's own `CODIGO_S` is
preferred for the sub-basin because it is the field used by the simulation
scenarios; the crosswalk resolves the water body.
"""
function site_level_vectors(sites::AbstractVector, site_df, crosswalk)
    lookup = Dict{String,NamedTuple}()
    for row in eachrow(crosswalk)
        code = string(row.CODIGO)
        isempty(code) && continue
        lookup[code] = (
            subcatchment=string(row.CODIGO_S),
            water_body=string(row.ID_masa),
            water_body_name=string(row.water_body_name),
            zone=string(row.zone),
            subzone=string(row.subzone),
        )
    end

    site_to_subcatchment = Dict{String,String}()
    if site_df !== nothing
        for row in eachrow(site_df)
            site_to_subcatchment[string(row.CODIGO)] = string(row.CODIGO_S)
        end
    end

    n = length(sites)
    subcatchment = Vector{String}(undef, n)
    water_body = Vector{String}(undef, n)
    water_body_name = Vector{String}(undef, n)
    zone = Vector{String}(undef, n)
    subzone = Vector{String}(undef, n)
    for (i, site) in enumerate(sites)
        code = string(site)
        info = get(lookup, code, nothing)
        sc = get(site_to_subcatchment, code, "")
        if isempty(sc)
            sc = info === nothing ? "" : info.subcatchment
        end
        subcatchment[i] = isempty(sc) ? GUADEX_UNKNOWN_SUBCATCHMENT : sc
        water_body[i] = info === nothing || isempty(info.water_body) ?
            GUADEX_UNASSIGNED_WATER_BODY : info.water_body
        water_body_name[i] = info === nothing ? "" : info.water_body_name
        zone[i] = info === nothing ? "" : info.zone
        subzone[i] = info === nothing ? "" : info.subzone
    end
    return (
        subcatchment=subcatchment,
        water_body=water_body,
        water_body_name=water_body_name,
        zone=zone,
        subzone=subzone,
        basin=fill(GUADEX_BASIN_ID, n),
    )
end

# =============================================================================
# --- Site metrics ---
# =============================================================================

"""
    species_indices(species, codes)

Indices of `codes` inside `species`, preserving the order of `species`, and
skipping codes that are not modelled.
"""
function species_indices(species::AbstractVector, codes::AbstractVector)
    wanted = Set(lowercase(string(c)) for c in codes)
    return [i for (i, sp) in enumerate(species) if lowercase(string(sp)) in wanted]
end

_snapshot_index(times::AbstractVector, target::Real) =
    clamp(searchsortedfirst(times, target), 1, length(times))

function _report_offsets(times::AbstractVector, days_per_year::Real)
    last_year = floor(Int, times[end] / days_per_year)
    return collect(0:last_year)
end

"""
    site_connectivity_metrics(sites, dams, distance_matrix)

Per-site connectivity summary derived from the effective passability matrix:
number of restricted inbound links, mean inbound/outbound passability, and the
number of barriers that restrict upstream movement.
"""
function site_connectivity_metrics(sites::AbstractVector, dams, distance_matrix)
    n = length(sites)
    barriers_in = zeros(Int, n)
    barriers_out = zeros(Int, n)
    sum_pass_in = zeros(Float64, n)
    sum_pass_out = zeros(Float64, n)
    connected_in = zeros(Int, n)
    connected_out = zeros(Int, n)

    I, J, _ = findnz(sparse(distance_matrix))
    for k in 1:length(I)
        i, j = I[k], J[k]
        (i == j) && continue
        # Dispersal entries are keyed [destination, origin]: movement j -> i
        # uses dams[i, j], so the same entry is the inbound passability of i and
        # the outbound passability of the origin j.
        inbound = dams[i, j]
        connected_in[i] += 1
        sum_pass_in[i] += inbound
        inbound < 1.0 && (barriers_in[i] += 1)

        connected_out[j] += 1
        sum_pass_out[j] += inbound
        inbound < 1.0 && (barriers_out[j] += 1)
    end

    mean_in = [connected_in[i] > 0 ? sum_pass_in[i] / connected_in[i] : 1.0 for i in 1:n]
    mean_out = [connected_out[i] > 0 ? sum_pass_out[i] / connected_out[i] : 1.0 for i in 1:n]
    return (
        barrier_links_in=barriers_in,
        barrier_links_out=barriers_out,
        mean_in_passability=mean_in,
        mean_out_passability=mean_out,
    )
end

"""
    compute_site_metrics(times, states, levels, species, native_codes, invasive_codes;
                         days_per_year, threshold, year_offsets, year_labels,
                         temperature_baseline, warming, ...)

Build the long per-site metric table.  `states` is the ODE solution vector
(`sol.u`), each entry a flattened `n_sites * n_species` state.

Optional site attributes (`connectivity`, `habitat_suitability`) are repeated on
every year row.  `warming` is an `n_sites × n_year` matrix of temperature
anomalies (relative to `temperature_baseline`).
"""
function compute_site_metrics(times::AbstractVector, states::AbstractVector,
        sites::AbstractVector, levels, species::AbstractVector,
        native_codes::AbstractVector, invasive_codes::AbstractVector;
        days_per_year::Real=365.0, threshold::Real=0.1,
        year_offsets::Union{Nothing,AbstractVector}=nothing,
        year_labels::Union{Nothing,AbstractVector}=nothing,
        temperature_baseline::Union{Nothing,AbstractVector}=nothing,
        warming::Union{Nothing,AbstractMatrix}=nothing,
        connectivity::Union{Nothing,NamedTuple}=nothing,
        habitat_suitability::Union{Nothing,AbstractVector}=nothing,
        upstream_cost::Union{Nothing,Real}=nothing)

    n_sites = length(levels.subcatchment)
    n_species = length(species)
    offsets = year_offsets === nothing ? _report_offsets(times, days_per_year) : collect(year_offsets)
    labels = year_labels === nothing ? offsets : collect(year_labels)
    length(labels) == length(offsets) ||
        error("year_labels and year_offsets must have the same length")

    native_idx = species_indices(species, native_codes)
    invasive_idx = species_indices(species, invasive_codes)

    baseline_idx = _snapshot_index(times, 0.0)
    baseline = reshape(states[baseline_idx], n_sites, n_species)
    native_baseline = [sum(baseline[i, native_idx] .> threshold) for i in 1:n_sites]

    barrier_count = connectivity === nothing ? nothing : connectivity.barrier_links_in
    barrier_count_out = connectivity === nothing ? nothing : connectivity.barrier_links_out
    pass_in = connectivity === nothing ? nothing : connectivity.mean_in_passability
    pass_out = connectivity === nothing ? nothing : connectivity.mean_out_passability

    rows = Vector{NamedTuple}(undef, n_sites * length(offsets))
    row = 0
    for (k, offset) in enumerate(offsets)
        idx = _snapshot_index(times, Float64(offset) * days_per_year)
        mat = reshape(states[idx], n_sites, n_species)
        label = labels[k]
        for i in 1:n_sites
            row += 1
            native_rich = sum(mat[i, native_idx] .> threshold)
            invasive_rich = sum(mat[i, invasive_idx] .> threshold)
            total_rich = sum(mat[i, :] .> threshold)
            native_biomass = sum(mat[i, native_idx])
            invasive_biomass = sum(mat[i, invasive_idx])
            total_biomass = sum(mat[i, :])
            relative = native_baseline[i] > 0 ? native_rich / native_baseline[i] : 1.0
            temp_base = temperature_baseline === nothing ? NaN : Float64(temperature_baseline[i])
            # A missing warming matrix means "no projected change", so report
            # the baseline temperature rather than an unknown (NaN) value.
            delta = warming === nothing ? 0.0 : Float64(warming[i, k])
            temp_proj = isnan(temp_base) || isnan(delta) ? NaN : temp_base + delta
            habitat = habitat_suitability === nothing ? NaN : Float64(habitat_suitability[i])
            rows[row] = (
                year=label,
                CODIGO=string(sites[i]),
                subcatchment=levels.subcatchment[i],
                water_body=levels.water_body[i],
                water_body_name=levels.water_body_name[i],
                zone=levels.zone[i],
                subzone=levels.subzone[i],
                basin=levels.basin[i],
                native_richness=native_rich,
                invasive_richness=invasive_rich,
                total_richness=total_rich,
                native_biomass=native_biomass,
                invasive_biomass=invasive_biomass,
                total_biomass=total_biomass,
                native_richness_relative=relative,
                native_extinction_risk=max(0.0, 1.0 - relative),
                temperature_c=temp_proj,
                delta_temperature_c=delta,
                habitat_suitability=habitat,
                habitat_degradation_index=isnan(habitat) ? NaN : 1.0 - habitat,
                barrier_links_in=barrier_count === nothing ? -1 : barrier_count[i],
                barrier_links_out=barrier_count_out === nothing ? -1 : barrier_count_out[i],
                mean_in_passability=pass_in === nothing ? NaN : pass_in[i],
                mean_out_passability=pass_out === nothing ? NaN : pass_out[i],
                upstream_cost=upstream_cost === nothing ? NaN : Float64(upstream_cost),
            )
        end
    end
    return DataFrame(rows)
end

# =============================================================================
# --- Four-level aggregation ---
# =============================================================================

const LEVEL_METRIC_EXCLUDE = Set(["year", "CODIGO", "subcatchment", "water_body",
    "water_body_name", "zone", "subzone", "basin"])

_metric_columns(df::DataFrame) =
    [name for name in names(df) if !(name in LEVEL_METRIC_EXCLUDE)]

"""
    aggregate_metrics(site_metrics; level)

Aggregate the per-site metric table to one level by averaging every numeric
metric across the sites of each group and year.  `level` is one of
`:subcatchment`, `:water_body`, `:basin`.
"""
function aggregate_metrics(site_metrics::DataFrame; level::Symbol)
    if nrow(site_metrics) == 0
        return DataFrame()
    end
    group_col = level
    group_col in propertynames(site_metrics) ||
        error("site_metrics has no column $group_col; rebuild with site_level_vectors")

    metric_cols = _metric_columns(site_metrics)
    grouped = groupby(site_metrics, [group_col, :year])
    group_ids = String[]
    years = Int[]
    n_sites = Int[]
    water_body_names = String[]
    means = Dict{Symbol,Vector{Float64}}(
        Symbol("mean_", col) => Float64[] for col in metric_cols)

    for sub in grouped
        first_row = sub[1, :]
        push!(group_ids, string(first_row[group_col]))
        push!(years, Int(first_row.year))
        push!(n_sites, nrow(sub))
        push!(water_body_names, group_col == :water_body ? string(first_row.water_body_name) : "")
        for col in metric_cols
            values = Float64[]
            for v in sub[!, col]
                (v isa Real && isfinite(float(v))) && push!(values, Float64(v))
            end
            push!(means[Symbol("mean_", col)], isempty(values) ? NaN : mean(values))
        end
    end

    order = sortperm(group_ids)
    out = DataFrame()
    out[!, group_col] = group_ids[order]
    group_col == :water_body && (out[!, :water_body_name] = water_body_names[order])
    out[!, :year] = years[order]
    out[!, :n_sites] = n_sites[order]
    for col in metric_cols
        out[!, Symbol("mean_", col)] = means[Symbol("mean_", col)][order]
    end
    return out
end

# =============================================================================
# --- Writers ---
# =============================================================================

function _write_level_tables(site_metrics::DataFrame, output_dir::AbstractString)
    mkpath(output_dir)
    tables = Dict{Symbol,DataFrame}()

    CSV.write(joinpath(output_dir, "level_sampling_point.csv"), site_metrics)

    for (level, filename) in ((:subcatchment, "level_subcatchment.csv"),
                              (:water_body, "level_water_body.csv"),
                              (:basin, "level_basin.csv"))
        agg = aggregate_metrics(site_metrics; level=level)
        if nrow(agg) > 0
            sort!(agg, [level, :year])
        end
        tables[level] = agg
        CSV.write(joinpath(output_dir, filename), agg)
    end
    return tables
end

const VIEWER_METRICS = [
    (:native_richness, "Native richness", "species"),
    (:invasive_richness, "Invasive richness", "species"),
    (:total_richness, "Total richness", "species"),
    (:native_extinction_risk, "Native richness loss (relative)", "fraction"),
    (:native_biomass, "Native biomass", "density"),
    (:invasive_biomass, "Invasive biomass", "density"),
    (:total_biomass, "Total biomass", "density"),
    (:temperature_c, "Projected water temperature", "degC"),
    (:delta_temperature_c, "Water-temperature change", "degC"),
]

function _viewer_steps(site_metrics::DataFrame)
    return sort(unique(site_metrics.year))
end

function _viewer_site_series(site_metrics::DataFrame, metric::Symbol)
    steps = _viewer_steps(site_metrics)
    data = Dict{String,Vector{Any}}()
    lookup = Dict{Tuple{String,Int},Float64}()
    for row in eachrow(site_metrics)
        v = row[metric]
        lookup[(string(row.CODIGO), Int(row.year))] =
            (v isa Real && isfinite(float(v))) ? Float64(v) : NaN
    end
    for code in sort(unique(string.(site_metrics.CODIGO)))
        data[code] = Any[get(lookup, (code, Int(s)), NaN) for s in steps]
    end
    return steps, data
end

"""
    write_viewer_outputs(site_metrics, level_tables, output_dir; name, primary_metric)

Write the files the `viz/` explorer consumes.  A wide long CSV carries every
viewer metric at every step (the explorer detects the `step` column and creates
one variable per metric), plus one canonical time-series JSON for the primary
metric and a per-site metric JSON at the final step.

When `level_tables` is supplied, per-level time-series JSON files are also
written so non-site levels can be attached to the viewer later.
"""
function write_viewer_outputs(site_metrics::DataFrame, output_dir::AbstractString;
        name::AbstractString="GuadeX simulation output",
        primary_metric::Symbol=:native_extinction_risk,
        level_tables::Union{Nothing,Dict}=nothing)
    mkpath(output_dir)
    steps = _viewer_steps(site_metrics)
    codes = sort(unique(string.(site_metrics.CODIGO)))
    lookup = Dict{Tuple{String,Int,Symbol},Any}()
    metric_syms = [m[1] for m in VIEWER_METRICS]
    for row in eachrow(site_metrics)
        for metric in metric_syms
            v = row[metric]
            lookup[(string(row.CODIGO), Int(row.year), metric)] =
                (v isa Real && isfinite(float(v))) ? Float64(v) : NaN
        end
    end

    # 1. Wide long CSV keyed by CODIGO + step (viewer CSV parser).  Use CSV.jl
    # so site codes and values are quoted/escaped correctly.
    csv_path = joinpath(output_dir, "guadex_results_timeseries.csv")
    n_rows = length(codes) * length(steps)
    codigo_col = Vector{String}(undef, n_rows)
    step_col = Vector{Int}(undef, n_rows)
    value_cols = Dict(metric => Vector{Float64}(undef, n_rows) for metric in metric_syms)
    r = 0
    for code in codes, step in steps
        r += 1
        codigo_col[r] = code
        step_col[r] = Int(step)
        for metric in metric_syms
            value_cols[metric][r] = lookup[(code, Int(step), metric)]
        end
    end
    csv_table = DataFrame(CODIGO=codigo_col, step=step_col)
    for metric in metric_syms
        csv_table[!, metric] = value_cols[metric]
    end
    CSV.write(csv_path, csv_table)

    # 2. Canonical single-variable time series (matches viz/README example).
    primary_label = name
    primary_unit = ""
    for (metric, label, unit) in VIEWER_METRICS
        if metric == primary_metric
            primary_label = label
            primary_unit = unit
        end
    end
    _, primary_data = _viewer_site_series(site_metrics, primary_metric)
    _write_timeseries_json(joinpath(output_dir, "guadex_results_$(primary_metric).json"),
        primary_label, primary_unit, string.(steps), primary_data)

    # 3. Per-site metrics at the final step.
    final_step = steps[end]
    metrics_data = Dict{String,Dict{String,Any}}()
    for code in codes
        rec = Dict{String,Any}()
        for (metric, _, _) in VIEWER_METRICS
            v = lookup[(code, Int(final_step), metric)]
            isfinite(v) && (rec[string(metric)] = v)
        end
        metrics_data[code] = rec
    end
    _write_metrics_json(joinpath(output_dir, "guadex_results_metrics.json"), name, metrics_data)

    # 4. Optional per-level time series (JSON keyed by group id), one file per
    # metric so the viewer can attach non-site levels later.
    if level_tables !== nothing
        for (level, table) in level_tables
            nrow(table) == 0 && continue
            id_col = level == :basin ? nothing : level
            for metric in (:native_richness, :native_extinction_risk, :total_biomass)
                col = Symbol("mean_", metric)
                hasproperty(table, col) || continue
                data = Dict{String,Vector{Any}}()
                for row in eachrow(table)
                    key = id_col === nothing ? GUADEX_BASIN_ID : string(row[id_col])
                    series = get!(data, key, fill(NaN, length(steps)))
                    pos = findfirst(==(row.year), steps)
                    pos === nothing && continue
                    v = row[col]
                    series[pos] = (v isa Real && isfinite(float(v))) ? Float64(v) : NaN
                end
                _write_timeseries_json(joinpath(output_dir,
                        "level_$(level)_$(col)_timeseries.json"),
                    "GuadeX $(level) $(metric)", "", string.(steps), data)
            end
        end
    end

    return output_dir
end

# =============================================================================
# --- Temperature projections (guadex_tw) ---
# =============================================================================

"""
    load_temperature_projections(path)

Load `guadex_tw/outputs/tables/water_temp_future_2045.csv`.
"""
function load_temperature_projections(path::AbstractString)
    isfile(path) || error("temperature projection table not found: $path")
    return CSV.read(path, DataFrame)
end

function _linear_interp(xs::AbstractVector, ys::AbstractVector, x::Real)
    length(xs) == length(ys) || error("interpolation vectors differ in length")
    length(xs) == 0 && return NaN
    x <= xs[1] && return Float64(ys[1])
    x >= xs[end] && return Float64(ys[end])
    for i in 1:(length(xs) - 1)
        if xs[i] <= x <= xs[i + 1]
            w = (x - xs[i]) / (xs[i + 1] - xs[i])
            return Float64(ys[i]) * (1 - w) + Float64(ys[i + 1]) * w
        end
    end
    return Float64(ys[end])
end

"""
    basin_warming_curve(projections, scenario, gcm, start_year, end_year)

Basin-level year-by-year warming curve for one GCM/scenario.  The across-site
median `delta_tw_mean` of each projection window is anchored at its window
midpoint, with zero warming at `start_year`, and linearly interpolated.
"""
function basin_warming_curve(projections::DataFrame, scenario::AbstractString,
        gcm::AbstractString, start_year::Int, end_year::Int)
    sub = filter(row -> string(row.scenario) == scenario && string(row.gcm) == gcm, projections)
    nrow(sub) == 0 && error("no projections for scenario=$scenario gcm=$gcm")

    windows = sort(unique([(Int(r.period_start), Int(r.period_end)) for r in eachrow(sub)]))
    nodes = Float64[start_year]
    values = Float64[0.0]
    for (p0, p1) in windows
        midpoint = (p0 + p1) / 2.0
        midpoint <= start_year && continue
        deltas = [Float64(r.delta_tw_mean) for r in eachrow(sub)
                  if Int(r.period_start) == p0 && Int(r.period_end) == p1]
        isempty(deltas) && continue
        push!(nodes, midpoint)
        push!(values, median(deltas))
    end
    order = sortperm(nodes)
    nodes = nodes[order]
    values = values[order]

    years = collect(start_year:end_year)
    curve = [_linear_interp(nodes, values, Float64(y)) for y in years]
    return years, curve
end

"""
    warming_matrix(n_sites, years, curve; elevations, b1, b3, z_ref, elevation_scaling)

Expand a basin warming curve to an `n_sites × length(years)` matrix.  When
`elevation_scaling` is true the anomaly is scaled by the calibrated
water-to-air response `(b1 + b3·z)`, matching the elevation dependence reported
by the guadex_tw pipeline; otherwise every site receives the basin curve.
"""
function warming_matrix(n_sites::Int, years::AbstractVector, curve::AbstractVector;
        elevations::Union{Nothing,AbstractVector}=nothing,
        b1::Real=0.51, b3::Real=-0.163, z_ref::Real=0.0,
        elevation_scaling::Bool=false)
    length(years) == length(curve) || error("years and curve must align")
    matrix = Matrix{Float64}(undef, n_sites, length(years))
    for i in 1:n_sites
        factor = 1.0
        if elevation_scaling
            elevations === nothing && error("elevation_scaling requires elevations")
            z = Float64(elevations[i]) / 1000.0
            reference = b1 + b3 * z_ref
            factor = reference == 0 ? 1.0 : (b1 + b3 * z) / reference
        end
        for (k, value) in enumerate(curve)
            matrix[i, k] = factor * Float64(value)
        end
    end
    return matrix
end

# =============================================================================
# --- Orchestration ---
# =============================================================================

"""
    export_run_outputs(output_dir; sol_t, sol_u, sites, species, site_df, ...)

Write the four-level CSV tables, viewer files and metadata for one run.
"""
function export_run_outputs(output_dir::AbstractString;
        sol_t::AbstractVector, sol_u::AbstractVector,
        sites::AbstractVector, species::AbstractVector,
        site_df=nothing,
        crosswalk_path::AbstractString=joinpath("data", "site_waterbody_crosswalk.csv"),
        native_species::AbstractVector=String[],
        invasive_species::AbstractVector=String[],
        days_per_year::Real=365.0,
        threshold::Real=0.1,
        report_year_offsets::Union{Nothing,AbstractVector}=nothing,
        report_year_labels::Union{Nothing,AbstractVector}=nothing,
        temperature_baseline::Union{Nothing,AbstractVector}=nothing,
        warming::Union{Nothing,AbstractMatrix}=nothing,
        dams=nothing,
        distance_matrix=nothing,
        habitat_suitability::Union{Nothing,AbstractVector}=nothing,
        upstream_cost::Union{Nothing,Real}=nothing,
        run_metadata::AbstractDict=Dict{String,Any}(),
        primary_metric::Symbol=:native_extinction_risk,
        viewer_name::AbstractString="GuadeX simulation output",
        require_crosswalk::Bool=false)

    crosswalk = load_site_level_crosswalk(crosswalk_path)
    levels = site_level_vectors(sites, site_df, crosswalk)
    assigned_water_bodies = count(!=(GUADEX_UNASSIGNED_WATER_BODY), levels.water_body)
    if require_crosswalk && (isempty(levels.water_body) || assigned_water_bodies == 0)
        error("no sites resolved to a water body; refusing to write an " *
              "all-UNASSIGNED water-body level (crosswalk: $crosswalk_path)")
    end

    connectivity = nothing
    if dams !== nothing && distance_matrix !== nothing
        connectivity = site_connectivity_metrics(sites, dams, distance_matrix)
    end

    site_metrics = compute_site_metrics(sol_t, sol_u, sites, levels, species,
        native_species, invasive_species;
        days_per_year=days_per_year,
        threshold=threshold,
        year_offsets=report_year_offsets,
        year_labels=report_year_labels,
        temperature_baseline=temperature_baseline,
        warming=warming,
        connectivity=connectivity,
        habitat_suitability=habitat_suitability,
        upstream_cost=upstream_cost)

    mkpath(output_dir)
    level_tables = _write_level_tables(site_metrics, joinpath(output_dir, "levels"))

    write_viewer_outputs(site_metrics, joinpath(output_dir, "viewer");
        name=viewer_name, primary_metric=primary_metric, level_tables=level_tables)

    metadata = Dict{String,Any}()
    for (k, v) in run_metadata
        metadata[string(k)] = v
    end
    metadata["n_sites"] = length(sites)
    metadata["n_species"] = length(species)
    metadata["levels"] = ["sampling_point", "subcatchment", "water_body", "basin"]
    metadata["crosswalk_file"] = crosswalk_path
    metadata["crosswalk_water_bodies"] = assigned_water_bodies
    metadata["outputs"] = [
        "levels/level_sampling_point.csv",
        "levels/level_subcatchment.csv",
        "levels/level_water_body.csv",
        "levels/level_basin.csv",
        "viewer/guadex_results_timeseries.csv",
        "viewer/guadex_results_metrics.json",
        "viewer/guadex_results_$(primary_metric).json",
    ]
    _write_value_json(joinpath(output_dir, "run_metadata.json"), metadata)

    return (
        site_metrics=site_metrics,
        level_tables=level_tables,
        levels=levels,
        output_dir=output_dir,
    )
end
