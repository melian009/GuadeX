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
    <run_dir>/viewer/guadex_results_<metric>_timeseries.csv   # one metric per file
    <run_dir>/viewer/guadex_results_metrics.json
    <run_dir>/viewer/guadex_results_<metric>.json
    <run_dir>/viewer/level_<level>_mean_<metric>_timeseries.json          # group-keyed
    <run_dir>/viewer/level_<level>_mean_<metric>_by_site_timeseries.json  # renderable
    <run_dir>/viewer/species_<sp>_<metric>_timeseries.json                # species-level
    <run_dir>/run_metadata.json

The viewer files are keyed exactly as `viz/src/data/ResultsModel.js` expects:
`{name, unit, steps, data}` for time series and `{name, data}` for per-site
metrics, with site keys equal to `CODIGO`.  Each `viewer/*_timeseries.csv` file
holds a single metric (`CODIGO,step,value`) so the viewer's time slider spans the
years of one metric rather than one entry per metric × year (minor #14).
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
        xf = Float64(x)
        # Only collapse to an integer when the conversion is exact and safe:
        # |x| <= 2^53 is exactly representable and far below typemax(Int64), so a
        # ratio that has blown up (e.g. division by a near-zero burn-in baseline)
        # is written as a JSON float instead of throwing InexactError.
        return (isinteger(xf) && abs(xf) <= 9.007199254740992e15) ?
            string(Int(xf)) : string(xf)
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

function _snapshot_index(times::AbstractVector, target::Real; atol::Real=1.0e-6)
    idx = clamp(searchsortedfirst(times, target), 1, length(times))
    # The `saveat` grid can place the intended point an ulp below `target`
    # (e.g. `365/12` accumulated twelve times, or a rounded interval like
    # 30.4167).  `searchsortedfirst` would then skip it and report the next
    # grid point, shifting the snapshot by a step.  Snap to the preceding grid
    # point only when it is, to within round-off, exactly the requested time, so
    # alignment stays exact rather than nearest-neighbour.
    if Float64(times[idx]) != Float64(target) && idx > 1 &&
            abs(Float64(times[idx - 1]) - Float64(target)) <= atol
        idx -= 1
    end
    return idx
end

# End-of-year reporting offsets (E5).  Offset `k` is the snapshot at
# `t = k · days_per_year`, i.e. the state at the end of year `k`, so the default
# grid is 1-based (`1:last_year`) and matches the drivers, which were moved from
# 0-based to 1-based offsets in E5.  The previous 0-based `0:last_year` grid
# reported the pre-forcing initial state as the first row.
function _report_offsets(times::AbstractVector, days_per_year::Real)
    last_year = floor(Int, times[end] / days_per_year)
    return collect(1:last_year)
end

"""
    quasi_extinction_flags(density_series; threshold, baseline_density, q=0.1,
                           persistence=3, established=baseline_density > threshold)

Boolean flags for one site × species series: `true` from the moment the density
has been below `max(presence_threshold, q · baseline_density)` for `persistence`
consecutive annual snapshots.  `baseline_density` is the spun-up (or t = 0)
reference density used to define a collapse relative to the species' own local
abundance, so the 0.1 presence threshold does not hide sub-threshold declines
(WP5).

A site where the species was **not established** in the baseline — baseline
density at or below the presence threshold — is never flagged: a site that was
never occupied cannot be quasi-extinct (E6).  `established` defaults to that
test and can be passed explicitly when the caller has already computed the mask.
"""
function quasi_extinction_flags(density_series::AbstractVector;
        threshold::Real, baseline_density::Real,
        established::Bool=Float64(baseline_density) > Float64(threshold),
        q::Real=0.1, persistence::Int=3)
    established || return falses(length(density_series))
    cutoff = max(Float64(threshold), q * Float64(baseline_density))
    flags = falses(length(density_series))
    run = 0
    for (k, d) in enumerate(density_series)
        run = (isfinite(d) && d < cutoff) ? run + 1 : 0
        flags[k] = run >= persistence
    end
    return flags
end

"""
    time_to_quasi_extinction(flags; year_labels)

First year label at which a quasi-extinction flag turns true, or `missing`.
"""
function time_to_quasi_extinction(flags::AbstractVector{Bool}; year_labels::AbstractVector)
    idx = findfirst(flags)
    return idx === nothing ? missing : Int(year_labels[idx])
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
        upstream_cost::Union{Nothing,Real}=nothing,
        baseline_species_density::Union{Nothing,AbstractMatrix}=nothing,
        quasi_extinction_q::Real=0.1,
        quasi_extinction_persistence::Int=3)

    n_sites = length(levels.subcatchment)
    n_species = length(species)
    offsets = year_offsets === nothing ? _report_offsets(times, days_per_year) : collect(year_offsets)
    labels = year_labels === nothing ? offsets : collect(year_labels)
    length(labels) == length(offsets) ||
        error("year_labels and year_offsets must have the same length")

    native_idx = species_indices(species, native_codes)
    invasive_idx = species_indices(species, invasive_codes)

    # The reference state is the spun-up equilibrium when supplied (WP4),
    # otherwise the t = 0 snapshot.  It defines relative abundance and the
    # quasi-extinction collapse threshold.
    if baseline_species_density === nothing
        baseline_idx = _snapshot_index(times, 0.0)
        baseline = reshape(states[baseline_idx], n_sites, n_species)
    else
        size(baseline_species_density) == (n_sites, n_species) ||
            error("baseline_species_density must be n_sites × n_species")
        baseline = Float64.(baseline_species_density)
    end
    native_baseline = [sum(baseline[i, native_idx] .> threshold) for i in 1:n_sites]
    native_biomass_baseline = [sum(baseline[i, native_idx]) for i in 1:n_sites]
    total_biomass_baseline = [sum(baseline[i, :]) for i in 1:n_sites]

    # Annual snapshot matrices, used for the quasi-extinction run-length flags.
    snapshot_mats = [reshape(states[_snapshot_index(times, Float64(offset) * days_per_year)],
        n_sites, n_species) for offset in offsets]
    quasi = falses(n_sites, n_species, length(offsets))
    for i in 1:n_sites, s in 1:n_species
        series = [snapshot_mats[k][i, s] for k in 1:length(offsets)]
        # E6: only sites where the species was established in the baseline can
        # go quasi-extinct; never-occupied sites are excluded.
        quasi[i, s, :] .= quasi_extinction_flags(series;
            threshold=threshold, baseline_density=baseline[i, s],
            established=baseline[i, s] > threshold,
            q=quasi_extinction_q, persistence=quasi_extinction_persistence)
    end

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
            native_occupancy = isempty(native_idx) ? NaN :
                count(>(threshold), mat[i, native_idx]) / length(native_idx)
            invasive_occupancy = isempty(invasive_idx) ? NaN :
                count(>(threshold), mat[i, invasive_idx]) / length(invasive_idx)
            total_occupancy = n_species == 0 ? NaN : count(>(threshold), mat[i, :]) / n_species
            native_biomass_relative = native_biomass_baseline[i] > 0 ?
                native_biomass / native_biomass_baseline[i] : 1.0
            total_biomass_relative = total_biomass_baseline[i] > 0 ?
                total_biomass / total_biomass_baseline[i] : 1.0
            native_quasi_extinct = isempty(native_idx) ? 0 : count(quasi[i, native_idx, k])
            # E6: divide by the native species established at this site in the
            # baseline, not by every modelled native species.
            native_established = native_baseline[i]
            native_quasi_extinct_fraction = native_established > 0 ?
                native_quasi_extinct / native_established : NaN
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
                native_baseline_richness=native_baseline[i],
                invasive_richness=invasive_rich,
                total_richness=total_rich,
                native_biomass=native_biomass,
                invasive_biomass=invasive_biomass,
                total_biomass=total_biomass,
                native_richness_relative=relative,
                realised_richness_loss=max(0.0, 1.0 - relative),
                native_biomass_relative=native_biomass_relative,
                total_biomass_relative=total_biomass_relative,
                native_occupancy=native_occupancy,
                invasive_occupancy=invasive_occupancy,
                total_occupancy=total_occupancy,
                native_quasi_extinct=native_quasi_extinct,
                native_quasi_extinct_fraction=native_quasi_extinct_fraction,
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

"""
    compute_species_metrics(times, states, sites, species; days_per_year, threshold,
        year_offsets, year_labels, baseline_species_density, q, persistence)

Per site × species annual table (WP5): density, thresholded presence, relative
abundance against the spun-up baseline, whether the species was established in
the baseline (`baseline_established`), the quasi-extinction flag and the
time-to-quasi-extinction for each site/species.  `export_run_outputs` writes this
as `levels/species_timeseries.csv`.
"""
function compute_species_metrics(times::AbstractVector, states::AbstractVector,
        sites::AbstractVector, species::AbstractVector;
        days_per_year::Real=365.0, threshold::Real=0.1,
        year_offsets::Union{Nothing,AbstractVector}=nothing,
        year_labels::Union{Nothing,AbstractVector}=nothing,
        baseline_species_density::Union{Nothing,AbstractMatrix}=nothing,
        quasi_extinction_q::Real=0.1, quasi_extinction_persistence::Int=3)
    n_sites = length(sites)
    n_species = length(species)
    offsets = year_offsets === nothing ? _report_offsets(times, days_per_year) : collect(year_offsets)
    labels = year_labels === nothing ? offsets : collect(year_labels)
    length(labels) == length(offsets) ||
        error("year_labels and year_offsets must have the same length")

    if baseline_species_density === nothing
        baseline = reshape(states[_snapshot_index(times, 0.0)], n_sites, n_species)
    else
        size(baseline_species_density) == (n_sites, n_species) ||
            error("baseline_species_density must be n_sites × n_species")
        baseline = Float64.(baseline_species_density)
    end

    snapshot_mats = [reshape(states[_snapshot_index(times, Float64(offset) * days_per_year)],
        n_sites, n_species) for offset in offsets]

    # E6: baseline-established mask, threaded into the flags and stored in the
    # table so the denominator used by `quasi_extinction_summary` is auditable.
    baseline_established = baseline .> threshold

    qe = falses(n_sites, n_species, length(offsets))
    time_to_qe = Matrix{Union{Missing,Int}}(missing, n_sites, n_species)
    for i in 1:n_sites, s in 1:n_species
        series = [snapshot_mats[k][i, s] for k in 1:length(offsets)]
        flags = quasi_extinction_flags(series; threshold=threshold,
            baseline_density=baseline[i, s],
            established=baseline_established[i, s],
            q=quasi_extinction_q,
            persistence=quasi_extinction_persistence)
        qe[i, s, :] .= flags
        time_to_qe[i, s] = time_to_quasi_extinction(flags; year_labels=labels)
    end

    rows = NamedTuple[]
    for (k, label) in enumerate(labels)
        mat = snapshot_mats[k]
        for (i, site) in enumerate(sites), (s, sp) in enumerate(species)
            density = mat[i, s]
            base = baseline[i, s]
            push!(rows, (
                year=Int(label),
                CODIGO=string(site),
                species=string(sp),
                density=density,
                baseline_density=base,
                baseline_established=baseline_established[i, s],
                relative_density=base > 0 ? density / base : NaN,
                present=density > threshold,
                quasi_extinct=qe[i, s, k],
                time_to_quasi_extinction=time_to_qe[i, s],
            ))
        end
    end
    return DataFrame(rows)
end

"""
    quasi_extinction_summary(species_metrics; threshold=0.1)

One row per species: final-year occupancy and quasi-extinct site fraction,
biomass change against the baseline, and the median time to quasi-extinction
(over sites that reached it).

`quasi_extinct_fraction` is computed over the **baseline-established** sites
only, i.e. the sites where the species' baseline density is above the presence
threshold (E6).  `occupancy_final` and `n_sites` still cover every modelled
site, while `n_sites_established` exposes the denominator.  When a species was
established at no site in the baseline there is no site at which it could go
quasi-extinct, and the fraction is reported as `0.0` (never `NaN` and never
`1.0`); the established mask is taken from the `baseline_established` column
when present and otherwise recomputed as `baseline_density > threshold`.
"""
function quasi_extinction_summary(species_metrics::DataFrame; threshold::Real=0.1)
    nrow(species_metrics) == 0 && return DataFrame()
    years = sort(unique(species_metrics.year))
    final_year = years[end]
    rows = NamedTuple[]
    for (sp, sub) in pairs(groupby(species_metrics, :species))
        final_rows = sub[sub.year .== final_year, :]
        base = sub[sub.year .== years[1], :]
        times = Float64[Float64(t) for t in unique(sub.time_to_quasi_extinction) if !ismissing(t)]
        times = filter(isfinite, times)
        base_biomass = sum(base.baseline_density)
        final_biomass = sum(final_rows.density)
        n_sites = nrow(final_rows)
        established = hasproperty(final_rows, :baseline_established) ?
            Bool.(final_rows.baseline_established) :
            (Float64.(final_rows.baseline_density) .> Float64(threshold))
        n_sites_established = count(established)
        push!(rows, (
            species=string(first(sub.species)),
            n_sites=n_sites,
            n_sites_established=n_sites_established,
            baseline_biomass=base_biomass,
            final_biomass=final_biomass,
            relative_biomass_change=base_biomass > 0 ? final_biomass / base_biomass - 1.0 : NaN,
            occupancy_final=n_sites == 0 ? NaN :
                count(final_rows.present) / n_sites,
            quasi_extinct_fraction=n_sites_established == 0 ? 0.0 :
                count(final_rows.quasi_extinct .& established) / n_sites_established,
            median_time_to_quasi_extinction=isempty(times) ? missing : median(times),
        ))
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
    (:realised_richness_loss, "Realised richness loss (relative)", "fraction"),
    (:native_biomass, "Native biomass", "density"),
    (:invasive_biomass, "Invasive biomass", "density"),
    (:total_biomass, "Total biomass", "density"),
    (:native_biomass_relative, "Native biomass (relative to baseline)", "fraction"),
    (:total_biomass_relative, "Total biomass (relative to baseline)", "fraction"),
    (:native_occupancy, "Native occupancy", "fraction"),
    (:invasive_occupancy, "Invasive occupancy", "fraction"),
    (:total_occupancy, "Total occupancy", "fraction"),
    (:native_quasi_extinct_fraction, "Native quasi-extinct fraction", "fraction"),
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

const VIEWER_LEVEL_METRICS = (:native_richness, :realised_richness_loss, :total_biomass)
const VIEWER_SPECIES_METRICS = (:density, :relative_density, :present, :quasi_extinct)

# A4: human-readable words for the `name` field the viewer shows in its legend.
# Machine filenames are unchanged; only the embedded description is spelled out.
const VIEWER_LEVEL_LABELS = Dict(
    :subcatchment => "sub-catchment mean",
    :water_body => "water-body mean",
    :basin => "basin mean",
)
const VIEWER_SPECIES_LABELS = Dict(
    :density => "density",
    :relative_density => "relative density",
    :present => "presence",
    :quasi_extinct => "quasi-extinct flag",
)

_species_metric_value(v) = v isa Bool ? (v ? 1.0 : 0.0) :
    (v isa Real && isfinite(float(v)) ? Float64(v) : NaN)

"""
    write_viewer_outputs(site_metrics, output_dir; name, primary_metric,
        level_tables, species_metrics, site_levels)

Write the files the `viz/` explorer consumes.

* one `guadex_results_<metric>_timeseries.csv` (`CODIGO,step,value`) per metric,
  so each viewer slider spans the years of a single metric instead of mixing
  every metric × year in one control (minor #14);
* the canonical single-variable JSON for `primary_metric` and the per-site
  metric JSON at the final step;
* per-level time-series JSON keyed by group id (`level_<level>_mean_<metric>_...`);
* `level_<level>_mean_<metric>_by_site_timeseries.json`, the same aggregate
  repeated onto every member site, so the existing site layer can *render* the
  sub-catchment / water-body / basin aggregate (minor #14);
* `species_<sp>_<metric>_timeseries.json`, one site-keyed time series per metric
  per species, built from the WP5 `species_metrics` table (minor #14).

All time-series JSON files use the `{name, unit, steps, data}` shape expected by
`viz/src/data/ResultsModel.js`, with `data` keyed by `CODIGO`.
"""
function write_viewer_outputs(site_metrics::DataFrame, output_dir::AbstractString;
        name::AbstractString="GuadeX simulation output",
        primary_metric::Symbol=:realised_richness_loss,
        level_tables::Union{Nothing,Dict}=nothing,
        species_metrics::Union{Nothing,DataFrame}=nothing,
        site_levels=nothing)
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

    # 1. One single-variable CSV per metric, keyed by CODIGO + step.  A single
    # wide CSV with 15 metric columns and a `step` column is read by the viewer
    # as one key per (metric, year) pair — 15 metrics × 20 years in one slider.
    # Splitting it keeps each file a proper single-series time series.
    for metric in metric_syms
        n_rows = length(codes) * length(steps)
        codigo_col = Vector{String}(undef, n_rows)
        step_col = Vector{Int}(undef, n_rows)
        value_col = Vector{Float64}(undef, n_rows)
        r = 0
        for code in codes, step in steps
            r += 1
            codigo_col[r] = code
            step_col[r] = Int(step)
            value_col[r] = lookup[(code, Int(step), metric)]
        end
        CSV.write(joinpath(output_dir, "guadex_results_$(metric)_timeseries.csv"),
            DataFrame(CODIGO=codigo_col, step=step_col, value=value_col))
    end

    # 2. Canonical single-variable time series (matches viz/README example).
    # A4: the embedded `name` is the legend text, so spell out run + metric.
    metric_labels = Dict(m => l for (m, l, _) in VIEWER_METRICS)
    primary_label = name
    primary_unit = ""
    for (metric, label, unit) in VIEWER_METRICS
        if metric == primary_metric
            primary_label = "$(name) · $(label)"
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

    # 4. Per-level aggregates: keep the group-keyed JSON (canonical) and add a
    # site-expanded copy so the site layer can actually render them.
    if level_tables !== nothing
        for (level, table) in level_tables
            nrow(table) == 0 && continue
            id_col = level == :basin ? nothing : level
            group_of_site = if level == :basin
                fill(GUADEX_BASIN_ID, length(codes))
            elseif site_levels === nothing || !hasproperty(site_levels, level)
                nothing
            else
                string.(getproperty(site_levels, level))
            end
            for metric in VIEWER_LEVEL_METRICS
                col = Symbol("mean_", metric)
                hasproperty(table, col) || continue
                group_lookup = Dict{Tuple{String,Int},Float64}()
                for row in eachrow(table)
                    key = id_col === nothing ? GUADEX_BASIN_ID : string(row[id_col])
                    v = row[col]
                    group_lookup[(key, Int(row.year))] =
                        (v isa Real && isfinite(float(v))) ? Float64(v) : NaN
                end
                data = Dict{String,Vector{Any}}()
                for (key, year) in keys(group_lookup)
                    data[key] = Any[get(group_lookup, (key, Int(s)), NaN) for s in steps]
                end
                _write_timeseries_json(joinpath(output_dir,
                        "level_$(level)_mean_$(metric)_timeseries.json"),
                    "$(name) · $(get(VIEWER_LEVEL_LABELS, level, string(level))) " *
                    "$(get(metric_labels, metric, string(metric)))",
                    "", string.(steps), data)

                group_of_site === nothing && continue
                site_data = Dict{String,Vector{Any}}()
                for (i, code) in enumerate(codes)
                    g = group_of_site[i]
                    site_data[code] = Any[get(group_lookup, (g, Int(s)), NaN) for s in steps]
                end
                _write_timeseries_json(joinpath(output_dir,
                        "level_$(level)_mean_$(metric)_by_site_timeseries.json"),
                    "$(name) · $(get(VIEWER_LEVEL_LABELS, level, string(level))) " *
                    "$(get(metric_labels, metric, string(metric))) " *
                    "(same aggregate repeated for every member site)",
                    "", string.(steps), site_data)
            end
        end
    end

    # 5. Species-level projections: one site-keyed time series per metric per
    # species (minor #14), from the WP5 species_metrics table.
    if species_metrics !== nothing && nrow(species_metrics) > 0
        series_steps = sort(unique(Int.(species_metrics.year)))
        species_order = sort(unique(string.(species_metrics.species)))
        for metric in VIEWER_SPECIES_METRICS
            hasproperty(species_metrics, metric) || continue
            sp_lookup = Dict{Tuple{String,String,Int},Float64}()
            for row in eachrow(species_metrics)
                sp_lookup[(string(row.CODIGO), string(row.species), Int(row.year))] =
                    _species_metric_value(row[metric])
            end
            for sp in species_order
                data = Dict{String,Vector{Any}}()
                for code in codes
                    data[code] = Any[get(sp_lookup, (code, sp, Int(s)), NaN)
                                     for s in series_steps]
                end
                _write_timeseries_json(joinpath(output_dir,
                        "species_$(sp)_$(metric)_timeseries.json"),
                    "$(name) · $(sp) $(get(VIEWER_SPECIES_LABELS, metric, string(metric)))",
                    "", string.(series_steps), data)
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
        primary_metric::Symbol=:realised_richness_loss,
        viewer_name::AbstractString="GuadeX simulation output",
        require_crosswalk::Bool=false,
        migratory_species::AbstractVector=String[],
        baseline_species_density::Union{Nothing,AbstractMatrix}=nothing,
        quasi_extinction_q::Real=0.1,
        quasi_extinction_persistence::Int=3,
        species_upper_limits::Union{Nothing,AbstractVector}=nothing,
        daily_forcing::Union{Nothing,NamedTuple}=nothing)

    crosswalk = load_site_level_crosswalk(crosswalk_path)
    levels = site_level_vectors(sites, site_df, crosswalk)
    # `levels.water_body` has one entry per site, so this counts SITES with an
    # assigned water body, not distinct water bodies.
    sites_with_assigned_water_body =
        count(!=(GUADEX_UNASSIGNED_WATER_BODY), levels.water_body)
    water_body_ids = Set(levels.water_body)
    n_assigned_water_bodies = count(!=(GUADEX_UNASSIGNED_WATER_BODY), water_body_ids)
    # True number of water-body levels: assigned ids plus the UNASSIGNED category.
    n_water_bodies = n_assigned_water_bodies +
        (GUADEX_UNASSIGNED_WATER_BODY in water_body_ids ? 1 : 0)
    if require_crosswalk && (isempty(levels.water_body) || sites_with_assigned_water_body == 0)
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
        upstream_cost=upstream_cost,
        baseline_species_density=baseline_species_density,
        quasi_extinction_q=quasi_extinction_q,
        quasi_extinction_persistence=quasi_extinction_persistence)

    mkpath(output_dir)
    level_tables = _write_level_tables(site_metrics, joinpath(output_dir, "levels"))

    # WP5: per-species abundance / occupancy / quasi-extinction table.
    species_metrics = compute_species_metrics(sol_t, sol_u, sites, species;
        days_per_year=days_per_year,
        threshold=threshold,
        year_offsets=report_year_offsets,
        year_labels=report_year_labels,
        baseline_species_density=baseline_species_density,
        quasi_extinction_q=quasi_extinction_q,
        quasi_extinction_persistence=quasi_extinction_persistence)
    CSV.write(joinpath(output_dir, "levels", "species_timeseries.csv"), species_metrics)
    quasi_summary = quasi_extinction_summary(species_metrics)
    CSV.write(joinpath(output_dir, "levels", "quasi_extinction_summary.csv"), quasi_summary)

    # WP5: exposure diagnostics (days above each species' empirical upper limit).
    exposure = nothing
    exposure_summary = nothing
    if daily_forcing !== nothing && species_upper_limits !== nothing
        exposure = exposure_table(sites, species, daily_forcing.temps, daily_forcing.dates;
            upper_limits=species_upper_limits)
        CSV.write(joinpath(output_dir, "levels", "exposure_sites.csv"), exposure)
        # C7: restrict the reported exposure to the sites where each species was
        # baseline established (E6), and report mean/median/max there rather than
        # the basin maximum.  The `site_set` column names the sites used.
        exposure_summary = established_exposure_summary(exposure, species_metrics)
        CSV.write(joinpath(output_dir, "levels", "exposure_summary.csv"), exposure_summary)
    end

    # A4: when no explicit viewer name was supplied, derive a descriptive one
    # from the run metadata (scenario / stage / case / GCM) so the `name` field
    # inside every viewer JSON identifies the experiment, not just "simulation".
    effective_viewer_name = viewer_name
    if viewer_name == "GuadeX simulation output" && !isempty(run_metadata)
        detail = String[]
        for key in ("scenario", "stage", "case", "exploitation_scenario",
                    "passability_scenario", "gcm")
            haskey(run_metadata, key) || continue
            value = run_metadata[key]
            value === nothing && continue
            text = string(value)
            (isempty(text) || text in detail) && continue
            push!(detail, text)
        end
        isempty(detail) || (effective_viewer_name = "GuadeX · " * join(detail, " · "))
    end

    write_viewer_outputs(site_metrics, joinpath(output_dir, "viewer");
        name=effective_viewer_name, primary_metric=primary_metric, level_tables=level_tables,
        species_metrics=species_metrics, site_levels=levels)

    metadata = Dict{String,Any}()
    for (k, v) in run_metadata
        metadata[string(k)] = v
    end
    metadata["n_sites"] = length(sites)
    metadata["n_species"] = length(species)
    metadata["levels"] = ["sampling_point", "subcatchment", "water_body", "basin"]
    metadata["crosswalk_file"] = crosswalk_path
    metadata["sites_with_assigned_water_body"] = sites_with_assigned_water_body
    metadata["n_assigned_water_bodies"] = n_assigned_water_bodies
    metadata["n_water_bodies"] = n_water_bodies
    metadata["migratory_species"] = String.(migratory_species)
    metadata["quasi_extinction_q"] = Float64(quasi_extinction_q)
    metadata["quasi_extinction_persistence"] = quasi_extinction_persistence
    # E6: only baseline-established sites enter the quasi-extinction fraction.
    metadata["quasi_extinction_denominator"] = "baseline_established_sites"
    metadata["quasi_extinction_zero_established_fraction"] = 0.0
    # C7: the reported exposure is over the baseline-established sites of each
    # species, not the basin maximum; `exposure_summary.csv` names the site set.
    if exposure_summary !== nothing
        metadata["exposure_site_set"] = isempty(exposure_summary) ? "all_sites" :
            string(first(exposure_summary.site_set))
        metadata["exposure_summary_statistic"] =
            "mean/median/max days above the species thermal limit over the exposure site set"
    end
    metadata["baseline_reference"] = baseline_species_density === nothing ? "t0" : "supplied"
    outputs = [
        "levels/level_sampling_point.csv",
        "levels/level_subcatchment.csv",
        "levels/level_water_body.csv",
        "levels/level_basin.csv",
        "levels/species_timeseries.csv",
        "levels/quasi_extinction_summary.csv",
        "viewer/guadex_results_<metric>_timeseries.csv",
        "viewer/guadex_results_metrics.json",
        "viewer/guadex_results_$(primary_metric).json",
        "viewer/level_<level>_mean_<metric>_by_site_timeseries.json",
        "viewer/species_<sp>_<metric>_timeseries.json",
    ]
    exposure === nothing || push!(outputs, "levels/exposure_sites.csv")
    exposure_summary === nothing || push!(outputs, "levels/exposure_summary.csv")
    metadata["outputs"] = outputs
    _write_value_json(joinpath(output_dir, "run_metadata.json"), metadata)

    return (
        site_metrics=site_metrics,
        species_metrics=species_metrics,
        quasi_extinction_summary=quasi_summary,
        exposure=exposure,
        exposure_summary=exposure_summary,
        level_tables=level_tables,
        levels=levels,
        output_dir=output_dir,
    )
end
