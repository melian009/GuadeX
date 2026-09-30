"""
    Data Preparation Module for ODE Model

This module provides functions to prepare data from various CSV files
for use in the fish metacommunity ODE model.

All functions are designed to handle large files efficiently by reading
them in chunks or using streaming approaches.
"""

# Species codes mapping (from density matrix columns)
const SPECIES_CODES = [
  "SA", "LS", "ST", "SP", "IL", "PW", "CP", "AA", "AH", "LR", "MC", "AB", "IO",  # Native
  "OM", "LG", "GH", "AA", "CG", "CC", "MS", "AM", "TT", "EL", "GL"  # Exotic + others
]

# Full species names mapping
const SPECIES_NAMES = Dict(
  "AB" => "Aphanius baeticus",
  "SP" => "Squalius pyrenaicus",
  "CP" => "Cobitis paludica",
  "IO" => "Iberochondrostoma oretanum",
  "MC" => "Mugil cephalus",
  "CC" => "Cyprinus carpio",
  "GH" => "Gambusia holbrooki",
  "AM" => "Ameiurus melas",
  "SA" => "Squalius alburnoides",
  "PW" => "Pseudochondrostoma willkommii",
  "LG" => "Lepomis gibbosus",
  "LS" => "Luciobarbus sclateri",
  "OM" => "Oncorhynchus mykiss",
  "CG" => "Carassius gibelio",
  "AH" => "Anaecypris hispanica",
  "LR" => "Liza ramada",
  "ST" => "Salmo trutta",
  "TT" => "Tinca tinca",
  "AA" => "Anguilla anguilla",
  "GL" => "Gobio lozanoi",
  "MS" => "Micropterus salmoides",
  "IL" => "Iberochondrostoma lemmingii",
  "EL" => "Esox lucius",
)

"""
    parse_temperature_range(temp_str::AbstractString)

Parse temperature range string like "8 to 30" and return the midpoint as thermal optimum.

Input file is "data/ABIOTIC/caracteristicas_peces_Guadalquivir_03-04-2018.csv"
"""
function parse_temperature_range(temp_str::AbstractString)
    if isempty(temp_str) || temp_str == ""
        return 15.0  # Default fallback
    end

    # Handle "X to Y" format
    if occursin(" to ", temp_str)
        parts = split(temp_str, " to ")
        if length(parts) == 2
            try
                t_min = parse(Float64, strip(parts[1]))
                t_max = parse(Float64, strip(parts[2]))
                return (t_min + t_max) / 2.0
            catch
                return 15.0
            end
        end
    end

    # Try to parse as single number
    try
        return parse(Float64, temp_str)
    catch
        return 15.0
    end
end

"""
    parse_temperature_range_and_sigma(temp_str::AbstractString)

Parse temperature range string like "8 to 30" and return the thermal optimum
(midpoint), thermal breadth (sigma) and the empirical lower/upper limits.  Sigma
is derived from the temperature range using the approximation sigma ≈ range / 6.

The returned `lower`/`upper` bounds are the *empirical* limits from the trait
table (WP3) and are used for the heat-stress term; when the input is a single
value or unparseable, they widen to ±Inf so no stress is applied.

Input file is "data/ABIOTIC/caracteristicas_peces_Guadalquivir_03-04-2018.csv"
"""
function parse_temperature_range_and_sigma(temp_str::AbstractString)
    # Default values
    default_optimum = 15.0
    default_sigma = 3.0

    if isempty(temp_str) || temp_str == ""
        return (optimum=default_optimum, sigma=default_sigma,
            lower=-Inf, upper=Inf)
    end

    # Handle "X to Y" format
    if occursin(" to ", temp_str)
        parts = split(temp_str, " to ")
        if length(parts) == 2
            try
                t_min = parse(Float64, strip(parts[1]))
                t_max = parse(Float64, strip(parts[2]))
                thermal_range = t_max - t_min
                # Sigma is roughly 1/6 of the temperature range
                # This ensures ~95% of the thermal niche falls within ±2σ
                sigma = thermal_range / 6.0
                optimum = (t_min + t_max) / 2.0
                return (optimum=optimum, sigma=sigma, lower=t_min, upper=t_max)
            catch
                return (optimum=default_optimum, sigma=default_sigma,
                    lower=-Inf, upper=Inf)
            end
        end
    end

    # Try to parse as single number - no range info, use default sigma
    try
        optimum = parse(Float64, temp_str)
        return (optimum=optimum, sigma=default_sigma, lower=-Inf, upper=Inf)
    catch
        return (optimum=default_optimum, sigma=default_sigma,
            lower=-Inf, upper=Inf)
    end
end

"""
    parse_elevation(elev_str::AbstractString)

Parse elevation string, handling ranges and missing values.
Returns a single elevation value (midpoint if range, default if missing).

Input file is "data/ABIOTIC/caracteristicas_peces_Guadalquivir_03-04-2018.csv"
"""
function parse_elevation(elev_str::AbstractString)
  if isempty(elev_str) || elev_str == ""
    return 500.0
  end
  if occursin(" to ", elev_str)
    parts = split(elev_str, " to ")
    if length(parts) == 2
      try
        e_min = parse(Float64, strip(parts[1]))
        e_max = parse(Float64, strip(parts[2]))
        return (e_min + e_max) / 2.0
      catch
        return 500.0
      end
    end
  end
  try
    return parse(Float64, elev_str)
  catch
    return 500.0
  end
end

"""
    load_species_characteristics(file::String)

Load species ecological characteristics from the species traits file.
Returns a DataFrame with species codes and their thermal parameters.

NB: This function takes midpoint values for temperature and elevation ranges, and parses max size for growth rate scaling. The thermal breadth (sigma) is derived from the temperature range (sigma ≈ range/6).
"""
function load_species_characteristics(file::String)
    println("Loading species characteristics from: $file")

    # Read the semicolon-delimited file
    df = CSV.read(file, DataFrame; delim=';')

    # Parse temperature ranges to get thermal optima, sigma (thermal breadth)
    # and the empirical lower/upper limits used by the heat-stress term.
    thermal_params = parse_temperature_range_and_sigma.(string.(df.TEMPERATURE_C))
    df.thermal_optimum = [p.optimum for p in thermal_params]
    df.thermal_sigma = [p.sigma for p in thermal_params]
    df.thermal_lower = [p.lower for p in thermal_params]
    df.thermal_upper = [p.upper for p in thermal_params]

    df.elevation_optimum = parse_elevation.(string.(df.ELEVATION_m))

    # Parse max size for scaling growth rates
    df.max_size_mm = Float64.(df.MAX_SIZE_mm)

    println("Loaded characteristics for $(nrow(df)) species")
    return df
end

"""
    load_site_data(connectivity_file::String, environmental_file::String)

Load site-level data from connectivity and environmental files.
Returns a DataFrame with site information including coordinates and elevation.
"""
function load_site_data(connectivity_file::String, environmental_file::String)
    println("Loading site data from connectivity file: $connectivity_file")

    # Load connectivity data (small file)
    connectivity_df = CSV.read(connectivity_file, DataFrame)

    # Filter out site 1.30.20 which lacks reticular distance info
    connectivity_df = filter(row -> row.CODIGO != "1.30.20", connectivity_df)

    println("Loaded $(nrow(connectivity_df)) sites from connectivity data")

    # Load environmental data (medium file)
    println("Loading environmental data from: $environmental_file")
    environmental_df = CSV.read(environmental_file, DataFrame)
    println("Loaded $(nrow(environmental_df)) sites from environmental data")

    # Merge on CODIGO
    site_df = innerjoin(connectivity_df, environmental_df, on=:CODIGO, makeunique=true)

    # Handle "No existe" values in dam distance columns
    site_df.Demb_arr_m = replace(site_df."Demb arr.(m)", "No existe" => "0")
    site_df.Demb_ab_m = replace(site_df."Demb ab.(m)", "No existe" => "0")

    # Convert to numeric
    site_df.Demb_arr_m = parse.(Float64, site_df.Demb_arr_m)
    site_df.Demb_ab_m = parse.(Float64, site_df.Demb_ab_m)

    println("Merged data for $(nrow(site_df)) sites")

    return site_df
end

"""
    load_species_density_data(density_file::String)

Load species density data from the fish density matrix.
Returns a DataFrame with sites as rows and species densities as columns.
"""
function load_species_density_data(density_file::String)
    println("Loading species density data from: $density_file")

    # Read density matrix
    density_df = CSV.read(density_file, DataFrame)

    # Get species density columns (ending with _DEN)
    density_cols = [c for c in names(density_df) if endswith(c, "_DEN")]

    # Extract just the species codes
    species_codes = [replace(c, "_DEN" => "") for c in density_cols]

    println("Loaded density data for $(length(species_codes)) species at $(nrow(density_df)) sites")

    return density_df, species_codes
end

"""
    _normalise_cedex_scenario(value)

Return a comparable scenario code while preserving the original labels in the
loaded data.  The two CEDEX files currently use labels such as `SSP 585` and
`SSP285`, so whitespace is removed but no spelling correction is attempted.
"""
function _normalise_cedex_scenario(value)
    return uppercase(replace(strip(string(value)), r"\s+" => ""))
end

"""
    load_cedex_var(file::AbstractString)

Load the 2045 CEDEX basin-level projection table.  The source contains rows
for precipitation (`PRE`), potential and actual evapotranspiration (`ETP`,
`ETR`), recharge (`REC`), and runoff (`ESC`) under two scenarios.  The columns
`A`--`U` are source region/classes and are intentionally not mapped to model
sites here.  Duplicate source headers are retained by CSV.jl as unique names
(e.g. `MAX` and `MAX_1`).
"""
function load_cedex_var(file::AbstractString)
    println("Loading CEDEX basin projections from: $file")
    df = CSV.read(file, DataFrame; normalizenames=false)
    required = ["VAR", "ESCEN"]
    missing_columns = setdiff(required, names(df))
    isempty(missing_columns) || error("CEDEX VAR file is missing columns: $(join(missing_columns, ", "))")

    # Keep source values and add a stable key for joins/lookups.  Do not infer
    # whether these values are percentages, millimetres, or temperatures.
    df.cedex_variable = uppercase.(string.(df[!, "VAR"]))
    df.cedex_scenario = _normalise_cedex_scenario.(df[!, "ESCEN"])
    numeric_columns = setdiff(names(df), ["VAR", "ESCEN", "cedex_variable", "cedex_scenario"])
    for column in numeric_columns
        df[!, column] = passmissing(Float64).(df[!, column])
    end

    println("Loaded $(nrow(df)) CEDEX basin projection rows")
    return df
end

"""
    load_cedex_esc_uts(file::AbstractString)

Load the two-row-header CEDEX seasonal table in tidy form.  Each output row
contains one UTS, scenario, measure (`percent` or `mm`), season, and value.
`Media` and `Mediana` rows are retained with `row_type == \"summary\"`; callers
can exclude them when constructing spatial inputs.
"""
function load_cedex_esc_uts(file::AbstractString)
    println("Loading CEDEX seasonal UTS projections from: $file")
    raw = CSV.read(file, DataFrame; header=false, types=String,
        normalizenames=false, silencewarnings=true)
    ncol(raw) >= 17 || error("CEDEX UTS file must contain at least 17 columns")

    seasons = ["OND", "EFM", "AMJ", "JAS"]
    # Source columns 2:17 are four groups of four seasonal values.  The first
    # header row carries the scenario and measure for each group.
    group_headers = [string(raw[1, 2 + 4 * group]) for group in 0:3]
    out = DataFrame(
        uts=String[], scenario=String[], measure=String[], season=String[],
        value=Union{Missing,Float64}[], row_type=String[], source_row=Int[]
    )

    for row_number in 3:nrow(raw)
        uts_value = raw[row_number, 1]
        ismissing(uts_value) && continue
        uts = strip(string(uts_value))
        isempty(uts) && continue
        row_type = uppercase(uts) in ("MEDIA", "MEDIANA") ? "summary" : "spatial"

        for group in 0:3
            header = strip(group_headers[group + 1])
            isempty(header) && continue
            tokens = split(header)
            length(tokens) >= 3 || error("Unexpected CEDEX UTS header '$header'")
            scenario = _normalise_cedex_scenario(tokens[2])
            measure = lowercase(tokens[3]) == "%" ? "percent" : lowercase(tokens[3])

            for season_index in 1:4
                cell = raw[row_number, 2 + group * 4 + season_index - 1]
                value = missing
                if !ismissing(cell) && !isempty(strip(string(cell)))
                    value = try
                        parse(Float64, strip(string(cell)))
                    catch
                        missing
                    end
                end
                push!(out, (uts, scenario, measure, seasons[season_index],
                    value, row_type, row_number))
            end
        end
    end

    println("Loaded $(nrow(out)) tidy CEDEX UTS projection values")
    return out
end

"""
    select_cedex_uts(uts_df; uts, scenario, measure, season,
                     include_summary=false)

Select a single CEDEX UTS value.  This helper deliberately requires an
explicit UTS, scenario, measure, and season so that a seasonal projection is
never silently applied to all sites.
"""
function select_cedex_uts(uts_df::DataFrame; uts::AbstractString,
                          scenario::AbstractString, measure::AbstractString,
                          season::AbstractString, include_summary::Bool=false)
    scenario_code = _normalise_cedex_scenario(scenario)
    measure_code = lowercase(strip(measure))
    season_code = uppercase(strip(season))
    selected = filter(row -> row.uts == uts &&
        row.scenario == scenario_code && row.measure == measure_code &&
        row.season == season_code && (include_summary || row.row_type == "spatial"), uts_df)
    nrow(selected) == 1 || error("Expected one CEDEX UTS value, found $(nrow(selected))")
    return selected[1, :]
end

"""
    load_obstacles(file::AbstractString)

Load the 2026 inventory of obstacles that are not completely passable.  Raw
Spanish categorical fields are retained.  Parsed convenience columns are
added for coordinates, height, and the source `IF` index.  `IF` is not
interpreted as a passability probability by this loader.
"""
function load_obstacles(file::AbstractString)
    println("Loading obstacle inventory from: $file")
    df = CSV.read(file, DataFrame; normalizenames=false)
    required = ["ID_CLAVE_MGM", "COORD_X", "COORD_Y", "TRAMO_COD"]
    missing_columns = setdiff(required, names(df))
    isempty(missing_columns) || error("Obstacle file is missing columns: $(join(missing_columns, ", "))")

    parse_optional_float(value) = begin
        if ismissing(value) || isempty(strip(string(value)))
            missing
        else
            text = strip(string(value))
            uppercase(text) in ("NA", "SD", "DE", "NO", "N/D", "-") ? missing :
                try
                    parse(Float64, replace(text, ',' => '.'))
                catch
                    missing
                end
        end
    end

    df.coord_x_m = parse_optional_float.(df[!, "COORD_X"])
    df.coord_y_m = parse_optional_float.(df[!, "COORD_Y"])
    df.height_m = parse_optional_float.(df[!, "Altura"])
    df.if_index = parse_optional_float.(df[!, "IF"])
    df.has_valid_coordinates = .!(ismissing.(df.coord_x_m) .| ismissing.(df.coord_y_m))

    # E12 classification inputs: reservoir capacity (de-duplication key) and the
    # normalised structure status.  `IF` is still retained uninterpreted.
    if hasproperty(df, :Cap_emba)
        df.cap_emba_m3 = parse_optional_float.(df[!, "Cap_emba"])
    end
    if hasproperty(df, :ESTADO)
        df.status_code = [_normalise_obstacle_status(v) for v in df[!, "ESTADO"]]
        df.status_class = [_classify_obstacle_status(v)[1] for v in df[!, "ESTADO"]]
    end

    println("Loaded $(nrow(df)) obstacles; $(count(df.has_valid_coordinates)) have valid coordinates")
    return df
end

"""
    build_site_coordinate_matrix(site_df, sites)

Return an `n_sites × 2` matrix of UTM X/Y coordinates in model site order.
"""
function build_site_coordinate_matrix(site_df::DataFrame, sites::Vector{String})
    coordinate_lookup = Dict{String,Tuple{Float64,Float64}}()
    for row in eachrow(site_df)
        x = hasproperty(row, :UTMX) ? row.UTMX : missing
        y = hasproperty(row, :UTMY) ? row.UTMY : missing
        if !ismissing(x) && !ismissing(y)
            coordinate_lookup[string(row.CODIGO)] = (Float64(x), Float64(y))
        end
    end

    coordinates = fill(NaN, length(sites), 2)
    for (index, site) in enumerate(sites)
        if haskey(coordinate_lookup, site)
            coordinates[index, :] .= coordinate_lookup[site]
        end
    end
    return coordinates
end

# =============================================================================
# E12: obstacle passability on the selected (legacy or corrected) network.
#
# Status classification (documented criteria; no MGM codebook ships with the
# inventory, so the code meanings below are inferred from the source's own
# free-text descriptions and recorded explicitly):
#   * non-operational / removed -> excluded from the barrier overlay:
#       AB = abandoned, DE = demolished, AR = ruined, FS = out of service, plus
#       the spelled-out "Abandonado en buen estado" / "Abandonado en ruinas".
#   * operational -> kept: EX (existing / in service) and "En explotación".
#   * genuinely ambiguous -> annotated and (by default) excluded, never
#       silently dropped: SC, DM, ND, OT and a missing/blank or unknown code.
# Out-of-basin structures are those whose water-body code (CODMAS) does not
# start with the Guadalquivir basin-district prefix ES050.
#
# De-duplication: the legacy dam layer is built from per-site distances to the
# nearest upstream/downstream embalse and exposes no barrier identifiers, so an
# exact ID join is impossible.  The documented de-duplication key is therefore
# "reservoir-class obstacle on a directional link the legacy dam layer already
# restricts": an obstacle is a legacy duplicate only when its `Cap_emba`
# reservoir capacity parses as a finite number AND `legacy_dams` already
# restricts that directional link.  Such obstacles are excluded from the
# overlay so the same physical dam is not applied twice.
#
# The inventory's `IF` passability index is deliberately NOT translated to a
# passability: its scale/direction is unconfirmed and 547/1,658 rows are
# non-numeric.  Uniform `obstacle_passability` /
# `obstacle_downstream_passability` are used instead (documented limitation).
# =============================================================================

const OBSTACLE_NONOPERATIONAL_STATUS = Set([
    "AB", "DE", "AR", "FS", "ABANDONADOENBUENESTADO", "ABANDONADOENRUINAS",
])
const OBSTACLE_AMBIGUOUS_STATUS = Set(["SC", "DM", "ND", "OT"])
const OBSTACLE_ACTIVE_STATUS = Set(["EX", "ENEXPLOTACION"])
const GUADALQUIVIR_BASIN_PREFIX = "ES050"

"""
    _normalise_obstacle_status(value) -> String

Uppercase, strip and de-accent an `ESTADO` cell so classification does not
depend on the source's casing or accents.  Returns `""` for a missing/blank
cell.
"""
function _normalise_obstacle_status(value)
    (value === missing || isempty(strip(string(value)))) && return ""
    raw = string(value)
    # The inventory CSV is Latin-1, so a raw accented byte is not valid UTF-8
    # and `uppercase` raises `InvalidCharError`.  Decode through Latin-1 in
    # that case; valid UTF-8/ASCII inputs are left untouched.
    text = try
        uppercase(strip(raw))
    catch
        uppercase(strip(String([Char(b) for b in codeunits(raw)])))
    end
    text = replace(text, "Á" => "A", "É" => "E", "Í" => "I", "Ó" => "O",
        "Ú" => "U", "Ü" => "U", "Ñ" => "N")
    return replace(text, r"\s+" => "")
end

"""
    _classify_obstacle_status(value) -> (class::String, code::String)

Deterministically classify a raw `ESTADO` value into `"active"`,
`"nonoperational"` or `"ambiguous"`.  Unknown, blank and missing codes are
reported as `"ambiguous"` so they are annotated rather than silently dropped.
"""
function _classify_obstacle_status(value)
    code = _normalise_obstacle_status(value)
    code in OBSTACLE_ACTIVE_STATUS && return ("active", code)
    code in OBSTACLE_NONOPERATIONAL_STATUS && return ("nonoperational", code)
    return ("ambiguous", code)
end

function _obstacle_status_class(row)
    hasproperty(row, :status_class) && return string(row.status_class)
    hasproperty(row, :ESTADO) && return _classify_obstacle_status(row.ESTADO)[1]
    return "active"
end

function _obstacle_status_code(row)
    hasproperty(row, :status_code) && return string(row.status_code)
    hasproperty(row, :ESTADO) && return _classify_obstacle_status(row.ESTADO)[2]
    return ""
end

function _obstacle_capacity(row)
    hasproperty(row, :cap_emba_m3) && return row.cap_emba_m3
    if hasproperty(row, :Cap_emba)
        text = strip(string(row.Cap_emba))
        return try
            parse(Float64, replace(text, ',' => '.'))
        catch
            missing
        end
    end
    return missing
end

"""
    build_obstacle_passability_matrix(obstacles, site_df, sites, distances, elevations;
                                      matching_tolerance=2000.0,
                                      obstacle_passability=0.1,
                                      obstacle_downstream_passability=0.5,
                                      tree_edges=Tuple{Int,Int,Float64}[],
                                      legacy_dams=nothing,
                                      include_ambiguous_status=false,
                                      basin_prefix="ES050")

Create an obstacle overlay for the selected site network.  Each obstacle is
matched to the nearest straight line segment between connected site pairs,
using UTM coordinates and the supplied tolerance.  On the corrected
(`connectivity_method = :on_path`) network pass `tree_edges`, the directed
`(upstream, downstream, metres)` interface returned by
[`build_on_path_distance_matrix`](@ref); obstacles are then snapped only to
those links and the tree orientation defines the movement direction.  The
legacy call (no `tree_edges`) keeps the original straight-line matching and the
elevation-based direction.

Barriers ACCUMULATE along a link: the effective passability is the PRODUCT of
the individual obstacle passabilities on it (two 0.1 barriers give 0.01),
separately for the upstream and downstream directions.  A link's upstream
factor is applied to movement towards the higher/tree-upstream site and the
downstream factor to the reverse movement (downstream stays more accessible).
Equal-elevation legacy edges are left unchanged because their direction cannot
be inferred.

Non-operational structures (demolished/abandoned), out-of-basin structures,
legacy duplicates and (by default) ambiguous-status structures are excluded
from the overlay; every obstacle is still annotated in the diagnostics with
its `status_class` and `outcome`, and the counts are returned in `metadata`.
Ties between equidistant links are broken deterministically by site code, and
obstacles are processed in `ID_CLAVE_MGM` order so results never depend on file
or `Dict` order.  The `IF` index is not interpreted (see the section comment).

Returns `(passability=overlay, diagnostics=DataFrame, metadata=NamedTuple)`.
"""
function build_obstacle_passability_matrix(obstacles::DataFrame, site_df::DataFrame,
                                            sites::Vector{String}, distances,
                                            elevations::AbstractVector;
                                            matching_tolerance::Float64=2000.0,
                                            obstacle_passability::Float64=0.1,
                                            obstacle_downstream_passability::Float64=0.5,
                                            tree_edges::Vector{Tuple{Int,Int,Float64}}=Tuple{Int,Int,Float64}[],
                                            legacy_dams::Union{Nothing,AbstractMatrix}=nothing,
                                            include_ambiguous_status::Bool=false,
                                            basin_prefix::AbstractString=GUADALQUIVIR_BASIN_PREFIX)
    matching_tolerance >= 0 || error("matching_tolerance must be non-negative")
    0.0 <= obstacle_passability <= 1.0 || error("obstacle_passability must be between 0 and 1")
    0.0 <= obstacle_downstream_passability <= 1.0 || error("obstacle_downstream_passability must be between 0 and 1")

    n_sites = length(sites)
    length(elevations) == n_sites || error("elevations must contain one value per site")
    coordinates = build_site_coordinate_matrix(site_df, sites)

    # Candidate links in a deterministic order.  On the corrected network the
    # caller supplies the directed tree edges; the legacy path derives the
    # links from the sparse distance matrix exactly as before.
    edges = Tuple{Int,Int}[]
    tree_direction = Dict{Tuple{Int,Int},Tuple{Int,Int}}()   # (a,b) -> (up, down)
    if !isempty(tree_edges)
        edge_set = Set{Tuple{Int,Int}}()
        for (u, v, _) in tree_edges
            (1 <= u <= n_sites && 1 <= v <= n_sites) ||
                error("tree edge ($u, $v) is out of range for $n_sites sites")
            u == v && error("tree edge ($u, $v) is a self-loop")
            a, b = minmax(u, v)
            push!(edge_set, (a, b))
            tree_direction[(a, b)] = (u, v)   # u upstream, v downstream
        end
        edges = sort!(collect(edge_set))
    else
        for i in 1:n_sites, j in (i + 1):n_sites
            if (isfinite(coordinates[i, 1]) && isfinite(coordinates[i, 2]) &&
                isfinite(coordinates[j, 1]) && isfinite(coordinates[j, 2]) &&
                (distances[i, j] > 0 || distances[j, i] > 0))
                push!(edges, (i, j))
            end
        end
    end

    edge_tie_key(i, j) = sites[i] <= sites[j] ? (sites[i], sites[j]) : (sites[j], sites[i])

    overlay = ones(Float64, n_sites, n_sites)
    obstacle_ids = String[]
    matched = Bool[]
    edge_from = String[]
    edge_to = String[]
    restricted_origin = String[]
    restricted_destination = String[]
    match_distance_m = Union{Missing,Float64}[]
    status = String[]
    status_class = String[]
    outcome = String[]
    duplicate_flag = Bool[]
    cap_values = Union{Missing,Float64}[]
    link_obstacle_counts = Dict{Tuple{String,String},Int}()

    n_obs = nrow(obstacles)
    id_column = hasproperty(obstacles, :ID_CLAVE_MGM) ?
        string.(obstacles[!, :ID_CLAVE_MGM]) : string.(1:n_obs)
    order = sortperm(id_column)

    for index in order
        row = obstacles[index, :]
        obstacle_id = id_column[index]
        push!(obstacle_ids, obstacle_id)

        class = _obstacle_status_class(row)
        code = _obstacle_status_code(row)
        cap = _obstacle_capacity(row)
        push!(status, code)
        push!(status_class, class)
        push!(cap_values, cap)

        raw_codmas = hasproperty(row, :CODMAS) ? string(row.CODMAS) : ""
        out_of_basin = !isempty(strip(raw_codmas)) &&
            !startswith(uppercase(strip(raw_codmas)), uppercase(string(basin_prefix)))

        x = hasproperty(row, :coord_x_m) ? row.coord_x_m : missing
        y = hasproperty(row, :coord_y_m) ? row.coord_y_m : missing

        best_edge = nothing
        best_distance = Inf
        best_key = (Inf, "", "")
        if !ismissing(x) && !ismissing(y) && isfinite(Float64(x)) && isfinite(Float64(y))
            for (i, j) in edges
                x1, y1 = coordinates[i, 1], coordinates[i, 2]
                x2, y2 = coordinates[j, 1], coordinates[j, 2]
                (isfinite(x1) && isfinite(y1) && isfinite(x2) && isfinite(y2)) || continue
                dx, dy = x2 - x1, y2 - y1
                length_squared = dx * dx + dy * dy
                length_squared == 0 && continue
                projection = ((Float64(x) - x1) * dx + (Float64(y) - y1) * dy) / length_squared
                # Do not attach an obstacle to an unrelated edge through an endpoint.
                (projection < 0.0 || projection > 1.0) && continue
                closest_x = x1 + projection * dx
                closest_y = y1 + projection * dy
                distance_to_edge = hypot(Float64(x) - closest_x, Float64(y) - closest_y)
                # Deterministic tie-break: nearer edge first, then site codes.
                key = (distance_to_edge, edge_tie_key(i, j)...)
                if key < best_key
                    best_key = key
                    best_distance = distance_to_edge
                    best_edge = (i, j)
                end
            end
        end

        is_matched = best_edge !== nothing && best_distance <= matching_tolerance
        push!(matched, is_matched)

        up = 0
        down = 0
        is_duplicate = false
        if is_matched
            i, j = best_edge
            push!(edge_from, sites[i]); push!(edge_to, sites[j])
            push!(match_distance_m, best_distance)

            # Direction: corrected tree orientation when available, otherwise
            # the legacy elevation rule.
            if haskey(tree_direction, (i, j))
                up, down = tree_direction[(i, j)]
            elseif elevations[i] > elevations[j]
                up, down = i, j
            elseif elevations[j] > elevations[i]
                up, down = j, i
            end

            if cap !== missing && isfinite(Float64(cap)) && legacy_dams !== nothing && up != 0
                is_duplicate = legacy_dams[up, down] < 1.0 || legacy_dams[down, up] < 1.0
            end
        else
            push!(edge_from, ""); push!(edge_to, "")
            push!(match_distance_m, best_edge === nothing ? missing : best_distance)
        end
        push!(duplicate_flag, is_duplicate)

        # Outcome precedence is fixed so the counts are additive and
        # independent of whether the structure also snapped to a link.
        if out_of_basin
            origin, destination = "", ""
            push!(outcome, "excluded_out_of_basin")
        elseif class == "nonoperational"
            origin, destination = "", ""
            push!(outcome, "excluded_nonoperational")
        elseif is_duplicate
            origin, destination = "", ""
            push!(outcome, "excluded_duplicate_legacy")
        elseif class == "ambiguous" && !include_ambiguous_status
            origin, destination = "", ""
            push!(outcome, "excluded_ambiguous_status")
        elseif !is_matched
            origin, destination = "", ""
            push!(outcome, "unmatched")
        elseif up == 0
            origin, destination = "", ""
            push!(outcome, "no_direction")
        else
            # Accumulate barriers: the effective passability is the product of
            # the individual obstacle passabilities in each direction.
            overlay[up, down] *= obstacle_passability
            overlay[down, up] *= obstacle_downstream_passability
            origin, destination = sites[down], sites[up]
            push!(outcome, "applied")
            key = edge_tie_key(i, j)
            link_obstacle_counts[key] = get(link_obstacle_counts, key, 0) + 1
        end
        push!(restricted_origin, origin)
        push!(restricted_destination, destination)
    end

    diagnostics = DataFrame(
        obstacle_id=obstacle_ids, matched=matched, edge_from=edge_from,
        edge_to=edge_to, restricted_origin=restricted_origin,
        restricted_destination=restricted_destination,
        match_distance_m=match_distance_m, status=status,
        status_class=status_class, outcome=outcome,
        duplicate_of_legacy_dam=duplicate_flag, cap_emba_m3=cap_values
    )

    tally(name) = count(==(name), outcome)
    ambiguous_statuses = sort(unique([status[k] for k in eachindex(status)
        if status_class[k] == "ambiguous"]))
    ambiguous_statuses = [isempty(code) ? "<blank>" : code for code in ambiguous_statuses]
    nonoperational_statuses = sort(unique([status[k] for k in eachindex(status)
        if status_class[k] == "nonoperational"]))
    metadata = (
        n_total = n_obs,
        n_matched = count(matched),
        n_unmatched = tally("unmatched"),
        n_applied = tally("applied"),
        n_no_direction = tally("no_direction"),
        n_excluded_out_of_basin = tally("excluded_out_of_basin"),
        n_excluded_nonoperational = tally("excluded_nonoperational"),
        n_excluded_ambiguous_status = tally("excluded_ambiguous_status"),
        n_excluded_duplicate_legacy = tally("excluded_duplicate_legacy"),
        n_excluded_total = tally("excluded_out_of_basin") +
            tally("excluded_nonoperational") + tally("excluded_ambiguous_status") +
            tally("excluded_duplicate_legacy"),
        n_status_active = count(==("active"), status_class),
        n_status_nonoperational = count(==("nonoperational"), status_class),
        n_status_ambiguous = count(==("ambiguous"), status_class),
        ambiguous_statuses = ambiguous_statuses,
        nonoperational_statuses = nonoperational_statuses,
        include_ambiguous_status = include_ambiguous_status,
        basin_prefix = string(basin_prefix),
        duplicate_legacy_key = "reservoir capacity (Cap_emba numeric) on a directional link already restricted by legacy_dams",
        if_index_used_for_passability = false,
        use_tree_edges = !isempty(tree_edges),
        link_obstacle_counts = link_obstacle_counts
    )
    return (passability=overlay, diagnostics=diagnostics, metadata=metadata)
end

# =============================================================================
# Interaction parsing (issue C1)
#
# Conventions implemented here:
#  * The model equation reads `α[s, j]` as the effect of species `j` (source)
#    on species `s` (target), stored at `α[target, source]`.
#  * A CSV cell at (row A, column B) describes the relationship between the two
#    species A and B.  The species named in the cell text is the AFFECTED one
#    (the target); the other species of the pair is the source.
#  * Undirected competition/interference (no named species) is symmetric: both
#    directions receive the mechanism value.
#  * `"No coexist"` denotes allopatry, i.e. ABSENCE of interaction (α = 0):
#    co-occurrence is already governed by the thermal and habitat filters.
#  * Cells that combine an undirected clause and a directed clause are applied
#    as both; the more specific directed clause defines the coefficient in its
#    own direction, and every directed effect of such a cell is flagged
#    `ambiguous = true` in the long-format table.
# =============================================================================

const _DIRECTED_VERB_REGEX =
    r"(?:affects|displaces|interferes|interfere|interfiere)\s+([A-Za-z]{1,4})\s+through"
const _BARE_TARGET_REGEX = r"\band\s+([A-Za-z]{1,4})\s+through"
const _DEPREDA_REGEX = r"([A-Za-z]{1,4})\s+depreda\s+([A-Za-z]{1,4})"
const _UNDIRECTED_REGEX =
    r"(?:affects|displaces|interferes|interfere|interfiere)\s+through"

"""
    _interaction_mechanism_name(text::AbstractString)

Return the qualitative mechanism category implied by `text`, using the same
priority order as [`parse_interaction_string`](@ref).
"""
function _interaction_mechanism_name(text::AbstractString)
    lower = lowercase(string(text))
    occursin("no coexist", lower) && return "none"
    occursin("displaces", lower) && return "displacement"
    occursin("predation", lower) && return "predation"
    if occursin("competition", lower) || occursin("interfere", lower) ||
       occursin("interfiere", lower)
        return "competition"
    end
    occursin("affects", lower) && return "effect"
    return "none"
end

"""
    parse_interaction_cell(text, row_code, col_code) -> Vector{NamedTuple}

Parse a single matrix cell into one directed effect per `(source, target)` pair.

`row_code` and `col_code` are the species of the cell's pair.  A species named
in the text is the target (the affected species); the other member of the pair
is the source.  Returned rows are `(source, target, mechanism, alpha,
ambiguous, raw_text)` with species codes uppercased.

Rules:
- `"<verb> <Code> through <mechanism>"` (verb ∈ affects/displaces/interfere/
  interferes/interfiere) names the target.
- `"<Code> through <mechanism>"` after `"and"` (no verb) names the target.
- `"<CodeA> depreda <CodeB>"` is predation of CodeA (source) on CodeB (target).
- An undirected clause (verb directly followed by `through`) applies to both
  directions.
- `"No coexist"` yields no effect at all (allopatry → α = 0).
"""
function parse_interaction_cell(text, row_code, col_code)
    raw = ismissing(text) ? "" : string(text)
    stripped = strip(raw)
    effects = NamedTuple[]
    isempty(stripped) && return effects

    lower = lowercase(stripped)
    occursin("no coexist", lower) && return effects

    row = lowercase(string(row_code))
    col = lowercase(string(col_code))

    # (source, target, value, mechanism, specific?)
    clauses = Tuple{String,String,Float64,String,Bool}[]

    for m in eachmatch(_DIRECTED_VERB_REGEX, lower)
        named = lowercase(m.captures[1])
        (named == row || named == col) || continue
        source = named == row ? col : row
        push!(clauses, (source, named, parse_interaction_string(lower),
                        _interaction_mechanism_name(lower), true))
    end

    for m in eachmatch(_BARE_TARGET_REGEX, lower)
        named = lowercase(m.captures[1])
        (named == row || named == col) || continue
        source = named == row ? col : row
        clause = lower[m.offset:end]
        push!(clauses, (source, named, parse_interaction_string(clause),
                        _interaction_mechanism_name(clause), true))
    end

    for m in eachmatch(_DEPREDA_REGEX, lower)
        source = lowercase(m.captures[1])
        target = lowercase(m.captures[2])
        (target == row || target == col) || continue
        push!(clauses, (source, target, -0.5, "predation", true))
    end

    has_undirected = occursin(_UNDIRECTED_REGEX, lower)
    if has_undirected
        m = match(_UNDIRECTED_REGEX, lower)
        clause = lower[m.offset:end]
        value = parse_interaction_string(clause)
        mechanism = _interaction_mechanism_name(clause)
        push!(clauses, (row, col, value, mechanism, false))
        push!(clauses, (col, row, value, mechanism, false))
    end

    ambiguous = has_undirected && any(c -> c[5], clauses)

    resolved = Dict{Tuple{String,String}, Tuple{Float64,String,Bool}}()
    for (source, target, value, mechanism, specific) in clauses
        key = (source, target)
        if !haskey(resolved, key) || (specific && !resolved[key][3])
            resolved[key] = (value, mechanism, specific)
        end
    end

    for ((source, target), (value, mechanism, _)) in resolved
        push!(effects, (source=uppercase(source), target=uppercase(target),
                        mechanism=mechanism, alpha=value, ambiguous=ambiguous,
                        raw_text=stripped))
    end
    return effects
end

"""
    build_interaction_long_table(interaction_file; species_codes=nothing) -> DataFrame

Parse a semicolon-delimited interaction matrix into an explicit long-format
table with columns `source, target, mechanism, alpha, ambiguous, raw_text`.

One row is emitted per directed effect (i.e. per `(source, target)` pair), so an
undirected symmetric relationship contributes two rows.  Rows belonging to a
cell that mixes an undirected clause with a directed clause are flagged
`ambiguous = true`.

The loader [`load_interaction_matrix`](@ref) reads this table when a sibling
`*_long.csv` (or `interaction_matrix_long.csv`) exists.
"""
function build_interaction_long_table(interaction_file::String;
        species_codes::Union{Nothing,AbstractVector{<:AbstractString}}=nothing)
    interaction_df = CSV.read(interaction_file, DataFrame; delim=';')
    rename!(interaction_df, 1 => :Species)

    row_labels = string.(interaction_df.Species)
    column_codes = String.(names(interaction_df)[2:end])

    records = NamedTuple[]
    for (i, row) in enumerate(eachrow(interaction_df))
        for col_code in column_codes
            for effect in parse_interaction_cell(row[Symbol(col_code)], row_labels[i], col_code)
                push!(records, effect)
            end
        end
    end

    if isempty(records)
        return DataFrame(source=String[], target=String[], mechanism=String[],
                         alpha=Float64[], ambiguous=Bool[], raw_text=String[])
    end

    table = DataFrame(records)
    if species_codes === nothing
        order = lowercase.(unique(vcat(row_labels, column_codes)))
    else
        order = lowercase.(collect(species_codes))
    end
    rank = Dict(code => i for (i, code) in enumerate(order))
    ranks = [(get(rank, lowercase(r.source), typemax(Int)),
              get(rank, lowercase(r.target), typemax(Int))) for r in eachrow(table)]
    perm = sortperm(ranks)
    return table[perm, :]
end

"""
    _interaction_long_path(interaction_file::String) -> Union{String,Nothing}

Return the path of the long-format table to read for `interaction_file`, or
`nothing` when none is present.  Prefers a sibling with the same basename plus
`_long.csv`, then the canonical `interaction_matrix_long.csv` in the same
directory.
"""
function _interaction_long_path(interaction_file::String)
    same_basename = splitext(interaction_file)[1] * "_long.csv"
    isfile(same_basename) && return same_basename

    canonical = joinpath(dirname(interaction_file), "interaction_matrix_long.csv")
    isfile(canonical) && return canonical

    return nothing
end

"""
    load_interaction_matrix(interaction_file::String, species_codes::Vector{String})

Load and parse the species interaction matrix.

Returns a numeric matrix `α` where `α[target, source]` is the effect of
`species_codes[source]` on `species_codes[target]`.  When a long-format table
is available next to `interaction_file` it is read directly; otherwise the
semicolon matrix is parsed from its cell text (see [`parse_interaction_cell`](@ref)).
`"No coexist"` cells are treated as no interaction (α = 0).  Row labels and
column names are matched case-insensitively and their union is used, so a row
without a matching column (e.g. *Lepomis gibbosus*) is not silently dropped.
"""
function load_interaction_matrix(interaction_file::String, species_codes::Vector{String})
    println("Loading interaction matrix from: $interaction_file")

    n_species = length(species_codes)
    lower_codes = lowercase.(species_codes)
    code_to_idx = Dict(lower_codes[i] => i for i in 1:n_species)
    interaction_matrix = zeros(n_species, n_species)

    long_file = _interaction_long_path(interaction_file)
    if long_file !== nothing
        println("Reading long-format interaction table: $long_file")
        long_df = CSV.read(long_file, DataFrame)
        for row in eachrow(long_df)
            source = lowercase(string(row.source))
            target = lowercase(string(row.target))
            (haskey(code_to_idx, source) && haskey(code_to_idx, target)) || continue
            interaction_matrix[code_to_idx[target], code_to_idx[source]] = Float64(row.alpha)
        end
        println("Created $(n_species)x$(n_species) interaction matrix from long table")
        return interaction_matrix
    end

    # Fall back to parsing the semicolon-delimited matrix directly.
    interaction_df = CSV.read(interaction_file, DataFrame; delim=';')
    rename!(interaction_df, 1 => :Species)

    row_labels = string.(interaction_df.Species)
    column_codes = String.(names(interaction_df)[2:end])

    for (i, row) in enumerate(eachrow(interaction_df))
        row_code = row_labels[i]
        for col_code in column_codes
            for effect in parse_interaction_cell(row[Symbol(col_code)], row_code, col_code)
                source = lowercase(effect.source)
                target = lowercase(effect.target)
                (haskey(code_to_idx, source) && haskey(code_to_idx, target)) || continue
                interaction_matrix[code_to_idx[target], code_to_idx[source]] = effect.alpha
            end
        end
    end

    println("Created $(n_species)x$(n_species) interaction matrix")
    return interaction_matrix
end

"""
    parse_interaction_string(interaction_str::String)

Parse an interaction string and return a numeric value.
- "No coexist" => 0.0 (allopatry: absence of interaction)
- "displaces" => strong negative (-0.8)
- "predation" or "affects ... predation" => negative (-0.5)
- "competition", "interfere", "interfiere", "affects ... competition" => negative (-0.3)
- "affects" (without specific mechanism) => moderate negative (-0.2)
- "coexist, neutral" => zero (0.0)
- "coexist" without negative qualifier => zero (0.0)

The function checks patterns in order of priority (most negative first).
"""
function parse_interaction_string(interaction_str::Union{String, Missing})
    # Handle missing or empty values
    if ismissing(interaction_str) || isempty(interaction_str)
        return 0.0
    end

    interaction_str = strip(string(interaction_str))

    # Empty after stripping or just semicolon
    if interaction_str == "" || interaction_str == ";"
        return 0.0
    end

    # Convert to lowercase for case-insensitive matching
    interaction_lower = lowercase(interaction_str)

    # 1. "No coexist" describes allopatry; the thermal and habitat filters
    #    already decide whether two species can co-occur, so this is no
    #    measurable interaction.
    if occursin("no coexist", interaction_lower)
        return 0.0
    end

    # 2. Strong negative: displaces (complete displacement)
    if occursin("displaces", interaction_lower)
        return -0.8
    end

    # 3. Predation (direct predation effect) - including "affects ... predation"
    if occursin("predation", interaction_lower)
        return -0.5
    end

    # 4. Competition or interference - including "affects ... competition"
    if occursin("competition", interaction_lower) ||
       occursin("interfere", interaction_lower) ||
       occursin("interfiere", interaction_lower)
        return -0.3
    end

    # 5. Moderate negative: affects (some effect but not complete displacement)
    # This catches "affects" when not followed by predation or competition
    if occursin("affects", interaction_lower)
        return -0.2
    end

    # 6. Neutral coexistence - "coexist" with neutral qualifier or alone
    if occursin("coexist", interaction_lower) || occursin("neutral", interaction_lower)
        return 0.0
    end

    # Default: no interaction (neutral)
    return 0.0
end

"""
    build_distance_matrix(distance_file::String, sites::Vector{String},
                         site_to_subcatchment::Dict{String, String},
                         site_to_river_distance::Dict{String, Float64},
                         site_to_elevation::Dict{String, Float64})

Distances between sites.
Build a sparse distance matrix from the distance file.
Only includes connections between ADJACENT sites within the same subcatchment,
plus connections between outlet sites (closest to main river) based on elevation.

This models the main river implicitly by connecting subcatchment outlets
based on their elevation, allowing fish to disperse between subcatchments.
"""
function build_distance_matrix(distance_file::String, sites::Vector{String},
                               site_to_subcatchment::Dict{T, T2},
                               site_to_river_distance::Dict{T3, Float64},
                               site_to_elevation::Dict{T4, Float64}) where T <: AbstractString where T2 <: AbstractString where T3 <: AbstractString where T4 <: AbstractString
    println("Building distance matrix from: $distance_file")
    println("Only including connections between adjacent sites in the same subcatchment")

    # Create site to index mapping
    site_to_idx = Dict(s => i for (i, s) in enumerate(sites))
    n_sites = length(sites)

    # Group sites by subcatchment
    subcatchment_to_sites = Dict{String, Vector{String}}()
    for site in sites
        if haskey(site_to_subcatchment, site)
            sc = site_to_subcatchment[site]
            if !haskey(subcatchment_to_sites, sc)
                subcatchment_to_sites[sc] = String[]
            end
            push!(subcatchment_to_sites[sc], site)
        end
    end

    # For each subcatchment, sort sites by distance to river and find adjacent pairs
    adjacent_pairs = Set{Tuple{String, String}}()  # (upstream, downstream) pairs

    # Track outlet sites (closest to main river) for each subcatchment
    outlet_sites = String[]

    for (sc, sites_in_sc) in subcatchment_to_sites
        # Sort sites by distance to river (ascending = downstream first)
        sorted_sites = sort(sites_in_sc, by=s -> get(site_to_river_distance, s, Inf))

        # Connect each site to its adjacent neighbor (upstream <-> downstream)
        for i in 1:(length(sorted_sites)-1)
            downstream = sorted_sites[i]      # Closer to river
            upstream = sorted_sites[i+1]       # Farther from river
            push!(adjacent_pairs, (upstream, downstream))
        end

        # The outlet site (closest to main river) connects to the main river network
        if length(sorted_sites) > 0
            outlet = sorted_sites[1]  # First = closest to river
            push!(outlet_sites, outlet)
        end
    end

    # Connect outlet sites to each other based on elevation
    # This models the main river implicitly - lower elevation outlets are "downstream"
    if length(outlet_sites) > 1
        # Sort outlets by elevation (ascending = downstream first)
        sorted_outlets = sort(outlet_sites, by=s -> get(site_to_elevation, s, Inf))

        for i in 1:(length(sorted_outlets)-1)
            downstream_outlet = sorted_outlets[i]
            upstream_outlet = sorted_outlets[i+1]
            push!(adjacent_pairs, (upstream_outlet, downstream_outlet))
        end
    end

    println("Found $(length(outlet_sites)) outlet sites connected to form main river network")

    println("Found $(length(adjacent_pairs)) adjacent site pairs in the network")

    # Initialize sparse matrix components
    I = Int[]
    J = Int[]
    V = Float64[]

    # Stream the distance file
    total_rows = 0
    valid_connections = 0

    reader = CSV.File(distance_file; delim=';')

    for row in reader
        total_rows += 1

        origin = row.ID_ORIGIN
        dest = row.ID_DESTINATION

        # Skip if origin or destination not in our site list
        if !haskey(site_to_idx, origin) || !haskey(site_to_idx, dest)
            continue
        end

        # Skip self-connections
        dist = row.RETICULAR_DIST
        if origin == dest || dist <= 0
            continue
        end

        # CRITICAL: Only include connections between adjacent sites in the same subcatchment
        # Check if this is an adjacent pair (either direction)
        if (origin, dest) ∈ adjacent_pairs || (dest, origin) ∈ adjacent_pairs
            i = site_to_idx[dest]
            j = site_to_idx[origin]

            push!(I, i)
            push!(J, j)
            push!(V, dist)
            valid_connections += 1
        end

        if total_rows % 1000000 == 0
            println("Processed $total_rows rows...")
        end
    end

    println("Processed $total_rows total distance records")
    println("Found $valid_connections valid adjacent connections")

    # Create sparse matrix
    distance_matrix = sparse(I, J, V, n_sites, n_sites)

    println("Distance matrix: $(nnz(distance_matrix)) non-zero entries out of $(n_sites*n_sites) possible")

    return distance_matrix
end

# ---------------------------------------------------------------------------
# C4: dendritic graph reconstructed from the network-distance matrix.
#
# The preferred route (snap sites to `SW_Line_4C_.shp`, orient the river lines
# and follow the channels) is not feasible with this layer: it is a WISE
# water-body *inventory* with no flow-direction attribute (no FromNode/ToNode,
# no order), no elevation (the geometry is 2D EPSG:25830, `has_z == false`) and
# no clean noded reach topology (3821 digitised parts for 360 multi-line
# features), so flow direction is ambiguous.  We therefore use the documented
# fallback: the on-path parent rule plus an explicit minimum-spanning-tree join
# of the remaining roots.  See `docs/GuadeX_Correction_Plan_Sept2026.md` (C4).
# ---------------------------------------------------------------------------

# Default tolerance of the on-path test `D(i,k) ≈ d(i) − d(k)`.  The two data
# sources (ConnectivityUTM and the reticular matrix) are not perfectly
# consistent, so a small relative + absolute slack is required.  These are the
# values of the reviewer's prototype (issue C4).
const ON_PATH_DEFAULT_RTOL = 0.02   # relative slack (fraction of D)
const ON_PATH_DEFAULT_ATOL = 200.0  # absolute slack (metres)

"""
    _consider_on_path_parent!(parent, parent_distance, i, k, dist, dg, sites,
                              rtol, atol)

Update the on-path parent of site `i` with candidate `k` (which must be
downstream, `dg[k] < dg[i]`) when the pair lies on one flow path
`|dist − (dg[i] − dg[k])| ≤ rtol·dist + atol` and the candidate is nearer (smaller
`dist`).  Exact ties are broken by site code so the result never depends on
`Dict` iteration order.
"""
function _consider_on_path_parent!(parent::Vector{Int}, parent_distance::Vector{Float64},
        i::Int, k::Int, dist::Float64, dg::Vector{Float64},
        sites::Vector{String}, rtol::Float64, atol::Float64)
    on_path = abs(dist - (dg[i] - dg[k])) <= rtol * dist + atol
    on_path || return nothing

    best = parent[i]
    if best == 0 || dist < parent_distance[i] ||
       (dist == parent_distance[i] && sites[k] < sites[best])
        parent[i] = k
        parent_distance[i] = dist
    end
    return nothing
end

"""
    build_on_path_distance_matrix(distance_file, sites, site_to_river_distance,
                                  site_to_elevation;
                                  rtol=0.02, atol=200.0,
                                  store_both_directions=true) -> NamedTuple

Reconstruct the dendritic dispersal graph from the network-distance matrix
(issue C4, fallback route).

Algorithm
1. Stream `distance_file` once and keep the reticular (network) distance of
   every unordered pair of model sites.
2. On-path parent rule: site `i` is connected to the *nearest* site `k` such
   that `dg[k] < dg[i]` (downstream, with `dg` the distance to the Guadalquivir)
   and `|D(i,k) − (dg[i] − dg[k])| ≤ rtol·D(i,k) + atol`, i.e. `k` lies on the
   same flow path from `i` to the river mouth.  Sites without such a candidate
   are roots.
3. The roots are joined by a minimum spanning tree on the network distance
   (Kruskal; edges ordered by `(weight, site-code pair)` for reproducibility).
   Each root edge is oriented from the endpoint farther from the mouth
   (larger `dg`) to the nearer one.

The result is a spanning tree with exactly `n − 1` edges.  Its parent edges are
flow-path edges; the `n_roots − 1` root-join edges connect distinct flow paths
(they are necessarily *not* on-path edges).  The returned sparse matrix stores
both directions of every tree edge (2·(n−1) non-zeros) so it is a drop-in
replacement for the symmetric legacy matrix used by the ODE dispersal kernel.

Returns a NamedTuple with `distance_matrix`, `parent` (`0` = root), `roots`,
`tree_edges` (directed `(upstream, downstream, distance)` triples) and
`diagnostics`.
"""
function build_on_path_distance_matrix(distance_file::String, sites::Vector{String},
        site_to_river_distance::Dict{T3, Float64},
        site_to_elevation::Dict{T4, Float64};
        rtol::Float64 = ON_PATH_DEFAULT_RTOL,
        atol::Float64 = ON_PATH_DEFAULT_ATOL,
        store_both_directions::Bool = true) where {T3 <: AbstractString} where {T4 <: AbstractString}
    println("Building on-path (river-network) distance matrix from: $distance_file")
    println("  on-path tolerance: rtol=$rtol, atol=$(atol) m")

    site_to_idx = Dict(s => i for (i, s) in enumerate(sites))
    n_sites = length(sites)
    dg = Float64[get(site_to_river_distance, s, Inf) for s in sites]
    dg_site = Dict(sites[i] => dg[i] for i in 1:n_sites)
    elevation = Float64[get(site_to_elevation, s, NaN) for s in sites]

    # --- streaming pass: network distance of every model-site pair ----------
    pair_distance = Dict{Tuple{Int, Int}, Float64}()
    total_rows = 0
    reader = CSV.File(distance_file; delim=';')
    for row in reader
        total_rows += 1
        i = get(site_to_idx, row.ID_ORIGIN, 0)
        k = get(site_to_idx, row.ID_DESTINATION, 0)
        (i == 0 || k == 0 || i == k) && continue
        dist = Float64(row.RETICULAR_DIST)
        (isfinite(dist) && dist > 0) || continue
        a, b = i < k ? (i, k) : (k, i)
        key = (a, b)
        # The matrix is symmetric but read twice; keeping the smaller value
        # makes the stored distance independent of row order.
        pair_distance[key] = haskey(pair_distance, key) ?
            min(pair_distance[key], dist) : dist
    end
    println("  read $total_rows distance records; kept $(length(pair_distance)) model-site pairs")

    # --- on-path parent rule ------------------------------------------------
    parent = zeros(Int, n_sites)
    parent_distance = fill(Inf, n_sites)
    for ((a, b), dist) in pair_distance
        if dg[b] < dg[a]
            _consider_on_path_parent!(parent, parent_distance, a, b, dist, dg, sites, rtol, atol)
        elseif dg[a] < dg[b]
            _consider_on_path_parent!(parent, parent_distance, b, a, dist, dg, sites, rtol, atol)
        end
    end

    parent_edges = Tuple{Int, Int, Float64}[]
    parent_length = 0.0
    for i in 1:n_sites
        if parent[i] != 0
            push!(parent_edges, (i, parent[i], parent_distance[i]))
            parent_length += parent_distance[i]
        end
    end
    roots = Int[i for i in 1:n_sites if parent[i] == 0]
    println("  $(length(parent_edges)) parent edges, $(length(roots)) roots")

    # --- join the roots with an MST on network distance ---------------------
    root_set = Set(roots)
    root_pos = Dict(r => p for (p, r) in enumerate(roots))
    mst_candidates = Tuple{Float64, Int, Int}[]
    for ((a, b), dist) in pair_distance
        (a in root_set && b in root_set) && push!(mst_candidates, (dist, a, b))
    end
    sort!(mst_candidates, by = e -> (e[1],
        min(sites[e[2]], sites[e[3]]), max(sites[e[2]], sites[e[3]])))

    uf = collect(1:length(roots))
    function _find(x)
        while uf[x] != x
            uf[x] = uf[uf[x]]
            x = uf[x]
        end
        return x
    end

    root_edges = Tuple{Int, Int, Float64}[]
    root_length = 0.0
    for (dist, a, b) in mst_candidates
        x, y = _find(root_pos[a]), _find(root_pos[b])
        x == y && continue
        uf[x] = y
        root_length += dist
        # Orient upstream -> downstream: larger distance-to-mouth first.
        if dg[a] > dg[b] || (dg[a] == dg[b] && sites[a] < sites[b])
            push!(root_edges, (a, b, dist))
        else
            push!(root_edges, (b, a, dist))
        end
    end
    println("  $(length(root_edges)) root-join MST edges " *
            "($(round(root_length / 1000, digits = 1)) km)")
    length(root_edges) == length(roots) - 1 ||
        error("root MST is disconnected: $(length(root_edges)) edges for $(length(roots)) roots")

    # --- combine and build the sparse matrix --------------------------------
    tree_edges = vcat(parent_edges, root_edges)
    isempty(tree_edges) && n_sites > 1 &&
        error("on-path graph has no edges for $(n_sites) sites")

    I = Int[]
    J = Int[]
    V = Float64[]
    for (u, v, w) in tree_edges
        push!(I, v); push!(J, u); push!(V, w)   # v (downstream) <- u (upstream)
        if store_both_directions
            push!(I, u); push!(J, v); push!(V, w)
        end
    end
    distance_matrix = sparse(I, J, V, n_sites, n_sites)

    parent_inversions = Tuple{String, String, Float64, Float64}[]
    for (u, v, _) in parent_edges
        if isfinite(elevation[u]) && isfinite(elevation[v]) && elevation[u] < elevation[v]
            push!(parent_inversions, (sites[u], sites[v], elevation[u], elevation[v]))
        end
    end
    root_inversions = Tuple{String, String, Float64, Float64}[]
    for (u, v, _) in root_edges
        if isfinite(elevation[u]) && isfinite(elevation[v]) && elevation[u] < elevation[v]
            push!(root_inversions, (sites[u], sites[v], elevation[u], elevation[v]))
        end
    end

    diagnostics = (
        method = :on_path,
        rtol = rtol,
        atol = atol,
        n_sites = n_sites,
        n_edges = length(tree_edges),
        n_parent_edges = length(parent_edges),
        n_root_join_edges = length(root_edges),
        n_roots = length(roots),
        nnz = nnz(distance_matrix),
        parent_length_km = parent_length / 1000.0,
        root_join_length_km = root_length / 1000.0,
        total_length_km = (parent_length + root_length) / 1000.0,
        parent_elevation_inversions = parent_inversions,
        root_join_elevation_inversions = root_inversions
    )
    println("  tree: $(diagnostics.n_edges) edges, total length " *
            "$(round(diagnostics.total_length_km, digits = 1)) km, " *
            "$(length(parent_inversions)) parent elevation inversion(s)")

    return (
        distance_matrix = distance_matrix,
        parent = parent,
        roots = roots,
        tree_edges = tree_edges,
        diagnostics = diagnostics
    )
end

"""
    build_elevation_vector(site_df::DataFrame, sites::Vector{String})

Build a vector of elevations for each site.
"""
function build_elevation_vector(site_df::DataFrame, sites::Vector{String})
    site_to_elevation = Dict(row.CODIGO => row.ALTITUD for row in eachrow(site_df))

    elevations = Float64[]
    for site in sites
        if haskey(site_to_elevation, site)
            push!(elevations, site_to_elevation[site])
        else
            push!(elevations, 500.0)  # Default elevation
        end
    end

    return elevations
end

"""
    build_dam_passability_matrix(site_df::DataFrame, sites::Vector{String}, distances, elevations)

Build a dam passability matrix based on dam distances, directional river flow,
and the physical location of dams relative to site pairs.

A dam reduces passability (to 0.1) for a connection j → i only when the dam
lies between the two sites along the river path. This is determined by:
- Direction: whether i is downstream (e_i < e_j) or upstream (e_i > e_j) of j
- Distance: whether the dam-to-site distance is less than the total river distance between sites
"""
function build_dam_passability_matrix(site_df::DataFrame, sites::Vector{String}, distances, elevations)
    n_sites = length(sites)

    site_dam_info = Dict{String, NamedTuple{(:dist_upstream, :dist_downstream), Tuple{Float64, Float64}}}()

    for row in eachrow(site_df)
        codigo = row.CODIGO
        dist_up = row.Demb_arr_m
        dist_down = row.Demb_ab_m

        if dist_up == 0
            dist_up = 1000_000.0
        end
        if dist_down == 0
            dist_down = 1000_000.0
        end

        site_dam_info[codigo] = (dist_upstream=dist_up, dist_downstream=dist_down)
    end

    dams = ones(n_sites, n_sites)

    for j in 1:n_sites
        origin = sites[j]
        if !haskey(site_dam_info, origin)
            continue
        end
        origin_info = site_dam_info[origin]
        e_j = elevations[j]

        for i in 1:n_sites
            if i == j
                continue
            end
            dest = sites[i]
            if !haskey(site_dam_info, dest)
                continue
            end

            d_ij = distances[i, j]
            if d_ij == 0 || isinf(d_ij)
                continue
            end

            dest_info = site_dam_info[dest]
            e_i = elevations[i]

            if e_i <= e_j
                # downstream or flat: origin's downstream dam or dest's upstream dam between them
                if (origin_info.dist_downstream < 1000_000 && origin_info.dist_downstream < d_ij) ||
                   (dest_info.dist_upstream < 1000_000 && dest_info.dist_upstream < d_ij)
                    dams[i, j] = 0.1
                end
            end
            if e_i >= e_j
                # upstream or flat: origin's upstream dam or dest's downstream dam between them
                if (origin_info.dist_upstream < 1000_000 && origin_info.dist_upstream < d_ij) ||
                   (dest_info.dist_downstream < 1000_000 && dest_info.dist_downstream < d_ij)
                    dams[i, j] = 0.1
                end
            end
        end
    end

    return dams
end

# C5/E1: the model's site temperature level is the corrected per-site
# water-temperature baseline (1986-2005) produced by the `guadex_tw` pipeline.
# It is selected by its EXACT column name, never by a "TEMP" substring match, so
# the legacy sub-catchment air-temperature climatology (`TEMP_MEDIA_SC`) cannot
# be picked up accidentally.
const DEFAULT_WATER_TEMPERATURE_FILE =
    "guadex_tw/outputs/tables/water_temp_baseline_guadex_sites.csv"
const WATER_TEMPERATURE_COLUMN = "tw_baseline_mean"

"""
    load_site_water_temperature_baseline(path; value_col=WATER_TEMPERATURE_COLUMN)

Load the corrected per-site water-temperature baseline (1986-2005) that defines
the model's site temperature level (C5/E1).  The table must contain a `site_id`
column and the exact `value_col` column (`tw_baseline_mean`); the value is joined
to model sites by site code elsewhere.

Returns a `Dict{String,Float64}`.  Fails loudly on a missing file, a missing
column, a duplicate site or a non-finite value, so a mis-selected column (e.g.
the legacy air-temperature `TEMP_MEDIA_SC`) or a truncated file cannot silently
produce a wrong temperature.
"""
function load_site_water_temperature_baseline(path::AbstractString;
        value_col::AbstractString=WATER_TEMPERATURE_COLUMN)
    isfile(path) || error(
        "water-temperature baseline file not found: $path. The model site " *
        "temperature level is the corrected per-site water-temperature baseline " *
        "(column '$WATER_TEMPERATURE_COLUMN'). Pass the correct " *
        "`water_temperature_file`, or opt out explicitly with " *
        "`water_temperature_file=nothing, require_water_temperature=false` to use " *
        "the legacy air-temperature source.")
    df = CSV.read(path, DataFrame)
    hasproperty(df, :site_id) || error(
        "water-temperature baseline $path has no `site_id` column " *
        "(found: $(join(names(df), ", ")))")
    hasproperty(df, Symbol(value_col)) || error(
        "water-temperature baseline $path has no '$value_col' column " *
        "(found: $(join(names(df), ", "))). Refusing to substitute another " *
        "temperature column such as TEMP_MEDIA_SC.")

    lookup = Dict{String,Float64}()
    for row in eachrow(df)
        site = string(row.site_id)
        haskey(lookup, site) &&
            error("water-temperature baseline $path contains duplicate site '$site'")
        raw = row[Symbol(value_col)]
        raw === missing &&
            error("water-temperature baseline $path has a missing '$value_col' for site '$site'")
        value = Float64(raw)
        isfinite(value) || error(
            "water-temperature baseline $path has a non-finite '$value_col' for site '$site'")
        lookup[site] = value
    end
    return lookup
end

"""
    _legacy_site_temperatures(site_df, sites, legacy_temperature_column=nothing)

Explicit legacy temperature source: the named column of `site_df` when
`legacy_temperature_column` is given (e.g. `"TEMP_MEDIA_SC"`), otherwise the
elevation-based estimate (≈6.5 °C/km lapse rate).  Reached only via an explicit
opt-out from the corrected water-temperature baseline.
"""
function _legacy_site_temperatures(site_df::DataFrame, sites::Vector{String},
        legacy_temperature_column::Union{Nothing,AbstractString}=nothing)
    if legacy_temperature_column !== nothing
        hasproperty(site_df, Symbol(legacy_temperature_column)) || error(
            "legacy temperature column '$legacy_temperature_column' not found in " *
            "the site data (found: $(join(names(site_df), ", ")))")
        col = Symbol(legacy_temperature_column)
        site_to_temp = Dict(string(row.CODIGO) => Float64(row[col]) for row in eachrow(site_df))
        return [get(site_to_temp, string(s), 15.0) for s in sites]
    end
    # Rough approximation: temperature decreases ~6.5°C per 1000m
    elevations = build_elevation_vector(site_df, sites)
    return 20.0 .- (elevations ./ 1000.0 .* 6.5)
end

"""
    extract_site_temperatures(site_df::DataFrame, sites::Vector{String}; kwargs...)

Return the site temperature level used by the model (C5/E1).

By default the level is the corrected per-site water-temperature baseline
(1986-2005) read from `water_temperature_file` and selected by the exact column
`tw_baseline_mean`, joined to `sites` by site code.  The baseline file contains
776 sites while the model has 775 (`1.30.20` is excluded from the model), so
extra rows are ignored; a model site missing from the baseline is an error by
default rather than a silent fallback.

# Keyword arguments
- `water_temperature_file`: corrected baseline CSV
  (default [`DEFAULT_WATER_TEMPERATURE_FILE`](@ref)).  `nothing` disables it.
- `require_water_temperature`: when `true` (default) the file, the exact column
  and every model site are mandatory; a missing site raises an error.  Set to
  `false` only for synthetic use, which falls back per missing site to the
  explicit legacy source.
- `legacy_temperature_column`: exact `site_df` column to use as the legacy
  source (e.g. `"TEMP_MEDIA_SC"`); `nothing` uses the elevation estimate.
- `value_col`: baseline column name; defaults to `tw_baseline_mean`.

Legacy behaviour is therefore available only on an explicit opt-out, and
`TEMP_MEDIA_SC` is never selected implicitly.
"""
function extract_site_temperatures(site_df::DataFrame, sites::Vector{String};
        water_temperature_file::Union{Nothing,AbstractString}=DEFAULT_WATER_TEMPERATURE_FILE,
        require_water_temperature::Bool=true,
        legacy_temperature_column::Union{Nothing,AbstractString}=nothing,
        value_col::AbstractString=WATER_TEMPERATURE_COLUMN)

    if water_temperature_file === nothing
        require_water_temperature && error(
            "no water-temperature source configured. Pass `water_temperature_file` " *
            "(default '$DEFAULT_WATER_TEMPERATURE_FILE'), or opt out explicitly with " *
            "`water_temperature_file=nothing, require_water_temperature=false` to use " *
            "the legacy temperature source.")
        println("Site temperature level: corrected water-temperature baseline disabled; " *
                "using explicit legacy source " *
                (legacy_temperature_column === nothing ? "(elevation estimate)" :
                 "(column '$legacy_temperature_column')"))
        return _legacy_site_temperatures(site_df, sites, legacy_temperature_column)
    end

    lookup = load_site_water_temperature_baseline(water_temperature_file; value_col=value_col)
    temperatures = Vector{Float64}(undef, length(sites))
    missing_sites = String[]
    for (k, site) in enumerate(sites)
        key = string(site)
        if haskey(lookup, key)
            temperatures[k] = lookup[key]
        else
            push!(missing_sites, key)
            temperatures[k] = NaN
        end
    end

    if !isempty(missing_sites)
        require_water_temperature && error(
            "water-temperature baseline $water_temperature_file does not cover " *
            "$(length(missing_sites)) model site(s): $(join(missing_sites, ", ")). " *
            "A missing site must not silently receive another temperature; pass " *
            "`require_water_temperature=false` only to allow an explicit legacy fallback.")
        fallback = _legacy_site_temperatures(site_df, sites, legacy_temperature_column)
        for k in eachindex(temperatures)
            isnan(temperatures[k]) && (temperatures[k] = fallback[k])
        end
        @warn "water-temperature baseline is missing model sites; used the explicit " *
              "legacy fallback for: $(join(missing_sites, ", "))"
    end

    return temperatures
end

"""
    extract_habitat_suitability(site_df::DataFrame, sites::Vector{String})

Extract habitat suitability index for each site.
Uses IET (Índice de Estado Trófico) as habitat quality indicator.
Return a suitability score where higher IET = lower suitability (normalized).
"""
function extract_habitat_suitability(site_df::DataFrame, sites::Vector{String})
    # Try to find habitat quality column
    # IET is a good indicator (lower is better - oligotrophic). Think of it as a "health check" that tells you how much organic matter (mostly algae) is growing in the water.
    # We'll use a simple transformation: higher IET = lower suitability

    site_to_iet = Dict(row.CODIGO => row.IET for row in eachrow(site_df))

    # Transform IET to suitability (simple inverse, normalized)
    iet_values = collect(values(site_to_iet))
    iet_min, iet_max = minimum(iet_values), maximum(iet_values)

    suitability = Float64[]
    for site in sites
        if haskey(site_to_iet, site)
            iet = site_to_iet[site]
            # Normalize: lower IET = higher suitability
            suit = 1.0 - (iet - iet_min) / (iet_max - iet_min + 1e-6)
            push!(suitability, max(0.1, suit))  # Minimum suitability of 0.1
        else
            push!(suitability, 0.5)  # Default
        end
    end

    return suitability
end

# =============================================================================
# Explicit, opt-in biological-assumption options (E13-E16, E18).
#
# Every option below defaults to the behaviour the model had before the option
# existed, so the legacy outputs and the test suite are reproduced exactly.  The
# options are read from the `[biological_options]` section of the parameter file
# (see `parameters.jl`) and exposed through [`prepare_ode_data`](@ref).
# =============================================================================

# E16: species documented with `REPRODU_WITHIN_THE_BASIN = no` yet given a
# local literature growth rate.  With `nonreproducing_local_growth = :zero`
# their local r is set to zero; recruitment is to be represented later as
# estuarine immigration.
const NONREPRODUCING_SPECIES = Set(["AA", "LR", "MC"])

# E16: the Guadiato cluster (subcatchment 14) holds the four documented eel
# records from a fish farm that the density README says must not be used for the
# native/alien calculations.  `exclude_fishfarm_eel_records = true` zeroes their
# AA_DEN before any downstream builder reads the density table.
const FISHFARM_EEL_SUBCATCHMENT = "14.0"

# E15: capacity assigned to a `SIN_PECES = DIFICIL` site under
# `fishless_dificil_capacity = :near_zero`.  Matches the ODE's own internal
# capacity floor (`max(K, 1e-6)`), so the site behaves as a transit-only node:
# the logistic term saturates immediately and no local population can establish.
const NEAR_ZERO_CARRYING_CAPACITY = 1.0e-6

# E18: minimum site conductivity (uS/cm) for a strictly brackish species under
# `salinity_envelope = true`.  1000 uS/cm is the conventional freshwater/brackish
# boundary; the site `CONDUCTIVIDAD` column is used as the salinity proxy.
const SALINITY_BRACKISH_MIN_CONDUCTIVITY_US_CM = 1000.0

"""
    build_intrinsic_growth_rates(density_df::DataFrame, species_codes::Vector{String}, sites::Vector{String}, species_chars_df::DataFrame;
        seasonal_rates::AbstractDict=Dict{String,Float64}(),
        absence_growth_fraction::Real=0.1,
        nonreproducing_local_growth::Symbol=:legacy)

Build intrinsic growth rates matrix from density data.
Uses literature-derived daily intrinsic growth rates for each species.

# Rate conversion (minor #4)
The default keeps the historical, tested behaviour: **every** species' annual
rate is converted as `r_daily = r_annual / 365`; absent species receive
`absence_growth_fraction · r_daily`.  The previously documented seasonal
conversion `r_daily = log1p(r_seasonal) / 183` is now offered as an *explicit,
opt-in* option: pass a species→rate mapping as `seasonal_rates` (e.g.
`Dict("GH" => 4.0)`) and the listed species use the seasonal formula while all
others keep `r_annual / 365`.  The default (`seasonal_rates` empty) is
unchanged, so existing runs and results are reproduced exactly.

# E13 — absence growth fraction
`absence_growth_fraction` is the fraction of the local daily rate given to a
species observed absent at a site (default `0.1`, the historical value).  Setting
it to `1.0` removes the presence-based penalty entirely, so temperature and
habitat alone decide where a species can establish (the review's recommendation);
the default is deliberately **not** changed.  Set it to `0.0` to forbid local
growth at absent sites.

# E16 — non-reproducing species
`nonreproducing_local_growth = :zero` sets the local `r` of the species with
`REPRODU_WITHIN_THE_BASIN = no` (`AA`, `LR`, `MC`) to zero; the default
(`:legacy`) keeps their literature rate.  This is opt-in and reversible.

# Literature Sources
- References are based on empirical studies from the Guadalquivir River Basin

# Species and Citations
Native species:
- AB (Aphanius baeticus): r_annual ≈ 1.0-1.5 /year, seasonal r ≈ 0.03-0.05 (Ref 7:成熟&葡萄)
- AH (Anaecypris hispanica): r_annual ≈ 0.8-1.2 /year (Ref 11: PMC)
- SP (Squalius pyrenaicus): r_annual ≈ 0.4-0.7 /year (Ref 11: PMC)
- PW (Pseudochondrostoma willkommii): r_annual ≈ 0.2-0.4 /year (Ref 6:研究)
- LS (Luciobarbus sclateri): r_annual ≈ 0.15-0.25 /year (Ref 6:研究)
- SA (Squalius alburnoides): r_annual ≈ 0.5-0.9 /year (Ref 11: PMC)
- IL (Iberochondrostoma lemmingii): r_annual ≈ 0.4-0.7 /year (Ref 11: PMC)
- CP (Cobitis paludica): r_annual ≈ 0.5-0.8 /year (Ref 19)
- IO (Iberochondrostoma oretanum): r_annual ≈ 0.3-0.6 /year (Ref 5)

Invasive species:
- GH (Gambusia holbrooki): seasonal r ≈ 0.029/day, annual r can exceed 4.0 /year (Ref 4,7)
- MS (Micropterus salmoides): r_annual ≈ 0.3-0.5 /year (Ref 2)
- LG (Lepomis gibbosus): r_annual ≈ 0.4-0.6 /year (Ref 23)
- CC (Cyprinus carpio): r_annual ≈ 0.3-0.5 /year (Ref 1)
- CG (Carassius gibelio): r_annual ≈ 0.5-1.1 /year (Ref 18)
- AM (Ameiurus melas): r_annual ≈ 0.25-0.45 /year (Ref 1)
- OM (Oncorhynchus mykiss): r_annual ≈ 0.3-0.5 /year (Ref 1)
- EL (Esox lucius): r_annual ≈ 0.2-0.4 /year (Ref 1)
- GL (Gobio lozanoi): r_annual ≈ 0.6-1.0 /year (Ref 18)
- TT (Tinca tinca): r_annual ≈ 0.2-0.4 /year (Ref 1)

Other species:
- AA (Anguilla anguilla): r_annual ≈ 0.05-0.15 /year (Ref 9)
- MC (Mugil cephalus): r_annual ≈ 0.2-0.4 /year (Ref 18)
- LR (Liza ramada): r_annual ≈ 0.2-0.4 /year (Ref 18)
- ST (Salmo trutta): r_annual ≈ 0.3-0.5 /year (typical for salmonids)
"""
function build_intrinsic_growth_rates(density_df::DataFrame, species_codes::Vector{String},
                                       sites::Vector{String}, species_chars_df::DataFrame;
                                       seasonal_rates::AbstractDict=Dict{String,Float64}(),
                                       absence_growth_fraction::Real=0.1,
                                       nonreproducing_local_growth::Symbol=:legacy)

    absence_growth_fraction >= 0 ||
        error("absence_growth_fraction must be non-negative (got $absence_growth_fraction)")
    nonreproducing_local_growth in (:legacy, :zero) ||
        error("nonreproducing_local_growth must be :legacy or :zero (got :$(nonreproducing_local_growth))")

    n_sites = length(sites)
    n_species = length(species_codes)

    # Create site to density row mapping
    site_to_row = Dict(row.CODIGO => rownum for (rownum, row) in enumerate(eachrow(density_df)))

    # Literature-derived annual intrinsic growth rates (r) for each species
    # These are the maximum per capita rates of increase under ideal conditions
    # Converted to daily rates: r_daily = r_annual / 365
    # Citations correspond to references in docs/Fish Growth and Dispersal Data Request.md
    annual_growth_rates = Dict{String, Float64}(
        # Native Endemics
        "AB" => 1.2,   # Aphanius baeticus - high growth, short lifespan (Ref 7)
        "AH" => 1.0,   # Anaecypris hispanica - high growth, short lifespan (Ref 11)
        "SP" => 0.55,  # Squalius pyrenaicus - medium growth (Ref 11)
        "PW" => 0.3,   # Pseudochondrostoma willkommii - medium growth (Ref 6)
        "LS" => 0.2,   # Luciobarbus sclateri - low growth, late maturity (Ref 6)
        "SA" => 0.7,   # Squalius alburnoides - medium-high growth (Ref 11)
        "IL" => 0.55,  # Iberochondrostoma lemmingii - medium growth (Ref 11)
        "CP" => 0.65,  # Cobitis paludica - medium growth (Ref 19)
        "IO" => 0.45,  # Iberochondrostoma oretanum - medium growth (Ref 5)

        # Invasive Species
        "GH" => 4.0,   # Gambusia holbrooki - extremely high growth (Ref 4,7)
        "MS" => 0.4,    # Micropterus salmoides - medium growth (Ref 2)
        "LG" => 0.5,   # Lepomis gibbosus - medium growth (Ref 23)
        "CC" => 0.4,   # Cyprinus carpio - medium growth (Ref 1)
        "CG" => 0.8,   # Carassius gibelio - medium-high growth (Ref 18)
        "AM" => 0.35,  # Ameiurus melas - medium-low growth (Ref 1)
        "OM" => 0.4,   # Oncorhynchus mykiss - medium growth (Ref 1)
        "EL" => 0.3,   # Esox lucius - medium-low growth (Ref 1)
        "GL" => 0.8,   # Gobio lozanoi - medium-high growth (Ref 18)
        "TT" => 0.3,   # Tinca tinca - medium-low growth (Ref 1)

        # Diadromous/Marine Species
        "AA" => 0.1,   # Anguilla anguilla - very low growth, long lifespan (Ref 9)
        "MC" => 0.3,   # Mugil cephalus - medium growth (Ref 18)
        "LR" => 0.3,   # Liza ramada - medium growth (Ref 18)

        # Salmonids
        "ST" => 0.4,   # Salmo trutta - medium growth (typical for salmonids)
    )

    growth_rates = zeros(n_sites, n_species)

    for (s_idx, sp_code) in enumerate(species_codes)
        # Get annual growth rate from literature (use midpoint of ranges)
        r_annual = get(annual_growth_rates, sp_code, 0.5)  # default 0.5 if unknown

        # Convert to a daily instantaneous rate.  Default (minor #4): the
        # historical `r_annual / 365` for every species.  Opt-in seasonal
        # conversion for species explicitly listed in `seasonal_rates`:
        # `log1p(r_seasonal) / 183`, matching the documented Gambusia formula.
        r_daily = if haskey(seasonal_rates, sp_code)
            log1p(Float64(seasonal_rates[sp_code])) / 183.0
        else
            r_annual / 365.0
        end

        # E16: opt-in zeroing of the local rate for species that do not reproduce
        # within the basin.  Applied before the absence branch so both present and
        # absent cells are zero.
        if nonreproducing_local_growth == :zero && sp_code in NONREPRODUCING_SPECIES
            r_daily = 0.0
        end

        for (site_idx, site) in enumerate(sites)
            if haskey(site_to_row, site)
                row_idx = site_to_row[site]
                density_col = Symbol("$(sp_code)_DEN")

                if hasproperty(density_df, density_col)
                    density = density_df[row_idx, density_col]

                    # If species is present (density > 0), use full growth rate
                    # If absent, use reduced rate (potential colonization from nearby)
                    if density > 0
                        growth_rates[site_idx, s_idx] = r_daily
                    else
                        # E13: fraction of the local rate at an observed absence
                        # (legacy default 0.1).
                        growth_rates[site_idx, s_idx] = r_daily * Float64(absence_growth_fraction)
                    end
                end
            end
        end
    end

    return growth_rates
end

"""
    build_carrying_capacity(density_df::DataFrame, site_df::DataFrame, sites::Vector{String}, species_codes::Vector{String}; scaling=10.0)

Build site-specific carrying capacities from observed fish density data.

The carrying capacity K_i for each site is derived from the observed total fish density,
scaled by `scaling` to represent the maximum sustainable biomass the site can support.
The scaling factor accounts for:
- Natural fluctuations around observed densities
- Additional habitat not sampled during surveys
- Density-dependent regulation allowing populations to exceed observed levels

`scaling = 1.0` sets K_i to the observed total density (with the floor below), i.e.
the observed snapshot is treated as the equilibrium.  `scaling = 10.0` is the legacy
convention.  A minimum-capacity floor of 5.0 observed-density units (10th percentile
of non-zero observations when larger) is multiplied by the same `scaling`, so the
effective multiplier against observed density is `scaling` for every site.

# E14 — isolated-pool capacity
`pool_capacity_mode` controls how sites flagged `EN_POZAS = Si` (isolated pools)
enter the capacity estimate.  `:legacy` (default) uses their raw observed density;
`:cap` clips a pool site's capacity at the **non-pool median**; `:exclude`
replaces it with that median, so a single-visit pool sample cannot set the
capacity at all.  The 84 pool sites hold ~56% of total density, so this choice
strongly affects the reported high-K tail.

# E15 — fishless DIFICIL capacity
`fishless_dificil_capacity = :near_zero` sets K to
`NEAR_ZERO_CARRYING_CAPACITY` at sites with `SIN_PECES = DIFICIL` (the field team
judged them unable to hold fish), making them transit-only nodes.  The default
`:legacy` keeps their observed capacity.  This interacts with E13: raising
`absence_growth_fraction` while zeroing K at these sites pulls establishment in
opposite directions.

# Arguments
- `density_df`: DataFrame with species density data (from load_species_density_data)
- `site_df`: DataFrame with site data (must carry `EN_POZAS` / `SIN_PECES` for the
  corresponding non-legacy options)
- `sites`: Vector of site codes in order
- `species_codes`: Vector of species codes
- `scaling`: factor applied to the observed total density (default 10.0)
- `pool_capacity_mode`: `:legacy` (default), `:exclude` or `:cap`
- `fishless_dificil_capacity`: `:legacy` (default) or `:near_zero`

# Returns
- Vector of carrying capacities for each site
"""
function build_carrying_capacity(density_df::DataFrame, site_df::DataFrame, sites::Vector{String}, species_codes::Vector{String}; scaling::Real=10.0,
        pool_capacity_mode::Symbol=:legacy, fishless_dificil_capacity::Symbol=:legacy)
    println("Building site-specific carrying capacities from density data...")

    pool_capacity_mode in (:legacy, :exclude, :cap) ||
        error("pool_capacity_mode must be :legacy, :exclude or :cap (got :$(pool_capacity_mode))")
    fishless_dificil_capacity in (:legacy, :near_zero) ||
        error("fishless_dificil_capacity must be :legacy or :near_zero (got :$(fishless_dificil_capacity))")

    site_to_idx = Dict{String, Int}()
    for (rownum, row) in enumerate(eachrow(density_df))
        site_to_idx[row.CODIGO] = rownum
    end

    density_cols = [Symbol("$(sp)_DEN") for sp in species_codes]

    K_scaling = Float64(scaling)
    K_scaling > 0 || error("carrying-capacity scaling must be > 0")

    # Site-attribute lookups for the E14/E15 options.  A missing column is only an
    # error when the corresponding option is actually requested, so the legacy
    # call works with a bare `site_df`.
    if pool_capacity_mode != :legacy
        hasproperty(site_df, :EN_POZAS) ||
            error("pool_capacity_mode=:$(pool_capacity_mode) requires an EN_POZAS column in site_df")
    end
    if fishless_dificil_capacity != :legacy
        hasproperty(site_df, :SIN_PECES) ||
            error("fishless_dificil_capacity=:$(fishless_dificil_capacity) requires a SIN_PECES column in site_df")
    end
    row_lookup = Dict{String, Int}()
    for (rownum, row) in enumerate(eachrow(site_df))
        row_lookup[string(row.CODIGO)] = rownum
    end
    function site_flag(site, column)
        hasproperty(site_df, column) || return ""
        idx = get(row_lookup, string(site), 0)
        idx == 0 && return ""
        return _normalise_obstacle_status(site_df[idx, column])
    end
    is_pool = pool_capacity_mode == :legacy ? falses(length(sites)) :
        [site_flag(site, :EN_POZAS) == "SI" for site in sites]
    is_dificil = fishless_dificil_capacity == :legacy ? falses(length(sites)) :
        [site_flag(site, :SIN_PECES) == "DIFICIL" for site in sites]

    raw_capacities = Float64[]
    for site in sites
        if haskey(site_to_idx, site)
            row_idx = site_to_idx[site]
            row = density_df[row_idx, :]
            total_density = 0.0
            for col in density_cols
                if hasproperty(row, col)
                    val = row[col]
                    if !ismissing(val) && !isnan(val) && val > 0
                        total_density += val
                    end
                end
            end
            push!(raw_capacities, total_density * K_scaling)
        else
            push!(raw_capacities, 0.0)
        end
    end

    # The minimum-capacity floor is expressed in *observed-density* units and then
    # multiplied by the same factor, so the effective multiplier against observed
    # density is `K_scaling` for every site.  The legacy absolute floor of 50.0
    # was defined at the legacy 10x base, i.e. 5.0 observed-density units, which
    # reproduces the historical K floor exactly when `scaling == 10.0`.
    K_min_observed = 5.0
    nonzero_observed = filter(c -> c > 0, raw_capacities ./ K_scaling)
    floor_observed = isempty(nonzero_observed) ? K_min_observed :
        max(K_min_observed, quantile(nonzero_observed, 0.1))
    K_floor = K_scaling * floor_observed

    carrying_capacity = Float64[]
    for cap in raw_capacities
        push!(carrying_capacity, max(cap, K_floor))
    end

    # E14: post-process ONLY the isolated-pool sites against the non-pool median,
    # so every non-pool site keeps the legacy capacity exactly.
    if pool_capacity_mode != :legacy
        nonpool = [carrying_capacity[i] for i in eachindex(carrying_capacity) if !is_pool[i]]
        isempty(nonpool) &&
            error("pool_capacity_mode=:$(pool_capacity_mode) requires at least one non-pool site")
        nonpool_reference = median(nonpool)
        for i in eachindex(carrying_capacity)
            is_pool[i] || continue
            carrying_capacity[i] = pool_capacity_mode == :cap ?
                min(carrying_capacity[i], nonpool_reference) : nonpool_reference
        end
        println("Isolated-pool capacity mode :$(pool_capacity_mode): $(count(is_pool)) pool " *
                "site(s) adjusted to the non-pool median $(nonpool_reference)")
    end

    # E15: post-process the fishless DIFICIL sites after the floor so a
    # transit-only node is genuinely near-zero rather than raised to K_floor.
    if fishless_dificil_capacity == :near_zero
        for i in eachindex(carrying_capacity)
            is_dificil[i] && (carrying_capacity[i] = NEAR_ZERO_CARRYING_CAPACITY)
        end
        println("Fishless DIFICIL capacity mode :near_zero: $(count(is_dificil)) site(s) set to " *
                "K = $(NEAR_ZERO_CARRYING_CAPACITY)")
    end

    println("Carrying capacity range: $(minimum(carrying_capacity)) - $(maximum(carrying_capacity))")
    println("Mean carrying capacity: $(mean(carrying_capacity))")
    println("K floor: $K_floor ($(floor_observed) observed-density units x $(K_scaling))")

    return carrying_capacity
end

"""
    build_dispersal_scaling(species_codes::Vector{String})

Build species-specific dispersal scaling factors based on literature values.
These scaling factors convert the base dispersal matrix to species-specific rates.

# Literature Sources
Dispersal rates are reported as km/year in the literature and converted to relative scaling
factors (normalized so that the median = 1.0). The base dispersal matrix uses the
dispersal_intensity parameter, and these scaling factors adjust per species.

# Species Dispersal Rates (km/year) and References:
Native Species:
- AB (Aphanius baeticus): < 0.5 km/year - highly fragmented populations (Ref 3)
- AH (Anaecypris hispanica): 0.1-0.3 km/year - highly sedentary (Ref 11, 15, 16)
- SP (Squalius pyrenaicus): 1-5 km/year (Ref 11, 12)
- PW (Pseudochondrostoma willkommii): 15-40 km/year - migratory/potadromous (Ref 9, 10)
- LS (Luciobarbus sclateri): 10-30 km/year - migratory (Ref 6, 11)
- SA (Squalius alburnoides): 2-8 km/year (Ref 11)
- IL (Iberochondrostoma lemmingii): 0.5-2 km/year (Ref 11)
- CP (Cobitis paludica): 0.5-2 km/year (Ref 19)
- IO (Iberochondrostoma oretanum): < 1 km/year - fragmented (Ref 5)

Invasive Species:
- GH (Gambusia holbrooki): 8-42 km/year - highly dispersive (Ref 7, 10)
- MS (Micropterus salmoides): 10-25 km/year (Ref 2, 10)
- LG (Lepomis gibbosus): 8-42 km/year (Ref 5, 10)
- CC (Cyprinus carpio): 10-30 km/year (Ref 1)
- CG (Carassius gibelio): 5-15 km/year (Ref 18)
- AM (Ameiurus melas): 5-15 km/year (Ref 1)
- OM (Oncorhynchus mykiss): 5-20 km/year (Ref 1)
- EL (Esox lucius): 5-20 km/year (Ref 1)
- GL (Gobio lozanoi): 1-5 km/year (Ref 18)
- TT (Tinca tinca): 2-10 km/year (Ref 1)

Other Species:
- AA (Anguilla anguilla): < 10 km/year - dam restricted (Ref 9)
- MC (Mugil cephalus): 50-100 km/year - euryhaline, high mobility (Ref 18)
- LR (Liza ramada): 50-100 km/year - euryhaline, high mobility (Ref 18)
- ST (Salmo trutta): 5-20 km/year - typical for salmonids

# References (from docs/Fish Growth and Dispersal Data Request.md):
1. Freshwater Fish Biodiversity in a Large Mediterranean Basin (Guadalquivir River)
2. Conservation status of freshwater fish in the Guadalquivir River Basin
3. Persistence despite isolation: Temporal genomic structure in Aphanius baeticus
4. Spatio-temporal and transmission dynamics of sarcoptic mange (for Gambusia birth rates)
5. Threatened Freshwater Fishes of the Mediterranean Basin
6. Age, growth and reproduction of the barbel, Barbus sclateri
7. Age, growth and reproduction of Aphanius iberus in the lower Guadalquivir
9. Why and when do freshwater fish migrate? (Iberian Peninsula)
10. Forensic reconstruction of Ictalurus punctatus invasion routes
11. Broad-scale sampling of primary freshwater fish populations (PMC/PeerJ)
12. Broad-scale sampling of primary freshwater fish populations (PeerJ)
15. Microsatellite analysis of genetic population structure of Anaecypris hispanica
16. Spatial and temporal variation in population genetic structure of Nile tilapia
18. A Long-Term Spatiotemporal Analysis of the Fish Community
19. Study on Invasive Alien Species – Development of Risk Assessments
"""
function build_dispersal_scaling(species_codes::Vector{String})
    # Annual dispersal rates in km/year from literature
    annual_dispersal_rates = Dict{String, Float64}(
        # Native Endemics
        "AB" => 0.3,   # < 0.5 km/year - highly fragmented (Ref 3)
        "AH" => 0.2,   # 0.1-0.3 km/year - highly sedentary (Ref 11, 15)
        "SP" => 3.0,   # 1-5 km/year (Ref 11, 12)
        "PW" => 27.5,  # 15-40 km/year - migratory (Ref 9, 10)
        "LS" => 20.0,  # 10-30 km/year - migratory (Ref 6, 11)
        "SA" => 5.0,   # 2-8 km/year (Ref 11)
        "IL" => 1.25,  # 0.5-2 km/year (Ref 11)
        "CP" => 1.25,  # 0.5-2 km/year (Ref 19)
        "IO" => 0.5,   # < 1 km/year - fragmented (Ref 5)

        # Invasive Species
        "GH" => 25.0,  # 8-42 km/year - highly dispersive (Ref 7, 10)
        "MS" => 17.5,  # 10-25 km/year (Ref 2, 10)
        "LG" => 25.0,  # 8-42 km/year (Ref 5, 10)
        "CC" => 20.0,  # 10-30 km/year (Ref 1)
        "CG" => 10.0,  # 5-15 km/year (Ref 18)
        "AM" => 10.0,  # 5-15 km/year (Ref 1)
        "OM" => 12.5,  # 5-20 km/year (Ref 1)
        "EL" => 12.5,  # 5-20 km/year (Ref 1)
        "GL" => 3.0,   # 1-5 km/year (Ref 18)
        "TT" => 6.0,   # 2-10 km/year (Ref 1)

        # Diadromous/Marine Species
        "AA" => 5.0,   # < 10 km/year - dam restricted (Ref 9)
        "MC" => 75.0,  # 50-100 km/year - euryhaline (Ref 18)
        "LR" => 75.0,  # 50-100 km/year - euryhaline (Ref 18)

        # Salmonids
        "ST" => 12.5,  # 5-20 km/year - typical for salmonids
    )

    # Calculate scaling factors relative to median dispersal rate
    # All rates are normalized so that the median species has scale = 1.0
    rates = [get(annual_dispersal_rates, sp, 5.0) for sp in species_codes]  # default 5.0 km/year
    median_rate = median(rates)

    scaling = Float64[]
    for sp in species_codes
        rate = get(annual_dispersal_rates, sp, 5.0)
        push!(scaling, rate / median_rate)
    end

    return scaling
end

"""
    observed_density_matrix(density_df, sites, species)

Build the observed initial-state matrix (`n_sites × n_species`, in `sites`
order) from the observed density table, aligning rows to model sites **by site
code** (`CODIGO`) rather than by row position (minor #3).

This is deliberately robust to a re-ordered or shuffled density table and to a
density table that carries extra sites: an `innerjoin` on `CODIGO` with
`order=:left` reindexes the density rows into `sites` order.  It errors loudly
when a modelled site has no density row (rather than silently shifting states
between sites), when duplicate `CODIGO` rows are present, or when a required
`<sp>_DEN` column is missing.  `NaN` entries are replaced by `0.0` and negative
densities are clamped to zero.
"""
function observed_density_matrix(density_df::DataFrame, sites::AbstractVector,
        species::AbstractVector)
    density_cols = [Symbol("$(sp)_DEN") for sp in species]
    missing_cols = [c for c in density_cols if !hasproperty(density_df, c)]
    isempty(missing_cols) ||
        error("density table is missing required column(s): $(missing_cols)")

    keyed = DataFrame(CODIGO=String.(density_df.CODIGO))
    for col in density_cols
        values = Float64[]
        for v in density_df[!, col]
            push!(values, ismissing(v) ? NaN : Float64(v))
        end
        keyed[!, col] = values
    end

    code_counts = Dict{String,Int}()
    for code in keyed.CODIGO
        code_counts[code] = get(code_counts, code, 0) + 1
    end
    duplicated = [code for (code, n) in code_counts if n > 1]
    isempty(duplicated) ||
        error("density table has duplicate CODIGO row(s): $(sort(duplicated))")

    target = DataFrame(CODIGO=String.(sites))
    joined = innerjoin(target, keyed; on=:CODIGO, order=:left)
    nrow(joined) == length(sites) ||
        error("observed density table has no row for $(length(sites) - nrow(joined)) " *
              "model site(s): $(sort(setdiff(String.(sites), String.(joined.CODIGO))))")

    matrix = Matrix{Float64}(joined[:, density_cols])
    replace!(matrix, NaN => 0.0)
    return max.(matrix, 0.0)
end

"""
    drop_fishfarm_eel_records(density_df, site_df;
                              subcatchment=FISHFARM_EEL_SUBCATCHMENT)
        -> (density_df, dropped_sites)

E16: return a copy of `density_df` with `AA_DEN` zeroed at the sites in the
Guadiato fish-farm cluster (subcatchment `14.0`).  The density README documents
four eel records there as fish-farm escapes that must not enter the native/alien
calculations.  Records at every other site are untouched; the function is a no-op
(and returns an empty `dropped_sites`) when there is nothing to drop.  Opt-in via
`exclude_fishfarm_eel_records = true`; reversible because the input table is not
mutated.
"""
function drop_fishfarm_eel_records(density_df::DataFrame, site_df::DataFrame;
        subcatchment::AbstractString=FISHFARM_EEL_SUBCATCHMENT)
    hasproperty(density_df, :AA_DEN) ||
        error("exclude_fishfarm_eel_records requires an AA_DEN column in the density table")

    site_subcatchment = Dict{String, String}()
    if hasproperty(site_df, :CODIGO_S)
        for row in eachrow(site_df)
            site_subcatchment[string(row.CODIGO)] = string(row.CODIGO_S)
        end
    end

    out = copy(density_df)
    dropped = String[]
    for i in axes(out, 1)
        site = string(out.CODIGO[i])
        value = out[i, :AA_DEN]
        (value === missing || !(value isa Real) || !(Float64(value) > 0)) && continue
        get(site_subcatchment, site, "") == string(subcatchment) || continue
        out[i, :AA_DEN] = 0.0
        push!(dropped, site)
    end
    return (out, sort!(dropped))
end

"""
    apply_salinity_envelope(growth_rates, species_codes, species_chars_df, site_df,
                            sites; min_conductivity=SALINITY_BRACKISH_MIN_CONDUCTIVITY_US_CM)
        -> (growth_rates, n_filtered)

E18: multiply the local growth of **strictly brackish** species by zero at sites
whose `CONDUCTIVIDAD` (the salinity proxy) is below `min_conductivity`.  Species
whose trait is `brackish/freshwater` (euryhaline) and freshwater species are left
unchanged, because the observed data place them across the conductivity range.
Sites without a finite conductivity value are left unchanged ("where the data
support it").  The envelope is a local growth multiplier, mathematically
equivalent to a per-species habitat-suitability factor, so the ODE is unchanged.

Returns the filtered growth matrix (a copy) and the number of site × species cells
set to zero.
"""
function apply_salinity_envelope(growth_rates, species_codes, species_chars_df,
        site_df::DataFrame, sites::Vector{String};
        min_conductivity::Real=SALINITY_BRACKISH_MIN_CONDUCTIVITY_US_CM)
    salinity_trait = Dict{String, String}()
    if hasproperty(species_chars_df, :SP) && hasproperty(species_chars_df, :SALINITY)
        for row in eachrow(species_chars_df)
            salinity_trait[lowercase(string(row.SP))] = lowercase(string(row.SALINITY))
        end
    end
    function strictly_brackish(code)
        trait = get(salinity_trait, lowercase(string(code)), "")
        return occursin("brackish", trait) && !occursin("freshwater", trait)
    end

    conductivity = Dict{String, Float64}()
    if hasproperty(site_df, :CONDUCTIVIDAD)
        for row in eachrow(site_df)
            value = row.CONDUCTIVIDAD
            (value === missing || !(value isa Real) || isnan(Float64(value))) && continue
            conductivity[string(row.CODIGO)] = Float64(value)
        end
    end

    filtered = copy(growth_rates)
    n_filtered = 0
    for (s_idx, code) in enumerate(species_codes)
        strictly_brackish(code) || continue
        for (i, site) in enumerate(sites)
            cond = get(conductivity, string(site), NaN)
            if !isnan(cond) && cond < Float64(min_conductivity) && filtered[i, s_idx] != 0.0
                filtered[i, s_idx] = 0.0
                n_filtered += 1
            end
        end
    end
    return (filtered, n_filtered)
end

"""
    prepare_ode_data(;
        connectivity_file::String = "data/ConnectivityUTM.csv",
        density_file::String = "data/BIOTIC/FishDensity_and_Juveniles_Matrix.csv",
        species_chars_file::String = "data/ABIOTIC/caracteristicas_peces_Guadalquivir_03-04-2018.csv",
        environmental_file::String = "data/ABIOTIC/Matriz_Ambiental_Data.csv",
        distance_file::String = "data/Matrix_distances_1037puntos_BRUTO_FINAL.csv",
        interaction_file::String = "data/BIOTIC/Interacciones_peces_Guadalquivir_03-04-2018_ENG.csv",
        upstream_cost::Float64 = 0.01,
        cedex_var_file::Union{Nothing,String} = nothing,
        cedex_esc_uts_file::Union{Nothing,String} = nothing,
        obstacles_file::Union{Nothing,String} = nothing,
        obstacle_mode::Symbol = :legacy,
        obstacle_matching_tolerance::Float64 = 2000.0,
        obstacle_passability::Float64 = 0.1,
        obstacle_downstream_passability::Float64 = 0.5,
        connectivity_method::Symbol = :legacy,
        on_path_rtol::Float64 = 0.02,
        on_path_atol::Float64 = 200.0
    )

Prepare all data needed for the ODE metacommunity model.
## Arguments
- `connectivity_file`: Path to site connectivity data
- `density_file`: Path to species density data
- `species_chars_file`: Path to species characteristics data
- `environmental_file`: Path to environmental data
- `distance_file`: Path to distance matrix data
- `interaction_file`: Path to species interaction data
- `upstream_cost`: Additional cost factor for upstream dispersal
- `connectivity_method`: `:legacy` (default) keeps the original
  sub-catchment chaining in [`build_distance_matrix`](@ref); `:on_path` builds the
  dendritic tree from the network-distance matrix with
  [`build_on_path_distance_matrix`](@ref) (issue C4).
- `on_path_rtol`, `on_path_atol`: tolerance used by the `:on_path` rule
  `|D(i,k) − (d(i) − d(k))| ≤ rtol·D(i,k) + atol` (metres).
- `obstacles_file`: Optional obstacle inventory; loaded but ignored unless `obstacle_mode = :overlay`
- `obstacle_mode`: `:legacy` (default, no obstacle overlay) or `:overlay`
- `obstacle_matching_tolerance`: Matching tolerance in metres for the obstacle overlay
- `obstacle_passability`: Uniform passability applied to upstream movement on
  matched edges (the inventory's `IF` index is deliberately not interpreted)
- `obstacle_downstream_passability`: Uniform passability applied to downstream
  movement on matched edges.  Barriers accumulate by multiplication along a
  link, so two matched obstacles of 0.1 give 0.01 upstream.  With
  `connectivity_method = :on_path` obstacles are snapped only to the corrected
  tree edges; non-operational, out-of-basin, legacy-duplicate and (by default)
  ambiguous-status structures are excluded and annotated (counts in
  `obstacle_metadata`).
- `obstacle_include_ambiguous_status`: when `true` the structures with an
  unconfirmed `ESTADO` code (SC/DM/ND/OT/blank) are kept as barriers instead of
  being excluded; they remain flagged `status_class = "ambiguous"`.
- `water_temperature_file`: corrected per-site water-temperature baseline (1986-2005)
  used as the model's site temperature level (C5/E1); selected by the exact column
  `tw_baseline_mean`. Pass `nothing` to opt out.
- `require_water_temperature`: when `true` (default) the baseline file, its exact
  column and every model site are mandatory (missing data raises an error).
- `legacy_temperature_column`: explicit legacy `site_df` column (e.g.
  `"TEMP_MEDIA_SC"`) used only when the corrected baseline is disabled or (with
  `require_water_temperature=false`) incomplete; `nothing` uses the elevation estimate.

## Explicit biological-assumption options (E13-E16, E18)
All default to the legacy behaviour; none is applied unless requested.  See
`docs/climate_scenarios.md`.
- `absence_growth_fraction` (E13, default `0.1`): fraction of the local rate
  given to a species observed absent at a site.  `1.0` lets temperature/habitat
  decide establishment; `0.0` forbids local growth at absent sites.
- `pool_capacity_mode` (E14, default `:legacy`): `:exclude` / `:cap` adjust the
  carrying capacity of `EN_POZAS = Si` isolated-pool sites against the non-pool
  median (a post-processing step inside [`build_carrying_capacity`](@ref)).
- `fishless_dificil_capacity` (E15, default `:legacy`): `:near_zero` sets K to
  `NEAR_ZERO_CARRYING_CAPACITY` at `SIN_PECES = DIFICIL` sites.
- `nonreproducing_local_growth` (E16, default `:legacy`): `:zero` sets the local
  `r` of `AA`, `LR`, `MC` to zero.
- `exclude_fishfarm_eel_records` (E16, default `false`): drop the four documented
  Guadiato fish-farm eel records before any builder reads the density table.
- `salinity_envelope` (E18, default `false`): when `true`, zero the local growth
  of strictly brackish species at sites below
  `SALINITY_BRACKISH_MIN_CONDUCTIVITY_US_CM`.  No elevation envelope is provided,
  deliberately: a static elevation envelope would block upslope range shifts.

## Returns a NamedTuple with:
- params: MetacommunityParams struct
- sites: Vector of site codes
- species: Vector of species codes
- distance_matrix: Sparse distance matrix
- elevations: Vector of elevations
- dams: Effective dam/obstacle passability matrix
- legacy_dams: Passability matrix from the legacy connectivity fields
- obstacle_overlay: Passability matrix derived from the optional obstacle inventory
- obstacle_mapping_diagnostics: Per-obstacle network matching diagnostics
  (including `status_class`, `outcome`, `duplicate_of_legacy_dam`)
- obstacle_metadata: E12 classification/accumulation counts (matched, applied,
  per-category exclusions, ambiguous/non-operational status codes, per-link
  obstacle counts); `nothing` when no overlay was applied
- connectivity_method: the method actually used (`:legacy` or `:on_path`)
- connectivity_diagnostics: `nothing` for `:legacy`, otherwise the
  `build_on_path_distance_matrix` diagnostics (roots, edge counts, elevation
  inversions, total link length)
- connectivity_tree_edges, connectivity_parent: directed dendritic edges and
  downstream-parent index of the selected graph (empty for `:legacy`).  These
  are the E12 interface for rebuilding dam/obstacle passability on the
  corrected network.
- cedex_var_df, cedex_esc_uts_df, obstacles_df: Optional loaded updated inputs
- biological_options: the resolved E13-E16/E18 option values actually applied
- fishfarm_eel_records_dropped: site codes whose AA record was dropped (E16)
- salinity_envelope_filtered_cells: site x species growth cells zeroed (E18)
"""
function prepare_ode_data(;
    connectivity_file::String = "data/ConnectivityUTM.csv",
    density_file::String = "data/BIOTIC/FishDensity_and_Juveniles_Matrix.csv",
    species_chars_file::String = "data/ABIOTIC/caracteristicas_peces_Guadalquivir_03-04-2018.csv",
    environmental_file::String = "data/ABIOTIC/Matriz_Ambiental_Data.csv",
    distance_file::String = "data/Matrix_distances_1037puntos_BRUTO_FINAL.csv",
    interaction_file::String = "data/BIOTIC/Interacciones_peces_Guadalquivir_03-04-2018_ENG.csv",
    upstream_cost::Float64 = 0.01,
    cedex_var_file::Union{Nothing,String} = nothing,
    cedex_esc_uts_file::Union{Nothing,String} = nothing,
    obstacles_file::Union{Nothing,String} = nothing,
    obstacle_mode::Symbol = :legacy,
    obstacle_matching_tolerance::Float64 = 2000.0,
    obstacle_passability::Float64 = 0.1,
    obstacle_downstream_passability::Float64 = 0.5,
    obstacle_include_ambiguous_status::Bool = false,
    connectivity_method::Symbol = :legacy,
    on_path_rtol::Float64 = ON_PATH_DEFAULT_RTOL,
    on_path_atol::Float64 = ON_PATH_DEFAULT_ATOL,
    heat_stress_rate::Float64 = 0.0,
    carrying_capacity_base_scaling::Float64 = 10.0,
    carrying_capacity_scaling::Float64 = 1.0,
    thermal_optima_override::Union{Nothing,AbstractVector} = nothing,
    water_temperature_file::Union{Nothing,String} = DEFAULT_WATER_TEMPERATURE_FILE,
    require_water_temperature::Bool = true,
    legacy_temperature_column::Union{Nothing,String} = nothing,
    absence_growth_fraction::Float64 = 0.1,
    nonreproducing_local_growth::Symbol = :legacy,
    exclude_fishfarm_eel_records::Bool = false,
    pool_capacity_mode::Symbol = :legacy,
    fishless_dificil_capacity::Symbol = :legacy,
    salinity_envelope::Bool = false
)
    println("="^60)
    println("Preparing data for ODE metacommunity model")
    println("="^60)

    # 1. Load site data
    println("\n[1/12] Loading site data...")
    site_df = load_site_data(connectivity_file, environmental_file)
    sites = String.(site_df.CODIGO)
    n_sites = length(sites)
    println("Found $n_sites sites")

    obstacle_mode in (:legacy, :overlay) || error("obstacle_mode must be :legacy or :overlay")
    connectivity_method in (:legacy, :on_path) ||
        error("connectivity_method must be :legacy or :on_path (got :$(connectivity_method))")
    # E13-E16/E18 opts.  `build_intrinsic_growth_rates` /
    # `build_carrying_capacity` re-validate, but checking here fails before any
    # expensive data preparation.
    absence_growth_fraction >= 0 ||
        error("absence_growth_fraction must be non-negative (got $absence_growth_fraction)")
    nonreproducing_local_growth in (:legacy, :zero) ||
        error("nonreproducing_local_growth must be :legacy or :zero (got :$(nonreproducing_local_growth))")
    pool_capacity_mode in (:legacy, :exclude, :cap) ||
        error("pool_capacity_mode must be :legacy, :exclude or :cap (got :$(pool_capacity_mode))")
    fishless_dificil_capacity in (:legacy, :near_zero) ||
        error("fishless_dificil_capacity must be :legacy or :near_zero (got :$(fishless_dificil_capacity))")

    # Optional updated inputs are loaded explicitly and kept separate from the
    # legacy site/environment tables.  This avoids silently changing the model
    # when a file is present but no crosswalk or transformation was selected.
    cedex_var_df = cedex_var_file === nothing ? nothing : load_cedex_var(cedex_var_file)
    cedex_esc_uts_df = cedex_esc_uts_file === nothing ? nothing : load_cedex_esc_uts(cedex_esc_uts_file)
    obstacles_df = obstacles_file === nothing ? nothing : load_obstacles(obstacles_file)

    # 2. Load species density data
    println("\n[2/12] Loading species density data...")
    density_df, species_codes = load_species_density_data(density_file)
    n_species = length(species_codes)
    println("Found $n_species species: $species_codes")

    # E16: optionally drop the four documented Guadiato fish-farm eel records
    # before ANY downstream builder (growth, capacity, initial state) sees them.
    fishfarm_eel_records_dropped = String[]
    if exclude_fishfarm_eel_records
        density_df, fishfarm_eel_records_dropped = drop_fishfarm_eel_records(density_df, site_df)
        println("E16 exclude_fishfarm_eel_records: dropped $(length(fishfarm_eel_records_dropped)) " *
                "fish-farm eel record(s): $(join(fishfarm_eel_records_dropped, ", "))")
    end

    # 3. Load species characteristics
    println("\n[3/12] Loading species characteristics...")
    species_chars_df = load_species_characteristics(species_chars_file)

    # Get thermal parameters for our species
    thermal_optima = Float64[]
    thermal_sigmas = Float64[]
    thermal_lower_limits = Float64[]
    thermal_upper_limits = Float64[]

    # Use thermal optimum, sigma (thermal breadth) and empirical limits from
    # species characteristics.  Sigma is derived from the temperature range in
    # the data (sigma ≈ range/6); the limits are the observed range bounds (WP3).
    sp_lookup = Dict(lowercase(r.SP) => (opt=r.thermal_optimum, sig=r.thermal_sigma,
        lower=hasproperty(r, :thermal_lower) ? r.thermal_lower : -Inf,
        upper=hasproperty(r, :thermal_upper) ? r.thermal_upper : Inf)
        for r in eachrow(species_chars_df))

    for sp in species_codes
        key = lowercase(sp)
        vals = get(sp_lookup, key, nothing)
        if vals !== nothing
            push!(thermal_optima, vals.opt)
            push!(thermal_sigmas, vals.sig)
            push!(thermal_lower_limits, vals.lower)
            push!(thermal_upper_limits, vals.upper)
        else
            push!(thermal_optima, 15.0)
            push!(thermal_sigmas, 3.0)
            push!(thermal_lower_limits, -Inf)
            push!(thermal_upper_limits, Inf)
        end
    end
    if thermal_optima_override !== nothing
        length(thermal_optima_override) == n_species ||
            error("thermal_optima_override has $(length(thermal_optima_override)) entries, expected $n_species")
        thermal_optima = Float64.(thermal_optima_override)
        println("Thermal optima overridden (WP3 optimum sweep): $thermal_optima")
    end
    println("Thermal optima: $thermal_optima")
    println("Thermal sigmas: $thermal_sigmas")
    println("Thermal limits (lower/upper): $thermal_lower_limits / $thermal_upper_limits")

    # 4. Load interaction matrix
    println("\n[4/12] Loading interaction matrix...")
    interaction_matrix = load_interaction_matrix(interaction_file, species_codes)

    # 5. Build distance matrix
    # Create subcatchment mapping from site data
    site_to_subcatchment = Dict{String, String}(string(row.CODIGO) => string(row.CODIGO_S) for row in eachrow(site_df))
    # Create distance to river mapping (Dist.Guadalq.(m))
    site_to_river_distance = Dict{String, Float64}(string(row.CODIGO) => Float64(coalesce(row."Dist.Guadalq.(m)", 0.0)) for row in eachrow(site_df))
    # Create elevation mapping
    site_to_elevation = Dict{String, Float64}(string(row.CODIGO) => Float64(coalesce(row.ALTITUD, 500.0)) for row in eachrow(site_df))

    println("\n[5/12] Building distance matrix (method=:$(connectivity_method))...")
    connectivity_diagnostics = nothing
    # E12 interface: the directed dendritic edges (upstream, downstream, metres)
    # and the downstream-parent index (0 = root) of the selected graph.  Empty
    # for the legacy method, which has no flow-oriented topology.
    connectivity_tree_edges = Tuple{Int, Int, Float64}[]
    connectivity_parent = Int[]
    if connectivity_method == :legacy
        distance_matrix = build_distance_matrix(distance_file, sites, site_to_subcatchment,
            site_to_river_distance, site_to_elevation)
    else
        on_path_result = build_on_path_distance_matrix(distance_file, sites,
            site_to_river_distance, site_to_elevation;
            rtol=on_path_rtol, atol=on_path_atol)
        distance_matrix = on_path_result.distance_matrix
        connectivity_diagnostics = on_path_result.diagnostics
        connectivity_tree_edges = on_path_result.tree_edges
        connectivity_parent = on_path_result.parent
        println("On-path graph: $(connectivity_diagnostics.n_edges) edges " *
                "($(connectivity_diagnostics.n_parent_edges) parent + " *
                "$(connectivity_diagnostics.n_root_join_edges) root-join), " *
                "$(connectivity_diagnostics.n_roots) roots, total " *
                "$(round(connectivity_diagnostics.total_length_km, digits=1)) km")
    end

    # 6. Build elevation vector
    println("\n[6/12] Extracting elevations...")
    elevations = build_elevation_vector(site_df, sites)
    println("Elevation range: $(minimum(elevations)) - $(maximum(elevations)) m")

    # 7. Build dam passability matrix
    println("\n[7/12] Building dam passability matrix...")
    legacy_dams = build_dam_passability_matrix(site_df, sites, distance_matrix, elevations)
    dams = copy(legacy_dams)
    obstacle_overlay = ones(Float64, n_sites, n_sites)
    obstacle_metadata = nothing
    obstacle_mapping_diagnostics = DataFrame(
        obstacle_id=String[], matched=Bool[], edge_from=String[], edge_to=String[],
        restricted_origin=String[], restricted_destination=String[],
        match_distance_m=Union{Missing,Float64}[], status=String[],
        status_class=String[], outcome=String[],
        duplicate_of_legacy_dam=Bool[], cap_emba_m3=Union{Missing,Float64}[]
    )
    if obstacles_df !== nothing && obstacle_mode == :overlay
        println("Applying obstacle passability overlay (tolerance=$(obstacle_matching_tolerance)m, upstream passability=$(obstacle_passability), downstream passability=$(obstacle_downstream_passability))...")
        # On the corrected network, snap obstacles only to the C4 tree edges and
        # use the tree's upstream/downstream orientation (E12); the legacy path
        # keeps the straight-line matching and elevation-based direction.
        obstacle_result = build_obstacle_passability_matrix(
            obstacles_df, site_df, sites, distance_matrix, elevations;
            matching_tolerance=obstacle_matching_tolerance,
            obstacle_passability=obstacle_passability,
            obstacle_downstream_passability=obstacle_downstream_passability,
            tree_edges=connectivity_tree_edges,
            legacy_dams=legacy_dams,
            include_ambiguous_status=obstacle_include_ambiguous_status
        )
        obstacle_overlay = obstacle_result.passability
        obstacle_mapping_diagnostics = obstacle_result.diagnostics
        obstacle_metadata = obstacle_result.metadata
        # Accumulate the distinct legacy-dam and obstacle layers rather than
        # taking the (non-compounding) elementwise minimum.
        dams = legacy_dams .* obstacle_overlay
        println("Obstacles: matched $(obstacle_metadata.n_matched)/$(obstacle_metadata.n_total), " *
                "applied $(obstacle_metadata.n_applied), excluded $(obstacle_metadata.n_excluded_total) " *
                "(non-operational $(obstacle_metadata.n_excluded_nonoperational), " *
                "out-of-basin $(obstacle_metadata.n_excluded_out_of_basin), " *
                "ambiguous $(obstacle_metadata.n_excluded_ambiguous_status), " *
                "legacy duplicates $(obstacle_metadata.n_excluded_duplicate_legacy))")
    end

    # 8. Extract environmental parameters
    println("\n[8/12] Extracting environmental parameters...")
    temperatures = extract_site_temperatures(site_df, sites;
        water_temperature_file=water_temperature_file,
        require_water_temperature=require_water_temperature,
        legacy_temperature_column=legacy_temperature_column)
    habitat_suitability = extract_habitat_suitability(site_df, sites)
    println("Site temperature level source: " *
            (water_temperature_file === nothing ? "explicit legacy / elevation" :
             "corrected water-temperature baseline $water_temperature_file"))
    println("Temperature range: $(minimum(temperatures)) - $(maximum(temperatures))")

    # 9. Build intrinsic growth rates
    println("\n[9/12] Building intrinsic growth rates...")
    intrinsic_growth_rates = build_intrinsic_growth_rates(density_df, species_codes, sites, species_chars_df;
        absence_growth_fraction=absence_growth_fraction,
        nonreproducing_local_growth=nonreproducing_local_growth)

    # E18: opt-in salinity envelope.  Applied multiplicatively to the local growth
    # matrix (equivalent to a per-species habitat factor), so the ODE is untouched.
    salinity_envelope_filtered_cells = 0
    if salinity_envelope
        intrinsic_growth_rates, salinity_envelope_filtered_cells = apply_salinity_envelope(
            intrinsic_growth_rates, species_codes, species_chars_df, site_df, sites)
        println("E18 salinity envelope: zeroed $(salinity_envelope_filtered_cells) site x species " *
                "growth cell(s) for strictly brackish species below " *
                "$(SALINITY_BRACKISH_MIN_CONDUCTIVITY_US_CM) uS/cm")
    end

    # 10. Build species dispersal scaling factors
    println("\n[10/12] Building species dispersal scaling factors...")
    dispersal_scaling = build_dispersal_scaling(species_codes)
    println("Dispersal scaling range: $(minimum(dispersal_scaling)) - $(maximum(dispersal_scaling))")
    println("Median-normalized scaling factors (median = 1.0)")

    # 11. Build carrying capacities from observed density data
    println("\n[11/12] Building carrying capacities...")
    carrying_capacity = build_carrying_capacity(density_df, site_df, sites, species_codes;
        scaling=carrying_capacity_base_scaling,
        pool_capacity_mode=pool_capacity_mode,
        fishless_dificil_capacity=fishless_dificil_capacity)
    println("Carrying capacity base scaling (observed-density multiplier): $(carrying_capacity_base_scaling)x")
    if carrying_capacity_scaling != 1.0
        carrying_capacity = carrying_capacity .* carrying_capacity_scaling
        println("Carrying capacity scaled by $(carrying_capacity_scaling) (WP4 sensitivity): " *
                "range $(minimum(carrying_capacity)) - $(maximum(carrying_capacity))")
    end

    # 12. Precompute dispersal matrix (using species-specific dispersal coefficients)
    println("\n[12/12] Precomputing dispersal matrix...")
    dispersal_matrix = precompute_dispersal_matrix(
        n_sites,
        Matrix(distance_matrix),
        elevations,
        upstream_cost,
        dams,
        species_codes
    )
    println("Dispersal matrix: $(Guadex.nnz(dispersal_matrix)) non-zero entries")

    # Create MetacommunityParams
    println("\n" * "="^60)
    println("Creating MetacommunityParams...")
    println("="^60)

    params = MetacommunityParams(
        n_sites,
        n_species,
        interaction_matrix,
        dispersal_matrix,
        dispersal_scaling,
        intrinsic_growth_rates,
        temperatures,
        habitat_suitability,
        thermal_optima,
        thermal_sigmas,
        carrying_capacity,
        thermal_lower_limits,
        thermal_upper_limits,
        heat_stress_rate
    )
    heat_stress_rate > 0 &&
        println("Heat-stress mortality enabled: k=$(heat_stress_rate) 1/day/degC^2")

    return (
        params = params,
        sites = sites,
        species = species_codes,
        distance_matrix = distance_matrix,
        elevations = elevations,
        dams = dams,
        legacy_dams = legacy_dams,
        obstacle_overlay = obstacle_overlay,
        obstacle_mapping_diagnostics = obstacle_mapping_diagnostics,
        # E12: obstacle classification/accumulation counts and per-link obstacle
        # counts (nothing when no overlay was applied).
        obstacle_metadata = obstacle_metadata,
        connectivity_method = connectivity_method,
        connectivity_diagnostics = connectivity_diagnostics,
        connectivity_tree_edges = connectivity_tree_edges,
        connectivity_parent = connectivity_parent,
        obstacles_df = obstacles_df,
        cedex_var_df = cedex_var_df,
        cedex_esc_uts_df = cedex_esc_uts_df,
        site_df = site_df,
        species_chars_df = species_chars_df,
        density_df = density_df,
        # E13-E16/E18: the resolved option values actually applied, for metadata
        # and the parameter digest.
        biological_options = (
            absence_growth_fraction = absence_growth_fraction,
            pool_capacity_mode = pool_capacity_mode,
            fishless_dificil_capacity = fishless_dificil_capacity,
            nonreproducing_local_growth = nonreproducing_local_growth,
            exclude_fishfarm_eel_records = exclude_fishfarm_eel_records,
            salinity_envelope = salinity_envelope
        ),
        fishfarm_eel_records_dropped = fishfarm_eel_records_dropped,
        salinity_envelope_filtered_cells = salinity_envelope_filtered_cells
    )
end

"""
    save_ode_data(data::NamedTuple, output_dir::String)

Save prepared ODE data to files for later use.
"""
function save_ode_data(data::NamedTuple, output_dir::String)
    println("Saving ODE data to $output_dir")

    # Save sites
    CSV.write(joinpath(output_dir, "sites.csv"), DataFrame(site=data.sites))

    # Save species
    CSV.write(joinpath(output_dir, "species.csv"), DataFrame(species=data.species))

    # Save distance matrix (as COO triplets)
    distance_df = DataFrame(
        i = Int[],
        j = Int[],
        distance = Float64[]
    )

    I, J, V = findnz(data.distance_matrix)
    for (i, j, v) in zip(I, J, V)
        push!(distance_df, (i, j, v))
    end
    CSV.write(joinpath(output_dir, "distance_matrix.csv"), distance_df)

    # Save elevations
    CSV.write(joinpath(output_dir, "elevations.csv"),
              DataFrame(site=data.sites, elevation=data.elevations))

    # Save temperatures
    CSV.write(joinpath(output_dir, "temperatures.csv"),
              DataFrame(site=data.sites, temperature=data.params.temperatures))

    # Save habitat suitability
    CSV.write(joinpath(output_dir, "habitat_suitability.csv"),
              DataFrame(site=data.sites, suitability=data.params.habitat_suitability))

    # Save interaction matrix
    CSV.write(joinpath(output_dir, "interaction_matrix.csv"),
              DataFrame(data.params.interaction_matrix, :auto))

    # Save growth rates
    CSV.write(joinpath(output_dir, "growth_rates.csv"),
              DataFrame(data.params.intrinsic_growth_rates, :auto))

    println("Data saved successfully!")
end
