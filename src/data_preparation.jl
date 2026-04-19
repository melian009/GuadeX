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

Parse temperature range string like "8 to 30" and return both the thermal optimum
(midpoint) and thermal breadth (sigma). Sigma is derived from the temperature range
using the approximation sigma ≈ range / 6.

Input file is "data/ABIOTIC/caracteristicas_peces_Guadalquivir_03-04-2018.csv"
"""
function parse_temperature_range_and_sigma(temp_str::AbstractString)
    # Default values
    default_optimum = 15.0
    default_sigma = 3.0

    if isempty(temp_str) || temp_str == ""
        return (optimum=default_optimum, sigma=default_sigma)
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
                return (optimum=optimum, sigma=sigma)
            catch
                return (optimum=default_optimum, sigma=default_sigma)
            end
        end
    end

    # Try to parse as single number - no range info, use default sigma
    try
        optimum = parse(Float64, temp_str)
        return (optimum=optimum, sigma=default_sigma)
    catch
        return (optimum=default_optimum, sigma=default_sigma)
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

    # Parse temperature ranges to get thermal optima and sigma (thermal breadth)
    thermal_params = parse_temperature_range_and_sigma.(string.(df.TEMPERATURE_C))
    df.thermal_optimum = [p.optimum for p in thermal_params]
    df.thermal_sigma = [p.sigma for p in thermal_params]

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
    load_interaction_matrix(interaction_file::String, species_codes::Vector{String})

Load and parse the species interaction matrix.
Returns a numeric interaction matrix where:
- Positive values indicate facilitation/positive effect
- Negative values indicate competition/predation
- Zero indicates neutral coexistence
"""
function load_interaction_matrix(interaction_file::String, species_codes::Vector{String})
    println("Loading interaction matrix from: $interaction_file")

    # Read the interaction matrix (semicolon delimited, with row names in first column)
    interaction_df = CSV.read(interaction_file, DataFrame; delim=';')

    # The first column contains species codes (row names)
    rename!(interaction_df, 1 => :Species)

    # Get species in the matrix (same order as columns without the first row)
    matrix_species = names(interaction_df)[2:end]

    # Create a mapping from species code to column index
    species_to_idx = Dict(sp => i for (i, sp) in enumerate(matrix_species))

    # Initialize interaction matrix
    n_species = length(species_codes)
    interaction_matrix = zeros(n_species, n_species)

    # Parse interaction values
    for row in eachrow(interaction_df)
        sp1 = row.Species
        sp1_lower = lowercase(sp1)
        if !haskey(species_to_idx, sp1)
            continue
        end

        # Find this species in our target species list
        lower_species_codes = lowercase.(species_codes)
        if sp1_lower ∈ lower_species_codes
            target_idx1 = findfirst(==(sp1_lower), lower_species_codes)

            for sp2 in matrix_species
                if hasproperty(row, Symbol(sp2))
                    interaction_str = row[Symbol(sp2)]

                    # Parse interaction string
                    value = parse_interaction_string(interaction_str)

                    sp2lower = lowercase(sp2)
                    if sp2lower ∈ lower_species_codes
                        target_idx2 = findfirst(==(sp2lower), lower_species_codes)
                        interaction_matrix[target_idx1, target_idx2] = value
                    end
                end
            end
        end
    end

    println("Created $(n_species)x$(n_species) interaction matrix")
    return interaction_matrix
end

"""
    parse_interaction_string(interaction_str::String)

Parse an interaction string and return a numeric value.
- "No coexist" => strong negative (-1.0)
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

    # 1. Strongest negative: No coexistence
    if occursin("no coexist", interaction_lower)
        return -1.0
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
    build_dam_passability_matrix(site_df::DataFrame, sites::Vector{String})

Build a dam passability matrix based on dam distances.
Lower passability if there are dams nearby.
"""
function build_dam_passability_matrix(site_df::DataFrame, sites::Vector{String})
    n_sites = length(sites)

    # Create site index mapping
    # site_to_idx = Dict(s => i for (i, s) in enumerate(sites))

    # Get dam information per site
    site_dam_info = Dict{String, NamedTuple{(:dist_upstream, :dist_downstream), Tuple{Float64, Float64}}}()

    for row in eachrow(site_df)
        codigo = row.CODIGO
        dist_up = row.Demb_arr_m
        dist_down = row.Demb_ab_m

        # If "No existe" was replaced with 0, it means no dam
        # We use a large value to indicate no dam effect
        if dist_up == 0
            dist_up = 1000_000.0  # No upstream dam
        end
        if dist_down == 0
            dist_down = 1000_000.0  # No downstream dam
        end

        site_dam_info[codigo] = (dist_upstream=dist_up, dist_downstream=dist_down)
    end

    # Build passability matrix
    # Passability = 1.0 if no dam, decreases with dam proximity
    dams = ones(n_sites, n_sites)

    for j in 1:n_sites  # origin
        for i in 1:n_sites  # destination
            if i == j
                continue
            end

            origin = sites[j]
            dest = sites[i]

            if haskey(site_dam_info, origin) && haskey(site_dam_info, dest)
                # Check if there's a dam between origin and destination
                # For now, use simple heuristic: if origin has downstream dam at all, reduce passability
                dam_info = site_dam_info[origin]
                if dam_info.dist_downstream < 1000_000
                    # Dam downstream of origin - reduce passability
                    dams[i, j] = 0.1
                end
            end
        end
    end

    return dams
end

"""
    extract_site_temperatures(site_df::DataFrame, sites::Vector{String})

Extract temperature data for each site.
Uses TEMP_MEDIA_SC (mean temperature of subcatchment) if available.
"""
function extract_site_temperatures(site_df::DataFrame, sites::Vector{String})
    # Try to find temperature column
    temp_col = nothing

    # Check for various temperature column names
    for col in names(site_df)
        if occursin("TEMP", uppercase(col)) || occursin("TEMPERATURA", uppercase(col))
            temp_col = col
            break
        end
    end

    if temp_col === nothing
        # Use elevation-based temperature estimate
        # Rough approximation: temperature decreases ~6.5°C per 1000m
        println("No temperature column found, using elevation-based estimate")
        elevations = build_elevation_vector(site_df, sites)
        temps = 20.0 .- (elevations ./ 1000.0 .* 6.5)
        return temps
    end

    # Extract temperatures
    site_to_temp = Dict(row.CODIGO => row[Symbol(temp_col)] for row in eachrow(site_df))

    temperatures = Float64[]
    for site in sites
        if haskey(site_to_temp, site)
            push!(temperatures, site_to_temp[site])
        else
            push!(temperatures, 15.0)  # Default temperature
        end
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

"""
    build_intrinsic_growth_rates(density_df::DataFrame, species_codes::Vector{String}, sites::Vector{String})

Build intrinsic growth rates matrix from density data.
Uses literature-derived daily intrinsic growth rates for each species.

# Literature Sources
- References are based on empirical studies from the Guadalquivir River Basin
- Annual rates are converted to daily rates using: r_daily = r_annual / 365
- For seasonal rates (e.g., Gambusia), converted using: r_daily = log(1 + r_seasonal) / 183

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
                                       sites::Vector{String}, species_chars_df::DataFrame)

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

        # Convert annual rate to daily instantaneous rate
        # r_daily = r_annual / 365
        # This assumes continuous exponential growth, which is appropriate for
        # fish populations with overlapping generations
        r_daily = r_annual / 365.0

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
                        growth_rates[site_idx, s_idx] = r_daily * 0.1
                    end
                end
            end
        end
    end

    return growth_rates
end

"""
    build_carrying_capacity(density_df::DataFrame, site_df::DataFrame, sites::Vector{String}, species_codes::Vector{String})

Build site-specific carrying capacities from observed fish density data.

The carrying capacity K_i for each site is derived from the observed total fish density,
scaled by a factor to represent the maximum sustainable biomass the site can support.
The scaling factor accounts for:
- Natural fluctuations around observed densities
- Additional habitat not sampled during surveys
- Density-dependent regulation allowing populations to exceed observed levels

# Arguments
- `density_df`: DataFrame with species density data (from load_species_density_data)
- `site_df`: DataFrame with site data
- `sites`: Vector of site codes in order
- `species_codes`: Vector of species codes

# Returns
- Vector of carrying capacities for each site
"""
function build_carrying_capacity(density_df::DataFrame, site_df::DataFrame, sites::Vector{String}, species_codes::Vector{String})
    println("Building site-specific carrying capacities from density data...")

    site_to_idx = Dict{String, Int}()
    for (rownum, row) in enumerate(eachrow(density_df))
        site_to_idx[row.CODIGO] = rownum
    end

    density_cols = [Symbol("$(sp)_DEN") for sp in species_codes]

    K_scaling = 1.5

    carrying_capacity = Float64[]
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
            push!(carrying_capacity, total_density * K_scaling)
        else
            push!(carrying_capacity, 10.0)
        end
    end

    println("Carrying capacity range: $(minimum(carrying_capacity)) - $(maximum(carrying_capacity))")
    println("Mean carrying capacity: $(mean(carrying_capacity))")

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
    prepare_ode_data(;
        connectivity_file::String = "data/ConnectivityUTM.csv",
        density_file::String = "data/BIOTIC/FishDensity_and_Juveniles_Matrix.csv",
        species_chars_file::String = "data/ABIOTIC/caracteristicas_peces_Guadalquivir_03-04-2018.csv",
        environmental_file::String = "data/ABIOTIC/Matriz_Ambiental_Data.csv",
        distance_file::String = "data/Matrix_distances_1037puntos_BRUTO_FINAL.csv",
        interaction_file::String = "data/BIOTIC/Interacciones_peces_Guadalquivir_03-04-2018_ENG.csv",
        upstream_cost::Float64 = 0.01,
        dispersal_intensity::Float64 = 0.1
    )

Prepare all data needed for the ODE metacommunity model.
## Arguments
- `connectivity_file`: Path to site connectivity data
- `density_file`: Path to species density data
- `species_chars_file`: Path to species characteristics data
- `environmental_file`: Path to environmental data
- `distance_file`: Path to distance matrix data
- `interaction_file`: Path to species interaction data
- `max_distance`: Maximum distance for connections in the distance matrix
- `upstream_cost`: Additional cost factor for upstream dispersal
- `dispersal_intensity`: Scaling factor for dispersal rates in the model

## Returns a NamedTuple with:
- params: MetacommunityParams struct
- sites: Vector of site codes
- species: Vector of species codes
- distance_matrix: Sparse distance matrix
- elevations: Vector of elevations
- dams: Dam passability matrix
"""
function prepare_ode_data(;
    connectivity_file::String = "data/ConnectivityUTM.csv",
    density_file::String = "data/BIOTIC/FishDensity_and_Juveniles_Matrix.csv",
    species_chars_file::String = "data/ABIOTIC/caracteristicas_peces_Guadalquivir_03-04-2018.csv",
    environmental_file::String = "data/ABIOTIC/Matriz_Ambiental_Data.csv",
    distance_file::String = "data/Matrix_distances_1037puntos_BRUTO_FINAL.csv",
    interaction_file::String = "data/BIOTIC/Interacciones_peces_Guadalquivir_03-04-2018_ENG.csv",
    upstream_cost::Float64 = 0.01,
    dispersal_intensity::Float64 = 0.1
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

    # 2. Load species density data
    println("\n[2/12] Loading species density data...")
    density_df, species_codes = load_species_density_data(density_file)
    n_species = length(species_codes)
    println("Found $n_species species: $species_codes")

    # 3. Load species characteristics
    println("\n[3/12] Loading species characteristics...")
    species_chars_df = load_species_characteristics(species_chars_file)

    # Get thermal parameters for our species
    thermal_optima = Float64[]
    thermal_sigmas = Float64[]

    # Use thermal optimum and sigma (thermal breadth) from species characteristics
    # Sigma is derived from the temperature range in the data (sigma ≈ range/6)
    for sp in species_codes
        sp_row = findfirst(row -> row.SP == sp, eachrow(species_chars_df))
        if sp_row !== nothing
            push!(thermal_optima, species_chars_df.thermal_optimum[sp_row])
            push!(thermal_sigmas, species_chars_df.thermal_sigma[sp_row])
        else
            push!(thermal_optima, 15.0)
            push!(thermal_sigmas, 3.0)
        end
    end
    println("Thermal optima: $thermal_optima")
    println("Thermal sigmas: $thermal_sigmas")

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

    println("\n[5/12] Building distance matrix...")
    distance_matrix = build_distance_matrix(distance_file, sites, site_to_subcatchment, site_to_river_distance, site_to_elevation)

    # 6. Build elevation vector
    println("\n[6/12] Extracting elevations...")
    elevations = build_elevation_vector(site_df, sites)
    println("Elevation range: $(minimum(elevations)) - $(maximum(elevations)) m")

    # 7. Build dam passability matrix
    println("\n[7/12] Building dam passability matrix...")
    dams = build_dam_passability_matrix(site_df, sites)

    # 8. Extract environmental parameters
    println("\n[8/12] Extracting environmental parameters...")
    temperatures = extract_site_temperatures(site_df, sites)
    habitat_suitability = extract_habitat_suitability(site_df, sites)
    println("Temperature range: $(minimum(temperatures)) - $(maximum(temperatures))")

    # 9. Build intrinsic growth rates
    println("\n[9/12] Building intrinsic growth rates...")
    intrinsic_growth_rates = build_intrinsic_growth_rates(density_df, species_codes, sites, species_chars_df)

    # 10. Build species dispersal scaling factors
    println("\n[10/12] Building species dispersal scaling factors...")
    dispersal_scaling = build_dispersal_scaling(species_codes)
    println("Dispersal scaling range: $(minimum(dispersal_scaling)) - $(maximum(dispersal_scaling))")
    println("Median-normalized scaling factors (median = 1.0)")

    # 11. Build carrying capacities from observed density data
    println("\n[11/12] Building carrying capacities...")
    carrying_capacity = build_carrying_capacity(density_df, site_df, sites, species_codes)

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
        carrying_capacity
    )

    return (
        params = params,
        sites = sites,
        species = species_codes,
        distance_matrix = distance_matrix,
        elevations = elevations,
        dams = dams,
        site_df = site_df,
        species_chars_df = species_chars_df,
        density_df = density_df
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