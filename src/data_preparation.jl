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
                         site_to_river_distance::Dict{String, Float64})

Distances between sites.
Build a sparse distance matrix from the distance file.
Only includes connections between ADJACENT sites within the same subcatchment.

This is critical because:
1. The river network is dendritic - sites in different subcatchments are not connected
2. Within each subcatchment, sites form a linear chain where each site only connects
   to its immediate upstream and downstream neighbors (like a river reach)
"""
function build_distance_matrix(distance_file::String, sites::Vector{String},
                               site_to_subcatchment::Dict{T, T2},
                               site_to_river_distance::Dict{T3, Float64}) where T <: AbstractString where T2 <: AbstractString where T3 <: AbstractString
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

    for (sc, sites_in_sc) in subcatchment_to_sites
        # Sort sites by distance to river (ascending = downstream first)
        sorted_sites = sort(sites_in_sc, by=s -> get(site_to_river_distance, s, Inf))

        # Connect each site to its adjacent neighbor (upstream <-> downstream)
        for i in 1:(length(sorted_sites)-1)
            downstream = sorted_sites[i]      # Closer to river
            upstream = sorted_sites[i+1]       # Farther from river
            push!(adjacent_pairs, (upstream, downstream))
        end
    end

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
Uses species max size and presence/absence to estimate base growth rates.
"""
function build_intrinsic_growth_rates(density_df::DataFrame, species_codes::Vector{String},
                                       sites::Vector{String}, species_chars_df::DataFrame)

    n_sites = length(sites)
    n_species = length(species_codes)

    # Create site to density row mapping
    site_to_row = Dict(row.CODIGO => rownum for (rownum, row) in enumerate(eachrow(density_df)))

    # Get max sizes for each species
    species_to_max_size = Dict(row.SP => row.max_size_mm for row in eachrow(species_chars_df))

    # Base growth rate is proportional to 1/max_size (smaller species grow faster)
    # Then scale by density presence (sites with species get positive growth)

    growth_rates = zeros(n_sites, n_species)

    for (s_idx, sp_code) in enumerate(species_codes)
        # Get base growth rate from max size
        max_size = get(species_to_max_size, sp_code, 100.0)
        base_rate = 1.0 / max_size * 10.0  # Scale factor

        for (site_idx, site) in enumerate(sites)
            if haskey(site_to_row, site)
                row_idx = site_to_row[site]
                density_col = Symbol("$(sp_code)_DEN")

                if hasproperty(density_df, density_col)
                    density = density_df[row_idx, density_col]

                    # If species is present (density > 0), use base rate
                    # If absent, use very low rate (potential colonization)
                    if density > 0
                        growth_rates[site_idx, s_idx] = base_rate
                    else
                        growth_rates[site_idx, s_idx] = base_rate * 0.1
                    end
                end
            end
        end
    end

    return growth_rates
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
    println("\n[1/8] Loading site data...")
    site_df = load_site_data(connectivity_file, environmental_file)
    sites = String.(site_df.CODIGO)
    n_sites = length(sites)
    println("Found $n_sites sites")

    # 2. Load species density data
    println("\n[2/8] Loading species density data...")
    density_df, species_codes = load_species_density_data(density_file)
    n_species = length(species_codes)
    println("Found $n_species species: $species_codes")

    # 3. Load species characteristics
    println("\n[3/8] Loading species characteristics...")
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
    println("\n[4/8] Loading interaction matrix...")
    interaction_matrix = load_interaction_matrix(interaction_file, species_codes)

    # 5. Build distance matrix
    # Create subcatchment mapping from site data
    site_to_subcatchment = Dict(row.CODIGO => string(row.CODIGO_S) for row in eachrow(site_df))
    # Create distance to river mapping (Dist.Guadalq.(m))
    site_to_river_distance = Dict(row.CODIGO => row."Dist.Guadalq.(m)" for row in eachrow(site_df))

    println("\n[5/8] Building distance matrix...")
    distance_matrix = build_distance_matrix(distance_file, sites, site_to_subcatchment, site_to_river_distance)

    # 6. Build elevation vector
    println("\n[6/8] Extracting elevations...")
    elevations = build_elevation_vector(site_df, sites)
    println("Elevation range: $(minimum(elevations)) - $(maximum(elevations)) m")

    # 7. Build dam passability matrix
    println("\n[7/8] Building dam passability matrix...")
    dams = build_dam_passability_matrix(site_df, sites)

    # 8. Extract environmental parameters
    println("\n[8/8] Extracting environmental parameters...")
    temperatures = extract_site_temperatures(site_df, sites)
    habitat_suitability = extract_habitat_suitability(site_df, sites)
    println("Temperature range: $(minimum(temperatures)) - $(maximum(temperatures))")

    # 9. Build intrinsic growth rates
    println("\n[9/8] Building intrinsic growth rates...")
    intrinsic_growth_rates = build_intrinsic_growth_rates(density_df, species_codes, sites, species_chars_df)

    # 10. Precompute dispersal matrix
    println("\n[10/8] Precomputing dispersal matrix...")
    dispersal_matrix = precompute_dispersal_matrix(
        n_sites,
        Matrix(distance_matrix),
        elevations,
        upstream_cost,
        dispersal_intensity,
        dams
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
        intrinsic_growth_rates,
        temperatures,
        habitat_suitability,
        thermal_optima,
        thermal_sigmas
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