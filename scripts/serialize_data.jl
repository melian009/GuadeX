# ============================================================================
# Data Serialization Script
#
# This script generates the serialized data file (data.jld2) required by
# the sensitivity analysis. Run this script once before executing run_sensitivity.jl
#
# Setup: First activate the local environment
#   julia> using Pkg
#   julia> Pkg.activate(".")
#   julia> Pkg.instantiate()
#   julia> include("serialize_data.jl")
# ============================================================================

using Pkg
Pkg.activate(".")

using JLD2
using DataFrames
using CSV
using SparseArrays
using LinearAlgebra
using DifferentialEquations
using Statistics

include("src/OdeModel.jl")
include("src/Visualization.jl")

using .OdeModel
using .Visualization

const SPECIES_CODES = [
  "SA", "LS", "ST", "SP", "IL", "PW", "CP", "AA", "AH", "LR", "MC", "AB", "IO",
  "OM", "LG", "GH", "AA", "CG", "CC", "MS", "AM", "TT", "EL", "GL"
]

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

function parse_temperature_range_and_sigma(temp_str::AbstractString)
    default_optimum = 15.0
    default_sigma = 3.0

    if isempty(temp_str) || temp_str == ""
        return (optimum=default_optimum, sigma=default_sigma)
    end

    if occursin(" to ", temp_str)
        parts = split(temp_str, " to ")
        if length(parts) == 2
            try
                t_min = parse(Float64, strip(parts[1]))
                t_max = parse(Float64, strip(parts[2]))
                thermal_range = t_max - t_min
                sigma = thermal_range / 6.0
                optimum = (t_min + t_max) / 2.0
                return (optimum=optimum, sigma=sigma)
            catch
                return (optimum=default_optimum, sigma=default_sigma)
            end
        end
    end

    try
        optimum = parse(Float64, temp_str)
        return (optimum=optimum, sigma=default_sigma)
    catch
        return (optimum=default_optimum, sigma=default_sigma)
    end
end

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

function load_species_characteristics(file::String)
    println("Loading species characteristics from: $file")
    df = CSV.read(file, DataFrame; delim=';')
    thermal_params = parse_temperature_range_and_sigma.(string.(df.TEMPERATURE_C))
    df.thermal_optimum = [p.optimum for p in thermal_params]
    df.thermal_sigma = [p.sigma for p in thermal_params]
    df.elevation_optimum = parse_elevation.(string.(df.ELEVATION_m))
    df.max_size_mm = Float64.(df.MAX_SIZE_mm)
    println("Loaded characteristics for $(nrow(df)) species")
    return df
end

function load_site_data(connectivity_file::String, environmental_file::String)
    println("Loading site data from connectivity file: $connectivity_file")
    connectivity_df = CSV.read(connectivity_file, DataFrame)
    connectivity_df = filter(row -> row.CODIGO != "1.30.20", connectivity_df)
    println("Loaded $(nrow(connectivity_df)) sites from connectivity data")
    println("Loading environmental data from: $environmental_file")
    environmental_df = CSV.read(environmental_file, DataFrame)
    println("Loaded $(nrow(environmental_df)) sites from environmental data")
    site_df = innerjoin(connectivity_df, environmental_df, on=:CODIGO, makeunique=true)
    site_df.Demb_arr_m = replace(site_df."Demb arr.(m)", "No existe" => "0")
    site_df.Demb_ab_m = replace(site_df."Demb ab.(m)", "No existe" => "0")
    site_df.Demb_arr_m = parse.(Float64, site_df.Demb_arr_m)
    site_df.Demb_ab_m = parse.(Float64, site_df.Demb_ab_m)
    println("Merged data for $(nrow(site_df)) sites")
    return site_df
end

function load_species_density_data(density_file::String)
    println("Loading species density data from: $density_file")
    density_df = CSV.read(density_file, DataFrame)
    density_cols = [c for c in names(density_df) if endswith(c, "_DEN")]
    species_codes = [replace(c, "_DEN" => "") for c in density_cols]
    println("Loaded density data for $(length(species_codes)) species at $(nrow(density_df)) sites")
    return density_df, species_codes
end

function load_interaction_matrix(interaction_file::String, species_codes::Vector{String})
    println("Loading interaction matrix from: $interaction_file")
    interaction_df = CSV.read(interaction_file, DataFrame; delim=';')
    rename!(interaction_df, 1 => :Species)
    matrix_species = names(interaction_df)[2:end]
    species_to_idx = Dict(sp => i for (i, sp) in enumerate(matrix_species))
    n_species = length(species_codes)
    interaction_matrix = zeros(n_species, n_species)

    function parse_interaction_string(interaction_str::Union{String, Missing})
        if ismissing(interaction_str) || isempty(interaction_str)
            return 0.0
        end
        interaction_str = strip(string(interaction_str))
        if interaction_str == "" || interaction_str == ";"
            return 0.0
        end
        interaction_lower = lowercase(interaction_str)
        if occursin("no coexist", interaction_lower)
            return -1.0
        end
        if occursin("displaces", interaction_lower)
            return -0.8
        end
        if occursin("predation", interaction_lower)
            return -0.5
        end
        if occursin("competition", interaction_lower) ||
           occursin("interfere", interaction_lower) ||
           occursin("interfiere", interaction_lower)
            return -0.3
        end
        if occursin("affects", interaction_lower)
            return -0.2
        end
        if occursin("coexist", interaction_lower) || occursin("neutral", interaction_lower)
            return 0.0
        end
        return 0.0
    end

    for row in eachrow(interaction_df)
        sp1 = row.Species
        sp1_lower = lowercase(sp1)
        if !haskey(species_to_idx, sp1)
            continue
        end
        lower_species_codes = lowercase.(species_codes)
        if sp1_lower ∈ lower_species_codes
            target_idx1 = findfirst(==(sp1_lower), lower_species_codes)
            for sp2 in matrix_species
                if hasproperty(row, Symbol(sp2))
                    interaction_str = row[Symbol(sp2)]
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

function build_distance_matrix(distance_file::String, sites::Vector{String},
                               site_to_subcatchment::Dict{T, T2},
                               site_to_river_distance::Dict{T3, Float64},
                               site_to_elevation::Dict{T4, Float64}) where T <: AbstractString where T2 <: AbstractString where T3 <: AbstractString where T4 <: AbstractString
    println("Building distance matrix from: $distance_file")
    site_to_idx = Dict(s => i for (i, s) in enumerate(sites))
    n_sites = length(sites)

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

    adjacent_pairs = Set{Tuple{String, String}}()
    outlet_sites = String[]

    for (sc, sites_in_sc) in subcatchment_to_sites
        sorted_sites = sort(sites_in_sc, by=s -> get(site_to_river_distance, s, Inf))
        for i in 1:(length(sorted_sites)-1)
            downstream = sorted_sites[i]
            upstream = sorted_sites[i+1]
            push!(adjacent_pairs, (upstream, downstream))
        end
        if length(sorted_sites) > 0
            outlet = sorted_sites[1]
            push!(outlet_sites, outlet)
        end
    end

    if length(outlet_sites) > 1
        sorted_outlets = sort(outlet_sites, by=s -> get(site_to_elevation, s, Inf))
        for i in 1:(length(sorted_outlets)-1)
            downstream_outlet = sorted_outlets[i]
            upstream_outlet = sorted_outlets[i+1]
            push!(adjacent_pairs, (upstream_outlet, downstream_outlet))
        end
    end

    println("Found $(length(outlet_sites)) outlet sites connected to form main river network")
    println("Found $(length(adjacent_pairs)) adjacent site pairs in the network")

    I = Int[]
    J = Int[]
    V = Float64[]
    total_rows = 0
    valid_connections = 0

    reader = CSV.File(distance_file; delim=';')
    for row in reader
        total_rows += 1
        origin = row.ID_ORIGIN
        dest = row.ID_DESTINATION
        if !haskey(site_to_idx, origin) || !haskey(site_to_idx, dest)
            continue
        end
        dist = row.RETICULAR_DIST
        if origin == dest || dist <= 0
            continue
        end
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

    distance_matrix = sparse(I, J, V, n_sites, n_sites)
    println("Distance matrix: $(nnz(distance_matrix)) non-zero entries out of $(n_sites*n_sites) possible")

    return distance_matrix
end

function build_elevation_vector(site_df::DataFrame, sites::Vector{String})
    site_to_elevation = Dict(row.CODIGO => row.ALTITUD for row in eachrow(site_df))
    elevations = Float64[]
    for site in sites
        if haskey(site_to_elevation, site)
            push!(elevations, site_to_elevation[site])
        else
            push!(elevations, 500.0)
        end
    end
    return elevations
end

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
                if (origin_info.dist_downstream < 1000_000 && origin_info.dist_downstream < d_ij) ||
                   (dest_info.dist_upstream < 1000_000 && dest_info.dist_upstream < d_ij)
                    dams[i, j] = 0.1
                end
            end
            if e_i >= e_j
                if (origin_info.dist_upstream < 1000_000 && origin_info.dist_upstream < d_ij) ||
                   (dest_info.dist_downstream < 1000_000 && dest_info.dist_downstream < d_ij)
                    dams[i, j] = 0.1
                end
            end
        end
    end
    return dams
end

function extract_site_temperatures(site_df::DataFrame, sites::Vector{String})
    temp_col = nothing
    for col in names(site_df)
        if occursin("TEMP", uppercase(col)) || occursin("TEMPERATURA", uppercase(col))
            temp_col = col
            break
        end
    end

    if temp_col === nothing
        println("No temperature column found, using elevation-based estimate")
        elevations = build_elevation_vector(site_df, sites)
        temps = 20.0 .- (elevations ./ 1000.0 .* 6.5)
        return temps
    end

    site_to_temp = Dict(row.CODIGO => row[Symbol(temp_col)] for row in eachrow(site_df))
    temperatures = Float64[]
    for site in sites
        if haskey(site_to_temp, site)
            push!(temperatures, site_to_temp[site])
        else
            push!(temperatures, 15.0)
        end
    end
    return temperatures
end

function extract_habitat_suitability(site_df::DataFrame, sites::Vector{String})
    site_to_iet = Dict(row.CODIGO => row.IET for row in eachrow(site_df))
    iet_values = collect(values(site_to_iet))
    iet_min, iet_max = minimum(iet_values), maximum(iet_values)

    suitability = Float64[]
    for site in sites
        if haskey(site_to_iet, site)
            iet = site_to_iet[site]
            suit = 1.0 - (iet - iet_min) / (iet_max - iet_min + 1e-6)
            push!(suitability, max(0.1, suit))
        else
            push!(suitability, 0.5)
        end
    end
    return suitability
end

function build_intrinsic_growth_rates(density_df::DataFrame, species_codes::Vector{String},
                                       sites::Vector{String}, species_chars_df::DataFrame)
    n_sites = length(sites)
    n_species = length(species_codes)

    site_to_row = Dict(row.CODIGO => rownum for (rownum, row) in enumerate(eachrow(density_df)))

    annual_growth_rates = Dict{String, Float64}(
        "AB" => 1.2,
        "AH" => 1.0,
        "SP" => 0.55,
        "PW" => 0.3,
        "LS" => 0.2,
        "SA" => 0.7,
        "IL" => 0.55,
        "CP" => 0.65,
        "IO" => 0.45,
        "GH" => 4.0,
        "MS" => 0.4,
        "LG" => 0.5,
        "CC" => 0.4,
        "CG" => 0.8,
        "AM" => 0.35,
        "OM" => 0.4,
        "EL" => 0.3,
        "GL" => 0.8,
        "TT" => 0.3,
        "AA" => 0.1,
        "MC" => 0.3,
        "LR" => 0.3,
        "ST" => 0.4,
    )

    growth_rates = zeros(n_sites, n_species)

    for (s_idx, sp_code) in enumerate(species_codes)
        r_annual = get(annual_growth_rates, sp_code, 0.5)
        r_daily = r_annual / 365.0

        for (site_idx, site) in enumerate(sites)
            if haskey(site_to_row, site)
                row_idx = site_to_row[site]
                density_col = Symbol("$(sp_code)_DEN")

                if hasproperty(density_df, density_col)
                    density = density_df[row_idx, density_col]
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

function build_dispersal_scaling(species_codes::Vector{String})
    annual_dispersal_rates = Dict{String, Float64}(
        "AB" => 0.3,
        "AH" => 0.2,
        "SP" => 3.0,
        "PW" => 27.5,
        "LS" => 20.0,
        "SA" => 5.0,
        "IL" => 1.25,
        "CP" => 1.25,
        "IO" => 0.5,
        "GH" => 25.0,
        "MS" => 17.5,
        "LG" => 25.0,
        "CC" => 20.0,
        "CG" => 10.0,
        "AM" => 10.0,
        "OM" => 12.5,
        "EL" => 12.5,
        "GL" => 3.0,
        "TT" => 6.0,
        "AA" => 5.0,
        "MC" => 75.0,
        "LR" => 75.0,
        "ST" => 12.5,
    )

    rates = [get(annual_dispersal_rates, sp, 5.0) for sp in species_codes]
    median_rate = median(rates)

    scaling = Float64[]
    for sp in species_codes
        rate = get(annual_dispersal_rates, sp, 5.0)
        push!(scaling, rate / median_rate)
    end

    return scaling
end

function build_carrying_capacity(density_df::DataFrame, site_df::DataFrame, sites::Vector{String}, species_codes::Vector{String})
    println("Building site-specific carrying capacities from density data...")

    site_to_idx = Dict{String, Int}()
    for (rownum, row) in enumerate(eachrow(density_df))
        site_to_idx[row.CODIGO] = rownum
    end

    density_cols = [Symbol("$(sp)_DEN") for sp in species_codes]
    K_scaling = 10.0

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

    nonzero_caps = filter(c -> c > 0, raw_capacities)
    K_min = 50.0
    K_floor = isempty(nonzero_caps) ? K_min : max(K_min, quantile(nonzero_caps, 0.1))

    carrying_capacity = Float64[]
    for cap in raw_capacities
        push!(carrying_capacity, max(cap, K_floor))
    end

    println("Carrying capacity range: $(minimum(carrying_capacity)) - $(maximum(carrying_capacity))")
    println("Mean carrying capacity: $(mean(carrying_capacity))")
    println("K floor (10th percentile of non-zero K, min $K_min): $K_floor")

    return carrying_capacity
end

function prepare_data(;
    connectivity_file::String = "data/ConnectivityUTM.csv",
    density_file::String = "data/BIOTIC/FishDensity_and_Juveniles_Matrix.csv",
    species_chars_file::String = "data/ABIOTIC/caracteristicas_peces_Guadalquivir_03-04-2018.csv",
    environmental_file::String = "data/ABIOTIC/Matriz_Ambiental_Data.csv",
    distance_file::String = "data/Matrix_distances_1037puntos_BRUTO_FINAL.csv",
    interaction_file::String = "data/BIOTIC/Interacciones_peces_Guadalquivir_03-04-2018_ENG.csv",
    upstream_cost::Float64 = 0.01
)
    println("="^60)
    println("Preparing data for sensitivity analysis")
    println("="^60)

    println("\n[1/12] Loading site data...")
    site_df = load_site_data(connectivity_file, environmental_file)
    sites = String.(site_df.CODIGO)
    n_sites = length(sites)
    println("Found $n_sites sites")

    println("\n[2/12] Loading species density data...")
    density_df, species_codes = load_species_density_data(density_file)
    n_species = length(species_codes)
    println("Found $n_species species: $species_codes")

    println("\n[3/12] Loading species characteristics...")
    species_chars_df = load_species_characteristics(species_chars_file)

    thermal_optima = Float64[]
    thermal_sigmas = Float64[]
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

    println("\n[4/12] Loading interaction matrix...")
    interaction_matrix = load_interaction_matrix(interaction_file, species_codes)

    site_to_subcatchment = Dict{String, String}(string(row.CODIGO) => string(row.CODIGO_S) for row in eachrow(site_df))
    site_to_river_distance = Dict{String, Float64}(string(row.CODIGO) => Float64(coalesce(row."Dist.Guadalq.(m)", 0.0)) for row in eachrow(site_df))
    site_to_elevation = Dict{String, Float64}(string(row.CODIGO) => Float64(coalesce(row.ALTITUD, 500.0)) for row in eachrow(site_df))

    println("\n[5/12] Building distance matrix...")
    distance_matrix = build_distance_matrix(distance_file, sites, site_to_subcatchment, site_to_river_distance, site_to_elevation)

    println("\n[6/12] Extracting elevations...")
    elevations = build_elevation_vector(site_df, sites)
    println("Elevation range: $(minimum(elevations)) - $(maximum(elevations)) m")

    println("\n[7/12] Building dam passability matrix...")
    dams = build_dam_passability_matrix(site_df, sites, distance_matrix, elevations)

    println("\n[8/12] Extracting environmental parameters...")
    temperatures = extract_site_temperatures(site_df, sites)
    habitat_suitability = extract_habitat_suitability(site_df, sites)
    println("Temperature range: $(minimum(temperatures)) - $(maximum(temperatures))")

    println("\n[9/12] Building intrinsic growth rates...")
    intrinsic_growth_rates = build_intrinsic_growth_rates(density_df, species_codes, sites, species_chars_df)

    println("\n[10/12] Building species dispersal scaling factors...")
    dispersal_scaling = build_dispersal_scaling(species_codes)
    println("Dispersal scaling range: $(minimum(dispersal_scaling)) - $(maximum(dispersal_scaling))")
    println("Median-normalized scaling factors (median = 1.0)")

    println("\n[11/12] Building carrying capacities...")
    carrying_capacity = build_carrying_capacity(density_df, site_df, sites, species_codes)

    println("\n[12/12] Precomputing dispersal matrix...")
    dispersal_matrix = precompute_dispersal_matrix(
        n_sites,
        Matrix(distance_matrix),
        elevations,
        upstream_cost,
        dams,
        species_codes
    )
    println("Dispersal matrix: $(nnz(dispersal_matrix)) non-zero entries")

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

println("\n" * "="^60)
println("Serializing model data for sensitivity analysis")
println("="^60)

data = prepare_data(upstream_cost = 0.05)

println("\n[Computing u0] Extracting initial conditions from density data...")
density_cols = [Symbol("$(sp)_DEN") for sp in data.species]
density_df_filtered = filter(row -> row.CODIGO in data.sites, data.density_df)
u0 = Matrix(density_df_filtered[:, density_cols])
replace!(u0, NaN => 0.0)
u0 = max.(u0, 0.0)
u0_flat = vec(u0)

println("u0 dimensions: $(length(u0_flat)) (n_sites * n_species = $(data.params.n_sites) * $(data.params.n_species))")
println("u0 range: $(minimum(u0_flat)) - $(maximum(u0_flat))")

output_file = joinpath("data", "data.jld2")
println("\nSaving serialized data to: $output_file")

mkpath("data")

jldsave(output_file;
    params = data.params,
    sites = data.sites,
    species = data.species,
    distance_matrix = data.distance_matrix,
    elevations = data.elevations,
    dams = data.dams,
    site_df = data.site_df,
    density_df = data.density_df,
    u0 = u0_flat,
    n_sites = data.params.n_sites,
    n_species = data.params.n_species
)

println("\n" * "="^60)
println("Data serialization complete!")
println("="^60)
println("\nYou can now run the sensitivity analysis:")
println("  julia> include(\"run_sensitivity.jl\")")