using Pkg; Pkg.activate(".")
using JLD2
using CSV
using DataFrames
using Guadex

# =============================================================================
# --- Standalone exporter: JLD2 simulation output -> four-level + viewer ---
#
# Usage:
#   julia --project=. scripts/export_run_outputs.jl <simulation_output.jld2> [output_dir]
#
# Reads a simulation saved by the entry scripts (sol_t / sol_u / sites / species,
# optionally warming / years / temperature_baseline / scenario / gcm) and writes
# the four reporting levels and the viz/ viewer files next to it (or to
# `output_dir`).  Native/invasive species groups come from parameters.toml.
# =============================================================================

const _GUADEX_ROOT = isfile(joinpath(@__DIR__, "..", "parameters.jl")) ? dirname(@__DIR__) : @__DIR__
include(joinpath(_GUADEX_ROOT, "parameters.jl"))
using .SimulationParameters
const GUADEX_PARAMS = SimulationParameters.load()
SimulationParameters.require_sections(GUADEX_PARAMS, "general", "species")

const DAYS_PER_YEAR = Int(GUADEX_PARAMS["general"]["days_per_year"])
const NATIVE_SPECIES = String.(GUADEX_PARAMS["species"]["native"])
const INVASIVE_SPECIES = String.(GUADEX_PARAMS["species"]["invasive"])
const CROSSWALK_PATH = joinpath(_GUADEX_ROOT, "data", "site_waterbody_crosswalk.csv")

if isempty(ARGS) || ARGS[1] in ("-h", "--help")
    println("Usage: julia --project=. scripts/export_run_outputs.jl <simulation_output.jld2> [output_dir]")
    exit(isempty(ARGS) ? 1 : 0)
end

input_path = ARGS[1]
isfile(input_path) || error("simulation file not found: $input_path")
output_dir = length(ARGS) >= 2 ? ARGS[2] : joinpath(dirname(input_path), "export")

data = JLD2.load(input_path)
println("Loaded $input_path")
println("  keys: $(join(sort(collect(keys(data))), ", "))")

for required in ("sol_t", "sol_u", "sites", "species")
    haskey(data, required) || error("simulation file is missing '$required'")
end

sol_t = data["sol_t"]
sol_u = data["sol_u"]
sites = String.(data["sites"])
species = String.(data["species"])

warming = haskey(data, "warming") ? data["warming"] : nothing
temperature_baseline = haskey(data, "temperature_baseline") ? data["temperature_baseline"] : nothing
years = haskey(data, "years") ? Int.(data["years"]) : nothing
start_year = haskey(data, "start_year") ? Int(data["start_year"]) : nothing

year_offsets = nothing
year_labels = nothing
if years !== nothing
    year_labels = years
    year_offsets = 0:(length(years) - 1)
else
    last_year = floor(Int, sol_t[end] / DAYS_PER_YEAR)
    year_offsets = collect(0:last_year)
    year_labels = start_year === nothing ? year_offsets : collect(start_year:(start_year + last_year))
end

if warming !== nothing && temperature_baseline !== nothing
    size(warming, 2) == length(year_offsets) ||
        error("warming has $(size(warming, 2)) columns but $(length(year_offsets)) report years were found")
end

metadata = Dict{String,Any}("script" => "scripts/export_run_outputs.jl", "source" => input_path)
for key in ("scenario", "gcm", "start_year", "end_year", "save_interval_days",
    "elevation_scaling", "temperature_increase", "upstream_cost",
    "passability_scenario", "interaction_matrix_type", "obstacle_mode",
    "obstacle_matched_count", "obstacle_total_count", "simulation_years")
    haskey(data, key) && (metadata[key] = data[key])
end

result = export_run_outputs(output_dir;
    sol_t=sol_t, sol_u=sol_u,
    sites=sites, species=species,
    site_df=nothing,
    crosswalk_path=CROSSWALK_PATH,
    native_species=NATIVE_SPECIES,
    invasive_species=INVASIVE_SPECIES,
    days_per_year=DAYS_PER_YEAR,
    report_year_offsets=year_offsets,
    report_year_labels=year_labels,
    temperature_baseline=temperature_baseline,
    warming=warming,
    dams=haskey(data, "dams") ? data["dams"] : nothing,
    habitat_suitability=haskey(data, "habitat_suitability") ? data["habitat_suitability"] : nothing,
    upstream_cost=haskey(data, "upstream_cost") ? Float64(data["upstream_cost"]) : nothing,
    run_metadata=metadata)

println("Wrote four-level and viewer outputs to: $(result.output_dir)")
