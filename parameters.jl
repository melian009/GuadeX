# =============================================================================
# SimulationParameters: loader for the GuadeX run settings file.
#
# Each maintained entry script includes this file and then reads values from
# `GUADEX_PARAMS` (a raw TOML Dict).  The parameter file defaults to
# `parameters.toml` next to this file and can be overridden with the
# GUADEX_PARAMETERS_FILE environment variable.
#
# Library defaults in src/data_preparation.jl are intentionally unchanged.
# =============================================================================

module SimulationParameters

using TOML

export load, float_vector, scenario_dict, scenario_library, require_sections

"""
    default_file()

Return the parameter file path: the GUADEX_PARAMETERS_FILE environment
variable if set, otherwise `parameters.toml` in the project root.
"""
function default_file()
    return get(ENV, "GUADEX_PARAMETERS_FILE", joinpath(@__DIR__, "parameters.toml"))
end

"""
    load(path=nothing)

Parse the TOML parameter file and return it as a nested Dict.  Prints a
provenance line with the resolved file path.
"""
function load(path::Union{Nothing,AbstractString}=nothing)
    file = path === nothing ? default_file() : String(path)
    isfile(file) || error("parameter file not found: $file")
    raw = TOML.parsefile(file)
    println("[SimulationParameters] using parameter file: $file")
    return raw
end

"""
    require_sections(raw, sections...)

Error unless every requested top-level section exists.
"""
function require_sections(raw::AbstractDict, sections::AbstractString...)
    missing_sections = [section for section in sections if !haskey(raw, section)]
    isempty(missing_sections) ||
        error("parameter file is missing sections: $(join(missing_sections, ", "))")
    return nothing
end

"""
    float_vector(values)

Return `values` as a `Vector{Float64}`.
"""
function float_vector(values)
    return Float64.(values)
end

"""
    expand_scenario(subcatchments, spec)

Build a `Dict{Float64,Float64}` from a scenario spec: `factor` applied to every
subcatchment, overridden per subcatchment by the optional `overrides` table
(keys are stringified subcatchment ids).
"""
function expand_scenario(subcatchments::Vector{Float64}, spec::AbstractDict)
    haskey(spec, "factor") || error("scenario spec requires a 'factor': $(keys(spec))")
    factor = Float64(spec["factor"])
    scenario = Dict{Float64,Float64}(subcatchment => factor for subcatchment in subcatchments)
    if haskey(spec, "overrides")
        for (subcatchment_key, value) in spec["overrides"]
            scenario[parse(Float64, string(subcatchment_key))] = Float64(value)
        end
    end
    return scenario
end

"""
    scenario_dict(subcatchments, spec)

Return a single expanded scenario as a `Dict{Float64,Float64}`.  Used by
entry scripts that apply exactly one scenario (e.g. `run_model.jl`).
"""
function scenario_dict(subcatchments, spec::AbstractDict)
    return expand_scenario(float_vector(subcatchments), spec)
end

"""
    scenario_library(subcatchments, library, names)

Expand the named scenarios from `library` (the parsed `[scenarios.<kind>]`
section) into `Dict{String,Dict{Float64,Float64}}`, in the order given by
`names`.  Errors on unknown scenario names so typos fail loudly.
"""
function scenario_library(subcatchments, library::AbstractDict, names)
    all_subcatchments = float_vector(subcatchments)
    scenarios = Dict{String,Dict{Float64,Float64}}()
    for name in names
        key = string(name)
        haskey(library, key) || error("unknown scenario '$key' in the scenario library")
        scenarios[key] = expand_scenario(all_subcatchments, library[key])
    end
    return scenarios
end

end # module SimulationParameters
