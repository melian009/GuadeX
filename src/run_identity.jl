# =============================================================================
# Run identity: stable digests for reproducible, auditable runs (E9/E10).
#
# A run's results depend on the resolved model state (interaction matrix,
# growth/temperature/carrying-capacity/thermal arrays, heat-stress slope,
# dispersal scaling, the interaction-form flag), on the resolved run options
# (forcing mode, baseline period, burn-in controls, K scaling, upstream cost,
# obstacle/passability configuration) and on the code version.  All of these are
# folded into one deterministic 16-hex-character FNV-1a digest by
# [`parameter_digest`](@ref).  The digest is included in the burn-in cache key
# and the run fingerprints, and is written into `run_metadata.json` and a
# `run_digest.txt` sidecar so a stale output can never be silently reused.
# =============================================================================

const _CODE_VERSION_CACHE = Ref{Union{Nothing,String}}(nothing)

"""
    stable_digest(s)

Deterministic 16-hex-character FNV-1a digest.  Unlike `Base.hash`, it is stable
across sessions and Julia versions, which is what makes it usable for cache keys
and run fingerprints.
"""
function stable_digest(s::AbstractString)
    h = 0xcbf29ce484222325
    for b in codeunits(s)
        h = (h ⊻ UInt64(b)) * 0x00000100000001b3
    end
    return string(h, base=16, pad=16)
end

"""
    code_version()

Best-effort short git commit (`git rev-parse --short HEAD`) of the checkout
running the model, memoized per session.  Returns `"unknown"` when git is
unavailable or the command fails, so the digest is always computable.
"""
function code_version()
    if _CODE_VERSION_CACHE[] === nothing
        _CODE_VERSION_CACHE[] = try
            strip(read(`git rev-parse --short HEAD`, String))
        catch
            "unknown"
        end
    end
    return _CODE_VERSION_CACHE[]
end

"""
    file_fingerprint(path)

Size/mtime fingerprint of an input file (`"missing"` when absent).  Used to
invalidate a cached result when the forcing or projection file changes without
changing its name.
"""
function file_fingerprint(path::AbstractString)
    isfile(path) || return "missing"
    info = stat(path)
    return "sz$(info.size)_mt$(round(Int, info.mtime))"
end

# Digest of one array's contents, preserving the stored order (deterministic
# because the model arrays themselves are built in a deterministic order).
function _array_digest(values)
    io = IOBuffer()
    print(io, size(values), ":")
    for value in vec(values)
        print(io, repr(value), ",")
    end
    return stable_digest(String(take!(io)))
end

# Digest of one resolved run-option value.  Arrays/dicts are handled so a nested
# config (e.g. a list of climate models) serializes deterministically.
function _value_digest(value)
    if value isa AbstractArray
        return _array_digest(value)
    elseif value isa AbstractDict
        keys_sorted = sort(collect(keys(value)), by=string)
        return "{" * join(("$(key)=$(_value_digest(value[key]))" for key in keys_sorted), ";") * "}"
    elseif value isa Symbol || value isa AbstractString
        return string(value)
    elseif value isa Bool
        return value ? "true" : "false"
    elseif value isa Integer
        return string(value)
    elseif value isa Real
        return repr(Float64(value))
    end
    return string(value)
end

"""
    parameter_digest(params; options=Dict(), code=code_version())

Stable digest over the resolved model state that affects results — the
interaction matrix, dispersal matrix and scaling, intrinsic growth rates,
temperatures, habitat suitability, carrying capacity, thermal optima/sigmas and
lower/upper limits, the heat-stress rate and the `interaction_inside_growth`
flag — together with the resolved run `options` (forcing mode, baseline period,
burn-in years/tolerance/criterion, carrying-capacity scaling, upstream cost,
passability/obstacle configuration, ...) and the code version.  Option keys are
sorted, so the digest is deterministic and independent of dict insertion order.
"""
function parameter_digest(params::MetacommunityParams;
        options::AbstractDict=Dict{String,Any}(),
        code::AbstractString=code_version())
    pieces = String[
        "n_sites=$(params.n_sites)",
        "n_species=$(params.n_species)",
        "interaction_matrix=$(_array_digest(params.interaction_matrix))",
        "dispersal_matrix=$(_array_digest(params.dispersal_matrix))",
        "dispersal_scaling=$(_array_digest(params.dispersal_scaling))",
        "intrinsic_growth_rates=$(_array_digest(params.intrinsic_growth_rates))",
        "temperatures=$(_array_digest(params.temperatures))",
        "habitat_suitability=$(_array_digest(params.habitat_suitability))",
        "thermal_optima=$(_array_digest(params.thermal_optima))",
        "thermal_sigmas=$(_array_digest(params.thermal_sigmas))",
        "carrying_capacity=$(_array_digest(params.carrying_capacity))",
        "thermal_lower_limits=$(_array_digest(params.thermal_lower_limits))",
        "thermal_upper_limits=$(_array_digest(params.thermal_upper_limits))",
        "heat_stress_rate=$(_value_digest(params.heat_stress_rate))",
        "interaction_inside_growth=$(params.interaction_inside_growth)",
    ]
    for key in sort(collect(keys(options)), by=string)
        push!(pieces, "option.$key=$(_value_digest(options[key]))")
    end
    push!(pieces, "code_version=$code")
    return stable_digest(join(pieces, "|"))
end
