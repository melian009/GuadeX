# =============================================================================
# E7 — Global sensitivity toolkit (Latin hypercube + Sobol indices).
#
# This file adds a *global* sensitivity design alongside the existing
# within-design factor ranking (`scripts/plot_sensitivity_effects.jl`,
# `parameter_effect_rank.csv`, `variance_decomposition.csv`).  The two answer
# different questions and are deliberately kept side by side:
#
#   * within-design ranking: "which swept factor moves the response most inside
#     one balanced full-factorial sweep?" — valid only for that design and
#     relative to its chosen levels.
#   * this global design: "how much of the response variance is attributable to
#     each parameter (and its interactions) over the agreed *ranges*?" — a
#     variance decomposition of a continuous parameter space (Sobol indices),
#     with bootstrap confidence intervals.
#
# It is written in pure Julia with no new dependencies: the deterministic
# PRNG (SplitMix64 + Fisher–Yates) and the Saltelli A/B/AB/BA design are
# implemented here.  The analytic Ishigami function is used in the unit tests
# to validate the estimators against known S1/ST.
#
# Nothing here changes the model equations, the run offsets or any existing
# default; it only samples parameters that the model already exposes and calls
# the existing ODE/set_* helpers.
# =============================================================================

# ---------------------------------------------------------------------------
# Deterministic PRNG (SplitMix64) and helpers
# ---------------------------------------------------------------------------

"""
    GlobalSensitivityRNG(seed)

Small, self-contained deterministic 64-bit generator (SplitMix64) used by the
Latin-hypercube sampler and the bootstrap.  Implemented locally so the toolkit
needs no `Random` dependency.  The same seed always yields the same stream, so
every Sobol result is reproducible.
"""
mutable struct GlobalSensitivityRNG
    state::UInt64
end

GlobalSensitivityRNG(seed::Integer) =
    GlobalSensitivityRNG(reinterpret(UInt64, Int64(seed) % Int64))

const _SPLITMIX_GOLDEN = 0x9e3779b97f4a7c15

function _next_u64(rng::GlobalSensitivityRNG)
    rng.state += _SPLITMIX_GOLDEN
    z = rng.state
    z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9
    z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
    return z ⊻ (z >> 31)
end

# Uniform in [0, 1) with 53 random bits.
function _next_float(rng::GlobalSensitivityRNG)
    return (_next_u64(rng) >> 11) * (1.0 / 9007199254740992.0)
end

_rand_index(rng::GlobalSensitivityRNG, k::Integer) = Int(_next_u64(rng) % UInt64(k)) + 1

"""
    randperm(rng, n)

Fisher–Yates permutation of `1:n` drawn from `rng`.  Used by the Latin-hypercube
sampler, where one sample per stratum requires a permutation (no replacement).
"""
function randperm(rng::GlobalSensitivityRNG, n::Integer)
    perm = collect(1:Int(n))
    for i in Int(n):-1:2
        j = _rand_index(rng, i)
        perm[i], perm[j] = perm[j], perm[i]
    end
    return perm
end

"""
    bootstrap_indices(rng, n)

Bootstrap row indices: `n` independent draws from `1:n` **with replacement**.
This is deliberately not a permutation — a permutation would keep every
bootstrap sample's mean unchanged and collapse the confidence intervals.
"""
function bootstrap_indices(rng::GlobalSensitivityRNG, n::Integer)
    n >= 1 || error("n must be >= 1 (got $n)")
    return [Int(_next_u64(rng) % UInt64(n)) + 1 for _ in 1:Int(n)]
end

_resolve_rng(rng::Union{Nothing,GlobalSensitivityRNG}, seed) =
    rng === nothing ? GlobalSensitivityRNG(seed === nothing ? 0 : seed) : rng

# ---------------------------------------------------------------------------
# Latin hypercube sampler
# ---------------------------------------------------------------------------

"""
    latin_hypercube(n, d; seed=nothing, rng=nothing) -> Matrix{Float64}

`n × d` Latin-hypercube sample in the unit cube `[0, 1)^d`.  Each of the `d`
columns is a permutation of the `n` equiprobable strata, with one uniform draw
inside each stratum; the columns are independently permuted, so the design is
marginally uniform over `[0, 1)` by construction and has one sample per stratum
in every dimension.

Pass a fixed `seed` (or an explicit `rng`) for a deterministic design.
"""
function latin_hypercube(n::Integer, d::Integer; seed=nothing, rng=nothing)
    n >= 1 || error("n must be >= 1 (got $n)")
    d >= 1 || error("d must be >= 1 (got $d)")
    generator = _resolve_rng(rng, seed)
    U = Matrix{Float64}(undef, Int(n), Int(d))
    for j in 1:d
        perm = randperm(generator, n)
        for i in 1:n
            U[i, j] = (perm[i] - 1 + _next_float(generator)) / n
        end
    end
    return U
end

"""
    lhs_samples(n, lower, upper; seed=nothing, rng=nothing) -> Matrix{Float64}

Latin-hypercube sample mapped onto the parameter box `[lower, upper]` (both
length-`d` vectors).  Equivalent to `latin_hypercube(n, d)` then rescaling each
column to its range.
"""
function lhs_samples(n::Integer, lower::AbstractVector, upper::AbstractVector;
        seed=nothing, rng=nothing)
    length(lower) == length(upper) ||
        error("lower and upper must have the same length")
    return rescale_columns(latin_hypercube(n, length(lower); seed=seed, rng=rng),
        lower, upper)
end

"""
    rescale_columns(U, lower, upper)

Map every column of the unit-cube matrix `U` linearly onto its `[lower_j,
upper_j]` range.
"""
function rescale_columns(U::AbstractMatrix, lower::AbstractVector, upper::AbstractVector)
    size(U, 2) == length(lower) == length(upper) ||
        error("U must have one column per parameter (got $(size(U, 2)) columns, " *
              "$(length(lower)) bounds)")
    out = Matrix{Float64}(undef, size(U, 1), size(U, 2))
    for j in 1:length(lower)
        lo = Float64(lower[j]); hi = Float64(upper[j])
        hi >= lo || error("upper[$j] < lower[$j]")
        for i in 1:size(U, 1)
            out[i, j] = lo + Float64(U[i, j]) * (hi - lo)
        end
    end
    return out
end

"""
    rescale_row(u, lower, upper)

Map a single unit-cube row onto `[lower, upper]`.
"""
function rescale_row(u::AbstractVector, lower::AbstractVector, upper::AbstractVector)
    length(u) == length(lower) == length(upper) ||
        error("row and bounds must have the same length")
    return [Float64(lower[j]) + Float64(u[j]) * (Float64(upper[j]) - Float64(lower[j]))
            for j in eachindex(u)]
end

# ---------------------------------------------------------------------------
# Saltelli A/B/AB/BA design
# ---------------------------------------------------------------------------

"""
    SaltelliDesign

Saltelli et al. (2010) design for first-order (`S1`) and total-order (`ST`)
Sobol indices over `d` parameters and `n` base samples.  It holds the unit-cube
matrices `A`, `B` and, for every parameter `i`, `AB[i]` (A with column `i` taken
from B) and `BA[i]` (B with column `i` taken from A).  The full design therefore
requires `n * (2d + 2)` model evaluations.

Use [`rescale_columns`](@ref) (or [`saltelli_design`](@ref)'s `lower`/`upper`
arguments) to map it onto a parameter box.
"""
struct SaltelliDesign
    A::Matrix{Float64}
    B::Matrix{Float64}
    AB::Vector{Matrix{Float64}}
    BA::Vector{Matrix{Float64}}
end

"""
    saltelli_design(n, d; seed=nothing, rng=nothing) -> SaltelliDesign

Build the Saltelli A/B/AB/BA design in the unit cube.  `A` and `B` are
independent uniform samples; `AB[i]` replaces column `i` of `A` with column `i`
of `B`, and `BA[i]` the reverse.  Total model runs: `n * (2d + 2)`.
"""
function saltelli_design(n::Integer, d::Integer; seed=nothing, rng=nothing)
    n >= 2 || error("n must be >= 2 for a meaningful variance estimate (got $n)")
    d >= 1 || error("d must be >= 1 (got $d)")
    generator = _resolve_rng(rng, seed)
    A = _uniform_matrix(generator, n, d)
    B = _uniform_matrix(generator, n, d)
    AB = Matrix{Float64}[copy(A) for _ in 1:d]
    BA = Matrix{Float64}[copy(B) for _ in 1:d]
    for i in 1:d
        AB[i][:, i] .= view(B, :, i)
        BA[i][:, i] .= view(A, :, i)
    end
    return SaltelliDesign(A, B, AB, BA)
end

function _uniform_matrix(generator::GlobalSensitivityRNG, n, d)
    U = Matrix{Float64}(undef, Int(n), Int(d))
    for j in 1:d, i in 1:n
        U[i, j] = _next_float(generator)
    end
    return U
end

"""
    n_model_runs(design) -> Int

Number of model evaluations the Saltelli design requires: `n * (2d + 2)`.
"""
n_model_runs(design::SaltelliDesign) = size(design.A, 1) * (2 * size(design.A, 2) + 2)

"""
    saltelli_sample_matrix(design; lower=nothing, upper=nothing, seed=nothing, rng=nothing)
        -> Matrix{Float64}

All `n * (2d + 2)` sample rows stacked in a fixed order (`A`, `B`, `AB_1…AB_d`,
`BA_1…BA_d`).  When `lower`/`upper` are supplied the rows are mapped onto the
parameter box.  Useful for archiving the exact evaluated design.
"""
function saltelli_sample_matrix(design::SaltelliDesign;
        lower::Union{Nothing,AbstractVector}=nothing,
        upper::Union{Nothing,AbstractVector}=nothing,
        seed=nothing, rng=nothing)
    d = size(design.A, 2)
    rows = Matrix{Float64}[design.A, design.B]
    append!(rows, design.AB)
    append!(rows, design.BA)
    stacked = reduce(vcat, rows)
    (lower === nothing && upper === nothing) && return stacked
    (lower === nothing || upper === nothing) &&
        error("supply both lower and upper, or neither")
    return reduce(vcat, [rescale_columns(m, lower, upper) for m in rows])
end

# ---------------------------------------------------------------------------
# Sobol estimators (Saltelli/Jansen) and bootstrap confidence intervals
# ---------------------------------------------------------------------------

# Percentile of the finite bootstrap replicates.  A metric that is constant
# across the design (zero variance) yields NaN indices; those are excluded so
# the confidence interval still works for the other metrics.
function _finite_quantile(values, p::Real)
    finite = Float64[v for v in values if isfinite(v)]
    isempty(finite) && return NaN
    return quantile(finite, p)
end

"""
    sobol_indices_from_outputs(design, Y_A, Y_B, Y_AB, Y_BA;
        bootstrap=1000, confidence=0.95, seed=nothing, rng=nothing,
        parameter_names=nothing, metric_names=nothing)

Estimate first-order (`S1`) and total-order (`ST`) Sobol indices from already
evaluated outputs.

`Y_A` and `Y_B` are `n × m` matrices (rows = base samples, columns = response
metrics); `Y_AB` and `Y_BA` are length-`d` vectors of `n × m` matrices.

Estimators (Saltelli et al. 2010; Jansen 1999), averaged over the AB and BA
orientations to reduce bias:

    S1_i = mean( Y_B[:,k] .* (Y_AB[i][:,k] .- Y_A[:,k]) ) / V_k
    ST_i = mean( (Y_A[:,k] .- Y_AB[i][:,k]).^2 ) / (2 V_k)

with the symmetric BA forms averaged in, and `V_k` the variance of the pooled
A∪B outputs for metric `k`.  Percentile bootstrap confidence intervals are
computed by resampling the `n` rows with replacement `bootstrap` times.

Returns a NamedTuple with `S1`, `ST` (`d × m`), `S1_ci`, `ST_ci`
(`d × m × 2`), `variance` (length `m`), `n`, `d`, `n_model_runs` and the
optional names.
"""
function sobol_indices_from_outputs(design::SaltelliDesign,
        Y_A::AbstractMatrix, Y_B::AbstractMatrix,
        Y_AB::AbstractVector, Y_BA::AbstractVector;
        bootstrap::Integer=1000, confidence::Real=0.95,
        seed=nothing, rng=nothing,
        parameter_names::Union{Nothing,AbstractVector}=nothing,
        metric_names::Union{Nothing,AbstractVector}=nothing)
    n = size(design.A, 1)
    d = size(design.A, 2)
    size(Y_A) == size(Y_B) || error("Y_A and Y_B must have the same size")
    size(Y_A, 1) == n || error("Y_A must have n rows")
    length(Y_AB) == d || error("Y_AB must have one matrix per parameter")
    length(Y_BA) == d || error("Y_BA must have one matrix per parameter")
    m = size(Y_A, 2)
    all(size(Y_AB[i]) == (n, m) for i in 1:d) ||
        error("every Y_AB[i] must be n × m")
    all(size(Y_BA[i]) == (n, m) for i in 1:d) ||
        error("every Y_BA[i] must be n × m")
    0.0 < confidence < 1.0 || error("confidence must be in (0, 1)")

    pooled = vcat(Y_A, Y_B)
    V = [var(view(pooled, :, k)) for k in 1:m]

    S1 = Matrix{Float64}(undef, d, m)
    ST = Matrix{Float64}(undef, d, m)
    _sobol_core!(S1, ST, V, Y_A, Y_B, Y_AB, Y_BA)

    generator = _resolve_rng(rng, seed)
    S1_boot = Array{Float64}(undef, d, m, bootstrap)
    ST_boot = Array{Float64}(undef, d, m, bootstrap)
    for b in 1:bootstrap
        idx = bootstrap_indices(generator, n)
        _sobol_core!(view(S1_boot, :, :, b), view(ST_boot, :, :, b), V,
            Y_A[idx, :], Y_B[idx, :],
            [Y_AB[i][idx, :] for i in 1:d], [Y_BA[i][idx, :] for i in 1:d])
    end

    α = (1 - Float64(confidence)) / 2
    S1_ci = Array{Float64}(undef, d, m, 2)
    ST_ci = Array{Float64}(undef, d, m, 2)
    for i in 1:d, k in 1:m
        S1_ci[i, k, 1] = _finite_quantile(view(S1_boot, i, k, :), α)
        S1_ci[i, k, 2] = _finite_quantile(view(S1_boot, i, k, :), 1 - α)
        ST_ci[i, k, 1] = _finite_quantile(view(ST_boot, i, k, :), α)
        ST_ci[i, k, 2] = _finite_quantile(view(ST_boot, i, k, :), 1 - α)
    end

    return (
        S1 = S1, ST = ST, S1_ci = S1_ci, ST_ci = ST_ci,
        variance = V, n = n, d = d, n_model_runs = n_model_runs(design),
        parameter_names = parameter_names === nothing ? nothing : collect(parameter_names),
        metric_names = metric_names === nothing ? nothing : collect(metric_names),
        confidence = Float64(confidence), bootstrap = Int(bootstrap),
    )
end

# Core estimators writing into preallocated S1/ST (used for the point estimate
# and for every bootstrap replicate).
function _sobol_core!(S1, ST, V, Y_A, Y_B, Y_AB, Y_BA)
    d = length(Y_AB)
    for i in 1:d
        for k in 1:size(Y_A, 2)
            vk = V[k]
            if !(isfinite(vk) && vk > 0)
                S1[i, k] = NaN
                ST[i, k] = NaN
                continue
            end
            s1_ab = sum(Y_B[:, k] .* (Y_AB[i][:, k] .- Y_A[:, k])) / size(Y_A, 1)
            s1_ba = sum(Y_A[:, k] .* (Y_BA[i][:, k] .- Y_B[:, k])) / size(Y_A, 1)
            S1[i, k] = 0.5 * (s1_ab + s1_ba) / vk

            st_ab = sum((Y_A[:, k] .- Y_AB[i][:, k]).^2) / (2 * size(Y_A, 1))
            st_ba = sum((Y_B[:, k] .- Y_BA[i][:, k]).^2) / (2 * size(Y_A, 1))
            ST[i, k] = 0.5 * (st_ab + st_ba) / vk
        end
    end
    return nothing
end

"""
    sobol_analyze(evaluate, design; lower, upper, bootstrap=1000,
                  confidence=0.95, seed=nothing, rng=nothing,
                  parameter_names=nothing, metric_names=nothing)

Driver: evaluate `evaluate(x)` for every sample row of `design` (mapped onto the
parameter box `[lower, upper]`) and return the Sobol analysis.

`evaluate` receives a length-`d` parameter vector and must return a scalar or a
vector of response metrics.  A scalar is reported as a single metric; a vector
as `m` metrics sharing one index table per metric.
"""
function sobol_analyze(evaluate, design::SaltelliDesign;
        lower::AbstractVector, upper::AbstractVector,
        bootstrap::Integer=1000, confidence::Real=0.95,
        seed=nothing, rng=nothing,
        parameter_names::Union{Nothing,AbstractVector}=nothing,
        metric_names::Union{Nothing,AbstractVector}=nothing)
    d = size(design.A, 2)
    length(lower) == length(upper) == d ||
        error("lower/upper must have one bound per parameter")

    Y_A = _evaluate_matrix(evaluate, design.A, lower, upper)
    Y_B = _evaluate_matrix(evaluate, design.B, lower, upper)
    Y_AB = [_evaluate_matrix(evaluate, design.AB[i], lower, upper) for i in 1:d]
    Y_BA = [_evaluate_matrix(evaluate, design.BA[i], lower, upper) for i in 1:d]

    return sobol_indices_from_outputs(design, Y_A, Y_B, Y_AB, Y_BA;
        bootstrap=bootstrap, confidence=confidence, seed=seed, rng=rng,
        parameter_names=parameter_names, metric_names=metric_names)
end

function _evaluate_matrix(evaluate, U::AbstractMatrix, lower, upper)
    values = [evaluate(rescale_row(view(U, i, :), lower, upper)) for i in 1:size(U, 1)]
    m = length(values[1])
    Y = Matrix{Float64}(undef, size(U, 1), m)
    for i in eachindex(values)
        length(values[i]) == m || error("evaluate returned inconsistent lengths")
        Y[i, :] .= Float64.(values[i])
    end
    return Y
end

# ---------------------------------------------------------------------------
# Parameter transforms used by the driver (model-equivalent, no equation change)
# ---------------------------------------------------------------------------

"""
    scale_intrinsic_growth_rate(params, factor)

Copy of `params` with every local intrinsic growth rate multiplied by `factor`
(the global growth multiplier).  `r` enters the ODE linearly inside the growth
bracket, so this is a pure parameter change.
"""
function scale_intrinsic_growth_rate(params::MetacommunityParams, factor::Real)
    factor > 0 || error("growth multiplier must be > 0 (got $factor)")
    return MetacommunityParams(
        params.n_sites, params.n_species, params.interaction_matrix,
        params.dispersal_matrix, params.dispersal_scaling,
        params.intrinsic_growth_rates .* Float64(factor), params.temperatures,
        params.habitat_suitability, params.thermal_optima, params.thermal_sigmas,
        params.carrying_capacity, params.thermal_lower_limits,
        params.thermal_upper_limits, params.heat_stress_rate,
        params.interaction_inside_growth)
end

"""
    scale_dispersal(params, factor)

Copy of `params` with the precomputed dispersal matrix multiplied by `factor`
(the dispersal-rate multiplier).  Scales emigration and immigration together,
so it changes turnover speed but not the network topology.
"""
function scale_dispersal(params::MetacommunityParams, factor::Real)
    factor > 0 || error("dispersal multiplier must be > 0 (got $factor)")
    return set_dispersal_matrix(params, params.dispersal_matrix .* Float64(factor))
end

# ---------------------------------------------------------------------------
# Basin response metrics from a state
# ---------------------------------------------------------------------------

"""
    basin_response_metrics(final_state, baseline_state, species; native, threshold=0.1)

Basin-level response summarised directly from two flattened states (no file
round-trip), for use inside a Sobol driver.  Returns a NamedTuple with:

- `native_biomass`: site-mean native biomass in the final state,
- `native_richness`: site-mean number of native species above `threshold` in the
  final state,
- `realised_richness_loss`: `max(0, 1 - final_richness / baseline_richness)`,
  the deterministic richness contraction against the supplied baseline (E8),
  clipped at 0.

`species` is the model species order; `native` lists the native codes (a
non-native code in `native` is ignored).  The baseline is required so the loss
metric is defined consistently with `compute_site_metrics`.
"""
function basin_response_metrics(final_state::AbstractVector, baseline_state::AbstractVector,
        species::AbstractVector; native::AbstractVector, threshold::Real=0.1)
    length(species) > 0 || error("species must not be empty")
    length(final_state) == length(baseline_state) ||
        error("final_state and baseline_state must have the same length")
    n_species = length(species)
    length(final_state) % n_species == 0 ||
        error("state length is not a multiple of the species count")
    n_sites = length(final_state) ÷ n_species

    native_idx = [i for i in 1:n_species if string(species[i]) in string.(native)]
    isempty(native_idx) && error("none of the native codes are in the species list")

    Uf = reshape(final_state, n_sites, n_species)
    Ub = reshape(baseline_state, n_sites, n_species)
    final_biomass = [sum(Uf[i, native_idx]) for i in 1:n_sites]
    final_richness = [count(>(Float64(threshold)), view(Uf, i, native_idx)) for i in 1:n_sites]
    base_richness = [count(>(Float64(threshold)), view(Ub, i, native_idx)) for i in 1:n_sites]

    mean_final_richness = mean(final_richness)
    mean_base_richness = mean(base_richness)
    loss = mean_base_richness > 0 ?
        max(0.0, 1.0 - mean_final_richness / mean_base_richness) : 0.0
    return (
        native_biomass = mean(final_biomass),
        native_richness = mean_final_richness,
        realised_richness_loss = loss,
    )
end
