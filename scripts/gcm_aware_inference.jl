using Pkg; Pkg.activate(".")
using CSV
using DataFrames
using Statistics
using LinearAlgebra
using Guadex

# GCM-aware inference for the aggregate dose-response slopes (Table 5 of
# docs/FinalReportOctober2026.tex).
#
# Pooled OLS and the GCM fixed-effect fit come from Guadex.dose_response.  With
# MixedModels deliberately absent, the standard error is hand-rolled as a
# GCM-clustered (sandwich) SE on y ~ x, and the two-sided p value uses a
# Student-t reference with (number of GCM clusters - 1) degrees of freedom.
#
#   julia --project=. scripts/gcm_aware_inference.jl
#
# Override the root with GUADEX_CLIMATE_OUTPUT_DIR (default:
# results/climate_scenarios_phase2).

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "scripts", "plot_paths.jl"))

const CLIMATE_ROOT = report_output_root("GUADEX_CLIMATE_OUTPUT_DIR",
    joinpath("results", "climate_scenarios_phase2"); base=ROOT)

# --- Student-t two-sided p value (regularized incomplete beta) ---------------

# Lanczos log-gamma (SpecialFunctions is not a project dependency).
const _LANCZOS = [676.5203681218851, -1259.1392167224028,
    771.32342877765313, -176.61502916214059, 12.507343278686905,
    -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7]

function _loggamma(z)
    z < 0.5 && return log(pi / sin(pi * z)) - _loggamma(1.0 - z)
    z -= 1.0
    x = 0.99999999999980993
    for (i, g) in enumerate(_LANCZOS)
        x += g / (z + i)
    end
    t = z + length(_LANCZOS) - 0.5
    return 0.5 * log(2pi) + (z + 0.5) * log(t) - t + log(x)
end

function _betacf(a, b, x)
    MAXIT, EPS, FPMIN = 200, 3.0e-14, 1.0e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    abs(d) < FPMIN && (d = FPMIN)
    d = 1.0 / d
    h = d
    for m in 1:MAXIT
        m2 = 2m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d; abs(d) < FPMIN && (d = FPMIN)
        c = 1.0 + aa / c; abs(c) < FPMIN && (c = FPMIN)
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d; abs(d) < FPMIN && (d = FPMIN)
        c = 1.0 + aa / c; abs(c) < FPMIN && (c = FPMIN)
        d = 1.0 / d
        del = d * c
        h *= del
        abs(del - 1.0) < EPS && break
    end
    return h
end

function _betai(a, b, x)
    (x <= 0.0 || x >= 1.0) && return x <= 0.0 ? 0.0 : 1.0
    bt = exp(_loggamma(a + b) - _loggamma(a) - _loggamma(b) +
             a * log(x) + b * log1p(-x))
    return x < (a + 1.0) / (a + b + 2.0) ?
           bt * _betacf(a, b, x) / a :
           1.0 - bt * _betacf(b, a, 1.0 - x) / b
end

t_two_sided_p(t, df) = _betai(df / 2, 0.5, df / (df + t^2))

# --- GCM-clustered sandwich SE for y ~ 1 + x --------------------------------

function cluster_robust_slope(x, y, groups)
    n = length(x)
    X = hcat(ones(n), x)
    XtX_inv = inv(X' * X)
    beta = XtX_inv * (X' * y)
    resid = y - X * beta
    labels = sort(unique(groups))
    meat = zeros(2, 2)
    for g in labels
        idx = findall(==(g), groups)
        s = X[idx, :]' * resid[idx]
        meat .+= s * s'
    end
    G = length(labels)
    scale = (G / (G - 1)) * ((n - 1) / (n - 2))
    V = XtX_inv * meat * XtX_inv * scale
    se = sqrt(V[2, 2])
    t = beta[2] / se
    return (slope=beta[2], se=se, t=t, p=t_two_sided_p(t, G - 1), groups=G, n=n)
end

# GCM fixed-effect design y ~ x + GCM, with a GCM-clustered sandwich SE.
function cluster_robust_fixed_effect(x, y, groups)
    n = length(x)
    labels = sort(unique(groups))
    k = 2 + length(labels) - 1
    X = zeros(n, k)
    X[:, 1] .= 1.0
    X[:, 2] .= x
    for (j, g) in enumerate(labels)
        j == 1 && continue
        X[:, 1 + j] .= Float64.(groups .== g)
    end
    XtX_inv = inv(X' * X)
    beta = XtX_inv * (X' * y)
    resid = y - X * beta
    meat = zeros(k, k)
    for g in labels
        idx = findall(==(g), groups)
        s = X[idx, :]' * resid[idx]
        meat .+= s * s'
    end
    G = length(labels)
    scale = (G / (G - 1)) * ((n - 1) / (n - k))
    V = XtX_inv * meat * XtX_inv * scale
    se = sqrt(V[2, 2])
    t = beta[2] / se
    return (slope=beta[2], se=se, t=t, p=t_two_sided_p(t, G - 1), groups=G, n=n)
end

function main()
    index = CSV.read(joinpath(CLIMATE_ROOT, "runs_index.csv"), DataFrame)
    scenario = index[index.scenario .!= "control", :]
    x = Float64.(scenario.realised_warming_2026_2045_mean_degc)
    groups = String.(scenario.gcm)
    metrics = [
        ("total_density", "basin_total_biomass"),
        ("native_richness", "basin_native_richness"),
        ("realised_richness_loss", "basin_realised_richness_loss"),
    ]
    out = DataFrame(metric=String[], pooled_slope=Float64[], pooled_r2=Float64[],
        pooled_cluster_se=Float64[], pooled_p=Float64[],
        fixed_effect_slope=Float64[], fixed_effect_r2=Float64[],
        fe_cluster_se=Float64[], fe_p=Float64[],
        per_gcm_sd=Float64[], between_gcm_sd=Float64[], n=Int[], n_groups=Int[])
    for (name, col) in metrics
        y = Float64.(scenario[!, col])
        d = Guadex.dose_response(x, y, groups)
        cp = cluster_robust_slope(x, y, groups)
        cf = cluster_robust_fixed_effect(x, y, groups)
        push!(out, (name, d.slope, d.r2, cp.se, cp.p,
            d.fixed_effect_slope, d.fixed_effect_r2, cf.se, cf.p,
            d.slope_spread, d.between_group_sd, d.n, d.n_groups))
        println(rpad(name, 26),
            " pooled=$(round(d.slope, sigdigits=3)) (p=$(round(cp.p, digits=4)))",
            " FE=$(round(d.fixed_effect_slope, sigdigits=3)) (p=$(round(cf.p, digits=4)))",
            " R2=$(round(d.r2, digits=3))/$(round(d.fixed_effect_r2, digits=3))",
            " per_gcm_sd=$(round(d.slope_spread, sigdigits=3))",
            " between_gcm_sd=$(round(d.between_group_sd, sigdigits=3))")
    end
    path = joinpath(CLIMATE_ROOT, "gcm_aware_inference.csv")
    CSV.write(path, out)
    println("wrote $path")
end

main()
