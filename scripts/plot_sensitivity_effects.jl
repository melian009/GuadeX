using Pkg; Pkg.activate(".")
using DataFrames
using CSV
using CairoMakie
using Statistics
using Printf
using Guadex

# =============================================================================
# Parameter-sensitivity comparison for the obstacle / interaction sweeps.
#
#   julia --project=. scripts/plot_sensitivity_effects.jl
#
# Reads the already-completed (read-only) sweeps:
#   results/sensitivity_obstacles/runs_index.csv               (32 runs)
#   results/sensitivity_obstacles/alt_interactions/runs_index.csv (48 runs)
#
# IMPORTANT design caveat (stated in the report as well): the combined design
# is NOT orthogonal.  Within each sweep the grid is a balanced full factorial
# (obstacle sweep: model x upstream_cost x passability = 2x4x4 = 32;
# alt-interaction sweep: matrix x model x upstream_cost x passability =
# 3x2x2x4 = 48), so a standard sum-of-squares decomposition is valid *inside*
# a sweep.  Across the two sweeps the thermal-sigma multiplier (0.3) is fully
# confounded with the sweep AND with a shorter burn-in (21 vs 373 years), so
# the sigma effect is reported as a confounded contrast, not a clean effect.
#
# Writes:
#   results/sensitivity_obstacles/report_plots/parameter_effect_rank.csv
#   results/sensitivity_obstacles/report_plots/variance_decomposition.csv
#   results/sensitivity_obstacles/report_plots/fig_tornado.png
#   results/sensitivity_obstacles/report_plots/fig_parameter_importance_ranked.png
#   results/sensitivity_obstacles/report_plots/fig_parameter_importance_dashboard.png
#   results/sensitivity_obstacles/report_plots/fig_interaction_effects.png
# =============================================================================

const ROOT = joinpath("results", "sensitivity_obstacles")
const MAIN_INDEX = joinpath(ROOT, "runs_index.csv")
const ALT_INDEX = joinpath(ROOT, "alt_interactions", "runs_index.csv")
const REPORT_DIR = joinpath(ROOT, "report_plots")
mkpath(REPORT_DIR)

const METRICS = [
    "basin_total_biomass",
    "basin_native_biomass",
    "basin_invasive_biomass",
    "basin_native_richness",
    "basin_native_extinction_risk",
    "st_relative_biomass_change",
    "st_final_biomass",
    "st_quasi_extinct_fraction",
]

const METRIC_LABELS = Dict(
    "basin_total_biomass" => "Basin total biomass",
    "basin_native_biomass" => "Basin native biomass",
    "basin_invasive_biomass" => "Basin invasive biomass",
    "basin_native_richness" => "Basin native richness",
    "basin_native_extinction_risk" => "Native extinction risk",
    "st_relative_biomass_change" => "ST relative biomass change",
    "st_final_biomass" => "ST final biomass",
    "st_quasi_extinct_fraction" => "ST quasi-extinct fraction",
)

# Level encodings.  The first entry is the reference level for % effects.
const PASS_LEVELS = ["baseline", "improved_passability", "reduced_passability", "blocked"]
const UC_LEVELS = [0.01, 0.05, 0.1, 0.5]
const MATRIX_LEVELS = ["original", "random", "invasive_favoring"]
const SIGMA_LEVELS = [1.0, 0.3]   # reference = main sweep (multiplier 1.0)

isfile(MAIN_INDEX) || error("missing $MAIN_INDEX")
isfile(ALT_INDEX) || error("missing $ALT_INDEX")

main = CSV.read(MAIN_INDEX, DataFrame)
alt = CSV.read(ALT_INDEX, DataFrame)

# Composite model key and tidy columns.
main.model = string.(main.climate_scenario, "/", main.gcm)
alt.model = string.(alt.climate_scenario, "/", alt.gcm)
main.sigma = fill(1.0, nrow(main))
alt.sigma = Float64.(alt.thermal_sigma_multiplier)
# The obstacle sweep did not record an interaction-matrix field; it runs the
# default ("original") matrix.  Tag it so the two sweeps can be column-aligned.
main.interaction_matrix = fill("original", nrow(main))

const COOL_MODEL = "ssp126/IITM-ESM"      # coolest sensitivity end
const WARM_MODEL = "ssp245/UKESM1-0-LL"   # warmest sensitivity end

# ---------------------------------------------------------------------------
# Generic block-matched effect helpers
# ---------------------------------------------------------------------------

"""
    effect_spans(rows, metric, factor_col, levels, block_cols) -> (records, sigma)

For every matched block (all `block_cols` held fixed) that contains >=2 of the
factor `levels`, compute the span (max-min) of `metric` across those levels.
`span_std` standardises by the run-set standard deviation of the metric; the
`span_pct` additionally divides by the absolute reference-level value (first
entry of `levels`), guarded against near-zero baselines.
"""
function effect_spans(rows, metric, factor_col, levels, block_cols)
    y = Float64.(rows[!, Symbol(metric)])
    σ = std(y)
    σ = (isfinite(σ) && σ > 0) ? σ : 1.0
    table = Dict{Any,Dict{Any,Float64}}()
    for r in eachrow(rows)
        key = Tuple(r[c] for c in block_cols)
        table[key] = get!(table, key, Dict{Any,Float64}())
        table[key][r[factor_col]] = Float64(r[Symbol(metric)])
    end
    recs = NamedTuple[]
    for (key, lvmap) in table
        present = [l for l in levels if haskey(lvmap, l)]
        length(present) < 2 && continue
        vals = [lvmap[l] for l in present]
        span = maximum(vals) - minimum(vals)
        ref = lvmap[first(levels)]
        denom = abs(ref) > 1e-8 ? abs(ref) : maximum(abs.(vals))
        denom > 0 || (denom = 1.0)
        push!(recs, (
            block=key,
            span_std=span / σ,
            span_pct=100.0 * span / denom,
            signed=[lv => (lvmap[lv] - ref) / σ for lv in present if lv != first(levels)],
            ref=ref,
        ))
    end
    return recs, σ
end

medianspan(recs, field) = isempty(recs) ? NaN : median(r[field] for r in recs)
maxspan(recs, field) = isempty(recs) ? NaN : maximum(r[field] for r in recs)

# ---------------------------------------------------------------------------
# Build factor-level effect records for the ranking
# ---------------------------------------------------------------------------

rank_rows = NamedTuple[]
const PCT_CAP = 500.0
clamp_pct(x) = isnan(x) ? NaN : min(x, PCT_CAP)
push_rank!(factor, metric, recs; analysis) = begin
    isempty(recs) && return
    push!(rank_rows, (
        analysis=analysis,
        factor=factor,
        metric=metric,
        n_blocks=length(recs),
        median_span_std=medianspan(recs, :span_std),
        max_span_std=maxspan(recs, :span_std),
        median_span_pct=clamp_pct(medianspan(recs, :span_pct)),
        max_span_pct=clamp_pct(maxspan(recs, :span_pct)),
    ))
end

for metric in METRICS
    # passability (blocks: model x upstream_cost)
    recs, _ = effect_spans(main, metric, :passability_scenario, PASS_LEVELS, [:model, :upstream_cost])
    push_rank!("passability", metric, recs; analysis="obstacle")

    # upstream cost (blocks: model x passability)
    recs, _ = effect_spans(main, metric, :upstream_cost, UC_LEVELS, [:model, :passability_scenario])
    push_rank!("upstream_cost", metric, recs; analysis="obstacle")

    # climate model / temperature (blocks: upstream_cost x passability)
    recs, _ = effect_spans(main, metric, :model, [COOL_MODEL, WARM_MODEL],
        [:upstream_cost, :passability_scenario])
    push_rank!("climate_model", metric, recs; analysis="obstacle")

    # interaction matrix (blocks: model x upstream_cost x passability)
    recs, _ = effect_spans(alt, metric, :interaction_matrix, MATRIX_LEVELS,
        [:model, :upstream_cost, :passability_scenario])
    push_rank!("interaction_matrix", metric, recs; analysis="alt_interactions")

    # thermal sigma (blocks: model x upstream_cost in {0.01,0.5} x passability)
    combined = vcat(
        select(main, Not([:interaction_matrix]), :sigma),
        select(alt[alt.interaction_matrix .== "original", :], Not([:interaction_matrix]), :sigma);
        cols=:union)
    recs, _ = effect_spans(combined, metric, :sigma, SIGMA_LEVELS,
        [:model, :upstream_cost, :passability_scenario])
    push_rank!("thermal_sigma", metric, recs; analysis="sigma_confounded")

    # --- interaction effects -------------------------------------------------
    # passability x upstream cost: for each passability scenario, how much does
    # its effect (vs baseline passability) vary across the four upstream costs?
    # A non-zero spread is the interaction; a pure main effect would be constant.
    inter = Float64[]
    σm = std(Float64.(main[!, Symbol(metric)]))
    σm = (isfinite(σm) && σm > 0) ? σm : 1.0
    for model in unique(main.model)
        for p in ["blocked", "reduced_passability", "improved_passability"]
            eff_uc = Float64[]
            for uc in UC_LEVELS
                sel = main[(main.model .== model) .&
                           (abs.(main.upstream_cost .- uc) .< 1e-9) .&
                           (main.passability_scenario .== p), :]
                bas = main[(main.model .== model) .&
                           (abs.(main.upstream_cost .- uc) .< 1e-9) .&
                           (main.passability_scenario .== "baseline"), :]
                (isempty(sel) || isempty(bas)) && continue
                push!(eff_uc, Float64(sel[1, Symbol(metric)]) - Float64(bas[1, Symbol(metric)]))
            end
            length(eff_uc) >= 2 && push!(inter, maximum(eff_uc) - minimum(eff_uc))
        end
    end
    if !isempty(inter)
        push!(rank_rows, (
            analysis="obstacle_interaction", factor="passability_x_upstream_cost", metric=metric,
            n_blocks=length(inter),
            median_span_std=median(inter) / σm, max_span_std=maximum(inter) / σm,
            median_span_pct=NaN, max_span_pct=NaN))
    end

    # passability x temperature: difference of the passability effect between
    # the cool and warm climate-model sensitivity ends.
    inter2 = Float64[]
    σm = std(Float64.(main[!, Symbol(metric)])); σm = σm > 0 ? σm : 1.0
    for uc in UC_LEVELS
        for p in ["blocked", "reduced_passability", "improved_passability"]
            vals = Float64[]
            for model in [COOL_MODEL, WARM_MODEL]
                sel = main[(main.model .== model) .&
                           (abs.(main.upstream_cost .- uc) .< 1e-9) .&
                           (main.passability_scenario .== p), :]
                bas = main[(main.model .== model) .&
                           (abs.(main.upstream_cost .- uc) .< 1e-9) .&
                           (main.passability_scenario .== "baseline"), :]
                (isempty(sel) || isempty(bas)) && continue
                push!(vals, Float64(sel[1, Symbol(metric)]) - Float64(bas[1, Symbol(metric)]))
            end
            length(vals) == 2 && push!(inter2, abs(vals[2] - vals[1]))
        end
    end
    if !isempty(inter2)
        push!(rank_rows, (
            analysis="obstacle_interaction", factor="passability_x_temperature", metric=metric,
            n_blocks=length(inter2),
            median_span_std=median(inter2) / σm, max_span_std=maximum(inter2) / σm,
            median_span_pct=NaN, max_span_pct=NaN))
    end
end

rank = DataFrame(rank_rows)
CSV.write(joinpath(REPORT_DIR, "parameter_effect_rank.csv"), rank)

# Normalised within-metric importance and mean rank across metrics.
const RANK_FACTORS = ["passability", "upstream_cost", "climate_model",
                      "interaction_matrix", "thermal_sigma",
                      "passability_x_upstream_cost", "passability_x_temperature"]
norm = DataFrame(factor=String[], metric=String[], importance=Float64[], rank=Int[])
for metric in METRICS
    sub = rank[rank.metric .== metric, :]
    for f in RANK_FACTORS
        row = sub[sub.factor .== f, :]
        isempty(row) && continue
        push!(norm, (factor=f, metric=metric, importance=row.median_span_std[1], rank=0))
    end
    m = maximum(norm.importance[norm.metric .== metric])
    m = (isfinite(m) && m > 0) ? m : 1.0
    idx = findall(norm.metric .== metric)
    sort!(idx, by=i -> -norm.importance[i])
    for (r, i) in enumerate(idx); norm.rank[i] = r; end
    norm.importance[idx] .= norm.importance[idx] ./ m
end
overall = combine(groupby(norm, :factor),
    :importance => mean => :mean_importance,
    :rank => mean => :mean_rank)
sort!(overall, :mean_rank)

# ---------------------------------------------------------------------------
# Variance decomposition (within-sweep, balanced full factorial)
# ---------------------------------------------------------------------------

function ss_group(y, keys)
    gm = mean(y)
    d = Dict{Any,Vector{Float64}}()
    for (k, v) in zip(keys, y); push!(get!(d, k, Float64[]), v); end
    s = 0.0
    for (_, vs) in d; s += length(vs) * (mean(vs) - gm)^2; end
    return s
end

function decompose(rows, metric, factors, labels)
    y = Float64.(rows[!, Symbol(metric)]); y = y .- mean(y)
    total = sum(abs2, y)
    total > 0 || return NamedTuple[]
    main_ss = Dict{String,Float64}()
    for f in factors
        main_ss[f] = ss_group(y, [r[f] for r in eachrow(rows)])
    end
    out = NamedTuple[]
    for f in factors
        push!(out, (term=labels[f], ss_fraction=main_ss[f] / total))
    end
    for (i, a) in enumerate(factors), b in factors[(i+1):end]
        sab = ss_group(y, [(r[a], r[b]) for r in eachrow(rows)])
        push!(out, (term="$(labels[a]) × $(labels[b])",
            ss_fraction=(sab - main_ss[a] - main_ss[b]) / total))
    end
    explained = sum(x.ss_fraction for x in out)
    push!(out, (term="higher-order + residual", ss_fraction=1.0 - explained))
    return out
end

vdec = NamedTuple[]
const OBSTACLE_FACTORS = Dict(
    "model" => "climate model",
    "upstream_cost" => "upstream cost",
    "passability_scenario" => "passability")
const ALT_FACTORS = Dict(
    "interaction_matrix" => "interaction matrix",
    "model" => "climate model",
    "upstream_cost" => "upstream cost",
    "passability_scenario" => "passability")

for metric in METRICS
    for d in decompose(main, metric, collect(keys(OBSTACLE_FACTORS)), OBSTACLE_FACTORS)
        push!(vdec, (sweep="obstacle", metric=metric, term=d.term, ss_fraction=d.ss_fraction))
    end
    for d in decompose(alt, metric, collect(keys(ALT_FACTORS)), ALT_FACTORS)
        push!(vdec, (sweep="alt_interactions", metric=metric, term=d.term, ss_fraction=d.ss_fraction))
    end
end
CSV.write(joinpath(REPORT_DIR, "variance_decomposition.csv"), DataFrame(vdec))

# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------

function hbar!(ax, y, x0, x1; color=RGBAf(0.27, 0.45, 0.77, 0.85), height=0.62)
    xlo, xhi = min(x0, x1), max(x0, x1)
    poly!(ax, [Point2f(xlo, y - height / 2), Point2f(xhi, y - height / 2),
               Point2f(xhi, y + height / 2), Point2f(xlo, y + height / 2)], color=color)
end

const COL_PASS = RGBAf(0.27, 0.51, 0.71, 0.9)
const COL_UC = RGBAf(1.00, 0.55, 0.00, 0.9)
const COL_CLIM = RGBAf(0.70, 0.13, 0.13, 0.9)
const COL_MAT = RGBAf(0.18, 0.55, 0.34, 0.9)
const COL_SIG = RGBAf(0.50, 0.00, 0.50, 0.9)
const COL_PXU = RGBAf(0.00, 0.50, 0.50, 0.9)
const COL_PXT = RGBAf(0.86, 0.08, 0.24, 0.9)
const COL_GRAY = RGBAf(0.5, 0.5, 0.5, 0.9)

factor_color = Dict(
    "passability" => COL_PASS,
    "upstream_cost" => COL_UC,
    "climate_model" => COL_CLIM,
    "interaction_matrix" => COL_MAT,
    "thermal_sigma" => COL_SIG,
    "passability_x_upstream_cost" => COL_PXU,
    "passability_x_temperature" => COL_PXT,
)

# --- Tornado: signed median level contrasts per metric ---------------------
function tornado_panels!(ax, metric)
    rows = NamedTuple[]
    σ_main = std(Float64.(main[!, Symbol(metric)])); σ_main = σ_main > 0 ? σ_main : 1.0
    for (factor, col, levels, blocks) in [
        ("passability", :passability_scenario, PASS_LEVELS, [:model, :upstream_cost]),
        ("upstream_cost", :upstream_cost, UC_LEVELS, [:model, :passability_scenario]),
        ("climate_model", :model, [COOL_MODEL, WARM_MODEL], [:upstream_cost, :passability_scenario]),
    ]
        recs, _ = effect_spans(main, metric, col, levels, blocks)
        for lv in levels[2:end]
            vs = [v for r in recs for (l, v) in r.signed if l == lv && isfinite(v)]
            isempty(vs) && continue
            push!(rows, (factor=factor, level=string(lv), val=median(vs)))
        end
    end
    for (factor, col, levels, blocks) in [
        ("interaction_matrix", :interaction_matrix, MATRIX_LEVELS,
            [:model, :upstream_cost, :passability_scenario]),
    ]
        recs, σ_alt = effect_spans(alt, metric, col, levels, blocks)
        for lv in levels[2:end]
            vs = [v for r in recs for (l, v) in r.signed if l == lv && isfinite(v)]
            isempty(vs) && continue
            push!(rows, (factor=factor, level=string(lv), val=median(vs) * σ_alt / σ_main))
        end
    end
    # thermal sigma (confounded): combined main (1.0) vs alt original (0.3)
    combined = vcat(select(main, Not([:interaction_matrix]), :sigma),
        select(alt[alt.interaction_matrix .== "original", :], Not([:interaction_matrix]), :sigma); cols=:union)
    recs, σ_c = effect_spans(combined, metric, :sigma, SIGMA_LEVELS,
        [:model, :upstream_cost, :passability_scenario])
    vs = [v for r in recs for (l, v) in r.signed if l == 0.3 && isfinite(v)]
    if !isempty(vs)
        push!(rows, (factor="thermal_sigma", level="0.3 vs 1.0", val=median(vs) * σ_c / σ_main))
    end

    sort!(rows, by=r -> r.val)
    labels = [r.level for r in rows]
    ys = collect(1:length(rows))
    for (i, r) in enumerate(rows)
        hbar!(ax, i, 0.0, r.val; color=get(factor_color, r.factor, COL_GRAY))
    end
    vlines!(ax, 0.0; color=:black, linewidth=0.8)
    ax.yticks = (ys, labels)
    ax.yticklabelsize = 8
    ax.xlabel = "signed median effect (σ units)"
    ax.title = METRIC_LABELS[metric]
    xmax = maximum(abs.(getfield.(rows, :val))) * 1.25
    xmax = xmax > 0 ? xmax : 1.0
    xlims!(ax, -xmax, xmax)
end

fig = Figure(size=(2200, 2600))
Label(fig[0, :], "Tornado: signed factor-level effects (median across matched blocks)", fontsize=18, font=:bold)
for (i, metric) in enumerate(METRICS)
    r, c = divrem(i - 1, 2)
    ax = Axis(fig[r + 1, c + 1])
    tornado_panels!(ax, metric)
end
legend_factors = ["passability", "upstream_cost", "climate_model", "interaction_matrix",
                  "thermal_sigma"]
Legend(fig[1:4, 3],
    [PolyElement(color=get(factor_color, f, COL_GRAY)) for f in legend_factors],
    replace.(legend_factors, "_" => " ");
    framevisible=false, title="factor", titlefontsize=11)
save_figure(fig, joinpath(REPORT_DIR, "fig_tornado.png"); size=(2200, 2600))

# --- Ranked contribution chart --------------------------------------------
fig2 = Figure(size=(1900, 1000))
Label(fig2[0, :], "Ranked parameter contribution (normalised within-metric effect size)",
    fontsize=18, font=:bold)

ax_rank = Axis(fig2[1, 1];
    xlabel="mean normalised effect (1 = largest within a metric)",
    title="Overall ranking (lower mean rank = more important)")
ys = collect(1:nrow(overall))
for (i, r) in enumerate(eachrow(overall))
    hbar!(ax_rank, i, 0, r.mean_importance;
        color=get(factor_color, r.factor, COL_GRAY))
    text!(ax_rank, r.mean_importance + 0.02, i;
        text=@sprintf("%.2f  (rank %.1f)", r.mean_importance, r.mean_rank),
        align=(:left, :center), fontsize=10)
end
ax_rank.yticks = (ys, replace.(overall.factor, "_" => " "))
ax_rank.yticklabelsize = 11
xlims!(ax_rank, 0, 1.25)

ax_hm = Axis(fig2[1, 2];
    xlabel="factor", ylabel="metric",
    xticks=(1:length(RANK_FACTORS), replace.(RANK_FACTORS, "_" => " ")),
    xticklabelrotation=π / 4,
    yticks=(1:length(METRICS), [METRIC_LABELS[m] for m in METRICS]),
    title="Normalised effect size per metric")
Mmat = fill(NaN, length(METRICS), length(RANK_FACTORS))
for (j, f) in enumerate(RANK_FACTORS), (i, m) in enumerate(METRICS)
    sub = norm[(norm.factor .== f) .& (norm.metric .== m), :]
    isempty(sub) || (Mmat[i, j] = sub.importance[1])
end
hm = heatmap!(ax_hm, 1:length(RANK_FACTORS), 1:length(METRICS), Mmat;
    colormap=:viridis, colorrange=(0, 1))
for i in 1:length(METRICS), j in 1:length(RANK_FACTORS)
    isnan(Mmat[i, j]) && continue
    text!(ax_hm, j, i; text=@sprintf("%.2f", Mmat[i, j]),
        align=(:center, :center), fontsize=8,
        color=Mmat[i, j] > 0.6 ? :white : :black)
end
Colorbar(fig2[1, 3], hm; label="normalised effect")
save_figure(fig2, joinpath(REPORT_DIR, "fig_parameter_importance_ranked.png"); size=(1900, 1000))

# --- Interaction-effect magnitude chart -----------------------------------
inter_rank = rank[rank.factor .∈ Ref(["passability_x_upstream_cost", "passability_x_temperature"]), :]
fig3 = Figure(size=(1600, 800))
Label(fig3[0, :], "Interaction effects on parameter sensitivity", fontsize=18, font=:bold)
ax_i = Axis(fig3[1, 1];
    xlabel="median standardised effect (σ units)",
    yticks=(1:length(METRICS), [METRIC_LABELS[m] for m in METRICS]),
    title="Spread of the passability effect across upstream cost / climate model")
n = length(METRICS)
for (i, m) in enumerate(METRICS)
    sub = inter_rank[inter_rank.metric .== m, :]
    p_uc = sub[sub.factor .== "passability_x_upstream_cost", :]
    p_t = sub[sub.factor .== "passability_x_temperature", :]
    y = i
    isempty(p_uc) || hbar!(ax_i, y + 0.17, 0, p_uc.median_span_std[1];
        color=COL_PXU, height=0.3)
    isempty(p_t) || hbar!(ax_i, y - 0.17, 0, p_t.median_span_std[1];
        color=COL_PXT, height=0.3)
end
Legend(fig3[1, 2], [PolyElement(color=COL_PXU), PolyElement(color=COL_PXT)],
    ["passability × upstream cost", "passability × temperature"], framevisible=false)
save_figure(fig3, joinpath(REPORT_DIR, "fig_interaction_effects.png"); size=(1600, 800))

# --- Dashboard ------------------------------------------------------------
dash_metric = "basin_native_biomass"
fig4 = Figure(size=(2000, 900))
Label(fig4[0, :], "Parameter importance dashboard — GuadeX obstacle sensitivity (2026–2045)",
    fontsize=18, font=:bold)
ax_d = Axis(fig4[1, 1])
tornado_panels!(ax_d, dash_metric)
ax_r = Axis(fig4[1, 2]; xlabel="mean normalised effect",
    title="Overall ranking")
for (i, r) in enumerate(eachrow(overall))
    hbar!(ax_r, i, 0, r.mean_importance; color=get(factor_color, r.factor, COL_GRAY))
end
ax_r.yticks = (collect(1:nrow(overall)), replace.(overall.factor, "_" => " "))
ax_r.yticklabelsize = 10
xlims!(ax_r, 0, 1.15)
save_figure(fig4, joinpath(REPORT_DIR, "fig_parameter_importance_dashboard.png"); size=(2000, 900))

println("Parameter ranking: $(joinpath(REPORT_DIR, "parameter_effect_rank.csv"))")
println("Variance decomposition: $(joinpath(REPORT_DIR, "variance_decomposition.csv"))")
println("Overall ranking (mean_rank, mean_importance):")
for r in eachrow(overall)
    @printf("  %-32s %6.2f  %6.3f\n", r.factor, r.mean_rank, r.mean_importance)
end
