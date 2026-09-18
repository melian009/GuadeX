using Pkg; Pkg.activate(".")
using DataFrames
using CSV
using CairoMakie
using Statistics
using Printf
using Guadex

# =============================================================================
# Report-ready figures for the GuadeX obstacle / interaction sensitivity study.
#
#   julia --project=. scripts/plot_report_figures.jl
#
# Reads (read-only):
#   results/sensitivity_obstacles/runs_index.csv                (32 runs)
#   results/sensitivity_obstacles/alt_interactions/runs_index.csv (48 runs)
#   results/climate_scenarios_k1x_burnin/runs_index.csv         (44 scenario runs)
#   <run_dir>/export/levels/level_subcatchment.csv
#   <run_dir>/export/levels/level_water_body.csv
#
# Writes PNGs to:
#   results/sensitivity_obstacles/report_plots/
#   results/sensitivity_obstacles/alt_interactions/report_plots/
#
# Runs 2026-2045, per-site daily water temperature, spin-up 373 yr, K=1x,
# heat-stress k=3.7817692640400324e-5 (verified from run_metadata.json).
# =============================================================================

const ROOT = joinpath("results", "sensitivity_obstacles")
const ALT_ROOT = joinpath(ROOT, "alt_interactions")
const CLIMATE_INDEX = joinpath("results", "climate_scenarios_k1x_burnin", "runs_index.csv")
const OUT = joinpath(ROOT, "report_plots")
const ALT_OUT = joinpath(ALT_ROOT, "report_plots")
mkpath(OUT); mkpath(ALT_OUT)

const PASS_LEVELS = ["baseline", "improved_passability", "reduced_passability", "blocked"]
const PASS_COLORS = Dict(
    "baseline" => RGBAf(0.20, 0.20, 0.20, 0.9),
    "improved_passability" => RGBAf(0.13, 0.55, 0.75, 0.95),
    "reduced_passability" => RGBAf(0.95, 0.60, 0.10, 0.95),
    "blocked" => RGBAf(0.75, 0.15, 0.15, 0.95))
const PASS_LABELS = Dict(
    "baseline" => "Baseline", "improved_passability" => "Improved",
    "reduced_passability" => "Reduced", "blocked" => "Blocked")

const COOL_MODEL = "ssp126/IITM-ESM"
const WARM_MODEL = "ssp245/UKESM1-0-LL"
const METRIC_LABELS = Dict(
    "basin_total_biomass" => "Basin total biomass",
    "basin_native_biomass" => "Basin native biomass",
    "basin_invasive_biomass" => "Basin invasive biomass",
    "basin_native_richness" => "Basin native richness",
    "basin_native_extinction_risk" => "Native extinction risk",
    "st_final_biomass" => "ST final biomass",
    "st_relative_biomass_change" => "ST relative biomass change")

main = CSV.read(joinpath(ROOT, "runs_index.csv"), DataFrame)
alt = CSV.read(joinpath(ALT_ROOT, "runs_index.csv"), DataFrame)
clim = CSV.read(CLIMATE_INDEX, DataFrame)

main.model = string.(main.climate_scenario, "/", main.gcm)
alt.model = string.(alt.climate_scenario, "/", alt.gcm)
clim.model = string.(clim.scenario, "/", clim.gcm)

function hbar!(ax, y, x0, x1; color=RGBAf(0.27, 0.45, 0.77, 0.9), height=0.62)
    xlo, xhi = min(x0, x1), max(x0, x1)
    poly!(ax, [Point2f(xlo, y - height / 2), Point2f(xhi, y - height / 2),
               Point2f(xhi, y + height / 2), Point2f(xlo, y + height / 2)], color=color)
end

# =============================================================================
# Figure 1 — ST (cold-water keystone) response vs upstream cost
# =============================================================================
function st_vs_uc_figure!(path; data=main, title="ST response to upstream cost and passability")
    uc_levels = sort(unique(Float64.(data.upstream_cost)))
    models = sort(unique(String.(data.model)))
    metrics = ["st_final_biomass", "st_relative_biomass_change"]
    fig = Figure(size=(1500, 1150))
    Label(fig[0, :], title, fontsize=17, font=:bold)
    for (mi, metric) in enumerate(metrics), (ci, model) in enumerate(models)
        ax = Axis(fig[mi, ci];
            xlabel="Upstream cost", ylabel=mi == 1 ? METRIC_LABELS[metric] : "",
            title=ci == 1 ? "$model" : "$model",
            xticks=(1:length(uc_levels), string.(uc_levels)))
        for p in PASS_LEVELS
            xs = Float64[]; ys = Float64[]
            for (i, uc) in enumerate(uc_levels)
                sel = data[(data.model .== model) .&
                           (abs.(data.upstream_cost .- uc) .< 1e-9) .&
                           (data.passability_scenario .== p), :]
                isempty(sel) && continue
                push!(xs, i); push!(ys, Float64(sel[1, Symbol(metric)]))
            end
            isempty(xs) && continue
            scatterlines!(ax, xs, ys; color=PASS_COLORS[p], linewidth=2.2,
                marker=:circle, markersize=10, label=PASS_LABELS[p])
        end
        mi == 1 && ci == 1 && axislegend(ax; position=:rt, framevisible=false, labelsize=10)
    end
    save_figure(fig, path; size=(1500, 1150))
end

st_vs_uc_figure!(joinpath(OUT, "fig_st_vs_upstream_cost.png"))

# =============================================================================
# Figure 2 — passability x upstream-cost heatmaps (Delta vs baseline), per model
# =============================================================================
function heatmap_effects!(path; data=main, metrics, title)
    models = sort(unique(String.(data.model)))
    uc_levels = sort(unique(Float64.(data.upstream_cost)))
    nrow_ = length(models); ncol = length(metrics)
    fig = Figure(size=(430 * ncol + 260, 330 * nrow_ + 90))
    Label(fig[0, 1:ncol], title, fontsize=17, font=:bold)
    hm = nothing
    for (ri, model) in enumerate(models), (ci, metric) in enumerate(metrics)
        M = fill(NaN, length(PASS_LEVELS), length(uc_levels))
        for (pi, p) in enumerate(PASS_LEVELS), (ui, uc) in enumerate(uc_levels)
            sel = data[(data.model .== model) .& (data.passability_scenario .== p) .&
                       (abs.(data.upstream_cost .- uc) .< 1e-9), :]
            bas = data[(data.model .== model) .& (data.passability_scenario .== "baseline") .&
                       (abs.(data.upstream_cost .- uc) .< 1e-9), :]
            (isempty(sel) || isempty(bas)) && continue
            M[pi, ui] = Float64(sel[1, Symbol(metric)]) - Float64(bas[1, Symbol(metric)])
        end
        finite = filter(isfinite, vec(M))
        clim = isempty(finite) ? 1.0 : max(quantile(abs.(finite), 0.95), 1e-12)
        ax = Axis(fig[ri, ci];
            title=ci == 1 ? "$model" : "",
            xlabel=ri == nrow_ ? "Upstream cost" : "",
            ylabel=ci == 1 ? "Passability" : "",
            xticks=(1:length(uc_levels), string.(uc_levels)),
            yticks=(1:length(PASS_LEVELS), [PASS_LABELS[p] for p in PASS_LEVELS]))
        hm = heatmap!(ax, 1:length(uc_levels), 1:length(PASS_LEVELS), M ./ clim;
            colormap=:balance, colorrange=(-1, 1))
        for pi in 1:length(PASS_LEVELS), ui in 1:length(uc_levels)
            isnan(M[pi, ui]) && continue
            text!(ax, ui, pi; text=@sprintf("%.3g", M[pi, ui]),
                align=(:center, :center), fontsize=8, color=:black,
                strokecolor=:white, strokewidth=0.7)
        end
    end
    Colorbar(fig[nrow_, ncol + 1], hm;
        label="Δ vs baseline (scaled by panel 95th pct)", width=16)
    save_figure(fig, path; size=(430 * ncol + 260, 330 * nrow_ + 90))
end

heatmap_effects!(joinpath(OUT, "fig_passability_uc_heatmaps.png");
    metrics=["basin_total_biomass", "basin_native_biomass", "basin_native_richness", "st_final_biomass"],
    title="Passability × upstream cost — change vs baseline passability (obstacle sweep)")

"""
Facet passability x upstream-cost change grids by interaction matrix x model x
metric.  The `random` matrix is kept for completeness but flagged: its
near-zero burn-in baselines make the *absolute* changes large and the colour
scaling per panel is robust (95th percentile).
"""
function matrix_heatmaps!(path)
    matrices = ["original", "random", "invasive_favoring"]
    models = sort(unique(String.(alt.model)))
    metrics = ["st_final_biomass", "basin_native_biomass"]
    uc_levels = sort(unique(Float64.(alt.upstream_cost)))
    cols = [(m, met) for m in models for met in metrics]
    fig = Figure(size=(380 * length(cols) + 200, 300 * length(matrices) + 220))
    Label(fig[0, 1:length(cols)],
        "Passability × upstream cost change by interaction matrix (σ=0.3 sweep) — " *
        "colours scaled per panel (95th pct); numbers are raw Δ",
        fontsize=14, font=:bold)
    for (ci, (model, metric)) in enumerate(cols)
        Label(fig[1, ci], "$model\n$(METRIC_LABELS[metric])"; fontsize=11, font=:bold)
    end
    hm = nothing
    for (ri, mat) in enumerate(matrices), (ci, (model, metric)) in enumerate(cols)
        M = fill(NaN, length(PASS_LEVELS), length(uc_levels))
        for (pi, p) in enumerate(PASS_LEVELS), (ui, uc) in enumerate(uc_levels)
            sel = alt[(alt.interaction_matrix .== mat) .& (alt.model .== model) .&
                      (alt.passability_scenario .== p) .&
                      (abs.(alt.upstream_cost .- uc) .< 1e-9), :]
            bas = alt[(alt.interaction_matrix .== mat) .& (alt.model .== model) .&
                      (alt.passability_scenario .== "baseline") .&
                      (abs.(alt.upstream_cost .- uc) .< 1e-9), :]
            (isempty(sel) || isempty(bas)) && continue
            M[pi, ui] = Float64(sel[1, Symbol(metric)]) - Float64(bas[1, Symbol(metric)])
        end
        finite = filter(isfinite, vec(M))
        clim = isempty(finite) ? 1.0 : max(quantile(abs.(finite), 0.95), 1e-12)
        ax = Axis(fig[ri + 1, ci];
            xlabel=ri == length(matrices) ? "Upstream cost" : "",
            ylabel=ci == 1 ? replace(mat, "_" => " ") : "",
            xticks=(1:length(uc_levels), string.(uc_levels)),
            yticks=(1:length(PASS_LEVELS), [PASS_LABELS[p] for p in PASS_LEVELS]))
        hm = heatmap!(ax, 1:length(uc_levels), 1:length(PASS_LEVELS), M ./ clim;
            colormap=:balance, colorrange=(-1, 1))
        for pi in 1:length(PASS_LEVELS), ui in 1:length(uc_levels)
            isnan(M[pi, ui]) && continue
            text!(ax, ui, pi; text=@sprintf("%.3g", M[pi, ui]),
                align=(:center, :center), fontsize=7, color=:black,
                strokecolor=:white, strokewidth=0.7)
        end
    end
    Colorbar(fig[length(matrices) + 1, length(cols) + 1], hm;
        label="Δ vs baseline (scaled by panel 95th pct)", width=16)
    save_figure(fig, path; size=(380 * length(cols) + 200, 300 * length(matrices) + 220))
end

matrix_heatmaps!(joinpath(ALT_OUT, "fig_matrix_uc_heatmaps.png"))

# =============================================================================
# Figure 3 — temperature effect across the 44-run climate ensemble
# =============================================================================
function temperature_figure!(path)
    d = clim[clim.scenario .!= "control", :]
    metrics = ["basin_total_biomass", "basin_native_richness", "basin_native_extinction_risk"]
    scen_colors = Dict("ssp126" => RGBAf(0.11, 0.36, 0.60, 0.9),
        "ssp245" => RGBAf(0.20, 0.60, 0.35, 0.9),
        "ssp370" => RGBAf(0.90, 0.55, 0.10, 0.9),
        "ssp585" => RGBAf(0.70, 0.15, 0.15, 0.9))
    fig = Figure(size=(2000, 560))
    Label(fig[0, :], "Basin response vs end-of-century warming across 44 climate-model runs " *
        "(k=1×, spin-up 373 yr, 2026–2045; sensitivity-sweep extremes marked)",
        fontsize=16, font=:bold)
    for (mi, metric) in enumerate(metrics)
        ax = Axis(fig[1, mi];
            xlabel="End-of-century warming (°C)", ylabel=METRIC_LABELS[metric],
            title=METRIC_LABELS[metric])
        for sc in unique(String.(d.scenario))
            sub = d[d.scenario .== sc, :]
            scatter!(ax, Float64.(sub.warming_end_degc), Float64.(sub[!, Symbol(metric)]);
                color=scen_colors[sc], markersize=11, label=sc, strokewidth=0)
        end
        x = Float64.(d.warming_end_degc); y = Float64.(d[!, Symbol(metric)])
        X = hcat(ones(length(x)), x)
        β = X \ y
        yhat = X * β
        ss_tot = sum(abs2, y .- mean(y)); ss_res = sum(abs2, y .- yhat)
        r2 = ss_tot > 0 ? 1 - ss_res / ss_tot : 0.0
        xx = range(minimum(x), maximum(x); length=50)
        lines!(ax, xx, β[1] .+ β[2] .* xx; color=:black, linewidth=1.6, linestyle=:dash)
        text!(ax, 0.03, 0.96; space=:relative, align=(:left, :top), fontsize=11,
            text=@sprintf("slope = %.3g /°C\nR² = %.2f", β[2], r2))
        # Mark the two sensitivity-sweep extremes.
        for (m, lab) in [(COOL_MODEL, "SSP126 / IITM-ESM (coolest)"),
                         (WARM_MODEL, "SSP245 / UKESM1 (warmest)")]
            sub = d[d.model .== m, :]
            isempty(sub) && continue
            scatter!(ax, [Float64(sub.warming_end_degc[1])], [Float64(sub[1, Symbol(metric)])];
                color=:magenta, marker=:star5, markersize=20, strokecolor=:black,
                strokewidth=1.0, label=lab)
        end
        mi == 1 && axislegend(ax; position=:rb, framevisible=false, labelsize=9)
    end
    save_figure(fig, path; size=(2000, 560))
end

temperature_figure!(joinpath(OUT, "fig_temperature_response.png"))

# =============================================================================
# Figure 4 — interaction-matrix robustness of the passability effect
# =============================================================================
function matrix_robustness2!(path)
    metrics = ["st_final_biomass", "basin_native_biomass"]
    matrices = ["original", "random", "invasive_favoring"]
    mcolors = Dict("original" => RGBAf(0.20, 0.45, 0.75, 0.9),
                   "random" => RGBAf(0.55, 0.55, 0.55, 0.9),
                   "invasive_favoring" => RGBAf(0.75, 0.30, 0.10, 0.9))
    models = sort(unique(String.(alt.model)))
    fig = Figure(size=(2000, 700))
    Label(fig[0, 1:2], "Robustness of the passability effect across interaction matrices " *
        "(span across passability levels; σ=0.3 sweep, two upstream costs per matrix)",
        fontsize=15, font=:bold)
    for (ri, metric) in enumerate(metrics), (ci, model) in enumerate(models)
        ax = Axis(fig[ri, ci];
            title="$(METRIC_LABELS[metric]) — $model",
            xlabel=ri == length(metrics) ? "span across passability levels" : "",
            ylabel=ci == 1 ? "interaction matrix" : "",
            yticks=(1:length(matrices), replace.(matrices, "_" => " ")))
        for (mi, mat) in enumerate(matrices)
            for uc in sort(unique(Float64.(alt.upstream_cost)))
                sub = alt[(alt.interaction_matrix .== mat) .& (alt.model .== model) .&
                          (abs.(alt.upstream_cost .- uc) .< 1e-9), :]
                vals = [Float64(sub[sub.passability_scenario .== p, Symbol(metric)][1])
                        for p in PASS_LEVELS if !isempty(sub[sub.passability_scenario .== p, :])]
                length(vals) < 2 && continue
                span = maximum(vals) - minimum(vals)
                off = uc < 0.3 ? -0.18 : 0.18
                hbar!(ax, mi + off, 0, span; color=mcolors[mat], height=0.28)
            end
        end
    end
    Legend(fig[1:2, 3], [PolyElement(color=mcolors[m]) for m in matrices],
        [replace(m, "_" => " ") for m in matrices]; framevisible=false, title="matrix")
    save_figure(fig, path; size=(2000, 700))
end

matrix_robustness2!(joinpath(OUT, "fig_interaction_matrix_robustness.png"))
matrix_robustness2!(joinpath(ALT_OUT, "fig_interaction_matrix_robustness.png"))

# =============================================================================
# Figure 5 — per-subcatchment / per-water-body effect for the strongest signal
# =============================================================================
function read_final(run_dir, level, year)
    df = CSV.read(joinpath(run_dir, "export", "levels", "level_$(level).csv"), DataFrame)
    return df[df.year .== year, :]
end

function subcatchment_waterbody_figure!(path; model=WARM_MODEL, uc=0.5, year=2045)
    sub = main[(main.model .== model) .& (abs.(main.upstream_cost .- uc) .< 1e-9), :]
    isempty(sub) && (println("skip subcatchment figure: no run for $model uc=$uc"); return)
    base_dir = sub[sub.passability_scenario .== "baseline", :run_dir][1]
    block_dir = sub[sub.passability_scenario .== "blocked", :run_dir][1]

    sc_b = read_final(base_dir, "subcatchment", year)
    sc_k = read_final(block_dir, "subcatchment", year)
    wb_b = read_final(base_dir, "water_body", year)
    wb_k = read_final(block_dir, "water_body", year)

    sc = innerjoin(select(sc_b, :subcatchment, :mean_native_biomass, :mean_total_biomass),
        select(sc_k, :subcatchment, :mean_native_biomass, :mean_total_biomass);
        on=:subcatchment, makeunique=true)
    sc.dnative = sc.mean_native_biomass_1 .- sc.mean_native_biomass
    sc.dtotal = sc.mean_total_biomass_1 .- sc.mean_total_biomass

    wb = innerjoin(select(wb_b, :water_body, :mean_native_biomass),
        select(wb_k, :water_body, :mean_native_biomass);
        on=:water_body, makeunique=true)
    wb.dnative = wb.mean_native_biomass_1 .- wb.mean_native_biomass

    fig = Figure(size=(2000, 700))
    Label(fig[0, :], "Obstacle-overlay effect: blocked minus baseline passability " *
        "($model, upstream cost=$uc, year $year)", fontsize=16, font=:bold)

    ax1 = Axis(fig[1, 1]; title="Per subcatchment — Δ native biomass",
        ylabel="Δ native biomass")
    s = sort(sc, :dnative)
    yy = collect(1:nrow(s))
    for (i, r) in enumerate(eachrow(s))
        hbar!(ax1, i, 0, r.dnative; color=r.dnative < 0 ? RGBAf(0.75, 0.15, 0.15, 0.9) :
            RGBAf(0.13, 0.55, 0.30, 0.9), height=0.72)
    end
    vlines!(ax1, 0.0; color=:black, linewidth=0.8)
    ax1.yticks = (yy, string.(s.subcatchment)); ax1.yticklabelsize = 7

    ax2 = Axis(fig[1, 2]; title="Per water body — Δ native biomass",
        ylabel="Δ native biomass")
    w = sort(wb, :dnative)
    # show the 25 most negative and 25 most positive
    w = vcat(first(w, 25), last(w, 25)) |> unique
    yy2 = collect(1:nrow(w))
    for (i, r) in enumerate(eachrow(w))
        hbar!(ax2, i, 0, r.dnative; color=r.dnative < 0 ? RGBAf(0.75, 0.15, 0.15, 0.9) :
            RGBAf(0.13, 0.55, 0.30, 0.9), height=0.72)
    end
    vlines!(ax2, 0.0; color=:black, linewidth=0.8)
    ax2.yticks = (yy2, [length(string(x)) > 22 ? string(x)[1:22] * "…" : string(x)
                        for x in w.water_body]); ax2.yticklabelsize = 6

    ax3 = Axis(fig[1, 3]; title="Distribution of Δ native biomass",
        ylabel="Δ native biomass", xticks=(1:2, ["Subcatchments", "Water bodies"]))
    for (i, v) in enumerate([sc.dnative, wb.dnative])
        vals = filter(isfinite, Float64.(v))
        isempty(vals) && continue
        scatter!(ax3, fill(i, length(vals)), vals; color=RGBAf(0.3, 0.3, 0.3, 0.35),
            markersize=7, strokewidth=0)
        lines!(ax3, [i - 0.22, i + 0.22], fill(median(vals), 2);
            color=:black, linewidth=2.5)
    end
    hlines!(ax3, 0.0; color=:red, linestyle=:dash, linewidth=1.0)
    save_figure(fig, path; size=(2000, 700))
end

subcatchment_waterbody_figure!(joinpath(OUT, "fig_subcatchment_waterbody_effects.png"))
subcatchment_waterbody_figure!(joinpath(OUT, "fig_subcatchment_waterbody_effects_cool.png");
    model=COOL_MODEL, uc=0.5)

println("Report figures written to:\n  $OUT\n  $ALT_OUT")
for f in sort(readdir(OUT)), g in [OUT]
    isfile(joinpath(g, f)) && endswith(f, ".png") &&
        println(@sprintf("  %-60s %8.1f KB", f, filesize(joinpath(g, f)) / 1024))
end
