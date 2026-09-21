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
# All figures use the shared report style (large fonts relative to the canvas,
# no descriptive titles, panels lettered a, b, c ...).  The explanation of every
# element lives in the LaTeX caption, not on the canvas.
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

# Display names.  Units are part of the label so no axis is unitless.
const METRIC_LABELS = Dict(
    "basin_total_biomass" => "Total biomass (units)",
    "basin_native_biomass" => "Native biomass (units)",
    "basin_invasive_biomass" => "Invasive biomass (units)",
    "basin_native_richness" => "Native richness (sp.)",
    "basin_native_extinction_risk" => "Richness loss (fraction)",
    "st_final_biomass" => "Brown trout biomass (units)",
    "st_relative_biomass_change" => "Brown trout change (fraction)")

const COST_COLORS = Dict(0.01 => RGBAf(0.27, 0.45, 0.77, 0.95),
    0.5 => RGBAf(0.90, 0.55, 0.10, 0.95))

function save_report(fig, path; size)
    Makie.save(path, fig; size=size, px_per_unit=2)
    trim_figure!(path)
    println("wrote $(basename(path))")
end

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

function horizontal_legend(fig, row, cols, handles, labels; title="")
    Legend(fig[row, cols], handles, labels; orientation=:horizontal,
        framevisible=false, title=title, titlefontsize=REPORT_LEGEND_FONTSIZE)
end

# =============================================================================
# Figure 9 — Brown trout response vs upstream cost (2 x 2, equal panel size)
# =============================================================================
function st_vs_uc_figure!(path; data=main)
    report_theme!()
    uc_levels = sort(unique(Float64.(data.upstream_cost)))
    models = sort(unique(String.(data.model)))
    metrics = ["st_final_biomass", "st_relative_biomass_change"]
    fig = Figure(size=(1500, 1600))
    for (ci, model) in enumerate(models)
        Label(fig[0, ci + 1], model; fontsize=REPORT_AXIS_LABEL_FONTSIZE, font=:bold)
    end
    for (mi, metric) in enumerate(metrics)
        Label(fig[mi, 1], METRIC_LABELS[metric]; rotation=π / 2,
            fontsize=REPORT_AXIS_LABEL_FONTSIZE, font=:bold, halign=:center)
    end
    for (mi, metric) in enumerate(metrics), (ci, model) in enumerate(models)
        ax = Axis(fig[mi, ci + 1]; width=540, height=480,
            xlabel = mi == length(metrics) ? "Upstream movement cost" : "",
            xticks=(1:length(uc_levels), string.(uc_levels)))
        panel_letter!(ax, "(" * string(Char(96 + (mi - 1) * length(models) + ci)) * ")")
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
            scatterlines!(ax, xs, ys; color=PASS_COLORS[p], linewidth=2.8,
                marker=:circle, markersize=13)
        end
    end
    equal_panel_columns!(fig, 0.45, 1, 1)
    horizontal_legend(fig, 3, 1:3,
        [LineElement(color=PASS_COLORS[p], linewidth=3.0) for p in PASS_LEVELS],
        [PASS_LABELS[p] for p in PASS_LEVELS]; title="Passability state")
    save_report(fig, path; size=(1500, 1600))
end

st_vs_uc_figure!(joinpath(OUT, "fig_st_vs_upstream_cost.png"))

# =============================================================================
# Figure 8 — passability x upstream-cost heatmaps (Delta vs baseline)
# =============================================================================
function heatmap_effects!(path; data=main, metrics)
    report_theme!()
    models = sort(unique(String.(data.model)))
    uc_levels = sort(unique(Float64.(data.upstream_cost)))
    nrow_ = length(metrics); ncol = length(models)
    fig = Figure(size=(1700, 2400))
    for (ci, model) in enumerate(models)
        Label(fig[0, ci + 1], model; fontsize=REPORT_AXIS_LABEL_FONTSIZE, font=:bold)
    end
    for (ri, metric) in enumerate(metrics)
        Label(fig[ri, 1], METRIC_LABELS[metric]; rotation=π / 2,
            fontsize=REPORT_AXIS_LABEL_FONTSIZE, font=:bold, halign=:center)
    end
    hm = nothing
    for (ri, metric) in enumerate(metrics), (ci, model) in enumerate(models)
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
        ax = Axis(fig[ri, ci + 1]; width=560, height=460,
            xlabel = ri == nrow_ ? "Upstream movement cost" : "",
            ylabel = ci == 1 ? "Passability state" : "",
            xticks=(1:length(uc_levels), string.(uc_levels)),
            yticks=(1:length(PASS_LEVELS), [PASS_LABELS[p] for p in PASS_LEVELS]))
        panel_letter!(ax, "(" * string(Char(96 + (ri - 1) * ncol + ci)) * ")")
        hm = heatmap!(ax, 1:length(uc_levels), 1:length(PASS_LEVELS), M ./ clim;
            colormap=:balance, colorrange=(-1, 1))
        for pi in 1:length(PASS_LEVELS), ui in 1:length(uc_levels)
            isnan(M[pi, ui]) && continue
            text!(ax, ui, pi; text=@sprintf("%.3g", M[pi, ui]),
                align=(:center, :center), fontsize=REPORT_ANNOTATION_FONTSIZE,
                color=:black, strokecolor=:white, strokewidth=0.8)
        end
    end
    Colorbar(fig[1:nrow_, ncol + 2], hm; label="Change vs baseline (scaled)")
    equal_panel_columns!(fig, 0.45, 1, 1, 0.3)
    save_report(fig, path; size=(1700, 2400))
end

heatmap_effects!(joinpath(OUT, "fig_passability_uc_heatmaps.png");
    metrics=["basin_total_biomass", "basin_native_biomass",
             "basin_native_richness", "st_final_biomass"])

# =============================================================================
# Figure 5 — temperature effect across the 44-run climate ensemble
# =============================================================================
function temperature_figure!(path)
    report_theme!()
    d = clim[clim.scenario .!= "control", :]
    metrics = ["basin_total_biomass", "basin_native_richness",
               "basin_native_extinction_risk"]
    scen_colors = Dict("ssp126" => RGBAf(0.11, 0.36, 0.60, 0.9),
        "ssp245" => RGBAf(0.20, 0.60, 0.35, 0.9),
        "ssp370" => RGBAf(0.90, 0.55, 0.10, 0.9),
        "ssp585" => RGBAf(0.70, 0.15, 0.15, 0.9))
    fig = Figure(size=(1500, 1600))
    for (mi, metric) in enumerate(metrics)
        r, c = divrem(mi - 1, 2)
        ax = Axis(fig[r + 1, c + 1]; width=560, height=480,
            xlabel="Warming by 2045 (°C)",
            ylabel=METRIC_LABELS[metric])
        panel_letter!(ax, "(" * string(Char(96 + mi)) * ")")
        for sc in unique(String.(d.scenario))
            sub = d[d.scenario .== sc, :]
            scatter!(ax, Float64.(sub.warming_end_degc), Float64.(sub[!, Symbol(metric)]);
                color=scen_colors[sc], markersize=13, strokewidth=0)
        end
        x = Float64.(d.warming_end_degc); y = Float64.(d[!, Symbol(metric)])
        X = hcat(ones(length(x)), x)
        β = X \ y
        yhat = X * β
        ss_tot = sum(abs2, y .- mean(y)); ss_res = sum(abs2, y .- yhat)
        r2 = ss_tot > 0 ? 1 - ss_res / ss_tot : 0.0
        xx = range(minimum(x), maximum(x); length=50)
        lines!(ax, xx, β[1] .+ β[2] .* xx; color=:black, linewidth=2.0, linestyle=:dash)
        report_annotation!(ax, @sprintf("OLS slope = %.3g per °C\nR² = %.2f", β[2], r2),
            position=:lt)
        # Mark the two sensitivity-sweep extremes used by the obstacle sweep.
        for (m, lab) in [(COOL_MODEL, "Cool endpoint (SSP1-2.6 / IITM-ESM)"),
                         (WARM_MODEL, "Warm endpoint (SSP2-4.5 / UKESM1)")]
            sub = d[d.model .== m, :]
            isempty(sub) && continue
            scatter!(ax, [Float64(sub.warming_end_degc[1])], [Float64(sub[1, Symbol(metric)])];
                color=:magenta, marker=:star5, markersize=24, strokecolor=:black,
                strokewidth=1.0)
        end
    end
    handles = [MarkerElement(color=scen_colors[s], marker=:circle, markersize=14)
               for s in ["ssp126", "ssp245", "ssp370", "ssp585"]]
    push!(handles, MarkerElement(color=:magenta, marker=:star5, markersize=18,
        strokecolor=:black, strokewidth=1.0))
    labels = ["ssp126", "ssp245", "ssp370", "ssp585",
              "Selected cool/warm\nendpoint (IITM-ESM, UKESM1)"]
    equal_panel_columns!(fig, 1, 1)
    Legend(fig[2, 2], handles, labels; framevisible=false,
        title="Point = one GCM × SSP run", titlefontsize=REPORT_LEGEND_FONTSIZE)
    save_report(fig, path; size=(1500, 1600))
end

temperature_figure!(joinpath(OUT, "fig_temperature_response.png"))

# =============================================================================
# Figure 10 — interaction-matrix robustness of the passability effect
# =============================================================================
function matrix_robustness2!(path)
    report_theme!()
    metrics = ["st_final_biomass", "basin_native_biomass"]
    matrices = ["original", "random", "invasive_favoring"]
    matrix_labels = Dict("original" => "Original", "random" => "Random placebo",
        "invasive_favoring" => "Invasive-favouring")
    models = sort(unique(String.(alt.model)))
    fig = Figure(size=(1600, 1500))
    for (ci, model) in enumerate(models)
        Label(fig[0, ci + 1], model; fontsize=REPORT_AXIS_LABEL_FONTSIZE, font=:bold)
    end
    for (ri, metric) in enumerate(metrics)
        Label(fig[ri, 1], METRIC_LABELS[metric]; rotation=π / 2,
            fontsize=REPORT_AXIS_LABEL_FONTSIZE, font=:bold, halign=:center)
    end
    for (ri, metric) in enumerate(metrics), (ci, model) in enumerate(models)
        ax = Axis(fig[ri, ci + 1]; width=540, height=500,
            xlabel = ri == length(metrics) ? "Span across passability states (model units)" : "",
            ylabel = "Interaction matrix",
            yticks=(1:length(matrices), [matrix_labels[m] for m in matrices]))
        panel_letter!(ax, "(" * string(Char(96 + (ri - 1) * length(models) + ci)) * ")")
        for (mi, mat) in enumerate(matrices)
            for uc in sort(unique(Float64.(alt.upstream_cost)))
                sub = alt[(alt.interaction_matrix .== mat) .& (alt.model .== model) .&
                          (abs.(alt.upstream_cost .- uc) .< 1e-9), :]
                vals = [Float64(sub[sub.passability_scenario .== p, Symbol(metric)][1])
                        for p in PASS_LEVELS if !isempty(sub[sub.passability_scenario .== p, :])]
                length(vals) < 2 && continue
                span = maximum(vals) - minimum(vals)
                off = uc < 0.3 ? -0.18 : 0.18
                hbar!(ax, mi + off, 0, span; color=get(COST_COLORS, uc, RGBAf(0.4, 0.4, 0.4, 0.9)),
                    height=0.30)
            end
        end
    end
    # No legend: the two bars per matrix are described in the caption (blue =
    # low upstream cost 0.01, orange = high upstream cost 0.50).
    equal_panel_columns!(fig, 0.45, 1, 1)
    save_report(fig, path; size=(1600, 1500))
end

matrix_robustness2!(joinpath(OUT, "fig_interaction_matrix_robustness.png"))
matrix_robustness2!(joinpath(ALT_OUT, "fig_interaction_matrix_robustness.png"))

# =============================================================================
# Figure 11 — per-subcatchment / per-water-body blockage effect
# =============================================================================
function read_final(run_dir, level, year)
    df = CSV.read(joinpath(run_dir, "export", "levels", "level_$(level).csv"), DataFrame)
    return df[df.year .== year, :]
end

function subcatchment_waterbody_figure!(path; model=WARM_MODEL, uc=0.5, year=2045)
    report_theme!()
    sub = main[(main.model .== model) .& (abs.(main.upstream_cost .- uc) .< 1e-9), :]
    isempty(sub) && (println("skip subcatchment figure: no run for $model uc=$uc"); return)
    base_dir = sub[sub.passability_scenario .== "baseline", :run_dir][1]
    block_dir = sub[sub.passability_scenario .== "blocked", :run_dir][1]

    sc_b = read_final(base_dir, "subcatchment", year)
    sc_k = read_final(block_dir, "subcatchment", year)
    wb_b = read_final(base_dir, "water_body", year)
    wb_k = read_final(block_dir, "water_body", year)

    sc = innerjoin(select(sc_b, :subcatchment, :mean_native_biomass),
        select(sc_k, :subcatchment, :mean_native_biomass);
        on=:subcatchment, makeunique=true)
    sc.dnative = sc.mean_native_biomass_1 .- sc.mean_native_biomass

    wb = innerjoin(select(wb_b, :water_body, :mean_native_biomass),
        select(wb_k, :water_body, :mean_native_biomass);
        on=:water_body, makeunique=true)
    wb.dnative = wb.mean_native_biomass_1 .- wb.mean_native_biomass

    # Both panels show the 10 most negative and the 10 most positive units (20
    # bars), so the two panels are the same size and the labels stay readable.
    # The full sub-catchment / water-body distributions are summarised in the
    # report text.
    fig = Figure(size=(1400, 1900))
    ax1 = Axis(fig[1, 1]; width=900, height=780, ylabel="Δ native biomass (units)")
    panel_letter!(ax1, "(a)")
    s = vcat(first(sort(sc, :dnative), 10), last(sort(sc, :dnative), 10)) |> unique
    yy = collect(1:nrow(s))
    for (i, r) in enumerate(eachrow(s))
        hbar!(ax1, i, 0, r.dnative; color=r.dnative < 0 ? RGBAf(0.75, 0.15, 0.15, 0.9) :
            RGBAf(0.13, 0.55, 0.30, 0.9), height=0.72)
    end
    vlines!(ax1, 0.0; color=:black, linewidth=1.0)
    ax1.yticks = (yy, string.(s.subcatchment))

    ax2 = Axis(fig[2, 1]; width=900, height=780,
        xlabel="Δ native biomass (units)", ylabel="Water body")
    panel_letter!(ax2, "(b)")
    w = vcat(first(sort(wb, :dnative), 10), last(sort(wb, :dnative), 10)) |> unique
    yy2 = collect(1:nrow(w))
    for (i, r) in enumerate(eachrow(w))
        hbar!(ax2, i, 0, r.dnative; color=r.dnative < 0 ? RGBAf(0.75, 0.15, 0.15, 0.9) :
            RGBAf(0.13, 0.55, 0.30, 0.9), height=0.72)
    end
    vlines!(ax2, 0.0; color=:black, linewidth=1.0)
    ax2.yticks = (yy2, [length(string(x)) > 26 ? string(x)[1:26] * "…" : string(x)
                        for x in w.water_body])
    save_report(fig, path; size=(1400, 1900))
    return (sc=sc, wb=wb)
end

subcatchment_waterbody_figure!(joinpath(OUT, "fig_subcatchment_waterbody_effects.png"))
subcatchment_waterbody_figure!(joinpath(OUT, "fig_subcatchment_waterbody_effects_cool.png");
    model=COOL_MODEL, uc=0.5)

println("Report figures written to:\n  $OUT\n  $ALT_OUT")
for f in sort(readdir(OUT)), g in [OUT]
    isfile(joinpath(g, f)) && endswith(f, ".png") &&
        println(@sprintf("  %-60s %8.1f KB", f, filesize(joinpath(g, f)) / 1024))
end
