using Pkg; Pkg.activate(".")
using DataFrames
using CSV
using CairoMakie
using Statistics
using Guadex

# =============================================================================
# Summary of the obstacle / upstream-cost sensitivity sweeps.
#
# Reads the `runs_index.csv` written by run_sensitivity_report.jl and
# run_alt_interactions.jl and plots, for every interaction matrix x climate
# model, the effect of each passability scenario relative to the baseline
# passability run at the same upstream cost.
#
#   julia --project=. scripts/plot_obstacle_sensitivity.jl [results_root] [figures_dir]
#
# Defaults: results_root = results/sensitivity_obstacles,
#           figures_dir = <results_root>/summary_plots
# =============================================================================

const RESULTS_ROOT = length(ARGS) >= 1 ? ARGS[1] : joinpath("results", "sensitivity_obstacles")
const FIGURES_DIR = length(ARGS) >= 2 ? ARGS[2] : joinpath(RESULTS_ROOT, "summary_plots")

const INDEX_PATH = joinpath(RESULTS_ROOT, "runs_index.csv")
const PASSABILITY_ORDER = ["baseline", "improved_passability", "reduced_passability", "blocked"]
const PASSABILITY_LABELS = Dict(
    "baseline" => "Baseline", "improved_passability" => "Improved",
    "reduced_passability" => "Reduced", "blocked" => "Blocked")

# (column, display label, palette)
const EFFECT_METRICS = [
    ("basin_total_biomass", "Δ total biomass", :balance),
    ("basin_native_biomass", "Δ native biomass", :balance),
    ("basin_native_richness", "Δ native richness", :balance),
    ("st_relative_biomass_change", "ST relative biomass change", :balance),
]

isfile(INDEX_PATH) || error("runs index not found: $INDEX_PATH (run the sweep first)")
raw = CSV.read(INDEX_PATH, DataFrame)
isempty(raw) && error("runs index is empty: $INDEX_PATH")

has_matrix = hasproperty(raw, :interaction_matrix)

function passability_rank(name)
    idx = findfirst(==(String(name)), PASSABILITY_ORDER)
    return idx === nothing ? length(PASSABILITY_ORDER) + 1 : idx
end

# Group keys: interaction matrix (if present) x climate model.
groups = Vector{NamedTuple}()
if has_matrix
    for m in unique(String.(raw.interaction_matrix))
        sub = raw[raw.interaction_matrix .== m, :]
        for key in unique([(String(r.climate_scenario), String(r.gcm)) for r in eachrow(sub)])
            push!(groups, (matrix=m, scenario=key[1], gcm=key[2]))
        end
    end
else
    for key in unique([(String(r.climate_scenario), String(r.gcm)) for r in eachrow(raw)])
        push!(groups, (matrix="", scenario=key[1], gcm=key[2]))
    end
end

function group_rows(group)
    sub = raw
    has_matrix && (sub = sub[sub.interaction_matrix .== group.matrix, :])
    return sub[(sub.climate_scenario .== group.scenario) .& (sub.gcm .== group.gcm), :]
end

"""
    effect_matrix(rows, metric, upstream_costs, passabilities)

`n_pass × n_uc` matrix of the metric *difference* from the baseline-passability
run at the same upstream cost (NaN where a run is missing).
"""
function effect_matrix(rows, metric, upstream_costs, passabilities)
    M = fill(NaN, length(passabilities), length(upstream_costs))
    for (pi, pass) in enumerate(passabilities)
        for (ui, uc) in enumerate(upstream_costs)
            sel = rows[(rows.passability_scenario .== pass) .&
                       (abs.(rows.upstream_cost .- uc) .< 1e-9), :]
            isempty(sel) && continue
            value = Float64(sel[1, Symbol(metric)])
            base = rows[(rows.passability_scenario .== "baseline") .&
                        (abs.(rows.upstream_cost .- uc) .< 1e-9), :]
            reference = isempty(base) ? value : Float64(base[1, Symbol(metric)])
            M[pi, ui] = value - reference
        end
    end
    return M
end

mkpath(FIGURES_DIR)
effect_rows = NamedTuple[]

for group in groups
    rows = group_rows(group)
    nrow(rows) == 0 && continue
    upstream_costs = sort(unique(Float64.(rows.upstream_cost)))
    passabilities = sort(unique(String.(rows.passability_scenario)); by=passability_rank)

    group_label = "obstacle sensitivity"
    has_matrix && (group_label *= " — $(group.matrix)")
    group_label *= " — $(group.scenario) / $(group.gcm)"

    n_metrics = length(EFFECT_METRICS)
    fig = Figure(size=(430 * n_metrics + 60, 440))

    for (mi, (metric, label, _)) in enumerate(EFFECT_METRICS)
        M = effect_matrix(rows, metric, upstream_costs, passabilities)
        finite = filter(isfinite, vec(M))
        # A near-zero burn-in baseline (possible under the null interaction
        # models) can make a ratio explode; use a robust colour limit so one
        # outlier cannot wash out the rest of the panel.
        if isempty(finite)
            clim = 1.0
        else
            q = quantile(abs.(finite), 0.95)
            clim = q > 0 ? q : max(maximum(abs.(finite)), 1e-12)
        end

        ax = Axis(fig[1, mi];
            title=label,
            xlabel="Upstream cost",
            ylabel=mi == 1 ? "Passability scenario" : "",
            xticks=(1:length(upstream_costs), string.(upstream_costs)),
            yticks=(1:length(passabilities),
                [get(PASSABILITY_LABELS, p, p) for p in passabilities]))
        hm = heatmap!(ax, 1:length(upstream_costs), 1:length(passabilities), M;
            colormap=:balance, colorrange=(-clim, clim))
        for pi in 1:length(passabilities), ui in 1:length(upstream_costs)
            isnan(M[pi, ui]) && continue
            text!(ax, ui, pi; text=string(round(M[pi, ui], digits=3)),
                align=(:center, :center), fontsize=9,
                color=abs(M[pi, ui]) > 0.55 * clim ? :white : :black)
        end
        Colorbar(fig[1, mi + n_metrics], hm; label="Δ vs baseline passability", width=18)
    end

    Label(fig[0, 1:n_metrics], group_label, fontsize=15, font=:bold)

    tag = has_matrix ? "$(group.matrix)__$(group.scenario)__$(group.gcm)" :
        "$(group.scenario)__$(group.gcm)"
    tag = replace(tag, r"[^A-Za-z0-9._-]" => "_")
    save_figure(fig, joinpath(FIGURES_DIR, "obstacle_effects_$tag.png"); size=(430 * n_metrics, 420))

    # Long-form effect table for downstream reporting.
    for (metric, label, _) in EFFECT_METRICS, (pi, pass) in enumerate(passabilities)
        for (ui, uc) in enumerate(upstream_costs)
            sel = rows[(rows.passability_scenario .== pass) .&
                       (abs.(rows.upstream_cost .- uc) .< 1e-9), :]
            isempty(sel) && continue
            value = Float64(sel[1, Symbol(metric)])
            base = rows[(rows.passability_scenario .== "baseline") .&
                        (abs.(rows.upstream_cost .- uc) .< 1e-9), :]
            reference = isempty(base) ? value : Float64(base[1, Symbol(metric)])
            push!(effect_rows, (
                interaction_matrix=group.matrix, climate_scenario=group.scenario,
                gcm=group.gcm, upstream_cost=uc, passability_scenario=pass,
                metric=metric, value=value, baseline_value=reference,
                effect=value - reference))
        end
    end
end

if !isempty(effect_rows)
    CSV.write(joinpath(RESULTS_ROOT, "sensitivity_effects.csv"), DataFrame(effect_rows))
    println("Effects table: $(joinpath(RESULTS_ROOT, "sensitivity_effects.csv"))")
end
println("Summary figures: $FIGURES_DIR")
