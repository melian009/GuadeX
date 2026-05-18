using Pkg; Pkg.activate(".")
using JLD2
using CairoMakie
using Statistics
using Guadex

"""
# USAGE

## Original behavior (flat structure)
julia scripts/plot_sensitivity_summary.jl results/sensitivity_temp_passability

## Alt interactions (matrix-nested structure)
julia scripts/plot_sensitivity_summary.jl results/sensitivity_alt_interactions

## Default (no arg)
julia scripts/plot_sensitivity_summary.jl
"""

const DAYS_PER_YEAR = 365
const SIMULATION_YEARS = 3
const NATIVE_SPECIES = ["AB", "AH", "SP", "PW", "LS", "SA", "IL", "CP", "IO"]
const INVASIVE_SPECIES = ["GH", "MS", "LG", "CC", "CG", "AM", "OM", "EL", "GL", "TT"]

# --- Configuration ---
if length(ARGS) >= 1
    BASE_DIR = ARGS[1]
else
    BASE_DIR = "results/sensitivity_temp_passability"
end

# Auto-detect matrix-level nesting
base_entries = readdir(BASE_DIR; join=false, sort=false)
matrix_names = String[]
has_matrices = false
for entry in base_entries
    subdir = joinpath(BASE_DIR, entry)
    if isdir(subdir) && !(entry in ["summary_plots", "sensitivity_summary_plots"])
        nested = readdir(subdir; join=false, sort=false)
        if any(startswith("dT_"), nested)
            push!(matrix_names, entry)
        end
    end
end
sort!(matrix_names)
has_matrices = !isempty(matrix_names)

if has_matrices
    println("Detected $(length(matrix_names)) interaction matrix variants: $(join(matrix_names, ", "))")
else
    println("No matrix-level nesting detected – working with flat structure.")
end

temperature_increases = [0.0, 1.0, 2.0, 3.0]
upstream_costs = [0.01, 0.05, 0.1, 0.5]
passability_scenarios = ["baseline", "improved_passability", "reduced_passability", "blocked"]
pass_labels = Dict(
    "baseline" => "Baseline",
    "improved_passability" => "Improved Passability",
    "reduced_passability" => "Reduced Passability",
    "blocked" => "Blocked",
)

function classify_species_indices(all_species, target_group)
    return [findfirst(==(sp), all_species) for sp in target_group if sp in all_species]
end

function jld2_path(base_dir, prefix, dt, uc, pass_name)
    if isempty(prefix)
        return joinpath(base_dir, "dT_$(dt)C", "uc_$(uc)", pass_name, "simulation_output.jld2")
    else
        return joinpath(base_dir, prefix, "dT_$(dt)C", "uc_$(uc)", pass_name, "simulation_output.jld2")
    end
end

# =============================================================================
# --- Load Scalar Metrics From All Runs ---
# =============================================================================

function load_all_metrics(base_dir; path_prefix="")
    target = isempty(path_prefix) ? base_dir : joinpath(base_dir, path_prefix)
    println("Loading metrics from $target ...")
    metrics = []

    for dt in temperature_increases
        for uc in upstream_costs
            for pass_name in passability_scenarios
                jld2_file = jld2_path(base_dir, path_prefix, dt, uc, pass_name)
                if !isfile(jld2_file)
                    continue
                end

                try
                    jldopen(jld2_file, "r") do f
                        sol_u = f["sol_u"]
                        species = f["species"]
                        sites = f["sites"]
                        n_sites = length(sites)
                        n_species = length(species)

                        native_idx = classify_species_indices(species, NATIVE_SPECIES)
                        invasive_idx = classify_species_indices(species, INVASIVE_SPECIES)

                        u0 = reshape(sol_u[1], n_sites, n_species)
                        uf = reshape(sol_u[end], n_sites, n_species)

                        total_biomass_0 = sum(u0)
                        total_biomass_T = sum(uf)

                        native_richness_0 = mean(Float64[sum(u0[i, native_idx] .> 0.1) for i in 1:n_sites])
                        native_richness_T = mean(Float64[sum(uf[i, native_idx] .> 0.1) for i in 1:n_sites])

                        invasive_richness_0 = mean(Float64[sum(u0[i, invasive_idx] .> 0.1) for i in 1:n_sites])
                        invasive_richness_T = mean(Float64[sum(uf[i, invasive_idx] .> 0.1) for i in 1:n_sites])

                        push!(metrics, (
                            dt = dt,
                            uc = uc,
                            pass = pass_name,
                            total_biomass_0 = total_biomass_0,
                            total_biomass_T = total_biomass_T,
                            biomass_change = total_biomass_T - total_biomass_0,
                            native_richness_0 = native_richness_0,
                            native_richness_T = native_richness_T,
                            native_richness_change = native_richness_T - native_richness_0,
                            invasive_richness_0 = invasive_richness_0,
                            invasive_richness_T = invasive_richness_T,
                            invasive_richness_change = invasive_richness_T - invasive_richness_0,
                        ))
                    end
                catch e
                    println("  ERROR loading $jld2_file: $e")
                    continue
                end
            end
        end
    end

    println("  Loaded $(length(metrics)) simulation results.")
    return metrics
end

function get_metric_matrix(metrics, pass_name; field=:native_richness_change)
    nT = length(temperature_increases)
    nU = length(upstream_costs)
    M = fill(NaN, nU, nT)
    for m in metrics
        if m.pass == pass_name
            ti = findfirst(==(m.dt), temperature_increases)
            ui = findfirst(==(m.uc), upstream_costs)
            if ti !== nothing && ui !== nothing
                M[ui, ti] = getfield(m, field)
            end
        end
    end
    return M, temperature_increases, upstream_costs
end

# Detect which passability scenarios actually have data
function active_passability_scenarios(metrics)
    pass_with_data = unique([m.pass for m in metrics])
    return [p for p in passability_scenarios if p in pass_with_data]
end

# =============================================================================
# --- Figure 1: Heatmap Grid — ΔRichness: Native & Invasive ---
# =============================================================================

function plot_richness_change_heatmaps(metrics, output_dir; title_suffix="")
    println("Generating richness change heatmaps...")

    diverging_cmap = cgrad([RGBf(0.129, 0.4, 0.675), RGBf(0.97, 0.97, 0.97), RGBf(0.698, 0.094, 0.169)])

    active_pass = active_passability_scenarios(metrics)
    n_pass = length(active_pass)
    if n_pass == 0
        println("  No data — skipping heatmaps.")
        return
    end

    n_rows = 2
    n_cols = n_pass

    all_native = Float64[getfield(m, :native_richness_change) for m in metrics]
    all_invasive = Float64[getfield(m, :invasive_richness_change) for m in metrics]
    all_vals = vcat(all_native, all_invasive)
    valid_vals = filter(!isnan, all_vals)
    if isempty(valid_vals)
        println("  No valid data — skipping heatmaps.")
        return
    end
    clim_max = maximum(abs.(valid_vals))
    clim_max = clim_max == 0.0 ? 1.0 : clim_max
    clims = (-clim_max, clim_max)

    fig = Figure(size = (280 * n_cols + 80, 280 * n_rows + 30))

    for (ri, (label, field)) in enumerate([("native_richness_change", :native_richness_change), ("invasive_richness_change", :invasive_richness_change)])
        for (ci, pass_name) in enumerate(active_pass)
            M, temps, costs = get_metric_matrix(metrics, pass_name; field=field)

            ax = Axis(fig[ri, ci];
                title = ri == 1 ? get(pass_labels, pass_name, pass_name) : "",
                xlabel = ri == n_rows ? "ΔT (°C)" : "",
                ylabel = ci == 1 ? "Upstream Cost" : "",
                xticks = (1:length(temps), string.(temps)),
                yticks = (1:length(costs), string.(costs)),
                xticklabelrotation = 0,
            )

            hm = heatmap!(ax, 1:length(temps), 1:length(costs), M;
                colormap = diverging_cmap,
                colorrange = clims)

            for ui in 1:size(M, 1), ti in 1:size(M, 2)
                if !isnan(M[ui, ti])
                    text!(ax, ti, ui; text = string(round(M[ui, ti]; digits=2)),
                        color = abs(M[ui, ti]) > clim_max * 0.5 ? :white : :black,
                        fontsize = 9, align = (:center, :center))
                end
            end

            if ri == n_rows && ci == n_cols
                Colorbar(fig[1:n_rows, n_cols + 1], hm; label = "Δ Richness", vertical = true, width = 25)
            end
        end
    end

    title_text = "Δ Species Richness: Native (top) vs Invasive (bottom)"
    if !isempty(title_suffix)
        title_text *= "  [$title_suffix]"
    end
    Label(fig[0, 1:n_cols], title_text, fontsize = 15, font = :bold)

    save_figure(fig, joinpath(output_dir, "richness_change_heatmaps.png"); size = (280 * n_cols + 80, 280 * n_rows + 30))
end

# =============================================================================
# --- Figure 2: Final-State Native vs Invasive Richness Scatter ---
# =============================================================================

function plot_native_vs_invasive_scatter(metrics, output_dir; title_suffix="")
    println("Generating native vs invasive richness scatter...")

    markers = Dict("baseline" => :circle,
                   "improved_passability" => :utriangle,
                   "reduced_passability" => :dtriangle,
                   "blocked" => :rect)
    colors_map = Dict("baseline" => :dodgerblue,
                      "improved_passability" => :darkgreen,
                      "reduced_passability" => :orange,
                      "blocked" => :firebrick)

    active_pass = active_passability_scenarios(metrics)
    if isempty(active_pass)
        println("  No data — skipping scatter.")
        return
    end

    fig = Figure(size = (900, 750))

    ax = Axis(fig[1, 1];
        title = "Final Native vs Invasive Richness",
        xlabel = "Mean Native Richness (Year $SIMULATION_YEARS)",
        ylabel = "Mean Invasive Richness (Year $SIMULATION_YEARS)")

    for pass_name in active_pass
        x_vals = Float64[getfield(m, :native_richness_T) for m in metrics if m.pass == pass_name]
        y_vals = Float64[getfield(m, :invasive_richness_T) for m in metrics if m.pass == pass_name]
        if !isempty(x_vals)
            scatter!(ax, x_vals, y_vals;
                marker = get(markers, pass_name, :circle),
                color = get(colors_map, pass_name, :black),
                markersize = 12,
                label = get(pass_labels, pass_name, pass_name))
        end
    end

    Legend(fig[2, 1], ax; orientation = :horizontal, fontsize = 10)

    title_text = "Trade-off Between Native and Invasive Species Richness"
    if !isempty(title_suffix)
        title_text *= "  [$title_suffix]"
    end
    Label(fig[0, :], title_text, fontsize = 14, font = :bold)

    save_figure(fig, joinpath(output_dir, "native_vs_invasive_scatter.png"))
end

# =============================================================================
# --- Figure 3: Line Plots — Richness × Temperature (by passability & upstream cost) ---
# =============================================================================

function plot_sensitivity_lines(metrics, output_dir; title_suffix="")
    println("Generating sensitivity line plots...")

    active_pass = active_passability_scenarios(metrics)
    n_pass = length(active_pass)
    if n_pass == 0
        println("  No data — skipping line plots.")
        return
    end

    uc_colors = [:black, :dodgerblue, :darkorange, :crimson]
    uc_markers = [:circle, :utriangle, :dtriangle, :diamond]

    rows = ceil(Int, sqrt(n_pass))
    cols = ceil(Int, n_pass / rows)

    fig = Figure(size = (280 * cols + 60, 260 * rows + 40))
    axes_refs = []

    for (ci, pass_name) in enumerate(active_pass)
        ri = (ci - 1) ÷ cols + 1
        cj = (ci - 1) % cols + 1

        ax = Axis(fig[ri, cj];
            title = get(pass_labels, pass_name, pass_name),
            xlabel = ri == rows ? "ΔT (°C)" : "",
            ylabel = cj == 1 ? "Mean Richness" : "",
            xticks = (0:3, string.(0:3)))

        sub_metrics = [m for m in metrics if m.pass == pass_name]

        for (ui, uc) in enumerate(upstream_costs)
            xs = Float64[]
            ys_nat = Float64[]
            ys_inv = Float64[]
            for dt in temperature_increases
                match = findfirst(m -> m.dt == dt && m.uc == uc, sub_metrics)
                if match !== nothing
                    push!(xs, dt)
                    push!(ys_nat, sub_metrics[match].native_richness_T)
                    push!(ys_inv, sub_metrics[match].invasive_richness_T)
                end
            end
            if !isempty(xs)
                lines!(ax, xs, ys_nat; color = uc_colors[ui], linewidth = 2, linestyle = :solid)
                scatter!(ax, xs, ys_nat; color = uc_colors[ui], marker = uc_markers[ui],
                    markersize = 10, label = "uc=$(uc) (native)")
                lines!(ax, xs, ys_inv; color = uc_colors[ui], linewidth = 2, linestyle = :dash)
                scatter!(ax, xs, ys_inv; color = uc_colors[ui], marker = uc_markers[ui],
                    markersize = 10, label = "uc=$(uc) (invasive)")
            end
        end

        push!(axes_refs, ax)
    end

    Legend(fig[rows + 1, :], axes_refs[1]; orientation = :horizontal, fontsize = 8, nbanks = 2)

    title_text = "Native (solid) & Invasive (dashed) Richness vs Temperature — by Passability & Upstream Cost"
    if !isempty(title_suffix)
        title_text *= "  [$title_suffix]"
    end
    Label(fig[0, :], title_text, fontsize = 13, font = :bold)

    save_figure(fig, joinpath(output_dir, "sensitivity_lines.png"))
end

# =============================================================================
# --- Figure 4: Time Series for Key Scenarios ---
# =============================================================================

function load_timeseries(base_dir, path_prefix, dt, uc, pass_name)
    jld2_file = jld2_path(base_dir, path_prefix, dt, uc, pass_name)
    jldopen(jld2_file, "r") do f
        return f["sol_t"], f["sol_u"], f["species"], f["sites"]
    end
end

function compute_richness_timeseries(sol_t, sol_u, species, sites, species_indices; threshold=0.1)
    n_sites = length(sites)
    n_species = length(species)
    n_t = length(sol_t)
    richness = zeros(Float64, n_t)
    for t in 1:n_t
        u_mat = reshape(sol_u[t], n_sites, n_species)
        site_vals = Float64[sum(u_mat[i, species_indices] .> threshold) for i in 1:n_sites]
        richness[t] = mean(site_vals)
    end
    return richness
end

function compute_biomass_timeseries(sol_t, sol_u, species, sites)
    n_sites = length(sites)
    n_species = length(species)
    n_t = length(sol_t)
    biomass = zeros(Float64, n_t)
    for t in 1:n_t
        u_mat = reshape(sol_u[t], n_sites, n_species)
        biomass[t] = sum(u_mat)
    end
    return biomass
end

function plot_timeseries_comparison(metrics, base_dir, output_dir; path_prefix="", title_suffix="")
    println("Generating time series comparison...")

    active_pass = active_passability_scenarios(metrics)
    if isempty(active_pass)
        println("  No data — skipping timeseries.")
        return
    end

    # Use uc values actually present in the loaded metrics
    actual_ucs = unique([m.uc for m in metrics])
    ref_uc = 0.05 in actual_ucs ? 0.05 : actual_ucs[1]

    # Use temperature values actually present
    actual_temps = unique([m.dt for m in metrics])
    ref_dt = temperature_increases[1] in actual_temps ? temperature_increases[1] : actual_temps[1]

    # Check that the timeseries files exist for the reference combination
    test_file = jld2_path(base_dir, path_prefix, ref_dt, ref_uc, active_pass[1])
    if !isfile(test_file)
        println("  Timeseries files not found — skipping timeseries comparison.")
        return
    end

    fig = Figure(size = (1650, 800))

    temp_colors = [:dodgerblue, :darkorange, :crimson, :darkviolet]
    pass_colors_default = Dict(zip(["baseline", "improved_passability", "reduced_passability", "blocked"],
                                    [:dodgerblue, :darkgreen, :orange, :crimson]))

    # Panel A: Temperature effect on biomass & richness (baseline passability or first active)
    ref_pass = "baseline" in active_pass ? "baseline" : active_pass[1]
    uc = ref_uc
    ax_bio = Axis(fig[1, 1];
        title = "Total Biomass — Temperature Effect ($ref_pass, uc=$uc)",
        xlabel = "Time (years)", ylabel = "Total Biomass")
    ax_rich = Axis(fig[2, 1];
        title = "Species Richness — Temperature Effect ($ref_pass, uc=$uc)",
        xlabel = "Time (years)", ylabel = "Mean Richness")

    for (ti, dt) in enumerate(temperature_increases)
        ts_file = jld2_path(base_dir, path_prefix, dt, uc, ref_pass)
        if !isfile(ts_file)
            continue
        end
        sol_t, sol_u, species, sites = load_timeseries(base_dir, path_prefix, dt, uc, ref_pass)
        time_years = sol_t ./ DAYS_PER_YEAR

        native_idx = classify_species_indices(species, NATIVE_SPECIES)
        invasive_idx = classify_species_indices(species, INVASIVE_SPECIES)

        bio = compute_biomass_timeseries(sol_t, sol_u, species, sites)
        nat_rich = compute_richness_timeseries(sol_t, sol_u, species, sites, native_idx)
        inv_rich = compute_richness_timeseries(sol_t, sol_u, species, sites, invasive_idx)

        lines!(ax_bio, time_years, bio; color = temp_colors[ti], linewidth = 2,
            label = "ΔT=$(dt)°C")

        lines!(ax_rich, time_years, nat_rich; color = temp_colors[ti], linewidth = 2,
            linestyle = :solid, label = "ΔT=$(dt)°C (native)")
        lines!(ax_rich, time_years, inv_rich; color = temp_colors[ti], linewidth = 2,
            linestyle = :dash, label = "ΔT=$(dt)°C (invasive)")
    end

    legend_col = GridLayout(fig[1:2, 3])
    Legend(legend_col[1, 1], ax_bio; fontsize = 8, tellheight = false)
    Legend(legend_col[2, 1], ax_rich; fontsize = 7, nbanks = 2, tellheight = false)

    # Panel B: Passability effect on biomass & richness
    dt = 1.0
    if !(dt in temperature_increases)
        dt = temperature_increases[1]
    end
    ax_bio2 = Axis(fig[1, 2];
        title = "Total Biomass — Passability Effect (ΔT=$(dt)°C, uc=$uc)",
        xlabel = "Time (years)", ylabel = "Total Biomass")
    ax_rich2 = Axis(fig[2, 2];
        title = "Species Richness — Passability Effect (ΔT=$(dt)°C, uc=$uc)",
        xlabel = "Time (years)", ylabel = "Mean Richness")

    for pass_name in active_pass
        ts_file = jld2_path(base_dir, path_prefix, dt, uc, pass_name)
        if !isfile(ts_file)
            continue
        end
        sol_t, sol_u, species, sites = load_timeseries(base_dir, path_prefix, dt, uc, pass_name)
        time_years = sol_t ./ DAYS_PER_YEAR

        native_idx = classify_species_indices(species, NATIVE_SPECIES)
        invasive_idx = classify_species_indices(species, INVASIVE_SPECIES)

        bio = compute_biomass_timeseries(sol_t, sol_u, species, sites)
        nat_rich = compute_richness_timeseries(sol_t, sol_u, species, sites, native_idx)
        inv_rich = compute_richness_timeseries(sol_t, sol_u, species, sites, invasive_idx)

        lines!(ax_bio2, time_years, bio; color = get(pass_colors_default, pass_name, :black), linewidth = 2,
            label = get(pass_labels, pass_name, pass_name))

        lines!(ax_rich2, time_years, nat_rich; color = get(pass_colors_default, pass_name, :black), linewidth = 2,
            linestyle = :solid, label = "$(get(pass_labels, pass_name, pass_name)) (native)")
        lines!(ax_rich2, time_years, inv_rich; color = get(pass_colors_default, pass_name, :black), linewidth = 2,
            linestyle = :dash, label = "$(get(pass_labels, pass_name, pass_name)) (invasive)")
    end

    Legend(legend_col[3, 1], ax_bio2; fontsize = 8, tellheight = false)
    Legend(legend_col[4, 1], ax_rich2; fontsize = 7, nbanks = 2, tellheight = false)

    title_text = "Time Series Comparison — Temperature & Passability Effects"
    if !isempty(title_suffix)
        title_text *= "  [$title_suffix]"
    end
    Label(fig[0, :], title_text, fontsize = 15, font = :bold)

    save_figure(fig, joinpath(output_dir, "timeseries_comparison.png"))
end

# =============================================================================
# --- Figure 5: Biomass Change Bar Chart ---
# =============================================================================

function plot_biomass_change_bars(metrics, output_dir; title_suffix="")
    println("Generating biomass change bar chart...")

    active_pass = active_passability_scenarios(metrics)
    n_pass = length(active_pass)
    if n_pass == 0
        println("  No data — skipping bar chart.")
        return
    end

    fig = Figure(size = (220 * n_pass + 220, 680))

    ax = Axis(fig[1, 1];
        title = "Total Biomass Change (Year 0 → Year $SIMULATION_YEARS) by Passability Scenario",
        xlabel = "Passability Scenario",
        ylabel = "Δ Total Biomass",
        xticklabelrotation = 0.3)

    bar_width = 0.18
    offsets = range(-1.5 * bar_width, 1.5 * bar_width; length=4)
    pass_colors_bar = [:dodgerblue, :darkgreen, :orange, :crimson]

    for (ci, pass_name) in enumerate(active_pass)
        vals = [getfield(m, :biomass_change) for m in metrics if m.pass == pass_name]
        if isempty(vals)
            continue
        end
        color_idx = findfirst(==(pass_name), passability_scenarios)
        color_idx = color_idx === nothing ? 1 : color_idx
        xpos = ci .+ offsets[ci]
        barplot!(ax, fill(xpos, length(vals)), vals;
            color = (pass_colors_bar[mod1(color_idx, 4)], 0.6),
            width = bar_width,
            label = get(pass_labels, pass_name, pass_name),
            strokecolor = pass_colors_bar[mod1(color_idx, 4)],
            strokewidth = 1)

        mean_val = mean(vals)
        scatter!(ax, [ci], [mean_val];
            color = pass_colors_bar[mod1(color_idx, 4)],
            marker = :diamond,
            markersize = 14,
            strokecolor = :black,
            strokewidth = 1)
    end

    ax.xticks = (1:n_pass, [get(pass_labels, p, p) for p in active_pass])

    Legend(fig[2, 1], ax; orientation = :horizontal, fontsize = 9)

    title_text = "Biomass Change Across All Parameter Combinations"
    if !isempty(title_suffix)
        title_text *= "  [$title_suffix]"
    end
    Label(fig[0, :], title_text, fontsize = 14, font = :bold)

    save_figure(fig, joinpath(output_dir, "biomass_change_bars.png"))
end

# =============================================================================
# --- Run All Plots For a Given Metrics Collection ---
# =============================================================================

function run_all_plots(metrics, base_dir, output_dir; path_prefix="", title_suffix="")
    mkpath(output_dir)
    if isempty(metrics)
        println("No metrics to plot — skipping.")
        return
    end
    plot_richness_change_heatmaps(metrics, output_dir; title_suffix=title_suffix)
    plot_native_vs_invasive_scatter(metrics, output_dir; title_suffix=title_suffix)
    plot_sensitivity_lines(metrics, output_dir; title_suffix=title_suffix)
    plot_timeseries_comparison(metrics, base_dir, output_dir; path_prefix=path_prefix, title_suffix=title_suffix)
    plot_biomass_change_bars(metrics, output_dir; title_suffix=title_suffix)
    println("  Plots saved to: $output_dir")
end

# =============================================================================
# --- Main ---
# =============================================================================

println("="^60)
println("Sensitivity Summary Plots")
println("  Base directory: $BASE_DIR")
println("="^60)

if has_matrices
    for matrix_name in matrix_names
        println("\n--- Matrix: $matrix_name ---")
        local metrics = load_all_metrics(BASE_DIR; path_prefix=matrix_name)
        local out_dir = joinpath(BASE_DIR, matrix_name, "summary_plots")
        run_all_plots(metrics, BASE_DIR, out_dir; path_prefix=matrix_name, title_suffix=matrix_name)
    end
else
    metrics = load_all_metrics(BASE_DIR)
    out_dir = joinpath(BASE_DIR, "summary_plots")
    run_all_plots(metrics, BASE_DIR, out_dir)
end

println("\nAll summary plots complete.")
println("="^60)
