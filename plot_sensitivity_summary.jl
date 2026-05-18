using Pkg; Pkg.activate(".")
using JLD2
using CairoMakie
using Statistics
using Guadex

const DAYS_PER_YEAR = 365
const SIMULATION_YEARS = 3
const NATIVE_SPECIES = ["AB", "AH", "SP", "PW", "LS", "SA", "IL", "CP", "IO"]
const INVASIVE_SPECIES = ["GH", "MS", "LG", "CC", "CG", "AM", "OM", "EL", "GL", "TT"]

const BASE_DIR = "results/sensitivity_temp_passability"
const OUTPUT_DIR = "results/sensitivity_summary_plots"
mkpath(OUTPUT_DIR)

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

# =============================================================================
# --- Load Scalar Metrics From All Runs ---
# =============================================================================

function load_all_metrics()
    println("Loading metrics from all simulation runs...")
    metrics = []

    for dt in temperature_increases
        for uc in upstream_costs
            for pass_name in passability_scenarios
                run_dir = joinpath(BASE_DIR, "dT_$(dt)C", "uc_$(uc)", pass_name)
                jld2_file = joinpath(run_dir, "simulation_output.jld2")

                if !isfile(jld2_file)
                    println("  WARNING: Missing $jld2_file")
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

# =============================================================================
# --- Figure 1: Heatmap Grid — ΔRichness: Native & Invasive ---
# =============================================================================

function plot_richness_change_heatmaps(metrics)
    println("Generating richness change heatmaps...")

    diverging_cmap = cgrad([RGBf(0.129, 0.4, 0.675), RGBf(0.97, 0.97, 0.97), RGBf(0.698, 0.094, 0.169)])

    n_pass = length(passability_scenarios)
    n_rows = 2
    n_cols = n_pass

    all_native = Float64[getfield(m, :native_richness_change) for m in metrics]
    all_invasive = Float64[getfield(m, :invasive_richness_change) for m in metrics]
    all_vals = vcat(all_native, all_invasive)
    clim_max = maximum(abs.(all_vals))
    clim_max = clim_max == 0.0 ? 1.0 : clim_max
    clims = (-clim_max, clim_max)

    fig = Figure(size = (280 * n_cols + 80, 280 * n_rows + 30))

    row_labels = ["Native ΔRichness", "Invasive ΔRichness"]

    for (ri, (label, field)) in enumerate([("native_richness_change", :native_richness_change), ("invasive_richness_change", :invasive_richness_change)])
        for (ci, pass_name) in enumerate(passability_scenarios)
            M, temps, costs = get_metric_matrix(metrics, pass_name; field=field)

            ax = Axis(fig[ri, ci];
                title = ri == 1 ? pass_labels[pass_name] : "",
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

    Label(fig[0, 1:n_cols], "Δ Species Richness: Native (top) vs Invasive (bottom)", fontsize = 15, font = :bold)

    save_figure(fig, joinpath(OUTPUT_DIR, "richness_change_heatmaps.png"); size = (280 * n_cols + 80, 280 * n_rows + 30))
end

# =============================================================================
# --- Figure 2: Final-State Native vs Invasive Richness Scatter ---
# =============================================================================

function plot_native_vs_invasive_scatter(metrics)
    println("Generating native vs invasive richness scatter...")

    markers = Dict("baseline" => :circle,
                   "improved_passability" => :utriangle,
                   "reduced_passability" => :dtriangle,
                   "blocked" => :rect)
    colors_map = Dict("baseline" => :dodgerblue,
                      "improved_passability" => :darkgreen,
                      "reduced_passability" => :orange,
                      "blocked" => :firebrick)

    fig = Figure(size = (900, 700))

    ax = Axis(fig[1, 1];
        title = "Final Native vs Invasive Richness",
        xlabel = "Mean Native Richness (Year 3)",
        ylabel = "Mean Invasive Richness (Year 3)")

    for pass_name in passability_scenarios
        x_vals = Float64[getfield(m, :native_richness_T) for m in metrics if m.pass == pass_name]
        y_vals = Float64[getfield(m, :invasive_richness_T) for m in metrics if m.pass == pass_name]
        scatter!(ax, x_vals, y_vals;
            marker = markers[pass_name],
            color = colors_map[pass_name],
            markersize = 12,
            label = pass_labels[pass_name])
    end

    axislegend(ax; position = :rb, fontsize = 10)

    Label(fig[0, :], "Trade-off Between Native and Invasive Species Richness", fontsize = 14, font = :bold)

    save_figure(fig, joinpath(OUTPUT_DIR, "native_vs_invasive_scatter.png"))
end

# =============================================================================
# --- Figure 3: Line Plots — Richness × Temperature (by passability & upstream cost) ---
# =============================================================================

function plot_sensitivity_lines(metrics)
    println("Generating sensitivity line plots...")

    fig = Figure(size = (1000, 700))

    n_uc = length(upstream_costs)
    uc_colors = [:black, :dodgerblue, :darkorange, :crimson]
    uc_markers = [:circle, :utriangle, :dtriangle, :diamond]

    pass_titles = pass_labels
    rows = 2
    cols = 2

    for (ci, pass_name) in enumerate(passability_scenarios)
        ri = (ci - 1) ÷ cols + 1
        cj = (ci - 1) % cols + 1

        ax_nat = Axis(fig[ri, cj];
            title = pass_titles[pass_name],
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
            lines!(ax_nat, xs, ys_nat; color = uc_colors[ui], linewidth = 2, linestyle = :solid)
            scatter!(ax_nat, xs, ys_nat; color = uc_colors[ui], marker = uc_markers[ui],
                markersize = 10, label = "uc=$(upstream_costs[ui]) (native)")
            lines!(ax_nat, xs, ys_inv; color = uc_colors[ui], linewidth = 2, linestyle = :dash)
            scatter!(ax_nat, xs, ys_inv; color = uc_colors[ui], marker = uc_markers[ui],
                markersize = 10, label = "uc=$(upstream_costs[ui]) (invasive)")
        end

        if ri == 1 && cj == 1
            axislegend(ax_nat; position = :lt, fontsize = 7, nbanks = 2)
        end
    end

    Label(fig[0, :],
        "Native (solid) & Invasive (dashed) Richness vs Temperature — by Passability & Upstream Cost",
        fontsize = 13, font = :bold)

    save_figure(fig, joinpath(OUTPUT_DIR, "sensitivity_lines.png"))
end

# =============================================================================
# --- Figure 4: Time Series for Key Scenarios ---
# =============================================================================

function load_timeseries(dt, uc, pass_name)
    run_dir = joinpath(BASE_DIR, "dT_$(dt)C", "uc_$(uc)", pass_name)
    jld2_file = joinpath(run_dir, "simulation_output.jld2")
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

function plot_timeseries_comparison(metrics)
    println("Generating time series comparison...")

    ref_uc = 0.05

    fig = Figure(size = (1400, 800))

    pass_names_all = ["baseline", "improved_passability", "reduced_passability", "blocked"]
    temp_colors = [:dodgerblue, :darkorange, :crimson, :darkviolet]
    pass_colors = Dict(zip(pass_names_all, [:dodgerblue, :darkgreen, :orange, :crimson]))

    # Panel A: Temperature effect on biomass & richness (baseline passability)
    pass_name = "baseline"
    uc = ref_uc
    ax_bio = Axis(fig[1, 1];
        title = "Total Biomass — Temperature Effect (baseline, uc=$uc)",
        xlabel = "Time (years)", ylabel = "Total Biomass")
    ax_rich = Axis(fig[2, 1];
        title = "Species Richness — Temperature Effect (baseline, uc=$uc)",
        xlabel = "Time (years)", ylabel = "Mean Richness")

    for (ti, dt) in enumerate([0.0, 1.0, 2.0, 3.0])
        sol_t, sol_u, species, sites = load_timeseries(dt, uc, pass_name)
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

    axislegend(ax_bio; position = :rb, fontsize = 8)
    axislegend(ax_rich; position = :rb, fontsize = 7, nbanks = 2)

    # Panel B: Passability effect on biomass & richness (ΔT=1°C, uc=0.05)
    dt = 1.0
    ax_bio2 = Axis(fig[1, 2];
        title = "Total Biomass — Passability Effect (ΔT=$(dt)°C, uc=$uc)",
        xlabel = "Time (years)", ylabel = "Total Biomass")
    ax_rich2 = Axis(fig[2, 2];
        title = "Species Richness — Passability Effect (ΔT=$(dt)°C, uc=$uc)",
        xlabel = "Time (years)", ylabel = "Mean Richness")

    for pass_name in pass_names_all
        sol_t, sol_u, species, sites = load_timeseries(dt, uc, pass_name)
        time_years = sol_t ./ DAYS_PER_YEAR

        native_idx = classify_species_indices(species, NATIVE_SPECIES)
        invasive_idx = classify_species_indices(species, INVASIVE_SPECIES)

        bio = compute_biomass_timeseries(sol_t, sol_u, species, sites)
        nat_rich = compute_richness_timeseries(sol_t, sol_u, species, sites, native_idx)
        inv_rich = compute_richness_timeseries(sol_t, sol_u, species, sites, invasive_idx)

        lines!(ax_bio2, time_years, bio; color = pass_colors[pass_name], linewidth = 2,
            label = pass_labels[pass_name])

        lines!(ax_rich2, time_years, nat_rich; color = pass_colors[pass_name], linewidth = 2,
            linestyle = :solid, label = "$(pass_labels[pass_name]) (native)")
        lines!(ax_rich2, time_years, inv_rich; color = pass_colors[pass_name], linewidth = 2,
            linestyle = :dash, label = "$(pass_labels[pass_name]) (invasive)")
    end

    axislegend(ax_bio2; position = :rb, fontsize = 8)
    axislegend(ax_rich2; position = :rb, fontsize = 7, nbanks = 2)

    Label(fig[0, :], "Time Series Comparison — Temperature & Passability Effects",
        fontsize = 15, font = :bold)

    save_figure(fig, joinpath(OUTPUT_DIR, "timeseries_comparison.png"))
end

# =============================================================================
# --- Figure 5: Biomass Change Bar Chart ---
# =============================================================================

function plot_biomass_change_bars(metrics)
    println("Generating biomass change bar chart...")

    fig = Figure(size = (1100, 600))

    akw = Axis(fig[1, 1];
        title = "Total Biomass Change (Year 0 → Year 3) by Passability Scenario",
        xlabel = "Passability Scenario",
        ylabel = "Δ Total Biomass",
        xticklabelrotation = 0.3)

    n_pass = length(passability_scenarios)
    n_combos = length(temperature_increases) * length(upstream_costs)
    bar_width = 0.18
    offsets = range(-1.5 * bar_width, 1.5 * bar_width; length=4)

    pass_colors_bar = [:dodgerblue, :darkgreen, :orange, :crimson]

    for (ci, pass_name) in enumerate(passability_scenarios)
        vals = [getfield(m, :biomass_change) for m in metrics if m.pass == pass_name]
        xpos = ci .+ offsets[ci]
        barplot!(akw, fill(xpos, length(vals)), vals;
            color = (pass_colors_bar[ci], 0.6),
            width = bar_width,
            label = pass_labels[pass_name],
            strokecolor = pass_colors_bar[ci],
            strokewidth = 1)

        mean_val = mean(vals)
        scatter!(akw, [ci], [mean_val];
            color = pass_colors_bar[ci],
            marker = :diamond,
            markersize = 14,
            strokecolor = :black,
            strokewidth = 1)
    end

    akw.xticks = (1:n_pass, [pass_labels[p] for p in passability_scenarios])

    axislegend(akw; position = :lt, fontsize = 9)

    Label(fig[0, :], "Biomass Change Across All Parameter Combinations", fontsize = 14, font = :bold)

    save_figure(fig, joinpath(OUTPUT_DIR, "biomass_change_bars.png"))
end

# =============================================================================
# --- Main ---
# =============================================================================

println("="^60)
println("Sensitivity Summary Plots")
println("="^60)

metrics = load_all_metrics()

plot_richness_change_heatmaps(metrics)
plot_native_vs_invasive_scatter(metrics)
plot_sensitivity_lines(metrics)
plot_timeseries_comparison(metrics)
plot_biomass_change_bars(metrics)

println("\nAll summary plots saved to: $OUTPUT_DIR")
println("="^60)
