# =============================================================================
# alt_interactions_comparison.jl
#
# Generates comparison plots across the three alternative interaction matrix
# setups from results/sensitivity_alt_interactions/.
#
# Usage:
#   julia --project=. scripts/alt_interactions_comparison.jl
#
#   Optional arguments:
#     arg1 (Float64): upstream_cost (default 0.01; also plots 0.5)
#     arg2 (String):  passability scenario (default "baseline")
#
#   Examples:
#     julia --project=. scripts/alt_interactions_comparison.jl
#     julia --project=. scripts/alt_interactions_comparison.jl 0.05 blocked
#
# Output:
#   results/sensitivity_alt_interactions/comparison_plots/
#
# Plot 1: Multi-panel figure showing overall native richness retention (S/S₀)
#         vs ΔT across all uc × passability combinations. Each panel has three
#         lines (original/random/invasive-favoured) with markers coloured by ΔT.
#
# Plot 2: Per-subcatchment heatmaps (S/S₀ by subcatchment × ΔT) for the
#         specified uc × passability combination, one panel per matrix type.
#
# Metric:
#   S/S₀ = mean native species richness at end of year 3 ÷ mean native species
#   richness at simulation start (t=0), computed per-site and averaged across
#   sites where at least one native species was present at t=0.
# =============================================================================

using Pkg; Pkg.activate(".")
using JLD2
using DataFrames
using Guadex
using CairoMakie
using Statistics
using Printf

const DAYS_PER_YEAR = 365
const FINAL_YEAR = 3
const RICHNESS_THRESHOLD = 0.1

const BASE_RESULTS_DIR = "results/sensitivity_alt_interactions"
const OUTPUT_DIR = joinpath(BASE_RESULTS_DIR, "comparison_plots")

const NATIVE_SPECIES_CODES = ["AB", "AH", "SP", "PW", "LS", "SA", "IL", "CP", "IO"]

const MATRIX_NAMES = ["original", "random", "invasive_favoring"]
const MATRIX_DISPLAY = Dict(
    "original"           => "Original",
    "random"             => "Random",
    "invasive_favoring"  => "Invasive-Favoured",
)
const MATRIX_COLORS = Dict(
    "original"           => :royalblue,
    "random"             => :darkorange,
    "invasive_favoring"  => :firebrick,
)
const MATRIX_MARKERS = Dict(
    "original"           => :circle,
    "random"             => :rect,
    "invasive_favoring"  => :utriangle,
)

const TEMPERATURE_INCREASES = [0.0, 1.0, 2.0, 3.0]
const DT_CMAP = :YlOrRd
const UPSTREAM_COSTS = [0.01, 0.5]
const PASS_SCENARIOS = ["baseline", "reduced_passability", "blocked"]

# =============================================================================
# Data Loading & Computation
# =============================================================================

function classify_species_indices(all_species, target_group)
    return [findfirst(==(sp), all_species) for sp in target_group if sp in all_species]
end

function find_time_index(sol_t, target_time)
    for (i, t) in enumerate(sol_t)
        if t >= target_time
            return i
        end
    end
    return length(sol_t)
end

function compute_native_richness_per_site(u_vec, n_sites, n_species, native_idx)
    u_mat = reshape(u_vec, n_sites, n_species)
    richness = zeros(Float64, n_sites)
    for i in 1:n_sites
        richness[i] = sum(u_mat[i, native_idx] .> RICHNESS_THRESHOLD)
    end
    return richness
end

function load_single_result(matrix_name, dT, uc, pass_name)
    jld2_path = joinpath(BASE_RESULTS_DIR, matrix_name, "dT_$(dT)C", "uc_$(uc)", pass_name, "simulation_output.jld2")
    if !isfile(jld2_path)
        return nothing
    end
    data = JLD2.load(jld2_path)
    sol_t = data["sol_t"]
    sol_u = data["sol_u"]
    species = data["species"]
    sites = data["sites"]
    n_sites = length(sites)
    n_species = length(species)
    native_idx = classify_species_indices(species, NATIVE_SPECIES_CODES)
    t_start = 1
    t_end = find_time_index(sol_t, FINAL_YEAR * DAYS_PER_YEAR)
    richness_start = compute_native_richness_per_site(sol_u[t_start], n_sites, n_species, native_idx)
    richness_end   = compute_native_richness_per_site(sol_u[t_end],   n_sites, n_species, native_idx)
    return (sites=sites, species=species, richness_start=richness_start, richness_end=richness_end,
            n_sites=n_sites, n_species=n_species)
end

function compute_mean_end_richness(r)
    valid = r.richness_start .> 0.0
    if !any(valid)
        return NaN
    end
    return mean(r.richness_end[valid])
end

function build_subcatchment_map(site_df, sites)
    sc_set = Set{String}()
    site_to_sc = Dict{String, String}()
    filtered = filter(row -> row.CODIGO in sites, site_df)
    for row in eachrow(filtered)
        sc = string(row.CODIGO_S)
        push!(sc_set, sc)
        site_to_sc[string(row.CODIGO)] = sc
    end
    return site_to_sc, sort(collect(sc_set))
end

# =============================================================================
# Plot 1: Overall comparison — multi-panel (uc × pass)
# =============================================================================

function plot_overall_comparison_grid(all_data)
    n_uc = length(UPSTREAM_COSTS)
    n_pass = length(PASS_SCENARIOS)

    fig = Figure(size = (380 * n_pass + 80, 320 * n_uc + 60))

    for (ui, uc) in enumerate(UPSTREAM_COSTS)
        for (pi, pass_name) in enumerate(PASS_SCENARIOS)
            key = (uc, pass_name)
            if !haskey(all_data, key)
                continue
            end
            data = all_data[key]

            ax = Axis(fig[ui, pi],
                xlabel = ui == n_uc ? "ΔT (°C)" : "",
                ylabel = pi == 1 ? "S / S₀ (ref: Original ΔT=0)" : "",
                xticks = (0:1:3, ["0", "1", "2", "3"]),
                yminorticksvisible = true,
                yminorticks = IntervalsBetween(2),
                title = "uc=$(uc)  $(pass_name)",
                titlefont = :regular,
            )
            xlims!(ax, -0.3, 3.3)

            ref_val = NaN
            if haskey(data, "original") && haskey(data["original"], 0.0)
                ref_val = compute_mean_end_richness(data["original"][0.0])
            end

            for matrix_name in MATRIX_NAMES
                if !haskey(data, matrix_name)
                    continue
                end
                ratios = Float64[]
                dts = Float64[]
                marker_colors = RGBAf[]
                for dT in TEMPERATURE_INCREASES
                    if !haskey(data[matrix_name], dT)
                        continue
                    end
                    end_mean = compute_mean_end_richness(data[matrix_name][dT])
                    if isnan(end_mean) || isnan(ref_val) || ref_val <= 0.0
                        continue
                    end
                    push!(ratios, end_mean / ref_val)
                    push!(dts, dT)
                    frac = (dT - minimum(TEMPERATURE_INCREASES)) / max(1e-9, maximum(TEMPERATURE_INCREASES) - minimum(TEMPERATURE_INCREASES))
                    push!(marker_colors, cgrad(DT_CMAP)[frac])
                end
                if length(dts) < 2
                    continue
                end
                lines!(ax, dts, ratios, color = MATRIX_COLORS[matrix_name], linewidth = 2,
                       label = MATRIX_DISPLAY[matrix_name])
                scatter!(ax, dts, ratios, color = marker_colors, marker = MATRIX_MARKERS[matrix_name],
                         markersize = 10, strokecolor = MATRIX_COLORS[matrix_name], strokewidth = 1.2)
            end
        end
    end

    legend_axis = Axis(fig[1:n_uc, n_pass + 1], tellwidth = false, tellheight = false)
    hidedecorations!(legend_axis)
    hidespines!(legend_axis)
    elements = [LineElement(color = MATRIX_COLORS[m], linewidth = 2) for m in MATRIX_NAMES]
    Legend(fig[1, n_pass + 1], elements, [MATRIX_DISPLAY[m] for m in MATRIX_NAMES],
           "Matrix Type", orientation = :vertical, framevisible = false,
           tellwidth = false, tellheight = false)

    Colorbar(fig[2:3, n_pass + 1], colormap = DT_CMAP, limits = (0, 3),
             label = "ΔT (°C)", ticks = 0:1:3, vertical = true, width = 25,
             tellwidth = false, tellheight = true)

    Label(fig[0, :], "Native Richness Retention (S/S₀) vs Warming Across Parameter Combinations",
          fontsize = 14, font = :bold)
    Label(fig[n_uc + 1, 1:n_pass],
          "Climate Warming Stress Gradient (ΔT in °C)",
          fontsize = 11)

    out_path = joinpath(OUTPUT_DIR, "comparison_overall_grid.png")
    save(out_path, fig, size = (380 * n_pass + 80, 320 * n_uc + 60))
    println("Saved: $out_path")
    return fig
end

# =============================================================================
# Plot 2: Per-subcatchment heatmap — single (uc, pass) combo
# =============================================================================

function plot_subcatchment_heatmaps(all_data, all_sc_data, site_to_sc, uc, pass_name)
    key = (uc, pass_name)
    if !haskey(all_data, key) || !haskey(all_sc_data, key)
        println("No data for uc=$(uc), pass=$(pass_name)")
        return nothing
    end
    data = all_data[key]
    sc_data = all_sc_data[key]

    all_scs_all = Set{String}()
    for matrix_name in MATRIX_NAMES
        if haskey(sc_data, matrix_name)
            union!(all_scs_all, keys(sc_data[matrix_name]))
        end
    end
    all_scs = sort(collect(all_scs_all))
    n_sc = length(all_scs)
    if n_sc == 0
        println("No subcatchment data available.")
        return nothing
    end

    sc_labels = String[]
    tick_positions = Int[]
    for (i, sc) in enumerate(all_scs)
        if i % 5 == 1 || i == n_sc
            push!(sc_labels, sc)
        else
            push!(sc_labels, "")
        end
        push!(tick_positions, i)
    end

    fig = Figure(size = (1600, 500))
    clims = (0.2, 2.0)

    for (mi, matrix_name) in enumerate(MATRIX_NAMES)
        if !haskey(data, matrix_name) || !haskey(sc_data, matrix_name)
            continue
        end
        heatmap_data = fill(NaN, n_sc, length(TEMPERATURE_INCREASES))
        for (ri, sc) in enumerate(all_scs)
            if haskey(sc_data[matrix_name], sc)
                for (ci, dT) in enumerate(TEMPERATURE_INCREASES)
                    if haskey(sc_data[matrix_name][sc], dT)
                        heatmap_data[ri, ci] = sc_data[matrix_name][sc][dT]
                    end
                end
            end
        end

        ax = Axis(fig[1, mi],
            title = MATRIX_DISPLAY[matrix_name],
            titlefont = :bold,
            xlabel = mi == 2 ? "ΔT (°C)" : "",
            ylabel = mi == 1 ? "Subcatchment" : "",
            xticks = (1:length(TEMPERATURE_INCREASES), string.(Int.(TEMPERATURE_INCREASES))),
            yticks = (tick_positions, sc_labels),
            yticklabelsize = 5,
        )

        heatmap!(ax, heatmap_data', colormap = Reverse(:RdYlGn), colorrange = clims)

        for ri in 1:n_sc, ci in 1:length(TEMPERATURE_INCREASES)
            v = heatmap_data[ri, ci]
            if !isnan(v)
                text_color = (0.5 < v < 1.5) ? :black : :white
                text!(ax, ci, ri, text = @sprintf("%.2f", v),
                      color = text_color, fontsize = 6, align = (:center, :center))
            end
        end
    end

    Colorbar(fig[1, 4], colormap = Reverse(:RdYlGn), limits = clims,
             label = "S/S₀", vertical = true, width = 25, labelsize = 12)

    Label(fig[0, :], "Per-Subcatchment Native Richness Retention (uc=$(uc), $(pass_name))",
          fontsize = 14, font = :bold)

    colsize!(fig.layout, 1, Relative(0.28))
    colsize!(fig.layout, 2, Relative(0.28))
    colsize!(fig.layout, 3, Relative(0.28))
    colsize!(fig.layout, 4, Relative(0.05))
    colgap!(fig.layout, 10)

    out_path = joinpath(OUTPUT_DIR, "comparison_subcatchment_heatmap_uc$(uc)_$(pass_name).png")
    save(out_path, fig, size = (1600, 500))
    println("Saved: $out_path")
    return fig
end

# =============================================================================
# Data Collection
# =============================================================================

function collect_all_results()
    all_data = Dict{Tuple{Float64, String}, Dict{String, Dict{Float64, Any}}}()
    all_sc_data = Dict{Tuple{Float64, String}, Dict{String, Dict{String, Dict{Float64, Float64}}}}()

    for uc in UPSTREAM_COSTS
        for pass_name in PASS_SCENARIOS
            key = (uc, pass_name)
            has_any = false
            for matrix_name in MATRIX_NAMES
                jld2_path = joinpath(BASE_RESULTS_DIR, matrix_name, "dT_0.0C", "uc_$(uc)", pass_name, "simulation_output.jld2")
                if isfile(jld2_path)
                    has_any = true
                    break
                end
            end
            if !has_any
                println("  Skipping uc=$(uc), pass=$(pass_name) (no data)")
                continue
            end

            data = Dict{String, Dict{Float64, Any}}()
            sc_data_matrix = Dict{String, Dict{String, Dict{Float64, Float64}}}()
            for matrix_name in MATRIX_NAMES
                data[matrix_name] = Dict{Float64, Any}()
                sc_data_matrix[matrix_name] = Dict{String, Dict{Float64, Float64}}()
                for dT in TEMPERATURE_INCREASES
                    r = load_single_result(matrix_name, dT, uc, pass_name)
                    if r === nothing
                        continue
                    end
                    data[matrix_name][dT] = r
                end
            end
            all_data[key] = data
            all_sc_data[key] = sc_data_matrix
        end
    end
    return all_data, all_sc_data
end

function compute_subcatchment_ratios(all_data, all_sc_data, site_to_sc)
    for (key, data) in all_data
        uc, pass_name = key
        sc_data_matrix = all_sc_data[key]
        for matrix_name in MATRIX_NAMES
            if !haskey(data, matrix_name)
                continue
            end
            for dT in TEMPERATURE_INCREASES
                if !haskey(data[matrix_name], dT)
                    continue
                end
                r = data[matrix_name][dT]
                sc_ratios = Dict{String, Vector{Float64}}()
                for (si, site) in enumerate(r.sites)
                    site_str = string(site)
                    if !haskey(site_to_sc, site_str)
                        continue
                    end
                    sc = site_to_sc[site_str]
                    if !haskey(sc_ratios, sc)
                        sc_ratios[sc] = Float64[]
                    end
                    if r.richness_start[si] > 0.0
                        push!(sc_ratios[sc], r.richness_end[si] / r.richness_start[si])
                    end
                end
                for (sc, ratios) in sc_ratios
                    if !isempty(ratios)
                        sc_data_matrix[matrix_name][sc] = get(sc_data_matrix[matrix_name], sc, Dict{Float64, Float64}())
                        sc_data_matrix[matrix_name][sc][dT] = mean(ratios)
                    end
                end
            end
        end
    end
    return all_sc_data
end

# =============================================================================
# Main
# =============================================================================

function main()
    println("="^60)
    println("Alternative Interaction Matrix Comparison Plots")
    println("="^60)

    uc_target = 0.01
    pass_target = "baseline"
    if length(ARGS) >= 1
        uc_target = parse(Float64, ARGS[1])
    end
    if length(ARGS) >= 2
        pass_target = ARGS[2]
    end

    println("\nSettings: default uc = $(uc_target), passability = $(pass_target)")
    mkpath(OUTPUT_DIR)

    println("\nLoading site_df for subcatchment mapping...")
    data_base = prepare_ode_data(upstream_cost = 0.05)
    site_to_sc, all_scs = build_subcatchment_map(data_base.site_df, data_base.sites)
    println("  $(length(all_scs)) subcatchments found.")

    println("\nLoading simulation results across all uc × pass combos...")
    all_data, all_sc_data = collect_all_results()
    n_total = sum(sum(length(keys(data[m])) for m in keys(data)) for (_, data) in all_data)
    println("  Loaded $(n_total) total simulation result sets across $(length(all_data)) parameter combos.")

    println("\nComputing per-subcatchment ratios...")
    all_sc_data = compute_subcatchment_ratios(all_data, all_sc_data, site_to_sc)

    println("\n=== Plot 1: Overall comparison grid (all uc × pass) ===")
    plot_overall_comparison_grid(all_data)

    println("\n=== Plot 2: Per-subcatchment heatmap (uc=$(uc_target), $(pass_target)) ===")
    plot_subcatchment_heatmaps(all_data, all_sc_data, site_to_sc, uc_target, pass_target)

    println("\nDone. Output saved to: $(abspath(OUTPUT_DIR))")
end

main()
