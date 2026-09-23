using Pkg; Pkg.activate(".")
using DataFrames
using CSV
using CairoMakie
using Statistics

# =============================================================================
# Figures for the integrated climate / obstacle reports.
#
# Reads ONLY already-completed exported tables (no simulation) and writes the
# report figures next to the results they visualise:
#   results/climate_scenarios_k1x_burnin/figures/report_plots/fig01_thermal_niches.png
#   results/climate_scenarios_k1x_burnin/figures/report_plots/fig12_st_climate_response.png
#   results/sensitivity_obstacles/alt_interactions/report_plots/fig13_interaction_matrix_effects.png
#
# All figures use the shared report style (fonts >= 18 pt, no descriptive
# titles, panels lettered a, b, c ...); the captions in the LaTeX reports
# explain every visual element.
# =============================================================================

# Reuse the report theme / panel-letter helpers defined for the Guadex package.
include(joinpath(@__DIR__, "..", "src", "report_style.jl"))

const ROOT = normpath(joinpath(@__DIR__, ".."))
# Report figures live with the results they visualise, not in docs/.
const OUT_CLIMATE = joinpath(ROOT, "results", "climate_scenarios_k1x_burnin", "figures", "report_plots")
const OUT_ALT = joinpath(ROOT, "results", "sensitivity_obstacles", "alt_interactions", "report_plots")
mkpath(OUT_CLIMATE); mkpath(OUT_ALT)

const CLIMATE_ROOT = joinpath(ROOT, "results", "climate_scenarios_k1x_burnin")
const ALT_INDEX = joinpath(ROOT, "results", "sensitivity_obstacles", "alt_interactions", "runs_index.csv")
const SPECIES_CHARS = joinpath(ROOT, "data", "ABIOTIC", "caracteristicas_peces_Guadalquivir_03-04-2018.csv")
const SITE_TABLE = joinpath(CLIMATE_ROOT, "control", "baseline", "export", "levels", "level_sampling_point.csv")

const NATIVE_CODES = ["AB", "AH", "SP", "PW", "LS", "SA", "IL", "CP", "IO", "ST"]
const SSP_COLORS = Dict("ssp126" => RGBAf(0.11, 0.36, 0.60, 0.9),
    "ssp245" => RGBAf(0.20, 0.60, 0.35, 0.9),
    "ssp370" => RGBAf(0.90, 0.55, 0.10, 0.9),
    "ssp585" => RGBAf(0.70, 0.15, 0.15, 0.9),
    "control" => RGBAf(0.0, 0.0, 0.0, 0.9))
const SSP_ORDER = ["ssp126", "ssp245", "ssp370", "ssp585", "control"]

run_path(row) = normpath(joinpath(ROOT, replace(String(row.run_dir), "\\" => "/")))

function read_st_metrics(run_dir)
    path = joinpath(run_dir, "export", "levels", "quasi_extinction_summary.csv")
    isfile(path) || return nothing
    df = CSV.read(path, DataFrame)
    sub = df[df.species .== "ST", :]
    isempty(sub) && return nothing
    return (relative = Float64(sub.relative_biomass_change[1]),
        final = Float64(sub.final_biomass[1]),
        baseline = Float64(sub.baseline_biomass[1]))
end

function read_st_exposure(run_dir)
    path = joinpath(run_dir, "export", "levels", "exposure_sites.csv")
    isfile(path) || return nothing
    df = CSV.read(path, DataFrame)
    sub = df[df.species .== "ST", :]
    isempty(sub) && return nothing
    return (max_days = maximum(Int.(sub.exposure_days)),
        max_energy = maximum(Float64.(sub.exceedance_energy)))
end

# Axes here use explicit width/height, so the requested Figure size can be
# smaller than the space the panels plus their labels, legends and titles
# occupy; Makie then clips the overflow at the canvas edge.  `resize_to_layout!`
# grows the scene to the tight bounding box of the layout, including every
# protrusion, so no content is cut off while font sizes and panel dimensions
# stay as authored.  `size` is kept for call-site compatibility but no longer
# overrides the fitted size.
function save_report(fig, path; size=nothing)
    resize_to_layout!(fig)
    Makie.save(path, fig; px_per_unit=2)
    trim_figure!(path)
    println("wrote $(path)")
end

# ---------------------------------------------------------------------------
# Figure 1: native thermal niches vs site temperature.
# ---------------------------------------------------------------------------
function parse_range(s)
    m = match(r"([0-9]+(?:\.[0-9]+)?)\s*to\s*([0-9]+(?:\.[0-9]+)?)", lowercase(strip(String(s))))
    m === nothing && return nothing
    return (parse(Float64, m.captures[1]), parse(Float64, m.captures[2]))
end

function thermal_curve(opt, sigma, t)
    return exp.(-(t .- opt) .^ 2 ./ (2 * sigma^2))
end

function fig01_thermal_niches()
    report_theme!()
    chars = CSV.read(SPECIES_CHARS, DataFrame; delim=';')
    lookup = Dict{String,Tuple{Float64,Float64}}()
    for r in eachrow(chars)
        rng = parse_range(r.TEMPERATURE_C)
        rng === nothing && continue
        lo, hi = rng
        lookup[lowercase(String(r.SP))] = (0.5 * (lo + hi), (hi - lo) / 6)
    end

    site = CSV.read(SITE_TABLE, DataFrame)
    temps = Float64.(site.temperature_c[site.year .== 2026])
    baseline = mean(temps)

    idx = CSV.read(joinpath(CLIMATE_ROOT, "runs_index.csv"), DataFrame)
    warming_delta = median(Float64.(idx.warming_end_degc[idx.scenario .!= "control"]))

    fig = Figure(size=(1700, 1050))
    ax = Axis(fig[1, 1]; width=760, height=650,
        xlabel = "Water temperature (°C)",
        ylabel = "Thermal suitability / rel. frequency")

    lo, hi = 8.0, 26.0
    edges = collect(lo:1.0:hi)
    counts = zeros(Int, length(edges) - 1)
    for t in temps
        (lo <= t < hi) || continue
        counts[clamp(floor(Int, t - lo) + 1, 1, length(counts))] += 1
    end
    max_count = maximum(counts; init = 0)
    max_count == 0 || barplot!(ax, edges[1:(end - 1)] .+ 0.5, counts ./ max_count;
        width = 0.95, color = (:gray, 0.18))

    groups = Dict{Tuple{Float64,Float64},Vector{String}}()
    for code in NATIVE_CODES
        haskey(lookup, lowercase(code)) || continue
        opt, sigma = lookup[lowercase(code)]
        push!(get!(groups, (round(opt, digits = 4), round(sigma, digits = 4)), String[]), code)
    end

    t_grid = range(lo - 0.5, hi + 1.0; length = 300)
    palette = [:dodgerblue, :seagreen, :darkorange, :firebrick, :mediumpurple,
        :teal, :goldenrod, :slateblue, :olive]
    handles = Any[PolyElement(color = (:gray, 0.35))]
    labels = String["Site-temperature distribution (relative frequency)"]
    for (i, (key, codes)) in enumerate(sort(collect(groups); by = first))
        opt, sigma = key
        lines!(ax, t_grid, thermal_curve(opt, sigma, t_grid);
            color = palette[mod1(i, length(palette))], linewidth = 2.6)
        push!(handles, LineElement(color = palette[mod1(i, length(palette))], linewidth = 2.6))
        push!(labels, length(codes) == 1 ? "$(only(codes)): optimum $(round(opt, digits = 1)) °C" :
            "optimum $(round(opt, digits = 1)) °C ×$(length(codes)) ($(join(codes, ",")))")
    end

    vlines!(ax, [baseline]; color = :black, linewidth = 2.0, linestyle = :dot)
    vlines!(ax, [baseline + warming_delta]; color = :crimson, linewidth = 2.0, linestyle = :dot)
    text!(ax, baseline - 0.15, 1.06; text = "Baseline mean $(round(baseline, digits = 1)) °C",
        align = (:right, :bottom), fontsize = REPORT_ANNOTATION_FONTSIZE)
    text!(ax, baseline + warming_delta + 0.15, 1.06;
        text = "+$(round(warming_delta, digits = 2)) °C → $(round(baseline + warming_delta, digits = 1)) °C",
        align = (:left, :bottom), fontsize = REPORT_ANNOTATION_FONTSIZE, color = :crimson)
    ylims!(ax, 0.0, 1.12)

    Legend(fig[1, 2], handles, labels; framevisible = false,
        title = "Native species (optimum)", titlefontsize = REPORT_LEGEND_FONTSIZE,
        nbanks = 1)
    equal_panel_columns!(fig, 1, 0.55)
    save_report(fig, joinpath(OUT_CLIMATE, "fig01_thermal_niches.png"); size=(1700, 1050))
end

# ---------------------------------------------------------------------------
# Figure 12: ST climate response and thermal exposure across the 44-run ensemble
# ---------------------------------------------------------------------------
function fig12_st_climate_response()
    report_theme!()
    idx = CSV.read(joinpath(CLIMATE_ROOT, "runs_index.csv"), DataFrame)

    scenario_rows = DataFrame(scenario = String[], gcm = String[], warming = Float64[],
        rel = Float64[], final = Float64[], max_days = Int[], max_energy = Float64[])
    control = nothing
    for r in eachrow(idx)
        dir = run_path(r)
        st = read_st_metrics(dir)
        ex = read_st_exposure(dir)
        st === nothing && continue
        ex === nothing && continue
        if String(r.scenario) == "control"
            control = (warming = 0.0, rel = st.relative, final = st.final,
                max_days = ex.max_days, max_energy = ex.max_energy)
            continue
        end
        push!(scenario_rows, (String(r.scenario), String(r.gcm), Float64(r.warming_end_degc),
            st.relative, st.final, ex.max_days, ex.max_energy))
    end

    warm = scenario_rows.warming
    rel_pct = 100 .* scenario_rows.rel
    n = length(warm)
    X = hcat(ones(n), warm)
    beta = X \ rel_pct
    yhat = X * beta
    r2 = 1 - sum(abs2, rel_pct .- yhat) / sum(abs2, rel_pct .- mean(rel_pct))

    fig = Figure(size=(1700, 1000))
    ax1 = Axis(fig[1, 1]; width=420, height=540,
        xlabel = "Warming by 2045 (°C)",
        ylabel = "Brown trout change (%)")
    panel_letter!(ax1, "(a)")
    hlines!(ax1, [0.0]; color = (:black, 0.4), linewidth = 1)
    for sc in unique(scenario_rows.scenario)
        s = scenario_rows[scenario_rows.scenario .== sc, :]
        scatter!(ax1, s.warming, 100 .* s.rel; color = SSP_COLORS[sc], markersize = 14,
            strokewidth = 0)
    end
    xx = range(minimum(warm), maximum(warm); length = 50)
    lines!(ax1, xx, beta[1] .+ beta[2] .* xx; color = :black, linewidth = 2.0, linestyle = :dash)
    control === nothing || scatter!(ax1, [0.0], [100 * control.rel]; color = :black,
        marker = :diamond, markersize = 18)
    text!(ax1, 0.03, 0.06; space = :relative,
        text = "OLS slope = $(round(beta[2], digits = 2)) % per °C\nR² = $(round(r2, digits = 2))",
        align = (:left, :bottom), fontsize = REPORT_ANNOTATION_FONTSIZE)

    ax2 = Axis(fig[1, 2]; width=420, height=540,
        xlabel = "Warming by 2045 (°C)",
        ylabel = "Brown trout biomass (units)")
    panel_letter!(ax2, "(b)")
    for sc in unique(scenario_rows.scenario)
        s = scenario_rows[scenario_rows.scenario .== sc, :]
        scatter!(ax2, s.warming, s.final; color = SSP_COLORS[sc], markersize = 14, strokewidth = 0)
    end
    control === nothing || scatter!(ax2, [0.0], [control.final]; color = :black,
        marker = :diamond, markersize = 18)

    ax3 = Axis(fig[1, 3]; width=420, height=540,
        xlabel = "Warming by 2045 (°C)",
        ylabel = "Days per year above 20 °C")
    panel_letter!(ax3, "(c)")
    for sc in unique(scenario_rows.scenario)
        s = scenario_rows[scenario_rows.scenario .== sc, :]
        scatter!(ax3, s.warming, s.max_days; color = SSP_COLORS[sc], markersize = 14, strokewidth = 0)
    end
    control === nothing || begin
        scatter!(ax3, [0.0], [control.max_days]; color = :black, marker = :diamond, markersize = 18)
        hlines!(ax3, [control.max_days]; color = (:black, 0.5), linewidth = 1.4, linestyle = :dot)
        text!(ax3, 0.03, 0.06; space = :relative,
            text = "No-warming control = $(control.max_days) days per year",
            align = (:left, :bottom), fontsize = REPORT_ANNOTATION_FONTSIZE)
    end

    handles = Any[]
    labels = String[]
    for sc in SSP_ORDER
        haskey(SSP_COLORS, sc) || continue
        is_control = sc == "control"
        push!(handles, is_control ?
            MarkerElement(color = :black, marker = :diamond, markersize = 18) :
            MarkerElement(color = SSP_COLORS[sc], marker = :circle, markersize = 15))
        push!(labels, is_control ? "No-warming control (0 °C)" : sc)
    end
    equal_panel_columns!(fig, 1, 1, 1)
    Legend(fig[2, 1:3], handles, labels; orientation = :horizontal,
        framevisible = false, title = "Point = one GCM × SSP run",
        titlefontsize = REPORT_LEGEND_FONTSIZE)
    save_report(fig, joinpath(OUT_CLIMATE, "fig12_st_climate_response.png"); size=(1700, 1000))
    println("fig12 uses n=$(n) scenario runs")
end

# ---------------------------------------------------------------------------
# Figure 13: interaction-matrix main effects (the -90 % headline)
# ---------------------------------------------------------------------------
function fig13_matrix_effects()
    report_theme!()
    alt = CSV.read(ALT_INDEX, DataFrame)
    matrices = ["original", "random", "invasive_favoring"]
    matrix_labels = ["Original", "Random\n(placebo)", "Invasive-\nfavouring"]
    mcolor = Dict("original" => RGBAf(0.20, 0.45, 0.75, 0.9),
        "random" => RGBAf(0.55, 0.55, 0.55, 0.9),
        "invasive_favoring" => RGBAf(0.75, 0.30, 0.10, 0.9))
    metrics = ["basin_native_biomass", "basin_invasive_biomass", "st_final_biomass", "basin_native_richness"]
    ylabels = ["Native biomass (units)", "Invasive biomass (units)",
        "Brown trout biomass (units)", "Native richness (sp.)"]

    fig = Figure(size=(1600, 1500))
    for (i, (metric, ylab)) in enumerate(zip(metrics, ylabels))
        row, col = divrem(i - 1, 2)
        ax = Axis(fig[row + 1, col + 1]; width=520, height=520,
            ylabel = ylab,
            xticks = (1:3, matrix_labels))
        panel_letter!(ax, "(" * string(Char(96 + i)) * ")")
        for (mi, m) in enumerate(matrices)
            sub = Float64.(alt[alt.interaction_matrix .== m, Symbol(metric)])
            jitter = 0.10 .* sin.(1:length(sub))
            scatter!(ax, fill(mi, length(sub)) .+ jitter, sub; color = mcolor[m],
                markersize = 12, strokewidth = 0)
            mu = mean(sub)
            lines!(ax, [mi - 0.22, mi + 0.22], [mu, mu]; color = :black, linewidth = 3.0)
        end
        orig = mean(Float64.(alt[alt.interaction_matrix .== "original", Symbol(metric)]))
        inv = mean(Float64.(alt[alt.interaction_matrix .== "invasive_favoring", Symbol(metric)]))
        report_annotation!(ax,
            "Original → invasive-favouring:\n$(round(orig, digits = 1)) → $(round(inv, digits = 1)) " *
            "($(round(100 * (inv - orig) / orig, digits = 0))%)", position = :rt)
    end

    equal_panel_columns!(fig, 1, 1)
    save_report(fig, joinpath(OUT_ALT, "fig13_interaction_matrix_effects.png"); size=(1600, 1500))
end

fig01_thermal_niches()
fig12_st_climate_response()
fig13_matrix_effects()
println("done; figures in $OUT_CLIMATE and $OUT_ALT")
