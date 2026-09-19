using Pkg; Pkg.activate(".")
using DataFrames
using CSV
using CairoMakie
using Statistics

# =============================================================================
# Figures for docs/integrated_climate_obstacle_report.md
#
# Reads ONLY already-completed exported tables (no simulation, results/ is
# read-only) and writes new report figures to docs/figures/.
#
# Produces:
#   docs/figures/fig01_thermal_niches.png        corrected Fig. 1 (legend outside)
#   docs/figures/fig12_st_climate_response.png   ST climate + exposure response
#   docs/figures/fig13_interaction_matrix_effects.png  matrix main-effect levels
# =============================================================================

const ROOT = normpath(joinpath(@__DIR__, ".."))
const OUT = joinpath(ROOT, "docs", "figures")
mkpath(OUT)

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

# ---------------------------------------------------------------------------
# Figure 1 (corrected): native thermal niches vs site temperature.
# The only change from the original diagnostic is that the legend is placed in
# its own figure column so it cannot overlap the curves or the histogram.
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

    fig = Figure(size=(1450, 720))
    ax = Axis(fig[1, 1];
        title = "Native thermal niches vs site temperature",
        xlabel = "Water temperature (°C)",
        ylabel = "Thermal suitability / relative frequency")

    lo, hi = 8.0, 26.0
    edges = collect(lo:1.0:hi)
    counts = zeros(Int, length(edges) - 1)
    for t in temps
        (lo <= t < hi) || continue
        counts[clamp(floor(Int, t - lo) + 1, 1, length(counts))] += 1
    end
    max_count = maximum(counts; init = 0)
    max_count == 0 || barplot!(ax, edges[1:(end - 1)] .+ 0.5, counts ./ max_count;
        width = 0.95, color = (:gray, 0.18), label = "site temperatures (rel. freq.)")

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
    labels = String["site temperatures (rel. freq.)"]
    for (i, (key, codes)) in enumerate(sort(collect(groups); by = first))
        opt, sigma = key
        lines!(ax, t_grid, thermal_curve(opt, sigma, t_grid);
            color = palette[mod1(i, length(palette))], linewidth = 2.2)
        push!(handles, LineElement(color = palette[mod1(i, length(palette))], linewidth = 2.2))
        push!(labels, length(codes) == 1 ? "$(only(codes)): opt $(round(opt, digits = 1))" :
            "opt $(round(opt, digits = 1)) ×$(length(codes)) ($(join(codes, ",")))")
    end

    vlines!(ax, [baseline]; color = :black, linewidth = 1.6, linestyle = :dot)
    vlines!(ax, [baseline + warming_delta]; color = :crimson, linewidth = 1.6, linestyle = :dot)
    text!(ax, baseline - 0.15, 1.06; text = "baseline $(round(baseline, digits = 1))°C",
        align = (:right, :bottom), fontsize = 10)
    text!(ax, baseline + warming_delta + 0.15, 1.06;
        text = "+$(round(warming_delta, digits = 2))°C → $(round(baseline + warming_delta, digits = 1))°C",
        align = (:left, :bottom), fontsize = 10, color = :crimson)
    ylims!(ax, 0.0, 1.12)

    # Legend in its own column, fully outside the plotting area.
    Legend(fig[1, 2], handles, labels; framevisible = false, labelsize = 11,
        title = "Native species (optimum °C)", titlesize = 11, nbanks = 1)
    Label(fig[0, 1:2],
        "Thermal niches explain the warming response: most native optima lie above current water temperature",
        fontsize = 14, font = :bold)

    save(joinpath(OUT, "fig01_thermal_niches.png"), fig; px_per_unit = 2)
    println("wrote fig01_thermal_niches.png")
end

# ---------------------------------------------------------------------------
# Figure 12: ST climate response and thermal exposure across the 44-run ensemble
# ---------------------------------------------------------------------------
function fig12_st_climate_response()
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

    fig = Figure(size=(2000, 640))
    Label(fig[0, 1:3],
        "Cold-water sentinel: ST response and thermal exposure across 44 climate-model runs " *
        "(2026–2045, control shown black)", fontsize = 16, font = :bold)

    ax1 = Axis(fig[1, 1]; xlabel = "Warming by 2045 (°C, basin mean)",
        ylabel = "ST relative biomass change (%)",
        title = "ST basin biomass change")
    hlines!(ax1, [0.0]; color = (:black, 0.4), linewidth = 1)
    for sc in unique(scenario_rows.scenario)
        s = scenario_rows[scenario_rows.scenario .== sc, :]
        scatter!(ax1, s.warming, 100 .* s.rel; color = SSP_COLORS[sc], markersize = 11,
            label = sc, strokewidth = 0)
    end
    xx = range(minimum(warm), maximum(warm); length = 50)
    lines!(ax1, xx, beta[1] .+ beta[2] .* xx; color = :black, linewidth = 1.6, linestyle = :dash)
    control === nothing || scatter!(ax1, [0.0], [100 * control.rel]; color = :black,
        marker = :diamond, markersize = 15, label = "control")
    text!(ax1, 0.03, 0.06; space = :relative,
        text = "slope = $(round(beta[2], digits = 2)) %/°C\nR² = $(round(r2, digits = 2))",
        align = (:left, :bottom), fontsize = 11)
    axislegend(ax1; position = :rt, framevisible = false, labelsize = 9)

    ax2 = Axis(fig[1, 2]; xlabel = "Warming by 2045 (°C, basin mean)",
        ylabel = "ST final biomass (model units)", title = "ST final biomass")
    for sc in unique(scenario_rows.scenario)
        s = scenario_rows[scenario_rows.scenario .== sc, :]
        scatter!(ax2, s.warming, s.final; color = SSP_COLORS[sc], markersize = 11, strokewidth = 0)
    end
    control === nothing || scatter!(ax2, [0.0], [control.final]; color = :black,
        marker = :diamond, markersize = 15)

    ax3 = Axis(fig[1, 3]; xlabel = "Warming by 2045 (°C, basin mean)",
        ylabel = "Max days/yr above 20 °C", title = "ST thermal exposure")
    for sc in unique(scenario_rows.scenario)
        s = scenario_rows[scenario_rows.scenario .== sc, :]
        scatter!(ax3, s.warming, s.max_days; color = SSP_COLORS[sc], markersize = 11, strokewidth = 0)
    end
    control === nothing || begin
        scatter!(ax3, [0.0], [control.max_days]; color = :black, marker = :diamond, markersize = 15)
        hlines!(ax3, [control.max_days]; color = (:black, 0.5), linewidth = 1.2, linestyle = :dot)
        text!(ax3, 0.03, 0.06; space = :relative,
            text = "control = $(control.max_days) days/yr", align = (:left, :bottom), fontsize = 11)
    end

    save(joinpath(OUT, "fig12_st_climate_response.png"), fig; px_per_unit = 2)
    println("wrote fig12_st_climate_response.png (n=$(n) scenario runs)")
end

# ---------------------------------------------------------------------------
# Figure 13: interaction-matrix main effects (the -90 % / -45 % headline)
# ---------------------------------------------------------------------------
function fig13_matrix_effects()
    alt = CSV.read(ALT_INDEX, DataFrame)
    matrices = ["original", "random", "invasive_favoring"]
    mcolor = Dict("original" => RGBAf(0.20, 0.45, 0.75, 0.9),
        "random" => RGBAf(0.55, 0.55, 0.55, 0.9),
        "invasive_favoring" => RGBAf(0.75, 0.30, 0.10, 0.9))
    metrics = ["basin_native_biomass", "basin_invasive_biomass", "st_final_biomass", "basin_native_richness"]
    titles = ["Basin native biomass", "Basin invasive biomass", "ST final biomass", "Basin native richness"]

    fig = Figure(size=(1750, 980))
    Label(fig[0, 1:2],
        "Interaction-matrix main effects (σ = 0.3 sweep; 16 runs per matrix, both climate ends and costs)",
        fontsize = 16, font = :bold)

    for (i, (metric, title)) in enumerate(zip(metrics, titles))
        row, col = divrem(i - 1, 2)
        ax = Axis(fig[row + 1, col + 1]; title = title,
            xticks = (1:3, ["original", "random", "invasive\nfavoring"]))
        for (mi, m) in enumerate(matrices)
            sub = Float64.(alt[alt.interaction_matrix .== m, Symbol(metric)])
            jitter = 0.10 .* sin.(1:length(sub))
            scatter!(ax, fill(mi, length(sub)) .+ jitter, sub; color = mcolor[m],
                markersize = 9, strokewidth = 0)
            mu = mean(sub)
            lines!(ax, [mi - 0.22, mi + 0.22], [mu, mu]; color = :black, linewidth = 2.6)
        end
        if i == 1
            orig = mean(Float64.(alt[alt.interaction_matrix .== "original", Symbol(metric)]))
            inv = mean(Float64.(alt[alt.interaction_matrix .== "invasive_favoring", Symbol(metric)]))
            text!(ax, 0.02, 0.30; space = :relative,
                text = "original → invasive-favoring:\n$(round(orig, digits = 1)) → $(round(inv, digits = 1)) " *
                       "($(round(100 * (inv - orig) / orig, digits = 0))%)",
                align = (:left, :top), fontsize = 11)
        end
        if i == 3
            orig = mean(Float64.(alt[alt.interaction_matrix .== "original", Symbol(metric)]))
            inv = mean(Float64.(alt[alt.interaction_matrix .== "invasive_favoring", Symbol(metric)]))
            text!(ax, 0.02, 0.16; space = :relative,
                text = "original → invasive-favoring:\n$(round(orig, digits = 0)) → $(round(inv, digits = 0)) " *
                       "($(round(100 * (inv - orig) / orig, digits = 0))%)",
                align = (:left, :top), fontsize = 11)
        end
    end

    save(joinpath(OUT, "fig13_interaction_matrix_effects.png"), fig; px_per_unit = 2)
    println("wrote fig13_interaction_matrix_effects.png")
end

fig01_thermal_niches()
fig12_st_climate_response()
fig13_matrix_effects()
println("done; figures in $OUT")
