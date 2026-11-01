# =============================================================================
# network_root_join_diagnostics.jl
#
# Quantitative diagnostics for the reconstructed on-path river network and its
# obstacle overlay (second-review items 4.7 / C4 / E12).
#
# Reports:
#   * parent edges vs root joins: count, length, share, median/max length;
#   * obstacles matched and applied on root-join edges;
#   * compounded passability per link: max obstacle count, min factor, links
#     with factor < 1e-6, links with >= 5 obstacles;
#   * elevation inversions for parent edges vs root joins;
#   * root counts at the three tolerances (0.01, 100 m), (0.02, 200 m),
#     (0.05, 500 m).
#
# Usage:
#   julia --project=. scripts/network_root_join_diagnostics.jl
#
# Env overrides:
#   GUADEX_PARAMETERS_FILE  default parameters_climate_scenarios_corrected.toml
#   GUADEX_NETWORK_OUT      output CSV prefix (default under results/)
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using Guadex
using DataFrames
using CSV
using Statistics

const ROOT = dirname(@__DIR__)
include(joinpath(ROOT, "parameters.jl"))
using .SimulationParameters
get!(ENV, "GUADEX_PARAMETERS_FILE", "parameters_climate_scenarios_corrected.toml")
const PARAMS = SimulationParameters.load()
SimulationParameters.require_sections(PARAMS, "general", "inputs", "obstacles",
    "species", "subcatchments", "run_climate_scenarios")

const RC = PARAMS["run_climate_scenarios"]
const DATA = Guadex.prepare_ode_data(
    upstream_cost = Float64(get(RC, "upstream_cost", 0.05)),
    cedex_var_file = get(PARAMS["inputs"], "cedex_var_file", nothing),
    cedex_esc_uts_file = get(PARAMS["inputs"], "cedex_esc_uts_file", nothing),
    obstacles_file = get(PARAMS["inputs"], "obstacles_file", nothing),
    obstacle_mode = Symbol(get(PARAMS["obstacles"], "mode", "legacy")),
    obstacle_matching_tolerance = Float64(get(PARAMS["obstacles"], "matching_tolerance_m", 2000.0)),
    obstacle_passability = Float64(get(PARAMS["obstacles"], "upstream_passability", 0.1)),
    obstacle_downstream_passability = Float64(get(PARAMS["obstacles"], "downstream_passability", 0.5)),
    connectivity_method = SimulationParameters.connectivity_method(PARAMS),
    carrying_capacity_base_scaling = Float64(get(RC, "carrying_capacity_base_scaling", 1.0))
)

diags = DATA.connectivity_diagnostics
diags === nothing && error("the corrected config did not select the on-path graph")
parent = DATA.connectivity_parent
sites = DATA.sites

is_parent_edge(u, v) = parent[u] == v
parent_lengths = Float64[]
root_lengths = Float64[]
root_pairs = Set{Tuple{String,String}}()
for (u, v, w) in DATA.connectivity_tree_edges
    if is_parent_edge(u, v)
        push!(parent_lengths, w)
    else
        push!(root_lengths, w)
        push!(root_pairs, sites[u] <= sites[v] ? (sites[u], sites[v]) : (sites[v], sites[u]))
    end
end
parent_km = sum(parent_lengths) / 1000
root_km = sum(root_lengths) / 1000
total_km = parent_km + root_km

println("\n=== Parent edges vs root joins ===")
println("parent edges        : $(length(parent_lengths))  $(round(parent_km, digits=1)) km  " *
        "($(round(100 * parent_km / total_km, digits=1)) %)")
println("root-join edges     : $(length(root_lengths))  $(round(root_km, digits=1)) km  " *
        "($(round(100 * root_km / total_km, digits=1)) %)")
println("roots               : $(diags.n_roots) of $(length(sites)) sites")
println("parent edge length  : median $(round(median(parent_lengths) / 1000, digits=1)) km, " *
        "max $(round(maximum(parent_lengths) / 1000, digits=1)) km")
println("root-join length    : median $(round(median(root_lengths) / 1000, digits=1)) km, " *
        "max $(round(maximum(root_lengths) / 1000, digits=1)) km")
println("elevation inversions: parent $(length(diags.parent_elevation_inversions)), " *
        "root $(length(diags.root_join_elevation_inversions))")

md = DATA.obstacle_metadata
mapping = DATA.obstacle_mapping_diagnostics
root_sites = Set(sites[i] for i in eachindex(parent) if parent[i] == 0)
on_root_matched = 0
on_root_applied = 0
incident_matched = 0
incident_applied = 0
if mapping !== nothing
    for r in eachrow(mapping)
        (String(r.edge_from) == "" || String(r.edge_to) == "") && continue
        key = String(r.edge_from) <= String(r.edge_to) ?
            (String(r.edge_from), String(r.edge_to)) :
            (String(r.edge_to), String(r.edge_from))
        is_applied = String(r.outcome) == "applied"
        if key in root_pairs
            global on_root_matched += r.matched ? 1 : 0
            global on_root_applied += is_applied ? 1 : 0
        end
        if String(r.edge_from) in root_sites || String(r.edge_to) in root_sites
            global incident_matched += r.matched ? 1 : 0
            global incident_applied += is_applied ? 1 : 0
        end
    end
end
println("\n=== Obstacles on root-join edges ===")
println("matched on roots    : $on_root_matched of $(md.n_matched) matched")
println("applied on roots    : $on_root_applied of $(md.n_applied) applied")
println("matched incident to a root site: $incident_matched (applied $incident_applied)")

counts = [n for (_, n) in md.link_obstacle_counts]
println("\n=== Compounded passability ===")
if !isempty(counts)
    up_factor(n) = Float64(get(PARAMS["obstacles"], "upstream_passability", 0.1))^n
    factors = up_factor.(counts)
    println("links with obstacles: $(length(counts)) of $(diags.n_edges)")
    println("max obstacles/link  : $(maximum(counts)) (upstream factor $(minimum(factors)))")
    println("links >= 5 obstacles: $(count(x -> x >= 5, counts))")
    println("links factor < 1e-6 : $(count(x -> x < 1e-6, factors))")
    println("min upstream factor : $(minimum(factors)) (link with $(maximum(counts)) obstacles)")
end

# Root count sensitivity to the on-path tolerance.
site_df = DATA.site_df
site_codes = String.(site_df.CODIGO)
d_river = Dict{String,Float64}(string(r.CODIGO) => Float64(coalesce(r."Dist.Guadalq.(m)", 0.0))
                               for r in eachrow(site_df))
elevation = Dict{String,Float64}(string(r.CODIGO) => Float64(coalesce(r.ALTITUD, 500.0))
                                 for r in eachrow(site_df))
println("\n=== Root count vs on-path tolerance ===")
root_sensitivity = NamedTuple[]
const DISTANCE_FILE = joinpath(ROOT, "data", "Matrix_distances_1037puntos_BRUTO_FINAL.csv")
for (rtol, atol) in ((0.01, 100.0), (0.02, 200.0), (0.05, 500.0))
    g = Guadex.build_on_path_distance_matrix(DISTANCE_FILE, site_codes, d_river, elevation;
        rtol=rtol, atol=atol)
    push!(root_sensitivity, (rtol=rtol, atol_m=atol, n_roots=g.diagnostics.n_roots,
        n_parent_edges=g.diagnostics.n_parent_edges, total_length_km=g.diagnostics.total_length_km))
    println("rtol=$(rtol), atol=$(atol) m -> $(g.diagnostics.n_roots) roots, " *
            "$(g.diagnostics.n_parent_edges) parent edges, " *
            "$(round(g.diagnostics.total_length_km, digits=1)) km")
end

OUT_DIR = get(ENV, "GUADEX_NETWORK_OUT", joinpath(ROOT, "results", "network_diagnostics_corrected"))
mkpath(OUT_DIR)
CSV.write(joinpath(OUT_DIR, "root_join_summary.csv"), DataFrame([
    (metric="n_parent_edges", value=length(parent_lengths)),
    (metric="n_root_join_edges", value=length(root_lengths)),
    (metric="n_roots", value=diags.n_roots),
    (metric="n_sites", value=length(sites)),
    (metric="parent_length_km", value=parent_km),
    (metric="root_join_length_km", value=root_km),
    (metric="total_length_km", value=total_km),
    (metric="root_join_share_pct", value=100 * root_km / total_km),
    (metric="root_join_median_km", value=median(root_lengths) / 1000),
    (metric="root_join_max_km", value=maximum(root_lengths) / 1000),
    (metric="obstacles_matched_on_roots", value=on_root_matched),
    (metric="obstacles_applied_on_roots", value=on_root_applied),
    (metric="obstacles_matched_incident_to_root_site", value=incident_matched),
    (metric="obstacles_applied_incident_to_root_site", value=incident_applied),
    (metric="parent_elevation_inversions", value=length(diags.parent_elevation_inversions)),
    (metric="root_elevation_inversions", value=length(diags.root_join_elevation_inversions)),
    (metric="links_with_obstacles", value=length(counts)),
    (metric="max_obstacles_per_link", value=isempty(counts) ? 0 : maximum(counts)),
    (metric="min_link_passability", value=isempty(counts) ? 1.0 :
        Float64(get(PARAMS["obstacles"], "upstream_passability", 0.1))^maximum(counts)),
    (metric="links_ge5_obstacles", value=count(x -> x >= 5, counts)),
    (metric="links_factor_below_1e-6", value=isempty(counts) ? 0 :
        count(x -> x < 1e-6,
            Float64(get(PARAMS["obstacles"], "upstream_passability", 0.1)).^counts)),
]))
CSV.write(joinpath(OUT_DIR, "root_count_tolerance.csv"), DataFrame(root_sensitivity))
println("\nWrote $OUT_DIR")
