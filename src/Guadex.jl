module Guadex

using Statistics
using DataFrames
using CSV
using Graphs
using GraphMakie
using CairoMakie
using NetworkLayout
using Colors
using DifferentialEquations
using LinearAlgebra
using SparseArrays
using Makie
using GeoMakie

# Export main functions
export build_stream_graph, get_graph_statistics, find_upstream_sites, find_downstream_sites, visualize_stream_graph, analyze_stream_connectivity, build_nearest_neighbor_connections!, build_threshold_connections!, build_mst_connections!, analyze_stream_connectivity
export plot_catchment_network
export metacommunity_ode!, MetacommunityParams, precompute_dispersal_matrix

# Export time-varying temperature schedule model
export TemperatureSchedule, ScheduledMetacommunityParams
export metacommunity_ode_scheduled!, temperature_delta

# Export daily/seasonal forcing, optimum sweep and heat-stress calibration (WP1-WP3)
export load_daily_temperature_forcing, daily_forcing_matrix
export load_daily_forcing_any, is_wide_daily_forcing, wide_forcing_matrix, wide_baseline_means
export daily_temperature_schedule, annual_mean_schedule, annual_mean_deltas
export annual_mean_deltas_by_year, baseline_climatology_schedule
export optimum_sweep_optima, exceedance_energy, calibrate_heat_stress_rate
export exposure_days, exposure_table

# Export equilibrium spin-up and K sensitivity (WP4)
export spin_up, scale_carrying_capacity, set_thermal_optima, set_heat_stress_rate
export set_dispersal_matrix, set_interaction_matrix, set_thermal_sigma_multiplier
export with_temperature_baseline, site_totals

# Export four-level reporting / viewer output functions
export load_site_level_crosswalk, site_level_vectors, species_indices
export compute_site_metrics, compute_species_metrics, quasi_extinction_flags
export quasi_extinction_summary, time_to_quasi_extinction
export aggregate_metrics, site_connectivity_metrics
export export_run_outputs, write_viewer_outputs
export load_temperature_projections, basin_warming_curve, warming_matrix

# Export climate-scenario figure functions
export plot_climate_scenario_figures, plot_climate_run_figure
export plot_ensemble_level_figure, plot_across_scenarios_figure, plot_final_year_summary
export discover_climate_runs, load_climate_run, run_level_stats, ensemble_level_stats
export climate_metric_column

# Export climate diagnostic figures and their data readers
export ClimateBasinSeries, read_climate_basin_series
export plot_climate_forcing_response, plot_climate_thermal_niches
export plot_climate_community_filling, plot_climate_diagnostics

# Export data preparation functions
export prepare_ode_data, save_ode_data
export load_species_characteristics, load_site_data, load_species_density_data
export load_interaction_matrix, build_distance_matrix
export load_cedex_var, load_cedex_esc_uts, select_cedex_uts, load_obstacles
export build_elevation_vector, build_dam_passability_matrix
export build_site_coordinate_matrix, build_obstacle_passability_matrix
export extract_site_temperatures, extract_habitat_suitability
export build_intrinsic_growth_rates

# Export visualization functions
export plot_ode_solution, plot_total_biomass, plot_species_richness
export plot_avg_total_biomass, plot_avg_species_richness
export plot_sites_map, plot_site_connectivity_map, plot_subcatchment_network
export plot_combined_analysis, save_figure

# Export sensitivity visualization functions
export plot_richness_change_per_site, plot_richness_change_per_subcatchment
export plot_richness_timeseries_grid

include("graph_construction.jl")
include("visualize_graph.jl")
include("ode_model.jl")
include("temperature_forcing.jl")
include("spin_up.jl")
include("data_preparation.jl")
include("outputs.jl")
include("visualization.jl")
include("climate_figures.jl")
include("climate_diagnostics.jl")

end # module Guadex
