# Climate-driven simulations and four-level outputs

This note describes the climate-scenario runner, the four-level reporting
outputs requested by the GuadeX team (`docs/Reporting/README.NewSimus_Sept26.md`),
and the files consumed by the `viz/` 3-D explorer (`viz/README.md`).

## Reporting levels

Every simulation output can be reported at four nested levels:

| level | id | source |
| :--- | :--- | :--- |
| sampling point | `CODIGO` (e.g. `1.1.2`) | model site table |
| sub-basin | `CODIGO_S` (e.g. `1.1`) | model site table (used by the scenarios) |
| water body | `ID_masa` / `EUMASCod` (e.g. `ES050MSPF011002001`) | `data/site_waterbody_crosswalk.csv` |
| whole basin | `ES050` | all modelled sites |

The crosswalk is generated from the GIS archive (no geometry is needed, only the
attribute tables):

```bash
cd viz
npm run crosswalk     # writes ../data/site_waterbody_crosswalk.csv (1037 rows)
```

For each level the pipeline reports the mean of the per-site metrics
(richness, biomass, relative native-richness loss / extinction risk,
temperature, habitat suitability, connectivity), plus `n_sites`. See
`src/outputs.jl` for the exact column set.

## Temperature scenarios

`guadex_tw/outputs/tables/water_temp_future_2045.csv` holds per-GCM x SSP
projections (11 GCMs x 4 SSPs) of `delta_tw_mean` for three windows
(2021-2040, 2036-2055, 2041-2070) at 6 sites. The runner turns these into a
year-by-year basin warming curve:

1. For each window, take the across-site median `delta_tw_mean` of the selected
   GCM/SSP.
2. Anchor zero warming at `start_year` and the window medians at their midpoints
   (2030.5, 2045.5, 2055.5).
3. Linearly interpolate the curve at each simulation year from `start_year` to
   `end_year` (default 2026-2045).

The curve is then expanded to every site. By default the same basin anomaly is
applied to all sites. With `elevation_scaling = true` the anomaly is scaled by
the calibrated water-to-air response `(b1 + b3·z)` using the coefficients
reported by the guadex_tw pipeline (`b1 ≈ 0.51`, `b3 ≈ -0.163`), so high reaches
warm slightly less. This is a model-based refinement, not an independent spatial
projection, and it is therefore off by default.

## Running the scenarios

```bash
julia --project=. run_climate_scenarios.jl
```

Settings live in `[run_climate_scenarios]` in `parameters.toml`:

```toml
[run_climate_scenarios]
start_year = 2026
end_year = 2045
save_interval_days = 30.4167   # 365/12, so monthly saves land on year boundaries
elevation_scaling = false
presence_threshold = 0.1
upstream_cost = 0.05
make_figures = true            # generate the four-level figures after the run
temperature_projections_file = "guadex_tw/outputs/tables/water_temp_future_2045.csv"
scenarios = ["ssp126", "ssp245", "ssp370", "ssp585"]
gcms = ["ACCESS-CM2", "CMCC-CM2-SR5", ...]
```

The time step is **daily** (1 model time unit = 1 day); states are saved
**monthly** to keep outputs small while retaining annual and seasonal metrics.
Each GCM x SSP run is written to
`results/climate_scenarios/<ssp>/<gcm>/` and a `runs_index.csv` summarises the
runs.

Environment overrides (useful for smoke tests):

| variable | meaning |
| :--- | :--- |
| `GUADEX_CLIMATE_START_YEAR`, `GUADEX_CLIMATE_END_YEAR` | horizon |
| `GUADEX_CLIMATE_GCMS` | comma-separated GCM subset |
| `GUADEX_CLIMATE_SAVE_DAYS` | save interval (days) |
| `GUADEX_CLIMATE_ELEVATION_SCALING` | `1`/`true` enables elevation scaling |
| `GUADEX_CLIMATE_MAX_RUNS` | cap the number of runs (testing aid) |
| `GUADEX_CLIMATE_FORCE` | `1`/`true` recomputes runs whose `simulation_output.jld2` already exists |
| `GUADEX_CLIMATE_PLOT` | `0`/`false` skips the post-run figure generation (`make_figures`) |
| `GUADEX_CLIMATE_PROJECTIONS_FILE` | alternative projection table |

Completed runs are skipped by default (set `GUADEX_CLIMATE_FORCE=1` to recompute),
and `runs_index.csv` is rewritten after each run so an interrupted ensemble can
be resumed without losing the summary of finished runs.

A quick smoke run:

```bash
GUADEX_CLIMATE_MAX_RUNS=1 GUADEX_CLIMATE_END_YEAR=2027 GUADEX_CLIMATE_GCMS=ACCESS-CM2 \
  julia --project=. run_climate_scenarios.jl
```

## Output layout

For one run directory `<run_dir>`:

```
<run_dir>/simulation_output.jld2
<run_dir>/export/run_metadata.json
<run_dir>/export/levels/level_sampling_point.csv
<run_dir>/export/levels/level_subcatchment.csv
<run_dir>/export/levels/level_water_body.csv
<run_dir>/export/levels/level_basin.csv
<run_dir>/export/viewer/guadex_results_timeseries.csv
<run_dir>/export/viewer/guadex_results_metrics.json
<run_dir>/export/viewer/guadex_results_native_extinction_risk.json
<run_dir>/export/viewer/level_subcatchment_mean_<metric>_timeseries.json
<run_dir>/export/viewer/level_water_body_mean_<metric>_timeseries.json
<run_dir>/export/viewer/level_basin_mean_<metric>_timeseries.json
```

### Viewer formats

`guadex_results_timeseries.csv` is keyed by `CODIGO` with a `step` column
(year); the explorer's CSV parser turns each metric into a time series
(`"native_richness @ 2027"`, ...). `guadex_results_<metric>.json` matches the
single-variable time-series format in `viz/README.md`:

```json
{ "name": "...", "unit": "fraction", "steps": ["2026", "2027"],
  "data": { "1.1.2": [0.0, 0.12] } }
```

`guadex_results_metrics.json` matches the per-site metrics format:

```json
{ "name": "...", "data": { "1.1.2": { "native_richness": 3, "total_biomass": 63.1 } } }
```

The per-level time-series JSON files use the same shape but are keyed by
subcatchment, water body or `ES050`; the current viewer renders only site-level
data, so these are ready for a future level selector.

## Climate figures

`src/climate_figures.jl` turns the four-level CSVs into report figures. It reads
only `export/levels/*.csv` (plus `export/run_metadata.json` for discovery) and
never re-runs the model, so figures can be regenerated at any time:

```bash
julia --project=. scripts/plot_climate_scenarios.jl [results_root] [figures_dir]
```

Defaults: `results_root = results/climate_scenarios` and
`figures_dir = <results_root>/figures`. `results_root` must contain either a
`runs_index.csv` (preferred) or a `<scenario>/<gcm>/export/levels/level_basin.csv`
tree. `GUADEX_CLIMATE_END_YEAR` overrides the expected final year.

The same plotting runs automatically at the end of `run_climate_scenarios.jl`
when `make_figures = true` (the default); disable it with
`make_figures = false` or `GUADEX_CLIMATE_PLOT=0`.

### Figure inventory

| path | content |
| :--- | :--- |
| `per_run/<scenario>__<gcm>.png` | one run: 5 metrics (rows) x 4 levels (columns); nested levels show the across-unit median with a 10-90% band (site spread at the sampling point) |
| `ensemble/ensemble_<scenario>_<level>.png` | one scenario + level: 5 metric panels with the median across the 11 GCMs and 25-75 / 10-90% bands |
| `comparison/across_scenarios_<level>.png` | one level: 5 metric panels, one ensemble-median line per SSP (with a light 10-90% band) |
| `summary_2045.png` | final-year summary: 4 levels x 5 metrics, per-scenario median with 25-75% and 10-90% spreads |
| `figure_inventory.csv` | manifest (category, scenario, gcm, metric, level, path, status) |

The shipping ensemble yields 65 PNGs: 44 `per_run` (4 SSPs x 11 GCMs), 16
`ensemble` (4 SSPs x 4 levels), 4 `comparison` and 1 final-year summary. Files
are written at ~200 dpi (`px_per_unit = 2`).

### Diagnostic figures

`src/climate_diagnostics.jl` adds three explanatory figures (again from the
exported CSVs only) that answer *why* the ensemble behaves as it does:

```bash
julia --project=. scripts/plot_climate_diagnostics.jl [results_root] [figures_dir]
```

| path | content |
| :--- | :--- |
| `diagnostics/forcing_and_response.png` | warming forcing vs richness/biomass response, final-year warming-response scatter, extinction-risk trajectory |
| `diagnostics/thermal_niches.png` | native species thermal optima vs the site-temperature distribution and the end-of-horizon warming |
| `diagnostics/community_filling.png` | initial vs final per-site native-richness distribution and the site-by-site change |

They are generated automatically by `run_climate_scenarios.jl` when
`make_figures = true`. See
`docs/climate_scenarios_results_report.md` for their interpretation.

Run discovery and completeness detection mirror `runs_index.csv`: a run whose
year range does not cover the reference horizon (for example the 2-year smoke
run `GUADEX_CLIMATE_END_YEAR=2027`) is flagged `[INCOMPLETE]` in its per-run
title and excluded from every ensemble band and summary, so a truncated run can
never bias the across-GCM statistics.

## Other entry points

The four-level + viewer export is also wired into:

* `scripts/run_model.jl` - single management scenario.
* `run_sensitivity_report.jl` - temperature x upstream cost x passability sweep.
* `run_alt_interactions.jl` - alternative interaction matrices.
* `run_climate_scenarios.jl` - all climate scenarios.

For an already-saved JLD2 file:

```bash
julia --project=. scripts/export_run_outputs.jl <simulation_output.jld2> [output_dir]
```

This re-export uses the crosswalk for water bodies; connectivity metrics are only
included when the JLD2 stores `dams` (the in-memory exports during a run always
include them).

## Assumptions and limitations

1. **Uniform basin warming by default.** The projections cover 6 sites; the
   runner uses the across-site median anomaly. Per-site spatial structure is
   only represented when `elevation_scaling` is enabled.
2. **Projection windows, not daily GCM series.** The shipped per-GCM table holds
   window-mean anomalies, so the annual curve linearly interpolates window
   midpoints. The ensemble-median daily series lives in
   `guadex_tw/data/interim/tw_future_ensemble_daily.parquet` and is not used
   per GCM here.
3. **Stationarity** of the calibrated Tw-Ta relation is assumed, as in the
   guadex_tw report.
4. **Habitat deterioration** is a placeholder: `habitat_degradation_index =
   1 - habitat_suitability`, where suitability is the IET-derived index already
   used by the model. A dedicated deterioration input has not been defined.
5. **Monthly saves** are used for reporting. The integration itself remains
   daily; set `save_interval_days` lower if higher temporal resolution is needed.

## Modelling-improvement features (WP0-WP6)

The features below implement the modelling-improvement plan
(`.kilo/plans/model-improvement-plan.md`). Every one is opt-in and the defaults
reproduce the previous behaviour, so old runs remain reproducible.

### WP0 — species classification

`parameters.toml` classifies `ST` (Salmo trutta) as **native** and adds a third
`migratory` group (`AA`, `AAL`, `LR`, `MC`). ST is the cold-water keystone whose
empirical optimum (12 °C, range 4-20 °C) is below current basin temperatures, so
it is the native species through which warming can appear as a decline. The
native richness metric therefore counts 10 species; migratory species are neither
native nor invasive and are reported separately.

### WP1 — per-site daily water temperature

`guadex_tw/scripts/13_project_guadex_sites.py` projects the calibrated spatial Tw
model at every GuadeX site coordinate (read from `data/ConnectivityUTM.csv`,
ETRS89/UTM 30N) for all 11 GCMs × 4 SSPs. It reads the Guadalquivir bounding box
of each 2 GB PNACC file remotely once per (GCM, experiment) — ~50× cheaper than
one read per site — caches the per-(GCM, experiment) arrays as `.npy`, and fills
short NaN gaps. It writes:

* `water_temp_daily_guadex_sites_wide.csv` — `date, scenario, <site columns>`,
  ensemble median across GCMs, for `historical` (1986-2005) and each SSP
  (2026-2045). This is what the runner reads.
* `water_temp_daily_guadex_sites_wide_<scenario>_<gcm>.csv` (with `--per-gcm`) —
  the same layout per GCM, so the ODE ensemble carries the GCM spread instead of
  the ensemble median.
* `water_temp_baseline_guadex_sites.csv` — per-site 1986-2005 baseline mean.
* `water_temp_daily_guadex_sites.parquet` (with `--long`) — the long
  `site_id, scenario, date, tw_ensemble_median` archival product.

The Julia loaders auto-detect the layout: `load_daily_forcing_any` /
`is_wide_daily_forcing` / `wide_forcing_matrix` for the wide format,
`load_daily_temperature_forcing` / `daily_forcing_matrix` for the long format.

### WP2 — daily vs annual-mean forcing

`[run_climate_scenarios] forcing_mode` selects:

* `"annual_mean"` (default) — the historical annual-node warming curve;
* `"daily"` — per-site daily forcing from WP1, with anomalies taken against the
  fixed `baseline_period_start`..`baseline_period_end` window instead of anchoring
  zero at the start year.

Set `daily_forcing_per_gcm = true` to read the per-GCM files
(`<base>_<scenario>_<gcm>.csv`) so each GCM run is differenced against its own
historical baseline; otherwise all GCMs in a scenario share the ensemble-median
forcing and produce identical trajectories.

`daily_temperature_schedule`, `baseline_climatology_schedule`,
`annual_mean_deltas` and `annual_mean_deltas_by_year` in
`src/temperature_forcing.jl` build the schedules; `TemperatureSchedule` already
interpolates arbitrary node spacing, so a daily schedule is just
`days_per_year = 1`. The export still reports annual-mean anomalies, and a daily
schedule with zero anomaly reproduces the static model exactly (unit-tested).

### WP3 — heat-stress mortality

The ODE gains a per-capita loss active only above each species' empirical upper
thermal limit:

```
m_heat,s(t) = k · max(0, T_i(t) − T_upper,s)²      # 1/day
dU[i,s] = N_is·(r_eff·logistic + interaction) − N_is·m_heat + dispersal
```

`T_upper,s` is the trait-table upper bound (parsed by
`parse_temperature_range_and_sigma`), and `k` is a **single shared slope**, not 24
parameters. Configure it in `[temperature_stress]`:

```toml
[temperature_stress]
enabled = false        # default: heat stress off (pre-WP3 model)
k = 0.0                # 1/day/degC^2
calibrate = false      # derive k from the baseline exceedance energy
max_annual_loss = 0.05 # baseline-viability constraint used by `calibrate`
limits_source = "trait_table"
```

`calibrate_heat_stress_rate` sets `k` so the worst baseline species/site keeps at
least `1 - max_annual_loss` annual survival. The thermal optimum (E5) is swept
separately with `thermal_optima_fraction` (0 = cold edge, 0.5 = midpoint,
1 = warm edge) via `optimum_sweep_optima`; this is not a fit. Narrowing σ is out
of scope.

### WP4 — equilibrium start and rebased metrics

`[run_climate_scenarios] spin_up` integrates to steady state under baseline
forcing before any scenario and reuses the state for every run;
`spin_up_max_years` / `spin_up_tol` control the stop criterion and
`spin_up_criterion` selects the measure (`"basin"` = basin-total relative annual
change, the robust default; `"q95"` = 95th percentile of the per-site change;
`"max"` = strict legacy max-over-sites, which a single near-empty site can make
unreachable). The run reports `converged = false` when the cap is hit;
`spin_up_progress_every` prints a progress line every N year-blocks. In daily mode
the spin-up uses the **seasonal** baseline climatology
(`spin_up(...; schedule=...)`), because the annual-mean and seasonal equilibria
differ substantially (the E1 gate: end biomass 759 vs 678) and only the seasonal
one matches the runs.

Carrying capacity has **two** multipliers:
`carrying_capacity_base_scaling` maps observed total density to the base K inside
`build_carrying_capacity` (legacy default `10.0`, i.e. K = 10× observed), and
`carrying_capacity_scaling` is the WP4 sensitivity multiplier applied on top
(`1×`, `3×`, `10×`). The effective multiplier against observed density is their
product, and both are recorded in `run_metadata.json`. Setting
`carrying_capacity_base_scaling = 1.0` treats the observed snapshot as the
capacity level. When spin-up is on, `native_richness_relative`,
`native_extinction_risk` and the new biomass ratios are rebased on the spun-up
state rather than the t = 0 snapshot, which was not an equilibrium.

### WP5 — abundance, occupancy and quasi-extinction

`compute_species_metrics` writes `levels/species_timeseries.csv` (per site ×
species density, presence, relative density, quasi-extinction flag and
time-to-quasi-extinction), and `quasi_extinction_summary` writes
`levels/quasi_extinction_summary.csv`. A species is **quasi-extinct** at a site
when its density is below `max(presence_threshold, q · baseline_density)` for
`quasi_extinction_persistence` consecutive annual snapshots (`q =
quasi_extinction_q`, default 0.1). Site metrics add `native_occupancy`,
`invasive_occupancy`, `total_occupancy`, `native_biomass_relative`,
`total_biomass_relative`, `native_quasi_extinct` and
`native_quasi_extinct_fraction`. When daily forcing is available,
`levels/exposure_sites.csv` reports days above each species' upper limit and the
annual squared exceedance energy (`exposure_table`).

### WP6 — staged experiments

`run_climate_experiments.jl` runs the staged design (E0 spin-up realism, E1
seasonality, E2 heat stress, E3 K sensitivity, E4 optimum sweep) and writes one
case per directory plus `results/climate_experiments/stage_index.csv`; E5 (the
44-run ensemble plus control) is delegated to `run_climate_scenarios.jl` with the
daily/heat-stress/spin-up flags. Settings live in
`[run_climate_experiments]`; `GUADEX_EXPERIMENTS_STAGES=E0,E2` selects stages.

### New environment overrides

| variable | meaning |
| :--- | :--- |
| `GUADEX_CLIMATE_FORCING_MODE` | `annual_mean` (default) or `daily` |
| `GUADEX_CLIMATE_DAILY_FILE` | per-site daily forcing CSV (WP1) |
| `GUADEX_CLIMATE_BASELINE_START` / `_END` | daily anomaly baseline window |
| `GUADEX_CLIMATE_SPIN_UP` | `1` to start from the spun-up equilibrium |
| `GUADEX_CLIMATE_SPIN_UP_YEARS` / `_TOL` | spin-up stop criterion |
| `GUADEX_CLIMATE_K_SCALING` | carrying-capacity multiplier |
| `GUADEX_CLIMATE_OPTIMA_FRACTION` | thermal-optimum sweep position |
| `GUADEX_CLIMATE_HEAT_STRESS` | `1` to enable WP3 |
| `GUADEX_CLIMATE_HEAT_STRESS_K` | shared slope `k` |
| `GUADEX_CLIMATE_HEAT_STRESS_CALIBRATE` | calibrate `k` from baseline forcing |
| `GUADEX_CLIMATE_QE_Q` / `_PERSISTENCE` | quasi-extinction definition |
| `GUADEX_CLIMATE_CONTROL` | `1` adds a no-warming control run (`control`/`baseline`) |
| `GUADEX_CLIMATE_OUTPUT_DIR` | results root (default `results/climate_scenarios`) |
| `GUADEX_CLIMATE_DAILY_PER_GCM` | `1` uses the per-GCM daily files |

## Validation

* `test/test_outputs.jl` covers the crosswalk mapping, per-site metrics,
  four-level aggregation, the temperature schedule, the warming curve, the
  end-to-end file export, and the WP5 species/quasi-extinction metrics.
* `test/test_temperature_forcing.jl` covers the empirical limits, the optimum
  sweep, heat-stress calibration, exposure diagnostics, the daily/annual-mean
  schedule builders, the WP1 loaders and the WP4 spin-up/K helpers.
* `test/test_ode.jl` checks that the scheduled ODE reproduces the static model
  when the anomalies are zero, responds correctly to warming, and that the
  heat-stress term is zero below the upper limit and monotone above it.
* `test/test_parameters.jl` checks the WP0 classification.
* `test/test_climate_figures.jl` covers run discovery, incomplete-run detection,
  per-run level statistics, across-GCM ensemble quantiles and the diagnostic
  `read_climate_basin_series` reader from synthetic CSVs (no figures are
  rendered in tests).
* `GUADEX_CONFIG_ONLY=1` on any entry script validates the configuration without
  loading data.
