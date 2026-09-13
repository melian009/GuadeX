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

## Validation

* `test/test_outputs.jl` covers the crosswalk mapping, per-site metrics,
  four-level aggregation, the temperature schedule, the warming curve and the
  end-to-end file export.
* `test/test_climate_figures.jl` covers run discovery, incomplete-run detection,
  per-run level statistics and across-GCM ensemble quantiles from synthetic CSVs
  (no figures are rendered in tests).
* `test/test_ode.jl` checks that the scheduled ODE reproduces the static model
  when the anomalies are zero and responds correctly to warming.
* `GUADEX_CONFIG_ONLY=1` on any entry script validates the configuration without
  loading data.
