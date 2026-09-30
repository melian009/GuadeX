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
(richness, biomass, relative native-richness loss (a deterministic realised
contraction, not a probability),
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

### Realised forcing vs the warming-end proxy (C7)

`warming_end_degc` in the run index and `run_metadata.json` is a **proxy**: it is
the last point of `basin_warming_curve`, computed from the 6-gauge
`water_temp_future_2045.csv` window means. The simulations, however, are forced
with the daily per-GCM series at the 775 sites, so the applied warming differs.

Each run therefore also persists the **realised** mean applied warming anomaly
over every site, derived from the annual-mean `warming` matrix the ODE is forced
with:

| field | meaning |
| :--- | :--- |
| `realised_warming_2026_2045_mean_degc` | site-and-window mean over 2026-2045 |
| `realised_warming_2036_2045_mean_degc` | site-and-window mean over 2036-2045 |
| `realised_warming_n_sites` | number of sites the mean covers (775) |
| `realised_warming_source` | `daily_schedule` or `annual_mean_schedule` |
| `warming_end_degc` | proxy, retained for continuity (`warming_end_degc_source = "proxy_basin_warming_curve"`) |

These fields are written to `<run_dir>/export/run_metadata.json` and to the
`runs_index.csv` columns. The dose-response analyses use the realised column as
the regressor (`Guadex.realised_forcing_axis`) and fall back to
`warming_end_degc` only when it is absent, labelling the fallback explicitly.
Runs written before the C7 change have `NaN` in the realised columns until they
are re-run.

### Re-selecting the cool/warm climate endpoints (C7)

The two endpoints of the 44-run ensemble used as the background climate for the
obstacle/upstream-cost sweeps live in `[obstacle_sensitivity].climate_models`.
They must not be hard-coded from memory: `Guadex.select_climate_endpoints`
derives them from the realised forcing (falling back to `warming_end_degc`, then
to the configured pair) and reports which path it used. On the next re-run:

```julia
using DataFrames, CSV, Guadex
idx = CSV.read(joinpath("results", "climate_scenarios", "runs_index.csv"), DataFrame)
endpoints = Guadex.select_climate_endpoints(idx)
println(endpoints.source)          # "realised" once the columns are populated
println(endpoints.cool.scenario, "/", endpoints.cool.gcm)
println(endpoints.warm.scenario, "/", endpoints.warm.gcm)
```

Put those `scenario`/`gcm` pairs into `[obstacle_sensitivity].climate_models`
(or the `GUADEX_SENSITIVITY_CLIMATE_MODELS` override). For a sensitivity index
the scenario column is `climate_scenario`:
`select_climate_endpoints(idx; scenario_col=:climate_scenario)`.

### Dose-response with GCM structure (C7)

The ensemble is clustered by GCM, so a pooled OLS slope over all runs is
optimistic. `Guadex.dose_response(x, y, groups)` fits the pooled slope and R²,
then adds GCM as a categorical fixed effect (`y ~ x + GCM`, which absorbs the
between-GCM level shifts, the fixed-effect analogue of a random intercept) and
reports the per-GCM slope spread (standard deviation, min/max) and the
between-GCM spread of the mean response. `MixedModels` is deliberately **not**
added: it is not a dependency, and the fixed-effect slope plus the per-GCM spread
carries the same information for this diagnostic. The report figures annotate the
pooled slope, the GCM spread and R².

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

### Corrected re-run configuration

The corrected re-run (post the C1/C2/C5/C7/E4 fixes) is driven by the committed
configuration `parameters_climate_scenarios_corrected.toml`:

```bash
GUADEX_PARAMETERS_FILE=parameters_climate_scenarios_corrected.toml \
  julia --project=. run_climate_scenarios.jl
# the matching climate-consistent obstacle/interaction sweeps:
GUADEX_PARAMETERS_FILE=parameters_climate_scenarios_corrected.toml \
  julia --project=. run_sensitivity_report.jl
```

It differs from the legacy September reproduction in exactly the settings that
the corrections changed: the E4 two-criteria spin-up stop rule
(`spin_up_criterion = "both"`, `spin_up_composition_tol = 1e-2`,
`spin_up_min_years = 10`) in `[run_climate_scenarios]` and
`[obstacle_sensitivity]`, a heat-stress slope derived from the corrected
baseline (`[temperature_stress].calibrate = true`), and the corrected output
root `results/climate_scenarios_corrected` (a new directory, so the
pre-correction September outputs under `results/climate_scenarios` are preserved
for comparison). The E4 composition statistic is the **robust** 95th percentile
of `|Δ log N|` over the site × species cells whose density exceeds the 0.1-unit
presence floor in at least one of the two compared years; the raw all-cells
statistic is recorded alongside it as `spin_up_composition_q95_change_all_cells`
(a diagnostic only — it plateaus near 1.4e-2 on real data and never reaches the
old 1e-3 tolerance). The corrected per-site water-temperature
source (`guadex_tw/outputs/tables/water_temp_baseline_guadex_sites.csv`, exact
column `tw_baseline_mean`) and the corrected interaction form
(`interaction_inside_growth = true`) are code defaults with no configuration
opt-in, so this file uses them by construction. The generic `parameters.toml`
defaults are left unchanged for the other scripts and the test suite.

The corrected file also makes the **team's decisions** explicit:

* **C3 interim projection route** — `projection_route = "interim_observed"` with
  `interim_spin_up_years = 3` (see below);
* **the six biological-option values** (E13–E16, E18) — `absence_growth_fraction
  = 1.0`, `pool_capacity_mode = "cap"`, `fishless_dificil_capacity =
  "near_zero"`, `nonreproducing_local_growth = "zero"`,
  `exclude_fishfarm_eel_records = true`, `salinity_envelope = true` (see
  [Explicit biological-assumption options](#explicit-biological-assumption-options-e13-e16-e18));
* **the uniform obstacle-passability assumption** (the IF note below).

### C3 interim projection route

The team chose the **interim** route over full simulation-based calibration
(Correction Plan decision #1): start projections from the **observed** community
with only a **short** spin-up and report every scenario relative to the matched
no-warming control, instead of starting from a long model-generated burn-in.
The route is an explicit, opt-in setting in `[run_climate_scenarios]`:

```toml
projection_route = "interim_observed"   # or "full_burnin" (default when absent)
interim_spin_up_years = 3               # allowed 0-5; N = 3 here
```

* `"full_burnin"` (the default, and the only behaviour when the key is absent)
  integrates to the E4 two-criteria stop rule described above. **The E4
  convergence machinery is not disabled for other configurations**; the interim
  route bypasses it only for the config that selects it.
* `"interim_observed"` starts from the observed density matrix and integrates
  **exactly** `interim_spin_up_years` one-year baseline blocks
  (`Guadex.interim_observed_spin_up`, `src/spin_up.jl`). `N = 3` is a defensible
  short relaxation: it lets the fastest local/compositional transients (e.g. from
  the presence penalty and the pool-cap change) relax, while the state stays
  anchored to the observations. Any residual transient common to a scenario and
  the control cancels because results are reported as `scenario − control`.
  `N = 0` is allowed and starts exactly from the observed snapshot.
* The route **requires `control_run = true`** (the runner errors otherwise), and
  the route and the number of spin-up years are written to
  `export/run_metadata.json` (`projection_route`, `interim_spin_up_years`,
  `spin_up_years`).
* Results under this route are **provisional** (pending full calibration): they
  are reported as per-metric deltas via `Guadex.scenario_minus_control` /
  `Guadex.scenario_minus_control_table` (`src/climate_diagnostics.jl`), which
  writes `figures/diagnostics/scenario_minus_control.csv` from the figures path.

Environment overrides: `GUADEX_CLIMATE_PROJECTION_ROUTE`,
`GUADEX_CLIMATE_INTERIM_SPIN_UP_YEARS`.

> The interim route's status as an uncalibrated, conditional projection is stated for the
> October report in `docs/GuadeX_Limitations_and_Scope_Oct2026.md` (item 2).

### Obstacle passability and the `IF` index

Barrier passability is currently modelled using a uniform structural assumption
(0.1 upstream, 0.5 downstream) pending full standardisation of the regional
obstacle inventory's qualitative passability index (IF).

> The standing IF/passability caveat is carried into the October report in
> `docs/GuadeX_Limitations_and_Scope_Oct2026.md` (item 4).

### Committed September configuration (LEGACY reproduction)

The September 2026 ensemble under `results/climate_scenarios_k1x_burnin/`
(11 GCMs x 4 SSPs plus a no-warming control) is reproduced by
`legacy/parameters_climate_scenarios_k1x_burnin.toml`. **This file is a legacy
reproduction only**: it deliberately pins the pre-E4 `spin_up_criterion =
"basin"` rule and the September heat-stress slope, and must **not** be used for
the corrected results (use `parameters_climate_scenarios_corrected.toml`). It is
a copy of `parameters.toml` with only two sections changed, so the default
configuration (and the tests that depend on it) is untouched:

```bash
GUADEX_PARAMETERS_FILE=legacy/parameters_climate_scenarios_k1x_burnin.toml \
  julia --project=. run_climate_scenarios.jl
```

The committed settings (read back from the September
`export/run_metadata.json`) are:

| setting | value | section |
| :--- | :--- | :--- |
| `forcing_mode` | `"daily"` | `[run_climate_scenarios]` |
| `daily_forcing_per_gcm` | `true` (per-GCM series) | `[run_climate_scenarios]` |
| `baseline_period_start` / `_end` | `1986` / `2005` | `[run_climate_scenarios]` |
| `spin_up` | `true` | `[run_climate_scenarios]` |
| `spin_up_max_years` / `_tol` / `_criterion` | `600` / `5e-5` / `"basin"` (converges at **373 years**) | `[run_climate_scenarios]` |
| `carrying_capacity_base_scaling` / `_scaling` | `1.0` / `1.0` (K = 1x observed) | `[run_climate_scenarios]` |
| `control_run` | `true` | `[run_climate_scenarios]` |
| `start_year` / `end_year` | `2026` / `2045` (20 years) | `[run_climate_scenarios]` |
| `output_dir` | `"results/climate_scenarios_k1x_burnin"` | `[run_climate_scenarios]` |
| `enabled` | `true` | `[temperature_stress]` |
| `k` | `3.7817692640400324e-5` (`calibrate = false`, pinned) | `[temperature_stress]` |

`k` is pinned with `calibrate = false` to the value the September run recorded,
so the configuration reproduces it exactly without depending on the calibration
data. The default `[temperature_stress]` (`enabled = false`) is unchanged.

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

A completed run is **reused only when its stored digest matches** the current
one; otherwise it is recomputed. The digest (`parameter_digest`) covers the
resolved model state (interaction matrix, growth/temperature/carrying-capacity/
thermal arrays, heat-stress slope, dispersal matrix and scaling, the
`interaction_inside_growth` flag), the resolved run options (forcing mode,
baseline period, burn-in years/tolerance/criterion, K scaling, upstream cost,
obstacle/passability configuration, scenario, GCM and the warming curve) and the
code version (`git rev-parse --short HEAD`, or `"unknown"` when git is
unavailable). It is written to `export/run_metadata.json`
(`parameter_digest`, `code_version`, and the full `resolved_config`) and to the
`run_digest.txt` sidecar next to `simulation_output.jld2`. Outputs produced
before this change carry no digest and are therefore re-run once. Set
`GUADEX_CLIMATE_FORCE=1` to recompute regardless, and note that
`runs_index.csv` is rewritten after each run so an interrupted ensemble can be
resumed without losing the summary of finished runs.

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
<run_dir>/export/levels/exposure_sites.csv
<run_dir>/export/levels/exposure_summary.csv
<run_dir>/export/viewer/guadex_results_<metric>_timeseries.csv
<run_dir>/export/viewer/guadex_results_metrics.json
<run_dir>/export/viewer/guadex_results_realised_richness_loss.json
<run_dir>/export/viewer/level_subcatchment_mean_<metric>_timeseries.json
<run_dir>/export/viewer/level_water_body_mean_<metric>_timeseries.json
<run_dir>/export/viewer/level_basin_mean_<metric>_timeseries.json
<run_dir>/export/viewer/level_<level>_mean_<metric>_by_site_timeseries.json
<run_dir>/export/viewer/species_<sp>_<metric>_timeseries.json
```

### Viewer formats

Each `guadex_results_<metric>_timeseries.csv` is keyed by `CODIGO` with a
`step` column (year) and a single `value` column, so the explorer builds one
proper time series per metric. (An all-metrics wide CSV would instead be read as
one key per metric × year, mixing the two axes in a single slider.)
`guadex_results_<metric>.json` matches the single-variable time-series format in
`viz/README.md`:

```json
{ "name": "...", "unit": "fraction", "steps": ["2026", "2027"],
  "data": { "1.1.2": [0.0, 0.12] } }
```

`guadex_results_metrics.json` matches the per-site metrics format:

```json
{ "name": "...", "data": { "1.1.2": { "native_richness": 3, "total_biomass": 63.1 } } }
```

The per-level time-series JSON files use the same shape but are keyed by
subcatchment, water body or `ES050`.  `level_<level>_mean_<metric>_by_site_timeseries.json`
repeats each aggregate onto its member sites, so the existing site layer renders
it directly; `species_<sp>_<metric>_timeseries.json` provides one site-keyed
time series per metric per species.

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
| `per_run/<scenario>__<gcm>.png` | one run: 5 metrics (rows) x 4 levels (columns); nested levels show the across-unit mean with a 10-90% band (site spread at the sampling point) |
| `ensemble/ensemble_<scenario>_<level>.png` | one scenario + level: 5 metric panels with the mean across GCMs of the level mean and 25-75 / 10-90% bands across GCMs |
| `comparison/across_scenarios_<level>.png` | one level: 5 metric panels, one mean line per SSP (with a light 10-90% across-GCM band) |
| `summary_2045.png` | final-year summary: 4 levels x 5 metrics, per-scenario mean with 25-75% and 10-90% spreads across GCMs |
| `figure_inventory.csv` | manifest (category, scenario, gcm, metric, level, path, status) |

Ensemble statistics aggregate **per run before pooling**: each run is first
reduced to its mean across units (sites, sub-basins, water bodies or the single
basin) and those run-level means are then summarised across runs. The line is the
mean across GCMs and the bands are the 25-75 / 10-90% across GCMs, so nested
levels show the aggregate trend with a tight, readable scale. Pooling raw units
instead pins the central value to the modal unit for discrete or zero-inflated
metrics (richness, realised richness loss), which is why sub-basin / water-body / site
panels previously rendered flat; `run_level_stats` still exposes the within-run
across-unit spread used by the per-run figures.

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
| `diagnostics/forcing_and_response.png` | warming forcing vs richness/biomass response, final-year warming-response scatter, realised richness-loss trajectory |
| `diagnostics/thermal_niches.png` | native species thermal optima vs the site-temperature distribution and the end-of-horizon warming |
| `diagnostics/community_filling.png` | initial vs final per-site native-richness distribution and the site-by-site change |

They are generated automatically by `run_climate_scenarios.jl` when
`make_figures = true`. See
`legacy/climate_scenarios_results_report.md` for their interpretation.

Run discovery and completeness detection mirror `runs_index.csv`: a run whose
year range does not cover the reference horizon (for example the 2-year smoke
run `GUADEX_CLIMATE_END_YEAR=2027`) is flagged `[INCOMPLETE]` in its per-run
title and excluded from every ensemble band and summary, so a truncated run can
never bias the across-GCM statistics.

## Other entry points

The four-level + viewer export is also wired into:

* `scripts/run_model.jl` - single management scenario.
* `run_sensitivity_report.jl` - obstacle/passability x upstream-cost sweep under two extreme climate projections (see below).
* `run_alt_interactions.jl` - the same sweep with alternative interaction matrices and a thermal-sigma multiplier.
* `run_climate_scenarios.jl` - all climate scenarios.

For an already-saved JLD2 file:

```bash
julia --project=. scripts/export_run_outputs.jl <simulation_output.jld2> [output_dir]
```

This re-export uses the crosswalk for water bodies; connectivity metrics are only
included when the JLD2 stores `dams` (the in-memory exports during a run always
include them).

## Obstacle / upstream-cost sensitivity (daily forcing)

`run_sensitivity_report.jl` and `run_alt_interactions.jl` answer *"what do
obstacles and upstream dispersal cost do to the fish community?"* under the same
climate-consistent setup as the `climate_scenarios_k1x_burnin` ensemble:

* per-site **daily** water temperature from the two extreme projections of the
  44-run ensemble (`ssp126/IITM-ESM` and `ssp245/UKESM1-0-LL`; the climate
  ensemble index reports 0.28 / 1.37 °C basin-median warming by 2045);
* a seasonal baseline **burn-in** (cached under `results/sensitivity_spinup_cache`
  and reused across runs);
* calibrated heat stress, K = 1x observed, WP0 classification (ST native,
  migratory group reported separately).

Shared settings live in `[obstacle_sensitivity]`; the swept grids in
`[run_sensitivity_report]` (upstream cost × passability) and
`[run_alt_interactions]` (adds interaction matrix × thermal sigma). Each run is
written to
`results/sensitivity_obstacles[/alt_interactions]/<matrix>/[sig_<sigma>/]<model>/uc_<cost>/<pass>/`
with the four-level + viewer export, plus a `runs_index.csv` summary and a
`sensitivity_effects.csv` / `summary_plots/` produced by:

```bash
julia --project=. scripts/plot_obstacle_sensitivity.jl results/sensitivity_obstacles
```

In the sensitivity index, `warming_end_degc` is the **maximum site-level**
annual-mean anomaly of the forcing actually delivered to the ODE (it differs
from the basin-median curve reported by the climate ensemble index).

Resuming is safe: a run counts as complete only when its full export is present
**and** a `.sensitivity_complete.jld2` marker (or, for runs written before
markers existed, the stored JLD2) matches the current science-affecting settings
— horizon, upstream cost, passability, K, optima fraction, interaction matrix and
thermal sigma. A rerun under changed settings is recomputed rather than reused or
relabelled. `[obstacle_sensitivity].spin_up = false` skips the burn-in and starts
from observed densities.

Configuration smoke test (no data loading):

```bash
GUADEX_CONFIG_ONLY=1 julia --project=. run_sensitivity_report.jl
```

Useful overrides: `GUADEX_SENSITIVITY_END_YEAR`, `GUADEX_SENSITIVITY_MAX_RUNS`,
`GUADEX_SENSITIVITY_FORCE`, `GUADEX_SENSITIVITY_OUTPUT_DIR`,
`GUADEX_SENSITIVITY_CLIMATE_MODELS` (`scenario:gcm,scenario:gcm`),
`GUADEX_SENSITIVITY_K_BASE`, `GUADEX_SENSITIVITY_SPIN_UP`, `GUADEX_SPINUP_REUSE`,
`GUADEX_SPINUP_FORCE`. Completed runs are skipped unless
`GUADEX_SENSITIVITY_FORCE=1`, so an interrupted sweep resumes where it stopped.
The burn-in cache key covers the upstream cost, obstacle configuration, forcing
file size/mtime, initial state **and the full parameter digest** (model arrays,
heat-stress slope, run options and code version), so changing any model
parameter, the code version or the configuration recomputes the burn-in instead
of silently reusing a stale equilibrium. Completed sensitivity runs are stored
with a completion marker containing the run fingerprint, which now also carries
that digest.

### Within-design ranking vs global sensitivity (E7)

`scripts/plot_sensitivity_effects.jl` produces `parameter_effect_rank.csv` and
`variance_decomposition.csv`.  These are **within-design** summaries: a balanced
full-factorial sweep lets a sum-of-squares decomposition answer *"which swept
factor moved the response most at the chosen levels?"* inside that sweep.  They do
**not** generalise across the full parameter ranges, and the thermal-sigma axis is
confounded across the two sweeps.

The complementary **global** design lives in `src/global_sensitivity.jl` and is
run by `scripts/global_sensitivity_sobol.jl` (Latin hypercube + Saltelli
A/B/AB/BA, first- and total-order Sobol indices with bootstrap confidence
intervals) over the agreed ranges in `[global_sensitivity]`.  It answers *"how
much of the response variance is attributable to each parameter and its
interactions over its range?"*.  Both are kept: the within-design ranking remains
valid for the balanced sweep; the Sobol index is the global answer.  See
`docs/global_sensitivity_methods_note.md` for the ranges, the Ishigami estimator
validation, the pilot and the production compute estimate.

## Assumptions and limitations

> **Consolidated limitations (October 2026 report).** The report's Limitations and Scope of
> Inference are consolidated in **`docs/GuadeX_Limitations_and_Scope_Oct2026.md`**. That
> document reflects the team's standing decisions on **C6** (distributional thermal niches),
> **C3** (the interim observed-start projection route), **E21** (the same-day air–water
> calibration and its lagged uncertainty band), the obstacle **IF**/passability assumption and
> **E8** (the deterministic realised richness-loss metric). The items below are the
> runner-specific subset; the consolidated document is authoritative for the report text.

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

## Explicit biological-assumption options (E13-E16, E18)

The review flagged a set of biological assumptions that were previously baked
into the code.  They are now explicit, **opt-in** options under
`[biological_options]` in `parameters.toml`, read through
`SimulationParameters.biological_options` and exposed by `prepare_ode_data`.
**Every option defaults to the current behaviour**, so the shipped configuration
and the test suite reproduce the historical outputs exactly (a real-data
regression test pins the default growth and carrying-capacity arrays).  The
resolved values are written to `run_metadata.json` (and the flat
`biological_options` keys of `resolved_config`) and folded into the parameter
digest, so changing any option invalidates cached burn-ins and completed runs.

| option | default | alternative | where implemented | trade-off |
| :--- | :--- | :--- | :--- | :--- |
| **E13** `absence_growth_fraction` | `0.1` | `1.0`, or `0.0` | `build_intrinsic_growth_rates` | Fraction of the local rate given to a species observed absent. `1.0` removes the presence penalty so temperature/habitat decide establishment (the review's recommendation); `0.0` forbids local growth at absent sites. Raising it makes rare/invasive spread easier; the default was kept. |
| **E14** `pool_capacity_mode` | `"legacy"` | `"cap"` / `"exclude"` | `build_carrying_capacity` | 84 `EN_POZAS = Si` isolated-pool sites hold ~56% of total density and drive the high-K tail (max K ≈ 95,333 → ≈ 8,057 under either alternative). `cap` clips a pool site's K at the non-pool median; `exclude` replaces it with that median. Both only touch pool sites. |
| **E15** `fishless_dificil_capacity` | `"legacy"` | `"near_zero"` | `build_carrying_capacity` | At the 56 `SIN_PECES = DIFICIL` sites the field team judged unable to hold fish, set K = 1e-6 (the ODE's own floor), making them transit-only nodes. **Interacts with E13**: raising `absence_growth_fraction` while zeroing K pulls establishment in opposite directions. |
| **E16** `nonreproducing_local_growth` | `"legacy"` | `"zero"` | `build_intrinsic_growth_rates` | `AA`, `LR`, `MC` have `REPRODU_WITHIN_THE_BASIN = no` yet receive a local r. `zero` sets local r = 0; recruitment is to be represented later as estuarine immigration. Without that term the species only decays locally. |
| **E16** `exclude_fishfarm_eel_records` | `false` | `true` | `drop_fishfarm_eel_records` (called in `prepare_ode_data`) | Drops the four documented Guadiato fish-farm eel records (subcatchment 14.0; AA presences 9 → 5) before any builder reads the density table. Reversible: the input table is copied, not mutated. |
| **E18** `salinity_envelope` | `false` | `true` | `apply_salinity_envelope` (called in `prepare_ode_data`) | Zeros the local growth of the strictly brackish species (`AB`, `LR`, `MC`) at sites below 1000 uS/cm `CONDUCTIVIDAD` (1,584 cells on the current data). The euryhaline `brackish/freshwater` species (`AA`, `GH`, `CC`, `CG`) are left unchanged because the observations place them across the conductivity range. The envelope is a local growth multiplier, mathematically equivalent to a per-species habitat factor, so the ODE is untouched. |

> **No elevation envelope.** Only a salinity/brackish envelope is offered.  A
> static elevation envelope is deliberately **not** implemented: it would block
> upslope range shifts under warming (the review's own caveat, recorded in
> `docs/GuadeX_Correction_Plan_Sept2026.md` §0.1).

**These are decision options, not calibrated values.** The team's chosen values
for the corrected re-run (Correction Plan decision #6) are set in
`parameters_climate_scenarios_corrected.toml`:

| option | team value | effect |
| :--- | :--- | :--- |
| **E13** `absence_growth_fraction` | `1.0` | temperature/habitat decide establishment (no absent-site penalty) |
| **E14** `pool_capacity_mode` | `"cap"` | clip the 84 isolated-pool sites at the non-pool median capacity |
| **E15** `fishless_dificil_capacity` | `"near_zero"` | the 56 `SIN_PECES = DIFICIL` sites become transit-only (K = 1e-6) |
| **E16** `nonreproducing_local_growth` | `"zero"` | `AA`, `LR`, `MC` get local r = 0 |
| **E16** `exclude_fishfarm_eel_records` | `true` | drop the four documented Guadiato fish-farm eel records |
| **E18** `salinity_envelope` | `true` | zero strictly brackish `AB`/`LR`/`MC` below the conductivity threshold (no elevation envelope) |

The generic `parameters.toml` **retains the safe legacy defaults** (`0.1`,
`"legacy"`, `"legacy"`, `"legacy"`, `false`, `false`), so the historical/legacy
reproduction and the default test suite are unchanged; only a config that sets
these keys (the corrected file) opts in. Every value is still recorded in
`run_metadata.json` and folded into the parameter digest.

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
`realised_richness_loss` and the new biomass ratios are rebased on the spun-up
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
  per-run level statistics, across-GCM ensemble quantiles, the per-run-then-pool
  ensemble reduction (including a regression test that a skewed sub-basin
  distribution yields the level aggregate rather than the modal unit) and the
  diagnostic `read_climate_basin_series` reader from synthetic CSVs (no figures
  are rendered in tests).
* `GUADEX_CONFIG_ONLY=1` on any entry script validates the configuration without
  loading data.
