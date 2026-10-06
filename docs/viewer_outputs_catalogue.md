# Viewer output catalogue

This document catalogues every file class written by `Guadex.write_viewer_outputs`
(`src/outputs.jl`) into a run's `export/viewer/` directory, maps the current
machine filenames to their content, key namespace, time steps and viewer
behaviour, and records the reviewer's proposed descriptive naming as an
old -> proposed mapping. It exists to answer the annex item "Issues Found in the
Viewer and in Its Output Files" (Part A: A4, A5, A7, A8) with facts taken from
the corrected files on disk under `results/climate_scenarios_corrected/`,
`results/sensitivity_obstacles_corrected/` and `results/alt_interactions_corrected/`.

## (a) Naming convention, and why machine names are retained

Viewer filenames are machine-consumed keys: the explorer discovers them by
pattern and deep links point at them, so they are **not** renamed (A4). They
follow `guadex_results_<metric>_timeseries.csv`, `level_<level>_mean_<metric>_<variant>_timeseries.json`
and `species_<sp>_<metric>_timeseries.json`. Human-readable description is
carried by the `name` field **inside** each JSON (shown in the viewer legend);
as of the A4 fix that field is composed as
`"<run label> · <level/species> <human metric words>"`, and the run label is
derived from `run_metadata` (e.g. `scenario`, `stage`, `case`,
`exploitation_scenario`, `passability_scenario`, `gcm`) when no explicit
`viewer_name` is supplied. Filenames stay stable; descriptions become readable.

## (b) File patterns written per run (`export/viewer/`)

Counts are per run; totals below are on-disk counts across the corrected corpora.

| # | File pattern | Content | Key | Steps | Viewer use |
|---|---|---|---|---|---|
| 1 | `guadex_results_<metric>_timeseries.csv` | One metric as a long table `CODIGO,step,value`; 15 files, one per `VIEWER_METRICS` entry | `CODIGO` | 20 years (2026-2045) | Drawable (site layer) |
| 2 | `guadex_results_<primary_metric>.json` | Canonical single-variable series for the run's primary metric (`realised_richness_loss` by default) | `CODIGO` | 20 | Drawable |
| 3 | `guadex_results_metrics.json` | All 15 metrics at the final step, per site | `CODIGO` | final step only | Drawable (per-site metrics panel) |
| 4 | `level_<level>_mean_<metric>_timeseries.json` | Basin / sub-catchment / water-body **mean** of `native_richness`, `realised_richness_loss`, `total_biomass`; 9 files | group id: `ES050` (basin), `1.1` (sub-catchment), `ES050MSPF...` (water body) | 20 | **Not drawable by the site layer** (keys are groups, not sites); canonical aggregate |
| 5 | `level_<level>_mean_<metric>_by_site_timeseries.json` | The same group mean **repeated onto every member site**; 9 files | `CODIGO` | 20 | Drawable (A5 fix) |
| 6 | `species_<sp>_<metric>_timeseries.json` | Per-species `density`, `relative_density`, `present`, `quasi_extinct`; 24 species x 4 = 96 files | `CODIGO` | 20 | Drawable (A8 fix) |

Per run: `15 + 2 + 9 + 9 + 96 = 131` files (15 CSV + 116 JSON).

On-disk evidence: 45 viewer dirs / 5895 files in `results/climate_scenarios_corrected`,
32 / 4192 in `results/sensitivity_obstacles_corrected`, 48 / 6288 in
`results/alt_interactions_corrected` (each run = 131 files; species = 96/run,
`_by_site` = 9/run). No file named `*native_extinction_risk*` exists under any
`*_corrected` corpus.

## (c) Reviewer-proposed descriptive naming (old -> proposed)

Machine names are kept; the proposed scheme is recorded here so a *distribution*
copy (for humans) could be generated without breaking the pipeline. The reviewer's
example prefix `08_` is an ordering index within the numbered catalogue of this
document. Example run directory (real): `ssp126__IITM-ESM/uc_0.05/blocked`.

| Ending / class | Current machine name | Reviewer-style proposed descriptive name |
|---|---|---|
| primary JSON `.json` | `guadex_results_realised_richness_loss.json` | `08_ssp126_IITM-ESM_uc-0.05_barriers-blocked_richness-loss_2026-2045.json` |
| per-metric CSV `_timeseries.csv` | `guadex_results_native_richness_timeseries.csv` | `08_ssp126_IITM-ESM_uc-0.05_barriers-blocked_native-richness_2026-2045.csv` |
| per-site metrics `.json` | `guadex_results_metrics.json` | `08_ssp126_IITM-ESM_uc-0.05_barriers-blocked_all-metrics-final-step.json` |
| level aggregate `.json` | `level_water_body_mean_total_biomass_timeseries.json` | `08_ssp126_IITM-ESM_uc-0.05_barriers-blocked_water-body-mean_total-biomass_2026-2045.json` |
| level by-site `.json` | `level_water_body_mean_total_biomass_by_site_timeseries.json` | `08_..._water-body-mean_total-biomass_by-site_2026-2045.json` |
| species `.json` | `species_AA_density_timeseries.json` | `08_..._species-AA_density_2026-2045.json` |

The "20-run example" is the corrected sensitivity/obstacle corpus
(`results/sensitivity_obstacles_corrected/`, 32 run dirs = 8 scenario/GCM x
upstream-cost x 4 obstacle modes); each run's 131 files would receive the same
prefix scheme with its own index. The reviewer's proposed name is illustrative
and is **not** applied to machine output.

## (d) A5 separation guidance: aggregate vs site-keyed vs per-species

- **Viewer-intended (drawable):** classes 1, 2, 3, 5 and 6 above. Their `data`
  keys are site `CODIGO` values, so the 3-D site layer can match geometry and
  draw. The `_by_site` files (class 5) are the conversion of the aggregate
  levels onto site keys and are the reason the September "nothing drawn" symptom
  is resolved.
- **Canonical aggregate, not viewer-intended (class 4):** `level_<level>_mean_<metric>_timeseries.json`
  remains keyed by group id (`ES050`, `1.1`, `ES050MSPF...`). It is retained for
  provenance/plots and *cannot* be drawn by the site layer; the drawable version
  is class 5. These 9 files per run are the only non-site-keyed files in
  `viewer/`; the viewer counts them as "sites" but draws nothing, so a future
  viewer should either exclude class 4 or warn "0 of N keys match site codes"
  (a `viz/src` change, out of scope here).
- **Per-species (class 6):** keyed by `CODIGO`, drawable; no per-species totals
  are missing — `density`, `relative_density`, `present` and `quasi_extinct`
  are all emitted for every modelled species.

## Note on absent / stale files

The corrected pipeline does **not** emit any `native_extinction_risk` file or
metric (A7). Pre-correction artefacts with that name still exist on disk under
`results/climate_experiments/` (e.g. `.../E1/annual_mean/export/viewer/guadex_results_native_extinction_risk_timeseries.json`);
these are stale outputs of the earlier code, not regenerated by the current
pipeline, and are separate from the `*_corrected` corpora. The viewer's own demo
fixture `viz/public/data/results.demo-metrics.json` carries a neutral demo
field, `richness_loss_2100` (renamed from the former `extinction_risk_2100`);
it is a demo asset, not a model output. The demo time-series fixture's label was
likewise renamed to "realised richness loss".
