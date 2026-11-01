# GuadeX Phase 2 — P2.1a: Parser-Fix Effect on the Corrected Climate Ensemble

**Scope:** re-run the 45-run corrected climate ensemble against the regenerated interaction
table and quantify, not assert, what the interaction-parser fix changed.

**Status:** companion evidence for Phase 2 item P2.1 (`§4` of
`.kilo/plans/1793526105792-second-review-triage-plan.md`). This document does not revise the
report or the status record; that is P2.6.

---

## 1. Provenance and commands

All work was done in `E:\projects\GuadeX` (Windows, PowerShell 7, Julia 1.12.6) on 2026-11-01.

Pre-fix audit baseline (read-only, untouched):

- `results/climate_scenarios_corrected/runs_index.csv` — **46 lines** (header + 45 runs).
- Run metadata `code_version = "9702169"`; `parameter_digest = "f5f95873698951df"` for
  `ssp126/ACCESS-CM2` (`results/climate_scenarios_corrected/ssp126/ACCESS-CM2/export/run_metadata.json`).

New run:

```powershell
$env:GUADEX_PARAMETERS_FILE='parameters_climate_scenarios_corrected.toml'
$env:GUADEX_CLIMATE_OUTPUT_DIR='results/climate_scenarios_phase2'
$env:GUADEX_CLIMATE_PLOT='0'
julia --project=. run_climate_scenarios.jl *> results/phase2_climate_run.log
```

- Wall time: **≈23.1 min** (log `CreationTime` 13:25:56 → `LastWriteTime` 13:49:03 local).
- Exit status: exited normally; final log line **“Climate scenarios complete: 45 run(s)
  written under 'results/climate_scenarios_phase2'”**; no `ERROR`/`LoadError`/`failed`
  matches in `results/phase2_climate_run.log`.

Interaction-table diff and experiment comparison were computed with Julia
`DataFrames`/`CSV` over `data/BIOTIC/interaction_matrix_long_prefix.csv` and
`data/BIOTIC/interaction_matrix_long.csv`.

---

## 2. New-run verification

| check | result |
|---|---|
| `results/climate_scenarios_phase2/runs_index.csv` | **46 lines = header + 45 runs** |
| Scenario coverage | `ssp126` 11, `ssp245` 11, `ssp370` 11, `ssp585` 11, `control` 1 |
| GCMs | 11 GCMs + `baseline` control (`ACCESS-CM2`, `CMCC-CM2-SR5`, `CNRM-ESM2-1`, `EC-Earth3-Veg`, `IITM-ESM`, `KACE-1-0-G`, `MIROC6`, `MPI-ESM1-2-HR`, `MRI-ESM2-0`, `NorESM2-MM`, `UKESM1-0-LL`) |
| Required artifacts per run | 45/45 have `run_digest.txt`, `export/run_metadata.json`, `simulation_output.jld2` |
| `code_version` (all 45) | `c340398` (= `git rev-parse --short HEAD`) |
| `resolved_config.competition_structure` (all 45) | `shared_k` |
| `parameter_digest` | 45/45 **different** from the same run’s pre-fix digest; e.g. `ssp126/ACCESS-CM2` `f5f95873698951df` → `96a0c5621104ebb9` |
| `run_digest.txt` vs metadata `parameter_digest` | identical for all 45 runs |

The requested top-level key `competition_structure` is not emitted at the top level of
`run_metadata.json`; the value is recorded in `resolved_config.competition_structure`
(nested key, present in all 45 runs).

---

## 3. α-matrix before/after

Inputs: pre-fix `data/BIOTIC/interaction_matrix_long_prefix.csv` and corrected
`data/BIOTIC/interaction_matrix_long.csv`. Both have **227 directed rows** and the same
227 unique `(source,target)` keys. Columns compared: `source,target,mechanism,alpha,ambiguous,raw_text`.

### 3.1 The 16 reverse-direction cells

Each row below is the reverse direction of a directed interaction whose named direction stays
−0.5 predation. Old = pre-fix, new = corrected. All 16 were already `ambiguous = true`
before the fix; the fix corrects the value and mechanism.

| direction (source→target) | old α | old mechanism | new α | new mechanism |
|---|---|---|---|---|
| AAL→AM | −0.5 | predation | −0.3 | competition |
| AAL→CC | −0.5 | predation | −0.3 | competition |
| AAL→OM | −0.5 | predation | −0.3 | competition |
| CG→AAL | −0.5 | predation | −0.3 | competition |
| CG→AM | −0.5 | predation | −0.3 | competition |
| CG→CC | −0.5 | predation | −0.3 | competition |
| CG→EL | −0.5 | predation | −0.3 | competition |
| CG→LG | −0.5 | predation | −0.3 | competition |
| GH→AAL | −0.5 | predation | −0.3 | competition |
| GH→AM | −0.5 | predation | −0.3 | competition |
| GH→CC | −0.5 | predation | −0.3 | competition |
| GH→LG | −0.5 | predation | −0.3 | competition |
| LG→OM | −0.5 | predation | −0.3 | competition |
| TT→AM | −0.5 | predation | −0.3 | competition |
| TT→CC | −0.5 | predation | −0.3 | competition |
| TT→LG | −0.5 | predation | −0.3 | competition |

The mirror (named) direction of each pair (e.g. `AM→AAL`, `AAL→GH`) is unchanged at
−0.5 predation, as intended: the directed clause remains authoritative for the named
source→target.

### 3.2 The 9 newly flagged ambiguous pairs

These are the mixed “interfere through competition and predation” cells with no nominated
direction. The fix keeps −0.5 predation both ways as a documented judgement but sets
`ambiguous = true` on both directions. 9 unordered pairs → **18 directed rows**.

| pair | both directions | old α / mechanism / ambiguous | new α / mechanism / ambiguous |
|---|---|---|---|
| AAL–LG | AAL→LG, LG→AAL | −0.5 / predation / false | −0.5 / predation / **true** |
| AM–CC | AM→CC, CC→AM | −0.5 / predation / false | −0.5 / predation / **true** |
| AM–EL | AM→EL, EL→AM | −0.5 / predation / false | −0.5 / predation / **true** |
| AM–LG | AM→LG, LG→AM | −0.5 / predation / false | −0.5 / predation / **true** |
| CC–EL | CC→EL, EL→CC | −0.5 / predation / false | −0.5 / predation / **true** |
| CC–LG | CC→LG, LG→CC | −0.5 / predation / false | −0.5 / predation / **true** |
| EL–LG | EL→LG, LG→EL | −0.5 / predation / false | −0.5 / predation / **true** |
| EL–MS | EL→MS, MS→EL | −0.5 / predation / false | −0.5 / predation / **true** |
| EL–OM | EL→OM, OM→EL | −0.5 / predation / false | −0.5 / predation / **true** |

### 3.3 Counts

| quantity | pre-fix | corrected | change |
|---|---|---|---|
| directed rows | 227 | 227 | 0 |
| changed directed rows | — | — | **34** = 16 reverse + 18 ambiguous |
| `ambiguous = true` directed rows | 44 | **62** | +18 |
| `ambiguous = true` unordered pairs | 22 | **31** | +9 |

The plan’s “31 cells / 62 rows” refers to the corrected table’s **total** ambiguous set:
31 unordered pairs / 62 directed rows. It decomposes as 22 pairs already flagged
(including the 16 reverse-direction pairs) plus the 9 newly flagged mixed pairs.
The full changed-cell list is saved as `results/phase2/interaction_changed_cells.csv`.

---

## 4. Experiment-level effect

### 4.1 Definitions

Both `runs_index.csv` files were joined on `scenario,gcm` (45 rows each; the control row
`control/baseline` is included and is trivially unchanged). For each run `r` and metric `m`:

- `SC_pre(r) = pre_value(r) − pre_control(m)`
- `SC_post(r) = post_value(r) − post_control(m)`
- `delta_change(r) = SC_post(r) − SC_pre(r)`

“Sign changed” means `sign(SC_pre) ≠ sign(SC_post)` (a zero pre-delta therefore counts as a
change when the post-delta is non-zero). “Changed by more than 1 %” is reported twice:
relative to the pre-fix scenario-minus-control delta (`|delta_change| > 0.01·|SC_pre|`) and
relative to the pre-fix value itself.

### 4.2 Summary by metric (45 runs each)

| metric (runs_index column) | mean \|Δ(SC)\| | max \|Δ(SC)\| | max \|SC_pre\| | max \|SC_post\| | sign changes | >1 % of SC delta | >1 % of value |
|---|---|---|---|---|---|---|---|
| `basin_native_richness` | 0.000545 | 0.002581 | 0.011613 | 0.011613 | 3 | 16 | 0 |
| `basin_realised_richness_loss` | 0.000180 | 0.001290 | 0.005247 | 0.004925 | 3 | 15 | 3 |
| `basin_total_biomass` | 0.000364 | 0.001074 | 0.127763 | 0.127058 | 0 | 22 | 0 |
| `realised_warming_2026_2045_mean_degc` | **0.0** | **0.0** | 1.09069 | 1.09069 | 0 | 0 | 0 |
| `realised_warming_2036_2045_mean_degc` | **0.0** | **0.0** | 1.29214 | 1.29214 | 0 | 0 | 0 |
| `realised_warming_n_sites` | **0.0** | **0.0** | 0.0 | 0.0 | 0 | 0 | 0 |

The realised-warming columns are **bit-identical** before and after the fix, so the comparison
isolates the interaction matrix; the forcing was not re-derived.

The underlying scenario-minus-control magnitudes are likewise small and barely move:
`mean |SC_pre| → mean |SC_post|` is 0.00356 → 0.00376 for native richness,
0.00169 → 0.00169 for richness loss, and 0.03493 → 0.03490 for total biomass
(per-run values are in `results/phase2/climate_prefix_vs_phase2.csv`).

### 4.3 Sign changes

`basin_native_richness` (3):

| scenario / GCM | SC_pre | SC_post |
|---|---|---|
| ssp126 / IITM-ESM | −0.001290 | 0.0 |
| ssp245 / CMCC-CM2-SR5 | 0.0 | −0.001290 |
| ssp245 / IITM-ESM | 0.0 | −0.001290 |

`basin_realised_richness_loss` (3):

| scenario / GCM | SC_pre | SC_post |
|---|---|---|
| ssp126 / CMCC-CM2-SR5 | −0.000387 | +0.000903 |
| ssp245 / CMCC-CM2-SR5 | −0.000387 | +0.000903 |
| ssp245 / NorESM2-MM | −0.000323 | +0.000968 |

`basin_total_biomass`: **0 sign changes.**

All six flips involve deltas of at most 0.00129 in the metric’s own units (1/775 for the
site-mean richness) and are rounding-scale, not direction reversals of the warming response.

### 4.4 Basin-level aggregates not in `runs_index.csv`

The requested `st_final_biomass`, `st_relative_biomass_change` and `st_quasi_extinct_fraction`
columns do **not** exist. The corresponding available final-year columns in
`export/levels/level_basin.csv` are used instead (45 runs):

| level_basin.csv column | mean \|Δ\| | max \|Δ\| | sign changes | >1 % relative |
|---|---|---|---|---|
| `mean_total_biomass` | 0.001206 | 0.001951 | 0 | 0 |
| `mean_total_biomass_relative` | 2.049e8 | 1.078e9 | 0 | 1 (but see §6) |
| `mean_native_quasi_extinct_fraction` | 0.000187 | 0.000891 | 0 | 0 |

Saved as `results/phase2/basin_aggregate_prefix_vs_phase2.csv`.

### 4.5 Largest changes

Top 5 by absolute change in the scenario-minus-control delta across the three biological
basin metrics:

| rank | scenario / GCM | metric | pre | post | Δ | % of SC delta |
|---|---|---|---|---|---|---|
| 1 | ssp585 / MRI-ESM2-0 | `basin_native_richness` | 2.650323 | 2.647742 | **−0.002581** | −25.0 |
| 2 | ssp370 / EC-Earth3-Veg | `basin_native_richness` | 2.642581 | 2.643871 | +0.001290 | +50.0 |
| 3 | ssp245 / IITM-ESM | `basin_native_richness` | 2.640000 | 2.638710 | −0.001290 | −100.0 |
| 4 | ssp245 / CMCC-CM2-SR5 | `basin_native_richness` | 2.640000 | 2.638710 | −0.001290 | −100.0 |
| 5 | ssp370 / NorESM2-MM | `basin_native_richness` | 2.646452 | 2.645161 | −0.001290 | −50.0 |

Largest changes in the other metrics:

| metric | scenario / GCM | pre | post | Δ |
|---|---|---|---|---|
| `basin_realised_richness_loss` | ssp126 / CMCC-CM2-SR5 | 0.071806 | 0.073097 | +0.001290 |
| `basin_realised_richness_loss` | ssp245 / NorESM2-MM | 0.071871 | 0.073161 | +0.001290 |
| `basin_total_biomass` | ssp370 / MPI-ESM1-2-HR | 30.274394 | 30.274512 | −0.001074 |
| `basin_total_biomass` | ssp245 / UKESM1-0-LL | 30.366615 | 30.366865 | −0.000943 |

`Δ` = post − pre. 1/775 = 0.0012903; basin richness changes are therefore 1–2 site-presence
flips out of 775.

---

## 5. Plain-English interpretation

**The parser fix was not material to the weak-warming conclusion.** The corrected interaction
matrix changes 34 directed entries (16 reverse directions that wrongly inherited −0.5 predation
instead of −0.3 competition, plus 18 ambiguity flags on mixed cells) but moves the 2026–2045
climate-endpoint metrics by at most **0.0026 species** of basin mean
native richness (2 of 775 sites), **0.0013** of realised richness loss, and **0.0011** biomass
units on a ≈30-unit mean (≈0.004 %). The applied forcing is bit-identical, so this is a pure
interaction-matrix effect.

Under weak warming the scenario-minus-control deltas are themselves tiny (max |SC| ≈ 0.012
species, ≈0.005 richness loss, ≈0.13 biomass). The fix preserves the **sign** of the delta in
42/45 richness runs, 42/45 richness-loss runs and **45/45** biomass runs; the three richness
and three richness-loss flips are near-zero rounding-scale changes, not a reversal of the
warming response. No run’s final basin richness or total biomass changed by more than 1 % of
its own value; the three richness-loss runs that moved by >1 % of value moved by ≈0.0013
absolute (≈0.0718 → ≈0.0731).

The fix is therefore important for the **correctness and auditability of the α matrix** (and
for the `c = 1 − α` competition coefficients implied by it), but it does not rescue or
overturn the climate conclusion. The Phase 2 narrative can retain the weak-warming
disposition with this effect quantified; the material structural sensitivity remains the
competition-form switch (P2.3), not this parser repair.

---

## 6. Caveats

- The comparison is a **delta of deltas** (`SC_post − SC_pre`); the control row is included
  and is trivially zero.
- The “>1 % of SC delta” column is inflated by near-zero pre-fix deltas (e.g. a 0.00032
  species change is 100 % of a 0.00032 species delta). The absolute columns are the
  interpretable scale.
- `mean_total_biomass_relative` is numerically pathological (≈1e11, division by a near-zero
  baseline); its absolute deltas are not interpretable and only the relative change is
  meaningful (largest 1.008 %, `ssp126/EC-Earth3-Veg`). The sign never changes.
- The requested `st_*` aggregate columns do not exist under those names; `level_basin.csv`
  equivalents were used and are labelled.
- This task ran the climate ensemble only. The obstacle and interaction sweeps are out of
  scope here (P2.1b/P2.4) and were not run.
- `results/climate_scenarios_corrected`, `results/sensitivity_obstacles_corrected`,
  `results/alt_interactions_corrected` and
  `data/BIOTIC/interaction_matrix_long_prefix.csv` were read only.

---

## 7. Output files

| path | content |
|---|---|
| `results/climate_scenarios_phase2/runs_index.csv` | new 45-run index (46 lines) |
| `results/phase2_climate_run.log` | run log, ends “Climate scenarios complete” |
| `results/phase2/climate_prefix_vs_phase2.csv` | 270 rows = 45 runs × 6 metrics, long format, with pre/post value, control, deltas and change |
| `results/phase2/parser_fix_effect_summary.csv` | per-metric summary (13 columns) |
| `results/phase2/basin_aggregate_prefix_vs_phase2.csv` | 45 × 3 final-year basin aggregates |
| `results/phase2/interaction_changed_cells.csv` | 34 changed directed cells (16 reverse + 18 ambiguity flags) |

---

## 8. Reproduction

The comparison is fully specified by §4.1: read both `runs_index.csv` files with
`DataFrames`/`CSV`, inner-join on `scenario,gcm`, subtract each dataset’s `control/baseline`
row to form `SC_pre`/`SC_post`, and report the statistics in §4.2–§4.4. The interaction-table
diff is defined in §3 (join the two long tables on `source,target` and compare
`alpha,mechanism,ambiguous,raw_text`). The derived artefacts are the six files in §7; no
model code was changed for this task.
