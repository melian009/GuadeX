# GuadeX Phase 2 — P2.1b: Parser-Fix Effect on the Obstacle / Upstream-Cost Sweep

**Scope:** re-run the 32-run obstacle/passability × upstream-cost sweep against the regenerated
interaction table and quantify, not assert, what the interaction-parser fix changed in this
experiment.

**Status:** companion evidence for Phase 2 item P2.1 (`§4` of
`.kilo/plans/1793526105792-second-review-triage-plan.md`). This document does not revise the
report or the status record; that is P2.6. It complements
`docs/GuadeX_Phase2_ParserFix_Effect.md`, which covers the climate ensemble.

---

## 1. Provenance and commands

All work was done in `E:\projects\GuadeX` (Windows, PowerShell 7, Julia 1.12.6) on 2026-11-01.
`git rev-parse --short HEAD` = **`c340398`**; the parser fix and the other Phase 1 edits are
uncommitted working-tree changes, as in the handover.

Pre-fix audit baseline (read-only, untouched — last write 29 Oct 2026):

- `results/sensitivity_obstacles_corrected/runs_index.csv` — **33 lines** (header + 32 runs);
  32 unique `(climate_scenario, gcm, upstream_cost, passability_scenario)` keys.
- One run's metadata,
  `results/sensitivity_obstacles_corrected/ssp126__IITM-ESM/uc_0.01/baseline/export/run_metadata.json`,
  records **`"parameter_digest": "6b97f30399ade16c"`**. It does **not** contain a
  `code_version` key (the obstacle-sweep metadata writer has never emitted one; see §3.3).
- The 32 pre-fix runs each carry `simulation_output.jld2`, `export/run_metadata.json` and the
  completion marker `.sensitivity_complete.jld2`; they do **not** carry `run_digest.txt`.
  The pre-fix code version is documented as `9702169` in the Phase 2 handover and in
  `docs/GuadeX_Phase2_ParserFix_Effect.md` §1 (the climate metadata from the same correction
  series records `code_version = "9702169"`).

New run:

```powershell
$env:GUADEX_PARAMETERS_FILE='parameters_climate_scenarios_corrected.toml'
$env:GUADEX_SENSITIVITY_OUTPUT_DIR='results/sensitivity_obstacles_phase2'
julia --project=. run_sensitivity_report.jl *> results/phase2_obstacle_run.log
```

- Wall time: the script reports **`831.8 s` ≈ 13 min 52 s** (log `CreationTime` 13:56:36 →
  `LastWriteTime` 14:11:34 local, which includes Julia/package start-up). The configured
  3-year interim route (`projection_route = "interim_observed"`) makes this much faster than
  the September full-burn-in equivalent (~82 min quoted in the task brief).
- Exit status: the background process exited normally; the closing status block reports
  **“Obstacle sensitivity complete: 32 run(s) under 'results/sensitivity_obstacles_phase2' in 831.8s”**
  (followed only by the `Run index:` line and a closing rule);
  zero matches for `ERROR`, `LoadError`, `failed` or `Stacktrace` in the log.
- A `GUADEX_CONFIG_ONLY=1` smoke test before the run resolved the expected grid: 2 climate
  models (`ssp126/IITM-ESM` cool, `ssp245/UKESM1-0-LL` warm) × 4 upstream costs
  (`0.01, 0.05, 0.1, 0.5`) × 4 passability scenarios (`baseline`, `improved_passability`,
  `reduced_passability`, `blocked`), `competition: :shared_k`,
  `passability: floor=0.0, on_path_only=false`, output `results/sensitivity_obstacles_phase2`.

---

## 2. New-run verification

| check | result |
|---|---|
| `results/sensitivity_obstacles_phase2/runs_index.csv` | **33 lines = header + 32 runs** |
| Grid | 2 climate models × 4 upstream costs × 4 passability scenarios = **32 unique keys**, no duplicates |
| Required per-run artifacts | **32/32** have `simulation_output.jld2`, `export/run_metadata.json`, `export/levels/level_basin.csv`, `export/levels/species_timeseries.csv`, `export/levels/quasi_extinction_summary.csv` and `.sensitivity_complete.jld2` |
| `run_digest.txt` | **0/32** — the obstacle sweeps do not write this file; identity is the completion marker + metadata digest (see §3.3) |
| `competition_structure` (all 32) | `shared_k` |
| `passability_floor` (all 32) | `0.0` |
| `obstacles_on_path_only` (all 32) | `false` |
| `projection_route` / `interim_spin_up_years` | `interim_observed` / `3` |
| `parameter_digest` | present in 32/32 |
| `code_version` key in `run_metadata.json` | **absent** in 32/32 (writer limitation, not a mis-run); see §3.3 |
| `Guadex.code_version()` in this tree | **`c340398`** (the digest folds this value; `src/run_identity.jl:132`) |
| `run_digest.txt` vs metadata | n/a (no `run_digest.txt`) |

A representative metadata file is
`results/sensitivity_obstacles_phase2/ssp126__IITM-ESM/uc_0.01/baseline/export/run_metadata.json`.

### 2.1 Digest comparison against the pre-fix sweep

The run key `(climate_scenario, gcm, upstream_cost, passability_scenario)` was joined across the
two sweeps and the stored `parameter_digest` compared (32 runs in common):

| quantity | result |
|---|---|
| runs with a changed digest | **32/32** |
| distinct pre-fix digests | 32 |
| distinct post-fix digests | 32 |
| example | `ssp126/IITM-ESM uc=0.01 baseline`: `6b97f30399ade16c` → **`580fba8364bccfe4`** |

Full list: `results/phase2/obstacle_digest_comparison.csv`. The digests differ because the
interaction matrix is a digest input (`src/run_identity.jl:112`) and the code version changed
from `9702169` to `c340398`.

### 2.2 Forcing unchanged

To show that the comparison isolates the interaction matrix, the climate/forcing columns of the
two `runs_index.csv` files were compared run by run:

| column | rows differing (of 32) |
|---|---|
| `realised_warming_2026_2045_mean_degc` | **0** |
| `realised_warming_2036_2045_mean_degc` | **0** |
| `realised_warming_n_sites` | **0** |
| `warming_end_degc` | **0** |
| `basin_mean_temperature_c` | **0** |

The applied forcing and the basin mean temperature are bit-identical, exactly as in the climate
re-run (§2 of `GuadeX_Phase2_ParserFix_Effect.md`).

---

## 3. Effect on the sweep metrics

### 3.1 Method

The two `runs_index.csv` files were inner-joined on
`climate_scenario,gcm,upstream_cost,passability_scenario` (32 matches). The obstacle sweep
reports **absolute final-year basin metrics**, not scenario-minus-control deltas, so the
comparison is a direct value difference:

- `delta_abs = post − pre`
- `delta_rel = delta_abs / |pre|` (undefined for a zero baseline; none occur in these metrics)

Metrics analysed: `basin_native_richness`, `basin_realised_richness_loss`,
`basin_total_biomass`, `basin_native_biomass`, `basin_invasive_biomass`, and the three `st_*`
columns present in the obstacle index (`st_final_biomass`, `st_relative_biomass_change`,
`st_quasi_extinct_fraction`).

### 3.2 Summary by metric (32 runs each)

| metric | mean Δ (signed) | mean \|Δ\| | max \|Δ\| | mean \|Δ\| % | max \|Δ\| % | runs >1 % of value | sign changes |
|---|---|---|---|---|---|---|---|
| `basin_native_richness` | −0.000282 | 0.000524 | 0.001290 | 0.0208 | 0.0637 | 0 | 0 |
| `basin_realised_richness_loss` | +0.000057 | 0.000138 | 0.001290 | 0.1705 | 1.5819 | 2 | 0 |
| `basin_total_biomass` | −0.000599 | 0.001824 | 0.006304 | 0.00596 | 0.01902 | 0 | 0 |
| `basin_native_biomass` | −0.003102 | 0.003102 | 0.004592 | 0.01243 | 0.01817 | 0 | 0 |
| `basin_invasive_biomass` | +0.002507 | 0.002945 | 0.005023 | 0.06035 | 0.10658 | 0 | 0 |
| `st_final_biomass` | −0.000007 | 0.012920 | 0.056972 | 0.00339 | 0.01605 | 0 | 0 |
| `st_relative_biomass_change` | +0.000037 | 0.000040 | 0.000135 | 0.01191 | 0.03582 | 0 | 0 |
| `st_quasi_extinct_fraction` | +0.000434 | 0.000434 | 0.013889 | 0.2404 | 7.6923 | 1 | 0 |

“Sign changes” counts runs where `sign(pre) ≠ sign(post)` with `delta ≠ 0`; the count is **0 for
every metric**. The mean signed deltas are all ≤ 0.0031 in their own units and alternate in
sign across metrics, i.e. rounding-scale drift rather than a directional shift.

`basin_realised_richness_loss` has two runs moving >1 % of their own value, but both by exactly
**0.0012903** absolute (one site presence, 1/775): `ssp245/UKESM1-0-LL uc=0.1 baseline`
0.103914 → 0.105204 (+1.24 %) and `ssp245/UKESM1-0-LL uc=0.1 blocked` 0.081570 → 0.080280
(−1.58 %).

### 3.3 Identity-artifact note

Neither the pre-fix nor the new obstacle run metadata contains a `code_version` key. The sweep
writer (`sensitivity_run_metadata`, `sensitivity_core.jl:756-814`) emits `parameter_digest`,
`competition_structure`, `passability_floor` and `obstacles_on_path_only` but not the raw code
version string. The code version still enters every run's identity through
`parameter_digest(params; code=code_version())` (`src/run_identity.jl:106-133`) and is written
into the `.sensitivity_complete.jld2` fingerprint as part of that digest. In this tree
`Guadex.code_version()` returns **`c340398`**, and 32/32 digests differ from the pre-fix sweep
(§2.1), so the new runs are provably not the pre-fix ones.

---

## 4. Factor ranking before vs after

### 4.1 Method and validation

The ranking replicates exactly the obstacle-analysis part of
`scripts/plot_sensitivity_effects.jl`:

- For each metric and factor, `effect_spans()` holds the other two design columns fixed
  (blocks), computes the span `max − min` of the metric across the factor's levels in each
  block, and standardises it by the standard deviation of the metric over the run set
  (`span_std`).
- The factor score is `median_span_std` over its blocks; factors are ranked descending.
- Blocks: `passability_scenario` over `model × upstream_cost` (8 blocks);
  `upstream_cost` over `model × passability_scenario` (8 blocks); `climate_model` over
  `upstream_cost × passability_scenario` (16 blocks).

The replication is validated against the existing pre-fix file: for
`basin_native_richness`, the recomputed pre-fix values are **identical** to
`results/sensitivity_obstacles_corrected/report_plots/parameter_effect_rank.csv`
(`upstream_cost` 2.7584184769, `passability` 0.4891381602, `climate_model` 0.0080186584).

### 4.2 Result

Rank order is **identical before and after the parser fix for all eight metrics**:
Spearman ρ = **1.00** (n = 3 factors; 6·Σd² = 0), and the top two factors overlap in 8/8 metrics.

| metric | pre-fix order | post-fix order | rank changes | ρ |
|---|---|---|---|---|
| `basin_native_richness` | upstream_cost > passability > climate_model | same | 0 | 1.00 |
| `basin_realised_richness_loss` | upstream_cost > passability > climate_model | same | 0 | 1.00 |
| `basin_total_biomass` | upstream_cost > passability > climate_model | same | 0 | 1.00 |
| `basin_native_biomass` | upstream_cost > passability > climate_model | same | 0 | 1.00 |
| `basin_invasive_biomass` | passability > upstream_cost > climate_model | same | 0 | 1.00 |
| `st_final_biomass` | passability > upstream_cost > climate_model | same | 0 | 1.00 |
| `st_relative_biomass_change` | passability > upstream_cost > climate_model | same | 0 | 1.00 |
| `st_quasi_extinct_fraction` | upstream_cost > passability > climate_model | same | 0 | 1.00 |

Score magnitudes move by < 1 % of their own value in almost every case
(`median_span_std`, pre → post):

| metric | upstream_cost | passability | climate_model |
|---|---|---|---|
| `basin_native_richness` | 2.7584 → 2.7574 | 0.4891 → 0.4873 | 0.008019 → 0.008021 |
| `basin_realised_richness_loss` | 2.6572 → 2.6576 | 0.3682 → 0.3748 | 0.0063 → 0.0077 |
| `basin_total_biomass` | 2.1372 → 2.1367 | 1.4716 → 1.4719 | 0.0418 → 0.0419 |
| `basin_native_biomass` | 2.2951 → 2.2956 | 1.3301 → 1.3295 | 0.0444 → 0.0443 |
| `basin_invasive_biomass` | 1.2219 → 1.2175 | 2.1487 → 2.1563 | 0.0427 → 0.0407 |
| `st_final_biomass` | 0.4146 → 0.4148 | 2.3507 → 2.3502 | 0.2522 → 0.2513 |
| `st_relative_biomass_change` | 0.4146 → 0.4148 | 2.3507 → 2.3502 | 0.2522 → 0.2513 |
| `st_quasi_extinct_fraction` | 2.1683 → 2.1721 | 1.3423 → 1.3447 | **0.0516 → 0.0000** |

The only non-trivial score move is `st_quasi_extinct_fraction` × `climate_model`, and it is
explained by the single changed run below: in the pre-fix data the cool/warm difference at
`uc=0.05, reduced_passability` was 1/72 of established sites; after the fix the two climate
models coincide there, so the **median** climate span over the 16 blocks becomes 0 while the
**max** climate span is unchanged (0.2065 → 0.2069). The factor remains third and no ordering
changes.

### 4.3 Signed factor contrasts

To test conclusions rather than just ordering, the signed marginal contrasts were computed per
matched block and averaged:

- passability: `improved − baseline`, `reduced − baseline`, `blocked − baseline`
  (blocks = model × upstream_cost, 8 each)
- upstream cost: `0.5 − 0.01` (blocks = model × passability, 8)
- climate: `warm − cool` (blocks = upstream_cost × passability, 16)

**40/40 contrasts keep their sign** (0 flips) and the magnitudes move marginally. For the
primary `basin_native_richness`:

| contrast | pre | post |
|---|---|---|
| improved_passability − baseline | +0.028548 | +0.028548 |
| reduced_passability − baseline | −0.043065 | −0.043065 |
| blocked − baseline | −0.140645 | −0.140484 |
| upstream_cost 0.5 − 0.01 | −0.847258 | −0.847097 |
| warm − cool | −0.002258 | −0.002500 |

The same qualitative pattern holds for the other metrics: higher upstream cost and reduced/blocked
passability lower native richness and biomass; passability dominates invasive biomass and trout
(ST) metrics; climate is the smallest factor. The fix does not reverse or weaken any of these
directions.

---

## 5. Largest changes

Per-metric largest absolute change (post − pre):

| metric | run | pre | post | Δ | relative |
|---|---|---|---|---|---|
| `basin_native_richness` | ssp126/IITM-ESM, uc=0.05, baseline | 2.638710 | 2.640000 | **+0.001290** | +0.049 % |
| `basin_native_richness` | ssp245/UKESM1-0-LL, uc=0.05, reduced | 2.588387 | 2.587097 | −0.001290 | −0.050 % |
| `basin_realised_richness_loss` | ssp245/UKESM1-0-LL, uc=0.1, baseline | 0.103914 | 0.105204 | **+0.001290** | +1.24 % |
| `basin_realised_richness_loss` | ssp245/UKESM1-0-LL, uc=0.1, blocked | 0.081570 | 0.080280 | −0.001290 | −1.58 % |
| `basin_total_biomass` | ssp245/UKESM1-0-LL, uc=0.01, blocked | 33.148673 | 33.142369 | **−0.006304** | −0.019 % |
| `basin_native_biomass` | ssp126/IITM-ESM, uc=0.01, blocked | 27.962560 | 27.957969 | **−0.004592** | −0.016 % |
| `basin_invasive_biomass` | ssp245/UKESM1-0-LL, uc=0.5, baseline | 4.713242 | 4.718265 | **+0.005023** | +0.107 % |
| `st_final_biomass` | ssp245/UKESM1-0-LL, uc=0.5, reduced | 355.023729 | 355.080701 | **+0.056972** | +0.016 % |
| `st_relative_biomass_change` | ssp245/UKESM1-0-LL, uc=0.5, reduced | −0.378001 | −0.377865 | **+0.000135** | +0.036 % |
| `st_quasi_extinct_fraction` | ssp126/IITM-ESM, uc=0.05, reduced | 0.180556 | 0.194444 | **+0.013889** | +7.69 % |

Scales: `0.0012903 = 1/775` (one basin site-presence flip on a 775-site mean);
`0.013889 = 1/72` (one site on the 72 baseline-established ST sites). Several runs are tied at
exactly ±1/775 for native richness and richness loss, so the “largest” row is one of a tied set.
The largest biomass change is 0.006 unit on a ≈33-unit basin mean (0.019 %); the largest trout
endpoint change is 0.057 biomass units on ≈355 (0.016 %). Every other run moves less than these
values.

---

## 6. Plain-English interpretation

**The parser fix did not change the obstacle/passability conclusions.** The corrected interaction
matrix changes 34 directed entries (16 reverse directions that wrongly inherited −0.5 predation
instead of −0.3 competition, plus 18 ambiguity flags on mixed cells), but across the 32-run
obstacle × upstream-cost sweep it moves the final-year metrics by at most:

- **0.0013 native species** of basin mean richness (1 site of 775; ≤ 0.064 % of value),
- **0.0013** of realised richness loss,
- **0.0063 biomass units** on a ≈33-unit basin mean (0.019 %),
- **0.057** ST final biomass on ≈355 (0.016 %),
- **0.0139** of the ST quasi-extinct fraction (1 site of 72; 0.18 → 0.19).

The applied forcing and basin mean temperature are **bit-identical** before and after, so this is
a pure interaction-matrix effect. The factor ranking is **identical for all eight metrics**
(Spearman ρ = 1.00) and **all 40 signed factor-level contrasts keep their sign**. The qualitative
story of this experiment therefore stands: upstream dispersal cost is the dominant factor for
native richness, native biomass, total biomass, richness loss and ST quasi-extinction; passability
is the dominant factor for invasive biomass and the trout (ST) endpoints; climate model is the
smallest factor. The single visible score move
(`st_quasi_extinct_fraction` × `climate_model`) is a one-site change that makes the two climate
models coincide in one block; it does not reorder factors or reverse a direction.

The parser fix is material for the **correctness and auditability** of the α matrix (and the
implied `c = 1 − α` competition coefficients), not for the obstacle/passability conclusions. As in
the climate ensemble, the structural sensitivity that matters is the competition-form switch
(P2.3) and the passability-floor/on-path-only variants (P2.4), not this repair.

---

## 7. Caveats

- The obstacle sweep publishes **absolute final-year basin metrics**, not scenario-minus-control
  deltas. The deltas here are direct post − pre differences, unlike the climate comparison in
  `docs/GuadeX_Phase2_ParserFix_Effect.md` §4.1.
- Neither the pre-fix nor the new obstacle `run_metadata.json` records a raw `code_version`
  string; provenance is the `parameter_digest` (which folds the code version) plus the
  `.sensitivity_complete.jld2` fingerprint. `run_digest.txt` is not written by the sweeps.
- The pre-fix obstacle metadata also omits `competition_structure`, `passability_floor` and
  `obstacles_on_path_only` (they were added to the writer in the Phase 1 edits). The new metadata
  records all three explicitly.
- Relative deltas near small baselines are inflated: 2 of 32 `basin_realised_richness_loss` runs
  exceed 1 % of value, and 1 of 32 `st_quasi_extinct_fraction` runs exceeds 1 % of value (the
  1/72 site change). The absolute columns are the interpretable scale.
- The ranking was recomputed with the exact `effect_spans`/`median_span_std` definition from
  `scripts/plot_sensitivity_effects.jl` and validated to reproduce the pre-fix
  `parameter_effect_rank.csv` obstacle rows exactly. The sweep itself does not emit
  `report_plots/`; no figures were regenerated for this task.
- `results/sensitivity_obstacles_corrected`, `results/climate_scenarios_corrected`,
  `results/alt_interactions_corrected` and `data/BIOTIC/interaction_matrix_long_prefix.csv`
  were read only. The pre-fix obstacle index and metadata retain their 29 Oct timestamps.
- This task ran the obstacle sweep only. The climate re-run is P2.1a
  (`docs/GuadeX_Phase2_ParserFix_Effect.md`); the interaction sweep and the P2.4 floor /
  on-path-only variants are out of scope here.

---

## 8. Output files

| path | content |
|---|---|
| `results/sensitivity_obstacles_phase2/runs_index.csv` | new 32-run index (33 lines) |
| `results/sensitivity_obstacles_phase2/<model>/uc_<cost>/<pass>/` | 32 complete runs: `simulation_output.jld2`, `export/run_metadata.json`, full exports, `.sensitivity_complete.jld2` |
| `results/phase2_obstacle_run.log` | run log, ends “Obstacle sensitivity complete: 32 run(s)” |
| `results/phase2/obstacle_prefix_vs_phase2.csv` | 256 rows = 32 runs × 8 metrics; columns `climate_scenario,gcm,upstream_cost,passability_scenario,metric,pre,post,delta_abs,delta_rel` |
| `results/phase2/obstacle_effect_summary.csv` | per-metric delta summary (12 columns) |
| `results/phase2/obstacle_rank_prefix_vs_phase2.csv` | per-metric factor scores (`median_span_std`, etc.) for pre-fix and phase2 |
| `results/phase2/obstacle_rank_comparison.csv` | factor rank pre vs post per metric |
| `results/phase2/obstacle_factor_contrasts.csv` | 40 signed factor contrasts, pre vs post, with sign-flip flag |
| `results/phase2/obstacle_digest_comparison.csv` | 32 run-dir digest pairs, `differs` flag |
| `results/phase2/analyze_obstacle_effect.jl` | the analysis script that produced the above (build artefact, gitignored) |

---

## 9. Reproduction

1. Pre-fix baseline: read `results/sensitivity_obstacles_corrected/runs_index.csv`
   (33 lines) and the run metadata (e.g. `parameter_digest = 6b97f30399ade16c`).
2. New run: set `GUADEX_PARAMETERS_FILE=parameters_climate_scenarios_corrected.toml` and
   `GUADEX_SENSITIVITY_OUTPUT_DIR=results/sensitivity_obstacles_phase2`, then
   `julia --project=. run_sensitivity_report.jl`. The configured grid is
   2 models × 4 upstream costs × 4 passability scenarios. Success is the log line
   “Obstacle sensitivity complete: 32 run(s)”.
3. Deltas: inner-join the two indices on
   `climate_scenario,gcm,upstream_cost,passability_scenario`; compute `post − pre` and
   `(post − pre)/|pre|` for the eight metrics. `results/phase2/analyze_obstacle_effect.jl`
   performs steps 2–4 and writes the five CSVs in §8.
4. Ranking: apply `effect_spans()` / `median_span_std` from
   `scripts/plot_sensitivity_effects.jl` separately to each index (blocks as in §4.1) and sort
   factors descending. The script reproduces the pre-fix
   `report_plots/parameter_effect_rank.csv` obstacle rows exactly, which validates the
   replication.

No model code was changed for this task.
