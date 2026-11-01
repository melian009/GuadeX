# GuadeX Phase 2 — P2.4: Passability Floor and On-Path-Only Obstacles

**Scope:** run the 32-run obstacle/passability × upstream-cost sweep under (A) a
compounded-passability floor of `1e-3` and (B) obstacle snapping restricted to
on-path parent edges, then compare each variant with the corrected P2.1b baseline
(`results/sensitivity_obstacles_phase2`).

**Status:** companion evidence for Phase 2 item P2.4 (`§4` of
`.kilo/plans/1793526105792-second-review-triage-plan.md`). This document does not
revise the report or the status record; that is P2.6. It complements
`docs/GuadeX_Phase2_ObstacleSweep_Update.md` (P2.1b).

---

## 1. Provenance, commands and elapsed time

All work was done in `E:\projects\GuadeX` (Windows, PowerShell 7, Julia 1.12.6) on
2026-11-01. `git rev-parse --short HEAD` = **`c340398`**; the Phase 1 edits are
uncommitted working-tree changes, as in the handover. No source file was changed
for this task.

Baseline (read only, not overwritten):
`results/sensitivity_obstacles_phase2` — 32 runs, wall time **831.8 s** as recorded
in `docs/GuadeX_Phase2_ObstacleSweep_Update.md`.

Variant A — floor `1e-3`:

```powershell
$env:GUADEX_PARAMETERS_FILE='parameters_climate_scenarios_corrected.toml'
$env:GUADEX_OBSTACLE_FLOOR='0.001'
$env:GUADEX_OBSTACLES_PARENT_ONLY='0'
$env:GUADEX_SENSITIVITY_OUTPUT_DIR='results/sensitivity_obstacles_floor1e3'
julia --project=. run_sensitivity_report.jl *> results/phase2_obstacle_floor1e3_run.log
```

- Wall time: **789.3 s** (13 min 09 s), normal exit.
- Closing log line:
  `Obstacle sensitivity complete: 32 run(s) under 'results/sensitivity_obstacles_floor1e3' in 789.3s`.

Variant B — on-path parent edges only:

```powershell
$env:GUADEX_PARAMETERS_FILE='parameters_climate_scenarios_corrected.toml'
$env:GUADEX_OBSTACLE_FLOOR='0'
$env:GUADEX_OBSTACLES_PARENT_ONLY='1'
$env:GUADEX_SENSITIVITY_OUTPUT_DIR='results/sensitivity_obstacles_parentonly'
julia --project=. run_sensitivity_report.jl *> results/phase2_obstacle_parentonly_run.log
```

- Wall time: **806.2 s** (13 min 26 s), normal exit.
- Closing log line:
  `Obstacle sensitivity complete: 32 run(s) under 'results/sensitivity_obstacles_parentonly' in 806.2s`.

A `GUADEX_CONFIG_ONLY=1` smoke test was run for both variants before the sweeps and
resolved the expected grid: 2 climate models × 4 upstream costs × 4 passability
scenarios, `competition: :shared_k`, and the intended `floor` / `on_path_only`
values.

---

## 2. Run and digest inventory

| variant | output directory | runs | `passability_floor` | `obstacles_on_path_only` | `competition_structure` | representative digest | elapsed |
|---|---|---|---|---|---|---|---|
| baseline (P2.1b) | `results/sensitivity_obstacles_phase2` | 32 | 0.0 | false | `shared_k` | `580fba8364bccfe4` | 831.8 s |
| A floor `1e-3` | `results/sensitivity_obstacles_floor1e3` | 32 | 0.001 | false | `shared_k` | `f2d862ed4f4eeea2` | 789.3 s |
| B parent-only | `results/sensitivity_obstacles_parentonly` | 32 | 0.0 | true | `shared_k` | `39b9faa914d04dc7` | 806.2 s |

Verification:

| check | A (floor) | B (parent-only) |
|---|---|---|
| `runs_index.csv` | 33 lines = header + 32 runs | 33 lines = header + 32 runs |
| unique `(scenario,gcm,uc,passability)` keys | 32/32 | 32/32 |
| metadata files present | 32/32 | 32/32 |
| metadata `passability_floor` | 0.001 in 32/32 | 0.0 in 32/32 |
| metadata `obstacles_on_path_only` | false in 32/32 | true in 32/32 |
| metadata `competition_structure` | `shared_k` in 32/32 | `shared_k` in 32/32 |
| run digests differing from baseline | **32/32** | **32/32** |
| forcing columns differing from baseline | **0/32** for all five columns | **0/32** for all five columns |
| log final line | 32 runs | 32 runs |

The five forcing columns checked were `realised_warming_2026_2045_mean_degc`,
`realised_warming_2036_2045_mean_degc`, `realised_warming_n_sites`,
`warming_end_degc` and `basin_mean_temperature_c`; all were bit-identical to the
baseline. The sweeps write no `run_digest.txt` and no raw `code_version` key, as
in the baseline; identity is the `parameter_digest` (which folds the code version)
plus the `.sensitivity_complete.jld2` fingerprint. Each sweep has 16 distinct
digest values because the digest does not include the climate model; the two GCMs
that share `(upstream_cost, passability_scenario)` share a digest but have
different fingerprints (the fingerprint includes scenario/GCM), so no run is
reused.

---

## 3. Why these two variants

### 3.1 Floor value and justification

The obstacle overlay applies a uniform upstream passability of **0.1** per
matched structure, compounded multiplicatively along a link
(`src/data_preparation.jl:749-750`). A link with `n` matched structures therefore
has an upstream compounded passability of `0.1^n`:

| barriers `n` | compounded upstream factor |
|---|---|
| 1 | 1e-1 |
| 2 | 1e-2 |
| 3 | 1e-3 |
| 4 | 1e-4 |
| 33 | 1e-33 |

The chosen floor is **1e-3 = 0.1³**, i.e. the passability of three standard
upstream barriers in series. It is a defensible minimum passable link because:

- the uniform 0.1 coefficient is a structural assumption, not a calibrated
  per-barrier value. A product of four or more barriers (`1e-4` to `1e-33`) is
  dominated by the assumed coefficient and by floating-point scale, not by an
  observed hydraulic property. Below `1e-3` the exponent is effectively numerical
  noise for the dispersal kernel.
- the floor is still a 99.9 % reduction relative to an unobstructed link. It does
  not declare the link passable in any practical sense; it only stops a
  33-barrier product from collapsing to `1e-33` and makes the result insensitive
  to the arbitrary coefficient beyond the third barrier.
- the floor applies only to restricted links (`0 < p < floor`); unrestricted
  links (passability 1.0) are untouched (`src/ode_model.jl:410-417`).
- because `run_sensitivity_scenario` applies the passability-scenario factor to
  the restricted dam entries first and then calls `precompute_dispersal_matrix`
  (`sensitivity_core.jl:562-570`), the floor also limits links pushed below
  `1e-3` by the `reduced` (×0.5) and `blocked` (×0.1) scenarios.

### 3.2 Parent-only rationale

The corrected network has **508 parent edges and 266 reconstructed MST root-join
edges** (774 total; 13,871.5 km, of which 5,234.5 km = 37.7 % are root joins).
Root joins are a reconstruction device, not observed river reaches. In the
baseline overlay, **341 obstacles match a root-join edge and 245 are applied
there**. If the obstacle/passability conclusions depend on barriers assigned to
those reconstructed edges, they are fragile.

The parent-only switch (`obstacle_scope_edges`, `src/data_preparation.jl:2406-2419`)
restricts obstacle snapping to the 508 parent edges, so root joins carry no
obstacles. It is a topological robustness check, not a pure "delete 245 barriers"
experiment: an obstacle that no longer lies within 2 km of a parent edge becomes
unmatched, while one that does can re-snap to a nearby parent edge and stay
applied. §5 reports both effects.

---

## 4. Compounded-passability diagnostics (variant A)

Baseline overlay, recomputed from `data_base.obstacle_metadata` and the directed
dam matrix:

| quantity | baseline value |
|---|---|
| obstacles total / matched / applied / excluded | 1658 / 1036 / 796 / 387 |
| links carrying at least one applied obstacle | 330 |
| maximum obstacles on one link | **33** |
| minimum obstacle-only upstream factor (`0.1^n`) | **1e-33** |
| minimum combined dam factor (legacy dam × obstacle overlay) | **1e-34** |
| links with obstacle factor < 1e-6 | **16** |
| links with ≥ 5 obstacles | **28** |
| links with ≥ 4 obstacles (obstacle factor < 1e-3) | 49 |

### 4.1 The 33-obstacle link

The extreme link is **`1.33.2`–`1.33.58`**:

| property | value |
|---|---|
| obstacles on the link | 33 |
| root join? | **no — it is a parent edge** |
| obstacle-only upstream factor | **1e-33** |
| combined dam upstream / downstream factor | **1e-34** / 1.16e-11 |
| floor effect | both directions raised to **1e-3** |
| parent-only effect | **unchanged: 33 obstacles, factor still 1e-33** |

The combined upstream factor is one order below the obstacle-only product because
a legacy dam restriction also multiplies the same directed entry. The floor
therefore raises the link from `1e-34` to `1e-3`, a 31-order change. Parent-only
does **not** touch this link: it is on the observed parent topology, so restricting
snapping to parent edges leaves all 33 obstacles in place.

### 4.2 How many links the floor changes

For the baseline passability scenario, the floor changes **69 directed dam
entries**, covering **56 unordered links** (49 of them already have obstacle
factors below `1e-3`; the remaining 7 fall below the floor only after the legacy
dam restriction is multiplied in). The scenario factor is applied before the
floor, so the affected set grows under reduced/blocked passability:

| passability scenario | scenario factor | directed dam entries below 1e-3 | unordered links below 1e-3 | dispersal entries changed |
|---|---|---|---|---|
| baseline | 1.0 | 69 | 56 | 69 |
| improved | 1.5 | 67 | 56 | 67 |
| reduced | 0.5 | 126 | 112 | 126 |
| blocked | 0.1 | 140 | 112 | 140 |

### 4.3 Does the 33-link dominate the floor's input change?

For each scenario, the floor adds positive dispersal rate to the 69–140 changed
entries. The share attributable to the two directed entries of the 33-link is:

| passability scenario | total added dispersal rate | 33-link added rate | 33-link share |
|---|---|---|---|
| baseline | 2.315e-5 | 2.396e-7 | **1.03 %** |
| improved | 2.138e-5 | 2.396e-7 | 1.12 % |
| reduced | 3.910e-5 | 2.396e-7 | 0.61 % |
| blocked | 5.708e-5 | 2.396e-7 | **0.42 %** |

The 33-link has the largest per-entry ratio (1e31), but it accounts for only
about **1 %** of the floor's added dispersal capacity because the added rate is
also scaled by `D_daily / distance`. The floor changes 56 other links in the
baseline scenario and 112 in reduced/blocked.

---

## 5. Parent-only matching diagnostics (variant B)

| quantity | baseline | parent-only | change |
|---|---|---|---|
| obstacles matched | 1036 | **759** | −277 |
| obstacles applied | 796 | **599** | −197 |
| obstacles excluded (status/out-of-basin/duplicate) | 387 | 387 | 0 |
| applied on root-join edges | 245 | 0 (by construction) | −245 |
| root-edge applied obstacles that became unmatched | — | **197** | — |
| root-edge applied obstacles re-snapped to a parent edge | — | **48** | — |
| applied on parent edges | 551 | 599 | +48 |
| links carrying obstacles | 330 | **230** | −100 |
| max obstacles on one link | 33 | 33 | 0 |
| minimum obstacle factor | 1e-33 | 1e-33 | 0 |

Of the 245 obstacles that the baseline applied on reconstructed root-join edges,
**197 become unmatched** and **48 are re-snapped onto a parent edge and remain
applied**. The 111 root-join links that carried obstacles lose all 245; the 48
re-snapped obstacles land on 23 parent links (12 already had obstacles, 11 are
new), so total links with obstacles fall from 330 to 230.

The 33-obstacle link is not affected because it is a parent edge. The
second-most extreme link, **`1.33.41`–`1.33.71` (32 obstacles, 1e-32), is a root
join and loses all 32**. Other large removals include `1.24.1`–`1.24.3`
(12 obstacles) and `1.35.19`–`1.35.25` (10). Re-snapping can also concentrate
obstacles: `1.33.2`–`1.33.96` rises from 9 to 28 obstacles. The full per-link
table is `results/phase2/passability_link_diagnostics.csv`.

---

## 6. Paired deltas against the baseline

Method: inner join on `climate_scenario,gcm,upstream_cost,passability_scenario`
(32 matches); `delta = variant − baseline`. The sweep publishes absolute
final-year basin metrics, not scenario-minus-control deltas.

### 6.1 Floor `1e-3`

| metric | mean Δ | mean \|Δ\| | max \|Δ\| | mean \|Δ\| % | max \|Δ\| % | sign flips |
|---|---|---|---|---|---|---|
| `basin_native_richness` | +0.01331 | 0.01331 | 0.02968 | 0.523 | 1.115 | 0 |
| `basin_realised_richness_loss` | +0.000134 | 0.000376 | 0.001935 | 0.242 | 1.268 | 0 |
| `basin_total_biomass` | +0.02732 | 0.02732 | 0.04472 | 0.0912 | 0.1349 | 0 |
| `basin_native_biomass` | +0.005582 | 0.005582 | 0.012744 | 0.0219 | 0.0469 | 0 |
| `basin_invasive_biomass` | +0.02173 | 0.02173 | 0.03495 | 0.445 | 0.692 | 0 |
| `st_final_biomass` | +0.1733 | 0.1733 | 0.4980 | 0.0444 | 0.1165 | 0 |
| `st_relative_biomass_change` | +0.0002913 | 0.0002913 | 0.0008583 | 0.0984 | 0.3416 | 0 |
| `st_quasi_extinct_fraction` | −0.0004340 | 0.0004340 | 0.013889 | 0.223 | 7.143 | 0 |

The largest absolute changes occur in the **blocked** scenario: native richness
+0.02968 at `ssp126/IITM-ESM uc=0.01 blocked` (2.66581 → 2.69548), ST final
biomass +0.498 at `ssp126/IITM-ESM uc=0.05 blocked` (427.324 → 427.822), and
total biomass +0.0447 (33.142 → 33.187). The `st_quasi_extinct_fraction` maximum
is a one-site change (1/72) in `ssp126/IITM-ESM uc=0.05 reduced` (0.1944 →
0.1806). Mean native-richness delta rises from +0.0102 (baseline passability) to
+0.0198 (blocked), consistent with the floor changing more links after the
blocked scaling.

### 6.2 Parent-only

| metric | mean Δ | mean \|Δ\| | max \|Δ\| | mean \|Δ\| % | max \|Δ\| % | sign flips |
|---|---|---|---|---|---|---|
| `basin_native_richness` | +0.2135 | 0.2135 | 0.2968 | 8.457 | 10.519 | 0 |
| `basin_realised_richness_loss` | +0.000491 | 0.006955 | 0.01733 | 9.702 | 27.763 | 0 |
| `basin_total_biomass` | −0.7161 | 0.7161 | 1.1690 | 2.335 | 3.527 | 0 |
| `basin_native_biomass` | −0.3639 | 0.3639 | 0.6001 | 1.433 | 2.443 | 0 |
| `basin_invasive_biomass` | −0.3522 | 0.3528 | 0.6480 | 7.099 | 12.834 | 0 |
| `st_final_biomass` | −0.3432 | 5.677 | 15.682 | 1.557 | 4.660 | 0 |
| `st_relative_biomass_change` | −0.001593 | 0.009845 | 0.02833 | 2.825 | 6.904 | 0 |
| `st_quasi_extinct_fraction` | +0.01833 | 0.03964 | 0.11420 | 28.904 | 137.037 | 0 |

Parent-only is a materially larger perturbation than the floor. The largest
native-richness change is +0.2968 at `ssp126/IITM-ESM uc=0.01 baseline`
(2.9574 → 3.2542); the largest ST final-biomass change is −15.68 at
`ssp245/UKESM1-0-LL uc=0.5 improved` (336.52 → 320.84); the largest ST
quasi-extinction change is +0.1142 at `ssp126/IITM-ESM uc=0.05 blocked`
(0.0833 → 0.1975). The ST biomass response is scenario-dependent in sign
(negative in improved/blocked, positive in baseline/reduced), which is why the
mean signed delta is small relative to the mean absolute delta.

---

## 7. Factor ranking

The ranking replicates the `effect_spans()` / `median_span_std` method used in
`results/phase2/analyze_obstacle_effect.jl` (validated there against the pre-fix
`report_plots/parameter_effect_rank.csv`).

**The factor order is identical in the baseline and both variants for all eight
metrics: Spearman ρ = 1.00, top-two factors identical 8/8, no sign flips in any
metric.**

| metric | order (baseline = floor = parent-only) |
|---|---|
| `basin_native_richness` | upstream_cost > passability > climate_model |
| `basin_realised_richness_loss` | upstream_cost > passability > climate_model |
| `basin_total_biomass` | upstream_cost > passability > climate_model |
| `basin_native_biomass` | upstream_cost > passability > climate_model |
| `basin_invasive_biomass` | passability > upstream_cost > climate_model |
| `st_final_biomass` | passability > upstream_cost > climate_model |
| `st_relative_biomass_change` | passability > upstream_cost > climate_model |
| `st_quasi_extinct_fraction` | upstream_cost > passability > climate_model |

Score movements (`median_span_std`, baseline → floor → parent-only) for the
metrics that move most:

| metric | upstream_cost | passability | climate_model |
|---|---|---|---|
| `basin_native_richness` | 2.7574 → 2.7572 → 2.7544 | 0.4873 → 0.4477 → 0.4536 | 0.00802 → 0.00990 → 0.00859 |
| `basin_invasive_biomass` | 1.2175 → 1.2243 → 1.4440 | 2.1563 → 2.1509 → 2.2317 | 0.0407 → 0.0584 → 0.0653 |
| `st_final_biomass` | 0.4148 → 0.4145 → 0.9372 | 2.3502 → 2.3498 → 2.3428 | 0.2513 → 0.2511 → 0.2614 |
| `st_quasi_extinct_fraction` | 2.1721 → 2.1683 → 2.7452 | 1.3447 → 1.3423 → 0.7445 | 0.0000 → 0.0516 → 0.0931 |

The only notable score move is the parent-only ST quasi-extinction split:
upstream cost strengthens and passability weakens, but upstream cost remains
first and passability second. The baseline `climate_model` score of 0.0 is the
one-block coincidence already documented in the P2.1b update (two GCMs coincide
in one block; the median span over 16 blocks is zero).

---

## 8. The 33-obstacle link: what the evidence does and does not show

**Finding.** The `1e-33` link is `1.33.2`–`1.33.58`, a **parent edge with 33
obstacles**. It is not one of the reconstructed root joins. Therefore:

- **The floor acts on it directly.** The combined upstream dam factor is
  `1e-34`; the floor raises it (and its downstream direction) to `1e-3`.
- **Parent-only does not act on it.** Restricting snapping to parent edges leaves
  all 33 obstacles and the `1e-33` obstacle factor in place. Parent-only removes
  the *second*-most extreme link (`1.33.41`–`1.33.71`, 32 obstacles, root join)
  and the other 213 root-edge-applied barriers, but not the 33-link.

**Does it drive the ranking?** No, on the aggregate evidence:

1. The floor changes 56 links in the baseline scenario and 112 under
   reduced/blocked, not just the 33-link. The 33-link contributes only
   **0.42–1.12 %** of the total dispersal rate added by the floor.
2. The floor response is small (native richness +0.013 basin mean, maximum
   +0.030; total biomass maximum +0.045 on ≈33; ST final biomass maximum +0.498
   on ≈427) and leaves every factor order and every signed factor contrast intact.
3. The much larger parent-only response is driven by removing 245 obstacles from
   111 root-join links (including a 32-obstacle link), not by the 33-link, which
   parent-only leaves untouched.

**What the evidence does not show.** No single-link experiment was run that
removes or restores only `1.33.2`–`1.33.58`. The results therefore cannot rule
out a local, non-linear effect at the few sites adjacent to that link that cancels
in the basin mean, and parent-only cannot isolate the link because it keeps it.
The supported statement is the aggregate one: the 33-link is not necessary to
reproduce the factor ranking, and raising it by 31 orders of magnitude moves the
basin metrics by at most the small floor deltas above.

---

## 9. Plain-English conclusion

**The `1e-33` link does not drive the passability ranking, and the conclusions
survive both variants.**

- **Ranking and direction are stable.** For all eight metrics the factor order is
  identical in the baseline, the `1e-3` floor and the parent-only overlay
  (Spearman ρ = 1.00; top two factors unchanged 8/8; zero sign flips in 256
  metric-run pairs per variant). Upstream dispersal cost remains dominant for
  native richness, native biomass, total biomass, richness loss and ST
  quasi-extinction; passability remains dominant for invasive biomass and the ST
  biomass endpoints; climate model remains the smallest factor.
- **The floor is a small perturbation.** It changes at most 0.030 native species
  of basin richness (1.1 %), 0.045 biomass units on ≈33 (0.13 %), 0.50 ST biomass
  units on ≈427 (0.12 %), and one ST quasi-extinction site (1/72). It makes the
  model insensitive to the arbitrary 0.1-per-barrier coefficient beyond three
  barriers and removes the `1e-33`/`1e-34` numerical extreme.
- **Parent-only is a material structural sensitivity for magnitudes, not for
  ranking.** Removing obstacles from reconstructed root-join edges raises native
  richness by 8.5 % on average (up to 10.5 %), lowers total biomass by 2.3 %,
  lowers invasive biomass by 7.1 %, and changes ST quasi-extinction by up to
  0.114 absolute (137 % of a small baseline fraction). The direction of the
  factor effects and the factor order are unchanged, but the absolute magnitudes
  of native richness, invasive biomass and ST quasi-extinction carry a structural
  uncertainty from the root-join obstacle overlay.
- **The 33-link is not the cause of either result.** The floor changes it but its
  share of the added dispersal capacity is about 1 %; parent-only leaves it in
  place and still moves the results more, because it removes 245 root-edge
  barriers. The `1e-33` value is an extreme but not a load-bearing number.

---

## 10. Output files

| path | content |
|---|---|
| `results/sensitivity_obstacles_floor1e3/` | 32 complete floor runs (`runs_index.csv`, `simulation_output.jld2`, exports, completion markers) |
| `results/sensitivity_obstacles_parentonly/` | 32 complete parent-only runs |
| `results/phase2_obstacle_floor1e3_run.log` | floor run log (ends with 32 runs) |
| `results/phase2_obstacle_parentonly_run.log` | parent-only run log (ends with 32 runs) |
| `results/phase2/obstacle_floor_vs_baseline.csv` | 256 rows = 32 runs × 8 metrics, `baseline`, `variant_value`, `delta_abs`, `delta_rel` |
| `results/phase2/obstacle_parentonly_vs_baseline.csv` | same for parent-only |
| `results/phase2/passability_variant_summary.csv` | per-variant per-metric delta summary, factor order, Spearman ρ |
| `results/phase2/passability_link_diagnostics.csv` | per-link obstacle counts, root/parent flag, compounded factors, floor and parent-only effects |
| `results/phase2/passability_matching_diagnostics.csv` | matched/applied/root-join change counts |
| `results/phase2/passability_floor_scenario_counts.csv` | floor burden and 33-link share by passability scenario |
| `results/phase2/passability_variant_by_scenario.csv` | paired deltas by passability scenario |
| `results/phase2/passability_rank_detail.csv` | `median_span_std` scores for baseline and both variants |
| `results/phase2/passability_variant_verification.csv` | index/digest/metadata checks |
| `results/phase2/analyze_passability_variants.jl` | analysis script that produced the above |

---

## 11. Caveats

- The obstacle sweep publishes **absolute final-year basin metrics**, not
  scenario-minus-control deltas; the deltas here are direct variant − baseline
  differences, exactly as in the P2.1b update.
- Parent-only is **not** a clean "remove 245 barriers" experiment: 197 root-edge
  applied obstacles become unmatched but 48 re-snap to parent edges, raising some
  links (e.g. `1.33.2`–`1.33.96` from 9 to 28 obstacles). Its delta therefore
  combines removal and re-snapping.
- The floor is applied **after** the passability-scenario scaling, so its affected
  set grows from 56 unordered links (baseline factor 1.0) to 112 (reduced/blocked
  factors 0.5/0.1). This is why the floor deltas are largest in the blocked
  scenario.
- The 33-link isolation is **aggregate**, not a single-link counterfactual; see
  §8 for the exact limits of the claim.
- The sweep metadata carries no raw `code_version` string (writer limitation,
  unchanged from the baseline). Provenance is the `parameter_digest` (which folds
  `Guadex.code_version() = c340398`) and the completion-marker fingerprint.
- `results/sensitivity_obstacles_phase2`, `results/sensitivity_obstacles_corrected`,
  `results/climate_scenarios_corrected` and the other `*_corrected` directories
  were read only and were not overwritten. No source code was changed.
- This task ran the two P2.4 variants only. The climate ensemble and the
  interaction sweep were not re-run.
