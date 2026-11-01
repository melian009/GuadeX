# GuadeX Phase 2 — P2.3: Competition-Structure Sensitivity

**Scope:** run the full 44-scenario climate ensemble (plus the shared no-warming control)
under the reviewer's species-specific competition structure (`c = −α`) and quantify how much
the sign and magnitude of the scenario-minus-control climate response depend on the choice,
relative to the primary shared-K structure (`c = 1 − α`).

**Status:** evidence for Phase 2 item P2.3 (`§4` of
`.kilo/plans/1793526105792-second-review-triage-plan.md`) and decision 1 of
`docs/GuadeX_Phase2_Decisions.md` (keep `c = 1 − α` primary; run `c = −α` as a sensitivity).
This document does not revise the report or the status record; that is P2.6.

**Non-negotiable framing:** the primary specification remains the shared-K structure
(`c = 1 − α`). The species-specific structure is a sensitivity that shows the dependence; it
is **not** a replacement and does not change any primary result directory.

---

## 1. Provenance and commands

All work was done in `E:\projects\GuadeX` (Windows, PowerShell 7, Julia 1.12) on 2026-11-01.
`git rev-parse --short HEAD` = **`c340398`**; the Phase 1 edits remain uncommitted
working-tree changes, as in the handover. No source file was changed for this task.

Config-only smoke test (resolves the switch before the long run):

```powershell
$env:GUADEX_PARAMETERS_FILE='parameters_climate_scenarios_corrected.toml'
$env:GUADEX_COMPETITION_STRUCTURE='species_specific'
$env:GUADEX_CLIMATE_OUTPUT_DIR='results/climate_scenarios_competition_alt'
$env:GUADEX_CLIMATE_PLOT='0'
$env:GUADEX_CONFIG_ONLY='1'
julia --project=. run_climate_scenarios.jl
```

It printed `projection route: :interim_observed (interim_spin_up_years=3; control_run=true)`
and **`competition structure: :species_specific`**.

The ensemble run:

```powershell
$env:GUADEX_PARAMETERS_FILE='parameters_climate_scenarios_corrected.toml'
$env:GUADEX_COMPETITION_STRUCTURE='species_specific'
$env:GUADEX_CLIMATE_OUTPUT_DIR='results/climate_scenarios_competition_alt'
$env:GUADEX_CLIMATE_PLOT='0'
julia --project=. run_climate_scenarios.jl *> results/phase2_competition_alt_run.log
```

- Wall time: **30 min 08 s** (log `CreationTime` 14:50:18 → `LastWriteTime` 15:20:26 local).
- Exit status: exited normally; the closing line is
  **“Climate scenarios complete: 45 run(s) written under 'results/climate_scenarios_competition_alt'”**.
  No `ERROR`, `LoadError`, `failed` or `Stacktrace` match in the log. The single
  “parsed expected 10 columns …” CSV-widening warning is also present in
  `results/phase2_climate_run.log`, so it is a pre-existing data-file quirk, not new here.
- The main banner records **`competition:   :species_specific; passability floor: 0.0;
  on-path-only obstacles: false`**; everything else (config, forcing, 3-year interim route,
  six biological options, on-path network) is the corrected primary configuration.

The delta analysis is `results/phase2/analyze_competition_structure.jl` (build artefact,
gitignored), run as `julia --project=. results/phase2/analyze_competition_structure.jl`;
its log is `results/phase2_competition_analysis.log`.

---

## 2. Run inventory and identity evidence

| check | result |
|---|---|
| `results/climate_scenarios_competition_alt/runs_index.csv` | **46 lines = header + 45 runs** |
| Scenario coverage | `ssp126` 11, `ssp245` 11, `ssp370` 11, `ssp585` 11, `control` 1 |
| GCMs | 11 GCMs + the `baseline` control |
| `simulation_output.jld2` / `export/run_metadata.json` / `run_digest.txt` | **45/45** each |
| `resolved_config.competition_structure` | **`species_specific` in 45/45** |
| `projection_route` / `interim_spin_up_years` | `interim_observed` / `3` in 45/45 |
| `code_version` | **`c340398`** in 45/45 (= current `HEAD`) |
| `run_digest.txt` vs metadata `parameter_digest` | identical in 45/45 |
| Digest differs from the corresponding primary run | **45/45** |
| Realised forcing vs primary | **bit-identical**: max abs difference 0.0 for `warming_end_degc`, `realised_warming_2026_2045_mean_degc`, `realised_warming_2036_2045_mean_degc` |

The digest differs only because the competition structure is folded into the run identity
(`resolved_config.competition_structure`, which enters `parameter_digest`); the forcing and
every other option are unchanged. Examples:

| run | primary digest | alternative digest |
|---|---|---|
| `ssp126/ACCESS-CM2` | `96a0c5621104ebb9` | `3c48bfabfa8b2ed5` |
| `ssp245/UKESM1-0-LL` | `dd0490977b04dc79` | `b0a50abb7f3153bd` |
| `control/baseline` | `96ceb1dbce8cbbe6` | `c9fe913c8a32b3ca` |

As in the primary re-run report, the writer emits the structure at
`resolved_config.competition_structure`, not as a top-level `run_metadata.json` key.

---

## 3. The two formulations (LaTeX-ready plain text)

Let `i` index sites, `s,j` species, `K_i` the per-site carrying capacity, `α_sj` the
interaction coefficient of `j` on `s`, `r_eff` the environment-filtered intrinsic growth rate,
and `heat = k · max(0, T_i(t) − T_s^max)^2` the heat-stress term.

**Primary — shared per-site carrying capacity (`:shared_k`, `c = 1 − α`):**

```latex
\frac{dN_{is}}{dt} = N_{is}\left[ r_{\mathrm{eff}}\left(1 - \sum_j \frac{N_{ij}}{K_i}
  + \sum_j \alpha_{sj}\frac{N_{ij}}{K_i}\right) - heat \right],
\qquad c_{sj} = 1 - \alpha_{sj}
```

Because the empirical matrix has a zero diagonal, `α_ss = 0` and `c_ss = 1`. The bracket is
equivalently `1 − Σ_j c_sj N_ij / K_i`: the same `K_i` caps the whole local community, and
heterospecific competition is `c_sj = 1 − α_sj` with `α_sj < 0`.

**Alternative — species-specific self-limitation, no shared capacity (`:species_specific`, `c = −α`):**

```latex
\frac{dN_{is}}{dt} = N_{is}\left[ r_{\mathrm{eff}}\left(1 - \frac{N_{is}}{K_i}
  + \sum_{j\neq s} \alpha_{sj}\frac{N_{ij}}{K_i}\right) - heat \right],
\qquad c_{ss} = 1,\quad c_{sj} = -\alpha_{sj}\ (j \neq s)
```

The bracket is equivalently `1 − c_ss N_is/K_i − Σ_{j≠s} c_sj N_ij/K_i`. The species' own
density is the only self-limitation term; heterospecific interactions enter with `c_sj = −α_sj`.

**Coefficient ranges in the corrected table.** Negative α takes three values by mechanism:
`−0.3` (competition, 90 directed rows), `−0.5` (predation, 73), `−0.8` (displacement, 64).
Therefore:

- shared-K: `c = 1 − α ∈ {1.3, 1.5, 1.8}`, i.e. **c ∈ [1.3, 1.8]**;
- species-specific: `c = −α ∈ {0.3, 0.5, 0.8}`, i.e. **c ∈ [0.3, 0.8]**; the competition and
  predation rows alone give **[0.3, 0.5]**, and the `−0.8` displacement rows extend it to 0.8.

The structural consequence: under shared-K every heterospecific coefficient exceeds
self-limitation (`c_sj ≥ 1.3 > c_ss = 1`), so priority effects and exclusion are possible; under
the alternative every heterospecific coefficient is below self-limitation
(`c_sj ≤ 0.8 < c_ss = 1`), which favours stable coexistence. The two structures therefore differ
qualitatively, not just in a scaling factor.

---

## 4. How the comparison is defined

`results/phase2/competition_structure_deltas.csv` (88 rows = 44 scenarios × 2 structures) and
`results/phase2/competition_endpoint_summary.csv` (24 rows = 6 metrics × 2 endpoints × 2
structures) use:

- `δ_shared(r) = value_shared(r) − control_shared` and
  `δ_alt(r) = value_alt(r) − control_alt`, each from **its own `runs_index.csv` control row**
  (`control/baseline`), as required for the interim route.
- `Δδ(r) = δ_alt(r) − δ_shared(r)`; the reported “Δ” statistics are over the 44 matched runs.
- Metrics: the three index columns `basin_native_richness`, `basin_realised_richness_loss`,
  `basin_total_biomass`, plus the trout (ST) endpoints. `runs_index.csv` has no `st_*` columns
  and `level_basin.csv` is basin-wide only, so the ST metrics are read from each run's
  `export/levels/quasi_extinction_summary.csv` row `species == "ST"`:
  `final_biomass` → `st_final_biomass`, `relative_biomass_change` →
  `st_relative_biomass_change`, `quasi_extinct_fraction` → `st_quasi_extinct_fraction`.
- A “sign flip” is `sign(δ_alt) ≠ sign(δ_shared)` (a zero delta counts as a sign); the strictly
  opposite non-zero count is reported separately.

Control (no-warming, final-year) values differ between the structures, because the structure
changes the community state itself:

| metric | shared-K control | species-specific control |
|---|---|---|
| `basin_native_richness` | 2.640000 | 3.251613 |
| `basin_realised_richness_loss` | 0.0721935 | 0.0273763 |
| `basin_total_biomass` | 30.278003 | 63.566854 |
| `st_final_biomass` | 357.987906 | 758.864616 |
| `st_relative_biomass_change` | −0.372772 | +0.257020 |
| `st_quasi_extinct_fraction` | 0.291667 | 0.184211 |

Only scenario-minus-control deltas are compared across structures; the level differences above
are not interpreted as a model improvement or deterioration.

---

## 5. Delta dependence across all 44 scenarios

| metric | mean \|Δδ\| | max \|Δδ\| | mean signed Δδ | flips (`sign()`) | strictly opposite non-zero | mean \|δ_shared\| | mean \|δ_alt\| | ratio |
|---|---|---|---|---|---|---|---|---|
| `basin_native_richness` | 0.009150 | 0.021935 | −0.009150 | 30 | 26 | 0.003842 | 0.007126 | 1.86× |
| `basin_realised_richness_loss` | 0.001705 | 0.004925 | +0.001514 | 37 | 6 | 0.001733 | 0.0000748 | 0.04× |
| `basin_total_biomass` | 0.139127 | 0.457481 | +0.044738 | 3 | 3 | 0.035696 | 0.172867 | 4.84× |
| `st_final_biomass` | 52.262 | 77.457 | −52.262 | 0 | 0 | 7.171 | 59.434 | 8.29× |
| `st_relative_biomass_change` | 0.085884 | 0.127520 | −0.085884 | 0 | 0 | 0.012565 | 0.098449 | 7.83× |
| `st_quasi_extinct_fraction` | 0.011148 | 0.027047 | +0.011148 | 32 | 3 | 0.002210 | 0.009569 | 4.33× |

Per-structure delta distributions (positive / zero / negative runs, and mean):

| metric | shared-K | species-specific |
|---|---|---|
| `basin_native_richness` | 26 / 4 / 14, mean **+0.002023** | **0 / 0 / 44**, mean **−0.007126** |
| `basin_realised_richness_loss` | 4 / 0 / 40, mean **−0.001551** | 5 / 31 / 8, mean **−0.000037** |
| `basin_total_biomass` | 32 / 0 / 12, mean **+0.022244** | 31 / 0 / 13, mean **+0.066982** |
| `st_final_biomass` | 0 / 0 / 44, mean **−7.171** | 0 / 0 / 44, mean **−59.434** |
| `st_relative_biomass_change` | 0 / 0 / 44, mean **−0.012565** | 0 / 0 / 44, mean **−0.098449** |
| `st_quasi_extinct_fraction` | 1 / 37 / 6, mean **−0.001578** | 30 / 14 / 0, mean **+0.009569** |

### 5.1 Sign-flip runs

`basin_native_richness` — 30 flips (26 strictly opposite non-zero):
`ssp126/ACCESS-CM2, ssp126/CNRM-ESM2-1, ssp126/IITM-ESM, ssp126/KACE-1-0-G, ssp126/MIROC6,
ssp126/MPI-ESM1-2-HR, ssp126/MRI-ESM2-0, ssp126/UKESM1-0-LL, ssp245/CNRM-ESM2-1,
ssp245/EC-Earth3-Veg, ssp245/MIROC6, ssp245/MPI-ESM1-2-HR, ssp245/MRI-ESM2-0,
ssp370/ACCESS-CM2, ssp370/CNRM-ESM2-1, ssp370/EC-Earth3-Veg, ssp370/KACE-1-0-G,
ssp370/MIROC6, ssp370/MPI-ESM1-2-HR, ssp370/MRI-ESM2-0, ssp370/NorESM2-MM,
ssp585/CMCC-CM2-SR5, ssp585/CNRM-ESM2-1, ssp585/EC-Earth3-Veg, ssp585/IITM-ESM,
ssp585/KACE-1-0-G, ssp585/MIROC6, ssp585/MPI-ESM1-2-HR, ssp585/MRI-ESM2-0,
ssp585/NorESM2-MM`.

`basin_realised_richness_loss` — 37 flips (6 strictly opposite):
`ssp126/ACCESS-CM2, ssp126/CMCC-CM2-SR5, ssp126/CNRM-ESM2-1, ssp126/EC-Earth3-Veg,
ssp126/KACE-1-0-G, ssp126/MIROC6, ssp126/MRI-ESM2-0, ssp126/UKESM1-0-LL, ssp245/ACCESS-CM2,
ssp245/CNRM-ESM2-1, ssp245/EC-Earth3-Veg, ssp245/KACE-1-0-G, ssp245/MIROC6,
ssp245/MPI-ESM1-2-HR, ssp245/MRI-ESM2-0, ssp245/NorESM2-MM, ssp245/UKESM1-0-LL,
ssp370/ACCESS-CM2, ssp370/CMCC-CM2-SR5, ssp370/CNRM-ESM2-1, ssp370/EC-Earth3-Veg,
ssp370/KACE-1-0-G, ssp370/MIROC6, ssp370/MPI-ESM1-2-HR, ssp370/MRI-ESM2-0,
ssp370/NorESM2-MM, ssp370/UKESM1-0-LL, ssp585/ACCESS-CM2, ssp585/CMCC-CM2-SR5,
ssp585/CNRM-ESM2-1, ssp585/EC-Earth3-Veg, ssp585/IITM-ESM, ssp585/KACE-1-0-G,
ssp585/MIROC6, ssp585/MRI-ESM2-0, ssp585/NorESM2-MM, ssp585/UKESM1-0-LL`.
Nearly all of these are shared-K non-zero → alternative exactly zero (31 of 44 alternative
deltas are zero).

`basin_total_biomass` — 3 flips (all strictly opposite):
`ssp245/CNRM-ESM2-1, ssp370/CNRM-ESM2-1, ssp370/MPI-ESM1-2-HR`.

`st_final_biomass` — **0**; `st_relative_biomass_change` — **0** (all 44 negative under both).

`st_quasi_extinct_fraction` — 32 flips (3 strictly opposite):
`ssp126/ACCESS-CM2, ssp126/CNRM-ESM2-1, ssp126/EC-Earth3-Veg, ssp126/KACE-1-0-G,
ssp126/MRI-ESM2-0, ssp126/UKESM1-0-LL, ssp245/ACCESS-CM2, ssp245/CNRM-ESM2-1,
ssp245/EC-Earth3-Veg, ssp245/KACE-1-0-G, ssp245/MIROC6, ssp245/MRI-ESM2-0,
ssp245/NorESM2-MM, ssp245/UKESM1-0-LL, ssp370/ACCESS-CM2, ssp370/CNRM-ESM2-1,
ssp370/EC-Earth3-Veg, ssp370/KACE-1-0-G, ssp370/MIROC6, ssp370/MPI-ESM1-2-HR,
ssp370/MRI-ESM2-0, ssp370/NorESM2-MM, ssp370/UKESM1-0-LL, ssp585/ACCESS-CM2,
ssp585/CMCC-CM2-SR5, ssp585/CNRM-ESM2-1, ssp585/EC-Earth3-Veg, ssp585/IITM-ESM,
ssp585/KACE-1-0-G, ssp585/MIROC6, ssp585/NorESM2-MM, ssp585/UKESM1-0-LL`.
Most are shared-K zero → alternative small positive (14 alternative zeros remain).

Per-run values behind every flag are in `competition_structure_deltas.csv` columns `delta_<metric>`, `alternative_delta_<metric>`, `delta_diff_<metric>`, `sign_flip_<metric>`, and in the long-form `competition_sign_flips.csv` (102 rows).

---

## 6. Endpoints

`Guadex.select_climate_endpoints` on the primary index selects, from the realised-forcing
column `realised_warming_2036_2045_mean_degc`:

| endpoint | scenario / GCM | realised ΔT 2036–2045 (°C) |
|---|---|---|
| cool | `ssp126 / IITM-ESM` | 0.180427 |
| warm | `ssp245 / UKESM1-0-LL` | 1.292140 |

The alternative index selects exactly the same two runs (same scenario, same GCM, same
warming); `source = "realised"` in both cases. Selection provenance is saved as
`results/phase2/competition_endpoint_selection.txt`.

Endpoint scenario-minus-control deltas:

| metric | cool shared-K | cool species-specific | warm shared-K | warm species-specific |
|---|---|---|---|---|
| `basin_native_richness` | 0.000000 | −0.006452 | −0.006452 | −0.006452 |
| `basin_realised_richness_loss` | −0.000645 | −0.000323 | −0.0000645 | 0.000000 |
| `basin_total_biomass` | +0.001338 | +0.016809 | +0.088862 | +0.450550 |
| `st_final_biomass` | −2.595 | −28.973 | −10.692 | −74.563 |
| `st_relative_biomass_change` | −0.004546 | −0.047992 | −0.018734 | −0.123510 |
| `st_quasi_extinct_fraction` | 0.000000 | 0.000000 | 0.000000 | +0.013158 |

At the **warm endpoint the basin native-richness delta is numerically identical under both
structures: −0.006452 species** (= −5/775 basin sites). At the cool endpoint the primary is
exactly zero while the alternative is −0.006452. The endpoint total-biomass and trout signals
are larger under the alternative (warm: biomass +0.089 → +0.451; ST relative change
−0.019 → −0.124; ST final biomass −10.7 → −74.6).

---

## 7. Does the weak-warming conclusion survive, reverse, or amplify?

**Magnitude: it survives, and remains small in absolute terms.** Across all 44 scenarios the
basin native-richness delta is at most 0.0129 species (alternative) or 0.0116 (shared-K), i.e.
≤ 10 of 775 basin sites and < 0.5 % of the ≈3-species basin mean. The warm-endpoint response is
identical under both structures. No scenario under either structure produces a large basin
richness loss.

**Sign of the average richness response: it reverses.** Shared-K gives a mixed, slightly
positive mean delta (+0.0020; 26 positive, 4 zero, 14 negative), i.e. on average no native
richness loss. Species-specific gives a uniformly negative response (44/44 negative, mean
−0.0071). The *typical* run therefore changes from “no loss” to “a small but consistent loss”,
and 30/44 runs flip sign. The headline “weak warming” is preserved; the stronger claim “warming
does not reduce native richness” is not robust to the competition structure.

**Realised richness-loss response: it amplifies down to zero.** Shared-K has a systematic
negative delta (40/44 negative, mean −0.00155); under species-specific 31/44 deltas are exactly
zero and the mean is −3.7×10⁻⁵. The small richness-loss signal seen under the primary does not
survive the alternative; it is a property of the shared-capacity coupling, not a robust warming
effect.

**Total biomass: direction survives, magnitude amplifies.** Both structures give a mean
positive delta; the alternative is ~4.8× larger in mean absolute terms and 3/44 runs flip sign.

**Trout: direction survives, magnitude amplifies strongly.** All 44 runs lose trout biomass
under both structures (0 flips), but the alternative's losses are ~8× larger:
`st_relative_biomass_change` mean delta −0.0126 → −0.0984 (≈ −1.3 → −9.8 percentage points of
the observed baseline), and `st_final_biomass` mean delta −7.2 → −59.4. Part of the absolute
difference reflects the much higher alternative control level (ST control final biomass 759 vs
358); the relative-change comparison is the level-free one and still amplifies ~7.8×.

**Net reading.** The competition structure does not overturn the weak-warming headline for
basin native richness: the endpoint change is tiny and, at the warm end, numerically identical.
It does change the sign of the ensemble-average richness response, removes the small
richness-loss signal, and materially amplifies the biomass and trout responses. The report
should therefore state the structural dependence explicitly, with the shared-K result as the
primary number.

---

## 8. Caveats

- **Primary remains `c = 1 − α`.** The species-specific structure is a sensitivity only. Every
  primary result directory (`results/climate_scenarios_phase2`, `results/climate_scenarios_corrected`,
  the obstacle/interaction sweeps) was read only; this task wrote only to
  `results/climate_scenarios_competition_alt` and `results/phase2/`.
- **Shared controls.** Per decision 2, the matched per-GCM controls (P2.2) are deferred; each
  structure is differenced against its own shared `control/baseline` row. The P2.2 offset
  applies equally to both structures and is not removed here.
- **Interim route.** Both ensembles use the corrected interim route: observed start, 3 baseline
  years, scenario-minus-control reporting, six biological options, on-path network, daily
  per-GCM forcing. The forcing is bit-identical across the two structures.
- **ST metric source.** `st_final_biomass`, `st_relative_biomass_change` and
  `st_quasi_extinct_fraction` are read from `export/levels/quasi_extinction_summary.csv`
  (species `ST`). They are not columns of the climate `runs_index.csv`, and `level_basin.csv`
  has only basin-wide aggregates. The control values for these metrics are also read the same
  way.
- **Sign-flip definition.** A flip is `sign(δ_alt) ≠ sign(δ_shared)`, so a shared-K non-zero
  turning into an alternative zero counts as a flip (this inflates the raw counts for
  `basin_realised_richness_loss`, where 31/44 alternative deltas are zero, and for
  `st_quasi_extinct_fraction`, where 37/44 shared-K deltas are zero). The strictly opposite
  non-zero counts are 26, 6, 3, 0, 0 and 3 respectively.
- **Control states differ.** The structure changes the no-warming community (control richness
  2.64 → 3.25; ST control biomass 358 → 759). Deltas are the valid comparison quantity;
  cross-structure levels are not.
- **Coefficient range.** The competition + predation entries alone give the alternative
  `c ∈ [0.3, 0.5]`. The corrected table also carries `−0.8` displacement rows, so the full
  alternative range is `[0.3, 0.8]` (shared-K `[1.3, 1.8]`).
- No source code was changed; no obstacle or interaction sweep was run for this item.

---

## 9. Output files

| path | content |
|---|---|
| `results/climate_scenarios_competition_alt/runs_index.csv` | 45-run alternative index (46 lines) |
| `results/climate_scenarios_competition_alt/<scenario>/<gcm>/` | 45 complete runs: `simulation_output.jld2`, `export/run_metadata.json`, full exports, `run_digest.txt` |
| `results/phase2_competition_alt_run.log` | ensemble log, ends “Climate scenarios complete: 45 run(s)” |
| `results/phase2/competition_structure_deltas.csv` | 88 rows = 44 scenarios × 2 structures; values, deltas against each structure's control, Δδ and sign-flip flags |
| `results/phase2/competition_endpoint_summary.csv` | 24 rows = 6 metrics × 2 endpoints × 2 structures; values, controls, deltas and the per-metric Δδ statistics |
| `results/phase2/competition_sign_flips.csv` | 102 rows = every sign-flip run per metric, with both structures' deltas and a `strictly_opposite` flag |
| `results/phase2/competition_endpoint_selection.txt` | endpoint selection from both indices |
| `results/phase2/analyze_competition_structure.jl` | analysis script (build artefact) |
| `results/phase2_competition_analysis.log` | analysis log |

---

## 10. Reproduction

1. Set the four environment variables in §1 and run `julia --project=. run_climate_scenarios.jl`;
   success is “Climate scenarios complete: 45 run(s)”.
2. Verify 45 runs with the three artifacts, `resolved_config.competition_structure =
   "species_specific"`, `projection_route = "interim_observed"`,
   `interim_spin_up_years = 3`, `code_version = c340398`, and a digest different from the
   primary run at the same relative path.
3. Run `julia --project=. results/phase2/analyze_competition_structure.jl`. It recomputes each
   structure's control from its own index, reads the ST rows, joins on `scenario,gcm`, and calls
   `Guadex.select_climate_endpoints` on both indices.
4. The CSVs in §9 fully specify the numbers in §5–§6; the script is deterministic.
