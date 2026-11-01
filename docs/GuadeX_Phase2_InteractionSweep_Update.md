# GuadeX Phase 2 — P2.1c: Parser-Fix Effect on the Alternative-Interaction Sweep

**Scope:** re-run the 48-run alternative-interaction sweep (`original` / `random` /
`invasive_favoring` matrices × cool/warm climate × upstream cost × passability) against the
regenerated interaction table, and quantify — not assert — what the interaction-parser fix
changed in it.

**Status:** companion evidence for Phase 2 item P2.1 (§4 of
`.kilo/plans/1793526105792-second-review-triage-plan.md`). This document does not revise the
report or the status record; that is P2.6. It completes the P2.1 series:
`docs/GuadeX_Phase2_ParserFix_Effect.md` (climate ensemble) and
`docs/GuadeX_Phase2_ObstacleSweep_Update.md` (obstacle sweep).

---

## 1. Provenance and commands

All work was done in `E:\projects\GuadeX` (Windows, PowerShell 7, Julia 1.12.6) on 2026-11-01.
`git rev-parse --short HEAD` = **`c340398`**; the parser fix and the other Phase 1 edits are
uncommitted working-tree changes, as in the handover.

Pre-fix audit baseline (read-only, untouched — run directories 30 Sep 2026, index 29 Oct 2026):

- `results/alt_interactions_corrected/runs_index.csv` — **49 lines** (header + 48 runs),
  48 unique run keys.
- `results/alt_interactions_corrected/burnin_summary.csv` — 3 rows, `years = 3`,
  `criteria = "basin"` for all three matrices.
- Run metadata carries `parameter_digest` (24 distinct values) but **no `code_version` key**
  and none of the Phase-1 identity keys (`competition_structure`, `passability_floor`,
  `obstacles_on_path_only`). The pre-fix code version is documented as `9702169` in the Phase 2
  handover and in the climate metadata of the same correction series.
- Example: `original/sig_0.3/ssp126__IITM-ESM/uc_0.01/baseline` →
  `parameter_digest = "314e1de67e715462"`.

New run (started as a tracked background process; no early stop):

```powershell
$env:GUADEX_PARAMETERS_FILE='parameters_climate_scenarios_corrected.toml'
$env:GUADEX_ALT_OUTPUT_DIR='results/alt_interactions_phase2'
julia --project=. run_alt_interactions.jl *> results/phase2_alt_interactions_run.log
```

- Background process `bgp_197426a23001X1jv63EmThDIJd`, exited normally.
- Wall time **21 min 47 s** (log `CreationTime` 14:17:54 → `LastWriteTime` 14:39:42 local); the
  script's own counter reports **1234.3 s** (20 min 34 s), the difference being Julia/package
  start-up and data loading.
- Final log line: **"Alternative-interaction sensitivity complete: 48 run(s) under
  'results/alt_interactions_phase2' in 1234.3s"**. Log has 226 lines and zero matches for
  `ERROR`, `LoadError`, `failed` or `Stacktrace`.
- One benign warning: `CSV ... parsed expected 10 columns, but didn't reach end of line around
  data row: 1. Parsing extra columns and widening final columnset` while reading
  `data/Matrix_distances_1037puntos_BRUTO_FINAL.csv`. The same warning occurs in the obstacle
  sweep log (`results/phase2_obstacle_run.log`), i.e. it is a pre-existing property of that
  input file, not of this run.
- A `GUADEX_CONFIG_ONLY=1` smoke test resolved the expected grid before the run: 2 climate
  models (`ssp126/IITM-ESM` cool, `ssp245/UKESM1-0-LL` warm) × 2 upstream costs
  (`0.01, 0.5`) × 4 passability scenarios (`baseline`, `improved_passability`,
  `reduced_passability`, `blocked`), sigma multiplier `0.3`, route `:interim_observed`
  (3 baseline years), output `results/alt_interactions_phase2`. (The banner's `output:` line
  prints the config default `results/alt_interactions_corrected`; the env override is applied
  after the banner, at `run_alt_interactions.jl:290`, and the closing line is authoritative.)

---

## 2. New-run verification

| check | result |
|---|---|
| `results/alt_interactions_phase2/runs_index.csv` | **49 lines = header + 48 runs** |
| Grid | 3 matrices × 1 sigma (`0.3`) × 2 climate models × 2 upstream costs × 4 passability scenarios = **48 unique keys**, no duplicates |
| Per-matrix split | `original` 16, `random` 16, `invasive_favoring` 16 |
| Per-model split | `ssp126/IITM-ESM` 24, `ssp245/UKESM1-0-LL` 24 |
| Per-cost split | `0.01` 24, `0.5` 24 |
| Per-passability split | `baseline` 12, `improved_passability` 12, `reduced_passability` 12, `blocked` 12 |
| `burnin_summary.csv` | **3 rows**, `years = 3`, `criteria = "basin"` (interim route), `converged = false` for all three (expected: the E4 stop rule is bypassed on this route) |
| Required per-run artifacts | **48/48** have `simulation_output.jld2`, `export/run_metadata.json` and `.sensitivity_complete.jld2` |
| `run_digest.txt` | **0/48** — the sweeps do not write this sidecar (only `run_climate_scenarios.jl` does) |
| `competition_structure` (all 48) | `shared_k` |
| `passability_floor` (all 48) | `0.0` |
| `obstacles_on_path_only` (all 48) | `false` |
| `projection_route` / `interim_spin_up_years` (all 48) | `interim_observed` / `3` |
| `parameter_digest` | present in 48/48 (24 distinct values) |
| `code_version` key in `run_metadata.json` | **absent** in 48/48 (writer limitation, as in the obstacle sweep — see §3.3 of `GuadeX_Phase2_ObstacleSweep_Update.md`) |

### 2.1 Burn-in comparison

| matrix | pre-fix `total_biomass_change` | new | pre-fix `composition_q95` | new |
|---|---|---|---|---|
| `original` | 0.07142970466789052 | 0.07132353222849723 | 0.6731975124153126 | 0.6731325298053666 |
| `random` | 0.07404444746411394 | **0.07404444746411394** | 0.6711429784751932 | **0.6711429784751932** |
| `invasive_favoring` | 0.07175624478198782 | 0.07183791637233125 | 0.674552434974807 | 0.6731407286065982 |

The `random` burn-in is bit-identical, because the random matrix is generated from a fixed seed
and its scale `abs(minimum(original α))` is still `0.8` (the minimum α is −0.8 in both tables);
its diagnostics are therefore immune to the parser fix. The `original` drift is the α change and
the `invasive_favoring` drift is the AAL group change (§5.3), not the parser.

`criteria = "basin"` is by construction on this route: `interim_observed_spin_up`
(`src/spin_up.jl`) calls the spin-up with `criterion = :basin`, `tol = 0`, `composition_tol = 0`
and `max_years = interim_spin_up_years`, so the E4 two-criteria stop rule is bypassed and the run
integrates exactly 3 baseline years (`converged = false` is the expected flag). The composition
diagnostics are still recorded for information.

---

## 3. Digest evidence

The `parameter_digest` was read from every run's `export/run_metadata.json` in both sweeps and
joined on the run directory:

| quantity | result |
|---|---|
| runs in common | **48** |
| runs with a changed digest | **48/48** |
| distinct pre-fix digests | 24 |
| distinct post-fix digests | 24 |
| example (pre → post) | `original/sig_0.3/ssp126__IITM-ESM/uc_0.01/baseline`: `314e1de67e715462` → **`c071ab6d7ada9a9a`** |

There are 24 distinct digests per sweep, not 48, because the digest folds the interacting model
state (which includes the matrix), the two upstream costs and the four passability levels but
**not** `climate_scenario`/`gcm`; the two climate models therefore share a digest for the same
`(matrix, cost, passability)` cell. This is a pre-existing property of
`sensitivity_run_fingerprint` (`sensitivity_core.jl`), and it holds identically on both sides, so
it does not weaken the comparison. Full list:
`results/phase2/alt_interactions_digest_comparison.csv`.

The digests differ because the interaction matrix is a digest input (`src/run_identity.jl:112`)
and the code version moved from `9702169` to `c340398`.

### 3.1 Forcing unchanged

To show that the comparison isolates the interaction structure, the climate/forcing columns of
the two indices were compared run by run:

| column | rows differing (of 48) |
|---|---|
| `realised_warming_2026_2045_mean_degc` | **0** |
| `realised_warming_2036_2045_mean_degc` | **0** |
| `realised_warming_n_sites` | **0** |
| `warming_end_degc` | **0** |
| `basin_mean_temperature_c` | **0** |

---

## 4. Effect on the sweep metrics

### 4.1 Method

The two `runs_index.csv` files were inner-joined on
`interaction_matrix, thermal_sigma_multiplier, climate_scenario, gcm, upstream_cost,
passability_scenario` (48 matches). The alt-interaction sweep publishes **absolute final-year
basin metrics**, not scenario-minus-control deltas, so the comparison is a direct value
difference `delta_abs = post − pre`; `delta_rel = delta_abs / |pre|`.

Metrics analysed: `basin_native_richness`, `basin_realised_richness_loss`,
`basin_total_biomass`, `basin_native_biomass`, `basin_invasive_biomass`,
`st_final_biomass`, `st_relative_biomass_change`, `st_quasi_extinct_fraction`.

### 4.2 Summary by metric (48 runs each)

| metric | mean \|Δ\| | max \|Δ\| | mean \|Δ\| % | max \|Δ\| % | runs with Δ ≠ 0 | runs >1 % of value | sign changes |
|---|---|---|---|---|---|---|---|
| `basin_native_richness` | 0.000215 | 0.001290 | 0.0078 | 0.0605 | 8 | 0 | 0 |
| `basin_realised_richness_loss` | 0.0000067 | 0.000323 | 0.0028 | 0.135 | 1 | 0 | 0 |
| `basin_total_biomass` | 0.009017 | 0.033088 | 0.0253 | 0.0883 | 32 | 0 | 0 |
| `basin_native_biomass` | 0.006512 | 0.026560 | 0.0203 | 0.0793 | 32 | 0 | 0 |
| `basin_invasive_biomass` | 0.004099 | 0.013044 | 0.137 | 0.374 | 44 | 0 | 0 |
| `st_final_biomass` | 0.114132 | 0.507790 | 0.0259 | 0.113 | 32 | 0 | 0 |
| `st_relative_biomass_change` | 0.000153 | 0.000855 | 0.0487 | 0.276 | 32 | 0 | 0 |
| `st_quasi_extinct_fraction` | **0.0** | **0.0** | 0 | 0 | **0** | 0 | 0 |

No metric changes sign in any run, no run moves by more than 1 % of its own value, and
`st_quasi_extinct_fraction` is bit-identical in all 48 runs (the 72-site ST denominator never
gains or loses a site). Scales: `0.0012903 = 1/775` (one site-presence flip on the 775-site
basin mean); `0.013889 = 1/72` (one site on the baseline-established ST set).

### 4.3 Decomposition by interaction matrix

This is the key table: the three matrices answer different questions about what changed.

| metric | `original` mean \|Δ\| / max | `random` mean \|Δ\| / max | `invasive_favoring` mean \|Δ\| / max |
|---|---|---|---|
| `basin_native_richness` | 0.000403 / 0.001290 (5/16 nonzero) | **0 / 0** | 0.000242 / 0.001290 (3/16 nonzero) |
| `basin_realised_richness_loss` | **0 / 0** | **0 / 0** | 0.0000202 / 0.000323 (1/16) |
| `basin_total_biomass` | 0.007224 / 0.016307 | **0 / 0** | 0.019826 / 0.033088 |
| `basin_native_biomass` | 0.002825 / 0.006891 | **0 / 0** | 0.016710 / 0.026560 |
| `basin_invasive_biomass` | 0.009179 / 0.013044 | 4.2e-16 / **1.3e-15** | 0.003118 / 0.006946 |
| `st_final_biomass` | 0.184367 / 0.443185 | **0 / 0** | 0.158029 / 0.507790 |
| `st_relative_biomass_change` | 0.000225 / 0.000582 | **0 / 0** | 0.000233 / 0.000855 |
| `st_quasi_extinct_fraction` | **0 / 0** | **0 / 0** | **0 / 0** |

- **`random` is numerically identical** (the only non-zero cells are 1e-15 floating-point noise
  in invasive biomass). The random matrix is drawn from a fixed seed and its scale is unchanged,
  so this row also proves that the Phase-1 code changes between `9702169` and `c340398` are
  **numerically inert** for this sweep: any change in the other two rows is a change in the
  *matrix*, not in the solver, the burn-in or the forcing.
- **`original` differs only through the α table.** This is the direct parser-fix effect: the 16
  reverse directions moved from −0.5 predation to −0.3 competition (less heterospecific
  suppression), which raises invasive biomass (max +0.0130) and lowers native biomass and ST
  biomass slightly (max −0.0069 and −0.44 on a ≈448 scale), with one site of 775 flipping in
  native richness in 5 of 16 runs.
- **`invasive_favoring` is not a parser effect at all.** Its matrix is *constructed* from the
  native/invasive species groups, not parsed from the table, so the α fix cannot enter it. Its
  drift is the AAL reclassification (commit `099a87d`, 29 Oct 2026: `AAL` moved from the
  `migratory` group to the `invasive` group) — see §5.3 for the exact reproduction.

---

## 5. Largest changes

### 5.1 `original` matrix — the parser signal

| metric | run | pre | post | Δ |
|---|---|---|---|---|
| `basin_native_richness` | ssp126/IITM-ESM, uc=0.01, baseline | 3.178065 | 3.176774 | **−0.001290** (= −1/775) |
| `basin_native_richness` | ssp245/UKESM1-0-LL, uc=0.01, baseline | 3.185806 | 3.184516 | −0.001290 |
| `basin_total_biomass` | ssp245/UKESM1-0-LL, uc=0.01, blocked | 39.125141 | 39.141448 | **+0.016307** |
| `basin_native_biomass` | ssp245/UKESM1-0-LL, uc=0.5, improved_passability | 30.804194 | 30.797303 | −0.006891 |
| `basin_invasive_biomass` | ssp245/UKESM1-0-LL, uc=0.01, blocked | 3.681588 | 3.694632 | **+0.013044** |
| `st_final_biomass` | ssp126/IITM-ESM, uc=0.01, baseline | 468.056562 | 467.613377 | **−0.443185** |
| `st_relative_biomass_change` | ssp126/IITM-ESM, uc=0.01, baseline | −0.281254 | −0.281836 | −0.000582 |

The invasive-biomass increase is the expected direction: the corrected reverse directions are
weaker competition (−0.3 instead of −0.5), so the 11 changed non-AAL directed cells
(`CG→AM/CC/EL/LG`, `GH→AM/CC/LG`, `LG→OM`, `TT→AM/CC/LG`) suppress natives and other invasives
less. The invasive-biomass delta is positive in 16/16 `original` runs and the total-biomass delta
in 16/16; the native-biomass delta is negative in 11/16 (the other 5 are small positives).

### 5.2 `invasive_favoring` matrix — the AAL signal

| metric | run | pre | post | Δ |
|---|---|---|---|---|
| `basin_total_biomass` | ssp126/IITM-ESM, uc=0.01, blocked | 37.467951 | 37.434862 | **−0.033088** |
| `basin_native_biomass` | ssp126/IITM-ESM, uc=0.01, blocked | 33.487648 | 33.461088 | **−0.026560** |
| `basin_invasive_biomass` | ssp245/UKESM1-0-LL, uc=0.01, blocked | 3.814061 | 3.807114 | −0.006946 |
| `st_final_biomass` | ssp245/UKESM1-0-LL, uc=0.01, improved_passability | 448.386207 | 447.878417 | **−0.507790** |
| `st_relative_biomass_change` | ssp245/UKESM1-0-LL, uc=0.01, improved_passability | −0.310445 | −0.311300 | −0.000855 |
| `basin_realised_richness_loss` | ssp245/UKESM1-0-LL, uc=0.5, improved_passability | 0.238538 | 0.238860 | +0.000323 |

All 16 `invasive_favoring` runs move in the same direction (native and total biomass down), i.e.
a systematic matrix change rather than noise.

### 5.3 Resolving the `invasive_favoring` residual

The pre-fix sweep predates commit `099a87d`: at the pre-fix code state
(`git show d58424f:parameters_climate_scenarios_corrected.toml`) the corrected configuration
listed `migratory = ["AA", "AAL", "LR", "MC"]`; the current file lists
`invasive = [..., "TT", "AAL"]` and `migratory = ["AA", "LR", "MC"]`.
`make_invasive_favoring_matrix` builds −0.8 native-suppression entries from the native/invasive
index sets, so adding AAL to the invasive group changes the constructed matrix even though no
parsed α value is involved.

This was verified by rebuilding the 3-year interim burn-in twice in the current tree
(`results/phase2/diagnose_invasive_residual.jl`):

| variant | years | `total_biomass_change` | `composition_q95_change` | active cells |
|---|---|---|---|---|
| current groups (AAL invasive) | 3 | 0.07183791637233125 | 0.6731407286065982 | 2681 |
| AAL excluded from the invasive group | 3 | **0.07175624478198782** | **0.674552434974807** | **2682** |
| *pre-fix recorded* | 3 | *0.07175624478198782* | *0.674552434974807* | *2682* |
| *new run recorded* | 3 | *0.07183791637233125* | *0.6731407286065982* | *2681* |

The "AAL excluded" variant reproduces the pre-fix burn-in diagnostics **to the last digit**, and
the current variant reproduces the new run's. The `invasive_favoring` row of §4.3 is therefore
the AAL group change, not the parser fix. (AAL is effectively absent from the runs — its basin
densities are ~1e-35 — which is why this group change does not move `basin_invasive_biomass`
materially and why the `random` matrix's invasive biomass is unchanged at 1e-15.)

---

## 6. Matrix ordering: empirical vs random vs invasive-favouring

### 6.1 Metric-mean ordering (16 runs per matrix)

| metric | pre-fix order | new order | rank changes | Spearman ρ |
|---|---|---|---|---|
| `basin_native_richness` | original > invasive_favoring > random | same | 0 | 1.00 |
| `basin_realised_richness_loss` | random > invasive_favoring > original | same | 0 | 1.00 |
| `basin_total_biomass` | original > invasive_favoring > random | same | 0 | 1.00 |
| `basin_native_biomass` | original > invasive_favoring > random | same | 0 | 1.00 |
| `basin_invasive_biomass` | invasive_favoring > original > random | same | 0 | 1.00 |
| `st_final_biomass` | original > invasive_favoring > random | same | 0 | 1.00 |
| `st_relative_biomass_change` | original > invasive_favoring > random | same | 0 | 1.00 |
| `st_quasi_extinct_fraction` | random > invasive_favoring > original | same | 0 | 1.00 |

The ordering is identical for all eight metrics (ρ = 1.00, n = 3 factors). The mean values move
only marginally:

| metric | original (pre → post) | random (pre → post) | invasive_favoring (pre → post) |
|---|---|---|---|
| `basin_native_richness` | 2.605645 → 2.605403 | 2.563306 → 2.563306 | 2.592419 → 2.592177 |
| `basin_total_biomass` | 35.384988 → 35.392212 | 33.482405 → 33.482405 | 34.974510 → 34.954684 |
| `basin_native_biomass` | 32.276657 → 32.274715 | 30.699318 → 30.699318 | 31.747007 → 31.730297 |
| `basin_invasive_biomass` | 2.959381 → 2.968560 | 2.634148 → 2.634148 | 3.078559 → 3.075442 |
| `st_final_biomass` | 445.877461 → 445.693570 | 432.497975 → 432.497975 | 438.409577 → 438.433738 |
| `st_relative_biomass_change` | −0.315312 → −0.315500 | −0.333604 → −0.333604 | −0.325787 → −0.325823 |
| `st_quasi_extinct_fraction` | 0.2578125 → 0.2578125 | 0.2682292 → 0.2682292 | 0.2612847 → 0.2612847 |

The largest ordering-relevant move is `basin_invasive_biomass`: the empirical matrix gains
(+0.0092) and the invasive-favouring matrix loses (−0.0031), narrowing the gap between them from
0.1192 to 0.1069 — but `invasive_favoring` still has the most invasive biomass in every block.

### 6.2 Per-block ordering

Within each of the 16 `(climate model × upstream cost × passability)` blocks, the three matrices
were ranked by each metric. The ordering is **identical in 16/16 blocks for all eight metrics**
(0 changes). The α fix therefore does not reorder the matrices anywhere in the design, not merely
on the pooled mean.

### 6.3 What the ordering means for the interaction-robustness conclusion

The robustness argument used in the report is that the qualitative climate/obstacle story
survives *three structurally different interaction matrices*: the empirical one, a random signed
matrix, and one deliberately built to favour invasives. That argument depends on the **ordering
and the response directions**, not on the third decimal of any metric. Both are preserved
exactly:

- empirical > invasive-favouring > random for native richness, native biomass, total biomass,
  `st_final_biomass` and `st_relative_biomass_change` — 16/16 blocks;
- invasive-favouring > empirical > random for invasive biomass — 16/16 blocks;
- for realised richness loss and ST quasi-extinction the *pooled-mean* order is
  random > invasive-favouring > empirical, but it is not uniform across blocks (richness loss:
  11/16 blocks in that order, 5/16 with the last two swapped; ST quasi-extinction: three
  patterns at 4/5/7 blocks). What matters is that **each block keeps its own pre-fix ordering**
  (16/16 for all eight metrics), so neither the pooled order nor any block's order is disturbed;
- zero sign flips in any of the 48 runs on any of the eight metrics.

---

## 7. Plain-English interpretation

**The parser fix did not change the interaction-robustness conclusion.**

The corrected interaction table fixes 16 reverse-direction entries (−0.5 predation → −0.3
competition) and flags 9 mixed cells (34 changed directed cells in total). In the 48-run
alternative-interaction sweep that moves the final-year basin metrics by at most:

- **0.0013 native species** of basin mean richness — exactly 1 site of 775, in 8 of 48 runs;
- **0.00032** of realised richness loss — 1 run, and 0 in the 16 `original`-matrix runs;
- **0.033 biomass units** on a ≈35-unit basin mean (0.09 %);
- **0.027** native biomass, **0.013** invasive biomass;
- **0.51** ST final biomass on a ≈448 scale (0.11 %);
- **exactly 0.0** for the ST quasi-extinct fraction in every run.

The applied forcing and basin mean temperature are **bit-identical** before and after (0
differing rows in five forcing columns), and the **`random`-matrix runs are bit-identical**
(1e-15 excepted), which proves the Phase-1 code changes are numerically inert here. The change in
the empirical (`original`) matrix is therefore the parser effect, and it is one order of magnitude
below the design's own factor effects (upstream cost moves native richness by ≈0.85 species in
the obstacle sweep).

The third matrix (`invasive_favoring`) moves more than the parser effect would suggest, but that
is **not** a parser effect: the matrix is constructed from the species groups, and the pre-fix
baseline predates the AAL reclassification (099a87d). Rebuilding the burn-in with AAL excluded
reproduces the pre-fix diagnostics exactly (§5.3). This matters for reporting: the pre-fix
interaction sweep is an audit baseline for the *parser*, but it also carries the older AAL
grouping, so its `invasive_favoring` numbers must not be presented as "the same model minus the
parser fix".

The three-matrix ordering — the actual robustness test — is unchanged for all eight metrics
(Spearman ρ = 1.00) and in all 16 design blocks. The material structural sensitivity for this
experiment remains the competition form (P2.3), the sigma multiplier and the passability floor
(P2.4), not this table repair.

---

## 8. Caveats

- The sweeps publish **absolute final-year basin metrics**, not scenario-minus-control deltas,
  so the comparison is a direct `post − pre` difference, as in
  `GuadeX_Phase2_ObstacleSweep_Update.md` §3.1.
- Neither sweep's `run_metadata.json` records a raw `code_version` string; provenance is the
  `parameter_digest` (which folds the code version and the model arrays) plus the
  `.sensitivity_complete.jld2` fingerprint. `run_digest.txt` is not written by the sweeps.
- The digest is 24-distinct rather than 48-distinct because it excludes `climate_scenario`/`gcm`
  (pre-existing `sensitivity_run_fingerprint` behaviour, identical on both sides).
- The pre-fix baseline carries the **pre-099a87d species grouping** (AAL migratory). For
  `original` and `random` this is immaterial (AAL biomass ~1e-35); for `invasive_favoring` it is
  the entire delta (§5.3). Only the `original`-matrix row is a clean parser-fix measurement.
- `basin_realised_richness_loss` is unchanged in all 16 `original` runs and changes in exactly
  one `invasive_favoring` run; a relative percentage on a 0.0003 change is not interpretable,
  the absolute column is.
- The run banner prints the config-default output root (`results/alt_interactions_corrected`)
  rather than the `GUADEX_ALT_OUTPUT_DIR` override; the closing line names
  `results/alt_interactions_phase2` and every run directory is under it.
- `results/alt_interactions_corrected`, `results/sensitivity_obstacles_corrected`,
  `results/climate_scenarios_corrected` and `data/BIOTIC/interaction_matrix_long_prefix.csv`
  were read only; their timestamps are unchanged (run directories 30 Sep 2026;
  `runs_index.csv` and `report_plots/` last written 29 Oct 2026).
- No model source was changed for this task. The two scripts written to `results/phase2/`
  (`analyze_alt_interactions_effect.jl`, `diagnose_invasive_residual.jl`) are build artefacts
  under the gitignored `results/` tree.

---

## 9. Output files

| path | content |
|---|---|
| `results/alt_interactions_phase2/runs_index.csv` | new 48-run index (49 lines) |
| `results/alt_interactions_phase2/burnin_summary.csv` | 3 burn-in rows, `years = 3`, `criteria = basin` |
| `results/alt_interactions_phase2/<matrix>/sig_0.3/<model>/uc_<cost>/<pass>/` | 48 complete runs (`simulation_output.jld2`, `export/run_metadata.json`, full exports, `.sensitivity_complete.jld2`) |
| `results/phase2_alt_interactions_run.log` | run log, ends "Alternative-interaction sensitivity complete: 48 run(s)" |
| `results/phase2/alt_interactions_prefix_vs_phase2.csv` | 384 rows = 48 runs × 8 metrics; `interaction_matrix, thermal_sigma_multiplier, climate_scenario, gcm, upstream_cost, passability_scenario, metric, pre, post, delta_abs, delta_rel` |
| `results/phase2/alt_interactions_effect_summary.csv` | per-metric delta summary (13 columns) |
| `results/phase2/alt_interactions_matrix_order.csv` | per-metric matrix means and ranks, pre and post |
| `results/phase2/alt_interactions_matrix_order_comparison.csv` | per-metric ordering, rank changes, Spearman |
| `results/phase2/alt_interactions_matrix_block_order.csv` | per-metric identical/changed block counts |
| `results/phase2/alt_interactions_digest_comparison.csv` | 48 run-dir digest pairs with a `differs` flag |
| `results/phase2/analyze_alt_interactions_effect.jl` | the analysis script that produced the above |
| `results/phase2/analyze_alt_interactions_effect.log` | analysis stdout |
| `results/phase2/diagnose_invasive_residual.jl` | the §5.3 AAL-group reproduction |

---

## 10. Reproduction

1. Pre-fix baseline: read `results/alt_interactions_corrected/runs_index.csv` (49 lines) and
   `burnin_summary.csv` (3 rows, `years = 3`); the run metadata carries `parameter_digest` but no
   `code_version`.
2. New run: `GUADEX_PARAMETERS_FILE=parameters_climate_scenarios_corrected.toml`,
   `GUADEX_ALT_OUTPUT_DIR=results/alt_interactions_phase2`, then
   `julia --project=. run_alt_interactions.jl`. Success is the log line
   "Alternative-interaction sensitivity complete: 48 run(s)".
3. Deltas: `julia --project=. results/phase2/analyze_alt_interactions_effect.jl` — inner-joins
   the two indices on the six key columns, writes the six CSVs, and prints the digest,
   ordering and largest-change tables used above.
4. AAL attribution: `julia --project=. results/phase2/diagnose_invasive_residual.jl` rebuilds the
   `invasive_favoring` matrix with and without AAL in the invasive group and compares the 3-year
   interim burn-in diagnostics against both baselines.
