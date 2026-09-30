# GuadeX Review Assessment

**Adjudication of `docs/GuadeX_Review_and_Corrections_Sept2026.md` against the repository at commit `9702169`.**

This document records, point by point, which of the reviewer's claims are supported by the
code, data and published outputs and which are not. It is deliberately not a rubber stamp:
the reviewer is right about the substance of most of the critical issues, but several
specific figures are wrong, two claims contradict each other, and several of the
reviewer's own proposed fixes are technically incorrect. The companion document
`docs/GuadeX_Correction_Plan_Sept2026.md` turns the fair points into an ordered work plan.

## How this assessment was produced

Each block of the review (C1–C7, E1–E10, E11–E19, E20–E24, minor issues, the report
statements table, "What stands", and Part II) was checked independently against primary
evidence: the Julia source (`src/`), the run scripts, `parameters.toml`, the raw input
CSVs under `data/`, the water-temperature sub-project `guadex_tw/`, and the committed
simulation outputs under `results/`. Counts and aggregates were recomputed from the raw
files rather than taken from the review or the report.

Verdicts used throughout:

- **FAIR** — the claim is correct as stated (minor wording aside).
- **PARTLY FAIR** — the underlying concern is real, but a specific number, attribution or
  framing is wrong or overstated; or the claim is valid for one configuration but not
  another.
- **NOT FAIR** — the claim is incorrect, or is presented as an error when the code/outputs
  actually support the report.
- **UNVERIFIABLE** — the figure depends on data or scripts not present in the repository
  (e.g. the reviewer's own GIS/heavy computations), so it cannot be confirmed or refuted
  from the repo alone.

The strongest sub-claims were reproduced numerically; where that was not possible the
document says so rather than repeating the number.

---

## 1. Verdict summary

### Critical issues

| ID | Issue | Verdict | The part that is not accurate |
|----|-------|---------|-------------------------------|
| C1 | Interaction matrix loaded inverted and one-way | **FAIR (core); PARTLY (figures)** | "114 column-named" includes the skipped *Lepomis* row (101 loaded); row-named cells are 13, not 7; "53 competition" entries not reproduced (48–52); **"89 pairs, Σα = −65.5" contradicts the review's own *Lepomis* claim — the loaded matrix gives 79 pairs, Σα = −58.1**. |
| C2 | Interaction term not scaled by growth | **FAIR (code + magnitude); PARTLY (framing/fix)** | The interaction term genuinely sits outside `r·W`. But it is a documented formulation, not a code/equation mismatch; and the reviewer's proposed coefficient `c = −α` is **wrong** (should be `c = 1 − α`). |
| C3 | Initial state does not match observed community | **FAIR** | Model native richness is **2.51**, not 2.52; four observed-site cells in the table use raw presence (`>0`) rather than the stated `>0.1` threshold. |
| C4 | Dispersal graph does not follow the river network | **FAIR (mechanism + most numbers); PARTLY (79%; fix)** | The "79%" is UNVERIFIABLE (an on-path proxy gives 86.7%); "cannot be reproduced deterministically" is overstated (arbitrary but stable within an environment). The reviewer's fallback tree length **8,640 km is the parent-edge sum only — the full tree is 13,871 km**, so it is not "close to the 8,783 km MST". |
| C5 | Site temperature is 1961–1990 sub-catchment air climatology | **FAIR** | All numbers reproduced. Minor: the README column is `Tºm_SC`, not `TEMP_MEDIA_SC`; the README never literally says "air"; one sub-catchment (`32.1`) has two values. The fix needs plumbing the reviewer does not mention. |
| C6 | Thermal niches from distributional ranges, not physiology | **PARTLY FAIR** | Mechanism correct, but **15 species share "8 to 30", not 16**. Literature claims (growth optima, CTmax) are external and not verifiable from the repo. |
| C7 | Warming axis / exposure do not reflect applied forcing | **FAIR (points 1,2,4,5,6); PARTLY (point 3)** | The correlation **r = 0.81 is not reproducible** (annual ≈ 0.77; decadal 0.96–0.97). "Six gauging stations / 20-year windows" is loose (one window is 30 years). Concrete examples and the SSP3 < SSP1-2.6 inversion are correct. |

### Major issues E1–E24

| ID | Verdict | Key qualification |
|----|---------|-------------------|
| E1 | **FAIR** | Three temperature levels/ series confirmed; basin offset 1.27 °C, trout 2.60 °C — both exact. |
| E2 | **PARTLY FAIR** | Like-for-like concern valid; "halves variability" wrong (removes ~68%); "0.34–0.52 °C" not reproduced. |
| E3 | **FAIR** | One burn-in reused across all 32 runs; SS 0.96–0.97 confirmed. |
| E4 | **PARTLY FAIR** | Weak total-biomass criterion confirmed; but **"21-year interaction burn-in" is wrong — matrix-specific: 21 / 494 / 532 years**. |
| E5 | **FAIR** | Offsets 0–19 ⇒ "2045" is 1 Jan 2045; only 19 years of forcing. |
| E6 | **FAIR** | Denominator is all 775 sites; trout fraction 0.9561 = 1 − occupancy exactly. |
| E7 | **FAIR** | Cross-design standardisation and confounding confirmed. |
| E8 | **FAIR** | Deterministic, no replicates, metric is realised richness loss. |
| E9 | **FAIR** | Cache key/fingerprint omit parameter values and code version; existing-file skip never invalidates. |
| E10 | **PARTLY FAIR** | TOML genuinely does not reproduce the report. But "undocumented environment variables" is **wrong** — documented in `docs/climate_scenarios.md` and recorded in run metadata; the real gap is a committed scenario config. |
| E11 | **FAIR** | `IET` is bank stability, not trophic state; outlier 27.67 confirmed. |
| E12 | **PARTLY FAIR** | Uniform non-accumulating 0.1/0.5 and 51.7% legacy restriction confirmed; "99 demolished" and "797/968" not reproducible from flat files. |
| E13 | **FAIR** | 10% absence factor and 38.9% empty sites confirmed; "prevents expansion by construction" is overstated (0.1·r is still positive). |
| E14 | **FAIR** | 84 pool sites, 56.2% of density, 30 m² pool with 9,533 *Gambusia*. |
| E15 | **FAIR** | 56 `DIFÍCIL` sites exact; fishless count is 293, not 294. |
| E16 | **FAIR** | Eel/mullets `REPRODU...=no` yet local r; 4 of 9 eel records in the Guadiato fish-farm cluster. |
| E17 | **FAIR** | Mange/catfish citations confirmed; sedentary species 6–25 km/yr; *Alburnus* (AAL) absent from both dicts. |
| E18 | **FAIR** | Salinity/elevation read but unused; 477 sub-500 m trout-eligible sites exact. (Review's own caveat: adding an elevation envelope would block upslope shifts.) |
| E19 | **PARTLY FAIR** | 261 excluded sites and unused CEDEX confirmed; "dry or effluent-only" attribution unverifiable; **"17% runoff decline" lower bound not reproduced** (some increases). |
| E20 | **FAIR** | 257 coded 1/6–1/3 °C values confirmed; QC flags unused. But the SD drop 2.65→2.34 is the **Spain** model (the Guadalquivir calibration contains none of the 257). The Part II claim that flag `06a` "should be two-sided" is wrong — it already uses `abs()`. |
| E21 | **PARTLY FAIR** | Coefficients 0.51/−0.163 are the **national** model (Guadalquivir is 0.436/−0.299); lagged slopes 0.60/0.65 reproduce for Spain, but fit is better only for d3/d7, worse for d30; only 3 in-sample sites. "15–30%" is the reviewer's own inference. |
| E22 | **FAIR** | LOBO bias −1.106 °C, Feb −2.22 / Mar −2.29; script 13 omits script 10's corrections. |
| E23 | **FAIR** | 64.1% of Guadalquivir obs from 2024; 2 single-observation sites; 17 model sites above the 1,127 m calibration ceiling. |
| E24 | **FAIR (substance)** | "All sites same anomaly" is indeed wrong for the daily runs (per-site anomalies, r ≈ −0.8 with elevation); `z_ref = 0` confirmed. Exact 0.82–1.08 / r = −0.79 not reproduced. |

### Minor issues

| # | Verdict | Key qualification |
|---|---------|-------------------|
| 1 Rare-species site counts | **FAIR** | Exact. |
| 2 Juvenile-only presences | **PARTLY FAIR** | "25 cases" not reproduced (23 cells across 8 sites); IO has no juvenile column — confirmed. |
| 3 Row-order alignment | **FAIR** | The named function `initial_state_from_density` does not exist; the fragility is in three run scripts. |
| 4 Rate conversion | **FAIR** | Comment promises `log(1+r)/183` for *Gambusia*; code applies `r/365`. |
| 5 Calibration of k | **PARTLY FAIR** | Mean-year calibration confirmed; "5.9%" not independently reproduced. |
| 6 Calendar | **FAIR** | 7,305 forcing days vs 7,300 simulated; day 366 dropped. |
| 7 Monthly steps | **PARTLY FAIR** | AIC equivalence confirmed (3654.5 vs 3654.9); "1.55 °C step" not reproduced. |
| 8 Hot tail | **PARTLY FAIR** | Curvature exists but is ≈+0.09 °C (25–30 °C) and −0.43 °C (30–50 °C), not "≈0.7 °C above 25 °C". |
| 9 Mixed-model details | **FAIR** | REML LRT across fixed effects, ICC at 0 °C, in-sample models in CV table confirmed; "strengthens with ML" unverified. |
| 10 Numerical details | **FAIR** | Clipping, unclipped immigration vs clipped emigration, and the positivity callback creating mass all confirmed. |
| 11 Obstacle count / coordinates | **PARTLY FAIR** | ED50 vs ETRS89 split confirmed (site 1.14.10 shifts ~234 m, not 200 m); **the 966–986 count variation is not reproduced (stable 973)** — tie order changes which edge is assigned, not the matched count. |
| 12 Auxiliary data | **PARTLY FAIR** | `PERIMETRO_MOJADO`/precipitation and the bogus "SSP285" confirmed; "ten headers" and "863 m coordinate disagreement" not reproduced (coordinates agree). |
| 13 Code and tests | **FAIR** | Duplicate `AA`/missing `AAL`, vacuous tests, and missing mass/topology/sign/temperature tests all confirmed. |
| 14 Viewer | **FAIR** | Mixed 15-metric/20-year slider, no species export, README 1,037 vs 775 confirmed. |

### Report statements table and "What stands"

Most rows are **FAIR**: the reviewer's characterisation of the report is accurate for the
water-body count (289, not 774), the temperature level, the same-anomaly text, the habitat
label, the 16.0 °C mean, the warming axis, invasive richness 0.2168, the un-operable
blocking mechanism, the dispersal-multiple error (two migratory taxa are 0.5×, not >5×),
the quasi-extinction denominator, and the README counts.

Two rows need qualification:

- **"K = 1 × observed"** — **NOT FAIR as a needed revision.** The September climate runs
  really did use `carrying_capacity_base_scaling = 1` (confirmed in run metadata, mean
  K ≈ 86.6 vs observed ≈ 84.4). This is a reproducibility/`parameters.toml` issue (E10),
  not a report error.
- ***Alburnus* "coded Aa"** — **PARTLY FAIR.** The report uses full species names and
  contains no "Aa" code; the real, report-acknowledged defect is the *misclassification*
  as migratory. The "coded Aa" wording is the review's own construction.

The **"What stands"** positives are largely supported (water-temperature statistics,
Table-3 slopes/R², rare-species ranges, quasi-extinction arithmetic, dispersal median).
Some are marked only **UNVERIFIABLE** (matrix symmetry/completeness, sum-of-squares
formulas, Tsit5 settings, mass conservation), not verified.

### Part II proposed fixes

| Fix | Verdict | Problem |
|-----|---------|---------|
| C2 Lotka–Volterra rewrite (`c = 1 − α`) | **YES** | Correct and small; the earlier `c = −α` in the diagnosis is the error. |
| E5 offsets 1…20 | **YES** | Must be applied in all three run scripts; watch `searchsortedfirst` alignment. |
| E9 parameter digest | **YES** | Reuses existing digests in `sensitivity_core.jl`. |
| E22 monthly bias in script 13 | **YES** | Correct; residual is pooled over basins, so validate. |
| E10 IET rename, report text | **YES** | Cosmetic/labelling. |
| Mixed-model LRT `reml=False` | **YES** | Trivial. |
| C1 direction parsing | **PARTLY** | Regex bugs: two different patterns; Spanish `interfiere`/`depreda` and "…and X through predation" cells are silently mislabelled; mixed-mechanism cells are promised to be flagged but are not; **"iterate over all rows" does not fix *Lepomis*** (no column exists). |
| C4 `dendritic_parents` fallback | **PARTLY** | The snippet yields a forest of 267 roots; the claimed ~8,640 km counts only parent edges — the complete tree is **13,871 km**, and main-stem topology is left to an undirected MST. |
| C5 water-temperature level | **PARTLY** | Conceptually right, but `extract_site_temperatures` has no access to the guadex_tw table; needs new plumbing, a site-code join (776 vs 775), and the E22 correction applied consistently. |
| C7 realised forcing / exposure at occupied sites | **PARTLY** | Needs a persisted realised-anomaly field and a join to baseline-established cells. |
| E6 established-only quasi-extinction | **PARTLY** | The summary function has no `baseline_established`/denominator argument; that must be added. |
| Obstacles from inventory IF index | **NO** | IF semantics unconfirmed, max 7.97, and 445/1,658 rows are non-numeric (NA/DE/SD); cannot be implemented as `IF/10`. |
| River-line graph rebuild | **PARTLY** | `SW_Line_4C_.shp` is in-repo and geopandas is available, but no snap/flow code exists and overlays must be rebuilt. |
| Per-GCM controls, Sobol, input rebuilds | **PARTLY** | Feasible but compute- or decision-heavy; per-GCM controls need a historical series in each forcing file. |

---

## 2. Points that are not accurate or are over-stated

This is the material the user asked for: **not every correction should be accepted.** The
following are wrong, contradicted, or cannot be reproduced.

1. **C1, the "89 pairs, Σα = −65.5" figure.** The loaded matrix (which skips *Lepomis*)
   has 79 non-zero invasive→native effects summing to −58.1. The 89 / −65.7 values require
   the *Lepomis* row that C1 itself says is skipped. The review contradicts itself. The
   headline claim — invasive→native effects are zero in all 100 pairs — is exact.
2. **C1, "iterate over all rows" fixes *Lepomis*.** False: the CSV has a `Lg` row but no
   `Lg` column, and the guard tests membership in the column list, so the row is skipped
   regardless. A column (or a long-format table) is required.
3. **C2, numerical coefficient.** The recommended `c_sj = −α_sj` is wrong; the correct
   competitive-LV mapping is `c_sj = 1 − α_sj` (with `c_ss = 1`). The review's own
   verification identity, `K(1−c₁₂)/(1−c₁₂c₂₁)`, only holds for the corrected mapping.
4. **C4, "total link length 31,215 km against 8,783 km MST" combined with a prototyped
   tree "of about 8,640 km".** The 31,215 km and 8,783 km figures reproduce. But the
   proposed tree's parent edges sum to 8,637 km and the **full tree is 13,871 km** once the
   267 roots are joined; a spanning tree cannot be lighter than the MST, so the review's
   "close to MST" conclusion is internally inconsistent. The fallback needs an explicit,
   documented root-joining rule.
5. **C4, "79% of within-sub-catchment links jump between different tributaries".**
   UNVERIFIABLE from the repo (no tributary id in the flat files); an independent on-path
   proxy gives 86.7%. Treat as directional, not exact.
6. **C7, "r = 0.81" and "SSP3-7.0 below SSP1-2.6" as a firm correlation.** The inversion in
   the 2045 single-year means is real (SSP3 0.696 < SSP1 0.725) and it reverses on the
   2036–2045 mean. The correlation is ≈0.77 for the 2045 snapshot, not 0.81; decadal means
   correlate 0.96–0.97. The substance (axis ≠ forcing) stands; the coefficient does not.
7. **E4 / report, "21-year interaction burn-in".** Wrong: the alternative-interaction
   matrices used 21, 494 and 532 years. This error is in the report too, and the review
   inherited it.
8. **E10, "settings existed only as undocumented environment variables".** Wrong: the
   overrides are documented in `docs/climate_scenarios.md` and recorded in each run's
   `run_metadata.json`. The valid complaint is the absence of a committed climate-scenario
   configuration.
9. **E12, "99 matched structures are demolished or abandoned".** Not reproducible from the
   flat files without the MGM codebook (a generous reading gives 84–119). The matched-count
   variability (966–986) is real for tie permutations, but the single observed value is 973.
10. **E19, "runoff declines of 17–55%".** The upper decline (−55%) is real; the lower bound
    is not reproduced — some projections increase.
11. **E20, "the residual SD falls from 2.65 to 2.34".** That is the **Spain** model. The
    Guadalquivir calibration (the one relevant to the report) has SD 2.54 and contains none
    of the 257 coded values, so the quoted improvement does not apply to this basin.
12. **E20, "the flag in 06a should also be two-sided, as in 02".** `06a` already computes
    `abs(resid) > 4·MAD`; it is not one-sided. The proposed rule is also partly unsound
    (it would delete plausible cold groundwater-fed summer readings).
13. **E21, slope/coefficient framing.** The verified 0.510 / −0.163 are the **national**
    coefficients; the Guadalquivir model is 0.436 / −0.299. Lagged air improves fit for
    3- and 7-day means but not 30-day, on only 3 in-sample sites. "The model runs 0.5 °C
    too cold" is really E22; "15–30%" is the reviewer's own inference.
14. **Minor #11, "number of matched obstacles varies between 966 and 986".** Across tie
    permutations yes, but the tie changes *which edge* an obstacle is assigned to, not the
    count; the observed count is stably 973. The coordinate shift is ~234 m, not ~200 m.
15. **Minor #2, "25 cases" of juvenile-only presence.** 23 value-cells across 8 site rows;
    the count needs a stated definition.
16. **Minor #12, "ten headers mistranslated" and "863 m".** The precipitation/perimeter
    example is confirmed, but the count of ten and the 863 m coordinate discrepancy are
    not reproduced; the two coordinate sources agree for site 1.14.10.
17. **"K = 1 × observed" as a report error.** Not an error: the runs used 1× and metadata
    confirms it. It is a configuration-reproducibility issue only.
18. **Part II obstacle passability `IF/10`.** Not implementable as stated (unconfirmed
    direction/scale, 445 non-numeric rows).
19. **Several UNVERIFIABLE positives in "What stands"** (matrix symmetry/completeness,
    sum-of-squares formulas, Tsit5 settings). These should not be repeated as established
    facts without a check.

---

## 3. Net assessment

The reviewer's **central diagnosis is sound and well evidenced**: the interaction matrix is
loaded in the wrong direction and only one-way, the interaction term is not scaled by
growth, the initial state is model-generated and does not match the survey, the dispersal
graph is not the river tree, the site temperature level is a 1961–1990 sub-catchment air
climatology, the thermal niches are range midpoints, and the warming axis/exposure are not
the applied forcing. These are independently reproducible and they do affect the report's
qualitative conclusions. The review is also right that the model's own outputs contradict
several sentences in the report.

The reviewer's weaknesses are in the details and in the proposed remedies:

- **Numbers:** the C1 pair count/Σ, C4 tree length, C7 correlation, E4 burn-in years,
  E12 demolished count, E19 lower runoff bound, E20 basin model and E21 model identity are
  each wrong or unverifiable. The report's "21-year burn-in" and the review's copy of it are
  both wrong.
- **Internal contradictions:** C1's 89 pairs requires the *Lepomis* row that C1 says is
  skipped; C4's 8,640 km tree contradicts its own 8,783 km MST.
- **Remedies:** `c = −α` should be `c = 1 − α`; "iterate rows" does not fix *Lepomis*; the
  dendritic-parents snippet is not a tree; the obstacle `IF` mapping is not usable.

None of this rescues the September results. It means the review should be adopted as a
**diagnostic agenda with corrections**, not applied literally. The companion plan separates
the fixes that are safe to implement now from those that need redesign or a team decision.
