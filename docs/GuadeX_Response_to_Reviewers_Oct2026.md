# GuadeX — Response to the Reviewers

**Subject:** Reviewer report *"GuadeX: Technical Review and Proposed Corrections"*
(`docs/GuadeX_Review_and_Corrections_Sept2026.md`), reviewed at repository commit `9702169`.

**Response date:** October 2026.
**Corrected report:** `docs/FinalReportOctober2026.tex` (compiled PDF `docs/FinalReportOctober2026.pdf`),
superseding `legacy/FinalReportSeptember2026.tex`.
**Adjudication of the review:** `docs/GuadeX_Review_Assessment_Sept2026.md`.
**Implementation record:** `docs/GuadeX_Correction_Status_Sept2026.md`.
**Consolidated limitations:** `docs/GuadeX_Limitations_and_Scope_Oct2026.md`.

We thank the reviewers for a detailed, technically serious and largely accurate review. The
central diagnosis is correct and independently reproducible: the interaction matrix was loaded
with the wrong direction and only one-way, the interaction term was not scaled by growth, the
initial state was model-generated, the dispersal graph was not the river tree, the site
temperature level was a sub-catchment air climatology, the thermal niches were range midpoints,
and the warming axis and exposure metric did not reflect the applied forcing. Those defects are
real, they affected the September conclusions, and we have corrected them.

In the interest of a fair and evidence-based record, however, this response also sets out, in
detail, the places where a specific reviewer number is wrong, internally contradictory, or not
reproducible from the repository, and where a proposed remedy is technically incorrect or not
implemented as stated. We do **not** claim any fix we did not make, and we do **not** overstate
the corrected results: the corrected projections are interim, uncalibrated and reported as
relative scenario-minus-control differences, not forecasts.

---

## 1. Summary of actions

1. **The model was corrected and re-run, and the report reissued.** The corrections
   (C1–C7, E1–E22, and most minor issues) were implemented in the working tree, exercised by the
   test suite, and the full corrected experiment design was re-executed:
   **45 climate runs** (11 GCMs × 4 SSPs + a matched no-warming control), **32 obstacle-design
   runs** and **48 alternative-interaction runs**. The corrected results are reported in
   **`docs/FinalReportOctober2026.tex`** (compiled `FinalReportOctober2026.pdf`), which states
   explicitly that its values supersede the September 2026 release.
2. **Corrected configuration and reproducibility.** A corrected configuration
   (`parameters_climate_scenarios_corrected.toml`) and a legacy-reproduction configuration
   (`legacy/parameters_climate_scenarios_k1x_burnin.toml`) now exist; the full resolved configuration
   plus a model/code digest is written to each run's `run_metadata.json`. Nothing is committed:
   all correction work is in the working tree on base commit `9702169`.
3. **Scope.** Where a correction would add a new process or data layer (full ABC-SMC
   calibration, stochastic extinction, physiological niches, hydrology/drought, alien
   introductions, obstacle-`IF` standardisation), the change is either implemented as an explicit
   opt-in sensitivity or deliberately deferred and documented as a scope decision (§4).
4. **Honesty statement.** The corrected projections follow the reviewers' own permitted interim
   route: observed start plus a three-year baseline spin-up, no calibration, reported as
   scenario-minus-control. They are **conditional, relative, deterministic sensitivity results**,
   not absolute forecasts and not extinction probabilities.

**Status key used below:** *Accepted* · *Partly accepted* · *Not accepted* · *Deferred* ·
*Not verifiable*.

---

## 2. Point-by-point responses

### 2.1 Critical issues (C1–C7)

**C1 — The interaction matrix is loaded inverted and only one-way.**
*Reviewer's point:* `α[A,B]` stores the effect of the column species on the row species,
opposite to the model equation; only the lower triangle is filled, so invasive→native effects are
zero.
*Assessment:* **Accepted in substance; partly accepted on the reviewer's figures** (see §3).
*Action and evidence:* The matrix is now parsed from the cell text in the correct direction:
the **named target** is the affected species and the partner is the source, stored as
`α[target, source]`. Undirected competition is made symmetric, Spanish *"depreda"/"interfiere"*
is handled, *"No coexist"* is treated as absence of interaction, mixed-mechanism cells are
flagged ambiguous, and the missing *Lepomis* (`LG`) row **and column** are now present. The
result is exported as an explicit long-format table. Evidence: `src/data_preparation.jl`
(`parse_interaction_cell`, `build_interaction_long_table`, long-table loader),
`scripts/build_interaction_long_table.jl`, `data/BIOTIC/interaction_matrix_long.csv`
(**227 data rows**, **22 ambiguous source cells / 44 ambiguous rows**), a **24×24** loaded matrix
including `LG` and `AAL`, and `test/test_data_preparation.jl` (parsing and loader testsets). The
corrected report adopts the corrected direction in Eq. (1) and §3.2, and the corrected
interaction-structure sweep now shows only a small native→invasive redistribution (§3.4/§4.4)
rather than the September ~90 % collapse.
*Reviewer figures not adopted:* the "89 pairs / Σα = −65.5" and the "iterate over all rows fixes
*Lepomis*" claim — see §3, items 1–2.

**C2 — The interaction term is not scaled by the growth rate.**
*Reviewer's point:* Σα·N/K enters outside the `r·W` bracket, making interactions 1–3 orders of
magnitude stronger than growth.
*Assessment:* **Accepted** (the formulation defect is real); **partly accepted** for the proposed
coefficient.
*Action and evidence:* The interaction term is now **inside** the intrinsic-growth bracket as a
competitive Lotka–Volterra form, `clamp(1 − ΣN/K + interaction, −1, 2)`, with
`c_sj = 1 − α_sj` and `c_ss = 1`; the legacy outside-growth form is opt-in via
`interaction_inside_growth = false`. Evidence: `src/ode_model.jl`, and
`test/test_ode.jl` (testset *interaction form: competitive Lotka–Volterra (C2)*, 13 tests)
asserting the two-species equilibrium `K(1 − c₁₂)/(1 − c₁₂c₂₁)` and that per-capita growth
doubles with `r`. Report Eq. (1).
*Reviewer remedy not adopted:* the proposed `c = −α` is **incorrect**; the correct mapping is
`c = 1 − α` — see §3, item 3.

**C3 — The initial state after burn-in does not match the observed community.**
*Reviewer's point:* after a 373-year burn-in the 2026 state differs profoundly from the survey
(invasives gone except *Lepomis*; *Gambusia* zero; native richness 1.32 → ~2.5).
*Assessment:* **Accepted.**
*Action and evidence:* We adopted the reviewers' **own permitted interim measure**, not the full
calibration: `projection_route = "interim_observed"` starts from the observed 2006–2009 community
and integrates a fixed three-year baseline spin-up, with a **mandatory matched control**
(`control_run = true`), and reports every scenario as per-metric `scenario − control` deltas.
Evidence: `src/spin_up.jl` (`interim_observed_spin_up`), `run_climate_scenarios.jl`,
`src/climate_diagnostics.jl` (`scenario_minus_control`); tests *C3 interim projection route* and
*scenario_minus_control (C3 interim route)*. Report §3.4 (interim route) and Limitations 2. The
corrected run reproduces none of the September artefacts: invasive richness is ≈0.71 species per
site (not the flat 0.217), and the residual transient is explicit (control total biomass −35.5 %,
native biomass −40.8 %, report §3.1/§4.1). Full ABC-SMC calibration remains deferred (§4).
*Minor:* the reviewer's model native richness 2.52 is 2.51 in our recomputation (assessment §1).

**C4 — The dispersal graph does not follow the river network.**
*Reviewer's point:* within each sub-catchment, sites are chained by distance to the Guadalquivir
and outlets by elevation, producing artificial links (31,215 km vs an 8,783 km MST).
*Assessment:* **Accepted in mechanism and most numbers; partly accepted** on the "79 %" and on
the reviewer's fallback tree length (see §3, items 4–5).
*Action and evidence:* We implemented the reviewers' **immediate fallback** (on-path parent rule
plus an explicit minimum-spanning-tree root join) because the regional channel layer
`SW_Line_4C_.shp` carries no flow-direction attribute or elevation; deterministic tie-breaking
is used. Verified on real data: **774 = n − 1 edges**, **267 roots joined by 266 MST edges**,
**total 13,871.5 km** (legacy 31,206.3 km), **4 parent-elevation inversions** (legacy 238),
acyclic and deterministic. Evidence: `src/data_preparation.jl`
(`build_on_path_distance_matrix`), `test/test_on_path_graph.jl` (testset 2,647 tests). Report
§1.2/§2.3 (13,871.5 km) and Limitations 8. The preferred flow-directed GIS build remains
deferred (§4).

**C5 — The site temperature level is sub-catchment air temperature for 1961–1990.**
*Reviewer's point:* `extract_site_temperatures` takes the first column containing "TEMP"
(`TEMP_MEDIA_SC`, a 1961–1990 sub-catchment climatology), and the level and anomaly refer to
different periods.
*Assessment:* **Accepted.**
*Action and evidence:* The model level is now the corrected `guadex_tw` baseline
`tw_baseline_mean` (1986–2005), selected by **exact column name** and joined by site code with a
loud fallback; `level + anomaly` reconstructs the corrected daily series, and the heat-stress
calibration and exposure diagnostics use the same array. Evidence: `src/data_preparation.jl`
(`WATER_TEMPERATURE_COLUMN`), `src/temperature_forcing.jl`; tests *corrected site
water-temperature baseline (C5/E1)*, *model temperature equals corrected daily series (C5/E1)*
(16 tests) and *driver daily schedules reconstruct the corrected series* (18 tests). The
corrected trout-site mean is **12.83 °C** (test over the 36 field trout sites); the corrected
baseline range is 9.79–18.06 °C; the report uses a basin-mean of 15.82 °C and a mean of 13.72 °C
over the 72 baseline-established model trout sites (§3.4/§4.2).
*Reviewer details not adopted:* the README column is `Tºm_SC`, not `TEMP_MEDIA_SC`, and the
README does not literally call it "air"; one sub-catchment (`32.1`) has two values. The fix also
needed new plumbing (a site-code join across 776 vs 775 records) that the review did not mention
(assessment §1, Part II C5).

**C6 — Thermal niches are derived from distributional ranges rather than physiology.**
*Reviewer's point:* optima are range midpoints with σ = range/6; warm-adapted species are given a
19 °C optimum; an asymmetric physiological curve is needed.
*Assessment:* **Partly accepted** — the mechanism is correct, but "sixteen of the 24 species"
share "8 to 30 °C"; the actual count is **15**, and the literature claims (growth optima, CTmax)
are external and not verifiable from the repository.
*Action and evidence:* The symmetric empirical-range Gaussian is retained as an explicitly
labelled **first-order sensitivity baseline**, not a calibrated physiological niche, and this is
stated in the report (§2.2) and in Limitations 1. No physiological replacement was implemented in
this version — it requires per-species literature values (optimum for growth, CTmax/UILT) that
are not in the repository — and is deferred (§4).

**C7 — The warming axis and trout exposure do not reflect the applied forcing.**
*Reviewer's point:* `warming_end_degc` comes from a six-station basin warming curve, not the
daily per-GCM forcing actually applied; the exposure metric is a basin maximum at a troutless
lowland site.
*Assessment:* **Accepted on points 1, 2, 4, 5 and 6; partly accepted on point 3** — the quoted
`r = 0.81` is not reproducible (see §3, item 6).
*Action and evidence:* The **realised applied warming** is now persisted per run (2026–2045 and
2036–2045 means over sites) and used as the dose–response regressor; exposure is computed only
over baseline-established sites; the dose–response is GCM-aware (fixed effect GCM, no new
dependency); and endpoints are re-selected from the realised forcing
(`select_climate_endpoints`). Evidence: `src/temperature_forcing.jl`
(`realised_warming_anomaly`), `src/climate_diagnostics.jl`, `run_climate_scenarios.jl`,
`docs/climate_scenarios.md`; tests *realised warming anomaly (C7)*, *established_exposure_summary
(C7)*, *dose_response handles GCM structure (C7)*, *select_climate_endpoints / realised_forcing_axis
(C7)*. Report §3.4/§3.6, Figures 5–6 and 12; the wrong-axis slopes in the September Table 3 are
superseded. The corrected endpoints are the coolest SSP1-2.6/IITM-ESM (realised 0.208 °C) and
warmest SSP2-4.5/UKESM1-0-LL (realised 1.091 °C) runs (report §3.4/§3.5).

### 2.2 Major issues (E1–E24)

**E1 — Three different temperature series in one simulation.**
*Reviewer's point:* the equation, the heat-stress constant *k* and the exposure table use
different temperature series, differing by ~1.3 °C (2.6 °C at trout sites).
*Assessment:* **Accepted.**
*Action and evidence:* A single corrected water-temperature series is now used by the model
equation, the *k* calibration and the exposure diagnostics (see C5). Evidence: the C5/E1 testsets
and `src/temperature_forcing.jl`; the report states the shared series in §3.4. The reviewer's
offsets (basin 1.27 °C, trout 2.60 °C) reproduce exactly (assessment §1).

**E2 — The control is not a like-for-like counterfactual.**
*Reviewer's point:* control and burn-in used a smoothed ensemble median while scenarios used
per-GCM daily series, and each scenario stepped into already-realised warming.
*Assessment:* **Partly accepted** — the concern is valid, but "halves day-to-day variability" is
wrong (it removes ≈68 %), and the "0.34–0.52 °C" band was not reproduced.
*Action and evidence:* The corrected design pairs every scenario with a matched no-warming
control for the **same GCM** and reports scenario-minus-control differences. Evidence: report
§3.4/§3.5, `src/climate_diagnostics.jl` (`scenario_minus_control`), `control_run = true`. The
pre-2026 offset between the projection products and the 1986–2005 baseline persists and is now
explicitly documented (report Limitations 11); a present-day baseline window driven by each GCM's
own historical series was **not** implemented and is deferred (§4).

**E3 — The obstacle sweep starts from a single equilibrium.**
*Reviewer's point:* all 32 runs started from one burn-in at cost 0.05 / baseline passability, so
the sums of squares measure transient re-equilibration.
*Assessment:* **Partly accepted.** The specific defect (a single equilibrium tied to one
cost/passability pair) is removed by the interim route and by within-design reporting; a separate
burn-in for every combination was not run.
*Action and evidence:* Under the corrected configuration every sweep run starts from the observed
community with a three-year baseline spin-up, and the corrected sweep is analysed in matched
factorial blocks with balanced **within-sweep sums of squares** and blocked-minus-baseline
contrasts (report §3.3/§3.5; `results/sensitivity_obstacles_corrected/report_plots/`). A dedicated
per-combination burn-in remains deferred (§4).

**E4 — The burn-in convergence criterion is weak.**
*Reviewer's point:* the criterion only bounds total biomass while composition drifts, so the
"21-year" interaction burn-in is probably premature.
*Assessment:* **Partly accepted** — the weak criterion is confirmed, but the "21-year" figure is
**wrong**: the alternative-interaction matrices used **21, 494 and 532 years** (and the error is
also present in the September report; review §7, item 7).
*Action and evidence:* A robust composition stop rule was added — the 95th percentile of
`|Δ log N|` over **active** site×species cells (density floor 0.1) in at least one compared year —
combined with the basin-total rule (`spin_up_criterion = "both"`, tol 1e-2, min 10 years); the raw
all-cells statistic is retained as a diagnostic, and the actual matrix-specific burn-in years are
recorded in `results/alt_interactions_corrected/burnin_summary.csv`. Evidence: `src/spin_up.jl`
(`composition_q95_change`, `composition_diagnostics`), tests *spin_up convergence rule (E4)* (29
tests) and *burn-in diagnostics flow into run metadata and index*. See §3, item 7.

**E5 — "2045" is the state on 1 January 2045.**
*Reviewer's point:* offsets 0–19 mean only 19 years are reported and the 2026 row precedes any
forcing.
*Assessment:* **Accepted.**
*Action and evidence:* Reporting offsets are now `1:20`, so the row labelled 2045 is the
end-of-year state (`t = 7300`). Evidence: `run_climate_scenarios.jl`
(`YEAR_OFFSETS = collect(1:SIMULATION_YEARS)`), `test/test_outputs.jl` (*annual snapshot
alignment*, 30 tests); applied consistently across the three drivers. Report §3.4.

**E6 — The quasi-extinct fraction counts sites that were never occupied.**
*Reviewer's point:* with a baseline density of 0 the site is flagged, so for trout the fraction is
≈1 − occupancy (the reported 0.956 "saturation").
*Assessment:* **Accepted.**
*Action and evidence:* Quasi-extinction is evaluated only over sites where the species was
**baseline-established**; a species absent everywhere in the baseline yields fraction 0.0; the
established denominator is exposed. Evidence: `src/outputs.jl` (`baseline_established`,
`established` argument, `n_sites_established`); tests *quasi_extinction_flags never flags
non-established sites (E6)*, *quasi-extinction fraction over established sites (E6)*, *species
absent everywhere in the baseline → fraction 0 (E6)*, *trout-style fraction is not 1 − occupancy
(E6)*. Report §3.6/§4.5: the corrected trout fraction is **0.278–0.306** (mean 0.290). The
reviewer's diagnosis (fraction 0.9561 = 1 − occupancy exactly, over all 775 sites) is confirmed.

**E7 — The factor ranking is not a valid comparison.**
*Reviewer's point:* effects are level spans standardised across different run sets, and factors
are confounded.
*Assessment:* **Accepted.**
*Action and evidence:* The ranking is replaced by balanced **within-design sums of squares**
(`variance_decomposition.csv`, `parameter_effect_rank.csv`), and a global sensitivity toolkit plus
a calibration pilot (Latin-hypercube sampling with Sobol indices) was implemented. Evidence:
report §3.3/§3.5, `test/test_global_sensitivity.jl`, `results/global_sensitivity_pilot`. The
interaction matrix remains confounded with a thermal-breadth multiplier of 0.3 and a shorter
spin-up; this is now documented (report §2.5, Limitations 11). A **full** global Sobol production
run is deferred pending agreed parameter ranges (§4).

**E8 — The model does not estimate extinction probability.**
*Reviewer's point:* the model is deterministic with no replicates; the exported
`realised_richness_loss` is a realised threshold contraction, not an extinction probability.
*Assessment:* **Accepted.**
*Action and evidence:* The metric is renamed and documented as a deterministic **realised
richness loss** and is explicitly stated not to be an extinction probability; no stochasticity
was added. Evidence: report Eq. (4)/§3.6 and Limitations 5. Adding demographic/environmental
stochasticity with replicates and species-specific quasi-extinction thresholds is deferred (§4).

**E9 — Old results can be silently reused.**
*Reviewer's point:* cache keys and fingerprints omit parameter values and code version, and an
existing output file is always skipped.
*Assessment:* **Accepted.**
*Action and evidence:* A parameter **and code digest** is folded into burn-in cache keys and
completed-run fingerprints; an existing run is reused **only when its stored digest matches**;
`run_metadata.json`/`run_digest.txt` sidecars are written. Evidence: `src/run_identity.jl`,
`sensitivity_core.jl`, `run_climate_scenarios.jl`; tests *parameter digest*, *burn-in cache key
includes the parameter digest*, *run fingerprint includes the parameter digest*. `GUADEX_CLIMATE_FORCE=1`
recomputes regardless. Report §3.4.

**E10 — `parameters.toml` does not reproduce the report.**
*Reviewer's point:* the TOML sets different settings from those actually used, which existed only
as undocumented environment variables.
*Assessment:* **Partly accepted** — the TOML genuinely did not reproduce the report, but the
environment overrides were **documented**, not undocumented (see §3, item 8). The real gap was the
absence of a committed scenario configuration.
*Action and evidence:* Two configurations now exist — the corrected run config
(`parameters_climate_scenarios_corrected.toml`) and the legacy-reproduction config
(`legacy/parameters_climate_scenarios_k1x_burnin.toml`) — and the full resolved settings plus digest are
written to `run_metadata.json`. Evidence: those files, `docs/climate_scenarios.md`, and the run
metadata/`run_digest.txt` sidecars.

**E11 — The habitat index is not a trophic index.**
*Reviewer's point:* `IET` is the bank-stability index (*índice de estabilidad del talud*), not a
trophic-state index; one outlier (27.67) sets the scale.
*Assessment:* **Accepted.**
*Action and evidence:* `IET` is relabelled the **bank-stability** index in the code docstring and
in the report, which now describes it as a rescaled bank-stability proxy held constant (§2.2).
Evidence: `src/data_preparation.jl` docstring; report §2.2; report and decision brief text
corrected.

**E12 — Obstacles are matched poorly and given uniform passability.**
*Reviewer's point:* obstacles fall on artificial links; all receive 0.1/0.5; multiple obstacles do
not accumulate; the legacy dam layer already restricts ~52 % of links; 99 matched structures are
demolished or abandoned.
*Assessment:* **Partly accepted** — the uniform 0.1/0.5 and the 51.7 % legacy restriction are
confirmed; "99 demolished" and "797/968" are not reproducible from the flat files (see §3,
item 9).
*Action and evidence:* Obstacles are snapped to the corrected tree; passabilities **multiply**
along a link (0.1·0.1 = 0.01) instead of taking the element-wise minimum; non-operational and
ambiguous structures are classified and excluded/annotated; legacy de-duplication and
deterministic assignment are applied. Evidence: `src/data_preparation.jl`
(`build_obstacle_passability_matrix`, status classification), `test/test_data_preparation.jl`
(`n_status_active = 1274`, `n_status_nonoperational = 185`, `n_status_ambiguous = 199`),
`test/test_on_path_graph.jl` overlay counts. Verified overlay: matched **973/1658** (legacy) or
**1036/1658** (on-path), applied **727/796**, excluded **387**. Report §1.2/§2.3 (1,036 matched
within 2,000 m). The `IF` index is deliberately **not** used (deferred; §4).

**E13 — The 2006–2009 presences fix species distributions (10 % absence factor).**
*Reviewer's point:* absences receive 10 % of *r* and 39 % of sites had no fish, preventing
expansion by construction.
*Assessment:* **Accepted** (the 10 % factor and 38.9 % empty-site fraction are confirmed; the
"prevents expansion by construction" wording is overstated, since 0.1·r is still positive).
*Action and evidence:* `absence_growth_fraction` is now an explicit option; the corrected run sets
it to **1.0** so temperature and habitat decide establishment. Evidence:
`parameters.toml`/`parameters_climate_scenarios_corrected.toml` `[biological_options]`,
`SimulationParameters.biological_options`, and `test/test_data_preparation.jl` *Explicit
biological options (E13–E16, E18)* (46 tests). Report §3.2.

**E14 — Carrying capacity is driven by isolated pools.**
*Reviewer's point:* 84 pool sites hold 56 % of density; a 30 m² pool contains 9,533 units of
*Gambusia*; sampling was a single visit in different seasons.
*Assessment:* **Accepted.**
*Action and evidence:* `pool_capacity_mode` is an explicit option; the corrected run uses
`"cap"`, clipping the 84 isolated-pool capacities at the **non-pool median**. Evidence: tests
above; report §3.2. The non-pool-median choice is arbitrary but documented (status §5.3).

**E15 — Fishless sites are treated as habitat that can be colonised.**
*Reviewer's point:* of ~294 fishless sites, 56 were flagged `SIN_PECES = DIFÍCIL`.
*Assessment:* **Accepted** (the fishless count is 293, not 294).
*Action and evidence:* `fishless_dificil_capacity` is an explicit option; the corrected run uses
`"near_zero"` (`K = 1e-6`) at the 56 `SIN_PECES = DIFÍCIL` sites, making them transit-only.
Evidence: tests above; report §3.2.

**E16 — Species that do not reproduce in the basin are modelled as local populations.**
*Reviewer's point:* eel and both mullets have `REPRODU_WITHIN_THE_BASIN = no` but receive a local
*r*; four of nine eel records are fish-farm escapees in the Guadiato.
*Assessment:* **Accepted.**
*Action and evidence:* `nonreproducing_local_growth = "zero"` is applied to `AA`, `LR` and `MC`,
and `exclude_fishfarm_eel_records = true` drops the four documented Guadiato fish-farm eel records
(AA presences 9 → 5). Evidence: tests above; report §3.2.

**E17 — Growth and dispersal rates lack reliable sources.**
*Reviewer's point:* the cited sources include a kit-fox mange study (for *Gambusia* growth) and a
channel-catfish study (for its dispersal); sedentary species receive 6–25 km/yr; *Alburnus* is
missing from both dictionaries.
*Assessment:* **Partly accepted.** The source-quality and *Alburnus* concerns are valid; the main
rate-table rebuild is not implemented in this version.
*Action and evidence:* The corrected 24-species matrix now includes `AAL`, and the report states
plainly that *Alburnus alburnus* (code `AAL`) is an exotic species still **misclassified as
migratory**, with missing interaction and dispersal trait data and a placeholder dispersal scalar,
so invasive richness and biomass are reported as **lower bounds** (report §2.1). A single
source-annotated `data/species_parameters.csv` rebuilding all growth and dispersal rates is
deferred (§4).

**E18 — Available habitat filters are not used.**
*Reviewer's point:* salinity (brackish for AB, LR, MC) and elevation ranges are read but not
applied; brown trout is eligible at 477 sub-500 m sites.
*Assessment:* **Accepted** (477 sub-500 m trout-eligible sites confirmed exactly).
*Action and evidence:* `salinity_envelope` is an explicit option; the corrected run applies it
(1,584 cells of strictly brackish AB/LR/MC below the conductivity threshold zeroed). An
**elevation** envelope is deliberately **not** offered, because the review itself notes it would
block upslope range shifts. Evidence: tests above; report §3.2.

**E19 — Hydrology is absent.**
*Reviewer's point:* the 261 dry or effluent-only sites are removed, so fish cross them at no cost,
and the loaded CEDEX 2045 runoff files (declines of 17–55 %) are unused.
*Assessment:* **Partly accepted** — the 261 excluded sites and the unused CEDEX files are
confirmed; the "dry or effluent-only" attribution is not verifiable from the repository, and the
lower runoff bound ("17 %") is **not reproduced** (some projections increase; see §3, item 10).
*Action and evidence:* No hydrological driver was added. This is recorded as a deliberate scope
limitation (Limitations 6: no discharge, flow intermittency, drought or abstraction), and the
runoff projections remain unused. A hydrology/drought module is deferred (§4).

**E20 — Coded values contaminate the water-temperature calibration.**
*Reviewer's point:* 257 Waterbase observations are exactly 1/6 or 1/3 °C, they pass quality
control, and the outlier flag is never applied; removing them (and refitting) lowers the residual
SD from 2.65 to 2.34 °C.
*Assessment:* **Accepted that the coded values are real** — but the cited SD improvement is the
**Spain** model, not the Guadalquivir calibration (see §3, item 11), and the claim that flag
`06a` is one-sided is **wrong** (`resid` is already defined as an absolute value).
*Action and evidence / not changed:* The coded-value filter and refit were **not** implemented in
this version. The corrected pipeline regenerated the temperature products with the E22 bias
correction only; the same-day calibration still uses the Waterbase input. We did not claim this
fix, because adopting it well requires a physical-consistency rule and a decision on the affected
cold-groundwater readings. It is listed as deferred (§4).

**E21 — The sensitivity of water to air warming is probably underestimated.**
*Reviewer's point:* the slope is from same-day air temperature (0.51); 7- and 30-day means give
0.60 and 0.64 and fit better; the size of the bias is the reviewer's own inference (15–30 %).
*Assessment:* **Partly accepted.** The like-for-like concern and the lagged slopes are valid, but
the verified 0.510 / −0.163 coefficients are the **national** model (the Guadalquivir model is
0.436 / −0.299); the lagged fit is better only for 3- and 7-day means, not 30-day; and only ~3
in-sample lag-eligible sites exist (see §3, item 12).
*Action and evidence:* The same-day calibration is retained as primary, and the lagged
relationship is reported as an **uncertainty band** on projected warming. Evidence:
`guadex_tw/outputs/tables/air_water_slope_uncertainty_band.csv` (ratios ×0.634 to ×1.276, i.e.
×0.63–×1.28); report §4.1 (*E21 warming uncertainty band*) and Limitations 3. A lagged refit
cannot be a primary calibration on three sites and is deferred (§4).

**E22 — Cold bias in the Guadalquivir is not corrected in the series used by the fish model.**
*Reviewer's point:* the held-out-basin bias is −1.1 °C (−2.2 °C in February–March), and script 13
omits the bias and elevation correction that script 10 applies.
*Assessment:* **Accepted.**
*Action and evidence:* A per-month held-out-basin (LOBO) bias correction was applied in the
`guadex_tw` pipeline, and the 48 model-facing temperature products were rewritten. Evidence:
`guadex_tw/scripts/13b_apply_bias_correction.py`,
`guadex_tw/models/lobo_monthly_bias_correction.json`,
`guadex_tw/logs/bias_correction_13b_run.log`; verified pooled annual **+1.1059 °C**, Feb
**+2.2237 °C**, Mar **+2.2867 °C**, 48 files processed. Anomalies are unchanged by design, as the
review anticipated. Report §3.4.

**E23 — Limited empirical support and extrapolation in elevation.**
*Reviewer's point:* 64 % of Guadalquivir observations are from 2024, two sites have one
observation, and 17 model sites lie above the highest calibration site (1,127 m).
*Assessment:* **Accepted.**
*Action and evidence / not changed:* The high-elevation extraction requires regional climate
projections not present in the repository, so the elevation ceiling was **not** completed; the
corrected baseline is documented (range 9.79–18.06 °C, report §3.4) and the extrapolation
remains a stated limitation. Deferred (§4).

**E24 — The report misdescribes the forcing.**
*Reviewer's point:* §3.4 stated all sites received the same anomaly, but it varies by site
(0.82–1.08 °C in one run, r = −0.79 with elevation); the elevation scaling would use a hard-coded
reference of 0 m.
*Assessment:* **Accepted in substance** — the "all sites same anomaly" text was indeed wrong; the
exact 0.82–1.08 range and r = −0.79 were not reproduced.
*Action and evidence:* The report now describes per-site daily anomalies and the persisted
realised warming (report §3.4/§4.1); the corrected runs use the delivered per-site daily series,
so the legacy `warming_matrix` helper's default `z_ref = 0.0` is not used by the corrected route.
Evidence: `src/temperature_forcing.jl`; report §3.4.

### 2.3 Minor issues

1. **Rare-species site counts.** *Accepted* (exact). No model change; the E6 correction removes the
   denominator artefact, but statistics for species on one or a few sites remain fragile and are
   noted as a limitation.
2. **Juvenile-only presences.** *Partly accepted.* The absence of a juvenile column for IO is
   confirmed; the "25 cases" count is not reproduced from the flat files (**23 value-cells across
   8 site rows**; a stated definition is needed). Not changed in this version.
3. **Row-order alignment.** *Accepted.* The initial state is now aligned by **site code** rather
   than row position. Evidence: `run_climate_scenarios.jl`, `scripts/serialize_data.jl` (imports
   `Guadex` and calls the shared extractor), and the alignment tests. The named function
   `initial_state_from_density` does not exist; the fragility was in three run scripts
   (assessment §1).
4. **Rate conversion.** *Partly accepted.* `r/365` is the default; a documented **opt-in**
   seasonal conversion `log1p(r_seasonal)/183` now exists (matching the documented *Gambusia*
   formula). Evidence: `src/data_preparation.jl` (`seasonal_rates`). The option is code-only and is
   not yet exposed through `parameters.toml` (status §5.8).
5. **Calibration of *k*.** *Partly accepted.* *k* is now recalibrated from the corrected baseline
   so that the most exposed species–site combination under baseline climate loses no more than
   5 % of annual survival (report §3.2); the reviewer's "5.9 %" was not independently reproduced.
6. **Calendar.** *Accepted.* The corrected pipeline uses a single explicit 365-day-year
   convention (`DAYS_PER_YEAR = 365`, end-of-year at `t = 7300`), which removes the day-366
   climatology drop; the leap-day handling is documented.
7. **Monthly steps.** *Partly accepted.* The harmonic alternative was fitted and compared but the
   monthly fixed effects are retained because the AIC is equivalent (3654.5 vs 3654.9); the
   reviewer's "1.55 °C step" was not reproduced. Not changed.
8. **Hot tail.** *Partly accepted.* Curvature exists but is smaller than stated (≈ +0.09 °C over
   25–30 °C and −0.43 °C over 30–50 °C, not "≈0.7 °C above 25 °C"); a logistic form was not
   adopted. Deferred.
9. **Mixed-model details.** *Accepted.* The mixed-model LRT is refitted with ML (`reml=False`) and
   the invalid REML comparison is kept for audit only. Evidence: `guadex_tw/scripts/08_models.py`.
   The REML LRT across fixed effects, the ICC at 0 °C and the in-sample models in the CV table are
   confirmed; the quoted ML statistics are not persisted (status §6) and "strengthens with ML"
   remains unverified in the tree.
10. **Numerical details.** *Accepted.* The clipping is now documented and exercised (report
    Eq. (1)); mass-conservation, topology and sign tests were added (see C1/C4/E1 tests). The
    unclipped-immigration/clipped-emigration asymmetry and the positivity-callback mass are
    documented residuals.
11. **Obstacle count and coordinates.** *Partly accepted.* The ED50/ETRS89 split is confirmed
    (site 1.14.10 shifts **~234 m**, not ~200 m), but the "966–986" variation is **not
    reproduced**: tie order changes *which edge* an obstacle is assigned to, not the count, and
    the observed matched count is stably **973** (see §3, item 13).
12. **Auxiliary data.** *Partly accepted.* The `PERIMETRO_MOJADO`/precipitation example and the
    bogus "SSP285" scenario are confirmed; the "ten headers" count and the "863 m" coordinate
    discrepancy are **not** reproduced (the two coordinate sources agree for site 1.14.10). The
    unused columns were not rebuilt in this version.
13. **Code and tests.** *Accepted.* Vacuous tests were replaced and mass-conservation, topology,
    α-sign and temperature-level tests were added. Evidence: `test/test_ode.jl`,
    `test/test_on_path_graph.jl`, `test/test_data_preparation.jl`,
    `test/test_temperature_forcing.jl`. Residual: the dead `SPECIES_CODES` constant in
    `scripts/serialize_data.jl` still duplicates `"AA"` and omits `"AAL"`; the live species list
    is resolved elsewhere and does include `AAL` (status §5.7).
14. **Viewer.** *Accepted.* Species-level time-series export and aggregate rendering were added
    (e.g. `results/.../export/viewer/species_<sp>_*_timeseries.json`), and the site count is
    corrected to **775** in the report. Remaining viewer cosmetic issues (mixed metric/year
    slider) are noted.

### 2.4 Statements in the reports that would benefit from revision

All rows of the review's statements table were reviewed and, where they were correct, acted on:

| Statement (reviewer) | Assessment | Action / corrected reading |
|---|---|---|
| 775 sites grouped into 774 water bodies | Accepted | Report now reports **289** water-body reporting units (288 assigned + 1 unassigned) and **774 network edges**, explicitly not a water-body count (§1.2/§2.1). |
| *Alburnus* coded "Aa" and classed migratory | Partly accepted | The report uses full species names and contains no "Aa"; the real, acknowledged defect is the **misclassification as migratory** (`AAL`), now stated as a lower-bound caveat (§2.1). |
| Daily absolute water temperature = sub-catchment climatology | Accepted | Report describes the corrected `guadex_tw` water-temperature series (`tw_baseline_mean` + anomaly) (§3.4). |
| All sites receive the same basin-scale anomaly | Accepted | Text corrected to **per-site daily anomalies**; realised per-run warming is tabulated (§3.4/§4.1). |
| Habitat index from the trophic-state index | Accepted | *IET* is now the **bank-stability** index (§2.2). |
| Basin mean 16.0 °C; trout 1.5 σ warm | Accepted | Corrected basin-mean **15.82 °C**; trout **0.65 σ** warm of optimum (§3.4/§4.2). |
| Max days > 20 °C rising 131→139–174 | Accepted | Replaced by exposure at the **72 baseline-established trout sites** (baseline max 90; 2045 80–135; mean 23→45 days) (§4.2, Fig. 12). |
| Basin-mean warming 0.74–0.97 °C | Accepted | Replaced by **realised** warming: pathway means 0.54/0.64/0.62/0.70 °C; run span 0.21–1.09 °C (§4.1). |
| Invasive richness 0.2168 because no invader crosses threshold | Accepted | Corrected value **≈0.71** species per site; the flat September value was a C1/C3 artefact (§4.5). |
| No warming-driven invasive expansion detected | Accepted | The corrected run supports the *statement* (invasive richness essentially invariant), but the September *reason* was an artefact; now reported with its caveat (§4.5). |
| Blocking helps by restricting competitively effective invasives | Accepted | Statement retired: blocking raised basin native biomass monotonically; the stated mechanism is unsupported (§4.3). |
| Invasive-favouring matrix is a deliberate extreme | Accepted | Direction now correct; the matrix remains a **bounding scenario** with a small contrast (−1.6 % native, +4.0 % invasive) (§4.4). |
| 3-D dendritic network connecting adjacent sites | Accepted | Replaced by the corrected **on-path river-network tree** (§1.2). |
| Migratory taxa have dispersal > 5× median | Accepted | Corrected: only MC and LR exceed the median; AA and AAL are about **0.5×** (§2.2). |
| Final-year (2045) metrics | Accepted | Now **end-of-year** (E5) (§3.4). |
| Trout quasi-extinction near saturation | Accepted | Corrected fraction **0.278–0.306**; the 0.94–0.97 values were a denominator artefact (E6) (§4.5). |
| Carrying capacity equal to 1× observed | **Not accepted as a report error** | The runs *did* use 1×; this is a configuration-reproducibility issue (E10), not a report error (§3, item 16). |
| 1,037 sites; 12 native / 11 exotic | Accepted | Report uses **775 sites**; **10 native, 10 invasive, 4 migratory** (§2.1). |
| Observed state used as starting point | Accepted | Survey dates (2006–2009) and the interim observed start are stated (§3.4). |

### 2.5 Part II proposed fixes

| Proposed fix | Assessment | Outcome |
|---|---|---|
| C2 Lotka–Volterra rewrite `c = 1 − α` | Adopted | Implemented and tested; the earlier `c = −α` is corrected (§2.1 C2). |
| E5 offsets 1…20 | Adopted | Applied in all three drivers (§2.2 E5). |
| E9 parameter digest | Adopted | Cache keys and fingerprints (§2.2 E9). |
| E22 monthly bias, script 13 | Adopted | Implemented in `13b` (§2.2 E22). |
| E10/E11 IET rename and report text | Adopted | Bank-stability label (E11). |
| Mixed-model LRT `reml=False` | Adopted | `08_models.py` (§2.3 #9). |
| C1 direction parsing | Partly adopted | Text parsing and long table implemented, but the reviewer's regex has two bugs and its "iterate over all rows" does **not** fix *Lepomis* (no column); we added the column explicitly (§3, items 1–2). |
| C4 `dendritic_parents` fallback | Partly adopted | The on-path rule was kept but required an explicit MST root join; the reviewer's snippet yields 267 roots and its ~8,640 km is the parent-edge sum only; the full tree is 13,871 km (§3, item 4). |
| C5 water-temperature level | Partly adopted | Conceptually right; implemented with new plumbing (site-code join 776 vs 775) and the E22 correction applied consistently, none of which the review mentioned (§2.1 C5). |
| C7 realised forcing / exposure at occupied sites | Partly adopted | Persisted realised-anomaly field plus a join to baseline-established cells implemented (§2.1 C7). |
| E6 established-only quasi-extinction | Adopted | `baseline_established`/denominator argument added (§2.2 E6). |
| Obstacles from the inventory `IF` index | **Not adopted** | `IF` semantics are unconfirmed (max 7.97; **547/1658** non-numeric rows in the repository, vs the review's 445); `IF/10` cannot be implemented safely (§3, item 17). |
| River-line graph rebuild | Partly adopted | `SW_Line_4C_.shp` has no flow-direction attribute or elevation; the on-path fallback was used and the preferred flow-directed build is deferred (§4). |
| Per-GCM controls, Sobol, input rebuilds | Partly adopted | Matched per-design control is in the corrected design; global Sobol toolkit + pilot implemented, production run deferred; E14/E17/E18 options implemented, E17 rate-table rebuild, E19 hydrology and E21 lagged refit deferred (§4). |

---

## 3. Corrections we did not adopt because they were not accurate

These are the specific claims we judged incorrect, internally contradictory or not reproducible
from the repository. For each we give the correct value or reading and its evidence. The full
adjudication is in `docs/GuadeX_Review_Assessment_Sept2026.md` §2.

1. **C1 — "89 pairs, Σα = −65.5".** The loaded matrix (which skips *Lepomis*, as C1 itself says)
   has **79** non-zero invasive→native effects summing to **−58.1**. The 89 / −65.7 figures
   require the *Lepomis* row that C1 says is skipped, so the review contradicts itself. The
   headline claim — invasive→native effects are zero in all 100 pairs — is exact.
2. **C1 — "iterate over all rows" fixes *Lepomis*.** False: the CSV has an `Lg` **row** but no
   `Lg` **column**, and the guard tests membership in the column list, so the row is skipped
   regardless. A column (or a long-format table) is required; we added both.
3. **C2 — `c = −α`.** The correct competitive Lotka–Volterra mapping is **`c = 1 − α`** (with
   `c_ss = 1`). The review's own verification identity, `K(1 − c₁₂)/(1 − c₁₂c₂₁)`, holds only for
   the corrected mapping. Evidence: `test/test_ode.jl`.
4. **C4 — the prototyped tree "of about 8,640 km".** The parent edges sum to 8,637 km, but the
   **full tree is 13,871 km** once the 267 roots are joined. 31,215 km and the 8,783 km MST
   reproduce, but a spanning tree cannot be lighter than the MST, so "close to the MST" is
   internally inconsistent. Evidence: report §1.2 (13,871.5 km); `test/test_on_path_graph.jl`.
5. **C4 — "79 % of within-sub-catchment links jump between different tributaries".**
   Not verifiable from the repository (no tributary id in the flat files); an independent on-path
   proxy gives **86.7 %**. Treat as directional, not exact.
6. **C7 — "r = 0.81".** Not reproducible: the 2045 single-year snapshot correlates at ≈**0.77**,
   and decadal means at 0.96–0.97. The SSP3-7.0 < SSP1-2.6 inversion is real in the 2045 single-year
   means but reverses on the 2036–2045 mean. The substance (axis ≠ forcing) stands; the coefficient
   does not.
7. **E4 — "21-year interaction burn-in".** Wrong (and copied into the September report): the
   alternative-interaction matrices used **21, 494 and 532 years**. Evidence:
   `results/alt_interactions_corrected/burnin_summary.csv`.
8. **E10 — "settings existed only as undocumented environment variables".** Wrong: the overrides
   are documented in `docs/climate_scenarios.md` and recorded in each run's `run_metadata.json`.
   The valid complaint was the absence of a committed scenario configuration, which is now
   provided.
9. **E12 — "99 matched structures are demolished or abandoned".** Not reproducible from the flat
   files without the manager's codebook (a generous reading gives 84–119). The matched-count
   "variation between 966 and 986" is real for tie permutations, but the tie changes *which edge*
   an obstacle is assigned to, not the count; the observed value is stably **973**.
10. **E19 — "runoff declines of 17–55 %".** The upper decline (−55 %) is real; the **lower bound is
    not reproduced**, because some projections increase. The "dry or effluent-only" attribution
    is also not verifiable from the repository.
11. **E20 — "the residual SD falls from 2.65 to 2.34 °C".** That is the **Spain** model. The
    Guadalquivir calibration relevant to the report has SD **2.54** and contains none of the 257
    coded values, so the quoted improvement does not apply to this basin. In addition, the claim
    that flag `06a` is one-sided is wrong: `resid` is defined as `(tw_obs − med).abs()` before the
    `resid > 4·MAD` test (`guadex_tw/scripts/06a_qc_waterbase_direct.py`, lines 84–86), so the
    flag is already two-sided.
12. **E21 — model identity.** The verified **0.510 / −0.163** coefficients are the **national**
    model; the Guadalquivir model is **0.436 / −0.299**. The lagged fit is better only for 3- and
    7-day means (not 30-day), on only ~3 in-sample lag-eligible sites; "the model runs 0.5 °C too
    cold" is really E22, and the "15–30 %" size is the reviewer's own inference, not a measured
    quantity.
13. **Minor #11 — "the number of matched obstacles varies between 966 and 986".** Across tie
    permutations this can occur, but the tie changes edge assignment, not the count; the observed
    count is stably **973**. The coordinate shift is **~234 m**, not ~200 m.
14. **Minor #2 — "25 cases" of juvenile-only presence.** Our count from the flat files is **23
    value-cells across 8 site rows**; the number needs a stated definition.
15. **Minor #12 — "ten headers mistranslated" and "863 m".** The precipitation/perimeter example
    and the "SSP285" scenario are confirmed, but the count of ten and the 863 m coordinate
    discrepancy are not reproduced; the two coordinate sources agree for site 1.14.10.
16. **"K = 1 × observed" as a report error.** Not an error: the September climate runs really did
    use `carrying_capacity_base_scaling = 1` (confirmed in run metadata; mean K ≈ 86.6 vs observed
    ≈ 84.4). This is a configuration-reproducibility issue, now fixed by committing the configs
    (E10), not a false statement in the report.
17. **Part II obstacle passability `IF/10`.** Not implementable as stated: the `IF` direction and
    scale are unconfirmed (max 7.97), and the repository reports **547/1658** non-numeric rows
    (na/DE/SD), not 445. Decoding `IF` prematurely would silently mis-scale every barrier, so the
    uniform 0.1/0.5 assumption is retained and documented (Limitations 4).
18. **Unverifiable positives in "What stands".** The matrix symmetry/completeness, the
    sum-of-squares formulas and the Tsit5 settings were marked UNVERIFIABLE in the assessment;
    we do not repeat them as established facts.

Additional small numeric corrections, for completeness: C3 native richness is **2.51** (not 2.52);
C6 **15** species share "8–30 °C" (not 16); E2 removes **~68 %** of day-to-day variability (not
half); E15 gives **293** fishless sites (not 294); C1 has **13** row-named cells (not 7) and the
"53 competition" entries are **48–52**.

---

## 4. Changes deferred or rejected as out of scope for this version

Each of the following is a deliberate scope decision for this iteration, not an oversight. The
deferral reasons and consequences are recorded in `docs/GuadeX_Correction_Status_Sept2026.md` §4
and `docs/GuadeX_Limitations_and_Scope_Oct2026.md`.

- **Full ABC-SMC calibration (C3).** We used the interim route the review explicitly permits
  (observed start + short spin-up + scenario-minus-control). Simulation-based inference on *r*,
  *K*, interactions and dispersal against pre-registered TSS/occupancy criteria, with spatial
  block cross-validation and posterior predictive checks, needs a computing campaign and
  pre-registered acceptance criteria. Consequence: the corrected results are **provisional and
  uncalibrated** and are labelled as such.
- **Stochasticity and extinction probability (E8).** Rather than add demographic/environmental
  stochasticity and replicates, we renamed and re-documented the metric as a deterministic
  **realised richness loss**. Adding replicates multiplies the run cost; the naming fix removes the
  over-claim without it.
- **Physiological thermal niches (C6).** The symmetric empirical-range Gaussian is retained as an
  explicitly labelled first-order sensitivity baseline. Replacement requires per-species
  literature values (optimum-for-growth, CTmax/UILT) and an asymmetric performance curve that are
  not in the repository; the approximation has little leverage under the mild, cool-baseline
  forcing but matters under strong warming.
- **Lagged air–water refit (E21).** A lagged refit cannot be a primary calibration on ~3
  in-sample lag-eligible sites; the lagged response is instead expressed as a **×0.63–×1.28
  uncertainty band** and treated as a conservative sensitivity envelope.
- **Global Sobol production run (E7).** The toolkit and a pilot were implemented and the corrected
  ranking uses balanced within-design sums of squares; a full global design requires **agreed
  parameter ranges and dependencies**, without which the indices would be uninterpretable.
- **Hydrology, drought and abstraction (E19).** No dynamic hydrological driver was added. The 261
  excluded reaches are not represented; the CEDEX runoff projections remain unused. This is
  documented as a standing limitation (Limitations 6).
- **Alien-species introduction scenarios.** The model has no explicit introduction events; which
  species to add and where/when are **team decisions**, and adding them would change
  invasive-richness projections.
- **Obstacle `IF` standardisation (E12).** The uniform structural assumption (0.1 upstream /
  0.5 downstream) is retained and documented because the `IF` index cannot yet be decoded safely.
- **E20 coded-value filter and refit.** Not adopted in this version: it needs a
  physical-consistency rule and affects plausible cold groundwater-fed summer readings; the
  correction was applied only as the E22 basin bias.
- **E23 high-elevation extraction.** Requires regional climate projections not present in the
  repository; the elevation ceiling remains a documented limitation.
- **E17 single source-annotated rate table, E14 seasonal standardisation, E2 per-GCM
  historical-series present-day baseline, and a per-combination obstacle burn-in (E3).** These need
  data or decisions beyond this iteration; the corrected runs use explicit, opt-in options where
  possible and matched controls where not.
- **Preferred flow-directed GIS network build (C4).** The regional channel layer carries no
  flow-direction attribute or elevation, so the validated on-path fallback is used; the
  flow-directed build is future work.
- **Elevation envelope (E18).** Deliberately **not** offered, because it would block upslope range
  shifts — the reviewers themselves make this point.

---

## 5. Verification

- **Test suite.** The full Julia suite passes: **3,693 / 3,693 tests, 67 testsets, 0 failures**
  (re-run 2026-09-30). It includes the C4 on-path graph testset (2,647 tests), the C2
  interaction-form testset (13 tests), the E4 spin-up-rule testset (29 tests), the C5/E1
  temperature-reconstruction testsets (16 + 18 tests), the biological-option testset (46 tests),
  and the new C3 interim-route testsets, together with mass-conservation, topology, α-sign and
  temperature-level tests.
- **Corrected re-runs.** All three corrected experiments were executed and are present:
  **45 climate rows** (`results/climate_scenarios_corrected/runs_index.csv`: 44 GCM×SSP scenarios
  plus one no-warming control), **32 obstacle rows**
  (`results/sensitivity_obstacles_corrected/runs_index.csv`), and **48 interaction rows**
  (`results/alt_interactions_corrected/runs_index.csv`). The design is summarised in report
  Table 1 (§3.5).
- **E21 uncertainty band.** `guadex_tw/outputs/tables/air_water_slope_uncertainty_band.csv`
  tabulates central/low/high warming by pathway; the slope ratios are **0.6344–1.2765**
  (×0.63–×1.28). For the GuadeX 2026–2045 horizon the central/low/high ensemble-median water
  warming is 0.580/0.368/0.741 °C (SSP1-2.6) to 0.730/0.463/0.931 °C (SSP5-8.5); the report notes
  the upper end may not apply to the Guadalquivir subset (§4.1, Limitations 3).
- **Limitations document.** `docs/GuadeX_Limitations_and_Scope_Oct2026.md` consolidates the
  standing caveats into **11 numbered limitations and a scope statement**, and records the C6, C3,
  E21, `IF` and E8 decisions. The report's Limitations section (§6) is drawn from it.
- **Invariants.** `level + anomaly` reconstructs the corrected daily series (tests assert to
  `atol = 1e-9`, exactly for the driver schedules); the corrected trout-site mean is 12.83 °C;
  the on-path tree is acyclic, deterministic and has 774 edges / 13,871.5 km / 4 parent
  inversions; obstacle passabilities multiply; a mismatched digest forces recomputation; and the
  default biological options reproduce the historical arrays byte-for-byte.
- **Not persisted / unverified.** The ML (`reml=False`) LRT statistics are computed but not
  committed; the smoke-run wall time and the 8.9e-16 residual are not recorded; the E4 real-data
  convergence year (~423) is asserted in comments and exercised only with synthetic data.

---

## 6. Closing statement

The reviewers identified genuine, consequential defects, and the corrected model confirms that the
September conclusions changed materially: the passability sign reversal is not recovered, the
invasive-favouring contrast is small rather than catastrophic, the flat "invasive richness 0.2168"
becomes ≈0.71, and the aggregate climate response is weak and dominated by a transient common to
the no-warming control. We are grateful for the diagnosis and have adopted its substance.

We have also been explicit, as a matter of record, about the places where the review's specific
numbers and remedies were wrong, and about the corrections we deliberately deferred. The corrected
results remain **interim and uncalibrated**: they are conditional, relative, deterministic
sensitivity results bounded by the interim projection route, the symmetric distributional thermal
niches, the uniform barrier-passability assumption, the on-path network reconstruction and the
same-day air–water calibration. They should be read as internally consistent
scenario-minus-control differences, not absolute forecasts or extinction probabilities. We would
welcome the reviewers' further comments on the corrected release and on the deferred items.

*Prepared on base commit `9702169`; all correction work, the corrected configuration, the new
tests, the corrected outputs and this response are in the working tree only.*
