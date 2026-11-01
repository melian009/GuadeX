# Limitations and Scope of Inference

**GuadeX — October 2026 report.** This note consolidates the standing limitations of the
current model and its projections into a single numbered statement, written to be read as the
Limitations section of the forthcoming October 2026 report. It supersedes the scattered
caveats in `docs/climate_scenarios.md`, `docs/GuadeX_Correction_Status_Sept2026.md` and
`docs/GuadeX_Correction_Plan_Sept2026.md` for the items below, and it records how the team's
**C6, C3, E21, IF and E8** decisions are represented. No model code or `.tex` report is changed
by this document.

**Phase 2 update (2026-11-01).** The post-parser-fix experiments have been re-run, the validation
has been repeated against the Phase 2 control, and the competition-structure and
passability-floor/on-path-only sensitivities are now reported. Items 13–15 are new; items 2, 3
and 8 have been updated. The primary post-fix outputs are the `results/*_phase2` directories;
the `results/*_corrected` directories remain the pre-fix audit baseline.

---

## Limitations

1. **Thermal niches are distributional, not physiological (C6).** The thermal response of each
   species is represented by a symmetric Gaussian whose optimum is the midpoint of the reported
   thermal range and whose breadth (σ) is fixed at one sixth of that range. This construction
   is a first-order sensitivity baseline, as the review explicitly accepts, and it is retained
   here rather than treated as a calibrated physiological niche. Because every baseline site in
   the basin lies at or below 18.06 °C, and because the warming delivered over the simulation
   horizon is mild (approximately 0.6 °C in the central projections), the bluntness of the
   symmetric curves has little leverage on the present thermal results. That insensitivity is
   not a general property of the model: it is a consequence of the mild forcing and the cool
   baseline. Under stronger warming, and for warm-adapted species whose true growth optima lie
   above the range midpoint, the symmetric-Gaussian approximation becomes a materially limiting
   assumption. Replacement of the range midpoint with literature-based optimum-for-growth and
   critical thermal maximum (CTmax / UILT) values, and replacement of the symmetric Gaussian
   with an asymmetric physiological performance curve (for example the Sharpe–Schoolfield
   formulation), are therefore reserved as future work and are required before the thermal
   results can be interpreted under strong-warming scenarios.

2. **The projections are interim and uncalibrated (C3).** The projections follow the review's
   interim route: they are initialised from the observed 2006–2009 community and advanced
   through a short three-year baseline spin-up, rather than from a long model-generated burn-in.
   Every scenario is reported as a scenario-minus-control difference against **one shared
   no-warming control**, not a per-GCM control: the control is driven by the day-of-year
   climatology of the ensemble-median historical series, whereas each scenario uses its own
   GCM's daily series. The transient common to all runs cancels, but a GCM-specific offset
   remains. The team's Phase 2 decision (2026-11-01, decision 2) retained the shared control as an
   interim counterfactual and deferred the 11 matched per-GCM controls (P2.2); the residual
   offset is therefore stated, not removed. The intrinsic growth rates (r), carrying capacities
   (K), interaction coefficients and dispersal parameters are **not** calibrated, and the model
   does **not** reproduce the observed community to a pre-registered standard (the measured
   start-state skill is in item 13). The results are
   consequently conditional and relative, not absolute forecasts: they should be read as
   internally consistent differences between a warming pathway and the shared control, and not
   as predictions of the absolute state of the community. Full simulation-based calibration
   (approximate Bayesian computation with sequential Monte Carlo, ABC-SMC) against
   pre-registered true-skill-statistic (TSS) and occupancy criteria, with spatial block
   cross-validation and posterior predictive checks, remains future work; the acceptance
   thresholds are fixed in item 13.

3. **The air–water warming response is a same-day calibration from the national model, with a
   lagged uncertainty band (E20/E21).** The projections are driven by the **national all-Spain**
   air–water model (Stage-B fit: `ta = 0.5097`, `n = 39,944`, 980 sites), not by a
   Guadalquivir-only fit (`ta = 0.436`). The 257 coded Waterbase values are therefore in the
   calibration that reaches the GuadeX forcing. The primary relationship is calibrated on same-day
   air temperature, giving a slope of approximately 0.510 across Spain and 0.436 for the
   Guadalquivir. The team's Phase 2 decision (2026-11-01, decision 3) retained the national model
   as an explicit, documented choice and **deferred the coded-value refit and forcing
   regeneration (P2.5)**; E20 remains a stated limitation, not a resolved defect. Lagged air fits
   (3-, 7- and 30-day means) yield higher slopes, up to approximately 0.65 Spain-wide, but they
   rest on only about three in-sample lag-eligible sites (n = 104, one year of data); a lagged
   refit therefore cannot serve as a valid primary calibration. The lagged relationship is instead
   expressed as an uncertainty band on the projected warming, spanning ×0.63 to ×1.28 of the
   central value (`guadex_tw/outputs/tables/air_water_slope_uncertainty_band.csv`). Within the
   Guadalquivir subset the lagged slopes are lower than the same-day slope, so the upper end of the
   band may not apply to the basin, and the band should be treated as a conservative sensitivity
   envelope rather than as a calibrated range for this system.

4. **Barrier passability rests on a uniform structural assumption (IF).** Barrier passability is
   currently modelled using a uniform structural assumption (0.1 upstream, 0.5 downstream)
   pending full standardisation of the regional obstacle inventory's qualitative passability
   index (IF). The IF index is not decoded in the present runs; using it prematurely would
   silently mis-scale every barrier. This assumption is a deliberate, documented stand-in for
   the qualitative index and is a first-order control on connectivity that should be revisited
   once the inventory's passability semantics are standardised.

5. **The richness-loss metric is a deterministic, realised contraction, not an extinction
   probability (E8).** The richness-loss metric reports a realised, threshold-based contraction
   of the community under the deterministic trajectory; it is not a probability of extinction
   for any species. The model contains no demographic or environmental stochasticity and is
   therefore run without replicates. Consequently, the metric describes the single deterministic
   trajectory implied by the parameters and forcing, and it cannot be interpreted as, or
   converted into, a quasi-extinction or extinction probability. Adding demographic and
   environmental stochasticity with replicate simulations and species-specific quasi-extinction
   thresholds remains future work, and until it is added the reported metric should be described
   only as realised richness loss.

6. **Absence of hydrological dynamics, drought and water abstraction.** The model does not
   represent river discharge, flow intermittency, drought or water abstraction. Habitat
   suitability enters through a static environmental index and through the site temperature
   series, but there is no dynamic hydrological driver and no representation of dry reaches,
   reduced runoff or abstraction-driven habitat loss. Warming is therefore the only
   climate-linked process that varies over the projection, and the results should not be read as
   projections of the combined effect of climate and water management.

7. **No age, stage or body-size structure and no evolution or acclimation.** Individuals are
   represented by a single aggregated abundance or biomass variable per species and site; there
   is no age, stage or body-size structure and no representation of ontogenetic shifts in
   habitat or thermal tolerance. There is likewise no evolution or acclimation: thermal traits
   are fixed, so populations cannot adapt to warming, and the model cannot represent
   range shifts that depend on adaptation. Both are plausible buffering mechanisms under
   sustained warming and their omission means the projections may overstate thermal stress
   relative to a model that allowed adaptation.

8. **The river network is an on-path reconstruction, not a flow-directed channel layer.** The
   connectivity graph is reconstructed by an on-path parent rule with an explicit
   minimum-spanning-tree root join, because the regional channel layer carries no flow-direction
   attribute or elevation. It is a topology defined by on-path proximity rather than a
   flow-directed channel network derived from the hydrography. The reconstruction is acyclic,
   deterministic and validated by topology, but its edges do not encode true flow direction, and
   the root-joining step produces flow-undirected connections. The root joins are 37.7 % of the
   total tree length (5,234.5 km of 13,871.5 km; median 15.8 km, maximum 206.6 km), cover 267 of
   the 775 sites, and add 41 elevation inversions to the 4 on parent edges. Of the 1,036 matched
   obstacles, 341 (245 applied) lie on a root-join edge; a broader count of edges incident to a
   root site gives 531 (392), but that also includes parent edges entering a root. Because
   obstacle passability multiplies along a link, one **parent edge** (`1.33.2`–`1.33.58`) carries
   33 obstacles and an upstream factor of 1e-33; 28 links carry five or more obstacles and 16 fall
   below 1e-6. The Phase 2 floor (`1e-3`) and parent-only sensitivities are now reported: the
   factor ranking is identical under the primary, floor and parent-only variants (Spearman
   ρ = 1.00 for all eight metrics), the 1e-33 link is unaffected by parent-only, the floor raises
   it to 1e-3, and it contributes only 0.42–1.12 % of the dispersal rate added by the floor.
   Parent-only removes 245 applied obstacles from 111 root-join links (197 become unmatched,
   48 re-snap to parent edges). **No passability floor is applied in the primary run.** Connectivity
   results should be interpreted relative to this reconstructed topology and to the unattested
   extreme passabilities.

9. **Biomass is reported in model units without an absolute scale.** Abundance and biomass are
   expressed in model units. They are internally consistent across sites, species and scenarios
   and are suitable for relative comparison, but they carry no absolute physical scale and
   cannot be read as real standing stocks or densities without an external calibration.

10. **Short projection horizon (2026–2045).** The projections span 2026–2045, a twenty-year
    horizon. This is short relative to the generation times of the long-lived species and to the
    timescales of community reorganisation, and it truncates the analysis before the
    differentiating effects of the higher-emission pathways are fully expressed. Statements
    about scenario divergence therefore refer to differences that are emerging by 2045, not to
    long-run outcomes.

11. **Climate-model spread overlaps pathway differences.** The scenario differences are small
    relative to the spread among the individual GCMs within any single pathway, and the two
    sources of variation overlap. The central warming delivered by the lower pathways cannot be
   separated from pathway differences on the strength of a single GCM, which is why the results
   are reported as ensemble statistics and as scenario-minus-control differences against the
   shared control, and why no ordering of pathways should be inferred where the GCM spread
   overlaps it.

12. **Alien-species introductions are not represented.** The model contains no explicit
    introduction events. Which alien species to add, and where and when they arrive, remain team
    decisions. Invasive richness and biomass are therefore counted from the current species pool;
    the reported invasive metrics are conditional on that pool and do not project colonisation by
    new arrivals, so they are lower bounds on any future that includes introductions.

13. **The start-state validation is measured and marginal (4.6).** Against the post-spin-up state
    of the Phase 2 shared control, pooled sensitivity is **0.975** and pooled specificity **0.920**
    over 24 species × 775 sites; the widespread native LS has TSS **0.417** with observed occupancy
    45.3 % against modelled 77.2 % (a **+31.9-point over-prediction**). The per-species table is
    `results/climate_scenarios_phase2/validation/initial_state_validation.csv`; the values are
    identical at reported precision to the pre-fix validation. The model over-predicts occupancy,
    especially for the dominant native cyprinid. The Phase 2 acceptance criteria for the future
    calibration campaign — **per-species TSS ≥ 0.40 where ≥ 10 sites are occupied**, and
    **|modelled − observed occupancy| ≤ 0.10 for LS, SA, SP, PW and CP** — are calibration targets,
    not gates for this release. The current LS passes the TSS criterion but fails the occupancy
    criterion.

14. **The interaction-parser fix is quantified and immaterial to the reported conclusions (4.4).**
    The corrected long table changes 16 reverse-direction entries from −0.5 predation to −0.3
    competition and flags 9 mixed "competition and predation" cells (34 changed directed cells;
    227 rows, 31 ambiguous source cells / 62 rows). Re-running the three experiments against the
    fixed table moves the 2026–2045 basin metrics by at most **0.0026 native species** of basin
    mean richness (2 of 775 sites), **0.0013** of realised richness loss and **0.0011** biomass
    units on a ≈30-unit mean in the climate ensemble; at most 0.0013 species, 0.0063 biomass
    units and 0.057 ST biomass in the obstacle sweep; and at most 0.0013 species, 0.033 biomass
    units and 0.51 ST biomass in the interaction sweep. There are no aggregate sign changes and
    the factor/matrix orderings are unchanged (Spearman ρ = 1.00). The fix is required for the
    correctness and auditability of the α matrix and the implied `c = 1 − α` coefficients, but it
    does not drive the climate, obstacle or interaction conclusions. Evidence:
    `docs/GuadeX_Phase2_ParserFix_Effect.md`,
    `docs/GuadeX_Phase2_ObstacleSweep_Update.md` and
    `docs/GuadeX_Phase2_InteractionSweep_Update.md`.

15. **The competition structure materially changes the interpretation of the climate response
    (4.5).** The primary specification remains the shared-capacity `c = 1 − α`; the species-specific
    `c = −α` alternative was run over all 44 scenarios as a sensitivity. The **magnitude** of the
    weak-warming result is robust: basin native richness changes by at most 0.0129 species (≤ 10 of
    775 sites) under either structure, and the warm-endpoint richness delta is numerically
    identical (−0.006452). The **sign of the ensemble-mean richness response depends on the
    structure**: shared-K gives a mixed, slightly positive mean (+0.0020; 26 positive / 4 zero /
    14 negative), whereas the alternative is uniformly negative (44/44; mean −0.0071), and 30 of
    44 runs flip sign. The small shared-K richness-loss signal collapses to approximately zero
    under the alternative (mean −0.001551 → −0.000037; 31/44 deltas exactly zero), while the total
    biomass and trout responses amplify (~4.8× and ~7.8–8.3×). The claim "warming does not reduce
    native richness" is therefore not robust to the competition structure; the weak-warming
    headline and every primary number are unchanged. Evidence:
    `docs/GuadeX_Phase2_Competition_Sensitivity.md`.

---

## Scope of inference

Taken together, these limitations define a narrow and explicit scope. The present results are a
conditional, internally consistent, deterministic sensitivity assessment of the effect of a
mild, prescribed thermal forcing on a lightly structured and uncalibrated metacommunity model,
reported against a shared no-warming control. They are informative about the sign, relative
magnitude and spatial pattern of the modelled thermal response, and about how that response
behaves within the modelled assumptions. They are not absolute forecasts, not extinction
probabilities, and not integrated assessments of climate together with hydrology, management or
adaptation. Each of the limitations above names the specific mechanism that would have to be
added, calibrated or standardised before a wider interpretation could be supported.

---

## Source artifacts and decision record

* Thermal niches (C6), the interim projection route (C3), the air–water uncertainty band (E21),
  barrier passability (IF) and the richness-loss metric (E8) are decided and documented in
  `docs/GuadeX_Correction_Status_Sept2026.md` (§4 Team decisions, §5 Known residuals) and
  `docs/climate_scenarios.md`.
* The Phase 2 decisions (competition structure, shared control, national air–water model,
  validation acceptance thresholds) are recorded in `docs/GuadeX_Phase2_Decisions.md`; the run
  provenance is in `docs/GuadeX_Phase2_Run_Record.md`.
* The Phase 2 evidence for items 8 and 13–15 is in
  `docs/GuadeX_Phase2_ParserFix_Effect.md`,
  `docs/GuadeX_Phase2_ObstacleSweep_Update.md`,
  `docs/GuadeX_Phase2_InteractionSweep_Update.md`,
  `docs/GuadeX_Phase2_Competition_Sensitivity.md` and
  `docs/GuadeX_Phase2_PassabilityFloor_Sensitivity.md`.
* The air–water slope uncertainty band is tabulated in
  `guadex_tw/outputs/tables/air_water_slope_uncertainty_band.csv`.
* The corrected-run configuration that selects the C3 interim route and the team's biological
  options is `parameters_climate_scenarios_corrected.toml`.
