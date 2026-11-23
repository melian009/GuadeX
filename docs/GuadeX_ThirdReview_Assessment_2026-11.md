# Third review assessment — 2026-11-23

**Reviewer:** Carlos F. (third review, `docs/comentarios_technical_report_CarlosF.md`).
**Reviewed snapshot:** `b8ee15a`. **Current tree:** `5a8077c` plus this triage. The only
commit between the reviewed snapshot and the current tree is `5a8077c`, which rewrote the
manuscript's terminology from "biomass" to "density" (and the same wording in scripts and
docs). The reviewer's quoted "biomass" sentences therefore read "density" in the current tree;
the substantive reviewer points are unchanged.

This note records the disposition of the eleven major points (R1–R11) and twelve suggested items
(S1–S12) and the evidence for each. It is the companion to the edited manuscript
(`docs/FinalReportOctober2026.tex`, recompiled to a 42-page PDF) and to the updated
`docs/GuadeX_Limitations_and_Scope_Oct2026.md`.

---

## 1. Major points (R1–R11)

| Item | Reviewer's core claim | Disposition | Action / evidence |
|---|---|---|---|
| **R1** | References are inaccurate or do not exist (Altermatt, Perkin & Gido, Tonkin, Herrera). | **Accepted** | Corrected the three metadata errors in `docs/final_report_references.bib`. `Herrera2022` (DOI `10.3390/w14152345`) resolves to an unrelated canyoning paper and was removed; the dam/potamodromous sentence is carried by the already-cited `RamosMerchante2018`. The four corrected DOIs resolve: `10.1007/s10452-013-9450-3`, `10.1890/12-0318.1`, `10.1038/d41586-019-01877-1`, `10.1111/j.1095-8649.2010.02762.x`. Added `Elliott2010` for the trout sublethal-threshold statement. |
| **R2** | Eq. (1) does not describe the implemented model (extra `+Σc·N/K` term; `11.5` vs `1.5` units). | **Accepted for Eq. (1)** | The `11.5`/`1.5` reading was a formatting artefact; no such text remains. Eq. (1) was the only incorrect expression: it carried `+Σ c_sj N/K` alongside `−ΣN/K`, which cancels under `c = 1 − α`. Eq. (1) now carries `+Σ α_sj N/K`, and the following sentence states the identity `1 − ΣN/K + ΣαN/K ≡ 1 − ΣcN/K` (same as Eq. 7 and `src/ode_model.jl:198-228`). A substitution check is written into the text. `docs/GuadeX_Phase2_Competition_Sensitivity.md:104-122` agrees. |
| **R3** | "No coexistence" is impossible; the `[−1,2]` clamp is undocumented; `blocked = 0.1` is unclear. | **Accepted** | (a) The parser maps "no coexistence" to allopatry (`α = 0`), not to `−1.0`; the mapping is now listed in §3.5 (`:814,837,874,1078,1104-1109` in `src/data_preparation.jl`). (b) The clamp is documented in §3.2 immediately after Eq. (1). (c) "Blocked (0.1)" now reads "blocked (0.1; strong attenuation, not complete closure)". |
| **R4** | Abstract claims "rose from 23 to 45 days", but 42.7 d in 2026 already exceeds 45. | **Accepted** | The 42.7 d is the single `ssp126/ACCESS-CM2` 2026 value. The abstract and §4.3 now give the ensemble mean: 37.30 d in 2026 (23.28–51.72) and 45.35 d in 2045 (30.19–65.72). Of the 21.9-d rise from the 23.43-d 1986–2005 baseline, 13.87 d are a pre-2026 product offset and 8.05 d arise within 2026–2045. Evidence: `results/climate_scenarios_phase2/*/*/export/levels/exposure_summary.csv`. |
| **R5** | Section 5.3's caveat is too vague; the mechanisms of the "blocking helps natives" result are not stated; the trout "same direction" sentence is imprecise. | **Accepted** | §5.3 now names the three mechanisms: destination-independent dispersal, near-zero-capacity `DIFICIL` sites, and the absence of local extinction/recolonisation. §5.6 now says warming and hydraulic structure act through different mechanisms and magnitudes: warming lowers trout density, reduced passability raises it, and only higher upstream cost acts weakly alongside warming. |
| **R6** | "Climate endpoint is a minor factor" contradicts the 0.142 climate SS for interaction-sweep native density. | **Accepted** | The claim is now scoped to the obstacle sweep (`SS ≤ 0.015`) and the interaction-sweep value (0.142) is stated in the same sentence, the abstract, and the results. Table 4 (`tab:ss`) already carries it. |
| **R7** | Figure 13 mixes 2036–2045 and 2026–2045 windows; panel (c) shows an all-year maximum; §4.3 cites the wrong figure; Figure 10's trout caption overstates passability. | **Accepted** | `scripts/plot_integrated_report_figures.jl` now passes `window=:early` to `realised_forcing_axis` and filters exposure to 2045. `scripts/plot_report_figures.jl` does the same and relabels the fitted line "Pooled OLS slope" (it was labelled "GCM-fixed" while plotting pooled fields). §4.3 now cites Figure 13; the caption gives the 2026–2045 regressor and the 2045 exposure range. The two rows of Figure 10 differ in visual weight because they show different quantities (final density and relative change), and the caption states this and gives the effect size in units. |
| **R8** | Trout text uses 20 °C but the reviewer said 25 °C; E21 and realised warming are inconsistent; Figure 1's +0.62 does not match; Table 3 and Figure 1 disagree. | **Accepted, with correction** | The 20 °C limit is a sublethal feeding/growth threshold and is now stated as such, with the ~25 °C incipient lethal temperature cited (`Elliott2010`). E21 and the realised pathway means are now explicitly distinguished as different statistics (§4.1). Figure 1's reference line is regenerated to the same ensemble-median end-of-2045 statistic as Figure 2 (+0.74 °C → 16.6 °C; the 2026–2045 mean is +0.58 °C). Table 3/Figure 1 agreement is retained. One warming statistic is used per figure and each is named. |
| **R9** | Dose–response slopes have no GCM-aware inference; 44 runs are not 44 independent observations. | **Accepted** | `scripts/gcm_aware_inference.jl` now reports the pooled OLS slope, the GCM fixed-effect slope, GCM-clustered SE and two-sided $p$ (11 GCM clusters, 10 df), the per-GCM slope SD and the between-GCM SD (Table 5, `tab:inference`). The pooled slopes are unchanged; GCM-clustered $p$ is $0.110$, $0.323$ and $0.700$ for total density, native richness and richness loss. The fixed-effect fit gives $+0.0827$, $-0.00376$ and $-0.000825$ per °C; only total density reaches significance ($p=0.042$), on 11 between-GCM units, and the per-GCM slope SDs are 1.9×, 7.3× and 10× the fixed-effect slopes. The abstract, §4.1 and §7 report the weak/GCM-heterogeneous reading rather than a robust directional sign. Reported slope/$R^2$ values are rounded to three significant figures. |
| **R10** | The abstract's "all applied" overclaims the correction status. | **Accepted** | Reads "as implemented or deferred (see Section 6)" and retains the pointer to §4.7. |
| **R11** | Title date is stale; the 1e-33/1e-34 passability factors are inconsistent. | **Accepted (R11b), reviewer error (R11a)** | The title block now carries both "October 2026" and "Phase 2 outputs recomputed 1 November 2026". The 1e-33/1e-34 pair is correct: `0.1^33 = 10^-33` is the obstacle-only product; the legacy dam restriction multiplies the same directed entry by a further 0.1, giving the combined `10^-34` used in the primary run. The report now gives both in one clause. Changing the number would introduce an error; `docs/GuadeX_Phase2_PassabilityFloor_Sensitivity.md:150-183` documents the pair. |

---

## 2. Suggested items (S1–S12)

| Item | Disposition | Where / note |
|---|---|---|
| **S1** "Validation" overclaims a successful check | **Accepted** | Renamed the §4.2 heading and Table 2 caption to "consistency check"; the limitation and the status/limitations docs follow. The caveat that the comparison is marginal is retained. |
| **S2** "Blocked (0.1)" is opaque | **Accepted** | "blocked (0.1; strong attenuation, not complete closure)" in §3.4. |
| **S3** Trout 20 °C is sublethal | **Accepted** | Stated in §3.2 with `Elliott2010`; the ~25 °C incipient lethal temperature is named. |
| **S4** Factor ranking is range-dependent | **Accepted** | One sentence in §4.4: the ranking is conditional on the explored ranges (upstream cost spans 50×; climate has two endpoints over 20 years). |
| **S5** "Dominant cyprinids" is wrong (AB, CP are not cyprinids) | **Accepted** | Replaced with "dominant native species" in the manuscript, Figure 1 caption and `docs/climate_scenarios_biological_report.md`. |
| **S6** *Liza ramada* → *Chelon ramada* | **Accepted** | Changed in the manuscript, `docs/data_summary.md`, `docs/climate_scenarios_biological_report.md` and `docs/Fish Growth and Dispersal Data Request.md`; the code `LR` is unchanged. |
| **S7** Deliberate disconnection gives fixed basin totals | **Accepted** | Added to §2.2: zero local growth for the three migratory taxa holds their basin totals constant; only spatial redistribution is possible. |
| **S8** Species-pool provenance is unclear | **Accepted** | Added to §2.2: the pool is the 2006–2009 survey; exotics established but not recorded then (e.g. *Silurus glanis* in Iznájar Reservoir) are excluded. |
| **S9** Table 2's large control changes invite misreading | **Accepted** | §4.6 now states that those changes are the shared uncalibrated transient, not a warming signal, and that the scenario-minus-control differences are small. |
| **S10** The cancellation claim is untested | **Accepted, with new analysis** | Added `results/climate_scenarios_phase2/per_species_scenario_minus_control.csv` and Table 6 (`tab:per-species`), computed from the saved `species_timeseries.csv` exports. The conclusion now reads as opposing species responses (large gains in a few warm-adapted natives; losses in brown trout and several other natives), not as a zero net change. |
| **S11** Claimed analyses are not in the report | **Accepted** | The report now records the E7 Latin-hypercube/Saltelli–Sobol pilot (§4.7), with its size and bootstrap interval, and states production Sobol as future work (Limitations 17). The response document no longer cites §3.5/§4.4 as the evidence. |
| **S12** E23 high-elevation limitation is not stated | **Accepted** | Added as report §6, Limitations 16 and to `docs/GuadeX_Limitations_and_Scope_Oct2026.md` item 16. |

---

## 3. Reviewer claims not reproduced

- **R11a (1e-33 vs 1e-34).** Both are correct and describe different products (obstacle-only vs
  combined with the legacy dam restriction). The manuscript now distinguishes them; the number is
  unchanged.
- **R8 (42.7 d in 2026).** This is one run (`ssp126/ACCESS-CM2`), not the ensemble mean; the
  decomposition replaces the single value.
- **R2 (11.5 vs 1.5 units).** No such text exists in the current tree; the substantive
  equation error is Eq. (1) and is fixed.
- **R1 (Herrera2022).** Confirmed: the DOI resolves to an unrelated canyoning paper. The entry is
  removed rather than repaired.

---

## 4. Verification

- **PDF:** `docs/FinalReportOctober2026.pdf`, 42 pages, compiled with `latexmk` from
  `docs/FinalReportOctober2026.tex` (bibliography: `docs/final_report_references.bib`).
- **Figures:** regenerated from the Phase 2 outputs:
  `results/climate_scenarios_phase2/figures/report_plots/fig01_thermal_niches.png`,
  `results/sensitivity_obstacles_phase2/report_plots/fig_temperature_response.png`,
  `results/climate_scenarios_phase2/figures/report_plots/fig12_st_climate_response.png`, with
  Figure 2 (`thermal_niches.png`) and the other report plots checked against the same data. Hashes
  and regeneration commands are in `docs/figure_manifest.json`.
- **Figure numbering:** 15 `\begin{figure}` environments in the `.tex`. The reviewer's Figure 10 is
  the hydraulic-structure figure and its Figure 13 is the brown-trout dose–response
  (`fig12_st_climate_response.png`, PDF Figure 13). §4.3 cites Figure 13.
- **Tests:** `julia --project=. test/runtests.jl` — 3,941/3,941 tests pass (same suite as
  `docs/GuadeX_Phase2_Run_Record.md`; no failures in the run log, and the final
  `Global sensitivity toolkit (E7)` testset reports 70/70). The clean-checkout testset skips
  when the wide forcing file is absent, as expected.

---

## 5. Files touched

- `docs/FinalReportOctober2026.tex` (manuscript), `docs/FinalReportOctober2026.pdf` (recompiled),
  `docs/final_report_references.bib`.
- `docs/GuadeX_Limitations_and_Scope_Oct2026.md`, `docs/GuadeX_Response_to_Reviewers_Oct2026.md`,
  `docs/GuadeX_Correction_Status_Sept2026.md`.
- `docs/data_summary.md`, `docs/climate_scenarios_biological_report.md`,
  `docs/Fish Growth and Dispersal Data Request.md`.
- `scripts/plot_integrated_report_figures.jl`, `scripts/plot_report_figures.jl`,
  `scripts/per_species_scenario_minus_control.jl` (regenerates the Table 6 source CSV) and
  `scripts/gcm_aware_inference.jl` (regenerates the Table 5 source CSV).
- `docs/figure_manifest.json`, `docs/GuadeX_Release_Checklist.md`, this note.
- `results/climate_scenarios_phase2/per_species_scenario_minus_control.csv`, and the regenerated
  figures listed above.
