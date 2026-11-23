# GuadeX release checklist (manuscript gate)

Run this gate before publishing or sharing a report release. It is deliberately short and
mechanical: every item either passes or names the artefact that must be regenerated.

## 1. Provenance

- [ ] Record the commit (`git rev-parse HEAD`) and the primary output roots. The corrected
      October 2026 release uses `results/climate_scenarios_phase2` (45 runs),
      `results/sensitivity_obstacles_phase2` (32), and `results/alt_interactions_phase2` (48);
      `results/*_corrected` is the pre-fix audit baseline and must not be quoted as current.
- [ ] Confirm the run record (`docs/GuadeX_Phase2_Run_Record.md`) agrees with the output
      directories and the reported run counts.
- [ ] Confirm the title block carries both "October 2026" and "Phase 2 outputs recomputed
      1 November 2026" when any number comes from those outputs.

## 2. Code/decision-record to manuscript

- [ ] Every "implemented" or "documented" claim in
      `docs/GuadeX_Response_to_Reviewers_Oct2026.md` and
      `docs/GuadeX_Correction_Status_Sept2026.md` points to a location in the current
      `docs/FinalReportOctober2026.tex` (section, equation, table or limitation), not to an older
      draft.
- [ ] Every deferred item in the response has a matching limitation or future-work sentence in
      the report.
- [ ] `docs/GuadeX_Phase2_Decisions.md` and the report agree on the competition structure, the
      shared control, the national air–water model and the validation thresholds.
- [ ] No report sentence claims a statistical test, correction or sensitivity that the code and
      saved outputs do not support.

## 3. Figures

- [ ] `docs/figure_manifest.json` lists every `\includegraphics` file, in `\begin{figure}` order,
      with its SHA-256. Check the hashes:

      ```powershell
      python -c "import json,hashlib,os; m=json.load(open('docs/figure_manifest.json')); [print(f['figure'], hashlib.sha256(open(f['path'],'rb').read()).hexdigest()==f['sha256']) for f in m['figures']]"
      ```

- [ ] The ST dose-response figure is PDF Figure 13 (file `fig12_st_climate_response.png`) and
      §4.3 cites Figure 13. The obstacle hydraulic-structure figure is PDF Figure 10.
- [ ] Every panel annotation matches its caption: Figures 1 and 2 use the ensemble-median
      end-of-2045 warming (+0.74 °C); Figures 6 and 13 use the 2026–2045 mean realised warming
      (0.21–1.09 °C); any other window is labelled in the panel and the caption.
- [ ] Regenerate with the Phase 2 roots, not the legacy defaults:

      ```powershell
      $env:GUADEX_CLIMATE_OUTPUT_DIR="results/climate_scenarios_phase2"
      $env:GUADEX_SENSITIVITY_OUTPUT_DIR="results/sensitivity_obstacles_phase2"
      $env:GUADEX_ALT_OUTPUT_DIR="results/alt_interactions_phase2"
      julia --project=. scripts/plot_integrated_report_figures.jl
      julia --project=. scripts/plot_report_figures.jl
      ```

## 4. Numbers and rounding

- [ ] Sample size, slope and $R^2$ values are rounded to three significant figures throughout
      (for example `0.0751`, `−0.00526`, `−0.000673`, `R² = 0.143`). Tables and captions use the
      same rounding as the body.
- [ ] No number is quoted from a single run where the text claims an ensemble value (for example
      the 2026 exposure value is 37.30 d, not 42.7 d).
- [ ] Abstract, body, tables, captions, Conclusions and Limitations agree on the direction and
      significance of every aggregate slope.
- [ ] The E21 air–water values and the realised basin-mean warming are labelled as different
      statistics wherever both appear.
- [ ] The 1e-33 obstacle-only factor and the 1e-34 combined legacy-dam factor are both stated in
      the same clause when the extreme link is discussed.

## 5. References and terminology

- [ ] Every DOI resolves (spot-check the entries changed since the previous release). The
      Herrera entry must not return.
- [ ] *Chelon ramada*, not *Liza ramada*, in text and docs; code `LR` unchanged.
- [ ] AB (*Aphanius baeticus*) and CP (*Cobitis paludica*) are not called cyprinids.
- [ ] No stale strings: `-1.0 no coexistence`, `all results ... applied`, `dominant cyprinids`,
      `Liza ramada`, `blocked (0.1)` without the gloss, `Post-spin-up validation`, `18.37 °C`.

## 6. Build and tests

- [ ] `julia --project=. test/runtests.jl` passes; record the pass count and testsets. The
      clean-checkout testset may skip when the wide forcing file is absent.
- [ ] `latexmk -pdf -interaction=nonstopmode -halt-on-error docs/FinalReportOctober2026.tex`
      succeeds; record the page count in `docs/GuadeX_Correction_Status_Sept2026.md`.
- [ ] The PDF opens without unresolved references (`??`) or missing figures.
