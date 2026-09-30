# LEGACY — SUPERSEDED ARTIFACTS (DO NOT CITE)

Everything in this folder is **superseded by the corrected October 2026 release** and is kept
only for comparison and provenance. Do **not** treat any number, figure or conclusion here as a
current result. The current artifacts are:

- `docs/FinalReportOctober2026.tex` / `docs/FinalReportOctober2026.pdf` — current report.
- `docs/GuadeX_Response_to_Reviewers_Oct2026.md`, `docs/GuadeX_Review_Assessment_Sept2026.md`,
  `docs/GuadeX_Correction_Plan_Sept2026.md`, `docs/GuadeX_Correction_Status_Sept2026.md`,
  `docs/GuadeX_Limitations_and_Scope_Oct2026.md` — current correction record.
- `parameters_climate_scenarios_corrected.toml` — current run configuration.
- `results/climate_scenarios_corrected/`, `results/sensitivity_obstacles_corrected/`,
  `results/alt_interactions_corrected/` — current (corrected) results.

Superseded documents moved here (pre-correction, September 2026 and earlier):

| file | why legacy |
| :--- | :--- |
| `FinalReportSeptember2026.tex` (+ build artifacts) | Superseded by `FinalReportOctober2026.tex`. |
| `GuadalquivirDecisionBrief.tex` (+ artifacts) | The September brief; carries `\pendingrerun` numbers superseded by the corrected run. |
| `ReportingMay2026.tex` | Pre-correction reporting draft. |
| `integrated_climate_obstacle_report.md` (+ pdf) | Synthesises the pre-correction `k1x` climate ensemble and obstacle sweeps. |
| `publication_manuscript_climate_fragmentation.md` | September manuscript quoting pre-correction values (0.74–0.97 °C, −0.120 slope, 9–11 % trout decline, ~90 % invasive-favouring contrast). |
| `climate_scenarios_results_report.md` | Reviews the pre-plan `results/climate_scenarios` model (annual-mean forcing, K = 10×, no ST). |
| `parameters_climate_scenarios_k1x_burnin.toml` | Legacy reproduction config (pinned pre-E4 stop rule and September heat-stress slope). |

Legacy result directories are **left in place** (referenced by plotting-script defaults) and marked
in place with a `LEGACY_README.md`: `results/climate_scenarios`,
`results/climate_scenarios_k1x_burnin`, `results/climate_scenarios_improved`,
`results/sensitivity_obstacles` (including `alt_interactions/`).
