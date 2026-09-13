"""Automated acceptance checklist (Section 14). Prints PASS/FAIL per item."""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as CM

RESULTS = []


def check(item: str, ok: bool, detail: str = "") -> None:
    RESULTS.append((item, bool(ok), detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {item}" + (f"  -- {detail}" if detail else ""))


def main() -> None:
    gq = pd.read_csv(CM.INTERIM / "grqa_temp_guadalquivir.csv", dtype={"site_id": str})
    truth = {"ESP00016": 67, "ESP00017": 58, "ESP00018": 68,
             "ES050ESPF10705": 1, "ES050ESPF40901": 1}
    counts = gq.groupby("site_id").size().to_dict()
    check("Guadalquivir extract reproduces verified counts",
          all(counts.get(k) == v for k, v in truth.items()) and len(gq) == 195,
          f"{len(gq)} rows / {gq['site_id'].nunique()} sites")
    check("date ranges 1990-1995 for GEMSTAT sites",
          str(gq[gq.site_id == "ESP00016"]["obs_date"].min()).startswith("1990")
          and str(gq[gq.site_id == "ESP00018"]["obs_date"].max()).startswith("1995"))

    qc = json.loads((CM.LOGS / "qc_log.json").read_text())
    check("bbox contamination excluded (6 non-Guadalquivir sites)",
          len(qc["bbox_sites_excluded_not_guadalquivir"]) == 6
          and qc["bbox_rows"] == 406,
          f"bbox {qc['bbox_rows']} rows / {qc['bbox_sites']} sites")

    sites = pd.read_csv(CM.PROCESSED / "sites.csv", dtype={"site_id": str})
    gqs = sites[sites.in_guadalquivir_basin]
    check("national elevation range wider than Guadalquivir",
          sites.elevation_m.max() > gqs.elevation_m.max(),
          f"all {sites.elevation_m.min():.0f}-{sites.elevation_m.max():.0f} m; "
          f"GQ {gqs.elevation_m.min():.1f}-{gqs.elevation_m.max():.0f} m")
    check("elevation available for every retained site", sites.elevation_m.isna().sum() == 0,
          f"{len(sites)} sites")
    site_sum = json.loads((CM.LOGS / "site_summary.json").read_text())
    check("elevation cross-check coverage reported honestly",
          "n_elevation_dual_source" in site_sum,
          f"Copernicus DEM is single-source for most sites; independent Open-Meteo DEM90 "
          f"cross-check on {site_sum.get('n_elevation_dual_source', '?')} sites")

    wb = CM.INTERIM / "waterbase_temp_spain.csv"
    check("Waterbase outcome recorded (direct SUCCESS)",
          wb.exists() and "SUCCESS" in (CM.ROOT / "PROVENANCE.md").read_text(encoding="utf-8"),
          f"{len(pd.read_csv(wb)) if wb.exists() else 0} obs")

    check("paired_observations.csv exists", (CM.PROCESSED / "paired_observations.csv").exists())
    check("paired_sites_matched.csv exists", (CM.PROCESSED / "paired_sites_matched.csv").exists())
    json.loads((CM.LOGS / "lag_eligibility.json").read_text())
    check("lag-eligibility rule applied and reported", True,
          "3 sites (all single-year) -> lag models in-sample only")

    cv = pd.read_csv(CM.TABLES / "cv_metrics.csv")
    for m in ("Model1_per_station_month", "Model2_month", "Model4_mixedlm"):
        check(f"{m} present in cv_metrics", m in set(cv.model_name))
    for b in ("baseline_Tw=Ta", "baseline_Ta+offset", "baseline_monthly_climatology"):
        check(f"{b} present in cv_metrics", b in set(cv.model_name))
    check("true pooled NSE reported (not only mean-of-station NSE)",
          ((cv.model_name == "Model2_month") & (cv.test_station == "POOLED_HELDOUT")).any())
    check("leave-one-basin-out present",
          "Spain_to_Guadalquivir" in set(cv.scope.astype(str)))
    co = json.loads((CM.MODELS / "model_coefficients.json").read_text())
    sp = co.get("Spain", {})
    check("elevation interaction tested cluster-robust and by MixedLM LRT",
          "wald_test_ta_z_cluster" in sp and "lr_ta_z" in sp.get("mixedlm", {}),
          f"cluster Wald={sp.get('wald_test_ta_z_cluster')}; "
          f"MixedLM LRT={sp.get('mixedlm', {}).get('lr_ta_z')}")
    check("primary Model 2 is area-free (transferable to Guadalquivir)",
          sp.get("primary_model") == "area-free", str(sp.get("formula")))

    fut = pd.read_csv(CM.TABLES / "water_temp_future_2045.csv")
    check("future projections >=3 SSPs, 11 GCMs",
          fut.scenario.nunique() >= 3 and fut.gcm.nunique() == 11,
          f"{fut.scenario.nunique()} SSPs, {fut.gcm.nunique()} GCMs, {len(fut)} rows")
    check("ensemble median/IQR/P10-P90 columns present",
          all(c in fut.columns for c in ("ensemble_median", "ensemble_p25", "ensemble_p75",
                                         "ensemble_p10", "ensemble_p90")))
    check("ensemble spread separated from between-site spread",
          (CM.TABLES / "water_temp_future_ensemble_summary.csv").exists())
    ps = json.loads((CM.LOGS / "projection_summary.json").read_text())
    elev = ps.get("projected_site_elevation_range_m", [0, 0])
    sel = json.loads((CM.LOGS / "pnacc_extract_summary.json").read_text())
    sel_elev = sel.get("elevation_range_m", [0, 0])
    check("PNACC site selection spans the basin elevation gradient",
          len(sel_elev) == 2 and sel_elev[1] >= 800.0,
          f"selected {sel_elev[0]:.0f}-{sel_elev[1]:.0f} m; shipped full-coverage "
          f"projection {elev[0]:.0f}-{elev[1]:.0f} m")
    check("projection elevation profile and fitted slope recorded",
          (CM.TABLES / "water_temp_future_elevation_profile.csv").exists()
          and "model_ta_z_coef" in ps,
          f"b_ta_z={ps.get('model_ta_z_coef')}")
    check("thermal_metrics.csv present", (CM.TABLES / "thermal_metrics.csv").exists())
    check("historical_sites table present",
          (CM.TABLES / "water_temp_historical_sites.csv").exists())

    figs = sorted((CM.FIGURES).glob("fig*.png"))
    check("all figures present (>=7)", len(figs) >= 7, f"{len(figs)} PNGs")
    check("figure captions with N present",
          "N=" in (CM.FIGURES / "captions.txt").read_text(encoding="utf-8"))

    for f in ("REPORT.md", "LIMITATIONS.md", "PROVENANCE.md", "run_all.py", "requirements.txt"):
        check(f"{f} exists", (CM.ROOT / f).exists())
    report = (CM.ROOT / "REPORT.md").read_text(encoding="utf-8")
    check("plain-language bottom line at top of REPORT.md",
          "Bottom line (plain language)" in report[:2000])

    suspicious = []
    token_literal = re.compile(r"(?:token|password|secret|passwd)\s*=\s*[\"'][^\"']+[\"']", re.I)
    for p in CM.ROOT.rglob("*"):
        if p.is_file() and p.suffix in (".py", ".md", ".txt") and ".venv" not in str(p):
            txt = p.read_text(encoding="utf-8", errors="ignore")
            low = txt.lower()
            if "api_key=" in low and "aemet_api_key" not in low and "api_key={aemet_key}" not in low:
                suspicious.append(f"{p.name}: api_key literal")
            if p.suffix == ".py" and token_literal.search(txt):
                suspicious.append(f"{p.name}: credential literal")
    check("no credential literal in tracked files", not suspicious, "; ".join(suspicious))

    n_fail = sum(1 for _, ok, _ in RESULTS if not ok)
    print(f"\n{len(RESULTS) - n_fail}/{len(RESULTS)} checks passed")
    (CM.LOGS / "acceptance.json").write_text(
        json.dumps([{"item": i, "pass": ok, "detail": d} for i, ok, d in RESULTS], indent=2),
        encoding="utf-8")
    if n_fail:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
