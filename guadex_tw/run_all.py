"""Reproduce the whole GuadeX Tw pipeline from raw download to final tables.

Usage:
    python run_all.py                 # full pipeline
    python run_all.py --from 07       # resume from a step prefix
    python run_all.py --list          # show steps

Open-Meteo (step 05) is optional because its weighted rate limit blocks bulk
retrieval; NASA POWER (05b) is the consistent calibration predictor. Set
GUADEX_TW_OPENMETEO=1 to attempt it as well.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
PY = sys.executable

STEPS = [
    ("01", "scripts/01_download_grqa.py", "Download GRQA v1.3 TEMP table + metadata (HTTP range)"),
    ("02", "scripts/02_process_grqa.py", "QC + Guadalquivir/Spain extracts; verify ground truth"),
    ("03", "scripts/03_waterbase_direct.py", "Waterbase spatial metadata via EEA discodata"),
    ("04", "scripts/04_waterbase_sdi.py", "Waterbase WISE6 disaggregated CSV (Spanish river Tw)"),
    ("05", "scripts/05_fetch_openmeteo.py", "Open-Meteo ERA5 Ta + DEM (optional; rate-limited)"),
    ("06a", "scripts/06a_qc_waterbase_direct.py", "QC Waterbase direct + combine with GRQA"),
    ("05b", "scripts/05b_fetch_power.py", "NASA POWER daily Ta for all sites"),
    ("06", "scripts/06_build_sites.py", "Site table, dual elevation, classification"),
    ("07", "scripts/07_paired.py", "Paired Tw~Ta table with lapse correction + lags"),
    ("09", "scripts/09_extract_pnacc.py", "PNACC CMIP6 tmean at Guadalquivir sites"),
    ("08", "scripts/08_models.py", "Fit Models 1-4 + baselines; LOSO/time-split validation"),
    ("08b", "scripts/08b_model3_lags.py", "Model 3 lagged fits for cadence-eligible sites"),
    ("08c", "scripts/08c_leave_one_basin_out.py", "Leave-one-basin-out transfer test"),
    ("10", "scripts/10_project.py", "Project Tw to 2045 for 11 GCMs x 4 SSPs"),
    ("11", "scripts/11_figures.py", "GuadeX tables, thermal metrics and figures"),
    ("12", "scripts/12_verify.py", "Acceptance checklist (PASS/FAIL)"),
]


def run(step: str, script: str, desc: str) -> None:
    print(f"\n{'='*70}\n[{step}] {desc}\n    {script}\n{'='*70}", flush=True)
    t0 = time.time()
    r = subprocess.run([PY, str(HERE / script)], cwd=str(HERE.parent))
    dt = time.time() - t0
    if r.returncode != 0:
        print(f"[{step}] FAILED (exit {r.returncode}) after {dt:.0f}s", flush=True)
        raise SystemExit(r.returncode)
    print(f"[{step}] OK ({dt:.0f}s)", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--from", dest="start", default=None, help="resume from this step prefix")
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()
    if args.list:
        for s, sc, d in STEPS:
            print(f"{s:>4}  {sc:<40} {d}")
        return
    started = args.start is None
    for step, script, desc in STEPS:
        if step == "05" and os.environ.get("GUADEX_TW_OPENMETEO") != "1":
            print(f"[05] skipped (set GUADEX_TW_OPENMETEO=1 to run Open-Meteo)")
            continue
        if not started:
            if step == args.start:
                started = True
            else:
                continue
        run(step, script, desc)
    print("\nPipeline complete.")


if __name__ == "__main__":
    main()
