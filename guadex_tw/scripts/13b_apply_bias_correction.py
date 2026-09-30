"""E22: apply the held-out-basin monthly Tw correction to existing products.

``13_project_guadex_sites.py`` now adds the per-month leave-one-basin-out (LOBO)
residual ``obs - pred`` (see :mod:`tw_bias_correction`) when it projects.  This
standalone, idempotent post-processor applies the *same* correction to products
that were generated before that fix, without recomputing anything from raw
climate data:

* every wide daily product ``water_temp_daily_guadex_sites_wide*.csv`` (the
  ensemble-median file and the 44 per-GCM files) - all site values get
  ``resid[month]`` added;
* ``water_temp_baseline_guadex_sites.csv`` - the annual baseline mean gets the
  mean of the 12 monthly residuals (equal month weights);
* ``water_temp_future_2045.csv`` - additive level/seasonal columns are shifted,
  the ``delta_*`` anomaly columns are left untouched and ``tw_amplitude`` is
  recomputed from the shifted seasonal means;
* ``water_temp_historical_sites.csv`` - the modelled columns only.

Large files are streamed chunk-by-chunk, written to a sibling ``.tmp`` file and
atomically replaced, so memory stays bounded and a crash never truncates the
original.  Row counts, header, column order and float64 dtype are preserved.

Idempotency: each rewritten file gets a ``<file>.biascorr.json`` sidecar keyed
by the correction id.  A file whose sidecar matches is skipped (no double
correction); use ``--force`` to rewrite anyway.  ``13_project_guadex_sites.py``
writes the same sidecars, so a corrected file produced by either path is
recognised.

Usage:
    python scripts/13b_apply_bias_correction.py               # rewrite as needed
    python scripts/13b_apply_bias_correction.py --dry-run     # report only
    python scripts/13b_apply_bias_correction.py --force       # rewrite all
    python scripts/13b_apply_bias_correction.py --skip-validation
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402
import tw_bias_correction as BC  # noqa: E402

TABLES = C.TABLES
WIDE_GLOB = "water_temp_daily_guadex_sites_wide*.csv"
BASELINE_CSV = TABLES / "water_temp_baseline_guadex_sites.csv"
FUTURE_CSV = TABLES / "water_temp_future_2045.csv"
HISTORICAL_CSV = TABLES / "water_temp_historical_sites.csv"

GCMS = ["ACCESS-CM2", "CMCC-CM2-SR5", "CNRM-ESM2-1", "EC-Earth3-Veg", "IITM-ESM",
        "KACE-1-0-G", "MIROC6", "MPI-ESM1-2-HR", "MRI-ESM2-0", "NorESM2-MM",
        "UKESM1-0-LL"]
SCENARIOS = ["ssp126", "ssp245", "ssp370", "ssp585"]
BASELINE_YEARS = (1986, 2005)  # water_temp_baseline_guadex_sites.csv window

DATE_NAMES = ("date", "obs_date", "datetime", "timestamp", "time")
DATE_MARKERS = ("-", "/")
CHUNK_ROWS = 8192

# future_2045 additive level columns and the monthly residual that applies
FUTURE_ANNUAL_COLS = ["tw_mean", "ensemble_median", "ensemble_p25", "ensemble_p75",
                      "ensemble_p10", "ensemble_p90"]
FUTURE_SUMMER_COLS = ["tw_summer_mean"]
FUTURE_WINTER_COLS = ["tw_winter_mean"]
# columns that are model outputs in the historical validation table
HISTORICAL_MODELLED_COLS = ["tw_modelled", "ci90_low", "ci90_high"]


def _fmt(x) -> str:
    return "n/a" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.4f}"


def _parses_as_dates(series: pd.Series) -> bool:
    sample = series.dropna().astype(str).head(20)
    if sample.empty:
        return False
    texty = sample.str.contains("|".join(DATE_MARKERS), regex=True).mean() > 0.5
    if not texty:
        return False
    parsed = pd.to_datetime(sample, errors="coerce")
    if parsed.notna().mean() < 0.8:
        return False
    years = parsed.dt.year
    return bool(((years >= 1900) & (years <= 2200)).all())


def _pick_date_col(columns, sample: pd.DataFrame) -> str:
    lowered = {str(c).lower(): c for c in columns}
    for name in DATE_NAMES:
        if name in lowered and _parses_as_dates(sample[lowered[name]]):
            return lowered[name]
    for c in columns:
        if _parses_as_dates(sample[c]):
            return c
    raise ValueError(f"no date column found among {list(columns)}")


def _months_from(series: pd.Series) -> np.ndarray:
    d = pd.to_datetime(series, errors="coerce")
    return d.dt.month.to_numpy()


def _value_columns(columns, sample: pd.DataFrame, date_col: str) -> list[str]:
    """Numeric value columns: everything numeric except the date column."""
    cols = []
    for c in columns:
        if c == date_col or str(c).lower() == "scenario":
            continue
        if pd.api.types.is_numeric_dtype(sample[c]):
            cols.append(c)
        elif pd.to_numeric(sample[c], errors="coerce").notna().any():
            cols.append(c)
    return cols


def _atomic_replace(tmp: Path, path: Path) -> None:
    os.replace(str(tmp), str(path))


def correct_wide(path: Path, correction: dict, resid: np.ndarray) -> dict:
    """Stream one wide daily file and add ``resid[month]`` to every site value."""
    sample = pd.read_csv(path, nrows=5)
    date_col = _pick_date_col(list(sample.columns), sample)
    value_cols = _value_columns(list(sample.columns), sample, date_col)
    tmp = path.with_name(path.name + ".tmp")
    month_sum = np.zeros(13)
    month_cnt = np.zeros(13, dtype=np.int64)
    total_sum = 0.0
    total_cnt = 0
    n_rows = 0
    first = True
    try:
        with open(tmp, "w", newline="", encoding="utf-8") as out:
            for chunk in pd.read_csv(path, chunksize=CHUNK_ROWS):
                months = _months_from(chunk[date_col])
                vals = chunk[value_cols].to_numpy(dtype=np.float64)
                finite = np.isfinite(vals)
                for m in range(1, 13):
                    sel = months == m
                    if sel.any():
                        vv = vals[sel][finite[sel]]
                        month_sum[m] += float(np.nansum(vv))
                        month_cnt[m] += int(vv.size)
                total_sum += float(np.nansum(vals))
                total_cnt += int(finite.sum())
                vals = vals + resid[months - 1][:, None]
                chunk[value_cols] = vals
                chunk.to_csv(out, index=False, header=first)
                first = False
                n_rows += len(chunk)
        _atomic_replace(tmp, path)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise
    before = total_sum / total_cnt if total_cnt else float("nan")
    after = before + float(np.sum(month_cnt[1:] * resid) / total_cnt) if total_cnt else float("nan")
    feb_before = month_sum[2] / month_cnt[2] if month_cnt[2] else float("nan")
    feb_after = feb_before + resid[1]
    mar_before = month_sum[3] / month_cnt[3] if month_cnt[3] else float("nan")
    mar_after = mar_before + resid[2]
    return {"file": path.name, "kind": "wide_daily", "n_rows": n_rows,
            "before_mean": before, "after_mean": after,
            "feb_before": feb_before, "feb_after": feb_after,
            "mar_before": mar_before, "mar_after": mar_after}


def _read_write_csv(path: Path, transform) -> pd.DataFrame:
    """Read a small CSV, apply ``transform(df) -> df``, atomically replace."""
    df = pd.read_csv(path, dtype=str)
    raw = df.copy()
    out = transform(df)
    tmp = path.with_name(path.name + ".tmp")
    out.to_csv(tmp, index=False)
    _atomic_replace(tmp, path)
    return raw, out


def correct_baseline(path: Path, correction: dict, resid: np.ndarray | None = None) -> dict:
    # Exact day-weighted shift for the 1986-2005 window the baseline averages.
    annual = BC.period_weighted_mean(correction, *BASELINE_YEARS)

    def transform(df: pd.DataFrame) -> pd.DataFrame:
        col = "tw_baseline_mean"
        df[col] = pd.to_numeric(df[col], errors="coerce") + annual
        return df

    raw, out = _read_write_csv(path, transform)
    col = "tw_baseline_mean"
    before = pd.to_numeric(raw[col], errors="coerce").mean()
    after = pd.to_numeric(out[col], errors="coerce").mean()
    return {"file": path.name, "kind": "baseline_scalar", "n_rows": len(out),
            "before_mean": float(before), "after_mean": float(after),
            "shift": annual}


def correct_future(path: Path, correction: dict, resid: np.ndarray) -> dict:
    annual = BC.annual_mean(correction)
    jja = float(np.mean([resid[5], resid[6], resid[7]]))
    djf = float(np.mean([resid[11], resid[0], resid[1]]))

    def transform(df: pd.DataFrame) -> pd.DataFrame:
        # Exact day-weighted annual shift for each row's 20-year window.
        if {"period_start", "period_end"}.issubset(df.columns):
            shifts = np.array([BC.period_weighted_mean(correction, int(a), int(b))
                               for a, b in zip(pd.to_numeric(df["period_start"]),
                                               pd.to_numeric(df["period_end"]))])
        else:
            shifts = np.full(len(df), annual)
        for col in FUTURE_ANNUAL_COLS:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce") + shifts
        for col in FUTURE_SUMMER_COLS:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce") + jja
        for col in FUTURE_WINTER_COLS:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce") + djf
        if "tw_amplitude" in df.columns:
            df["tw_amplitude"] = (pd.to_numeric(df["tw_amplitude"], errors="coerce")
                                  + (jja - djf))
        return df

    raw, out = _read_write_csv(path, transform)
    before = pd.to_numeric(raw["tw_mean"], errors="coerce").mean()
    after = pd.to_numeric(out["tw_mean"], errors="coerce").mean()
    return {"file": path.name, "kind": "future_stats", "n_rows": len(out),
            "before_mean": float(before), "after_mean": float(after),
            "shift_annual": annual, "shift_jja": jja, "shift_djf": djf}


def correct_historical(path: Path, correction: dict, resid: np.ndarray) -> dict:
    def transform(df: pd.DataFrame) -> pd.DataFrame:
        months = _months_from(df["date"])
        for col in HISTORICAL_MODELLED_COLS:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce") + resid[months - 1]
        return df

    raw, out = _read_write_csv(path, transform)
    before = pd.to_numeric(raw["tw_modelled"], errors="coerce").mean()
    after = pd.to_numeric(out["tw_modelled"], errors="coerce").mean()
    return {"file": path.name, "kind": "historical_validation", "n_rows": len(out),
            "before_mean": float(before), "after_mean": float(after),
            "note": "modelled columns only; tw_obs and ta_mean untouched"}


def expected_wide_files() -> list[Path]:
    out = [TABLES / "water_temp_daily_guadex_sites_wide.csv"]
    for scen in SCENARIOS:
        for gcm in GCMS:
            out.append(TABLES / f"water_temp_daily_guadex_sites_wide_{scen}_{gcm}.csv")
    return out


def _disk_free_gb() -> float:
    try:
        return shutil.disk_usage(str(TABLES)).free / 1e9
    except OSError:
        return float("nan")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--force", action="store_true", help="rewrite even if already corrected")
    ap.add_argument("--dry-run", action="store_true", help="report without writing")
    ap.add_argument("--skip-validation", action="store_true",
                    help="do not touch water_temp_historical_sites.csv")
    ap.add_argument("--only", default="", help="comma-separated filename substrings to process")
    args = ap.parse_args()

    t0 = time.time()
    correction = BC.load_correction()
    resid = BC.correction_vector(correction)
    print(f"E22 correction id={correction['correction_id']} variant={correction['variant']}")
    print("monthly (obs-pred) C: " +
          ", ".join(f"{m}:{BC.month_value(correction, m):+.4f}" for m in BC.MONTHS))
    print(f"pooled annual +{BC.pooled_mean(correction):.4f} C | month-balanced "
          f"+{BC.annual_mean(correction):.4f} C | disk free {_disk_free_gb():.1f} GB")

    only = [s.strip() for s in args.only.split(",") if s.strip()]

    def selected(path: Path) -> bool:
        return not only or any(s in path.name for s in only)

    results = []
    errors = []
    skipped = []

    targets = expected_wide_files()
    present = {p.name for p in TABLES.glob(WIDE_GLOB)}
    missing = [p for p in targets if p.name not in present]
    for p in missing:
        print(f"WARNING: expected wide file missing: {p.name}")

    for path in targets:
        if not path.exists():
            continue
        if not selected(path):
            continue
        results.append(_process(path, correction, resid, args, correct_wide, errors, skipped))

    for path, fn in ((BASELINE_CSV, correct_baseline),
                     (FUTURE_CSV, correct_future)):
        if path.exists() and selected(path):
            results.append(_process(path, correction, resid, args, fn, errors, skipped))
        elif not path.exists():
            print(f"WARNING: expected product missing: {path.name}")

    if not args.skip_validation:
        if HISTORICAL_CSV.exists() and selected(HISTORICAL_CSV):
            print("WARNING: water_temp_historical_sites.csv is an in-sample validation table "
                  "(Spain-wide, not a fish-model input); correcting its modelled columns "
                  "only. Use --skip-validation to leave it untouched.")
            results.append(_process(HISTORICAL_CSV, correction, resid, args,
                                    correct_historical, errors, skipped))
    else:
        print("skipping water_temp_historical_sites.csv (--skip-validation)")

    elapsed = time.time() - t0
    print(f"\nprocessed {len(results)} file(s) in {elapsed:.1f}s; "
          f"skipped {len(skipped)} already-corrected; {len(errors)} error(s)")
    for r in results:
        n = r.get("n_rows")
        n_txt = f"{n:,}" if isinstance(n, int) else "n/a"
        print(f"  {r['kind']:20s} {r['file']}: n={n_txt} "
              f"mean {_fmt(r.get('before_mean'))} -> {_fmt(r.get('after_mean'))}")
    for s in skipped:
        print(f"  SKIP already corrected: {s}")

    summary = {"correction_id": correction["correction_id"],
               "variant": correction["variant"],
               "monthly_c": {str(m): BC.month_value(correction, m) for m in BC.MONTHS},
               "annual_mean_pooled_c": BC.pooled_mean(correction),
               "annual_mean_month_balanced_c": BC.annual_mean(correction),
               "elapsed_s": elapsed, "disk_free_gb": _disk_free_gb(),
               "n_processed": len(results), "n_skipped": len(skipped),
               "n_errors": len(errors), "n_missing": len(missing),
               "missing": [p.name for p in missing],
               "results": results, "errors": errors}
    C.LOGS.mkdir(parents=True, exist_ok=True)
    (C.LOGS / "bias_correction_13b_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")
    C.record("E22 applied held-out-basin bias correction to Tw products",
             status="FAILED" if errors else "SUCCESS", checksum=False,
             notes=f"{len(results)} files, {len(missing)} missing, "
                   f"{len(errors)} errors; correction {correction['correction_id']}")
    if errors:
        raise SystemExit(1)


def _process(path: Path, correction: dict, resid: np.ndarray, args, handler,
             errors: list, skipped: list) -> dict:
    if not args.force and BC.is_corrected(path, correction):
        skipped.append(path.name)
        marker = BC.read_marker(path)
        return {"file": path.name, "kind": "already_corrected",
                "n_rows": marker.get("n_rows"),
                "before_mean": marker.get("before_mean"),
                "after_mean": marker.get("after_mean")}
    if args.dry_run:
        print(f"DRY-RUN would correct {path.name}")
        return {"file": path.name, "kind": "dry_run",
                "n_rows": None, "before_mean": None, "after_mean": None}
    t = time.time()
    try:
        info = handler(path, correction, resid)
    except Exception as exc:  # noqa: BLE001
        print(f"ERROR correcting {path.name}: {type(exc).__name__}: {exc}")
        errors.append({"file": path.name, "error": f"{type(exc).__name__}: {exc}"})
        return {"file": path.name, "kind": "error", "n_rows": None,
                "before_mean": None, "after_mean": None}
    dt = time.time() - t
    info["seconds"] = round(dt, 1)
    BC.write_marker(path, correction, n_rows=info.get("n_rows"),
                    before=info.get("before_mean"), after=info.get("after_mean"))
    print(f"  corrected {path.name} ({dt:.1f}s) "
          f"mean {_fmt(info.get('before_mean'))} -> {_fmt(info.get('after_mean'))}")
    return info


if __name__ == "__main__":
    main()
