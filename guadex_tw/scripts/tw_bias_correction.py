"""E22: held-out-basin (LOBO) monthly cold-bias correction for GuadeX Tw.

The Guadalquivir water-temperature series produced by
``13_project_guadex_sites.py`` inherited the cold bias of the calibrated
``Model2`` (``10_project.py``): the elevation-aware leave-one-basin-out (LOBO)
transfer test (``08c_leave_one_basin_out.py``) reports a mean held-out residual
``obs - pred`` of about +1.1 degC, rising to +2.2 degC in February/March.

This module is the single source of truth for that correction:

* :func:`load_correction` returns the per-month ``obs - pred`` means.  It reads
  the canonical JSON snapshot ``models/lobo_monthly_bias_correction.json`` and
  recomputes (and rewrites) it from ``models/lobo_predictions.csv`` when the
  snapshot is missing or stale.
* :func:`correction_vector` returns a length-12 vector indexed by ``month - 1``.
* :func:`write_marker` / :func:`is_corrected` implement the idempotency sidecar
  used by both the projection script and the batch post-processor, so a file is
  never corrected twice.

Sign convention: the value is the amount to **add** to the model prediction
(``obs - pred``); it is positive because the model is too cold.
"""
from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
LOBO_CSV = ROOT / "models" / "lobo_predictions.csv"
CORRECTION_JSON = ROOT / "models" / "lobo_monthly_bias_correction.json"

#: The right LOBO variant: the national model includes elevation (z, ta_z), the
#: same specification used by ``10_project.py``/``13_project_guadex_sites.py``.
LOBO_VARIANT = "LOBO_Model2_elevation"
MONTHS = tuple(range(1, 13))
JJA = (6, 7, 8)
DJF = (12, 1, 2)
MARKER_SUFFIX = ".biascorr.json"


def compute_monthly_residuals(variant: str = LOBO_VARIANT) -> tuple[dict[int, float], dict[int, int], float]:
    """Mean held-out residual ``tw_obs - pred`` per calendar month.

    Reads ``models/lobo_predictions.csv`` (771 rows, ~50 kB) and returns
    ``(residual_by_month, n_by_month, pooled_mean)`` where ``pooled_mean`` is the
    mean over all observations (the ~+1.106 degC headline).  Only finite pairs
    are used.
    """
    if not LOBO_CSV.exists():
        raise FileNotFoundError(f"LOBO predictions not found: {LOBO_CSV}")
    sums = {m: 0.0 for m in MONTHS}
    counts = {m: 0 for m in MONTHS}
    pooled_sum = 0.0
    pooled_n = 0
    with open(LOBO_CSV, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            if row.get("model_name") != variant:
                continue
            obs = float(row["tw_obs"])
            pred = float(row["pred"])
            if not (np.isfinite(obs) and np.isfinite(pred)):
                continue
            date = row["obs_date"]
            month = int(date[5:7]) if len(date) >= 7 else int(date.split("-")[1])
            r = obs - pred
            sums[month] += r
            counts[month] += 1
            pooled_sum += r
            pooled_n += 1
    if pooled_n == 0:
        raise ValueError(f"no rows for LOBO variant {variant!r} in {LOBO_CSV}")
    resid = {m: sums[m] / counts[m] for m in MONTHS}
    return resid, counts, pooled_sum / pooled_n


def _correction_id(variant: str, resid: dict[int, float]) -> str:
    payload = variant + "|" + "|".join(f"{m}:{resid[m]:.12f}" for m in MONTHS)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:16]


def build_correction(variant: str = LOBO_VARIANT) -> dict:
    """Compute the full correction record (and its stable ``correction_id``)."""
    resid, counts, pooled = compute_monthly_residuals(variant)
    month_balanced = float(np.mean([resid[m] for m in MONTHS]))
    return {
        "variant": variant,
        "definition": "mean(tw_obs - pred) by calendar month over leave-one-basin-out "
                      "(Spain non-Guadalquivir train -> Guadalquivir test) folds",
        "sign": "add to model prediction (positive = model too cold)",
        "source": str(LOBO_CSV.relative_to(ROOT)).replace("\\", "/"),
        "months": {str(m): round(resid[m], 10) for m in MONTHS},
        "n_by_month": {str(m): counts[m] for m in MONTHS},
        "n_obs": int(sum(counts.values())),
        # Pooled over all held-out observations: the ~+1.106 degC headline bias.
        "annual_mean_pooled": float(pooled),
        # Mean of the 12 monthly means (equal month weights); what a full-year
        # daily series averages to once every day carries its month's residual.
        "annual_mean_month_balanced": month_balanced,
        "jja_mean": float(np.mean([resid[m] for m in JJA])),
        "djf_mean": float(np.mean([resid[m] for m in DJF])),
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "correction_id": _correction_id(variant, resid),
    }


def _as_float_months(raw: dict) -> dict[int, float]:
    return {int(k): float(v) for k, v in raw.items()}


def load_correction(refresh: bool = False) -> dict:
    """Load the canonical correction, recomputing the JSON snapshot if needed."""
    if not refresh and CORRECTION_JSON.exists():
        try:
            data = json.loads(CORRECTION_JSON.read_text(encoding="utf-8"))
            if data.get("variant") == LOBO_VARIANT and set(data.get("months", {})) == {
                    str(m) for m in MONTHS}:
                data["months"] = _as_float_months(data["months"])
                if "correction_id" not in data:
                    data["correction_id"] = _correction_id(data["variant"], data["months"])
                if "annual_mean_month_balanced" not in data:
                    data["annual_mean_month_balanced"] = annual_mean(data)
                if "annual_mean_pooled" not in data:
                    data["annual_mean_pooled"] = data["annual_mean_month_balanced"]
                return data
        except (json.JSONDecodeError, KeyError, TypeError, ValueError):
            pass
    data = build_correction()
    data["months"] = _as_float_months(data["months"])
    CORRECTION_JSON.write_text(json.dumps(data, indent=2), encoding="utf-8")
    return data


def month_value(correction: dict, month: int) -> float:
    """Correction for one month, accepting int or str keyed ``months``."""
    months = correction["months"]
    if month in months:
        return float(months[month])
    return float(months[str(month)])


def correction_vector(correction: dict) -> np.ndarray:
    """Length-12 vector of corrections, index ``month - 1`` (January = 0)."""
    return np.array([month_value(correction, m) for m in MONTHS], dtype=np.float64)


def annual_mean(correction: dict) -> float:
    """Equal-month mean of the 12 corrections (applied to a full-year series)."""
    if "annual_mean_month_balanced" in correction:
        return float(correction["annual_mean_month_balanced"])
    return float(np.mean(correction_vector(correction)))


def pooled_mean(correction: dict) -> float:
    """Observation-pooled mean held-out residual (the ~+1.106 degC headline)."""
    if "annual_mean_pooled" in correction:
        return float(correction["annual_mean_pooled"])
    return annual_mean(correction)


def period_weighted_mean(correction: dict, y0: int, y1: int) -> float:
    """Day-weighted mean correction over calendar years ``y0..y1`` inclusive.

    This is the exact mean shift a full daily series over ``y0..y1`` receives,
    accounting for the actual number of days per month (leap years included).
    """
    resid = correction_vector(correction)
    import pandas as pd
    dates = pd.date_range(f"{y0}-01-01", f"{y1}-12-31", freq="D")
    months = dates.month.to_numpy()
    return float(resid[months - 1].mean())


def marker_path(path: Path) -> Path:
    return path.with_name(path.name + MARKER_SUFFIX)


def write_marker(path: Path, correction: dict, *, n_rows: int | None = None,
                 before: float | None = None, after: float | None = None,
                 note: str = "") -> None:
    """Write the idempotency sidecar for one corrected file."""
    payload = {
        "file": path.name,
        "correction_id": correction["correction_id"],
        "variant": correction["variant"],
        "applied_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    }
    if n_rows is not None:
        payload["n_rows"] = int(n_rows)
    if before is not None:
        payload["before_mean"] = float(before)
    if after is not None:
        payload["after_mean"] = float(after)
    if note:
        payload["note"] = note
    marker_path(path).write_text(json.dumps(payload, indent=2), encoding="utf-8")


def read_marker(path: Path) -> dict:
    """Return the sidecar contents for ``path`` (empty dict when absent)."""
    marker = marker_path(path)
    if not marker.exists():
        return {}
    try:
        return json.loads(marker.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return {}


def is_corrected(path: Path, correction: dict) -> bool:
    """True when ``path`` already carries a marker for this exact correction."""
    return read_marker(path).get("correction_id") == correction["correction_id"]
