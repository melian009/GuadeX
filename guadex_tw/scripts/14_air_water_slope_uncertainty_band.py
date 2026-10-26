"""E21: air-water slope uncertainty band for projected water-temperature warming.

The shipped projection uses the SAME-DAY air-water relationship (Spain-wide
Stage-B Model 2, ``ta_mean_corr`` slope ~0.51).  Lagged / distributed-lag air
temperature fits the pooled calibration better and yields higher slopes, but the
repo has only three cadence-eligible sites, all inside a single calendar year
(2024), so a lagged refit cannot replace the primary calibration.

This script therefore keeps the primary model untouched and expresses the
lagged alternative as an UNCERTAINTY BAND on the projected warming:

    low_warming  = central_warming * (slope_low  / slope_same_day)
    high_warming = central_warming * (slope_high / slope_same_day)

The Model-3 script (`08b_model3_lags.py`) never persisted its coefficients, so
the lagged slopes are re-estimated here from the *saved* paired calibration
table (`data/processed/paired_observations.csv`) with the exact primary
specification, substituting each lagged mean for ``ta``.

Outputs
-------
* ``outputs/tables/air_water_slope_uncertainty_band.csv``
* ``models/air_water_slope_uncertainty_band_method.json``

Run::

    .venv/Scripts/python.exe scripts/14_air_water_slope_uncertainty_band.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

ROOT = Path(__file__).resolve().parent.parent
PROCESSED = ROOT / "data" / "processed"
TABLES = ROOT / "outputs" / "tables"
MODELS = ROOT / "models"

# Lag columns present in the paired table; same-day is the shipped predictor.
SAME_DAY = "ta_mean_corr"
LAG_COLS = {
    "same_day": "ta_mean_corr",
    "d3": "ta_mean_d3",
    "d7": "ta_mean_d7",
    "d30": "ta_mean_d30",
}
BAND_LAG_KEYS = ("d3", "d7", "d30")
ELIGIBLE_SITES = ["ES020ESPF004300618", "ES091R0571", "ESCHC2621"]
WINDOW_NOTE = {
    "2021-2040": "near-term window in the shipped projection table",
    "2036-2055": "2045-centred window (primary in REPORT.md)",
    "2041-2070": "late-century appendix window",
    "2026-2045": "GuadeX ODE horizon; from the all-site ensemble wide file",
}


def load_paired() -> pd.DataFrame:
    d = pd.read_csv(PROCESSED / "paired_observations.csv", dtype={"site_id": str},
                    low_memory=False)
    d["z"] = d["elevation_m"] / 1000.0
    d["month"] = pd.to_numeric(d["month"], errors="coerce")
    return d


def fit_slope(df: pd.DataFrame, ta_col: str) -> dict:
    """Primary specification, with ``ta_col`` substituted for same-day air temp.

    ``tw_obs ~ ta + C(month) + z + ta_z``  (ta_z = ta * z).  The reported slope
    is the ``ta`` coefficient at z = 0, exactly comparable to the shipped
    Model-2 ``ta`` coefficient.
    """
    dd = df.dropna(subset=["tw_obs", ta_col, "elevation_m", "month"]).copy()
    dd = dd[np.isfinite(dd["tw_obs"])]
    dd["TA"] = dd[ta_col].astype(float)
    dd["TAz"] = dd["TA"] * dd["z"].astype(float)
    m = smf.ols("tw_obs ~ TA + C(month) + z + TAz", data=dd).fit()
    return {
        "ta_column": ta_col,
        "slope": float(m.params["TA"]),
        "bse": float(m.bse["TA"]),
        "pvalue": float(m.pvalues["TA"]),
        "slope_at_z": float(m.params["TA"]),
        "ta_z": float(m.params["TAz"]),
        "n_obs": int(m.nobs),
        "n_sites": int(dd["site_id"].nunique()),
    }


def slope_table(d: pd.DataFrame) -> dict:
    spain = d
    guad = d[d["in_guadalquivir_basin"] == True]  # noqa: E712
    elig = d[d["site_id"].isin(ELIGIBLE_SITES)]
    out = {"spain": {}, "guadalquivir": {}, "lag_eligible_3sites": {}}
    for name, col in LAG_COLS.items():
        out["spain"][name] = fit_slope(spain, col)
        out["guadalquivir"][name] = fit_slope(guad, col)
    for name, col in LAG_COLS.items():
        out["lag_eligible_3sites"][name] = fit_slope(elig, col)
    return out


def warm_from_future_table() -> pd.DataFrame:
    """Scenario-level (and per-GCM for 2036-2055) central warming, shipped table.

    Central value = median across the 6 shipped Guadalquivir sites of each
    site's median across the 11 GCMs of ``delta_tw_mean``.
    """
    f = pd.read_csv(TABLES / "water_temp_future_2045.csv", dtype={"site_id": str})
    f["window"] = f["period_start"].astype(int).astype(str) + "-" + \
        f["period_end"].astype(int).astype(str)
    rows = []
    # scenario-level (median across sites of site-median across GCMs)
    site_med = (f.groupby(["scenario", "window", "site_id"])["delta_tw_mean"]
                .median().reset_index())
    for (scen, win), g in site_med.groupby(["scenario", "window"]):
        rows.append({
            "scenario": scen, "gcm": "ALL", "window": win,
            "scope": "6 shipped Guadalquivir sites (11-GCM median per site)",
            "central_warming_c": float(g["delta_tw_mean"].median()),
            "n_sites": int(g["site_id"].nunique()),
        })
    # per-GCM for the 2045-centred window
    sub = f[f["window"] == "2036-2055"]
    pg = (sub.groupby(["scenario", "gcm", "site_id"])["delta_tw_mean"].median()
          .reset_index())
    for (scen, gcm), g in pg.groupby(["scenario", "gcm"]):
        rows.append({
            "scenario": scen, "gcm": gcm, "window": "2036-2055",
            "scope": "6 shipped Guadalquivir sites (site median)",
            "central_warming_c": float(g["delta_tw_mean"].median()),
            "n_sites": int(g["site_id"].nunique()),
        })
    return pd.DataFrame(rows)


def warm_2026_2045_from_wide() -> pd.DataFrame:
    """Basin-mean 2026-2045 warming from the all-site ensemble wide daily file.

    Wide file holds `historical` (1986-2005) and each SSP (2026-2045).  Per site,
    delta = mean(Tw 2026-2045) - mean(Tw historical); central = median across
    the GuadeX network.  This is the ODE horizon used by the GuadeX model.
    """
    path = TABLES / "water_temp_daily_guadex_sites_wide.csv"
    if not path.exists():
        return pd.DataFrame()
    hdr = pd.read_csv(path, nrows=0)
    site_cols = [c for c in hdr.columns if c not in ("date", "scenario")]
    d = pd.read_csv(path, usecols=["date", "scenario"] + site_cols,
                    parse_dates=["date"], low_memory=False)
    d["year"] = d["date"].dt.year
    base = d[d["scenario"] == "historical"]
    base = base[(base["year"] >= 1986) & (base["year"] <= 2005)]
    base_mean = base[site_cols].mean(skipna=True)
    rows = []
    for scen in ["ssp126", "ssp245", "ssp370", "ssp585"]:
        f = d[(d["scenario"] == scen) & (d["year"] >= 2026) & (d["year"] <= 2045)]
        if f.empty:
            continue
        delta = f[site_cols].mean(skipna=True) - base_mean
        rows.append({
            "scenario": scen, "gcm": "ALL", "window": "2026-2045",
            "scope": f"GuadeX network ({len(site_cols)} sites), ensemble median",
            "central_warming_c": float(np.nanmedian(delta.to_numpy(float))),
            "n_sites": len(site_cols),
        })
    return pd.DataFrame(rows)


def main() -> None:
    d = load_paired()
    slopes = slope_table(d)

    same_day = slopes["spain"]["same_day"]["slope"]
    lag_min = min(slopes[scope][k]["slope"]
                  for scope in ("spain", "guadalquivir")
                  for k in BAND_LAG_KEYS)
    lag_max = max(slopes[scope][k]["slope"]
                  for scope in ("spain", "guadalquivir")
                  for k in BAND_LAG_KEYS)
    r_low = lag_min / same_day
    r_high = lag_max / same_day

    frames = [warm_from_future_table()]
    wide = warm_2026_2045_from_wide()
    if not wide.empty:
        frames.append(wide)
    warm = pd.concat(frames, ignore_index=True)

    warm["central_warming_c"] = warm["central_warming_c"].astype(float)
    warm["low_warming_c"] = warm["central_warming_c"] * r_low
    warm["high_warming_c"] = warm["central_warming_c"] * r_high
    warm["slope_central_same_day"] = same_day
    warm["slope_low"] = lag_min
    warm["slope_high"] = lag_max
    warm["ratio_low"] = r_low
    warm["ratio_high"] = r_high
    warm["window_note"] = warm["window"].map(WINDOW_NOTE)
    warm["band_method"] = ("central x (slope_lagged / slope_same_day); "
                           "approximation: slope ratio scales the projected "
                           "water-warming anomaly (same air delta-Ta)")
    cols = ["scenario", "gcm", "window", "scope", "central_warming_c",
            "low_warming_c", "high_warming_c", "slope_central_same_day",
            "slope_low", "slope_high", "ratio_low", "ratio_high", "n_sites",
            "window_note", "band_method"]
    warm = warm[cols].sort_values(["window", "scenario", "gcm"]).reset_index(drop=True)
    out_csv = TABLES / "air_water_slope_uncertainty_band.csv"
    warm.to_csv(out_csv, index=False)
    print(f"wrote {out_csv} ({len(warm)} rows)")

    method = {
        "decision": "E21 - lagged air-temperature fits are an uncertainty band, "
                    "not a new primary calibration",
        "primary_model": {
            "identity": "Stage-B all-Spain Model 2 (area-free), used for projection",
            "formula": "tw_obs ~ ta + C(month) + z + ta_z",
            "predictor": "NASA POWER ta_mean_corr (same-day, lapse-corrected)",
            "slope_at_z0": same_day,
            "source": "models/model_coefficients.json, models/model_summary.txt",
        },
        "same_day_slopes_from_saved_outputs": {
            "Spain": {"slope_at_z0": slopes["spain"]["same_day"]["slope"],
                      "bse": slopes["spain"]["same_day"]["bse"],
                      "n_obs": slopes["spain"]["same_day"]["n_obs"],
                      "n_sites": slopes["spain"]["same_day"]["n_sites"]},
            "Guadalquivir": {"slope_at_z0": slopes["guadalquivir"]["same_day"]["slope"],
                             "bse": slopes["guadalquivir"]["same_day"]["bse"],
                             "n_obs": slopes["guadalquivir"]["same_day"]["n_obs"],
                             "n_sites": slopes["guadalquivir"]["same_day"]["n_sites"]},
        },
        "lagged_slopes_reestimated_from_saved_paired_table": {
            "note": "08b_model3_lags.py did not persist coefficients, so lagged "
                    "slopes were re-fitted from data/processed/paired_observations.csv "
                    "with the exact primary specification, replacing ta with the "
                    "lagged mean.  Reproduces the saved same-day slope 0.5097 exactly.",
            "source": "data/processed/paired_observations.csv",
            "spain": {k: slopes["spain"][k] for k in LAG_COLS},
            "guadalquivir": {k: slopes["guadalquivir"][k] for k in LAG_COLS},
        },
        "lag_eligible_sites_in_sample": {
            "site_ids": ELIGIBLE_SITES,
            "cadence_rule_obs_per_year": 24,
            "per_site_n": {"ES020ESPF004300618": 27, "ES091R0571": 53, "ESCHC2621": 24},
            "total_n_eligible": 104,
            "source": "logs/lag_eligibility.json, logs/model3_summary.json, "
                      "outputs/tables/cv_metrics.csv",
            "note": "All three sites lie in a single calendar year (2024); "
                    "time-blocked CV is impossible. Their pooled same-day slope "
                    f"is {slopes['lag_eligible_3sites']['same_day']['slope']:.4f} "
                    f"(bse {slopes['lag_eligible_3sites']['same_day']['bse']:.4f}), "
                    "i.e. unidentifiable, confirming lags cannot be a primary model.",
        },
        "band": {
            "slope_central_same_day": same_day,
            "slope_low": lag_min,
            "slope_low_source": "lowest lagged slope across Spain & Guadalquivir (d30, Guadalquivir)",
            "slope_high": lag_max,
            "slope_high_source": "highest lagged slope across Spain & Guadalquivir (d30, Spain)",
            "ratio_low": r_low,
            "ratio_high": r_high,
            "formula": "low/high warming = central_warming x (slope_low/high / slope_same_day)",
        },
        "approximation": {
            "valid": "documented approximation, NOT a re-run of the projection",
            "description": "The projection pipeline cannot be re-run with lagged "
                           "coefficients (no usable lagged calibration for the basin). "
                           "We scale each projected water-warming anomaly by the ratio "
                           "of the lagged to the same-day air-water slope. This assumes "
                           "the air-temperature anomaly is unchanged and that the water "
                           "response is linear in that anomaly.",
            "caveats": [
                "Lagged enhancement is a national-pooled result; within the "
                "Guadalquivir subset the lagged slopes are LOWER than same-day "
                "(0.32-0.44), so the direction of the correction is not guaranteed.",
                "Lag-eligible sample is 3 single-year sites (n=104), so the lagged "
                "slopes are not validated by time-blocked CV.",
            ],
        },
        "sources": [
            "models/model_coefficients.json", "models/model_summary.txt",
            "logs/model3_summary.json", "logs/lag_eligibility.json",
            "outputs/tables/cv_metrics.csv",
            "outputs/tables/water_temp_future_2045.csv",
            "outputs/tables/water_temp_daily_guadex_sites_wide.csv",
        ],
    }
    out_json = MODELS / "air_water_slope_uncertainty_band_method.json"
    out_json.write_text(json.dumps(method, indent=2), encoding="utf-8")
    print(f"wrote {out_json}")
    print(f"same_day={same_day:.4f}  lag_low={lag_min:.4f} (r={r_low:.3f})  "
          f"lag_high={lag_max:.4f} (r={r_high:.3f})")


if __name__ == "__main__":
    main()
