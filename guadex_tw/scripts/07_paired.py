"""Build the paired Tw~Ta calibration table (Section 7).

Joins each QC daily water-temperature observation with NASA POWER daily air
temperature at the site, lapse-corrected to the site elevation, and derives the
backward-looking lag means. Missing values in a lag window are never filled.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C

GAMMA = -0.0065  # K/m, environmental lapse rate (parameter, not truth)
LAGS = {
    "ta_mean_d1": 1, "ta_mean_d3": 3, "ta_mean_d7": 7,
    "ta_mean_d14": 14, "ta_mean_d30": 30,
    "ta_mean_d0_d2": 3, "ta_mean_d0_d6": 7, "ta_mean_d0_d13": 14, "ta_mean_d0_d29": 30,
}
POWER_CACHE = C.RAW / "nasa_power"


def _safe(s: str) -> str:
    import re
    return re.sub(r"[^0-9A-Za-z_.-]+", "_", str(s))[:120]


def load_sites() -> pd.DataFrame:
    p = C.PROCESSED / "sites.csv"
    if p.exists():
        return pd.read_csv(p, dtype={"site_id": str})
    raise FileNotFoundError("run 06_build_sites.py first")


def power_series(sid: str) -> tuple[pd.DataFrame, dict]:
    path = POWER_CACHE / f"{_safe(sid)}.json"
    if not path.exists():
        return pd.DataFrame(), {}
    j = json.loads(path.read_text(encoding="utf-8"))
    if "properties" not in j:
        return pd.DataFrame(), {}
    params = j["properties"]["parameter"]
    coords = (j.get("geometry", {}).get("coordinates") or [None, None, None])
    dates = sorted(params["T2M"].keys())
    df = pd.DataFrame({
        "date": pd.to_datetime(dates, format="%Y%m%d"),
        "ta_mean_raw": [params["T2M"].get(d, np.nan) for d in dates],
        "ta_max_raw": [params["T2M_MAX"].get(d, np.nan) for d in dates],
        "ta_min_raw": [params["T2M_MIN"].get(d, np.nan) for d in dates],
    })
    for c in ("ta_mean_raw", "ta_max_raw", "ta_min_raw"):
        df[c] = pd.to_numeric(df[c], errors="coerce").replace(-999.0, np.nan)
    meta = {"power_lon": coords[0], "power_lat": coords[1], "power_elevation_m": coords[2]}
    return df, meta


def main() -> None:
    sites = load_sites()
    obs = pd.read_csv(C.PROCESSED / "water_temp_all_qc.csv", parse_dates=["obs_date"],
                      dtype={"site_id": str})
    obs = obs.merge(sites[["site_id", "elevation_m", "site_class", "in_guadalquivir_basin",
                           "source"]].rename(columns={"source": "site_source"}),
                    on="site_id", how="left")
    print(f"{len(obs):,} QC observations across {obs['site_id'].nunique()} sites")

    paired_frames = []
    matched_rows = []
    for sid, g in obs.groupby("site_id"):
        s = sites[sites["site_id"] == sid]
        z_site = float(s["elevation_m"].iloc[0]) if len(s) and pd.notna(s["elevation_m"].iloc[0]) else np.nan
        ser, meta = power_series(sid)
        if ser.empty:
            continue
        z_cell = meta.get("power_elevation_m")
        if z_site is not np.nan and z_cell is not None and np.isfinite(z_cell):
            lapse = GAMMA * (z_site - z_cell)
        else:
            lapse = 0.0
        ser = ser.sort_values("date").set_index("date")
        ser["ta_mean"] = ser["ta_mean_raw"] + lapse
        ser["ta_max"] = ser["ta_max_raw"] + lapse
        ser["ta_min"] = ser["ta_min_raw"] + lapse
        ser["ta_mean_corr"] = ser["ta_mean"]
        ser["ta_mean_same_day"] = ser["ta_mean"]
        for col, w in LAGS.items():
            ser[col] = ser["ta_mean"].rolling(w, min_periods=w).mean()
        g = g.copy()
        g = g.merge(ser.reset_index(), left_on="obs_date", right_on="date", how="left")
        g["lapse_correction_c"] = lapse
        paired_frames.append(g)
        matched_rows.append({
            "site_id": sid, "site_name": s["site_name"].iloc[0] if len(s) else "",
            "lat": s["lat"].iloc[0] if len(s) else np.nan,
            "lon": s["lon"].iloc[0] if len(s) else np.nan,
            "elevation_m": z_site, "site_class": s["site_class"].iloc[0] if len(s) else "",
            "power_cell_lat": meta.get("power_lat"), "power_cell_lon": meta.get("power_lon"),
            "power_cell_elevation_m": z_cell,
            "power_distance_km": round(math.hypot(
                (meta.get("power_lat", np.nan) - s["lat"].iloc[0]) * 111.0,
                (meta.get("power_lon", np.nan) - s["lon"].iloc[0]) * 111.0
                * math.cos(math.radians(s["lat"].iloc[0]))), 2) if len(s) else np.nan,
            "lapse_correction_c": round(float(lapse), 3),
            "nearest_rocio_cell_lat": s["nearest_rocio_cell_lat"].iloc[0] if len(s) else np.nan,
            "nearest_rocio_cell_lon": s["nearest_rocio_cell_lon"].iloc[0] if len(s) else np.nan,
            "nearest_era5_land_cell_lat": s["nearest_era5_land_cell_lat"].iloc[0] if len(s) else np.nan,
            "nearest_era5_land_cell_lon": s["nearest_era5_land_cell_lon"].iloc[0] if len(s) else np.nan,
            "nearest_aemet_idema": s["nearest_aemet_idema"].iloc[0] if len(s) else "",
            "nearest_aemet_distance_km": s["nearest_aemet_distance_km"].iloc[0] if len(s) else np.nan,
        })

    paired = pd.concat(paired_frames, ignore_index=True)
    paired["source"] = paired["site_source"]
    keep = ["site_id", "site_name", "obs_date", "tw_obs", "n_samples", "obs_time",
            "site_class", "elevation_m", "source", "lat", "lon", "in_guadalquivir_basin",
            "upstream_basin_area",
            "ta_mean", "ta_min", "ta_max", "ta_mean_corr", "ta_mean_same_day",
            "ta_mean_d1", "ta_mean_d3", "ta_mean_d7", "ta_mean_d14", "ta_mean_d30",
            "ta_mean_d0_d2", "ta_mean_d0_d6", "ta_mean_d0_d13", "ta_mean_d0_d29",
            "lapse_correction_c", "iqr_outlier_flag", "mad_outlier_flag"]
    keep = [c for c in keep if c in paired.columns]
    paired = paired[keep].sort_values(["site_id", "obs_date"])
    paired["doy"] = paired["obs_date"].dt.dayofyear
    paired["month"] = paired["obs_date"].dt.month
    paired["year"] = paired["obs_date"].dt.year
    paired.to_csv(C.PROCESSED / "paired_observations.csv", index=False)

    matched = pd.DataFrame(matched_rows).sort_values("site_id")
    matched.to_csv(C.PROCESSED / "paired_sites_matched.csv", index=False)

    n_with_ta = int(paired["ta_mean_corr"].notna().sum())
    print(f"paired rows: {len(paired):,}; with Ta: {n_with_ta:,}")
    print(f"wrote {C.PROCESSED/'paired_observations.csv'} and paired_sites_matched.csv")
    C.record("Generated paired_observations.csv / paired_sites_matched.csv", status="SUCCESS",
             path=str(C.PROCESSED / "paired_observations.csv"), checksum=False,
             notes=f"NASA POWER Ta, GAMMA={GAMMA} K/m; {n_with_ta} paired daily obs; "
                   f"lag windows require full coverage (no interpolation)")


if __name__ == "__main__":
    main()
