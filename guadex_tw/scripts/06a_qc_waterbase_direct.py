"""QC the direct Waterbase extract and combine it with the GRQA QC table.

Applies the same Section-5 rules to WATERBASE_DIRECT, de-duplicates against GRQA
by site+date and by rounded location+date, and writes
data/processed/water_temp_all_qc.csv (GRQA + WATERBASE_DIRECT).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C

CELSIUS = {"deg c", "°c", "c", "cel", "degc", "degree c", "celsius"}
MIN_OBS = 12


def main() -> None:
    log: dict = {}
    grqa = pd.read_csv(C.INTERIM / "water_temp_qc.csv", parse_dates=["obs_date"],
                       dtype={"site_id": str})
    grqa["source"] = grqa["source"].astype(str)
    if "iqr_outlier_flag" not in grqa:
        grqa["iqr_outlier_flag"] = False
    log["grqa_qc_rows"] = int(len(grqa))
    log["grqa_sites"] = int(grqa["site_id"].nunique())

    wb = pd.read_csv(C.INTERIM / "waterbase_temp_spain.csv", parse_dates=["obs_date"],
                     dtype={"site_id": str})
    n0 = len(wb)
    wb["unit_norm"] = wb["unit"].astype(str).str.strip().str.lower()
    wb = wb[wb["unit_norm"].isin(CELSIUS)]
    n_units = n0 - len(wb)
    wb["obs_value"] = pd.to_numeric(wb["obs_value"], errors="coerce")
    bad = ~wb["obs_value"].between(-1.0, 45.0)
    n_phys = int(bad.sum())
    wb = wb[~bad]
    dl = wb["detection_limit_flag"].notna() & wb["detection_limit_flag"].astype(str).str.lower().ne("")
    n_dl = int(dl.sum())
    wb = wb[~dl]
    wb["iqr_outlier_flag"] = False
    log["direct_units_removed"] = int(n_units)
    log["direct_plausibility_removed"] = n_phys
    log["direct_detection_limit_removed"] = n_dl

    # daily aggregation
    wb["_d"] = wb["obs_date"].dt.normalize()
    agg = wb.groupby(["site_id", "_d"], as_index=False).agg(
        tw_obs=("obs_value", "mean"), n_samples=("obs_value", "size"),
        site_name=("site_name", "first"), lat=("lat", "first"), lon=("lon", "first"),
        source=("source", "first"), drainage_region_name=("drainage_region_name", "first"),
        site_class=("site_class", "first"),
    ).rename(columns={"_d": "obs_date"})
    log["direct_daily_rows"] = int(len(agg))

    # de-dup against GRQA by site+date and by rounded location+date
    key_site = set(zip(grqa["site_id"], grqa["obs_date"].dt.normalize()))
    key_loc = set(zip(grqa["lat"].round(2), grqa["lon"].round(2),
                      grqa["obs_date"].dt.normalize()))
    before = len(agg)
    mask = [
        not ((r.site_id, r.obs_date) in key_site
             or (round(r.lat, 2), round(r.lon, 2), r.obs_date) in key_loc)
        for r in agg.itertuples()
    ]
    agg = agg[mask]
    log["direct_rows_removed_as_grqa_duplicates"] = before - len(agg)

    # >=12 obs per site
    counts = agg.groupby("site_id").size()
    small = set(counts[counts < MIN_OBS].index)
    n_small = int(agg["site_id"].isin(small).sum())
    agg = agg[~agg["site_id"].isin(small)]
    log["direct_sites_dropped_lt12"] = int(len(small))
    log["direct_obs_removed_lt12"] = n_small

    # MAD outlier flag per site-month
    med = agg.groupby(["site_id", agg["obs_date"].dt.month])["tw_obs"].transform("median")
    resid = (agg["tw_obs"] - med).abs()
    mad = resid.groupby([agg["site_id"], agg["obs_date"].dt.month]).transform("median").replace(0, np.nan)
    agg["mad_outlier_flag"] = (resid > 4 * mad).fillna(False)
    agg["site_country"] = "Spain"
    agg["obs_time"] = np.nan
    agg["upstream_basin_area"] = np.nan
    agg["site_ts_availability"] = np.nan
    agg["site_ts_continuity"] = np.nan
    log["direct_sites_retained"] = int(agg["site_id"].nunique())
    log["direct_obs_retained"] = int(len(agg))

    common_cols = ["site_id", "site_name", "site_country", "lat", "lon", "obs_date",
                   "tw_obs", "n_samples", "obs_time", "source", "drainage_region_name",
                   "upstream_basin_area", "iqr_outlier_flag", "mad_outlier_flag",
                   "site_ts_availability", "site_ts_continuity"]
    for c in common_cols:
        if c not in grqa:
            grqa[c] = np.nan
        if c not in agg:
            agg[c] = np.nan
    # site_class only from direct
    grqa["site_class_ext"] = np.nan
    agg["site_class_ext"] = agg["site_class"]
    combined = pd.concat([grqa[common_cols + ["site_class_ext"]],
                          agg[common_cols + ["site_class_ext"]]], ignore_index=True)
    combined["source"] = combined["source"].astype(str)
    combined.to_csv(C.PROCESSED / "water_temp_all_qc.csv", index=False)

    log["combined_rows"] = int(len(combined))
    log["combined_sites"] = int(combined["site_id"].nunique())
    log["combined_by_source"] = combined["source"].value_counts().to_dict()
    log["combined_span"] = [str(combined["obs_date"].min().date()),
                            str(combined["obs_date"].max().date())]
    (C.LOGS / "waterbase_direct_qc.json").write_text(json.dumps(log, indent=2, default=str),
                                                     encoding="utf-8")
    C.record("Generated water_temp_all_qc.csv (GRQA + Waterbase direct)", status="SUCCESS",
             path=str(C.PROCESSED / "water_temp_all_qc.csv"), checksum=False,
             notes=json.dumps(log, default=str))
    print(json.dumps(log, indent=2, default=str))


if __name__ == "__main__":
    main()
