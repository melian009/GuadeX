"""Process GRQA v1.3 TEMP_GRQA.csv into raw / QC / Spain / Guadalquivir extracts.

Streams the 1.05 GB semicolon CSV in chunks (never holds it in memory), applies
the Section 5 QC rules, and clips the Guadalquivir subset with the official CHG
ES050 water-body catchment polygons.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Point
from shapely.ops import unary_union

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C

SRC = C.RAW / "TEMP_GRQA.csv"
DUP = C.RAW / "TEMP_GRQA_dup_obs.csv"
BASIN_SHP = C.ROOT.parent / "data" / "Version_02-01-2026-Masas" / "Cuencas_masas_agua_4c.shp"

USECOLS = [
    "obs_id", "lat_wgs84", "lon_wgs84", "obs_date", "obs_time", "site_id", "site_name",
    "site_country", "upstream_basin_area", "upstream_basin_area_unit",
    "drainage_region_name", "param_code", "param_name", "obs_value", "unit", "source",
    "obs_percentile", "obs_iqr_outlier", "site_ts_availability", "site_ts_continuity",
    "detection_limit_flag", "param_form",
]

IBERIA = dict(lat_min=36.0, lat_max=44.0, lon_min=-10.0, lon_max=3.0)
GQ_BBOX = dict(lat_min=36.5, lat_max=39.5, lon_min=-7.5, lon_max=-2.0)

SPAIN_COUNTRY = {"spain", "es", "esp", "españa", "espana"}

GQ_NAME_RE = (
    r"Guadalquivir|Guadalbull|Genil|Guadalimar|Jándula|Jandula|Bembézar|Bembezar|Viar|"
    r"Guadajoz|Corbones|Guadaira|Guadiana Menor|Rivera de Huelva"
)

CELSIUS = {
    "deg c", "°c", "c", "cel", "degc", "degree c", "degrees c", "degrees celsius",
    "degree celsius", "celsius", "oc", "ºc", "deg. c", "degrees centigrade",
    "centigrade", "degree centigrade",
}

QC_LOG: dict[str, object] = {}


def normalise_unit(u: object) -> str:
    if u is None or (isinstance(u, float) and np.isnan(u)):
        return ""
    return str(u).strip().lower().replace("\u00b0", "°")


def load_basin_polygon():
    g = gpd.read_file(BASIN_SHP)
    es050 = g[g["EUMASCod"].astype(str).str.startswith("ES050")]
    poly = unary_union(es050.geometry.values)
    return g, es050, poly


def main() -> None:
    gdf, es050, basin_ll = load_basin_polygon()
    # project basin polygon to WGS84 for lon/lat point tests
    basin_wgs = gpd.GeoSeries([basin_ll], crs=25830).to_crs(4326).iloc[0]

    print("Scanning TEMP_GRQA.csv ...")
    unit_counts: dict[str, int] = {}
    keep_parts = []
    bbox_parts = []
    n_total = 0
    reader = pd.read_csv(
        SRC, sep=";", usecols=USECOLS, dtype=str, chunksize=400_000,
        encoding="utf-8", encoding_errors="replace", low_memory=False,
    )
    for ci, ch in enumerate(reader):
        n_total += len(ch)
        ch = ch[ch["param_code"].astype(str).str.upper() == "TEMP"]
        if ch.empty:
            continue
        ch["_u"] = ch["unit"].map(normalise_unit)
        lat = pd.to_numeric(ch["lat_wgs84"], errors="coerce")
        lon = pd.to_numeric(ch["lon_wgs84"], errors="coerce")
        country = ch["site_country"].astype(str).str.strip().str.lower()
        in_iberia = (
            lat.between(IBERIA["lat_min"], IBERIA["lat_max"])
            & lon.between(IBERIA["lon_min"], IBERIA["lon_max"])
        )
        # keep a bbox-of-any-country snapshot for the contamination cross-check
        bb = ch[
            lat.between(GQ_BBOX["lat_min"], GQ_BBOX["lat_max"])
            & lon.between(GQ_BBOX["lon_min"], GQ_BBOX["lon_max"])
        ]
        if not bb.empty:
            bbox_parts.append(bb[["site_id", "site_name", "source", "lat_wgs84", "lon_wgs84"]])
        is_spain = country.isin(SPAIN_COUNTRY)
        is_empty_country = country.isin(["", "nan", "none", "null"])
        sel = is_spain | (is_empty_country & in_iberia)
        ch = ch[sel]
        if ch.empty:
            continue
        for k, v in ch["_u"].value_counts().items():
            unit_counts[k] = unit_counts.get(k, 0) + int(v)
        ch = ch[ch["_u"].isin(CELSIUS)]
        if ch.empty:
            continue
        keep_parts.append(ch)
        if (ci + 1) % 10 == 0:
            print(f"  chunk {ci+1}: total={n_total:,} kept so far="
                  f"{sum(len(p) for p in keep_parts):,}")

    raw = pd.concat(keep_parts, ignore_index=True)
    del keep_parts
    print(f"TEMP rows scanned: {n_total:,}; Iberian/Spain TEMP kept: {len(raw):,}")

    QC_LOG["units_all_iberia_temp"] = dict(sorted(unit_counts.items(), key=lambda x: -x[1])[:15])
    QC_LOG["n_iberia_temp_rows_rule1_2"] = int(len(raw))

    raw = raw.rename(columns={"lat_wgs84": "lat", "lon_wgs84": "lon"})
    for c in ("lat", "lon", "obs_value", "upstream_basin_area"):
        raw[c] = pd.to_numeric(raw[c], errors="coerce")
    raw["obs_date"] = pd.to_datetime(raw["obs_date"], errors="coerce")
    raw["obs_iqr_outlier"] = raw["obs_iqr_outlier"].astype("string").str.lower().str.strip()
    raw = raw.dropna(subset=["lat", "lon", "obs_date", "obs_value"])

    raw["site_id"] = raw["site_id"].astype(str).str.strip()
    raw["source"] = raw["source"].astype(str).str.strip()
    raw.to_csv(C.INTERIM / "water_temp_raw.csv", index=False)

    # ---------- QC ----------
    qc = raw.copy()
    n0 = len(qc)

    # rule 3: physical plausibility
    bad_phys = ~qc["obs_value"].between(-1.0, 45.0)
    n_bad_phys = int(bad_phys.sum())
    qc = qc[~bad_phys]
    QC_LOG["rule3_plausibility_removed"] = n_bad_phys

    # rule 4: cross-source duplicates -> drop obs_id_2 members
    dup = pd.read_csv(DUP, sep=";", usecols=["obs_id_1", "obs_id_2"], dtype=str)
    dup2 = set(dup["obs_id_2"].dropna().astype(str))
    n_dup = int(qc["obs_id"].isin(dup2).sum())
    qc = qc[~qc["obs_id"].isin(dup2)]
    QC_LOG["rule4_duplicate_pairs_in_meta"] = int(len(dup))
    QC_LOG["rule4_duplicate_rows_removed"] = n_dup

    # rule 6: detection limits
    dl_raw = qc["detection_limit_flag"]
    dl = dl_raw.astype("string").str.strip().str.lower()
    has_dl = dl_raw.notna() & ~dl.isin(["", "nan", "none", "null", "0", "false", "no"])
    n_dl = int(has_dl.sum())
    qc = qc[~has_dl]
    QC_LOG["rule6_detection_limit_removed"] = n_dl

    # rule 5: flag IQR outliers (do not delete)
    qc["iqr_outlier_flag"] = qc["obs_iqr_outlier"].isin(["yes", "1", "true", "y"])

    # rule 8: daily aggregation
    qc["obs_time"] = qc["obs_time"].astype(str).replace({"nan": None})
    qc["_day"] = qc["obs_date"].dt.normalize()
    agg = (
        qc.sort_values(["site_id", "_day"])
        .groupby(["site_id", "_day"], as_index=False)
        .agg(
            tw_obs=("obs_value", "mean"),
            n_samples=("obs_value", "size"),
            obs_value_sd=("obs_value", "std"),
            obs_time=("obs_time", "first"),
            site_name=("site_name", "first"),
            site_country=("site_country", "first"),
            lat=("lat", "first"),
            lon=("lon", "first"),
            drainage_region_name=("drainage_region_name", "first"),
            upstream_basin_area=("upstream_basin_area", "first"),
            upstream_basin_area_unit=("upstream_basin_area_unit", "first"),
            source=("source", lambda s: "|".join(sorted(set(s.dropna().astype(str))))),
            param_form=("param_form", "first"),
            obs_iqr_outlier=("obs_iqr_outlier", "first"),
            iqr_outlier_flag=("iqr_outlier_flag", "any"),
            site_ts_availability=("site_ts_availability", "first"),
            site_ts_continuity=("site_ts_continuity", "first"),
        )
        .rename(columns={"_day": "obs_date"})
    )
    QC_LOG["n_daily_obs_after_aggregation"] = int(len(agg))
    QC_LOG["n_rows_collapsed_by_aggregation"] = int(len(qc) - len(agg))

    # rule 7: per-site >= 12 usable observations
    counts = agg.groupby("site_id").size()
    small = set(counts[counts < 12].index)
    n_small_rows = int(agg["site_id"].isin(small).sum())
    agg = agg[~agg["site_id"].isin(small)]
    QC_LOG["rule7_sites_dropped_lt12"] = int(len(small))
    QC_LOG["rule7_obs_removed_lt12"] = int(n_small_rows)
    QC_LOG["sites_dropped_lt12_ids"] = sorted(small)

    # rule 9: MAD outliers vs site-month median (flag only)
    med = agg.groupby(["site_id", agg["obs_date"].dt.month])["tw_obs"].transform("median")
    resid = agg["tw_obs"] - med
    mad = resid.abs().groupby([agg["site_id"], agg["obs_date"].dt.month]).transform("median")
    mad = mad.replace(0, np.nan)
    agg["mad_outlier_flag"] = (resid.abs() > 4 * mad).fillna(False)
    QC_LOG["rule9_mad_outlier_flagged"] = int(agg["mad_outlier_flag"].sum())

    # rule 10: continuity / availability distribution
    QC_LOG["rule10_site_ts_continuity"] = (
        agg.groupby("site_id")["site_ts_continuity"].first().astype(str)
        .value_counts().head(15).to_dict()
    )
    QC_LOG["rule10_site_ts_availability"] = (
        agg.groupby("site_id")["site_ts_availability"].first().astype(str)
        .value_counts().head(15).to_dict()
    )

    agg["month"] = agg["obs_date"].dt.month
    agg["year"] = agg["obs_date"].dt.year
    agg.to_csv(C.INTERIM / "water_temp_qc.csv", index=False)

    # ---------- domain extracts (raw counts, rules 1-2 only) ----------
    def raw_extract(df: pd.DataFrame) -> pd.DataFrame:
        return df.copy()

    # Spain (all retained raw rows; bbox is Iberia + Balearics via country)
    spain = raw_extract(raw)
    spain.to_csv(C.INTERIM / "grqa_temp_spain.csv", index=False)

    # Guadalquivir: clip by official ES050 water-body catchments, plus the
    # ES050 (Guadalquivir district) Waterbase site-id prefix, plus name rules.
    gpts = gpd.GeoDataFrame(
        raw.copy(), geometry=gpd.points_from_xy(raw["lon"], raw["lat"]), crs=4326
    )
    inside = gpts.within(basin_wgs).values
    sid = raw["site_id"].astype(str)
    name_hit = raw["site_name"].astype(str).str.contains(GQ_NAME_RE, case=False, na=False)
    drain_hit = raw["drainage_region_name"].astype(str).str.contains(
        "guadalquivir", case=False, na=False
    )
    is_gq = inside | sid.str.upper().str.startswith("ES050") | name_hit | drain_hit
    gq = raw[is_gq].copy()
    gq.to_csv(C.INTERIM / "grqa_temp_guadalquivir.csv", index=False)
    QC_LOG["gq_match_polygon"] = int(inside.sum())
    QC_LOG["gq_match_siteid_es050"] = int(sid.str.upper().str.startswith("ES050").sum())
    QC_LOG["gq_match_name"] = int(name_hit.sum())
    QC_LOG["gq_match_drainage"] = int(drain_hit.sum())

    # ---------- verification against the known ground truth ----------
    truth = {
        "ESP00016": 67, "ESP00017": 58, "ESP00018": 68,
        "ES050ESPF10705": 1, "ES050ESPF40901": 1,
    }
    summ = (
        gq.groupby(["site_id", "site_name", "source"]).agg(
            n_obs=("obs_id", "size"),
            first=("obs_date", "min"), last=("obs_date", "max"),
            lat=("lat", "first"), lon=("lon", "first"),
        ).reset_index()
    )
    print("\n=== Guadalquivir raw extract ===")
    print(summ.to_string(index=False))
    total_gq = int(len(gq))
    ok = total_gq == sum(truth.values()) and set(summ["site_id"]) == set(truth)
    print(f"ground-truth total rows: {sum(truth.values())}, got {total_gq}, match={ok}")

    QC_LOG["gq_raw_rows"] = total_gq
    QC_LOG["gq_raw_sites"] = int(summ["site_id"].nunique())
    QC_LOG["gq_matches_ground_truth"] = bool(ok)
    QC_LOG["spain_raw_rows"] = int(len(spain))
    QC_LOG["spain_raw_sites"] = int(spain["site_id"].nunique())

    # bbox cross-check: how many sites/rows the naive bbox pulls in (all countries)
    bbox_all = pd.concat(bbox_parts, ignore_index=True) if bbox_parts else pd.DataFrame()
    bbox_all.to_csv(C.INTERIM / "grqa_bbox_anycountry_temp.csv", index=False)
    QC_LOG["bbox_rows"] = int(len(bbox_all))
    QC_LOG["bbox_sites"] = int(bbox_all["site_id"].nunique())
    QC_LOG["bbox_sites_excluded_not_guadalquivir"] = sorted(
        set(bbox_all["site_id"]) - set(gq["site_id"])
    )

    (C.LOGS / "qc_log.json").write_text(json.dumps(QC_LOG, indent=2, default=str), encoding="utf-8")
    # save verification table
    summ.to_csv(C.TABLES / "grqa_guadalquivir_verification.csv", index=False)

    C.record("Generated: water_temp_raw/qc, grqa_temp_spain/guadalquivir", status="SUCCESS",
             path=str(C.INTERIM / "water_temp_qc.csv"), checksum=False,
             notes="generated by 02_process_grqa.py from GRQA v1.3 TEMP_GRQA.csv")
    C.record("CHG ES050 water-body catchments (repo copy)",
             path=str(BASIN_SHP),
             notes="used to clip Guadalquivir; CHG WFS GetFeature returned HTTP 401",
             license="CHG IDE / public", checksum=True)
    print("\nQC log:", json.dumps(QC_LOG, indent=2, default=str)[:2500])


if __name__ == "__main__":
    main()
