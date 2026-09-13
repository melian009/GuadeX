"""Build data/processed/sites.csv: one row per water-temperature site.

Joins GRQA / Waterbase sites with Open-Meteo elevation and Copernicus DEM GLO-30
(3x3 median) cross-check, basin membership, river/reservoir classification and
observation cadence statistics.
"""
from __future__ import annotations

import json
import math
import os
import re
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import Window
from shapely.geometry import Point
from shapely.ops import unary_union

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C

BASIN_SHP = C.ROOT.parent / "data" / "Version_02-01-2026-Masas" / "Cuencas_masas_agua_4c.shp"
DEM_CACHE = C.INTERIM / "dem_samples.csv"

RESERVOIR_RE = re.compile(
    r"embalse|presa|reservoir|azud|pantano|izn[aá]jar|tranco|negrat[ií]n|rumblar|"
    r"puente nuevo|guadal[eé]n|giribaile|bemb[eé]zar|j[aá]ndula|"
    r"embassament|barragem", re.I)
CANAL_RE = re.compile(r"canal|acequia|canalizad", re.I)
ESTUARY_RE = re.compile(r"estuario|estuary|desembocadura|marisma|do[nñ]ana|pluma|"
                        r"bah[ií]a|port|puerto", re.I)


def dem_tile(lat: float, lon: float) -> str:
    la, lo = math.floor(lat), math.floor(lon)
    ns = "N" if la >= 0 else "S"
    ew = "E" if lo >= 0 else "W"
    return f"Copernicus_DSM_COG_10_{ns}{abs(la):02d}_00_{ew}{abs(lo):03d}_00_DEM"


def dem_sample(lat: float, lon: float) -> tuple[float, float]:
    """Return (centre_cell, 3x3 median) elevation from Copernicus DEM GLO-30."""
    t = dem_tile(lat, lon)
    url = f"https://copernicus-dem-30m.s3.amazonaws.com/{t}/{t}.tif"
    with rasterio.open("/vsicurl/" + url) as src:
        r, c = src.index(lon, lat)
        win = Window(c - 1, r - 1, 3, 3)
        arr = src.read(1, window=win, boundless=True, fill_value=np.nan)
    arr = arr.astype(float)
    if arr.size == 0:
        return np.nan, np.nan
    center = float(arr[arr.shape[0] // 2, arr.shape[1] // 2])
    med = float(np.nanmedian(arr)) if np.isfinite(arr).any() else np.nan
    return center, med


def load_sites() -> pd.DataFrame:
    frames = []
    qc = C.PROCESSED / "water_temp_all_qc.csv"
    if qc.exists():
        d = pd.read_csv(qc, parse_dates=["obs_date"], dtype={"site_id": str})
        if "site_class_ext" not in d:
            d["site_class_ext"] = np.nan
        g = d.groupby("site_id", as_index=False).agg(
            site_name=("site_name", "first"), lat=("lat", "first"), lon=("lon", "first"),
            source=("source", "first"), n_obs_qc=("obs_date", "size"),
            first_date=("obs_date", "min"), last_date=("obs_date", "max"),
            wb_reservoir=("site_class_ext", "first"),
        )
        frames.append(g)
    raw = C.INTERIM / "water_temp_raw.csv"
    if raw.exists():
        d = pd.read_csv(raw, parse_dates=["obs_date"], dtype={"site_id": str},
                        usecols=["site_id", "site_name", "lat", "lon", "source", "obs_date"])
        g = d.groupby("site_id", as_index=False).agg(
            site_name=("site_name", "first"), lat=("lat", "first"), lon=("lon", "first"),
            source=("source", "first"), n_obs_raw=("obs_date", "size"),
        )
        frames.append(g)
    # raw Guadalquivir extract so the single-obs Waterbase points appear in the table
    gqr = C.INTERIM / "grqa_temp_guadalquivir.csv"
    if gqr.exists():
        d = pd.read_csv(gqr, parse_dates=["obs_date"], dtype={"site_id": str})
        g = d.groupby("site_id", as_index=False).agg(
            site_name=("site_name", "first"), lat=("lat", "first"), lon=("lon", "first"),
            source=("source", "first"), n_obs_raw=("obs_date", "size"),
        )
        frames.append(g)
    s = pd.concat(frames, ignore_index=True)
    s = s.groupby("site_id", as_index=False).agg(
        site_name=("site_name", "first"), lat=("lat", "first"), lon=("lon", "first"),
        source=("source", "first"),
        n_obs_raw=("n_obs_raw", "max"), n_obs_qc=("n_obs_qc", "max"),
        first_date=("first_date", "min"), last_date=("last_date", "max"),
        wb_reservoir=("wb_reservoir", "first"),
    )
    for c in ("n_obs_raw", "n_obs_qc"):
        if c not in s.columns:
            s[c] = np.nan
    return s


def main() -> None:
    print("Loading site list ...")
    sites = load_sites()
    print(f"{len(sites)} unique sites")

    # cadence
    dailies = []
    p = C.PROCESSED / "water_temp_all_qc.csv"
    if p.exists():
        dailies.append(pd.read_csv(p, parse_dates=["obs_date"], dtype={"site_id": str}))
    d = pd.concat(dailies, ignore_index=True)
    by_year = d.groupby(["site_id", d["obs_date"].dt.year]).size()
    med_per_year = by_year.groupby("site_id").median().rename("median_obs_per_year")
    cad = d.sort_values(["site_id", "obs_date"]).groupby("site_id")["obs_date"].diff().dt.days
    cadence = cad.groupby(d["site_id"]).median().rename("obs_cadence_days_median")
    sites = sites.merge(med_per_year, on="site_id", how="left")
    sites = sites.merge(cadence, on="site_id", how="left")

    # Open-Meteo elevation
    oe = C.INTERIM / "openmeteo_elevation.csv"
    if oe.exists():
        sites = sites.merge(pd.read_csv(oe, dtype={"site_id": str})[
                                ["site_id", "openmeteo_elevation_m"]],
                            on="site_id", how="left")
    else:
        sites["openmeteo_elevation_m"] = np.nan
    sites["elevation_m"] = sites["openmeteo_elevation_m"]

    # Copernicus DEM sampling (cached, resumable)
    if DEM_CACHE.exists():
        cache = pd.read_csv(DEM_CACHE, dtype={"site_id": str})
    else:
        cache = pd.DataFrame(columns=["site_id", "dem_elevation_m", "dem_center_m"])
    done = set(cache["site_id"])
    new_rows = []
    t = C.Throttle(0.0)
    for _, row in sites.iterrows():
        sid = str(row["site_id"])
        if sid in done:
            continue
        try:
            center, med = dem_sample(float(row["lat"]), float(row["lon"]))
            new_rows.append({"site_id": sid, "dem_elevation_m": med, "dem_center_m": center})
        except Exception as exc:  # noqa: BLE001
            new_rows.append({"site_id": sid, "dem_elevation_m": np.nan, "dem_center_m": np.nan,
                             "dem_error": repr(exc)})
        if len(new_rows) % 20 == 0 and new_rows:
            print(f"  DEM {len(new_rows)} sampled ...", flush=True)
            pd.concat([cache, pd.DataFrame(new_rows)], ignore_index=True).to_csv(DEM_CACHE, index=False)
    if new_rows:
        cache = pd.concat([cache, pd.DataFrame(new_rows)], ignore_index=True)
        cache.to_csv(DEM_CACHE, index=False)
    sites = sites.merge(cache[["site_id", "dem_elevation_m", "dem_center_m"]],
                        on="site_id", how="left")

    # dual-source reconciliation
    sites["elevation_diff_m"] = sites["elevation_m"] - sites["dem_elevation_m"]
    sites["elevation_discrepancy_flag"] = sites["elevation_diff_m"].abs() > 50
    sites["elevation_source"] = np.where(
        sites["dem_elevation_m"].notna() & sites["openmeteo_elevation_m"].notna(),
        "Copernicus_DEM_GLO30(median)|OpenMeteo_DEM90(site)",
        np.where(sites["dem_elevation_m"].notna(), "Copernicus_DEM_GLO30(median)",
                 np.where(sites["openmeteo_elevation_m"].notna(), "OpenMeteo_DEM90(site)", "FAILED")),
    )
    # prefer DEM median as the elevation_m used by models, else Open-Meteo
    sites["elevation_m"] = np.where(
        sites["dem_elevation_m"].notna(), sites["dem_elevation_m"],
        sites["openmeteo_elevation_m"])

    # basin membership
    g = gpd.read_file(BASIN_SHP)
    es050 = g[g["EUMASCod"].astype(str).str.startswith("ES050")]
    poly = unary_union(es050.geometry.values)
    pts = gpd.GeoSeries([Point(x, y) for x, y in zip(sites["lon"], sites["lat"])], crs=4326)
    pts_utm = pts.to_crs(25830)
    sites["in_guadalquivir_basin"] = [bool(poly.contains(p)) for p in pts_utm]
    sites["basin_name"] = np.where(sites["in_guadalquivir_basin"], "Guadalquivir (ES050)", "")

    # classification
    def classify(row) -> str:
        name = str(row.get("site_name", ""))
        src = str(row.get("source", ""))
        if ESTUARY_RE.search(name):
            return "estuary"
        if RESERVOIR_RE.search(name):
            return "reservoir"
        if CANAL_RE.search(name):
            return "canal"
        if src == "WATERBASE_DIRECT" and str(row.get("wb_reservoir", "")).lower() in ("reservoir", "true", "1"):
            return "reservoir"
        if src == "WATERBASE_DIRECT" and str(row.get("wb_reservoir", "")).lower() in ("river", "false", "0"):
            return "river"
        return "river" if name else "unknown"

    sites["site_class"] = sites.apply(classify, axis=1)

    # approximate nearest grid cells (documented as index-cell approximations)
    sites["nearest_rocio_cell_lat"] = (np.round(sites["lat"] / 0.05) * 0.05).round(4)
    sites["nearest_rocio_cell_lon"] = (np.round(sites["lon"] / 0.05) * 0.05).round(4)
    sites["nearest_era5_land_cell_lat"] = (np.round(sites["lat"] / 0.1) * 0.1).round(3)
    sites["nearest_era5_land_cell_lon"] = (np.round(sites["lon"] / 0.1) * 0.1).round(3)
    sites["gridcell_elevation_m"] = sites["dem_elevation_m"]

    # nearest AEMET station (only if API key available)
    aemet_key = C.env_key("AEMET_API_KEY")
    sites["nearest_aemet_idema"] = ""
    sites["nearest_aemet_name"] = ""
    sites["nearest_aemet_distance_km"] = np.nan
    sites["aemet_elevation_m"] = np.nan
    if aemet_key:
        try:
            base = "https://opendata.aemet.es/opendata/api"
            r = requests_get_json(f"{base}/valores/climatologicos/inventarioestaciones/todasestaciones?api_key={aemet_key}")
            inv = pd.DataFrame(r)
            inv["lat"] = pd.to_numeric(inv["latitud"], errors="coerce")
            inv["lon"] = pd.to_numeric(inv["longitud"], errors="coerce")
            inv["alt"] = pd.to_numeric(inv["altitud"], errors="coerce")
            for i, row in sites.iterrows():
                dd = np.hypot((inv["lat"] - row["lat"]) * 111.0,
                              (inv["lon"] - row["lon"]) * 111.0 * math.cos(math.radians(row["lat"])))
                j = int(np.nanargmin(dd.values))
                sites.at[i, "nearest_aemet_idema"] = inv.iloc[j]["indicativo"]
                sites.at[i, "nearest_aemet_name"] = inv.iloc[j]["nombre"]
                sites.at[i, "nearest_aemet_distance_km"] = round(float(dd.iloc[j]), 2)
                sites.at[i, "aemet_elevation_m"] = inv.iloc[j]["alt"]
            C.record("AEMET station inventory", url=f"{base}/valores/climatologicos/inventarioestaciones/todasestaciones",
                     status="SUCCESS", license="AEMET OpenData",
                     notes=f"{len(inv)} stations, used for nearest-station column")
        except Exception as exc:  # noqa: BLE001
            print(f"AEMET inventory failed: {exc!r}")
            C.record("AEMET station inventory", url="https://opendata.aemet.es/", status="FAILED",
                     error=repr(exc), license="AEMET OpenData")
    else:
        print("AEMET_API_KEY not set; nearest_aemet_* left blank")
        C.record("AEMET station inventory", url="https://opendata.aemet.es/",
                 status="SKIPPED", error="AEMET_API_KEY not set",
                 notes="no-auth path sufficient; AEMET requires a key")

    cols = ["site_id", "site_name", "source", "lat", "lon", "elevation_m", "dem_elevation_m",
            "openmeteo_elevation_m", "elevation_source", "elevation_diff_m",
            "elevation_discrepancy_flag", "basin_name", "in_guadalquivir_basin", "site_class",
            "n_obs_raw", "n_obs_qc", "first_date", "last_date", "median_obs_per_year",
            "obs_cadence_days_median", "nearest_aemet_idema", "nearest_aemet_name",
            "nearest_aemet_distance_km", "aemet_elevation_m", "nearest_rocio_cell_lat",
            "nearest_rocio_cell_lon", "nearest_era5_land_cell_lat",
            "nearest_era5_land_cell_lon", "gridcell_elevation_m"]
    out = sites[cols].sort_values("site_id")
    out.to_csv(C.PROCESSED / "sites.csv", index=False)
    print(f"wrote {C.PROCESSED / 'sites.csv'} ({len(out)} sites)")

    gq = out[out["in_guadalquivir_basin"]]
    summary = {
        "n_sites": int(len(out)),
        "n_sites_guadalquivir": int(len(gq)),
        "elevation_range_all_m": [float(out["elevation_m"].min()), float(out["elevation_m"].max())],
        "elevation_range_guadalquivir_m": [float(gq["elevation_m"].min()), float(gq["elevation_m"].max())],
        "n_elevation_discrepancy_gt50m": int(out["elevation_discrepancy_flag"].sum()),
        "n_missing_elevation": int(out["elevation_m"].isna().sum()),
        "site_class_counts": out["site_class"].value_counts().to_dict(),
        "source_counts": out["source"].value_counts().to_dict(),
    }
    (C.LOGS / "site_summary.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    C.record("Generated sites.csv", path=str(C.PROCESSED / "sites.csv"), status="SUCCESS",
             notes=json.dumps(summary, default=str))
    print(json.dumps(summary, indent=2, default=str))


def requests_get_json(url: str):
    import requests
    r = requests.get(url, timeout=120, headers={"User-Agent": C.USER_AGENT})
    r.raise_for_status()
    j = r.json()
    if isinstance(j, dict) and "datos" in j:
        r2 = requests.get(j["datos"], timeout=120, headers={"User-Agent": C.USER_AGENT})
        r2.raise_for_status()
        return r2.json()
    return j


if __name__ == "__main__":
    main()
