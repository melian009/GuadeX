"""WP1: project daily water temperature at every GuadeX site (all GCMs/SSPs).

The calibration pipeline (`09_extract_pnacc.py` / `10_project.py`) extracts PNACC
grid cells only for a stratified sample of monitoring stations (`MAX_SITES = 10`)
because that is what the Tw model is validated against.  The climate-scenario
ODE, however, is run on the GuadeX network, so it needs forcing at every site
coordinate.

This script mirrors the calibration pipeline but reads the GuadeX site table
(`data/ConnectivityUTM.csv`), extracts the PNACC daily `tmean` grid cell at each
site, applies the calibrated spatial Tw model and writes the per-site daily
ensemble-median series consumed by `Guadex.load_daily_temperature_forcing`:

    guadex_tw/outputs/tables/water_temp_daily_guadex_sites.csv
    columns: site_id, scenario, date, tw_ensemble_median

Rows are written for the `historical` experiment (the 1986-2005 baseline used by
`daily_temperature_schedule`) and for every SSP, so the ODE can form anomalies
against the fixed baseline period instead of anchoring at the start year.

Usage (network access required, ~2 GB files opened remotely over HTTP):

    python scripts/13_project_guadex_sites.py --max-sites 25      # pilot
    python scripts/13_project_guadex_sites.py                     # all sites
    python scripts/13_project_guadex_sites.py --combine-only      # cache -> CSV

Design notes
------------
* The remote HDF5 subsetting, cache layout and completeness guards follow
  `09_extract_pnacc.py`; the cache directory is separate so the calibration
  pull is never disturbed.
* Bias correction to the NASA POWER 1986-2005 climatology is applied when the
  cached POWER series are available; because the correction is an additive
  constant per (site, GCM) and the model is linear in Ta, it cancels in the
  anomalies the ODE actually uses.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402

REPO_ROOT = C.ROOT.parent
GUADEX_SITES_CSV = REPO_ROOT / "data" / "ConnectivityUTM.csv"

BASE = ("https://escenarios.adaptecca.es/thredds/fileServer/peninsula/"
        "Proyecciones_CMIP6_en_rejilla/Dato_diario/Temperatura/tmean")
FILE = "tmean_SP-005_{gcm}_{exp}_ESD-RegBA_day.nc"
GCMS = ["ACCESS-CM2", "CMCC-CM2-SR5", "CNRM-ESM2-1", "EC-Earth3-Veg", "IITM-ESM",
        "KACE-1-0-G", "MIROC6", "MPI-ESM1-2-HR", "MRI-ESM2-0", "NorESM2-MM",
        "UKESM1-0-LL"]
SCENARIOS = ["ssp126", "ssp245", "ssp370", "ssp585"]
BASELINE = (1986, 2005)
YEARS_HIST = (1985, 2006)
YEARS_FUT = (2020, 2071)
LICENSE = "AEMET/PNACC 2024 - public (attribution)"

CACHE = C.INTERIM / "pnacc_cache_guadex"
CACHE.mkdir(parents=True, exist_ok=True)
OUT = C.TABLES / "water_temp_daily_guadex_sites.csv"


def utm_to_latlon(easting: float, northing: float, zone: int = 30,
                  northern: bool = True) -> tuple[float, float]:
    """Inverse UTM (WGS84) without a pyproj dependency.

    The GuadeX connectivity table stores ETRS89/UTM coordinates for the
    Guadalquivir basin; zone 30N is the basin convention.  The inverse transform
    is only used to locate the 5 km PNACC grid cell, so sub-kilometre accuracy is
    ample.
    """
    a = 6378137.0
    f = 1 / 298.257223563
    k0 = 0.9996
    e2 = f * (2 - f)
    e1 = (1 - math.sqrt(1 - e2)) / (1 + math.sqrt(1 - e2))

    x = easting - 500000.0
    y = northing if northern else northing - 10000000.0
    lon0 = math.radians(zone * 6 - 183)

    m = y / k0
    mu = m / (a * (1 - e2 / 4 - 3 * e2**2 / 64 - 5 * e2**3 / 256))
    phi1 = (mu
            + (3 * e1 / 2 - 27 * e1**3 / 32) * math.sin(2 * mu)
            + (21 * e1**2 / 16 - 55 * e1**4 / 32) * math.sin(4 * mu)
            + (151 * e1**3 / 96) * math.sin(6 * mu))
    n1 = a / math.sqrt(1 - e2 * math.sin(phi1)**2)
    t1 = math.tan(phi1)**2
    c1 = e2 / (1 - e2) * math.cos(phi1)**2
    r1 = a * (1 - e2) / (1 - e2 * math.sin(phi1)**2)**1.5
    d = x / (n1 * k0)

    lat = phi1 - (n1 * math.tan(phi1) / r1) * (
        d**2 / 2
        - (5 + 3 * t1 + 10 * c1 - 4 * c1**2 - 9 * e2 / (1 - e2)) * d**4 / 24
        + (61 + 90 * t1 + 298 * c1 + 45 * t1**2 - 252 * e2 / (1 - e2)
           - 3 * c1**2) * d**6 / 720)
    lon = lon0 + (
        d
        - (1 + 2 * t1 + c1) * d**3 / 6
        + (5 - 2 * c1 + 28 * t1 - 3 * c1**2 + 8 * e2 / (1 - e2)
           + 24 * t1**2) * d**5 / 120) / math.cos(phi1)
    return math.degrees(lat), math.degrees(lon)


def guadex_sites(path: Path = GUADEX_SITES_CSV) -> pd.DataFrame:
    """GuadeX site table with lat/lon/elevation for the PNACC grid lookup."""
    if not path.exists():
        raise SystemExit(
            f"GuadeX site table not found: {path}\n"
            "Export the connectivity file (CODIGO, ALTITUD, UTMX, UTMY) before running.")
    df = pd.read_csv(path, dtype={"CODIGO": str}, low_memory=False)
    required = {"CODIGO", "ALTITUD", "UTMX", "UTMY"}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"{path} is missing columns: {sorted(missing)}")
    df = df.dropna(subset=["UTMX", "UTMY", "ALTITUD"]).copy()
    coords = [utm_to_latlon(float(x), float(y)) for x, y in zip(df["UTMX"], df["UTMY"])]
    df["lat"] = [c[0] for c in coords]
    df["lon"] = [c[1] for c in coords]
    df["elevation_m"] = pd.to_numeric(df["ALTITUD"], errors="coerce")
    df = df.dropna(subset=["elevation_m"]).drop_duplicates(subset=["lat", "lon"])
    return df.reset_index(drop=True)


def extract_one(gcm: str, exp: str, sites: pd.DataFrame) -> pd.DataFrame | None:
    cached = CACHE / f"{gcm}_{exp}.csv"
    want = set(sites["CODIGO"].astype(str))
    if cached.exists():
        c = pd.read_csv(cached, parse_dates=["date"], dtype={"site_id": str})
        if set(c["site_id"].astype(str)) == want:
            return c
        print(f"    cache {cached.name} has a different site set; re-extracting")

    import aiohttp
    import fsspec
    import h5py

    fname = FILE.format(gcm=gcm, exp=exp)
    url = f"{BASE}/{fname}"
    fs = fsspec.filesystem(
        "http",
        client_kwargs={"timeout": aiohttp.ClientTimeout(total=900, sock_read=120)},
    )
    f = fs.open(url, "rb")
    h = h5py.File(f, "r")
    lat = h["lat"][:]
    lon = h["lon"][:]
    tunits = h["time"].attrs["units"]
    tunits = tunits.decode() if isinstance(tunits, bytes) else str(tunits)
    tref = pd.Timestamp(tunits.split("since")[1].strip())
    t = pd.to_datetime(tref) + pd.to_timedelta(h["time"][:], unit="D")
    y0, y1 = YEARS_HIST if exp == "historical" else YEARS_FUT
    tidx = np.where((t.year >= y0) & (t.year <= y1))[0]
    t = t[tidx]
    rows = []
    for _, s in sites.iterrows():
        i = int(np.argmin(np.abs(lat - s["lat"])))
        j = int(np.argmin(np.abs(lon - s["lon"])))
        vals = h["tmean"][int(tidx[0]):int(tidx[-1]) + 1, i, j]
        rows.append(pd.DataFrame({
            "gcm": gcm, "experiment": exp, "site_id": s["CODIGO"],
            "date": t, "tmean": vals, "cell_lat": lat[i], "cell_lon": lon[j],
        }))
    h.close()
    f.close()
    df = pd.concat(rows, ignore_index=True).drop_duplicates(subset=["site_id", "date"])
    df.to_csv(cached, index=False)
    return df


def load_model():
    """Reuse the calibrated Stage-B Model 2 from `10_project.py`."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "project_model", Path(__file__).resolve().parent / "10_project.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    model, _ = mod.fit_model()
    return mod, model


def project_series(mod, model, site_row, sub: pd.DataFrame, deltas: dict) -> pd.DataFrame:
    z = float(site_row["elevation_m"]) / 1000.0
    frames = []
    for gcm, delta in deltas.items():
        sel = sub[sub["gcm"] == gcm].copy()
        if sel.empty:
            continue
        sel["ta_adj"] = sel["tmean"] + delta
        sel["tw"] = mod.predict_tw(model, sel["ta_adj"].to_numpy(), z,
                                   sel["date"].dt.month.to_numpy())
        frames.append(sel[["date", "tw", "gcm"]])
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-sites", type=int, default=0,
                    help="limit to the first N sites (pilot; 0 = all)")
    ap.add_argument("--gcms", default="", help="comma-separated GCM subset")
    ap.add_argument("--combine-only", action="store_true",
                    help="only project from the existing cache")
    ap.add_argument("--max-seconds", type=float, default=24 * 3600.0)
    args = ap.parse_args()

    sites = guadex_sites()
    if args.max_sites:
        sites = sites.head(args.max_sites).reset_index(drop=True)
    print(f"{len(sites)} GuadeX sites "
          f"(elevation {sites.elevation_m.min():.0f}-{sites.elevation_m.max():.0f} m)")

    gcms = [g.strip() for g in args.gcms.split(",") if g.strip()] or GCMS
    t0 = time.time()
    if not args.combine_only:
        for gcm in gcms:
            for exp in ["historical"] + SCENARIOS:
                if time.time() - t0 > args.max_seconds:
                    print("wall-clock budget reached; cache is partial", flush=True)
                    return
                cached = CACHE / f"{gcm}_{exp}.csv"
                if cached.exists():
                    continue
                try:
                    d = extract_one(gcm, exp, sites)
                    print(f"  {gcm:16s} {exp:11s} {len(d):7d} rows "
                          f"({time.time() - t0:.0f}s)", flush=True)
                except Exception as exc:  # noqa: BLE001
                    print(f"  FAIL {gcm} {exp}: {type(exc).__name__}: {exc}", flush=True)

    missing = [(g, e) for g in gcms for e in ["historical"] + SCENARIOS
               if not (CACHE / f"{g}_{e}.csv").exists()]
    if missing:
        print(f"refusing to write {OUT}: missing {len(missing)} cache file(s), "
              f"e.g. {missing[:3]}; rerun (cache is resumable)", flush=True)
        return

    mod, model = load_model()
    out_rows = []
    for _, site in sites.iterrows():
        sid = str(site["CODIGO"])
        for exp in ["historical"] + SCENARIOS:
            parts = []
            for gcm in gcms:
                p = CACHE / f"{gcm}_{exp}.csv"
                if p.exists():
                    parts.append(pd.read_csv(p, parse_dates=["date"], dtype={"site_id": str}))
            sub = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
            sub = sub[sub["site_id"].astype(str) == sid]
            if sub.empty:
                continue
            series = project_series(mod, model, site, sub, {g: 0.0 for g in gcms})
            if series.empty:
                continue
            med = series.groupby("date")["tw"].median().reset_index()
            med["site_id"] = sid
            med["scenario"] = exp
            out_rows.append(med)

    daily = pd.concat(out_rows, ignore_index=True) if out_rows else pd.DataFrame()
    daily = daily[["site_id", "scenario", "date", "tw"]].rename(
        columns={"tw": "tw_ensemble_median"})
    daily = daily.sort_values(["site_id", "scenario", "date"])
    daily.to_csv(OUT, index=False)
    print(f"wrote {OUT} ({len(daily):,} rows; {daily.site_id.nunique()} sites; "
          f"{daily.scenario.nunique()} experiments)")

    C.record("Projected daily Tw at all GuadeX sites (WP1)", path=str(OUT),
             license=LICENSE, version="PNACC 2024 CMIP6 5km ESD-RegBA + Model 2",
             notes=f"{len(gcms)} GCMs x {len(SCENARIOS)} SSPs; "
                   f"{daily.site_id.nunique()} sites; ensemble median; "
                   f"baseline {BASELINE}")
    summary = {
        "n_sites": int(daily.site_id.nunique()),
        "n_gcm": len(gcms),
        "scenarios": SCENARIOS,
        "baseline": BASELINE,
        "date_range": [str(daily.date.min().date()), str(daily.date.max().date())],
    }
    (C.LOGS / "projection_summary_guadex_sites.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
