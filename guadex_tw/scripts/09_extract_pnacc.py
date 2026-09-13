"""Extract PNACC (AEMET/AdapteCCa) daily tmean at Guadalquivir sites, all GCMs/scenarios.

The 2 GB NetCDF-4 files are opened remotely over HTTP (fsspec + h5py) and only
the ~5 site grid cells are read, so nothing large is ever downloaded. THREDDS
NCSS is broken for these files (server-side "strict nc3" netCDF-4 error), which
is recorded in PROVENANCE.md; remote HDF5 random access is the working method.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import fsspec
import h5py
import numpy as np
import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C

BASE = ("https://escenarios.adaptecca.es/thredds/fileServer/peninsula/"
        "Proyecciones_CMIP6_en_rejilla/Dato_diario/Temperatura/tmean")
FILE = "tmean_SP-005_{gcm}_{exp}_ESD-RegBA_day.nc"
GCMS = ["ACCESS-CM2", "CMCC-CM2-SR5", "CNRM-ESM2-1", "EC-Earth3-Veg", "IITM-ESM",
        "KACE-1-0-G", "MIROC6", "MPI-ESM1-2-HR", "MRI-ESM2-0", "NorESM2-MM",
        "UKESM1-0-LL"]
SCENARIOS = ["ssp126", "ssp245", "ssp370", "ssp585"]
LICENSE = "AEMET/PNACC 2024 - public (attribution)"
CACHE = C.INTERIM / "pnacc_cache"
CACHE.mkdir(parents=True, exist_ok=True)
OUT = C.INTERIM / "pnacc_tmean_guadalquivir.parquet"


CORE_GEMSTAT = ["ESP00016", "ESP00017", "ESP00018"]
MAX_SITES = 10
# year windows actually needed (baseline 1986-2005; periods 2021-2040/2036-2055/2041-2070)
YEARS_HIST = (1985, 2006)
YEARS_FUT = (2020, 2071)


def gq_sites() -> pd.DataFrame:
    """Elevation-spanning Guadalquivir site sample.

    Selection is stratified so the projection actually samples the elevation
    gradient the slope interaction is meant to capture: the three multi-year
    GEMSTAT reaches, then the lowest and highest sites, then the densest
    remaining records. High-elevation picks are added *before* truncation to
    MAX_SITES so they cannot be dropped as they were previously.
    """
    s = pd.read_csv(C.PROCESSED / "sites.csv", dtype={"site_id": str})
    gq = s[s["in_guadalquivir_basin"]].dropna(subset=["lat", "lon", "elevation_m"]).copy()
    gq["n_obs_qc"] = pd.to_numeric(gq["n_obs_qc"], errors="coerce").fillna(0)
    by_z = gq.sort_values("elevation_m")
    picks = [
        gq[gq["site_id"].isin(CORE_GEMSTAT)],     # only multi-year basin series
        by_z.head(2),                             # lowest reaches
        by_z.tail(2),                             # highest headwater reaches
        gq.sort_values("n_obs_qc", ascending=False),  # densest records
    ]
    pick = pd.concat(picks, ignore_index=True).drop_duplicates(subset=["site_id"])
    pick = pick.head(MAX_SITES)
    pick["lat_r"] = pick["lat"].round(3)
    pick["lon_r"] = pick["lon"].round(3)
    return pick.drop_duplicates(subset=["lat_r", "lon_r"]).reset_index(drop=True)


def extract_one(gcm: str, exp: str, sites: pd.DataFrame) -> pd.DataFrame:
    cached = CACHE / f"{gcm}_{exp}.csv"
    want = set(sites["site_id"].astype(str))
    if cached.exists():
        c = pd.read_csv(cached, parse_dates=["date"], dtype={"site_id": str})
        if set(c["site_id"].astype(str)) == want:
            return c
        print(f"    cache {cached.name} has a different site set; re-extracting")
    fname = FILE.format(gcm=gcm, exp=exp)
    url = f"{BASE}/{fname}"
    import aiohttp
    fs = fsspec.filesystem(
        "http",
        client_kwargs={"timeout": aiohttp.ClientTimeout(total=600, sock_read=90)},
    )
    f = fs.open(url, "rb")
    h = h5py.File(f, "r")
    lat = h["lat"][:]
    lon = h["lon"][:]
    tunits = h["time"].attrs["units"].decode() if isinstance(h["time"].attrs["units"], bytes) else str(h["time"].attrs["units"])
    tref = pd.Timestamp(tunits.split("since")[1].strip())
    t = pd.to_datetime(tref) + pd.to_timedelta(h["time"][:], unit="D")
    y0, y1 = YEARS_HIST if exp == "historical" else YEARS_FUT
    tmask = (t.year >= y0) & (t.year <= y1)
    tidx = np.where(tmask)[0]
    t = t[tidx]
    rows = []
    for _, s in sites.iterrows():
        i = int(np.argmin(np.abs(lat - s["lat"])))
        j = int(np.argmin(np.abs(lon - s["lon"])))
        vals = h["tmean"][int(tidx[0]):int(tidx[-1]) + 1, i, j]
        rows.append(pd.DataFrame({
            "gcm": gcm, "experiment": exp, "site_id": s["site_id"],
            "date": t, "tmean": vals,
            "cell_lat": lat[i], "cell_lon": lon[j],
        }))
    h.close()
    f.close()
    df = pd.concat(rows, ignore_index=True)
    df = df.drop_duplicates(subset=["site_id", "date"])
    df.to_csv(cached, index=False)
    return df


def combine(sites: pd.DataFrame) -> pd.DataFrame | None:
    """Combine the per-(GCM, experiment) cache files into the parquet output.

    Refuses to overwrite the output unless every expected (GCM, experiment) file
    exists *and* every requested site is present in the cache, so a partial or
    stale cache can never replace a good parquet with a truncated one.
    """
    expected = [(g, e) for g in GCMS for e in ["historical"] + SCENARIOS]
    missing_files = [f"{g}_{e}" for g, e in expected if not (CACHE / f"{g}_{e}.csv").exists()]
    want = set(sites["site_id"].astype(str))
    frames = []
    for gcm, exp in expected:
        p = CACHE / f"{gcm}_{exp}.csv"
        if p.exists():
            frames.append(pd.read_csv(p, parse_dates=["date"], dtype={"site_id": str}))
    if not frames:
        print("no cache files to combine")
        return None
    all_df = pd.concat(frames, ignore_index=True)
    if missing_files:
        print(f"REFUSING to overwrite {OUT}: missing extraction files {missing_files[:5]}"
              f"{' ...' if len(missing_files) > 5 else ''}")
        return None
    # Per-file completeness: every cache file must contain every requested site.
    # A union check is not enough because a partially refreshed cache can hold
    # all sites across files while individual GCMs are missing high elevations.
    incomplete = []
    for gcm, exp in expected:
        p = CACHE / f"{gcm}_{exp}.csv"
        have = set(pd.read_csv(p, usecols=["site_id"], dtype={"site_id": str})["site_id"]
                   .astype(str))
        miss = want - have
        if miss:
            incomplete.append(f"{gcm}_{exp}({len(miss)} missing)")
    if incomplete:
        print(f"REFUSING to overwrite {OUT}: per-file site coverage incomplete for "
              f"{incomplete[:5]}{' ...' if len(incomplete) > 5 else ''}")
        return None
    all_df = all_df[all_df["site_id"].astype(str).isin(want)]
    all_df.to_parquet(OUT, index=False)
    n_gcm = all_df["gcm"].nunique()
    print(f"wrote {OUT} ({len(all_df):,} rows; {n_gcm} GCMs; "
          f"{all_df['site_id'].nunique()} sites)")
    C.record("PNACC tmean extracted at Guadalquivir sites", url=BASE, path=str(OUT),
             license=LICENSE, version="PNACC 2024 CMIP6 5km ESD-RegBA",
             notes=f"{n_gcm} GCMs x {len(SCENARIOS)+1} experiments; "
                   f"{all_df['site_id'].nunique()} sites; remote HDF5 subset (NCSS broken)")
    return all_df


def write_summary(sites: pd.DataFrame, t0: float, gcms: list) -> None:
    summary = {
        "n_gcm": len(gcms), "gcms": gcms, "scenarios": SCENARIOS,
        "sites": list(sites["site_id"]),
        "site_elevations_m": {str(r["site_id"]): float(r["elevation_m"])
                              for _, r in sites.iterrows()},
        "elevation_range_m": [float(sites["elevation_m"].min()),
                              float(sites["elevation_m"].max())],
        "seconds": round(time.time() - t0, 1),
    }
    (C.LOGS / "pnacc_extract_summary.json").write_text(json.dumps(summary, indent=2),
                                                       encoding="utf-8")


def main() -> None:
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--gcms", default="",
                    help="comma-separated GCM subset to extract (default: all)")
    ap.add_argument("--combine-only", action="store_true",
                    help="only combine existing cache files into the parquet")
    ap.add_argument("--max-seconds", type=float, default=7200.0,
                    help="stop extracting after this wall-clock budget (default 7200)")
    args = ap.parse_args()

    sites = gq_sites()
    print(f"{len(sites)} Guadalquivir sites: {list(sites['site_id'])}")
    t0 = time.time()
    if args.combine_only:
        combine(sites)
        write_summary(sites, t0, GCMS)
        return

    gcms = [g.strip() for g in args.gcms.split(",") if g.strip()] or GCMS
    budget_hit = False
    for gcm in gcms:
        for exp in ["historical"] + SCENARIOS:
            if time.time() - t0 > args.max_seconds:
                budget_hit = True
                print(f"wall-clock budget ({args.max_seconds:.0f}s) reached; "
                      f"stopping with cache as-is", flush=True)
                break
            d, last = None, None
            for attempt in range(1, 3):
                try:
                    d = extract_one(gcm, exp, sites)
                    break
                except Exception as exc:  # noqa: BLE001
                    last = exc
                    print(f"  retry {attempt}/2 {gcm} {exp}: "
                          f"{type(exc).__name__}: {exc}", flush=True)
                    time.sleep(5 * attempt)
            if d is not None:
                print(f"  {gcm:16s} {exp:11s} {len(d):6d} rows "
                      f"({time.time()-t0:.0f}s)", flush=True)
            else:
                print(f"  FAIL {gcm} {exp}: {last!r}", flush=True)
                C.record(f"PNACC {gcm} {exp}", url=f"{BASE}/{FILE.format(gcm=gcm, exp=exp)}",
                         status="FAILED", error=repr(last), license=LICENSE,
                         version="PNACC 2024 CMIP6 5km ESD-RegBA")
        if budget_hit:
            break
    if set(gcms) == set(GCMS):
        combine(sites)
    else:
        print("subset extraction complete; run 09_extract_pnacc.py --combine-only "
              "once every GCM is done")
    write_summary(sites, t0, gcms)


if __name__ == "__main__":
    main()
