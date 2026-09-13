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


CORE = ["ESP00016", "ESP00017", "ESP00018", "ES050ESPF10705", "ES050ESPF40901"]
MAX_SITES = 8
# year windows actually needed (baseline 1986-2005; periods 2021-2040/2036-2055/2041-2070)
YEARS_HIST = (1985, 2006)
YEARS_FUT = (2020, 2071)


def gq_sites() -> pd.DataFrame:
    """Representative Guadalquivir sites: the 5 core + highest-N + highest-elevation."""
    s = pd.read_csv(C.PROCESSED / "sites.csv", dtype={"site_id": str})
    gq = s[s["in_guadalquivir_basin"]].copy()
    core = gq[gq["site_id"].isin(CORE)]
    rest = gq[~gq["site_id"].isin(CORE)].copy()
    rest["n_obs_qc"] = pd.to_numeric(rest["n_obs_qc"], errors="coerce").fillna(0)
    by_n = rest.sort_values("n_obs_qc", ascending=False).head(5)
    by_z = rest.sort_values("elevation_m", ascending=False).head(5)
    pick = pd.concat([core, by_n, by_z], ignore_index=True)
    pick = pick.drop_duplicates(subset=["site_id"]).dropna(subset=["lat", "lon"])
    pick = pick.head(MAX_SITES)
    pick["lat_r"] = pick["lat"].round(3)
    pick["lon_r"] = pick["lon"].round(3)
    return pick.drop_duplicates(subset=["lat_r", "lon_r"]).reset_index(drop=True)


def extract_one(gcm: str, exp: str, sites: pd.DataFrame) -> pd.DataFrame:
    cached = CACHE / f"{gcm}_{exp}.csv"
    if cached.exists():
        return pd.read_csv(cached, parse_dates=["date"])
    fname = FILE.format(gcm=gcm, exp=exp)
    url = f"{BASE}/{fname}"
    fs = fsspec.filesystem("http")
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
    df = pd.concat(rows, ignore_index=True)
    df = df.drop_duplicates(subset=["site_id", "date"])
    df.to_csv(cached, index=False)
    return df


def main() -> None:
    sites = gq_sites()
    print(f"{len(sites)} Guadalquivir sites: {list(sites['site_id'])}")
    frames = []
    t0 = time.time()
    for gcm in GCMS:
        for exp in ["historical"] + SCENARIOS:
            try:
                d = extract_one(gcm, exp, sites)
                frames.append(d)
                print(f"  {gcm:16s} {exp:11s} {len(d):6d} rows "
                      f"({time.time()-t0:.0f}s)", flush=True)
            except Exception as exc:  # noqa: BLE001
                print(f"  FAIL {gcm} {exp}: {exc!r}")
                C.record(f"PNACC {gcm} {exp}", url=f"{BASE}/{FILE.format(gcm=gcm, exp=exp)}",
                         status="FAILED", error=repr(exc), license=LICENSE,
                         version="PNACC 2024 CMIP6 5km ESD-RegBA")
    if frames:
        all_df = pd.concat(frames, ignore_index=True)
        all_df.to_parquet(OUT, index=False)
        print(f"wrote {OUT} ({len(all_df):,} rows)")
        C.record("PNACC tmean extracted at Guadalquivir sites", url=BASE, path=str(OUT),
                 license=LICENSE, version="PNACC 2024 CMIP6 5km ESD-RegBA",
                 notes=f"{len(GCMS)} GCMs x {len(SCENARIOS)+1} experiments; "
                       f"{len(sites)} sites; remote HDF5 subset (NCSS broken)")
    summary = {
        "n_gcm": len(GCMS), "scenarios": SCENARIOS, "sites": list(sites["site_id"]),
        "seconds": round(time.time() - t0, 1),
    }
    (C.LOGS / "pnacc_extract_summary.json").write_text(json.dumps(summary, indent=2),
                                                       encoding="utf-8")


if __name__ == "__main__":
    main()
