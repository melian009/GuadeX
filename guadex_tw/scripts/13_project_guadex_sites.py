"""WP1: project daily water temperature at every GuadeX site (all GCMs/SSPs).

The calibration pipeline (`09_extract_pnacc.py` / `10_project.py`) extracts PNACC
grid cells only for a stratified sample of monitoring stations (`MAX_SITES = 10`)
because that is what the Tw model is validated against.  The climate-scenario ODE
runs on the GuadeX network, so it needs forcing at every site coordinate.

This script reads the GuadeX site table (`data/ConnectivityUTM.csv`, ETRS89/UTM
zone 30N), extracts the PNACC daily `tmean` at each site, applies the calibrated
spatial Tw model and writes the per-site daily ensemble-median series.

Outputs (in `guadex_tw/outputs/tables/`):

* `water_temp_daily_guadex_sites_wide.csv` - the product consumed by
  `Guadex.load_wide_daily_forcing`.  Columns: `date, scenario, <site codes...>`;
  rows for `historical` (1986-2005 baseline) and each SSP (2026-2045).
* `water_temp_baseline_guadex_sites.csv` - per-site baseline mean Tw
  (`site_id, tw_baseline_mean`) for the 1986-2005 window.
* `water_temp_daily_guadex_sites.parquet` - optional long format
  (`site_id, scenario, date, tw_ensemble_median`) with `--long`.

Efficiency
----------
PNACC files are 2 GB NetCDF-4 opened remotely over HTTP.  Reading one grid cell
at a time costs ~2 s per read (~20 min per file at 776 sites); reading the
Guadalquivir bounding box once per file is ~50x more bandwidth-efficient, so the
extraction reads the bbox `[t0:t1, lat0:lat1, lon0:lon1]` and then slices the
nearest cell for every site.  Per-(GCM, experiment) arrays are cached as
`.npy` (float32) and the run is resumable.

Usage:
    python scripts/13_project_guadex_sites.py                 # all GCMs
    python scripts/13_project_guadex_sites.py --gcms ACCESS-CM2 --max-seconds 600
    python scripts/13_project_guadex_sites.py --max-sites 30   # pilot
    python scripts/13_project_guadex_sites.py --combine-only   # cache -> outputs
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
FUTURE = (2026, 2045)   # ODE horizon (run_climate_scenarios end_year)

CACHE = C.INTERIM / "pnacc_cache_guadex"
CACHE.mkdir(parents=True, exist_ok=True)
WIDE = C.TABLES / "water_temp_daily_guadex_sites_wide.csv"
LONG = C.TABLES / "water_temp_daily_guadex_sites.parquet"
BASELINE_CSV = C.TABLES / "water_temp_baseline_guadex_sites.csv"
CELL_CACHE = CACHE / "cell_index.csv"
LICENSE = "AEMET/PNACC 2024 - public (attribution)"


def guadex_sites(path: Path = GUADEX_SITES_CSV) -> pd.DataFrame:
    """GuadeX site table with lat/lon/elevation (ETRS89/UTM 30N -> WGS84)."""
    if not path.exists():
        raise SystemExit(f"GuadeX site table not found: {path}")
    df = pd.read_csv(path, dtype={"CODIGO": str}, low_memory=False)
    missing = {"CODIGO", "ALTITUD", "UTMX", "UTMY"} - set(df.columns)
    if missing:
        raise SystemExit(f"{path} is missing columns: {sorted(missing)}")
    df = df.dropna(subset=["UTMX", "UTMY", "ALTITUD"]).copy()
    try:
        from pyproj import Transformer
        t = Transformer.from_crs("EPSG:25830", "EPSG:4326", always_xy=True)
        lon, lat = t.transform(df["UTMX"].to_numpy(float), df["UTMY"].to_numpy(float))
        df["lat"], df["lon"] = lat, lon
    except Exception:  # pragma: no cover - fallback to the embedded transform
        coords = [_utm_to_latlon(float(x), float(y)) for x, y in zip(df["UTMX"], df["UTMY"])]
        df["lat"] = [c[0] for c in coords]
        df["lon"] = [c[1] for c in coords]
    df["elevation_m"] = pd.to_numeric(df["ALTITUD"], errors="coerce")
    df = df.dropna(subset=["elevation_m", "lat", "lon"]).reset_index(drop=True)
    return df


def _utm_to_latlon(easting: float, northing: float, zone: int = 30) -> tuple[float, float]:
    a, f, k0 = 6378137.0, 1 / 298.257223563, 0.9996
    e2 = f * (2 - f)
    e1 = (1 - math.sqrt(1 - e2)) / (1 + math.sqrt(1 - e2))
    x, y = easting - 500000.0, northing
    lon0 = math.radians(zone * 6 - 183)
    m = y / k0
    mu = m / (a * (1 - e2 / 4 - 3 * e2**2 / 64 - 5 * e2**3 / 256))
    phi1 = (mu + (3 * e1 / 2 - 27 * e1**3 / 32) * math.sin(2 * mu)
            + (21 * e1**2 / 16 - 55 * e1**4 / 32) * math.sin(4 * mu)
            + (151 * e1**3 / 96) * math.sin(6 * mu))
    n1 = a / math.sqrt(1 - e2 * math.sin(phi1)**2)
    t1 = math.tan(phi1)**2
    c1 = e2 / (1 - e2) * math.cos(phi1)**2
    r1 = a * (1 - e2) / (1 - e2 * math.sin(phi1)**2)**1.5
    d = x / (n1 * k0)
    lat = phi1 - (n1 * math.tan(phi1) / r1) * (
        d**2 / 2 - (5 + 3 * t1 + 10 * c1 - 4 * c1**2 - 9 * e2 / (1 - e2)) * d**4 / 24
        + (61 + 90 * t1 + 298 * c1 + 45 * t1**2 - 252 * e2 / (1 - e2) - 3 * c1**2) * d**6 / 720)
    lon = lon0 + (d - (1 + 2 * t1 + c1) * d**3 / 6
                  + (5 - 2 * c1 + 28 * t1 - 3 * c1**2 + 8 * e2 / (1 - e2) + 24 * t1**2) * d**5 / 120) / math.cos(phi1)
    return math.degrees(lat), math.degrees(lon)


def _open_remote(gcm: str, exp: str):
    import aiohttp
    import fsspec
    import h5py
    fs = fsspec.filesystem("http",
        client_kwargs={"timeout": aiohttp.ClientTimeout(total=1800, sock_read=180)})
    f = fs.open(f"{BASE}/{FILE.format(gcm=gcm, exp=exp)}", "rb")
    return f, h5py.File(f, "r")


def _time_dates(h) -> pd.DatetimeIndex:
    tunits = h["time"].attrs["units"]
    tunits = tunits.decode() if isinstance(tunits, bytes) else str(tunits)
    tref = pd.Timestamp(tunits.split("since")[1].strip())
    return pd.to_datetime(tref) + pd.to_timedelta(h["time"][:], unit="D")


def cell_index(sites: pd.DataFrame) -> pd.DataFrame:
    """Nearest PNACC grid cell (0-based lat/lon index) for every site."""
    if CELL_CACHE.exists():
        cached = pd.read_csv(CELL_CACHE, dtype={"site_id": str})
        if set(cached.site_id) == set(sites.CODIGO.astype(str)):
            return cached
    f, h = _open_remote(GCMS[0], SCENARIOS[0])
    lat, lon = h["lat"][:], h["lon"][:]
    h.close()
    f.close()
    i = np.abs(lat[None, :] - sites["lat"].to_numpy()[:, None]).argmin(axis=1)
    j = np.abs(lon[None, :] - sites["lon"].to_numpy()[:, None]).argmin(axis=1)
    out = pd.DataFrame({"site_id": sites.CODIGO.astype(str), "i": i, "j": j})
    out.to_csv(CELL_CACHE, index=False)
    return out


def dates_for(exp: str) -> pd.DatetimeIndex:
    """Cache the date vector for an experiment from the first GCM file."""
    path = CACHE / f"dates_{exp}.npy"
    if path.exists():
        return pd.DatetimeIndex(np.load(path, allow_pickle=True))
    f, h = _open_remote(GCMS[0], exp)
    dates = _time_dates(h)
    h.close()
    f.close()
    years = BASELINE if exp == "historical" else FUTURE
    sel = (dates.year >= years[0]) & (dates.year <= years[1])
    np.save(path, dates[sel].to_numpy())
    return dates[sel]


def extract_one(gcm: str, exp: str, sites: pd.DataFrame, cells: pd.DataFrame,
                dates: pd.DatetimeIndex) -> np.ndarray:
    """Nearest-cell daily tmean (n_sites x n_days) for one (GCM, experiment)."""
    cached = CACHE / f"{gcm}_{exp}_{len(sites)}_ta.npy"
    if cached.exists():
        return np.load(cached)
    f, h = _open_remote(gcm, exp)
    full_dates = _time_dates(h)
    years = BASELINE if exp == "historical" else FUTURE
    tidx = np.where((full_dates.year >= years[0]) & (full_dates.year <= years[1]))[0]
    i = cells["i"].to_numpy()
    j = cells["j"].to_numpy()
    i0, i1 = int(i.min()), int(i.max())
    j0, j1 = int(j.min()), int(j.max())
    # One bounding-box read over the whole horizon is far cheaper over HTTP
    # than one read per site.
    slab = np.asarray(h["tmean"][int(tidx[0]):int(tidx[-1]) + 1, i0:i1 + 1, j0:j1 + 1],
                      dtype=np.float32)
    h.close()
    f.close()
    n = min(slab.shape[0], len(dates))
    out = np.empty((len(sites), n), dtype=np.float32)
    for k in range(len(sites)):
        out[k, :] = slab[:n, i[k] - i0, j[k] - j0]
    np.save(cached, out)
    return out


def _fill_gaps(frame: pd.DataFrame) -> pd.DataFrame:
    """Linearly interpolate short NaN gaps along time (some GCM files miss a few
    days).  Edge NaNs are back/forward filled so every site has a daily series."""
    cols = [c for c in frame.columns if c not in ("date", "scenario")]
    frame[cols] = frame[cols].astype("float64").interpolate(axis=0, limit_direction="both")
    return frame


def project_tw(mod, model, ta: np.ndarray, z: np.ndarray, months: np.ndarray) -> np.ndarray:
    tw = np.empty_like(ta, dtype=np.float64)
    for k in range(ta.shape[0]):
        tw[k, :] = mod.predict_tw(model, ta[k, :].astype(float), float(z[k]), months)
    return tw


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gcms", default="", help="comma-separated GCM subset")
    ap.add_argument("--max-sites", type=int, default=0, help="pilot: first N sites")
    ap.add_argument("--max-seconds", type=float, default=1e12)
    ap.add_argument("--combine-only", action="store_true")
    ap.add_argument("--long", action="store_true", help="also write the long parquet")
    ap.add_argument("--per-gcm", action="store_true",
                    help="also write one wide file per (scenario, GCM) so the ODE "
                         "ensemble carries the GCM spread instead of the ensemble median")
    args = ap.parse_args()

    sites = guadex_sites()
    if args.max_sites:
        sites = sites.head(args.max_sites).reset_index(drop=True)
    gcms = [g.strip() for g in args.gcms.split(",") if g.strip()] or GCMS
    print(f"{len(sites)} sites; {len(gcms)} GCMs; elevation "
          f"{sites.elevation_m.min():.0f}-{sites.elevation_m.max():.0f} m", flush=True)

    t0 = time.time()
    if not args.combine_only:
        cells = cell_index(sites)
        for exp in ["historical"] + SCENARIOS:
            dates = dates_for(exp)
            for gcm in gcms:
                if (CACHE / f"{gcm}_{exp}_{len(sites)}_ta.npy").exists():
                    continue
                if time.time() - t0 > args.max_seconds:
                    print("wall-clock budget reached; cache is partial "
                          "(rerun to resume)", flush=True)
                    return
                try:
                    arr = extract_one(gcm, exp, sites, cells, dates)
                    print(f"  {gcm:16s} {exp:11s} {arr.shape} ({time.time()-t0:.0f}s)",
                          flush=True)
                except Exception as exc:  # noqa: BLE001
                    print(f"  FAIL {gcm} {exp}: {type(exc).__name__}: {exc}", flush=True)

    missing = [(g, e) for g in gcms for e in ["historical"] + SCENARIOS
               if not (CACHE / f"{g}_{e}_{len(sites)}_ta.npy").exists()]
    if missing:
        print(f"refusing to write outputs: {len(missing)} cache file(s) missing, "
              f"e.g. {missing[:3]}", flush=True)
        return

    # --- project Tw and ensemble-median per experiment ---
    mod, model = _load_model()
    z = (sites.elevation_m.to_numpy(float) / 1000.0)
    site_cols = sites.CODIGO.astype(str).tolist()

    per_exp: dict[str, pd.DataFrame] = {}
    baseline_rows = []
    for exp in ["historical"] + SCENARIOS:
        dates = dates_for(exp)
        months = dates.month.to_numpy()
        stack = []
        for gcm in gcms:
            ta = np.load(CACHE / f"{gcm}_{exp}_{len(sites)}_ta.npy")
            stack.append(project_tw(mod, model, ta, z, months))
        tw_med = np.nanmedian(np.stack(stack), axis=0)  # (n_sites, n_days)
        frame = pd.DataFrame(tw_med.T, columns=site_cols)
        frame.insert(0, "date", dates.to_numpy())
        frame.insert(1, "scenario", exp)
        frame = _fill_gaps(frame)
        per_exp[exp] = frame
        if exp == "historical":
            base_mean = np.nanmean(tw_med, axis=1)
            baseline_rows.append(pd.DataFrame({"site_id": site_cols,
                                               "tw_baseline_mean": base_mean}))
        print(f"  projected {exp}: {frame.shape}", flush=True)

    wide = pd.concat(per_exp.values(), ignore_index=True)
    wide = wide.sort_values(["scenario", "date"])
    wide.to_csv(WIDE, index=False)
    print(f"wrote {WIDE} ({len(wide):,} rows x {wide.shape[1]} cols)")

    baseline = pd.concat(baseline_rows, ignore_index=True)
    baseline.to_csv(BASELINE_CSV, index=False)
    print(f"wrote {BASELINE_CSV} ({len(baseline)} sites)")

    if args.per_gcm:
        hist_dates = dates_for("historical")
        hist_months = hist_dates.month.to_numpy()
        for scen in SCENARIOS:
            f_dates = dates_for(scen)
            f_months = f_dates.month.to_numpy()
            for gcm in gcms:
                ta_h = np.load(CACHE / f"{gcm}_historical_{len(sites)}_ta.npy")
                ta_f = np.load(CACHE / f"{gcm}_{scen}_{len(sites)}_ta.npy")
                tw_h = project_tw(mod, model, ta_h, z, hist_months)
                tw_f = project_tw(mod, model, ta_f, z, f_months)
                frame_h = pd.DataFrame(tw_h.T, columns=site_cols)
                frame_h.insert(0, "date", hist_dates.to_numpy())
                frame_h.insert(1, "scenario", "historical")
                frame_f = pd.DataFrame(tw_f.T, columns=site_cols)
                frame_f.insert(0, "date", f_dates.to_numpy())
                frame_f.insert(1, "scenario", scen)
                out = _fill_gaps(pd.concat([frame_h, frame_f], ignore_index=True))
                path = C.TABLES / f"water_temp_daily_guadex_sites_wide_{scen}_{gcm}.csv"
                out.to_csv(path, index=False)
            print(f"  per-GCM wide files written for {scen} ({len(gcms)} GCMs)", flush=True)

    if args.long:
        long = wide.melt(id_vars=["date", "scenario"], var_name="site_id",
                         value_name="tw_ensemble_median")
        long.to_parquet(LONG, index=False)
        print(f"wrote {LONG} ({len(long):,} rows)")

    C.record("Projected daily Tw at all GuadeX sites (WP1)", path=str(WIDE),
             license=LICENSE, version="PNACC 2024 CMIP6 5km ESD-RegBA + Model 2",
             notes=f"{len(gcms)} GCMs x {len(SCENARIOS)} SSPs; {len(sites)} sites; "
                   f"ensemble median; baseline {BASELINE}; horizon {FUTURE}")
    (C.LOGS / "projection_summary_guadex_sites.json").write_text(json.dumps({
        "n_sites": len(sites), "n_gcm": len(gcms), "scenarios": SCENARIOS,
        "baseline": BASELINE, "future": FUTURE,
    }, indent=2), encoding="utf-8")


def _load_model():
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "project_model", Path(__file__).resolve().parent / "10_project.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    model, _ = mod.fit_model()
    return mod, model


if __name__ == "__main__":
    main()
