"""Fetch daily air temperature and DEM elevation from Open-Meteo (no auth).

Every response is cached under data/raw/openmeteo/ and never re-requested.
Air temperature is fetched for the site's observation window plus a 31-day
lead-in so that backward-looking 30-day lags are fully populated.
"""
from __future__ import annotations

import json
import re
import sys
import time
from pathlib import Path

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C

CACHE = C.RAW / "openmeteo"
CACHE.mkdir(parents=True, exist_ok=True)
ARCHIVE = "https://archive-api.open-meteo.com/v1/archive"
ELEV = "https://api.open-meteo.com/v1/elevation"
THROTTLE = C.Throttle(2.5)


def request_json(url: str, params: dict, max_attempts: int = 8) -> tuple[dict, str]:
    """GET with exponential backoff on HTTP 429/5xx; honours Retry-After."""
    delay = 60.0
    last = None
    for attempt in range(1, max_attempts + 1):
        THROTTLE.wait()
        try:
            r = requests.get(url, params=params, timeout=120,
                             headers={"User-Agent": C.USER_AGENT})
        except requests.RequestException as exc:
            last = repr(exc)
            time.sleep(min(delay, 300))
            delay = min(delay * 1.6, 600)
            continue
        if r.status_code == 429 or r.status_code >= 500:
            ra = r.headers.get("Retry-After")
            wait = float(ra) if (ra and str(ra).strip().isdigit()) else delay
            print(f"    HTTP {r.status_code}; backing off {wait:.0f}s "
                  f"(attempt {attempt}/{max_attempts})", flush=True)
            time.sleep(wait)
            delay = min(delay * 1.6, 600)
            last = f"HTTP {r.status_code}"
            continue
        r.raise_for_status()
        return r.json(), r.url
    raise RuntimeError(f"gave up after {max_attempts} attempts: {last}")

SAFE = re.compile(r"[^0-9A-Za-z_.-]+")


def safe(s: str) -> str:
    return SAFE.sub("_", str(s))[:120]


def cached_get(url: str, logical_key: str, params: dict) -> dict:
    path = CACHE / f"{safe(logical_key)}.json"
    if path.exists():
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except Exception:  # noqa: BLE001
            path.unlink()
    j, req_url = request_json(url, params)
    j["_request_url"] = req_url
    j["_retrieved_utc"] = C.utcnow_iso()
    path.write_text(json.dumps(j), encoding="utf-8")
    return j


def site_list() -> pd.DataFrame:
    frames = []
    qc = C.INTERIM / "water_temp_qc.csv"
    if qc.exists():
        d = pd.read_csv(qc, parse_dates=["obs_date"])
        g = d.groupby("site_id").agg(
            lat=("lat", "first"), lon=("lon", "first"), site_name=("site_name", "first"),
            source=("source", "first"), first=("obs_date", "min"), last=("obs_date", "max"),
            n_obs=("obs_date", "size"),
        ).reset_index()
        frames.append(g)
    # also the raw Guadalquivir extract (includes single-obs Waterbase points)
    gq = C.INTERIM / "grqa_temp_guadalquivir.csv"
    if gq.exists():
        d = pd.read_csv(gq, parse_dates=["obs_date"]).rename(columns={"lat_wgs84": "lat", "lon_wgs84": "lon"})
        g = d.groupby("site_id").agg(
            lat=("lat", "first"), lon=("lon", "first"), site_name=("site_name", "first"),
            source=("source", "first"), first=("obs_date", "min"), last=("obs_date", "max"),
            n_obs=("obs_date", "size"),
        ).reset_index()
        frames.append(g)
    wb = C.INTERIM / "waterbase_temp_spain.csv"
    if wb.exists():
        d = pd.read_csv(wb, parse_dates=["obs_date"])
        g = d.groupby("site_id").agg(
            lat=("lat", "first"), lon=("lon", "first"), site_name=("site_name", "first"),
            source=("source", "first"), first=("obs_date", "min"), last=("obs_date", "max"),
            n_obs=("obs_date", "size"),
        ).reset_index()
        frames.append(g)
    allsites = pd.concat(frames, ignore_index=True).drop_duplicates(subset=["site_id"])
    return allsites


def main() -> None:
    sites = site_list()
    print(f"{len(sites)} sites to fetch")
    elev_rows = []
    n_new, n_cached, n_fail = 0, 0, 0
    for i, row in sites.iterrows():
        sid = str(row["site_id"])
        lat, lon = float(row["lat"]), float(row["lon"])
        start = (pd.Timestamp(row["first"]) - pd.Timedelta(days=31)).strftime("%Y-%m-%d")
        end = pd.Timestamp(row["last"]).strftime("%Y-%m-%d")
        klat, klon = round(lat, 4), round(lon, 4)
        key = f"ta_{sid}_{klat}_{klon}"
        existed = (CACHE / f"{safe(key)}.json").exists()
        try:
            cached_get(ARCHIVE, key, {
                "latitude": lat, "longitude": lon,
                "start_date": start, "end_date": end,
                "daily": "temperature_2m_mean,temperature_2m_max,temperature_2m_min",
                "timezone": "UTC",
            })
            ekey = f"elev_{sid}_{klat}_{klon}"
            ej = cached_get(ELEV, ekey, {"latitude": lat, "longitude": lon})
            elev_rows.append({
                "site_id": sid, "lat": lat, "lon": lon,
                "openmeteo_elevation_m": (ej.get("elevation") or [None])[0],
            })
            if existed:
                n_cached += 1
            else:
                n_new += 1
        except Exception as exc:  # noqa: BLE001
            n_fail += 1
            print(f"  FAIL {sid}: {exc!r}")
            C.record(f"Open-Meteo {sid}", url=ARCHIVE, status="FAILED", error=repr(exc),
                     notes=f"lat={lat} lon={lon} {start}..{end}")
        if (i + 1) % 25 == 0:
            print(f"  {i+1}/{len(sites)} ({n_new} new, {n_cached} cached, {n_fail} fail)",
                  flush=True)
    pd.DataFrame(elev_rows).to_csv(C.INTERIM / "openmeteo_elevation.csv", index=False)
    print(f"done: {n_new} new, {n_cached} cached, {n_fail} failed")
    C.record("Open-Meteo ERA5 archive + DEM elevation (bulk)", url=ARCHIVE,
             path=str(CACHE), status="SUCCESS" if n_fail == 0 else "PARTIAL",
             license="CC-BY-4.0 (Open-Meteo, ERA5/Copernicus DEM)",
             notes=f"{len(sites)} sites, {n_new} new, {n_cached} cached, {n_fail} failed")


if __name__ == "__main__":
    main()
