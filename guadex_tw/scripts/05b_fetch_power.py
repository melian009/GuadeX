"""Fetch daily air temperature from NASA POWER (no auth) for every site.

NASA POWER (MERRA-2, 0.5 x 0.625 deg) returns T2M / T2M_MAX / T2M_MIN and the
grid-cell elevation, which lets us apply the lapse-rate correction
Ta_corr = Ta_gridcell + GAMMA * (z_site - z_gridcell).

Open-Meteo ERA5 was the intended workhorse but its weighted rate limit blocked
bulk retrieval (HTTP 429); NASA POWER is used as the consistent calibration
predictor and Open-Meteo cache is retained as a partial cross-check.
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

CACHE = C.RAW / "nasa_power"
CACHE.mkdir(parents=True, exist_ok=True)
URL = "https://power.larc.nasa.gov/api/temporal/daily/point"
THROTTLE = C.Throttle(1.0)
POWER_MIN_DATE = pd.Timestamp("1981-01-01")
SAFE = re.compile(r"[^0-9A-Za-z_.-]+")


def safe(s: str) -> str:
    return SAFE.sub("_", str(s))[:120]


def site_list() -> pd.DataFrame:
    frames = []
    for p in (C.PROCESSED / "water_temp_all_qc.csv", C.INTERIM / "grqa_temp_guadalquivir.csv"):
        if p.exists():
            d = pd.read_csv(p, parse_dates=["obs_date"], dtype={"site_id": str})
            g = d.groupby("site_id", as_index=False).agg(
                site_name=("site_name", "first"), lat=("lat", "first"),
                lon=("lon", "first"), first=("obs_date", "min"), last=("obs_date", "max"),
            )
            frames.append(g)
    s = pd.concat(frames, ignore_index=True)
    s = s.groupby("site_id", as_index=False).agg(
        site_name=("site_name", "first"), lat=("lat", "first"), lon=("lon", "first"),
        first=("first", "min"), last=("last", "max"))
    return s


def fetch_power(lat: float, lon: float, start: str, end: str, retries: int = 4) -> dict:
    params = {
        "parameters": "T2M,T2M_MAX,T2M_MIN", "community": "AG",
        "longitude": lon, "latitude": lat, "start": start, "end": end,
        "format": "JSON", "time-standard": "UTC",
    }
    delay = 10.0
    last = None
    for attempt in range(1, retries + 1):
        THROTTLE.wait()
        try:
            r = requests.get(URL, params=params, timeout=120,
                             headers={"User-Agent": C.USER_AGENT})
            if r.status_code in (429, 503) or r.status_code >= 500:
                time.sleep(delay)
                delay = min(delay * 1.8, 300)
                last = f"HTTP {r.status_code}"
                continue
            r.raise_for_status()
            j = r.json()
            if "properties" not in j or "parameter" not in j.get("properties", {}):
                msg = j.get("messages") or j.get("header") or "unknown POWER error"
                raise ValueError(f"POWER error response: {msg}")
            return j
        except ValueError:
            raise
        except Exception as exc:  # noqa: BLE001
            last = repr(exc)
            time.sleep(delay)
            delay = min(delay * 2, 60)
    raise RuntimeError(f"NASA POWER failed after {retries}: {last}")


def main() -> None:
    sites = site_list()
    print(f"{len(sites)} sites", flush=True)
    rows = []
    n_new = n_cached = n_fail = n_skip = 0
    for i, r in sites.iterrows():
        sid = str(r["site_id"])
        path = CACHE / f"{safe(sid)}.json"
        first = pd.Timestamp(r["first"])
        last = pd.Timestamp(r["last"])
        start = max(first - pd.Timedelta(days=31), POWER_MIN_DATE).strftime("%Y%m%d")
        end = last.strftime("%Y%m%d")
        if last < POWER_MIN_DATE:
            n_skip += 1
            C.record(f"NASA POWER {sid}", url=URL, status="FAILED",
                     error=f"observation period {first.date()}..{last.date()} precedes "
                           f"NASA POWER start 1981-01-01",
                     notes="no POWER data available for this site")
            continue
        try:
            if path.exists():
                j = json.loads(path.read_text(encoding="utf-8"))
                if "properties" not in j or "parameter" not in j.get("properties", {}):
                    path.unlink()
                    raise ValueError(f"cached response invalid: {j.get('messages')}")
                n_cached += 1
            else:
                j = fetch_power(float(r["lat"]), float(r["lon"]), start, end)
                j["_retrieved_utc"] = C.utcnow_iso()
                path.write_text(json.dumps(j), encoding="utf-8")
                n_new += 1
            elev = (j.get("geometry", {}).get("coordinates") or [None, None, None])[2]
            rows.append({"site_id": sid, "nasa_power_elevation_m": elev,
                         "power_start": start, "power_end": end})
        except Exception as exc:  # noqa: BLE001
            n_fail += 1
            print(f"  FAIL {sid}: {exc!r}", flush=True)
            C.record(f"NASA POWER {sid}", url=URL, status="FAILED", error=repr(exc),
                     notes=f"lat={r['lat']} lon={r['lon']} {start}..{end}")
        if (i + 1) % 20 == 0:
            print(f"  {i+1}/{len(sites)} ({n_new} new, {n_cached} cached, "
                  f"{n_fail} fail, {n_skip} skipped)", flush=True)
    pd.DataFrame(rows).to_csv(C.INTERIM / "nasa_power_elevation.csv", index=False)
    print(f"done: {n_new} new, {n_cached} cached, {n_fail} failed, {n_skip} skipped")
    C.record("NASA POWER daily T2M/T2M_MAX/T2M_MIN (bulk)", url=URL, path=str(CACHE),
             status="SUCCESS" if n_fail == 0 else "PARTIAL",
             license="NASA open data (MERRA-2)",
             notes=f"{len(sites)} sites, {n_new} new, {n_cached} cached, {n_fail} failed; "
                   f"0.5x0.625 deg, gridcell elevation returned for lapse correction")


if __name__ == "__main__":
    main()
