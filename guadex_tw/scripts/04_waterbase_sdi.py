"""Download the EEA SDI WISE6 disaggregated CSV (Deflate64 zip) and keep Spanish river Tw.

THREDDS-like large HTTP member ranges stall on this Nextcloud share, so the
1.64 GB zip is downloaded in full (checksummed) and then decompressed locally as
a stream; only Spanish river temperature rows are retained. The member uses
compression method 9 (Deflate64), handled by the `inflate64` package.
"""
from __future__ import annotations

import csv
import struct
import sys
import time
from pathlib import Path

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C

TOKEN = "ZcJ37WK6oaWonFJ"
ZIP_URL = ("https://sdi.eea.europa.eu/datashare/public.php/webdav/"
           "WISE6_DisaggregatedData-csv.zip")
LICENSE = "EEA open data"
VERSION = "WISE6 v2025 (Waterbase WQ ICM 2026, 1900-2025)"
AUTH = (TOKEN, "")
ZIP = C.RAW / "WISE6_DisaggregatedData-csv.zip"

OUT_RAW = C.INTERIM / "waterbase_obs_raw_spain.csv"
OUT_SPAIN = C.INTERIM / "waterbase_temp_spain.csv"
OUT_GQ = C.INTERIM / "waterbase_temp_guadalquivir.csv"
TEMP_CODE = "EEA_3121-01-5"
TEMP_LABEL = "Water temperature"


def download_zip() -> None:
    if ZIP.exists() and ZIP.stat().st_size > 1_600_000_000:
        print(f"zip already present ({ZIP.stat().st_size:,} bytes)")
        return
    print("Downloading WISE6 disaggregated zip ...")
    t0 = time.time()
    with requests.get(ZIP_URL, auth=AUTH, stream=True, timeout=(60, 300)) as r:
        r.raise_for_status()
        total = int(r.headers.get("content-length", 0))
        n = 0
        with open(ZIP, "wb") as fh:
            for chunk in r.iter_content(1 << 22):
                fh.write(chunk)
                n += len(chunk)
                if n // (256 << 20) != (n - len(chunk)) // (256 << 20):
                    print(f"  {n/1e6:.0f}/{total/1e6:.0f} MB "
                          f"({time.time()-t0:.0f}s)", flush=True)
    print(f"  done {n:,} bytes in {time.time()-t0:.0f}s")
    C.record("WISE6_DisaggregatedData-csv.zip", url=ZIP_URL, path=str(ZIP),
             license=LICENSE, version=VERSION,
             notes="EEA SDI Nextcloud public share; Deflate64 zip")


def member_info() -> tuple[int, int]:
    """Return (data_start, csize) for the single CSV member."""
    tail = requests.get(ZIP_URL, headers={"Range": "bytes=-300000"}, auth=AUTH, timeout=180).content
    eocd = tail.rfind(b"PK\x05\x06")
    cd_size = struct.unpack_from("<I", tail, eocd + 12)[0]
    cd_off = struct.unpack_from("<I", tail, eocd + 16)[0]
    cd = requests.get(ZIP_URL, headers={"Range": f"bytes={cd_off}-{cd_off + cd_size - 1}"},
                      auth=AUTH, timeout=180).content
    csize = struct.unpack_from("<I", cd, 20)[0]
    lhoff = struct.unpack_from("<I", cd, 42)[0]
    # read local header from the local copy
    with open(ZIP, "rb") as fh:
        fh.seek(lhoff)
        hdr = fh.read(300)
    fn_len, ex_len = struct.unpack_from("<2H", hdr, 26)
    return lhoff + 30 + fn_len + ex_len, csize


def stream_local(data_start: int, csize: int):
    from inflate64 import Inflater

    inf = Inflater()
    buf = b""
    t0 = time.time()
    with open(ZIP, "rb") as fh:
        fh.seek(data_start)
        remaining = csize
        read_total = 0
        while remaining > 0:
            comp = fh.read(min(1 << 22, remaining))
            if not comp:
                break
            remaining -= len(comp)
            read_total += len(comp)
            buf += inf.inflate(comp)
            parts = buf.split(b"\n")
            buf = parts.pop()  # keep the trailing partial line
            for line in parts:
                if line:
                    yield line
            if read_total // (256 << 20) != (read_total - len(comp)) // (256 << 20):
                print(f"  ... {read_total/1e6:.0f} MB compressed in "
                      f"({time.time()-t0:.0f}s)", flush=True)
        try:
            buf += inf.flush()
        except AttributeError:
            pass  # inflate64>=1.0 returns all output from inflate()
        for line in buf.split(b"\n"):
            if line:
                yield line


def main() -> None:
    t0 = time.time()
    download_zip()
    if not ZIP.exists():
        print("zip missing; aborting")
        return
    data_start, csize = member_info()
    print(f"member data_start={data_start} csize={csize:,}")

    n_lines = n_es = n_temp = 0
    with open(OUT_RAW, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["uid", "site_id", "obs_date", "obs_value", "unit", "below_loq",
                    "category", "label", "code"])
        first = True
        for raw in stream_local(data_start, csize):
            n_lines += 1
            if first:
                first = False
                continue
            if not raw.startswith(b"ES,"):
                continue
            parts = raw.decode("utf-8", "replace").split(",")
            if len(parts) < 28:
                continue
            n_es += 1
            if parts[4] != TEMP_CODE and parts[5] != TEMP_LABEL:
                continue
            if parts[3] != "RW":
                continue
            n_temp += 1
            w.writerow([parts[27], parts[1], parts[8], parts[10], parts[7],
                        parts[11], parts[3], parts[5], parts[4]])
    print(f"CSV lines: {n_lines:,}; ES rows: {n_es:,}; ES river Tw rows: {n_temp:,}")
    C.record("WISE6 disaggregated CSV (EEA SDI local stream)", url=ZIP_URL,
             path=str(OUT_RAW), license=LICENSE, version=VERSION,
             notes=f"inflate64; {n_temp} ES river Tw rows kept")

    if not (C.INTERIM / "waterbase_spatial_spain.csv").exists():
        print("spatial metadata missing; run 03_waterbase_direct.py first")
        return
    spat = pd.read_csv(C.INTERIM / "waterbase_spatial_spain.csv", dtype={"site_id": str})
    obs = pd.read_csv(OUT_RAW, dtype={"site_id": str, "uid": str})
    obs["obs_date"] = pd.to_datetime(obs["obs_date"], format="%Y%m%d", errors="coerce")
    obs["obs_value"] = pd.to_numeric(obs["obs_value"], errors="coerce")
    df = obs.merge(spat, on="site_id", how="left")
    df = df.dropna(subset=["obs_date", "obs_value", "lat", "lon"])
    df["obs_id"] = "WB_" + df["uid"].astype(str)
    df["unit"] = "Deg C"
    df["param_code"] = "TEMP"
    df["param_name"] = "Water Temperature"
    df["source_param_code"] = TEMP_CODE
    df["source"] = "WATERBASE_DIRECT"
    df["site_country"] = "Spain"
    df["drainage_region_name"] = df["rbd_name"]
    df["detection_limit_flag"] = df["below_loq"].map(
        {True: "below_loq", False: None, 1: "below_loq", 0: None,
         "1": "below_loq", "0": None, "True": "below_loq", "False": None})
    df["site_class"] = df["reservoir"].map(
        {True: "reservoir", False: "river", "True": "reservoir",
         "true": "reservoir", "False": "river", "false": "river"}).fillna("unknown")
    cols = ["obs_id", "site_id", "site_name", "site_country", "lat", "lon", "obs_date",
            "obs_value", "unit", "param_code", "param_name", "source_param_code", "source",
            "drainage_region_name", "wb_id", "wb_name", "rbd_id", "site_class",
            "detection_limit_flag", "category", "label", "code"]
    out = df[cols].copy()
    out.to_csv(OUT_SPAIN, index=False)
    gq = out[out["rbd_id"].astype(str).str.startswith("ES050")]
    gq.to_csv(OUT_GQ, index=False)
    print(f"Spain: {len(out):,} obs / {out['site_id'].nunique()} sites; "
          f"Guadalquivir ES050: {len(gq):,} obs / {gq['site_id'].nunique() if len(gq) else 0} sites")
    print(f"period {out['obs_date'].min().date()} .. {out['obs_date'].max().date()}")
    C.record("Waterbase direct (WISE6 SDI)", url=ZIP_URL, path=str(OUT_SPAIN),
             license=LICENSE, version=VERSION, status="SUCCESS",
             notes=f"{len(out)} ES river Tw obs, {out['site_id'].nunique()} sites; "
                   f"{len(gq)} in Guadalquivir ES050; {round(time.time()-t0)}s")


if __name__ == "__main__":
    main()
