"""Download GRQA v1.3 water-temperature table and metadata.

Uses HTTP range requests to fetch only TEMP_GRQA.csv out of the 1.147 GB zip,
decompressing as a stream so the ~1 GB table is never held fully in memory.
Falls back to a full download + local extract only if range requests fail.
"""
from __future__ import annotations

import io
import struct
import sys
import zipfile
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C

URL = "https://zenodo.org/api/records/7056647/files/GRQA_data_v1.3.zip/content"
META_URL = "https://zenodo.org/api/records/7056647/files/GRQA_meta.zip/content"
RECORD_API = "https://zenodo.org/api/records/7056647"
DOI = "10.5281/zenodo.7056647"
LICENSE = "CC-BY-4.0"
VERSION = "1.3"
UA = {"User-Agent": C.USER_AGENT}
TIMEOUT = (60, 600)

OUT_CSV = C.RAW / "TEMP_GRQA.csv"
META_ZIP = C.RAW / "GRQA_meta.zip"


def get(url: str, headers: dict | None = None, **kw) -> requests.Response:
    h = dict(UA)
    if headers:
        h.update(headers)
    r = requests.get(url, headers=h, timeout=TIMEOUT, **kw)
    r.raise_for_status()
    return r


def stream_decompress(r: requests.Response, csize: int, out: Path, method: int) -> int:
    """Stream a zip member to disk. method 8=deflate, 0=stored."""
    dec = None
    if method == 8:
        import zlib

        dec = zlib.decompressobj(-15)
    written = 0
    with open(out, "wb") as fh:
        for chunk in r.iter_content(chunk_size=1 << 23):
            if not chunk:
                continue
            blob = dec.decompress(chunk) if dec is not None else chunk
            fh.write(blob)
            written += len(blob)
        if dec is not None:
            fh.write(dec.flush())
    return written


def range_extract_csv() -> bool:
    print("Probing zip central directory via HTTP range ...")
    tail = get(URL, {"Range": "bytes=-300000"}).content
    head = get(URL, {"Range": "bytes=0-1"})
    total = int(head.headers["content-range"].split()[1].split("-")[0])
    base = total - len(tail)
    print(f"  archive size: {total} bytes")

    eocd = tail.rfind(b"PK\x05\x06")
    if eocd < 0:
        print("  no EOCD found")
        return False
    cd_size = struct.unpack_from("<I", tail, eocd + 12)[0]
    cd_off = struct.unpack_from("<I", tail, eocd + 16)[0]

    cd = get(URL, {"Range": f"bytes={cd_off}-{cd_off + cd_size - 1}"}).content
    entries = {}
    p = 0
    while p < len(cd) and cd[p:p + 4] == b"PK\x01\x02":
        method = struct.unpack_from("<H", cd, p + 10)[0]
        csize, usize = struct.unpack_from("<2I", cd, p + 20)
        name_len, extra_len, comment_len = struct.unpack_from("<3H", cd, p + 28)
        lhoff = struct.unpack_from("<I", cd, p + 42)[0]
        name = cd[p + 46:p + 46 + name_len].decode("utf-8", "replace")
        entries[name] = dict(csize=csize, usize=usize, lhoff=lhoff, method=method)
        p += 46 + name_len + extra_len + comment_len
    print(f"  {len(entries)} central-directory entries")

    target = [k for k in entries if k.endswith("TEMP_GRQA.csv")]
    if not target:
        print("  TEMP_GRQA.csv not present in archive")
        return False
    name = target[0]
    e = entries[name]
    print(f"  member: {name} csize={e['csize']} usize={e['usize']} method={e['method']}")

    hdr = get(URL, {"Range": f"bytes={e['lhoff']}-{e['lhoff'] + 300}"}).content
    fn_len, ex_len = struct.unpack_from("<2H", hdr, 26)
    data_start = e["lhoff"] + 30 + fn_len + ex_len

    r = requests.get(
        URL,
        headers={**UA, "Range": f"bytes={data_start}-{data_start + e['csize'] - 1}"},
        timeout=TIMEOUT,
        stream=True,
    )
    r.raise_for_status()
    written = stream_decompress(r, e["csize"], OUT_CSV, e["method"])
    print(f"  wrote {written} bytes (expected {e['usize']})")
    if written != e["usize"]:
        print("  WARNING: uncompressed size mismatch")
    return True


def full_download_fallback() -> bool:
    print("Falling back to full archive download ...")
    zpath = C.RAW / "GRQA_data_v1.3.zip"
    with requests.get(URL, headers=UA, timeout=TIMEOUT, stream=True) as r:
        r.raise_for_status()
        with open(zpath, "wb") as fh:
            for chunk in r.iter_content(1 << 23):
                fh.write(chunk)
    with zipfile.ZipFile(zpath) as zf:
        member = [n for n in zf.namelist() if n.endswith("TEMP_GRQA.csv")][0]
        with zf.open(member) as src, open(OUT_CSV, "wb") as dst:
            while True:
                b = src.read(1 << 23)
                if not b:
                    break
                dst.write(b)
    return True


def download_meta() -> None:
    print("Downloading GRQA_meta.zip ...")
    with requests.get(META_URL, headers=UA, timeout=TIMEOUT, stream=True) as r:
        r.raise_for_status()
        with open(META_ZIP, "wb") as fh:
            for chunk in r.iter_content(1 << 20):
                fh.write(chunk)
    with zipfile.ZipFile(META_ZIP) as zf:
        for m in zf.namelist():
            if m.lower().endswith((".csv", ".txt", ".md")):
                dest = C.RAW / Path(m).name
                with zf.open(m) as src, open(dest, "wb") as dst:
                    dst.write(src.read())
                print(f"  extracted {dest.name}")
    C.record(
        "GRQA_meta.zip", url=META_URL, path=str(META_ZIP), license=LICENSE,
        doi=DOI, version=VERSION, notes="cross-source duplicate pairs + param stats",
    )


def main() -> None:
    if OUT_CSV.exists() and OUT_CSV.stat().st_size > 1_000_000_000:
        print(f"Already have {OUT_CSV} ({OUT_CSV.stat().st_size} bytes)")
    else:
        try:
            ok = range_extract_csv()
        except Exception as exc:  # noqa: BLE001
            print(f"Range extraction failed: {exc!r}")
            C.record("GRQA TEMP_GRQA.csv", url=URL, status="FAILED",
                     error=f"range extraction: {exc!r}", license=LICENSE, doi=DOI,
                     version=VERSION, notes="will retry full download")
            ok = full_download_fallback()
        if ok:
            C.record(
                "GRQA TEMP_GRQA.csv", url=URL, path=str(OUT_CSV), license=LICENSE,
                doi=DOI, version=VERSION,
                notes="extracted from GRQA_data_v1.3.zip via HTTP range + streaming zlib",
            )
    download_meta()
    print("Done.")


if __name__ == "__main__":
    main()
