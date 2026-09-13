"""Shared paths, provenance logging and small IO helpers for the GuadeX Tw pipeline.

Every download or generated artefact is logged to data/raw/provenance.jsonl.
PROVENANCE.md is rendered from that JSONL so the log cannot drift from the files.
"""
from __future__ import annotations

import hashlib
import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RAW = DATA / "raw"
INTERIM = DATA / "interim"
PROCESSED = DATA / "processed"
MODELS = ROOT / "models"
TABLES = ROOT / "outputs" / "tables"
FIGURES = ROOT / "outputs" / "figures"
LOGS = ROOT / "logs"
PROV_JSONL = RAW / "provenance.jsonl"

for _d in (RAW, INTERIM, PROCESSED, MODELS, TABLES, FIGURES, LOGS):
    _d.mkdir(parents=True, exist_ok=True)


def utcnow_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


USER_AGENT = "GuadeX-Tw-pipeline/1.0 (hydrology research; contact: GuadeX project)"


def checksums(path: Path) -> dict:
    """Return md5, sha256 and size for a file (streamed, safe for large files)."""
    md5 = hashlib.md5()
    sha = hashlib.sha256()
    n = 0
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            md5.update(chunk)
            sha.update(chunk)
            n += len(chunk)
    return {"md5": md5.hexdigest(), "sha256": sha.hexdigest(), "bytes": n}


def record(
    source: str,
    *,
    url: str = "",
    path: str = "",
    license: str = "",
    doi: str = "",
    version: str = "",
    status: str = "SUCCESS",
    http_status: object = "",
    error: str = "",
    notes: str = "",
    checksum: bool = True,
) -> None:
    """Append one provenance row and (re)render PROVENANCE.md."""
    row = {
        "timestamp_utc": utcnow_iso(),
        "source": source,
        "url": url,
        "path": str(path).replace("\\", "/"),
        "license": license,
        "doi": doi,
        "version": version,
        "status": status,
        "http_status": http_status,
        "error": error,
        "notes": notes,
    }
    if checksum and path and Path(path).exists() and Path(path).is_file():
        try:
            row.update(checksums(Path(path)))
        except Exception as exc:  # pragma: no cover
            row["checksum_error"] = str(exc)
    with open(PROV_JSONL, "a", encoding="utf-8") as fh:
        fh.write(json.dumps(row, ensure_ascii=False) + "\n")
    render_provenance_md()


def render_provenance_md() -> None:
    rows = []
    if PROV_JSONL.exists():
        for line in PROV_JSONL.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    # last write wins for a given (source, path)
    dedup: dict = {}
    for r in rows:
        dedup[(r.get("source", ""), r.get("path", "") or r.get("url", ""))] = r
    rows = list(dedup.values())

    cols = [
        "timestamp_utc", "source", "version", "status", "http_status",
        "bytes", "md5", "sha256", "license", "doi", "url", "path", "notes", "error",
    ]
    lines = [
        "# PROVENANCE",
        "",
        "Every remote file or generated artefact that feeds the GuadeX Tw pipeline.",
        "Rows are appended automatically by `scripts/common.record()` and de-duplicated by",
        "(source, path/url) with last-write-wins. Checksums are streamed (MD5 + SHA256).",
        "",
        f"For download failures the `status` column is `FAILED` and `error` carries the exact",
        f"message; `http_status` the HTTP code. Generated at {utcnow_iso()}.",
        "",
    ]
    header = "| " + " | ".join(cols) + " |"
    sep = "|" + "|".join(["---"] * len(cols)) + "|"
    lines += [header, sep]
    for r in rows:
        cells = []
        for c in cols:
            v = r.get(c, "")
            v = "" if v is None else str(v)
            v = v.replace("|", "\\|").replace("\n", " ")
            cells.append(v)
        lines.append("| " + " | ".join(cells) + " |")
    lines.append("")
    (ROOT / "PROVENANCE.md").write_text("\n".join(lines), encoding="utf-8")


class Throttle:
    """Simple politeness throttle (minimum seconds between calls)."""

    def __init__(self, min_interval: float = 1.0):
        self.min_interval = min_interval
        self._last = 0.0

    def wait(self) -> None:
        dt = time.time() - self._last
        if dt < self.min_interval:
            time.sleep(self.min_interval - dt)
        self._last = time.time()


def env_key(name: str) -> str | None:
    v = os.environ.get(name)
    return v.strip() if v else None
