"""Scan the `input/` tree, classify runs as new / partial / done / invalid."""

from __future__ import annotations

import json
from collections import Counter
from dataclasses import dataclass
from pathlib import Path


@dataclass
class Run:
    name: str
    refs_dir: Path
    fastq_dir: Path
    output_dir: Path
    status: str  # "new" | "done" | "partial" | "invalid"
    summary: str  # e.g. "12 GREEN · 4 flagged" or "missing fastq/"


def _has_refs(refs_dir: Path) -> bool:
    if not refs_dir.is_dir():
        return False
    for p in refs_dir.iterdir():
        if p.suffix.lower() in (".gb", ".dna"):
            return True
    return False


def _has_fastqs(fastq_dir: Path) -> bool:
    if not fastq_dir.is_dir():
        return False
    for _ in fastq_dir.glob("*.fastq.gz"):
        return True
    return False


def _summarize_done(plate_map_json: Path) -> str:
    """Read plate_map.json and return 'N GREEN · M flagged' style summary."""
    try:
        data = json.loads(plate_map_json.read_text())
    except (OSError, json.JSONDecodeError):
        return "done (unreadable summary)"

    wells = data if isinstance(data, list) else data.get("wells", [])
    if not isinstance(wells, list):
        return "done"

    verdicts = Counter()
    for w in wells:
        if isinstance(w, dict):
            v = w.get("verdict")
            if v:
                verdicts[str(v)] += 1

    if not verdicts:
        return "done"

    green = verdicts.get("GREEN", 0)
    flagged = sum(c for v, c in verdicts.items() if v != "GREEN" and v != "NO_DATA")
    parts = [f"{green} GREEN"]
    if flagged:
        parts.append(f"{flagged} flagged")
    return " · ".join(parts)


def scan_input_dir(input_root: Path, output_root: Path) -> list[Run]:
    """Walk input_root/*/, return one Run per subdirectory.

    Status is derived from presence of output_root/<name>/plate_map.json:
      - missing refs/ or fastq/ → "invalid"
      - output/<name>/plate_map.json exists → "done" (summary from JSON)
      - output/<name>/ exists but no plate_map.json → "partial"
      - otherwise → "new"
    """
    input_root = Path(input_root)
    output_root = Path(output_root)
    runs: list[Run] = []

    if not input_root.is_dir():
        return runs

    for entry in sorted(input_root.iterdir()):
        if not entry.is_dir() or entry.name.startswith("."):
            continue

        refs_dir = entry / "refs"
        fastq_dir = entry / "fastq"
        output_dir = output_root / entry.name
        plate_map_json = output_dir / "plate_map.json"

        if not _has_refs(refs_dir) or not _has_fastqs(fastq_dir):
            missing = []
            if not _has_refs(refs_dir):
                missing.append("refs/")
            if not _has_fastqs(fastq_dir):
                missing.append("fastq/")
            status = "invalid"
            summary = "missing " + " or ".join(missing)
        elif plate_map_json.is_file():
            status = "done"
            summary = _summarize_done(plate_map_json)
        elif output_dir.is_dir():
            status = "partial"
            summary = "output/ exists but no plate_map.json"
        else:
            status = "new"
            summary = "ready to run"

        runs.append(
            Run(
                name=entry.name,
                refs_dir=refs_dir,
                fastq_dir=fastq_dir,
                output_dir=output_dir,
                status=status,
                summary=summary,
            )
        )

    return runs
