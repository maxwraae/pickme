"""Headless pipeline runner — shared by the TUI and the `pickme run` CLI.

Mirrors the per-well loop from `scripts/identify.py` but with a progress
callback so the TUI can show live updates.
"""

from __future__ import annotations

import os
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Callable, Optional

from tn5verify.refs import load_folder
from tn5verify.grouping import build_groups
from tn5verify.characterize import annotate
from tn5verify.target import build_multi_fasta
from tn5verify.align import run_bwa
from tn5verify.identify import classify
from tn5verify.integrity import evaluate
from tn5verify.features import annotate_gaps
from tn5verify.insert_scan import scan as insert_scan
from tn5verify.render.json_dump import write as write_json
from tn5verify.render.well_report import write as write_report
from tn5verify.render.plate import write_xlsx

from pickme.runs import Run


# Progress signature: (well_id, stage, payload)
# stages: "start", "align", "identify", "done", "summary"
ProgressCallback = Callable[[str, str, dict], None]


_WELL_RE = re.compile(r"([A-H][0-9]{1,2})")


def _discover_wells(fastq_dir: Path) -> list[tuple[str, str, str]]:
    """Return list of (well_id, r1, r2) from fastq_dir."""
    wells: list[tuple[str, str, str]] = []
    seen: dict[str, str] = {}

    r1_files = sorted(
        list(fastq_dir.glob("*_R1*.fastq.gz"))
        + list(fastq_dir.glob("*_R1*.fastq"))
    )

    for r1_path in r1_files:
        r2_name = r1_path.name.replace("_R1", "_R2")
        r2_path = r1_path.parent / r2_name
        if not r2_path.exists():
            print(
                f"  WARNING: R2 not found for {r1_path.name}, skipping.",
                file=sys.stderr,
            )
            continue

        matches = _WELL_RE.findall(r1_path.stem)
        if not matches:
            print(
                f"  WARNING: no well_id found in {r1_path.name}, skipping.",
                file=sys.stderr,
            )
            continue
        well_id = matches[-1]

        if well_id in seen:
            print(
                f"  WARNING: well_id {well_id} already seen ({seen[well_id]}), "
                f"skipping {r1_path.name}.",
                file=sys.stderr,
            )
            continue

        seen[well_id] = r1_path.name
        wells.append((well_id, str(r1_path), str(r2_path)))

    return wells


def _noop(well_id: str, stage: str, payload: dict) -> None:
    pass


def run_pipeline(
    run: Run,
    scan_mode: str = "quick",
    insert_region: Optional[tuple[int, int]] = None,
    threads: Optional[int] = None,
    on_progress: Optional[ProgressCallback] = None,
) -> dict:
    """Run the full Tn5 verification pipeline on one Run.

    Writes plate_map.{xlsx,txt,json} into run.output_dir.
    Returns the verdict histogram as a dict.
    """
    cb: ProgressCallback = on_progress or _noop

    if threads is None:
        threads = max(1, (os.cpu_count() or 2) // 2)

    run.output_dir.mkdir(parents=True, exist_ok=True)
    bam_root = run.output_dir / "bam"
    bam_root.mkdir(parents=True, exist_ok=True)

    wells = _discover_wells(run.fastq_dir)
    if not wells:
        raise RuntimeError(
            f"No wells discovered under {run.fastq_dir}. "
            "Check filenames match *_R1*.fastq.gz with a well ID."
        )

    cb("", "start", {"n_wells": len(wells), "fastq_dir": str(run.fastq_dir)})

    # Stage 1: Load references
    refs = load_folder(run.refs_dir)

    # Stage 2: Group references
    groups = build_groups(refs)

    # Stage 3: Characterize groups (backbone + variant regions)
    groups = annotate(groups, refs)

    # Stage 4: Build alignment target
    target_fa = build_multi_fasta(groups, refs, bam_root)

    # Stage 5: Align all wells
    all_stats = {}
    for well_id, r1, r2 in wells:
        well_dir = bam_root / well_id
        well_dir.mkdir(exist_ok=True)
        stats = run_bwa(Path(r1), Path(r2), target_fa, well_dir, threads=threads)
        all_stats[well_id] = stats
        cb(
            well_id,
            "align",
            {
                "mapped_reads": stats.mapped_reads,
                "total_reads": stats.total_reads,
            },
        )

    # Parse insert region
    scan_start: Optional[int] = None
    scan_end: Optional[int] = None
    if scan_mode == "insert":
        if not insert_region:
            raise ValueError(
                "scan_mode='insert' requires insert_region=(start, end)."
            )
        scan_start, scan_end = insert_region

    ref_by_name = {r.name: r for r in refs}

    # Stage 6: Identify and evaluate
    results = []
    for well_id, _, _ in wells:
        stats = all_stats[well_id]
        result = classify(stats.bam_path, groups, refs)
        result.well_id = well_id
        result = evaluate(result, stats.bam_path, groups, refs)
        annotate_gaps(result, refs)

        if (
            result.called_member
            and result.called_member in ref_by_name
            and scan_mode != "quick"
        ):
            result.insert_scan = insert_scan(
                stats.bam_path,
                result.called_member,
                ref_by_name[result.called_member],
                mode=scan_mode,
                scan_start=scan_start,
                scan_end=scan_end,
            )

        results.append(result)
        flagged = (
            len(result.insert_scan.mutations) if result.insert_scan else 0
        )
        cb(
            well_id,
            "identify",
            {
                "verdict": str(result.verdict),
                "called_member": result.called_member or "unknown",
                "flagged": flagged,
            },
        )

    # Stage 7: Write outputs
    out_xlsx = run.output_dir / "plate_map.xlsx"
    out_txt = run.output_dir / "plate_map.txt"
    out_json = run.output_dir / "plate_map.json"
    write_xlsx(results, out_xlsx)
    write_report(results, out_txt)
    write_json(results, out_json)

    verdict_counts = Counter(str(r.verdict) for r in results)
    cb(
        "",
        "summary",
        {
            "verdicts": dict(verdict_counts),
            "xlsx": str(out_xlsx),
            "txt": str(out_txt),
            "json": str(out_json),
        },
    )

    return dict(verdict_counts)
