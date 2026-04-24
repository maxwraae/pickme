#!/usr/bin/env python3
"""Tn5 construct verification pipeline."""

import argparse
import os
import csv
import re
import sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

# Add project root to sys.path so tn5verify package is importable
sys.path.insert(0, str(Path(__file__).parent.parent))

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


def _parse_args():
    parser = argparse.ArgumentParser(
        description="Tn5 construct verification pipeline."
    )
    parser.add_argument(
        "--fastq", required=True,
        help="Path to folder containing paired-end FASTQ files."
    )
    parser.add_argument(
        "--refs", required=True,
        help="Path to folder containing .gb/.dna reference files."
    )
    parser.add_argument(
        "--output", required=True,
        help="Path for output .xlsx plate map."
    )
    parser.add_argument(
        "--well-map",
        help="CSV file with columns well_id,fastq_r1,fastq_r2 (overrides filename inference)."
    )
    parser.add_argument(
        "--threads", type=int,
        default=max(1, os.cpu_count() // 2),
        help="Alignment threads (default: half of CPU count)."
    )
    parser.add_argument(
        "--run-dir",
        help="Directory for BAM outputs; defaults to run_{fastq_folder_name}/ next to the fastq folder."
    )
    parser.add_argument(
        "--scan-mode",
        choices=["quick", "insert", "whole"],
        default="quick",
        help=(
            "Mutation scan scope. "
            "'quick' (default): only signature positions that distinguish group members. "
            "'insert': scan a region specified with --insert-region. "
            "'whole': scan the entire called reference. "
            "All modes apply the same ≥3-read / ≥1%% noise floor. "
            "Whole adds ~1-2 min per 96-well plate."
        ),
    )
    parser.add_argument(
        "--insert-region",
        help=(
            "Native-coord region to scan when --scan-mode=insert, format 'START-END' "
            "(0-based, half-open). Example: --insert-region 500-3500."
        ),
    )
    return parser.parse_args()


def _discover_wells_from_folder(fastq_dir: Path) -> list[tuple[str, str, str]]:
    """
    Scan fastq_dir for *_R1*.fastq.gz (or *_R1*.fastq), pair with R2,
    extract well_id using the LAST match of ([A-H][0-9]{1,2}) in the filename.
    Returns list of (well_id, r1_path, r2_path).
    """
    wells: list[tuple[str, str, str]] = []
    seen_well_ids: dict[str, str] = {}  # well_id -> r1 filename (for dup detection)

    # Collect all R1 files
    r1_files = sorted(
        list(fastq_dir.glob("*_R1*.fastq.gz")) +
        list(fastq_dir.glob("*_R1*.fastq"))
    )

    well_pattern = re.compile(r"([A-H][0-9]{1,2})")

    for r1_path in r1_files:
        # Find R2 partner
        r2_name = r1_path.name.replace("_R1", "_R2")
        r2_path = r1_path.parent / r2_name
        if not r2_path.exists():
            print(f"  WARNING: R2 not found for {r1_path.name}, skipping.", file=sys.stderr)
            continue

        # Extract well_id: last match in filename stem
        matches = well_pattern.findall(r1_path.stem)
        if not matches:
            print(f"  WARNING: no well_id found in {r1_path.name}, skipping.", file=sys.stderr)
            continue
        well_id = matches[-1]  # last match

        # Duplicate detection
        if well_id in seen_well_ids:
            print(
                f"  WARNING: well_id {well_id} already seen ({seen_well_ids[well_id]}), "
                f"skipping {r1_path.name}.",
                file=sys.stderr,
            )
            continue

        seen_well_ids[well_id] = r1_path.name
        wells.append((well_id, str(r1_path), str(r2_path)))

    return wells


def _discover_wells_from_csv(csv_path: Path) -> list[tuple[str, str, str]]:
    """Read well map CSV with columns well_id,fastq_r1,fastq_r2."""
    wells = []
    with open(csv_path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            wells.append((row["well_id"], row["fastq_r1"], row["fastq_r2"]))
    return wells


def main():
    args = _parse_args()

    fastq_dir = Path(args.fastq)
    refs_dir = Path(args.refs)
    threads = args.threads

    # Determine run_dir
    if args.run_dir:
        run_dir = Path(args.run_dir)
    else:
        run_dir = fastq_dir.parent / f"run_{fastq_dir.name}"
    run_dir.mkdir(parents=True, exist_ok=True)

    # Well discovery
    if args.well_map:
        wells = _discover_wells_from_csv(Path(args.well_map))
    else:
        wells = _discover_wells_from_folder(fastq_dir)

    if not wells:
        print("ERROR: no wells discovered. Check --fastq folder or --well-map.", file=sys.stderr)
        sys.exit(1)

    # ── Stage 1: Load references ────────────────────────────────────────────────
    print(f"Loading references from {refs_dir}...")
    refs = load_folder(refs_dir)
    print(f"  {len(refs)} references loaded")

    # ── Stage 2: Group references ───────────────────────────────────────────────
    print("Grouping references...")
    groups = build_groups(refs)
    print(f"  {len(groups)} groups ({sum(len(g.members) > 1 for g in groups)} multi-member)")

    # ── Stage 3: Characterize groups ────────────────────────────────────────────
    print("Characterizing groups (backbone + variant regions)...")
    groups = annotate(groups, refs)

    # ── Stage 4: Build alignment target ─────────────────────────────────────────
    print(f"Building combined alignment target in {run_dir}...")
    target_fa = build_multi_fasta(groups, refs, run_dir)

    # ── Stage 5: Align all wells ─────────────────────────────────────────────────
    print(f"Aligning {len(wells)} wells (threads={threads})...")
    all_stats = {}
    for well_id, r1, r2 in wells:
        well_dir = run_dir / well_id
        well_dir.mkdir(exist_ok=True)
        print(f"  {well_id}...", end=" ", flush=True)
        stats = run_bwa(Path(r1), Path(r2), target_fa, well_dir, threads=threads)
        all_stats[well_id] = stats
        print(f"{stats.mapped_reads}/{stats.total_reads} mapped")

    # Parse --insert-region once up front if scan-mode=insert
    scan_start: int | None = None
    scan_end: int | None = None
    if args.scan_mode == "insert":
        if not args.insert_region:
            print("ERROR: --scan-mode=insert requires --insert-region START-END", file=sys.stderr)
            sys.exit(1)
        try:
            s, e = args.insert_region.split("-", 1)
            scan_start, scan_end = int(s), int(e)
        except ValueError:
            print(f"ERROR: --insert-region must be 'START-END', got '{args.insert_region}'", file=sys.stderr)
            sys.exit(1)

    ref_by_name = {r.name: r for r in refs}

    # ── Stage 6: Identify and evaluate ─────────────────────────────────────────
    print(f"Identifying and evaluating wells (scan-mode={args.scan_mode})...")
    results = []
    for well_id, _, _ in wells:
        stats = all_stats[well_id]
        result = classify(stats.bam_path, groups, refs)
        result.well_id = well_id
        result = evaluate(result, stats.bam_path, groups, refs)
        annotate_gaps(result, refs)

        # Full-reference mutation scan on the called member (quick mode = no-op)
        if result.called_member and result.called_member in ref_by_name and args.scan_mode != "quick":
            result.insert_scan = insert_scan(
                stats.bam_path,
                result.called_member,
                ref_by_name[result.called_member],
                mode=args.scan_mode,
                scan_start=scan_start,
                scan_end=scan_end,
            )

        results.append(result)
        flagged = len(result.insert_scan.mutations) if result.insert_scan else 0
        extra = f"  scan:{flagged} flagged" if args.scan_mode != "quick" else ""
        print(f"  {well_id}: {result.verdict} ({result.called_member or 'unknown'}){extra}")

    # ── Stage 7: Write outputs ───────────────────────────────────────────────────
    print("Writing outputs...")
    out_path = Path(args.output)
    write_xlsx(results, out_path)
    write_report(results, out_path.with_suffix(".txt"))
    write_json(results, out_path.with_suffix(".json"))
    print(f"Done. Plate map: {out_path}")

    # Summary
    from collections import Counter
    verdict_counts = Counter(r.verdict for r in results)
    print("\nVerdict summary:")
    for verdict, count in sorted(verdict_counts.items()):
        print(f"  {verdict}: {count}")


if __name__ == "__main__":
    main()
