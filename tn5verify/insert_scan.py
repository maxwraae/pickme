"""
insert_scan.py — full-reference mutation scan.

Applies the same noise-floor algorithm used at discriminating positions
(≥3 reads AND ≥1% disagreeing with the reference base), but across a
user-chosen range: the whole plasmid or an arbitrary insert region.

Shares behavior with integrity.py's pileup: direct read iteration that
does NOT filter mapq=0 reads (competitive BWA assigns mapq=0 to reads
that tie-break between sibling contigs; those carry real signal).
"""
from __future__ import annotations

from pathlib import Path
from typing import Literal

import pysam

from .types import InsertScanResult, MutationHit, Reference

_PAD = 300  # bases prepended for refs >= 500 bp
_BASES = "ACGT"

# Noise floor — identical to per-position classification in integrity.py
_MIN_READS = 3
_MIN_FRACTION = 0.01
# We only trust a call if total coverage at the position is high enough
_MIN_COVERAGE = 10


def _is_padded(ref: Reference) -> bool:
    return len(ref.sequence) >= 500


def _count_bases_at_range(
    bam: pysam.AlignmentFile,
    contig: str,
    padded_start: int,
    padded_end: int,
) -> tuple[list[int], list[int], list[int], list[int]]:
    """Return per-position A, C, G, T counts across [padded_start, padded_end).

    Uses pysam.count_coverage which is vectorised in C. Accepts mapq=0 reads
    via read_callback='nofilter' — safe because the BAM is already deduped and
    BWA doesn't emit unmapped/secondary/supplementary by default at this stage.
    Base-quality threshold matches our single-position pileup (≥20).
    """
    acgt = bam.count_coverage(
        contig,
        padded_start,
        padded_end,
        quality_threshold=20,
        read_callback="nofilter",
    )
    # Order returned by pysam: A, C, G, T
    return list(acgt[0]), list(acgt[1]), list(acgt[2]), list(acgt[3])


def scan(
    bam_path: Path,
    member: str,
    ref: Reference,
    mode: Literal["quick", "insert", "whole"],
    scan_start: int | None = None,
    scan_end: int | None = None,
    backbone_intervals: tuple[tuple[int, int], ...] = (),
) -> InsertScanResult:
    """Pileup every position in the chosen range; flag disagreements with ref.

    Parameters
    ----------
    mode : 'quick' | 'insert' | 'whole'
        'quick'   — returns an empty scan (signature-position detail is the scan).
        'insert'  — scan [scan_start, scan_end) only.
        'whole'   — scan the entire reference, native coords [0, len(ref)).
    scan_start, scan_end : int | None
        Required when mode='insert'. Native-coord half-open interval.
    backbone_intervals : reserved for future use (not consulted today).

    Returns
    -------
    InsertScanResult with the list of flagged positions and scan statistics.
    """
    seq = ref.sequence.upper()
    seq_len = len(seq)
    padded = _is_padded(ref)
    pad_offset = _PAD if padded else 0

    if mode == "quick":
        return InsertScanResult(
            mode="quick",
            scan_start=0,
            scan_end=0,
            positions_scanned=0,
            mean_coverage=0.0,
        )

    if mode == "whole":
        scan_start = 0
        scan_end = seq_len
    elif mode == "insert":
        if scan_start is None or scan_end is None:
            raise ValueError("scan_start and scan_end required for mode='insert'")
        scan_start = max(0, scan_start)
        scan_end = min(seq_len, scan_end)
    else:
        raise ValueError(f"unknown scan mode: {mode}")

    if scan_end <= scan_start:
        return InsertScanResult(
            mode=mode,
            scan_start=scan_start,
            scan_end=scan_end,
            positions_scanned=0,
            mean_coverage=0.0,
        )

    padded_start = scan_start + pad_offset
    padded_end = scan_end + pad_offset

    bam = pysam.AlignmentFile(str(bam_path), "rb")
    try:
        a_counts, c_counts, g_counts, t_counts = _count_bases_at_range(
            bam, member, padded_start, padded_end
        )
    finally:
        bam.close()

    mutations: list[MutationHit] = []
    total_depth = 0
    n_positions = scan_end - scan_start

    for i in range(n_positions):
        native_pos = scan_start + i
        expected = seq[native_pos]
        counts = {
            "A": a_counts[i],
            "C": c_counts[i],
            "G": g_counts[i],
            "T": t_counts[i],
        }
        depth = counts["A"] + counts["C"] + counts["G"] + counts["T"]
        total_depth += depth

        if depth < _MIN_COVERAGE:
            continue  # not enough evidence to call a mutation here
        if expected not in _BASES:
            continue  # skip N, ambiguity codes

        # Look for any real non-expected allele (majority or minor)
        flagged_base = None
        for base, count in counts.items():
            if base == expected:
                continue
            if count < _MIN_READS:
                continue
            if count / depth < _MIN_FRACTION:
                continue
            # Real disagreement — flag this position
            flagged_base = base
            break

        if flagged_base is None:
            continue

        # Determine majority base for the report (can be expected if heterozygous)
        majority_base = max(counts, key=lambda b: counts[b])
        majority_fraction = counts[majority_base] / depth

        # Feature overlap — look up .gb features touching this native position
        affected: list[str] = []
        for feature in ref.features:
            if feature.start <= native_pos < feature.end:
                affected.append(f"{feature.label} ({feature.kind})")

        mutations.append(MutationHit(
            pos=native_pos,
            expected=expected,
            A=counts["A"],
            C=counts["C"],
            G=counts["G"],
            T=counts["T"],
            majority_base=majority_base,
            majority_fraction=majority_fraction,
            affected_features=affected,
        ))

    mean_coverage = total_depth / n_positions if n_positions > 0 else 0.0

    return InsertScanResult(
        mode=mode,
        scan_start=scan_start,
        scan_end=scan_end,
        positions_scanned=n_positions,
        mean_coverage=mean_coverage,
        mutations=mutations,
    )
