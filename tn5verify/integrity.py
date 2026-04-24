"""
integrity.py — M4: per-well integrity evaluation.

Takes a WellResult from classify() and fills in:
  - variant_regions (list[RegionResult])
  - backbone_windows (list[WindowResult])
  - verdict / verdict_reason
"""
from __future__ import annotations

import subprocess
import tempfile
from collections import Counter
from pathlib import Path
from typing import Literal

import pandas as pd
import pysam

from .types import Group, Reference, RegionResult, WellResult, WindowResult

_PAD = 300  # prepended bases for refs >= 500 bp


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _find_tool(name: str) -> str:
    import shutil
    found = shutil.which(name)
    if found:
        return found
    for prefix in ["/opt/homebrew/bin", "/usr/local/bin"]:
        candidate = Path(prefix) / name
        if candidate.exists():
            return str(candidate)
    raise FileNotFoundError(f"Could not find {name} in PATH or common locations")


def _is_padded(ref_by_name: dict[str, Reference], member: str) -> bool:
    return len(ref_by_name[member].sequence) >= 500


def _padded_pos(native: int, padded: bool) -> int:
    return native + _PAD if padded else native


# ---------------------------------------------------------------------------
# §2.1–2.4  Per-position pileup
# ---------------------------------------------------------------------------

_BASES = frozenset("ACGT")


def _pileup_position(
    bam: pysam.AlignmentFile,
    contig: str,
    padded_pos: int,
) -> dict:
    """Return {A, T, C, G, pos} counts at a single padded position.

    Uses direct read iteration rather than pysam.pileup() because pysam.pileup
    filters out mapq=0 reads regardless of min_mapping_quality. Competitive
    BWA assigns mapq=0 to reads that tie-break between sibling contigs, and
    those are exactly the reads that carry the discriminating signal at a
    variant position. Base quality ≥20 still screens sequencing noise.
    """
    counts: Counter[str] = Counter()
    for read in bam.fetch(contig, padded_pos, padded_pos + 1):
        if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate:
            continue
        if read.query_sequence is None:
            continue
        # Skip if read doesn't actually span padded_pos (fetch may return overlapping mates)
        if read.reference_start > padded_pos or read.reference_end is None or read.reference_end <= padded_pos:
            continue
        # Find the query-position that aligns to padded_pos
        for qp, rp in read.get_aligned_pairs(matches_only=True):
            if rp == padded_pos:
                b = read.query_sequence[qp].upper()
                if b in _BASES:
                    bq = read.query_qualities[qp] if read.query_qualities is not None else 40
                    if bq >= 20:
                        counts[b] += 1
                break
    return {"pos": padded_pos, "A": counts["A"], "T": counts["T"],
            "C": counts["C"], "G": counts["G"]}


def _classify_position(
    pileup_dict: dict,
    called_member: str,
    native_pos: int,
    region: object,  # VariantRegion
) -> str:
    """
    Classify a single position given pileup counts and the VariantRegion context.

    Returns one of: CLEAN, THIN, UNCOVERED, MIXED, MUTATED
    """
    total = pileup_dict["A"] + pileup_dict["T"] + pileup_dict["C"] + pileup_dict["G"]

    if total < 3:
        return "UNCOVERED"

    base_counts = {b: pileup_dict[b] for b in "ACGT"}
    expected_base = region.discriminating_positions.get(called_member, {}).get(native_pos)
    alternatives = set(region.alternative_bases.get(called_member, {}).get(native_pos, ()))

    def classify_nonexpected(base: str) -> str:
        # Minor/majority allele that isn't the expected base at this discriminating position.
        # If another member expects it (same MSA column), it's MIXED (identity confusion).
        # Otherwise it's a novel base nobody ordered → MUTATED.
        return "MIXED" if base in alternatives else "MUTATED"

    # Homozygous mutation: majority differs from expected
    majority_base = max(base_counts, key=lambda b: base_counts[b])
    if expected_base is not None and majority_base != expected_base and base_counts[majority_base] > 0:
        return classify_nonexpected(majority_base)

    # Real minor alleles (consensus matches expected, but a sub-population of different base is present)
    worst_minor: str | None = None
    for base, count in base_counts.items():
        if expected_base is not None and base == expected_base:
            continue
        if count == 0:
            continue
        # Noise-floor: real iff count >= 3 AND fraction >= 1 %
        if count < 3 or count / total < 0.01:
            continue
        candidate = classify_nonexpected(base)
        if worst_minor is None or candidate == "MUTATED":
            worst_minor = candidate

    if worst_minor is not None:
        return worst_minor

    if total >= 10:
        return "CLEAN"
    return "THIN"


def _region_status_from_positions(
    statuses: list[str],
) -> Literal["trusted_clean", "thin_clean", "partial", "mixed", "mutated"]:
    """Aggregate per-position statuses into a region-level status."""
    if not statuses:
        return "partial"  # no positions sampled → treat as partial
    status_set = set(statuses)
    if "MUTATED" in status_set:
        return "mutated"
    if "MIXED" in status_set:
        return "mixed"
    if "UNCOVERED" in status_set:
        return "partial"
    if "THIN" in status_set:
        return "thin_clean"
    return "trusted_clean"


def _variant_verdict_from_regions(
    region_statuses: list[str],
) -> str:
    if not region_statuses:
        return "CLEAN"  # no regions → singleton → clean
    if "mutated" in region_statuses:
        return "MUT"
    if "mixed" in region_statuses:
        return "MIX"
    if "partial" in region_statuses:
        return "PARTIAL"
    if "thin_clean" in region_statuses:
        return "THIN"
    return "CLEAN"


# ---------------------------------------------------------------------------
# §2.5  Build RegionResults
# ---------------------------------------------------------------------------

def _multi_member_pileup(
    bam: pysam.AlignmentFile,
    group_members: tuple[str, ...],
    ref_by_name: dict[str, Reference],
    native_pos: int,
) -> dict:
    """Sum pileup counts at the same native position across every member's contig.

    Needed because competitive BWA alignment spreads reads spanning a discriminating
    position across multiple sibling contigs. A mutation at such a position
    would leave the called_member's contig empty there, hiding the mutation.
    Safe only when members share a coordinate system within the region
    (non-insertion size classes).
    """
    counts = {"A": 0, "T": 0, "C": 0, "G": 0, "pos": native_pos}
    for member in group_members:
        padded = _is_padded(ref_by_name, member)
        pp = _padded_pos(native_pos, padded)
        d = _pileup_position(bam, member, pp)
        for b in "ACGT":
            counts[b] += d[b]
    return counts


def _evaluate_region(
    bam: pysam.AlignmentFile,
    region_index: int,
    region,  # VariantRegion
    called_member: str,
    group_members: tuple[str, ...],
    ref_by_name: dict[str, Reference],
) -> RegionResult:
    padded = _is_padded(ref_by_name, called_member)
    disc_positions = region.discriminating_positions.get(called_member, {})

    if not disc_positions:
        # No discriminating positions for this member → clean by default
        return RegionResult(
            region_index=region_index,
            positions=[],
            region_status="trusted_clean",
        )

    # Select positions to sample
    if region.size_class == "point_or_short":
        sample_positions = list(disc_positions.keys())
    else:
        # long / insertion: sample up to 20 evenly-spaced positions
        positions_list = list(disc_positions.keys())
        n = len(positions_list)
        if n <= 20:
            sample_positions = positions_list
        else:
            indices = [int(i * (n - 1) / 19) for i in range(20)]
            sample_positions = [positions_list[i] for i in indices]

    # Multi-member pileup when coordinate systems align across the group.
    # Insertion regions have per-member coordinate offsets, so fall back
    # to single-member pileup on the called_member's contig.
    use_multi_member = region.size_class != "insertion" and len(group_members) > 1

    position_records = []
    pos_statuses = []

    for native_pos in sample_positions:
        if use_multi_member:
            pdict = _multi_member_pileup(bam, group_members, ref_by_name, native_pos)
        else:
            pp = _padded_pos(native_pos, padded)
            pdict = _pileup_position(bam, called_member, pp)
        status = _classify_position(pdict, called_member, native_pos, region)
        rec = dict(pdict)
        rec["pos"] = native_pos  # store as native coord in the record
        rec["status"] = status
        rec["expected"] = region.discriminating_positions.get(called_member, {}).get(native_pos, "?")
        position_records.append(rec)
        pos_statuses.append(status)

    region_status = _region_status_from_positions(pos_statuses)

    return RegionResult(
        region_index=region_index,
        positions=position_records,
        region_status=region_status,
    )


# ---------------------------------------------------------------------------
# §2.7  Backbone windows via mosdepth
# ---------------------------------------------------------------------------

def _run_mosdepth(bam_path: Path, tmp_dir: Path) -> pd.DataFrame | None:
    """
    Run mosdepth with --by 100 and return the regions BED as a DataFrame.
    Returns None if mosdepth is unavailable or output is missing.
    Tolerates non-zero exit (mosdepth sometimes crashes with SIGABRT on macOS
    but still writes valid output).
    """
    try:
        mosdepth = _find_tool("mosdepth")
    except FileNotFoundError:
        return None

    prefix = tmp_dir / "mosdepth"
    bed_path = tmp_dir / "mosdepth.regions.bed.gz"

    result = subprocess.run(
        [mosdepth, "--by", "100", "--no-per-base", "--mapq", "20",
         str(prefix), str(bam_path)],
        capture_output=True,
        text=True,
    )

    if not bed_path.exists():
        return None

    try:
        df = pd.read_csv(
            str(bed_path),
            sep="\t",
            header=None,
            names=["chrom", "start", "end", "coverage"],
            compression="gzip",
        )
        return df
    except Exception:
        return None


def _build_backbone_windows(
    bam_path: Path,
    called_member: str,
    group: Group,
    ref_by_name: dict[str, Reference],
) -> list[WindowResult]:
    """Run mosdepth and classify backbone windows for the called member."""
    with tempfile.TemporaryDirectory() as tmpdir:
        df = _run_mosdepth(bam_path, Path(tmpdir))

    if df is None or df.empty:
        return []

    # Filter to called member's contig
    member_df = df[df["chrom"] == called_member].copy()
    if member_df.empty:
        return []

    # Determine backbone intervals in padded coords
    padded = _is_padded(ref_by_name, called_member)
    raw_intervals = group.backbone_intervals.get(called_member, ())
    backbone_ranges: list[tuple[int, int]] = [
        (_padded_pos(s, padded), _padded_pos(e, padded))
        for s, e in raw_intervals
    ]

    # Filter windows to those overlapping any backbone interval
    def overlaps_backbone(row: pd.Series) -> bool:
        w_start, w_end = int(row["start"]), int(row["end"])
        for b_start, b_end in backbone_ranges:
            if w_start < b_end and w_end > b_start:
                return True
        return False

    backbone_df = member_df[member_df.apply(overlaps_backbone, axis=1)].copy()
    if backbone_df.empty:
        # No windows overlap backbone — return all windows on this contig
        backbone_df = member_df.copy()

    coverages = backbone_df["coverage"].values
    median_cov = float(pd.Series(coverages).median())

    windows: list[WindowResult] = []
    for _, row in backbone_df.iterrows():
        cov = float(row["coverage"])
        w_start = int(row["start"])
        w_end = int(row["end"])

        # Classification
        if median_cov == 0:
            # If median is 0, any window with ≥ 0 is "intact"
            status: Literal["intact", "warn", "broken"] = "intact"
        elif cov >= 0.20 * median_cov:
            status = "intact"
        else:
            status = "warn"  # final "broken" logic applied below

        windows.append(WindowResult(
            contig=called_member,
            start=w_start,
            end=w_end,
            mean_coverage=cov,
            status=status,
        ))

    # Post-process: mark "broken" for windows < 5% of median AND consecutive
    if median_cov > 0:
        broken_threshold = 0.05 * median_cov
        for i, w in enumerate(windows):
            if w.mean_coverage < broken_threshold:
                # Check if consecutive neighbor is also < threshold
                prev_broken = (i > 0 and windows[i - 1].mean_coverage < broken_threshold)
                next_broken = (i < len(windows) - 1 and windows[i + 1].mean_coverage < broken_threshold)
                if prev_broken or next_broken:
                    windows[i].status = "broken"

    return windows


def _backbone_verdict(windows: list[WindowResult]) -> str:
    if not windows:
        return "INTACT"
    statuses = {w.status for w in windows}
    if "broken" in statuses:
        return "BROKEN"
    warn_count = sum(1 for w in windows if w.status == "warn")
    if warn_count > 1:
        return "WARN"
    return "INTACT"


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def evaluate(
    result: WellResult,
    bam_path: Path,
    groups: list[Group],
    refs: list[Reference],
) -> WellResult:
    """
    Fill in variant_regions, backbone_windows, verdict, verdict_reason.
    Returns the same WellResult (mutated in place).
    """
    # Guards
    if result.verdict == "RED_UNKNOWN":
        return result
    if result.called_member is None:
        return result
    # Also skip if verdict is NO_DATA with no winning_group (truly no data)
    if result.verdict == "NO_DATA" and result.winning_group is None:
        return result

    ref_by_name = {r.name: r for r in refs}

    # Find the group corresponding to winning_group (stored as string index)
    group: Group | None = None
    if result.winning_group is not None:
        try:
            g_idx = int(result.winning_group)
            group = groups[g_idx]
        except (ValueError, IndexError):
            pass

    if group is None:
        return result

    # ------------------------------------------------------------------
    # §2.1–2.5  Variant region evaluation
    # ------------------------------------------------------------------
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    try:
        for i, vr in enumerate(group.variant_regions):
            rr = _evaluate_region(bam, i, vr, result.called_member, group.members, ref_by_name)
            result.variant_regions.append(rr)
    finally:
        bam.close()

    # §2.6 Variant verdict
    region_statuses = [rr.region_status for rr in result.variant_regions]
    variant_verdict = _variant_verdict_from_regions(region_statuses)

    # ------------------------------------------------------------------
    # §2.7  Backbone windows
    # ------------------------------------------------------------------
    result.backbone_windows = _build_backbone_windows(
        bam_path, result.called_member, group, ref_by_name
    )
    bb_verdict = _backbone_verdict(result.backbone_windows)

    # ------------------------------------------------------------------
    # §3  Final verdict (priority order)
    # ------------------------------------------------------------------
    incoming_verdict = result.verdict
    if variant_verdict == "MUT":
        result.verdict = "RED_MUT"
        result.verdict_reason = "mutated base in variant region"
    elif incoming_verdict == "RED_MIX":
        # Per-position analysis found no novel mutation → remain RED_MIX
        pass
    elif bb_verdict == "BROKEN":
        result.verdict = "RED_BROKEN"
        result.verdict_reason = "broken backbone coverage"
    elif variant_verdict == "PARTIAL":
        result.verdict = "YELLOW"
        result.verdict_reason = "partial coverage in variant region"
    elif variant_verdict == "THIN":
        result.verdict = "YELLOW"
        result.verdict_reason = "thin coverage"
    elif bb_verdict == "WARN":
        result.verdict = "YELLOW"
        result.verdict_reason = "backbone warning window"
    else:
        result.verdict = "GREEN"
        result.verdict_reason = ""

    return result
