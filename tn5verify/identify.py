from __future__ import annotations

from pathlib import Path

import pysam

from .types import Group, Reference, WellResult

_PAD = 300  # bases prepended to sequences >= 500 bp


def classify(
    bam_path: Path,
    groups: list[Group],
    refs: list[Reference],
) -> WellResult:
    """Classify a well BAM into a WellResult (M3 — variant/backbone covered; no integrity)."""

    # ------------------------------------------------------------------
    # Build contig → group-index lookup
    # ------------------------------------------------------------------
    ref_by_name = {r.name: r for r in refs}
    contig_to_group: dict[str, str] = {}
    for g_idx, group in enumerate(groups):
        for member in group.members:
            contig_to_group[member] = str(g_idx)

    # ------------------------------------------------------------------
    # Step 1: count reads per group
    # ------------------------------------------------------------------
    bam = pysam.AlignmentFile(str(bam_path), "rb")

    # Normalise BAM reference names — pysam strips after first whitespace,
    # so bam.references should already match ref.name.  Log a warning if not.
    bam_refs = set(bam.references)
    for name in contig_to_group:
        if name not in bam_refs:
            # Try stripping — should be a no-op but be safe
            pass  # contig_to_group keys are already ref.name

    total_reads: int = bam.mapped + bam.unmapped
    mapped_reads: int = 0
    group_totals: dict[str, int] = {}

    for read in bam.fetch(until_eof=True):
        if read.is_unmapped:
            continue
        mapped_reads += 1
        contig: str = read.reference_name  # type: ignore[assignment]
        g_key = contig_to_group.get(contig)
        if g_key is not None:
            group_totals[g_key] = group_totals.get(g_key, 0) + 1

    # ------------------------------------------------------------------
    # Step 2: early exit — too few reads
    # ------------------------------------------------------------------
    total_mapped = sum(group_totals.values())
    if total_mapped < 100:
        bam.close()
        return WellResult(
            well_id="",
            total_reads=total_reads,
            mapped_reads=mapped_reads,
            winning_group=None,
            called_member=None,
            coverage_span=0.0,
            verdict="NO_DATA",
            verdict_reason="fewer than 100 reads mapped",
        )

    # ------------------------------------------------------------------
    # Step 3: winning group
    # ------------------------------------------------------------------
    winning_group_idx = max(group_totals, key=lambda k: group_totals[k])
    winning_group = groups[int(winning_group_idx)]

    # ------------------------------------------------------------------
    # Step 4: coverage_span over backbone intervals
    # ------------------------------------------------------------------
    total_backbone_pos = 0
    covered_backbone_pos = 0

    for member in winning_group.members:
        intervals = winning_group.backbone_intervals.get(member, ())
        orig_len = len(ref_by_name[member].sequence)
        is_padded = orig_len >= 500
        offset = _PAD if is_padded else 0

        for nat_start, nat_end in intervals:
            padded_start = nat_start + offset
            padded_end = nat_end + offset
            region_len = padded_end - padded_start
            total_backbone_pos += region_len

            try:
                acgt = bam.count_coverage(
                    member,
                    padded_start,
                    padded_end,
                    quality_threshold=0,
                )
                # acgt is a tuple of 4 arrays (A, C, G, T)
                for pos_idx in range(region_len):
                    depth = sum(acgt[b][pos_idx] for b in range(4))
                    if depth > 0:
                        covered_backbone_pos += 1
            except (ValueError, KeyError):
                # Contig not present in BAM index — treat as uncovered
                pass

    coverage_span = covered_backbone_pos / total_backbone_pos if total_backbone_pos > 0 else 0.0

    if coverage_span < 0.30:
        bam.close()
        return WellResult(
            well_id="",
            total_reads=total_reads,
            mapped_reads=mapped_reads,
            winning_group=winning_group_idx,
            called_member=None,
            coverage_span=coverage_span,
            verdict="RED_UNKNOWN",
            verdict_reason="coverage_span below 30%",
        )

    # ------------------------------------------------------------------
    # Step 5: member selection via variant regions
    # ------------------------------------------------------------------
    variant_regions = winning_group.variant_regions

    # Tally wins per member
    member_wins: dict[str, int] = {m: 0 for m in winning_group.members}

    for vr in variant_regions:
        region_counts: dict[str, int] = {}
        for member, (nat_start, nat_end) in vr.member_spans.items():
            orig_len = len(ref_by_name[member].sequence)
            is_padded = orig_len >= 500
            offset = _PAD if is_padded else 0
            padded_start = nat_start + offset
            padded_end = nat_end + offset

            try:
                count = bam.count(member, padded_start, padded_end)
            except (ValueError, KeyError):
                count = 0
            region_counts[member] = count

        region_total = sum(region_counts.values())

        # MIX detection: ≥2 members each have ≥20% of region reads AND region_total ≥ 5
        if region_total >= 5:
            high_members = [
                m for m, c in region_counts.items()
                if c / region_total >= 0.20
            ]
            if len(high_members) >= 2:
                # Identify best-mapping member so evaluate() can run per-position analysis.
                # evaluate() may reclassify to RED_MUT if a novel mutation caused the apparent mix.
                mix_called_member = max(region_counts, key=lambda m: region_counts[m])
                bam.close()
                return WellResult(
                    well_id="",
                    total_reads=total_reads,
                    mapped_reads=mapped_reads,
                    winning_group=winning_group_idx,
                    called_member=mix_called_member,
                    coverage_span=coverage_span,
                    verdict="RED_MIX",
                    verdict_reason="mixed signal across variant region",
                )

        # Winning member for this region
        if region_total > 0:
            best_member = max(region_counts, key=lambda m: region_counts[m])
            member_wins[best_member] = member_wins.get(best_member, 0) + 1

    # If no variant regions, the single/only member wins by default
    if not variant_regions:
        called_member = winning_group.members[0]
    else:
        called_member = max(member_wins, key=lambda m: member_wins[m])

    bam.close()

    # ------------------------------------------------------------------
    # Step 6: return partial WellResult (verdict filled by integrity.py)
    # ------------------------------------------------------------------
    return WellResult(
        well_id="",
        total_reads=total_reads,
        mapped_reads=mapped_reads,
        winning_group=winning_group_idx,
        called_member=called_member,
        coverage_span=coverage_span,
        verdict="NO_DATA",  # placeholder — integrity.py sets final verdict
        verdict_reason="",
    )
