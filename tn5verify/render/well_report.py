from __future__ import annotations

from pathlib import Path

from ..types import WellResult

_DIVIDER = "\n" + "─" * 60 + "\n"

_ACTIONS = {
    "GREEN":       "Use — ready to transfect.",
    "YELLOW":      "Re-sequence deeper, or use with caution.",
    "RED_MUT":     "Discard — mutation detected.",
    "RED_MIX":     "Re-streak and re-sequence.",
    "RED_BROKEN":  "Discard — backbone deletion.",
    "RED_UNKNOWN": "Check references — construct not identified.",
    "NO_DATA":     "Re-sequence — insufficient reads.",
}

# Display labels shown in the header next to the verdict code
_VERDICT_LABELS = {
    "GREEN":       "CLEAN",
    "RED_MUT":     "RED — MUTATION",
    "RED_MIX":     "RED — MIXED WELL",
    "RED_BROKEN":  "RED — BROKEN BACKBONE",
    "RED_UNKNOWN": "RED — UNKNOWN CONSTRUCT",
    "NO_DATA":     "NO DATA",
}

_BAR_CHARS = {"intact": "▇", "warn": "▄", "broken": "▁"}


def _header_verdict(r: WellResult) -> str:
    """Human-readable verdict label for the block title, with a sub-reason for YELLOW."""
    if r.verdict == "YELLOW":
        reason = r.verdict_reason.lower()
        if "thin" in reason:
            suffix = " (low confidence — not enough reads)"
        elif "partial" in reason:
            suffix = " (partial coverage in variant region)"
        elif "backbone" in reason:
            suffix = " (backbone warning window)"
        else:
            suffix = ""
        return f"YELLOW{suffix}"
    return _VERDICT_LABELS.get(r.verdict, r.verdict)


def _pos_line(pos: dict) -> str:
    """One row of the per-position table, spec-style: expected base first, then others."""
    expected = pos.get("expected", "?")
    status = pos["status"]
    depth = pos["A"] + pos["T"] + pos["C"] + pos["G"]

    if status == "UNCOVERED" or depth == 0:
        return f"  pos {pos['pos']:<6}  expect {expected}    (no reads)              ⚠  (UNCOVERED)"

    # Expected base first, other bases after in fixed alphabetical order
    others = [b for b in "ACGT" if b != expected]
    chunks = [f"{expected}:{pos[expected]:<4}"]
    chunks.extend(f"{b}:{pos[b]:<4}" for b in others)
    bases_str = " ".join(chunks)

    tick = "✓" if status in ("CLEAN", "THIN") else "⚠"
    trailing = ""
    if status == "THIN":
        trailing = f"  ({depth} reads — verified but thin)"
    elif status == "MIXED":
        trailing = "  (minor allele matches another known construct — mixed well)"
    elif status == "MUTATED":
        trailing = "  (minor allele matches no known construct — real mutation)"
    elif status == "CLEAN":
        trailing = f"  ({depth}×)"

    return f"  pos {pos['pos']:<6}  expect {expected}    {bases_str}  {tick}{trailing}"


def _format_block(r: WellResult) -> str:
    called = r.called_member or "UNKNOWN"
    total = r.total_reads
    mapped = r.mapped_reads
    mapped_pct = (mapped / total * 100) if total > 0 else 0.0

    all_positions = [p for region in r.variant_regions for p in region.positions]
    position_depths = [p["A"] + p["T"] + p["C"] + p["G"] for p in all_positions]
    n_confirmed = sum(1 for p in all_positions if p["status"] in ("CLEAN", "THIN"))
    n_thin_or_uncov = sum(1 for p in all_positions if p["status"] in ("THIN", "UNCOVERED"))
    n_uncov = sum(1 for p in all_positions if p["status"] == "UNCOVERED")
    n_total = len(all_positions)
    reads_at_positions = sum(position_depths)
    expected_reads = sum(p.get(p.get("expected", ""), 0) for p in all_positions)
    disagreeing = reads_at_positions - expected_reads

    lines = []
    # ── Title ──────────────────────────────────────────────────────────────
    lines.append(f"{r.well_id} — {called}  ·  {_header_verdict(r)}")
    lines.append("")

    # ── "Why we called this clone X" summary ──────────────────────────────
    if r.verdict == "NO_DATA":
        lines.append("Why we couldn't call this well:")
    elif r.verdict in ("RED_UNKNOWN",):
        lines.append("Why we couldn't identify this well:")
    else:
        lines.append(f"Why we called this clone {called}:")

    lines.append(f"  Total reads in well:     {total:,}")
    lines.append(f"  Mapped to plasmid:       {mapped:,} ({mapped_pct:.0f}%)")
    lines.append(f"  Coverage span:           {r.coverage_span:.0%}")
    if n_total > 0 and position_depths:
        mean_d = sum(position_depths) / len(position_depths)
        lines.append(
            f"  Variant-region coverage: mean {mean_d:.0f}×  "
            f"min {min(position_depths)}×  max {max(position_depths)}×"
        )
        lines.append(f"  Variant positions:       {n_confirmed}/{n_total} confirmed")
    lines.append("")

    # ── Per-position pileup table ─────────────────────────────────────────
    if r.variant_regions:
        if r.verdict == "GREEN":
            lines.append(
                f"At each of {called}'s {n_total} signature positions, the reads"
            )
            lines.append("overwhelmingly show the expected base:")
        elif r.verdict == "YELLOW":
            lines.append(
                f"Signature-position detail ({n_confirmed}/{n_total} confirmed, "
                f"{n_thin_or_uncov} with thin or no coverage):"
            )
        elif r.verdict in ("RED_MUT", "RED_MIX"):
            lines.append(f"Signature-position detail ({n_total} positions inspected):")
        else:
            lines.append(f"Signature-position detail ({n_total} positions):")
        lines.append("")

        for region in r.variant_regions:
            for pos in region.positions:
                lines.append(_pos_line(pos))

        if reads_at_positions > 0:
            lines.append("")
            lines.append(
                f"  Disagreeing reads across all positions: "
                f"{disagreeing} out of {reads_at_positions}."
            )
            # Only claim "clean" when we actually have enough coverage to back the claim.
            if disagreeing == 0 and r.verdict == "GREEN":
                lines.append(
                    "  (Sequencing noise typically gives ~1 per 1,000 reads. This well is clean.)"
                )
            elif disagreeing > 0:
                noise_rate = disagreeing / reads_at_positions * 1000
                lines.append(
                    f"  ({noise_rate:.1f} per 1,000 — sequencing noise typically gives ~1 per 1,000.)"
                )

        if r.verdict == "YELLOW" and n_thin_or_uncov > 0:
            lines.append("")
            lines.append(
                f"  {n_thin_or_uncov} of {n_total} positions have <10 reads "
                f"({n_uncov} with <3). We can't distinguish between"
            )
            lines.append(
                "  \"correct clone with low coverage\" and \"something else we can't see\"."
            )

    lines.append("")

    # ── Backbone block ────────────────────────────────────────────────────
    lines.append("Backbone:")
    if not r.backbone_windows:
        lines.append("  (not evaluated)")
    else:
        windows = r.backbone_windows
        statuses = {w.status for w in windows}
        if statuses == {"intact"}:
            lines.append("  Intact. All windows covered, no drop-outs.")
        else:
            # Render the coverage bar (28-char summary)
            n = len(windows)
            if n >= 28:
                indices = [int(i * (n - 1) / 27) for i in range(28)]
                bar_windows = [windows[i] for i in indices]
            else:
                bar_windows = windows
            bar = "".join(_BAR_CHARS.get(w.status, "▁") for w in bar_windows)
            lines.append(f"  Coverage profile (100 bp windows):")
            lines.append(f"  |{bar}|")

            # Find broken stretches and annotate under the bar
            broken = [w for w in windows if w.status == "broken"]
            warn = [w for w in windows if w.status == "warn"]

            if broken:
                start = broken[0].start
                end = broken[-1].end
                bp_lost = end - start
                lines.append(
                    f"  Broken stretch: {bp_lost} bp with zero coverage, position {start}–{end}."
                )
                # .gb feature overlap
                feats: list[str] = []
                for w in broken:
                    for f in w.affected_features:
                        if f not in feats:
                            feats.append(f)
                if feats:
                    lines.append("")
                    lines.append("  From .gb annotation, the deleted region overlaps:")
                    for f in feats:
                        lines.append(f"     • {f}")
                    lines.append("")
                    lines.append(
                        "  Interpretation: a chunk of the plasmid is missing. "
                        "If a selection marker (AmpR, KanR, etc.) is affected,"
                    )
                    lines.append(
                        "  this clone couldn't have survived antibiotic selection — "
                        "likely contamination or a mini-prep artifact."
                    )
            if warn and not broken:
                lines.append(
                    f"  {len(warn)} window(s) with coverage below 20 % of median — flagged but not broken."
                )

    lines.append("")

    # ── Insert mutation scan (whole plasmid or user-specified region) ─────
    if r.insert_scan is not None and r.insert_scan.mode != "quick":
        scan = r.insert_scan
        scope = (
            f"whole plasmid [{scan.scan_start}–{scan.scan_end})"
            if scan.mode == "whole" else
            f"insert region [{scan.scan_start}–{scan.scan_end})"
        )
        lines.append(
            f"Insert mutation scan ({scan.mode}): "
            f"{len(scan.mutations)} position(s) flagged across "
            f"{scan.positions_scanned:,} bp at mean {scan.mean_coverage:.0f}× coverage."
        )
        lines.append(f"  Scope: {scope}")
        if not scan.mutations:
            lines.append("  ✓ No mutations detected above noise floor (≥3 reads AND ≥1 %).")
        else:
            lines.append("")
            for m in scan.mutations:
                depth = m.A + m.C + m.G + m.T
                # Show expected base first, then others alphabetically
                others = " ".join(
                    f"{b}:{getattr(m, b)}" for b in "ACGT" if b != m.expected
                )
                pct = m.majority_fraction * 100
                lines.append(
                    f"  pos {m.pos:<6}  expect {m.expected}  "
                    f"{m.expected}:{getattr(m, m.expected)}  {others}  "
                    f"→ {m.majority_base} ({pct:.0f}% majority, {depth}× total)"
                )
                if m.affected_features:
                    # Collapse to unique, keep order
                    seen = []
                    for f in m.affected_features:
                        if f not in seen:
                            seen.append(f)
                    lines.append(f"     overlaps: {', '.join(seen)}")
            lines.append("")
            lines.append(
                "  Note: filter by feature context — mutations in selection markers "
                "(AmpR etc.) or origins usually don't affect your experiment;"
            )
            lines.append(
                "  mutations inside your expressed insert do."
            )
        lines.append("")

    lines.append(f"Verdict: {r.verdict}")
    if r.verdict_reason:
        lines.append(f"Reason:  {r.verdict_reason}")
    lines.append("")
    action = _ACTIONS.get(r.verdict, "Unknown verdict.")
    lines.append(f"Suggested action: {action}")

    return "\n".join(lines)


def write(results: list[WellResult], path: Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    blocks = [_format_block(r) for r in results]
    text = _DIVIDER.join(blocks)
    path.write_text(text, encoding="utf-8")
