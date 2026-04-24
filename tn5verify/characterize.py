from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

from .types import Group, Reference, VariantRegion


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


def annotate(groups: list[Group], refs: list[Reference]) -> list[Group]:
    ref_by_name = {r.name: r for r in refs}
    result: list[Group] = []

    for group in groups:
        if len(group.members) == 1:
            name = group.members[0]
            seq_len = len(ref_by_name[name].sequence)
            result.append(Group(
                members=group.members,
                backbone_intervals={name: ((0, seq_len),)},
                variant_regions=(),
            ))
        else:
            result.append(_annotate_multi(group, ref_by_name))

    return result


def _parse_fasta(text: str) -> dict[str, str]:
    seqs: dict[str, str] = {}
    current_name = None
    parts: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if current_name is not None:
                seqs[current_name] = "".join(parts)
            current_name = line[1:].split()[0]
            parts = []
        else:
            parts.append(line.strip())
    if current_name is not None:
        seqs[current_name] = "".join(parts)
    return seqs


def _annotate_multi(group: Group, ref_by_name: dict[str, Reference]) -> Group:
    mafft = _find_tool("mafft")
    members = list(group.members)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        input_fa = tmp / "input.fa"
        output_aln = tmp / "output.aln"

        with input_fa.open("w") as f:
            for name in members:
                seq = ref_by_name[name].sequence
                f.write(f">{name}\n{seq}\n")

        try:
            proc = subprocess.run(
                [mafft, "--auto", "--preservecase", str(input_fa)],
                capture_output=True,
                text=True,
                check=True,
            )
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"mafft failed: {e.stderr}") from e

        output_aln.write_text(proc.stdout)
        aligned = _parse_fasta(proc.stdout)

    # Ensure order matches members
    aln_seqs = [aligned[name] for name in members]
    aln_len = len(aln_seqs[0])

    # Track native positions
    native_pos = [0] * len(members)
    # Per-column: is_backbone
    backbone_cols: list[bool] = []

    for col in range(aln_len):
        chars = [aln_seqs[i][col] for i in range(len(members))]
        # Backbone: all identical and none is a gap
        if len(set(chars)) == 1 and chars[0] != "-":
            backbone_cols.append(True)
        else:
            backbone_cols.append(False)

    # Build backbone intervals and variant regions in native coords
    # We need to track native positions per column
    # native_pos[i] = native position at the START of column col (before processing col)
    native_positions_at_col: list[list[int]] = []
    nat = [0] * len(members)
    for col in range(aln_len):
        native_positions_at_col.append(list(nat))
        for i in range(len(members)):
            if aln_seqs[i][col] != "-":
                nat[i] += 1

    # Group contiguous backbone vs variant columns
    # Backbone intervals in native coords per member
    backbone_intervals: dict[str, list[tuple[int, int]]] = {name: [] for name in members}
    variant_regions: list[VariantRegion] = []

    col = 0
    while col < aln_len:
        if backbone_cols[col]:
            # Find run of backbone columns
            start_col = col
            while col < aln_len and backbone_cols[col]:
                col += 1
            end_col = col  # exclusive
            # Native coords for each member
            for i, name in enumerate(members):
                nat_start = native_positions_at_col[start_col][i]
                nat_end = native_positions_at_col[end_col][i] if end_col < aln_len else nat[i]
                # Only add if non-empty
                if nat_end > nat_start:
                    backbone_intervals[name].append((nat_start, nat_end))
        else:
            # Find run of variant columns
            start_col = col
            while col < aln_len and not backbone_cols[col]:
                col += 1
            end_col = col  # exclusive

            # Compute member spans and discriminating positions
            member_spans: dict[str, tuple[int, int]] = {}
            discriminating_positions: dict[str, dict[int, str]] = {}
            alternative_bases: dict[str, dict[int, tuple[str, ...]]] = {}
            has_insertion = False

            for i, name in enumerate(members):
                nat_start = native_positions_at_col[start_col][i]
                nat_end = native_positions_at_col[end_col][i] if end_col < aln_len else nat[i]
                member_spans[name] = (nat_start, nat_end)

                # Discriminating positions: native coords where this member differs from at least one other
                disc_pos: dict[int, str] = {}
                alt_bases: dict[int, tuple[str, ...]] = {}
                local_nat = native_positions_at_col[start_col][i]
                for c in range(start_col, end_col):
                    ch = aln_seqs[i][c]
                    chars = [aln_seqs[j][c] for j in range(len(members))]
                    # Check for gaps (insertion/deletion)
                    if "-" in chars:
                        has_insertion = True
                    if ch != "-":
                        # Is this position discriminating?
                        if len(set(chars)) > 1:
                            ch_up = ch.upper()
                            disc_pos[local_nat] = ch_up  # record expected base for this member
                            # Alternative bases: other members' non-gap, non-self bases at the same MSA column
                            others = {
                                chars[j].upper()
                                for j in range(len(members))
                                if j != i and chars[j] != "-" and chars[j].upper() != ch_up
                            }
                            alt_bases[local_nat] = tuple(sorted(others))
                        local_nat += 1
                    # else: gap for this member, no native position increment

                discriminating_positions[name] = disc_pos
                alternative_bases[name] = alt_bases

            # Size class
            if has_insertion:
                size_class = "insertion"
            else:
                # All native region lengths
                max_native_len = max(
                    member_spans[name][1] - member_spans[name][0]
                    for name in members
                )
                if max_native_len <= 50:
                    size_class = "point_or_short"
                else:
                    size_class = "long"

            variant_regions.append(VariantRegion(
                member_spans=member_spans,
                discriminating_positions=discriminating_positions,
                alternative_bases=alternative_bases,
                size_class=size_class,
            ))

    # Convert backbone intervals to tuples
    bb_tuples: dict[str, tuple[tuple[int, int], ...]] = {
        name: tuple(intervals) for name, intervals in backbone_intervals.items()
    }

    return Group(
        members=group.members,
        backbone_intervals=bb_tuples,
        variant_regions=tuple(variant_regions),
    )
