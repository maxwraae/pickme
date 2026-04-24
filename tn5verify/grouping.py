from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

from .types import Group, Reference


def _find_tool(name: str) -> str:
    """Return full path for a tool, checking common homebrew locations."""
    import shutil
    found = shutil.which(name)
    if found:
        return found
    for prefix in ["/opt/homebrew/bin", "/usr/local/bin"]:
        candidate = Path(prefix) / name
        if candidate.exists():
            return str(candidate)
    raise FileNotFoundError(f"Could not find {name} in PATH or common locations")


def build_groups(refs: list[Reference]) -> list[Group]:
    minimap2 = _find_tool("minimap2")

    same_group_pairs: list[tuple[str, str]] = []

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)

        # Write individual FASTA files
        fasta_paths: dict[str, Path] = {}
        for ref in refs:
            fa = tmp / f"{ref.name}.fa"
            fa.write_text(f">{ref.name}\n{ref.sequence}\n")
            fasta_paths[ref.name] = fa

        ref_by_name = {r.name: r for r in refs}

        # Compare every pair
        for i, a in enumerate(refs):
            for b in refs[i + 1:]:
                if a.name >= b.name:
                    name_a, name_b = b.name, a.name
                else:
                    name_a, name_b = a.name, b.name

                fa_a = fasta_paths[name_a]
                fa_b = fasta_paths[name_b]

                try:
                    result = subprocess.run(
                        [minimap2, "-x", "asm5", "-c", "--cs", str(fa_a), str(fa_b)],
                        capture_output=True,
                        text=True,
                        check=True,
                    )
                except subprocess.CalledProcessError as e:
                    raise RuntimeError(
                        f"minimap2 failed for {name_a} vs {name_b}: {e.stderr}"
                    ) from e

                stdout = result.stdout.strip()
                if not stdout:
                    continue

                best_mlen = 0
                for line in stdout.splitlines():
                    cols = line.split("\t")
                    if len(cols) < 11:
                        continue
                    try:
                        nmatch = int(cols[9])
                        alen = int(cols[10])
                        identity = nmatch / alen if alen > 0 else 0
                        if identity >= 0.98:
                            if nmatch > best_mlen:
                                best_mlen = nmatch
                    except (ValueError, IndexError):
                        continue

                min_len = min(
                    len(ref_by_name[name_a].sequence),
                    len(ref_by_name[name_b].sequence),
                )
                if best_mlen >= 0.5 * min_len:
                    same_group_pairs.append((name_a, name_b))

    # Union-Find
    parent: dict[str, str] = {r.name: r.name for r in refs}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x: str, y: str) -> None:
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[ry] = rx

    for a, b in same_group_pairs:
        union(a, b)

    # Build groups
    groups_dict: dict[str, list[str]] = {}
    for ref in refs:
        root = find(ref.name)
        groups_dict.setdefault(root, []).append(ref.name)

    groups: list[Group] = []
    for members_list in groups_dict.values():
        groups.append(Group(
            members=tuple(sorted(members_list)),
            backbone_intervals={},
            variant_regions=(),
        ))

    return groups
