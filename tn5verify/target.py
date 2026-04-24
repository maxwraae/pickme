from __future__ import annotations

import subprocess
from pathlib import Path

from .types import Group, Reference


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


def build_multi_fasta(
    groups: list[Group],
    refs: list[Reference],
    out_dir: Path,
) -> Path:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    fa_path = out_dir / "combined_padded.fa"

    with fa_path.open("w") as fh:
        for ref in refs:
            seq = ref.sequence
            orig_len = len(seq)
            if orig_len >= 500:
                padded = seq[-300:] + seq
            else:
                padded = seq
            fh.write(f">{ref.name}  original_length={orig_len}\n")
            # Write in 60-char lines
            for i in range(0, len(padded), 60):
                fh.write(padded[i:i + 60] + "\n")

    # bwa index
    bwa = _find_tool("bwa")
    result = subprocess.run(
        [bwa, "index", str(fa_path)],
        capture_output=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"bwa index failed: {result.stderr.decode()}")

    # samtools faidx
    samtools = _find_tool("samtools")
    result = subprocess.run(
        [samtools, "faidx", str(fa_path)],
        capture_output=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"samtools faidx failed: {result.stderr.decode()}")

    return fa_path
