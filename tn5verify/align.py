from __future__ import annotations

import re
import subprocess
from dataclasses import dataclass
from pathlib import Path


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


@dataclass
class AlignmentStats:
    total_reads: int
    mapped_reads: int
    bam_path: Path


def run_bwa(
    r1: Path,
    r2: Path,
    target_fa: Path,
    out_dir: Path,
    threads: int = 4,
) -> AlignmentStats:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    bwa = _find_tool("bwa")
    samtools = _find_tool("samtools")

    namesort_bam = out_dir / "well.namesort.bam"
    fixmate_bam = out_dir / "well.fixmate.bam"
    well_bam = out_dir / "well.bam"
    dedup_bam = out_dir / "well.dedup.bam"

    # Step 1: bwa mem | samtools sort -n (name sort for fixmate)
    # Quote all paths to handle spaces in directory names
    pipeline_cmd = (
        f'"{bwa}" mem -M -t {threads} "{target_fa}" "{r1}" "{r2}" '
        f'| "{samtools}" sort -n -@ 2 -o "{namesort_bam}" -'
    )
    result = subprocess.run(
        pipeline_cmd,
        shell=True,
        capture_output=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"bwa mem | samtools sort failed: {result.stderr.decode()}"
        )

    # Step 2: samtools fixmate -m (adds ms score tag required by markdup)
    result = subprocess.run(
        [samtools, "fixmate", "-m", str(namesort_bam), str(fixmate_bam)],
        capture_output=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"samtools fixmate failed: {result.stderr.decode()}"
        )

    # Step 3: samtools sort (coordinate sort for markdup)
    result = subprocess.run(
        [samtools, "sort", "-@", "2", "-o", str(well_bam), str(fixmate_bam)],
        capture_output=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"samtools sort (coord) failed: {result.stderr.decode()}"
        )

    # Step 4: samtools index well.bam
    result = subprocess.run(
        [samtools, "index", str(well_bam)],
        capture_output=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"samtools index (well.bam) failed: {result.stderr.decode()}"
        )

    # Step 5: samtools markdup -r well.bam well.dedup.bam
    result = subprocess.run(
        [samtools, "markdup", "-r", str(well_bam), str(dedup_bam)],
        capture_output=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"samtools markdup failed: {result.stderr.decode()}"
        )

    # Step 6: samtools index well.dedup.bam
    result = subprocess.run(
        [samtools, "index", str(dedup_bam)],
        capture_output=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"samtools index (well.dedup.bam) failed: {result.stderr.decode()}"
        )

    # Step 7: flagstat to get read counts
    result = subprocess.run(
        [samtools, "flagstat", str(dedup_bam)],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"samtools flagstat failed: {result.stderr}"
        )

    total_reads = 0
    mapped_reads = 0
    for line in result.stdout.splitlines():
        if "+ 0 in total" in line:
            m = re.match(r"(\d+)", line)
            if m:
                total_reads = int(m.group(1))
        elif "+ 0 mapped" in line:
            m = re.match(r"(\d+)", line)
            if m:
                mapped_reads = int(m.group(1))

    return AlignmentStats(
        total_reads=total_reads,
        mapped_reads=mapped_reads,
        bam_path=dedup_bam,
    )
