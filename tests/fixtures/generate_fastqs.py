"""
Generate synthetic paired-end FASTQ fixtures for M2 tests.

Reads are generated from the rev340 reference sequence (positions 400-850),
producing valid FASTQ pairs that bwa mem can align at >= 80% mapping rate.
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

# --------------------------------------------------------------------------- #
# rev340 sequence (4000 bp, from .gb file)                                    #
# --------------------------------------------------------------------------- #
REV340_SEQ = (
    "GCCCGCCGCATGAACAGAATGGTAACGCCACCGAAGGCCACGTGTGCAAGACTGGAGGCT"
    "AAAAATTCCAACTCGTGGCCAAATCTCGAAGTTCACAAAGTAACCAGATATCCCCGTGCT"
    "CAATGTTATGTTTGCCTCTGCATTTGTCTAGATAGGTAAAGTAATGAAGTAAAGTTCCAT"
    "AATGAAAACATCGAAGTCACGGCCCGACTGCGCGAGTTTAG CGCCTGCACCAGAAGATCG"
    "CTCCGCACATCTCCTAAGCAAGCGGTGGGTGTTACCGCTGACTCCAGCCAGCCGTATCAA"
    "GGGTGCTTATGCCATACTGTAATGGCAGATTATGTCCCATAGGGCATCCCCCCCAAAGGG"
    "TTGCCCTAAAGGTCTTATTTACCTAAGTTTAAGCATTCACAACGTCCAAGTCTACGGCGG"
    "AACCGCGTACCACCAGATCATCGGCCGAAATAATTGCGGGTGGGCAGCTGACCAATGGAT"
    "GCATTATCGTTCCAGACGATTTAAGCCTGAAATAGCAGAGATTATCAACCGAGCGTCGTA"
    "GCACGAAAAGCGATCAATAGGAATCACACAACAAAGAGGGCTACCTTTCACTGAGAGTAA"
    "CCTGCTGTATAGGATAT CTGAATTTGGGCGTTGGATAACGGGACCGTAAAGTATTCCCCC"  # noqa: it has a space but we'll strip
    # ... we'll load it properly from the .gb file using BioPython
).replace(" ", "")

FIXTURES = Path(__file__).parent
REFS_DIR = FIXTURES / "refs"
FASTQS_DIR = FIXTURES / "fastqs"


def reverse_complement(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def load_rev340_sequence() -> str:
    """Load the rev340 sequence from the .gb fixture file."""
    gb_path = REFS_DIR / "rev340.gb"
    try:
        from Bio import SeqIO
        record = SeqIO.read(str(gb_path), "genbank")
        return str(record.seq).upper()
    except ImportError:
        # Fallback: parse ORIGIN section manually
        seq_lines = []
        in_origin = False
        with gb_path.open() as fh:
            for line in fh:
                if line.startswith("ORIGIN"):
                    in_origin = True
                    continue
                if in_origin:
                    if line.startswith("//"):
                        break
                    # Strip line numbers and spaces
                    seq_lines.append("".join(line.split()[1:]))
        return "".join(seq_lines).upper()


def generate_reads(seq: str, read_len: int = 150, n_pairs: int = 50) -> tuple[list[str], list[str]]:
    """
    Generate n_pairs of paired-end reads from seq.

    R1: forward read starting at positions distributed across the sequence
    R2: reverse-complement of the insert end (simulating paired-end)
    Insert size: ~300 bp
    """
    r1_reads: list[str] = []
    r2_reads: list[str] = []

    seq_len = len(seq)
    insert_size = 300

    # Distribute starting positions across the sequence to maximise coverage
    # Start well within the padded zone so reads align
    step = max(1, (seq_len - read_len - insert_size) // n_pairs)

    for i in range(n_pairs):
        start = (i * step) % (seq_len - read_len - insert_size)
        start = max(0, start)

        # R1: forward read
        r1_end = start + read_len
        if r1_end > seq_len:
            r1_end = seq_len
            start = r1_end - read_len
        r1_seq = seq[start:r1_end]

        # R2: reverse complement of the end of the insert
        r2_start = start + insert_size
        r2_end = r2_start + read_len
        if r2_end > seq_len:
            r2_end = seq_len
            r2_start = r2_end - read_len
        r2_seq = reverse_complement(seq[r2_start:r2_end])

        r1_reads.append(r1_seq)
        r2_reads.append(r2_seq)

    return r1_reads, r2_reads


def write_fastq_gz(reads: list[str], path: Path, read_name_prefix: str = "read") -> None:
    with gzip.open(path, "wt") as fh:
        for i, seq in enumerate(reads):
            quality = "I" * len(seq)  # Phred 40 — best quality
            fh.write(f"@{read_name_prefix}_{i + 1}\n")
            fh.write(f"{seq}\n")
            fh.write("+\n")
            fh.write(f"{quality}\n")


def generate_random_reads(read_len: int = 150, n_pairs: int = 50, seed: int = 42) -> tuple[list[str], list[str]]:
    """Generate n_pairs of purely random 150 bp paired-end reads (should not map to any reference)."""
    import random
    rng = random.Random(seed)
    bases = "ACGT"
    r1_reads: list[str] = []
    r2_reads: list[str] = []
    for _ in range(n_pairs):
        r1 = "".join(rng.choice(bases) for _ in range(read_len))
        r2 = "".join(rng.choice(bases) for _ in range(read_len))
        r1_reads.append(r1)
        r2_reads.append(r2)
    return r1_reads, r2_reads


def generate_mutated_reads() -> None:
    """
    Generate F2 fixture: rev340 with a novel mutation at position 508 (0-based).

    Position 508 (0-based) = position 509 (1-based):
      rev340 = G, rev341 = C, rev342 = G
    We mutate to 'T', which is not expected by any member → should yield RED_MUT.
    """
    FASTQS_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading rev340 sequence for F2 mutated fixture...")
    seq = load_rev340_sequence()

    # Verify the bases at position 508 (0-based) match expectations
    mut_pos = 508  # 0-based
    print(f"  Base at pos {mut_pos} (0-based) in rev340: '{seq[mut_pos]}'")
    assert seq[mut_pos] == "G", f"Expected G at pos {mut_pos}, got '{seq[mut_pos]}'"

    # Load rev341 and rev342 to double-check the mutation is novel
    rev341_gb = REFS_DIR / "rev341.gb"
    rev342_gb = REFS_DIR / "rev342.gb"

    def load_gb_seq(path: Path) -> str:
        try:
            from Bio import SeqIO
            record = SeqIO.read(str(path), "genbank")
            return str(record.seq).upper()
        except ImportError:
            seq_lines = []
            in_origin = False
            with path.open() as fh:
                for line in fh:
                    if line.startswith("ORIGIN"):
                        in_origin = True
                        continue
                    if in_origin:
                        if line.startswith("//"):
                            break
                        seq_lines.append("".join(line.split()[1:]))
            return "".join(seq_lines).upper()

    seq341 = load_gb_seq(rev341_gb)
    seq342 = load_gb_seq(rev342_gb)
    print(f"  Base at pos {mut_pos} in rev341: '{seq341[mut_pos]}'")
    print(f"  Base at pos {mut_pos} in rev342: '{seq342[mut_pos]}'")

    # All member bases at this position: G (rev340), C (rev341), G (rev342)
    # 'T' is not present in any member → mutation is novel
    all_bases = {seq[mut_pos], seq341[mut_pos], seq342[mut_pos]}
    novel_base = "T"
    assert novel_base not in all_bases, f"'{novel_base}' is not novel; present in {all_bases}"

    # Introduce the mutation
    mutated_seq = seq[:mut_pos] + novel_base + seq[mut_pos + 1:]
    print(f"  Mutated base at pos {mut_pos}: '{seq[mut_pos]}' → '{novel_base}'")

    print("Generating 50 paired-end read pairs from mutated rev340 sequence...")
    r1_reads, r2_reads = generate_reads(mutated_seq, read_len=150, n_pairs=50)

    r1_path = FASTQS_DIR / "F2_R1.fastq.gz"
    r2_path = FASTQS_DIR / "F2_R2.fastq.gz"

    write_fastq_gz(r1_reads, r1_path, read_name_prefix="F2_read")
    write_fastq_gz(r2_reads, r2_path, read_name_prefix="F2_read")

    print(f"  Written: {r1_path}")
    print(f"  Written: {r2_path}")


def main() -> None:
    FASTQS_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading rev340 sequence...")
    seq = load_rev340_sequence()
    print(f"  Sequence length: {len(seq)} bp")

    # The FASTA written by target.py prepends last 300 bp as padding for seqs >= 500 bp
    # So aligned reads should come from the padded sequence.
    # We generate reads from the native rev340 sequence; they will align to the padded ref.
    print("Generating 50 paired-end read pairs (150 bp, insert ~300 bp)...")
    r1_reads, r2_reads = generate_reads(seq, read_len=150, n_pairs=50)

    r1_path = FASTQS_DIR / "B4_R1.fastq.gz"
    r2_path = FASTQS_DIR / "B4_R2.fastq.gz"

    write_fastq_gz(r1_reads, r1_path, read_name_prefix="B4_read")
    write_fastq_gz(r2_reads, r2_path, read_name_prefix="B4_read")

    print(f"  Written: {r1_path}")
    print(f"  Written: {r2_path}")

    # G5: random reads that should not map to any reference → NO_DATA / RED_UNKNOWN
    print("Generating G5 random reads (50 pairs, 150 bp) ...")
    g5_r1, g5_r2 = generate_random_reads(read_len=150, n_pairs=50, seed=42)
    g5_r1_path = FASTQS_DIR / "G5_R1.fastq.gz"
    g5_r2_path = FASTQS_DIR / "G5_R2.fastq.gz"
    write_fastq_gz(g5_r1, g5_r1_path, read_name_prefix="G5_read")
    write_fastq_gz(g5_r2, g5_r2_path, read_name_prefix="G5_read")
    print(f"  Written: {g5_r1_path}")
    print(f"  Written: {g5_r2_path}")

    # F2: mutated reads for RED_MUT test
    generate_mutated_reads()

    print("Done.")


if __name__ == "__main__":
    main()
