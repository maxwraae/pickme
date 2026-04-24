# pickme — developer/agent reference (Tn5 construct verification)

> User-facing docs are in `README.md`. This file is for contributors and AI agents working on the codebase.

Pipeline that takes Tn5 sequencing FASTQs from a 96-well plate and tells you,
per well: **what construct is in it, and is it what you expected?**

Design spec: see `README.md` for quickstart; module docstrings explain algorithm details.

---

## What it does

For each well it answers two questions:

1. **What construct is in this well?** — competitive BWA alignment against
   every reference, winning-group + member selection via variant regions.
2. **Is it what you expected?** — per-position pileup at signature positions,
   backbone coverage check, optional whole-insert mutation scan.

Output is a single `.xlsx` with three sheets — a colored 96-well plate grid,
a "best pick per construct" table, and a per-well analysis block you can
hand to anyone at the bench without explanation.

---

## Dependencies

External tools (install via conda, brew, etc.):

- `bwa` ≥ 0.7.17
- `samtools` ≥ 1.17
- `minimap2` ≥ 2.26
- `mafft` ≥ 7.5
- `mosdepth` ≥ 0.3
- Python ≥ 3.11 with `pysam`, `biopython`, `snapgene-reader`, `openpyxl`, `pandas`, `numpy`

See `env.yml` and `pyproject.toml` for pinned versions.

---

## Input layout

```
refs_folder/
    construct_a.gb          # GenBank
    construct_b.gb
    construct_c.dna         # SnapGene also supported
    ...

fastq_folder/
    plate_name_A1_S1_R1_001.fastq.gz
    plate_name_A1_S1_R2_001.fastq.gz
    plate_name_A2_S2_R1_001.fastq.gz
    ...
```

**Well detection**: the last `[A-H][0-9]{1,2}` match in each filename's stem
is used as the well ID. So `construct_A5_S197_R1_001.fastq.gz`
becomes well `A5`. If your naming doesn't fit, pass a `--well-map` CSV with
columns `well_id,fastq_r1,fastq_r2`.

---

## Running it

### The simplest case (preferred)

Drop your run into `~/Documents/pickme/input/<run_name>/{refs,fastq}/`, then:

```bash
pickme                    # TUI — select the run and hit enter
pickme run <run_name>     # headless equivalent
```

Reports are written to `~/Documents/pickme/output/<run_name>/plate_map.{xlsx,txt,json}`
and BAMs to `output/<run_name>/bam/`.

### Legacy CLI (full path control)

```bash
python3 scripts/identify.py \
    --fastq "path/to/fastq_folder" \
    --refs "path/to/refs_folder" \
    --output plate_map.xlsx
```

Outputs land in a new `run_<fastq_folder_name>/` directory next to the fastq
folder. The `.xlsx`, `.txt`, and `.json` reports are written at `--output`.

### The three scan modes (how deep to check for mutations)

Same noise-floor algorithm (≥3 reads AND ≥1 % to count a disagreement as
real), different scan range.

| Mode | What it checks | When to use | Cost on a 96-well plate |
|------|----------------|-------------|-------------------------|
| `quick` *(default)* | Only the ~20 signature positions that distinguish siblings | Fast identity/mix/mut check at the variant sites only | Baseline (~60 s) |
| `insert` | A user-specified region `START-END` | You know where your insert is and only care about that range | +2 s |
| `whole` | The entire called reference | You want to catch any mutation anywhere in the plasmid | +5 s |

```bash
# Default — only signature positions
pickme run <run_name>

# Scan the whole plasmid for mutations (recommended — cost is negligible)
pickme run <run_name> --scan-mode whole

# Scan only a specific insert region (0-based, half-open)
pickme run <run_name> --scan-mode insert --insert-region 500-3500

# Legacy CLI (full path control) — same flags
python3 scripts/identify.py --fastq ... --refs ... --output plate.xlsx --scan-mode whole
```

### All CLI flags (legacy `scripts/identify.py`)

- `--fastq DIR` — folder with paired-end FASTQs (required)
- `--refs DIR` — folder with `.gb` / `.dna` reference maps (required)
- `--output PATH` — output `.xlsx` path; `.txt` and `.json` written alongside (required)
- `--well-map CSV` — override well-ID inference; columns `well_id,fastq_r1,fastq_r2`
- `--threads N` — BWA threads per well (default: half of CPU count)
- `--run-dir DIR` — BAM output dir (default: `run_<fastq_folder_name>/`)
- `--scan-mode {quick,insert,whole}` — mutation scan range (default: `quick`)
- `--insert-region START-END` — native coords, required when `--scan-mode=insert`

---

## Output

### `plate_map.xlsx`

Three sheets:

1. **Plate Map** — 96-well colored grid. Cell shows
   `construct_name / 62× / 20/20` (coverage · confirmed positions).
   - Green = CLEAN
   - Yellow = YELLOW (with sub-reason in the Well Reports sheet)
   - Red = MIX, MUT, BROKEN, or UNKNOWN
   - Grey = NO_DATA
   - Intensity scales with coverage; gold border marks the best pick per construct
2. **Picks** — table mapping construct → best well
3. **Well Reports** — per-well narrative block. For each well:
   - "Why we called this clone X:" header with read counts
   - Per-position pileup with expected base first, A/C/G/T counts
   - Disagreeing-reads summary
   - YELLOW wells get an explanation of what's uncertain
   - Backbone summary (intact / bp deletion with feature overlap)
   - Insert mutation scan results (if `--scan-mode=insert` or `whole`)
   - Verdict, reason, suggested action

### `plate_map.txt`

Same per-well blocks as the "Well Reports" sheet, dumped as plain text.
Handy for diff / grep.

### `plate_map.json`

Full `WellResult` data — use for custom downstream analysis or re-rendering.

### `run_<folder>/<well>/`

Per-well BAMs (`well.bam`, `well.dedup.bam`, indexes). Keep if you want to
re-pileup or inspect in IGV.

---

## How to read the verdicts

| Verdict | Meaning | Action |
|---------|---------|--------|
| `GREEN` (CLEAN) | Identity confirmed, all positions clean, backbone intact | Use — ready to transfect |
| `YELLOW` | Low coverage, partial coverage, or minor backbone warning | Re-sequence deeper, or use with caution |
| `RED_MIX` | Two constructs detected in the same well | Re-streak and re-sequence |
| `RED_MUT` | A mutation was found that matches no known reference | Discard — synthesis or cloning error |
| `RED_BROKEN` | Chunk of the plasmid is missing (often a selection marker) | Discard |
| `RED_UNKNOWN` | Reads don't match any reference you ordered | Check references folder or contamination |
| `NO_DATA` | Fewer than 100 reads mapped | Re-sequence — insufficient reads |

**YELLOW has a sub-reason** shown in parentheses in the report title:
- `(low confidence — not enough reads)` — variant positions have thin coverage
- `(partial coverage in variant region)` — some variant positions uncovered
- `(backbone warning window)` — coverage dropped in ≥2 backbone windows

---

## Tuning knobs (edit in code)

The noise-floor and trust thresholds match clinical-genotyping defaults.
If you need to change them, they live in:

- `tn5verify/identify.py` — winning-group threshold (`< 100` reads → NO_DATA),
  coverage-span threshold (`< 30 %` → RED_UNKNOWN), MIX detection (≥20 % each
  member, region_total ≥ 5).
- `tn5verify/integrity.py` — per-position noise floor (`≥3` reads AND `≥1 %`),
  trust tier (`≥10` reads → CLEAN; 3–9 → THIN; `<3` → UNCOVERED).
- `tn5verify/integrity.py::_build_backbone_windows` — backbone window size
  (100 bp), intact threshold (≥20 % of median), broken threshold (<5 % of
  median with consecutive neighbor).
- `tn5verify/insert_scan.py` — insert-scan minimum coverage to call a
  mutation (`_MIN_COVERAGE = 10`).

---

## When things go wrong

### "No wells discovered"

Filename pattern doesn't match `[A-H][0-9]{1,2}`. Use `--well-map`.

### "Could not find <tool> in PATH"

External binary not installed or not in `PATH`. The code checks
`/opt/homebrew/bin` and `/usr/local/bin` as fallbacks. Install via conda/brew.

### All wells are RED_UNKNOWN

Either contamination (reads don't match your refs) or your refs folder is
wrong. Spot-check: does `run_<folder>/<well>/well.dedup.bam` have reads
mapped to any contig? (`samtools flagstat`)

### Whole-plasmid scan shows many flagged positions on every well

Likely systematic — Tn5 tagmentation has edge artifacts in the first/last
~10 bp of each fragment, and some reference-specific low-complexity regions
can cause persistent miscalls. If flagged positions cluster at the same
coordinates across every well, treat as noise and investigate (trim reads,
mask reference, or raise `_MIN_COVERAGE`).

### BWA mapq = 0 reads at variant positions

Expected — competitive mapping tie-breaks reads that match multiple siblings
equally. The pileup code accepts mapq=0 reads on purpose (see
`integrity._pileup_position` docstring). Base-quality ≥20 filter is still
enforced to drop sequencing errors.

---

## Tests

```bash
python3 -m pytest tests/ -q
```

Six fixtures cover every verdict path (GREEN, YELLOW, RED_MUT, RED_MIX,
RED_BROKEN, RED_UNKNOWN, NO_DATA, singleton group). Regenerate fixtures with
`python3 tests/fixtures/generate_fastqs.py`.

---

## File map

```
pickme/
  __main__.py                   — `pickme` and `pickme run <name>` entrypoint
  app.py                        — Textual TUI
  runs.py                       — scan input/, status from output/<name>/plate_map.json
  runner.py                     — shared headless pipeline runner (used by TUI + CLI)
scripts/identify.py             — legacy CLI entry point (full path control)
tn5verify/
  refs.py                       — load .gb / .dna
  grouping.py                   — pairwise minimap2 → union-find groups
  characterize.py               — MAFFT → backbone + variant regions
  target.py                     — build padded multi-FASTA + bwa index
  align.py                      — bwa mem → fixmate → markdup per well
  identify.py                   — winning group + member selection + MIX
  integrity.py                  — pileup at variant positions, backbone windows
  features.py                   — .gb feature overlap for broken windows
  insert_scan.py                — whole-plasmid / insert-region mutation scan
  render/
    plate.py                    — Excel plate map, picks, well reports sheet
    well_report.py              — per-well narrative block formatter
    json_dump.py                — raw WellResult JSON
  types.py                      — dataclasses (Reference, Group, VariantRegion, WellResult, ...)
tests/
  test_m1.py ... test_m6.py     — per-milestone fixtures
  fixtures/                     — refs + simulated FASTQs for each verdict path
```
