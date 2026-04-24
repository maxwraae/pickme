# pickme

Tn5 construct verification for 96-well plates — drop a run in, get a colored plate map out.

## What it does

For each well on a 96-well plate it answers two questions:

1. **What construct is in this well?** — competitive BWA alignment against every reference you ordered, then winning-group + member selection using variant regions.
2. **Is it what you expected?** — per-position pileup at signature positions, backbone coverage check, optional whole-insert mutation scan.

Output is a single `.xlsx` with three sheets (plate grid, best-pick-per-construct table, per-well narrative) plus matching `.txt` and `.json` reports. Cells are colored by verdict so you can hand the sheet to anyone at the bench without explanation.

## Install

```bash
git clone https://github.com/maxwraae/pickme
cd pickme
pip install -e .
```

You also need five external binaries on `PATH`. On macOS:

```bash
brew install bwa samtools minimap2 mafft mosdepth
```

Or via conda (see `env.yml`):

```bash
conda env create -f env.yml
conda activate tn5verify
pip install -e .
```

## Use it

The ritual is: drop a run into `input/`, then run `pickme`.

1. Make a folder under `input/` named after your run (e.g. `input/260422_tn5/`).
2. Inside it, make two subfolders: `refs/` and `fastq/`.
3. Put your reference maps (`.gb` or `.dna`) in `refs/`.
4. Put the paired-end FASTQs (`*_R1*.fastq.gz` / `*_R2*.fastq.gz`) in `fastq/`.
5. Run `pickme` (TUI) or `pickme run <run_name>` (headless).

Results land at `output/<run_name>/plate_map.{xlsx,txt,json}`. BAMs go to `output/<run_name>/bam/`.

**Well detection:** the last `[A-H][0-9]{1,2}` match in each filename's stem is the well ID. So `construct_A5_S197_R1_001.fastq.gz` → well `A5`.

### Headless flags

```bash
pickme run <run_name> --scan-mode whole        # scan the whole plasmid
pickme run <run_name> --scan-mode insert --insert-region 500-3500
pickme run <run_name> --threads 8
```

## Directory layout

```
pickme/
├── input/                          # you drop runs here
│   └── <run_name>/
│       ├── refs/    (*.gb, *.dna)
│       └── fastq/   (*_R1*.fastq.gz, *_R2*.fastq.gz)
├── output/                         # results land here
│   └── <run_name>/
│       ├── plate_map.xlsx
│       ├── plate_map.txt
│       ├── plate_map.json
│       └── bam/<well_id>/well.dedup.bam
├── pickme/                         # TUI + CLI
├── tn5verify/                      # pipeline package
├── scripts/identify.py             # legacy CLI (full path control)
└── tests/
```

## Reading the verdicts

| Verdict | Meaning | Action |
|---------|---------|--------|
| `GREEN` (CLEAN) | Identity confirmed, all positions clean, backbone intact | Use — ready to transfect |
| `YELLOW` | Low coverage, partial coverage, or minor backbone warning | Re-sequence deeper, or use with caution |
| `RED_MIX` | Two constructs detected in the same well | Re-streak and re-sequence |
| `RED_MUT` | A mutation was found that matches no known reference | Discard — synthesis or cloning error |
| `RED_BROKEN` | Chunk of the plasmid is missing (often a selection marker) | Discard |
| `RED_UNKNOWN` | Reads don't match any reference you ordered | Check references folder or contamination |
| `NO_DATA` | Fewer than 100 reads mapped | Re-sequence — insufficient reads |

YELLOW has a sub-reason in the report title — `(low confidence)`, `(partial coverage)`, or `(backbone warning window)`.

## For developers / AI agents

See [`CLAUDE.md`](CLAUDE.md) for architecture, tuning knobs, and per-module orientation. Run the test suite with `pytest tests/ -q`.

## License

MIT — see [`LICENSE`](LICENSE).
