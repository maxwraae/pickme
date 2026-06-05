# pickme

**Verify which plasmid clone is in each well of a sequenced plate — "can I pick this one?" — with nothing but Python.**

You cloned something, picked colonies into a 96-well plate, miniprepped, tagmented (Tn5) and sequenced. `pickme` reads the FASTQs, compares each well to the reference maps you ordered, and hands you a colored plate map that answers, per well:

- **is it the construct I expected** (clean match),
- **is it a mix** of two constructs (contaminated colony),
- **does it carry a real mutation**, or
- **can't I match it to anything you gave me** (unknown).

No `bwa`, no `samtools`, no `conda`, no compiler. It installs and runs the same on **macOS (Intel & Apple Silicon)** and **Windows**. The only thing a collaborator needs is Python 3.10+.

---

## Quick start

```bash
python pickme.py /path/to/your/run_folder
```

That's it. On the first run it `pip install`s its two pure-wheel dependencies (`openpyxl`, `sequence-align`) automatically, then writes results into `run_folder/pickme_results/`.

The folder just needs to contain:

```
run_folder/
├── *_R1*.fastq.gz   +   *_R2*.fastq.gz     # paired-end reads, one pair per well
└── *.gb                                     # GenBank maps of the constructs you ordered
```

`pickme` finds them itself (references and FASTQs may sit in sub-folders). From each FASTQ filename it reads the **plate** (an `aRY####`-style token) and the **well** (the last `A–H` + number, e.g. `..._A5_S101_R1_001.fastq.gz` → well `A5`).

> **Give it every map you put on the plate.** Identity is decided by competition between the references you supply. A well whose true plasmid you *didn't* give a map for is reported as `UNKNOWN` rather than force-fit to the nearest map — but the more complete your reference set, the sharper every call.

---

## What you get

Everything lands in `pickme_results/`:

| file | what it is |
|---|---|
| `plate_map.xlsx` | one **colored 96-well grid per plate** (green = clean match, red = problem), a **Summary** table, and a **Reports** sheet with a plain-English paragraph per well |
| `plate_map.txt`  | the same per-well reports as plain text |
| `plate_map.json` | every number, for your own scripts |

Each well gets one of these verdicts:

| verdict | colour | meaning | what to do |
|---|---|---|---|
| **GREEN** | green | clean clone, matches the map | **pick it** |
| **RED_MUT** | red | a real, well-specific mutation | discard / re-pick |
| **RED_HET** | pink | two substantial alleles — not a pure clone | re-streak and re-pick |
| **RED_MIX** | orange | two different constructs in one well | re-streak and re-sequence |
| **RED_BROKEN** | red | a large deletion (≥150 bp) | discard |
| **RED_UNKNOWN** | grey | reads don't match any map you supplied | check references / contamination |
| **YELLOW** | yellow | partial or thin coverage — can't be sure | re-sequence deeper |
| **NO_DATA** | grey | too few reads (empty well) | re-sequence |

A green well's report still lists, for transparency, any **systematic map mismatches** (positions where your `.gb` map disagrees with *every* clone — i.e. the map is stale, not the clone) and any **low-level variants** (<30%, real but a minor subpopulation). Neither condemns the clone.

---

## How it works (exactly)

`pickme` is two layers: **identity** (which construct is this?) then **quality** (is it true?). It's a pure-Python reimplementation of the standard "competitive read mapping + pileup variant calling" workflow (`bwa`/`minimap2` + `samtools`), reproduced base-for-base so it can run anywhere.

### Setup (once)
1. **Read the references.** Each `.gb` → a DNA string + its annotated features. (A small built-in GenBank parser — no Biopython.)
2. **Index k-mers.** Every 15-mer of every reference, both strands → a lookup `kmer → (reference, strand, position)`.
3. **Group siblings.** References sharing ≥50% of their 15-mers are one *family* (e.g. near-identical clone variants). This shrinks the search.
4. **Map the distinctive regions.** With a real global aligner (`sequence-align`, Needleman–Wunsch), find, for each construct, the positions that distinguish it from its siblings — its "fingerprint" / insert region.

### Per well — identity
5. **Competitive read mapping.** Each read is assigned to the reference it matches best (most exact 15-mers). Reads on shared backbone match every sibling equally and **abstain**; a read that spans a distinguishing position matches its true construct and casts the deciding vote.
6. **Call it.** The construct with the most clearly-best reads is the identity. If no construct wins a clear majority of the deciding reads → **mix / ambiguous** (`RED_MIX`).
7. **Absolute-fit check.** A well is only confirmed as construct *X* if its reads actually cover *X*'s **distinctive region** — not just the shared vector. A well that covers the backbone but 0% of *X*'s distinctive positions is **something else** (`RED_UNKNOWN`), not *X*. This is what stops an off-target plasmid from being force-fit to the nearest map.

### Per well — quality
8. **Pile up.** Lay every read onto the called construct (seed-chaining, which handles indels of any size by letting a gap between matched blocks *be* the indel) → at each position, a count of A/C/G/T and a depth.
9. **Judge each nucleotide.** At each position with depth ≥10, read the consensus and test every non-reference base against **sequencing noise**: a base is a *real* variant only if its read count exceeds what the error rate (default 1%, `--error-rate`) explains at that depth — a binomial test, so 500 reads at 5% is real but 4 reads at 13% is noise. A real variant is then sorted by fraction:
   - **≥70%** of reads → homozygous **mutation**,
   - **30–70%** → **heterogeneous** (two real alleles, not a pure clone),
   - **<30%** → a real but minor **subpopulation** (reported, not disqualifying).
10. **Separate map errors from mutations.** Across the whole plate, any difference that recurs in most wells of a construct is the **map** being wrong (siblings share a backbone, so they share its export errors) — it's demoted and never condemns a clone. Only *well-specific* mutations / heterozygosity make a well red.

### Roll-up
11. **Verdict** (first match wins): `NO_DATA` → `RED_UNKNOWN` → `RED_MIX` → `RED_BROKEN` → `RED_MUT` → `RED_HET` → `YELLOW` → `GREEN`.

---

## Options

```
python pickme.py [folder]
    --refs DIR          reference folder (default: auto-detect under folder)
    --fastq DIR         FASTQ folder   (default: auto-detect under folder)
    --out DIR           output folder  (default: folder/pickme_results)
    --cap N             max reads per R1/R2 file (speed cap; default 20000)
    --error-rate F      assumed per-base sequencing error for the mutation
                        test (default 0.01 = 1%); lower it for high-Q data
```

## Requirements

- **Python 3.10+** — the only hard requirement.
- Auto-installed on first run: [`openpyxl`](https://pypi.org/project/openpyxl/) (pure Python) and [`sequence-align`](https://pypi.org/project/sequence-align/) (ships pre-built wheels for macOS Intel/ARM, Linux x86/ARM, and Windows 32/64). No build tools needed.
- A 96-well plate of ~5–10k reads/well runs in a couple of minutes (it's pure Python by design — that's the trade for "runs anywhere").

## License

MIT — see [`LICENSE`](LICENSE).
