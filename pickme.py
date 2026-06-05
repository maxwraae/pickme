#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pickme.py  —  "Can I pick this clone?" Standalone construct verification.

WHAT YOU DO
-----------
Put this file in a folder that contains:
  * paired-end FASTQs   (*_R1*.fastq.gz  +  *_R2*.fastq.gz)
  * reference maps      (*.gb GenBank files; *.dna SnapGene optional)
...then run:

    python pickme.py          # (Windows: py pickme.py)

It finds the FASTQs and the .gb maps by itself and writes a colored Excel
plate map telling you, for every well: which construct is in it, and whether
it matches the map you ordered.

WORKS FOR ANYONE
----------------
The only thing needed is Python 3.10+. On first run the script pip-installs its
two dependencies (openpyxl + sequence-align), both shipping cross-platform
wheels. No bwa, no samtools, no conda, no compiler — it installs and runs the
same on macOS (Intel & Apple Silicon), Linux, and Windows. It's a pure-Python
reimplementation of the standard "competitive read mapping + pileup variant
calling" workflow, so it reproduces bwa/minimap2 + samtools without them.

HOW IT WORKS (two layers: identity, then quality)
-------------------------------------------------
SETUP   index every reference's 15-mers; group near-identical siblings into
        families; with a real global aligner (sequence-align) map each
        construct's DISTINCTIVE positions (what separates it from its siblings).
IDENTITY  competitive read mapping — each read is assigned to the reference it
        matches best; backbone reads tie and abstain, distinguishing reads
        decide. Most clearly-best reads = the construct. No clear winner = MIX.
        A well is only confirmed as X if its reads cover X's distinctive region
        (not just the shared vector); otherwise UNKNOWN (something else).
QUALITY  pile reads onto the called construct (seed-chaining handles indels of
        any size). At each position, a non-reference base is a REAL variant only
        if its read count exceeds sequencing error at that depth (binomial test,
        --error-rate). Real variants sort by fraction: >=70% = mutation,
        30-70% = heterogeneous (not a pure clone), <30% = minor subpopulation.
        Differences recurring across most wells of a construct = stale-map
        error, demoted; only well-specific events condemn a clone.
VERDICT  NO_DATA -> RED_UNKNOWN -> RED_MIX -> RED_BROKEN -> RED_MUT -> RED_HET
        -> YELLOW -> GREEN.

Output lands in ./pickme_results/ :  plate_map.xlsx  +  .txt  +  .json
See README.md for the full description. Run `python pickme.py --help` for options.

Author: Max Wraae.  MIT License.
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import re
import subprocess
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field, asdict
from pathlib import Path


# ─────────────────────────────────────────────────────────────────────────────
# 0.  Dependency bootstrap — make "just run it" true (openpyxl is pure-Python)
# ─────────────────────────────────────────────────────────────────────────────

def _ensure_deps() -> None:
    need = []
    for mod, pkg in (("openpyxl", "openpyxl"), ("sequence_align", "sequence-align")):
        try:
            __import__(mod)
        except ImportError:
            need.append(pkg)
    if not need:
        return
    print(f"[setup] Installing {', '.join(need)} (one-time) ...", flush=True)
    cmd = [sys.executable, "-m", "pip", "install", "--quiet", *need]
    for extra in ([], ["--user"]):
        try:
            subprocess.check_call(cmd + extra)
            print("[setup] Done.\n", flush=True)
            return
        except subprocess.CalledProcessError:
            continue
    print(f"\nERROR: could not auto-install {', '.join(need)}. Run once yourself:\n"
          f"    {sys.executable} -m pip install {' '.join(need)}\n", file=sys.stderr)
    sys.exit(1)


# ─────────────────────────────────────────────────────────────────────────────
# 1.  Data model
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class Feature:
    kind: str
    label: str
    start: int
    end: int


@dataclass
class Reference:
    name: str
    sequence: str
    features: list = field(default_factory=list)
    group: int = -1


@dataclass
class Diff:
    pos: int          # native 0-based reference position
    expected: str
    observed: str
    depth: int
    fraction: float
    kind: str         # "mutation" | "heterogeneous" | "systematic"
    features: list = field(default_factory=list)
    note: str = ""    # e.g. "matches sibling AOPE_50_5"
    ref_count: int = 0
    alt_count: int = 0


@dataclass
class WellResult:
    sample: str
    plate: str
    well: str
    total_reads: int = 0
    mapped_reads: int = 0
    called: str | None = None
    called_fraction: float = 0.0
    runner_up: str | None = None
    runner_up_fraction: float = 0.0
    breadth: float = 0.0
    mean_depth: float = 0.0
    distinctive_cov: float = 1.0  # fraction of the called construct's distinctive
    has_distinctive: bool = False # positions actually covered by reads
    member_support: int = 0      # reads that uniquely fingerprint the called clone
    mix_flag: bool = False
    n_mut: int = 0               # clone-specific homozygous mutations (disqualifying)
    n_het: int = 0               # heterogeneous positions (two substantial alleles)
    n_minor: int = 0             # real but low-level subpopulation positions
    n_systematic: int = 0        # differences shared across the plate = map error
    differences: list = field(default_factory=list)
    deletion: tuple | None = None      # (start, end, bp) suspected internal deletion
    deletion_features: list = field(default_factory=list)
    verdict: str = "NO_DATA"
    reason: str = ""


# ─────────────────────────────────────────────────────────────────────────────
# 2.  Sequence helpers
# ─────────────────────────────────────────────────────────────────────────────

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


# ─────────────────────────────────────────────────────────────────────────────
# 3.  Reference loading  (pure-Python GenBank parser — no Biopython needed)
# ─────────────────────────────────────────────────────────────────────────────

def parse_genbank(path: Path) -> Reference:
    text = path.read_text(errors="replace")
    seq_chars, feature_lines = [], []
    in_origin = in_features = False
    for line in text.splitlines():
        if line.startswith("ORIGIN"):
            in_origin, in_features = True, False
            continue
        if line.startswith("FEATURES"):
            in_features = True
            continue
        if in_origin:
            if line.startswith("//"):
                in_origin = False
                continue
            seq_chars.append(re.sub(r"[^A-Za-z]", "", line))
        elif in_features:
            if line and not line.startswith(" "):
                in_features = False
            else:
                feature_lines.append(line)
    sequence = "".join(seq_chars).upper()

    features: list[Feature] = []
    cur_kind = cur_label = cur_span = None
    loc_re = re.compile(r"(\d+)\.\.(\d+)")

    def flush():
        if cur_kind and cur_span and cur_kind != "source":
            features.append(Feature(cur_kind, cur_label or cur_kind,
                                    cur_span[0], cur_span[1]))

    for line in feature_lines:
        m = re.match(r"\s{5}(\S+)\s+(.*)", line)
        if m and not line.lstrip().startswith("/"):
            flush()
            cur_kind, cur_label = m.group(1), None
            nums = loc_re.findall(m.group(2))
            if nums:
                cur_span = (min(int(a) for a, _ in nums) - 1,
                            max(int(b) for _, b in nums))
            else:
                cur_span = None
        else:
            qm = re.search(r'/(?:label|gene|product|note)="?([^"]+)"?', line)
            if qm and cur_label is None:
                cur_label = qm.group(1).strip()
    flush()
    return Reference(name=path.stem, sequence=sequence, features=features)


def try_parse_snapgene(path: Path) -> Reference | None:
    try:
        import snapgene_reader  # type: ignore
    except Exception:
        return None
    try:
        data = snapgene_reader.parse(str(path))
        seq = (data.get("seq") or data.get("sequence") or "").upper()
        feats = [Feature(f.get("type", "misc_feature"),
                         f.get("label", f.get("name", "feature")),
                         int(f.get("start", 0)), int(f.get("end", 0)))
                 for f in (data.get("features", []) or [])]
        return Reference(path.stem, seq, feats) if seq else None
    except Exception:
        return None


def load_references(refs: list[Path]) -> list[Reference]:
    out: list[Reference] = []
    for p in sorted(refs):
        if p.suffix.lower() == ".gb":
            try:
                r = parse_genbank(p)
                if r.sequence:
                    out.append(r)
                else:
                    print(f"  ! {p.name}: no sequence, skipped")
            except Exception as e:
                print(f"  ! {p.name}: parse failed ({e}), skipped")
        elif p.suffix.lower() == ".dna":
            r = try_parse_snapgene(p)
            if r:
                out.append(r)
            else:
                print(f"  ! {p.name}: .dna needs `pip install snapgene_reader` — skipped")
    seen = {}
    for r in out:
        seen.setdefault(r.name, r)
    return list(seen.values())


# ─────────────────────────────────────────────────────────────────────────────
# 4.  k-mer index
# ─────────────────────────────────────────────────────────────────────────────

K = 15


def build_index(refs: list[Reference]):
    """kmer -> list of (ref_idx, strand, forward_strand_pos)."""
    index: dict[str, list] = defaultdict(list)
    for ri, ref in enumerate(refs):
        seq, n = ref.sequence, len(ref.sequence)
        for i in range(n - K + 1):
            index[seq[i:i + K]].append((ri, 1, i))
        rc = revcomp(seq)
        for i in range(n - K + 1):
            index[rc[i:i + K]].append((ri, -1, n - i - K))
    return index


def build_unique(index, refs):
    """kmer -> ref_idx for k-mers that occur in exactly ONE reference.
    These 'member-unique' k-mers are what distinguish near-identical siblings;
    shared backbone k-mers cannot."""
    uniq: dict[str, int] = {}
    for kmer, post in index.items():
        rset = {p[0] for p in post}
        if len(rset) == 1:
            uniq[kmer] = next(iter(rset))
    return uniq


def classify_read(read: str, index, uniq, refs):
    """One pass over a read's k-mers. Returns (group, member_idx_or_None,
    mapped_bool). `group` is which construct family the read belongs to (robust,
    from shared k-mers). `member` is the specific sibling it matches, decided
    ONLY by member-unique k-mers; None means a backbone read that can't pick a
    sibling (it still counts toward the group)."""
    gc = Counter()   # group -> hits
    mc = Counter()   # ref_idx -> unique-kmer hits
    seeds = 0
    for off in range(0, len(read) - K + 1):
        post = index.get(read[off:off + K])
        if not post:
            continue
        seeds += 1
        rset = {p[0] for p in post}
        for g in {refs[r].group for r in rset}:
            gc[g] += 1
        u = uniq.get(read[off:off + K])
        if u is not None:
            mc[u] += 1
    if seeds < 3 or not gc:
        return (None, None, False)
    # member-unique k-mers discriminate BOTH the sibling and the family, so let
    # them choose the group; fall back to shared-k-mer group only for backbone-
    # only reads that carry no unique k-mer at all.
    if mc:
        member = mc.most_common(1)[0][0]
        group = refs[member].group
    else:
        member = None
        group = gc.most_common(1)[0][0]
    return (group, member, True)


# ─────────────────────────────────────────────────────────────────────────────
# 5.  Seed-chaining pileup  (the engine)
# ─────────────────────────────────────────────────────────────────────────────

def seedchain_tally(read: str, target_idx: int, ref_len: int,
                    index, depth: list, base_counts):
    """Place `read` onto reference `target_idx` by chaining its matching seed
    blocks, and tally per-position base counts. A jump between colinear blocks
    is an indel (handled, not smeared). Returns True if the read was placed."""
    # gather seeds for this target on both strands
    fwd, rev = [], []
    L = len(read)
    for off in range(0, L - K + 1):
        for (ri, st, pos) in index.get(read[off:off + K], ()):
            if ri != target_idx:
                continue
            (fwd if st == 1 else rev).append((off, pos))
    if len(fwd) + len(rev) < 3:
        return False
    use_rev = len(rev) > len(fwd)
    seeds = rev if use_rev else fwd
    oriented = revcomp(read) if use_rev else read
    # convert to oriented-read offset + diagonal (ref_pos - oriented_off)
    pts = []
    for off, pos in seeds:
        ooff = (L - K - off) if use_rev else off
        pts.append((ooff, pos - ooff))
    pts.sort()
    # group seeds into blocks of constant diagonal, split when offsets jump >40
    bydiag = defaultdict(list)
    for ooff, diag in pts:
        bydiag[diag].append(ooff)
    blocks = []
    for diag, offs in bydiag.items():
        offs.sort()
        s = p = offs[0]
        for o in offs[1:]:
            if o - p > 40:
                blocks.append((s, p + K, diag))
                s = o
            p = o
        blocks.append((s, p + K, diag))
    # tally; assign each oriented position once, largest blocks win overlaps
    covered = bytearray(L)
    for (omin, omax, diag) in sorted(blocks, key=lambda b: -(b[1] - b[0])):
        for o in range(max(0, omin), min(L, omax)):
            if covered[o]:
                continue
            covered[o] = 1
            a = oriented[o]
            if a in "ACGT":
                rp = (o + diag) % ref_len
                depth[rp] += 1
                base_counts[rp][a] += 1
    return True


def seedchain_map(query: str, target: str):
    """Align `query` onto `target` by seed-chaining and return {target_pos:
    query_base} for every aligned position. Used to compare two REFERENCES so we
    learn exactly which positions distinguish sibling clones (SNPs and indels),
    in the target's coordinate frame."""
    tidx = defaultdict(list)
    for i in range(len(target) - K + 1):
        tidx[target[i:i + K]].append(i)
    seeds = []
    for off in range(len(query) - K + 1):
        for pos in tidx.get(query[off:off + K], ()):
            seeds.append((off, pos - off))
    if not seeds:
        return {}
    bydiag = defaultdict(list)
    for off, diag in seeds:
        bydiag[diag].append(off)
    blocks = []
    for diag, offs in bydiag.items():
        offs.sort()
        s = p = offs[0]
        for o in offs[1:]:
            if o - p > 40:
                blocks.append((s, p + K, diag))
                s = o
            p = o
        blocks.append((s, p + K, diag))
    mapping = {}
    covered = bytearray(len(query))
    for (omin, omax, diag) in sorted(blocks, key=lambda b: -(b[1] - b[0])):
        for o in range(max(0, omin), min(len(query), omax)):
            if covered[o]:
                continue
            covered[o] = 1
            tp = o + diag
            if 0 <= tp < len(target):
                mapping[tp] = query[o]
    return mapping


def build_distinctive(refs, ref_by_group):
    """Once, with a REAL global aligner: for every construct, find the positions
    that distinguish it from its same-family siblings (the insert / variant
    region). Seed-chaining can't do this — it skips exactly the differing
    positions — so we use sequence-align's Needleman-Wunsch (Hirschberg).
    Returns {ref_idx: {target_pos: {sibling_idx: base}}}."""
    from sequence_align.pairwise import hirschberg
    pair_cache = {}

    def align(qi, ti):
        key = (qi, ti)
        if key not in pair_cache:
            pair_cache[key] = hirschberg(
                list(refs[qi].sequence), list(refs[ti].sequence), gap="-",
                match_score=1, mismatch_score=-1, indel_score=-1)
        return pair_cache[key]

    out = {}
    for grp, members in ref_by_group.items():
        for ci in members:
            target = refs[ci].sequence
            disc = {}
            for m in members:
                if m == ci:
                    continue
                qa, qb = align(m, ci)        # qa=sibling, qb=target(called)
                tp = 0
                for x, y in zip(qa, qb):
                    if y != "-":             # consumes a target position
                        if x != "-" and x != target[tp] and x in "ACGT":
                            disc.setdefault(tp, {})[m] = x
                        tp += 1
            out[ci] = disc
    return out


def discriminating_map(target_idx, members, refs):
    """For the called construct `target_idx`, return {target_pos: {member_idx:
    base}} at every position where a sibling differs from it — the fingerprint
    positions used to confirm identity and detect mixes from the pileup."""
    target = refs[target_idx].sequence
    disc: dict[int, dict[int, str]] = {}
    for m in members:
        if m == target_idx:
            continue
        mp = seedchain_map(refs[m].sequence, target)
        for tp, qb in mp.items():
            if qb != target[tp] and qb in "ACGT":
                disc.setdefault(tp, {})[m] = qb
    return disc


# ─────────────────────────────────────────────────────────────────────────────
# 6.  FASTQ reading
# ─────────────────────────────────────────────────────────────────────────────

def read_fastq_seqs(path: Path, cap: int | None = None):
    op = gzip.open if path.suffix == ".gz" else open
    n = 0
    with op(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            seq = fh.readline().strip().upper()
            fh.readline()
            fh.readline()
            if seq:
                yield seq
                n += 1
                if cap and n >= cap:
                    break


# ─────────────────────────────────────────────────────────────────────────────
# 7.  Reference grouping (k-mer Jaccard -> union-find)
# ─────────────────────────────────────────────────────────────────────────────

def group_refs(refs: list[Reference]):
    sets = []
    for r in refs:
        s = set()
        for i in range(0, len(r.sequence) - K + 1, 7):
            s.add(r.sequence[i:i + K])
        sets.append(s)
    parent = list(range(len(refs)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i in range(len(refs)):
        for j in range(i + 1, len(refs)):
            if sets[i] and sets[j]:
                if len(sets[i] & sets[j]) / min(len(sets[i]), len(sets[j])) >= 0.5:
                    parent[find(j)] = find(i)
    groups = defaultdict(list)
    for i, r in enumerate(refs):
        r.group = find(i)
        groups[r.group].append(i)
    return groups


# ─────────────────────────────────────────────────────────────────────────────
# 8.  Per-well evaluation
# ─────────────────────────────────────────────────────────────────────────────

MIN_MAPPED = 100
ASSIGN_MARGIN = 5            # a read is assigned to a ref only if it beats 2nd by this
MIN_ASSIGNED = 10           # need this many discriminating reads to resolve the sibling
CONF_CLEAN = 0.80           # winner must hold this share of discriminating reads
DISTINCTIVE_MIN_POSITIONS = 5   # need this many distinctive positions to gate on them
DISTINCTIVE_COV_MIN = 0.30      # must cover this fraction of them, else "something else"
MIN_MEMBER_SUPPORT = 5        # need this many fingerprint reads to trust a call
MIX_FRACTION = 0.20
SIBLING_MIX_FRACTION = 0.15   # a 15% sibling contaminant already flags a mix
BREADTH_GREEN = 0.90
BREADTH_UNKNOWN = 0.80        # below this, the reads don't really cover any map
DIFF_MIN_DEPTH = 10           # need this depth to judge a position at all
# Quality layer: a base is a REAL variant only if its read count exceeds what
# sequencing error explains at that depth (a binomial test), not a flat %.
ERROR_RATE = 0.01             # assumed per-base sequencing error (set with --error-rate)
Z_SIGNIF = 5.0                # ~genome-wide significance for the error test
MIN_ALT_READS = 3             # absolute floor before any test
# once a base is statistically REAL, its fraction sorts it:
MUT_FRACTION = 0.70           # alt this dominant -> homozygous mutation
HET_FRACTION = 0.30           # both alleles this substantial -> heterogeneous (not pure)
THIN_DEPTH = 15
DELETION_MIN_BP = 30
SYSTEMATIC_FRACTION = 0.50    # a fixed diff in >=half a sibling's wells = map error
MINOR_YELLOW_COUNT = 3        # this many subclone-level positions = look closer
MIX_DISC_POSITIONS = 3        # sibling base present at >=this many disc positions = mix


def assign_read_competitive(read, index, refs):
    """Competitive read mapping (the standard identity method): count exact
    k-mer matches to every reference, on whichever strand the read maps. Return
    (best_ref, strict, group) where `strict` means the best reference beats the
    runner-up by a clear margin — i.e. the read spans a position that actually
    distinguishes the siblings. Backbone reads tie (strict=False) and abstain."""
    cnt = None
    for s in (read, revcomp(read)):
        c = Counter()
        for off in range(0, len(s) - K + 1):
            for (ri, _strand, _pos) in index.get(s[off:off + K], ()):
                c[ri] += 1
        if c:
            cnt = c
            break
    if not cnt:
        return None
    ranked = cnt.most_common()
    best, bestc = ranked[0]
    second = ranked[1][1] if len(ranked) > 1 else 0
    return (best, (bestc - second) >= ASSIGN_MARGIN, refs[best].group)


def _pileup_member_evidence(depth, base_counts, disc, called_idx, members, refs):
    """At every discriminating position, read the pileup and tally:
      major_for[m]  = positions where member m's base is the MAJORITY (identity)
      present_for[m]= positions where m's base is a real SECOND population (mix)
    This is the per-nucleotide agree/disagree logic that decides identity & mix
    from what the reads actually show, not from k-mer vote proportions."""
    callbase = refs[called_idx].sequence
    major_for, present_for = Counter(), Counter()
    ncov = 0
    for p, sib in disc.items():
        d = depth[p]
        if d < DIFF_MIN_DEPTH:
            continue
        ncov += 1
        c = base_counts[p]
        maj = c.most_common(1)[0][0]
        cb = callbase[p]
        if cb == maj:
            major_for[called_idx] += 1
        for m, mb in sib.items():
            if mb not in "ACGT":
                continue
            frac = c[mb] / d
            if mb == maj and mb != cb:
                major_for[m] += 1
            elif mb != cb and c[mb] >= DIFF_MIN_READS and frac >= SIBLING_MIX_FRACTION:
                present_for[m] += 1
    return ncov, major_for, present_for


def evaluate_well(res: WellResult, r1: Path, r2: Path, refs, index, uniq,
                  ref_by_group, distinct_maps, cap: int):
    reads = list(read_fastq_seqs(r1, cap)) + list(read_fastq_seqs(r2, cap))
    res.total_reads = len(reads)
    if not reads:
        res.verdict, res.reason = "NO_DATA", "no reads in file"
        return res

    # ── IDENTITY: competitive read assignment (the standard method) ──
    # Assign each read to the reference it matches best. Backbone reads tie and
    # abstain; reads that span a distinguishing position cast the deciding vote.
    # The construct with the most uniquely-best reads is the identity; if no one
    # construct dominates, the well is a mix / ambiguous.
    assigned = Counter()           # ref_idx -> reads that match it best (clearly)
    group_totals = Counter()       # family -> reads
    read_group = []                # per-read family (None if unmapped)
    mapped = 0
    for rd in reads:
        a = assign_read_competitive(rd, index, refs)
        if a is None:
            read_group.append(None)
            continue
        best, strict, grp = a
        mapped += 1
        group_totals[grp] += 1
        read_group.append(grp)
        if strict:
            assigned[best] += 1
    res.mapped_reads = mapped
    if mapped < MIN_MAPPED:
        res.verdict = "NO_DATA"
        res.reason = f"only {mapped} reads mapped (<{MIN_MAPPED})"
        return res

    total_assigned = sum(assigned.values())
    runner_idx, runner_conf = None, 0.0
    if total_assigned >= MIN_ASSIGNED:
        ranked = assigned.most_common()
        called_idx, conf = ranked[0][0], ranked[0][1] / total_assigned
        if len(ranked) > 1:
            runner_idx = ranked[1][0]
            runner_conf = ranked[1][1] / total_assigned
    else:
        # too few discriminating reads to resolve the sibling — best guess from
        # the busiest family, flagged low-confidence
        win_g = group_totals.most_common(1)[0][0]
        called_idx = max(ref_by_group[win_g], key=lambda m: assigned.get(m, 0))
        conf = 0.0
    called = refs[called_idx]
    win_group = called.group
    members = ref_by_group[win_group]
    res.called = called.name
    res.called_fraction = conf
    res.member_support = assigned.get(called_idx, 0)

    # mix / ambiguous: no single construct holds a clear majority of the
    # discriminating reads (A3-style 67/17/15 splits land here, not on a clean call)
    if total_assigned >= MIN_ASSIGNED and conf < CONF_CLEAN and runner_idx is not None:
        res.mix_flag = True
        res.runner_up = refs[runner_idx].name
        res.runner_up_fraction = runner_conf
        res.reason = (f"no single construct dominates — {res.called} "
                      f"{conf:.0%} / {res.runner_up} {runner_conf:.0%}")

    # ── pile reads onto the called construct ──
    # Pile EVERY read; seedchain_tally only places a read that has >=3 seeds
    # matching this construct, so shared backbone covers fully and truly foreign
    # reads are dropped automatically. (Restricting by best-ref family starved
    # APOE_80 wells, whose backbone reads tie-break onto APOE_50.)
    def pileup(idx):
        nn = len(refs[idx].sequence)
        dep = [0] * nn
        bc = [Counter() for _ in range(nn)]
        for rd in reads:
            seedchain_tally(rd, idx, nn, index, dep, bc)
        return dep, bc

    depth, base_counts = pileup(called_idx)
    disc = distinct_maps.get(called_idx, {})
    n = len(called.sequence)
    covered = sum(1 for d in depth if d > 0)
    res.breadth = covered / n if n else 0.0
    nz = [d for d in depth if d > 0]
    res.mean_depth = sum(nz) / len(nz) if nz else 0.0

    # ── absolute-fit check: does the well actually contain what makes THIS
    #    construct unique (its distinctive region), or just the shared vector? ──
    disc_positions = list(disc)
    res.has_distinctive = len(disc_positions) >= DISTINCTIVE_MIN_POSITIONS
    if res.has_distinctive:
        present = sum(1 for p in disc_positions if depth[p] >= 3)
        res.distinctive_cov = present / len(disc_positions)

    # ── QUALITY: per position, is each nucleotide true? ──
    # Stack the reads; at each position read the consensus. A non-reference base
    # counts only if it's a REAL variant (more reads than error explains). Then:
    #   - reference base no longer real, alt is        -> homozygous MUTATION
    #   - BOTH the reference base and the alt are real  -> HETEROGENEOUS (not pure)
    #   - only the reference is real                    -> clean (skip)
    diffs: list[Diff] = []
    for pos in range(n):
        d = depth[pos]
        if d < DIFF_MIN_DEPTH:
            continue
        exp = called.sequence[pos]
        if exp not in "ACGT":
            continue
        cc = base_counts[pos]
        ref_count = cc[exp]
        alt = max("ACGT", key=lambda b: -1 if b == exp else cc[b])
        alt_count = cc[alt]
        if not _is_real_variant(alt_count, d):
            continue                         # clean, or below the error floor
        f = alt_count / d
        if f >= MUT_FRACTION:
            kind = "mutation"                # alt dominant -> homozygous change
        elif f >= HET_FRACTION:
            kind = "heterogeneous"           # two substantial alleles -> not pure
        else:
            kind = "minor"                   # real, but a low-level subpopulation
        note = ""
        if alt in disc.get(pos, {}).values():
            note = "matches " + next(refs[m].name
                                     for m, mb in disc[pos].items() if mb == alt)
        diffs.append(Diff(
            pos=pos, expected=exp, observed=alt, depth=d, fraction=alt_count / d,
            kind=kind, features=_features_at(called, pos, pos + 1), note=note,
            ref_count=ref_count, alt_count=alt_count))
    res.differences = diffs

    # ── suspected internal deletion: a dip well below flanking coverage ──
    res.deletion = _deletion_dip(depth)
    if res.deletion:
        s, e, bp = res.deletion
        res.deletion_features = _features_at(called, s, e)

    _recount(res)
    decide_verdict(res)
    return res


def _is_real_variant(count, depth):
    """Is `count` reads of one base more than sequencing ERROR explains at this
    depth? Binomial test (normal approximation with continuity correction)
    against ERROR_RATE. This is the truth-test: it trusts 500 reads at 5% AND
    distrusts 4 reads at 13% — count and fraction together, no flat threshold."""
    if count < MIN_ALT_READS:
        return False
    mu = depth * ERROR_RATE
    sd = math.sqrt(depth * ERROR_RATE * (1.0 - ERROR_RATE))
    if sd <= 0:
        return count >= MIN_ALT_READS
    z = (count - 0.5 - mu) / sd
    return z >= Z_SIGNIF


def _recount(res: WellResult):
    # clone-specific homozygous mutations (the disqualifying kind)
    res.n_mut = sum(1 for x in res.differences
                    if x.kind == "mutation" and not x.note)
    # heterogeneous positions: two substantial alleles = not a pure clone
    res.n_het = sum(1 for x in res.differences if x.kind == "heterogeneous")
    res.n_minor = sum(1 for x in res.differences if x.kind == "minor")
    res.n_systematic = sum(1 for x in res.differences if x.kind == "systematic")


def decide_verdict(res: WellResult):
    """Turn the measurements into a PICK decision. Order = priority.
    Key principle: only WELL-SPECIFIC fixed differences (n_fixed) condemn a
    clone. Differences shared across the whole plate (n_systematic) are a map
    error, not a mutation, and never make a well red."""
    # genuinely not identifiable -> say UNKNOWN, never fake a MUT call
    if res.mapped_reads < MIN_MAPPED:
        res.verdict, res.reason = "NO_DATA", res.reason or "too few reads"
        return
    if res.breadth < BREADTH_UNKNOWN:
        res.verdict, res.called = "RED_UNKNOWN", None
        res.reason = (f"reads don't cover any reference you supplied "
                      f"(best breadth {res.breadth:.0%})")
        return
    # absolute fit: covers the shared vector but NOT this construct's distinctive
    # region -> it's something else (a construct whose map you didn't supply)
    if res.has_distinctive and res.distinctive_cov < DISTINCTIVE_COV_MIN:
        res.reason = (f"matches the shared backbone but not {res.called}'s "
                      f"distinctive region ({res.distinctive_cov:.0%}) — likely a "
                      f"construct not in your references")
        res.verdict, res.called = "RED_UNKNOWN", None
        return
    if res.mix_flag:
        res.verdict = "RED_MIX"
        return
    if res.deletion and res.deletion[2] >= 150:
        s, e, bp = res.deletion
        feat = ", ".join(res.deletion_features[:3]) or f"bp {s+1}-{e}"
        res.verdict, res.reason = "RED_BROKEN", f"{bp} bp deletion ({feat})"
        return
    if res.n_mut > 0:
        ex = next(x for x in res.differences
                  if x.kind == "mutation" and not x.note)
        res.verdict = "RED_MUT"
        res.reason = (f"{res.n_mut} clone-specific mutation(s), "
                      f"e.g. pos {ex.pos+1} {ex.expected}>{ex.observed} "
                      f"({ex.fraction:.0%} of {ex.depth}x)")
        return
    if res.n_het > 0:
        ex = next(x for x in res.differences if x.kind == "heterogeneous")
        res.verdict = "RED_HET"
        res.reason = (f"{res.n_het} heterogeneous position(s) — two alleles, not a "
                      f"pure clone, e.g. pos {ex.pos+1} "
                      f"{ex.expected}={ex.ref_count}/{ex.observed}={ex.alt_count}")
        return
    if res.breadth < BREADTH_GREEN:
        res.verdict = "YELLOW"
        res.reason = f"partial coverage (breadth {res.breadth:.0%})"
        return
    if res.mean_depth < THIN_DEPTH:
        res.verdict = "YELLOW"
        res.reason = f"thin coverage (mean {res.mean_depth:.0f}x)"
        return
    res.verdict = "GREEN"
    notes = []
    if res.n_systematic:
        notes.append(f"{res.n_systematic} systematic map mismatch(es)")
    if res.n_minor:
        notes.append(f"{res.n_minor} low-level variant(s) <30%")
    res.reason = "matches the map" + (f"  ({'; '.join(notes)} ignored)"
                                      if notes else "")


def mark_systematic(results, refs):
    """Second pass across the whole plate: a fixed difference that shows up in
    most wells of a construct FAMILY is the MAP being wrong (the .gb export is
    stale), not the clones being mutated — siblings share the same backbone, so
    they share the same map errors. Pool by family, not by individual sibling,
    or a construct with only 2-3 wells never reaches a stable threshold. Re-label
    those positions 'systematic' so they stop condemning good clones."""
    # pool per SIBLING (one coordinate frame); a difference seen in most wells of
    # that sibling is the map being wrong, not the clones mutating.
    by_ref = defaultdict(list)
    for r in results:
        if r.called and not r.mix_flag:
            by_ref[r.called].append(r)
    for ref_name, wells in by_ref.items():
        n = len(wells)
        if n < 2:
            continue
        pos_count = Counter()
        for r in wells:
            for d in r.differences:
                # any real call recurring across the sibling's wells is a map /
                # sequence-context artifact, not a per-clone event
                if d.kind in ("mutation", "heterogeneous", "minor",
                              "systematic") and not d.note:
                    pos_count[(d.pos, d.observed)] += 1
        threshold = max(2, int(round(SYSTEMATIC_FRACTION * n)))
        systematic = {k for k, c in pos_count.items() if c >= threshold}
        for r in wells:
            for d in r.differences:
                if (d.pos, d.observed) in systematic and not d.note:
                    d.kind = "systematic"
            _recount(r)
            decide_verdict(r)


def _sibling_bases(refs, members, called_idx):
    """{member_idx: {called_pos: set(bases)}} where sibling differs from called.
    Lightweight: ungapped compare on shared k-mer anchors is unreliable, so we
    do a simple positional walk via difflib-free banded-free heuristic: compare
    sequences by anchoring on identical k-mers. For robustness we just compare
    the raw sequences position-by-position up to the shorter length AND via a
    second pass on the reverse to catch trailing differences."""
    out: dict[int, dict[int, set]] = {}
    called = refs[called_idx].sequence
    n = len(called)
    for m in members:
        if m == called_idx:
            continue
        sib = refs[m].sequence
        # anchor-based alignment using shared 15-mers (cheap, gap-tolerant)
        pos_map = _anchor_map(called, sib)
        d: dict[int, set] = {}
        for cp, sp in pos_map.items():
            if 0 <= cp < n and 0 <= sp < len(sib):
                if sib[sp] != called[cp] and sib[sp] in "ACGT":
                    d.setdefault(cp, set()).add(sib[sp])
        out[m] = d
    return out


def _anchor_map(a: str, b: str):
    """Map positions of a -> positions of b using shared 15-mer anchors and
    linear interpolation between anchors. Good enough to tell whether a
    sibling carries a given base near a position (MIX vs MUT)."""
    idx_b: dict[str, int] = {}
    for i in range(len(b) - K + 1):
        idx_b.setdefault(b[i:i + K], i)
    anchors = []
    for i in range(0, len(a) - K + 1, 5):
        j = idx_b.get(a[i:i + K])
        if j is not None:
            anchors.append((i, j))
    mapping: dict[int, int] = {}
    for (ai, bj) in anchors:
        for k in range(K):
            mapping[ai + k] = bj + k
    return mapping


def _features_at(ref: Reference, start: int, end: int):
    out = []
    for f in ref.features:
        if f.start < end and f.end > start:
            lab = f"{f.label} ({f.kind})"
            if lab not in out:
                out.append(lab)
    return out


def _deletion_dip(depth):
    """Find the longest internal run of NEAR-ZERO coverage — a real dropout, not
    a Tn5 coverage valley. Requires the run to be essentially empty (<5% of
    median AND <=2 reads) and flanked by well-covered sequence. Returns
    (start,end,bp) or None."""
    n = len(depth)
    nz = [d for d in depth if d > 0]
    if not nz:
        return None
    med = sorted(nz)[len(nz) // 2]
    if med < 20:                       # too shallow overall to judge a dropout
        return None
    thresh = min(2, 0.05 * med)        # near-zero, not merely "below average"
    best = None
    i = 0
    while i < n:
        if depth[i] <= thresh:
            j = i
            while j < n and depth[j] <= thresh:
                j += 1
            flank_ok = (i > 50 and j < n - 50
                        and depth[i - 1] >= 0.3 * med and depth[j] >= 0.3 * med)
            if flank_ok and (best is None or (j - i) > best[2]):
                best = (i, j, j - i)
            i = j
        else:
            i += 1
    return best if best and best[2] >= DELETION_MIN_BP else None


# ─────────────────────────────────────────────────────────────────────────────
# 9.  Discovery
# ─────────────────────────────────────────────────────────────────────────────

WELL_RE = re.compile(r"([A-H])(\d{1,2})")
PLATE_RE = re.compile(r"(aRY\d+|plate[_-]?\w+)", re.IGNORECASE)


def discover(root: Path, fastq_dir, refs_dir):
    rs = refs_dir or root
    ref_files = [p for p in rs.rglob("*") if p.suffix.lower() in (".gb", ".dna")]
    fs = fastq_dir or root
    r1_files = sorted(p for p in fs.rglob("*")
                      if p.name.endswith((".fastq.gz", ".fastq")) and "_R1" in p.name)
    return ref_files, r1_files


def parse_sample(r1: Path):
    stem = r1.name
    for ext in (".fastq.gz", ".fastq"):
        if stem.endswith(ext):
            stem = stem[:-len(ext)]
    stem = re.sub(r"_R1.*$", "", stem)
    wm = WELL_RE.findall(stem)
    well = (wm[-1][0] + wm[-1][1]) if wm else "?"
    pm = PLATE_RE.search(stem)
    plate = pm.group(1) if pm else "plate1"
    return stem, plate, well


def find_r2(r1: Path):
    c = r1.parent / r1.name.replace("_R1", "_R2")
    return c if c.exists() else None


# ─────────────────────────────────────────────────────────────────────────────
# 10.  Rendering
# ─────────────────────────────────────────────────────────────────────────────

_ACTION = {
    "GREEN": "Use — matches the map.",
    "YELLOW": "Check — partial/thin coverage.",
    "RED_MUT": "Clone has a real mutation — discard or re-pick.",
    "RED_HET": "Not a pure clone (two alleles) — re-streak and re-pick.",
    "RED_MIX": "Two constructs present — re-streak and re-sequence.",
    "RED_BROKEN": "Large deletion — discard.",
    "RED_UNKNOWN": "Not identified — check the reference folder / contamination.",
    "NO_DATA": "Re-sequence — not enough reads.",
}
_LABEL = {
    "GREEN": "MATCH", "YELLOW": "CHECK", "RED_MUT": "RED — MUTATION",
    "RED_HET": "RED — HETEROGENEOUS", "RED_MIX": "RED — MIXED WELL",
    "RED_BROKEN": "RED — DELETION", "RED_UNKNOWN": "RED — UNKNOWN",
    "NO_DATA": "NO DATA",
}


def text_block(r: WellResult) -> str:
    L = [f"{r.plate}  {r.well}  —  {r.called or 'UNKNOWN'}  ·  {_LABEL[r.verdict]}",
         f"  sample: {r.sample}",
         f"  reads: {r.total_reads:,} total, {r.mapped_reads:,} mapped "
         f"({(r.mapped_reads/r.total_reads*100 if r.total_reads else 0):.0f}%)"]
    if r.called:
        L.append(f"  construct: {r.called} ({r.called_fraction:.0%} of mapped reads)")
    if r.runner_up:
        L.append(f"  also present: {r.runner_up} ({r.runner_up_fraction:.0%})")
    L.append(f"  coverage: {r.breadth:.0%} of plasmid, mean {r.mean_depth:.0f}x")
    mut = [d for d in r.differences if d.kind == "mutation" and not d.note]
    het = [d for d in r.differences if d.kind == "heterogeneous"]
    systematic = [d for d in r.differences if d.kind == "systematic"]
    sibm = [d for d in r.differences if d.note and d.kind != "heterogeneous"]
    if mut:
        L.append(f"  MUTATIONS ({len(mut)}) — consensus differs from the map "
                 f"(homozygous, above the error floor):")
        for d in mut[:30]:
            ft = f"  [{', '.join(d.features)}]" if d.features else ""
            L.append(f"     pos {d.pos+1:<6} {d.expected}>{d.observed}  "
                     f"{d.alt_count}/{d.depth} reads ({d.fraction:.0%}){ft}")
        if len(mut) > 30:
            L.append(f"     ... and {len(mut)-30} more")
    if het:
        L.append(f"  HETEROGENEOUS positions ({len(het)}) — two real alleles, "
                 f"this colony isn't a single clone:")
        for d in het[:30]:
            ft = f"  [{', '.join(d.features)}]" if d.features else ""
            nb = f"  ({d.note})" if d.note else ""
            L.append(f"     pos {d.pos+1:<6} {d.expected}={d.ref_count} / "
                     f"{d.observed}={d.alt_count}  of {d.depth}x{ft}{nb}")
    minor = [d for d in r.differences if d.kind == "minor"]
    if minor:
        L.append(f"  low-level variants ({len(minor)}, <30% of reads) — real but "
                 f"a minor subpopulation; not disqualifying")
    if systematic:
        L.append(f"  systematic map mismatches ({len(systematic)}) — present in "
                 f"most wells, so the .gb map is off here, NOT this clone "
                 f"(ignored for the verdict)")
    if sibm:
        L.append(f"  positions matching a sibling construct ({len(sibm)}) — "
                 f"possible minor cross-contamination")
    if r.deletion:
        s, e, bp = r.deletion
        ft = f"  [{', '.join(r.deletion_features)}]" if r.deletion_features else ""
        L.append(f"  suspected deletion: ~{bp} bp of low coverage at "
                 f"{s+1}-{e}{ft}")
    L.append(f"  VERDICT: {r.verdict}" + (f" — {r.reason}" if r.reason else ""))
    L.append(f"  action: {_ACTION[r.verdict]}")
    return "\n".join(L)


def write_text(results, path):
    Path(path).write_text(("\n" + "─" * 64 + "\n").join(text_block(r) for r in results),
                          encoding="utf-8")


def write_json(results, path):
    Path(path).write_text(json.dumps([asdict(r) for r in results], indent=2),
                          encoding="utf-8")


def write_xlsx(results, path):
    from openpyxl import Workbook
    from openpyxl.styles import PatternFill, Alignment, Font, Border, Side
    from openpyxl.utils import get_column_letter

    wb = Workbook()
    wb.remove(wb.active)
    center = Alignment(horizontal="center", vertical="center", wrap_text=True)
    fills = {"GREEN": "C8E6C9", "YELLOW": "FFF59D", "RED_MIX": "FFCC80",
             "RED_MUT": "EF9A9A", "RED_HET": "F48FB1", "RED_BROKEN": "E57373",
             "RED_UNKNOWN": "E0E0E0", "NO_DATA": "EEEEEE"}

    by_plate = defaultdict(dict)
    for r in results:
        by_plate[r.plate][r.well] = r
    best = {}
    for r in results:
        if r.verdict == "GREEN" and r.called:
            if r.called not in best or r.mean_depth > best[r.called].mean_depth:
                best[r.called] = r
    gold = Side(style="thick", color="FFB300")
    gborder = Border(left=gold, right=gold, top=gold, bottom=gold)

    for plate in sorted(by_plate):
        ws = wb.create_sheet(plate[:31])
        for col in range(1, 13):
            ws.cell(row=1, column=col + 1, value=col).alignment = center
        for ri, letter in enumerate("ABCDEFGH"):
            ws.cell(row=ri + 2, column=1, value=letter).alignment = center
        for well, r in by_plate[plate].items():
            m = WELL_RE.fullmatch(well)
            if not m:
                continue
            row = ord(m.group(1)) - ord("A") + 2
            col = int(m.group(2)) + 1
            if not (2 <= row <= 9 and 2 <= col <= 13):
                continue
            tag = {"GREEN": "MATCH", "YELLOW": "CHECK", "RED_MUT": "MUT",
                   "RED_HET": "HET", "RED_MIX": "MIX", "RED_BROKEN": "DEL",
                   "RED_UNKNOWN": "UNK", "NO_DATA": "—"}[r.verdict]
            txt = f"{r.called or '?'}\n{tag}"
            if r.called:
                txt += f"\n{r.mean_depth:.0f}x {r.breadth:.0%}"
                if r.n_mut:
                    txt += f"\n{r.n_mut} mut"
                elif r.n_het:
                    txt += f"\n{r.n_het} het"
            cell = ws.cell(row=row, column=col, value=txt)
            cell.fill = PatternFill("solid", fgColor=fills.get(r.verdict, "EEEEEE"))
            cell.alignment = center
            cell.font = Font(size=8)
            if best.get(r.called) is r:
                cell.border = gborder
        for col in range(1, 14):
            ws.column_dimensions[get_column_letter(col)].width = 11
        for row in range(1, 10):
            ws.row_dimensions[row].height = 38

    ws2 = wb.create_sheet("Summary")
    ws2.append(["plate", "well", "sample", "construct", "verdict", "reason",
                "mapped_reads", "breadth", "mean_depth", "mutations",
                "heterogeneous", "also_present"])
    for r in sorted(results, key=lambda x: (x.plate, x.well)):
        ws2.append([r.plate, r.well, r.sample, r.called or "", r.verdict, r.reason,
                    r.mapped_reads, round(r.breadth, 3), round(r.mean_depth, 1),
                    r.n_mut, r.n_het, r.runner_up or ""])
    for i, w in enumerate([10, 6, 32, 14, 22, 42, 12, 9, 10, 10, 13, 14], 1):
        ws2.column_dimensions[get_column_letter(i)].width = w

    ws3 = wb.create_sheet("Reports")
    mono = Font(name="Menlo", size=10)
    rn = 1
    for r in sorted(results, key=lambda x: (x.plate, x.well)):
        for line in text_block(r).split("\n"):
            ws3.cell(row=rn, column=1, value=line).font = mono
            rn += 1
        rn += 1
    ws3.column_dimensions["A"].width = 100
    wb.save(str(path))


# ─────────────────────────────────────────────────────────────────────────────
# 11.  Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        prog="pickme.py",
        description="pickme — 'can I pick this clone?'. Drop this in a folder of "
                    "paired FASTQs + .gb reference maps and run; get a colored "
                    "plate map of which construct is in each well and whether "
                    "it's clean. Pure Python, runs on any Mac/Windows.")
    ap.add_argument("folder", nargs="?", default=".", help="Folder to scan.")
    ap.add_argument("--fastq", help="FASTQ folder (default: auto-detect).")
    ap.add_argument("--refs", help="Reference folder (default: auto-detect).")
    ap.add_argument("--out", help="Output folder (default: <folder>/pickme_results).")
    ap.add_argument("--cap", type=int, default=20000,
                    help="Max reads per R1/R2 file (speed cap; default 20000).")
    ap.add_argument("--error-rate", type=float, default=0.01,
                    help="Assumed per-base sequencing error for the mutation test "
                         "(default 0.01 = 1%%). Lower it for high-quality data.")
    args = ap.parse_args()

    global ERROR_RATE
    ERROR_RATE = args.error_rate
    _ensure_deps()

    root = Path(args.folder).expanduser().resolve()
    if not root.exists():
        print(f"ERROR: folder not found: {root}", file=sys.stderr)
        sys.exit(1)
    fastq_dir = Path(args.fastq).expanduser().resolve() if args.fastq else None
    refs_dir = Path(args.refs).expanduser().resolve() if args.refs else None

    print(f"Scanning: {root}")
    ref_files, r1_files = discover(root, fastq_dir, refs_dir)
    print(f"  {len(ref_files)} reference file(s), {len(r1_files)} FASTQ pair(s)")
    if not ref_files:
        print("ERROR: no .gb/.dna references found.", file=sys.stderr)
        sys.exit(1)
    if not r1_files:
        print("ERROR: no *_R1*.fastq(.gz) files found.", file=sys.stderr)
        sys.exit(1)

    print("Loading references ...")
    refs = load_references(ref_files)
    print(f"  {len(refs)} reference(s)")
    print("Indexing & grouping ...")
    index = build_index(refs)
    ref_by_group = group_refs(refs)
    uniq = build_unique(index, refs)
    distinct_maps = build_distinctive(refs, ref_by_group)
    print(f"  {len(ref_by_group)} construct group(s), "
          f"{len(uniq):,} member-unique k-mers")

    out_dir = (Path(args.out).expanduser().resolve() if args.out
               else root / "pickme_results")
    out_dir.mkdir(parents=True, exist_ok=True)

    results = []
    n = len(r1_files)
    for i, r1 in enumerate(sorted(r1_files), 1):
        r2 = find_r2(r1)
        sample, plate, well = parse_sample(r1)
        res = WellResult(sample=sample, plate=plate, well=well)
        if r2 is None:
            res.verdict, res.reason = "NO_DATA", "no R2 mate"
            results.append(res)
            print(f"  [{i}/{n}] {sample:38} no R2")
            continue
        evaluate_well(res, r1, r2, refs, index, uniq, ref_by_group,
                      distinct_maps, args.cap)
        results.append(res)

    # ── plate-wide pass: demote map errors so good clones aren't condemned ──
    print("Cross-checking the plate (separating map errors from real mutations) ...")
    mark_systematic(results, refs)

    for i, res in enumerate(results, 1):
        extra = (f"  {res.n_mut}mut/{res.n_het}het/{res.n_systematic}sys"
                 if res.called else "")
        if res.runner_up:
            extra += f"  +{res.runner_up}"
        print(f"  [{i}/{len(results)}] {res.sample:38} {res.verdict:11} "
              f"{res.called or '-':12} {res.breadth:.0%} {res.mean_depth:.0f}x{extra}")

    base = out_dir / "plate_map"
    write_xlsx(results, base.with_suffix(".xlsx"))
    write_text(results, base.with_suffix(".txt"))
    write_json(results, base.with_suffix(".json"))
    print(f"\nDone.\n  {base.with_suffix('.xlsx')}\n  {base.with_suffix('.txt')}"
          f"\n  {base.with_suffix('.json')}")
    counts = Counter(r.verdict for r in results)
    print("\nSummary:")
    for v in ("GREEN", "YELLOW", "RED_MIX", "RED_MUT", "RED_HET", "RED_BROKEN",
              "RED_UNKNOWN", "NO_DATA"):
        if counts.get(v):
            print(f"  {v:12} {counts[v]}")


if __name__ == "__main__":
    main()
