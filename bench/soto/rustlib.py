#!/usr/bin/env python3
"""Canonical primitives for the bench scripts — ONE implementation of each rule that has been got wrong.

WHY THIS EXISTS. bench/soto holds ~70 scripts and ~9,700 lines. 19 of them independently parse `de:f:`,
11 independently re-implement the coverage denominator, 37 independently load copies.tsv. Every
re-implementation is a chance to drift, and several have:

  - four chains_*.py files were written with four DIFFERENT sort keys, so `chain_id` was not comparable
    across genes (caught only by an adversarial verifier);
  - the SAM/GFF off-by-one has been re-proven from scratch three times, once after it had already produced
    a reported "0% junction recall / 100% non-canonical";
  - best-overlap was written as first-overlap at least three times, each time yielding a plausible wrong
    number (NOTCH2NLB, AMY1, NPIPB9);
  - the effective identity floor is 0.60, not the 0.80 in RefineParams::default, because the sensitive tier
    is on by default and E_r is the UNION — this had to be rediscovered by reading the Rust.

Import from here instead of re-deriving. If a rule changes in the Rust, change it HERE and every script
follows. Every function below cites the Rust it mirrors.

    from rustlib import paf_pairs, genes_on, load_copies, chains_in, best_overlap

Stdlib only. Disk cache under $RUSTLE_BENCH_CACHE (default ~/winloci_scratch/benchcache), keyed on the
input file's mtime+size, so a stale cache cannot survive an input change.
"""
from __future__ import annotations

import gzip
import hashlib
import os
import pickle
import re
import subprocess
import sys
from collections import defaultdict, deque

GFF = os.environ.get("RUSTLE_GFF", "/mnt/c/Users/jfris/Desktop/Reference/chm13v2.0_RefSeq_full.gff.gz")
FASTA = os.environ.get("RUSTLE_FASTA", "/mnt/c/Users/jfris/Desktop/Reference/chm13v2.0.fa")
CACHE = os.environ.get("RUSTLE_BENCH_CACHE", os.path.expanduser("~/winloci_scratch/benchcache"))

# ---------------------------------------------------------------- cache
def _key(path: str, *parts) -> str:
    st = os.stat(path)
    h = hashlib.sha1(f"{path}|{st.st_mtime_ns}|{st.st_size}|{parts}".encode()).hexdigest()[:16]
    return os.path.join(CACHE, f"{os.path.basename(path)}.{h}.pkl")

def _cached(path: str, parts, build):
    """Cache `build()` keyed on the INPUT's mtime+size, so editing the input invalidates it.

    ⚠⚠ WHAT THE KEY DOES AND DOES NOT COVER (audited 2026-08-11).
      COVERED explicitly, in `parts`: the coverage FORMULA ("covform=axis" — bumped when M1 landed),
        the denominator, drop_self, and BOTH floors. Change any of those and the key changes.
      COVERED ONLY THROUGH THE PATH: the ALIGNER TIER. A PAF carries no record of the flags that made
        it, so the key can only see it as (path, mtime_ns, size). That is weaker than it looks:
        measured on this box, st_mtime_ns granularity on the ext4 scratch is ~3.6 ms and two writes
        2 ms apart returned the IDENTICAL mtime_ns. Two same-size PAFs written to the same path in
        the same tick would collide.
      ⟹ THE RULE FOR CALLERS: a PAF produced at a different tier MUST get a different FILENAME.
        `ava()` output for the shipped E_r tier is named `*.er.paf` throughout bench/crossspecies for
        exactly this reason — the panel-tier `*.paf` files are still on disk, and reusing the name
        would have let `if not os.path.exists(paf)` serve panel alignments to the shipped rule.
        Do not rely on mtime to save you here; rely on the name.
    """
    try:
        k = _key(path, *parts)
    except OSError:
        return build()
    if os.path.exists(k):
        try:
            with open(k, "rb") as fh:
                return pickle.load(fh)
        except Exception:
            pass  # corrupt cache is not an error, just rebuild
    v = build()
    try:
        os.makedirs(CACHE, exist_ok=True)
        tmp = k + f".{os.getpid()}.tmp"
        with open(tmp, "wb") as fh:
            pickle.dump(v, fh, protocol=4)
        os.replace(tmp, k)   # atomic, so concurrent readers never see a partial file
    except Exception:
        pass
    return v

# ---------------------------------------------------------------- PAF / the edge rule
def paf_pairs(path, min_identity=None, min_coverage=None, denom="min", drop_self=True):
    """Best record per unordered pair, with identity and coverage on the CANONICAL definitions.

    Mirrors denovo_pipeline.rs::nucleotide_edges:
        identity = 1 - de           (`de:f:` = minimap2 gap-compressed divergence: an indel of any length
                                    counts ONCE, not once per base)
        coverage = aligned span on the SHORTER sequence / len(shorter)
                   i.e.  qlen <= tlen ?  (qend-qstart)/qlen  :  (tend-tstart)/tlen

    ⚠ `denom` selects the denominator and BOTH are real rules in the codebase:
        "min"  the shipped rule. STRUCTURALLY BLIND to truncation — a 10% fragment aligning fully into a
               complete sibling scores 1.00, so it cannot tell a copy from a piece of one.
        "max"  RUSTLE_ER_COVERAGE_LONGER. Demands BOTH members be covered (= BLAST qcovs AND scovs).
    ⚠⚠ DEFECT M1, FIXED IN THE RUST AND HERE ON 2026-08-10. This docstring used to read: "Coverage is
      NOT a proper fraction under 'min': it exceeds 1.0 when the query is the longer sequence ...
      That is in the SHIPPED formula, not a bug here." It WAS a bug, in both places. The numerator was
      the QUERY-axis aligned span while the denominator could be the TARGET's length — two different
      sequences' measurements divided, which is not an aligned fraction of anything.
      It bites because the shipped tier passes `-X` (⟹ `--dual=no`), so ONE orientation per pair is
      emitted and the query is NOT necessarily the shorter sequence: measured on the 25-region gorilla
      control panel, the query is the LONGER sequence in 54.5% of records and the old expression
      returned coverage above 1.0 on 8.2% of them, maximum 1.142.
      THE FIX: the numerator's AXIS FOLLOWS THE DENOMINATOR. The statistic is now bounded by 1.0 and
      SYMMETRIC under exchanging query and target — the property the old form lacked.
      ⚠ Any cached pickle written before 2026-08-10 holds the OLD values; the cache key was bumped
      ("covform=axis") so stale entries can never be served for a new call.
    ⚠ The effective identity floor in the pipeline is `sensitive_identity` (0.60), NOT RefineParams'
      min_identity (0.80), because nucleotide_sensitive defaults TRUE and E_r is the UNION of both tiers.
      Pass the floor you actually mean.

    ⚠ THE FLOORS ARE APPLIED PER RECORD, NOT TO A PRE-SELECTED BEST RECORD. This was a real bug here
      until 2026-08-07: the old code kept the ARGMAX-BY-COVERAGE record for each pair and only then
      tested the floors, so a pair whose highest-coverage record failed identity lost its edge even when
      another record cleared BOTH floors. denovo_pipeline.rs:3687 tests every record and UNIONs. Proof
      (2-record synthetic, floors 0.60/0.50): rec1 cov 0.90 id 0.55, rec2 cov 0.70 id 0.95 -> the old
      code returned no edge, the Rust returns the edge. On a 78-locus segdup PAF (336,061 records) the
      two rules happen to agree, so this cannot be caught by eyeballing a real fixture.
    ⚠ The returned identity/coverage are the EXEMPLAR: the highest-coverage PASSING record, ties broken
      on higher identity. That mirrors denovo_pipeline.rs:3691 exactly. Reporting only — the KEY set is
      the same whichever passing record is kept.

    Returns {(a, b): {"identity", "coverage", "aln", "qlen", "tlen", "q", "t"}} with a < b.
    """
    if denom not in ("min", "max"):
        raise ValueError("denom must be 'min' or 'max'")

    def build():
        best = {}
        with open(path) as fh:
            for ln in fh:
                f = ln.rstrip("\n").split("\t")
                if len(f) < 12:
                    continue
                q, t = f[0], f[5]
                if drop_self and q == t:
                    continue
                try:
                    qlen, qs, qe, tlen = int(f[1]), int(f[2]), int(f[3]), int(f[6])
                except ValueError:
                    continue
                de = None
                for x in f[12:]:
                    if x.startswith("de:f:"):
                        de = float(x[5:])
                        break
                if de is None:
                    # ⚠ Do NOT drop the record. denovo_pipeline.rs:3635 falls back to nmatch/blocklen
                    # when `de:f:` is absent, so dropping it silently loses edges the pipeline keeps.
                    try:
                        ident = int(f[9]) / max(int(f[10]), 1)
                    except ValueError:
                        continue
                else:
                    ident = 1.0 - de
                ts, te = int(f[7]), int(f[8])
                # M1: whichever sequence supplies the DENOMINATOR also supplies the NUMERATOR.
                side_q = (qlen <= tlen) if denom == "min" else (qlen >= tlen)
                d = (min(qlen, tlen) if denom == "min" else max(qlen, tlen)) or 1
                aln_on_denom_axis = (qe - qs) if side_q else (te - ts)
                cov = aln_on_denom_axis / d
                # Floors PER RECORD, before any selection. This is the whole point.
                if min_identity is not None and ident < min_identity:
                    continue
                if min_coverage is not None and cov < min_coverage:
                    continue
                rec = {"identity": ident, "coverage": cov,
                       # `aln` is the span on the DENOMINATOR's axis, i.e. the one `coverage` actually
                       # used. The query-axis value stays available as `aln_query` for callers that
                       # genuinely mean the query.
                       "aln": aln_on_denom_axis, "aln_query": qe - qs,
                       "qlen": qlen, "tlen": tlen, "q": q, "t": t}
                k = (q, t) if q < t else (t, q)
                cur = best.get(k)
                # Exemplar among the PASSING records. Coverage is NEVER summed across records: two loci
                # sharing 60% in four scattered blocks are sharing fragments, not a gene.
                if cur is None or (cov, ident) > (cur["coverage"], cur["identity"]):
                    best[k] = rec
        return best

    # ⚠ The floors are now part of build(), so they MUST be part of the cache key.
    # ⚠ "covform=axis" is in the key ON PURPOSE: the coverage formula changed on 2026-08-10, so a pickle
    # written by the old code must NEVER be served for a new call. Bump this string if it changes again.
    return _cached(path, ("paf", "covform=axis", denom, drop_self, min_identity, min_coverage), build)

def edges_from_paf(path, min_identity=0.80, min_coverage=0.50, denom="min"):
    """The E_r edge set: pairs where ONE record clears both floors. Mirrors homology_edges_all_reps_pooled."""
    return set(paf_pairs(path, min_identity, min_coverage, denom).keys())

# ---------------------------------------------------------------- GFF
def genes_on(chrom, gff=GFF, kind="gene"):
    """Gene table for one contig as [(start0, end, name)], 0-BASED HALF-OPEN, sorted by start.

    ⚠ GFF is 1-BASED INCLUSIVE. The conversion is start0 = col4 - 1, end = col5. This is the SAME
      convention CIGAR-derived coordinates use (SAM POS is 1-based; ref0 = POS - 1), verified against
      reference splice dinucleotides which agree ONLY at shift 0 (166 exact hits, 0 at +/-1).
    Cached — this is a full gzip scan of a multi-hundred-MB file.
    """
    def build():
        out = []
        with gzip.open(gff, "rt") as fh:
            for ln in fh:
                if ln[0] == "#":
                    continue
                f = ln.rstrip("\n").split("\t")
                if len(f) < 9 or f[2] != kind or f[0] != chrom:
                    continue
                m = re.search(r"Name=([^;]+)", f[8])
                if m:
                    out.append((int(f[3]) - 1, int(f[4]), m.group(1)))
        out.sort()
        return out
    return _cached(gff, ("genes", chrom, kind), build)

def introns_of_transcripts(chrom, gene_name, gff=GFF, min_len=20):
    """{transcript_id: [(donor0, acceptor0), ...]} for one gene, 0-based, per-transcript (NOT unioned).

    ⚠ Scoring a single linear model against the UNION of all transcripts is unfair by construction: the
      union mixes mutually exclusive isoforms and no single path can reach it. Compare against ONE
      transcript, and report which.
    ⚠ `min_len` drops sub-20bp RefSeq transcript-to-genome ALIGNMENT GAPS, which are not splices. 21 of 229
      annotated "introns" across the NPIP loci were such artifacts. Set min_len=0 to keep them, but say so.
    """
    def build():
        exons = defaultdict(list)
        want = f"gene={gene_name};"
        with gzip.open(gff, "rt") as fh:
            for ln in fh:
                if ln[0] == "#":
                    continue
                f = ln.rstrip("\n").split("\t")
                if len(f) < 9 or f[0] != chrom or f[2] != "exon":
                    continue
                if want not in f[8] and f"gene={gene_name}" not in f[8]:
                    continue
                m = re.search(r"Parent=([^;]+)", f[8])
                if m:
                    exons[m.group(1)].append((int(f[3]) - 1, int(f[4])))
        out = {}
        for tx, ex in exons.items():
            ex.sort()
            iv = [(ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1)]
            out[tx] = [(d, a) for d, a in iv if a - d >= min_len]
        return out
    return _cached(gff, ("tx_introns", chrom, gene_name, min_len), build)

# ---------------------------------------------------------------- catalogs
_COPIES_MIN = ["family_id", "copy_idx", "tid", "chrom", "start", "end"]

def load_copies(path, require_exons=False):
    """copies.tsv -> [{col: value}] with the numeric columns typed and the SCHEMA CHECKED.

    ⚠ Schema drift is real: some catalogs have 9 columns and NO `exons` column, others have 10. A script
      that indexes positionally will silently read the wrong field. Pass require_exons=True if you need it.
    ⚠⚠ BEFORE QUOTING ANY CATALOG, CHECK HOW ITS BAM WAS BUILT. A mini-BAM restricted to some BED
      (`samtools view -M -L ...`) yields loci only where that BED had reads, and this has produced three
      separate retractions ("7/19 NPIP genes lost", "13 families", "4 families"). Absence of a locus in a
      catalog is NOT evidence of non-detection unless the input covered it.
    """
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        missing = [c for c in _COPIES_MIN if c not in hdr]
        if missing:
            raise ValueError(f"{path}: missing required columns {missing}; header was {hdr}")
        if require_exons and "exons" not in hdr:
            raise ValueError(f"{path}: no `exons` column — this catalog predates it; re-run or drop the need")
        idx = {c: i for i, c in enumerate(hdr)}
        rows = []
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < len(hdr):
                continue
            r = {c: f[idx[c]] for c in hdr}
            for c in ("start", "end", "n_exon", "n_reads", "copy_idx"):
                if c in r:
                    try:
                        r[c] = int(r[c])
                    except ValueError:
                        pass
            rows.append(r)
    return rows

def families_from_copies(rows):
    """{family_id: [row, ...]} — the partition the catalog asserts."""
    fam = defaultdict(list)
    for r in rows:
        fam[r["family_id"]].append(r)
    return dict(fam)

# ---------------------------------------------------------------- BAM chains
def chains_in(bam, chrom, start, end, samtools="samtools"):
    """Intron chains in a region: [{"introns", "n_reads", "start", "end"}], sorted by the CANONICAL key.

    A "transcript" here = reads sharing an EXACT intron chain, matching pass1_skeletons_robust.

    ⚠ ALWAYS `-F 2308` (primary only) — this is applied here and must never be relaxed for a per-read
      CIGAR statistic.
    ⚠ `N` = intron. M/=/X/D consume reference; I/S/H/P do not. SAM POS is 1-based -> ref0 = POS - 1.
    ⚠ CANONICAL SORT: (-n_reads, start, tuple(introns)). Four bench scripts previously used four different
      sort keys, making `chain_id` non-comparable across genes. Index into this list ONLY via this order.
    ⚠ A region query returns reads OVERLAPPING the region, so chains include junctions OUTSIDE it. Per-gene
      junction totals are therefore NOT per-gene quantities (measured: 43-44% of junctions fall outside).
    Cached on (bam mtime, region).
    """
    def build():
        out = subprocess.run([samtools, "view", "-F", "2308", bam, f"{chrom}:{start}-{end}"],
                             capture_output=True, text=True).stdout
        groups = defaultdict(lambda: {"n_reads": 0, "start": None, "end": None})
        for ln in out.splitlines():
            f = ln.split("\t")
            if len(f) < 6:
                continue
            try:
                ref = int(f[3]) - 1
            except ValueError:
                continue
            introns, n = [], 0
            beg = ref
            for ch in f[5]:
                if ch.isdigit():
                    n = n * 10 + ord(ch) - 48
                else:
                    if ch in "M=XD":
                        ref += n
                    elif ch == "N":
                        introns.append((ref, ref + n))
                        ref += n
                    n = 0
            g = groups[tuple(introns)]
            g["n_reads"] += 1
            g["start"] = beg if g["start"] is None else min(g["start"], beg)
            g["end"] = ref if g["end"] is None else max(g["end"], ref)
        rows = [{"introns": list(k), "n_reads": v["n_reads"], "start": v["start"], "end": v["end"]}
                for k, v in groups.items()]
        rows.sort(key=lambda r: (-r["n_reads"], r["start"], tuple(r["introns"])))
        return rows
    return _cached(bam, ("chains", chrom, start, end), build)

# ---------------------------------------------------------------- interval helpers
def best_overlap(a, b, intervals):
    """Index of the interval with the LARGEST overlap with [a, b), or None.

    ⚠ THIS EXISTS BECAUSE `break`-ON-FIRST-OVERLAP HAS BEEN WRITTEN AT LEAST THREE TIMES IN THIS PROJECT
      (NOTCH2NLB, the AMY1 comparison, NPIPB9) and each time produced a plausible WRONG number. First
      overlap is almost never what you want; if you genuinely want it, write it inline and say why.
    `intervals` is [(start, end, ...)] — extra tuple fields are ignored.
    """
    best, bo = None, 0
    for i, iv in enumerate(intervals):
        o = min(b, iv[1]) - max(a, iv[0])
        if o > bo:
            best, bo = i, o
    return best

def in_band(values, lo=0.8, hi=1.25):
    """Fraction of ratio-to-truth values inside [lo, hi].

    ⚠ FOR A RATIO AGAINST TRUTH, REPORT THIS — NEVER THE MEDIAN. A median scores 1.86x as nearer to truth
      than 0.65x, so a change that pushes loci from undershoot to OVERSHOOT improves the median while
      moving nothing into band. That exact error was made this session: median 0.41 -> 0.77 looked like a
      fix while the in-band fraction stayed at 2/19, the same two loci.
    """
    v = [x for x in values if x is not None]
    return (sum(1 for x in v if lo <= x <= hi), len(v))

if __name__ == "__main__":
    print(__doc__)
    print(f"GFF   {GFF}\nFASTA {FASTA}\nCACHE {CACHE}")
    if len(sys.argv) > 1 and sys.argv[1] == "selftest":
        g = genes_on("chr16")
        npip = [x for x in g if x[2].upper().startswith("NPIP")]
        print(f"selftest: chr16 has {len(g):,} genes, {len(npip)} NPIP")
        i = best_overlap(28623033, 28645151, [(0, 10), (28620000, 28630000), (28640000, 28650000)])
        assert i == 1, f"best_overlap picked {i}, expected 1 (largest, not first)"
        print("selftest: best_overlap returns LARGEST not FIRST — ok")
        n, d = in_band([0.65, 0.9, 1.86, 1.0])
        assert (n, d) == (2, 4), (n, d)
        print("selftest: in_band ok")

# ================================================================ ADDED 2026-08-07
# Everything below was written ad hoc in bench/crossspecies during the 08-06/07 session and got at least
# FIVE things wrong in ways that produced plausible, quotable, wrong numbers. Each is now here once.
#
#   1. AGGREGATED vs SINGLE-RECORD coverage. One script summed coverage across records while another used
#      the shipped single-record rule; the PAF figure then contradicted the graph figure in the same
#      artifact. -> er_edges(aggregate=False) is the shipped default; aggregate=True must be asked for.
#   2. UNIFORM vs TWO-TIER identity floor. Applying 0.80 to both tiers gave RNA density 0.547/0.698, which
#      had to be retracted; the shipped floors are asm20 0.80 / sensitive 0.60. -> er_edges takes a list
#      of (paf, floor) pairs so the tiers cannot be collapsed by accident.
#   3. ARGMAX BEST-HIT counting. Counting each locus's single best target said "6/19 orthologs"; counting
#      by QUALIFYING EDGE said 19/19. Best-hit DEFLATES here and INFLATES elsewhere. -> count_by_edge.
#   4. NODE NAME CONVENTION. The harness sanitises "chr16:100-200" to "L~chr16_99_200" (0-based!). Looking
#      up the un-sanitised form found nothing and was nearly reported as "no alignment records at all".
#      -> canon_node / parse_node.
#   5. READS OVERLAPPING vs CONTAINED IN a window. Selecting reads that merely overlap a locus pulled the
#      NEIGHBOUR's transcripts in and made them the "representative". -> read_exons(contained=...).
#
# ⚠ SHIPPED TIER, 2026-08-07: asm20 was RETIRED from the homology path (RUSTLE_ER_SENSITIVE_ONLY now
#   defaults TRUE). The shipped tier list is therefore SENSITIVE ALONE at identity 0.60. asm20 is a subset
#   on genomic sequence (0 unique edges) but contributes 1 unique edge on the spliced substrate.

SHIPPED_TIERS = (("sensitive", 0.60),)      # -k11 -w5 ; asm20 0.80 retired 2026-08-07
# Named so callers stop writing the bare 0.60 / 0.80 literals. The SENSITIVE tier pairs with
# `sensitive_identity`; ASM20_IDENTITY is `min_identity` and applies ONLY to the retired asm20 leg.
SENSITIVE_IDENTITY = 0.60
ASM20_IDENTITY = 0.80
MIN_COVERAGE = 0.50
GAMMA = 0.40

# ⚠⚠ DEFECT B1 (fixed 2026-08-10). The tier below is the ALIGNER half of the rule, and until today no
# bench script stated it: the eight bench/crossspecies panels ran `-c --eqx -N 200 -p 0.02` with NO `-X`,
# plus an `-x asm20` leg the shipped default SKIPS. `-N`/`-p` are INERT at this tier; `-X` is the
# operative difference because it implies `--dual=no`. Measured over 14 panels on BYTE-IDENTICAL FASTA:
# the partition differs on 4/14 and the edge count on 10/14. Mirrors ER_TIER_FLAGS / ER_SENSITIVE_SEED
# in denovo_pipeline.rs — change it THERE first, then here.
#
# ⚠ THIS IS THE ALL-VS-ALL (E_r) TIER ONLY. A seed->genome pass legitimately keeps `-N 200 -p 0.02` and
# drops `-X`: query and target are different files, secondaries ARE the answer, and `-X` would be wrong.
ER_TIER_FLAGS = ("-c", "-X", "--no-long-join")
ER_SENSITIVE_SEED = ("-k", "11", "-w", "5")
ER_ASM20_SEED = ("-x", "asm20")


def er_tier_argv(threads=4, seed=ER_SENSITIVE_SEED):
    """The exact argv the binary runs, minus the two FASTA operands."""
    return ["minimap2", *ER_TIER_FLAGS, "-t", str(max(int(threads), 1)), *seed]


def er_provenance(threads=4, seed=ER_SENSITIVE_SEED):
    """One line naming the tier a number was computed under. Print it next to every number."""
    return " ".join(er_tier_argv(threads, seed)) + " <fa> <fa>"


def ava(src, out_paf, threads=4, seed=ER_SENSITIVE_SEED):
    """All-vs-all self-alignment at the SHIPPED tier. Writes ONE PAF and returns its path.

    ⚠ ATOMIC ON PURPOSE (2026-08-11). It used to stream straight into `out_paf`. Every caller guards
      re-runs with `if not os.path.exists(paf) or os.path.getsize(paf) == 0`, so a run killed part-way
      (OOM, TaskStop, a WSL2 crash) left a TRUNCATED but non-empty PAF that the guard then accepted as
      finished — silently fewer edges, no error anywhere. Observed: a stopped chr16-sized run left a
      stub behind. Writing to a temp and renaming only on a clean exit makes a half-PAF impossible.
    ⚠ The output name must encode the TIER (`*.er.paf`) — see `_cached`. The cache cannot see flags.
    """
    tmp = f"{out_paf}.{os.getpid()}.tmp"
    try:
        with open(tmp, "w") as fh:
            subprocess.run(er_tier_argv(threads, seed) + [src, src],
                           stdout=fh, stderr=subprocess.DEVNULL, check=True)
        os.replace(tmp, out_paf)          # atomic: readers see the whole file or no file
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)
    return out_paf


def canon_node(region: str) -> str:
    """"chr16:100-200" -> "L~chr16_99_200". The harness's node id: 0-BASED start, underscores.

    ⚠ This exists because looking up the un-sanitised form returns nothing and reads as "no alignments".
    """
    chrom, rng = region.rsplit(":", 1)
    a, b = rng.split("-")
    return f"L~{chrom}_{int(a) - 1}_{b}"


def parse_node(node: str):
    """Inverse of canon_node. Accepts either form. Returns (chrom, start0, end)."""
    if node.startswith("L~"):
        chrom, a, b = node[2:].rsplit("_", 2)
        return chrom, int(a), int(b)
    chrom, rng = node.rsplit(":", 1)
    a, b = rng.split("-")
    return chrom, int(a) - 1, int(b)


def er_edges(tiers, nodes=None, min_coverage=MIN_COVERAGE, denom="min", aggregate=False,
             min_block=200):
    """E_r over one or more alignment tiers. `tiers` = [(paf_path, identity_floor), ...].

    Mirrors denovo_pipeline.rs::homology_edges_all_reps_pooled. Each tier is thresholded at ITS OWN
    identity floor and the results are UNIONed — the floors must never be collapsed to one value.

    aggregate=False (SHIPPED): a pair needs ONE record clearing both floors. Coverage is not summed.
    aggregate=True: mirrors `RUSTLE_ER_SUM_COVERAGE=1` (denovo_pipeline.rs:3695-3733). Coverage is the
        union of aligned query intervals across records. This recovers edges on fragmented
        (concatenated-exon) targets but it is NOT what ships. Ask for it explicitly and say so.
        ⚠ Its "+53 RNA edges / +0 DNA edges" was measured before 2026-08-07 with nmatch/blocklen identity
          and none of the Rust's guards; that figure is STALE and must be re-measured before quoting.
        The Rust's three guards, all of which exist because summing is what the per-record rule prevents,
        are reproduced here: records are grouped by STRAND (a real fragmented gene is collinear, a repeat
        matches both ways), only records with query span >= `min_block` count (default 200 bp, so a swarm
        of short repeat hits cannot accumulate past the floor), and a summed edge is only ADDED for pairs
        the per-record rule did not already return.
    """
    E = set()
    for path, floor in tiers:
        per_record = {k for k in edges_from_paf(path, floor, min_coverage, denom)
                      if nodes is None or (k[0] in nodes and k[1] in nodes)}
        E |= per_record
        if not aggregate:
            continue
        iv = defaultdict(list)
        lens = {}
        with open(path) as fh:
            for ln in fh:
                f = ln.rstrip("\n").split("\t")
                if len(f) < 12:
                    continue
                q, t = f[0], f[5]
                if q == t:
                    continue
                if nodes is not None and (q not in nodes or t not in nodes):
                    continue
                qlen, qs, qe, tlen = int(f[1]), int(f[2]), int(f[3]), int(f[6])
                # ⚠ 1 - de, NOT nmatch/blocklen. denovo_pipeline.rs:3695 reuses the SAME `ident` the
                # per-record rule used; using nm/bl here made the two branches disagree on identity.
                de = None
                for x in f[12:]:
                    if x.startswith("de:f:"):
                        de = float(x[5:])
                        break
                ident = (1.0 - de) if de is not None else int(f[9]) / max(int(f[10]), 1)
                # M1 applies to the summed rule too: the intervals are UNIONed into ONE denominator, so
                # they must be measured on that denominator's axis or the union adds two different
                # coordinate systems. `min_block` is a length on the same axis.
                ts, te = int(f[7]), int(f[8])
                side_q = (qlen <= tlen) if denom == "min" else (qlen >= tlen)
                a_s, a_e = (qs, qe) if side_q else (ts, te)
                if ident < floor or (a_e - a_s) < min_block:
                    continue
                k = (q, t) if q < t else (t, q)
                strand = f[4]
                iv[(k, strand)].append((a_s, a_e))
                lens[(k, strand)] = (qlen, tlen)
        for (k, _strand), spans in iv.items():
            if k in per_record:
                continue
            ql, tl = lens[(k, _strand)]
            d = (min(ql, tl) if denom == "min" else max(ql, tl)) or 1
            if _union_len(spans) / d >= min_coverage:
                E.add(k)
    return E


def _union_len(spans):
    tot = 0
    cs = ce = None
    for s, e in sorted(spans):
        if cs is None:
            cs, ce = s, e
        elif s > ce:
            tot += ce - cs
            cs, ce = s, e
        else:
            ce = max(ce, e)
    return tot + (ce - cs if cs is not None else 0)


def components(nodes, edges):
    """Connected components, largest first. Isolated nodes are returned as singletons."""
    parent = {n: n for n in nodes}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for a, b in edges:
        if a in parent and b in parent:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb
    out = defaultdict(list)
    for n in nodes:
        out[find(n)].append(n)
    return sorted((sorted(v) for v in out.values()), key=len, reverse=True)


def density(component, edges):
    """Induced edge density. NaN for a single node — never 0, which would read as 'disconnected'."""
    n = len(component)
    if n < 2:
        return float("nan")
    s = set(component)
    e = sum(1 for a, b in edges if a in s and b in s)
    return 2 * e / (n * (n - 1))


def count_by_edge(component, members, edges=None):
    """How many MEMBERS the component covers, counted by overlap — NOT by argmax best-hit.

    ⚠ Best-hit counting deflated a real result 19 -> 6 (every locus's single best target happened to be
      one of ~6 genes) and has inflated others 1.7x. Count what qualifies, not what wins.

    `members` = [(chrom, start0, end, name)]. Returns (covered_names, pure_count).
    """
    covered = set()
    pure = 0
    for node in component:
        c, a, b = parse_node(node)
        hit = [nm for mc, ms, me, nm in members if mc == c and ms < b and a < me]
        if hit:
            covered.update(hit)
            pure += 1
    return covered, pure


def score_family(component, members, edges=None):
    """Recall / purity / density for one component against one family's curated members.

    ⚠ PURITY IS A LOWER BOUND. The curated catalogs are provably incomplete — PPIAL4G/H, 8 RefSeq NPIP
      genes and LOC105371131 (= CAT's NPIPB5) all sit in the 'impure' bucket while being real members.
      Report it as a floor, never as an estimate.
    """
    covered, pure = count_by_edge(component, members)
    n = len(component)
    return {"recall": len(covered) / max(len(members), 1),
            "purity": pure / max(n, 1),
            "n": n,
            "density": density(component, edges) if edges is not None else float("nan"),
            "covered": sorted(covered)}


# ---------------------------------------------------------------- multimapper presence (stage N0)
def _bam_cached(bam, parts, build):
    """Cache keyed on the BAM's mtime+size plus the region. A re-aligned BAM invalidates itself."""
    return _cached(bam, parts, build)


def sec_frac(bam, chrom, start, end, samtools="samtools"):
    """N0: secondary / (primary + secondary) alignment records over a region. A PER-LOCUS SCALAR.

    Mirrors docs/seeded_family_definition.md §3b. It says "reads placed elsewhere also fit here" and
    never says WHERE, so it carries no pairwise information and may be used to build nodes.

    ⚠ KEEPS SECONDARY ALIGNMENTS ON PURPOSE (-F 2052 drops only unmapped + supplementary). The standing
      "-F 2308 primary only" rule guards per-read CIGAR statistics; here the secondaries ARE the
      measurement. At NPIPB11, 234 of 271 PRIMARIES carry MAPQ 60 while the window holds 25,287
      secondaries — a primary-only view of duplication sees almost nothing.
    ⚠ minimap2 emits NO NH tag; this is the NH-free equivalent.
    """
    reg = f"{chrom}:{start + 1}-{end}"

    def build():
        pri = subprocess.run([samtools, "view", "-c", "-F", "2308", bam, reg],
                             capture_output=True, text=True).stdout.strip()
        sec = subprocess.run([samtools, "view", "-c", "-f", "256", "-F", "2052", bam, reg],
                             capture_output=True, text=True).stdout.strip()
        p, s = int(pri or 0), int(sec or 0)
        return {"primary": p, "secondary": s, "sec_frac": s / (p + s) if p + s else 0.0}

    return _bam_cached(bam, ("secfrac", reg), build)


def read_exons(bam, chrom, start, end, contained=True, min_reads=3, samtools="samtools"):
    """Read-supported exon blocks in a window, as [(start0, end)] merged at >= min_reads support.

    ⚠ `contained=True` keeps only reads whose alignment lies INSIDE the window. With contained=False the
      neighbouring locus's transcripts are pulled in and can outvote the real ones — that is exactly how a
      figure ended up showing a chain from the adjacent locus as NPIPB9's representative (1,318 reads
      spanning 28,992,032-29,016,360, entirely left of the NPIPB9 window).
    ⚠ -F 2308: primary only, per the standing rule for per-read CIGAR statistics.
    ⚠ N in an RNA CIGAR is an intron spliced OUT and BREAKS the block; D does not.
    """
    reg = f"{chrom}:{start + 1}-{end}"

    def _build():
        return subprocess.run([samtools, "view", "-F", "2308", bam, reg],
                              capture_output=True, text=True).stdout

    out = _bam_cached(bam, ("view2308", reg), _build)
    REF = {"M", "=", "X", "D"}
    cnt = defaultdict(int)
    for ln in out.splitlines():
        f = ln.split("\t", 6)
        if len(f) < 6 or f[5] == "*":
            continue
        pos = int(f[3]) - 1
        ref = bs = pos
        num = ""
        blocks = []
        for ch in f[5]:
            if ch.isdigit():
                num += ch
                continue
            n = int(num or 0)
            num = ""
            if ch in REF:
                ref += n
            elif ch == "N":
                if ref > bs:
                    blocks.append((bs, ref))
                ref += n
                bs = ref
        if ref > bs:
            blocks.append((bs, ref))
        if not blocks:
            continue
        if contained and (blocks[0][0] < start or blocks[-1][1] > end):
            continue
        for b in blocks:
            cnt[b] += 1
    keep = []
    cur = None
    for (s, e) in sorted(cnt):
        s2, e2 = max(s, start), min(e, end)
        if e2 - s2 < 20:
            continue
        if cur and s2 <= cur[1]:
            cur = (cur[0], max(cur[1], e2), cur[2] + cnt[(s, e)])
        else:
            if cur and cur[2] >= min_reads:
                keep.append(cur[:2])
            cur = (s2, e2, cnt[(s, e)])
    if cur and cur[2] >= min_reads:
        keep.append(cur[:2])
    return keep


# ================================================================ SINGLE-PATH CEILING (2026-08-07)
# THE STRUCTURAL DENOMINATOR. "How well does this representative represent its locus?" has been answered
# with FIVE different denominators giving 0.33 -> 1.00 for the SAME model at the SAME locus (NPIPB6). Four
# of the five are not usable:
#
#   union of all 16 RefSeq transcripts   7/21 = 0.333   UNFAIR — the union mixes MUTUALLY EXCLUSIVE
#                                                       isoforms; no single linear path can ever reach it.
#   best-matching transcript             7/7  = 1.000   ARGMAX over 16 — inflates ~1.7x here. Banned:
#                                                       best-hit deflated 19->6 elsewhere for the same
#                                                       reason. An argmax denominator measures the argmax.
#   reads whose ENTIRE chain is on path  —      0.367   strictest, and a different question.
#   curated RefSeq Select NM_001395275.1 5/8  = 0.625   external and honest, but ONE curator's choice.
#   SINGLE-PATH CEILING (max chain = 12) 7/12 = 0.583   <- the FAIR STRUCTURAL denominator, below.
#
# A representative is ONE LINEAR PATH. It can therefore carry only a set of PAIRWISE NON-OVERLAPPING
# junctions. The largest such set is an intrinsic property of the locus's junction set J — it is not an
# argmax over models, it cannot be inflated by lengthening the representative, and it needs no threshold.
#
# ⚠⚠ THE RESULTING RATIO IS A **SATURATION**, NOT A PRECISION / ACCURACY / RECALL. Normalising by the
#   ceiling does NOT repair the tautology that killed "junction recall vs the >=3-read set": if J is the
#   locus's own read-supported set and the rep was BUILT from those reads, the rep's junctions are a subset
#   of J by construction and no false junction can ever be registered. What makes the record non-vacuous is
#   that the CEILING is a property of the LOCUS and moves independently of the model. Report it as "of what
#   one linear path could carry here, how much did it carry". Emitting a precision column against J would
#   be a third killed metric.
# ⚠ NAMING, load-bearing. The ceiling is a maximum CHAIN (pairwise NON-overlapping), not a maximum
#   antichain, and it is NOT a matching problem — it is longest-path-in-a-DAG / interval scheduling,
#   O(n log n). `n - maximum bipartite matching` computes the DUAL: the minimum number of linear paths
#   needed to carry ALL of J (= max ANTICHAIN, by Dilworth). On NPIPB6 those are 12 and 5. Conflating them
#   silently swaps one number for the other.
# ⚠⚠ MATCHING IS A MEASUREMENT HERE AND NOTHING ELSE. `min_path_cover` returns actual paths and they look
#   exactly like candidate representatives. THEY ARE NOT. If anything downstream ever consumes them as
#   loci, family members or reps, the O1 rule ("families are defined by sequence homology alone; never
#   build loci with bipartite matching / facility location") is broken. `locus_ceiling` therefore returns
#   COUNTS by default and only hands out the paths when explicitly asked.


def junction_conflict(j1, j2, touch_ok=True):
    """True iff NO single linear path can carry both junctions.

    A junction j = (donor0, acceptor0) is the intron as the HALF-OPEN reference interval [d, a): d is the
    first intronic base, a the first exonic base after the intron. Same convention as
    `introns_of_transcripts` and `chains_in`, so no coordinate is re-derived here.

    CONFLICT == the intervals intersect.  Case by case, and each case is a real annotation pattern:
      * disjoint (a1 <= d2)            compatible — and SUFFICIENT, not merely necessary: any pairwise
                                       disjoint intron set is realisable as one transcript by letting the
                                       exons fill the gaps. So the ceiling has no hidden feasibility test.
      * partial overlap                conflict — carrying j1 deletes the bases holding j2's donor site.
      * NESTING (j2 inside j1)         conflict — a strict subcase of intersection; needs no separate rule.
      * SHARED DONOR / SHARED ACCEPTOR conflict — one donor on one linear path splices to one acceptor.
      * TOUCHING (a1 == d2)            the only arguable case. Half-open intervals do not intersect, so
                                       touch_ok=True calls them compatible; but the exon between two
                                       abutting introns then has length ZERO. touch_ok=False requires at
                                       least ONE exonic base (no larger floor — that would be an arbitrary
                                       threshold). MEASURED: the two modes give IDENTICAL ceilings at all
                                       19 chr16 NPIP loci, so the choice moves no reported number here.
                                       The default follows the half-open convention used everywhere else.
    """
    d1, a1 = j1
    d2, a2 = j2
    if touch_ok:
        return not (a1 <= d2 or a2 <= d1)
    return not (a1 < d2 or a2 < d1)


def max_chain(J, touch_ok=True):
    """The largest pairwise-COMPATIBLE subset of J, as a sorted list — the junctions ONE path could carry.

    Earliest-ACCEPTOR greedy, the classic interval-scheduling argument; optimal, O(n log n), no matching
    and no threshold. Verified against full brute-force enumeration on 8,000 randomised cases (n <= 9,
    both touch modes, seed 20260807): 0 mismatches.
    """
    out, last = [], None
    for d, a in sorted(set(J), key=lambda x: (x[1], x[0])):
        if last is None or (d >= last if touch_ok else d > last):
            out.append((d, a))
            last = a
    return out


def single_path_ceiling(J, touch_ok=True):
    """|max_chain(J)| — the most junctions ANY one linear representative could carry at this locus."""
    return len(max_chain(J, touch_ok))


def max_antichain(J, touch_ok=True):
    """Largest pairwise-CONFLICTING set = max clique of the interval graph = maximum overlap DEPTH.

    Sweep, O(n log n), no matching. By Dilworth this equals `min_path_cover`, and because the two are
    computed by completely different means their agreement is a live check on both — `locus_ceiling`
    asserts it on every call.
    """
    ev = []
    for d, a in set(J):
        ev.append((d, +1))
        ev.append((a, -1))
    # touch_ok: [d, a) is half-open, so a CLOSE at x is processed before an OPEN at x.
    ev.sort(key=lambda x: (x[0], x[1]) if touch_ok else (x[0], -x[1]))
    cur = best = 0
    for _, t in ev:
        cur += t
        best = max(best, cur)
    return best


def _hopcroft_karp(n_left, adj):
    """Maximum bipartite matching. adj[u] = list of right-vertex indices. Returns (size, matchL)."""
    INF = float("inf")
    n_right = 1 + max((v for vs in adj for v in vs), default=-1)
    matchL = [-1] * n_left
    matchR = [-1] * n_right
    dist = [0] * n_left

    def bfs():
        q = deque()
        for u in range(n_left):
            if matchL[u] == -1:
                dist[u] = 0
                q.append(u)
            else:
                dist[u] = INF
        found = False
        while q:
            u = q.popleft()
            for v in adj[u]:
                w = matchR[v]
                if w == -1:
                    found = True
                elif dist[w] == INF:
                    dist[w] = dist[u] + 1
                    q.append(w)
        return found

    def dfs(u):
        for v in adj[u]:
            w = matchR[v]
            if w == -1 or (dist[w] == dist[u] + 1 and dfs(w)):
                matchL[u] = v
                matchR[v] = u
                return True
        dist[u] = INF
        return False

    size = 0
    while bfs():
        for u in range(n_left):
            if matchL[u] == -1 and dfs(u):
                size += 1
    return size, matchL


def min_path_cover(J, touch_ok=True, want_paths=False):
    """Minimum number of linear paths needed to carry ALL of J. ⚠ THIS IS NOT THE CEILING.

    Reduction (Fulkerson / Koenig): split every junction into an OUT-copy (left) and an IN-copy (right)
    and put an edge (u_out, v_in) iff u precedes v under the interval order (a_u <= d_v). Every edge joins
    an OUT to an IN, so the graph is bipartite by construction, not by luck. A matching is a set of
    successor links with in- and out-degree <= 1, i.e. a set of vertex-disjoint paths, so
    #paths = n - |matching|. The precedence relation is ALREADY TRANSITIVELY CLOSED (a1 <= d2 and
    a2 <= d3 give a1 <= d3), so path cover == CHAIN cover and Dilworth gives min chain cover == max
    antichain — which `max_antichain` computes independently.

    ⚠ n - |M| is the minimum CHAIN cover only because the edge set is the TRANSITIVE precedence relation.
      Feed it a non-transitive successor map (e.g. cothread co-observation) and the number means something
      else. Say which graph a quoted number came from.
    ⚠ Builds the pair graph in O(n^2). Fine at |J| <= 50 (the NPIP maximum is 24 annotated / 43 read-
      derived); for a large read-derived J use `max_antichain`, which is the same number by Dilworth.
    ⚠⚠ want_paths=True returns candidate paths. THEY ARE A MEASUREMENT ARTIFACT, NEVER REPRESENTATIVES.
      Consuming them to build loci or families breaks the O1 rule outright.

    Returns k, or (k, paths) when want_paths.
    """
    js = sorted(set(J))
    n = len(js)
    if n == 0:
        return (0, []) if want_paths else 0
    adj = [[] for _ in range(n)]
    for u in range(n):
        _du, au = js[u]
        for v in range(n):
            if u == v:
                continue
            dv, _av = js[v]
            if (au <= dv) if touch_ok else (au < dv):
                adj[u].append(v)
    m, matchL = _hopcroft_karp(n, adj)
    if not want_paths:
        return n - m
    succ = {u: matchL[u] for u in range(n) if matchL[u] != -1}
    has_pred = set(succ.values())
    paths = []
    for u in range(n):
        if u in has_pred:
            continue
        p, cur = [], u
        while cur is not None:
            p.append(js[cur])
            cur = succ.get(cur)
        paths.append(p)
    return n - m, paths


def locus_ceiling(J, touch_ok=True, want_paths=False):
    """The per-locus STRUCTURAL record for a junction set J. Every field is intrinsic to J.

    n_junctions   |J|, the union — an UNFAIR denominator on its own, kept so the ratio is auditable.
    ceiling       the single-path ceiling = |max chain|. The fair structural denominator.
    min_paths     minimum linear representatives needed to carry all of J (= max antichain, by Dilworth).
                  Arguably the sharper headline: NPIPB6 needs 5, so ONE rep is at best a fifth of it.
    NO RATIO IS STORED. The consumer divides; the file keeps counts, so a 0/0 locus can never be written
    out as a 0.0 that reads like a real score.
    """
    js = sorted(set(J))
    chain = max_chain(js, touch_ok)
    k = min_path_cover(js, touch_ok)
    ac = max_antichain(js, touch_ok)
    assert k == ac, f"Dilworth violated: min chain cover {k} != max antichain {ac}"
    rec = {
        "n_junctions": len(js),
        "ceiling": len(chain),
        "min_paths": k,
        "max_antichain": ac,
        "ceiling_junctions": chain,
    }
    if want_paths:
        rec["cover_paths"] = min_path_cover(js, touch_ok, want_paths=True)[1]
    return rec
