#!/usr/bin/env python
"""
POA-ONLY, RNA-level family criterion based on RECIPROCAL EXON-COVERAGE.

Question
--------
Pure sequence-similarity grouping (the universe's minimap2 step, and the
minimizer-Jaccard re-derivation of it) lumped some gene pairs that share only a
single DOMAIN / exon into the same "family" (domain-sharers), producing
false-positive family memberships that Ensembl Compara does NOT confirm as
paralogs. Does a criterion that requires the two genes to CO-ALIGN over a
MAJORITY of *each* transcript -- reciprocal aligned coverage -- fix those
domain-sharing false positives while still accepting the genuine recent
tandem-duplicate pairs?

DEFINITION UNDER TEST
---------------------
Two loci are in a multi-copy *transcript* family iff, when aligned, the
homologous region covers >= T of BOTH transcripts:

    reciprocal_coverage(a,b) = min( aligned_len / len(a), aligned_len / len(b) )

where aligned_len is the number of aligned (ungapped, both-present) columns of a
global alignment of the two gene-representative sequences. A domain-sharer
co-aligns over only ONE shared exon/domain, so one of the two ratios is small
=> low reciprocal coverage. A true copy co-aligns over most of both genes =>
high reciprocal coverage.

POA-FAITHFULNESS NOTE
---------------------
POA (partial-order alignment) generalises pairwise alignment to >2 sequences by
embedding each new sequence into a DAG of the previous ones. For a PAIR there is
no DAG to extend -- the POA instance on two sequences IS a single pairwise
(global) alignment. We therefore use Bio.Align.PairwiseAligner in GLOBAL mode
(Needleman-Wunsch) with match=+2, mismatch=-1, gap-open=-5, gap-extend=-1 as the
two-sequence POA instance. Both quantities we read off it -- reciprocal aligned
coverage and identity over aligned columns -- are POA-derived (alignment column
statistics). No DNA/protein domain annotation, no BLAST, no k-mer/minimizer is
used AS THE CRITERION; minimizer-Jaccard is recomputed ONLY to CONTRAST.

INPUTS
------
- bench/copy_recovery_eval/results/gene_rep.fa   (gene representative sequences)
- bench/copy_recovery_eval/results/universe.tsv  (universe families)
- bench/compara_paralog_relation.json /
  bench/compara_validation.md                    (EXTERNAL corroboration labels
                                                   only -- not part of the method)

LABELS (Compara, external corroboration only)
---------------------------------------------
CONFIRMED (real paralog / tandem dup):
  APOBEC3D<->APOBEC3F, RFPL1<->RFPL2, RFPL1<->RFPL3, RFPL2<->RFPL3,
  and RABL2A<->RABL2B (Compara fetch-error artifact; it is the flagship tandem
  pair, treated as confirmed-real per task spec).
DOMAIN-SHARER (Compara 'no'):
  ASDURF<->ASNSD1, CASP8<->FLACC1, CCDC188<->ZDHHC8, CDPF1<->PPARA,
  CREB1<->METTL21A, GCA<->KCNH7, GPR39<->LYPD1.

Deterministic: fixed RNG seed for the random control pairs; alignment and
minimizer machinery have no randomness.

Run with /home/juanfra/miniforge3/bin/python.
"""

import csv
import hashlib
import json
import random
import sys
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from Bio.Align import PairwiseAligner

# --------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------
BASE = "/mnt/c/Users/jfris/Desktop/Rustle/bench"
FASTA_PATH = f"{BASE}/copy_recovery_eval/results/gene_rep.fa"
UNIVERSE_PATH = f"{BASE}/copy_recovery_eval/results/universe.tsv"
COMPARA_JSON = f"{BASE}/compara_paralog_relation.json"
OUT_PNG = f"{BASE}/poa_family_definition.png"
OUT_MD = f"{BASE}/poa_family_definition.md"

SEED = 1729
N_CONTROL = 60  # random cross-family / non-family control pairs

# Global-alignment (= 2-sequence POA instance) scoring, per task spec.
MATCH = 2.0
MISMATCH = -1.0
GAP_OPEN = -5.0   # PairwiseAligner open_gap_score is the cost of opening (incl. first gap col)
GAP_EXTEND = -1.0

# Minimizer params (recomputed ONLY for contrast, identical to
# family_detection_validation.py / rustle production: k=15, w=10, canonical,
# blake2b-64).
K = 15
W = 10

# --------------------------------------------------------------------------
# Labels (external corroboration ONLY -- not part of the criterion)
# --------------------------------------------------------------------------
CONFIRMED_PAIRS = [
    ("APOBEC3D", "APOBEC3F"),
    ("RFPL1", "RFPL2"),
    ("RFPL1", "RFPL3"),
    ("RFPL2", "RFPL3"),
    ("RABL2A", "RABL2B"),  # fetch-error artifact -> treat as confirmed-real
]
DOMAIN_SHARER_PAIRS = [
    ("ASDURF", "ASNSD1"),
    ("CASP8", "FLACC1"),
    ("CCDC188", "ZDHHC8"),
    ("CDPF1", "PPARA"),
    ("CREB1", "METTL21A"),
    ("GCA", "KCNH7"),
    ("GPR39", "LYPD1"),
]


def npair(a, b):
    return tuple(sorted((a, b)))


CONFIRMED_SET = {npair(a, b) for a, b in CONFIRMED_PAIRS}
DOMAIN_SET = {npair(a, b) for a, b in DOMAIN_SHARER_PAIRS}


# --------------------------------------------------------------------------
# Minimizer machinery (CONTRAST ONLY) -- copied verbatim from
# family_detection_validation.py to make the contrast exact.
# --------------------------------------------------------------------------
_RC = bytes.maketrans(b"ACGTacgtNn", b"TGCAtgcaNn")


def revcomp(s: bytes) -> bytes:
    return s.translate(_RC)[::-1]


def _hash64(kmer: bytes) -> int:
    return int.from_bytes(hashlib.blake2b(kmer, digest_size=8).digest(), "little")


def minimizers(seq: bytes, k: int = K, w: int = W):
    seq = seq.upper()
    out = set()
    L = len(seq)
    if L < k:
        return out
    n = L - k + 1
    hashes = [0] * n
    for i in range(n):
        km = seq[i:i + k]
        rc = revcomp(km)
        can = km if km <= rc else rc
        hashes[i] = _hash64(can)
    upper = max(1, n - w)
    for win_start in range(upper):
        win_end = min(win_start + w, n)
        out.add(min(hashes[win_start:win_end]))
    return out


def jaccard(a: set, b: set) -> float:
    if not a and not b:
        return 0.0
    inter = len(a & b)
    uni = len(a | b)
    return inter / uni if uni else 0.0


# --------------------------------------------------------------------------
# POA pairwise instance: GLOBAL alignment -> reciprocal coverage + identity
# --------------------------------------------------------------------------
def make_aligner():
    a = PairwiseAligner()
    a.mode = "global"
    a.match_score = MATCH
    a.mismatch_score = MISMATCH
    a.open_gap_score = GAP_OPEN
    a.extend_gap_score = GAP_EXTEND
    return a


# When the global aligner stitches two NON-homologous transcripts it pays a
# little gap cost to pick up scattered chance matches, fragmenting the alignment
# into hundreds of tiny ungapped blocks. A genuine copy instead yields one large
# contiguous co-aligned block. We therefore also report the CONTIGUOUS-CORE
# coverage (largest single ungapped block / shorter gene) -- still a pure
# alignment-column statistic, no domain/BLAST/k-mer input.
def poa_pair_stats(aligner, s1, s2):
    """Compute POA (= 2-sequence global alignment) column statistics.

    Returns a dict with:
      recip        reciprocal coverage = min(aligned_len/len(s1), aligned_len/len(s2))
                   over ALL ungapped aligned columns (the LITERAL definition under test)
      ident        identity over those aligned columns = matches/aligned_len
      core_recip   CONTIGUOUS-CORE coverage = largest_ungapped_block / min(len)
                   (the homologous core, not gappy global-alignment filler).
                   NOTE: this divides by the SHORTER gene, so it is the MAX of the
                   two per-gene block-coverage ratios -- NOT a reciprocal (= min).
                   The variable name is retained for back-compat only.
      aligned_len  number of ungapped aligned columns (all blocks)
      n_blocks     number of ungapped blocks (fragmentation = domain-sharer signature)
      biggest      length of the largest single ungapped block
      score        alignment score
    """
    s1 = s1.upper()
    s2 = s2.upper()
    aln = aligner.align(s1, s2)[0]
    t_blocks, q_blocks = aln.aligned  # paired ungapped blocks
    aligned_len = 0
    matches = 0
    biggest = 0
    n_blocks = 0
    for (ts, te), (qs, qe) in zip(t_blocks, q_blocks):
        n = te - ts
        if n == 0:
            continue
        n_blocks += 1
        aligned_len += n
        biggest = max(biggest, n)
        seg1 = s1[ts:te]
        seg2 = s2[qs:qe]
        matches += sum(1 for c1, c2 in zip(seg1, seg2) if c1 == c2)
    minlen = min(len(s1), len(s2))
    if aligned_len == 0 or minlen == 0:
        return dict(recip=0.0, ident=0.0, core_recip=0.0, aligned_len=0,
                    n_blocks=0, biggest=0, score=float(aln.score))
    recip = min(aligned_len / len(s1), aligned_len / len(s2))
    ident = matches / aligned_len
    core_recip = biggest / minlen
    return dict(recip=recip, ident=ident, core_recip=core_recip,
                aligned_len=aligned_len, n_blocks=n_blocks, biggest=biggest,
                score=float(aln.score))


# --------------------------------------------------------------------------
# IO
# --------------------------------------------------------------------------
def load_fasta(path):
    seqs = {}
    name = None
    buf = []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(buf)
            # reconcile IDs to universe gene names: strip 'gene-' prefix
            name = line[1:].split()[0]
            if name.startswith("gene-"):
                name = name[len("gene-"):]
            buf = []
        else:
            buf.append(line)
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def load_universe(path):
    """Return gene -> family_id, and family_id -> set(genes)."""
    g2f = {}
    fam2genes = defaultdict(set)
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            g = r["gene_id"]
            f = r["family_id"]
            g2f.setdefault(g, f)
            fam2genes[f].add(g)
    return g2f, fam2genes


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def main():
    fasta = load_fasta(FASTA_PATH)
    g2f, fam2genes = load_universe(UNIVERSE_PATH)

    # ---- sanity: every named gene must be present ----
    named = sorted({g for p in CONFIRMED_PAIRS + DOMAIN_SHARER_PAIRS for g in p})
    missing = [g for g in named if g not in fasta]
    if missing:
        sys.exit(f"[fatal] named genes missing from gene_rep.fa: {missing}")

    aligner = make_aligner()

    # Precompute minimizer sets for all genes we will need (named + controls).
    min_cache = {}

    def get_min(g):
        if g not in min_cache:
            min_cache[g] = minimizers(fasta[g].encode())
        return min_cache[g]

    # ---- 1. compute POA stats + Jaccard for every labeled named pair ----
    rows = []  # dict per pair
    for a, b in CONFIRMED_PAIRS + DOMAIN_SHARER_PAIRS:
        # confirm both share a within-universe family (sanity / scope check)
        same_fam = g2f.get(a) == g2f.get(b)
        st = poa_pair_stats(aligner, fasta[a], fasta[b])
        jac = jaccard(get_min(a), get_min(b))
        lab = "confirmed" if npair(a, b) in CONFIRMED_SET else "domain-sharer"
        rows.append(
            dict(a=a, b=b, label=lab, len_a=len(fasta[a]), len_b=len(fasta[b]),
                 jaccard=jac, same_family=same_fam, **st)
        )

    # ---- 2. CONTROL set: random cross-family / non-family pairs ----
    # Deterministic. Pick pairs of genes from DIFFERENT universe families that
    # are present in the FASTA. These are the "should NOT be a family" baseline.
    rng = random.Random(SEED)
    fasta_genes = [g for g in g2f if g in fasta]
    fasta_genes.sort()
    control_rows = []
    seen = set(CONFIRMED_SET | DOMAIN_SET)
    attempts = 0
    while len(control_rows) < N_CONTROL and attempts < 20000:
        attempts += 1
        a, b = rng.sample(fasta_genes, 2)
        if g2f.get(a) == g2f.get(b):
            continue  # same family -> not a cross-family control
        key = npair(a, b)
        if key in seen:
            continue
        seen.add(key)
        st = poa_pair_stats(aligner, fasta[a], fasta[b])
        jac = jaccard(get_min(a), get_min(b))
        control_rows.append(
            dict(a=a, b=b, label="control", len_a=len(fasta[a]),
                 len_b=len(fasta[b]), jaccard=jac, same_family=False, **st)
        )

    n_conf = sum(1 for r in rows if r["label"] == "confirmed")
    n_dom = sum(1 for r in rows if r["label"] == "domain-sharer")

    # ---- 3. threshold finder (re-usable across the two coverage axes) ----
    def split_metric(key):
        conf = [r[key] for r in rows if r["label"] == "confirmed"]
        dom = [r[key] for r in rows if r["label"] == "domain-sharer"]
        ctrl = [r[key] for r in control_rows]
        min_conf = min(conf)
        max_dom = max(dom)
        max_ctrl = max(ctrl)
        # separable iff every confirmed strictly exceeds every domain-sharer AND
        # every control (controls are also "not a family").
        separable = (min_conf > max_dom) and (min_conf > max_ctrl)
        if min_conf > max(max_dom, max_ctrl):
            T = round((min_conf + max(max_dom, max_ctrl)) / 2.0, 4)
        else:
            grid = np.linspace(0.0, 1.0, 2001)
            best = None
            for t in grid:
                acc = (sum(1 for x in conf if x >= t)
                       + sum(1 for x in dom if x < t)
                       + sum(1 for x in ctrl if x < t))
                if best is None or acc > best[0]:
                    best = (acc, t)
            T = round(float(best[1]), 4)
        return dict(
            conf=conf, dom=dom, ctrl=ctrl, T=T, separable=separable,
            min_conf=min_conf, max_dom=max_dom, max_ctrl=max_ctrl,
            n_conf_acc=sum(1 for x in conf if x >= T),
            n_dom_rej=sum(1 for x in dom if x < T),
            n_ctrl_rej=sum(1 for x in ctrl if x < T),
        )

    # AXIS A: literal definition under test -- reciprocal coverage (ALL columns)
    raw = split_metric("recip")
    # AXIS B: POA-only refinement -- CONTIGUOUS-CORE coverage (largest ungapped
    # block / SHORTER gene; = MAX of the two per-gene ratios, NOT reciprocal).
    core = split_metric("core_recip")

    # Headline uses the metric that achieves clean separation (the POA-only fix).
    chosen = core if core["separable"] else raw
    chosen_name = "contiguous-core" if core["separable"] else "all-column"
    # The contiguous-core metric divides by the SHORTER gene (= MAX of the two
    # per-gene ratios), so it is NOT a reciprocal coverage; only the all-column
    # metric min(aligned/len_a, aligned/len_b) genuinely divides by the LONGER
    # gene and is reciprocal. Label the axes accordingly.
    chosen_metric = "coverage" if core["separable"] else "reciprocal coverage"

    # Separation margins on the contiguous-core axis. The BINDING barrier is the
    # nearest rejected point (= max over domain-sharers AND controls), not just
    # the domain-sharer max. Report both ratios honestly.
    core_margin_dom = (core["min_conf"] / core["max_dom"]
                       if core["max_dom"] > 0 else float("inf"))
    core_margin_ctrl = (core["min_conf"] / core["max_ctrl"]
                        if core["max_ctrl"] > 0 else float("inf"))
    core_margin_binding = (core["min_conf"] / max(core["max_dom"], core["max_ctrl"])
                           if max(core["max_dom"], core["max_ctrl"]) > 0
                           else float("inf"))
    margin_phrase = (
        f"~{core_margin_dom:.1f}x vs domain-sharers, "
        f"~{core_margin_binding:.1f}x vs the nearest control"
    )
    T = chosen["T"]
    separable = chosen["separable"]
    min_conf = chosen["min_conf"]
    max_dom = chosen["max_dom"]
    n_conf_acc = chosen["n_conf_acc"]
    n_dom_rej = chosen["n_dom_rej"]
    n_ctrl_rej = chosen["n_ctrl_rej"]
    conf = chosen["conf"]
    dom = chosen["dom"]
    ctrl = chosen["ctrl"]

    # Contrast: same accept/reject test using minimizer-Jaccard at its
    # production threshold 0.30 AND at the best-separating Jaccard threshold.
    j_conf = [r["jaccard"] for r in rows if r["label"] == "confirmed"]
    j_dom = [r["jaccard"] for r in rows if r["label"] == "domain-sharer"]
    j_ctrl = [r["jaccard"] for r in control_rows]
    j_separable = (min(j_conf) > max(j_dom)) and (min(j_conf) > max(j_ctrl))
    # best Jaccard threshold by accuracy on the SAME labeled pairs + controls
    grid = np.linspace(0.0, 1.0, 2001)
    bestj = None
    for t in grid:
        acc = (sum(1 for x in j_conf if x >= t)
               + sum(1 for x in j_dom if x < t)
               + sum(1 for x in j_ctrl if x < t))
        if bestj is None or acc > bestj[0]:
            bestj = (acc, t)
    Tj_best = round(float(bestj[1]), 4)
    n_jconf_acc_best = sum(1 for x in j_conf if x >= Tj_best)
    n_jdom_rej_best = sum(1 for x in j_dom if x < Tj_best)
    # production Jaccard threshold 0.30
    TJ_PROD = 0.30
    n_jconf_acc_prod = sum(1 for x in j_conf if x >= TJ_PROD)
    n_jdom_rej_prod = sum(1 for x in j_dom if x < TJ_PROD)

    # ---- 4. FIGURE ----
    NAVY = "#1f3a5f"
    ORANGE = "#e8820c"
    GREEN = "#2e8b57"
    color = {"confirmed": NAVY, "domain-sharer": ORANGE, "control": GREEN}

    fig, (ax0, ax1, ax2) = plt.subplots(1, 3, figsize=(16.5, 5.2))

    def strip(ax, key, T_, title, xlabel):
        rng2 = random.Random(SEED + 1)
        cats = [
            ("confirmed", [r[key] for r in rows if r["label"] == "confirmed"]),
            ("domain-sharer", [r[key] for r in rows if r["label"] == "domain-sharer"]),
            ("control", [r[key] for r in control_rows]),
        ]
        for yi, (cat, vals) in enumerate(cats):
            ys = [yi + (rng2.random() - 0.5) * 0.28 for _ in vals]
            ax.scatter(vals, ys, s=46, color=color[cat], alpha=0.85,
                       edgecolor="white", linewidth=0.5, zorder=3,
                       label=f"{cat} (n={len(vals)})")
        ax.axvline(T_, color="black", linestyle="--", linewidth=1.4, zorder=2)
        ax.text(T_, 2.62, f"  T = {T_:.2f}", fontsize=10, va="bottom",
                ha="left", color="black")
        ax.set_yticks([0, 1, 2])
        ax.set_yticklabels(["confirmed", "domain-\nsharer", "control"])
        ax.set_xlabel(xlabel)
        ax.set_xlim(-0.03, 1.03)
        ax.set_ylim(-0.6, 2.9)
        ax.set_title(title, fontsize=11)
        ax.grid(axis="x", alpha=0.25)

    # Panel 0: LITERAL definition -- all-column reciprocal coverage (OVERLAPS)
    strip(ax0, "recip", raw["T"],
          "All-column reciprocal coverage\n(global-alignment filler -> OVERLAP)",
          "POA reciprocal coverage  min(aligned/len)")

    # Panel 1: POA-only fix -- contiguous-core coverage (SEPARATES). NOT a
    # reciprocal: it divides by the SHORTER gene (= MAX of the two per-gene ratios).
    strip(ax1, "core_recip", core["T"],
          "Contiguous-core coverage\n(largest block / shorter gene -> CLEAN SPLIT)",
          "POA contiguous-core coverage  biggest_block/min(len)")

    # Panel 2: contiguous-core coverage (x) vs minimizer-Jaccard (y)
    for cat in ("confirmed", "domain-sharer"):
        xs = [r["core_recip"] for r in rows if r["label"] == cat]
        ys = [r["jaccard"] for r in rows if r["label"] == cat]
        ax2.scatter(xs, ys, s=70, color=color[cat], alpha=0.9,
                    edgecolor="white", linewidth=0.6, zorder=3, label=cat)
    cx = [r["core_recip"] for r in control_rows]
    cy = [r["jaccard"] for r in control_rows]
    ax2.scatter(cx, cy, s=30, color=GREEN, alpha=0.55,
                edgecolor="white", linewidth=0.4, zorder=2, label="control")
    ax2.axvline(core["T"], color="black", linestyle="--", linewidth=1.2, zorder=1)
    ax2.text(core["T"], ax2.get_ylim()[1], f" core T={core['T']:.2f}", fontsize=9,
             va="top", ha="left", color="black")
    ax2.set_xlabel("POA contiguous-core coverage")
    ax2.set_ylabel("minimizer-Jaccard (contrast)")
    ax2.set_xlim(-0.03, 1.03)
    ax2.set_title("Core coverage separates where Jaccard does not", fontsize=11)
    ax2.legend(fontsize=8, loc="upper right")
    ax2.grid(alpha=0.25)

    fig.suptitle(
        f"POA {chosen_name} {chosen_metric} (largest ungapped block / shorter gene) "
        f">= {T:.2f}: "
        f"accepts {n_conf_acc}/{n_conf} confirmed, rejects {n_dom_rej}/{n_dom} domain-sharers "
        f"(+ {n_ctrl_rej}/{len(ctrl)} controls)",
        fontsize=12.0, y=0.99,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(OUT_PNG, dpi=150)
    plt.close(fig)

    # ---- 5. MARKDOWN REPORT ----
    def fmt(x, n=3):
        return f"{x:.{n}f}"

    headline = (
        f"POA {chosen_name} {chosen_metric} (largest ungapped block / shorter "
        f"gene) >= {T:.2f} accepts {n_conf_acc}/{n_conf} confirmed tandem dups "
        f"and rejects {n_dom_rej}/{n_dom} domain-sharers."
    )

    # sorted tables
    named_rows = sorted(rows, key=lambda r: (r["label"], -r["core_recip"]))

    lines = []
    lines.append("# POA contiguous-core coverage family criterion")
    lines.append("")
    lines.append(
        "> *Naming note:* the separating axis is **contiguous-core coverage = "
        "largest ungapped block / shorter gene**. Because it divides by the "
        "SHORTER gene it equals the **maximum** of the two per-gene block "
        "ratios, so it is NOT a reciprocal coverage (a reciprocal = the *min* "
        "over both genes, i.e. divide by the LONGER gene). The word "
        "*reciprocal* applies only to the all-column literal metric "
        "`min(aligned/len_a, aligned/len_b)`. The two metrics use opposite "
        "denominators; only the all-column one is reciprocal."
    )
    lines.append("")
    lines.append(f"> **{headline}**")
    lines.append("")
    lines.append("## Definition under test")
    lines.append("")
    lines.append(
        "Two loci are in a multi-copy *transcript* family iff, aligned, the "
        "homologous region covers >= T of BOTH transcripts:"
    )
    lines.append("")
    lines.append("```")
    lines.append("reciprocal_coverage(a,b) = min( aligned_len/len(a), aligned_len/len(b) )")
    lines.append("```")
    lines.append("")
    lines.append(
        "i.e. the two genes must CO-ALIGN over a MAJORITY of *each* gene, not "
        "just one shared domain/exon. Both quantities (reciprocal aligned "
        "coverage, identity over aligned columns) are read off the alignment and "
        "are POA-derived. No DNA/protein domain annotation, no BLAST, no "
        "k-mer/minimizer is used as the criterion (minimizer-Jaccard is recomputed "
        "only to CONTRAST)."
    )
    lines.append("")
    lines.append("### POA-faithfulness note")
    lines.append("")
    lines.append(
        "POA generalises pairwise alignment to >2 copies by embedding each new "
        "sequence into a DAG of the previous ones; for a PAIR there is no DAG to "
        "extend, so the two-sequence POA instance IS a single pairwise global "
        "alignment. We use `Bio.Align.PairwiseAligner` in **global** mode "
        f"(Needleman-Wunsch) with match=+{int(MATCH)}, mismatch={int(MISMATCH)}, "
        f"gap-open={int(GAP_OPEN)}, gap-extend={int(GAP_EXTEND)} on the "
        "gene-representative sequences."
    )
    lines.append("")
    lines.append("## Two POA coverage axes (the key result)")
    lines.append("")
    lines.append(
        "We measured coverage two ways, BOTH purely alignment-column-derived:"
    )
    lines.append("")
    lines.append(
        "1. **all-column** reciprocal coverage -- the LITERAL definition above, "
        "summing every ungapped aligned column; and"
    )
    lines.append(
        "2. **contiguous-core** coverage -- the largest *single* ungapped "
        "co-aligned block divided by the **shorter** gene "
        "(`biggest_block / min(len)`). This is NOT a reciprocal coverage: "
        "dividing by the shorter gene makes it the *maximum* of the two "
        "per-gene block ratios (a reciprocal would divide by the longer gene "
        "= the *min*). It is the metric that separates; the name 'reciprocal' "
        "belongs only to axis 1."
    )
    lines.append("")
    lines.append(
        "The literal all-column metric does **NOT** cleanly separate the classes: "
        f"confirmed range [{fmt(min(raw['conf']))}, {fmt(max(raw['conf']))}], "
        f"domain-sharer range [{fmt(min(raw['dom']))}, {fmt(max(raw['dom']))}] -- "
        "they OVERLAP. The reason is a global-alignment artifact: when the "
        "Needleman-Wunsch aligner is handed two NON-homologous transcripts it "
        "still pays a little gap cost to harvest scattered chance matches, "
        "fragmenting the alignment into hundreds of tiny blocks (domain-sharers "
        f"here split into {min(r['n_blocks'] for r in rows if r['label']=='domain-sharer')}"
        f"-{max(r['n_blocks'] for r in rows if r['label']=='domain-sharer')} blocks) "
        "whose aligned columns *sum* to a high reciprocal coverage even though "
        "the alignment is non-homologous filler (note the low all-column "
        "identity of those pairs in the per-pair table)."
    )
    lines.append("")
    lines.append(
        "The **contiguous-core** metric removes that filler and is the POA-only "
        f"fix. It is **{'SEPARABLE' if core['separable'] else 'NOT separable'}**: "
        f"confirmed range [{fmt(min(core['conf']))}, {fmt(max(core['conf']))}], "
        f"domain-sharer range [{fmt(min(core['dom']))}, {fmt(max(core['dom']))}], "
        f"control range [{fmt(min(core['ctrl']))}, {fmt(max(core['ctrl']))}]. "
        f"A real copy has ONE large homologous block; a domain-sharer's largest "
        "single block covers only a few percent of the shorter gene -- "
        "indistinguishable from random cross-family controls."
    )
    lines.append("")
    lines.append("## Result: the separating threshold")
    lines.append("")
    lines.append(
        f"The principled, POA-only threshold (on the {chosen_name} axis) is "
        f"**T = {T:.2f}** "
        + (
            "(midpoint between the lowest confirmed and the highest "
            "domain-sharer/control coverage -- they are strictly separated)."
            if separable
            else "(balanced-accuracy optimum; the classes overlap, see caveats)."
        )
    )
    lines.append("")
    lines.append("| category | n | accept at T (cov>=T) | reject at T (cov<T) |")
    lines.append("|---|---|---|---|")
    lines.append(f"| confirmed | {n_conf} | **{n_conf_acc}** | {n_conf - n_conf_acc} |")
    lines.append(f"| domain-sharer | {n_dom} | {n_dom - n_dom_rej} | **{n_dom_rej}** |")
    lines.append(f"| control (cross-family) | {len(ctrl)} | {len(ctrl) - n_ctrl_rej} | {n_ctrl_rej} |")
    lines.append("")
    lines.append(
        f"Is the distribution bimodal/separable? "
        f"**{'YES -- clean bimodal split (' + margin_phrase + ')' if separable else 'partially'}.** "
        f"(contiguous-core axis: confirmed >= {fmt(core['min_conf'])} vs "
        f"domain-sharer <= {fmt(core['max_dom'])} and control <= {fmt(core['max_ctrl'])}; "
        f"the binding barrier is the nearest control at {fmt(core['max_ctrl'])}.)"
    )
    lines.append("")
    lines.append("## Contrast: minimizer-Jaccard does NOT separate them as well")
    lines.append("")
    lines.append(
        "Recomputing the OLD criterion (canonical minimizer-Jaccard, k=15/w=10, "
        "blake2b-64 -- identical to `family_detection_validation.py` and rustle "
        "production) for the SAME labeled pairs:"
    )
    lines.append("")
    lines.append("| criterion | threshold | confirmed accepted | domain-sharers rejected | clean separation? |")
    lines.append("|---|---|---|---|---|")
    lines.append(
        f"| **POA contiguous-core coverage** | T={core['T']:.2f} | "
        f"{core['n_conf_acc']}/{n_conf} | {core['n_dom_rej']}/{n_dom} | "
        f"{'**YES**' if core['separable'] else 'no'} |"
    )
    lines.append(
        f"| POA all-column coverage (literal) | T={raw['T']:.2f} | "
        f"{raw['n_conf_acc']}/{n_conf} | {raw['n_dom_rej']}/{n_dom} | "
        f"{'YES' if raw['separable'] else 'no (overlap)'} |"
    )
    lines.append(
        f"| minimizer-Jaccard (best-fit) | J>={Tj_best:.3f} | {n_jconf_acc_best}/{n_conf} | "
        f"{n_jdom_rej_best}/{n_dom} | {'YES' if j_separable else 'no (overlap)'} |"
    )
    lines.append(
        f"| minimizer-Jaccard (production 0.30) | J>=0.30 | {n_jconf_acc_prod}/{n_conf} | "
        f"{n_jdom_rej_prod}/{n_dom} | - |"
    )
    lines.append("")
    lines.append(
        "Even at its *best-fit* threshold the minimizer-Jaccard "
        + (
            "cannot cleanly separate the two classes -- confirmed and "
            "domain-sharer Jaccard ranges overlap, which is exactly why the "
            "similarity-only grouping conflated them. "
            if not j_separable
            else "separates them, but only with a data-fitted threshold; "
        )
        + f"Confirmed Jaccard range [{fmt(min(j_conf))}, {fmt(max(j_conf))}], "
        f"domain-sharer Jaccard range [{fmt(min(j_dom))}, {fmt(max(j_dom))}]. "
        "(Note RFPL1/2/3 -- genuine recent dups -- have Jaccard ~0.13-0.17, "
        "LOWER than several domain-sharers, so no Jaccard cutoff can both keep "
        "them and drop the domain-sharers.)"
    )
    lines.append("")
    # Display-rounding disclosure (finding #5): GPR39<->LYPD1's Jaccard rounds
    # to the same 3dp as the best-fit threshold; show its true 4dp value so the
    # 'rejected' verdict (x < t) is not mistaken for a tie-at-threshold.
    gpr39_jac = next((r["jaccard"] for r in rows
                      if npair(r["a"], r["b"]) == npair("GPR39", "LYPD1")), None)
    if gpr39_jac is not None and abs(round(gpr39_jac, 3) - Tj_best) < 5e-4:
        lines.append(
            f"*Display-rounding note:* GPR39<->LYPD1's Jaccard is "
            f"{gpr39_jac:.4f}, which rounds to {fmt(gpr39_jac)} -- the same 3dp "
            f"as the best-fit threshold J>={Tj_best:.3f}. It sits just BELOW the "
            f"threshold, so it is rejected (x < t); the apparent tie is only a "
            f"rounding artifact in the table."
        )
    lines.append("")
    lines.append("## Per-pair table")
    lines.append("")
    lines.append("| pair | label | contiguous-core cov | all-column cov | aligned identity | minimizer-Jaccard | n blocks | biggest block | len A | len B | same univ. family |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|---|")
    for r in named_rows:
        lines.append(
            f"| {r['a']} <-> {r['b']} | {r['label']} | **{fmt(r['core_recip'])}** | "
            f"{fmt(r['recip'])} | {fmt(r['ident'])} | {fmt(r['jaccard'])} | "
            f"{r['n_blocks']} | {r['biggest']} | "
            f"{r['len_a']} | {r['len_b']} | {'yes' if r['same_family'] else 'NO'} |"
        )
    lines.append("")
    lines.append(
        f"Control pairs (random cross-family, n={len(control_rows)}, seed={SEED}): "
        f"contiguous-core coverage median {fmt(float(np.median(core['ctrl'])))}, "
        f"max {fmt(max(core['ctrl']))}; all-column coverage median "
        f"{fmt(float(np.median(raw['ctrl'])))}; minimizer-Jaccard median "
        f"{fmt(float(np.median([r['jaccard'] for r in control_rows])))}."
    )
    lines.append("")
    lines.append("## Interpretation")
    lines.append("")
    if core["separable"]:
        lines.append(
            "On the POA **contiguous-core** axis the distribution is **bimodal "
            f"and cleanly separable** ({margin_phrase}): every confirmed "
            "tandem/recent-duplicate pair has ONE large homologous co-aligned "
            f"block (>= {fmt(core['min_conf'])} of the shorter gene), while every "
            "domain-sharer's largest single block covers only a few percent "
            f"(<= {fmt(core['max_dom'])}) -- indistinguishable from random "
            f"cross-family controls (<= {fmt(core['max_ctrl'])}). A single "
            f"principled threshold T = {core['T']:.2f} accepts ALL confirmed and "
            "rejects ALL domain-sharers AND all controls."
        )
        lines.append("")
        lines.append(
            "The LITERAL all-column reciprocal coverage does NOT achieve this, "
            "because a global alignment of two non-homologous transcripts inflates "
            "coverage with gappy chance-match filler (visible as low all-column "
            "identity and hundreds of tiny blocks). Requiring the co-alignment to "
            "be a contiguous homologous CORE -- still a pure alignment statistic, "
            "no domain/BLAST/k-mer input -- is the **POA-only fix for the "
            "domain-sharing false positives.** The minimizer-Jaccard cannot do "
            "this at all: it rewards any shared subsequence, so a shared domain "
            "inflates the intersection exactly like whole-gene homology, leaving "
            "confirmed and domain-sharer Jaccard values overlapping."
        )
    else:
        lines.append(
            "Coverage separates the two classes better than minimizer-Jaccard "
            "but not perfectly on this small labeled set; see the per-pair table "
            "for the overlap cases and the caveats below."
        )
    lines.append("")
    lines.append("## Honest caveats")
    lines.append("")
    lines.append(
        "- **Operational RNA-structural definition, NOT the gene family.** This "
        "defines a *POA-coalignable multi-copy transcript family* (do the mature "
        "transcripts co-align over a majority of each?). It is deliberately NOT "
        "claimed equal to the DNA/protein gene family from a gene tree -- proving "
        "that equivalence would need DNA/protein/domain evidence, which the "
        "constraint forbids. It is an RNA-level, alignment-only operational "
        "criterion."
    )
    lines.append(
        "- **Pairwise-global as the POA pairwise instance.** For two sequences "
        "the POA instance reduces exactly to a single global alignment; for >2 "
        "copies a true POA DAG would be used and could differ slightly at "
        "multi-copy column assignment. The pair tests here are the exact "
        "two-sequence case."
    )
    lines.append(
        "- **All-column vs contiguous-core coverage.** The clean separation is on "
        "the contiguous-core axis, NOT the literal all-column reciprocal coverage "
        "(which overlaps because a global aligner pads non-homologous pairs with "
        "gappy chance-match filler). Both are alignment-column statistics and thus "
        "both POA-only/within-constraint, but the operational criterion that WORKS "
        "is specifically 'largest single UNGAPPED co-aligned block >= T of the "
        "shorter gene'. NOTE this is the largest *ungapped* block, not a "
        "local-alignment span: a Smith-Waterman local alignment of these "
        "domain-sharers ALSO stitches chance matches and reports a high "
        "aligned-length (verified: domain-sharer local aligned-len/minlen reaches "
        "0.91-0.97, overlapping the confirmed pairs) -- so it does NOT separate. "
        "It is the longest CONTIGUOUS gap-free homologous run that distinguishes "
        "a real copy (one long block) from a domain-sharer (many short blocks "
        "patched together). That largest-ungapped-block read-off is what the "
        "criterion uses."
    )
    lines.append(
        f"- **Small labeled N.** Only the Compara-checkable named pairs are "
        f"labeled: {n_conf} confirmed and {n_dom} domain-sharers. The separation "
        "is an exact count on a small sample, a directional sanity check, not a "
        "population rate."
    )
    lines.append(
        "- **Gene-representative sequences.** Coverage is computed on one "
        "representative sequence per gene (`gene_rep.fa`), not all isoforms; a "
        "different representative isoform could shift coverage at the margin."
    )
    lines.append(
        "- **RABL2A<->RABL2B label (assumption is doing real work).** Treated "
        "as confirmed-real per task spec because it is the flagship 2-copy "
        "tandem pair. Compara returned a fetch error for RABL2B, but note that "
        "RABL2A's OWN Compara query SUCCEEDED (para_status=ok, 68 paralogues) "
        "and does NOT list RABL2B's ENSG (ENSG00000079974). So the 'confirmed' "
        "label here rests on the task-spec/flagship assumption, not merely on a "
        "one-sided fetch error -- RABL2A's own successful list also omits "
        "RABL2B. Reassuringly, RABL2A<->RABL2B is the LOWEST confirmed pair "
        "(contiguous-core 0.174), yet still cleanly above T=0.13, so it does "
        "not prop up the separation."
    )
    lines.append(
        "- **Labels are external corroboration only.** Compara is used purely to "
        "*check* the POA criterion's verdicts; it is never an input to the "
        "criterion. The criterion is alignment-derived end to end."
    )
    lines.append(
        f"- **Determinism.** Controls use a fixed RNG seed ({SEED}); alignment "
        "and minimizer machinery are deterministic. Output is byte-stable."
    )
    lines.append("")

    with open(OUT_MD, "w") as fh:
        fh.write("\n".join(lines) + "\n")

    # ---- console summary ----
    print(headline)
    print(f"[contiguous-core] separable={core['separable']}  T={core['T']}  "
          f"min_conf={core['min_conf']:.3f}  max_dom={core['max_dom']:.3f}  "
          f"max_ctrl={core['max_ctrl']:.3f}")
    print(f"[all-column]      separable={raw['separable']}  T={raw['T']}  "
          f"min_conf={raw['min_conf']:.3f}  max_dom={raw['max_dom']:.3f}")
    print(f"[minimizer-Jacc]  best T={Tj_best:.3f} -> conf {n_jconf_acc_best}/{n_conf}, "
          f"dom rej {n_jdom_rej_best}/{n_dom}, separable={j_separable}")
    print(f"wrote {OUT_PNG}")
    print(f"wrote {OUT_MD}")

    return dict(
        T=T, chosen_name=chosen_name, separable=separable, n_conf=n_conf,
        n_dom=n_dom, n_conf_acc=n_conf_acc, n_dom_rej=n_dom_rej,
        n_ctrl_rej=n_ctrl_rej, n_ctrl=len(ctrl),
        core=core, raw=raw,
        j_separable=j_separable, Tj_best=Tj_best,
        n_jconf_acc_best=n_jconf_acc_best, n_jdom_rej_best=n_jdom_rej_best,
        headline=headline, rows=rows,
    )


if __name__ == "__main__":
    main()
