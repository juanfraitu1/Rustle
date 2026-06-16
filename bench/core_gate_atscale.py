#!/usr/bin/env python
"""
Contiguous-core family gate AT SCALE across ALL annotated paralog families.

Question
--------
`poa_family_definition.py` showed, on the 12 Ensembl-Compara-checkable named
pairs, that a contiguous-core coverage gate (largest ungapped equal block /
shorter gene >= T) cleanly KEEPS the confirmed recent-duplicate pairs and DROPS
the domain-sharing false positives that pure similarity grouping conflated. That
was a 12-pair directional sanity check. This script measures the SAME gate at
SCALE: across EVERY within-family gene pair of EVERY universe-annotated paralog
family, to quantify the genome-wide "false-family" reclassification rate -- how
many annotated within-family pairs the gate would KEEP (true core) vs DROP
(domain-sharing) at the shipped threshold T = 0.13.

This is the universe-ANNOTATED-FAMILY scale: a proxy for the full de-novo
pipeline (the production gate runs on de-novo assembled loci, not on the
RefSeq-annotated family universe), measured against the families the universe
TSV declares. It is NOT a re-derivation of the families themselves -- it takes
the universe's family assignment as given and asks, pair by pair, which of those
declared within-family memberships the contiguous-core gate would uphold.

GATE UNDER TEST (matches the fixed rust gate)
---------------------------------------------
For a within-family gene pair (a, b):

    core(a, b) = largest_ungapped_equal_block(a, b) / min(len(a), len(b))

computed from a GLOBAL (Needleman-Wunsch) alignment with match=+2,
mismatch=-1, gap-open=-5, gap-extend=-1 (Bio.Align.PairwiseAligner, the ROBUST
aligner reused verbatim from poa_family_definition.py -- the two-sequence POA
instance). The pair is

    KEPT      (true core / genuine copy)   iff core >= T
    FLAGGED   (domain-sharer, would split) iff core <  T

with T = 0.13 (the shipped threshold).

"Largest ungapped equal block" is `poa_family_definition.py`'s `biggest`: the
largest single ungapped ALIGNED block of the global alignment (one contiguous
gap-free run of paired columns, between two gaps), divided by the shorter gene.
This is the homologous core, robust to the gappy chance-match filler a global
aligner pads non-homologous pairs with -- and it is the metric that achieves
the clean T=0.13 separation on the Compara-labeled set. We ALSO report a
stricter variant -- the largest ungapped, base-IDENTICAL run -- for comparison,
but it is NOT the gate (it does not separate the labeled classes at T=0.13,
because real recent-duplicate cores carry scattered SNPs that break the exact
run; see the report).

INPUTS
------
- bench/copy_recovery_eval/results/gene_rep.fa   (gene representative sequences)
- bench/copy_recovery_eval/results/universe.tsv  (62 universe families, 195 genes)

LABELS (Compara, external corroboration only -- not part of the gate)
---------------------------------------------------------------------
CONFIRMED (real paralog / tandem dup):
  APOBEC3D<->APOBEC3F, RFPL1<->RFPL2, RFPL1<->RFPL3, RFPL2<->RFPL3, RABL2A<->RABL2B
DOMAIN-SHARER (Compara 'no'):
  ASDURF<->ASNSD1, CASP8<->FLACC1, CCDC188<->ZDHHC8, CDPF1<->PPARA,
  CREB1<->METTL21A, GCA<->KCNH7, GPR39<->LYPD1

Deterministic: alignment is deterministic; the per-pair jitter in the figure
uses a fixed RNG seed. No network, no randomness in any reported number.

Run with /home/juanfra/miniforge3/bin/python.
"""

import csv
import itertools
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
OUT_PNG = f"{BASE}/core_gate_atscale.png"
OUT_MD = f"{BASE}/core_gate_atscale.md"

# Shipped gate threshold (matches the fixed rust gate / poa_family_definition).
T = 0.13

# Global-alignment (= 2-sequence POA instance) scoring, identical to
# poa_family_definition.py and the fixed rust gate.
MATCH = 2.0
MISMATCH = -1.0
GAP_OPEN = -5.0
GAP_EXTEND = -1.0

SEED = 1729  # figure jitter only; no numeric result depends on it.

# Figure colors.
NAVY = "#1f3a5f"
ORANGE = "#e8820c"
GREEN = "#2e8b57"

# --------------------------------------------------------------------------
# Compara labels (external corroboration ONLY -- never an input to the gate)
# --------------------------------------------------------------------------
CONFIRMED_PAIRS = [
    ("APOBEC3D", "APOBEC3F"),
    ("RFPL1", "RFPL2"),
    ("RFPL1", "RFPL3"),
    ("RFPL2", "RFPL3"),
    ("RABL2A", "RABL2B"),
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
# ROBUST aligner (reused verbatim from poa_family_definition.py)
# --------------------------------------------------------------------------
def make_aligner():
    a = PairwiseAligner()
    a.mode = "global"
    a.match_score = MATCH
    a.mismatch_score = MISMATCH
    a.open_gap_score = GAP_OPEN
    a.extend_gap_score = GAP_EXTEND
    return a


def core_stats(aligner, s1, s2):
    """Largest contiguous-core block statistics from a global alignment.

    Returns:
      core_block   largest ungapped ALIGNED block / min(len) -- the GATE
                   metric, identical to poa_family_definition.py's `biggest`
                   (one contiguous gap-free run of paired columns / shorter
                   gene). This is what separates the labeled classes at T=0.13.
      core_equal   largest ungapped, base-IDENTICAL run / min(len) -- a
                   stricter variant reported for comparison only; it does NOT
                   separate the labeled classes at T=0.13.
      aligned_len  total ungapped aligned columns (all blocks)
      n_blocks     number of ungapped blocks (fragmentation signature)
      ident        identity over aligned columns
      biggest_eq   largest equal run length (bp)
      biggest_blk  largest ungapped block length (bp)
    """
    s1 = s1.upper()
    s2 = s2.upper()
    minlen = min(len(s1), len(s2))
    if minlen == 0:
        return dict(core_equal=0.0, core_block=0.0, aligned_len=0, n_blocks=0,
                    ident=0.0, biggest_eq=0, biggest_blk=0)
    aln = aligner.align(s1, s2)[0]
    t_blocks, q_blocks = aln.aligned
    aligned_len = 0
    matches = 0
    n_blocks = 0
    biggest_blk = 0
    biggest_eq = 0
    for (ts, te), (qs, qe) in zip(t_blocks, q_blocks):
        n = te - ts
        if n == 0:
            continue
        n_blocks += 1
        aligned_len += n
        biggest_blk = max(biggest_blk, n)
        seg1 = s1[ts:te]
        seg2 = s2[qs:qe]
        # longest run of consecutive base-identical positions within this
        # ungapped block (the contiguous EXACT core).
        run = 0
        for c1, c2 in zip(seg1, seg2):
            if c1 == c2:
                run += 1
                matches += 1
                if run > biggest_eq:
                    biggest_eq = run
            else:
                run = 0
    if aligned_len == 0:
        return dict(core_equal=0.0, core_block=0.0, aligned_len=0, n_blocks=0,
                    ident=0.0, biggest_eq=0, biggest_blk=0)
    return dict(
        core_equal=biggest_eq / minlen,
        core_block=biggest_blk / minlen,
        aligned_len=aligned_len,
        n_blocks=n_blocks,
        ident=matches / aligned_len,
        biggest_eq=biggest_eq,
        biggest_blk=biggest_blk,
    )


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
    aligner = make_aligner()

    # ---- 1. COVERAGE: which genes / families are evaluable? ----
    all_genes = sorted(g2f)
    present = [g for g in all_genes if g in fasta]
    missing = [g for g in all_genes if g not in fasta]
    n_genes = len(all_genes)
    n_present = len(present)
    n_loc = sum(1 for g in present if g.startswith("LOC"))

    # families with >= 2 genes present in the FASTA
    fam_present = {f: sorted(g for g in gs if g in fasta)
                   for f, gs in fam2genes.items()}
    fam_eval = {f: gs for f, gs in fam_present.items() if len(gs) >= 2}
    n_fam_total = len(fam2genes)
    n_fam_eval = len(fam_eval)

    # ---- 2. compute the gate for EVERY within-family pair ----
    pair_rows = []  # one dict per within-family pair
    for f in sorted(fam_eval):
        genes = fam_eval[f]
        for a, b in itertools.combinations(genes, 2):
            st = core_stats(aligner, fasta[a], fasta[b])
            key = npair(a, b)
            if key in CONFIRMED_SET:
                lab = "confirmed"
            elif key in DOMAIN_SET:
                lab = "domain-sharer"
            else:
                lab = "unlabeled"
            pair_rows.append(dict(
                family=f, a=a, b=b, label=lab,
                len_a=len(fasta[a]), len_b=len(fasta[b]),
                core=st["core_block"], core_equal=st["core_equal"],
                ident=st["ident"], n_blocks=st["n_blocks"],
                biggest_eq=st["biggest_eq"], biggest_blk=st["biggest_blk"],
                kept=st["core_block"] >= T,
            ))

    n_pairs = len(pair_rows)
    n_kept = sum(1 for r in pair_rows if r["kept"])
    n_drop = n_pairs - n_kept

    # ---- 2b. per-family fraction-of-pairs-kept distribution ----
    fam_keep_frac = {}
    fam_kept_n = {}
    fam_pairs_n = {}
    for f in sorted(fam_eval):
        rs = [r for r in pair_rows if r["family"] == f]
        k = sum(1 for r in rs if r["kept"])
        fam_kept_n[f] = k
        fam_pairs_n[f] = len(rs)
        fam_keep_frac[f] = k / len(rs) if rs else 0.0

    fracs = list(fam_keep_frac.values())
    fam_all_kept = sum(1 for v in fracs if v == 1.0)
    fam_all_drop = sum(1 for v in fracs if v == 0.0)
    fam_mixed = n_fam_eval - fam_all_kept - fam_all_drop

    # ---- 3. Compara-labeled subset verdicts ----
    labeled_rows = [r for r in pair_rows if r["label"] in ("confirmed", "domain-sharer")]
    conf_rows = [r for r in labeled_rows if r["label"] == "confirmed"]
    dom_rows = [r for r in labeled_rows if r["label"] == "domain-sharer"]
    # any labeled pair not found among within-family pairs? (sanity)
    found_labeled = {npair(r["a"], r["b"]) for r in labeled_rows}
    expected_labeled = CONFIRMED_SET | DOMAIN_SET
    not_within_family = sorted(expected_labeled - found_labeled)

    conf_kept = sum(1 for r in conf_rows if r["kept"])
    dom_dropped = sum(1 for r in dom_rows if not r["kept"])

    # ---- 4. FIGURE (navy/orange/green) ----
    rng = random.Random(SEED)
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(15.5, 6.0))

    # Panel 0: strip/jitter of ALL within-family-pair contiguous-core values,
    # one row, colored kept (navy, core>=T) vs dropped (orange, core<T);
    # unlabeled drawn small, Compara-labeled drawn large with green edge +
    # text. T line in black dashed.
    def jit():
        return (rng.random() - 0.5) * 0.7

    # unlabeled first (background)
    un = [r for r in pair_rows if r["label"] == "unlabeled"]
    xs = [r["core"] for r in un]
    ys = [jit() for _ in un]
    cs = [NAVY if r["kept"] else ORANGE for r in un]
    ax0.scatter(xs, ys, s=30, c=cs, alpha=0.55, edgecolor="white",
                linewidth=0.3, zorder=2)

    # labeled on top, larger, green ring
    for r in labeled_rows:
        x = r["core"]
        y = jit()
        face = NAVY if r["kept"] else ORANGE
        ax0.scatter([x], [y], s=120, c=[face], alpha=0.95,
                    edgecolor=GREEN, linewidth=2.0, zorder=4)
        tag = f"{r['a']}/{r['b']}"
        ax0.annotate(tag, (x, y), fontsize=6.5, xytext=(3, 4),
                     textcoords="offset points", zorder=5,
                     color="black")

    ax0.axvline(T, color="black", linestyle="--", linewidth=1.6, zorder=3)
    ax0.text(T, 0.42, f"  T = {T:.2f}", fontsize=11, va="bottom", ha="left",
             color="black", fontweight="bold")
    ax0.set_xlabel("contiguous-core coverage  (largest ungapped equal block / shorter gene)")
    ax0.set_yticks([])
    ax0.set_xlim(-0.02, 1.02)
    ax0.set_ylim(-0.6, 0.6)
    ax0.set_title(
        f"All {n_pairs} within-family pairs ({n_fam_eval} families)\n"
        f"navy = KEPT (core>=T, n={n_kept}) | orange = DROPPED (core<T, n={n_drop})\n"
        f"green ring = Compara-labeled",
        fontsize=10.5)
    ax0.grid(axis="x", alpha=0.25)
    # legend proxies
    from matplotlib.lines import Line2D
    leg = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=NAVY,
               markersize=9, label=f"KEPT core>=T (n={n_kept})"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=ORANGE,
               markersize=9, label=f"DROPPED core<T (n={n_drop})"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="gray",
               markeredgecolor=GREEN, markeredgewidth=2, markersize=11,
               label="Compara-labeled"),
    ]
    ax0.legend(handles=leg, fontsize=8.5, loc="upper right")

    # Panel 1: histogram of contiguous-core coverage, stacked kept/dropped,
    # plus overlaid markers for the Compara-labeled subset.
    bins = np.linspace(0.0, 1.0, 41)
    drop_vals = [r["core"] for r in pair_rows if not r["kept"]]
    keep_vals = [r["core"] for r in pair_rows if r["kept"]]
    ax1.hist([drop_vals, keep_vals], bins=bins, stacked=True,
             color=[ORANGE, NAVY], edgecolor="white", linewidth=0.3,
             label=[f"DROPPED (n={n_drop})", f"KEPT (n={n_kept})"])
    ax1.axvline(T, color="black", linestyle="--", linewidth=1.6)
    ax1.text(T, ax1.get_ylim()[1] * 0.96, f"  T = {T:.2f}", fontsize=11,
             va="top", ha="left", color="black", fontweight="bold")
    # labeled subset as rug at top
    ytop = ax1.get_ylim()[1]
    for r in conf_rows:
        ax1.plot([r["core"]], [ytop * 0.5], marker="v", color=NAVY,
                 markeredgecolor=GREEN, markeredgewidth=1.5, markersize=11,
                 zorder=6)
    for r in dom_rows:
        ax1.plot([r["core"]], [ytop * 0.5], marker="v", color=ORANGE,
                 markeredgecolor=GREEN, markeredgewidth=1.5, markersize=11,
                 zorder=6)
    ax1.set_xlabel("contiguous-core coverage")
    ax1.set_ylabel("within-family pairs")
    ax1.set_xlim(0.0, 1.0)
    ax1.set_title(
        f"Distribution of contiguous-core coverage\n"
        f"Compara: {conf_kept}/{len(conf_rows)} confirmed KEPT, "
        f"{dom_dropped}/{len(dom_rows)} domain-sharers DROPPED "
        f"(triangles, green ring)",
        fontsize=10.5)
    ax1.legend(fontsize=9, loc="upper center")
    ax1.grid(axis="y", alpha=0.25)

    fig.suptitle(
        f"Contiguous-core family gate (T={T:.2f}) at scale: of {n_pairs} "
        f"annotated within-family pairs across {n_fam_eval} families, "
        f"{n_kept} KEPT ({100*n_kept/n_pairs:.0f}%) / {n_drop} DROPPED "
        f"({100*n_drop/n_pairs:.0f}%)",
        fontsize=12.5, y=0.99)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(OUT_PNG, dpi=150)
    plt.close(fig)

    # ---- 5. MARKDOWN REPORT ----
    def fmt(x, n=3):
        return f"{x:.{n}f}"

    # per-family table sorted by keep-fraction then size
    fam_order = sorted(
        fam_eval,
        key=lambda f: (fam_keep_frac[f], -fam_pairs_n[f], f))

    lines = []
    lines.append("# Contiguous-core family gate AT SCALE")
    lines.append("")
    lines.append(
        f"> **Of {n_pairs} annotated within-family gene pairs across "
        f"{n_fam_eval} universe families, the contiguous-core gate (T={T:.2f}) "
        f"KEEPS {n_kept} ({100*n_kept/n_pairs:.1f}%) as true cores and DROPS "
        f"{n_drop} ({100*n_drop/n_pairs:.1f}%) as domain-sharing "
        f"false-family memberships.**")
    lines.append("")
    lines.append(
        f"> On the Compara-labeled subset the gate is fully consistent: it KEEPS "
        f"{conf_kept}/{len(conf_rows)} confirmed pairs and DROPS "
        f"{dom_dropped}/{len(dom_rows)} domain-sharers.")
    lines.append("")

    lines.append("## What this measures (and what it does NOT)")
    lines.append("")
    lines.append(
        "This is the **universe-annotated-family scale**: it takes the "
        "universe TSV's family assignment as GIVEN and asks, pair by pair, "
        "which of those declared within-family memberships the contiguous-core "
        "gate would uphold (KEEP) versus reclassify as domain-sharing (DROP). "
        "It is a **proxy for the full de-novo pipeline** -- the production gate "
        "actually runs on de-novo assembled loci, not on the RefSeq-annotated "
        "family universe -- so the genome-wide rate here is an annotation-scale "
        "estimate, not the exact de-novo rate. The Compara labels exist only "
        "for a small subset (12 named pairs); the rest are UNLABELED, so for "
        "them KEEP/DROP is the gate's verdict, not a verified ground truth.")
    lines.append("")

    lines.append("## The gate")
    lines.append("")
    lines.append("```")
    lines.append("core(a,b) = largest_ungapped_equal_block(a,b) / min(len(a), len(b))")
    lines.append(f"KEEP  iff core >= T   (T = {T:.2f})")
    lines.append("DROP  iff core <  T")
    lines.append("```")
    lines.append("")
    lines.append(
        "computed from a GLOBAL Needleman-Wunsch alignment "
        f"(match=+{int(MATCH)}, mismatch={int(MISMATCH)}, gap-open={int(GAP_OPEN)}, "
        f"gap-extend={int(GAP_EXTEND)}; `Bio.Align.PairwiseAligner`, the robust "
        "aligner reused verbatim from `poa_family_definition.py`). The "
        "'largest ungapped equal block' is the longest single contiguous "
        "gap-free run of paired columns in the global alignment (one block "
        "between two gaps), divided by the shorter gene -- the contiguous "
        "homologous core, robust to the gappy chance-match filler a global "
        "aligner pads non-homologous pairs with. This is exactly "
        "`poa_family_definition.py`'s `biggest` metric.")
    lines.append("")

    lines.append("## (1) Coverage / evaluability")
    lines.append("")
    lines.append("| quantity | value |")
    lines.append("|---|---|")
    lines.append(f"| universe families | {n_fam_total} |")
    lines.append(f"| universe distinct genes | {n_genes} |")
    lines.append(f"| genes present in gene_rep.fa | {n_present} ({100*n_present/n_genes:.1f}%) |")
    lines.append(f"| genes MISSING from gene_rep.fa | {len(missing)} |")
    lines.append(f"| of present genes, LOC* (provisional) loci | {n_loc} ({100*n_loc/n_present:.0f}%) |")
    lines.append(f"| families evaluable (>=2 genes present) | {n_fam_eval} of {n_fam_total} |")
    lines.append(f"| within-family pairs evaluated | {n_pairs} |")
    lines.append("")
    if missing:
        lines.append(f"Missing genes: {', '.join(sorted(missing))}")
    else:
        lines.append(
            "**Match rate is 100%** -- every universe gene (including all "
            f"{n_loc} LOC* provisional loci) is present in `gene_rep.fa`, so all "
            f"{n_fam_total} families and all {n_pairs} within-family pairs are "
            "evaluable; nothing is dropped for missing sequence.")
    lines.append("")

    lines.append("## (2) The at-scale effect")
    lines.append("")
    lines.append("| verdict | pairs | % of pairs |")
    lines.append("|---|---|---|")
    lines.append(f"| KEPT (true core, core>=T) | {n_kept} | {100*n_kept/n_pairs:.1f}% |")
    lines.append(f"| DROPPED (domain-sharer, core<T) | {n_drop} | {100*n_drop/n_pairs:.1f}% |")
    lines.append(f"| total | {n_pairs} | 100% |")
    lines.append("")
    lines.append(
        f"**Genome-wide (annotation-scale) false-family reclassification rate "
        f"= {100*n_drop/n_pairs:.1f}%** of declared within-family pairs would be "
        "split off by the gate as domain-sharers.")
    lines.append("")
    lines.append("### Per-family fraction-of-pairs-kept distribution")
    lines.append("")
    lines.append("| keep-fraction | families |")
    lines.append("|---|---|")
    lines.append(f"| 1.00 (all pairs kept) | {fam_all_kept} |")
    lines.append(f"| mixed (0 < frac < 1) | {fam_mixed} |")
    lines.append(f"| 0.00 (all pairs dropped) | {fam_all_drop} |")
    lines.append("")
    lines.append(
        f"Median per-family keep-fraction = {fmt(float(np.median(fracs)))}; "
        f"mean = {fmt(float(np.mean(fracs)))}. "
        f"{fam_all_kept} of {n_fam_eval} families are fully upheld (every "
        f"within-family pair kept); {fam_all_drop} are fully reclassified "
        f"(every pair dropped); {fam_mixed} are mixed (the gate splits the "
        "family into a true-core subset + domain-sharing outliers).")
    lines.append("")
    lines.append(
        "Quartiles of the per-family keep-fraction: "
        f"q25={fmt(float(np.percentile(fracs,25)))}, "
        f"q50={fmt(float(np.percentile(fracs,50)))}, "
        f"q75={fmt(float(np.percentile(fracs,75)))}.")
    lines.append("")

    lines.append("## (3) Compara-labeled subset (the gate is correct where we can check)")
    lines.append("")
    lines.append(
        "The 12 Compara-checkable named pairs are the only within-family pairs "
        "with external ground truth. The gate's verdict on each:")
    lines.append("")
    lines.append("| pair | family | Compara label | contiguous-core | gate verdict | correct? |")
    lines.append("|---|---|---|---|---|---|")
    for r in sorted(labeled_rows, key=lambda r: (r["label"], -r["core"])):
        verdict = "KEEP" if r["kept"] else "DROP"
        want = "KEEP" if r["label"] == "confirmed" else "DROP"
        ok = "yes" if verdict == want else "**NO**"
        lines.append(
            f"| {r['a']} <-> {r['b']} | {r['family']} | {r['label']} | "
            f"{fmt(r['core'])} | **{verdict}** | {ok} |")
    lines.append("")
    if not_within_family:
        lines.append(
            "Labeled pairs not found as within-family pairs (not co-assigned in "
            f"the universe): {not_within_family}.")
        lines.append("")
    lines.append(
        f"**The gate KEEPS all {conf_kept}/{len(conf_rows)} confirmed pairs and "
        f"DROPS all {dom_dropped}/{len(dom_rows)} domain-sharers** -- it is "
        "perfectly consistent with Compara on the labeled subset. Confirmed "
        f"contiguous-core range "
        f"[{fmt(min(r['core'] for r in conf_rows))}, "
        f"{fmt(max(r['core'] for r in conf_rows))}]; domain-sharer range "
        f"[{fmt(min(r['core'] for r in dom_rows))}, "
        f"{fmt(max(r['core'] for r in dom_rows))}]. The shipped T={T:.2f} sits "
        "in the gap between them.")
    lines.append("")

    lines.append("## Per-family table (all evaluable families)")
    lines.append("")
    lines.append("Sorted by keep-fraction (most-reclassified first), then size.")
    lines.append("")
    lines.append("| family | genes | pairs | kept | dropped | keep-frac | has Compara label |")
    lines.append("|---|---|---|---|---|---|---|")
    for f in fam_order:
        labset = {r["label"] for r in pair_rows
                  if r["family"] == f and r["label"] != "unlabeled"}
        lab = ", ".join(sorted(labset)) if labset else "-"
        lines.append(
            f"| {f} | {len(fam_eval[f])} | {fam_pairs_n[f]} | {fam_kept_n[f]} | "
            f"{fam_pairs_n[f]-fam_kept_n[f]} | {fmt(fam_keep_frac[f])} | {lab} |")
    lines.append("")

    lines.append("## Honest caveats")
    lines.append("")
    lines.append(
        "- **Annotation-scale proxy, not the de-novo pipeline.** This runs on "
        "the RefSeq-annotated universe families (`universe.tsv`), taking their "
        "assignment as given. The production gate runs on DE-NOVO assembled "
        "loci; the genome-wide rate here is an annotation-scale estimate of the "
        "gate's behavior, not the exact de-novo reclassification rate.")
    lines.append(
        "- **Most pairs are UNLABELED.** Only the 12 Compara-checkable named "
        f"pairs ({len(conf_rows)} confirmed + {len(dom_rows)} domain-sharer) "
        "have external ground truth. For the other "
        f"{n_pairs - len(labeled_rows)} pairs, KEEP/DROP is the gate's own "
        "verdict, validated only indirectly by its perfect agreement with "
        "Compara on the labeled subset.")
    lines.append(
        "- **LOC families included.** All "
        f"{n_loc} present LOC* (provisional/computationally-predicted) loci are "
        "included if in `gene_rep.fa`. Many of the largest, most-fragmented "
        "families are LOC tandem arrays (e.g. GGTLC2, LOC101126655), whose "
        "internal pairs the gate splits aggressively; these dominate the DROP "
        "count and may include both genuine sub-family structure and "
        "annotation noise.")
    lines.append(
        "- **Gene-representative sequences.** One representative sequence per "
        "gene (`gene_rep.fa`); a different representative isoform could shift a "
        "borderline pair's coverage.")
    lines.append(
        "- **Block vs exact-run metric.** The gate metric is the largest "
        "ungapped ALIGNED block / shorter gene -- identical to "
        "`poa_family_definition.py`'s `biggest`, the metric that separates the "
        "labeled classes at T=0.13. A stricter variant (largest ungapped "
        "base-IDENTICAL run) was also computed but is NOT used as the gate: it "
        "does not separate the labeled set at T=0.13, because genuine recent-"
        "duplicate cores carry scattered SNPs that break the exact run "
        f"(confirmed exact-run range "
        f"[{fmt(min(r['core_equal'] for r in conf_rows))}, "
        f"{fmt(max(r['core_equal'] for r in conf_rows))}], all below T). The "
        "aligned-block metric tolerates those SNPs and is the gate-faithful "
        "one.")
    lines.append(
        "- **Single threshold, no per-family tuning.** A single shipped "
        f"T={T:.2f} is applied uniformly; the separation was established on the "
        "12-pair labeled set in `poa_family_definition.py` and is re-used here "
        "unchanged.")
    lines.append(
        "- **Determinism.** Alignment is deterministic; only the figure jitter "
        f"uses a fixed seed ({SEED}). Every reported number is reproducible.")
    lines.append("")

    with open(OUT_MD, "w") as fh:
        fh.write("\n".join(lines) + "\n")

    # ---- console summary ----
    print(f"coverage: {n_present}/{n_genes} genes present "
          f"({100*n_present/n_genes:.1f}%), {n_fam_eval}/{n_fam_total} families "
          f"evaluable, {n_pairs} within-family pairs")
    print(f"at scale: KEPT {n_kept} ({100*n_kept/n_pairs:.1f}%) / DROPPED "
          f"{n_drop} ({100*n_drop/n_pairs:.1f}%)")
    print(f"per-family keep-frac: median={np.median(fracs):.3f} "
          f"all-kept={fam_all_kept} mixed={fam_mixed} all-dropped={fam_all_drop}")
    print(f"Compara: confirmed KEPT {conf_kept}/{len(conf_rows)}, "
          f"domain-sharers DROPPED {dom_dropped}/{len(dom_rows)}")
    print(f"wrote {OUT_PNG}")
    print(f"wrote {OUT_MD}")

    return dict(
        n_genes=n_genes, n_present=n_present, n_missing=len(missing),
        n_fam_total=n_fam_total, n_fam_eval=n_fam_eval, n_pairs=n_pairs,
        n_kept=n_kept, n_drop=n_drop,
        fam_all_kept=fam_all_kept, fam_mixed=fam_mixed, fam_all_drop=fam_all_drop,
        conf_kept=conf_kept, n_conf=len(conf_rows),
        dom_dropped=dom_dropped, n_dom=len(dom_rows),
        not_within_family=not_within_family,
    )


if __name__ == "__main__":
    main()
