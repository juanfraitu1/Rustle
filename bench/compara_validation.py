#!/usr/bin/env python
"""
compara_validation.py
=====================

NON-CIRCULAR validation of rustle's / the universe's paralog families against an
INDEPENDENT, human-curated, gene-tree-based truth: Ensembl Compara paralogy.

WHY THIS IS NON-CIRCULAR
------------------------
The earlier `family_detection_validation.py` validates one SEQUENCE-SIMILARITY
clustering (minimizer-Jaccard) against another (the minimap2-built
`universe.tsv`). Both reward shared subsequence, so that concordance is partly
methodological, not biological. This script swaps the truth set for Ensembl
Compara paralogy, which is built from MULTIPLE-SEQUENCE-ALIGNMENT gene trees +
species-tree reconciliation (`compara_paralog_relation.json`, downloaded from
the Ensembl REST API). Compara never sees our minimap2 alignments, so a "yes,
these are paralogs" from Compara is an INDEPENDENT confirmation.

KEY FRAMING (state it loudly)
-----------------------------
Ensembl Compara paralogy is GENE-TREE-BASED and COARSER than our families. A
Compara query for RABL2A returns the whole RAB / GTPase superfamily (~68
paralogues across the genome), whereas our family is the 2-copy tandem
RABL2A / RABL2B. So:

  * OUR families should be SUBSETS of Compara groups (we split ancient
    superfamilies into recent copy-clusters; that is correct granularity, not
    error).
  * Therefore PRECISION is the headline metric: of the gene PAIRS we put in the
    SAME family, what fraction does Compara independently also call paralogs?
  * RECALL vs Compara is intentionally NOT the metric: Compara lumping ancient
    paralogs we correctly separate is expected granularity. Reporting it as
    "missed" would be wrong.

WHAT IS COMPUTED (restricted to the mapped NAMED genes)
-------------------------------------------------------
  (1) PRECISION vs Compara of the UNIVERSE families: of all within-universe-
      family NAMED gene pairs (both genes mapped to an ENSG), the fraction
      Compara independently calls paralogs. De-circularizes the earlier work.
  (2) PRECISION vs Compara of RUSTLE's minimizer-Jaccard grouping at 0.30 and
      0.06 (re-using the EXACT criterion from family_detection_validation.py:
      canonical k=15/w=10 minimizers, blake2b-64, union-find), over the SAME
      named genes.
  (3) GRANULARITY: how many DISTINCT universe families fall inside ONE Compara
      paralog group (Compara is coarser / merges them) -- an observation, not an
      error. And the reverse: any universe family that SPANS >1 Compara group
      (a real over-merge by us).
  (4) COVERAGE up front: N named genes mapped, N families with >=2 mapped genes
      (the evaluable set), N evaluable within-family pairs.

HONEST CAVEATS (printed + in the .md)
-------------------------------------
  * HUMAN-AS-PROXY: the data are GORILLA genes; Compara paralogy is queried on
    the HUMAN ortholog's ENSG. Paralogy relationships are deeply conserved
    across great apes (the gorilla/human split is ~10 Mya, far younger than the
    duplications that created these families), but it is a proxy, not the gorilla
    gene tree.
  * COMPARA IS COARSER (superfamily-level): the central framing above. A
    NON-confirmation can mean EITHER a genuine universe false-grouping OR a
    Compara gap / different granularity; we name every non-confirmed pair so the
    reader can judge.
  * ONLY NAMED GENES MAPPED: 154 of 195 universe genes are LOC-only (no human
    symbol -> no ENSG -> not checkable here). Coverage is stated up front; this
    is a check on the NAMED-gene backbone of the families, not all 62 families.
  * SYMBOL -> ENSEMBL MAPPING: genes were resolved by symbol via the Ensembl
    REST lookup (see meta in the JSON); 40/41 got an ENSG, 1 lookup error
    (GP1BB), 4 paralogue-fetch errors. Mapping noise is possible.
  * SMALL N: only 10 universe families have >=2 mapped named genes -> 12
    within-family named pairs. This is a SMALL, HONEST sample; we report exact
    counts and do not over-generalize a percentage.

DETERMINISM: no RNG used at all (no sampling -- we evaluate the full mapped set).
The minimizer hash is blake2b (fixed). Output is byte-stable across runs.

Run:
  /home/juanfra/miniforge3/bin/python bench/compara_validation.py
"""

import os
import sys
import json
import csv
import hashlib
from collections import defaultdict
from itertools import combinations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --------------------------------------------------------------------------
# Fixed parameters (match family_detection_validation.py)
# --------------------------------------------------------------------------
K = 15                       # rustle minimizer k
W = 10                       # rustle minimizer w
T_RUSTLE = 0.30              # FAMILY_MERGE_JACCARD_DEFAULT
T_OPT = 0.06                 # data-optimal whole-gene threshold from prior work

BASE = "/mnt/c/Users/jfris/Desktop/Rustle/bench"
RESULTS = os.path.join(BASE, "copy_recovery_eval", "results")
UNIVERSE_TSV = os.path.join(RESULTS, "universe.tsv")
GENE_FA = os.path.join(RESULTS, "gene_rep.fa")
COMPARA_JSON = os.path.join(BASE, "compara_paralog_relation.json")

OUT_PNG = os.path.join(BASE, "compara_validation.png")
OUT_MD = os.path.join(BASE, "compara_validation.md")

NAVY = "#22313f"
ORANGE = "#e8590c"
GREEN = "#1e8449"

# --------------------------------------------------------------------------
# Minimizer machinery -- byte-for-byte the same scheme as
# family_detection_validation.py (canonical k-mer, blake2b-64, rustle windowing)
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
    uni = len(a | b)
    return len(a & b) / uni if uni else 0.0

# --------------------------------------------------------------------------
# Loaders
# --------------------------------------------------------------------------
def load_fasta(path):
    seqs = {}
    cur = None
    buf = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line[0] == ">":
                if cur is not None:
                    seqs[cur] = "".join(buf).encode()
                cur = line[1:].split()[0]
                if cur.startswith("gene-"):
                    cur = cur[len("gene-"):]
                buf = []
            else:
                buf.append(line)
        if cur is not None:
            seqs[cur] = "".join(buf).encode()
    return seqs

def load_universe(path):
    """gene_id -> family_id, and family_id -> set(gene_id)."""
    g2f = {}
    fams = defaultdict(set)
    with open(path) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            g = row["gene_id"].strip()
            f = row["family_id"].strip()
            if g.startswith("gene-"):
                g = g[len("gene-"):]
            g2f[g] = f
            fams[f].add(g)
    return g2f, fams


print(f"[load] reading inputs (deterministic, no RNG)", file=sys.stderr)
fasta = load_fasta(GENE_FA)
g2f, uni_fams = load_universe(UNIVERSE_TSV)
comp = json.load(open(COMPARA_JSON))
gene_info = comp["genes"]

# --------------------------------------------------------------------------
# Build the Compara relation among NAMED genes
# --------------------------------------------------------------------------
named = sorted(gene_info.keys())                       # 41 named genes attempted
ensg = {g: gene_info[g].get("ensg") for g in named}    # symbol -> ENSG (or None)
para = {g: set(gene_info[g].get("paralogues") or []) for g in named}

# A named gene is "mapped" iff it has an ENSG (so we can place it in Compara id
# space). It is "Compara-checkable" against a partner iff BOTH have an ENSG.
mapped_named = sorted(g for g in named if ensg.get(g))

# Reverse index: ENSG -> our gene symbol (for the granularity analysis).
ensg2sym = {ensg[g]: g for g in mapped_named}

def compara_paralog(a, b):
    """Symmetric: are a and b independently called paralogs by human Compara?
    True iff ENSG(b) appears in a's paralogue set OR ENSG(a) in b's set.
    Requires both mapped to an ENSG."""
    ea, eb = ensg.get(a), ensg.get(b)
    if not (ea and eb):
        return None  # un-checkable (missing ENSG)
    if eb in para.get(a, set()):
        return True
    if ea in para.get(b, set()):
        return True
    return False

# Sanity: reproduce the JSON's stated within-universe symmetric relation.
json_within = set(tuple(sorted(p)) for p in comp["within_universe_paralog_relation"])

# --------------------------------------------------------------------------
# (4) COVERAGE
# --------------------------------------------------------------------------
n_named = len(named)
n_named_mapped = len(mapped_named)                       # got ENSG
named_in_universe = [g for g in named if g in g2f]
# Universe families restricted to mapped named genes:
fam_named_mapped = {}                                    # family -> [mapped named genes]
for f, members in uni_fams.items():
    nm = sorted(g for g in members if g in ensg and ensg.get(g))
    if nm:
        fam_named_mapped[f] = nm
evaluable_fams = sorted(f for f, nm in fam_named_mapped.items() if len(nm) >= 2)
n_evaluable_fams = len(evaluable_fams)

# --------------------------------------------------------------------------
# (1) PRECISION vs Compara of the UNIVERSE families
#     Of all within-universe-family NAMED pairs (both mapped), fraction Compara
#     confirms as paralogs.
# --------------------------------------------------------------------------
uni_pairs = []          # (family, a, b, compara_bool)
for f in evaluable_fams:
    nm = fam_named_mapped[f]
    for a, b in combinations(nm, 2):
        uni_pairs.append((f, a, b, compara_paralog(a, b)))

uni_pairs_checkable = [p for p in uni_pairs if p[3] is not None]
uni_confirmed = [p for p in uni_pairs_checkable if p[3] is True]
uni_unconfirmed = [p for p in uni_pairs_checkable if p[3] is False]
uni_n = len(uni_pairs_checkable)
uni_conf_n = len(uni_confirmed)
uni_precision = (uni_conf_n / uni_n) if uni_n else float("nan")

print(f"[universe] within-family checkable pairs={uni_n} "
      f"confirmed={uni_conf_n} precision={uni_precision:.3f}", file=sys.stderr)
for f, a, b, ok in uni_pairs_checkable:
    print(f"    {'OK ' if ok else 'NO '} {f}: {a}<->{b}", file=sys.stderr)

# --------------------------------------------------------------------------
# (2) PRECISION vs Compara of RUSTLE's minimizer-Jaccard grouping
#     Cluster the MAPPED NAMED genes by union-find at T, then score every
#     within-predicted-cluster pair against Compara.
# --------------------------------------------------------------------------
# Restrict to mapped named genes that also have a sequence (all 41 do, but guard).
cluster_genes = [g for g in mapped_named if g in fasta]
missing_seq = [g for g in mapped_named if g not in fasta]
mins = {g: minimizers(fasta[g], K, W) for g in cluster_genes}
cluster_genes = [g for g in cluster_genes if mins[g]]   # drop empty (none expected)
ng = len(cluster_genes)
idx = {g: i for i, g in enumerate(cluster_genes)}

# Precompute pairwise Jaccard once (deterministic i<j walk).
pair_jac = {}
for i in range(ng):
    gi = cluster_genes[i]
    for j in range(i + 1, ng):
        gj = cluster_genes[j]
        pair_jac[(gi, gj)] = jaccard(mins[gi], mins[gj])

def cluster_at(threshold):
    parent = list(range(ng))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for i in range(ng):
        for j in range(i + 1, ng):
            if pair_jac[(cluster_genes[i], cluster_genes[j])] >= threshold:
                ra, rb = find(i), find(j)
                if ra != rb:
                    parent[ra] = rb
    clusters = defaultdict(list)
    for i in range(ng):
        clusters[find(i)].append(cluster_genes[i])
    # deterministic ordering
    return sorted((sorted(v) for v in clusters.values()), key=lambda c: c[0])

def rustle_precision(threshold):
    """Precision vs Compara of rustle's predicted clusters at `threshold`.
    Of every within-predicted-cluster NAMED pair (both mapped), fraction Compara
    confirms. Returns (precision, n_pairs, confirmed, pairs_list)."""
    clusters = cluster_at(threshold)
    pairs = []
    for c in clusters:
        if len(c) < 2:
            continue
        for a, b in combinations(sorted(c), 2):
            ok = compara_paralog(a, b)
            if ok is None:
                continue
            pairs.append((a, b, ok))
    conf = sum(1 for _a, _b, ok in pairs if ok)
    n = len(pairs)
    prec = (conf / n) if n else float("nan")
    return prec, n, conf, pairs, clusters

prec_030, n_030, conf_030, pairs_030, clusters_030 = rustle_precision(T_RUSTLE)
prec_006, n_006, conf_006, pairs_006, clusters_006 = rustle_precision(T_OPT)

print(f"[rustle@0.30] pairs={n_030} confirmed={conf_030} precision={prec_030 if n_030 else float('nan')}", file=sys.stderr)
print(f"[rustle@0.06] pairs={n_006} confirmed={conf_006} precision={prec_006 if n_006 else float('nan')}", file=sys.stderr)

def fmt_clusters(clusters):
    return [c for c in clusters if len(c) >= 2]

print(f"[rustle@0.30] non-trivial clusters: {fmt_clusters(clusters_030)}", file=sys.stderr)
print(f"[rustle@0.06] non-trivial clusters: {fmt_clusters(clusters_006)}", file=sys.stderr)

# --------------------------------------------------------------------------
# (3) GRANULARITY
#     a) How many DISTINCT universe families fall inside ONE Compara paralog
#        group (Compara coarser / merges them)?  -> observation, not error.
#     b) Reverse: any universe family that SPANS >1 Compara group? -> real
#        over-merge by us.
#
# We build Compara connected components over the MAPPED NAMED genes (union genes
# that are Compara paralogs of each other), then map universe families onto them.
# --------------------------------------------------------------------------
# Compara graph among mapped named genes.
cparent = {g: g for g in mapped_named}
def cfind(x):
    while cparent[x] != x:
        cparent[x] = cparent[cparent[x]]
        x = cparent[x]
    return x
for a, b in combinations(mapped_named, 2):
    if compara_paralog(a, b):
        ra, rb = cfind(a), cfind(b)
        if ra != rb:
            cparent[ra] = rb
compara_comp = defaultdict(list)
for g in mapped_named:
    compara_comp[cfind(g)].append(g)
compara_groups = [sorted(v) for v in compara_comp.values()]
# non-trivial Compara groups (>=2 of OUR mapped named genes)
compara_groups_nt = sorted((c for c in compara_groups if len(c) >= 2), key=lambda c: c[0])

# For each non-trivial Compara group, how many DISTINCT universe families do its
# member genes come from? (>1 => Compara merges multiple universe families =
# Compara is coarser.)
compara_merges_universe = []   # (compara_group_genes, set_of_universe_families)
for grp in compara_groups_nt:
    ufams = set(g2f.get(g) for g in grp)
    compara_merges_universe.append((grp, sorted(ufams)))
n_compara_merging_multi = sum(1 for _g, uf in compara_merges_universe if len(uf) > 1)

# Reverse: does any universe family's mapped named genes split across >1 Compara
# component? That would be a real over-merge by us. (Genes Compara does NOT link
# fall in their OWN singleton component, so a family with any non-confirmed
# internal pair will "span" multiple components -- this is exactly the precision
# story, re-expressed at component level.)
universe_spanning = []
for f in evaluable_fams:
    nm = fam_named_mapped[f]
    comps = set(cfind(g) for g in nm)
    if len(comps) > 1:
        universe_spanning.append((f, nm, len(comps)))

print(f"[granularity] non-trivial Compara groups (>=2 of our mapped named)="
      f"{len(compara_groups_nt)}; Compara groups merging >1 universe family="
      f"{n_compara_merging_multi}; universe families spanning >1 Compara comp="
      f"{len(universe_spanning)}", file=sys.stderr)
for grp, uf in compara_merges_universe:
    print(f"    COMPARA-GROUP {grp}  universe_families={uf}", file=sys.stderr)

# --------------------------------------------------------------------------
# FIGURE
# --------------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.4),
                               gridspec_kw={"width_ratios": [1.35, 1.0]})

# Left: precision-vs-Compara bars for the three groupings.
labels = ["universe\nfamilies",
          f"rustle\nJaccard@{T_RUSTLE:g}",
          f"rustle\nJaccard@{T_OPT:g}"]
precs = [uni_precision,
         prec_030 if n_030 else 0.0,
         prec_006 if n_006 else 0.0]
npairs = [uni_n, n_030, n_006]
nconf = [uni_conf_n, conf_030, conf_006]
colors = [NAVY, ORANGE, GREEN]
bars = ax1.bar(range(3), [p if not np.isnan(p) else 0 for p in precs],
               color=colors, width=0.62, edgecolor="white")
for i, (p, nc, npr) in enumerate(zip(precs, nconf, npairs)):
    if npr == 0:
        ax1.text(i, 0.02, "no\npairs", ha="center", va="bottom",
                 fontsize=8, color="#888888")
        continue
    ax1.text(i, max(p, 0.0) + 0.03, f"{p*100:.0f}%", ha="center", va="bottom",
             fontsize=11, fontweight="bold", color=colors[i])
    lab_y = 0.04 if p > 0.10 else (p + 0.10)
    lab_col = "white" if p > 0.10 else colors[i]
    ax1.text(i, lab_y, f"{nc}/{npr}\npairs", ha="center", va="bottom",
             fontsize=8.5, color=lab_col)
ax1.set_xticks(range(3))
ax1.set_xticklabels(labels, fontsize=9)
ax1.set_ylabel("precision vs Ensembl Compara\n(within-group pairs confirmed paralog)")
ax1.set_ylim(0, 1.08)
ax1.set_title("Independent (gene-tree) confirmation of within-family pairs",
              fontsize=10.5)
ax1.axhline(1.0, color="#cccccc", lw=0.8, ls=":")

# Right: coverage funnel (text-light bars).
cov_labels = ["named\ngenes", "got\nENSG", ">=2-mapped\nfamilies", "evaluable\npairs"]
cov_vals = [n_named, n_named_mapped, n_evaluable_fams, uni_n]
ax2.bar(range(4), cov_vals, color=[NAVY, NAVY, ORANGE, GREEN], width=0.62,
        edgecolor="white")
for i, v in enumerate(cov_vals):
    ax2.text(i, v + max(cov_vals) * 0.015, str(v), ha="center", va="bottom",
             fontsize=10, fontweight="bold", color=NAVY)
ax2.set_xticks(range(4))
ax2.set_xticklabels(cov_labels, fontsize=9)
ax2.set_ylabel("count")
ax2.set_ylim(0, max(cov_vals) * 1.15)
ax2.set_title("Coverage (small, honest sample)", fontsize=10.5)

fig.suptitle("Non-circular validation: human Ensembl Compara vs our families",
             fontsize=12.5)
fig.text(0.5, 0.005,
         "Compara is gene-tree-based & COARSER (superfamily-level) than our "
         "copy-clusters; PRECISION (subset confirmation) is the metric, not "
         "recall. Human-as-proxy for gorilla; named genes only.",
         ha="center", fontsize=7.3, color="#555555")
fig.tight_layout(rect=[0, 0.045, 1, 0.95])
fig.savefig(OUT_PNG, dpi=130)
print(f"[fig] wrote {OUT_PNG}", file=sys.stderr)

# --------------------------------------------------------------------------
# REPORT (.md)
# --------------------------------------------------------------------------
def pct(x):
    return "n/a" if (x is None or (isinstance(x, float) and np.isnan(x))) else f"{x*100:.0f}%"

# Headline number is the universe-family precision (the de-circularizing one).
headline = (
    f"On the **{n_evaluable_fams}** universe families ({uni_n} within-family "
    f"named gene pairs) that map to Ensembl, human Compara independently confirms "
    f"**{pct(uni_precision)}** ({uni_conf_n}/{uni_n}) of within-family pairs as "
    f"paralogs."
)

md = []
md.append("# Non-circular validation: Ensembl Compara vs our paralog families\n")
md.append(f"_Deterministic (no RNG). Minimizer k={K}, w={W}, canonical "
          f"blake2b-64. Generated by `bench/compara_validation.py`. Compara "
          f"release: `{comp['meta'].get('compara')}`, "
          f"`{comp['meta'].get('rest_endpoint')}`._\n")

md.append("## Headline (the non-circular number)\n")
md.append(f"> {headline}\n")
md.append(
    "This **de-circularizes** the earlier `family_detection_validation.py`, "
    "which compared one sequence-similarity clustering (minimizer-Jaccard) "
    "against another (the minimap2-built `universe.tsv`) -- both reward shared "
    "subsequence, so that agreement was partly methodological. Here the truth "
    "set is Ensembl Compara paralogy: **gene-tree + species-tree "
    "reconciliation**, which never sees our minimap2 alignments. A Compara "
    "'yes' is therefore an independent biological confirmation.\n")

md.append("## KEY FRAMING: Compara is COARSER, so PRECISION is the metric\n")
md.append(
    "Ensembl Compara paralogy is gene-tree-based and far **coarser** than our "
    "families. Querying Compara for `RABL2A` returns the entire RAB / small-GTPase "
    f"**superfamily** ({len(para.get('RABL2A', []))} paralogues genome-wide in "
    "this dump), whereas our family is the 2-copy tandem `RABL2A`/`RABL2B`. So "
    "**our families should be SUBSETS of Compara groups**: we split ancient "
    "superfamilies into recent copy-clusters, which is correct granularity, not "
    "error. The right question is therefore PRECISION -- *of the gene pairs we "
    "put in the same family, what fraction does Compara independently also call "
    "paralogs?* -- not recall. Compara lumping ancient paralogs we correctly "
    "separate is expected granularity and is **intentionally not scored as a "
    "miss**.\n")

md.append("## (4) Coverage (stated up front -- small, honest sample)\n")
md.append(f"- Universe: **{len(g2f)} genes / {len(uni_fams)} families** "
          f"(`universe.tsv`). Of these, **{comp['meta'].get('loc_only_genes')} "
          f"genes are LOC-only** (no human symbol) and cannot be mapped to "
          f"Ensembl -- they are out of scope here.\n")
md.append(f"- Named genes attempted: **{n_named}**. Got an ENSG id: "
          f"**{n_named_mapped}/{n_named}** (1 lookup error: `GP1BB`). "
          f"Non-empty paralogue set returned: "
          f"**{comp['meta'].get('named_genes_with_nonempty_paralogues')}** "
          f"(4 paralogue-fetch errors: `CREB1`, `NCAPH2`, `RABL2B`, `USP18` "
          f"returned an error and so contribute NO confirmations even if real).\n")
md.append(f"- **Evaluable set: {n_evaluable_fams} universe families with >=2 "
          f"mapped named genes -> {uni_n} checkable within-family pairs.** "
          f"This is the entire basis of the headline; it is a SMALL sample and "
          f"is reported as exact counts, not extrapolated.\n")
if missing_seq:
    md.append(f"- (Note: {len(missing_seq)} mapped named genes lacked a "
              f"sequence: {missing_seq}.)\n")

md.append("## (1) PRECISION vs Compara -- UNIVERSE families (the headline)\n")
md.append(f"Of the **{uni_n}** within-universe-family named pairs (both genes "
          f"mapped to an ENSG), Compara independently confirms **{uni_conf_n}** "
          f"as paralogs => **precision = {pct(uni_precision)}** "
          f"({uni_conf_n}/{uni_n}).\n")
md.append("\n| universe family | pair | Compara confirms? |\n")
md.append("|---|---|---|\n")
for f, a, b, ok in uni_pairs_checkable:
    md.append(f"| `{f}` | {a} <-> {b} | {'**yes**' if ok else 'no'} |\n")
md.append(
    "\n**Reading the non-confirmations.** Most of the `no` rows are pairs where "
    "the universe's minimap2 step grouped two DIFFERENTLY-NAMED genes (e.g. "
    "`GCA<->KCNH7`, `CDPF1<->PPARA`, `CREB1<->METTL21A`) that Compara's gene "
    "tree does NOT place in the same paralog group. A non-confirmation can mean "
    "EITHER (a) a genuine universe over-grouping (sequence similarity from a "
    "shared domain / repeat, not orthologous-block paralogy), OR (b) a Compara "
    "gap / coarser-but-different granularity, OR (c) one partner's paralogue "
    "fetch errored (`CREB1`, `RABL2B` did -- so `CREB1<->METTL21A` and "
    "`RABL2A<->RABL2B` could not be confirmed even if true). The confirmed "
    "pairs are the unambiguous tandem/recent-duplicate families: "
    + ", ".join(f"`{a}<->{b}`" for _f, a, b, _ok in uni_confirmed) + ".\n")

md.append("## (2) PRECISION vs Compara -- RUSTLE's minimizer-Jaccard grouping\n")
md.append("Re-using the EXACT criterion from `family_detection_validation.py` "
          "(canonical k=15/w=10 minimizers, blake2b-64, union-find) over the "
          "same mapped named genes:\n")
md.append("\n| grouping | within-cluster pairs | Compara-confirmed | precision |\n")
md.append("|---|---|---|---|\n")
md.append(f"| universe families | {uni_n} | {uni_conf_n} | {pct(uni_precision)} |\n")
md.append(f"| rustle Jaccard @ {T_RUSTLE:g} | {n_030} | {conf_030} | "
          f"{pct(prec_030) if n_030 else 'no pairs'} |\n")
md.append(f"| rustle Jaccard @ {T_OPT:g} | {n_006} | {conf_006} | "
          f"{pct(prec_006) if n_006 else 'no pairs'} |\n")
md.append("\nrustle's non-trivial predicted clusters (mapped named genes only):\n")
md.append(f"- @ {T_RUSTLE:g}: " +
          ("; ".join("{" + ", ".join(c) + "}" for c in fmt_clusters(clusters_030))
           or "none (no named gene pair reaches J>=0.30)") + "\n")
md.append(f"- @ {T_OPT:g}: " +
          ("; ".join("{" + ", ".join(c) + "}" for c in fmt_clusters(clusters_006))
           or "none") + "\n")
md.append(
    "\nThe restriction to the 40 mapped NAMED genes means rustle's "
    "minimizer-Jaccard rarely puts two of them in one cluster at the strict "
    "0.30 bar -- the named genes that share a universe family are often "
    "different-symbol genes whose whole-gene Jaccard is below 0.30 (the same "
    "broad within-family Jaccard distribution documented in the prior report). "
    "Where rustle DOES cluster a named pair, that pair is the relevant "
    "precision test; counts are reported exactly above and are small.\n")

md.append("## (3) GRANULARITY (observation, not error)\n")
md.append(f"- Non-trivial Compara groups containing >=2 of our mapped named "
          f"genes: **{len(compara_groups_nt)}**.\n")
md.append(f"- Compara groups that MERGE >1 distinct universe family "
          f"(Compara coarser): **{n_compara_merging_multi}**. This is EXPECTED "
          f"granularity -- Compara lumps at the superfamily level what we split "
          f"into copy-clusters -- and is reported as an observation, NOT an "
          f"error.\n")
if compara_merges_universe:
    md.append("\n| Compara group (our genes) | distinct universe families merged |\n")
    md.append("|---|---|\n")
    for grp, uf in compara_merges_universe:
        md.append(f"| {', '.join(grp)} | {len(uf)}: "
                  f"{', '.join('`'+x+'`' for x in uf)} |\n")
md.append(f"\n- REVERSE check (real over-merge by us): universe families whose "
          f"mapped named genes span >1 Compara component: "
          f"**{len(universe_spanning)}**.\n")
if universe_spanning:
    for f, nm, ncomp in universe_spanning:
        md.append(f"  - `{f}` ({', '.join(nm)}) spans {ncomp} Compara "
                  f"components -- i.e. Compara does not link all of these into "
                  f"one paralog group. NOTE this is the SAME signal as the "
                  f"non-confirmed pairs above, re-expressed at the component "
                  f"level: a universe family that pairs genes Compara does not "
                  f"relate. Whether it is a genuine over-merge or a Compara/"
                  f"granularity gap is exactly the ambiguity named in caveat 2.\n")

md.append("## Cross-check: JSON's stated within-universe paralog relation\n")
md.append(f"The JSON records **{len(json_within)}** within-universe symmetric "
          f"Compara paralog pairs (both genes mapped + linked): "
          + ", ".join(f"`{a}<->{b}`" for a, b in sorted(json_within)) + ". "
          "Note `GGT1<->GGTLC2` is in this list but the two genes sit in "
          "SEPARATE universe families (`GGT1` and `GGTLC2`), so Compara "
          "confirms a paralogy that our minimap2 universe **split** -- a recall "
          "observation (granularity), not a precision failure, and consistent "
          "with the framing that we cluster more finely than Compara.\n")

md.append("## Honest caveats\n")
md.append(
    "- **Human-as-proxy for gorilla.** The sequences are gorilla genes; "
    "Compara paralogy is queried on the HUMAN ortholog's ENSG. Paralogy is "
    "deeply conserved across great apes (these duplications predate the "
    "~10-Mya gorilla/human split), so human paralogy is a strong proxy -- but "
    "it is a proxy, not the gorilla gene tree.\n")
md.append(
    "- **Compara is COARSER (superfamily-level).** A NON-confirmation can mean "
    "a genuine universe false-grouping OR a Compara gap / different "
    "granularity. Every non-confirmed pair is named above so the reader can "
    "judge; we do not silently treat them all as errors.\n")
md.append(
    "- **Only NAMED genes mapped.** "
    f"{comp['meta'].get('loc_only_genes')} of {len(g2f)} universe genes are "
    "LOC-only (no symbol -> no ENSG) and are out of scope. This validates the "
    "NAMED-gene backbone of the families, not all 62 families.\n")
md.append(
    "- **Symbol -> Ensembl mapping noise.** Genes were resolved by symbol via "
    f"the Ensembl REST lookup; {n_named_mapped}/{n_named} got an ENSG, with 1 "
    "lookup error and 4 paralogue-fetch errors. An errored fetch yields an "
    "EMPTY paralogue set, which can only DEPRESS the confirmation rate (it "
    "cannot create a false confirmation) -- so the reported precision is a "
    "conservative lower bound for the affected pairs.\n")
md.append(
    "- **Recall vs Compara is intentionally NOT the metric.** Compara lumping "
    "ancient paralogs we correctly split is expected granularity, not error "
    "(see the framing section). We report granularity as an observation only.\n")
md.append(
    "- **Small N.** Only "
    f"{n_evaluable_fams} families / {uni_n} pairs are checkable. The headline is "
    "an exact count on a small sample, not a generalized rate; treat it as a "
    "directional, independent sanity check rather than a population estimate.\n")
md.append(f"- **Determinism.** No RNG anywhere (the full mapped set is "
          f"evaluated; no sampling). Minimizer hash is blake2b. Output is "
          f"byte-stable.\n")

with open(OUT_MD, "w") as fh:
    fh.write("".join(md))
print(f"[md] wrote {OUT_MD}", file=sys.stderr)

# --------------------------------------------------------------------------
# Machine-readable summary
# --------------------------------------------------------------------------
print("=== SUMMARY ===")
print(f"named_genes={n_named} mapped_with_ensg={n_named_mapped}")
print(f"evaluable_families={n_evaluable_fams} evaluable_pairs={uni_n}")
print(f"UNIVERSE_PRECISION={uni_precision:.4f} ({uni_conf_n}/{uni_n})")
print(f"confirmed_pairs=" + ";".join(f"{a}<->{b}" for _f, a, b, _ok in uni_confirmed))
print(f"unconfirmed_pairs=" + ";".join(f"{a}<->{b}" for _f, a, b, _ok in uni_unconfirmed))
print(f"RUSTLE@0.30 pairs={n_030} confirmed={conf_030} "
      f"precision={prec_030 if n_030 else float('nan')}")
print(f"RUSTLE@0.06 pairs={n_006} confirmed={conf_006} "
      f"precision={prec_006 if n_006 else float('nan')}")
print(f"compara_groups_nt={len(compara_groups_nt)} "
      f"compara_merging_multi_universe={n_compara_merging_multi} "
      f"universe_spanning_multi_compara={len(universe_spanning)}")
print(f"json_within_universe_pairs={len(json_within)}")
print(f"HEADLINE: {headline}")
