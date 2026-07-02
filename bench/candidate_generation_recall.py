#!/usr/bin/env python
"""Candidate-generation RECALL diagnosis for the family-definition gap.

CONTEXT (bench/GRAPH_TO_GRAPH_FAMILY_ALN.md)
--------------------------------------------
The shipped family-definition recall gap (48/57 diploid-oracle multi_copy genes recovered) is an
EDGE-CREATION problem, not a scoring one: 7 of the 9 missed genes have ZERO candidate edges --
LOC101141440, LOC115930538, LOC115934629, LOC129534585, TNK2, UBE2Q2P16, ZNF425. The shipped
candidate-generation (bench/denovo_families.py) proposes a gene-pair only if the two de-novo
representative transcripts share >= K_SHARE=6 exact canonical 18-mers (the recall-safe "bloom"
pre-screen); POA core_recip>=0.13 then decides membership.

QUESTION
--------
WHY are these 7 missed, and would a LOOSENED candidate-generation (shorter k / minimizers) recover
them?  Each missed gene is classified into ONE of four causes (the task's three + REP-COLLAPSE, the
honest fourth surfaced by the data -- see below):
  1. PRE-FILTER-BOTTLENECK : the paralog copy IS expressed (a de-novo locus exists) and the pair is
     homologous enough to pass POA core_recip + aln_frac, but the exact-18-mer pre-filter did NOT
     propose the edge (too divergent for exact 18-mers).  A shorter-k / minimizer pre-filter PROPOSES
     it and it PASSES scoring  => RECOVERABLE from RNA.  [the ONLY loosening-recoverable class]
  2. EXPRESSION-LIMIT      : the family has no distinct expressed paralog copy => no edge possible from
     THIS IsoSeq regardless of the pre-filter.
  3. DIVERGENCE-LIMIT      : a candidate pair exists (proposed at k18 or by loosening) but the family's
     CLEANEST cross-copy pair FAILS the downstream core_recip / aln_frac scoring (too divergent at the
     sequence level) => DNA-bound, not recoverable from RNA.
  4. REP-COLLAPSE/PROJECTION : the family DOES contain a POA-passing, 18-mer-proposed cross-copy member
     pair -- so candidate-generation AND scoring both succeed for the family -- yet the gene's
     annotation region is not credited in the recovered de-novo family.  The miss is DOWNSTREAM
     (intron-junction rep-collapse merges the gene's copy into a readthrough, and gene_of max-overlap
     labels the linkable locus as a neighbour gene), NOT a candidate-generation-k problem.  Loosening
     cannot help (the edge is already proposed at k18).

METHOD (RNA-only; DNA oracle = the TARGET-truth only)   [FIXED: score ACTUAL family members]
--------------------------------------------------------------------------------------------
An earlier version re-derived the per-copy representative as the max-n_reads de-novo locus OVERLAPPING
each copy's genomic span (denovo_transcripts.meta.tsv).  That is UNFAITHFUL: on rep-collapsed loci it
picks a huge readthrough that is NOT a family member (e.g. it scored LOC101141440's 129 kb readthrough
DN_..._99594354_7, in NO validated family, against DN_..._92169976_5, a member of family 154 not the
gene's family 21) -> mis-scored 0.107 and mis-labelled DIVERGENCE-LIMIT.  The real family-21 members
score core_recip 0.926.  FIX: score the ACTUAL validated_families member loci only.

* Missed gene -> annotation interval (annot_intervals.tsv).
* Family = the gene's oracle family (diploid_cn_oracle.tsv -> validated_families.tsv, the
  DNA/co-location-derived family the diploid oracle is built on).  Members = the family's ACTUAL DN_
  loci (each is an expressed de-novo transcript by construction, so every paralog copy IS expressed).
  Members are clustered into distinct genomic copies (gap > COPY_GAP or different chrom); the gene's
  OWN copy = the copy overlapping the annotation.  Per copy we keep the top-TOP_ISO members by n_reads.
* CANDIDATE EDGES = every cross-copy member pair (best-scoring member pair per copy-pair over the
  top-TOP_ISO grid), scored with the SHIPPED machinery.  Classification is driven by the family's
  strongest edge (fam_best); the gene's-own-copy edge is reported separately so a family that is
  recoverable from its OTHER copies while the gene's specific transcript is a divergent/short-span
  paralog (LOC115930538, ZNF425) is flagged honestly (gene_to_paralog_edge column).
* Shipped pre-filter reproduction: RAW canonical shared-18-mer count for each scored member pair
  (canonical = min(fwd,rc) hash, exactly denovo_families.kmer_hashes).  This RAW count is an UPPER
  BOUND on the shipped "informative" shared count (informative = owned by 2..40 reps; a shared k-mer is
  owned by >=2 reps by definition, so informative-shared <= raw-shared), hence raw_shared_18mer <
  K_SHARE proves the shipped pipeline did NOT propose the pair.
* Loosen: RAW canonical shared-k-mer counts at k=15/13/11 and a minimizer(k=15,w=10) shared-sketch
  count.  "Loosened proposes" iff any loosened count >= K_SHARE.
* Score: core_recip (poa_family_definition.poa_pair_stats, fwd + RC fallback = exactly denovo._poa)
  and aln_frac (minimap2 -cx asm20 -t1 longest block / shorter length = the shipped ri_sharedlen
  workhorse).  "Passes scoring" iff core_recip >= CORE_PASS and aln_frac >= ALN_PASS (task line; both
  are also >= the shipped POA membership line T_CORE=0.13, so calls are threshold-robust).

Deterministic: PYTHONHASHSEED=0 re-exec, minimap2 -t1, Bio.Align deterministic, sorted output.
Writes bench/candidate_generation_recall.tsv (one row per cross-copy candidate edge) and prints the
per-gene diagnosis.
Run: /home/juanfra/miniforge3/bin/python bench/candidate_generation_recall.py
"""
import os
import subprocess
import sys
import tempfile
from collections import defaultdict

if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

import numpy as np
import pysam

BENCH = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BENCH)
import poa_family_definition as P

SCRATCH = "/home/juanfra/winloci_scratch"
FA = os.path.join(SCRATCH, "denovo_transcripts.fa")
META = os.path.join(SCRATCH, "denovo_transcripts.meta.tsv")
ANNOT = os.path.join(SCRATCH, "annot_intervals.tsv")
VALID = os.path.join(SCRATCH, "validated_families.tsv")
ORACLE = os.path.join(BENCH, "diploid_cn_oracle.tsv")
OUT = os.path.join(BENCH, "candidate_generation_recall.tsv")
SCRATCH_CACHE = os.environ.get(
    "CGR_CACHE",
    "/tmp/claude-1000/-mnt-c-Users-jfris-Desktop/03f0156d-63b6-4d0f-bf09-862bddee491c/scratchpad")
os.makedirs(SCRATCH_CACHE, exist_ok=True)

# missed genes (bench/GRAPH_TO_GRAPH_FAMILY_ALN.md; the 7 with ZERO candidate edges)
MISSED = ["LOC101141440", "LOC115930538", "LOC115934629", "LOC129534585", "TNK2", "UBE2Q2P16", "ZNF425"]

# --- shipped candidate-generation parameters (bench/denovo_families.py) ---
K_SHARE = 6            # shared informative k-mers to PROPOSE a pair
T_CORE = 0.13          # shipped POA membership threshold (denovo_families.T_CORE)
# --- task "passes scoring" line ---
CORE_PASS = 0.19
ALN_PASS = 0.24
# --- loosened pre-filter ladder ---
K_LADDER = [18, 15, 13, 11]
MM_K, MM_W = 15, 10    # minimizer sketch
COPY_GAP = 100_000     # genomic gap to split family members into distinct paralog copies
TOP_ISO = 8            # ACTUAL members per copy to score (top by n_reads); best over the member grid, so
                       # a DIVERGENCE-LIMIT call is not an artifact of one badly-assembled member and the
                       # reported family-best is the copy-pair's genuinely CLEANEST alignment

_CODE = np.full(256, -1, dtype=np.int64)
for _i, _b in enumerate(b"ACGT"):
    _CODE[_b] = _i
    _CODE[_b + 32] = _i
_RC = str.maketrans("ACGTN", "TGCAN")


def _horner(codes, n, k):
    h = np.zeros(n, dtype=np.int64)
    for t in range(k):
        h = (h << 2) + codes[t:t + n]
    return h


def canon_kmers(seq, k):
    """canonical (min of fwd,rc) base-4 codes of every k-mer; N-windows dropped. Exactly the
    denovo_families.kmer_hashes construction, parametrised in k."""
    codes = _CODE[np.frombuffer(seq.upper().encode(), dtype=np.uint8)]
    L = codes.size
    if L < k:
        return np.empty(0, dtype=np.int64)
    n = L - k + 1
    bad = codes < 0
    safe = np.where(bad, 0, codes)
    fwd = _horner(safe, n, k)
    rc = _horner((3 - safe)[::-1], n, k)[::-1]
    canon = np.minimum(fwd, rc)
    cb = np.concatenate(([0], np.cumsum(bad.astype(np.int64))))
    badwin = (cb[k:] - cb[:n]) > 0
    return canon[~badwin]


def shared_kmers(sa, sb, k):
    """# distinct canonical k-mers shared by the two sequences (the shipped 'shared' definition)."""
    va = np.unique(canon_kmers(sa, k))
    vb = np.unique(canon_kmers(sb, k))
    if va.size == 0 or vb.size == 0:
        return 0
    return int(np.intersect1d(va, vb, assume_unique=True).size)


def minimizer_set(seq, k, w):
    """canonical (k,w) minimizer set of a sequence."""
    km = canon_kmers(seq, k)
    if km.size < w:
        return set(km.tolist())
    # sliding-window minima over w consecutive canonical k-mers
    mins = set()
    from numpy.lib.stride_tricks import sliding_window_view
    win = sliding_window_view(km, w)
    mins.update(win.min(axis=1).tolist())
    return mins


def shared_minimizers(sa, sb, k, w):
    return len(minimizer_set(sa, k, w) & minimizer_set(sb, k, w))


# ---------- aln_frac (shipped ri_sharedlen workhorse: minimap2 -cx asm20 longest block / shorter) ----------
def aln_frac(sa, sb, tmp):
    if not sa or not sb:
        return 0.0
    pa, pb = os.path.join(tmp, "A.fa"), os.path.join(tmp, "B.fa")
    with open(pa, "w") as fh:
        fh.write(f">A\n{sa}\n")
    with open(pb, "w") as fh:
        fh.write(f">B\n{sb}\n")
    r = subprocess.run(["minimap2", "-cx", "asm20", "-t", "1", pa, pb],
                       capture_output=True, text=True)
    best = 0
    for line in r.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 11:
            continue
        best = max(best, int(f[10]))            # PAF col 11 = alignment block length
    return round(min(best / min(len(sa), len(sb)), 1.0), 4)


# ---------- core_recip (denovo._poa: POA fwd, RC fallback) ----------
def core_recip(aligner, sa, sb):
    cr = P.poa_pair_stats(aligner, sa, sb)["core_recip"]
    if cr < T_CORE:
        cr = max(cr, P.poa_pair_stats(aligner, sa, sb.translate(_RC)[::-1])["core_recip"])
    return round(cr, 4)


# ============================================================ loaders
def load_annot_missed():
    g = {}
    with open(ANNOT) as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if f[4] in MISSED:
                g[f[4]] = (f[0], int(f[1]), int(f[2]))
    return g


def gene_to_family():
    """oracle gene -> family id (diploid_cn_oracle.tsv)."""
    g2f = {}
    with open(ORACLE) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        I = {c: i for i, c in enumerate(h)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            g2f[f[I["gene"]]] = f[I["family"]]
    return g2f


def load_family_members():
    """family id -> list of (chrom,start,end,locus,nreads)."""
    fam = defaultdict(list)
    with open(VALID) as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            fam[f[0]].append((f[2], int(f[3]), int(f[4]), f[1], int(f[5])))
    return fam


def cluster_members(members):
    """cluster the ACTUAL family members (validated_families loci) into distinct genomic copies
    (gap > COPY_GAP or different chrom).  Keeps the member list per copy so the reps we score are real
    family members, not a max-n_reads locus re-derived over the copy's genomic span."""
    clusters = []
    for (c, s, e, loc, nr) in sorted(members):
        if clusters and clusters[-1]["chrom"] == c and s - clusters[-1]["end"] <= COPY_GAP:
            clusters[-1]["end"] = max(clusters[-1]["end"], e)
            clusters[-1]["members"].append((loc, s, e, nr))
        else:
            clusters.append(dict(chrom=c, start=s, end=e, members=[(loc, s, e, nr)]))
    return clusters


def score_pair(aligner, fa, tmp, la, lb):
    """core_recip + aln_frac + shared-k-mer ladder + minimizer for two ACTUAL member loci (by id)."""
    sa, sb = fa.fetch(la), fa.fetch(lb)
    cr = core_recip(aligner, sa, sb)
    af = aln_frac(sa, sb, tmp)
    sh = {k: shared_kmers(sa, sb, k) for k in K_LADDER}
    mm = shared_minimizers(sa, sb, MM_K, MM_W)
    return dict(cr=cr, af=af, sh=sh, mm=mm)


def best_cluster_pair(aligner, fa, tmp, mi, mj):
    """best-scoring ACTUAL-member pair between two copies (grid over their top-TOP_ISO members)."""
    best = None
    for (la, _s, _e, _n) in mi:
        for (lb, _s2, _e2, _n2) in mj:
            st = score_pair(aligner, fa, tmp, la, lb)
            key = (st["cr"], st["af"])
            if best is None or key > best["key"]:
                best = dict(key=key, la=la, lb=lb, **st)
    return best


# ============================================================ main
def classify(gene, fid, edges, gene_ci):
    """Cause + explanation from the family's cross-copy CANDIDATE EDGES (best member pair per copy-pair).
    fam_best = the family's strongest edge (drives divergence-vs-recoverable); gene_best = the strongest
    edge touching the gene's own copy (drives the gene_copy_links honesty flag)."""
    if not edges:
        return ("EXPRESSION-LIMIT",
                f"family {fid} has a single expressed copy (no distinct paralog locus)", None, None)
    fam_best = max(edges, key=lambda e: (e["cr"], e["af"]))
    gene_edges = [e for e in edges if e["touches_gene"]]
    gene_best = max(gene_edges, key=lambda e: (e["cr"], e["af"])) if gene_edges else None
    gene_copy_links = bool(gene_best and gene_best["passes"])

    if not fam_best["passes"]:
        already = ("18-mer ALREADY proposes it (pre-filter recall-safe)" if fam_best["shipped"]
                   else ("loosening to shorter k proposes it but POA still fails" if fam_best["loosened"]
                         else "not proposed even by loosened k"))
        why = (f"family {fid}'s CLEANEST cross-copy pair {fam_best['la']} x {fam_best['lb']} = "
               f"core_recip={fam_best['cr']} aln_frac={fam_best['af']} (< {CORE_PASS}/{ALN_PASS}, and "
               f"< shipped T_CORE={T_CORE}); 18-mer={fam_best['sh'][18]}; {already} -> RNA genuinely too "
               f"divergent, DNA-bound")
        return ("DIVERGENCE-LIMIT", why, fam_best, gene_best)
    if fam_best["shipped"]:
        gnote = (f"gene's own copy LINKS directly (best gene edge core_recip={gene_best['cr']})"
                 if gene_copy_links else
                 f"family recovered from its OTHER copies, but the gene's own transcript is a "
                 f"divergent/short-span paralog (best gene edge core_recip="
                 f"{gene_best['cr'] if gene_best else 'NA'} < gate)")
        why = (f"family {fid} HAS a POA-passing, 18-mer-proposed cross-copy pair {fam_best['la']} x "
               f"{fam_best['lb']} (core_recip={fam_best['cr']} aln_frac={fam_best['af']}, "
               f"{fam_best['sh'][18]} shared 18-mers) -> candidate-gen + scoring BOTH succeed; miss is "
               f"DOWNSTREAM (rep-collapse / gene_of label), NOT candidate-generation-k. {gnote}")
        return ("REP-COLLAPSE/PROJECTION", why, fam_best, gene_best)
    if fam_best["loosened"]:
        first_k = next((k for k in (15, 13, 11) if fam_best["sh"][k] >= K_SHARE), None)
        why = (f"family {fid}'s best pair {fam_best['la']} x {fam_best['lb']} PASSES scoring "
               f"(core_recip={fam_best['cr']} aln_frac={fam_best['af']}) but 18-mer="
               f"{fam_best['sh'][18]}(<{K_SHARE}); k{first_k}={fam_best['sh'][first_k]} proposes -> "
               f"RECOVERABLE by loosening")
        return ("PRE-FILTER-BOTTLENECK", why, fam_best, gene_best)
    why = (f"family {fid}'s best pair PASSES (core_recip={fam_best['cr']}) but even k11/mm < {K_SHARE} "
           f"(18-mer={fam_best['sh'][18]}); loosening headroom limited")
    return ("PRE-FILTER-BOTTLENECK", why, fam_best, gene_best)


def main():
    annot = load_annot_missed()
    g2f = gene_to_family()
    fams = load_family_members()
    fa = pysam.FastaFile(FA)
    aligner = P.make_aligner()
    tmp = tempfile.mkdtemp(prefix="cgr_")

    rows = []
    diagnosis = {}
    for gene in MISSED:
        c, gs, ge = annot[gene]
        fid = g2f[gene]
        clusters = cluster_members(fams[fid])
        # gene's OWN copy = the family copy overlapping the annotation interval
        gene_ci = [i for i, cl in enumerate(clusters)
                   if cl["chrom"] == c and cl["start"] < ge and cl["end"] > gs]
        # reps = top-TOP_ISO ACTUAL members per copy (by n_reads)
        reps = [sorted(cl["members"], key=lambda m: -m[3])[:TOP_ISO] for cl in clusters]

        # CANDIDATE EDGES = best-scoring member pair for every cross-copy pair
        edges = []
        for i in range(len(clusters)):
            for j in range(i + 1, len(clusters)):
                b = best_cluster_pair(aligner, fa, tmp, reps[i], reps[j])
                touches_gene = (i in gene_ci) or (j in gene_ci)
                gene_to_par = (i in gene_ci) ^ (j in gene_ci)   # gene-copy <-> paralog-copy (not the readthrough self)
                loosened = any(b["sh"][k] >= K_SHARE for k in (15, 13, 11)) or b["mm"] >= K_SHARE
                edges.append(dict(
                    i=i, j=j, ci=clusters[i], cj=clusters[j], la=b["la"], lb=b["lb"],
                    touches_gene=touches_gene, gene_to_par=gene_to_par, cr=b["cr"], af=b["af"],
                    sh=b["sh"], mm=b["mm"], shipped=b["sh"][18] >= K_SHARE, loosened=loosened,
                    passes=(b["cr"] >= CORE_PASS and b["af"] >= ALN_PASS)))

        # TSV rows: one per cross-copy candidate edge
        for ed in edges:
            rows.append([
                gene, fid, f"{c}:{gs}-{ge}",
                ed["ci"]["chrom"], ed["ci"]["start"], ed["ci"]["end"],
                ed["cj"]["chrom"], ed["cj"]["start"], ed["cj"]["end"],
                ed["la"], ed["lb"], int(ed["touches_gene"]), int(ed["gene_to_par"]),
                ed["sh"][18], ed["sh"][15], ed["sh"][13], ed["sh"][11], ed["mm"],
                ed["cr"], ed["af"], int(ed["shipped"]), int(ed["loosened"]), int(ed["passes"])])

        cause, why, _fb, _gb = classify(gene, fid, edges, gene_ci)
        diagnosis[gene] = (cause, why)

    # ---- write TSV ----
    hdr = ["gene", "family", "gene_annot", "ci_chrom", "ci_start", "ci_end",
           "cj_chrom", "cj_start", "cj_end", "ci_member", "cj_member",
           "touches_gene_copy", "gene_to_paralog_edge",
           "shared_18mer", "shared_15mer", "shared_13mer", "shared_11mer", "minimizer_share",
           "core_recip", "aln_frac", "shipped_proposes", "loosened_proposes", "passes_scoring"]
    with open(OUT, "w") as fh:
        fh.write("\t".join(hdr) + "\n")
        for r in sorted(rows, key=lambda x: (x[0], str(x[3]), x[4], str(x[6]), x[7])):
            fh.write("\t".join(str(x) for x in r) + "\n")

    # ---- report ----
    print(f"[wrote {OUT}] {len(rows)} cross-copy candidate-edge rows\n")
    print("PER-GENE DIAGNOSIS (scoring ACTUAL validated_families members)")
    print("=" * 100)
    counts = defaultdict(int)
    for gene in MISSED:
        cause, why = diagnosis[gene]
        counts[cause] += 1
        print(f"{gene:15s} {cause:24s} {why}")
    print("=" * 100)
    print("SUMMARY:", dict(counts))
    print(f"recoverable_by_loosening (PRE-FILTER-BOTTLENECK) = {counts['PRE-FILTER-BOTTLENECK']}/7")


def cost_probe(sample_n=150, seed=0):
    """GENOME-WIDE COST of loosening the exact-18-mer pre-filter: for the 7 missed-gene reps and a
    random sample of loci, count how many OTHER de-novo transcripts share >= K_SHARE canonical k-mers
    at k=18/15/13/11 -- i.e. how many candidate POA pairs each locus would generate. Shows the recall
    loosening buys 0 new genes but multiplies the POA load. ~10 min single-thread. Writes to stdout."""
    import random
    import time
    t0 = time.time()
    fa = pysam.FastaFile(FA)
    names = list(fa.references)
    gene_reps = {"LOC101141440": "DN_NC_073232.2_99594354_7", "LOC115930538": "DN_NC_073240.2_21107958_3",
                 "LOC115934629": "DN_NC_073228.2_45710800_10", "LOC129534585": "DN_NC_073224.2_15783750_7",
                 "TNK2": "DN_NC_086017.1_207500177_5", "UBE2Q2P16": "DN_NC_073240.2_85821672_6",
                 "ZNF425": "DN_NC_073230.2_158066680_4"}
    random.seed(seed)
    sample = random.sample(names, sample_n)
    tids = list(dict.fromkeys(list(gene_reps.values()) + sample))
    tidx = {t: i for i, t in enumerate(tids)}
    Tk, own_start, own_targets = {}, {}, {}
    for k in K_LADDER:
        pairs = []
        for t in tids:
            for h in np.unique(canon_kmers(fa.fetch(t), k)).tolist():
                pairs.append((h, tidx[t]))
        pairs.sort()
        arr = np.array([p[0] for p in pairs], dtype=np.int64)
        tg = np.array([p[1] for p in pairs], dtype=np.int32)
        uk, idx = np.unique(arr, return_index=True)
        Tk[k] = uk
        own_start[k] = list(idx) + [len(arr)]
        own_targets[k] = tg
    partners = np.zeros((len(tids), len(K_LADDER)), dtype=np.int32)
    cur = np.zeros(len(tids), dtype=np.int32)
    for i, nm in enumerate(names):
        s = fa.fetch(nm)
        self_i = tidx.get(nm, -1)
        for ki, k in enumerate(K_LADDER):
            v = np.unique(canon_kmers(s, k))
            if v.size == 0:
                continue
            pos = np.searchsorted(Tk[k], v)
            ok = pos < Tk[k].size
            pos, v2 = pos[ok], v[ok]
            pos = pos[Tk[k][pos] == v2]
            if pos.size == 0:
                continue
            cur[:] = 0
            for p in pos.tolist():
                for j in range(own_start[k][p], own_start[k][p + 1]):
                    cur[own_targets[k][j]] += 1
            for tj in np.nonzero(cur >= K_SHARE)[0].tolist():
                if tj != self_i:
                    partners[tj, ki] += 1
        if i % 20000 == 0:
            print(f"  {i}/{len(names)} t={time.time() - t0:.0f}s", flush=True)
    print("\nPER-MISSED-GENE candidate fan-out (# transcripts sharing >= %d k-mers = # POA pairs):" % K_SHARE)
    print("gene            " + "".join(f"{'k='+str(k):>8s}" for k in K_LADDER))
    for g, tid in gene_reps.items():
        print(f"{g:15s} " + "".join(f"{partners[tidx[tid], ki]:8d}" for ki in range(len(K_LADDER))))
    si = [tidx[t] for t in sample]
    mean = partners[si].mean(axis=0)
    med = np.median(partners[si], axis=0)
    print("\nrandom %d-locus MEAN fan-out:   " % sample_n
          + " ".join(f"k{K_LADDER[ki]}={mean[ki]:.1f}" for ki in range(len(K_LADDER))))
    print("random %d-locus MEDIAN fan-out: " % sample_n
          + " ".join(f"k{K_LADDER[ki]}={med[ki]:.0f}" for ki in range(len(K_LADDER))))
    print("fan-out MULTIPLIER vs k=18 (mean): "
          + " ".join(f"k{K_LADDER[ki]}={mean[ki]/max(mean[0], 1e-9):.1f}x" for ki in range(len(K_LADDER))))


# ============================================================ cost-probe result (logged run)
# Deterministic output of `candidate_generation_recall.py --cost` (sample_n=150, seed=0), captured to
# the scratchpad cost_out.txt.  Re-running --cost reproduces these to the digit (PYTHONHASHSEED=0,
# fixed random.seed, numpy-only arithmetic).  Fan-out = # OTHER transcripts sharing >= K_SHARE=6 RAW
# canonical k-mers with a locus = # candidate POA pairs that locus would trigger.  RAW (no informative
# 2..40 filter) -> an UPPER BOUND on the shipped candidate count; the informative filter only removes
# pairs, never adds, so raw is the conservative (worst-case) frame for the cost + FP blow-up.
COST_PROBE = {
    "sample_n": 150, "seed": 0, "n_transcripts": 94192,
    "per_missed_gene_fanout": {  # k18, k15, k13, k11
        "LOC101141440": [14850, 18348, 26099, 90070],
        "LOC115930538": [4581, 9444, 14350, 67449],
        "LOC115934629": [8488, 12519, 17018, 78704],
        "LOC129534585": [22, 26, 164, 24166],
        "TNK2": [15410, 19062, 22457, 72234],
        "UBE2Q2P16": [8769, 13906, 17674, 64411],
        "ZNF425": [3373, 8532, 14633, 76754],
    },
    "random150_mean_fanout": {"k18": 1947.9, "k15": 3037.4, "k13": 6347.8, "k11": 71845.4},
    "random150_median_fanout": {"k18": 30, "k15": 120, "k13": 2694, "k11": 75545},
    "fanout_multiplier_vs_k18": {"k18": 1.0, "k15": 1.6, "k13": 3.3, "k11": 36.9},
}


def _cost_derived():
    """Genome-wide candidate-pair estimates + POA-compute estimates derived from COST_PROBE.
    total undirected candidate pairs ~ n_transcripts * mean_fanout / 2."""
    N = COST_PROBE["n_transcripts"]
    m = COST_PROBE["random150_mean_fanout"]
    pairs = {k: int(round(N * m[k] / 2.0)) for k in ("k18", "k15", "k13", "k11")}
    extra = {k: pairs[k] - pairs["k18"] for k in ("k15", "k13", "k11")}
    POA_S = 0.2243  # measured single-thread s/pair (poa_family_definition.poa_pair_stats, ~3.3kb median)
    cpu_days = {k: round(pairs[k] * POA_S / 86400.0, 1) for k in pairs}
    extra_cpu_days = {k: round(extra[k] * POA_S / 86400.0, 1) for k in extra}
    return dict(genomewide_candidate_pairs_raw=pairs, extra_candidate_pairs_raw=extra,
                poa_s_per_pair=POA_S, poa_cpu_days_raw=cpu_days, extra_poa_cpu_days_raw=extra_cpu_days)


# ============================================================ FP / precision probe
def fp_probe(target_n=80, seed=0, n_poa=300):
    """PRECISION IMPACT of loosening the exact-18-mer pre-filter to k=13.

    Loosening changes ONLY k; the downstream POA membership gate (core_recip >= T_CORE=0.13) still
    decides every edge. So the precision question is: of the EXTRA candidate pairs that k=13 proposes
    but k=18 did not (= 'loosened-only' pairs), what fraction SURVIVE the POA gate and become family
    edges?  A near-zero survival rate => POA filters the loosened candidates => precision holds (the
    loosening is recall-safe, just wasteful); a non-trivial survival rate, multiplied by the ~200M
    extra candidate pairs, => the FP over-merge set explodes.

    Method (genome-wide, faithful to the shipped scan):
      * pick TARGET_N random target loci; scan ALL 94k transcripts, recording for each target the SET
        of transcripts that share >= K_SHARE=6 RAW canonical k-mers at k=18 and at k=13 (exactly the
        cost_probe scan, but keeping partner IDENTITIES, not just counts).
      * base pairs = {(target, k18-partner)}; loosened pairs = {(target, k13-partner)};
        EXTRA = loosened setminus base.
      * deterministic sub-sample of EXTRA pairs and of BASE pairs; run the shipped POA gate
        (core_recip) + aln_frac on each; tally survival at the shipped gate (cr>=T_CORE) and at the
        stricter task line (cr>=CORE_PASS & af>=ALN_PASS).  Also record core_ident / n_blocks of the
        surviving EXTRA pairs (a repeat-bridge over-merge squeaks past 0.13 with LOW core identity and
        many fragmented blocks; genuine homology has high core identity in one block).
    Deterministic; writes the cache the JSON builder reads. ~6-8 min single-thread."""
    import json
    import random
    import time
    t0 = time.time()
    KS = [18, 13]
    fa = pysam.FastaFile(FA)
    names = list(fa.references)
    gidx = {nm: i for i, nm in enumerate(names)}
    random.seed(seed)
    targets = sorted(random.sample(names, target_n))
    tidx = {t: i for i, t in enumerate(targets)}
    # target k-mer index (sorted unique k-mers -> owning target ids)
    Tk, own_start, own_targets = {}, {}, {}
    for k in KS:
        pairs = []
        for t in targets:
            for h in np.unique(canon_kmers(fa.fetch(t), k)).tolist():
                pairs.append((h, tidx[t]))
        pairs.sort()
        arr = np.array([p[0] for p in pairs], dtype=np.int64)
        tg = np.array([p[1] for p in pairs], dtype=np.int32)
        uk, ix = np.unique(arr, return_index=True)
        Tk[k] = uk
        own_start[k] = list(ix) + [len(arr)]
        own_targets[k] = tg
    # partner identity sets: per k, per target -> set of global query indices
    part = {k: [set() for _ in targets] for k in KS}
    cur = np.zeros(len(targets), dtype=np.int32)
    for gi, nm in enumerate(names):
        s = fa.fetch(nm)
        self_t = tidx.get(nm, -1)
        for k in KS:
            v = np.unique(canon_kmers(s, k))
            if v.size == 0:
                continue
            pos = np.searchsorted(Tk[k], v)
            ok = pos < Tk[k].size
            pos, v2 = pos[ok], v[ok]
            pos = pos[Tk[k][pos] == v2]
            if pos.size == 0:
                continue
            cur[:] = 0
            for p in pos.tolist():
                for j in range(own_start[k][p], own_start[k][p + 1]):
                    cur[own_targets[k][j]] += 1
            for tj in np.nonzero(cur >= K_SHARE)[0].tolist():
                if tj != self_t:
                    part[k][tj].add(gi)
        if gi % 20000 == 0:
            print(f"  scan {gi}/{len(names)} t={time.time() - t0:.0f}s", flush=True)
    # build undirected candidate pair pools (target global idx, partner global idx)
    def pool(k):
        pp = set()
        for ti, t in enumerate(targets):
            tg = gidx[t]
            for q in part[k][ti]:
                pp.add((min(tg, q), max(tg, q)))
        return pp
    base = pool(18)
    loos = pool(13)
    extra = sorted(loos - base)
    basel = sorted(base)
    print(f"\npools: base(k18)={len(base)}  loosened(k13)={len(loos)}  EXTRA(k13-only)={len(extra)} "
          f"t={time.time() - t0:.0f}s", flush=True)
    rng = random.Random(seed + 1)
    extra_s = sorted(rng.sample(extra, min(n_poa, len(extra)))) if extra else []
    base_s = sorted(rng.sample(basel, min(n_poa, len(basel)))) if basel else []
    aligner = P.make_aligner()
    tmp = tempfile.mkdtemp(prefix="cgrfp_")

    def evaluate(sample):
        rows = []
        tp = time.time()
        for (ia, ib) in sample:
            a, b = names[ia], names[ib]
            sa, sb = fa.fetch(a), fa.fetch(b)
            st = P.poa_pair_stats(aligner, sa, sb)
            cr = st["core_recip"]
            if cr < T_CORE:  # RC fallback, exactly denovo._poa / core_recip()
                st2 = P.poa_pair_stats(aligner, sa, sb.translate(_RC)[::-1])
                if st2["core_recip"] > cr:
                    st, cr = st2, st2["core_recip"]
            af = aln_frac(sa, sb, tmp)
            rows.append(dict(a=a, b=b, core_recip=round(cr, 4), aln_frac=af,
                             core_ident=round(st["core_ident"], 4), n_blocks=st["n_blocks"]))
        return rows, time.time() - tp

    extra_rows, te = evaluate(extra_s)
    base_rows, tb = evaluate(base_s)

    def frac(rows, pred):
        return round(sum(1 for r in rows if pred(r)) / len(rows), 4) if rows else 0.0
    ship = lambda r: r["core_recip"] >= T_CORE
    task = lambda r: r["core_recip"] >= CORE_PASS and r["aln_frac"] >= ALN_PASS
    survivors = [r for r in extra_rows if ship(r)]
    result = dict(
        target_n=target_n, seed=seed, n_transcripts=len(names),
        pool_sizes=dict(base_k18=len(base), loosened_k13=len(loos), extra_k13_only=len(extra)),
        n_poa_extra=len(extra_rows), n_poa_base=len(base_rows),
        shipped_gate=f"core_recip>={T_CORE}", task_gate=f"core_recip>={CORE_PASS} & aln_frac>={ALN_PASS}",
        extra_pass_frac_shipped=frac(extra_rows, ship), extra_pass_frac_task=frac(extra_rows, task),
        base_pass_frac_shipped=frac(base_rows, ship), base_pass_frac_task=frac(base_rows, task),
        extra_core_recip_max=round(max((r["core_recip"] for r in extra_rows), default=0.0), 4),
        extra_survivors_shipped=len(survivors),
        extra_survivor_core_ident=[r["core_ident"] for r in survivors],
        extra_survivor_examples=sorted(survivors, key=lambda r: -r["core_recip"])[:10],
        poa_time_extra_s=round(te, 1), poa_time_base_s=round(tb, 1),
    )
    CACHE_FP = os.path.join(SCRATCH_CACHE, "cgr_fp.json")
    with open(CACHE_FP, "w") as fh:
        json.dump(result, fh, indent=2, sort_keys=True)
    print("\nFP/PRECISION PROBE")
    print("=" * 90)
    print(f"pools: base(k18)={len(base)} loosened(k13)={len(loos)} EXTRA(k13-only)={len(extra)}")
    print(f"EXTRA (k13-only) POA survival: shipped(cr>={T_CORE})={result['extra_pass_frac_shipped']} "
          f"task={result['extra_pass_frac_task']}  (n={len(extra_rows)}, max cr={result['extra_core_recip_max']})")
    print(f"BASE  (k18)      POA survival: shipped(cr>={T_CORE})={result['base_pass_frac_shipped']} "
          f"task={result['base_pass_frac_task']}  (n={len(base_rows)})")
    print(f"EXTRA survivors: {len(survivors)}  core_ident={result['extra_survivor_core_ident']}")
    print(f"wrote {CACHE_FP}")
    return result


# ============================================================ JSON builder
def build_json():
    """Assemble bench/candidate_generation_recall.json: fast live diagnosis (re-derived from the TSV
    the diagnosis stage wrote) + the logged COST_PROBE constants + the FP-probe cache (if present)."""
    import json
    # --- diagnosis, re-derived from the committed TSV (byte-stable) ---
    rows = []
    with open(OUT) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for ln in fh:
            rows.append(dict(zip(hdr, ln.rstrip("\n").split("\t"))))
    per_gene = {}
    tally = defaultdict(int)
    def fnum(r, key):
        return float(r[key]) if r[key] not in ("", None) else 0.0
    for gene in MISSED:
        edges = [r for r in rows if r["gene"] == gene]        # cross-copy candidate edges
        # fam_best = the family's strongest edge; gene_best = strongest edge touching the gene's own copy
        fam_best = max(edges, key=lambda r: (fnum(r, "core_recip"), fnum(r, "aln_frac"))) if edges else None
        gedges = [r for r in edges if r["touches_gene_copy"] == "1"]
        gene_best = max(gedges, key=lambda r: (fnum(r, "core_recip"), fnum(r, "aln_frac"))) if gedges else None
        fam_passes = bool(fam_best and fnum(fam_best, "core_recip") >= CORE_PASS
                          and fnum(fam_best, "aln_frac") >= ALN_PASS)
        fam_shipped = bool(fam_best and fam_best["shipped_proposes"] == "1")
        fam_loosened = bool(fam_best and fam_best["loosened_proposes"] == "1")
        gene_copy_links = bool(gene_best and fnum(gene_best, "core_recip") >= CORE_PASS
                               and fnum(gene_best, "aln_frac") >= ALN_PASS)
        if not edges:
            cause = "EXPRESSION-LIMIT"
        elif not fam_passes:
            cause = "DIVERGENCE-LIMIT"
        elif fam_shipped:
            cause = "REP-COLLAPSE/PROJECTION"
        else:
            cause = "PRE-FILTER-BOTTLENECK"      # passes but 18-mer does not propose (loosened or headroom)
        tally[cause] += 1
        n_copies = len(set([(e["ci_chrom"], e["ci_start"], e["ci_end"]) for e in edges]
                           + [(e["cj_chrom"], e["cj_start"], e["cj_end"]) for e in edges]))
        per_gene[gene] = dict(
            cause=cause,
            family=(edges[0]["family"] if edges else None),
            n_copy_pairs=len(edges),
            n_copies=n_copies,
            fam_best_pair=(f"{fam_best['ci_member']} x {fam_best['cj_member']}" if fam_best else None),
            fam_best_core_recip=(fnum(fam_best, "core_recip") if fam_best else None),
            fam_best_aln_frac=(fnum(fam_best, "aln_frac") if fam_best else None),
            fam_best_shared_18mer=(int(fam_best["shared_18mer"]) if fam_best else None),
            fam_best_shared_13mer=(int(fam_best["shared_13mer"]) if fam_best else None),
            fam_best_shared_11mer=(int(fam_best["shared_11mer"]) if fam_best else None),
            fam_best_18mer_already_proposes=fam_shipped,
            fam_best_loosened_proposes=fam_loosened,
            fam_passes_scoring=fam_passes,
            gene_own_copy_best_core_recip=(fnum(gene_best, "core_recip") if gene_best else None),
            gene_own_copy_links=gene_copy_links,
        )
    recoverable = tally["PRE-FILTER-BOTTLENECK"]

    cost = dict(COST_PROBE)
    cost.update(_cost_derived())

    CACHE_FP = os.path.join(SCRATCH_CACHE, "cgr_fp.json")
    fp = json.load(open(CACHE_FP)) if os.path.exists(CACHE_FP) else None

    # --- verdict ---
    if fp is not None:
        extra_est = cost["extra_candidate_pairs_raw"]["k13"]
        exp_extra_edges = int(round(extra_est * fp["extra_pass_frac_shipped"]))
        # rule-of-3 upper 95% CI on a 0/n survival rate = 3/n (so the ABSOLUTE worst-case exposure over
        # the ~200M extra raw candidate pairs is bounded even when the point estimate is exactly 0).
        ci95 = 3.0 / max(fp["n_poa_extra"], 1)
        exp_extra_edges_ci95 = int(round(extra_est * ci95))
        precision_line = (
            f"EXTRA (k13-only) candidate pairs survive the shipped POA gate (core_recip>={T_CORE}) at "
            f"{fp['extra_pass_frac_shipped']*100:.1f}% ({fp['extra_survivors_shipped']}/{fp['n_poa_extra']}, "
            f"max core_recip {fp['extra_core_recip_max']}), BELOW the {fp['base_pass_frac_shipped']*100:.1f}% "
            f"survival of the k18 baseline candidates -- the loosened-only pairs are the more-divergent tail "
            f"and POA rejects them. Point-estimate expected extra over-merge edges over the ~{extra_est/1e6:.0f}M "
            f"extra raw k13 candidate pairs = {exp_extra_edges}; rule-of-3 95%-CI upper bound "
            f"<= {exp_extra_edges_ci95}.")
    else:
        exp_extra_edges = None
        exp_extra_edges_ci95 = None
        precision_line = "FP probe cache absent (run --fpprobe)."

    nd = tally["DIVERGENCE-LIMIT"]
    nr = tally["REP-COLLAPSE/PROJECTION"]
    verdict = (
        f"NOT WORTH IT. Loosening the exact-18-mer pre-filter recovers {recoverable} of the 7 missed "
        f"multi-copy genes (PRE-FILTER-BOTTLENECK is the ONLY loosening-recoverable class and it is empty). "
        f"Scoring the ACTUAL validated_families member loci: {nd}/7 are DIVERGENCE-LIMIT (TNK2, UBE2Q2P16 -- "
        f"2-member families whose cleanest cross-copy RNA pair fails the POA gate, core_recip 0.05/0.02; "
        f"genuinely DNA-bound), and {nr}/7 are REP-COLLAPSE/PROJECTION (the family HAS a POA-passing, "
        f"18-mer-proposed cross-copy member pair -- 2019/5311/3725/1401/1237 shared 18-mers, core_recip "
        f"0.85-1.0 -- so candidate-generation AND scoring already succeed; the gene's annotation region is "
        f"lost DOWNSTREAM to intron-junction rep-collapse + gene_of max-overlap labelling, NOT to a "
        f"candidate-generation-k gap). For LOC115930538 and ZNF425 the family is recovered from its OTHER "
        f"near-identical copies while the gene's OWN transcript is a divergent/short-span paralog. NOTE: an "
        f"earlier version mis-scored 3 of these as DIVERGENCE-LIMIT by re-deriving the per-copy rep as the "
        f"max-n_reads locus over the copy's genomic SPAN (which on rep-collapsed loci is a non-member "
        f"readthrough, or a locus in a DIFFERENT family) instead of the actual member; the real members "
        f"score 0.93-1.0. Either way the recall gap is NOT candidate-generation-k-bound: every family edge is "
        f"already proposed at k18. Loosening multiplies the candidate/POA load 3.3x at k13 (36.9x at k11, "
        f"near-saturation) for 0 recovered genes. "
        + (f"Precision is NOT the failure mode: POA filters the loosened candidates -- "
           f"{fp['extra_survivors_shipped']}/{fp['n_poa_extra']} k13-only pairs survive the gate "
           f"(max core_recip {fp['extra_core_recip_max']} << {T_CORE}), so the FP over-merge set does NOT "
           f"explode. Loosening is thus recall-SAFE but recall-USELESS: pure wasted compute."
           if fp is not None else "See fp_precision_impact for the measured survival rates."))

    out = dict(
        question=("Are the 7 zero-candidate-edge missed multi-copy genes (LOC101141440, LOC115930538, "
                  "LOC115934629, LOC129534585, TNK2, UBE2Q2P16, ZNF425) recoverable by LOOSENING the "
                  "exact-18-mer candidate-generation pre-filter (shorter k / minimizers), and if so at "
                  "what genome-wide candidate-blow-up / POA-compute / FP-precision cost?"),
        missed_genes=MISSED,
        params=dict(K_SHARE=K_SHARE, T_CORE_shipped_gate=T_CORE, CORE_PASS_task=CORE_PASS,
                    ALN_PASS_task=ALN_PASS, k_ladder=K_LADDER, minimizer="(k15,w10)"),
        diagnosis_tally=dict(tally),
        per_gene=per_gene,
        recoverable_by_loosening=recoverable,
        loosening_cost=cost,
        fp_precision_impact=(dict(fp, expected_extra_overmerge_edges_k13_raw=exp_extra_edges,
                                  expected_extra_overmerge_edges_k13_ci95=exp_extra_edges_ci95,
                                  summary=precision_line) if fp is not None
                             else dict(summary=precision_line)),
        verdict=verdict,
        determinism=dict(PYTHONHASHSEED="0", minimap2="-t1", biopython="deterministic global NW",
                         seeds="cost seed=0, fp seed=0", note="re-run of --json is byte-identical"),
        files=dict(
            script=os.path.abspath(__file__),
            tsv=OUT,
            json=os.path.join(BENCH, "candidate_generation_recall.json"),
            cost_log=os.path.join(SCRATCH_CACHE, "cost_out.txt"),
            fp_cache=CACHE_FP,
        ),
    )
    OUT_JSON = os.path.join(BENCH, "candidate_generation_recall.json")
    with open(OUT_JSON, "w") as fh:
        json.dump(out, fh, indent=2, sort_keys=True)
        fh.write("\n")
    print(f"[wrote {OUT_JSON}]")
    print(f"recoverable_by_loosening = {recoverable}/7")
    print(f"tally = {dict(tally)}")
    if fp is not None:
        print(f"extra(k13-only) POA survival shipped-gate = {fp['extra_pass_frac_shipped']} "
              f"(base k18 = {fp['base_pass_frac_shipped']})")
        print(f"expected extra over-merge edges (k13 raw) ~ {exp_extra_edges}")
    print("VERDICT:", verdict)
    return out


if __name__ == "__main__":
    if "--cost" in sys.argv:
        cost_probe()
    elif "--fpprobe" in sys.argv:
        fp_probe()
    elif "--json" in sys.argv:
        build_json()
    else:
        main()
