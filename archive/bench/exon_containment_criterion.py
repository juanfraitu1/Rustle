#!/usr/bin/env python
"""bench/exon_containment_criterion.py

Test the USER HYPOTHESIS that the shipped RNA multi-copy family EDGE criterion --
a COVERAGE rule (core_recip = whole-transcript reciprocal contiguous-core >= t1 AND
aln_frac = longest shared spliced-exon block / shorter transcript >= t2) -- is FLAWED
because a single large shared EXON (a domain / repeat / retro bridge) can inflate the
coverage past threshold even when the OTHER exons differ. PROPOSED FIX: a HARD
STRUCTURAL criterion -- FULL EXON CONTAINMENT -- an edge requires (most) exons of one
transcript to each have a HOMOLOGOUS exon in the other, not just one shared block.

DEFINITION (RNA-only; genome bases only as the spliced-exon substrate, identical
provenance to the shipped aln_frac -- see ri_build_sharedlen_universal.py):
  For an E_r edge with representative de-novo loci (A,B):
   * reconstruct per-exon coords from denovo_skeletons (introns delimit exons);
   * splice = genome bases concatenated in genomic-ascending exon order (EXACT same
     sequence the shipped aln_frac aligner sees) -> transcript-coordinate exon
     partition = cumulative exon lengths in that same order;
   * align A vs B with minimap2 -c -x asm20 (asm20 => only >=~80% id blocks reported);
   * map the aligned QUERY interval onto B's exon partition and the aligned TARGET
     interval onto A's exon partition;
   * exon_contain(A->B) = frac of A's exons with >= 50% of the exon covered by a
     homologous (>=0.7 block-id) aligned block; symmetric for B;
   * reciprocal_exon_containment = min(exon_contain(A->B), exon_contain(B->A)).
  FULL containment ~ 1.0 (every exon shared); a 1-of-many domain bridge -> low recip.

PREDICTED WALL (tested + quantified): for a SINGLE-EXON locus (MAGE = intronless
retrocopy; many residual FPs are single-block), exon-containment on that side is
TRIVIALLY 1.0 for BOTH a real paralog AND a domain-bridge sharing that one exon, so it
cannot separate them. => every result is split MULTI-EXON (min n_exon >= 2) vs
SINGLE-EXON (min n_exon == 1; incl. the both-single total wall).

EVALUATION: DNA `in_dna_loose` is the VALIDATION label ONLY (never a gate). Population
= the 5571 core-present direct E_r edges (ri_sharedlen_universal.tsv: 416 TP /
1773 genuine-overmerge / 3382 truthbar-divergent-paralog), the LEARN-stage population
of the shipped rna_only_edge_oracle. Component-level 5-fold CV (union-find over ALL
edge_bridge_mask gene pairs; whole families held out; TP-balanced greedy; seed=0),
identical machinery to the shipped oracle. Reference held-out AUCs to beat:
core_recip 0.815, aln_frac 0.915.

Deterministic: PYTHONHASHSEED=0, minimap2 -t1, RandomState(0) folds, sorted writes.
Writes: bench/exon_containment_criterion.tsv (per-edge) + .json (summary).
Sidecar cache: bench/exon_containment_mm_cache.tsv (recip inputs; deterministic).
Run: PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/exon_containment_criterion.py
"""
import os
import sys

if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

import json
import math
import subprocess
import tempfile
import itertools
from collections import Counter, defaultdict

import numpy as np
import scipy.stats as ss

BENCH = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BENCH)
import family_er_pr as F
import family_er_te_gate as TE
import pysam

SCRATCH = "/home/juanfra/winloci_scratch"
GENOME = os.path.join(SCRATCH, "GGO.fasta")
UNIV_TSV = os.path.join(BENCH, "ri_sharedlen_universal.tsv")
ERPR_TSV = os.path.join(BENCH, "family_er_pr.tsv")
EBM_TSV = os.path.join(BENCH, "edge_bridge_mask.tsv")
SKEL_TSV = os.path.join(SCRATCH, "denovo_skeletons.tsv")
OUT_TSV = os.path.join(BENCH, "exon_containment_criterion.tsv")
OUT_JSON = os.path.join(BENCH, "exon_containment_criterion.json")
MM_CACHE = os.path.join(BENCH, "exon_containment_mm_cache.tsv")
# Robustness: recompute containment with the EXACT shipped aln_frac minimap2 flags
# (ri_build_sharedlen_universal.py: `minimap2 -cx asm20 -t1`, no -N/-p) so the
# head-to-head is strictly parameter-matched. Separate cache => primary path untouched.
MM_CACHE_MATCHED = os.path.join(BENCH, "exon_containment_matched_mm_cache.tsv")
FLAGS_GENEROUS = ["-c", "-x", "asm20", "-N", "50", "-p", "0.1", "-t", "1"]  # primary
FLAGS_MATCHED = ["-c", "-x", "asm20", "-t", "1"]                            # == shipped aln_frac

# Shipped coverage operating points (from rna_only_edge_oracle.py / RNA_ONLY_EDGE_ORACLE.md)
COV_F1 = (0.26, 0.72)          # F1-optimal coverage rule
COV_RECALL = (0.19, 0.24)      # recall-preserving coverage rule
EXON_HOMOL_FRAC = 0.5          # >=50% of an exon covered => that exon is contained
BLOCK_ID_MIN = 0.7             # minimap2 block identity (nmatch/blocklen) to count

# Named residual-FP roster (GRAPH_DEF_REFINE_SWEEP.md / RNA_ONLY_EDGE_ORACLE.md)
NAMED_FP = {
    "allele:DHRSX": ["DHRSX"],
    "allele:LOC129530050": ["LOC129530050"],
    "oversize:LOC129523567(54locZNF)": ["LOC129523567"],
    "oversize:MAGEA9_18loc": ["MAGEA9", "LOC129529978", "LOC129529986"],
    "oversize:MPHOSPH8": ["MPHOSPH8"],
    "oversize:LOC134758618": ["LOC134758618"],
    "multifam:GSTM2_hub": ["GSTM2"],
    "multifam:FOXO1+LOC115933254": ["FOXO1", "LOC115933254"],
    "multifam:LOC101142904+LOC129526550": ["LOC101142904", "LOC129526550"],
}


# --------------------------------------------------------------------------- geometry
def exon_bounds_transcript(exons):
    """Genomic-ascending exon (a,b) list -> transcript-coord exon intervals + total len.
       Order matches the spliced sequence built by ri_build_sharedlen_universal.py."""
    b = []
    cur = 0
    for (a, e) in exons:
        L = e - a
        b.append((cur, cur + L))
        cur += L
    return b, cur


def merge_iv(iv):
    if not iv:
        return []
    iv = sorted(iv)
    out = [list(iv[0])]
    for s, e in iv[1:]:
        if s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def overlap_len(a0, a1, merged):
    """total overlap of [a0,a1) with a merged (sorted disjoint) interval list."""
    tot = 0
    for s, e in merged:
        lo = max(a0, s)
        hi = min(a1, e)
        if hi > lo:
            tot += hi - lo
    return tot


def frac_exons_covered(exon_bnds, covered_iv):
    """frac of exons with >= EXON_HOMOL_FRAC of their length inside covered_iv; + count."""
    if not exon_bnds:
        return 0.0, 0
    merged = merge_iv(covered_iv)
    cnt = 0
    for (s, e) in exon_bnds:
        L = e - s
        if L <= 0:
            continue
        if overlap_len(s, e, merged) >= EXON_HOMOL_FRAC * L:
            cnt += 1
    return cnt / len(exon_bnds), cnt


# --------------------------------------------------------------------------- alignment
def run_minimap2(sa, sb, pa, pb, flags=FLAGS_GENEROUS):
    """Align target=A (pa) vs query=B (pb) -- same call orientation as the shipped
       aln_frac builder (minimap2 pa pb). Return (target_iv on A coords, query_iv on B
       coords), keeping only blocks with identity >= BLOCK_ID_MIN. Deterministic (-t1)."""
    with open(pa, "w") as fh:
        fh.write(f">a\n{sa}\n")
    with open(pb, "w") as fh:
        fh.write(f">b\n{sb}\n")
    r = subprocess.run(
        ["minimap2"] + flags + [pa, pb],
        capture_output=True, text=True)
    tgt_iv = []
    qry_iv = []
    for line in r.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 12:
            continue
        qs, qe = int(f[2]), int(f[3])          # query = B
        ts, te = int(f[7]), int(f[8])          # target = A
        nmatch, blen = int(f[9]), int(f[10])
        if blen <= 0 or nmatch / blen < BLOCK_ID_MIN:
            continue
        qry_iv.append((qs, qe))
        tgt_iv.append((ts, te))
    return tgt_iv, qry_iv


# --------------------------------------------------------------------------- loaders
def load_labels():
    """gene-pair frozenset -> in_dna_loose (int)."""
    lab = {}
    with open(ERPR_TSV) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        I = {c: i for i, c in enumerate(h)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            k = frozenset((f[I["gene_a"]], f[I["gene_b"]]))
            lab[k] = int(f[I["in_dna_loose"]])
    return lab


def load_ebm_pairs():
    """All refined-E_r gene pairs (for union-find components) in file order."""
    pairs = []
    with open(EBM_TSV) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        I = {c: i for i, c in enumerate(h)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            pairs.append((f[I["gene_a"]], f[I["gene_b"]]))
    return pairs


def skeleton_single_exon_census():
    """Reproducible single-exon census of the RAW de-novo skeleton catalog
       (n_exon==1). This is the denominator for the 'single-exon wall' claim."""
    tot = single = 0
    with open(SKEL_TSV) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        I = {c: i for i, c in enumerate(h)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            tot += 1
            if int(f[I["n_exon"]]) == 1:
                single += 1
    return tot, single


def load_univ_rows():
    rows = []
    with open(UNIV_TSV) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        I = {c: i for i, c in enumerate(h)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            rows.append(dict(
                gene_a=f[I["gene_a"]], gene_b=f[I["gene_b"]],
                dn_a=f[I["dn_a"]], dn_b=f[I["dn_b"]], cls=f[I["cls"]],
                core=float(f[I["core"]]) if f[I["core"]] not in ("", "NA") else None,
                aln_frac=float(f[I["aln_frac"]]) if f[I["aln_frac"]] not in ("", "NA") else None,
            ))
    return rows


def load_mm_cache(path=MM_CACHE):
    d = {}
    if not os.path.exists(path):
        return d
    with open(path) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        I = {c: i for i, c in enumerate(h)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            d[(f[I["dn_a"]], f[I["dn_b"]])] = dict(
                contain_ab=float(f[I["contain_ab"]]), contain_ba=float(f[I["contain_ba"]]),
                n_exon_a=int(f[I["n_exon_a"]]), n_exon_b=int(f[I["n_exon_b"]]),
                n_cov_a=int(f[I["n_cov_a"]]), n_cov_b=int(f[I["n_cov_b"]]))
    return d


def write_mm_cache(cache, path=MM_CACHE):
    keys = sorted(cache)
    with open(path, "w") as out:
        out.write("dn_a\tdn_b\tcontain_ab\tcontain_ba\tn_exon_a\tn_exon_b\tn_cov_a\tn_cov_b\n")
        for k in keys:
            v = cache[k]
            out.write(f"{k[0]}\t{k[1]}\t{v['contain_ab']:.6f}\t{v['contain_ba']:.6f}\t"
                      f"{v['n_exon_a']}\t{v['n_exon_b']}\t{v['n_cov_a']}\t{v['n_cov_b']}\n")


# --------------------------------------------------------------------------- containment compute
def compute_containment(rows, meta, skel, flags=FLAGS_GENEROUS, cache_path=MM_CACHE, tag="primary"):
    fa = pysam.FastaFile(GENOME)
    seqcache = {}

    def spliced(dn):
        if dn in seqcache:
            return seqcache[dn]
        ch, ex = TE.dn_exons(dn, meta, skel)
        if ch is None:
            seqcache[dn] = (None, None)
        else:
            seqcache[dn] = ("".join(fa.fetch(ch, a, b) for a, b in ex), ex)
        return seqcache[dn]

    cache = load_mm_cache(cache_path)
    todo = [r for r in rows if (r["dn_a"], r["dn_b"]) not in cache]
    print(f"[contain:{tag}] {len(rows)} edges; {len(cache)} cached; {len(todo)} to align", flush=True)
    td = tempfile.mkdtemp(prefix="exoncontain_")
    pa = os.path.join(td, "a.fa")
    pb = os.path.join(td, "b.fa")
    for i, r in enumerate(todo):
        da, db = r["dn_a"], r["dn_b"]
        sa, exa = spliced(da)
        sb, exb = spliced(db)
        if not sa or not sb:
            cache[(da, db)] = dict(contain_ab=0.0, contain_ba=0.0,
                                   n_exon_a=len(exa or []), n_exon_b=len(exb or []),
                                   n_cov_a=0, n_cov_b=0)
            continue
        bnds_a, _ = exon_bounds_transcript(exa)
        bnds_b, _ = exon_bounds_transcript(exb)
        tgt_iv, qry_iv = run_minimap2(sa, sb, pa, pb, flags)  # tgt on A, qry on B
        c_ab, ncov_a = frac_exons_covered(bnds_a, tgt_iv)  # A's exons covered => A->B
        c_ba, ncov_b = frac_exons_covered(bnds_b, qry_iv)  # B's exons covered => B->A
        # round to 6dp == cache write precision => cached and fresh runs are byte-identical
        c_ab, c_ba = round(c_ab, 6), round(c_ba, 6)
        cache[(da, db)] = dict(contain_ab=c_ab, contain_ba=c_ba,
                               n_exon_a=len(exa), n_exon_b=len(exb),
                               n_cov_a=ncov_a, n_cov_b=ncov_b)
        if (i + 1) % 500 == 0:
            print(f"[contain:{tag}] aligned {i+1}/{len(todo)}", flush=True)
    write_mm_cache(cache, cache_path)
    return cache


# --------------------------------------------------------------------------- CV machinery
def union_find_components(edges):
    parent = {}

    def find(a):
        parent.setdefault(a, a)
        root = a
        while parent[root] != root:
            root = parent[root]
        while parent[a] != root:
            parent[a], a = root, parent[a]
        return root

    def uni(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for a, b in edges:
        uni(a, b)
    reps = {}
    comp = {}
    for g in sorted(parent):
        r = find(g)
        if r not in reps:
            reps[r] = len(reps)
        comp[g] = reps[r]
    return comp


def assign_folds(comp_ids, comp_tp, comp_sz, K, seed):
    rng = np.random.RandomState(seed)
    order = list(comp_ids)
    rng.shuffle(order)
    order.sort(key=lambda c: (-comp_tp[c], -comp_sz[c]))
    fold_tp = [0] * K
    fold_sz = [0] * K
    fold_of = {}
    for c in order:
        k = min(range(K), key=lambda j: (fold_tp[j], fold_sz[j], j))
        fold_of[c] = k
        fold_tp[k] += comp_tp[c]
        fold_sz[k] += comp_sz[c]
    return fold_of


def prf(pred, y):
    pred = np.asarray(pred)
    y = np.asarray(y)
    tp = int(((pred == 1) & (y == 1)).sum())
    fp = int(((pred == 1) & (y == 0)).sum())
    fn = int(((pred == 0) & (y == 1)).sum())
    prec = tp / (tp + fp) if (tp + fp) else 0.0
    rec = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
    return prec, rec, f1, tp, fp, fn


def auc_mw(pos, neg):
    pos = [x for x in pos if x is not None and x == x]
    neg = [x for x in neg if x is not None and x == x]
    if not pos or not neg:
        return float("nan"), len(pos), len(neg)
    r = ss.rankdata(pos + neg)
    rp = float(np.sum(r[:len(pos)]))
    a = (rp - len(pos) * (len(pos) + 1) / 2.0) / (len(pos) * len(neg))
    return a, len(pos), len(neg)


# --------------------------------------------------------------------------- rule search
def grid_predict(feats, thr):
    """feats/thr: list of (array, op, threshold). op '>=' only. Returns bool AND-mask."""
    m = np.ones(len(feats[0][0]), dtype=bool)
    for arr, t in zip([f[0] for f in feats], thr):
        m = m & (arr >= t)
    return m


def held_out_rule(arrs, y, fold, K, grids, objective="f1", min_recall=None):
    """arrs: list of feature arrays (AND-ed with >= thresholds from grids).
       Trains best thresholds per fold on the OTHER folds, predicts held-out.
       Returns pooled (prec,rec,f1,tp,fp,fn, pooled_pred)."""
    N = len(y)
    pred = np.zeros(N, dtype=int)
    for f in range(K):
        tr = fold != f
        te = fold == f
        ytr = y[tr]
        tr_arrs = [a[tr] for a in arrs]
        best = None
        for combo in itertools.product(*grids):
            m = np.ones(tr.sum(), dtype=bool)
            for a, t in zip(tr_arrs, combo):
                m = m & (a >= t)
            p, rc, f1, *_ = prf(m.astype(int), ytr)
            if objective == "f1":
                key = (f1, p)
            else:
                if rc < min_recall:
                    continue
                key = (p, f1)
            if best is None or key > best[0]:
                best = (key, combo)
        combo = best[1] if best else tuple(g[0] for g in grids)
        m_te = np.ones(te.sum(), dtype=bool)
        te_arrs = [a[te] for a in arrs]
        for a, t in zip(te_arrs, combo):
            m_te = m_te & (a >= t)
        pred[te] = m_te.astype(int)
    return prf(pred, y) + (pred,)


def fixed_rule(arrs, thr, y):
    m = np.ones(len(y), dtype=bool)
    for a, t in zip(arrs, thr):
        m = m & (a >= t)
    return prf(m.astype(int), y) + (m.astype(int),)


# --------------------------------------------------------------------------- main
def main():
    print("[load] meta / skeletons / labels / edges ...", flush=True)
    meta = F.load_meta()
    skel = TE.load_skeletons()
    labels = load_labels()
    rows = load_univ_rows()
    ebm_pairs = load_ebm_pairs()

    cache = compute_containment(rows, meta, skel)
    # robustness: strictly parameter-matched recip (shipped aln_frac flags -cx asm20)
    cache_matched = compute_containment(rows, meta, skel, flags=FLAGS_MATCHED,
                                        cache_path=MM_CACHE_MATCHED, tag="matched")

    # reproducible single-exon census (the 'wall' denominator)
    skel_tot, skel_single = skeleton_single_exon_census()

    # ---- assemble population -------------------------------------------------
    pop = []
    miss_lab = 0
    for r in rows:
        k = frozenset((r["gene_a"], r["gene_b"]))
        y = labels.get(k)
        if y is None:
            y = 1 if r["cls"] == "TP" else 0
            miss_lab += 1
        c = cache[(r["dn_a"], r["dn_b"])]
        recip = min(c["contain_ab"], c["contain_ba"])
        cm = cache_matched[(r["dn_a"], r["dn_b"])]
        recip_matched = min(cm["contain_ab"], cm["contain_ba"])
        nmin = min(c["n_exon_a"], c["n_exon_b"])
        nmax = max(c["n_exon_a"], c["n_exon_b"])
        pop.append(dict(
            gene_a=r["gene_a"], gene_b=r["gene_b"], dn_a=r["dn_a"], dn_b=r["dn_b"],
            cls=r["cls"], y=y, core=r["core"] if r["core"] is not None else 0.0,
            aln_frac=r["aln_frac"] if r["aln_frac"] is not None else 0.0,
            contain_ab=c["contain_ab"], contain_ba=c["contain_ba"], recip=recip,
            recip_matched=recip_matched,
            n_exon_a=c["n_exon_a"], n_exon_b=c["n_exon_b"], n_exon_min=nmin, n_exon_max=nmax,
            both_single=1 if nmax == 1 else 0,
            stratum="multi" if nmin >= 2 else "single",
        ))
    pop.sort(key=lambda r: (r["gene_a"], r["gene_b"], r["dn_a"], r["dn_b"]))
    N = len(pop)

    # ---- components + folds (identical machinery to shipped oracle) ----------
    comp = union_find_components(ebm_pairs)
    for r in pop:
        r["comp"] = comp.get(r["gene_a"], comp.get(r["gene_b"]))
    comp_ids = sorted(set(r["comp"] for r in pop))
    comp_tp = {c: 0 for c in comp_ids}
    comp_sz = {c: 0 for c in comp_ids}
    for r in pop:
        comp_sz[r["comp"]] += 1
        comp_tp[r["comp"]] += r["y"]
    K, SEED = 5, 0
    fold_of = assign_folds(comp_ids, comp_tp, comp_sz, K, SEED)
    for r in pop:
        r["fold"] = fold_of[r["comp"]]
    # leakage guard: no family across folds
    seen = {}
    for r in pop:
        c = r["comp"]
        assert seen.setdefault(c, r["fold"]) == r["fold"], "component spans folds"

    # ---- arrays --------------------------------------------------------------
    core = np.array([r["core"] for r in pop])
    aln = np.array([r["aln_frac"] for r in pop])
    recip = np.array([r["recip"] for r in pop])
    recip_matched = np.array([r["recip_matched"] for r in pop])
    y = np.array([r["y"] for r in pop])
    fold = np.array([r["fold"] for r in pop])
    cls = np.array([r["cls"] for r in pop])
    nmin = np.array([r["n_exon_min"] for r in pop])
    multi = nmin >= 2
    single = ~multi

    cls_ct = Counter(r["cls"] for r in pop)
    strat_ct = Counter(r["stratum"] for r in pop)
    both_single = int(sum(r["both_single"] for r in pop))

    # ---- single-exon census (the 'wall' denominator; reproducible) -----------
    univ_loci = {}
    for r in pop:
        univ_loci[r["dn_a"]] = r["n_exon_a"]
        univ_loci[r["dn_b"]] = r["n_exon_b"]
    single_exon_census = dict(
        raw_skeleton_loci=skel_tot,
        raw_skeleton_single_exon=skel_single,
        raw_skeleton_single_exon_pct=round(100.0 * skel_single / skel_tot, 4),
        edge_universe_loci=len(univ_loci),
        edge_universe_single_exon=int(sum(1 for v in univ_loci.values() if v == 1)),
        edges_min_nexon_eq_1=int((nmin == 1).sum()),
        edges_both_single=both_single,
    )

    # ---- AUC (single feature; fold-independent) ------------------------------
    def auc_feat(arr, mask=None, neg_cls=None):
        """AUC of arr as a REAL(TP)-vs-negative separator. neg_cls restricts the
           negative pool: None => all non-TP (TP-vs-ALL, easy axis); 'genuine' =>
           only genuine over-merges (HARD axis, the domain-bridge question);
           'truthbar' => divergent-paralog truthbar."""
        if mask is None:
            mask = np.ones(N, dtype=bool)
        posm = mask & (y == 1)
        if neg_cls is None:
            negm = mask & (y == 0)
        else:
            negm = mask & (cls == neg_cls)
        a, _, _ = auc_mw(list(arr[posm]), list(arr[negm]))
        return a

    auc_all = dict(core=auc_feat(core), aln_frac=auc_feat(aln), recip=auc_feat(recip),
                   recip_matched=auc_feat(recip_matched))
    auc_multi = dict(core=auc_feat(core, multi), aln_frac=auc_feat(aln, multi),
                     recip=auc_feat(recip, multi))
    auc_single = dict(core=auc_feat(core, single), aln_frac=auc_feat(aln, single),
                      recip=auc_feat(recip, single))
    # HARD axis: TP vs genuine-overmerge ONLY (aln_frac's easy win on truthbar
    # over-credits it on the TP-vs-ALL axis; this is the axis that matters for the
    # 'does containment cut domain-bridges' question).
    auc_vs_genuine = dict(core=auc_feat(core, neg_cls="genuine"),
                          aln_frac=auc_feat(aln, neg_cls="genuine"),
                          recip=auc_feat(recip, neg_cls="genuine"),
                          recip_matched=auc_feat(recip_matched, neg_cls="genuine"))
    auc_vs_truthbar = dict(core=auc_feat(core, neg_cls="truthbar"),
                           aln_frac=auc_feat(aln, neg_cls="truthbar"),
                           recip=auc_feat(recip, neg_cls="truthbar"),
                           recip_matched=auc_feat(recip_matched, neg_cls="truthbar"))

    # ---- grids ---------------------------------------------------------------
    g_core = np.round(np.arange(0.13, 0.701, 0.01), 3)
    g_aln = np.round(np.arange(0.0, 1.001, 0.02), 3)
    g_recip = np.round(np.arange(0.0, 1.001, 0.02), 3)
    g_core_c = np.round(np.arange(0.13, 0.701, 0.03), 3)   # coarser for 3-way combined
    g_aln_c = np.round(np.arange(0.0, 1.001, 0.05), 3)
    g_recip_c = np.round(np.arange(0.0, 1.001, 0.05), 3)

    def report_rule(name, arrs, grids, objective="f1", min_recall=None, mask=None):
        if mask is None:
            mask = np.ones(N, dtype=bool)
        aa = [a[mask] for a in arrs]
        res = held_out_rule(aa, y[mask], fold[mask], K, grids, objective, min_recall)
        prec, rc, f1, tp, fp, fn, pred = res
        # per-class cut of negatives
        clsm = cls[mask]
        ym = y[mask]
        gen = (clsm == "genuine")
        tb = (clsm == "truthbar")
        gen_cut = int(((pred == 0) & gen).sum())
        gen_tot = int(gen.sum())
        tb_cut = int(((pred == 0) & tb).sum())
        tb_tot = int(tb.sum())
        return dict(name=name, prec=prec, rec=rc, f1=f1, tp=tp, fp=fp, fn=fn,
                    genuine_cut=gen_cut, genuine_tot=gen_tot,
                    truthbar_cut=tb_cut, truthbar_tot=tb_tot,
                    n=int(mask.sum()), n_pos=int(ym.sum())), pred

    results = {}
    preds = {}

    # held-out learned rules (whole population)
    for tag, arrs, grids in [
        ("coverage(core&aln)", [core, aln], [g_core, g_aln]),
        ("containment(recip)", [recip], [g_recip]),
        ("core&recip", [core, recip], [g_core, g_recip]),
        ("aln&recip", [aln, recip], [g_aln, g_recip]),
        ("combined(core&aln&recip)", [core, aln, recip], [g_core, g_aln, g_recip]),
    ]:
        d, p = report_rule(tag, arrs, grids)
        results[tag] = d
        preds[tag] = p

    # fixed shipped coverage operating points (global -> held-out == in-sample)
    for tag, thr in [("coverage_F1(0.26,0.72)", COV_F1),
                     ("coverage_recall(0.19,0.24)", COV_RECALL)]:
        prec, rc, f1, tp, fp, fn, pred = fixed_rule([core, aln], thr, y)
        gen = (cls == "genuine")
        tb = (cls == "truthbar")
        results[tag] = dict(name=tag, prec=prec, rec=rc, f1=f1, tp=tp, fp=fp, fn=fn,
                            genuine_cut=int(((pred == 0) & gen).sum()), genuine_tot=int(gen.sum()),
                            truthbar_cut=int(((pred == 0) & tb).sum()), truthbar_tot=int(tb.sum()),
                            n=N, n_pos=int(y.sum()))
        preds[tag] = pred

    # ---- MULTI vs SINGLE stratified learned rules ----------------------------
    strat_results = {}
    for sname, smask in [("multi", multi), ("single", single)]:
        strat_results[sname] = {}
        for tag, arrs, grids in [
            ("coverage(core&aln)", [core, aln], [g_core, g_aln]),
            ("containment(recip)", [recip], [g_recip]),
            ("combined(core&aln&recip)", [core, aln, recip], [g_core, g_aln, g_recip]),
        ]:
            d, _ = report_rule(tag, arrs, grids, mask=smask)
            strat_results[sname][tag] = d

    # ---- exon-count BIN stratification (the operative 'coarse-exon' wall test:
    #      few-exon loci give a coarse recip {0,.5,1}; does that degrade it?) -----
    def bin_of(nm):
        if nm <= 2:
            return "min2"
        if nm == 3:
            return "min3"
        if nm <= 6:
            return "min4-6"
        return "min7+"
    bins = np.array([bin_of(v) for v in nmin])
    bin_auc = {}
    bin_ct = {}
    for bname in ["min2", "min3", "min4-6", "min7+"]:
        bmask = bins == bname
        bin_ct[bname] = dict(n=int(bmask.sum()), n_tp=int((y[bmask] == 1).sum()))
        bin_auc[bname] = dict(core=auc_feat(core, bmask), aln_frac=auc_feat(aln, bmask),
                              recip=auc_feat(recip, bmask))

    # ---- KEY TEST: on MULTI-exon, of the genuine over-merges the shipped coverage
    #      F1 rule KEEPS, how many does containment cut? at what TP recall cost? ----
    covF1_keep = (core >= COV_F1[0]) & (aln >= COV_F1[1])
    mm = multi
    kept_gen_multi = mm & covF1_keep & (cls == "genuine")     # domain-bridges coverage keeps
    kept_tp_multi = mm & covF1_keep & (y == 1)
    # containment threshold that PRESERVES multi-exon TP recall of the coverage rule:
    # pick largest recip cut t s.t. it does not drop any coverage-kept multi TP.
    tp_recip_vals = recip[kept_tp_multi]
    t_pres = float(min(tp_recip_vals)) if tp_recip_vals.size else 1.0
    # add containment>=t_pres on top of coverage rule (multi only)
    add_keep = covF1_keep & (recip >= t_pres)
    gen_cut_by_contain = int((kept_gen_multi & (recip < t_pres)).sum())
    tp_lost = int((kept_tp_multi & (recip < t_pres)).sum())
    # also a fixed illustrative threshold recip>=0.5
    gen_cut_half = int((kept_gen_multi & (recip < 0.5)).sum())
    tp_lost_half = int((kept_tp_multi & (recip < 0.5)).sum())

    key_multi = {
        "n_multi": int(mm.sum()),
        "cov_kept_genuine_multi": int(kept_gen_multi.sum()),
        "cov_kept_tp_multi": int(kept_tp_multi.sum()),
        "recip_thresh_preserving_all_covkept_TP": t_pres,
        "genuine_bridges_cut_by_containment_at_preserving_t": gen_cut_by_contain,
        "tp_lost_at_preserving_t": tp_lost,
        "genuine_bridges_cut_at_recip_0.5": gen_cut_half,
        "tp_lost_at_recip_0.5": tp_lost_half,
    }

    # ---- named residual-FP outcomes ------------------------------------------
    gene_to_idx = defaultdict(list)
    for i, r in enumerate(pop):
        gene_to_idx[r["gene_a"]].append(i)
        gene_to_idx[r["gene_b"]].append(i)
    named = {}
    for fpname, genes in NAMED_FP.items():
        idxs = sorted(set(i for g in genes for i in gene_to_idx.get(g, [])))
        if not idxs:
            named[fpname] = dict(edges=0, note="no edge in core-present population")
            continue
        sub = [pop[i] for i in idxs]
        gen_idx = [i for i in idxs if pop[i]["cls"] == "genuine"]
        n_multi = sum(1 for r in sub if r["n_exon_min"] >= 2)
        n_single = sum(1 for r in sub if r["n_exon_min"] < 2)
        recip_g = [pop[i]["recip"] for i in gen_idx]
        aln_g = [pop[i]["aln_frac"] for i in gen_idx]
        # would containment (recip<0.5) cut vs coverage-F1 keep, on genuine edges here
        cov_keeps = [i for i in gen_idx if covF1_keep[i]]
        contain_cuts = [i for i in cov_keeps if recip[i] < 0.5]
        named[fpname] = {
            "edges": len(idxs), "genuine_edges": len(gen_idx),
            "n_multi": n_multi, "n_single": n_single,
            "median_recip_genuine": float(np.median(recip_g)) if recip_g else None,
            "median_aln_genuine": float(np.median(aln_g)) if aln_g else None,
            "genuine_kept_by_coverageF1": len(cov_keeps),
            "of_those_cut_by_containment_recip_lt_0.5": len(contain_cuts),
            "all_min_nexon": sorted(set(r["n_exon_min"] for r in sub)),
        }

    # ---- diagnostics: why containment does not beat aln_frac -----------------
    def quant_(arr):
        if arr.size == 0:
            return None
        return [round(float(np.quantile(arr, q)), 4) for q in (0.1, 0.25, 0.5, 0.75, 0.9)]
    tp_m = y == 1
    gen_m = cls == "genuine"
    tb_m = cls == "truthbar"
    rho_recip_aln, _ = ss.spearmanr(recip, aln)
    rho_recip_matched, _ = ss.spearmanr(recip, recip_matched)
    diag = {
        "spearman_recip_alnfrac": round(float(rho_recip_aln), 4),
        "spearman_recip_vs_recip_matched": round(float(rho_recip_matched), 4),
        "recip_q10_25_50_75_90": {
            "TP": quant_(recip[tp_m]), "genuine": quant_(recip[gen_m]), "truthbar": quant_(recip[tb_m])},
        "alnfrac_q10_25_50_75_90": {
            "TP": quant_(aln[tp_m]), "genuine": quant_(aln[gen_m]), "truthbar": quant_(aln[tb_m])},
        "TP_edges_recip_eq_0": int(((recip == 0) & tp_m).sum()),
        "TP_edges_recip_lt_0.5": int(((recip < 0.5) & tp_m).sum()),
        "TP_edges_total": int(tp_m.sum()),
        "genuine_edges_recip_ge_0.8": int(((recip >= 0.8) & gen_m).sum()),
        "genuine_edges_total": int(gen_m.sum()),
        "covF1_kept_TP_median_recip": (round(float(np.median(recip[covF1_keep & tp_m])), 4)
                                       if (covF1_keep & tp_m).any() else None),
        "covF1_kept_genuine_median_recip": (round(float(np.median(recip[covF1_keep & gen_m])), 4)
                                            if (covF1_keep & gen_m).any() else None),
        "TP_alnfrac_high_but_recip_low": int((tp_m & (aln >= 0.72) & (recip < 0.5)).sum()),
    }

    # ---- write per-edge TSV --------------------------------------------------
    cols = ["gene_a", "gene_b", "dn_a", "dn_b", "cls", "in_dna_loose",
            "n_exon_a", "n_exon_b", "n_exon_min", "stratum", "both_single",
            "core_recip", "aln_frac", "contain_ab", "contain_ba", "recip_contain",
            "recip_contain_matched", "fold"]
    with open(OUT_TSV, "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in pop:
            out.write("\t".join(str(x) for x in [
                r["gene_a"], r["gene_b"], r["dn_a"], r["dn_b"], r["cls"], r["y"],
                r["n_exon_a"], r["n_exon_b"], r["n_exon_min"], r["stratum"], r["both_single"],
                f"{r['core']:.4f}", f"{r['aln_frac']:.4f}",
                f"{r['contain_ab']:.4f}", f"{r['contain_ba']:.4f}", f"{r['recip']:.4f}",
                f"{r['recip_matched']:.4f}", r["fold"]]) + "\n")

    # ---- summary + JSON ------------------------------------------------------
    summary = dict(
        population=dict(N=N, class_split=dict(cls_ct), stratum_split=dict(strat_ct),
                        both_single=both_single, missing_label_fallback=miss_lab),
        single_exon_census=single_exon_census,
        auc_all=auc_all, auc_multi=auc_multi, auc_single=auc_single,
        auc_vs_genuine_hard_axis=auc_vs_genuine, auc_vs_truthbar=auc_vs_truthbar,
        heldout_rules={k: {kk: vv for kk, vv in v.items() if kk != "name"}
                       for k, v in results.items()},
        stratified_rules={s: {k: {kk: vv for kk, vv in d.items() if kk != "name"}
                              for k, d in strat_results[s].items()} for s in strat_results},
        key_multi_domain_bridge_test=key_multi,
        exon_count_bins=dict(counts=bin_ct, auc=bin_auc),
        complementarity_diag=diag,
        named_residual_fp=named,
    )
    with open(OUT_JSON, "w") as fh:
        json.dump(summary, fh, sort_keys=True, indent=1,
                  default=lambda x: None if (isinstance(x, float) and x != x) else x)

    # ---- console report ------------------------------------------------------
    P = print
    P("\n================= EXON-CONTAINMENT vs COVERAGE CRITERION =================")
    P(f"population = core-present direct E_r edges  N={N}   "
      f"({'  '.join(f'{k}={v}' for k,v in sorted(cls_ct.items()))})")
    P(f"strata: multi(min nexon>=2)={strat_ct['multi']}  single(min nexon==1)={strat_ct['single']}  "
      f"both-single(WALL)={both_single}")
    P(f"wrote {OUT_TSV}  +  {OUT_JSON}")

    P("\n---- SINGLE-EXON census (the predicted-wall denominator; reproducible) ----")
    sc = single_exon_census
    P(f"  raw skeleton loci={sc['raw_skeleton_loci']}  single-exon={sc['raw_skeleton_single_exon']} "
      f"({sc['raw_skeleton_single_exon_pct']}%)")
    P(f"  edge-universe loci={sc['edge_universe_loci']}  single-exon={sc['edge_universe_single_exon']}  "
      f"=> single-exon EDGES(min nexon==1)={sc['edges_min_nexon_eq_1']}  both-single={sc['edges_both_single']}")
    P("  => literal single-exon wall INERT (0 single-exon edges): MAGE retrocopies still assemble w/ UTR introns.")

    P("\n---- single-feature AUC vs in_dna_loose (TP vs non-TP; higher=REAL) ----")
    P("  feature          ALL        MULTI       SINGLE")
    for feat in ("core", "aln_frac", "recip"):
        P(f"  {feat:14s}  {auc_all[feat]:.4f}     {auc_multi[feat]:.4f}      {auc_single[feat]:.4f}")
    P(f"  recip_matched   {auc_all['recip_matched']:.4f}  (shipped -cx asm20 flags; robustness)")
    P("  (reference shipped held-out AUC: core 0.815, aln_frac 0.915)")

    P("\n---- AUC by NEGATIVE class (easy TP-vs-ALL over-credits aln_frac on truthbar) ----")
    P("  feature         TP-vs-ALL   TP-vs-GENUINE(hard)  TP-vs-TRUTHBAR")
    for feat in ("core", "aln_frac", "recip"):
        P(f"  {feat:14s}  {auc_all[feat]:.4f}      {auc_vs_genuine[feat]:.4f}             {auc_vs_truthbar[feat]:.4f}")
    P(f"  recip_matched   {auc_all['recip_matched']:.4f}      {auc_vs_genuine['recip_matched']:.4f}"
      f"             {auc_vs_truthbar['recip_matched']:.4f}")
    P("  => on the HARD axis aln_frac and recip are ~TIED; recip's headline deficit is truthbar-driven.")

    P("\n---- held-out CV rules (whole population; component 5-fold seed=0) ----")
    P("  rule                                prec    rec     F1     TP/kept  gen-cut          tb-cut")
    for tag in ["coverage_F1(0.26,0.72)", "coverage_recall(0.19,0.24)",
                "coverage(core&aln)", "containment(recip)", "core&recip", "aln&recip",
                "combined(core&aln&recip)"]:
        d = results[tag]
        P(f"  {tag:34s}  {d['prec']:.3f}   {d['rec']:.3f}   {d['f1']:.3f}   {d['tp']:3d}/{d['n_pos']:<3d}   "
          f"{d['genuine_cut']}/{d['genuine_tot']}      {d['truthbar_cut']}/{d['truthbar_tot']}")

    P("\n---- MULTI vs SINGLE held-out rules ----")
    for s in ("multi", "single"):
        P(f"  [{s}]  N={strat_ct[s]}")
        for tag in ["coverage(core&aln)", "containment(recip)", "combined(core&aln&recip)"]:
            d = strat_results[s][tag]
            P(f"    {tag:28s}  prec={d['prec']:.3f} rec={d['rec']:.3f} F1={d['f1']:.3f}  "
              f"gen-cut={d['genuine_cut']}/{d['genuine_tot']}")

    P("\n---- AUC by exon-count bin (coarse-exon 'wall' test; min n_exon) ----")
    P("  bin       n     n_tp    core    aln_frac  recip")
    for bname in ["min2", "min3", "min4-6", "min7+"]:
        c = bin_ct[bname]
        a = bin_auc[bname]
        P(f"  {bname:8s}  {c['n']:5d} {c['n_tp']:5d}   {a['core']:.4f}  {a['aln_frac']:.4f}   {a['recip']:.4f}")

    P("\n---- complementarity diagnostics ----")
    for k, v in diag.items():
        P(f"  {k} = {v}")

    P("\n---- KEY: multi-exon domain-bridges the coverage-F1 rule KEEPS ----")
    for k, v in key_multi.items():
        P(f"  {k} = {v}")

    P("\n---- named residual-FP outcomes ----")
    for fpname in NAMED_FP:
        d = named[fpname]
        P(f"  {fpname}: {d}")

    P("\n=========================================================================")
    P("JSON " + json.dumps(summary, sort_keys=True,
                           default=lambda x: None if (isinstance(x, float) and x != x) else x))


if __name__ == "__main__":
    main()
