#!/usr/bin/env python3
"""
vg_repeat_catalog.py  --  LIBRARY-FREE, read-derived REPEAT DEFINITION via the VARIATION-GRAPH framework.

A novel O1 object parallel to the family definition. No RepBase / Dfam / RepeatMasker enters the
DEFINITION (RepeatMasker soft-mask is used ONLY as external validation at the end).

VG construction  --  canonical-MINIMIZER MULTIPLICITY graph (bag-of-minimizers per gene)
----------------------------------------------------------------------------------------
HONEST NAMING (do not over-read "variation/segment graph"): the delivered object is the documented
ALTERNATIVE of the two the prompt allowed -- a canonical-minimizer / k-mer MULTIPLICITY graph, NOT
the maximal-shared-aligned-SEGMENT graph.  Each gene carries a SET of minimizer nodes (a bag of
minimizers); the family-edge oracle uses only node-SET intersection + minimum multiplicity.  Ordered
PATH structure is used ONLY for the cyclic/tandem flag.  So "path-multiplicity" here = node
set-membership across gene node-sets, and "shared node" = set intersection, not an aligned segment.
 * TRANSCRIPT = spliced exonic sequence of a de-novo locus (DN_<contig>_<start>_<nreads>), fetched
   soft-masked from GGO.fasta at the RNA exon coords (dn_exons).  Universe = the 3462 gene-assigned
   loci that participate in E_r families (F.build_genes_dict), i.e. the family-forming universe.
 * NODE  = a canonical (strand-agnostic) (k,w)-minimizer value  -- a shared sequence anchor.
   Canonical => an inverted paralog shares the same node.
 * A GENE's node SET = union of the minimizer values over its loci (order discarded => a bag/set).
 * PATH-MULTIPLICITY(node) = # DISTINCT GENES whose node set contains that node (set membership,
   NOT k-mer occurrence count; invariant multiplicity <= n_transcripts holds by construction).
 * CYCLE flag(node) = the node's minimizer recurs (>=2 runs) within a single transcript = tandem/cyclic
   (the only place ordered-path structure is consulted).

Library-free REPEAT node  :=  path-multiplicity >= m   (m swept)   OR   cyclic.
Repeat-AWARE edge oracle  :  keep a family edge (A,B) iff the node-SET intersection A&B contains a
                             node with multiplicity < m (a UNIQUE/low-mult minimizer); i.e. MASK
                             repeat nodes before edge formation (pure set intersection, no path).
                             score(A,B) = min multiplicity over shared nodes  (higher => over-merge).

Why this could beat the scalar read-degree: read/homology DEGREE flags the whole LOCUS as repeat-rich,
so gorilla LOC* paralogs (themselves repeat-rich) look like junk.  The graph keeps the pair iff they
share even ONE unique node -- the repeat nodes are masked but the private family node survives.

Evaluation (vs DNA truth = validation only)
 1. AUC(genuine-FP > TP) of the VG min-shared-multiplicity separator, vs core_recip / soft-mask /
    aln_frac baselines on the SAME rows (edge_bridge_mask.tsv + ri_sharedlen_cache.tsv).
 2. fam17 (27 genes / 16 protein families, the clearest repeat-bridge hub): does the VG flag the
    shared node high-multiplicity and cut the hub?  Do GSTM2 (protein DOMAIN) + MAGE (cardinality)
    correctly REMAIN (honest negative controls -- NOT flagged as repeats)?
 3. RepeatMasker concordance: do high-path-multiplicity nodes correspond to soft-masked genome?
 Component-level 5-fold CV (seed=0) for the masking operating point + recall cost.

Deterministic: 2-bit canonical encoding (no hashing), leftmost-tie minimizers, PYTHONHASHSEED=0,
component-sorted folds.

Outputs: bench/vg_repeat_catalog.tsv  (per-node: multiplicity,cyclic,is_repeat,n_genes,... ; then a
per-edge block: shares_only_repeat) and bench/vg_repeat_catalog.json (summary).
"""
import os, sys, json, math
from collections import defaultdict, deque

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pysam
import family_er_pr as F  # meta/annot/gene projection reuse

BENCH = os.path.dirname(os.path.abspath(__file__))
SCRATCH = "/home/juanfra/winloci_scratch"
GENOME = os.path.join(SCRATCH, "GGO.fasta")
SKEL = os.path.join(SCRATCH, "denovo_skeletons.tsv")
MASK_CACHE = os.path.join(BENCH, "edge_bridge_mask.tsv")     # per gene-pair: cls, core, bridge_mask
SHAREDLEN = os.path.join(BENCH, "ri_sharedlen_cache.tsv")    # per gene-pair: aln_frac (arbitration pop)
LEVELPR = os.path.join(BENCH, "family_level_pr_current.tsv")
OUT_TSV = os.path.join(BENCH, "vg_repeat_catalog.tsv")
OUT_JSON = os.path.join(BENCH, "vg_repeat_catalog.json")

K = 15          # minimizer k
W = 10          # minimizer window (w consecutive k-mers)
M_GRID = [2, 3, 4, 5, 6, 8, 10, 12, 15, 20]   # swept repeat-multiplicity thresholds
M_OP = 5        # headline operating threshold (reported explicitly; CV also selects one)
KFOLDS = 5
SEED = 0

_B = {65: 0, 67: 1, 71: 2, 84: 3, 97: 0, 99: 1, 103: 2, 116: 3}  # A C G T a c g t


# --------------------------------------------------------------- skeleton / exon extraction
def load_skeletons():
    skel = {}
    with open(SKEL) as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            skel[(f[0], int(f[1]), int(f[2]))] = f[5]
    return skel


def dn_exons(dn, meta, skel):
    mm = meta.get(dn)
    if mm is None:
        return None, None
    chrom, start, end = mm
    introns = skel.get((chrom, start, end))
    if introns is None:
        return chrom, [(start, end)]
    iv = []
    if introns.strip():
        for tok in introns.split(";"):
            a, b = tok.split("-")
            iv.append((int(a), int(b)))
    iv.sort()
    exons, cur = [], start
    for a, b in iv:
        if a > cur:
            exons.append((cur, a))
        cur = b
    if end > cur:
        exons.append((cur, end))
    return chrom, exons


# --------------------------------------------------------------- canonical minimizers
def minimizers(seq, k=K, w=W):
    """Yield collapsed runs of selected canonical minimizers over `seq` (case-carrying).
    Returns list of (canon_value, masked_bool); consecutive identical selections merged into one run.
    canon = min(fwd 2-bit code, revcomp code) => strand-agnostic.  N breaks the k-mer window."""
    n = len(seq)
    if n < k:
        return []
    kmask = (1 << (2 * k)) - 1
    shift = 2 * (k - 1)
    b = seq.encode("ascii", "replace")
    # rolling canonical code + rolling lowercase count over the k-window
    canon = []        # canon value per k-mer end position (or None if window invalid)
    masked = []       # majority-lowercase over the k-mer
    f = r = 0
    valid = 0
    lc = deque()      # lowercase flags of current window bases
    lc_sum = 0
    for i in range(n):
        c = b[i]
        code = _B.get(c, -1)
        is_lc = 1 if 97 <= c <= 122 else 0
        if code < 0:
            f = r = 0
            valid = 0
            lc.clear(); lc_sum = 0
            canon.append(None); masked.append(0)
            continue
        f = ((f << 2) | code) & kmask
        r = (r >> 2) | ((3 - code) << shift)
        lc.append(is_lc); lc_sum += is_lc
        if len(lc) > k:
            lc_sum -= lc.popleft()
        valid += 1
        if valid >= k:
            canon.append(f if f < r else r)
            masked.append(1 if lc_sum * 2 > k else 0)
        else:
            canon.append(None); masked.append(0)

    # window-min (monotonic deque) over canon positions that are not None
    out = []
    dq = deque()  # holds indices into canon, increasing value order, leftmost-tie kept
    last_val = None
    positions = [i for i in range(len(canon)) if canon[i] is not None]
    # simple O(n*w) with early skips is fine at this scale; use deque for O(n)
    for idx, p in enumerate(positions):
        v = canon[p]
        while dq and canon[dq[-1]] > v:   # strict => keep leftmost on ties
            dq.pop()
        dq.append(p)
        # window of the last w k-mers => positions[idx-w+1 .. idx]
        if idx >= w - 1:
            left_p = positions[idx - w + 1]
            while dq[0] < left_p:
                dq.popleft()
            sel = dq[0]
            val = canon[sel]
            if val != last_val:               # collapse consecutive identical selections into runs
                out.append((val, masked[sel]))
                last_val = val
    return out


# --------------------------------------------------------------- AUC + stats
def auc_gt(pos, neg):
    """Mann-Whitney AUC = P(pos > neg), 0.5 tie-corrected. None if either side empty."""
    if not pos or not neg:
        return None
    allv = sorted([(v, 0) for v in pos] + [(v, 1) for v in neg])
    # rank-sum
    ranks = {}
    i = 0
    n = len(allv)
    order = sorted(range(n), key=lambda j: allv[j][0])
    j = 0
    rank = [0.0] * n
    while j < n:
        k = j
        while k + 1 < n and allv[order[k + 1]][0] == allv[order[j]][0]:
            k += 1
        avg = (j + k) / 2.0 + 1.0
        for t in range(j, k + 1):
            rank[order[t]] = avg
        j = k + 1
    r_pos = sum(rank[i] for i in range(n) if allv[i][1] == 0)
    npos, nneg = len(pos), len(neg)
    u = r_pos - npos * (npos + 1) / 2.0
    return u / (npos * nneg)


def spearman(xs, ys):
    n = len(xs)
    if n < 3:
        return None
    def rk(a):
        order = sorted(range(n), key=lambda i: a[i])
        r = [0.0] * n
        j = 0
        while j < n:
            k = j
            while k + 1 < n and a[order[k + 1]] == a[order[j]]:
                k += 1
            avg = (j + k) / 2.0 + 1.0
            for t in range(j, k + 1):
                r[order[t]] = avg
            j = k + 1
        return r
    rx, ry = rk(xs), rk(ys)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((rx[i] - mx) * (ry[i] - my) for i in range(n))
    dx = math.sqrt(sum((rx[i] - mx) ** 2 for i in range(n)))
    dy = math.sqrt(sum((ry[i] - my) ** 2 for i in range(n)))
    return None if dx == 0 or dy == 0 else num / (dx * dy)


# --------------------------------------------------------------- main
def main():
    print("[1] gene projection ...", flush=True)
    meta = F.load_meta()
    annot = F.load_annot()
    gene_of = F.gene_of_factory(annot)
    raw_fams = F.load_raw_families()
    genes, gene_of_dn, _, _, _ = F.build_genes_dict(raw_fams, meta, gene_of)
    skel = load_skeletons()
    loci = sorted(gene_of_dn.keys())
    print(f"    {len(loci)} gene-assigned family-universe loci, "
          f"{len(set(gene_of_dn.values()))} distinct genes", flush=True)

    print("[2] build VG (fetch soft-masked spliced seq + canonical minimizers) ...", flush=True)
    fa = pysam.FastaFile(GENOME)
    node_genes = defaultdict(set)     # code -> set(gene)         (=> multiplicity)
    gene_nodes = defaultdict(set)     # gene -> set(code)
    cyc_max = defaultdict(int)        # code -> max #runs within a single transcript
    lc_occ = defaultdict(int)         # code -> masked run occurrences
    tot_occ = defaultdict(int)        # code -> run occurrences
    n_tx = defaultdict(int)           # code -> #transcripts carrying it
    done = 0
    for dn in loci:
        g = gene_of_dn.get(dn)
        if not g:
            continue
        chrom, exons = dn_exons(dn, meta, skel)
        if exons is None:
            continue
        try:
            s = "".join(fa.fetch(chrom, a, b) for a, b in exons)
        except Exception:
            continue
        if len(s) < K:
            continue
        runs = minimizers(s)
        # per-transcript run counts for cyclic flag
        loc_runcount = defaultdict(int)
        seen = set()
        for val, mk in runs:
            node_genes[val].add(g)
            seen.add(val)
            loc_runcount[val] += 1
            tot_occ[val] += 1
            lc_occ[val] += mk
        for val in seen:
            gene_nodes[g].add(val)
            n_tx[val] += 1
            if loc_runcount[val] > cyc_max[val]:
                cyc_max[val] = loc_runcount[val]
        done += 1
        if done % 500 == 0:
            print(f"    {done}/{len(loci)} loci, {len(node_genes)} nodes so far", flush=True)
    fa.close()

    mult = {c: len(gs) for c, gs in node_genes.items()}
    n_nodes = len(mult)
    n_shared = sum(1 for v in mult.values() if v >= 2)
    n_cyclic = sum(1 for c in mult if cyc_max[c] >= 2)
    print(f"    VG: {n_nodes} nodes, {n_shared} shared (mult>=2), {n_cyclic} cyclic", flush=True)

    # ---------------------------------------------------------- eval population
    print("[3] load eval population (edge_bridge_mask + sharedlen) ...", flush=True)
    CLS = {"TP": "TP", "genuine": "genuine", "truthbar": "truthbar"}
    pairs = {}   # frozenset(a,b) -> dict(cls, core, mask)
    with open(MASK_CACHE) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {h: i for i, h in enumerate(hdr)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            a, b = f[ix["gene_a"]], f[ix["gene_b"]]
            cls = f[ix["cls"]]
            core = f[ix["core"]]
            mask = f[ix["bridge_mask"]]
            pairs[frozenset((a, b))] = dict(
                a=a, b=b, cls=cls,
                core=(None if core == "" else float(core)),
                mask=(None if mask == "" else float(mask)))
    aln = {}
    with open(SHAREDLEN) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {h: i for i, h in enumerate(hdr)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            af = f[ix["aln_frac"]]
            aln[frozenset((f[ix["gene_a"]], f[ix["gene_b"]]))] = (None if af == "" else float(af))

    # VG score per eval pair
    for k, d in pairs.items():
        a, b = d["a"], d["b"]
        na, nb = gene_nodes.get(a), gene_nodes.get(b)
        if not na or not nb:
            d["shared"] = None
            continue
        shared = na & nb
        if not shared:
            d["shared"] = 0
            d["min_shared_mult"] = None
            d["max_shared_mult"] = None
            d["n_shared_repeat"] = None
            d["n_shared_unique"] = None
            continue
        sm = [mult[c] for c in shared]
        d["shared"] = len(shared)
        d["min_shared_mult"] = min(sm)
        d["max_shared_mult"] = max(sm)
        d["aln_frac"] = aln.get(k)

    # ---------------------------------------------------------- AUC comparison
    print("[4] AUC(genuine>TP) VG vs baselines ...", flush=True)

    def collect(feat, need_vg=True, cls_pos="genuine", cls_neg="TP"):
        pos, neg = [], []
        for d in pairs.values():
            if need_vg and d.get("shared") in (None, 0):
                continue
            v = feat(d)
            if v is None:
                continue
            if d["cls"] == cls_pos:
                pos.append(v)
            elif d["cls"] == cls_neg:
                neg.append(v)
        return pos, neg

    results = {}
    # VG min-shared-multiplicity: HIGHER => genuine over-merge
    p, n = collect(lambda d: d.get("min_shared_mult"))
    results["VG_min_shared_mult"] = dict(auc=auc_gt(p, n), n_genuine=len(p), n_TP=len(n),
                                         note="higher=genuine; primary VG separator")
    # VG max-shared-multiplicity
    p, n = collect(lambda d: d.get("max_shared_mult"))
    results["VG_max_shared_mult"] = dict(auc=auc_gt(p, n), n_genuine=len(p), n_TP=len(n))
    # baselines on the SAME VG-covered rows (matched anchor)
    p, n = collect(lambda d: (-d["core"]) if d.get("core") is not None else None)  # lower core=genuine
    results["core_recip_matched"] = dict(auc=auc_gt(p, n), n_genuine=len(p), n_TP=len(n),
                                         note="AUC on VG-covered rows; sign: lower core=genuine")
    p, n = collect(lambda d: d["mask"] if d.get("mask") is not None else None)     # higher mask=genuine
    results["softmask_matched"] = dict(auc=auc_gt(p, n), n_genuine=len(p), n_TP=len(n),
                                       note="higher softmask=genuine")
    p, n = collect(lambda d: (-d["aln_frac"]) if d.get("aln_frac") is not None else None)
    results["aln_frac_matched"] = dict(auc=auc_gt(p, n), n_genuine=len(p), n_TP=len(n),
                                       note="lower aln_frac=genuine (matched, sparse coverage)")
    # baselines on ALL rows that have the feature (unmatched, full population)
    p, n = collect(lambda d: (-d["core"]) if d.get("core") is not None else None, need_vg=False)
    results["core_recip_full"] = dict(auc=auc_gt(p, n), n_genuine=len(p), n_TP=len(n))
    p, n = collect(lambda d: d["mask"] if d.get("mask") is not None else None, need_vg=False)
    results["softmask_full"] = dict(auc=auc_gt(p, n), n_genuine=len(p), n_TP=len(n))

    # coverage
    cov = defaultdict(lambda: [0, 0])
    for d in pairs.values():
        has = d.get("shared") not in (None, 0)
        cov[d["cls"]][0] += 1
        cov[d["cls"]][1] += 1 if has else 0
    coverage = {c: dict(n=v[0], vg_covered=v[1]) for c, v in cov.items()}

    # ---------------------------------------------------------- masking oracle sweep + CV
    print("[5] repeat-masking oracle sweep + component CV ...", flush=True)

    def oracle_keep(d, m):
        """keep iff shares a node with multiplicity < m (a unique/low-mult segment)."""
        if d.get("shared") in (None, 0):
            return False  # no shared node => cut
        return d["min_shared_mult"] < m

    sweep = {}
    for m in M_GRID:
        tp_keep = tp_all = gen_keep = gen_all = tb_keep = tb_all = 0
        for d in pairs.values():
            keep = oracle_keep(d, m)
            if d["cls"] == "TP":
                tp_all += 1; tp_keep += keep
            elif d["cls"] == "genuine":
                gen_all += 1; gen_keep += keep
            elif d["cls"] == "truthbar":
                tb_all += 1; tb_keep += keep
        sweep[m] = dict(
            tp_recall=round(tp_keep / tp_all, 4) if tp_all else None,
            genuine_kept=gen_keep, genuine_cut=gen_all - gen_keep,
            genuine_cut_rate=round((gen_all - gen_keep) / gen_all, 4) if gen_all else None,
            truthbar_recall=round(tb_keep / tb_all, 4) if tb_all else None,
            tp_all=tp_all, genuine_all=gen_all, truthbar_all=tb_all)

    # component-level 5-fold CV (group genes into connected components of in_er edges)
    adj = defaultdict(set)
    for d in pairs.values():
        adj[d["a"]].add(d["b"]); adj[d["b"]].add(d["a"])
    seen = set(); comp_of = {}
    cid = 0
    for g in sorted(adj):
        if g in seen:
            continue
        stack = [g]; seen.add(g)
        while stack:
            x = stack.pop(); comp_of[x] = cid
            for y in adj[x]:
                if y not in seen:
                    seen.add(y); stack.append(y)
        cid += 1
    comps = sorted(set(comp_of.values()))
    fold_of = {c: (i % KFOLDS) for i, c in enumerate(comps)}  # seed=0 deterministic (sorted)

    def pair_fold(d):
        return fold_of[comp_of[d["a"]]]

    cv = []
    for test in range(KFOLDS):
        # select m on TRAIN maximizing genuine_cut_rate - (1 - tp_recall) (net over-merge removal)
        best_m, best_obj = None, -1e9
        for m in M_GRID:
            tpk = tpa = gk = ga = 0
            for d in pairs.values():
                if pair_fold(d) == test:
                    continue
                keep = oracle_keep(d, m)
                if d["cls"] == "TP":
                    tpa += 1; tpk += keep
                elif d["cls"] == "genuine":
                    ga += 1; gk += keep
            if tpa == 0 or ga == 0:
                continue
            tp_rec = tpk / tpa
            gen_cut = (ga - gk) / ga
            obj = gen_cut - (1 - tp_rec)
            if obj > best_obj:
                best_obj, best_m = obj, m
        # evaluate on TEST at best_m + test-fold AUC of VG separator
        tpk = tpa = gk = ga = 0
        pos, neg = [], []
        for d in pairs.values():
            if pair_fold(d) != test:
                continue
            keep = oracle_keep(d, best_m)
            if d["cls"] == "TP":
                tpa += 1; tpk += keep
                if d.get("min_shared_mult") is not None:
                    neg.append(d["min_shared_mult"])
            elif d["cls"] == "genuine":
                ga += 1; gk += keep
                if d.get("min_shared_mult") is not None:
                    pos.append(d["min_shared_mult"])
        cv.append(dict(fold=test, selected_m=best_m,
                       tp_recall=round(tpk / tpa, 4) if tpa else None,
                       genuine_cut_rate=round((ga - gk) / ga, 4) if ga else None,
                       test_auc=auc_gt(pos, neg), tp_all=tpa, genuine_all=ga))
    cv_aucs = [c["test_auc"] for c in cv if c["test_auc"] is not None]
    cv_tprec = [c["tp_recall"] for c in cv if c["tp_recall"] is not None]
    cv_gcut = [c["genuine_cut_rate"] for c in cv if c["genuine_cut_rate"] is not None]
    cv_summary = dict(
        folds=cv, n_components=len(comps),
        mean_test_auc=round(sum(cv_aucs) / len(cv_aucs), 4) if cv_aucs else None,
        mean_tp_recall=round(sum(cv_tprec) / len(cv_tprec), 4) if cv_tprec else None,
        mean_genuine_cut_rate=round(sum(cv_gcut) / len(cv_gcut), 4) if cv_gcut else None)

    # ---------------------------------------------------------- fam17 / GSTM2 / MAGE
    print("[6] fam17 / GSTM2 / MAGE case analysis ...", flush=True)
    fam17_genes = ["C9H11orf65", "COPS9", "DCAF16", "DCP1B", "GIMAP4", "KLHL33", "LIPA",
                   "LOC109024560", "LOC109025730", "LOC109028669", "LOC129530096", "LOC129530131",
                   "LOC134756589", "LRFN5", "NCALD", "NCR3LG1", "NFXL1", "OXCT1", "PDCD5",
                   "SLC11A2", "SMYD3", "SPATA33", "SPC25", "THNSL1", "TMEM38B", "TSPAN8", "UCHL1"]
    gstm2_genes = ["GSTM2", "LOC115930164", "LOC115930576", "LOC101126097", "LOC134756922",
                   "LOC101129940", "SEC22B", "LOC109023809", "TRIP13", "LOC129532045", "LOC134757399",
                   "LOC101131274", "LOC101135585", "LOC115932984", "LOC115933235", "LOC115933241",
                   "LOC129525330", "LOC129525331", "LOC129525599", "LOC134756861"]
    mage_genes = ["MAGEA1", "MAGEA4", "MAGEA12", "LOC129529976", "LOC129529978", "LOC129529983",
                  "LOC129529986", "MAGEA9", "LOC129530018", "MAGEB6", "MAGEB6B"]

    def case_edges(gene_list, label):
        """all VG-covered in_er edges internal to gene_list: min_shared_mult, cut@M_OP, softmask."""
        gs = set(gene_list)
        rows = []
        for d in pairs.values():
            if d["a"] in gs and d["b"] in gs:
                rows.append(d)
        edges = []
        cut = flagged_repeat = 0
        msms = []
        node_mults = []
        for d in rows:
            msm = d.get("min_shared_mult")
            keep = oracle_keep(d, M_OP)
            edges.append(dict(a=d["a"], b=d["b"], cls=d["cls"],
                              min_shared_mult=msm, max_shared_mult=d.get("max_shared_mult"),
                              softmask=d.get("mask"), cut_at_Mop=(not keep)))
            if not keep:
                cut += 1
            if msm is not None:
                msms.append(msm); node_mults.append(msm)
        return dict(label=label, n_genes_present=len(gs & set(gene_nodes.keys())),
                    n_edges=len(edges), n_cut_at_Mop=cut,
                    frac_cut=round(cut / len(edges), 4) if edges else None,
                    median_min_shared_mult=(sorted(msms)[len(msms) // 2] if msms else None),
                    max_min_shared_mult=(max(msms) if msms else None),
                    edges=sorted(edges, key=lambda e: -(e["min_shared_mult"] or 0))[:40])

    fam17 = case_edges(fam17_genes, "fam17_repeat_bridge_hub")
    gstm2 = case_edges(gstm2_genes, "GSTM2_protein_domain")
    mage = case_edges(mage_genes, "MAGE_cardinality")

    # softmask of the shared nodes in each case (repeat should be softmasked; domain should not)
    def case_node_softmask(gene_list):
        gs = set(gene_list)
        shared_nodes = set()
        gl = [g for g in gene_list if g in gene_nodes]
        for i in range(len(gl)):
            for j in range(i + 1, len(gl)):
                shared_nodes |= (gene_nodes[gl[i]] & gene_nodes[gl[j]])
        if not shared_nodes:
            return None
        # among shared nodes, the high-mult ones (the bridge) -- report their softmask
        hi = [c for c in shared_nodes if mult[c] >= M_OP]
        sm_hi = [lc_occ[c] / tot_occ[c] for c in hi if tot_occ[c] > 0]
        sm_all = [lc_occ[c] / tot_occ[c] for c in shared_nodes if tot_occ[c] > 0]
        return dict(n_shared_nodes=len(shared_nodes),
                    n_hi_mult=len(hi),
                    max_mult=max((mult[c] for c in shared_nodes), default=0),
                    mean_softmask_hi=round(sum(sm_hi) / len(sm_hi), 4) if sm_hi else None,
                    mean_softmask_all=round(sum(sm_all) / len(sm_all), 4) if sm_all else None)
    fam17["node_softmask"] = case_node_softmask(fam17_genes)
    gstm2["node_softmask"] = case_node_softmask(gstm2_genes)
    mage["node_softmask"] = case_node_softmask(mage_genes)

    # ---------------------------------------------------------- RepeatMasker concordance
    print("[7] RepeatMasker soft-mask concordance ...", flush=True)
    # HONEST CAVEAT (validation only; softmask NOT in the definition): the extreme high-mult tail is
    # a MIX -- some of the very top nodes are LOW-COMPLEXITY minimizers (2-bit code 0/2/8 = poly-A
    # runs AAAAAAAAAAAAAAA...) which are trivially both high-mult and soft-masked, while most others
    # are genuine interspersed-TE fragments (canonical AluY consensus motifs AAAGTGCTGGGATTA /
    # AATCCCAGCACTTTG / AGGCTGAGGCAGGAG).  So the 93% high-tail concordance reflects BOTH low-complexity
    # AND real interspersed TE rediscovered library-free -- not purely TE, not purely low-complexity.
    # node-level: does high multiplicity predict soft-masked?
    node_sm = {c: (lc_occ[c] / tot_occ[c]) for c in mult if tot_occ[c] > 0}
    # AUC: multiplicity predicting "node is soft-masked (>=0.5)"; on shared nodes (mult>=2)
    shared_codes = [c for c in mult if mult[c] >= 2]
    pos = [mult[c] for c in shared_codes if node_sm.get(c, 0) >= 0.5]   # softmasked nodes
    neg = [mult[c] for c in shared_codes if node_sm.get(c, 0) < 0.5]
    rm_auc = auc_gt(pos, neg)
    # of repeat nodes (mult>=M_OP), what fraction are majority soft-masked?
    for m in [2, 3, 5, 10]:
        pass
    concord = {}
    for m in [2, 3, M_OP, 10, 20]:
        rep = [c for c in mult if mult[c] >= m]
        if not rep:
            concord[m] = None
            continue
        smk = sum(1 for c in rep if node_sm.get(c, 0) >= 0.5)
        concord[m] = dict(n_repeat_nodes=len(rep),
                          frac_softmasked=round(smk / len(rep), 4),
                          mean_softmask=round(sum(node_sm.get(c, 0) for c in rep) / len(rep), 4))
    # baseline: fraction of ALL shared nodes soft-masked, and low-mult (mult==2) soft-masked
    low = [c for c in mult if mult[c] == 2]
    low_sm = round(sum(1 for c in low if node_sm.get(c, 0) >= 0.5) / len(low), 4) if low else None
    # spearman multiplicity vs softmask over shared nodes (subsample for speed if huge)
    sc = shared_codes
    if len(sc) > 40000:
        step = len(sc) // 40000 + 1
        sc = sc[::step]
    sp = spearman([mult[c] for c in sc], [node_sm.get(c, 0) for c in sc])
    rm = dict(auc_mult_predicts_softmasked=rm_auc,
              spearman_mult_softmask=(round(sp, 4) if sp is not None else None),
              concordance_by_threshold=concord,
              frac_softmasked_mult2_baseline=low_sm,
              n_shared_nodes=len(shared_codes))

    # ---------------------------------------------------------- write outputs
    print("[8] write outputs ...", flush=True)
    with open(OUT_TSV, "w") as fh:
        fh.write("# SECTION nodes (mult>=2)\n")
        fh.write("node_id\tmultiplicity\tn_genes\tcyclic\tis_repeat_m%d\tn_transcripts\tsoftmask_frac\n" % M_OP)
        for c in sorted((c for c in mult if mult[c] >= 2), key=lambda c: (-mult[c], c)):
            isrep = 1 if (mult[c] >= M_OP or cyc_max[c] >= 2) else 0
            fh.write(f"{c}\t{mult[c]}\t{mult[c]}\t{1 if cyc_max[c]>=2 else 0}\t{isrep}\t{n_tx[c]}\t"
                     f"{round(node_sm.get(c,0.0),4)}\n")
        fh.write("# SECTION edges (eval population)\n")
        fh.write("gene_a\tgene_b\tcls\tn_shared_nodes\tmin_shared_mult\tmax_shared_mult\t"
                 "shares_only_repeat_m%d\tkeep_m%d\tcore\tsoftmask\taln_frac\n" % (M_OP, M_OP))
        for k in sorted(pairs, key=lambda k: (pairs[k]["cls"], pairs[k]["a"], pairs[k]["b"])):
            d = pairs[k]
            sh = d.get("shared")
            msm = d.get("min_shared_mult")
            xsm = d.get("max_shared_mult")
            only_rep = "" if sh in (None, 0) else (1 if msm >= M_OP else 0)
            keep = 1 if oracle_keep(d, M_OP) else 0
            fh.write(f"{d['a']}\t{d['b']}\t{d['cls']}\t{sh if sh is not None else ''}\t"
                     f"{msm if msm is not None else ''}\t{xsm if xsm is not None else ''}\t"
                     f"{only_rep}\t{keep}\t{d.get('core','')}\t{d.get('mask','')}\t"
                     f"{d.get('aln_frac','')}\n")

    summary = dict(
        params=dict(K=K, W=W, M_GRID=M_GRID, M_OP=M_OP, KFOLDS=KFOLDS, SEED=SEED,
                    PYTHONHASHSEED=os.environ.get("PYTHONHASHSEED"),
                    canonical_minimizers=True, softmask_is_validation_only=True),
        universe=dict(n_loci=len(loci), n_genes=len(set(gene_of_dn.values())),
                      n_nodes=n_nodes, n_shared_nodes=n_shared, n_cyclic_nodes=n_cyclic),
        eval_coverage=coverage,
        auc_comparison=results,
        reference_baselines_prior_work=dict(
            aln_frac=0.915, core_recip=0.815, homology_graph_degree=0.83,
            read_multimap_degree="0.64-0.69", genome_softmask=0.69,
            note="prior bench/FAMILY_ER_*; matched anchors recomputed above on VG-covered rows"),
        masking_oracle_sweep=sweep,
        cv=cv_summary,
        fam17=fam17, gstm2=gstm2, mage=mage,
        repeatmasker_concordance=rm)
    with open(OUT_JSON, "w") as fh:
        json.dump(summary, fh, indent=1, default=str)

    # ---------------------------------------------------------- console report
    print("\n================= VG REPEAT CATALOG =================")
    print(f"universe: {len(loci)} loci / {len(set(gene_of_dn.values()))} genes; "
          f"nodes={n_nodes} shared={n_shared} cyclic={n_cyclic}")
    print(f"eval coverage (VG shared-node present): "
          + ", ".join(f"{c}={coverage[c]['vg_covered']}/{coverage[c]['n']}" for c in ("TP", "genuine", "truthbar")))
    print("\nAUC(genuine-FP > TP)  [higher => cleaner over-merge separator]")
    for name in ["VG_min_shared_mult", "VG_max_shared_mult", "core_recip_matched",
                 "softmask_matched", "aln_frac_matched", "core_recip_full", "softmask_full"]:
        r = results[name]
        print(f"  {name:22s} AUC={r['auc']}  (genuine={r['n_genuine']} TP={r['n_TP']})")
    print("  reference prior-work: aln_frac=0.915 core=0.815 homology_deg=0.83 readmm=0.64-0.69 softmask=0.69")
    print(f"\nCV(seed=0, {cv_summary['n_components']} comps, {KFOLDS}-fold): "
          f"mean_test_auc={cv_summary['mean_test_auc']} "
          f"tp_recall={cv_summary['mean_tp_recall']} genuine_cut={cv_summary['mean_genuine_cut_rate']}")
    print("\nmasking oracle sweep (keep iff shares node mult<m):")
    for m in M_GRID:
        s = sweep[m]
        print(f"  m={m:2d}  TP_recall={s['tp_recall']}  genuine_cut={s['genuine_cut_rate']}  "
              f"truthbar_recall={s['truthbar_recall']}")
    print(f"\nfam17 hub: {fam17['n_edges']} edges, cut@m{M_OP}={fam17['n_cut_at_Mop']} "
          f"(frac {fam17['frac_cut']}), median_min_shared_mult={fam17['median_min_shared_mult']}, "
          f"max_min_shared_mult={fam17['max_min_shared_mult']}, node_softmask={fam17['node_softmask']}")
    print(f"GSTM2    : {gstm2['n_edges']} edges, cut@m{M_OP}={gstm2['n_cut_at_Mop']} "
          f"(frac {gstm2['frac_cut']}), median_min_shared_mult={gstm2['median_min_shared_mult']}, "
          f"node_softmask={gstm2['node_softmask']}")
    print(f"MAGE     : {mage['n_edges']} edges, cut@m{M_OP}={mage['n_cut_at_Mop']} "
          f"(frac {mage['frac_cut']}), median_min_shared_mult={mage['median_min_shared_mult']}, "
          f"node_softmask={mage['node_softmask']}")
    print(f"\nRepeatMasker concordance: AUC(mult->softmasked)={rm['auc_mult_predicts_softmasked']}, "
          f"spearman={rm['spearman_mult_softmask']}, mult2_baseline_softmask={rm['frac_softmasked_mult2_baseline']}")
    for m in [2, 3, M_OP, 10, 20]:
        print(f"  repeat nodes mult>={m}: {rm['concordance_by_threshold'][m]}")
    print("====================================================")
    print(f"wrote {OUT_TSV}\nwrote {OUT_JSON}")


if __name__ == "__main__":
    main()
