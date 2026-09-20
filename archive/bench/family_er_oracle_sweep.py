#!/usr/bin/env python
"""Sweep the reframed RNA family (E_r) EDGE ORACLE to find the operating point that KEEPS the
divergent-paralog recall (the E_r win) while CUTTING the genuine over-merges (the AMY2A-ZNF91 type FP).

Context (bench/FAMILY_ER_PRECISION_RECALL.md, family_er_pr.py): the reframed E_r is the refined
gamma-quasi-clique gene-pair set (10,755 edges) = TP 450 (in cDNA-90% truth) + truthbar 6,142
(protein-homologous divergent paralogs BELOW the cDNA-90% bar -- the E_r WIN, the truth is incomplete)
+ genuine 4,163 (real over-merges: AMY2A-ZNF91, RPL14-ZNF669 -- unrelated genes bridged by the loose
0.13-core POA homology). This script gates those 10,755 edges by three oracles and reports the
precision/recall trade-off, to pick a concrete operating point.

TWO HONESTY CAVEATS baked into every number below (surfaced after three adversarial passes)
  (A) TAUTOLOGY. family_er_pr.classify_fp DEFINES truthbar-FP := protein_homologous AND id_ok AND cov_ok,
      so `truthbar` is a SUBSET of `in_ep` BY CONSTRUCTION (verified: 0 violations). Therefore "the protein
      route preserves 100% of the divergent-paralog recall at zero cost" is DEFINITIONALLY FORCED, not an
      empirical discovery. Read it as "zero cost to the PROTEIN-HOMOLOGOUS PROXY of the win." Any
      NON-protein-homologous divergent paralog is definitionally binned `genuine` and CUT; the 97
      protein-missed TP (real per the cDNA truth, yet in_ep=0) PROVE that category is non-empty, so the
      `genuine` count is an UPPER BOUND on true over-merges and the recall of non-protein paralogs is
      UNMEASURED.
  (B) BLOCK-RATE has two faces, both reported side-by-side, never just the flattering one:
        genuine_rate_block            = PRFAM-mixing impurity (a block spans >=2 distinct non-mega PRFAM);
                                        matches the FAMILY_ER_PRECISION_RECALL.md 13.4%/15.9% MD headline
                                        but is BLIND to LOC*/no-PRFAM bridges.
        genuine_rate_block_edgehost   = a block still co-hosts both endpoints of a surviving genuine edge;
                                        the DIRECT family-contamination measure. ~85-88% of residual
                                        genuine over-merges involve a LOC*/no-PRFAM gene invisible to the
                                        prfam-mixing metric, so edgehost is 3-4x higher. Cite BOTH.

GATE ORACLES (post-refine gene-pair filters over the fixed 10,755 refined edges)
  * none              : keep all (reference; reproduces the quick-check baseline).
  * protein (in_ep)   : keep iff the pair carries a reciprocal whole-protein mmseqs edge OR same
                        non-mega PRFAM (family_er_pr.protein_homologous). Crushes over-merge (genuine
                        4163->11) but COSTS 97 TP (450->353): non-coding / protein-missed real families.
  * id>=X             : keep iff cDNA nt id >= X (None fails). The WRONG knob: destroys divergent-paralog
                        recall (truthbar has no cDNA id row -> all dropped).
  * core>=t           : keep iff the pair's transcript-homology core_recip >= t. core_recip is a DN-LOCUS
                        edge property (bench/denovo_family_edges.tsv); a gene pair's core = the MAX core
                        over its direct co-blocked DN edges (transitive-only pairs have NO direct edge ->
                        core=None -> fail any t, they carry no direct homology evidence).
  * combined (in_ep OR core>=t_high)  : THE RECOMMENDED FAMILY. Protein confirms coding paralogs; a high
                        transcript-core confirms the non-coding / protein-missed real families the pure
                        protein gate drops; the unrelated bridges (no protein AND low/no core) are cut.

KNEE (task-literal): max (tp_retained + truthbar_retained) s.t. genuine_rate_block (PRFAM-mixing) <= a
small target (~1-3%). The block<=3% knee is combined t_high=0.70 (block 2.75%). We RECOMMEND that point.
A looser EDGE-genuine-rate<=3% budget lands on combined t=0.50 (block 6.15%, >2x the 3% target); it is
reported ONLY as an explicitly-labeled recovery-favoring ALTERNATIVE, never as "the knee."

For each operating point: n_edges, tp_retained, tp_noep_retained (of the 97 protein-missed TP -- the
recovery target), noncoding_tp_retained (of the 31 real non-coding TP, ALL protein-missed -- the honest
non-coding recall, replacing the inverted block-dominant count), truthbar_retained (the protein-homologous
divergent-paralog proxy = the E_r win), genuine_overmerge, genuine_rate_edge, BOTH block rates,
noncoding_dominant_blocks (+ how many host a REAL TP edge -- the rest are genuine-only over-merges), and
id-stratified recall on the reachable dna_loose subset.

Outputs: bench/family_er_oracle_sweep.tsv (one row per operating point), bench/family_er_oracle_sweep.json
(the full curve + the recommended point + its numbers + named survivors/casualties + the tautology note +
the disproportionate-non-coding-drop Fisher test).
Reproduce: PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/family_er_oracle_sweep.py
Byte-deterministic across launches (hash seed pinned in-script; gene-graph refine witness pinned seed=0).
"""
import json
import os
import sys
from collections import defaultdict, Counter

# Pin string-hash order so the louvain gene-graph re-refine yields a byte-identical witness partition
# across process launches (same discipline as family_er_pr.py).
if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

BENCH = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BENCH)
import family_er_pr as F           # loaders, projection, FP classifier, protein_homologous
import genome_family_def as G      # refine_families (gamma-quasi-clique), UF, GAMMA, SEED

DENOVO_EDGES = os.path.join(BENCH, "denovo_family_edges.tsv")
T_GRID = [0.13, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80]
BLOCK_TARGETS = [0.01, 0.02, 0.03]   # genuine block-over-merge (PRFAM-mixing) targets for the task-literal knee
NONCODING_BT = ("protein_coding", "NA", "?")   # a gene is "non-coding" iff biotype not in this set


def is_noncoding(g, g2bt):
    bt = g2bt.get(g)
    return bt is not None and bt not in NONCODING_BT


# ---------------------------------------------------------------- data prep
def prepare():
    meta = F.load_meta()
    annot = F.load_annot()
    gene_of = F.gene_of_factory(annot)
    dna_loose, ec_rna, rowinfo = F.load_truth()
    ep_direct = F.load_ep_direct()
    g2f, mega = F.load_prfam()
    raw_fams = F.load_raw_families()
    edges = F.load_edges()
    genes, gene_of_dn, gene_of_dn_raw, miss, ov = F.build_genes_dict(raw_fams, meta, gene_of)

    blocks0 = G.refine_families(raw_fams, edges, genes, gamma=G.GAMMA, seed=0)
    ref_pairs, ref_repr = F.blocks_to_gene_pairs(blocks0, gene_of_dn)

    # gene -> biotype (same source as family_er_pr: the DN-locus max-overlap projection over annot_intervals,
    # which is derived from GGO_genomic.gff)
    g2bt = {}
    for dn, g in gene_of_dn.items():
        g2bt.setdefault(g, genes[dn]["biotype"])

    # gene -> representative coord (longest annot interval) for the gene-graph gamma re-refine
    g2coord = {}
    for c, (starts, ivals, maxlen) in annot.items():
        for s, e, bt, g in ivals:
            if g in ref_repr:
                if g not in g2coord or (e - s) > (g2coord[g]["end"] - g2coord[g]["start"]):
                    g2coord[g] = dict(chrom=c, start=s, end=e, name=g, biotype=bt)

    # DN block membership -> gene-pair core_recip = MAX core over direct CO-BLOCKED DN edges
    dn2blk = {}
    for bi, blk in enumerate(blocks0):
        for dn in blk:
            dn2blk[dn] = bi
    gp_core = {}
    with open(DENOVO_EDGES) as fh:
        next(fh)
        for ln in fh:
            a, b, cc = ln.rstrip("\n").split("\t")
            cc = float(cc)
            ga, gb = gene_of_dn.get(a), gene_of_dn.get(b)
            if ga is None or gb is None or ga == gb:
                continue
            if dn2blk.get(a) is None or dn2blk.get(a) != dn2blk.get(b):
                continue
            k = frozenset((ga, gb))
            if k not in gp_core or cc > gp_core[k]:
                gp_core[k] = cc

    # per-pair attributes over the fixed refined set
    attr = {}
    for k in ref_pairs:
        if k in dna_loose:
            cls = "TP"
        else:
            c, _ = F.classify_fp(k, ep_direct, g2f, mega, rowinfo)
            cls = "truthbar" if c == "truthbar-FP" else "genuine"
        inep = F.protein_homologous(k, ep_direct, g2f, mega)[0]
        ri = rowinfo.get(k)
        idv = ri["id"] if ri else None
        attr[k] = dict(cls=cls, in_ep=inep, core=gp_core.get(k), id=idv)

    # (A) TAUTOLOGY verification: truthbar is a SUBSET of in_ep by construction (classify_fp requires
    # protein_homologous). Any violation would be a bug; 0 violations is DEFINITIONAL, not empirical.
    tautology_violations = [k for k in ref_pairs if attr[k]["cls"] == "truthbar" and not attr[k]["in_ep"]]
    assert not tautology_violations, f"truthbar NOT subset of in_ep: {len(tautology_violations)}"

    # REAL non-coding TP (pair-level): a cDNA-truth pair with >=1 non-coding gene. This REPLACES the
    # inverted block-dominant "noncoding_families_retained" count (which was ~75% genuine over-merge blocks).
    # ALL of these are protein-missed (in_ep=0) -> recoverable ONLY via the core route.
    noncoding_tp = {k for k in ref_pairs if attr[k]["cls"] == "TP"
                    and (is_noncoding(tuple(k)[0], g2bt) or is_noncoding(tuple(k)[1], g2bt))}
    coding_tp = {k for k in ref_pairs if attr[k]["cls"] == "TP" and k not in noncoding_tp}

    # reachable dna_loose subset (expressed & cross_map & both genes in E_r footprint) with id bin
    def strat(idv):
        if idv is None:
            return "na"
        if idv >= 0.99:
            return ">=0.99"
        if idv >= 0.95:
            return "0.95-0.99"
        if idv >= 0.90:
            return "0.90-0.95"
        return "<0.90"
    reach = []
    for k in dna_loose:
        ri = rowinfo.get(k)
        if ri is None or not (ri["expr_a"] == 1 and ri["expr_b"] == 1) or ri["cross_map"] != 1:
            continue
        a, b = tuple(k)
        if a in ref_repr and b in ref_repr:
            reach.append((k, strat(ri["id"])))

    return dict(ref_pairs=ref_pairs, ref_repr=ref_repr, attr=attr, g2bt=g2bt, g2coord=g2coord,
                reach=reach, dna_loose=dna_loose, blocks0=blocks0, g2f=g2f, mega=mega,
                noncoding_tp=noncoding_tp, coding_tp=coding_tp,
                tautology_violations=len(tautology_violations))


# ---------------------------------------------------------------- helpers on residual genuine edges
def prfam_invisible(k, D):
    """A genuine over-merge edge is INVISIBLE to the PRFAM-mixing block metric iff it cannot create a
    block spanning two distinct non-mega PRFAM -- i.e. at least one endpoint has no non-mega PRFAM
    (a LOC*/no-PRFAM/mega-blob gene). ~85-88% of residual over-merges are of this kind, which is why the
    prfam-mixing block rate is 3-4x lower than the direct edge-host block rate."""
    a, b = tuple(k)
    g2f, mega = D["g2f"], D["mega"]
    fa, fb = g2f.get(a), g2f.get(b)
    return not (fa is not None and fb is not None and fa not in mega and fb not in mega and fa != fb)


# ---------------------------------------------------------------- block-level (gamma re-refine)
def block_metrics(kept, attr, D):
    """Re-refine the SURVIVING gene graph into gamma-quasi-clique blocks (seed=0) and measure family-level
    over-merge TWO ways (both reported):
      PRIMARY (prfam_impure)  -- a block spans >=2 distinct NON-MEGA PRFAM (the MD 13.4%/15.9% headline;
                                 blind to LOC*/no-PRFAM bridges).
      DIRECT  (edge_host)     -- a block still co-hosts both endpoints of a surviving genuine edge (the
                                 honest family-contamination measure; 3-4x higher).
    Also returns the non-coding-DOMINANT block count and, of those, how many HOST A REAL TP EDGE (the rest
    are genuine-only over-merges -- this is what makes the old block-dominant 'noncoding_families' count
    inverted)."""
    g2coord, g2bt, g2f, mega = D["g2coord"], D["g2bt"], D["g2f"], D["mega"]
    if not kept:
        return dict(n_blocks=0, prfam_impure=0, edge_host=0, nc_dom=0, nc_dom_host_tp=0)
    uf = G.UF()
    for k in kept:
        a, b = tuple(k)
        uf.union(a, b)
    comps = [sorted(set(v)) for v in uf.groups().values() if len(set(v)) >= 2]
    ge = [(tuple(k)[0], tuple(k)[1]) for k in kept]

    # gamma re-refine (R): split each component into gamma-quasi-cliques, >=2 loci
    blocks = G.refine_families(comps, ge, g2coord, gamma=G.GAMMA, seed=0)
    gene_blk = {}
    for bi, blk in enumerate(blocks):
        for g in blk:
            gene_blk[g] = bi

    # blocks that host a retained TP edge (both endpoints in the same block)
    blk_tp = set()
    for k in kept:
        if attr[k]["cls"] == "TP":
            a, b = tuple(k)
            if gene_blk.get(a) is not None and gene_blk.get(a) == gene_blk.get(b):
                blk_tp.add(gene_blk[a])

    prfam_impure = 0
    nc_dom = 0
    nc_dom_host_tp = 0
    for bi, blk in enumerate(blocks):
        nonmega = {g2f[g] for g in blk if g in g2f and g2f[g] not in mega}
        if len(nonmega) >= 2:
            prfam_impure += 1
        bts = Counter(g2bt.get(g, "NA") for g in blk)
        dom = bts.most_common(1)[0][0]
        if dom not in NONCODING_BT:
            nc_dom += 1
            if bi in blk_tp:
                nc_dom_host_tp += 1
    edge_host = len({gene_blk[tuple(k)[0]] for k in kept if attr[k]["cls"] == "genuine"
                     and gene_blk.get(tuple(k)[0]) is not None
                     and gene_blk.get(tuple(k)[0]) == gene_blk.get(tuple(k)[1])})
    return dict(n_blocks=len(blocks), prfam_impure=prfam_impure, edge_host=edge_host,
                nc_dom=nc_dom, nc_dom_host_tp=nc_dom_host_tp)


# ---------------------------------------------------------------- one operating point
def compute_point(name, t, keep, D):
    ref_pairs, attr = D["ref_pairs"], D["attr"]
    kept = {k for k in ref_pairs if keep(k, attr[k])}
    n = len(kept)
    tp = sum(1 for k in kept if attr[k]["cls"] == "TP")
    tp_noep = sum(1 for k in kept if attr[k]["cls"] == "TP" and not attr[k]["in_ep"])
    truthbar = sum(1 for k in kept if attr[k]["cls"] == "truthbar")
    genuine = sum(1 for k in kept if attr[k]["cls"] == "genuine")
    genuine_kept = {k for k in kept if attr[k]["cls"] == "genuine"}
    nc_tp_ret = sum(1 for k in kept if k in D["noncoding_tp"])

    bm = block_metrics(kept, attr, D)
    nblk = bm["n_blocks"]
    inv = sum(1 for k in genuine_kept if prfam_invisible(k, D))

    # id-stratified recall on the reachable dna_loose subset (kept as the in_er set)
    bins = ["0.90-0.95", "0.95-0.99", ">=0.99"]
    den = Counter(s for _, s in D["reach"])
    num = Counter(s for k, s in D["reach"] if k in kept)
    rec = {b: (round(num[b] / den[b], 4) if den[b] else None) for b in bins}
    reach_tot_n = sum(den[b] for b in bins) + den["<0.90"] + den["na"]
    reach_tot_hit = sum(1 for k, s in D["reach"] if k in kept)

    return dict(
        t=(None if t is None else round(t, 3)), gate=name, n_edges=n,
        tp_retained=tp, tp_noep_retained=tp_noep,
        noncoding_tp_retained=nc_tp_ret, noncoding_tp_total=len(D["noncoding_tp"]),
        truthbar_retained=truthbar,
        genuine_overmerge=genuine,
        genuine_rate_edge=round(genuine / n, 4) if n else None,
        genuine_prfam_invisible=inv,
        n_blocks=nblk, blocks_overmerged=bm["prfam_impure"],
        genuine_rate_block=round(bm["prfam_impure"] / nblk, 4) if nblk else None,
        genuine_rate_block_edgehost=round(bm["edge_host"] / nblk, 4) if nblk else None,
        noncoding_dominant_blocks=bm["nc_dom"],
        noncoding_dominant_blocks_hosting_tp=bm["nc_dom_host_tp"],
        reach_recall_overall=round(reach_tot_hit / reach_tot_n, 4) if reach_tot_n else None,
        **{f"recall_id_{b}": rec[b] for b in bins},
        recall_den={b: den[b] for b in bins},
    )


# ---------------------------------------------------------------- gate predicates
def g_none(k, a):    return True
def g_protein(k, a): return a["in_ep"]
def core_ok(a, t):   return a["core"] is not None and a["core"] >= t


def predicate_for(gate, t):
    """Single source of truth for every gate -> keep(k, attr) predicate (used by the sweep loop, the knee
    reconstruction, the named survivors/casualties, and the determinism recheck)."""
    if gate == "none":
        return g_none
    if gate == "protein":
        return g_protein
    if gate.startswith("id>="):
        x = float(gate[4:])
        return lambda k, a, x=x: a["id"] is not None and a["id"] >= x
    if gate == "core":
        return lambda k, a, t=t: core_ok(a, t)
    if gate == "combined":
        return lambda k, a, t=t: a["in_ep"] or core_ok(a, t)
    if gate == "protein+core":
        return lambda k, a, t=t: a["in_ep"] and core_ok(a, t)
    raise ValueError(gate)


def main():
    print("[load] preparing data (reusing family_er_pr machinery) ...", flush=True)
    D = prepare()
    ref_pairs, attr = D["ref_pairs"], D["attr"]
    n_tp = sum(1 for k in ref_pairs if attr[k]["cls"] == "TP")
    n_tb = sum(1 for k in ref_pairs if attr[k]["cls"] == "truthbar")
    n_ge = sum(1 for k in ref_pairs if attr[k]["cls"] == "genuine")
    n_tp_noep = sum(1 for k in ref_pairs if attr[k]["cls"] == "TP" and not attr[k]["in_ep"])
    n_nc_tp = len(D["noncoding_tp"])
    print(f"[E_r]  refined gene-pairs={len(ref_pairs)}  TP={n_tp}  truthbar={n_tb}  genuine={n_ge}  "
          f"TP(in_ep=0, protein-missed)={n_tp_noep}  of which non-coding-TP={n_nc_tp} (all protein-missed)",
          flush=True)
    print(f"[tautology] truthbar subset of in_ep verified: {D['tautology_violations']} violations "
          f"(0 is DEFINITIONAL -- classify_fp requires protein_homologous; 'protein preserves 100% of the "
          f"E_r win' is forced, not discovered)", flush=True)

    # operating points: (gate, t)
    points = [("none", None), ("protein", None)]
    for x in (0.50, 0.80, 0.90):
        points.append((f"id>={x:.2f}", x))
    for t in T_GRID:
        points.append(("core", t))
    for t in T_GRID:
        points.append(("combined", t))       # in_ep OR core>=t_high  (the recommended family)
    for t in T_GRID:
        points.append(("protein+core", t))    # in_ep AND core>=t (extra: strict pruning reference)
    rows = [compute_point(gate, t, predicate_for(gate, t), D) for gate, t in points]

    # ---- sanity: reproduce the quick-check ----
    def find(gate, t=None):
        for r in rows:
            if r["gate"] == gate and r["t"] == (None if t is None else round(t, 3)):
                return r
        return None
    qc = find("none")
    assert (qc["n_edges"], qc["tp_retained"], qc["truthbar_retained"], qc["genuine_overmerge"]) == \
        (10755, 450, 6142, 4163), f"none mismatch: {qc}"
    qp = find("protein")
    assert (qp["n_edges"], qp["tp_retained"], qp["truthbar_retained"], qp["genuine_overmerge"]) == \
        (6506, 353, 6142, 11), f"protein mismatch: {qp}"
    qi = find("id>=0.80", 0.80)
    assert (qi["n_edges"], qi["tp_retained"], qi["truthbar_retained"], qi["genuine_overmerge"]) == \
        (478, 450, 9, 19), f"id>=0.80 mismatch: {qi}"
    print("[sanity] quick-check reproduced: none(10755/450/6142/4163) protein(6506/353/6142/11) "
          "id>=0.80(478/450/9/19)  OK", flush=True)

    # ---- KNEE selection (TASK-LITERAL): max (tp+truthbar) s.t. genuine_rate_block (PRFAM-mixing) <= target
    gated = [r for r in rows if r["gate"] not in ("none",) and r["genuine_rate_block"] is not None]
    knees = {}
    for tgt in BLOCK_TARGETS:
        cand = [r for r in gated if r["genuine_rate_block"] <= tgt]
        best = max(cand, key=lambda r: (r["tp_retained"] + r["truthbar_retained"], -r["genuine_overmerge"])) \
            if cand else None
        knees[f"block<= {tgt}"] = best

    # TRANSPARENCY: the combined-only knee per BLOCK target, and per EDGE-genuine-rate target. The EDGE
    # budget is a LOOSER precision constraint that lands on a more recovery-favorable point (combined t=0.50,
    # block 6.15% > 3%); it is reported as an alternative, NOT as the recommendation.
    combined_knees = {}
    for tgt in BLOCK_TARGETS:
        cand = [r for r in rows if r["gate"] == "combined" and r["genuine_rate_block"] is not None
                and r["genuine_rate_block"] <= tgt]
        combined_knees[f"block<= {tgt}"] = max(
            cand, key=lambda r: (r["tp_retained"] + r["truthbar_retained"], -r["genuine_overmerge"])) \
            if cand else None
    combined_edge_knees = {}
    for tgt in (0.01, 0.02, 0.03):
        cand = [r for r in rows if r["gate"] == "combined" and r["genuine_rate_edge"] is not None
                and r["genuine_rate_edge"] <= tgt]
        combined_edge_knees[f"edge<= {tgt}"] = max(
            cand, key=lambda r: (r["tp_retained"] + r["truthbar_retained"], -r["genuine_overmerge"])) \
            if cand else None

    # RECOMMENDATION = the TASK-LITERAL block<=3% knee (combined t_high=0.70). It keeps every
    # protein-homologous divergent paralog (truthbar 6142/6142, tautological), recovers 19/97 protein-missed
    # TP, and holds the PRFAM-mixing block over-merge to 2.75% (<=3% target). This is NOT the max-recovery
    # point; the trade-off (recover the 97 protein-missed TP <-> re-admit genuine over-merges) has no free
    # lunch because both live in overlapping core_recip.
    recommended = knees.get("block<= 0.03") or knees.get("block<= 0.02") or knees.get("block<= 0.01")

    # recovery-of-the-97 table: does 'in_ep OR core>=t_high' recover the 97 protein-missed TP, at what cost?
    protein_row = find("protein")
    recovery_of_97 = []
    for r in rows:
        if r["gate"] == "combined":
            recovery_of_97.append(dict(
                t_high=r["t"], tp_noep_recovered=r["tp_noep_retained"], of=n_tp_noep,
                noncoding_tp_recovered=r["noncoding_tp_retained"], of_noncoding=n_nc_tp,
                tp_total=r["tp_retained"],
                genuine_readmitted=r["genuine_overmerge"] - protein_row["genuine_overmerge"],
                genuine_total=r["genuine_overmerge"], genuine_rate_edge=r["genuine_rate_edge"],
                genuine_rate_block_prfam=r["genuine_rate_block"],
                genuine_rate_block_edgehost=r["genuine_rate_block_edgehost"],
                noncoding_dominant_blocks=r["noncoding_dominant_blocks"],
                noncoding_dominant_blocks_hosting_tp=r["noncoding_dominant_blocks_hosting_tp"]))

    # ---- named survivors / casualties for the recommended point ----
    named = {}
    disproportionate = {}
    if recommended:
        keep = predicate_for(recommended["gate"], recommended["t"])
        kept = {k for k in ref_pairs if keep(k, attr[k])}
        def nm(k):
            return " -- ".join(sorted(k))
        surv_tb = [nm(k) for k in kept if attr[k]["cls"] == "truthbar"][:15]
        surv_tp_core = [nm(k) for k in kept if attr[k]["cls"] == "TP" and not attr[k]["in_ep"]][:15]
        cas_gen = [nm(k) for k in ref_pairs if attr[k]["cls"] == "genuine" and k not in kept][:15]
        lost_tp = [nm(k) for k in ref_pairs if attr[k]["cls"] == "TP" and k not in kept][:20]
        kept_gen = [nm(k) for k in kept if attr[k]["cls"] == "genuine"][:20]
        nc_tp_kept = [nm(k) for k in D["noncoding_tp"] if k in kept]
        nc_tp_lost = [nm(k) for k in D["noncoding_tp"] if k not in kept]
        named = dict(survivors_divergent_paralog=surv_tb,
                     survivors_tp_recovered_via_core=surv_tp_core,
                     casualties_genuine_overmerge_cut=cas_gen,
                     casualties_tp_lost=lost_tp,
                     residual_genuine_kept=kept_gen,
                     noncoding_tp_kept=nc_tp_kept,
                     noncoding_tp_lost=nc_tp_lost)

        # (verdict fix) disproportionate NON-CODING TP drop at the recommended point, Fisher exact.
        from scipy.stats import fisher_exact
        nc_drop = sum(1 for k in D["noncoding_tp"] if k not in kept)
        nc_keep = n_nc_tp - nc_drop
        c_drop = sum(1 for k in D["coding_tp"] if k not in kept)
        c_keep = len(D["coding_tp"]) - c_drop
        odds, p = fisher_exact([[nc_drop, nc_keep], [c_drop, c_keep]])
        disproportionate = dict(
            noncoding_tp=n_nc_tp, noncoding_dropped=nc_drop, noncoding_retained=nc_keep,
            noncoding_drop_rate=round(nc_drop / max(1, n_nc_tp), 4),
            coding_tp=len(D["coding_tp"]), coding_dropped=c_drop, coding_retained=c_keep,
            coding_drop_rate=round(c_drop / max(1, len(D["coding_tp"])), 4),
            fisher_odds_ratio=round(float(odds), 3), fisher_p=float(p),
            fold=round((nc_drop / max(1, n_nc_tp)) / max(1e-9, c_drop / max(1, len(D["coding_tp"]))), 2),
            note=("the recommended precision-first point drops non-coding TP disproportionately; all "
                  "non-coding TP are protein-missed so the core route is their ONLY path and it recovers "
                  "few at high t_high. If non-coding families matter, use the recovery-favoring alternative "
                  "(combined t=0.30, 20/31 non-coding TP) at a higher block over-merge, OR a biotype-aware "
                  "core floor -- an explicit precision-vs-non-coding value judgment."))

    # ---- determinism: recompute the recommended point (incl. the louvain block re-refine) twice ----
    det = True
    if recommended:
        gate, t = recommended["gate"], recommended["t"]
        r2 = compute_point(gate, t, predicate_for(gate, t), D)
        det = (r2["genuine_rate_block"] == recommended["genuine_rate_block"]
               and r2["n_blocks"] == recommended["n_blocks"] and r2["n_edges"] == recommended["n_edges"]
               and r2["genuine_rate_block_edgehost"] == recommended["genuine_rate_block_edgehost"])

    # ============================== write TSV ==============================
    cols = ["t", "gate", "n_edges", "tp_retained", "tp_noep_retained", "noncoding_tp_retained",
            "truthbar_retained", "genuine_overmerge", "genuine_rate_edge", "genuine_prfam_invisible",
            "n_blocks", "blocks_overmerged", "genuine_rate_block", "genuine_rate_block_edgehost",
            "noncoding_dominant_blocks", "noncoding_dominant_blocks_hosting_tp",
            "reach_recall_overall", "recall_id_0.90-0.95", "recall_id_0.95-0.99", "recall_id_>=0.99"]
    with open(os.path.join(BENCH, "family_er_oracle_sweep.tsv"), "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in rows:
            out.write("\t".join("" if r.get(c) is None else str(r.get(c)) for c in cols) + "\n")

    # ============================== write JSON ==============================
    summary = dict(
        params=dict(gamma=G.GAMMA, seed=0, t_grid=T_GRID, block_targets=BLOCK_TARGETS,
                    core_recip="max over direct co-blocked DN edges; transitive-only pair -> None (fails)",
                    noncoding_biotype_excluded=list(NONCODING_BT)),
        baseline=dict(refined_gene_pairs=len(ref_pairs), TP=n_tp, truthbar=n_tb, genuine=n_ge,
                      TP_protein_missed=n_tp_noep, noncoding_TP=n_nc_tp,
                      note="truthbar = protein-homologous divergent paralogs below cDNA-90% = the E_r WIN; "
                           "genuine = real over-merges (unrelated-gene bridges) = the cost; noncoding_TP = "
                           "real cDNA-truth pairs with >=1 non-coding gene (ALL protein-missed)"),
        tautology=dict(
            truthbar_subset_of_in_ep_violations=D["tautology_violations"],
            note=("truthbar-FP := protein_homologous AND id_ok AND cov_ok, so truthbar is a SUBSET of in_ep "
                  "BY CONSTRUCTION. '100% of the E_r win preserved by the protein route at zero cost' is "
                  "therefore TAUTOLOGICAL (0 violations is definitional). Read it as 'zero cost to the "
                  "protein-homologous PROXY of the win.' Non-protein-homologous divergent paralogs are "
                  "definitionally binned 'genuine' and CUT; the 97 protein-missed TP prove that category is "
                  "non-empty, so the 'genuine' count is an UPPER BOUND on true over-merges and the recall of "
                  "non-protein paralogs is UNMEASURED.")),
        block_metric_note=("TWO block over-merge rates are reported per point and BOTH should be cited: "
                           "genuine_rate_block = PRFAM-mixing impurity (matches the MD 13.4%/15.9% headline "
                           "but is blind to LOC*/no-PRFAM bridges); genuine_rate_block_edgehost = a block "
                           "still co-hosts a surviving genuine edge (the direct measure, 3-4x higher). "
                           "genuine_prfam_invisible = the count of residual genuine edges invisible to the "
                           "prfam-mixing metric (~85-88% of them)."),
        curve=rows,
        knees=knees,
        recommended=recommended,
        recommended_note=("TASK-LITERAL knee: max(TP+truthbar) s.t. PRFAM-mixing block over-merge <= 3%. "
                          "= combined (in_ep OR core_recip>=0.70). Keeps every protein-homologous divergent "
                          "paralog, recovers 19/97 protein-missed TP (6/31 non-coding), block over-merge 2.75% "
                          "(prfam-mixing) / 11.2% (edge-host)."),
        combined_knees=combined_knees,
        combined_edge_knees=combined_edge_knees,
        recovery_of_97_protein_missed_TP=recovery_of_97,
        recommended_alternatives=dict(
            precision_first_pure_protein=protein_row,          # block 0.0025, 0 protein-missed / 0 non-coding recovered
            edge_budget_recovery_favoring_t050=find("combined", 0.50),   # recovers 33/97 (11/31 nc) but block 6.15% > 3% target
            noncoding_preserving_t030=find("combined", 0.30)),  # 20/31 non-coding TP but block ~11% (fails precision target)
        recommended_named=named,
        noncoding_disproportionate_drop=disproportionate,
        determinism=dict(recommended_reproducible=det, pythonhashseed=os.environ.get("PYTHONHASHSEED")),
    )
    with open(os.path.join(BENCH, "family_er_oracle_sweep.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

    # ============================== console report ==============================
    def fmt(r):
        return (f"{str(r['t'] if r['t'] is not None else ''):>5} {r['gate']:<13} "
                f"n={r['n_edges']:>6} TP={r['tp_retained']:>4} TPnoEp={r['tp_noep_retained']:>3} "
                f"ncTP={r['noncoding_tp_retained']:>2}/{r['noncoding_tp_total']:<2} "
                f"tbar={r['truthbar_retained']:>5} gen={r['genuine_overmerge']:>5} "
                f"g_edge={r['genuine_rate_edge']} g_blk[prfam/host]={r['genuine_rate_block']}/{r['genuine_rate_block_edgehost']} "
                f"rec[.90/.95/.99]={r['recall_id_0.90-0.95']}/{r['recall_id_0.95-0.99']}/{r['recall_id_>=0.99']}")
    print("\n================ E_r ORACLE SWEEP ================")
    print("  t   gate          n_edges TP  TPnoEp ncTP  tbar  gen   g_edge g_blk[prfam/edgehost] recall.90/.95/.99")
    for r in rows:
        print("  " + fmt(r))
    print("\n================ RECOVERY OF THE 97 PROTEIN-MISSED TP (combined = in_ep OR core>=t_high) ======")
    print("  t_high  TPnoEp_recov  ncTP_recov  genuine(+vs protein)  g_edge  g_blk[prfam/edgehost]  ncDomBlk(hostTP)")
    for x in recovery_of_97:
        print(f"   {x['t_high']:<5}  {x['tp_noep_recovered']:>3}/{x['of']:<3}     "
              f"{x['noncoding_tp_recovered']:>2}/{x['of_noncoding']:<2}     "
              f"{x['genuine_total']:>5} (+{x['genuine_readmitted']:<5})   {x['genuine_rate_edge']}   "
              f"{x['genuine_rate_block_prfam']}/{x['genuine_rate_block_edgehost']}   "
              f"{x['noncoding_dominant_blocks']}({x['noncoding_dominant_blocks_hosting_tp']})")
    print("\n================ KNEES (task-literal: PRFAM-mixing block over-merge <= target) ================")
    for tgt, r in knees.items():
        if r:
            print(f"   {tgt}: {r['gate']} t={r['t']}  TP={r['tp_retained']} genuine={r['genuine_overmerge']} "
                  f"g_blk_prfam={r['genuine_rate_block']} g_blk_edgehost={r['genuine_rate_block_edgehost']} "
                  f"recov97={r['tp_noep_retained']}/{n_tp_noep} recovNC={r['noncoding_tp_retained']}/{n_nc_tp}")
        else:
            print(f"   {tgt}: (none reach it)")
    print(" [transparency] looser EDGE-genuine-rate budget (NOT the recommendation):")
    for tgt, r in combined_edge_knees.items():
        if r:
            print(f"   {tgt}: combined t={r['t']}  TP={r['tp_retained']} genuine={r['genuine_overmerge']} "
                  f"g_edge={r['genuine_rate_edge']} g_blk_prfam={r['genuine_rate_block']} "
                  f"g_blk_edgehost={r['genuine_rate_block_edgehost']} recov97={r['tp_noep_retained']}/{n_tp_noep}")
    if recommended:
        print(f"\nRECOMMENDED (task-literal block<=3% knee): gate={recommended['gate']}  "
              f"t_high={recommended['t']}   [E_r edge = core_recip>=0.13 AND (in_ep OR core_recip>=0.70)]")
        print(f"  KEEPS: TP {recommended['tp_retained']}/450  (protein-missed recovered "
              f"{recommended['tp_noep_retained']}/{n_tp_noep}; of which NON-CODING "
              f"{recommended['noncoding_tp_retained']}/{n_nc_tp})  truthbar(protein-homologous divergent "
              f"paralog) {recommended['truthbar_retained']}/{n_tb} = 100% (TAUTOLOGICAL, truthbar subset in_ep)")
        print(f"  CUTS:  genuine over-merge to {recommended['genuine_overmerge']} edges "
              f"(from {n_ge}, {100*(1-recommended['genuine_overmerge']/n_ge):.1f}% cut); edge genuine rate "
              f"{recommended['genuine_rate_edge']} (from 0.3871); block over-merge PRFAM-mixing "
              f"{recommended['genuine_rate_block']} (from {rows[0]['genuine_rate_block']}) AND edge-host "
              f"{recommended['genuine_rate_block_edgehost']} (from {rows[0]['genuine_rate_block_edgehost']}) "
              f"-- {recommended['genuine_prfam_invisible']}/{recommended['genuine_overmerge']} residual "
              f"over-merges are LOC*/no-PRFAM bridges invisible to the prfam-mixing metric")
        print(f"  reachable recall id 0.90-0.95={recommended['recall_id_0.90-0.95']} "
              f"0.95-0.99={recommended['recall_id_0.95-0.99']} >=0.99={recommended['recall_id_>=0.99']}")
        if disproportionate:
            print(f"  RESIDUAL COST: loses {450-recommended['tp_retained']} TP; NON-CODING TP drop rate "
                  f"{disproportionate['noncoding_drop_rate']} vs coding {disproportionate['coding_drop_rate']} "
                  f"({disproportionate['fold']}x, Fisher OR={disproportionate['fisher_odds_ratio']} "
                  f"p={disproportionate['fisher_p']:.1e}); non-coding-dominant blocks "
                  f"{recommended['noncoding_dominant_blocks']} of which only "
                  f"{recommended['noncoding_dominant_blocks_hosting_tp']} host a REAL TP edge (rest are "
                  f"genuine-only over-merges -- why the old block-dominant 'noncoding_families' count was inverted).")
    print(f"\n[determinism] recommended point reproducible within-run={det}  (PYTHONHASHSEED=0 pinned)")
    print("wrote bench/family_er_oracle_sweep.tsv + bench/family_er_oracle_sweep.json")


if __name__ == "__main__":
    main()
