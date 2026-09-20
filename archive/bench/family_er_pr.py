#!/usr/bin/env python
"""Re-measure precision/recall for the REFRAMED RNA family definition E_r (transcript-homology
component, gamma-quasi-clique refined) and honestly stratify its residual false positives.

E_r pipeline recap
------------------
  * de-novo loci (DN_<contig>_<start>_<nreads>) are grouped by POA transcript homology
    (core_recip >= 0.13) into single-linkage components  -> RAW E_r  (bench/denovo_families.tsv)
  * each raw component is REFINED by the gamma-quasi-clique operator R
    (genome_family_def.refine_families, gamma=0.20, seed=0, >=2 distinct loci) -> REFINED E_r
  * each refined block is projected to GENE space by coordinate max-overlap (gene_of), then
    clique-expanded to unordered gene pairs -> the in_er edge set compared to the gene-keyed truths.

Two truths
----------
  cDNA-90%  : bench/family_def_dna_pr_edges.tsv, col in_dna_loose (all-vs-all cDNA id>=90% cov>=30%).
              STRICT -> under-counts divergent paralogs.  Rows = union(in_rna, in_dna_loose).
  protein Ep: reciprocal whole-protein homology, NO identity floor.
              Tier-1 = direct reciprocal mmseqs edge (aln.m8, evalue<=1e-5 & qcov,tcov>=0.5).
              Tier-2 = same PRFAM in protein_families_refined.tsv, EXCLUDING mega-families (>=500 genes,
                       GPCR/olfactory transitive-closure blobs) unless also backed by a Tier-1 edge.

The honest FP split (deliverable)
---------------------------------
  For every E_r FP edge (in_er & !in_dna_loose):
    protein_homologous := direct_ep(a,b)  OR  (cofam(a,b) & PRFAM not mega)
    id_ok    := id  is None  or id >= 0.50           (guards a repeat/Alu bridge; only if the pair is a truth row)
    cov_ok   := mincov is None or mincov >= 0.60     (guards a shared single-domain hit)
    REAL_DIVERGENT_PARALOG := protein_homologous & id_ok & cov_ok     -> cDNA-90% truth is incomplete, NOT an error
    GENUINE_OVER_MERGE     := not REAL_DIVERGENT_PARALOG              -> a real false positive (headline this)

Outputs: bench/family_er_pr.tsv (per gene-pair), bench/family_er_pr.json (summary).
Run: /home/juanfra/miniforge3/bin/python bench/family_er_pr.py
"""
import json
import os
import sys
from bisect import bisect_right
from collections import defaultdict, Counter
from itertools import combinations

# Pin string-hash iteration order so the louvain splitter (fed set-ordered nodes) yields a byte-identical
# witness partition across process launches. Without this, seed=0 is reproducible WITHIN a process but the
# blob-derived block SET wobbles ~1% blocks / ~5% pairs across launches (the gamma-quasi-clique certificate
# PROPERTY is seed/hash-invariant; the specific witness partition is not -- see genome_family_def.SEED).
if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

BENCH = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BENCH)
import genome_family_def as G   # refine_families, distinct_loci, GAMMA, SEED

SCRATCH = "/home/juanfra/winloci_scratch"
META = os.path.join(SCRATCH, "denovo_transcripts.meta.tsv")
ANNOT = os.path.join(SCRATCH, "annot_intervals.tsv")
ALN = os.path.join(SCRATCH, "aln.m8")
DENOVO_FAMS = os.path.join(BENCH, "denovo_families.tsv")
DENOVO_EDGES = os.path.join(BENCH, "denovo_family_edges.tsv")
TRUTH = os.path.join(BENCH, "family_def_dna_pr_edges.tsv")
PRFAM = os.path.join(BENCH, "protein_families_refined.tsv")

GAMMA = G.GAMMA          # 0.20
ALPHA_P = 1e-5           # protein_family_def.json
C_P = 0.50
MEGA_MIN = 500           # a PRFAM with >=500 genes is a transitive-closure mega-family (co-membership alone is weak)
ID_MIN = 0.50            # repeat/Alu-bridge guard (cDNA nt id)
COV_MIN = 0.60           # shared-single-domain guard (cDNA coverage)
OV_FLOOR = 0.20          # min-overlap floor for the DN-locus -> gene projection: keep a map iff the overlap
                         # covers >= OV_FLOOR of the DN span OR >= OV_FLOOR of the gene span. Drops incidental
                         # <20% clips (e.g. a DN locus grazing ALB by 0.18 of the gene) that otherwise ATTRIBUTE
                         # an unrelated gene into a block and INFLATE the genuine-over-merge (false-positive) count.
                         # The floor can only *remove* spurious attributions -> it tightens (improves) E_r precision;
                         # the pre-floor count is reported as a conservative upper bound (fp_split.genuine_over_merge_prefloor).


# ------------------------------------------------------------------ small helpers
def pct_or_dash(num, den):
    return f"{100*num/den:.1f}%" if den else "-"


# ------------------------------------------------------------------ loaders
def load_meta():
    """DN id -> (chrom, start, end)."""
    m = {}
    with open(META) as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            m[f[0]] = (f[1], int(f[2]), int(f[3]))
    return m


def load_annot():
    """contig -> (starts[], intervals[(start,end,biotype,gene)]) sorted by start; + max interval len per contig."""
    by_c = defaultdict(list)
    with open(ANNOT) as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            by_c[f[0]].append((int(f[1]), int(f[2]), f[3], f[4]))
    idx = {}
    for c, ivals in by_c.items():
        ivals.sort()
        starts = [x[0] for x in ivals]
        maxlen = max((e - s) for s, e, _, _ in ivals)
        idx[c] = (starts, ivals, maxlen)
    return idx


def gene_of_factory(annot):
    """max-overlap gene/biotype for a DN locus interval, plus the overlap fraction
    ovf = max(overlap/DN_span, overlap/gene_span). No floor is applied here (kept pure);
    the OV_FLOOR floor is applied by the caller so both the floored and unfloored maps
    are available for the sensitivity report. Returns (gene, biotype, ovf); (None,None,0.0)
    when the DN locus overlaps no annotation interval."""
    def gene_of(chrom, qs, qe):
        rec = annot.get(chrom)
        if rec is None:
            return None, None, 0.0
        starts, ivals, maxlen = rec
        hi = bisect_right(starts, qe)          # intervals with start <= qe
        best_ov = 0
        best = (None, None)
        best_span = 0
        i = hi - 1
        lo_bound = qs - maxlen
        while i >= 0:
            s, e, bt, g = ivals[i]
            if s < lo_bound:
                break
            if e > qs:                          # overlaps [qs,qe)
                ov = min(e, qe) - max(s, qs)
                if ov > best_ov:
                    best_ov = ov
                    best = (g, bt)
                    best_span = e - s
            i -= 1
        if best[0] is None:
            return None, None, 0.0
        dn_span = qe - qs
        frac_dn = (best_ov / dn_span) if dn_span > 0 else 0.0
        frac_gene = (best_ov / best_span) if best_span > 0 else 0.0
        return best[0], best[1], max(frac_dn, frac_gene)
    return gene_of


def load_raw_families():
    fams = []
    with open(DENOVO_FAMS) as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            fams.append(f[2].split(","))
    return fams


def load_edges():
    edges = []
    with open(DENOVO_EDGES) as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            edges.append((f[0], f[1]))
    return edges


def load_truth():
    """returns:
       dna_loose : set(frozenset({a,b}))            complete (rows>=1 => in_dna_loose complete)
       ec_rna    : set(frozenset({a,b}))            OLD read-conflict edges (in_rna=1)
       rowinfo   : frozenset -> dict(id,cov_a,cov_b,expr_a,expr_b,cross_map,in_dna_loose,in_rna)
    """
    dna_loose = set()
    ec_rna = set()
    rowinfo = {}
    with open(TRUTH) as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            a, b = f[0], f[1]
            key = frozenset((a, b))
            in_rna = int(f[2]); in_dl = int(f[3])
            idv = float(f[13]) if f[13] not in ("", "NA", "None") else None
            cova = float(f[14]) if f[14] not in ("", "NA", "None") else None
            covb = float(f[15]) if f[15] not in ("", "NA", "None") else None
            rowinfo[key] = dict(id=idv, cov_a=cova, cov_b=covb,
                                expr_a=int(f[6]), expr_b=int(f[7]), cross_map=int(f[8]),
                                in_dna_loose=in_dl, in_rna=in_rna)
            if in_dl == 1:
                dna_loose.add(key)
            if in_rna == 1:
                ec_rna.add(key)
    return dna_loose, ec_rna, rowinfo


def load_ep_direct():
    """reciprocal whole-protein edges from aln.m8, keyed by gene symbol; value = mean pident of the two dirs."""
    fwd = {}
    with open(ALN) as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            q, t = f[0], f[1]
            if q == t:
                continue
            try:
                ev, pid, qc, tc = float(f[2]), float(f[3]), float(f[4]), float(f[5])
            except ValueError:
                continue
            if ev > ALPHA_P or qc < C_P or tc < C_P:
                continue
            fwd[(q, t)] = pid
    edges = {}
    for (q, t), pid in fwd.items():
        if (t, q) in fwd:
            key = frozenset((q, t))
            edges[key] = (pid + fwd[(t, q)]) / 2.0
    return edges


def load_prfam():
    """gene -> PRFAM id; and set of MEGA PRFAM ids (>=MEGA_MIN genes)."""
    g2f = {}
    sizes = {}
    with open(PRFAM) as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            fid = f[0]
            genes = f[11].split(",")
            sizes[fid] = len(genes)
            for g in genes:
                g2f[g] = fid
    mega = {fid for fid, n in sizes.items() if n >= MEGA_MIN}
    return g2f, mega


# ------------------------------------------------------------------ E_r gene-pair projection
def blocks_to_gene_pairs(blocks, gene_of_dn):
    """clique-expand each DN-locus block into unordered distinct gene pairs; also return represented gene set."""
    pairs = set()
    represented = set()
    for block in blocks:
        gs = set()
        for dn in block:
            g = gene_of_dn.get(dn)
            if g:
                gs.add(g)
        represented |= gs
        for a, b in combinations(sorted(gs), 2):
            pairs.add(frozenset((a, b)))
    return pairs, represented


# ------------------------------------------------------------------ FP classifier
def protein_homologous(key, ep_direct, g2f, mega):
    a, b = tuple(key)
    if key in ep_direct:
        return True, "tier1"
    fa, fb = g2f.get(a), g2f.get(b)
    if fa is not None and fa == fb and fa not in mega:
        return True, "tier2_nonmega"
    return False, "none"


def classify_fp(key, ep_direct, g2f, mega, rowinfo):
    ph, tier = protein_homologous(key, ep_direct, g2f, mega)
    ri = rowinfo.get(key)
    idv = ri["id"] if ri else None
    if ri and ri["cov_a"] is not None and ri["cov_b"] is not None:
        mincov = min(ri["cov_a"], ri["cov_b"])
    else:
        mincov = None
    id_ok = (idv is None) or (idv >= ID_MIN)
    cov_ok = (mincov is None) or (mincov >= COV_MIN)
    real_paralog = ph and id_ok and cov_ok
    cls = "truthbar-FP" if real_paralog else "genuine-FP"
    return cls, dict(ph=ph, tier=tier, id=idv, mincov=mincov)


# ------------------------------------------------------------------ main
def build_genes_dict(raw_fams, meta, gene_of):
    """DN id -> {chrom,start,end,biotype,name} for every member of a raw family (needed by refine_families).
       Fallback (missing from meta): parse chrom/start from the id, zero-length locus.

       Returns (genes, gene_of_dn, gene_of_dn_raw, miss_meta, ovstats):
         gene_of_dn      : FLOORED map (ovf >= OV_FLOOR)  -- the headline projection.
         gene_of_dn_raw  : UNFLOORED max-overlap map      -- for the conservative pre-floor sensitivity line.
         ovstats         : overlap-fraction contamination counts (loci mapping on a <20% / <50% clip).
       Only the floored map feeds genes[dn].name/biotype so a spurious <20% clip does not seed a block name
       or inflate genuine-over-merge; block MEMBERSHIP is independent of the projection (distinct_loci + adjacency
       use meta coords / edges only), so the floor is precision-only and does not perturb the refined block set."""
    genes = {}
    gene_of_dn = {}
    gene_of_dn_raw = {}
    miss_meta = 0
    lt20 = lt50 = mapped_raw = 0
    for fam in raw_fams:
        for dn in fam:
            if dn in genes:
                continue
            mm = meta.get(dn)
            if mm is None:
                miss_meta += 1
                # DN_<contig...>_<start>_<nreads> ; contig may contain underscores
                parts = dn.split("_")
                start = int(parts[-2]); chrom = "_".join(parts[1:-2])
                end = start
            else:
                chrom, start, end = mm
            g, bt, ovf = gene_of(chrom, start, end)
            keep = bool(g) and ovf >= OV_FLOOR
            genes[dn] = dict(chrom=chrom, start=start, end=end,
                             biotype=(bt if keep else "NA"), name=(g if keep else dn), ov=ovf)
            if g:
                gene_of_dn_raw[dn] = g
                mapped_raw += 1
                if ovf < 0.20:
                    lt20 += 1
                if ovf < 0.50:
                    lt50 += 1
            if keep:
                gene_of_dn[dn] = g
    ovstats = dict(mapped_raw=mapped_raw, mapped_floored=len(gene_of_dn),
                   dropped_by_floor=mapped_raw - len(gene_of_dn),
                   lt20_overlap=lt20, lt50_overlap=lt50)
    return genes, gene_of_dn, gene_of_dn_raw, miss_meta, ovstats


def precision_block(pairs, dna_loose, ep_direct, g2f, mega, rowinfo):
    """precision of a gene-pair set vs cDNA truth and vs Ep (protein-homologous), + FP split."""
    n = len(pairs)
    tp_cdna = sum(1 for k in pairs if k in dna_loose)
    fp = [k for k in pairs if k not in dna_loose]
    ph_count = sum(1 for k in pairs
                   if (k in dna_loose) or protein_homologous(k, ep_direct, g2f, mega)[0])
    genuine = truthbar = 0
    for k in fp:
        cls, _ = classify_fp(k, ep_direct, g2f, mega, rowinfo)
        if cls == "genuine-FP":
            genuine += 1
        else:
            truthbar += 1
    return dict(n=n, tp_cdna=tp_cdna, fp=len(fp),
                prec_cdna=tp_cdna / n if n else 0.0,
                prec_ep=ph_count / n if n else 0.0,
                genuine=genuine, truthbar=truthbar,
                genuine_prec=1 - genuine / n if n else 0.0)


def main():
    print("[load] meta / annot / truth / protein-Ep ...", flush=True)
    meta = load_meta()
    annot = load_annot()
    gene_of = gene_of_factory(annot)
    dna_loose, ec_rna, rowinfo = load_truth()
    ep_direct = load_ep_direct()
    g2f, mega = load_prfam()
    print(f"       dna_loose={len(dna_loose)}  ec_rna={len(ec_rna)}  ep_direct_edges={len(ep_direct)}  "
          f"mega_prfam={sorted(mega)}", flush=True)

    raw_fams = load_raw_families()
    edges = load_edges()
    genes, gene_of_dn, gene_of_dn_raw, miss_meta, ovstats = build_genes_dict(raw_fams, meta, gene_of)
    n_members = sum(len(f) for f in raw_fams)
    print(f"[E_r]  raw_families={len(raw_fams)}  members={n_members}  distinct_dn={len(genes)}  "
          f"mapped_to_gene(floored>={OV_FLOOR})={len(gene_of_dn)}  (unfloored={ovstats['mapped_raw']}, "
          f"dropped_by_floor={ovstats['dropped_by_floor']})  missing_from_meta={miss_meta}", flush=True)
    print(f"[floor] pre-floor map contamination: {ovstats['lt50_overlap']} loci map on <50% overlap, "
          f"{ovstats['lt20_overlap']} on <20% -- these seed spurious gene attributions the floor removes", flush=True)

    # RAW E_r gene pairs (clique-expand raw single-linkage components)
    raw_pairs, raw_repr = blocks_to_gene_pairs(raw_fams, gene_of_dn)

    # REFINED E_r (gamma-quasi-clique) -- pin seed=0.  Block MEMBERSHIP is independent of the gene projection,
    # so the floored and unfloored maps share the SAME blocks0 (only the gene pairs projected out of it differ).
    blocks0 = G.refine_families(raw_fams, edges, genes, gamma=GAMMA, seed=0)
    ref_pairs, ref_repr = blocks_to_gene_pairs(blocks0, gene_of_dn)
    # UNFLOORED projection of the same blocks -> conservative pre-floor upper bound on genuine over-merge.
    ref_pairs_raw, _ = blocks_to_gene_pairs(blocks0, gene_of_dn_raw)
    print(f"[refine] seed0 blocks={len(blocks0)}  refined_gene_pairs(floored)={len(ref_pairs)}  "
          f"(pre-floor={len(ref_pairs_raw)})  raw_gene_pairs={len(raw_pairs)}", flush=True)

    # determinism at fixed seed + multi-seed wobble
    blocks0b = G.refine_families(raw_fams, edges, genes, gamma=GAMMA, seed=0)
    det = (sorted(map(tuple, blocks0)) == sorted(map(tuple, blocks0b)))
    seed_counts = {0: len(blocks0)}
    seed_pairs = {0: len(ref_pairs)}
    for s in (1, 2):
        bl = G.refine_families(raw_fams, edges, genes, gamma=GAMMA, seed=s)
        pr, _ = blocks_to_gene_pairs(bl, gene_of_dn)
        seed_counts[s] = len(bl); seed_pairs[s] = len(pr)
    print(f"[determ] seed0 reproducible={det}  block_counts={seed_counts}  pair_counts={seed_pairs}", flush=True)

    # ------- precision (RAW vs REFINED) -------
    prec_raw = precision_block(raw_pairs, dna_loose, ep_direct, g2f, mega, rowinfo)
    prec_ref = precision_block(ref_pairs, dna_loose, ep_direct, g2f, mega, rowinfo)
    # pre-floor (unfloored projection of the SAME blocks) -> conservative upper bound on the headline over-merge
    prec_ref_prefloor = precision_block(ref_pairs_raw, dna_loose, ep_direct, g2f, mega, rowinfo)
    floor_removed_genuine = prec_ref_prefloor["genuine"] - prec_ref["genuine"]
    print(f"[floor] genuine-over-merge  pre-floor={prec_ref_prefloor['genuine']}  ->  floored={prec_ref['genuine']}  "
          f"(floor removed {floor_removed_genuine} spurious over-merge edges, "
          f"{pct_or_dash(floor_removed_genuine, prec_ref_prefloor['genuine'])})", flush=True)

    # ------- FP sub-split (refined) with biotype detail + named examples -------
    def biotype_of(g):
        # any DN locus mapping to g -> its biotype; fallback via annot not needed (use genes dict names)
        for dn, gg in gene_of_dn.items():
            if gg == g:
                return genes[dn]["biotype"]
        return "NA"
    # cache gene->biotype
    g2bt = {}
    for dn, gg in gene_of_dn.items():
        g2bt.setdefault(gg, genes[dn]["biotype"])

    fp_refined = [k for k in ref_pairs if k not in dna_loose]
    genuine_examples = []
    clean_genuine_examples = []
    paralog_examples = []
    genuine_both_pc = genuine_one_pc = genuine_neither_pc = 0
    genuine_mega_cofam = 0            # genuine-FP where both genes share the SAME mega PRFAM (weak blob homology)
    genuine_znf = 0                   # genuine-FP both genes ZNF/ZFP symbols (KRAB-ZNF paralog sanity anchor)
    tier1_para = tier2_para = 0
    para_pident = []
    genuine_rows = []
    def is_znf(g):
        return g.startswith("ZNF") or g.startswith("ZFP") or g.startswith("ZSCAN") or g.startswith("ZKSCAN")
    for k in fp_refined:
        cls, info = classify_fp(k, ep_direct, g2f, mega, rowinfo)
        a, b = sorted(k)
        bta, btb = g2bt.get(a, "NA"), g2bt.get(b, "NA")
        pc = (bta == "protein_coding") + (btb == "protein_coding")
        fa, fb = g2f.get(a), g2f.get(b)
        mega_cofam = (fa is not None and fa == fb and fa in mega)
        genuine_rows.append((a, b, cls, info, bta, btb, pc))
        if cls == "genuine-FP":
            if pc == 2:
                genuine_both_pc += 1
            elif pc == 1:
                genuine_one_pc += 1
            else:
                genuine_neither_pc += 1
            if mega_cofam:
                genuine_mega_cofam += 1
            if is_znf(a) and is_znf(b):
                genuine_znf += 1
            if pc == 2 and len(genuine_examples) < 40:
                genuine_examples.append((a, b, bta, btb, info["id"], info["mincov"]))
            if pc == 2 and not mega_cofam and not (is_znf(a) and is_znf(b)) and len(clean_genuine_examples) < 20:
                clean_genuine_examples.append((a, b, bta, btb))
        else:
            if info["tier"] == "tier1":
                tier1_para += 1
                if info["id"] is not None:
                    para_pident.append(k)
            else:
                tier2_para += 1
            if len(paralog_examples) < 60:
                paralog_examples.append((a, b, info["tier"], ep_direct.get(k), info["id"]))

    # ------- block-level over-merge (Ep-impurity of refined blocks) -------
    n_blocks = len(blocks0)
    impure_blocks = 0
    blocks_ge2_gene = 0
    for block in blocks0:
        fams_here = set()
        for dn in block:
            g = gene_of_dn.get(dn)
            if g and g in g2f:
                fams_here.add(g2f[g])
        # count blocks that mix >=2 distinct NON-mega protein families
        nonmega = {f for f in fams_here if f not in mega}
        gs = {gene_of_dn[dn] for dn in block if dn in gene_of_dn}
        if len(gs) >= 2:
            blocks_ge2_gene += 1
        if len(nonmega) >= 2:
            impure_blocks += 1

    # ------- recall (identity-stratified) on expressed[/cross_map] subsets; E_r vs E_c -------
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

    def recall_table(require_crossmap, restrict_reachable):
        strata = ["0.90-0.95", "0.95-0.99", ">=0.99", "<0.90", "na"]
        agg = {s: dict(n=0, er=0, ec=0) for s in strata}
        tot = dict(n=0, er=0, ec=0)
        for key in dna_loose:
            ri = rowinfo.get(key)
            if ri is None:
                continue
            if not (ri["expr_a"] == 1 and ri["expr_b"] == 1):
                continue
            if require_crossmap and ri["cross_map"] != 1:
                continue
            a, b = tuple(key)
            if restrict_reachable and not (a in ref_repr and b in ref_repr):
                continue
            s = strat(ri["id"])
            in_er = 1 if key in ref_pairs else 0
            in_ec = 1 if ri["in_rna"] == 1 else 0
            agg[s]["n"] += 1; agg[s]["er"] += in_er; agg[s]["ec"] += in_ec
            tot["n"] += 1; tot["er"] += in_er; tot["ec"] += in_ec
        return agg, tot

    rec_shared_agg, rec_shared = recall_table(True, False)     # (a) expressed & cross_map
    rec_expr_agg, rec_expr = recall_table(False, False)        # (b) expressed only
    rec_reach_agg, rec_reach = recall_table(True, True)        # (c) expressed & cross_map & reachable
    rec_reach_all_agg, rec_reach_all = recall_table(False, True)

    # what the (a)->(c) reachability conditioning actually does (guardrail for the +gain claim):
    # it drops the pairs where E_r's de-novo footprint does not reach both genes. By construction in_er
    # requires both genes co-blocked => both represented => reachable, so E_r loses 0 TP to this cut, while
    # E_c can lose the pairs it correctly grouped that lie OUTSIDE E_r's footprint. Hence (c) isolates the
    # family-def CRITERION but its denominator is convenient; it is NEVER absolute recall superiority.
    cond_dropped_pairs = rec_shared["n"] - rec_reach["n"]
    cond_ec_lost = rec_shared["ec"] - rec_reach["ec"]
    cond_er_lost = rec_shared["er"] - rec_reach["er"]
    print(f"[cond] (a)->(c) reachability cut drops {cond_dropped_pairs} dna_loose pairs: "
          f"E_c loses {cond_ec_lost} correctly-grouped TP, E_r loses {cond_er_lost} "
          f"(0 by construction) -- (c) measures the criterion, not absolute recall", flush=True)

    def rate(d):
        return dict(n=d["n"],
                    er=round(d["er"] / d["n"], 4) if d["n"] else None,
                    ec=round(d["ec"] / d["n"], 4) if d["n"] else None,
                    er_ct=d["er"], ec_ct=d["ec"])

    # ------- represented-gene coverage of the cDNA truth -------
    dna_loose_genes = set()
    for k in dna_loose:
        dna_loose_genes |= set(k)
    covered = len(dna_loose_genes & ref_repr)
    # raw representation = any DN locus maps to the gene (superset of refined-block genes)
    raw_repr_genes = set(gene_of_dn.values())
    covered_raw = len(dna_loose_genes & raw_repr_genes)
    clean_genuine = genuine_both_pc - genuine_mega_cofam   # both-pc genuine that are NOT mega-blob paralogs
    print(f"[diag] raw_repr_genes={len(raw_repr_genes)} (covers {covered_raw} dna_loose genes)  "
          f"ref_repr_genes={len(ref_repr)} (covers {covered})", flush=True)
    print(f"[diag] genuine both_pc={genuine_both_pc}  of which mega-cofam(ZNF/GPCR blob)={genuine_mega_cofam}  "
          f"ZNF-ZNF={genuine_znf}  => clean both-pc genuine(non-mega)={clean_genuine}", flush=True)

    # ============================== write per-edge TSV ==============================
    with open(os.path.join(BENCH, "family_er_pr.tsv"), "w") as out:
        out.write("gene_a\tgene_b\tin_er_raw\tin_er_refined\tin_dna_loose\tin_ep\tep_tier\t"
                  "id\tmin_cov\tclass\n")
        all_keys = raw_pairs | ref_pairs | dna_loose
        for key in sorted(all_keys, key=lambda k: tuple(sorted(k))):
            a, b = sorted(key)
            in_raw = 1 if key in raw_pairs else 0
            in_ref = 1 if key in ref_pairs else 0
            in_dl = 1 if key in dna_loose else 0
            ph, tier = protein_homologous(key, ep_direct, g2f, mega)
            ri = rowinfo.get(key)
            idv = ri["id"] if ri else None
            mincov = (min(ri["cov_a"], ri["cov_b"]) if ri and ri["cov_a"] is not None
                      and ri["cov_b"] is not None else None)
            if in_ref:
                if in_dl:
                    cls = "TP"                                   # in refined E_r AND cDNA-loose truth
                else:
                    cls, _ = classify_fp(key, ep_direct, g2f, mega, rowinfo)  # genuine-FP / truthbar-FP
            else:
                cls = "FN" if in_dl else "not-er"               # cDNA-loose truth missed by refined E_r / neither
            out.write(f"{a}\t{b}\t{in_raw}\t{in_ref}\t{in_dl}\t{1 if ph else 0}\t{tier}\t"
                      f"{'' if idv is None else round(idv,4)}\t"
                      f"{'' if mincov is None else round(mincov,4)}\t{cls}\n")

    # ============================== summary JSON ==============================
    summary = dict(
        params=dict(gamma=GAMMA, seed=0, alpha_p=ALPHA_P, c_p=C_P, mega_min=MEGA_MIN,
                    id_min=ID_MIN, cov_min=COV_MIN, ov_floor=OV_FLOOR, mega_prfam=sorted(mega)),
        er=dict(raw_families=len(raw_fams), members=n_members, distinct_dn=len(genes),
                mapped_to_gene=len(gene_of_dn), mapped_to_gene_prefloor=ovstats["mapped_raw"],
                dropped_by_floor=ovstats["dropped_by_floor"], missing_from_meta=miss_meta,
                refined_blocks=n_blocks, blocks_ge2_gene=blocks_ge2_gene,
                raw_gene_pairs=len(raw_pairs), refined_gene_pairs=len(ref_pairs),
                refined_gene_pairs_prefloor=len(ref_pairs_raw),
                represented_genes=len(ref_repr)),
        overlap_floor=dict(ov_floor=OV_FLOOR,
                           mapped_prefloor=ovstats["mapped_raw"], mapped_floored=len(gene_of_dn),
                           dropped_by_floor=ovstats["dropped_by_floor"],
                           prefloor_maps_lt50pct=ovstats["lt50_overlap"],
                           prefloor_maps_lt20pct=ovstats["lt20_overlap"],
                           genuine_over_merge_prefloor=prec_ref_prefloor["genuine"],
                           genuine_over_merge_floored=prec_ref["genuine"],
                           floor_removed_genuine=floor_removed_genuine),
        determinism=dict(seed0_reproducible=det, block_counts=seed_counts, pair_counts=seed_pairs),
        precision_raw=prec_raw,
        precision_refined=prec_ref,
        fp_split_refined=dict(
            total_fp=len(fp_refined),
            genuine_over_merge=prec_ref["genuine"],
            genuine_rate=round(prec_ref["genuine"] / max(1, len(ref_pairs)), 4),
            real_divergent_paralog=prec_ref["truthbar"],
            paralog_rate=round(prec_ref["truthbar"] / max(1, len(fp_refined)), 4),
            genuine_both_pc=genuine_both_pc,
            genuine_one_pc=genuine_one_pc,
            genuine_neither_pc=genuine_neither_pc,
            genuine_mega_cofam=genuine_mega_cofam,
            genuine_znf_znf=genuine_znf,
            clean_both_pc_genuine_nonmega=clean_genuine,
            paralog_tier1_direct=tier1_para,
            paralog_tier2_nonmega=tier2_para),
        block_over_merge=dict(refined_blocks=n_blocks,
                              ep_impure_blocks=impure_blocks,
                              impure_rate=round(impure_blocks / max(1, n_blocks), 4)),
        recall=dict(
            shared_expr_crossmap=dict(overall=rate(rec_shared),
                                      strata={s: rate(rec_shared_agg[s]) for s in rec_shared_agg}),
            expressed_only=dict(overall=rate(rec_expr),
                                strata={s: rate(rec_expr_agg[s]) for s in rec_expr_agg}),
            reachable_expr_crossmap=dict(overall=rate(rec_reach),
                                         strata={s: rate(rec_reach_agg[s]) for s in rec_reach_agg}),
            reachable_expr_only=dict(overall=rate(rec_reach_all)),
            reachability_conditioning=dict(
                note=("(a)->(c) cut drops pairs outside E_r's de-novo footprint; E_r loses 0 TP by construction "
                      "(in_er => both genes represented => reachable), E_c loses the pairs it correctly grouped "
                      "outside that footprint. (c) isolates the CRITERION on a convenient denominator; it is not "
                      "absolute recall. Absolute recall (subsets a/b) has E_r BELOW E_c."),
                dropped_pairs=cond_dropped_pairs, ec_lost_tp=cond_ec_lost, er_lost_tp=cond_er_lost,
                representation_ceiling_refined=round(covered / max(1, len(dna_loose_genes)), 4),
                representation_ceiling_raw=round(covered_raw / max(1, len(dna_loose_genes)), 4))),
        truth_coverage=dict(dna_loose_genes=len(dna_loose_genes),
                            represented_genes_refined=len(ref_repr),
                            represented_genes_raw=len(raw_repr_genes),
                            dna_loose_genes_represented_refined=covered,
                            dna_loose_genes_represented_raw=covered_raw),
    )
    with open(os.path.join(BENCH, "family_er_pr.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

    # ============================== console report ==============================
    def pct(x):
        return f"{100*x:.1f}%"
    print("\n================ E_r PRECISION ================")
    print(f"RAW E_r     : gene_pairs={prec_raw['n']:>7}  TP_cdna={prec_raw['tp_cdna']:>6}  "
          f"prec_cdna={prec_raw['prec_cdna']:.3f}  prec_Ep={prec_raw['prec_ep']:.3f}  "
          f"genuine_over_merge={prec_raw['genuine']}  genuine_prec={prec_raw['genuine_prec']:.3f}")
    print(f"REFINED E_r : gene_pairs={prec_ref['n']:>7}  TP_cdna={prec_ref['tp_cdna']:>6}  "
          f"prec_cdna={prec_ref['prec_cdna']:.3f}  prec_Ep={prec_ref['prec_ep']:.3f}  "
          f"genuine_over_merge={prec_ref['genuine']}  genuine_prec={prec_ref['genuine_prec']:.3f}")
    print("\n================ FP SPLIT (refined) ================")
    print(f"total FP (in_er & !in_dna_loose) = {len(fp_refined)}")
    print(f"  REAL divergent paralog (truth-missed) = {prec_ref['truthbar']}  "
          f"({pct(prec_ref['truthbar']/max(1,len(fp_refined)))} of FP)  "
          f"[tier1_direct={tier1_para}  tier2_nonmega={tier2_para}]")
    print(f"  GENUINE over-merge (HEADLINE)         = {prec_ref['genuine']}  "
          f"({pct(prec_ref['genuine']/max(1,len(fp_refined)))} of FP)  "
          f"[both_pc={genuine_both_pc}  one_pc={genuine_one_pc}  neither_pc={genuine_neither_pc}]")
    print(f"     of which MEGA-cofam (ZNF/GPCR blob, weak-homology, arguably real paralog)={genuine_mega_cofam}"
          f"  [ZNF-ZNF={genuine_znf}]  => CLEAN both-pc genuine (non-blob)={clean_genuine}")
    print(f"  => genuine-over-merge RATE of all refined E_r edges = "
          f"{prec_ref['genuine']}/{len(ref_pairs)} = {pct(prec_ref['genuine']/max(1,len(ref_pairs)))}"
          f"   (genuine precision {prec_ref['genuine_prec']:.3f})")
    print(f"  min-overlap floor (>= {OV_FLOOR} of DN or gene span): genuine over-merge "
          f"{prec_ref_prefloor['genuine']} (pre-floor, conservative upper bound) -> {prec_ref['genuine']} (floored); "
          f"floor removed {floor_removed_genuine} spurious <20% clip attributions")
    print(f"  *** HEADLINE family-level over-merge = {impure_blocks}/{n_blocks} refined blocks Ep-impure "
          f"({pct(impure_blocks/max(1,n_blocks))}) -- lead with THIS, not the clique-inflated {pct(prec_ref['genuine']/max(1,len(ref_pairs)))} edge rate")
    print("\n  named GENUINE over-merge examples (both protein_coding, not Ep-cofamily):")
    for a, b, bta, btb, idv, mc in genuine_examples[:12]:
        print(f"    {a} -- {b}   [{bta}/{btb}]  id={idv}  mincov={mc}")
    print("  named CLEAN genuine over-merges (both-pc, NOT mega-blob, NOT ZNF-ZNF):")
    for a, b, bta, btb in clean_genuine_examples[:12]:
        print(f"    {a} -- {b}   [{bta}/{btb}]")
    print("\n  named REAL divergent-paralog 'FP' examples (Ep-homologous, below cDNA-90%):")
    for a, b, tier, pid, idv in paralog_examples[:15]:
        print(f"    {a} -- {b}   Ep={tier}  prot_pident={pid}  cdna_id={idv}")
    print("\n================ RECALL (E_r vs E_c=in_rna) ================")
    for label, agg, tot in [("(a) expressed & cross_map", rec_shared_agg, rec_shared),
                            ("(b) expressed only", rec_expr_agg, rec_expr),
                            ("(c) reachable & expressed & cross_map", rec_reach_agg, rec_reach)]:
        r = rate(tot)
        print(f"{label}: n={r['n']}  E_r={r['er']}  E_c={r['ec']}  gain={None if r['er'] is None else round(r['er']-r['ec'],4)}")
        for s in ["0.90-0.95", "0.95-0.99", ">=0.99"]:
            rr = rate(agg[s])
            print(f"      id {s:<10} n={rr['n']:>4}  E_r={rr['er']}  E_c={rr['ec']}")
    r = rate(rec_reach_all)
    print(f"(c') reachable & expressed (any cross_map): n={r['n']}  E_r={r['er']}  E_c={r['ec']}")
    print(f"  ^ CAVEAT: the +{round(rec_reach['er']/max(1,rec_reach['n'])-rec_reach['ec']/max(1,rec_reach['n']),3)} "
          f"gain in (c) is CONDITIONAL on reachability. The (a)->(c) cut drops {cond_dropped_pairs} dna_loose pairs; "
          f"E_c loses {cond_ec_lost} correctly-grouped TP, E_r loses {cond_er_lost} (0 by construction). "
          f"On the ABSOLUTE subsets (a)/(b) E_r recall is BELOW E_c. (c) isolates the CRITERION, not absolute recall.")
    print(f"\ntruth coverage: refined-block genes represent {covered}/{len(dna_loose_genes)} cDNA-loose genes "
          f"({pct(covered/max(1,len(dna_loose_genes)))} <- cite THIS for refined recall); ANY DN locus represents "
          f"{covered_raw}/{len(dna_loose_genes)} ({pct(covered_raw/max(1,len(dna_loose_genes)))}) "
          f"-- the representation ceiling that caps raw recall (a)/(b)")
    print(f"\nTRADE-OFF: E_r buys higher CONDITIONAL low-id recall (id0.90-0.95 reachable "
          f"{rate(rec_reach_agg['0.90-0.95'])['er']} vs E_c {rate(rec_reach_agg['0.90-0.95'])['ec']}) "
          f"at the cost of genuine precision {prec_ref['genuine_prec']:.3f} (vs E_c-refined 0.94). "
          f"E_r does NOT dominate E_c on both axes.")
    print(f"\ndeterminism: seed0 reproducible={det}  seed block counts={seed_counts}  pairs={seed_pairs}")
    print("\nwrote bench/family_er_pr.tsv + bench/family_er_pr.json")


if __name__ == "__main__":
    main()
