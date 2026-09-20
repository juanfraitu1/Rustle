#!/usr/bin/env python3
"""Post-refinement family merge by exon colinearity.

Insert this stage AFTER allele-demote and BEFORE validation/output in family_rna_refine.py.
It merges catalog families that share >= MIN_COLINEAR exons in colinear order between any
pair of loci, subject to the existing antisense/reciprocal-overlap and repeat-hub gates.

The goal is to recover known families that the gamma-quasi-clique refine splits at weak
edges (e.g. MAGEA fam510 + fam508) while keeping domain-sharer over-merges out (they
typically share only one isolated exon).
"""
import os
import sys
from collections import defaultdict, Counter

import pysam

BENCH = os.path.dirname(os.path.abspath(__file__))
SCRATCH = "/home/juanfra/winloci_scratch"
sys.path.insert(0, BENCH)
import recombination_bridge_detector as R
import colinear_multiexon_gate as CM

SKEL_TSV = os.path.join(SCRATCH, "denovo_skeletons.tsv")
META_TSV = os.path.join(SCRATCH, "denovo_transcripts.meta.tsv")
GENOME = os.path.join(SCRATCH, "GGO.fasta")

# Default exon identity threshold.  0.70 is the project-wide ID_THRESH used by
# recombination_bridge_detector / colinear_multiexon_gate.  It is strict enough to keep
# domain-sharers from chaining spurious exon matches, but still recovers MAGEA-style splits.
ID_THRESH = 0.70
MIN_COLINEAR = 2
WINDOW_BP = 5_000_000
# For window-only pairs (no shared gene symbol, no raw homology edge) we demand structural
# evidence stronger than two isolated shared exons: an adaptive adjacent-junction floor.
# The required floor is min(MIN_ADJACENT_JUNCTIONS, colinear_exons - 1), so a 2-exon hit
# needs 1 adjacent preserved splice junction and a >=3-exon hit needs 2.  This blocks
# domain-sharer neighbours (ANKRD18 + ANKRD36C) while keeping real split families whose raw
# edges or gene content already support them (MAGEA, GSTM).
MIN_ADJACENT_JUNCTIONS = 2

# Divergent same-chromosome duplicon merge (e.g. HERC2).  Recent segduplications can diverge
# below the 0.70 exon-identity floor while still retaining a colinear exon backbone.  We use
# a lower identity threshold but demand a longer colinear chain so domain-sharers stay out.
DIVERGENT_ID_THRESH = 0.55
DIVERGENT_MIN_COLINEAR = 4
DIVERGENT_MIN_ADJACENT_JUNCTIONS = 3
DIVERGENT_WINDOW_BP = 10_000_000

# Cross-chromosome domain-bridge split.  Same-chromosome duplication is the common case;
# cross-chromosome homology claims need stronger evidence or a curated annotation link.
CROSSCHROM_ID_THRESH = 0.80
CROSSCHROM_MIN_COLINEAR = 4
CROSSCHROM_MIN_ADJACENT_JUNCTIONS = 3


def load_meta():
    nreads = {}
    strand = {}
    with open(META_TSV) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        I = {c: i for i, c in enumerate(h)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            dn = f[I["id"]]
            nreads[dn] = int(f[I["n_reads"]])
            strand[dn] = f[I["strand"]] if len(f) > I["strand"] else "+"
    return nreads, strand


def _load_fasta_once():
    if not hasattr(_load_fasta_once, "fa"):
        _load_fasta_once.fa = pysam.FastaFile(GENOME)
    return _load_fasta_once.fa


def _load_skeletons_once():
    if not hasattr(_load_skeletons_once, "skel"):
        _load_skeletons_once.skel = R.load_skeletons()
    return _load_skeletons_once.skel


def block_chrom_span(block, genes):
    """Return (chrom, min_start, max_end) for a catalog block; uses dominant chromosome."""
    chroms = []
    starts = []
    ends = []
    for dn in block:
        L = genes.get(dn)
        if L:
            chroms.append(L["chrom"])
            starts.append(L["start"])
            ends.append(L["end"])
    if not chroms:
        return None, 0, 0
    ccount = Counter(chroms)
    dom_chrom = ccount.most_common(1)[0][0]
    idx = [i for i, c in enumerate(chroms) if c == dom_chrom]
    return dom_chrom, min(starts[i] for i in idx), max(ends[i] for i in idx)


def _block_recs(block, genes, skel, strand, fa, cache):
    key = tuple(block)
    if key not in cache:
        mem = [dict(dn=dn, gene=genes.get(dn, {}).get("gene", ""),
                    chrom=genes[dn]["chrom"], start=genes[dn]["start"], end=genes[dn]["end"])
               for dn in block if dn in genes]
        cache[key] = R.family_exons(mem, skel, strand, fa)
    return cache[key]


def _cross_exon_match(rec_a, rec_b):
    """Compute best cross matches only: best[(ia, ei, ib)] for a's exons against b."""
    best = {}
    ea = rec_a["exons"]
    eb = rec_b["exons"]
    for i, qa in enumerate(ea):
        if len(qa) < R.MIN_EXON:
            continue
        bid, bj = 0.0, -1
        for j, tb in enumerate(eb):
            if len(tb) < R.MIN_EXON:
                continue
            lo, hi = (len(qa), len(tb)) if len(qa) <= len(tb) else (len(tb), len(qa))
            if lo / hi < 0.5:
                continue
            if len(qa) <= len(tb):
                idv = R.aln_id(qa, tb)
            else:
                idv = R.aln_id(tb, qa)
            if idv > bid:
                bid, bj = idv, j
        best[(i, 0)] = (bid, bj)   # 0 is a dummy; we only need one direction per pair
    return best


def best_colinear_between_blocks(block_a, block_b, genes, skel, strand, fa,
                                 id_thresh=ID_THRESH, rec_cache=None):
    """Best strict-LIS colinear shared-exon count and adjacent-junction count between any locus in A and any in B.

    Returns (max_colinear_exons, max_adjacent_junctions, n_loci_a_with_skeleton,
             n_loci_b_with_skeleton).  Adjacent junctions = number of consecutive matched
    exon pairs (i,i+1) in A that map to consecutive exons (j,j+1) in B.
    """
    if rec_cache is None:
        rec_cache = {}
    recs_a = _block_recs(block_a, genes, skel, strand, fa, rec_cache)
    recs_b = _block_recs(block_b, genes, skel, strand, fa, rec_cache)
    if not recs_a or not recs_b:
        return 0, 0, len(recs_a), len(recs_b)
    max_col = 0
    max_junc = 0
    for ra in recs_a:
        for rb in recs_b:
            # Build a fake tensor for colinear_count and adjacent-junction counting
            best = {}
            for i, qa in enumerate(ra["exons"]):
                if len(qa) < R.MIN_EXON:
                    continue
                bid, bj = 0.0, -1
                for j, tb in enumerate(rb["exons"]):
                    if len(tb) < R.MIN_EXON:
                        continue
                    lo, hi = (len(qa), len(tb)) if len(qa) <= len(tb) else (len(tb), len(qa))
                    if lo / hi < 0.5:
                        continue
                    if len(qa) <= len(tb):
                        idv = R.aln_id(qa, tb)
                    else:
                        idv = R.aln_id(tb, qa)
                    if idv > bid:
                        bid, bj = idv, j
                best[(0, i, 1)] = (bid, bj)   # (a=0, i, b=1)
            col = CM.colinear_count(best, 0, 1, len(ra["exons"]), strict=True, thresh=id_thresh)
            if col > max_col:
                max_col = col
            # adjacent preserved junctions
            pairs = []
            for i in range(len(ra["exons"])):
                idv, j = best.get((0, i, 1), (0.0, -1))
                if idv >= id_thresh and j >= 0:
                    pairs.append((i, j))
            pairs.sort()
            junc = 0
            for k in range(1, len(pairs)):
                if pairs[k][0] == pairs[k - 1][0] + 1 and pairs[k][1] == pairs[k - 1][1] + 1:
                    junc += 1
            if junc > max_junc:
                max_junc = junc
    return max_col, max_junc, len(recs_a), len(recs_b)


def merge_catalog_colinear(catalog, genes, gene_of_dn, gene_strand=None, pair_repeat_mult=None,
                           raw_edge_pairs=None,
                           min_colinear=MIN_COLINEAR, min_adjacent_junctions=MIN_ADJACENT_JUNCTIONS,
                           id_thresh=ID_THRESH, window_bp=WINDOW_BP,
                           antisense_recip_min=0.50, mega_span_max=500_000,
                           repeat_mult_min=20, verbose=True):
    """Merge catalog blocks based on colinear shared exons.

    A pair of blocks is considered if:
      - they share a dominant chromosome, AND
      - either they share an annotated gene symbol, OR their spans are within window_bp.
    They merge if the best strict-LIS colinear shared-exon count between any two loci is
    >= min_colinear at id_thresh, and the bridge is not blocked by the antisense or repeat gates.

    Window-only pairs (no shared gene symbol and no raw homology edge from the pre-gate graph)
    must additionally show contiguous splice-junction preservation: >= min_adjacent_junctions
    adjacent matched-exon pairs, but never more than the chain itself can provide
    (max(0, colinear_exons - 1)).  This means a 2-exon colinear hit needs one adjacent pair,
    while a >=3-exon hit needs two.  The adaptive floor keeps domain-sharer neighbours
    (e.g. ANKRD18 + ANKRD36C: col=3, junc=1) from merging while still recovering short
    real split-family blocks (e.g. GSTM1/2/4/5: col=2, junc=1).
    """
    fa = _load_fasta_once()
    skel = _load_skeletons_once()
    _, strand = load_meta()
    rec_cache = {}

    n = len(catalog)
    block_meta = []
    for idx, block in enumerate(catalog):
        chrom, s, e = block_chrom_span(block, genes)
        block_genes = {gene_of_dn.get(dn) for dn in block if gene_of_dn.get(dn)}
        block_meta.append((idx, chrom, s, e, block, block_genes))

    def antisense_overlap(ga, gb):
        if not gene_strand:
            return False
        ia = gene_strand.get(ga)
        ib = gene_strand.get(gb)
        if ia is None or ib is None:
            return False
        ca, sa, ea, ta = ia
        cb, sb, eb, tb = ib
        if ca != cb or ta == tb:
            return False
        spa, spb = ea - sa, eb - sb
        if spa <= 0 or spb <= 0:
            return False
        if spa >= mega_span_max or spb >= mega_span_max:
            return False
        ov = min(ea, eb) - max(sa, sb)
        if ov <= 0:
            return False
        return (ov / min(spa, spb)) >= antisense_recip_min

    def repeat_hub(ga, gb):
        if not pair_repeat_mult:
            return False
        return pair_repeat_mult.get(frozenset((ga, gb)), 0) >= repeat_mult_min

    # Index raw edges by the catalog block indices they connect (if provided)
    raw_edges_between = defaultdict(set)
    if raw_edge_pairs:
        dn_to_idx = {}
        for idx, block in enumerate(catalog):
            for dn in block:
                dn_to_idx[dn] = idx
        for a, b in raw_edge_pairs:
            ia = dn_to_idx.get(a)
            ib = dn_to_idx.get(b)
            if ia is not None and ib is not None and ia != ib:
                raw_edges_between[min(ia, ib)].add(max(ia, ib))

    merge_edges = []
    chrom_blocks = defaultdict(list)
    for bm in block_meta:
        if bm[1] is not None:
            chrom_blocks[bm[1]].append(bm)

    for chrom, blist in sorted(chrom_blocks.items()):
        blist = sorted(blist, key=lambda x: x[2])
        m = len(blist)
        for i in range(m):
            idx_a, _, sa, ea, block_a, genes_a = blist[i]
            for j in range(i + 1, m):
                idx_b, _, sb, eb, block_b, genes_b = blist[j]
                within_window = (sb - ea) <= window_bp
                share_gene = bool(genes_a & genes_b)
                has_raw_edge = idx_b in raw_edges_between.get(idx_a, set())
                if not (share_gene or has_raw_edge or within_window):
                    if sb - ea > window_bp:
                        break
                    continue
                col, junc, na, nb = best_colinear_between_blocks(
                    block_a, block_b, genes, skel, strand, fa,
                    id_thresh=id_thresh, rec_cache=rec_cache)
                if col < min_colinear:
                    continue
                weak_support = within_window and not share_gene and not has_raw_edge
                # adaptive junction floor: a short 2-exon colinear block needs one adjacent
                # pair; a longer block needs the floor, capped by what the chain can provide.
                if weak_support and junc < min(min_adjacent_junctions, max(0, col - 1)):
                    continue
                skip = False
                for ga in sorted(genes_a):
                    for gb in sorted(genes_b):
                        if ga == gb:
                            continue
                        if antisense_overlap(ga, gb) or repeat_hub(ga, gb):
                            skip = True
                            break
                    if skip:
                        break
                if skip:
                    continue
                merge_edges.append(dict(
                    a=idx_a, b=idx_b, chrom=chrom,
                    colinear_exons=col, adjacent_junctions=junc, n_loci_a=na, n_loci_b=nb,
                    genes_a=sorted(genes_a)[:20], genes_b=sorted(genes_b)[:20],
                    within_window=within_window, share_gene=share_gene, has_raw_edge=has_raw_edge,
                    weak_support=weak_support))

    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry

    for e in merge_edges:
        union(e["a"], e["b"])

    groups = defaultdict(list)
    for idx, block in enumerate(catalog):
        groups[find(idx)].append(block)

    merged = []
    for root in sorted(groups.keys()):
        combined = []
        for block in groups[root]:
            combined.extend(block)
        merged.append(sorted(set(combined)))

    if verbose:
        print(f"[family_merge_colinear] {n} blocks -> {len(merged)} blocks; edges={len(merge_edges)} "
              f"(id_thresh={id_thresh}, min_colinear={min_colinear}, "
              f"min_adjacent_junctions={min_adjacent_junctions}, window_bp={window_bp})")
        for e in merge_edges:
            if e["share_gene"]:
                reason = "share_gene"
            elif e["has_raw_edge"]:
                reason = "raw_edge"
            elif e["weak_support"]:
                reason = "window_weak"
            else:
                reason = "window"
            need = min(min_adjacent_junctions, max(0, e["colinear_exons"] - 1))
            print(f"   merge fam{e['a']} + fam{e['b']}  chrom={e['chrom']} "
                  f"colinear={e['colinear_exons']} adjacent_junctions={e['adjacent_junctions']} "
                  f"need_junc={need} reason={reason} "
                  f"genes_a={','.join(e['genes_a'][:5])}... "
                  f"genes_b={','.join(e['genes_b'][:5])}...")

    return merged, merge_edges


def _union_find_components(n, edges):
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry

    for i, j in edges:
        union(i, j)
    groups = defaultdict(list)
    for idx in range(n):
        groups[find(idx)].append(idx)
    return list(groups.values())


def _gene_symbol_root(symbol):
    """Strip trailing letter/digit suffixes (e.g. RABL2A->RABL2, ZNF92->ZNF).

    Returns None for NA / empty / LOC symbols so different LOCs are not collapsed.
    """
    if not symbol or symbol == "NA" or symbol.startswith("LOC"):
        return None
    s = symbol.rstrip("0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ")
    return s if s else None


def merge_catalog_divergent_samechrom(catalog, genes, gene_of_dn, gene_strand=None,
                                      pair_repeat_mult=None, raw_edge_pairs=None,
                                      id_thresh=DIVERGENT_ID_THRESH,
                                      min_colinear=DIVERGENT_MIN_COLINEAR,
                                      min_adjacent_junctions=DIVERGENT_MIN_ADJACENT_JUNCTIONS,
                                      window_bp=DIVERGENT_WINDOW_BP,
                                      antisense_recip_min=0.50, mega_span_max=500_000,
                                      repeat_mult_min=20, verbose=True):
    """Second merge pass for divergent same-chromosome duplicons (e.g. HERC2).

    Only considers same-chromosome block pairs that share an annotated gene symbol OR are
    connected by a raw homology edge.  Merge if they share >= ``min_colinear`` exons in colinear
    order at ``id_thresh`` (lower than the standard pass) and the adaptive adjacent-junction
    floor is met.  The high colinear/junction bar (4 exons / 3 adjacent junctions by default)
    keeps short domain-sharers (e.g. the GSTM2 GST-domain hub) from merging.
    """
    fa = _load_fasta_once()
    skel = _load_skeletons_once()
    _, strand = load_meta()
    rec_cache = {}

    n = len(catalog)
    block_meta = []
    for idx, block in enumerate(catalog):
        chrom, s, e = block_chrom_span(block, genes)
        block_genes = {gene_of_dn.get(dn) for dn in block if gene_of_dn.get(dn)}
        block_meta.append((idx, chrom, s, e, block, block_genes))

    def antisense_overlap(ga, gb):
        if not gene_strand:
            return False
        ia = gene_strand.get(ga)
        ib = gene_strand.get(gb)
        if ia is None or ib is None:
            return False
        ca, sa, ea, ta = ia
        cb, sb, eb, tb = ib
        if ca != cb or ta == tb:
            return False
        spa, spb = ea - sa, eb - sb
        if spa <= 0 or spb <= 0:
            return False
        if spa >= mega_span_max or spb >= mega_span_max:
            return False
        ov = min(ea, eb) - max(sa, sb)
        if ov <= 0:
            return False
        return (ov / min(spa, spb)) >= antisense_recip_min

    def repeat_hub(ga, gb):
        if not pair_repeat_mult:
            return False
        return pair_repeat_mult.get(frozenset((ga, gb)), 0) >= repeat_mult_min

    raw_edges_between = defaultdict(set)
    if raw_edge_pairs:
        dn_to_idx = {}
        for idx, block in enumerate(catalog):
            for dn in block:
                dn_to_idx[dn] = idx
        for a, b in raw_edge_pairs:
            ia = dn_to_idx.get(a)
            ib = dn_to_idx.get(b)
            if ia is not None and ib is not None and ia != ib:
                raw_edges_between[min(ia, ib)].add(max(ia, ib))

    merge_edges = []
    chrom_blocks = defaultdict(list)
    for bm in block_meta:
        if bm[1] is not None:
            chrom_blocks[bm[1]].append(bm)

    for chrom, blist in sorted(chrom_blocks.items()):
        blist = sorted(blist, key=lambda x: x[2])
        m = len(blist)
        for i in range(m):
            idx_a, _, sa, ea, block_a, genes_a = blist[i]
            for j in range(i + 1, m):
                idx_b, _, sb, eb, block_b, genes_b = blist[j]
                share_gene = bool(genes_a & genes_b)
                # Require a shared annotated gene symbol and that the smaller block does not
                # introduce an unrelated named gene (e.g. a GSTM2 block must not pull in TRIP13).
                if not share_gene:
                    continue
                nonloc_a = {g for g in genes_a if g and g != "NA" and not g.startswith("LOC")}
                nonloc_b = {g for g in genes_b if g and g != "NA" and not g.startswith("LOC")}
                smaller_nonloc = nonloc_a if len(block_a) <= len(block_b) else nonloc_b
                larger_nonloc = nonloc_b if len(block_a) <= len(block_b) else nonloc_a
                if not (smaller_nonloc <= larger_nonloc):
                    continue
                # Block merging heterogeneous multi-named blocks (e.g. GSTM2 into GSTM1/2/4/5):
                # require the larger block to have a single non-LOC gene identity, so the merge is
                # a fragment-of-one-named-gene situation (e.g. HERC2), not a multi-paralog cluster.
                if len(larger_nonloc) != 1:
                    continue
                col, junc, na, nb = best_colinear_between_blocks(
                    block_a, block_b, genes, skel, strand, fa,
                    id_thresh=id_thresh, rec_cache=rec_cache)
                if col < min_colinear:
                    continue
                if junc < min(min_adjacent_junctions, max(0, col - 1)):
                    continue
                skip = False
                for ga in sorted(genes_a):
                    for gb in sorted(genes_b):
                        if ga == gb:
                            continue
                        if antisense_overlap(ga, gb) or repeat_hub(ga, gb):
                            skip = True
                            break
                    if skip:
                        break
                if skip:
                    continue
                merge_edges.append(dict(
                    a=idx_a, b=idx_b, chrom=chrom,
                    colinear_exons=col, adjacent_junctions=junc, n_loci_a=na, n_loci_b=nb,
                    genes_a=sorted(genes_a)[:20], genes_b=sorted(genes_b)[:20],
                    share_gene=share_gene,
                    divergent=True))

    components = _union_find_components(n, [(e["a"], e["b"]) for e in merge_edges])
    merged = []
    for comp in components:
        combined = []
        for idx in comp:
            combined.extend(block_meta[idx][4])
        merged.append(sorted(set(combined)))

    if verbose:
        print(f"[family_merge_colinear divergent] {n} blocks -> {len(merged)} blocks; "
              f"edges={len(merge_edges)} (id_thresh={id_thresh}, min_colinear={min_colinear}, "
              f"min_adjacent_junctions={min_adjacent_junctions}, window_bp={window_bp})")
        for e in merge_edges:
            need = min(min_adjacent_junctions, max(0, e["colinear_exons"] - 1))
            print(f"   merge fam{e['a']} + fam{e['b']}  chrom={e['chrom']} "
                  f"colinear={e['colinear_exons']} adjacent_junctions={e['adjacent_junctions']} "
                  f"need_junc={need} reason=share_gene "
                  f"genes_a={','.join(e['genes_a'][:5])}... "
                  f"genes_b={','.join(e['genes_b'][:5])}...")

    return merged, merge_edges


def split_crosschrom_domain_bridges(catalog, genes, gene_of_dn,
                                    id_thresh=CROSSCHROM_ID_THRESH,
                                    min_colinear=CROSSCHROM_MIN_COLINEAR,
                                    min_adjacent_junctions=CROSSCHROM_MIN_ADJACENT_JUNCTIONS,
                                    verbose=True):
    """Post-merge cross-chromosome domain-bridge split gate.

    For each family that spans more than one chromosome, keep cross-chromosome chromosome
    components connected only if there is strong colinear exon backbone support OR a shared
    (or related) annotated gene symbol.  Otherwise split into per-chromosome components.

    Strong support = best strict-LIS colinear shared-exon count >= ``min_colinear`` at
    ``id_thresh`` with adjacent junctions >= min(min_adjacent_junctions, col-1).  A higher
    identity threshold is required for cross-chrom claims because same-chrom duplication is
    the mechanistic default.
    """
    fa = _load_fasta_once()
    skel = _load_skeletons_once()
    _, strand = load_meta()
    rec_cache = {}

    new_catalog = []
    split_info = []
    for block in catalog:
        bychrom = defaultdict(list)
        for dn in block:
            L = genes.get(dn)
            if L:
                bychrom[L["chrom"]].append(dn)
        chroms = sorted(bychrom.keys())
        if len(chroms) <= 1:
            new_catalog.append(block)
            continue

        n = len(chroms)
        keep_edges = set()
        for i in range(n):
            for j in range(i + 1, n):
                c1, c2 = chroms[i], chroms[j]
                best_col, best_junc = 0, 0
                for dn1 in bychrom[c1]:
                    for dn2 in bychrom[c2]:
                        col, junc, _, _ = best_colinear_between_blocks(
                            [dn1], [dn2], genes, skel, strand, fa,
                            id_thresh=id_thresh, rec_cache=rec_cache)
                        if col > best_col or (col == best_col and junc > best_junc):
                            best_col, best_junc = col, junc

                syms1 = {gene_of_dn.get(dn) for dn in bychrom[c1] if gene_of_dn.get(dn)}
                syms2 = {gene_of_dn.get(dn) for dn in bychrom[c2] if gene_of_dn.get(dn)}
                nonloc1 = {s for s in syms1 if s and s != "NA" and not s.startswith("LOC")}
                nonloc2 = {s for s in syms2 if s and s != "NA" and not s.startswith("LOC")}

                sz1, sz2 = len(bychrom[c1]), len(bychrom[c2])
                # Multi-locus cross-chromosome clusters are kept; we only challenge links where
                # the smaller chromosome-side is a singleton (the classic domain-bridge pattern).
                singleton_bridge = min(sz1, sz2) == 1

                keep = False
                if not singleton_bridge:
                    keep = True
                else:
                    if nonloc1 & nonloc2:
                        keep = True
                    if (best_col >= min_colinear and
                            best_junc >= min(min_adjacent_junctions, max(0, best_col - 1))):
                        keep = True
                    if not keep:
                        roots1 = {_gene_symbol_root(s) for s in nonloc1}
                        roots2 = {_gene_symbol_root(s) for s in nonloc2}
                        if roots1 & roots2:
                            keep = True

                if keep:
                    keep_edges.add((i, j))

        components = _union_find_components(n, keep_edges)
        if len(components) == 1:
            new_catalog.append(block)
            continue

        for comp in components:
            new_block = []
            for ci in comp:
                new_block.extend(bychrom[chroms[ci]])
            new_catalog.append(sorted(set(new_block)))

        split_info.append(dict(
            original_n_loci=len(block),
            chroms=chroms,
            sizes={c: len(bychrom[c]) for c in chroms},
            kept_edges=[(chroms[i], chroms[j]) for i, j in keep_edges],
            split_into=[sorted(bychrom[chroms[ci]][0] for ci in comp) for comp in components]))

    if verbose:
        n_splits = len(split_info)
        n_loci_moved = sum(len(block) for block in new_catalog) - sum(len(block) for block in catalog)
        print(f"[family_merge_colinear crosschrom-split] {len(catalog)} blocks -> {len(new_catalog)} blocks; "
              f"splits={n_splits} (id_thresh={id_thresh}, min_colinear={min_colinear}, "
              f"min_adjacent_junctions={min_adjacent_junctions})")
        for si in split_info:
            print(f"   split n={si['original_n_loci']} chroms={si['chroms']} sizes={si['sizes']} "
                  f"kept={si['kept_edges']}")

    return new_catalog, split_info


if __name__ == "__main__":
    import csv
    CATALOG = os.path.join(BENCH, "family_rna_refine.tsv")

    def load_catalog():
        fam = defaultdict(list)
        with open(CATALOG) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                fam[row["family_id"]].append(row["member_dn"])
        return [sorted(v) for _, v in sorted(fam.items(), key=lambda x: int(x[0]))]

    import family_er_pr as FP
    import family_rna_refine as FR

    raw = FP.load_raw_families()
    meta = FP.load_meta()
    annot = FP.load_annot()
    gene_of = FP.gene_of_factory(annot)
    genes, gene_of_dn, *_ = FP.build_genes_dict(raw, meta, gene_of)
    gene_strand = FR.load_gene_strand() if hasattr(FR, "load_gene_strand") else {}
    pair_repeat_mult = FR.load_repeat_mult() if hasattr(FR, "load_repeat_mult") else {}

    catalog = load_catalog()
    merged, info = merge_catalog_colinear(
        catalog, genes, gene_of_dn, gene_strand, pair_repeat_mult,
        min_colinear=2, id_thresh=0.70, window_bp=5_000_000)
    print(f"\nmerged {len(catalog)} -> {len(merged)} families")
