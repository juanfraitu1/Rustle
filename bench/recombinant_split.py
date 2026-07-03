#!/usr/bin/env python
"""Recombinant-split gate -- split families held together by a RECOMBINANT/mosaic bridge.

WHAT THIS IS
------------
A RECOMBINANT bridge is a locus whose exon-chain splits COLINEARLY into >=2 DIFFERENT
sub-families' exons (different exons -> different sub-families) AND is an articulation point
whose best single-neighbour colinear coverage < COLINEAR_COV (0.60).  The shipped pairwise-edge
definition transitively merges G1-G3-G2 through such a bridge (a density-2/3 barbell survives the
gamma-quasi-clique); this gate SPLITS the family along the bridge into its colinear sub-families.

Validated in bench/RECOMBINATION_BRIDGE_DETECTOR.md (commit 439d2ae): the detector flagged 5
recombinant blocks (94, 110, 187, 210, 321), but only the HIGH-CONFIDENCE DISTRIBUTED-MOSAIC
subset is RECALL-SAFE to split.  The shipped gate additionally requires (per the detector's own
confidence tier) BOTH bridged sub-families ANNOTATED, minor block >= 2 exons, and NOT genomically
colocated -- which splits exactly the 2 clean distinct-gene mosaics (fid 210 GALNT17|LOC101126070
+ fid 187) at 0 oracle-recall cost, and correctly LEAVES ALONE the 3 unsafe blocks: fid 321 (a
READTHROUGH-chimera whose split would fragment RDH14x3) and fid 110 (an NA-dominated OVERSIZE blob
whose split would strand LOC134758618 -- a DNA-multi-copy-but-RNA-single-copy spurious co-membership
credit).  Two recall-safety guarantees (both RNA-only, no oracle): (1) never separate same-gene loci;
(2) HIGH-confidence-only.  Flagship = fid 210: bridge LOC115935264 -> LOC101126070 | GALNT17.

REUSE (nothing re-derived): bench/recombination_bridge_detector.py -- family_exons,
exon_match_tensor, build_graph, colinear_cov, ID_THRESH, COLINEAR_COV, and the per-bridge
RECOMBINANT classification (different-exons->different-subfams).  The colinear/exon logic is the
detector's exact, refuter-checked implementation.

RNA-ONLY: the split decision uses ONLY exon structure (denovo_skeletons introns) + genome-base
exon homology (the SAME provenance as the shipped aln_frac).  No DNA/protein/soft-mask column.

DETERMINISM: sorted articulation scan, sorted components, deterministic bridge/orphan assignment.
"""
import collections

import networkx as nx

import recombination_bridge_detector as RBD

COLINEAR_COV = RBD.COLINEAR_COV       # 0.60 -- RECOMBINATION_BRIDGE_DETECTOR.md; do NOT re-derive
ID_THRESH = RBD.ID_THRESH             # 0.70 exon-match identity
MAX_LOCI = 60                         # skip giant families (repeat hubs, already repeat-gated; never mosaics)
assert abs(COLINEAR_COV - 0.60) < 1e-9, "COLINEAR_COV drifted from the validated 0.60"

# --------------------------------------------------------------------------- cached loaders
_SKEL = _STRAND = _FA = None


def loaders():
    """Load (and cache) the RNA exon-structure + genome inputs the split needs.  Raises a clear
    error (mirroring the repeat gate) if an input is missing so the caller can --no-split-recombinants."""
    global _SKEL, _STRAND, _FA
    if _SKEL is None:
        try:
            import pysam
        except ImportError as e:
            raise RuntimeError("recombinant-split gate is DEFAULT-ON but pysam is missing; "
                               "install it or disable with --no-split-recombinants / "
                               "RUSTLE_NO_SPLIT_RECOMBINANTS=1") from e
        import os
        for p in (RBD.SKEL_TSV, RBD.META_TSV, RBD.GENOME):
            if not os.path.exists(p):
                raise FileNotFoundError(
                    f"recombinant-split gate is DEFAULT-ON but {p} is missing; "
                    f"disable with --no-split-recombinants / RUSTLE_NO_SPLIT_RECOMBINANTS=1")
        _SKEL = RBD.load_skeletons()
        _STRAND = RBD.load_strand()
        _FA = pysam.FastaFile(RBD.GENOME)
    return _SKEL, _STRAND, _FA


# --------------------------------------------------------------------------- recombinant bridge detection
def _members(block, genes, gene_of_dn):
    ms = []
    for dn in sorted(block):
        g = genes.get(dn)
        if g is None:
            continue
        ms.append(dict(dn=dn, gene=gene_of_dn.get(dn, "NA"),
                       chrom=g["chrom"], start=g["start"], end=g["end"]))
    return ms


def _recombinant_bridges(recs, best, G):
    """Indices of vg-split RECOMBINANT articulation bridges.  Replicates analyze_family's per-bridge
    test exactly: an articulation point that separates its neighbours into >=2 sub-families, where
    DIFFERENT exons hit DIFFERENT sub-families (>=2 singly-hit sub-families), and NO single neighbour
    colinearly covers >= COLINEAR_COV of the bridge (=> a mosaic, not one shared domain)."""
    if G.number_of_edges() == 0:
        return []
    out = []
    for b in sorted(nx.articulation_points(G)):
        Gb = G.copy()
        Gb.remove_node(b)
        nbrs = set(G.neighbors(b))
        bridged = sorted((sorted(c) for c in nx.connected_components(Gb) if set(c) & nbrs),
                         key=lambda c: c[0])
        if len(bridged) < 2:
            continue
        cid_of = {m: ci for ci, comp in enumerate(bridged) for m in comp}
        single = collections.Counter()
        for i in range(len(recs[b]["exons"])):
            hit = set()
            for m in nbrs:
                idv, _ = best.get((b, i, m), (0.0, -1))
                if idv >= ID_THRESH and m in cid_of:
                    hit.add(cid_of[m])
            if len(hit) == 1:
                single[next(iter(hit))] += 1
        if len(single) < 2:                    # NOT recombinant (needs different exons -> different subfams)
            continue
        best_single = max((RBD.colinear_cov(recs, best, b, m) for m in nbrs), default=0.0)
        if best_single >= COLINEAR_COV:        # some neighbour colinearly covers the bridge -> not a mosaic
            continue
        # HIGH-confidence DISTRIBUTED-MOSAIC only (RECOMBINATION_BRIDGE_DETECTOR.md confidence tier):
        # both bridged sub-families ANNOTATED + minor block >= 2 exons + NOT genomically colocated
        # (a colocated pair = READTHROUGH-CHIMERA assembly artifact; an NA-dominated sub-family =
        # an oversize blob, e.g. fid 110 LOC134758618/NLGN1/NA -- NOT a clean distinct-gene mosaic).
        # This is what separates the clean fid 210 (GALNT17|LOC101126070) from the messy oversize blobs
        # whose spurious co-membership credit is a DNA-multi-copy-but-RNA-single-copy recall trap.
        if RBD.subfams_colocated(recs, bridged):
            continue                           # readthrough chimera (assembly artifact), not distributed
        involved = sorted(single)              # sub-families each SOLELY hit by >=1 bridge exon
        minor = min(single[c] for c in involved)
        genes_involved = [RBD.subfam_dom_gene(recs, bridged[c]) for c in involved]
        both_annotated = all(g not in ("NA", "") for g in genes_involved)
        if not (both_annotated and minor >= 2):
            continue                           # LOW confidence (NA-dominant sub-family / 1-exon minor block)
        out.append(b)
    return out


def split_block(block, genes, gene_of_dn, skel, strand, fa):
    """Return a PARTITION of `block` (list of DN-blocks): split along vg-split recombinant bridges,
    or [sorted(block)] unchanged if none.  Never loses a locus.  Deterministic."""
    block = list(block)
    if len(block) < 3:
        return [sorted(block)]
    members = _members(block, genes, gene_of_dn)
    if len(members) < 3 or len(members) > MAX_LOCI:
        return [sorted(block)]
    recs = RBD.family_exons(members, skel, strand, fa)
    if len(recs) < 3:
        return [sorted(block)]
    best = RBD.exon_match_tensor(recs)
    G, _ = RBD.build_graph(recs, best)
    bridge_idx = _recombinant_bridges(recs, best, G)
    if not bridge_idx:
        return [sorted(block)]

    dn_of = {i: recs[i]["dn"] for i in range(len(recs))}
    bset = set(bridge_idx)
    # sub-families = connected components after removing ALL vg-split bridges
    Gs = G.copy()
    Gs.remove_nodes_from(bridge_idx)
    subfams = sorted((sorted(c) for c in nx.connected_components(Gs) if c), key=lambda c: c[0])
    subfam_dns = [set(dn_of[i] for i in sf) for sf in subfams]

    # assign each bridge to its best-colinear sub-family (>= COLINEAR_COV to a member) else own singleton
    for b in bridge_idx:
        bestc, bestsf = 0.0, -1
        for si, sf in enumerate(subfams):
            c = max((RBD.colinear_cov(recs, best, b, m) for m in sf), default=0.0)
            if c > bestc:
                bestc, bestsf = c, si
        if bestsf >= 0 and bestc >= COLINEAR_COV:
            subfam_dns[bestsf].add(dn_of[b])
        else:
            subfam_dns.append({dn_of[b]})       # mosaic bridge -> its own singleton

    # skeleton-less DNs (not in recs): attach to the first sub-family sharing a chromosome, else singleton
    placed = set()
    for s in subfam_dns:
        placed |= s
    for dn in sorted(block):
        if dn in placed:
            continue
        g = genes.get(dn)
        tgt = -1
        if g is not None:
            for si in range(len(subfam_dns)):
                if any(genes.get(x, {}).get("chrom") == g["chrom"] for x in sorted(subfam_dns[si])):
                    tgt = si
                    break
        if tgt >= 0:
            subfam_dns[tgt].add(dn)
        else:
            subfam_dns.append({dn})
        placed.add(dn)

    result = sorted((sorted(s) for s in subfam_dns if s), key=lambda c: c[0])

    # RECALL-SAFETY GUARD (provable): a genuine recombinant mosaic bridges DISTINCT genes; if the
    # partition would separate loci of the SAME gene into different sub-families, this is NOT a clean
    # recombinant over-merge -- it is breaking real same-gene copies (a recall loss, e.g. a readthrough
    # like fid 321 fragmenting RDH14x3, or IFITM1x2).  REJECT the whole split in that case (keep merged).
    # Because same-gene loci never get separated, no gene's max-loci-in-one-family can fall -> no
    # multi-copy family recovery is ever lost (0 recall cost, by construction, not just on controls).
    if len(result) > 1:
        gene_parts = collections.defaultdict(set)
        for pi, part in enumerate(result):
            for dn in part:
                g = gene_of_dn.get(dn)
                if g not in (None, "NA"):
                    gene_parts[g].add(pi)
        if any(len(pis) >= 2 for pis in gene_parts.values()):
            return [sorted(block)]        # would separate same-gene copies -> not a clean mosaic -> keep merged
    return result


def split_families(refined, genes, gene_of_dn):
    """Apply the recombinant-split gate to a list of refined family blocks.  Returns
    (new_refined, split_info) where split_info lists the families that were split (bridge -> subfams).
    Loads the RNA/genome inputs lazily (cached)."""
    skel, strand, fa = loaders()
    new_refined = []
    split_info = []
    for b in refined:
        parts = split_block(b, genes, gene_of_dn, skel, strand, fa)
        if len(parts) > 1:
            doms = []
            for p in parts:
                gs = [gene_of_dn.get(dn) for dn in p if gene_of_dn.get(dn) not in (None, "NA")]
                doms.append(collections.Counter(gs).most_common(1)[0][0] if gs else "NA")
            split_info.append(dict(n_loci=sum(len(p) for p in parts),
                                   n_parts=len(parts), subfam_genes=doms))
        new_refined.extend(parts)
    return new_refined, split_info
