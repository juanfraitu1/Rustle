#!/usr/bin/env python3
"""Candidates of the template-guided copy stage: projected template hits at positions no base locus covers.

Spec section 1, steps 2-3. A hit (minimap2 -x splice, PAF with `cg:Z:`) is kept when its identity (matches / block
length) is >= 0.80 and >= 50% of the template is aligned: the guided finder's registered thresholds (`guided_pipeline`
MIN_ID, MIN_COV; figure 8, G1). Its projected exons are its CIGAR's blocks split at `N` (`guided_pipeline.tx_exon_blocks`).
A hit is a candidate when every base locus covers < 50% of its projected exon union (figure 8's rule); the template's own
locus covers its own position, so self-hits drop out. Candidates on one contig and strand whose unions overlap by >= 50% of
the smaller union merge; the highest-identity hit is kept.
"""
import collections
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(HERE, "..", "..", "figures"))
import guided_pipeline as GP  # noqa: E402
import _liftoff as LF  # noqa: E402

Hit = collections.namedtuple("Hit", "template contig strand exons identity aligned")


def hits(paf_path, min_id=GP.MIN_ID, min_cov=GP.MIN_COV):
    """PAF -> [Hit] passing the identity and aligned-fraction thresholds (hits without a CIGAR are skipped)."""
    out = []
    for h in GP.parse_paf(paf_path):
        if not h["cg"] or not h["bl"] or not h["qlen"]:
            continue
        ident, aligned = h["nm"] / h["bl"], (h["qe"] - h["qs"]) / h["qlen"]
        if ident >= min_id and aligned >= min_cov:
            out.append(Hit(h["q"], h["chrom"], h["strand"], GP.tx_exon_blocks(h), ident, aligned))
    return out


def uncovered(hit_list, base_rows, min_cov=LF.MATCH_COV):
    """The hits whose projected exon union is covered < min_cov by every base locus."""
    idx = LF.Index(base_rows)
    return [h for h in hit_list if idx.best_cov(h.contig, h.exons)[0] < min_cov]


def _rank(h):
    return (-h.identity, -h.aligned, h.template, h.contig, h.exons[0][0])


def merge_candidates(hit_list):
    """Groups of hits on one contig and strand whose exon unions overlap by >= 50% of the smaller union (transitively).
    -> [(kept Hit, [the other Hits])], kept = the highest identity, then aligned fraction, then template id; sorted by
    contig, strand and start."""
    order = sorted(range(len(hit_list)), key=lambda i: (hit_list[i].contig, hit_list[i].strand,
                                                        hit_list[i].exons[0][0], hit_list[i].exons[-1][1]))
    parent = list(range(len(hit_list)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    active = []
    for i in order:
        h = hit_list[i]
        active = [j for j in active if hit_list[j].contig == h.contig and hit_list[j].strand == h.strand
                  and hit_list[j].exons[-1][1] > h.exons[0][0]]
        for j in active:
            o = hit_list[j]
            inter = LF.iv_inter(h.exons, o.exons)
            if inter and inter >= 0.5 * min(LF.iv_len(h.exons), LF.iv_len(o.exons)):
                parent[find(i)] = find(j)
        active.append(i)
    groups = collections.defaultdict(list)
    for i in range(len(hit_list)):
        groups[find(i)].append(hit_list[i])
    out = []
    for g in groups.values():
        g = sorted(g, key=_rank)
        out.append((g[0], g[1:]))
    return sorted(out, key=lambda kg: (kg[0].contig, kg[0].strand, kg[0].exons[0][0], kg[0].template))
