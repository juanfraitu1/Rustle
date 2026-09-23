#!/usr/bin/env python3
"""Score one node-construction MODE's clusters against Soto families, per
`docs/PREREG_semiguided_npip_2026-09-22.md` (§6x0).

The three modes (de novo / guided / semi-guided) differ ONLY in what an interval is. This scores any of
their `mcl_families` cluster files the same way so the comparison is about node construction and nothing
else.

  sensitivity = matched truth members / truth members
  precision   = matched truth members / members of the matched clusters
  F           = harmonic mean
  matching    = ONE-TO-ONE bipartite, maximising overlap (scipy linear_sum_assignment)
  collapse    = truth genes that lost their best locus to another truth gene (register 817's failure mode:
                a coarse interval swallows several genes, and one-to-one matching charges for every one
                beyond the first)

⚠ Locus -> gene is MAX-OVERLAP, ONE NAME PER LOCUS. Register 964/§6v9 showed multi-labelling is worse,
and §6v8 showed winner-take-all undercounts recall -- but the bias is identical in all three arms, which
is what a between-arm comparison needs.
⚠ Symbol-root truth is VOID (register 902). Truth is Soto `Family ID` only.

Usage:
  mode_family_score.py --clusters X.clusters.tsv --gff chr16.genes.gff --soto soto_famCN_S1C.tsv \
      --chrom chr16 [--family NPIP] [--label de-novo]
"""
import argparse
import collections
import csv
import re

import numpy as np
from scipy.optimize import linear_sum_assignment


def soto_truth(path):
    fam = {}
    for r in csv.DictReader(open(path), delimiter='\t'):
        f = (r.get('Family ID') or '').strip()
        g = (r.get('Gene Name') or '').strip()
        if f and f != 'N/A' and g:
            fam.setdefault(g, f)
    return fam


def gene_spans(gff, chrom):
    out = []
    for ln in open(gff):
        if ln.startswith('#'):
            continue
        f = ln.rstrip('\n').split('\t')
        if len(f) < 9 or f[0] != chrom or f[2] not in ('gene', 'pseudogene', 'ncRNA_gene'):
            continue
        m = re.search(r'Name=([^;]+)', f[8])
        if m:
            out.append((int(f[3]), int(f[4]), m.group(1)))
    out.sort()
    return out


def load_clusters(path, chrom):
    """cluster_id -> [(start, end)] for members on `chrom`."""
    rows = [ln.rstrip('\n').split('\t') for ln in open(path)]
    hdr = rows[0]
    ci = {h: i for i, h in enumerate(hdr)}
    need = ('cluster_id', 'chrom', 'start', 'end')
    if not all(k in ci for k in need):
        raise SystemExit(f"{path}: expected columns {need}, got {hdr[:8]}")
    out = collections.defaultdict(list)
    for r in rows[1:]:
        if len(r) < len(hdr) or r[ci['chrom']] != chrom:
            continue
        out[r[ci['cluster_id']]].append((int(r[ci['start']]), int(r[ci['end']])))
    return out


def main():
    ap = argparse.ArgumentParser()
    for a in ('--clusters', '--gff', '--soto'):
        ap.add_argument(a, required=True)
    ap.add_argument('--chrom', default='chr16')
    ap.add_argument('--family', default=None, help='restrict truth to Soto families containing this prefix')
    ap.add_argument('--label', default='arm')
    a = ap.parse_args()

    fam = soto_truth(a.soto)
    spans = gene_spans(a.gff, a.chrom)
    clusters = load_clusters(a.clusters, a.chrom)

    # locus -> ONE gene, by max overlap
    locus_gene = {}
    gene_claims = collections.defaultdict(list)
    for cid, members in clusters.items():
        for (s, e) in members:
            best = None
            for (gs, ge, g) in spans:
                if ge < s:
                    continue
                if gs > e:
                    break
                ov = min(e, ge) - max(s, gs)
                if ov > 0 and (best is None or ov > best[0]):
                    best = (ov, g)
            if best:
                locus_gene[(cid, s, e)] = best[1]
                gene_claims[best[1]].append((cid, s, e))

    # truth families on this chromosome
    on_chrom = {g for (_, _, g) in spans}
    truth = collections.defaultdict(set)
    for g, f in fam.items():
        if g in on_chrom:
            truth[f].add(g)
    if a.family:
        truth = {f: gs for f, gs in truth.items() if any(a.family in g for g in gs)}
    truth = {f: gs for f, gs in truth.items() if len(gs) >= 2}

    # predicted clusters as gene sets
    pred = {}
    for cid, members in clusters.items():
        gs = {locus_gene[(cid, s, e)] for (s, e) in members if (cid, s, e) in locus_gene}
        gs &= set().union(*truth.values()) if truth else set()
        if gs:
            pred[cid] = gs
    if not truth or not pred:
        print(f"{a.label}: no scoreable truth/prediction overlap")
        return

    T = list(truth); P = list(pred)
    M = np.zeros((len(T), len(P)))
    for i, tf in enumerate(T):
        for j, pc in enumerate(P):
            M[i, j] = len(truth[tf] & pred[pc])
    ri, cj = linear_sum_assignment(-M)
    matched = sum(M[i, j] for i, j in zip(ri, cj))
    tot_truth = sum(len(truth[t]) for t in T)
    tot_pred = sum(len(pred[P[j]]) for i, j in zip(ri, cj) if M[i, j] > 0)
    sens = matched / tot_truth if tot_truth else 0.0
    prec = matched / tot_pred if tot_pred else 0.0
    f1 = 2 * sens * prec / (sens + prec) if (sens + prec) else 0.0

    # collapse (register 817's failure mode): map each TRUTH GENE to the locus that best covers it, then
    # count genes that must share one locus. ⚠ It has to be computed gene->locus. Doing it locus->gene is
    # 0 BY CONSTRUCTION, because the max-overlap resolver gives every locus exactly one gene.
    truth_genes = set().union(*truth.values())
    all_loci = [(cid, s, e) for cid, ms in clusters.items() for (s, e) in ms]
    gene_span = {g: (gs, ge) for (gs, ge, g) in spans}
    best_locus = {}
    for g in truth_genes:
        if g not in gene_span:
            continue
        gs, ge = gene_span[g]
        best = None
        for (cid, s, e) in all_loci:
            ov = min(ge, e) - max(gs, s)
            if ov > 0 and (best is None or ov > best[0]):
                best = (ov, (cid, s, e))
        if best:
            best_locus[g] = best[1]
    share = collections.defaultdict(set)
    for g, L in best_locus.items():
        share[L].add(g)
    collapsed = sum(len(v) - 1 for v in share.values() if len(v) > 1)
    missing = len(truth_genes - set(best_locus))

    print(f"{a.label:>14} | truth {len(T)} fams / {tot_truth} genes | clusters {len(pred)} "
          f"| sens {sens:.3f} prec {prec:.3f} F {f1:.3f} | collapsed {collapsed} | no-locus {missing}")


if __name__ == '__main__':
    main()
