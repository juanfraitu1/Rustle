#!/usr/bin/env python3
"""Can an ALIGNMENT-FREE intron-chain certificate create edges minimap2 never proposes?
Per `docs/PREREG_trie_edge_construction_2026-09-22.md` (§6y2, md5 `619a1999`).

§6u4 split edge-construction loss into GATE REJECTION (22.7%, already addressed by `--min-cov-shorter`)
and NO ALIGNMENT (26.1%). No alignment-derived signal can reach the second by construction; this measures
whether intron structure can.

Rule: two genes are joined iff they share a 3-intron length shingle at 5% tolerance (log-binned) --
r1024's operating point (precision 0.979 / recall 13.5% overall, but 99.3% redundant with the aligner).
⚠ A gene with < 3 exons has no 3-shingle and is UNREACHABLE by this rule; that ceiling is reported.
"""
import argparse
import collections
import csv
import itertools
import math
import re


def chains(gff, chrom):
    """gene -> intron-length chain of its longest annotated transcript, oriented 5'->3'."""
    tx = collections.defaultdict(list)
    tx_gene = {}
    for ln in open(gff):
        if ln.startswith('#'):
            continue
        f = ln.rstrip('\n').split('\t')
        if len(f) < 9 or f[0] != chrom:
            continue
        if f[2] == 'exon':
            p = re.search(r'Parent=([^;]+)', f[8])
            if p:
                tx[p.group(1)].append((int(f[3]), int(f[4]), f[6]))
        elif f[2] in ('mRNA', 'transcript'):
            i = re.search(r'ID=([^;]+)', f[8])
            g = re.search(r'gene=([^;]+)', f[8]) or re.search(r'Parent=gene-([^;]+)', f[8])
            if i and g:
                tx_gene[i.group(1)] = g.group(1)
    out = {}
    for t, ex in tx.items():
        g = tx_gene.get(t)
        if not g or len(ex) < 3:
            continue
        ex.sort()
        iv = [ex[i + 1][0] - ex[i][1] - 1 for i in range(len(ex) - 1)]
        if ex[0][2] == '-':
            iv = iv[::-1]
        if g not in out or len(iv) > len(out[g]):
            out[g] = iv
    return out


def spans(gff, chrom):
    out = {}
    for ln in open(gff):
        if ln.startswith('#'):
            continue
        f = ln.rstrip('\n').split('\t')
        if len(f) < 9 or f[0] != chrom or f[2] not in ('gene', 'pseudogene', 'ncRNA_gene'):
            continue
        m = re.search(r'Name=([^;]+)', f[8])
        if m:
            out[m.group(1)] = (int(f[3]) - 1, int(f[4]))
    return out


def aligned_pairs(paf, gspan, chrom):
    """gene pairs minimap2 proposes ANY alignment for (node header = chrom:start-end)."""
    node = {}
    for g, (s, e) in gspan.items():
        node[f'{chrom}:{s}-{e}'] = g
        node[f'{chrom}:{s + 1}-{e}'] = g
    out = set()
    for ln in open(paf):
        f = ln.split('\t', 6)
        a, b = node.get(f[0]), node.get(f[5])
        if a and b and a != b:
            out.add((a, b) if a < b else (b, a))
    return out


def main():
    ap = argparse.ArgumentParser()
    for a in ('--gff', '--paf', '--truth'):
        ap.add_argument(a, required=True)
    ap.add_argument('--chrom', required=True)
    ap.add_argument('--tol', type=float, default=0.05)
    ap.add_argument('--k', type=int, default=3)
    a = ap.parse_args()

    ch = chains(a.gff, a.chrom)
    gsp = spans(a.gff, a.chrom)
    aln = aligned_pairs(a.paf, gsp, a.chrom)
    lab = {r['Gene Name']: r['Family ID']
           for r in csv.DictReader(open(a.truth), delimiter='\t') if r.get('Family ID')}

    # every referee-labelled gene, whether or not it has a usable chain -- the ceiling is part of the result
    scored = sorted(g for g in lab if g in gsp)
    withchain = [g for g in scored if g in ch and len(ch[g]) >= a.k]
    b = lambda x: int(math.log(max(x, 1)) / math.log(1 + a.tol))
    sh = {g: {tuple(b(ch[g][i + j]) for j in range(a.k)) for i in range(len(ch[g]) - a.k + 1)}
          for g in withchain}

    truth_pairs = {(x, y) for x, y in itertools.combinations(scored, 2) if lab[x] == lab[y]}
    noaln_truth = {p for p in truth_pairs if p not in aln}
    rec = {p for p in noaln_truth if p[0] in sh and p[1] in sh and sh[p[0]] & sh[p[1]]}
    # false positives among UNALIGNED pairs: different family, no alignment, certificate fires
    fp = 0
    for x, y in itertools.combinations(withchain, 2):
        if lab[x] == lab[y]:
            continue
        if (x, y) in aln:
            continue
        if sh[x] & sh[y]:
            fp += 1
    prec = len(rec) / (len(rec) + fp) if (len(rec) + fp) else float('nan')
    reachable = {p for p in noaln_truth if p[0] in sh and p[1] in sh}
    print(f"{a.chrom}: referee genes {len(scored)} (with a >={a.k}-intron chain: {len(withchain)})")
    print(f"  truth same-family pairs            {len(truth_pairs)}")
    print(f"  of which NO ALIGNMENT (the target) {len(noaln_truth)}  "
          f"[reachable by the rule: {len(reachable)} = {100*len(reachable)/max(len(noaln_truth),1):.1f}%]")
    print(f"  ⭐ recovered by the certificate     {len(rec)} = {100*len(rec)/max(len(noaln_truth),1):.1f}% of target")
    print(f"  false joins among UNALIGNED pairs  {fp}")
    print(f"  ⭐ precision on unaligned pairs     {prec:.3f}")


if __name__ == '__main__':
    main()
