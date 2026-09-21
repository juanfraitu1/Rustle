#!/usr/bin/env python3
"""Us vs Soto, scored against a NEUTRAL referee.

Every number this session has quoted used Soto as the truth, so it measures agreement with Soto, not
precision. To ask whether Soto is more precise ANYWHERE, both have to be scored against a third party.

Referee: **protein families** (§6ko's rule — longest CDS per gene, translated, all-vs-all blastp
e <= 1e-5, edge iff non-overlapping HSPs cover >= 0.30 of the longer protein, MCL I = 2.8, r2
exclusions). It is independent of BOTH comparators: it never sees our genomic alignment gate, and it
never sees Soto's SD/WSSD construction. It is amino-acid evidence about the product.

⚠ Register T15: "never consume the comparator's own files and call it replication" — nothing here reads
a Soto-derived file except the S1C family assignment being SCORED, which is the object under test.

Reports, pooled and stratified by referee-family size: pairwise precision / recall / F for each
comparator against the referee, where a "pair" is two genes the comparator places together.

Usage: soto_vs_us_referee.py --gff G --genome FA --soto S1C --clusters DIR --chroms chr2,chr8,chr10
"""
import argparse
import collections
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mcl_port  # noqa: E402
from protein_edge_gap import longest_cds, protein_edges, translate  # noqa: E402
from protein_families import excluded  # noqa: E402
from neighbourhood_jaccard import gene_names  # noqa: E402


def protein_referee(gff, genome, chrom, workdir):
    import pysam
    import re
    cds = longest_cds(gff, chrom)
    bt = {}
    for line in open(gff):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[0] != chrom or f[2] not in ('gene', 'pseudogene'):
            continue
        n = re.search(r'Name=([^;]+)', f[8]); b = re.search(r'gene_biotype=([^;]+)', f[8])
        if n:
            bt[n.group(1)] = b.group(1) if b else ''
    cds = {g: v for g, v in cds.items() if not excluded(bt.get(g, ''), 2)}
    fa = pysam.FastaFile(genome)
    out = f'{workdir}/{chrom}_ref'
    faa = out + '.proteins.faa'
    plen = {}
    if not os.path.exists(faa):
        with open(faa, 'w') as fh:
            for g, (st, segs) in cds.items():
                p = translate(fa, chrom, st, segs)
                if len(p) >= 10:
                    plen[g] = len(p); fh.write(f'>{g}\n{p}\n')
    if not plen:
        n = None
        for line in open(faa):
            if line.startswith('>'):
                n = line[1:].strip()
            elif n:
                plen[n] = len(line.strip()); n = None
    pe = protein_edges(faa, out, plen, '4')
    fams = mcl_port.mcl({(x, y): 1.0 for x, y in pe}, inflation=2.8)
    lab = {}
    for i, mem in enumerate(fams):
        mem = sorted(set(mem))
        if len(mem) >= 2:
            for g in mem:
                lab[g] = f'PF{i}'
    return lab


def pair_scores(label_of, ref):
    """precision/recall/F against the referee over a FIXED UNIVERSE: every referee-labelled gene is
    scored, and a gene the comparator never placed becomes its own singleton rather than being dropped.

    ⚠ Restricting the universe to genes the comparator labelled conditions the denominator on the
    prediction (register 770) and returns precision 1.000 by construction. Do not do that.
    """
    genes = list(ref)
    label_of = {g: label_of.get(g, f'__singleton__{g}') for g in genes}
    byc = collections.defaultdict(list)
    for g in genes:
        byc[label_of[g]].append(g)
    byr = collections.defaultdict(list)
    for g in genes:
        byr[ref[g]].append(g)
    pred = set()
    for v in byc.values():
        for i in range(len(v)):
            for j in range(i + 1, len(v)):
                pred.add(tuple(sorted((v[i], v[j]))))
    true = set()
    for v in byr.values():
        for i in range(len(v)):
            for j in range(i + 1, len(v)):
                true.add(tuple(sorted((v[i], v[j]))))
    tp = len(pred & true)
    p = tp / len(pred) if pred else 0.0
    r = tp / len(true) if true else 0.0
    f = 0.0 if p + r == 0 else 2 * p * r / (p + r)
    return p, r, f, len(pred), len(true), tp


def main():
    ap = argparse.ArgumentParser()
    for x in ('--gff', '--genome', '--soto', '--clusters', '--chroms'):
        ap.add_argument(x, required=True)
    ap.add_argument('--workdir', default='/mnt/linuxdisk/tmp/referee')
    a = ap.parse_args()
    os.makedirs(a.workdir, exist_ok=True)

    ids = collections.defaultdict(set)
    for r in csv.DictReader(open(a.soto), delimiter='\t'):
        f = (r.get('Family ID') or '').strip(); n = (r.get('Gene Name') or '').strip()
        if f and f != 'N/A' and n:
            ids[n].add(f)
    soto = {n: next(iter(v)) for n, v in ids.items() if len(v) == 1}

    ours, sot, ref = {}, {}, {}
    for chrom in a.chroms.split(','):
        names = gene_names(a.gff, chrom)
        lab = protein_referee(a.gff, a.genome, chrom, a.workdir)
        for g, f in lab.items():
            ref[f'{chrom}:{g}'] = f'{chrom}:{f}'
        p = f'{a.clusters}/{chrom}_fam.clusters.tsv'
        for line in open(p):
            if line.startswith('cluster_id'):
                continue
            q = line.rstrip('\n').split('\t')
            g = names.get(f'{q[5]}:{q[6]}-{q[7]}')
            if g:
                ours[f'{chrom}:{g}'] = f'{chrom}:{q[0]}'
        for sp, g in names.items():
            if g in soto:
                sot[f'{chrom}:{g}'] = soto[g]

    print(f"referee: protein families — {len(set(ref.values()))} families over {len(ref)} genes\n")
    print(f"  {'comparator':12s} {'precision':>10} {'recall':>8} {'F':>8} {'pairs called':>13} {'TP':>6}")
    for name, lab in (('OURS (MCL)', ours), ('SOTO', sot)):
        p, r, f, np_, nt, tp = pair_scores(lab, ref)
        print(f"  {name:12s} {p:>10.3f} {r:>8.3f} {f:>8.3f} {np_:>13} {tp:>6}")
    print(f"\n  (referee pairs available: {pair_scores(ours, ref)[4]})")

    print("\nstratified by REFEREE family size — 'in any part':")
    print(f"  {'ref fam size':>13} {'genes':>6} {'OURS prec':>10} {'SOTO prec':>10} {'OURS rec':>9} {'SOTO rec':>9}")
    byr = collections.defaultdict(list)
    for g, f in ref.items():
        byr[f].append(g)
    for lo, hi, lbl in ((2, 2, '2'), (3, 4, '3-4'), (5, 9, '5-9'), (10, 10**9, '>=10')):
        keep = {g for f, v in byr.items() if lo <= len(v) <= hi for g in v}
        sub = {g: f for g, f in ref.items() if g in keep}
        if not sub:
            continue
        po, ro = pair_scores({g: v for g, v in ours.items() if g in keep}, sub)[:2]
        ps, rs = pair_scores({g: v for g, v in sot.items() if g in keep}, sub)[:2]
        print(f"  {lbl:>13} {len(keep):>6} {po:>10.3f} {ps:>10.3f} {ro:>9.3f} {rs:>9.3f}")


if __name__ == '__main__':
    main()
