#!/usr/bin/env python3
"""Arm C of `docs/PREREG_jaccard_sweep_2026-09-20.md`: Jaccard threshold + the exon conjunct, then
connected components. No MCL.

§6t3 measured that the clustering operator buys +0.035 while the exon conjunct buys +0.062, so this arm
asks the obvious follow-up: keep the cheap biological conjunct, drop the expensive operator.

The conjunct, as `mcl_families` documents it: *the pair's best record must cover this FRACTION of the
SMALLER gene's exonic length with shared exon-to-exon evidence*.

⚠ APPROXIMATION, stated plainly: the faithful version needs per-base query→target correspondence (the
CIGAR) to require a base be exonic on BOTH sides. Here the aligned interval is projected onto each gene's
exons independently and the smaller of the two exonic coverages is used:

    f_ex ~= min(exonic bases of A under the alignment, exonic bases of B under the alignment)
            / min(exonic length of A, exonic length of B)

That is an UPPER bound on the true shared-exon fraction (it does not check that the same alignment
columns are exonic on both sides), so this arm is if anything flattered. Any adoption would have to
recompute it from the CIGAR inside `mcl_families`.

Usage: jaccard_plus_exon.py --paf P --gff G --chrom chrN --out PREFIX --t 0.50 --fex 0.60
"""
import argparse
import collections
import re


def exons_by_gene(gff, chrom):
    """(start1,end) span key -> merged exon intervals, for gene/pseudogene records with a Name."""
    span_of_name, name_of_span = {}, {}
    for line in open(gff):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[0] != chrom or f[2] not in ('gene', 'pseudogene'):
            continue
        n = re.search(r'Name=([^;]+)', f[8])
        if n:
            k = (int(f[3]), int(f[4]))
            span_of_name[n.group(1)] = k; name_of_span[k] = n.group(1)
    ex = collections.defaultdict(list)
    for line in open(gff):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[0] != chrom or f[2] != 'exon':
            continue
        g = re.search(r'gene=([^;]+)', f[8])
        if g and g.group(1) in span_of_name:
            ex[span_of_name[g.group(1)]].append((int(f[3]), int(f[4])))
    merged = {}
    for k, v in ex.items():
        v.sort(); out = []
        for s, e in v:
            if out and s <= out[-1][1] + 1:
                out[-1] = (out[-1][0], max(out[-1][1], e))
            else:
                out.append((s, e))
        merged[k] = out
    return merged


def overlap(iv, lo, hi):
    return sum(max(0, min(e, hi) - max(s, lo) + 1) for s, e in iv)


def main():
    ap = argparse.ArgumentParser()
    for x in ('--paf', '--gff', '--chrom', '--out'):
        ap.add_argument(x, required=True)
    ap.add_argument('--t', type=float, required=True)
    ap.add_argument('--fex', type=float, required=True)
    a = ap.parse_args()

    ex = exons_by_gene(a.gff, a.chrom)
    exlen = {k: sum(e - s + 1 for s, e in v) for k, v in ex.items()}

    def key(nm):
        try:
            ch, rng = nm.rsplit(':', 1); s, e = rng.split('-')
            return (int(s), int(e))
        except ValueError:
            return None

    best = {}
    for line in open(a.paf):
        f = line.rstrip('\n').split('\t')
        if len(f) < 11:
            continue
        ka, kb = key(f[0]), key(f[5])
        if not ka or not kb or ka == kb:
            continue
        la, lb, nm = int(f[1]), int(f[6]), int(f[9])
        den = la + lb - nm
        if den <= 0:
            continue
        j = nm / den
        if j < a.t:
            continue
        # aligned interval projected to genomic coordinates on each side
        qa0, qa1 = ka[0] + int(f[2]), ka[0] + int(f[3])
        qb0, qb1 = kb[0] + int(f[7]), kb[0] + int(f[8])
        ea = overlap(ex.get(ka, []), qa0, qa1)
        eb = overlap(ex.get(kb, []), qb0, qb1)
        den2 = min(exlen.get(ka, 0), exlen.get(kb, 0))
        fex = (min(ea, eb) / den2) if den2 else 0.0
        k = tuple(sorted((ka, kb)))
        if fex > best.get(k, (0, 0))[1]:
            best[k] = (j, fex)

    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x

    kept = 0
    for (x, y), (j, fex) in best.items():
        if fex >= a.fex:
            rx, ry = find(x), find(y)
            if rx != ry:
                parent[rx] = ry
            kept += 1

    comp = collections.defaultdict(list)
    for n in list(parent):
        comp[find(n)].append(n)
    clusters = {k: v for k, v in comp.items() if len(v) >= 2}
    out = a.out + '.clusters.tsv'
    with open(out, 'w') as fh:
        fh.write('cluster_id\tsize\tdensity\tfrac_in\tcorroborated\tchrom\tstart\tend\n')
        for i, (_, mem) in enumerate(sorted(clusters.items(), key=lambda kv: -len(kv[1]))):
            for (s, e) in sorted(mem):
                fh.write(f'JX{i}\t{len(mem)}\tNA\tNA\tNA\t{a.chrom}\t{s-1}\t{e}\n')
    sizes = sorted((len(v) for v in clusters.values()), reverse=True)
    print(f'{a.chrom} t={a.t} f_ex={a.fex}: {kept} pairs kept | {len(clusters)} components | '
          f'{sum(sizes)} members | largest {sizes[0] if sizes else 0}')


if __name__ == '__main__':
    main()
