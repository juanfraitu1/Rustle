#!/usr/bin/env python3
"""Keep a locus only when the reads still need it under the EM (register 1905; O2 extension, joint provenance test).

    em_filter.py ASSEMBLY.gtf BAM AS_TABLE OUT_PREFIX [--min-em 2] [--label] [--pooled] [--threads N]

Loci, read-to-locus rule and ties as tools/eq_classes.py (minimap2 ms, exact band). Reads tied between loci are shared by
its EM; by default every class has its own no-locus member (row 1877), --pooled = the tool's single `outside`.
Writes OUT_PREFIX.em.tsv (locus, placed, max, em, kept) and OUT_PREFIX.gtf = the assembly without the transcripts of
loci whose EM count is under --min-em (loci no read reaches are dropped too), the others tagged reads_em.
--label drops nothing (user 2026-11-22, the advisor's point: a family member that is found stays available): every
transcript is kept and tagged o2_status = "expressed" (reads placed here alone, EM >= --min-em), "shared" (no read placed
here alone but the EM still leaves >= --min-em: indistinguishable from its partners on these reads) or "homology_only"
(reads reach it, all better explained by other loci: the member exists, these reads do not show it expressed)."""
import sys, os, re, argparse, collections, multiprocessing
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'tools'))
import eq_classes as E

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('gtf'); ap.add_argument('bam'); ap.add_argument('tab'); ap.add_argument('out')
    ap.add_argument('--min-em', type=float, default=2); ap.add_argument('--pooled', action='store_true'); ap.add_argument('--label', action='store_true'); ap.add_argument('--threads', type=int, default=8)
    a = ap.parse_args()
    import pysam
    meta, E.G_GEX = E.load_loci(a.gtf)
    E.G_FLOOR = E.floor_fn(10); E.G_BAM, E.G_TAG, E.G_TIES = a.bam, 'ms', E.load_ties(a.tab, E.G_FLOOR, 'ms')
    with multiprocessing.get_context('fork').Pool(a.threads) as pool:
        res = pool.map(E.scan, list(pysam.AlignmentFile(a.bam).references), chunksize=1)
    placed, sets = collections.Counter(), collections.defaultdict(set)
    for pl, pairs, _ in res:
        placed.update(pl)
        for q, g, exact, *x in pairs:
            if exact: sets[q].add(g)
    placed.pop(E.OUTSIDE, None); classes = collections.Counter()
    for s in sets.values():
        if len(s) == 1: placed.update(s - {E.OUTSIDE})
        else:
            c = frozenset(s)
            if not a.pooled and E.OUTSIDE in c: c = frozenset((c - {E.OUTSIDE}) | {'outside#' + ','.join(sorted(c - {E.OUTSIDE}))})
            classes[c] += 1
    theta = E.em(classes, placed); mx = collections.Counter(placed)
    for c, n in classes.items():
        for l in c: mx[l] += n
    em = {l: theta.get(l, float(placed.get(l, 0))) for l in meta}
    keep = {l for l in meta if em[l] >= a.min_em}
    status = {l: ('expressed' if placed.get(l, 0) >= 1 else 'shared') if l in keep else 'homology_only' for l in meta}
    with open(a.out + '.em.tsv', 'w') as w:
        w.write('locus\tplaced\tmax\tem\tkept\tstatus\n')
        for l in sorted(meta, key=lambda l: (meta[l][0], meta[l][1], l)): w.write(f'{l}\t{placed.get(l, 0)}\t{mx[l]}\t{em[l]:.2f}\t{int(l in keep)}\t{status[l]}\n')
    n = 0
    with open(a.out + '.gtf.tmp', 'w') as w:
        for l in open(a.gtf):
            f = l.rstrip('\n').split('\t'); g = re.search(r'gene_id "([^"]+)"', f[8]) if len(f) > 8 else None
            if g and g.group(1) in meta and g.group(1) not in keep and not a.label: n += f[2] == 'transcript'; continue
            if g and g.group(1) in em and f[2] == 'transcript': l = l.rstrip('\n').rstrip(';') + f'; reads_em "{em[g.group(1)]:.1f}"; o2_status "{status[g.group(1)]}";\n'
            w.write(l)
    os.replace(a.out + '.gtf.tmp', a.out + '.gtf')
    c = collections.Counter(status.values())
    print(f"em_filter: {len(meta)} loci: expressed {c['expressed']}, shared {c['shared']}, homology_only {c['homology_only']} (under {a.min_em:g} reads by EM; {'kept and labelled' if a.label else str(n) + ' transcripts dropped'}) -> {a.out}.gtf")

if __name__ == '__main__':
    main()
