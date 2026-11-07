#!/usr/bin/env python3
"""Row 1571: rejoin locus pieces the family stage gives no sequence reason to keep apart (docs/PREREG_locus_pieces_2026-11-07.md).

Pieces (gene_id GENE.pK, from locus_pieces.py) of one original gene_id that land in the same family stay separate copies only
when their loci share a gated-PAF record of >= 300 bp at >= 0.70 identity (single linkage); each linked group keeps one row,
the piece with the most reads. Pieces in different families are kept as they are.

usage: piece_merge.py ARM_PREFIX OUT_COPIES_TSV     (reads ARM_PREFIX.fam.copies.tsv and ARM_PREFIX.fam.loci.gated.paf)
       piece_merge.py --selftest
"""
import collections, re, sys

MIN_ALEN, MIN_ID = 300, 0.70


def merge(rows, linked):
    """rows: dicts with fam, orig, key, reads, i. linked(a_key, b_key) -> bool. Returns the indices to keep."""
    groups = collections.defaultdict(list)
    for r in rows:
        groups[(r['fam'], r['orig'])].append(r)
    keep = set()
    for g in groups.values():
        parent = list(range(len(g)))
        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        for a in range(len(g)):
            for b in range(a + 1, len(g)):
                if linked(g[a]['key'], g[b]['key']):
                    parent[find(a)] = find(b)
        comps = collections.defaultdict(list)
        for a in range(len(g)):
            comps[find(a)].append(g[a])
        # a linked group is a set of pieces that align to each other = separate copies; unlinked pieces are one gene
        for c in comps.values():
            if len(c) == 1:
                continue
            for r in c:
                keep.add(r['i'])
        lone = [c[0] for c in comps.values() if len(c) == 1]
        if lone:
            keep.add(max(lone, key=lambda r: (r['reads'], -r['i']))['i'])
    return keep


def main(prefix, out):
    with open(prefix + '.fam.copies.tsv') as fh:
        hdr = fh.readline()
        cols = hdr.rstrip('\n').split('\t')
        raw = [l for l in fh]
    rows = []
    for i, l in enumerate(raw):
        r = dict(zip(cols, l.rstrip('\n').split('\t')))
        rows.append(dict(i=i, fam=r['family_id'], orig=re.sub(r'\.p\d+$', '', r['gene_id']), reads=float(r['n_reads']),
                         key=f"{r['chrom']}:{int(r['locus_start']) + 1}-{r['locus_end']}"))
    need = collections.Counter((r['fam'], r['orig']) for r in rows)
    keys = {r['key'] for r in rows if need[(r['fam'], r['orig'])] > 1}
    pairs = set()
    for line in open(prefix + '.fam.loci.gated.paf'):
        f = line.split('\t', 12)
        if f[0] in keys and f[5] in keys and f[0] != f[5] and int(f[10]) >= MIN_ALEN and int(f[9]) / int(f[10]) >= MIN_ID:
            pairs.add(frozenset((f[0], f[5])))
    keep = merge(rows, lambda a, b: frozenset((a, b)) in pairs)
    with open(out, 'w') as oh:
        oh.write(hdr)
        oh.writelines(raw[i] for i in sorted(keep))
    split = sum(1 for k, n in need.items() if n > 1)
    print(f"[piece_merge] copies {len(rows)} -> {len(keep)}; family x locus groups with >= 2 pieces {split}", file=sys.stderr)


def selftest():
    rows = [dict(i=0, fam='F', orig='G', key='a', reads=5), dict(i=1, fam='F', orig='G', key='b', reads=9),
            dict(i=2, fam='F', orig='H', key='c', reads=3), dict(i=3, fam='F', orig='H', key='d', reads=4),
            dict(i=4, fam='E', orig='G', key='e', reads=1)]
    linked = lambda x, y: {x, y} == {'c', 'd'}
    # G's two pieces in F do not align: one gene, keep its most-read piece (1); H's pieces align: tandem copies, keep both;
    # G's piece in another family E stays
    assert merge(rows, linked) == {1, 2, 3, 4}
    print('selftest ok')


if __name__ == '__main__':
    selftest() if sys.argv[1:] == ['--selftest'] else main(*sys.argv[1:3])
