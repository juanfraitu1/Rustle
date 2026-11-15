#!/usr/bin/env python3
"""Row 1676: add a `related` column (1 / 0 / NA) to an edgefeat.py table.   edgerelated.py SPECIES IN.tsv OUT.tsv
1 = Compara paralogs at any level or one HGNC gene group (pairscore's `related`); NA = a gene without a human Ensembl id."""
import sys, os, csv, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import famtruth as F


def main(species, src, dst):
    G = F.genes_from_argv(['--species', species], species)
    anyp = set()
    for r in csv.reader(F._open(F.PARALOGS), delimiter='\t'):
        if len(r) < 2 or r[0].startswith('Gene'): continue
        x, y = r[0].split('.')[0], r[1].split('.')[0]
        if x and y: anyp.add((min(x, y), max(x, y)))
    groups = collections.defaultdict(set)
    for r in csv.DictReader(open(os.environ['DENOVO_HGNC']), delimiter='\t'):
        for g in r['gene_group_id'].split('|'):
            if g and r['ensembl_gene_id']: groups[r['ensembl_gene_id']].add(g)
    n = collections.Counter()
    with open(dst, 'w') as w:
        for i, l in enumerate(open(src)):
            f = l.rstrip('\n').split('\t')
            if i == 0:
                w.write(l.rstrip('\n') + '\trelated\trel_kind\n'); continue
            ex, ey = ({e.split('.')[0] for e in (G[g].get('ens') or ())} for g in f[:2])
            cp = any((min(p, q), max(p, q)) in anyp for p in ex for q in ey); hg = any(groups[p] & groups[q] for p in ex for q in ey)
            v = 'NA' if not ex or not ey else str(int(cp or hg))
            kind = 'NA' if v == 'NA' else 'compara' if cp else 'hgnc_only' if hg else 'none'
            n[v] += 1; w.write(l.rstrip('\n') + f'\t{v}\t{kind}\n')
    print(src, dict(n))


if __name__ == '__main__':
    main(*sys.argv[1:])
