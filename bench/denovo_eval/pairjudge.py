"""LIFTMODE-GW3 asymmetric pair judge (docs/LEDGER_O1.md): recall at Compara Simiiformes (famtruth), precision against Compara at ANY level
+ HGNC gene groups + same symbol (PAR).

  pairjudge.py COPIES_TSV [COPIES_TSV...]   one line per catalog: true (Simiiformes) / older-or-HGNC / wrong pairs among judged genes
"""
import collections
import csv
import itertools
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import famtruth as F  # noqa: E402

HGNC = os.environ.get('DENOVO_HGNC', '/mnt/linuxdisk/tmp/hgnc/hgnc_complete_set.txt')


def main(paths):
    species = os.environ.get('PAIRJUDGE_SPECIES', 'gorilla')  # EDGE-HUMAN: 'human' judges CAT genes directly
    G = F.load_cat() if species == 'human' else F.load_ggo(F.load_cat())
    idx = F.exon_index(G)
    sim, universe = F.truth_pairs(G, F.PARALOGS, F.LEVELS['simiiformes'])
    sim = {tuple(sorted(p)) for p in sim}
    anyp = set()
    for r in csv.DictReader(open(F.PARALOGS), delimiter='\t'):
        a, b = r['Gene stable ID'].split('.')[0], r['Human paralogue gene stable ID'].split('.')[0]
        if a and b:
            anyp.add((min(a, b), max(a, b)))
    groups = collections.defaultdict(set)
    for r in csv.DictReader(open(HGNC), delimiter='\t'):
        for g in r['gene_group_id'].split('|'):
            if g and r['ensembl_gene_id']:
                groups[r['ensembl_gene_id']].add(g)

    def related(a, b):
        if G[a]['name'] == G[b]['name']:
            return 'same-symbol'
        ea, eb = {e.split('.')[0] for e in G[a]['ens']}, {e.split('.')[0] for e in G[b]['ens']}
        if any((min(x, y), max(x, y)) in anyp for x in ea for y in eb):
            return 'older'
        if any(groups[x] & groups[y] for x in ea for y in eb):
            return 'hgnc'
        return None

    for p in paths:
        fams, _, _ = F.load_catalog(p, G, idx)
        pred = {(a, b) for v in fams.values() for a, b in itertools.combinations(sorted(v), 2) if a in universe and b in universe}
        c = collections.Counter()
        wrong = []
        for a, b in pred:
            if (a, b) in sim:
                c['true'] += 1
            else:
                k = related(a, b)
                c[k or 'wrong'] += 1
                if not k:
                    wrong.append((G[a]['name'], G[b]['name']))
        print(f"{p}: pairs {len(pred)} | true(Simiiformes) {c['true']}/{len(sim)} "
              f"| older {c['older']} hgnc {c['hgnc']} same-symbol {c['same-symbol']} | WRONG {c['wrong']} {sorted(wrong)[:12]}")


if __name__ == '__main__':
    main(sys.argv[1:])
