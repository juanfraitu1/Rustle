#!/usr/bin/env python3
"""Per-exon shared-exon rule (2026-10-29, diagnosis after the column-exact arm): judge every locus pair of the gated all-vs-all
(table from exon_rules.py: best_col_frac = column rule, exon_frac = share of the smaller locus's exons with >= T of their bases joined
exon-to-exon, union over records at identity >= 0.70) against the truth: Compara at the registered levels (recent / simiiformes),
'older' (Compara at any level), 'hgnc' (shared HGNC group), 'same-symbol', else WRONG; LOC / unjudgeable pairs reported as 'unjudged'.
Rows: pairs that pass the column rule at 0.60 (the shipped set) and pairs the exon rule would ADD at each threshold (exon_frac >= X but
column < 0.60), by verdict.  Env: DENOVO_* as famtruth, DENOVO_HGNC.

  exon_rule_judge.py EXON_RULES_TSV LOCI_GFF3 [--species gorilla]
"""
import collections, csv, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import famtruth as F  # noqa: E402



class Judge:
    """Locus pair -> verdict: true-recent / true-simiiformes / older / hgnc / same-symbol / WRONG / same-gene / unjudged."""

    def __init__(self, species, loci_gff3):
        G = F.load_cat() if species == 'human' else F.load_ggo(F.load_cat())
        if species == 'hg002':
            sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'guided_eval'))
            import guided3 as g3
            g3.HG002_XREF = os.environ['DENOVO_HG002_XREF']
            F.GGO_GFF = os.environ['DENOVO_HG002_GFF']
            G = g3.load_genes('hg002')
        idx = F.exon_index(G)
        levels = {lv: {tuple(sorted(p)) for p in F.truth_pairs(G, F.PARALOGS, F.LEVELS[lv])[0]} for lv in ('recent', 'simiiformes')}
        _, universe = F.truth_pairs(G, F.PARALOGS, F.LEVELS['recent'])
        anyp = set()
        for r in csv.DictReader(open(F.PARALOGS), delimiter='\t'):
            a, b = r['Gene stable ID'].split('.')[0], r['Human paralogue gene stable ID'].split('.')[0]
            if a and b:
                anyp.add((min(a, b), max(a, b)))
        groups = collections.defaultdict(set)
        hg = os.environ.get('DENOVO_HGNC')
        if hg:
            for r in csv.DictReader(open(hg), delimiter='\t'):
                for g in r['gene_group_id'].split('|'):
                    if g and r['ensembl_gene_id']:
                        groups[r['ensembl_gene_id']].add(g)
        # locus -> gene via the exons of the loci GFF3
        ex = collections.defaultdict(list); meta = {}
        for line in open(loci_gff3):
            if line[0] == '#':
                continue
            f = line.rstrip('\n').split('\t')
            if f[2] == 'gene':
                meta[f'{f[0]}:{f[3]}-{f[4]}'] = (f[0], f[6], re.search(r'ID=gene-([^;]+)', f[8]).group(1))
            elif f[2] == 'exon':
                ex[re.search(r'gene=([^;]+)', f[8]).group(1)].append((int(f[3]) - 1, int(f[4])))
        gene_of = {k: F.map_copy(idx, c, s, sorted(ex[g])) for k, (c, s, g) in meta.items()}


        def _v(a, b):
            ga, gb = self.gene_of.get(a), self.gene_of.get(b)
            if ga is None or gb is None or ga == gb:
                return 'same-gene' if ga is not None and ga == gb else 'unjudged'
            if ga not in universe or gb not in universe:
                return 'unjudged'
            p = (min(ga, gb), max(ga, gb))
            for lv in ('recent', 'simiiformes'):
                if p in levels[lv]:
                    return 'true-' + lv
            if G[ga]['name'] == G[gb]['name']:
                return 'same-symbol'
            ea, eb = {e.split('.')[0] for e in G[ga]['ens']}, {e.split('.')[0] for e in G[gb]['ens']}
            if any((min(x, y), max(x, y)) in anyp for x in ea for y in eb):
                return 'older'
            if any(groups[x] & groups[y] for x in ea for y in eb):
                return 'hgnc'
            return 'WRONG'




        self.G, self.gene_of, self._v = G, gene_of, _v

    def verdict(self, a, b):
        return self._v(a, b)

    def name(self, a):
        g = self.gene_of.get(a)
        return self.G[g]['name'] if g else '-'


if __name__ == '__main__':
    J = Judge(sys.argv[3] if len(sys.argv) > 3 else 'gorilla', sys.argv[2])
    G, gene_of, verdict = J.G, J.gene_of, J.verdict
    rows = list(csv.DictReader(open(sys.argv[1]), delimiter='\t'))
    for r in rows:
        r['v'] = verdict(r['a'], r['b'])
        r['names'] = (G[gene_of[r['a']]]['name'] if gene_of.get(r['a']) else '-') + '/' + (G[gene_of[r['b']]]['name'] if gene_of.get(r['b']) else '-')
    order = ['true-recent', 'true-simiiformes', 'older', 'hgnc', 'same-symbol', 'WRONG', 'same-gene', 'unjudged']


    def table(name, sel):
        c = collections.Counter(r['v'] for r in sel)
        print(f'{name}\t{len(sel)}\t' + '\t'.join(str(c[k]) for k in order))


    print('set\tn\t' + '\t'.join(order))
    table('column>=0.60 (shipped column arm)', [r for r in rows if float(r['best_col_frac']) >= 0.60])
    for X in (0.5, 0.6, 0.7, 0.8, 0.9, 1.0):
        table(f'exon_frac>={X} added (column<0.60)', [r for r in rows if float(r['best_col_frac']) < 0.60 and float(r['exon_frac']) >= X])
    for X in (0.5, 0.6, 0.7, 0.8):
        table(f'union_col>={X} added (column<0.60)', [r for r in rows if float(r['best_col_frac']) < 0.60 and float(r['union_col_frac']) >= X])
    print('\n# true pairs (either level) with column < 0.60, by exon_frac:')
    for r in sorted((r for r in rows if r['v'].startswith('true') and float(r['best_col_frac']) < 0.60), key=lambda r: -float(r['exon_frac'])):
        print(f"{r['names']}\t{r['v']}\tcol {r['best_col_frac']}\tunion {r['union_col_frac']}\texon {r['exon_frac']} ({r['exons_joined']}/{r['n_exons_smaller']})\tid {r['max_id']}")
    print('\n# WRONG pairs the exon rule adds at 0.5 (column < 0.60):')
    for r in sorted((r for r in rows if r['v'] == 'WRONG' and float(r['best_col_frac']) < 0.60 and float(r['exon_frac']) >= 0.5), key=lambda r: -float(r['exon_frac'])):
        print(f"{r['names']}\tcol {r['best_col_frac']}\texon {r['exon_frac']} ({r['exons_joined']}/{r['n_exons_smaller']})\tid {r['max_id']}\t{r['a']} {r['b']}")

