#!/usr/bin/env python3
"""Pair precision / recall only, no F, no family matching (2026-10-29, user: "we don't need bipartite matching anymore, the only thing
that is useful is precision and recall", "we dont need F anymore, just precision and recall", "pseudogenes don't need their own category"). Truth = Compara paralog pairs at --level among protein-coding genes with a de
novo transcript in the run (famtruth's expressed universe). Predicted = gene pairs inside one family after mapping copies to genes
(famtruth.map_copy). A predicted pair where either gene is NOT protein-coding (pseudogene, lncRNA) is IGNORED: neither hit nor false.
Also reports 'related' (Compara at any level / same HGNC group / same symbol) among the non-true predicted pairs, so a false pair means
an unrelated protein-coding gene pair.

  pairscore.py --species human|gorilla|hg002 --gtf RUN.families.gtf [--level simiiformes] COPIES_TSV...
  env as famtruth (+ DENOVO_HG002_GFF / DENOVO_HG002_XREF for hg002, DENOVO_HGNC)
"""
import sys, os, csv, collections, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import famtruth as F
a = sys.argv[1:]
species = a[a.index('--species') + 1]; gtf = a[a.index('--gtf') + 1]
level = a[a.index('--level') + 1] if '--level' in a else 'simiiformes'
paths = [x for i, x in enumerate(a) if not x.startswith('--') and a[i - 1] not in ('--species', '--gtf', '--level')]
if species == 'hg002':
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'guided_eval'))
    import guided3 as g3
    g3.HG002_XREF = os.environ['DENOVO_HG002_XREF']; F.GGO_GFF = os.environ['DENOVO_HG002_GFF']
    G = g3.load_genes('hg002')
else:
    G = F.load_cat() if species == 'human' else F.load_ggo(F.load_cat())
idx = F.exon_index(G)
expressed = F.expressed_from_gtf(gtf, idx)
truth, universe = F.truth_pairs(G, F.PARALOGS, F.LEVELS[level])
truth = {tuple(sorted(p)) for p in truth if p[0] in expressed and p[1] in expressed}
anyp = set()
for r in csv.DictReader(open(F.PARALOGS), delimiter='\t'):
    x, y = r['Gene stable ID'].split('.')[0], r['Human paralogue gene stable ID'].split('.')[0]
    if x and y: anyp.add((min(x, y), max(x, y)))
groups = collections.defaultdict(set)
hg = os.environ.get('DENOVO_HGNC')
if hg:
    for r in csv.DictReader(open(hg), delimiter='\t'):
        for g in r['gene_group_id'].split('|'):
            if g and r['ensembl_gene_id']: groups[r['ensembl_gene_id']].add(g)
def related(x, y):
    if G[x]['name'] == G[y]['name']: return True
    ex, ey = {e.split('.')[0] for e in G[x]['ens']}, {e.split('.')[0] for e in G[y]['ens']}
    return any((min(p, q), max(p, q)) in anyp for p in ex for q in ey) or any(groups[p] & groups[q] for p in ex for q in ey)
print(f"# level {level}: truth pairs among expressed protein-coding genes {len(truth)} (universe {len(universe & expressed)})")
print("run\tpred_pairs\tignored_noncoding\thit\tfalse\tof_which_related\tprecision\trecall")
for p in paths:
    fams, _, _ = F.load_catalog(p, G, idx)
    pred = set()
    for v in fams.values():
        for x, y in itertools.combinations(sorted(v), 2): pred.add((x, y))
    coding = {(x, y) for x, y in pred if G[x]['biotype'] == 'protein_coding' and G[y]['biotype'] == 'protein_coding' and x in universe and y in universe}
    hit = len(coding & truth); false = [(x, y) for x, y in coding if (x, y) not in truth]
    rel = sum(1 for x, y in false if related(x, y))
    prec = hit / len(coding) if coding else 0.0; rec = hit / len(truth) if truth else 0.0
    print(f"{p}\t{len(coding)}\t{len(pred) - len(coding)}\t{hit}\t{len(false)}\t{rel}\t{prec:.3f}\t{rec:.3f}")
