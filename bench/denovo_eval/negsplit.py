#!/usr/bin/env python3
"""PREREG r1621 (docs/PREREG_families_negative_edges_2026-11-10.md): negative edges inside MCL families.

  negsplit.py step0 COPIES_TSV G0_TSV SHARE_TSV --gtf RUN.families.gtf [--species ...]   (env as pairscore)
  negsplit.py split COPIES_TSV G0_TSV SHARE_TSV LAMBDA OUT_COPIES_TSV
"""
import sys, os, collections, itertools
mode, copies, g0, share_p = sys.argv[1:5]
share = {}
for l in open(share_p):
    a, b, s = l.split('\t'); share[(a, b)] = float(s)
W = {}  # (a, b) sorted -> signed base weight: + when SHARE >= 0.5, - otherwise
for l in open(g0):
    a, b, w = l.rstrip('\n').split('\t'); k = (min(a, b), max(a, b))
    W[k] = float(w) if share.get(k, 0.0) >= 0.5 else -float(w)
rows = [l.rstrip('\n').split('\t') for l in open(copies)]; h = {k: i for i, k in enumerate(rows[0])}
node = lambda f: f"{f[h['chrom']]}:{int(f[h['locus_start']]) + 1}-{f[h['locus_end']]}"
if mode == 'step0':
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import famtruth as F
    a = sys.argv[5:]; species = a[a.index('--species') + 1] if '--species' in a else 'human'
    G = F.genes_from_argv(a, species); idx = F.exon_index(G)
    expressed = F.expressed_from_gtf(a[a.index('--gtf') + 1], idx)
    truth, universe = F.truth_pairs(G, F.need(F.PARALOGS, 'DENOVO_PARALOGS'), F.levels('simiiformes'))
    truth = {tuple(sorted(p)) for p in truth if p[0] in expressed and p[1] in expressed}
    fam = collections.defaultdict(lambda: collections.defaultdict(set))  # family -> gene -> nodes
    for f in rows[1:]:
        g = F.map_copy(idx, f[h['chrom']], f[h['strand']] if f[h['strand']] in '+-' else '+', F.blocks(f[h['exons']]))
        if g: fam[f[h['family_id']]][g].add(node(f))
    n = collections.Counter()
    for genes in fam.values():
        for x, y in itertools.combinations(sorted(genes), 2):
            if not (G[x]['biotype'] == G[y]['biotype'] == 'protein_coding' and x in universe and y in universe) or G[x]['name'] == G[y]['name']: continue
            ws = [W.get((min(p, q), max(p, q))) for p in genes[x] for q in genes[y]]
            kind = 'positive' if any(w is not None and w > 0 for w in ws) else 'negative' if any(w is not None for w in ws) else 'no edge'
            n[('hit' if (x, y) in truth else 'false', kind)] += 1
    for c in ('hit', 'false'):
        t = sum(v for (a_, _), v in n.items() if a_ == c)
        print(c, t, *(f"{k} {n[(c, k)]} ({n[(c, k)] / max(1, t):.3f})" for k in ('positive', 'negative', 'no edge')), sep='\t')
else:
    lam = float(sys.argv[5]); out = sys.argv[6]
    fam = collections.defaultdict(list)
    for i, f in enumerate(rows[1:], 1): fam[f[h['family_id']]].append(i)
    nsplit = 0
    for fid, idxs in fam.items():
        nodes = sorted({node(rows[i]) for i in idxs}); cl = {x: [x] for x in nodes}  # cluster id -> members
        net = collections.defaultdict(float)
        for x, y in itertools.combinations(nodes, 2):
            w = W.get((x, y))
            if w is not None: net[(x, y)] = w if w > 0 else lam * w
        while True:
            best = max(net.items(), key=lambda kv: (kv[1], kv[0]), default=None)
            if best is None or best[1] <= 0: break
            (x, y), _ = best; cl[x] += cl.pop(y)
            for z in list(cl):
                if z == x: continue
                v = net.pop((min(y, z), max(y, z)), 0.0)
                if v: net[(min(x, z), max(x, z))] = net.get((min(x, z), max(x, z)), 0.0) + v
            net.pop((x, y), None)
        parts = sorted(cl.values(), key=lambda m: (-len(m), m)); lab = {x: k for k, m in enumerate(parts) for x in m}
        nsplit += len(parts) > 1
        for i in idxs: rows[i][h['family_id']] = fid if lab[node(rows[i])] == 0 else f"{fid}.{lab[node(rows[i])]}"
    with open(out, 'w') as w:
        for f in rows: w.write('\t'.join(f) + '\n')
    print(f'lambda {lam}: {nsplit} of {len(fam)} families split', file=sys.stderr)
