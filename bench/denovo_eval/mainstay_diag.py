#!/usr/bin/env python3
"""Per named family (mainstays + YAGs), on one sketchfam arm: which family holds each expressed member, the foreign protein-coding
genes in those families and the edges that bring them in, and the best containment between the member sub-groups (why they split).
  mainstay_diag.py BASE.families.gtf ARM_PREFIX [NAME ...]    (ARM_PREFIX.fam.copies.tsv + ARM_PREFIX.edges.tsv; env as famtruth)"""
import sys, os, re, collections, csv
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import famtruth as F
gtf, arm = sys.argv[1:3]; only = set(sys.argv[3:])
G = F.load_cat(); idx = F.exon_index(G)
expressed = F.expressed_from_gtf(gtf, idx)
pats = dict(F.MAINSTAYS); pats.update(F.YAGS)
loc_gene = {}; loc_fam = {}
for r in csv.DictReader(open(arm + '.fam.copies.tsv'), delimiter='\t'):
    g = F.map_copy(idx, r['chrom'], r['strand'], F.blocks(r['exons']))
    loc_fam[r['tid']] = r['family_id']
    if g: loc_gene[r['tid']] = g
edges = collections.defaultdict(float)
for l in open(arm + '.edges.tsv'):
    a, b, n, c, m = l.split('\t')[:5]
    edges[(a, b)] = max(edges[(a, b)], float(c)); edges[(b, a)] = edges[(a, b)]
gene_locs = collections.defaultdict(set)
for t, g in loc_gene.items(): gene_locs[g].add(t)
fam_members = collections.defaultdict(set)
for t, f in loc_fam.items(): fam_members[f].add(t)
for name, pat in pats.items():
    if only and name not in only: continue
    rx = re.compile(pat)
    mem = {g for g, d in G.items() if rx.match(d['name']) and d['biotype'] == 'protein_coding' and g in expressed}
    if not mem: continue
    by_fam = collections.defaultdict(set)
    for g in mem:
        fs = {loc_fam[t] for t in gene_locs[g]} or {'-'}
        for f in fs: by_fam[f].add(G[g]['name'])
    print(f"\n## {name}: {len(mem)} expressed coding members in {len([f for f in by_fam if f != '-'])} families" + (f", unplaced {sorted(by_fam['-'])}" if '-' in by_fam else ''))
    for f, names in sorted(by_fam.items(), key=lambda x: -len(x[1])):
        if f == '-': continue
        foreign = collections.Counter()
        for t in fam_members[f]:
            g = loc_gene.get(t)
            if g and g not in mem and G[g]['biotype'] == 'protein_coding': foreign[G[g]['name']] += 1
        line = f"  {f}: {len(names)} members {sorted(names)}"
        if foreign:
            # strongest edge from each foreign gene into a member locus
            best = {}
            for t in fam_members[f]:
                g = loc_gene.get(t)
                if g and G[g]['name'] in foreign:
                    for u in fam_members[f]:
                        if loc_gene.get(u) in mem and (t, u) in edges:
                            best[G[g]['name']] = max(best.get(G[g]['name'], 0), edges[(t, u)])
            line += f"\n      foreign: " + ', '.join(f"{n} (edge to member {best.get(n, 0):.2f})" for n in sorted(foreign))
        print(line)
    fams = [f for f in by_fam if f != '-']
    if len(fams) > 1:
        for i in range(len(fams)):
            for j in range(i + 1, len(fams)):
                b = max((edges.get((t, u), 0) for t in fam_members[fams[i]] for u in fam_members[fams[j]]), default=0)
                print(f"  best edge {fams[i]} <-> {fams[j]}: {b:.2f}")
