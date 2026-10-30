#!/usr/bin/env python3
"""Mainstay / YAG family score (2026-10-29, user "focus on specific gene families ... ensure that all mainstays are correctly and
specifically recalled"). Per named family (famtruth MAINSTAYS + YAGS), human only:
  expressed    protein-coding members with a de novo locus
  recalled     members whose BEST locus (most exonic bp shared with the gene) is in the largest family holding the group
  recalled_annot  the same over ALL annotated members (semi-guided copies of unexpressed members count; 2026-10-30)
  families     number of families holding any member locus
  unrel_foreign  protein-coding genes in the largest family that are not members and not related to any member (Compara at any
                 level, HGNC group, same symbol); rel_foreign = the related ones (pseudogene-like siblings, unnamed copies)
  mainstay_score.py BASE.families.gtf COPIES_TSV...     (env as famtruth + DENOVO_HGNC)"""
import sys, os, re, csv, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import famtruth as F
a = sys.argv[1:]
species = a[a.index('--species') + 1] if '--species' in a else 'human'
if '--species' in a: i = a.index('--species'); del a[i:i + 2]
if species == 'hg002':  # Ensembl rapid-release genes via guided3; HGNC symbols from the xref so the family patterns match
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'guided_eval'))
    import guided3 as g3
    g3.HG002_XREF = os.environ['DENOVO_HG002_XREF']; F.GGO_GFF = os.environ['DENOVO_HG002_GFF']
    G = g3.load_genes('hg002')
    sym = {}
    for r in csv.DictReader(open(os.environ['DENOVO_HG002_XREF']), delimiter='\t'):
        if r['db_name'] == 'HGNC Symbol' and r['gene_stable_id']: sym.setdefault(r['gene_stable_id'], r['xref_label'])
    for g, d in G.items(): d['name'] = sym.get(g.split(':')[-1], d['name'])
elif species == 'gorilla':  # Yoo families: symbol or RefSeq description, any biotype (F.GGO_MAINSTAYS)
    G = F.load_ggo(F.load_cat())
else:
    G = F.load_cat()
idx = F.exon_index(G)
expressed = F.expressed_from_gtf(a[0], idx)
anyp = set()
for r in csv.DictReader(open(F.PARALOGS), delimiter='\t'):
    x, y = r['Gene stable ID'].split('.')[0], r['Human paralogue gene stable ID'].split('.')[0]
    if x and y: anyp.add((min(x, y), max(x, y)))
groups = collections.defaultdict(set)
if os.environ.get('DENOVO_HGNC'):
    for r in csv.DictReader(open(os.environ['DENOVO_HGNC']), delimiter='\t'):
        for g in r['gene_group_id'].split('|'):
            if g and r['ensembl_gene_id']: groups[r['ensembl_gene_id']].add(g)
def related(x, y):
    if G[x]['name'] == G[y]['name']: return True
    if species == 'gorilla' and G[x].get('desc') and G[x]['desc'].replace('-like', '') == G[y].get('desc', '').replace('-like', ''): return True
    ex, ey = {e.split('.')[0] for e in G[x]['ens']}, {e.split('.')[0] for e in G[y]['ens']}
    return any((min(p, q), max(p, q)) in anyp for p in ex for q in ey) or any(groups[p] & groups[q] for p in ex for q in ey)
def share(g, chrom, strand, ex):
    if G[g]['chrom'] != chrom or G[g]['strand'] != strand: return 0
    tot = 0
    for a, b in ex:
        for s, e in G[g]['exons']: tot += max(0, min(b, e) - max(a, s))
    return tot
pats = dict(F.MAINSTAYS); pats.update(F.YAGS)
if species == 'gorilla': pats = dict(F.GGO_MAINSTAYS)
def is_member(d, rx):  # gorilla: symbol or description, any biotype; human: symbol, protein-coding
    if species == 'gorilla': return bool(rx.search(d['name']) or rx.search(d.get('desc', '')))
    return bool(rx.match(d['name'])) and d['biotype'] == 'protein_coding'
print("run\tfamily\texpressed\trecalled\tfamilies\trel_foreign\tunrel_foreign\tunrel_names\tannotated\trecalled_annot")
for path in a[1:]:
    rows = [r for r in csv.DictReader(open(path), delimiter='\t')]
    loc_gene = {}; loc_fam = {}; loc_ex = {}
    for r in rows:
        ex = F.blocks(r['exons']); g = F.map_copy(idx, r['chrom'], r['strand'], ex)
        loc_fam[r['tid']] = r['family_id']; loc_ex[r['tid']] = (r['chrom'], r['strand'], ex)
        if g: loc_gene[r['tid']] = g
    gene_locs = collections.defaultdict(list)
    for t, g in loc_gene.items(): gene_locs[g].append(t)
    fam_genes = collections.defaultdict(set)
    for t, g in loc_gene.items(): fam_genes[loc_fam[t]].add(g)
    tot = collections.Counter()
    for name, pat in pats.items():
        rx = re.compile(pat)
        mem = {g for g, d in G.items() if is_member(d, rx)}
        exp = mem & expressed
        if not exp and species != 'gorilla': continue  # gorilla: Yoo families with no expressed copy are reported as 0 of N
        # read-through units (gorilla JMJD7-PLA2G4B, 2026-10-30): adjacent same-strand members within 20 kb are ONE unit; a unit
        # is expressed / recalled when any of its genes is, so a read-through locus scored against one gene counts the unit once
        unit_of = {g: g for g in mem}
        if species == 'gorilla' and '-' in name:
            ordered = sorted(mem, key=lambda g: (G[g]['chrom'], G[g]['strand'], G[g]['exons'][0][0]))
            for p, q_ in zip(ordered, ordered[1:]):
                if G[p]['chrom'] == G[q_]['chrom'] and G[p]['strand'] == G[q_]['strand'] and G[q_]['exons'][0][0] - G[p]['exons'][-1][1] <= 20_000:
                    unit_of[q_] = unit_of[p]
        best = {g: max(gene_locs[g], key=lambda t: share(g, *loc_ex[t])) for g in exp if gene_locs[g]}
        best_all = {g: max(gene_locs[g], key=lambda t: share(g, *loc_ex[t])) for g in mem if gene_locs[g]}  # semi-guided (2026-10-30): any annotated member with a copy, expressed or silent
        held = collections.Counter(loc_fam[t] for g in exp for t in gene_locs[g])
        top_by_best = collections.Counter(loc_fam[t] for t in best.values())  # largest = family holding most members' BEST loci
        if not held:
            print(f"{os.path.basename(path)}\t{name}\t{len(exp)}\t0\t0\t0\t0\t\t{len(mem)}\t0"); continue
        top = top_by_best.most_common(1)[0][0]
        rec = len({unit_of[g] for g, t in best.items() if loc_fam[t] == top})
        rec_all = len({unit_of[g] for g, t in best_all.items() if loc_fam[t] == top})
        n_exp = len({unit_of[g] for g in exp}); n_mem = len(set(unit_of.values()))
        foreign = {g for g in fam_genes[top] if g not in mem and G[g]['biotype'] == 'protein_coding'}
        unrel = sorted(G[g]['name'] for g in foreign if not any(related(g, m) for m in exp))
        print(f"{os.path.basename(path)}\t{name}\t{n_exp}\t{rec}\t{len(held)}\t{len(foreign) - len(unrel)}\t{len(unrel)}\t{','.join(unrel[:6])}\t{n_mem}\t{rec_all}")
        tot['exp'] += n_exp; tot['rec'] += rec; tot['unrel'] += len(unrel); tot['mem'] += n_mem; tot['rec_all'] += rec_all
    print(f"{os.path.basename(path)}\tALL\t{tot['exp']}\t{tot['rec']}\t\t\t{tot['unrel']}\t\t{tot['mem']}\t{tot['rec_all']}")
