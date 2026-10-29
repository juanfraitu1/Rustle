#!/usr/bin/env python3
"""Where the lift pipeline (U1: rep templates, reciprocal, retro-filtered, --extra-edges into MCL) loses each TRUE pair it misses.

  lift_misses.py RUN_PREFIX LIFT_DIR ARM_PREFIX [--levels simiiformes] [--species gorilla]
    RUN_PREFIX.families.gtf / RUN_PREFIX.fam.loci.gff3 / RUN_PREFIX.fam.clusters.tsv (base), LIFT_DIR/{placed,retro,edges}.tsv + keymap.tsv,
    ARM_PREFIX.fam.{copies,clusters}.tsv (the U1 arm). Stages: hit (arm) | edge-split (edge given, MCL kept them apart) | retro |
    one-way | no-placement; plus whether the gated all-vs-all has a record for the pair.
"""
import collections, itertools, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import famtruth as F  # noqa: E402

run, lift, arm = sys.argv[1:4]
levels = sys.argv[sys.argv.index('--levels') + 1] if '--levels' in sys.argv else 'simiiformes'
species = sys.argv[sys.argv.index('--species') + 1] if '--species' in sys.argv else 'gorilla'
G = F.load_cat() if species == 'human' else F.load_ggo(F.load_cat())
if species == 'hg002':
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'guided_eval'))
    import guided3 as g3
    g3.HG002_XREF = os.environ['DENOVO_HG002_XREF']
    F.GGO_GFF = os.environ['DENOVO_HG002_GFF']
    G = g3.load_genes('hg002')
idx = F.exon_index(G)
fams, _, _ = F.load_catalog(f'{arm}.fam.copies.tsv', G, idx)
in_fam = set().union(*fams.values()) if fams else set()
expressed = F.expressed_from_gtf(f'{run}.families.gtf', idx) | in_fam
pairs, universe = F.truth_pairs(G, F.PARALOGS, F.LEVELS[levels])
uexp = universe & expressed
true = sorted({p for p in pairs if p[0] in uexp and p[1] in uexp})
# loci -> gene
ex = collections.defaultdict(list); meta = {}
for line in open(f'{run}.fam.loci.gff3'):
    if line[0] == '#':
        continue
    f = line.rstrip('\n').split('\t')
    if f[2] == 'gene':
        meta[f'{f[0]}:{f[3]}-{f[4]}'] = (f[0], f[6], re.search(r'ID=gene-([^;]+)', f[8]).group(1))
    elif f[2] == 'exon':
        ex[re.search(r'gene=([^;]+)', f[8]).group(1)].append((int(f[3]) - 1, int(f[4])))
loci_of = collections.defaultdict(set)
for k, (c, s, g) in meta.items():
    gg = F.map_copy(idx, c, s, sorted(ex[g]))
    if gg:
        loci_of[gg].add(k)
key = dict(l.rstrip('\n').split('\t')[:2] for l in open(f'{lift}/keymap.tsv'))


def clusters(path):
    fam = {}
    for line in open(path):
        f = line.rstrip('\n').split('\t')
        if f[0] != 'cluster_id':
            fam[f'{f[5]}:{f[6]}-{f[7]}'] = f[0]
    return fam


base_c, arm_c = clusters(f'{run}.fam.clusters.tsv'), clusters(f'{arm}.fam.clusters.tsv')
placed = {(key[a], key[b]) for a, b in (l.rstrip('\n').split('\t')[:2] for l in open(f'{lift}/placed.tsv'))}
retro = {(key[a], key[b]) for a, b in (l.rstrip('\n').split('\t')[:2] for l in open(f'{lift}/retro.tsv'))}
edges = {tuple(sorted(l.rstrip('\n').split('\t')[:2])) for l in open(f'{lift}/edges.u1.tsv')}
recs = set()
for line in open(f'{run}.fam.loci.gated.paf'):
    q, _, _, _, _, t = line.split('\t', 6)[:6]
    if q != t:
        recs.add((min(q, t), max(q, t)))


def same(c, la, lb):
    return any(c.get(a) is not None and c.get(a) == c.get(b) for a in la for b in lb)


cnt = collections.Counter()
print('gene_a\tgene_b\tstage\tbase_hit\tpaf_record\tplaced_ab\tplaced_ba\tretro\tn_loci_a\tn_loci_b')
for ga, gb in true:
    la, lb = loci_of.get(ga, set()), loci_of.get(gb, set())
    ab = any((a, b) in placed for a in la for b in lb); ba = any((b, a) in placed for a in la for b in lb)
    rt = any((a, b) in retro or (b, a) in retro for a in la for b in lb)
    rec = any((min(a, b), max(a, b)) in recs for a in la for b in lb)
    if same(arm_c, la, lb):
        st = 'hit'
    elif any(tuple(sorted((a, b))) in edges for a in la for b in lb):
        st = 'edge-split'
    elif ab and ba:
        st = 'both-ways-no-edge'
    elif ab or ba:
        st = 'one-way' + ('+retro' if rt else '')
    elif rt:
        st = 'retro-only'
    elif not la or not lb:
        st = 'no-locus'
    else:
        st = 'no-placement'
    cnt[st] += 1
    print(f"{G[ga]['name']}\t{G[gb]['name']}\t{st}\t{int(same(base_c, la, lb))}\t{int(rec)}\t{int(ab)}\t{int(ba)}\t{int(rt)}\t{len(la)}\t{len(lb)}")
print('\n# summary:', dict(cnt), 'true pairs', len(true), file=sys.stderr)
