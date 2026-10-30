#!/usr/bin/env python3
"""Multi-copy loci: our assembled GTF vs StringTie / FLAIR GTFs on the mainstay + YAG families (2026-10-31, user "if we focus on
multi-copy loci can we already affirm that we are doing a better job than stringtie and flair?"). Per family member copy and tool:
  modeled      the tool has a spliced (>= 2 exon) transcript whose exons map to this copy (famtruth.map_copy, most shared exonic bp)
  junctions    distinct introns of those transcripts; `annot` = both ends on an exon-block boundary of the copy's annotated exon union
               (the CAT GFF3 is gene-level, so block boundaries stand in for annotated splice sites, both species alike)
  spliced_tx   number of spliced transcripts at the copy
Per family: copies expressed (modeled by ANY tool), per tool modeled copies, junction precision (annot / all), spliced transcripts.
  toolcompare.py --species human|gorilla LABEL=GTF ...     (env as famtruth)"""
import sys, os, re, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import famtruth as F
a = sys.argv[1:]
species = a[a.index('--species') + 1] if '--species' in a else 'human'
if '--species' in a: i = a.index('--species'); del a[i:i + 2]
G = F.load_ggo(F.load_cat()) if species == 'gorilla' else F.load_cat()
idx = F.exon_index(G)
pats = dict(F.GGO_MAINSTAYS) if species == 'gorilla' else {**F.MAINSTAYS, **F.YAGS}
def is_member(d, rx):
    if species == 'gorilla': return bool(rx.search(d['name']) or rx.search(d.get('desc', '')))
    return bool(rx.match(d['name'])) and d['biotype'] == 'protein_coding'
members = {}
for name, pat in pats.items():
    rx = re.compile(pat)
    for g, d in G.items():
        if is_member(d, rx): members[g] = name
bounds = {g: ({e for _, e in G[g]['exons']}, {s for s, _ in G[g]['exons']}) for g in members}  # 0-based half-open blocks
def scan(path):
    """gene -> dict(spliced, junc: Counter of (donor0, acceptor0))"""
    tx = collections.defaultdict(list); meta = {}
    for line in open(path):
        f = line.split('\t')
        if len(f) < 9 or f[2] != 'exon': continue
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        if not m: continue
        t = m.group(1); tx[t].append((int(f[3]) - 1, int(f[4]))); meta[t] = (f[0], f[6] if f[6] in '+-' else '+')
    out = collections.defaultdict(lambda: dict(spliced=0, junc=set()))
    for t, ex in tx.items():
        if len(ex) < 2: continue
        ex.sort()
        g = F.map_copy(idx, meta[t][0], meta[t][1], ex)
        if g not in members: continue
        out[g]['spliced'] += 1
        for (s0, e0), (s1, e1) in zip(ex, ex[1:]): out[g]['junc'].add((e0, s1))
    return out
tools = [(x.split('=', 1)[0], x.split('=', 1)[1]) for x in a]
res = {lab: scan(p) for lab, p in tools}
expressed = {g for r in res.values() for g in r}
fam_members = collections.defaultdict(set)
for g, n in members.items(): fam_members[n].add(g)
print("family\tannotated\texpressed(any tool)\t" + "\t".join(f"{lab}: modeled / junc_prec / spliced_tx" for lab, _ in tools))
tot = collections.defaultdict(lambda: [0, 0, 0, 0])
texp = tann = 0
for name in pats:
    mem = fam_members[name]; exp = mem & expressed
    if not exp: continue
    cells = []
    for lab, _ in tools:
        r = res[lab]; modeled = sum(1 for g in exp if g in r)
        ja = jn = sp = 0
        for g in exp:
            if g not in r: continue
            sp += r[g]['spliced']
            for d, ac in r[g]['junc']:
                jn += 1; ja += (d in bounds[g][0] and ac in bounds[g][1])
        cells.append(f"{modeled} / {ja / jn if jn else 0:.2f} / {sp}")
        t = tot[lab]; t[0] += modeled; t[1] += ja; t[2] += jn; t[3] += sp
    texp += len(exp); tann += len(mem)
    print(f"{name}\t{len(mem)}\t{len(exp)}\t" + "\t".join(cells))
print(f"ALL\t{tann}\t{texp}\t" + "\t".join(f"{tot[lab][0]} / {tot[lab][1] / tot[lab][2] if tot[lab][2] else 0:.2f} / {tot[lab][3]}" for lab, _ in tools))
