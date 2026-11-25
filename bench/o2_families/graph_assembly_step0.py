#!/usr/bin/env python3
"""Step 0 (register 1944, docs/STEP0_graph_assembly_2026-11-25.md): the one decision an assembly in the family graph's
coordinates makes differently: keep a one-read chain when another member has the same chain of exon columns.

    graph_assembly_step0.py BASE.gtf FLOOR1.gtf MULTICOPY.gtf U.fam.clusters.tsv EXONMAP.columns.tsv REFSEQ.gff OUT.tsv"""
import sys, re, bisect, collections
base, f1, mc, clusters, columns, gff, out = sys.argv[1:8]
TID = re.compile(r'transcript_id "([^"]+)"')

def transcripts(gtf, want=None):
    """{tid: (contig, strand, [(s, e)] 0-based half-open sorted, attributes of the transcript line)}"""
    ex, meta = collections.defaultdict(list), {}
    for l in open(gtf):
        f = l.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] not in ('exon', 'transcript'): continue
        t = TID.search(f[8]).group(1)
        if f[2] == 'exon': ex[t].append((int(f[3]) - 1, int(f[4]))); meta.setdefault(t, [f[0], f[6], ''])
        else: meta.setdefault(t, [f[0], f[6], ''])[2] = f[8]
    return {t: (meta[t][0], meta[t][1], sorted(e), meta[t][2]) for t, e in ex.items()}

def chain(c, s, ex): return (c, s, tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1)))

# members and columns
mem = collections.defaultdict(list); fam_of = {}; size = collections.Counter()
L = [l.rstrip('\n').split('\t') for l in open(clusters)]; h = {c: i for i, c in enumerate(L[0])}
for f in L[1:]:
    k = (f[h['chrom']], int(f[h['start']]), int(f[h['end']])); fam_of[k] = f[h['cluster_id']]; size[f[h['cluster_id']]] += 1
for (c, s, e), fam in fam_of.items():
    if size[fam] >= 2: mem[c].append((s, e, fam))
for v in mem.values(): v.sort()
col = collections.defaultdict(list)  # member key string -> [(s, e, column)]
for l in open(columns):
    f = l.rstrip('\n').split('\t')
    if f[0] == 'family': continue
    for x in f[8].split(','):
        k, r = x.rsplit('#', 1); a, b = r.split(':')[1].split('-'); col[k].append((int(a), int(b), f[2]))

def member(c, s, e):
    """the family member (of a family of >= 2) holding most of the span, or None; -> (key string, family)"""
    best = None
    for a, b, fam in mem.get(c, ()):
        if a >= e: break
        o = min(b, e) - max(a, s)
        if o > 0 and (best is None or o > best[0]): best = (o, f'{c}:{a}-{b}', fam)
    return best and best[1:]

def colchain(key, ex):
    """the exons as columns of their member: tuple of column ids, or None when an exon is in no column"""
    cs = col.get(key)
    if not cs: return None
    outc = []
    for s, e in ex:
        b = max(((min(e, y) - max(s, x), c) for x, y, c in cs), default=(0, None))
        if b[0] <= 0: return None
        outc.append(b[1])
    return tuple(outc)

# RefSeq intron chains (exon lines grouped by Parent)
rex = collections.defaultdict(list)
for l in open(gff):
    if l[0] == '#': continue
    f = l.split('\t')
    if len(f) > 8 and f[2] == 'exon':
        m = re.search(r'Parent=([^;\n]+)', f[8]) or TID.search(f[8])   # GFF3 or GTF
        if m: rex[(f[0], f[6], m.group(1))].append((int(f[3]) - 1, int(f[4])))
ref = {chain(c, s, sorted(e)) for (c, s, _), e in rex.items() if len(e) > 1}

B = transcripts(base); basechains = {chain(c, s, e) for c, s, e, _ in B.values() if len(e) > 1}
# the two-read transcripts as column chains, per family: {(family, columns): [(member, internal exon lengths)]}
have = collections.defaultdict(list)
for t, (c, s, e, _) in B.items():
    if len(e) < 2: continue
    m = member(c, e[0][0], e[-1][1])
    if not m: continue
    cc = colchain(m[0], e)
    if cc: have[(m[1], cc)].append((m[0], [b - a for a, b in e[1:-1]]))
D = {t: v for t, v in transcripts(mc).items() if 'multicopy_rescue "1"' in v[3] and len(v[2]) > 1}
Dchains = {chain(c, s, e) for c, s, e, _ in D.values()}
n = collections.Counter(); rows = []; seen = set()
for t, (c, s, e, at) in transcripts(f1).items():
    if len(e) < 2: continue
    ch = chain(c, s, e)
    if ch in basechains or ch in seen: continue
    seen.add(ch)
    m = member(c, e[0][0], e[-1][1])
    if not m: n['one-read chains outside a family member'] += 1; continue
    n['universe: new chains inside a family member'] += 1
    cc = colchain(m[0], e); inD = ch in Dchains; ann = ch in ref
    lens = [b - a for a, b in e[1:-1]]
    sib = [x for x in have.get((m[1], cc), ()) if x[0] != m[0]] if cc else []
    loose = bool(sib); strict = any(all(abs(p - q) <= 3 for p, q in zip(lens, x[1])) for x in sib)
    rows.append((t, m[1], m[0], len(e), int(inD), int(strict), int(loose), int(ann), 'no column' if cc is None else 'ok'))
    for name, flag in (('D', inD), ('G', strict), ('Gloose', loose)):
        if flag: n[name] += 1; n[name + ' annotated'] += ann
    if inD and strict: n['D and G'] += 1; n['D and G annotated'] += ann
    if inD and loose: n['D and Gloose annotated'] += ann
    if strict and not inD: n['G at an occupied place or outside D'] += 1; n['G at an occupied place or outside D, annotated'] += ann
    if cc is None: n['chains with an exon in no column'] += 1
with open(out, 'w') as w:
    w.write('transcript\tfamily\tmember\texons\tin_D\tG\tG_loose\tannotated\tcolumns\n'); w.writelines('\t'.join(map(str, r)) + '\n' for r in rows)
for k, v in sorted(n.items()): print(f'  {k}\t{v}')
sh = lambda a, b: f'{a} of {b} ({a / b:.3f})' if b else f'{a} of 0'
print(f"rescued transcripts tagged in the multi-copy GTF: {len(D)} ({len(Dchains)} chains); in the universe: {n['D']}")
print(f"(1) lost: D's annotated chains that G keeps {sh(n['D and G annotated'], n['D annotated'])}  bar >= 0.80   (loose arm {sh(n['D and Gloose annotated'], n['D annotated'])})")
print(f"(2) additions: annotated share of G {sh(n['G annotated'], n['G'])} against D {sh(n['D annotated'], n['D'])}   (loose arm {sh(n['Gloose annotated'], n['Gloose'])})")
