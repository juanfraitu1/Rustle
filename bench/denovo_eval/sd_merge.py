#!/usr/bin/env python3
"""SD-PROJECTION MERGE (2026-10-30, user "can we somehow overlap segmental duplications?"): join two sketch families when a
segmental-duplication pair (SEDEF table: CHM13 SD90 for human, GGO sedef for gorilla) PROJECTS a locus of one family onto a locus of
the other. Projection, not overlap: an SD block holds many genes (NPIP next to PKD1), so "both loci inside SD pairs" would glue
neighbours; instead each locus's exons are mapped through the SD pair (linear interpolation inside the pair, strand-aware) and the
edge is made only if the projected exons cover >= --cov of the exonic bp of a locus on the other side (and vice versa is not needed).
  sd_merge.py LOCI_GFF3 SD_BED(.gz) COPIES_TSV OUT_PREFIX [--min-id 0.9] [--cov 0.5] [--max-pairs 0] [--gate-edges SKETCH.edges.tsv --gate-cont 0.1]
Writes OUT.sd_edges.tsv (locus, locus, sd identity, cov_a, cov_b) and OUT.fam.copies.tsv (families merged by SD edges; loci without a
family that gain an SD edge to a family member are NOT added: this is a merge rule, members still come from the sketch)."""
import sys, os, gzip, csv, re, collections, bisect
a = sys.argv[1:]
opt = lambda k, d, f=float: f(a[a.index(k) + 1]) if k in a else d
MINID = opt('--min-id', 0.9); COV = opt('--cov', 0.5); MAXP = opt('--max-pairs', 0, int)
GATE = a[a.index('--gate-edges') + 1] if '--gate-edges' in a else None; GATEC = opt('--gate-cont', 0.1)
MINFRAC = opt('--min-frac', 0.0)  # merge two families only if SD edges touch >= this fraction of the SMALLER family's loci (one pseudogene edge must not merge two 30-member families)  # merge only SD edges that also have a sketch edge >= gate-cont
pos = [x for i, x in enumerate(a) if not x.startswith('--') and (i == 0 or not a[i - 1].startswith('--'))]
GFF, SD, COPIES, OUT = pos[:4]

loci = {}; exons = collections.defaultdict(list)
for line in open(GFF):
    if line[0] == '#': continue
    f = line.rstrip('\n').split('\t')
    if len(f) < 9: continue
    if f[2] == 'gene': loci[re.search(r'ID=gene-([^;]+)', f[8]).group(1)] = f[0]
    elif f[2] == 'exon': exons[re.search(r'gene=([^;]+)', f[8]).group(1)].append((int(f[3]) - 1, int(f[4])))
for t in exons: exons[t].sort()
bychrom = collections.defaultdict(list)  # chrom -> sorted (start, end, locus) by exon
for t, ex in exons.items():
    for s, e in ex: bychrom[loci[t]].append((s, e, t))
for c in bychrom: bychrom[c].sort()
starts = {c: [x[0] for x in v] for c, v in bychrom.items()}
exlen = {t: sum(e - s for s, e in ex) for t, ex in exons.items()}

def ident(f):  # first field after column 12 that parses as a fraction in (0, 1]
    for x in f[12:]:
        try:
            v = float(x)
            if 0 < v <= 1: return v
        except ValueError: pass
    return None

def overlapping(c, s, e):
    v = bychrom.get(c, [])
    i = bisect.bisect_left(starts.get(c, []), s - 1_000_000)
    out = []
    while i < len(v) and v[i][0] < e:
        if v[i][1] > s: out.append(v[i])
        i += 1
    return out

def cov_by(c, ivs, t):  # exonic bp of locus t covered by intervals ivs (on chrom c)
    if loci[t] != c: return 0.0
    cov = 0
    for s, e in exons[t]:
        for a_, b_ in ivs:
            lo, hi = max(s, a_), min(e, b_)
            if hi > lo: cov += hi - lo
    return cov / exlen[t]

edges = {}
n_sd = 0
op = gzip.open if SD.endswith('.gz') else open
for line in op(SD, 'rt'):
    if line[0] == '#': continue
    f = line.rstrip('\n').split('\t')
    if len(f) < 10: continue
    try: c1, s1, e1, c2, s2, e2 = f[0], int(f[1]), int(f[2]), f[3], int(f[4]), int(f[5])
    except ValueError: continue
    idn = ident(f)
    if idn is None or idn < MINID: continue
    rev = (f[8] != f[9]) if len(f) > 9 and f[8] in '+-' and f[9] in '+-' else False
    n_sd += 1
    if MAXP and n_sd > MAXP: break
    L1, L2 = e1 - s1, e2 - s2
    if L1 <= 0 or L2 <= 0: continue
    for (ca, sa, ea, cb, sb, eb) in ((c1, s1, e1, c2, s2, e2), (c2, s2, e2, c1, s1, e1)):
        La, Lb = ea - sa, eb - sb
        for t in {x[2] for x in overlapping(ca, sa, ea)}:
            # project exons of t through the pair
            proj = []
            for s, e in exons[t]:
                lo, hi = max(s, sa), min(e, ea)
                if hi <= lo: continue
                fa_, fb_ = (lo - sa) / La, (hi - sa) / La
                if rev: fa_, fb_ = 1 - fb_, 1 - fa_
                proj.append((sb + int(fa_ * Lb), sb + int(fb_ * Lb)))
            if not proj: continue
            ca_t = sum(hi - lo for lo, hi in proj) / exlen[t]  # share of t's exons inside the pair
            if ca_t < COV: continue
            cands = {x[2] for lo, hi in proj for x in overlapping(cb, lo, hi)}
            for u in cands:
                if u == t: continue
                cb_u = cov_by(cb, proj, u)
                if cb_u < COV: continue
                k = (min(t, u), max(t, u))
                if k not in edges or idn > edges[k][0]: edges[k] = (idn, ca_t, cb_u)
with open(OUT + '.sd_edges.tsv', 'w') as w:
    for (t, u), (idn, ca_t, cb_u) in sorted(edges.items()): w.write(f"{t}\t{u}\t{idn:.4f}\t{ca_t:.2f}\t{cb_u:.2f}\n")

rows = list(csv.DictReader(open(COPIES), delimiter='\t'))
fam = {r['tid']: r['family_id'] for r in rows}
gp = {f_: f_ for f_ in set(fam.values())}
def find(x):
    while gp[x] != x: gp[x] = gp[gp[x]]; x = gp[x]
    return x
sk = None
if GATE:
    sk = {}
    for line in open(GATE):
        x, y, _, c = line.split('\t')[:4]
        if float(c) >= GATEC: sk[(min(x, y), max(x, y))] = float(c)
    edges = {k: v for k, v in edges.items() if k in sk}
    print(f"gated by sketch edges >= {GATEC}: {len(edges)} SD edges remain", file=sys.stderr)
nm = 0
if MINFRAC > 0:
    members = collections.defaultdict(set)
    for t, f_ in fam.items(): members[f_].add(t)
    touch = collections.defaultdict(set)  # (famA, famB) -> loci of A and of B touched by SD edges
    for (t, u) in edges:
        if t in fam and u in fam and fam[t] != fam[u]:
            k = (min(fam[t], fam[u]), max(fam[t], fam[u])); touch[k].add(t); touch[k].add(u)
    cands = []
    for (fa_, fb_), s_ in touch.items():
        small = min(len(members[fa_]), len(members[fb_]))
        fr = min(len(s_ & members[fa_]) / len(members[fa_]), len(s_ & members[fb_]) / len(members[fb_]))
        if fr >= MINFRAC: cands.append((fa_, fb_))
    for fa_, fb_ in cands:
        if find(fa_) != find(fb_): gp[find(fa_)] = find(fb_); nm += 1
else:
    for (t, u) in edges:
        if t in fam and u in fam and find(fam[t]) != find(fam[u]): gp[find(fam[t])] = find(fam[u]); nm += 1
with open(OUT + '.fam.copies.tsv', 'w') as w:
    w.write('\t'.join(rows[0].keys()) + '\n')
    for r in rows:
        r = dict(r); r['family_id'] = find(r['family_id']); w.write('\t'.join(r.values()) + '\n')
print(f"SD pairs >= {MINID}: {n_sd}; locus-projection edges {len(edges)}; family merges {nm}; families {len(set(find(f_) for f_ in fam.values()))}", file=sys.stderr)
