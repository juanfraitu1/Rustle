"""usage: sdbio.py [SEDEF_BED] [GENOME.fai]   (env as famtruth: DENOVO_CAT_GFF, DENOVO_CAT_TSV, DENOVO_PARALOGS)
How do SDs relate to multi-copy gene families? Annotation-level (CHM13 CAT genes, human Compara paralog pairs, CHM13 SEDEF),
independent of Rustle's loci. Descriptive.
Per Compara age bin: share of protein-coding paralog pairs whose two genes sit in the two copies of one SD (span rule, as
--sd-damp), the identity of that SD, intra/inter-chromosomal, tandem (< 1 Mb), and for pairs with no SD: retrocopy signature
(one gene single-exon, the other multi-exon). Plus gene- and genome-level SD coverage."""
import sys, os, collections, bisect, random, statistics as st
sys.path.insert(0, os.path.expanduser('~/Rustle/bench/denovo_eval')); import famtruth as F
SEDEF = sys.argv[1] if len(sys.argv) > 1 else "/anvme/workspace/v138da13-rustle/sedef_dl/incoming/sedef_send/final_human.bed"  # CHM13 v2.0 SEDEF, 34 cols, identity col 21
FAI = sys.argv[2] if len(sys.argv) > 2 else os.path.expanduser("~/data/chm13/chm13v2.0.fa.fai")
BINS = [('young: human/Homininae/Hominoidea', {'Homo sapiens', 'Homininae', 'Hominoidea'}),
        ('Catarrhini/Simiiformes', {'Catarrhini', 'Simiiformes'}),
        ('Primates..Mammalia', {'Primates', 'Boreoeutheria', 'Eutheria', 'Theria', 'Mammalia'}),
        ('older than mammals', None)]
G = F.load_cat()
pc = {g for g, d in G.items() if d['biotype'] == 'protein_coding'}
def merged(bl):
    out = []
    for s, e in sorted(bl):
        if out and s <= out[-1][1]: out[-1][1] = max(out[-1][1], e)
        else: out.append([s, e])
    return out
EX = {g: merged(G[g]['exons']) for g in pc}
# SD records, both directions: by contig -> (s, e, pair id, side, identity)
by = collections.defaultdict(list); ident = []
for l in open(SEDEF):
    if l.startswith('#'): continue
    f = l.split('\t'); p = len(ident); ident.append(float(f[20]))
    by[f[0]].append((int(f[1]), int(f[2]), p, 0)); by[f[3]].append((int(f[4]), int(f[5]), p, 1))
for v in by.values(): v.sort()
st_ = {c: [x[0] for x in v] for c, v in by.items()}; mx = {c: max(x[1] - x[0] for x in v) for c, v in by.items()}
def holds(bl, s, e):
    tot = sum(y - x for x, y in bl); ins = sum(max(0, min(y, e) - max(x, s)) for x, y in bl)
    lo, hi = bl[0][0], bl[-1][1]
    return tot > 0 and 2 * ins >= tot and 2 * max(0, min(hi, e) - max(lo, s)) >= hi - lo
H = {}
for g in pc:
    c, bl = G[g]['chrom'], EX[g]; v = by.get(c, []); lo, hi = bl[0][0], bl[-1][1]; h = set()
    for i in range(bisect.bisect_left(st_.get(c, []), lo - mx.get(c, 0)), len(v)):
        s, e, p, side = v[i]
        if s >= hi: break
        if e > lo and holds(bl, s, e): h.add((p, side))
    H[g] = h
def link(a, b):
    ps = [p for p, side in H[a] if (p, 1 - side) in H[b]]
    return max((ident[p] for p in ps), default=None)
# Compara pairs by age bin
by_ens = collections.defaultdict(set)
for g in pc:
    for e in G[g]['ens']: by_ens[e].add(g)
pairs = collections.defaultdict(set)
for line in F._open(F.PARALOGS):
    f = line.rstrip('\n').split('\t')
    if len(f) < 4 or not f[1] or f[0].startswith('Gene'): continue
    lca = f[3]; b = next((n for n, s in BINS if s is None or lca in s), None)
    if b == BINS[-1][0] and any(lca in s for _, s in BINS[:-1]): continue
    for x in by_ens.get(f[0].split('.')[0], ()):
        for y in by_ens.get(f[1].split('.')[0], ()):
            if x != y and G[x]['name'] != G[y]['name']: pairs[b].add(tuple(sorted((x, y))))
random.seed(1)
print('bin\tpairs\tsampled\tSD-linked\tmedian SD identity\tSD id >= 0.98\tintra-chrom\ttandem<1Mb\tno-SD: retrocopy-like\tno-SD: intra-chrom')
for name, _ in BINS:
    ps = sorted(pairs[name]); smp = random.sample(ps, min(20000, len(ps)))
    lk, ids, intra, tandem, rc, nointra, nosd = 0, [], 0, 0, 0, 0, 0
    for a, b in smp:
        same = G[a]['chrom'] == G[b]['chrom']
        d = min(abs(EX[a][0][0] - EX[b][-1][1]), abs(EX[b][0][0] - EX[a][-1][1])) if same else None
        intra += same; tandem += same and d < 1_000_000
        i = link(a, b)
        if i is not None: lk += 1; ids.append(i)
        else:
            nosd += 1; nointra += same
            rc += (len(EX[a]) == 1) != (len(EX[b]) == 1)
    n = len(smp)
    print(f"{name}\t{len(ps)}\t{n}\t{lk / n:.3f}\t{st.median(ids) if ids else float('nan'):.4f}\t{sum(x >= 0.98 for x in ids) / max(1, len(ids)):.3f}\t"
          f"{intra / n:.3f}\t{tandem / n:.3f}\t{rc / max(1, nosd):.3f}\t{nointra / max(1, nosd):.3f}", flush=True)
# gene and genome level
young = {g for p in pairs[BINS[0][0]] | pairs[BINS[1][0]] for g in p}
in_sd = {g for g in pc if H[g]}
print(f"protein-coding genes {len(pc)}; inside an SD copy (span rule) {len(in_sd)} ({len(in_sd) / len(pc):.3f})")
print(f"  with a Simiiformes-or-younger paralog: {len(young)}, of which inside an SD {len(young & in_sd) / max(1, len(young)):.3f}")
print(f"  without one: {len(pc - young)}, of which inside an SD {len((pc - young) & in_sd) / max(1, len(pc - young)):.3f}")
glen = sum(int(l.split('\t')[1]) for l in open(FAI))
cov = 0
for c, v in by.items():
    end = -1
    for s, e, *_ in v:
        if e > end: cov += e - max(s, end); end = e
print(f"genome in SDs: {cov / glen:.3f} of {glen / 1e9:.2f} Gb")
