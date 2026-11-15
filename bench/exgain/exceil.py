"""EXCEIL (bench/exceil.py) on the rebuilt units + gate G3. Needs edlib (the repo's bench env).
usage: exceil.py DIR SEX.tsv BAM REF.fa MAT.fa PAT.fa ; writes DIR/exceil.tsv; exit 3 when G3 fails."""
import collections, random, sys
import edlib, pysam
from common import sexset, load_hits, merged
D, sexf, BAM = sys.argv[1:4]
FA = dict(zip(('ref', 'mat', 'pat'), sys.argv[4:7]))
COMP = str.maketrans('ACGTN', 'TGCAN')
def rc(s): return s.translate(COMP)[::-1]
H = load_hits(D, sexset(sexf))
uid = {}; gain = {}
for l in open(f'{D}/units.tsv'):
    f = l.rstrip('\n').split('\t')
    if f[0] != 'unit': gain[f[0]] = f[6] == '1'; uid[f[0]] = f[7]
nets = collections.defaultdict(list)
for l in open(f'{D}/net.tsv'):
    f = l.split('\t')
    if f[0] != 'family': nets[f[0]].append(f[1])
fa = {g: pysam.FastaFile(p) for g, p in FA.items()}
def seqs(g, ex):
    return [(rc if st == '-' else str)(fa[g].fetch(t, s, e).upper()) for t, s, e, st in merged(H[g][ex])]
bam = pysam.AlignmentFile(BAM)
loci_of = collections.defaultdict(list)
for l in open(f'{D}/copies.tsv'):
    f = l.rstrip('\n').split('\t')
    if f[0] != 'family_id': loci_of[f[0]].append((f[3], int(f[4]), int(f[5]), f[7]))
def level(u):
    ex = uid[u]; R = seqs('ref', ex); rs = set(R)
    hap_only = set()
    for g in ('mat', 'pat'):
        X = seqs(g, ex)
        if len(X) > len(R) or not gain[u]:
            hap_only |= {x for x in X if x not in rs}
    if not hap_only: return 'none', 0
    want = set(nets[u]); hits = 0
    for c, s, e, st in loci_of[u]:
        for r in bam.fetch(c, s, e):
            if r.query_name not in want or r.is_secondary or r.is_supplementary: continue
            want.discard(r.query_name)
            qp = [q for q, p in r.get_aligned_pairs(matches_only=True) if s - 20 <= p < e + 20]
            if not qp: continue
            seg = r.query_sequence[min(qp):max(qp) + 1]
            seg = rc(seg) if st == '-' else seg
            d = lambda x: edlib.align(x, seg, mode='HW')['editDistance']
            if min(d(x) for x in hap_only) < min(d(x) for x in rs): hits += 1
    return ('B' if hits >= 2 else 'A'), hits
G = sorted(u for u in uid if gain.get(u))
rng = random.Random(20261025); C = sorted(rng.sample(sorted(u for u in uid if u in gain and not gain[u]), 300))
share = {}
with open(f'{D}/exceil.tsv', 'w') as o:
    o.write('unit\tgroup\tlevel\thits\tnet_reads\n')
    for name, us in (('gain', G), ('control', C)):
        c = collections.Counter()
        for u in us:
            lv, h = level(u); c[lv] += 1; o.write(f'{u}\t{name}\t{lv}\t{h}\t{len(nets[u])}\n')
        n = len(us); share[name] = 100 * c['B'] / n
        print(f"{name} units n={n}: level A (sequence) {c['A'] + c['B']} ({100 * (c['A'] + c['B']) / n:.1f}%), level B (reads) {c['B']} "
              f"({share[name]:.1f}%), none {c['none']}")
ok = abs(share['gain'] - 41.9) <= 10 and abs(share['control'] - 14.0) <= 10
print(f"G3 level B gaining {share['gain']:.1f}% (41.9 +- 10), control {share['control']:.1f}% (14.0 +- 10): {'PASS' if ok else 'FAIL'}")
sys.exit(0 if ok else 3)
