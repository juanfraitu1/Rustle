# P19: abundance EM over the sample's reads; assign a near-tie read only at posterior >= 0.99
import sys, math, collections, bisect, pysam
SPLIT = sys.argv[1]; BAND = 10; POST = 0.99; ITERS = 2000
exec("import json\nR = json.load(open('errprof.json'))" )
src = open('score.py').read()
exec(src[src.index("LX = "):src.index("truth = ")])          # LX, LM, FA, run(), ll()
comp = str.maketrans('ACGTN', 'TGCAN'); rc = lambda s: s.translate(comp)[::-1]
# catalog copies -> units
spans = collections.defaultdict(list)
for l in open('/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_gwcat.copies.tsv'):
    if l.startswith('family_id'): continue
    f = l.split('\t'); spans[f[3]].append((int(f[4]), int(f[5]), f'{f[0]}|{f[1]}'))
def unit(r):
    for s, e, u in spans.get(r.reference_name, ()):
        if r.reference_start < e and r.reference_end > s: return u
    return f'{r.reference_name}:{r.reference_start // 10000}'
truth = {l.split('\t')[0]: l.rstrip('\n').split('\t') for l in open(f'p19.{SPLIT}.truth.tsv') if not l.startswith('read\t')}
recs = collections.defaultdict(list)
for r in pysam.AlignmentFile(f'p19.{SPLIT}.bam'):
    if not (r.is_unmapped or r.is_supplementary): recs[r.query_name].append(r)
reads = {}   # read -> {unit: likelihood (relative)}, primary unit, tie flag, truth unit
for name, rs in recs.items():
    p = next(r for r in rs if not r.is_secondary)
    band = [r for r in rs if r.get_tag('ms') >= p.get_tag('ms') - BAND]
    L = {}
    for r in band:
        if r.query_sequence is None:
            r.query_sequence = p.query_sequence if r.is_reverse == p.is_reverse else rc(p.query_sequence)
        u = unit(r); v = ll(r)[0]; L[u] = max(L.get(u, -1e18), v)
    m = max(L.values()); L = {u: math.exp(v - m) for u, v in L.items()}
    t = truth[name]
    reads[name] = (L, unit(p), sum(1 for v in L.values() if v > 1 - 1e-9) > 1, f'{t[1]}|{t[2]}')
units = sorted({u for L, *_ in reads.values() for u in L})
th = {u: 1 / len(units) for u in units}
for it in range(ITERS):
    acc = collections.Counter()
    for L, *_ in reads.values():
        z = sum(th[u] * v for u, v in L.items())
        for u, v in L.items(): acc[u] += th[u] * v / z
    new = {u: acc[u] / len(reads) for u in units}
    delta = max(abs(new[u] - th[u]) for u in units); th = new
    if delta < 1e-13: break
print(f'EM: {len(reads):,} reads, {len(units):,} units, {it + 1} iterations, last change {delta:.1e}')
c = collections.Counter()
for name, (L, pu, tie, tu) in reads.items():
    if len(L) < 2: continue
    z = sum(th[u] * v for u, v in L.items()); post = {u: th[u] * v / z for u, v in L.items()}
    best = max(post, key=post.get)
    g = 'LL-tied reads' if tie else 'decidable band reads'
    c[(g, 'reads')] += 1; c[(g, 'primary right')] += pu == tu; c[(g, 'truth among candidates')] += tu in L
    if post[best] >= POST:
        c[(g, 'EM assigned')] += 1; c[(g, 'EM assigned wrong')] += best != tu
        c[(g, 'EM assigned wrong | primary was right')] += best != tu and pu == tu
        c[(g, 'EM assigned right | primary was wrong')] += best == tu and pu != tu
    else: c[(g, 'EM abstains')] += 1; c[(g, 'EM abstains | primary right')] += pu == tu
for k in sorted(c): print(k, c[k])
