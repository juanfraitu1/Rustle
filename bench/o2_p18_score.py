# P18 scorer: minimap2's primary (rank by ms) vs an error-model log-likelihood on near-tie reads
import sys, json, math, collections, pysam
from scipy.stats import binomtest
SPLIT = sys.argv[1]; BAND = 10; POST = 0.99
R = json.load(open('errprof.json'))
LX = math.log(R['X']); P_TOT = sum(v for k, v in R.items() if k[:2] in ('X', 'I_', 'D_') or k == 'X'); LM = math.log(1 - P_TOT)
FA = pysam.FastaFile('/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta')
def run(ref, p, b):
    n = 0; i = p
    while i < len(ref) and ref[i] == b: n += 1; i += 1
    i = p - 1
    while i >= 0 and ref[i] == b: n += 1; i -= 1
    return n
def ll(r):
    ref = FA.fetch(r.reference_name, r.reference_start, r.reference_end).upper(); q = r.query_sequence.upper()
    rp = qp = 0; s = 0.0; introns = 0
    for op, L in r.cigartuples:
        if op == 7: s += L * LM; rp += L; qp += L
        elif op == 8: s += L * LX; rp += L; qp += L
        elif op == 1:
            x = q[qp:qp + L]; ctx = 'hp' if len(set(x)) == 1 and run(ref, rp, x[0]) >= 1 else 'nhp'
            s += math.log(R[f'I_{ctx}_{min(L, 3)}']); qp += L
        elif op == 2:
            x = ref[rp:rp + L]; ctx = 'hp' if len(set(x)) == 1 and run(ref, rp, x[0]) >= L + 1 else 'nhp'
            s += math.log(R[f'D_{ctx}_{min(L, 3)}']); rp += L
        elif op == 3: rp += L; introns += 1
        elif op == 4: s += L * LX; qp += L       # an unexplained (soft-clipped) base costs a substitution
    return s, introns
truth = {l.split('\t')[0]: l.rstrip('\n').split('\t') for l in open(f'{SPLIT}.truth.tsv') if not l.startswith('read\t')}
recs = collections.defaultdict(list)
for r in pysam.AlignmentFile(f'{SPLIT}.bam'):
    if r.is_unmapped or r.is_supplementary: continue
    recs[r.query_name].append(r)
def right(r, t):
    return r.reference_name == t[3] and r.reference_start < int(t[5]) and r.reference_end > int(t[4])
c = collections.Counter(); rows = []
for name, rs in recs.items():
    t = truth[name]; p = next(r for r in rs if not r.is_secondary)
    ms = lambda r: r.get_tag('ms')
    band = [r for r in rs if ms(r) >= ms(p) - BAND]
    c['reads'] += 1; c['primary right (all reads)'] += right(p, t)
    if len(band) < 2: continue
    # secondaries carry no SEQ (-Y keeps clips but SEQ is '*'): borrow the primary's sequence, oriented
    psq = p.query_sequence
    sc = []
    for r in band:
        if r.query_sequence is None:
            r.query_sequence = psq if r.is_reverse == p.is_reverse else psq.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]
        sc.append((ll(r), r))
    best = max(sc, key=lambda x: x[0][0]); lls = [x[0][0] for x in sc]
    m = max(lls); post = math.exp(best[0][0] - m) / sum(math.exp(v - m) for v in lls)
    nm = min(band, key=lambda r: (r.get_tag('NM'), r is not p))
    exact = sum(1 for v in lls if abs(v - m) < 1e-9) > 1
    c['band reads'] += 1; c['band: primary right'] += right(p, t)
    c['band: truth among candidates'] += any(right(r, t) for r in band)
    if exact: c['band: LL exact tie (left to primary)'] += 1; continue
    sw = best[1] is not p
    c['band: LL argmax right'] += right(best[1], t); c['band: fewest-edits right'] += right(nm, t)
    if sw:
        c['switches'] += 1; c['switch right'] += right(best[1], t); c['switch wrong (broke a right primary)'] += right(p, t)
    if post >= POST:
        c['posterior>=0.99 assigned'] += 1; c['posterior>=0.99 wrong'] += not right(best[1], t); c['posterior>=0.99 primary wrong on same reads'] += not right(p, t)
    rows.append((name, t[1], t[6], ms(p), len(band), right(p, t), right(best[1], t), sw, round(post, 4)))
n = c['band reads'] - c['band: LL exact tie (left to primary)']
for k in c: print(f'{k}: {c[k]:,}')
if c['switches']:
    print('H1 binomial p (switch right > 1/2):', binomtest(c['switch right'], c['switches'], 0.5, alternative='greater').pvalue)
open(f'{SPLIT}.rows.tsv', 'w').write('read\tfamily\tallele\tms\tn_band\tprimary_right\tll_right\tswitch\tposterior\n' + '\n'.join('\t'.join(map(str, x)) for x in rows) + '\n')
