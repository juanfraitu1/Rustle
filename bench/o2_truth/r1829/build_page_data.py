#!/usr/bin/env python3
"""r1829 page data: four JSON files {rows, meta} for the near-tie dashboard, each carrying the code that made it.
  build_page_data.py WORKDIR SEL.bam OUTDIR
WORKDIR holds census.tsv (layer 1) and choice.tsv (primary_choice.py); SEL.bam = the HG002 layer-2 selection (all records of
30,000 near-tied reads, residues 0-24)."""
import collections, json, os, sys
import numpy as np, pysam

W, SEL, OUT = sys.argv[1:4]
CODE = open(os.path.abspath(__file__)).read()
SAMPLES = ['A119b', 'HG002', 'OR6737', 'KB3781', 'PTR_testis', 'PTR_AG18354_lcl', 'PAB_testis', 'PAB_AG06213_fib',
           'PPA_testis', 'PPA_PR00251_fib', 'PPY_testis', 'PPY_AG05252_fib']
INFO = {'A119b': ('human', 'testis'), 'HG002': ('human', 'cell line (HG002)'), 'OR6737': ('gorilla', 'testis'),
        'KB3781': ('gorilla', 'fibroblast'), 'PTR_testis': ('chimpanzee', 'testis'), 'PTR_AG18354_lcl': ('chimpanzee', 'LCL'),
        'PAB_testis': ('Sumatran orangutan', 'testis'), 'PAB_AG06213_fib': ('Sumatran orangutan', 'fibroblast'),
        'PPA_testis': ('bonobo', 'testis'), 'PPA_PR00251_fib': ('bonobo', 'fibroblast'),
        'PPY_testis': ('Bornean orangutan', 'testis'), 'PPY_AG05252_fib': ('Bornean orangutan', 'fibroblast')}

def dump(name, rows, inputs):
    meta = [{'label': 'build_page_data.py', 'value': CODE, 'type': 'code', 'lang': 'python'}] + \
           [{'label': 'Input', 'value': x} for x in inputs]
    json.dump({'rows': rows, 'meta': meta}, open(os.path.join(OUT, name + '.json'), 'w'))
    print(name, len(rows))

# 1. census: tiers per sample
C = collections.defaultdict(dict)
for l in open(os.path.join(W, 'census.tsv')):
    f = l.rstrip('\n').split('\t')
    if f[0] != 'sample':
        C[f[0]][f[1]] = int(f[2])
rows = []
for s in SAMPLES:
    m = C[s]['molecules']
    rows.append({'sample': s, 'species': INFO[s][0], 'tissue': INFO[s][1], 'molecules': m,
                 'multimapper_share': round(C[s].get('two_places', 0) / m, 5),
                 'exact_share': round(C[s].get('exact', 0) / m, 5), 'ms10_share': round(C[s].get('ms10', 0) / m, 5),
                 'good_share': round(C[s].get('good', 0) / m, 5)})
dump('census', rows, [os.path.join(W, 'census.tsv') + ' (near_tie_census.py over each BAM\'s as_table)'])

# 2. primary choice: per sample, metric, gap band: primary higher / equal / rival higher
P = collections.defaultdict(lambda: collections.Counter())
tot = collections.defaultdict(dict)
for l in open(os.path.join(W, 'choice.tsv')):
    f = l.rstrip('\n').split('\t')
    if f[0] == 'sample':
        continue
    if f[1] == '-':
        tot[f[0]][f[3]] = int(f[4])
    else:
        P[(f[0], f[1], f[2])][f[3]] += int(f[4])
BANDS = ['0', '1-5', '6-10', '11-25', '26-50', '51-100', '>100']
rows = []
for s in SAMPLES:
    for m in ('AS', 'ms'):
        for b in BANDS:
            c = P[(s, m, b)]
            n = sum(c.values())
            rows.append({'sample': s, 'metric': m, 'gap': b, 'reads': n, 'primary_higher': c['primary higher'],
                         'equal': c['equal'], 'rival_higher': c['rival higher'],
                         'multimappers': tot[s].get('multimappers', 0)})
dump('primary_choice', rows, [os.path.join(W, 'choice.tsv') + ' (primary_choice.py over each BAM\'s as_table)'])

# 3. event costs: AS and ms fitted on the HG002 selection's records
def feats(r):
    c = collections.Counter()
    for op, n in r.cigartuples:
        if op == 7: c['match'] += n
        elif op == 8: c['mismatch'] += n
        elif op in (1, 2): c['gap_open'] += 1; c['gap_len'] += n
        elif op == 3: c['intron'] += 1
    return c
K = ['match', 'mismatch', 'gap_open', 'gap_len', 'intron']
X, Y, M, R = [], [], [], collections.defaultdict(list)
for r in pysam.AlignmentFile(SEL, check_sq=False):
    if r.is_unmapped or r.is_supplementary:
        continue
    f = feats(r)
    X.append([f[k] for k in K]); Y.append(r.get_tag('AS')); M.append(r.get_tag('ms'))
    R[r.query_name].append((r.get_tag('ms'), r.get_tag('AS'), r.reference_name, r.is_reverse, r.reference_start, f))
X = np.array(X, float); A = np.c_[X, np.ones(len(X))]
co = {n: np.linalg.lstsq(A, np.array(y, float), rcond=None)[0] for n, y in (('AS', Y), ('ms', M))}
label = {'match': 'match (per base)', 'mismatch': 'mismatched base (vs +1 for a match: a substitution costs 5)',
         'gap_open': 'indel opened', 'gap_len': 'indel, per base', 'intron': 'intron'}
rows = [{'event': label[k], 'AS_fitted': round(float(co['AS'][i]), 2), 'ms_fitted': round(float(co['ms'][i]), 2)}
        for i, k in enumerate(K)]
exact = float(np.mean(np.array(Y) == X[:, 0] - 4 * X[:, 1] - 6 * X[:, 2] - 1 * X[:, 3] - 24 * X[:, 4]))
rows.append({'event': 'records the integer model (+1, -4, -6 open, -1 per base, -24 per intron) reproduces exactly',
             'AS_fitted': round(exact, 3), 'ms_fitted': None})
rows.append({'event': 'records fitted', 'AS_fitted': len(Y), 'ms_fitted': len(M)})
dump('event_costs', rows, [SEL + ' (HG002 Kinnex on CHM13, all records of 30,000 near-tied reads, hash residues 0-24)'])

# 4. what differs between a read's two best places, by ms gap
G = collections.Counter(); per = collections.defaultdict(list)
for q, v in R.items():
    v.sort(key=lambda x: -x[0]); a = v[0]
    b = next((x for x in v[1:] if not (x[2] == a[2] and x[3] == a[3] and abs(x[4] - a[4]) < 5000)), None)
    if b is None:
        continue
    ea, eb = a[5], b[5]; dms = a[0] - b[0]
    dx = eb['mismatch'] - ea['mismatch']
    dgap = eb['gap_open'] - ea['gap_open']
    dint = ea['intron'] - eb['intron']
    dlen = (ea['match'] + ea['mismatch']) - (eb['match'] + eb['mismatch'])
    if dint: k = 'intron count differs'
    elif abs(dlen) > 20: k = 'aligned length differs'
    elif dgap: k = 'indels differ'
    elif dx: k = 'substitutions only'
    else: k = 'identical alignments'
    band = '0' if dms == 0 else '1-5' if dms <= 5 else '6-10' if dms <= 10 else '>10'
    G[(band, k)] += 1
    if k == 'substitutions only' and dx: per['ms per substitution'].append(dms / dx)
rows = [{'ms_gap': b, 'difference': k, 'reads': n} for (b, k), n in sorted(G.items())]
dump('gap_events', rows, [SEL])
print('median ms per substitution', np.median(per['ms per substitution']))
