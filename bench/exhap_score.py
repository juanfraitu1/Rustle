"""EXHAP scorer (LEDGER_O3). usage: exhap_score.py PREFIX SPLIT (in the exhap directory)"""
import csv, collections, sys
from scipy.stats import fisher_exact
pre, split = sys.argv[1], sys.argv[2]
U = {r['unit']: r for r in csv.DictReader(open('units.tsv'), delimiter='\t') if r['split'] == split}
T = {r['family']: r['gain'] for r in csv.DictReader(open('truth.tsv'), delimiter='\t') if r['split'] == split}
run = {r['family']: r for r in csv.DictReader(open(f'{pre}.families.tsv'), delimiter='\t')}
put = collections.Counter(r['family'] for r in csv.DictReader(open(f'{pre}.putative.tsv'), delimiter='\t'))
def rep(name, flagged, truth, pos):
    P = sum(v == pos for v in truth.values()); N = len(truth) - P
    tp = sum(truth[x] == pos for x in flagged); fp = len(flagged) - tp
    p = fisher_exact([[tp, fp], [P - tp, N - fp]], alternative='greater')[1]
    print(f'{name}: flagged {len(flagged)} ({pos} {tp}/{P}, other {fp}/{N}); sens {tp / max(1, P):.3f} '
          f'prec {tp / max(1, len(flagged)):.3f} base {P / len(truth):.3f} p={p:.3g}')
ut = {u: ('1' if r['gain'] == '1' else '0') for u, r in U.items() if u in run}
print(f'units run {len(ut)} of {len(U)}')
cu = {u for u in ut if int(run[u]['certified']) > 0}
pu = {u for u in ut if put[u] > 0}
rep('units, certified', cu, ut, '1'); rep('units, putative', pu, ut, '1'); rep('units, certified or putative', cu | pu, ut, '1')
fam = lambda us: {U[u]['family'] for u in us}
rep('families, certified', fam(cu), T, 'GAIN'); rep('families, certified or putative', fam(cu | pu), T, 'GAIN')
