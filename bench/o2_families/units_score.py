#!/usr/bin/env python3
"""Unit-level and divergence-class scores of assemblies of planted reads (register 1906).

    units_score.py PLANT_DIR SCENARIO LABEL=GTF [...] [--families F1,F2,...] [--out TABLE.tsv]

Needs PLANT_DIR/units.tsv and dist.tsv (units.py). A unit = copies with identical planted isoforms; expressed when any
member has planted reads. Per tool: unit recall and precision; silent copies with a transcript and expressed copies
without one, by the edits to the nearest expressed (for silent) / other expressed (for missed) copy of the family:
0, 1-5, 6-20, > 20 (no aligned partner counts as > 20). --families restricts everything to those families."""
import sys, os, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import score as S
a = sys.argv[1:]; fams = None; out = None
if '--families' in a: i = a.index('--families'); fams = set(a[i + 1].split(',')); del a[i:i + 2]
if '--out' in a: i = a.index('--out'); out = a[i + 1]; del a[i:i + 2]
d, sc = a[0], a[1]
truth, planted, fam = S.load_plan(d, sc)
creads = collections.Counter()
for t, n in planted.items(): creads[truth[t][3]] += n
U = {l.split('\t')[1]: l.split('\t')[3] for l in list(open(f'{d}/units.tsv'))[1:]}
D = collections.defaultdict(dict)
for l in list(open(f'{d}/dist.tsv'))[1:]:
    f = l.rstrip('\n').split('\t'); D[f[1]][f[2]] = int(f[3]); D[f[2]][f[1]] = int(f[3])
copies = [g for g in creads if fams is None or fam[g] in fams]
exp = {g for g in copies if creads[g] > 0}
def cls(g):
    e = min((D[g][x] for x in D[g] if x in exp and x != g), default=None)
    if U[g] in {U[x] for x in exp if x != g}: e = 0
    return '0' if e == 0 else '1-5' if e is not None and e <= 5 else '6-20' if e is not None and e <= 20 else '>20'
C = {g: cls(g) for g in copies}; K = ('0', '1-5', '6-20', '>20')
units = collections.defaultdict(set)
for g in copies: units[U[g]].add(g)
uexp = {u for u, gs in units.items() if gs & exp}
hdr = ['scenario', 'tool', 'units_expressed', 'units_found', 'unit_recall', 'units_silent', 'false_units', 'unit_precision'] + [f'silent_called_{k}' for k in K] + [f'silent_{k}' for k in K] + [f'missed_{k}' for k in K] + [f'expressed_{k}' for k in K]
rows = []
for x in a[2:]:
    lab, path = x.split('=', 1)
    if not os.path.exists(path): rows.append([sc, lab, 'MISSING']); continue
    called = S.score(truth, planted, fam, S.load_gtf(path))['_called'] & set(copies)
    ucall = {u for u, gs in units.items() if gs & called}
    r = lambda p, q: f'{p / q:.3f}' if q else 'NA'
    rows.append([sc, lab, len(uexp), len(ucall & uexp), r(len(ucall & uexp), len(uexp)), len(units) - len(uexp), len(ucall - uexp), r(len(ucall & uexp), len(ucall))]
                + [sum(1 for g in copies if g not in exp and g in called and C[g] == k) for k in K] + [sum(1 for g in copies if g not in exp and C[g] == k) for k in K]
                + [sum(1 for g in exp if g not in called and C[g] == k) for k in K] + [sum(1 for g in exp if C[g] == k) for k in K])
txt = '\t'.join(hdr) + '\n' + ''.join('\t'.join(map(str, r)) + '\n' for r in rows)
sys.stdout.write(txt)
if out: open(out, 'w').write(txt)
