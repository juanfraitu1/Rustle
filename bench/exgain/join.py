"""Row 1694 join: among EXCEIL-VISIBLE (level B) gaining and control units, the share on which EXHAP's certified tier
(PRE.families.tsv certified > 0) and putative tier (>= 1 row in PRE.putative.tsv) fire. Denominator = all VISIBLE units; a
VISIBLE unit absent from the copy_phasing output counts as not fired. usage: join.py DIR PRE ; writes DIR/join.tsv"""
import csv, collections, sys
D, pre = sys.argv[1], sys.argv[2]
rd = lambda p: list(csv.DictReader(open(p), delimiter='\t'))
run = {r['family']: int(r['certified']) for r in rd(f'{pre}.families.tsv')}
put = collections.Counter(r['family'] for r in rd(f'{pre}.putative.tsv'))
units = {r['unit']: r['gain'] for r in rd(f'{D}/units.tsv')}
for g in ('1', '0'):
    us = [u for u in units if units[u] == g and u in run]
    print(f"all units gain={g}: {len([u for u in units if units[u] == g])}, run {len(us)}, certified {sum(run[u] > 0 for u in us)}, putative {sum(put[u] > 0 for u in us)}")
ex = rd(f'{D}/exceil.tsv'); sh = {}
with open(f'{D}/join.tsv', 'w') as o:
    o.write('group\tunits\tvisible\tvisible_run\tcertified\tputative\tcertified_share\tputative_share\n')
    for grp in ('gain', 'control'):
        allu = [r['unit'] for r in ex if r['group'] == grp]; vis = [r['unit'] for r in ex if r['group'] == grp and r['level'] == 'B']
        c = sum(run.get(u, 0) > 0 for u in vis); p = sum(put[u] > 0 for u in vis); n = max(1, len(vis)); sh[grp] = (c / n, p / n)
        o.write(f"{grp}\t{len(allu)}\t{len(vis)}\t{sum(u in run for u in vis)}\t{c}\t{p}\t{c / n:.4f}\t{p / n:.4f}\n")
with open(f'{D}/join_units.tsv', 'w') as o:
    o.write('unit\tgroup\tlevel\thits\trun\tcertified\tputative\n')
    for r in ex: o.write(f"{r['unit']}\t{r['group']}\t{r['level']}\t{r['hits']}\t{int(r['unit'] in run)}\t{run.get(r['unit'], 0)}\t{put[r['unit']]}\n")
print(open(f'{D}/join.tsv').read())
# reading rule (prereg, written before the join)
best = max(sh['gain'])
print('reading: "the detector is the limit" ' + ('STANDS (best tier < 10% of VISIBLE gaining units)' if best < 0.10 else 'does not stand by the < 10% clause'))
for i, t in enumerate(('certified', 'putative')):
    fu = sh['gain'][i] >= 0.25 and sh['control'][i] <= sh['gain'][i] / 3
    print(f'  {t}: follow-up row {"WARRANTED" if fu else "not proposed"} (gaining {sh["gain"][i]:.3f}, control {sh["control"][i]:.3f})')
