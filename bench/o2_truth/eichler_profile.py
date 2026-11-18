#!/usr/bin/env python3
"""PREREG r1786 (docs/PREREG_truth_repair_eichler_profile_2026-11-17.md): how far the adapted Eichler filter carries, from neutral_truth's table.

  eichler_profile.py READS.tsv S0_TOTAL S1_TOTAL S2_TOTAL OUT.json     totals = the strata sizes before the caps"""
import sys, json, collections

BINS = [(0, 0, '0'), (1, 5, '1-5'), (6, 10, '6-10'), (11, 20, '11-20'), (21, 50, '21-50'), (51, 200, '51-200'), (201, 10 ** 9, '>200')]
SIZES = [(0, 999, '<1000'), (1000, 1999, '1000-1999'), (2000, 3999, '2000-3999'), (4000, 10 ** 9, '>=4000')]
def bin_of(x, B): return next(n for a, b, n in B if a <= x <= b)
assert bin_of(10, BINS) == '6-10' and bin_of(11, BINS) == '11-20' and bin_of(4000, SIZES) == '>=4000'

def main(tsv, n0, n1, n2, out):
    tot = {'S0': int(n0), 'S1': int(n1), 'S2': int(n2)}; took = collections.Counter(); R = []
    for l in open(tsv):
        f = l.rstrip('\n').split('\t')
        if f[0] == 'read': continue
        took[f[1]] += 1
        if f[1] == 'S0': continue
        R.append(dict(st=f[1], lead=int(f[2]), status=f[3], primary=int(f[4].split(';')[0].split('=')[1]), origin={int(i) for i in f[5].split(',') if i}, AS=int(f[6].split(';')[0].rsplit(':AS', 1)[1])))
    w = {k: tot[k] / took[k] for k in took}; n = sum(tot.values()); res = {}
    res['M1'] = {k: tot[k] / n for k in ('S0', 'S1', 'S2')}; res['totals'] = tot
    print(f"M1 reach: one placement {tot['S0']} ({tot['S0'] / n:.3f}), clear lead {tot['S1']} ({tot['S1'] / n:.3f}), near tie {tot['S2']} ({tot['S2'] / n:.3f})")
    C = collections.Counter()
    for r in R:
        if r['status'] != 'decidable': continue
        b = bin_of(r['lead'], BINS); C[(b, 'n')] += 1; C[(b, 'best')] += 0 in r['origin']; C[(b, 'prim')] += r['primary'] in r['origin']
        if r['primary'] != 0: C[(b, 'dis')] += 1; C[(b, 'dis_best')] += 0 in r['origin']; C[(b, 'dis_prim')] += r['primary'] in r['origin']
    print("\nM2 lead profile, decidable reads (sample counts)\nlead\treads\tbest_AS_right\tbest_wrong_share\tprimary_right\tprimary_wrong_share\tbest_is_not_primary\tthere_best_right\tthere_primary_right")
    res['M2'] = []
    for _, _, b in BINS:
        m = C[(b, 'n')]
        if not m: continue
        row = dict(lead=b, n=m, best_right=C[(b, 'best')], prim_right=C[(b, 'prim')], dis=C[(b, 'dis')], dis_best=C[(b, 'dis_best')], dis_prim=C[(b, 'dis_prim')]); res['M2'].append(row)
        print(f"{b}\t{m}\t{row['best_right']}\t{1 - row['best_right'] / m:.4f}\t{row['prim_right']}\t{1 - row['prim_right'] / m:.4f}\t{row['dis']}\t{row['dis_best']}\t{row['dis_prim']}")
    g = lambda bs, k: sum(C[(b, k)] for b in bs)
    hi = ['11-20', '21-50']; a = g(hi, 'n'); safe = 1 - g(hi, 'best') / a if a else None
    w610 = 1 - C[('6-10', 'best')] / C[('6-10', 'n')] if C[('6-10', 'n')] else None; w1120 = 1 - C[('11-20', 'best')] / C[('11-20', 'n')] if C[('11-20', 'n')] else None
    p610 = C[('6-10', 'prim')] / C[('6-10', 'n')] if C[('6-10', 'n')] else None
    res['readings'] = dict(best_wrong_11_50=safe, n_11_50=a, best_wrong_6_10=w610, best_wrong_11_20=w1120, primary_right_6_10=p610, n_6_10=C[('6-10', 'n')], n_11_20=C[('11-20', 'n')])
    print(f"\nreadings: best-AS wrong in lead 11-50: {safe} (n {a}); best-AS wrong 6-10 {w610} vs 11-20 {w1120}; primary right in 6-10 {p610}")
    allS1 = [r for r in R if r['st'] == 'S1']; d = [r for r in allS1 if r['primary'] != 0]; dd = [r for r in d if r['status'] == 'decidable']
    res['M3'] = dict(sample=len(allS1), disagree=len(d), share=len(d) / len(allS1), decidable=len(dd), best_right=sum(0 in r['origin'] for r in dd), prim_right=sum(r['primary'] in r['origin'] for r in dd))
    print(f"\nM3 clear-lead reads whose best-AS locus is not the primary: {len(d)} of {len(allS1)} sampled ({len(d) / len(allS1):.4f}); decidable {len(dd)}: best AS right {res['M3']['best_right']}, primary right {res['M3']['prim_right']}")
    S = collections.Counter()
    for r in R: S[(bin_of(r['AS'], SIZES), r['st'])] += w[r['st']]
    print("\nM4 near-tie share by the best alignment's AS (weighted to the full strata)\nAS\tclear_lead\tnear_tie\tnear_tie_share"); res['M4'] = []
    for _, _, b in SIZES:
        a1, a2 = S[(b, 'S1')], S[(b, 'S2')]
        if a1 + a2: res['M4'].append(dict(size=b, clear=round(a1), near=round(a2), share=a2 / (a1 + a2))); print(f"{b}\t{a1:.0f}\t{a2:.0f}\t{a2 / (a1 + a2):.4f}")
    json.dump(res, open(out, 'w'), indent=1)

if __name__ == '__main__': main(*sys.argv[1:6])
