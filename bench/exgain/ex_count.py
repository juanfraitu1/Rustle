"""Exon copy counts (bench/sdcn_count2.py, EX part) + validity gates G1, G2 of Amendment 1.
usage: ex_count.py DIR SEX.tsv   (DIR holds exons.bed, ex_sr_{ref,mat,pat}.paf); exit 3 when a gate fails."""
import collections, sys
from common import sexset
D, sex = sys.argv[1], sexset(sys.argv[2])
def loci(xs):
    iv = collections.defaultdict(list)
    for t, s, e in xs: iv[t].append((s, e))
    k = 0
    for t, v in iv.items():
        end = -1
        for s, e in sorted(v):
            if s >= end: k += 1
            end = max(end, e)
    return k
H = {g: collections.defaultdict(list) for g in ('ref', 'mat', 'pat')}
for g in H:
    for l in open(f'{D}/ex_sr_{g}.paf'):
        f = l.split('\t')
        H[g][f[0]].append(((int(f[3]) - int(f[2])) / int(f[1]), int(f[9]) / int(f[10]), f[5], int(f[7]), int(f[8])))
# G1 (frozen before any count; the original QC code is not in the repo): an exon's own locus is recovered when a reference hit at
# coverage >= 0.8 and identity >= 0.90 overlaps the exon's own coordinates.
exons = [l.split('\t') for l in open(f'{D}/exons.bed')]
own = sum(any(cov >= 0.8 and idn >= 0.90 and t == c and a < int(e) and b > int(s) for cov, idn, t, a, b in H['ref'][name.strip()])
          for c, s, e, name in exons)
g1 = own / len(exons)
U = sorted(set(H['ref']) | set(H['mat']) | set(H['pat']))
capped = [u for u in U if any(len(H[g][u]) >= 400 for g in H)]
U = [u for u in U if u not in set(capped)]
stable = None; out = open(f'{D}/ex_counts.tsv', 'w'); out.write('exon\tfamily\tn_ref\tn_mat\tn_pat\n')
for s in [(c, t) for c in (0.5, 0.8, 0.95) for t in (0.90, 0.95, 0.98)]:
    n = {u: [loci([(t, a, b) for cov, idn, t, a, b in H[g][u] if cov >= s[0] and idn >= s[1] and (g, t) not in sex])
             for g in ('ref', 'mat', 'pat')] for u in U}
    d = {u for u in U if n[u][1] != n[u][2]}; stable = d if stable is None else stable & d
    if s == (0.8, 0.90):
        for u in U: out.write(f"{u}\t{u.split('|')[1]}\t{n[u][0]}\t{n[u][1]}\t{n[u][2]}\n")
        fams = {u.split('|')[1] for u in U}; df = {u.split('|')[1] for u in d}; g2 = len(df) / len(fams)
        print(f'exons {len(exons)}, units counted {len(U)} (capped {len(capped)}); mat!=pat {len(d)} ({100 * len(d) / len(U):.1f}%); '
              f'ref!=mat {100 * sum(n[u][0] != n[u][1] for u in U) / len(U):.1f}%, ref!=pat {100 * sum(n[u][0] != n[u][2] for u in U) / len(U):.1f}%')
        print(f'families with >= 1 differing exon {len(df)}/{len(fams)} ({100 * g2:.1f}%)')
print(f'differ under all 9 settings: {len(stable)} ({100 * len(stable) / len(U):.1f}%)')
ok1, ok2 = g1 >= 0.99, abs(100 * g2 - 24.0) <= 8
print(f'G1 own locus recovered {own}/{len(exons)} ({100 * g1:.2f}%) bar >= 99%: {"PASS" if ok1 else "FAIL"}')
print(f'G2 families with a mat/pat-differing exon {100 * g2:.1f}% bar 24.0 +- 8: {"PASS" if ok2 else "FAIL"}')
sys.exit(0 if ok1 and ok2 else 3)
