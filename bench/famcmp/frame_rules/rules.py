"""rules.py TAG REF_GTF DEFAULT_GTF FINAL_GTF : frozen rules of register 1766 on one sample (needs TAG.new.tsv, TAG.dorf.tsv).
P1 drop default chains whose longest ORF stops > 50 nt before the last junction.
S1 admit a 1-read-floor chain: no such stop, ORF >= 900 nt, one read with exactly the chain, no longer read containing it, locus = none | extends.
S2 as S1 with ORF >= 300 nt and any locus class but 'inside', plus >= 1 sister-species read with exactly the lifted chain."""
import sys, collections
sys.path.insert(0, '/home/hpc/v138da/v138da13/Rustle/bench/famcmp'); from allmiss import chains
n, ref_p, def_p, fin_p = sys.argv[1:5]; key = lambda c, ch: (c, ','.join('%d-%d' % i for i in ch))
ref = {key(c, ch) for c, s, ch in chains(ref_p)}; fin = {key(c, ch) for c, s, ch in chains(fin_p)}
D = {}  # default chain -> (reads, ptc)
for l in open(n + '.dorf.tsv'):
    f = l.rstrip('\n').split('\t')
    if f[0] != 'chrom': k = (f[0], f[1]); D[k] = (max(int(f[3]), D.get(k, (0, 0))[0]), int(f[6]))
two = [k for k, v in D.items() if v[0] == 2]; bar = sum(k in ref for k in two) / len(two)
print('%s default: %d chains, %d correct, %.3f; 2-read chains %d, %.3f = bar for new chains; yield bar %d' % (n, len(D), sum(k in ref for k in D), sum(k in ref for k in D) / len(D), len(two), bar, round(.02 * sum(k in ref for k in D))))
def curve(S, target):  # precision of S cut by reads down to `target` correct chains
    c = t = 0
    for k in sorted(S, key=lambda k: -D[k][0]):
        t += 1; c += k in ref
        if c >= target: break
    return c / t, t
for name, S in (('ours.gtf', set(D)), ('ours.final.gtf', set(D) & fin)):
    drop = {k for k in S if D[k][1]}; keep = S - drop; ck = sum(k in ref for k in keep); cd = sum(k in ref for k in drop); cp, ct = curve(S, ck)
    print('P1 on %s: drops %d chains, %d correct (%.3f; %.3f of correct lost); kept %d, %.3f vs read cut at the same correct count %.3f (%d chains)' % (name, len(drop), cd, cd / len(drop), cd / (ck + cd), len(keep), ck / len(keep), cp, ct))
S1 = set(); S2 = set()
for l in open(n + '.new.tsv'):
    f = l.rstrip('\n').split('\t')
    if f[0] == 'chrom' or f[5] != '0' or f[9] == 'NA' or int(f[10]) > 0: continue
    k = (f[0], f[13]); o = int(f[4])
    if o >= 900 and f[9] in ('none', 'extends'): S1.add(k)
    if o >= 300 and f[9] != 'inside' and int(f[7]) >= 1: S2.add(k)
for name, S in (('S1', S1), ('S2', S2), ('S1 or S2', S1 | S2)):
    c = sum(k in ref for k in S); print('%s: %d new chains, %d correct, %.3f' % (name, len(S), c, c / len(S)))
for name, base in (('ours.gtf', set(D)), ('ours.final.gtf', set(D) & fin)):
    for lab, add in (('P1 + S1', S1), ('P1 + S1 + S2', S1 | S2)):
        R = {k for k in base if not D[k][1]} | add; c = sum(k in ref for k in R); cb = sum(k in ref for k in base)
        print('%s -> %s: %d / %d / %.3f  (was %d / %d / %.3f)' % (name, lab, len(R), c, c / len(R), len(base), cb, cb / len(base)))
pickle_out = n + '.rules.keys'
with open(pickle_out, 'w') as w:
    for k in S1: w.write('S1\t%s\t%s\n' % k)
    for k in S2: w.write('S2\t%s\t%s\n' % k)
