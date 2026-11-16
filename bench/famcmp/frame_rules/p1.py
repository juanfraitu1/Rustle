"""p1.py DORF.tsv REF_GTF FINAL_GTF : rule P1 alone on a default assembly (kept precision vs read cut at equal correct count) + RefSeq class of lost chains."""
import sys, re, collections
sys.path.insert(0, '/home/hpc/v138da/v138da13/Rustle/bench/famcmp'); from allmiss import chains
key = lambda c, ch: (c, ','.join('%d-%d' % i for i in ch)); ref = {key(c, ch) for c, s, ch in chains(sys.argv[2])}; fin = {key(c, ch) for c, s, ch in chains(sys.argv[3])}
D = {}
for l in open(sys.argv[1]):
    f = l.rstrip('\n').split('\t')
    if f[0] != 'chrom': k = (f[0], f[1]); D[k] = (max(int(f[3]), D.get(k, (0, 0))[0]), int(f[6]))
for name, S in (('ours.gtf', set(D)), ('final.gtf', set(D) & fin)):
    drop = {k for k in S if D[k][1]}; keep = S - drop; ck = sum(k in ref for k in keep); cd = sum(k in ref for k in drop); c = t = 0
    for k in sorted(S, key=lambda k: -D[k][0]):
        t += 1; c += k in ref
        if c >= ck: break
    print('%s: %d / %d / %.3f -> %d / %d / %.3f; dropped %d at %.3f, %.3f of correct lost; read cut at equal correct %.3f' % (name, len(S), ck + cd, (ck + cd) / len(S), len(keep), ck, ck / len(keep), len(drop), cd / len(drop), cd / (ck + cd), c / t))
ex = collections.defaultdict(list)
for l in open(sys.argv[2]):
    f = l.split('\t')
    if len(f) > 8 and f[2] == 'exon': ex[(f[0], re.search(r'transcript_id "([^"]+)"', f[8]).group(1))].append((int(f[3]), int(f[4])))
cls = collections.defaultdict(set)
for (c, t), e in ex.items():
    e.sort()
    if len(e) > 1: cls[(c, ','.join('%d-%d' % (e[i][1], e[i + 1][0] - 1) for i in range(len(e) - 1)))].add(re.sub(r'^rna-', '', t).split('_')[0])
n = collections.Counter(); d = collections.Counter()
for k in D:
    if k in ref:
        g = 'coding (NM/XM)' if cls[k] & {'NM', 'XM'} else 'non-coding (NR/XR)' if cls[k] & {'NR', 'XR'} else 'other'; n[g] += 1; d[g] += D[k][1]
print('  correct chains lost by RefSeq class: ' + '; '.join('%s %d of %d (%.3f)' % (g, d[g], n[g], d[g] / n[g]) for g in sorted(n)))
