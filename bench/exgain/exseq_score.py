"""Row 1711 scorer: the only script of the row that reads truth. usage: exseq_score.py DIR SEX.tsv REF.fa MAT.fa PAT.fa"""
import collections, glob, hashlib, sys
import pysam
from common import sexset, load_hits, merged
from exseq import rc
D, sexf = sys.argv[1:3]; fa = {g: pysam.FastaFile(p) for g, p in zip(('ref', 'mat', 'pat'), sys.argv[3:6])}
H = load_hits(D, sexset(sexf))
seqs = lambda g, ex: {(rc if st == '-' else str)(fa[g].fetch(t, s, e).upper()) for t, s, e, st in merged(H[g][ex])}
U = {f[0]: (f[1], f[6] == '1', f[7]) for f in (l.rstrip('\n').split('\t') for l in open(f'{D}/units.tsv')) if f[0] != 'unit'}
vis = {f[0] for f in (l.split('\t') for l in open(f'{D}/exceil.tsv')) if f[1] == 'gain' and f[2] == 'B'}
rows = [l.rstrip('\n').split('\t') for p in sorted(glob.glob(f'{D}/exseq_*.tsv')) for l in open(p)]
half = lambda u: int(hashlib.sha1(U[u][0].encode()).hexdigest(), 16) % 2
print(f'units scored {len(rows)} (gaining {sum(U[r[0]][1] for r in rows)}, with >= 1 exon-sequence read {sum(int(r[2]) > 0 for r in rows)}); VISIBLE gaining {len(vis)}')
for arm, col in (('H', 3), ('N', 4)):
    fired = [r for r in rows if r[col] == '1']; g = [r for r in fired if U[r[0]][1]]; right = 0
    for r in g:
        ex = U[r[0]][2]; hap_only = (seqs('mat', ex) | seqs('pat', ex)) - seqs('ref', ex)
        right += any(x.split(':')[0] in hap_only for x in r[5].split(','))
    v = sum(r[0] in vis for r in g)
    print(f"arm {arm}: fires gaining {len(g)} / non-gaining {len(fired) - len(g)}; precision {len(g) / max(1, len(fired)):.3f}; recall {len(g)}/469; "
          f"VISIBLE {v}/{len(vis)} ({v / len(vis):.3f}); right reason {right}/{len(g)}; "
          + '; '.join(f"half {h}: {sum(half(r[0]) == h for r in g)} / {sum(half(r[0]) == h for r in fired) - sum(half(r[0]) == h for r in g)}" for h in (0, 1))
          + f" -> {'PASS' if len(g) >= 0.5 * max(1, len(fired)) and v >= 0.25 * len(vis) and fired else 'FAIL'}")
