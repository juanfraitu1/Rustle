"""agree.py TAG : per new chain of TAG.new.tsv, the other primary reads at its junctions and whether they are colinear with it -> TAG.agree.tsv
A read sharing >= 1 intron with chain C is: agree (its chain is a contiguous part of C), longer (C is a contiguous part of it), else conflict.
Per junction: share of the reads having it that are colinear (agree, longer or the full read); jmin = the lowest share, jsolo = junctions only the full reads have; pairsolo = adjacent junction pairs no agree / longer read has."""
import sys, pickle, collections
n = sys.argv[1]; prim = pickle.load(open(n + '.t0.reads.pkl', 'rb')); has = collections.defaultdict(list)
for (c, ch), k in prim.items():
    for i in ch: has[(c, i)].append(ch)
def run(a, b):
    if len(a) > len(b) or a[0] not in b: return False
    k = b.index(a[0]); return b[k:k + len(a)] == a
H = None
with open(n + '.agree.tsv', 'w') as w:
    for l in open(n + '.new.tsv'):
        f = l.rstrip('\n').split('\t')
        if H is None: H = f; w.write(l.rstrip('\n') + '\tfull\tagree\tlonger\tconflict\tjmin\tjsolo\tpairsolo\n'); continue
        c = f[0]; C = tuple(tuple(map(int, x.split('-'))) for x in f[-1].split(',')); cls = {}
        for i in C:
            for r in has.get((c, i), ()):
                if r not in cls: cls[r] = 'full' if r == C else 'agree' if run(r, C) else 'longer' if run(C, r) else 'conflict'
        t = collections.Counter()
        for r, k in cls.items(): t[k] += prim[(c, r)]
        jmin = 1.0; solo = 0
        for i in C:
            tot = col = oth = 0
            for r in has.get((c, i), ()):
                k = prim[(c, r)]; tot += k
                if cls[r] != 'conflict': col += k
                if cls[r] != 'full': oth += k
            solo += oth == 0; jmin = min(jmin, col / tot if tot else 0)
        cov = set()
        for r, k in cls.items():
            if k in ('agree', 'longer'): cov.update(zip(r, r[1:]))
        psolo = sum(1 for p in zip(C, C[1:]) if p not in cov)
        w.write('\t'.join(f + [str(t['full']), str(t['agree']), str(t['longer']), str(t['conflict']), '%.3f' % jmin, str(solo), str(psolo)]) + '\n')
