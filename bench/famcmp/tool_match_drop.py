#!/usr/bin/env python3
"""r1633/r1634 drop side: what are the chains one assembly has and the other lacks, relative to the other's chains?

  tool_match_drop.py REFSEQ A_GTF B_GTF        classes of A-only chains against B's chains on the same contig and strand:
near = same intron count, every junction within 25 bp of a B chain's; sub = contiguous run of a B chain's introns;
super = contains a B chain as a contiguous run; other_shared = shares >= 1 intron with B; novel_locus = shares none."""
import sys, collections
sys.path.insert(0, __import__('os').path.dirname(__file__))
from tool_match import read
def classify(ch, Bset, by_first, by_intron, tol=25):
    n = len(ch)
    for d in range(-tol, tol + 1):
        for b in by_first.get((n, ch[0][0] + d), ()):
            if all(abs(x[0] - y[0]) <= tol and abs(x[1] - y[1]) <= tol for x, y in zip(ch, b)): return 'near'
    cand = set().union(*(by_intron.get(i, ()) for i in ch))
    if not cand: return 'novel_locus'
    run = lambda a, b: any(b[k:k + len(a)] == a for k in range(len(b) - len(a) + 1))
    if any(run(ch, b) for b in cand): return 'sub'
    if any(run(b, ch) for b in cand): return 'super'
    return 'other_shared'
if __name__ == '__main__':
    if sys.argv[1] == 'selftest':
        B = [((10, 20), (30, 40), (50, 60))]; bf = collections.defaultdict(list); bi = collections.defaultdict(set)
        for b in B:
            bf[(len(b), b[0][0])].append(b)
            for i in b: bi[i].add(b)
        f = lambda ch: classify(ch, set(B), bf, bi)
        assert [f(((12, 20), (30, 41), (50, 60))), f(((30, 40), (50, 60))), f(((1, 5), (10, 20), (30, 40), (50, 60))), f(((10, 20), (33, 99))), f(((70, 80),))] == ['near', 'sub', 'super', 'other_shared', 'novel_locus']
        print('ok'); sys.exit()
    ref, A, B = (read(p) for p in sys.argv[1:4]); tab = collections.defaultdict(lambda: [0, 0])
    idx = collections.defaultdict(lambda: (collections.defaultdict(list), collections.defaultdict(set)))
    for c, s, b in B:
        bf, bi = idx[(c, s)]; bf[(len(b), b[0][0])].append(b)
        for i in b: bi[i].add(b)
    for k in A - B:
        bf, bi = idx[k[:2]]; v = tab[classify(k[2], B, bf, bi)]; v[0] += 1; v[1] += k in ref
    print('class\tchains\tcorrect\tshare')
    for c, (n, ok) in sorted(tab.items(), key=lambda x: -x[1][0]): print(f'{c}\t{n}\t{ok}\t{ok / n:.3f}')
