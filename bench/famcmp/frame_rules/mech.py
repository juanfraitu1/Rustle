import sys, collections
sys.path.insert(0, '/home/hpc/v138da/v138da13/Rustle/bench/famcmp'); from allmiss import chains
tsv, ours_p = sys.argv[1:3]; ours = chains(ours_p)
byi = collections.defaultdict(set)
for c, s, ch in ours:
    for i in ch: byi[(c, s, i)].add(ch)
def run(a, b):  # a contiguous run of b -> offset or -1
    if len(a) > len(b) or a[0] not in b: return -1
    k = b.index(a[0]); return k if b[k:k + len(a)] == a else -1
cnt = collections.Counter(); H = None
for l in open(tsv):
    f = l.rstrip('\n').split('\t')
    if H is None: H = f; continue
    r = dict(zip(H, f))
    if r['tier'] not in ('full', 'linked') or any(int(r[k]) for k in H[9:-1]): continue
    c, s = r['chrom'], r['strand']; C = tuple(tuple(map(int, x.split('-'))) for x in r['chain'].split(','))
    O = set().union(*[byi.get((c, s, i), set()) for i in C])
    if not O: m = 'no model of ours shares an intron'
    else:
        sup = [(o, run(C, o)) for o in O]; sup = [(o, k) for o, k in sup if k >= 0]
        sub = [(o, run(o, C)) for o in O]; sub = [(o, k) for o, k in sub if k >= 0]
        if sup:
            o, k = max(sup, key=lambda x: len(x[0])); lo, hi = k > 0, k + len(C) < len(o)
            if s == '-': lo, hi = hi, lo
            m = 'inside a longer model of ours: ref lacks ' + ("both ends" if lo and hi else "5' introns" if lo else "3' introns")
        elif sub:
            o, k = max(sub, key=lambda x: len(x[0])); lo, hi = k > 0, k + len(o) < len(C)
            if s == '-': lo, hi = hi, lo
            m = 'ours is a part of it: ours lacks ' + ("both ends" if lo and hi else "5' introns" if lo else "3' introns")
        else:
            best = max(len(set(o) & set(C)) / len(set(o) | set(C)) for o in O)
            m = 'different combination, best intron Jaccard ' + ('>= 0.8' if best >= .8 else '0.5-0.8' if best >= .5 else '< 0.5')
    cnt[(r['tier'], m)] += 1
for k, v in sorted(cnt.items()): print(v, *k, sep='\t')
