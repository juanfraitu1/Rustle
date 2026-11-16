"""transplant.py REF_GTF READS.pkl OURS_GTF OUT : a read's local chain r set into a model M of ours that shares r's first and/or last intron.
J = M[:p] + r + M[q+1:]; kept when J is not a model of ours and takes >= 1 intron from M alone. Label = RefSeq chain, strand ignored."""
import sys, pickle, collections
sys.path.insert(0, '/home/hpc/v138da/v138da13/Rustle/bench/famcmp'); from allmiss import chains
ref_p, pkl, ours_p, out = sys.argv[1:5]
ref = {(c, ch) for c, s, ch in chains(ref_p)}; ours = {(c, ch) for c, s, ch in chains(ours_p)}; prim = pickle.load(open(pkl, 'rb'))
oi = collections.defaultdict(set)
for c, ch in ours:
    for i in ch: oi[(c, i)].add(ch)
J = {}
for (c, r), n in prim.items():
    if (c, r) in ours: continue
    L = oi.get((c, r[0]), set()); R = oi.get((c, r[-1]), set())
    for M in L | R:
        p = M.index(r[0]) if M in L else None; q = M.index(r[-1]) if M in R else None
        if p is not None and q is not None and q < p: continue
        if p is not None and M[p:p + len(r)] == r[:len(M) - p] and (q is not None or len(r) <= len(M) - p): continue  # r runs along M: no event
        if q is not None and p is None and len(r) <= q + 1 and M[q + 1 - len(r):q + 1] == r: continue
        j = (M[:p] if p else ()) + r + (M[q + 1:] if q is not None else ())
        ext = len(j) - len(r)
        if ext == 0 or (c, j) in ours or list(j) != sorted(set(j)) or any(j[k][1] >= j[k + 1][0] for k in range(len(j) - 1)): continue
        kind = 'both' if p is not None and q is not None else 'left' if p is not None else 'right'
        o = J.setdefault((c, j), [0, 0, 0, set(), set()]); mf = prim[(c, M)]
        if r not in o[3]: o[0] += n; o[3].add(r)
        o[1] = max(o[1], mf); o[2] = max(o[2], ext); o[4].add(kind)
with open(out, 'w') as w:
    w.write('chrom\tn_introns\tev_reads\tm_full\text\tkinds\tfull\tcorrect\tchain\n')
    for (c, j), o in J.items():
        w.write('\t'.join(map(str, [c, len(j), o[0], o[1], o[2], '+'.join(sorted(o[4])), prim[(c, j)], int((c, j) in ref), ','.join('%d-%d' % i for i in j)])) + '\n')
