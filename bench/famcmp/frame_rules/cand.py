"""cand.py REF_GTF READS.pkl OURS_GTF FASTA OUT : candidate chains ours does not output, with features and RefSeq label (strand ignored).
single = chain of exactly one primary record; stitch = join of two read chains a (upstream) + b sharing >= 1 intron, no record with the whole join."""
import sys, pickle, collections, pysam
sys.path.insert(0, '/home/hpc/v138da/v138da13/Rustle/bench/famcmp'); from allmiss import chains
ref_p, pkl, ours_p, fa_p, out = sys.argv[1:6]
ref = {(c, ch) for c, s, ch in chains(ref_p)}; ours = {(c, ch) for c, s, ch in chains(ours_p)}
prim = pickle.load(open(pkl, 'rb')); fa = pysam.FastaFile(fa_p)
intron = collections.Counter(); first = collections.defaultdict(list); has = collections.defaultdict(list)
for (c, ch), n in prim.items():
    first[(c, ch[0])].append(ch)
    for i in ch: intron[(c, i)] += n; has[(c, i)].append(ch)
oi = collections.defaultdict(list)
for c, ch in ours:
    for i in ch: oi[(c, i)].append(ch)
def run(a, b):
    if len(a) > len(b) or a[0] not in b: return False
    k = b.index(a[0]); return b[k:k + len(a)] == a
def canon(c, i):
    m = (fa.fetch(c, i[0], i[0] + 2) + fa.fetch(c, i[1] - 2, i[1])).upper(); return m in ('GTAG', 'CTAC', 'GCAG', 'CTGC', 'ATAC', 'GTAT')
def locus(c, ch):
    O = {o for i in ch for o in oi.get((c, i), ())}
    if not O: return 'none'
    if any(run(ch, o) for o in O): return 'inside'
    if any(run(o, ch) for o in O): return 'extends'
    return 'alt'
def feats(c, ch):
    longer = sum(prim[(c, b)] for b in set(has[(c, ch[0])]) if len(b) > len(ch) and run(ch, b))
    ooi = sum(1 for i in ch if (c, i) not in oi)  # introns in no model of ours
    return [len(ch), min(intron[(c, i)] for i in ch), locus(c, ch), longer, ooi, int(all(canon(c, i) for i in ch)), int((c, ch) in ref)]
with open(out, 'w') as w:
    w.write('kind\tchrom\tn_introns\tmin_intron\tlocus\tlonger_reads\tintrons_new\tcanon\tcorrect\tov\tsup\tnpairs\tna\tnb\tchain\n')
    J = {}
    for (c, a), na in prim.items():
        if na == 1 and (c, a) not in ours: w.write('\t'.join(map(str, ['single', c] + feats(c, a) + [0, 1, 0, 0, 0, ','.join('%d-%d' % i for i in a)])) + '\n')
        for j in range(1, len(a)):
            L = len(a) - j
            for b in first.get((c, a[j]), ()):
                if len(b) > L and b[:L] == a[j:]:
                    k = (c, a + b[L:]); nb = prim[(c, b)]; o = J.get(k)
                    J[k] = (max(o[0], L), max(o[1], min(na, nb)), o[2] + 1, max(o[3], na), max(o[4], nb)) if o else (L, min(na, nb), 1, na, nb)
    for (c, ch), v in J.items():
        if (c, ch) in prim or (c, ch) in ours: continue
        w.write('\t'.join(map(str, ['stitch', c] + feats(c, ch) + list(v) + [','.join('%d-%d' % i for i in ch)])) + '\n')
