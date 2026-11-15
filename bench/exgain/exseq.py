"""Row 1711 blind exon-gain detector (docs/PREREG_o3_exon_gain_newdetector_2026-11-14.md). Reads no truth.
usage: exseq.py DIR BAM REF.fa K N  -> DIR/exseq_K.tsv (batch K of N by sha1 of the unit) | exseq.py --selftest"""
import collections, hashlib, sys
MIN_READS, MIN_SHARE, MIN_ED = 3, 0.10, 2
COMP = str.maketrans('ACGTN', 'TGCAN')
def rc(s): return s.translate(COMP)[::-1]
def call(segs, refs, dist):
    """segs: read exon sequences; refs: reference locus sequences. -> (supported [(seq, n, min dist to refs)], H, N)"""
    c = collections.Counter(segs); n = len(segs)
    sup = [(s, k, min(dist(s, r) for r in refs)) for s, k in c.most_common() if k >= MIN_READS and k >= MIN_SHARE * n]
    return sup, len(sup) > 2 * len(refs), any(d >= MIN_ED for _, _, d in sup)
def exon_seq(r, s, e, strand):
    q = {p: q for q, p in r.get_aligned_pairs(matches_only=True) if p == s or p == e - 1}
    if len(q) < 2: return None
    seg = r.query_sequence[q[s]:q[e - 1] + 1].upper()
    return rc(seg) if strand == '-' else seg
if __name__ == '__main__':
    import edlib
    dist = lambda a, b: edlib.align(a, b)['editDistance']
    if sys.argv[1] == '--selftest':
        R = 'ACGTACGTAC'; A = 'ACGTACGTAT'; B = 'ACGAACGAAC'; C = 'TCGTACGTAC'
        assert call([R] * 5 + [A] * 5, [R], dist)[1:] == (False, False)            # two alleles, one edit: silent
        assert call([R] * 5 + [A] * 5 + [C] * 5, [R], dist)[1:] == (True, False)   # three haplotypes at one locus
        assert call([R] * 5 + [B] * 5, [R], dist)[1:] == (False, True)             # divergent sequence
        assert call([R] * 30 + [B] * 2, [R], dist)[1:] == (False, False)           # under the support floor
        print('selftest ok'); sys.exit()
    import pysam
    D, BAM, REF, K, N = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5])
    mine = lambda u: int(hashlib.sha1(u.encode()).hexdigest(), 16) % N == K
    loci = collections.defaultdict(list); nets = collections.defaultdict(set)
    for l in open(f'{D}/copies.tsv'):
        f = l.rstrip('\n').split('\t')
        if f[0] != 'family_id' and mine(f[0]): loci[f[0]].append((f[3], int(f[4]), int(f[5]), f[7]))
    for l in open(f'{D}/net.tsv'):
        f = l.split('\t')
        if f[0] in loci: nets[f[0]].add(f[1])
    bam = pysam.AlignmentFile(BAM); fa = pysam.FastaFile(REF)
    with open(f'{D}/exseq_{K}.tsv', 'w') as o:
        for u in sorted(nets):
            refs = [(rc if st == '-' else str)(fa.fetch(c, s, e).upper()) for c, s, e, st in loci[u]]
            segs = []; seen = set()
            for c, s, e, st in loci[u]:
                for r in bam.fetch(c, s, e):
                    if r.query_name in nets[u] and r.query_name not in seen and not (r.is_secondary or r.is_supplementary):
                        seen.add(r.query_name); x = exon_seq(r, s, e, st)
                        if x: segs.append(x)
            sup, H, Nf = call(segs, refs, dist) if segs else ([], False, False)
            o.write(f"{u}\t{len(loci[u])}\t{len(segs)}\t{int(H)}\t{int(Nf)}\t{','.join(f'{s}:{k}:{d}' for s, k, d in sup)}\n")
