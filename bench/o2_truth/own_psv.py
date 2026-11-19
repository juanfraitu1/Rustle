#!/usr/bin/env python3
"""PREREG r1783 (docs/PREREG_own_psv_2026-11-17.md): place near-tied reads by copy differences seen in the individual's own
MAPQ-60 reads, used only where the two copies' own allele sets are disjoint.

  own_psv.py READS.tsv SEL.bam REGION.bam REGIONS.bed
  OP_CALLS=FILE also writes one line per read: read, own call, loose call (locus index or NA); r1829, the rule unchanged."""
import sys, os, bisect, collections
from multiprocessing import Pool
import numpy as np, pysam
from allele_graph import table

BIT = np.zeros(256, dtype=np.uint8)
for k, b in enumerate('ACGT'): BIT[ord(b)] = 1 << k
CMASK = np.array([sum(((m >> k) & 1) << (3 - k) for k in range(4)) for m in range(16)], dtype=np.uint8)  # mask of the complement bases
COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
assert CMASK[1] == 8 and CMASK[2] == 4 and CMASK[0b0101] == 0b1010

def aligned(start, cigar):
    """(query index, reference position) of every match / mismatch column, clips counted in the query index"""
    q = []; r = []; p = start; i = 0
    for op, n in cigar:
        if op in (0, 7, 8): q.append(np.arange(i, i + n)); r.append(np.arange(p, p + n)); p += n; i += n
        elif op in (1, 4, 5): i += n
        elif op in (2, 3): p += n
    return (np.concatenate(q), np.concatenate(r)) if q else (np.zeros(0, int), np.zeros(0, int))
_q, _r = aligned(10, [(4, 2), (7, 2), (3, 100), (8, 1), (1, 1), (7, 1)]); assert _q.tolist() == [2, 3, 4, 6] and _r.tolist() == [10, 11, 112, 113]

def pile(a):
    """own-allele mask of one region from the MAPQ-60 primaries"""
    bam_path, c, s, e = a
    cov = np.array(pysam.AlignmentFile(bam_path).count_coverage(c, s, e, quality_threshold=0,
        read_callback=lambda r: not (r.is_secondary or r.is_supplementary or r.is_unmapped) and r.mapping_quality >= 60))
    tot = cov.sum(0); m = np.zeros(e - s, dtype=np.uint8)
    for k in range(4): m |= (((cov[k] >= 3) & (cov[k] >= 0.2 * tot)).astype(np.uint8) << k)
    return c, s, m

def votes(X, Y):
    """certified sites of the pair and the read's votes: (sites, for X, for Y); X, Y = (forward query index, forward mask, read bit)"""
    _, ix, iy = np.intersect1d(X[0], Y[0], assume_unique=True, return_indices=True)
    mx = X[1][ix]; my = Y[1][iy]; rb = X[2][ix]
    cert = (mx > 0) & (my > 0) & ((mx & my) == 0)
    return int(cert.sum()), int((cert & ((rb & mx) > 0)).sum()), int((cert & ((rb & my) > 0)).sum())
_a = (np.array([0, 1, 2, 3]), np.array([1, 1, 3, 0], np.uint8), np.array([1, 2, 1, 1], np.uint8))
_b = (np.array([1, 2, 3, 0]), np.array([2, 1, 2, 4], np.uint8), np.array([2, 1, 1, 1], np.uint8))
assert votes(_a, _b) == (2, 1, 1)  # site 0: A vs G, read A -> X; site 1: A vs C, read C -> Y; site 2 shares A; site 3 has no allele at X

def main(tsv, sel, region, bed):
    T = table(tsv)
    regs = [(l.split()[0], int(l.split()[1]), int(l.split()[2])) for l in open(bed)]
    with Pool(int(os.environ.get('OP_THREADS', 32))) as p: piles = p.map(pile, [(region,) + r for r in regs], chunksize=4)
    M = collections.defaultdict(list)
    for c, s, m in piles: M[c].append((s, m))
    for c in M: M[c].sort(key=lambda x: x[0])
    S = {c: [s for s, _ in v] for c, v in M.items()}
    print(f"regions {len(regs)}; positions with an own allele {sum(int((m > 0).sum()) for v in M.values() for _, m in v)}; with two or more {sum(int((np.bitwise_count(m) > 1).sum()) for v in M.values() for _, m in v)}")
    def mask_at(c, pos):
        out = np.zeros(len(pos), dtype=np.uint8)
        if c not in M or not len(pos): return out
        for j in range(bisect.bisect_right(S[c], int(pos[0])) - 1, len(M[c])):
            if j < 0: continue
            s, m = M[c][j]
            if s > pos[-1]: break
            k = (pos >= s) & (pos < s + len(m)); out[k] = m[pos[k] - s]
        return out
    recs = collections.defaultdict(list); seqs = {}
    for r in pysam.AlignmentFile(sel, check_sq=False):
        if r.is_unmapped or r.is_supplementary or r.query_name not in T: continue
        recs[r.query_name].append((r.reference_name, '-' if r.is_reverse else '+', r.reference_start, r.reference_end, r.get_tag('AS'), r.cigartuples))
        if not r.is_secondary: seqs[r.query_name] = r.query_sequence.translate(COMP)[::-1] if r.is_reverse else r.query_sequence
    C = collections.Counter()
    calls = open(os.environ['OP_CALLS'], 'w') if os.environ.get('OP_CALLS') else None
    for q, d in T.items():
        fw = np.frombuffer(seqs[q].encode(), dtype=np.uint8); L = len(fw); best = d['loci'][0][4]
        comp = [i for i, l in enumerate(d['loci']) if l[4] >= best - 10 or i == d['primary']] if d['st'] == 'S2' else [0, 1]
        V = {}
        for i in comp:
            c, s, e, st, a = d['loci'][i]
            rr = [x for x in recs[q] if x[0] == c and x[1] == st and x[2] < e and s < x[3]]
            if not rr: continue
            x = max(rr, key=lambda x: x[4]); qi, rp = aligned(x[2], x[5]); m = mask_at(c, rp)
            if st == '-': qi = L - 1 - qi; m = CMASK[m]; o = np.argsort(qi); qi = qi[o]; m = m[o]
            ok = qi < L; qi = qi[ok]; m = m[ok]
            V[i] = (qi, m, BIT[fw[qi]])
        res = {}
        for name, win in (('own', lambda vx, vy: vx >= 2 and vy == 0), ('own-loose', lambda vx, vy: vx - vy >= 2)):
            w = [i for i in V if len(V) > 1 and all(win(*votes(V[i], V[j])[1:]) for j in V if j != i)]
            res[name] = w[0] if len(w) == 1 else None
        anysite = any(votes(V[i], V[j])[0] for i in V for j in V if i < j)
        if calls: calls.write(f"{q}\t{'NA' if res['own'] is None else res['own']}\t{'NA' if res['own-loose'] is None else res['own-loose']}\t{int(anysite)}\n")
        k = (d['st'], d['status']); C[k + ('sites', anysite)] += 1
        for arm, call in [('primary', d['primary'])] + list(res.items()):
            if call is None: C[k + (arm, 'abstain')] += 1; fin = d['primary']
            else:
                fin = call; C[k + (arm, 'right' if call in d['origin'] else 'wrong')] += 1
                if call != d['primary']: C[k + (arm, 'moved', call in d['origin'])] += 1
            if arm != 'primary': C[k + (arm + '-or-primary', 'right' if fin in d['origin'] else 'wrong')] += 1
    for st in ('S2', 'S1'):
        n = sum(1 for d in T.values() if d['st'] == st and d['status'] == 'decidable')
        print(f"\n{st}: {n} decidable reads, with a certified site between competitors {C[(st, 'decidable', 'sites', True)]}\narm\tright\twrong\tabstain\twrong_share_of_calls\tmoved_right\tmoved_wrong")
        for arm in ('primary', 'own', 'own-or-primary', 'own-loose', 'own-loose-or-primary'):
            r, x, ab = (C[(st, 'decidable', arm, k)] for k in ('right', 'wrong', 'abstain'))
            print(f"{arm}\t{r}\t{x}\t{ab}\t{x / (r + x) if r + x else 0:.4f}\t{C[(st, 'decidable', arm, 'moved', True)]}\t{C[(st, 'decidable', arm, 'moved', False)]}")
        print("class\treads\twith_site\town_calls\town_moves\tloose_calls")
        for cl in ('outside', 'tie', 'no near-perfect hit'):
            n = sum(1 for d in T.values() if d['st'] == st and d['status'] == cl)
            print(f"{cl}\t{n}\t{C[(st, cl, 'sites', True)]}\t{n - C[(st, cl, 'own', 'abstain')]}\t{C[(st, cl, 'own', 'moved', True)] + C[(st, cl, 'own', 'moved', False)]}\t{n - C[(st, cl, 'own-loose', 'abstain')]}")

if __name__ == '__main__':
    main(*sys.argv[1:5])
