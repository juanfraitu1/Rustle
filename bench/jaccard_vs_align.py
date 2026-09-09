#!/usr/bin/env python3
"""Time MinHash-Jaccard against the pairwise DNA alignment O1 actually uses, on ONE input.

Advisor question: "You perform pairwise alignments of DNA sequences. Can you report runtimes of
this step? How does it compare to just using Jaccard indexes?"

Both methods see the same FASTA. minimap2's runtime is read from its own log (or re-timed);
the sketching and all-pairs stages here are timed separately, because a Jaccard pipeline is
sketch-once + compare-all-pairs and the two scale differently.

The comparison that matters is not speed alone but whether Jaccard can express O1's edge rule:
identity >= 0.70 AND coverage of the longer side >= 0.30 AND >= 300 shared bp.

usage: jaccard_vs_align.py <spans.fa> <aln.paf> [--k 21] [--sketch 1000]
"""
import sys, time, collections
import numpy as np
from scipy.sparse import csr_matrix

fa_p, paf_p = sys.argv[1], sys.argv[2]
def opt(k, d): return type(d)(sys.argv[sys.argv.index(k)+1]) if k in sys.argv else d
K = opt('--k', 21); S = opt('--sketch', 1000)

# ---------------------------------------------------------------- load
t0 = time.time()
names, seqs, cur, buf = [], [], None, []
for line in open(fa_p):
    if line[0] == '>':
        if cur is not None: seqs.append(''.join(buf))
        cur = line[1:].strip().split()[0]; names.append(cur); buf = []
    else: buf.append(line.strip())
seqs.append(''.join(buf))
t_load = time.time() - t0
print(f'loaded {len(names)} sequences, {sum(map(len,seqs)):,} bp in {t_load:.1f}s')

# ---------------------------------------------------------------- sketch
# Vectorised 2-bit k-mer sketching. A per-character Python loop over ~97 Mb is hours; this is
# seconds and computes the SAME bottom-S sketch. (Runtime caveat: a C implementation such as
# mash would still beat this — see the note in the write-up. The structural conclusion below
# does not depend on which implementation is timed.)
LUT = np.full(256, 255, dtype=np.uint8)
for i, c in enumerate('ACGT'):
    LUT[ord(c)] = i; LUT[ord(c.lower())] = i
MUL = np.uint64(0x9E3779B97F4A7C15)

def sketch(s):
    a = LUT[np.frombuffer(s.encode(), dtype=np.uint8)]
    n = a.size
    if n < K: return np.empty(0, dtype=np.uint64)
    ok = a != 255
    a = np.where(ok, a, 0).astype(np.uint64)
    # rolling forward and reverse-complement k-mers, k<=31 so 2k bits fit in uint64
    fwd = np.zeros(n - K + 1, dtype=np.uint64)
    rev = np.zeros(n - K + 1, dtype=np.uint64)
    two = np.uint64(2)
    for i in range(K):
        fwd = (fwd << two) | a[i:n - K + 1 + i]
        rev |= (np.uint64(3) - a[i:n - K + 1 + i]) << np.uint64(2 * i)
    # drop windows containing a non-ACGT base
    valid = np.convolve(ok.astype(np.int32), np.ones(K, dtype=np.int32), 'valid') == K
    canon = np.minimum(fwd, rev)[valid]
    if canon.size == 0: return np.empty(0, dtype=np.uint64)
    # splitmix-style finaliser so the bottom-S sketch is a uniform sample
    x = canon * MUL
    x ^= x >> np.uint64(29); x *= np.uint64(0xBF58476D1CE4E5B9); x ^= x >> np.uint64(32)
    u = np.unique(x)
    return u[:S]

t0 = time.time()
sk = [sketch(s) for s in seqs]
t_sketch = time.time() - t0
print(f'sketched  k={K} s={S}  in {t_sketch:.1f}s')

# ---------------------------------------------------------------- all pairs
t0 = time.time()
vocab = {}
rows, cols = [], []
for i, a in enumerate(sk):
    for v in a:
        j = vocab.setdefault(int(v), len(vocab))
        rows.append(i); cols.append(j)
M = csr_matrix((np.ones(len(rows), dtype=np.float32), (rows, cols)),
               shape=(len(sk), len(vocab)))
inter = (M @ M.T).tocoo()                            # shared sketch hashes per pair
sizes = np.array([len(a) for a in sk], dtype=np.float32)
t_pairs = time.time() - t0
print(f'all-pairs {len(sk)*(len(sk)-1)//2:,} comparisons in {t_pairs:.1f}s')
print(f'JACCARD TOTAL (sketch + all-pairs) = {t_sketch + t_pairs:.1f}s')

# Jaccard estimate for every pair the sketch says shares anything
jac = {}
for i, j, v in zip(inter.row, inter.col, inter.data):
    if i >= j: continue
    u = sizes[i] + sizes[j] - v
    if u > 0: jac[(int(i), int(j))] = v / u

# ---------------------------------------------------------------- O1's edges from the PAF
def base(n): return n
idx = {n: i for i, n in enumerate(names)}
edges = {}
LEN = {}
for line in open(paf_p):
    f = line.rstrip('\n').split('\t')
    q, ql, qs, qe, tn, tl, ts, te = f[0], int(f[1]), int(f[2]), int(f[3]), f[5], int(f[6]), int(f[7]), int(f[8])
    if q == tn: continue
    nm = int(f[9]); bl = int(f[10])
    if q not in idx or tn not in idx: continue
    a, b = idx[q], idx[tn]
    if a > b: a, b = b, a
    LEN[a] = ql if idx[q] == a else tl
    e = edges.setdefault((a, b), [0, 0, 0.0, max(ql, tl)])
    e[0] += qe - qs; e[1] += te - ts
    e[2] = max(e[2], nm / bl if bl else 0.0)
o1 = set()
for (a, b), (qcov, tcov, ident, longer) in edges.items():
    shared = max(qcov, tcov)
    if ident >= 0.70 and shared >= 300 and shared / longer >= 0.30:
        o1.add((a, b))
print(f'\nO1 edges from the alignment: {len(o1):,}   (identity>=0.70, cov_longer>=0.30, >=300bp)')

# ---------------------------------------------------------------- can Jaccard reproduce them?
jv = np.array([jac.get(e, 0.0) for e in o1])
print(f'  Jaccard on those SAME pairs: median {np.median(jv):.4f}  '
      f'q10 {np.quantile(jv,0.1):.4f}  q90 {np.quantile(jv,0.9):.4f}')
print(f'  O1 edges with Jaccard < 0.01: {int((jv < 0.01).sum()):,} '
      f'({(jv < 0.01).mean():.1%})  <- invisible to any Jaccard threshold')
allj = np.array(list(jac.values()))
print(f'  all sketched pairs with Jaccard >= 0.01: {int((allj>=0.01).sum()):,}')
for thr in (0.001, 0.005, 0.01, 0.05, 0.10):
    tp = int((jv >= thr).sum()); fp = int((allj >= thr).sum()) - tp
    print(f'   thr {thr:<6} recall of O1 edges {tp/len(o1):6.1%}   extra pairs above it {fp:,}')
