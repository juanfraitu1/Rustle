#!/usr/bin/env python3
"""IsoCon-style partition of a cluster into competing candidate transcripts (docs/PREREG_unmapped_rescue_2026-10-08.md Amendment 18).

significant_variants -> read_carry -> find_blocks (variants that travel together across reads) -> choose_split -> recurse on both sides.
partition(reads, align_fn, consensus_fn): align_fn(cons, reads) -> SAM lines of the reads on cons; consensus_fn([seq]) -> str."""
import os
import subprocess

import numpy as np
from scipy.stats import hypergeom

import polish as P

END_WINDOW = 15
MIN_LEAF = 3
MIN_SPLIT = 6
MAX_DEPTH = 6
ALPHA = 0.05


def significant_variants(cons, sam_lines, end=END_WINDOW, n_as_del=False):
    """[(kind, pos, allele)]: the significant alleles of Amendment 9's test, those the consensus lacks and those it carries as minority, without the first and last
    `end` consensus columns"""
    applied, minority, _e = P.corrections(cons, P.pileup(sam_lines, cons, n_as_del))
    L = len(cons)
    V = {(v[0], v[1], v[2]) for v in applied} | {(m[0], m[1], m[2]) for m in minority}
    return sorted((v for v in V if end <= v[1] < L - end), key=lambda v: (v[1], v[0], str(v[2])))


def read_carry(sam_lines, cons, V, n_as_del=False):
    """-> (read names, C, K): C[i, r] = read r carries variant i, K[i, r] = read r covers the variant's column (or gap)"""
    cols = {v[1] for v in V if v[0] in ("sub", "del")}
    names, per = [], []
    for ln in sam_lines:
        if ln[0] == "@":
            continue
        f = ln.rstrip("\n").split("\t")
        if int(f[1]) & 2308 or f[2] == "*":
            continue
        r, q, seq = int(f[3]) - 1, 0, f[9]
        start, allele, ins = r, {}, {}
        for n, op in P.CIG.findall(f[5]):
            n = int(n)
            if op in "=XM":
                for _ in range(n):
                    if r in cols:
                        allele[r] = cons[r] if op == "=" else seq[q]
                    r += 1
                    q += 1
            elif op == "I":
                ins[r] = seq[q:q + n]
                q += n
            elif op == "D" or (op == "N" and n_as_del):
                for _ in range(n):
                    if r in cols:
                        allele[r] = "-"
                    r += 1
            elif op == "N":
                r += n
            elif op == "S":
                q += n
        names.append(f[0])
        per.append((start, r, allele, ins))
    C = np.zeros((len(V), len(names)), dtype=bool)
    K = np.zeros_like(C)
    for j, (s, e, allele, ins) in enumerate(per):
        for i, (kind, pos, x) in enumerate(V):
            if kind == "ins":
                if s < pos < e:
                    K[i, j] = True
                    C[i, j] = ins.get(pos) == x
            elif s <= pos < e:
                K[i, j] = True
                C[i, j] = allele.get(pos) == ("-" if kind == "del" else x)
    return names, C, K


def find_blocks(C, K, cols, alpha=ALPHA, min_carry=MIN_LEAF):
    """connected components (>= 2 variants) of the graph whose edges are significantly positively associated variants of different columns: one-sided
    hypergeometric over the reads covering both, p < alpha / (number of pairs), both carried together by >= min_carry reads"""
    m = C.shape[0]
    if m < 2:
        return []
    Cf, Kf = C.astype(float), K.astype(float)
    a, N, K1, K2 = Cf @ Cf.T, Kf @ Kf.T, Cf @ Kf.T, Kf @ Cf.T
    iu = np.triu_indices(m, 1)
    a_, N_, K1_, K2_ = (np.rint(x[iu]).astype(int) for x in (a, N, K1, K2))
    diff = np.array(cols)[iu[0]] != np.array(cols)[iu[1]]
    ok = diff & (N_ > 0) & (a_ >= min_carry)
    thr = alpha / max(1, int(diff.sum()))
    p = np.ones(len(a_))
    p[ok] = hypergeom.sf(a_[ok] - 1, N_[ok], K1_[ok], K2_[ok])
    parent = list(range(m))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for i, j, pv, good in zip(iu[0], iu[1], p, ok):
        if good and pv < thr:
            parent[find(i)] = find(j)
    comp = {}
    for i in range(m):
        comp.setdefault(find(i), []).append(i)
    return [b for b in comp.values() if len(b) >= 2]


def choose_split(blocks, C, K, min_leaf=MIN_LEAF):
    """-> (smaller side, carriers boolean vector) of the block whose smaller side is largest (both sides >= min_leaf), or None"""
    n = C.shape[1]
    best = None
    for b in blocks:
        covered, carried = K[b].sum(0), C[b].sum(0)
        carriers = (covered > 0) & (2 * carried >= covered)
        small = min(int(carriers.sum()), n - int(carriers.sum()))
        if small >= min_leaf and (best is None or small > best[0]):
            best = (small, carriers)
    return best


def partition(reads, align_fn, consensus_fn, cons=None, depth=0, max_depth=MAX_DEPTH, min_leaf=MIN_LEAF, min_split=MIN_SPLIT, end=END_WINDOW, n_as_del=False):
    """-> [dict(reads=[names], cons=str)]: the leaves, i.e. the candidate transcripts and the reads that support them"""
    names = list(reads)
    if cons is None:
        cons = consensus_fn([reads[n] for n in names])
    leaf = [dict(reads=names, cons=cons)]
    if len(names) < min_split or depth >= max_depth:
        return leaf
    sam = align_fn(cons, reads)
    V = significant_variants(cons, sam, end, n_as_del)
    if len(V) < 2:
        return leaf
    nm, C, K = read_carry(sam, cons, V, n_as_del)
    blocks = find_blocks(C, K, [v[1] for v in V])
    best = choose_split(blocks, C, K, min_leaf) if blocks else None
    if best is None:
        return leaf
    side1 = {nm[i] for i in range(len(nm)) if best[1][i]}
    out = []
    for group in ([n for n in names if n in side1], [n for n in names if n not in side1]):
        out += partition({n: reads[n] for n in group}, align_fn, consensus_fn, None, depth + 1, max_depth, min_leaf, min_split, end, n_as_del)
    return out


def abpoa_consensus(seqs, k=100):
    """the frozen consensus: abPOA (heaviest bundle) of the k longest reads; needs pyabpoa (miniforge python)"""
    import pyabpoa
    pick = sorted(seqs, key=lambda s: -len(s))[:k]
    res = pyabpoa.msa_aligner(aln_mode="g", is_aa=False, cons_algrm="HB").msa(pick, out_cons=True, out_msa=False)
    return res.cons_seq[0] if res.cons_seq else pick[0]


def minimap_align_fn(workdir, threads=2, preset="map-hifi"):
    os.makedirs(workdir, exist_ok=True)

    def align(cons, reads):
        rf, cf = f"{workdir}/reads.fa", f"{workdir}/cons.fa"
        with open(rf, "w") as o:
            for n, s in reads.items():
                o.write(f">{n}\n{s}\n")
        open(cf, "w").write(f">c\n{cons}\n")
        return subprocess.run(f"minimap2 -ax {preset} --eqx -t {threads} {cf} {rf}", shell=True, capture_output=True, text=True, check=True).stdout.splitlines()
    return align


def edlib_align_fn():
    """Amendment 19: semi-global edit-distance alignment of each read inside the candidate (edlib HW, no distance cap), as SAM lines with an extended CIGAR. Keeps a
    long deletion (a skipped exon) as a deletion run, which a seed-chain-extend aligner turns into a soft clip. Needs edlib (miniforge python)."""
    import edlib

    def align(cons, reads):
        out = []
        for name, seq in reads.items():
            r = edlib.align(seq, cons, mode="HW", task="path", k=-1)
            if not r["locations"]:
                continue
            start = r["locations"][0][0]
            out.append("\t".join([name, "0", "cons", str(start + 1), "60", r["cigar"], "*", "0", "0", seq, "*"]))
        return out
    return align
