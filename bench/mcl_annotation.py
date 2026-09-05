#!/usr/bin/env python3
"""Define gene families by MCL clustering of the ANNOTATION, by SEQUENCE — never by gene symbol.

⭐ WHY MCL. `docs/NEGATIVE_RESULTS_REGISTER.md:474` already diagnosed the gap: transitive closure over the
pairwise homology graph gives superfamilies of 145 and 114 genes because it chains subfamilies through
domain hubs — "It is OrthoFinder without MCL, and skipping MCL is why."

⭐ WHY NO SYMBOLS. RFPL is 9/9 LOC-named and GOLGA6 14/14, so a `gene=NPIP*` grep is not a definition
(project-gene-naming-traps). FASTA headers here are bare coordinates; symbols are attached only when
REPORTING, never when clustering.

⚠ EDGE WEIGHT USES cov_longer, NOT cov_shorter. A ~300 bp Alu covers most of a short fragment and almost
none of a real gene, so a shorter-side weight lets shared repeat drive clustering (ledger §6cs/§6cz).
Weighting by the LONGER side asks how much of BOTH objects is shared, which is what paralogy means.
"""
import sys, collections, numpy as np
from scipy import sparse

def parse_paf(path, min_id=0.70, min_covlong=0.30, min_bp=300, cache=True, exonic=None):
    """One weight per unordered pair: the best identity x cov_longer over its records.

    ⭐ CACHED. The all-vs-all PAF is the expensive artifact (genome-wide: hours); parsing it into the
    weighted graph is the only step between it and MCL, and the graph does not change when you sweep
    inflation. The cache key is (path, size, mtime, thresholds) so an edited or regenerated PAF, or a
    different threshold, misses cleanly rather than silently serving a stale graph.
    """
    import hashlib, os, pickle
    key = None
    if cache:
        st = os.stat(path)
        key = hashlib.md5(
            f"{os.path.abspath(path)}|{st.st_size}|{int(st.st_mtime)}|{min_id}|{min_covlong}|{min_bp}|{'exonic' if exonic else 'span'}".encode()
        ).hexdigest()[:16]
        cf = os.path.join(os.path.dirname(os.path.abspath(path)), f".mclcache_{key}.pkl")
        if os.path.exists(cf):
            with open(cf, "rb") as fh:
                names, best = pickle.load(fh)
            print(f"[cache HIT] {cf}  {len(names)} nodes, {len(best)} edges", file=sys.stderr)
            return names, best
    names, best = _parse_paf_uncached(path, min_id, min_covlong, min_bp, exonic)
    if cache:
        with open(cf, "wb") as fh:
            pickle.dump((names, best), fh)
        print(f"[cache WRITE] {cf}", file=sys.stderr)
    return names, best


def _hkey(h):
    """FASTA header CONTIG:START-END -> the (chrom, start, end) key of the exonic-length map."""
    try:
        c, r = h.split(":"); a, b = r.split("-")
        return (c, int(a), int(b))
    except Exception:
        return None


def _parse_paf_uncached(path, min_id, min_covlong, min_bp, exonic=None):
    """`exonic`: optional {(chrom,start,end): exonic_bp}. ⭐ THE DENOMINATOR MATTERS ENORMOUSLY.
    The FASTA sequences are GENOMIC SPANS, and the median gorilla gene is only 23.15% exon, so a
    cov_longer floor of 0.30 measured against the SPAN demands more coverage than the median gene has
    exonic sequence at all. Measured: it removed 1,525/2,428 = 62.8% of genes that already passed
    identity and length, and switching to exonic bases took the graph 903 -> 1,850 nodes with 0 lost
    (ledger §6dc). Same class as the scale-free-denominator defects in §6co and the single-exon stubs."""
    best = {}
    names = {}
    def idx(n):
        if n not in names: names[n] = len(names)
        return names[n]
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12: continue
            q, ql, qs, qe, t, tl, ts, te = f[0], int(f[1]), int(f[2]), int(f[3]), f[5], int(f[6]), int(f[7]), int(f[8])
            if q == t: continue
            nm, bl = int(f[9]), int(f[10])
            if bl < min_bp: continue
            ident = nm / max(1, bl)
            if ident < min_id: continue
            aln = max(qe - qs, te - ts)
            if exonic:
                dq, dt = exonic.get(_hkey(q)), exonic.get(_hkey(t))
                den = max(dq or ql, dt or tl)
            else:
                den = max(ql, tl)
            cov_long = min(1.0, aln / max(1, den))
            if cov_long < min_covlong: continue
            w = ident * cov_long
            a, b = idx(q), idx(t)
            k = (a, b) if a < b else (b, a)
            if w > best.get(k, 0.0): best[k] = w
    return names, best

def mcl(M, inflation=2.0, max_iter=100, prune=1e-5, tol=1e-7):
    """Markov clustering: add self-loops, then alternate expansion (M@M) and inflation (elementwise
    power then column-renormalise) until the matrix stops changing."""
    n = M.shape[0]
    M = M + sparse.identity(n, format="csr")          # self-loops
    M = normalize_cols(M)
    for _ in range(max_iter):
        prev = M
        M = (M @ M).tocsr()                            # expand
        M.data = np.power(M.data, inflation)           # inflate
        M.data[M.data < prune] = 0.0
        M.eliminate_zeros()
        M = normalize_cols(M)
        if M.shape == prev.shape and abs(M - prev).max() < tol: break
    return M

def normalize_cols(M):
    s = np.asarray(M.sum(axis=0)).ravel()
    s[s == 0] = 1.0
    return M @ sparse.diags(1.0 / s)

def clusters_from(M, names_inv):
    """Attractors are rows with a non-zero diagonal-ish mass; group columns by their attractor row."""
    M = M.tocsc()
    lab = {}
    for j in range(M.shape[1]):
        s, e = M.indptr[j], M.indptr[j+1]
        if e == s: continue
        rows, vals = M.indices[s:e], M.data[s:e]
        lab[j] = int(rows[np.argmax(vals)])
    # union attractors that co-occur
    parent = {}
    def find(x):
        while parent.get(x, x) != x: parent[x] = parent.get(parent[x], parent[x]); x = parent[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb: parent[ra] = rb
    for j, r in lab.items(): union(j, r)
    out = collections.defaultdict(list)
    for j in lab: out[find(j)].append(names_inv[j])
    return list(out.values())

if __name__ == "__main__":
    paf = sys.argv[1]
    infl = float(sys.argv[2]) if len(sys.argv) > 2 else 2.0
    names, best = parse_paf(paf)
    inv = {v: k for k, v in names.items()}
    n = len(names)
    print(f"graph: {n} nodes, {len(best)} weighted edges (identity>=0.70, cov_longer>=0.30, >=300bp)", file=sys.stderr)
    if not best: sys.exit("no edges")
    r = np.array([k[0] for k in best] + [k[1] for k in best])
    c = np.array([k[1] for k in best] + [k[0] for k in best])
    v = np.array(list(best.values()) * 2)
    M = sparse.csr_matrix((v, (r, c)), shape=(n, n))
    cl = clusters_from(mcl(M, inflation=infl), inv)
    cl = [c for c in cl if len(c) >= 2]
    cl.sort(key=len, reverse=True)
    print(f"inflation={infl}  multi-member clusters={len(cl)}  "
          f"largest={len(cl[0]) if cl else 0}  members={sum(len(c) for c in cl)}", file=sys.stderr)
    for i, c in enumerate(cl):
        for m in sorted(c): print(f"CL{i}\t{len(c)}\t{m}")
