"""Executable verification of the copy-assignment theory note's combinatorial claims."""
import itertools
import networkx as nx

# A read is a frozenset of (column, allele) pairs.
# Helpers operate on lists of reads.


def observed(read):
    """Return a dict {column: allele} for a read."""
    return {j: a for (j, a) in read}


def conflict_graph(reads):
    """Build the conflict graph H: vertices = reads, edge iff they disagree at a co-observed column."""
    G = nx.Graph()
    G.add_nodes_from(range(len(reads)))
    for i, j in itertools.combinations(range(len(reads)), 2):
        oi, oj = observed(reads[i]), observed(reads[j])
        if any(oi[c] != oj[c] for c in (oi.keys() & oj.keys())):
            G.add_edge(i, j)
    return G


def jointly_consistent(read_idxs, reads):
    """True iff the reads share at least one allele-vector (no column observed with two alleles)."""
    col = {}
    for k in read_idxs:
        for c, a in observed(reads[k]).items():
            if col.setdefault(c, a) != a:
                return False
    return True


def mcc_bruteforce(reads):
    """Minimum Copy Cover size: smallest k s.t. reads partition into k jointly-consistent parts.

    Uses backtracking; exponential — intended for tiny instances only.
    """
    n = len(reads)
    idxs = list(range(n))
    for k in range(1, n + 1):
        if _can_partition(idxs, k, reads):
            return k
    return n


def _can_partition(idxs, k, reads):
    """Backtracking: assign each read a label in 0..k-1 such that every class is jointly consistent."""
    label = {}

    def bt(pos):
        if pos == len(idxs):
            return True
        for lab in range(k):
            label[idxs[pos]] = lab
            cls = [i for i in idxs[: pos + 1] if label[i] == lab]
            if jointly_consistent(cls, reads) and bt(pos + 1):
                return True
        del label[idxs[pos]]
        return False

    return bt(0)


def check_lemma_mcc_equals_chromatic():
    """Verify Lemma 1 (MCC = chi(H)) on the C5 gadget.

    Construction: 5 reads indexed 0..4 whose conflict graph is a 5-cycle C5.
    Each read i observes two columns (the edges incident to vertex i in C5)
    with opposite alleles, forcing a pairwise conflict between every adjacent pair.

    C5 is an odd cycle: chi(C5) = 3.  Lemma 1 predicts MCC = 3.
    DSATUR (greedy saturation-largest-first) is exact on C5.
    """
    edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]
    cols = {e: idx for idx, e in enumerate(edges)}

    reads = []
    for w in range(5):
        obs = []
        for e in edges:
            if w in e:
                obs.append((cols[e], 0 if w == e[0] else 1))
        reads.append(frozenset(obs))

    # Build conflict graph and verify it is exactly C5.
    H = conflict_graph(reads)
    assert H.number_of_nodes() == 5, "Expected 5 nodes"
    assert H.number_of_edges() == 5, f"Expected 5 edges (C5), got {H.number_of_edges()}"
    expected_edges = sorted((min(i, (i + 1) % 5), max(i, (i + 1) % 5)) for i in range(5))
    assert sorted(H.edges()) == expected_edges, (
        f"Conflict graph is not C5: {sorted(H.edges())}"
    )

    # Chromatic number via DSATUR (exact on C5: odd cycle -> chi = 3).
    chi = 1 + max(nx.algorithms.coloring.greedy_color(H, strategy="saturation_largest_first").values())
    assert chi == 3, f"Expected chi(C5) = 3, got {chi}"

    # Lemma 1: MCC must equal chi(H) = 3.
    mcc = mcc_bruteforce(reads)
    assert mcc == 3, f"MCC must equal chromatic number (3) on the C5 gadget, got {mcc}"

    return "Lemma 1 (MCC=chi) verified on C5: MCC=3=chi(C5)"


def reduction_instance(graph_edges, n):
    """Map graph Gamma=(V=range(n), edges) to an MCC instance: one binary column per edge, one read per vertex.

    Read w observes column e iff w in e, allele 0 if w==e[0] else 1.
    The conflict graph of the resulting instance is isomorphic to Gamma.
    """
    cols = {e: idx for idx, e in enumerate(graph_edges)}
    reads = []
    for w in range(n):
        obs = [(cols[e], 0 if w == e[0] else 1) for e in graph_edges if w in e]
        reads.append(frozenset(obs))
    return reads


def chromatic_number(G):
    """Brute-force chromatic number of a small graph."""
    n = G.number_of_nodes()
    nodes = list(G.nodes())
    for k in range(1, n + 1):
        for assign in itertools.product(range(k), repeat=n):
            coloring = {nodes[i]: assign[i] for i in range(n)}
            if all(coloring[u] != coloring[v] for u, v in G.edges()):
                return k
    return n


def check_thm1_reduction():
    """Verify Theorem 1's reduction on 4 small graphs.

    For each source graph Gamma, the reduction produces an MCC instance whose conflict graph
    is isomorphic to Gamma, and MCC == chi(Gamma).  Graphs tested:
      - triangle (K3):   chi = 3
      - 4-cycle  (C4):   chi = 2
      - 5-cycle  (C5):   chi = 3
      - star     (K1,3): chi = 2
    """
    cases = [
        ([(0, 1), (1, 2), (2, 0)], 3),                          # triangle  -> chi 3
        ([(0, 1), (1, 2), (2, 3), (3, 0)], 4),                  # C4        -> chi 2
        ([(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)], 5),          # C5        -> chi 3
        ([(0, 1), (0, 2), (0, 3)], 4),                          # star      -> chi 2
    ]
    for edges, n in cases:
        reads = reduction_instance(edges, n)
        H = conflict_graph(reads)
        Gamma = nx.Graph()
        Gamma.add_nodes_from(range(n))
        Gamma.add_edges_from(edges)
        assert nx.is_isomorphic(H, Gamma), (
            f"conflict graph not isomorphic to source graph for edges={edges}"
        )
        chi = chromatic_number(Gamma)
        mcc = mcc_bruteforce(reads)
        assert mcc == chi, (
            f"MCC ({mcc}) != chi ({chi}) for edges={edges}"
        )
    return "Thm 1: conflict graph == source graph and MCC == chi on 4 instances"


def reads_from_copies(copies, windows):
    """copies: list of allele-vectors (tuples in {0,1}^L). windows: list of (copy_index, set-of-columns)
    describing each read's origin copy and which columns it observes. Returns reads + true labels."""
    reads, labels = [], []
    for ci, cols in windows:
        reads.append(frozenset((j, copies[ci][j]) for j in cols))
        labels.append(ci)
    return reads, labels


def all_min_colorings(H, k):
    """All proper k-colorings of H (raw label tuples) -- tiny graphs only (k^n enumeration)."""
    n = H.number_of_nodes()
    out = []
    for assign in itertools.product(range(k), repeat=n):
        if all(assign[u] != assign[v] for u, v in H.edges()):
            out.append(assign)
    return out


def partition_of(labels):
    """Canonical partition (set of frozensets of read-indices) from a labeling."""
    parts = {}
    for i, l in enumerate(labels):
        parts.setdefault(l, set()).add(i)
    return frozenset(frozenset(p) for p in parts.values())


def check_thm2_recovery():
    """Verify Theorem 2 (identifiability) on a K=2 instance over L=3 columns satisfying condition C.

    Two copies differ at columns 0 and 2 (K_{ij} = 2, robust regime).  The reads tile the columns so
    that every read conflicts (in H) with at least one read of the foreign copy -- the per-read
    condition C2.  Theorem 2 then predicts: MCC = 2 AND the true partition is the UNIQUE minimum cover.
    """
    copies = [(0, 0, 0), (1, 0, 1)]  # differ at cols 0 and 2; agree at col 1
    # Each copy contributes 3 reads tiling all column-pairs {0,1},{1,2},{0,2}.
    windows = [
        (0, {0, 1}), (0, {1, 2}), (0, {0, 2}),
        (1, {0, 1}), (1, {1, 2}), (1, {0, 2}),
    ]
    reads, labels = reads_from_copies(copies, windows)
    H = conflict_graph(reads)

    # Per-read condition C2: every read conflicts with >=1 read of every foreign copy.
    adj = {node: set(H.neighbors(node)) for node in H.nodes()}
    for i, li in enumerate(labels):
        for fj in set(labels) - {li}:
            foreign = {k for k, lk in enumerate(labels) if lk == fj}
            assert adj[i] & foreign, f"C2 violated: read {i} (copy {li}) has no conflict with copy {fj}"

    k = mcc_bruteforce(reads)
    assert k == 2, f"MCC should be 2, got {k}"

    true_part = partition_of(labels)
    colorings = {partition_of(c) for c in all_min_colorings(H, 2)}
    assert colorings == {true_part}, (
        f"the true partition must be the UNIQUE minimum cover under C; got {colorings}"
    )
    return "Thm 2: K=2 instance under C -> unique minimum cover == true copies"


def check_thm2_K0_merge():
    """Verify the K=0 boundary of the K-bound corollary: identical copies are non-identifiable.

    When two copies are identical over every observed column (K_{ij} = 0), condition C1 fails: the
    reads produce no conflict edge, the minimum cover merges them into one part, and the true copies
    are provably unrecoverable (the MAGEA co-located regime, resolvable fraction 0/494).
    """
    copies = [(0, 1, 0), (0, 1, 0)]  # identical over all columns -> K_{ij} = 0
    windows = [(0, {0, 1}), (0, {1, 2}), (1, {0, 1}), (1, {1, 2})]
    reads, _ = reads_from_copies(copies, windows)
    H = conflict_graph(reads)
    assert H.number_of_edges() == 0, "identical copies produce no conflict edges"
    assert mcc_bruteforce(reads) == 1, (
        "K=0: minimum cover merges all reads into ONE copy (true copies unrecoverable)"
    )
    return "Thm 2 boundary: K=0 -> minimum cover = 1 (forced merge, non-identifiable)"


CHECKS = [check_lemma_mcc_equals_chromatic, check_thm1_reduction]
CHECKS.append(check_thm2_recovery)
CHECKS.append(check_thm2_K0_merge)


def main():
    for fn in CHECKS:
        print("OK  -", fn())


if __name__ == "__main__":
    main()
