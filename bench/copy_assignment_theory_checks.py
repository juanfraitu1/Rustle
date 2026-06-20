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


CHECKS = [check_lemma_mcc_equals_chromatic]


def main():
    for fn in CHECKS:
        print("OK  -", fn())


if __name__ == "__main__":
    main()
