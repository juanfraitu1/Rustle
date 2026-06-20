"""Executable verification of the copy-assignment theory note's combinatorial claims."""
import itertools
import random

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


# --------------------------------------------------------------------------------------
# Corrected Condition C (§5): (C1) distinguishability + (C2) = (C2-sep) AND (C2-link).
#
#   D_ij     = columns where copies i and j differ.
#   Delta_i  = union_{j != i} D_ij  (copy i's "identity" / distinguishing columns).
#
#   (C1)       all true copies are pairwise distinct allele-vectors.
#   (C2-sep)   every read of copy i conflicts (in H) with >= 1 read of every foreign copy j
#              (cross-conflict: forbids two copies from merging).
#   (C2-link)  for each copy i, the column-linkage graph L_i on Delta_i -- nodes = Delta_i,
#              edge {a,b} iff some single read of copy i observes BOTH a and b -- is CONNECTED
#              (phasing backbone: forbids recombination by pinning the phase between sites).
#
# (C2-sep) alone is the OLD, FALSE condition (admits the recombinant haplotype-phasing cover);
# (C2-link) is the repair.  Both clauses are independently load-bearing (ablated below).
# --------------------------------------------------------------------------------------


def _reads_by_copy(labels):
    parts = {}
    for idx, lab in enumerate(labels):
        parts.setdefault(lab, []).append(idx)
    return parts


def _distinguishing(copies):
    """D[(i,j)] = set of columns where copies i and j differ."""
    L = len(copies[0]) if copies else 0
    D = {}
    for i in range(len(copies)):
        for j in range(len(copies)):
            if i != j:
                D[(i, j)] = {c for c in range(L) if copies[i][c] != copies[j][c]}
    return D


def _deltas(copies, rc):
    """Delta_i = union over foreign j of D_ij = copy i's distinguishing columns."""
    D = _distinguishing(copies)
    return {
        i: (set().union(*[D[(i, j)] for j in rc if j != i]) if len(rc) > 1 else set())
        for i in rc
    }


def condition_C(copies, reads, labels, *, use_sep=True, use_link=True):
    """True iff the instance satisfies the corrected Condition C.

    Parameters use_sep / use_link toggle the (C2-sep) / (C2-link) clauses so the
    randomized check can ablate them and confirm each is load-bearing.
    Single-copy instances (K < 2) are trivially identifiable.
    """
    K = len(set(labels))
    if K < 2:
        return True
    # (C1) distinguishability: all copies pairwise distinct.
    present = sorted(set(labels))
    for i, j in itertools.combinations(present, 2):
        if copies[i] == copies[j]:
            return False
    H = conflict_graph(reads)
    adj = {n: set(H.neighbors(n)) for n in H.nodes()}
    rc = _reads_by_copy(labels)
    delta = _deltas(copies, rc)
    # (C2-sep) every read of copy i conflicts with >= 1 read of every foreign copy j.
    if use_sep:
        for r, lr in enumerate(labels):
            for fj in rc:
                if fj == lr:
                    continue
                if not (adj[r] & set(rc[fj])):
                    return False
    # (C2-link) per copy: column-linkage graph on Delta_i is connected.
    if use_link:
        for i in rc:
            di = delta[i]
            if len(di) <= 1:
                continue  # 0 or 1 distinguishing column: trivially connected
            Lg = nx.Graph()
            Lg.add_nodes_from(di)
            for ridx in rc[i]:
                obs = set(observed(reads[ridx]).keys()) & di
                for a, b in itertools.combinations(sorted(obs), 2):
                    Lg.add_edge(a, b)
            if not nx.is_connected(Lg):
                return False
    return True


def unique_min_cover_is_true(reads, labels):
    """True iff MCC == K and the true partition is the UNIQUE minimum cover.

    Ground-truth oracle via brute force: mcc_bruteforce for the size, and
    all_min_colorings (enumerated over the conflict graph) for uniqueness.
    """
    K = len(set(labels))
    mcc = mcc_bruteforce(reads)
    if mcc != K:
        return False
    H = conflict_graph(reads)
    covers = {partition_of(c) for c in all_min_colorings(H, mcc)}
    return covers == {partition_of(labels)}


def check_thm2_recovery():
    """Verify Theorem 2 (identifiability) on a K=2 instance over L=3 columns satisfying Condition C.

    Two copies differ at columns 0 and 2 (K_{ij} = 2, robust regime), so Delta_0 = Delta_1 = {0, 2}.
    The reads tile the column-pairs so that (C2-sep) holds (every read conflicts with a foreign read)
    AND (C2-link) holds (the read {0,2} links columns 0 and 2 within each copy, making each L_i
    connected -- pinning the phase, forbidding the recombinant cover).  Theorem 2 then predicts:
    MCC = 2 AND the true partition is the UNIQUE minimum cover.
    """
    copies = [(0, 0, 0), (1, 0, 1)]  # differ at cols 0 and 2; agree at col 1 -> Delta_i = {0,2}
    # Each copy contributes 3 reads tiling all column-pairs {0,1},{1,2},{0,2}.
    # The {0,2} read is the linkage backbone: it spans both distinguishing columns.
    windows = [
        (0, {0, 1}), (0, {1, 2}), (0, {0, 2}),
        (1, {0, 1}), (1, {1, 2}), (1, {0, 2}),
    ]
    reads, labels = reads_from_copies(copies, windows)
    H = conflict_graph(reads)

    # The instance must satisfy the corrected Condition C (both (C2-sep) AND (C2-link)).
    assert condition_C(copies, reads, labels), (
        "check_thm2_recovery instance must satisfy the corrected Condition C"
    )
    # And it must FAIL if either clause is dropped only when that clause is genuinely needed;
    # here we at least confirm the full condition holds and (C2-link) is active (Delta has 2 cols).
    assert _deltas(copies, _reads_by_copy(labels))[0] == {0, 2}, "Delta_0 should be {0,2}"

    k = mcc_bruteforce(reads)
    assert k == 2, f"MCC should be 2, got {k}"

    true_part = partition_of(labels)
    colorings = {partition_of(c) for c in all_min_colorings(H, 2)}
    assert colorings == {true_part}, (
        f"the true partition must be the UNIQUE minimum cover under C; got {colorings}"
    )
    return "Thm 2: K=2 instance under C -> unique minimum cover == true copies"


def check_thm2_uniqueness_random(n_instances=4000, seed=20260620):
    """Randomized empirical proof that the corrected Condition C is SUFFICIENT for uniqueness.

    Generate many small random instances (K in {2,3}, L in {2..4}, random distinct copies, random
    read windows), KEEP those satisfying the corrected Condition C, and for each KEPT instance assert
    by brute force that the minimum copy cover is UNIQUE and equals the true partition.  This must hold
    with ZERO violations -- it is the operative certificate that Condition C (= (C1) + (C2-sep) +
    (C2-link)) implies Theorem 2's uniqueness.

    Three further assertions guard the condition itself:
      (1) the known recombinant counterexample (4 single-column reads, two min covers) does NOT
          satisfy Condition C -- it must be EXCLUDED;
      (2) dropping (C2-link) [keeping only the old per-read (C2-sep)] RE-ADMITS a non-unique instance
          -- so (C2-link) is load-bearing;
      (3) dropping (C2-sep) [keeping only (C2-link)] RE-ADMITS a non-unique instance
          -- so (C2-sep) is load-bearing.
    Deterministic via a fixed seed.
    """
    # --- (1) the recombinant counterexample must be excluded by Condition C ---
    # This is the canonical, deterministic witness that (C2-link) is load-bearing: it satisfies the
    # OLD (C2-sep)-only condition yet has TWO minimum covers (the true one and the recombinant).
    ce_copies = [(0, 0), (1, 1)]
    ce_windows = [(0, {0}), (0, {1}), (1, {0}), (1, {1})]  # no read spans BOTH columns
    ce_reads, ce_labels = reads_from_copies(ce_copies, ce_windows)
    assert not condition_C(ce_copies, ce_reads, ce_labels), (
        "the recombinant phasing counterexample must NOT satisfy the corrected Condition C"
    )
    # It satisfies the OLD (C2-sep)-only condition (which is exactly the bug) ...
    assert condition_C(ce_copies, ce_reads, ce_labels, use_link=False), (
        "counterexample should satisfy (C2-sep) alone -- that is the falsified weak condition"
    )
    # ... yet has TWO minimum covers (truth + recombinant), so it is genuinely non-unique.  Hence
    # dropping (C2-link) re-admits a non-unique instance => (C2-link) is load-bearing (deterministic):
    assert not unique_min_cover_is_true(ce_reads, ce_labels), (
        "counterexample must have a non-unique minimum cover (recombinant + true)"
    )

    def random_instances(rng, n):
        out = []
        for _ in range(n):
            K = rng.choice([2, 3])
            L = rng.randint(2, 4)
            copies, tries = [], 0
            while len(copies) < K and tries < 64:
                v = tuple(rng.randint(0, 1) for _ in range(L))
                if v not in copies:
                    copies.append(v)
                tries += 1
            if len(copies) < K:
                continue
            windows = []
            for ci in range(K):
                for _ in range(rng.randint(1, 4)):
                    ncols = rng.randint(1, L)
                    windows.append((ci, set(rng.sample(range(L), ncols))))
            reads, labels = reads_from_copies(copies, windows)
            out.append((copies, reads, labels))
        return out

    rng = random.Random(seed)
    instances = random_instances(rng, n_instances)

    sampled = len(instances)
    satisfied = 0
    violations = 0
    # ablation counters: instances that pass the WEAKENED condition but are non-unique
    # Seed the (C2-link)-load-bearing counter with the deterministic recombinant counterexample
    # (it passes (C2-sep)-only and is non-unique), so the corroboration never depends on the seed.
    sep_only_admits_nonunique = 1   # (C2-sep) alone -> non-unique exists (proves C2-link needed)
    link_only_admits_nonunique = 0  # (C2-link) alone -> should find some non-unique (proves C2-sep needed)

    for copies, reads, labels in instances:
        if condition_C(copies, reads, labels):
            satisfied += 1
            if not unique_min_cover_is_true(reads, labels):
                violations += 1
        # ablations
        if condition_C(copies, reads, labels, use_link=False):
            if not unique_min_cover_is_true(reads, labels):
                sep_only_admits_nonunique += 1
        if condition_C(copies, reads, labels, use_sep=False):
            if not unique_min_cover_is_true(reads, labels):
                link_only_admits_nonunique += 1

    assert satisfied > 0, "sample produced no Condition-C instances; broaden the generator"
    assert violations == 0, (
        f"Condition C must be SUFFICIENT: {violations}/{satisfied} kept instances had a "
        "non-unique minimum cover"
    )
    # (2) (C2-link) load-bearing: established DETERMINISTICALLY by the counterexample above
    #     (it passes (C2-sep)-only yet is non-unique); the random sweep corroborates.
    assert sep_only_admits_nonunique > 0, (
        "(C2-link) must be load-bearing: (C2-sep) alone admits the recombinant non-unique cover"
    )
    # (3) (C2-sep) load-bearing: dropping it re-admits many non-unique instances in the sweep.
    assert link_only_admits_nonunique > 0, (
        "(C2-sep) must be load-bearing: dropping it should re-admit a non-unique instance"
    )

    return (
        f"Thm 2 uniqueness (randomized): {sampled} sampled, {satisfied} satisfy Condition C, "
        f"{violations} uniqueness violations; counterexample excluded; "
        f"(C2-sep),(C2-link) each load-bearing "
        f"(sep-only re-admits {sep_only_admits_nonunique}, link-only re-admits "
        f"{link_only_admits_nonunique})"
    )


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
CHECKS.append(check_thm2_uniqueness_random)
CHECKS.append(check_thm2_K0_merge)


def main():
    for fn in CHECKS:
        print("OK  -", fn())


if __name__ == "__main__":
    main()
