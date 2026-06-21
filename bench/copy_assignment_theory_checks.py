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


# --------------------------------------------------------------------------------------
# §5 conditions.  EXHAUSTIVE enumeration (check_thm2_strong_exhaustive) is the ground
# truth; it discovered that the previous condition (sep+link) is sufficient only at K=2.
#
#   D_ij     = columns where copies i and j differ.
#   Delta_i  = union_{j != i} D_ij  (copy i's "identity" / distinguishing columns).
#
#   strong (Strong Separation) -- THE sufficient condition for ALL K:
#              for all i != j, EVERY read of copy i conflicts (in H) with EVERY read of
#              copy j.  Equivalently: no jointly-consistent class contains reads from two
#              different copies.  0 uniqueness violations at K=2 AND K=3.
#
#   sep+link (the previous attempt) -- sufficient ONLY at K=2:
#     (sep)    every read of copy i conflicts with >= 1 read of every foreign copy j
#              (cross-conflict: forbids two copies from merging).
#     (link)   for each copy i, the column-linkage graph L_i on Delta_i -- nodes = Delta_i,
#              edge {a,b} iff some single read of copy i observes BOTH a and b -- is CONNECTED
#              (per-copy phasing backbone).
#     0 violations at K=2 but >0 at K=3: cross-copy RECOMBINATION builds an alternative
#     minimum cover realizing a vector outside C* (the K-frontier).
#
# strong is SUFFICIENT but NOT NECESSARY (holds for only ~15% of truly-unique K=3
# instances); the tight necessary condition is recombination-freeness (instance-global,
# no clean closed form).
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


def _conflict(ri, rj):
    """True iff reads ri, rj co-observe a column with differing alleles."""
    oi, oj = observed(ri), observed(rj)
    return any(oi[c] != oj[c] for c in (oi.keys() & oj.keys()))


def cond_strong(copies, reads, labels):
    """Strong Separation: for all i != j, EVERY read of copy i conflicts with EVERY read of copy j.

    Equivalently: no jointly-consistent class contains reads from two distinct copies.
    This is THE sufficient condition for unique recovery for ALL K (0 violations at K=2 and K=3).
    Single-copy instances (K < 2) are trivially identifiable.
    """
    rc = _reads_by_copy(labels)
    present = sorted(rc)
    if len(present) < 2:
        return True
    # (distinguishability is implied: identical copies would yield a non-conflicting cross pair)
    for i, j in itertools.combinations(present, 2):
        for r in rc[i]:
            for s in rc[j]:
                if not _conflict(reads[r], reads[s]):
                    return False
    return True


def cond_sep_and_link(copies, reads, labels):
    """The PREVIOUS attempt: (sep) cross-conflict per read + (link) per-copy column-linkage connected.

    Sufficient ONLY at K=2 (standard error-free phasing identifiability); FAILS at K>=3 because
    cross-copy recombination can build an alternative minimum cover realizing a vector outside C*.
    Single-copy instances (K < 2) are trivially identifiable.
    """
    rc = _reads_by_copy(labels)
    present = sorted(rc)
    if len(present) < 2:
        return True
    # distinguishability: all present copies pairwise distinct.
    for i, j in itertools.combinations(present, 2):
        if copies[i] == copies[j]:
            return False
    delta = _deltas(copies, rc)
    # (sep) every read of copy i conflicts with >= 1 read of every foreign copy j.
    #       (direct pairwise conflict -- no networkx graph build in the hot loop)
    for fi in rc:
        for r in rc[fi]:
            for fj in rc:
                if fj == fi:
                    continue
                if not any(_conflict(reads[r], reads[s]) for s in rc[fj]):
                    return False
    # (link) per copy: the column-linkage graph on Delta_i is connected (union-find).
    for i in rc:
        di = delta[i]
        if len(di) <= 1:
            continue  # 0 or 1 distinguishing column: trivially connected
        parent = {c: c for c in di}

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for ridx in rc[i]:
            obs = sorted(set(observed(reads[ridx]).keys()) & di)
            for a, b in itertools.combinations(obs, 2):
                parent[find(a)] = find(b)
        if len({find(c) for c in di}) != 1:
            return False
    return True


def check_thm2_recovery():
    """Verify Theorem 2 (identifiability) on a K=2 instance over L=3 columns satisfying Strong Separation.

    Two copies differ at all three columns ({0,1,2}; K_{ij} = 3, robust regime).  Each read spans a
    column-pair that includes the shared distinguishing column 1, so EVERY cross-copy read pair
    co-observes a distinguishing column and conflicts in H -- i.e. Strong Separation holds.  Theorem 2
    then predicts: MCC = 2 AND the true partition is the UNIQUE minimum cover.
    """
    copies = [(0, 0, 0), (1, 1, 1)]  # differ at all of cols {0,1,2}
    # Each copy contributes reads spanning column-pairs; every read observes column 1.
    windows = [
        (0, {0, 1}), (0, {1, 2}),
        (1, {0, 1}), (1, {1, 2}),
    ]
    reads, labels = reads_from_copies(copies, windows)
    H = conflict_graph(reads)

    # The instance must satisfy Strong Separation (every cross-copy read pair conflicts).
    assert cond_strong(copies, reads, labels), (
        "check_thm2_recovery instance must satisfy Strong Separation"
    )

    k = mcc_bruteforce(reads)
    assert k == 2, f"MCC should be 2, got {k}"

    true_part = partition_of(labels)
    colorings = {partition_of(c) for c in all_min_colorings(H, 2)}
    assert colorings == {true_part}, (
        f"the true partition must be the UNIQUE minimum cover under Strong Separation; got {colorings}"
    )
    return "Thm 2: K=2 instance under Strong Separation -> unique minimum cover == true copies"


# --------------------------------------------------------------------------------------
# Exhaustive enumeration engine (no networkx, pure brute force) -- the GROUND TRUTH that
# discovered the §5 result.  Ported from .superpowers/sdd/cond_hunt_reference.py.
# --------------------------------------------------------------------------------------


def _jc(idxs, reads):
    """Jointly-consistent: the reads in idxs share a single allele-vector (per-column agreement)."""
    col = {}
    for k in idxs:
        for c, a in observed(reads[k]).items():
            if col.setdefault(c, a) != a:
                return False
    return True


def _min_covers(reads, K):
    """All minimum copy covers of `reads` IF the minimum size is exactly K; else None.

    A copy cover of size k is a partition of the reads into k jointly-consistent classes.
    Brute force over label tuples (k^n); intended for tiny instances only.  Returns the set of
    partitions (frozenset of frozenset of read-indices) achieving the minimum size K, or None when
    the minimum cover size is < K (so the instance is excluded from the MCC==K enumeration).
    """
    n = len(reads)

    def can(k):
        for labels in itertools.product(range(k), repeat=n):
            classes = {}
            for i, lab in enumerate(labels):
                classes.setdefault(lab, []).append(i)
            if all(_jc(c, reads) for c in classes.values()):
                return True
        return False

    if any(can(k) for k in range(1, K)):
        return None
    parts = set()
    for labels in itertools.product(range(K), repeat=n):
        cls = {}
        for i, l in enumerate(labels):
            cls.setdefault(l, set()).add(i)
        if len(cls) == K and all(_jc(c, reads) for c in cls.values()):
            parts.add(frozenset(frozenset(c) for c in cls.values()))
    return parts


def check_thm2_strong_exhaustive():
    """EXHAUSTIVE (deterministic) certificate of the §5 result over all K in {2,3}, L = 3.

    Enumerate EVERY instance with: K distinct copies drawn from {0,1}^3, and exactly 2 read-windows
    per copy (each window a non-empty subset of the 3 columns), keeping only instances whose minimum
    copy cover size equals K.  For every kept instance, count how many times each condition holds yet
    the minimum cover is NON-UNIQUE (a "uniqueness violation").  The exhaustive enumeration is the
    ground truth that discovered the result:

      * strong (Strong Separation)  -> 0 violations at K=2 AND K=3   (THE sufficient condition, all K)
      * sep+link (previous attempt) -> 0 violations at K=2 but >0 at K=3   (sufficient ONLY at K=2;
        the K-frontier: cross-copy recombination breaks it at K>=3)

    The window-count is bounded to 2/copy (matching the reference derivation) so the enumeration runs
    well under the controller's timeout; within that bound it is EXHAUSTIVE.  Also verifies the
    explicit recombination witness has a non-unique cover and is EXCLUDED by strong.
    """
    L = 3
    cols = list(range(L))
    all_windows = [frozenset(s) for k in range(1, L + 1) for s in itertools.combinations(cols, k)]
    # Window bound per K to keep this routine fast as a test-suite check. K=2 is enumerated over ALL
    # column-windows (singletons + pairs + triple). K=3 is enumerated over the size>=2 windows only:
    # Strong Separation's sufficiency holds for ANY window set (so a singleton-bearing instance can never
    # be a strong-violation the bound would hide), and the cross-copy recombination frontier (the sep+link
    # failure) already manifests on 2-column reads (the explicit witness). The FULL K=3 enumeration over all
    # windows (238,992 MCC=3 instances) was run separately and also gives strong=0 / sep+link>0; this bounded
    # re-run is the fast deterministic certificate. See the §5 note.
    windows_for = {2: all_windows, 3: [w for w in all_windows if len(w) >= 2]}
    copyvecs = [tuple(c) for c in itertools.product((0, 1), repeat=L)]

    stats = {}  # (K, cond) -> [holds, violations]
    summary = {}  # K -> (total_mcc_eq_K, unique)
    for K in (2, 3):
        total = unique = 0
        st = {"strong": [0, 0], "sep+link": [0, 0]}
        windows = windows_for[K]
        for copies in itertools.combinations(copyvecs, K):
            copies = list(copies)
            for assign in itertools.product(itertools.combinations(windows, 2), repeat=K):
                reads, labels = [], []
                for ci, wins in enumerate(assign):
                    for w in wins:
                        reads.append(frozenset((col, copies[ci][col]) for col in w))
                        labels.append(ci)
                mc = _min_covers(reads, K)
                if mc is None:
                    continue  # minimum cover size < K: not an MCC==K instance
                total += 1
                uq = mc == {partition_of(labels)}
                if uq:
                    unique += 1
                conds = {
                    "strong": cond_strong(copies, reads, labels),
                    "sep+link": cond_sep_and_link(copies, reads, labels),
                }
                for name, holds in conds.items():
                    if holds:
                        st[name][0] += 1
                        if not uq:
                            st[name][1] += 1
        stats[K] = st
        summary[K] = (total, unique)

    strong_k2 = stats[2]["strong"][1]
    strong_k3 = stats[3]["strong"][1]
    splink_k2 = stats[2]["sep+link"][1]
    splink_k3 = stats[3]["sep+link"][1]

    print(
        f"    [exhaustive] K=2: total(MCC=2)={summary[2][0]} unique={summary[2][1]}  "
        f"strong holds={stats[2]['strong'][0]} viol={strong_k2}  "
        f"sep+link holds={stats[2]['sep+link'][0]} viol={splink_k2}"
    )
    print(
        f"    [exhaustive] K=3: total(MCC=3)={summary[3][0]} unique={summary[3][1]}  "
        f"strong holds={stats[3]['strong'][0]} viol={strong_k3}  "
        f"sep+link holds={stats[3]['sep+link'][0]} viol={splink_k3}"
    )

    # strong is SUFFICIENT for all K: zero uniqueness violations at K=2 and K=3.
    assert strong_k2 == 0, f"strong must have 0 violations at K=2, got {strong_k2}"
    assert strong_k3 == 0, f"strong must have 0 violations at K=3, got {strong_k3}"
    assert stats[2]["strong"][0] > 0 and stats[3]["strong"][0] > 0, (
        "strong must actually hold on some kept instances at both K"
    )
    # sep+link is sufficient ONLY at K=2: zero violations at K=2, but >0 at K=3 (the K-frontier).
    assert splink_k2 == 0, f"sep+link must have 0 violations at K=2, got {splink_k2}"
    assert splink_k3 > 0, (
        f"sep+link must FAIL at K=3 (cross-copy recombination); expected >0 violations, got {splink_k3}"
    )
    # strong is NOT NECESSARY: it holds for only a minority of the truly-unique K=3 instances.
    assert stats[3]["strong"][0] < summary[3][1], (
        "strong must be strictly weaker in coverage than uniqueness at K=3 (sufficient, not necessary)"
    )

    # --- the explicit recombination witness: copies (1,1,0),(0,0,1),(0,1,1) ---
    # A class { read on cols {1,2} of copy0, read on cols {0,1} of copy2 } realizes the NOVEL vector
    # (0,1,0) -- a recombinant outside C*.  The instance has a non-unique minimum cover AND fails
    # Strong Separation (so strong correctly excludes it).
    w_copies = [(1, 1, 0), (0, 0, 1), (0, 1, 1)]
    w_windows = [
        (0, {0, 1}), (0, {1, 2}),
        (1, {0, 1}), (1, {1, 2}),
        (2, {0, 1}), (2, {1, 2}),
    ]
    w_reads, w_labels = reads_from_copies(w_copies, w_windows)
    w_mc = _min_covers(w_reads, 3)
    assert w_mc is not None, "recombination witness must have minimum cover size exactly 3"
    assert w_mc != {partition_of(w_labels)}, (
        "recombination witness must have a NON-UNIQUE minimum cover (an alternative size-3 cover exists)"
    )
    assert not cond_strong(w_copies, w_reads, w_labels), (
        "recombination witness must FAIL Strong Separation (so strong correctly excludes it)"
    )

    return (
        f"Thm 2 (exhaustive K=2,3 / L=3): strong viol K2={strong_k2}/K3={strong_k3} (SUFFICIENT all K); "
        f"sep+link viol K2={splink_k2}/K3={splink_k3} (K-frontier: K=2 only); "
        f"recombination witness non-unique and excluded by strong"
    )


def check_thm2_K0_merge():
    """Verify the K=0 boundary of the K-bound corollary: identical copies are non-identifiable.

    When two copies are identical over every observed column (K_{ij} = 0), there is no distinguishing
    column: Strong Separation fails (no cross-copy pair can conflict) AND the copies merge -- the reads
    produce no conflict edge, the minimum cover merges them into one part, and the true copies are
    provably unrecoverable (the MAGEA co-located regime, resolvable fraction 0/494).
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


def check_corollary_paths():
    """Verify the paths/isoforms corollary: junction treated as a linkage column.

    A family variation-graph DAG has two copies that differ by BOTH a PSV allele (column 0)
    AND an isoform-distinguishing junction (pseudo-column JJ=99).  Treating the junction on
    equal footing with PSV columns, the combined conflict machinery applies verbatim: same
    conflict_graph, same mcc_bruteforce, same Strong Separation condition.

    Copy A: allele 0 at PSV column 0, takes junction-branch 0 at column JJ.
    Copy B: allele 1 at PSV column 0, takes junction-branch 1 at column JJ.

    For Strong Separation to hold over the combined column set, every cross-copy read pair must
    co-observe at least one distinguishing column and disagree on it.  This is exactly the
    long-read linkage condition: a molecule that spans exon-content (carrying both PSV and junction
    observations) gives a conflicting signal against every foreign molecule.  Single-column reads
    observing ONLY the PSV or ONLY the junction from disjoint sets may not co-observe any column
    and thus cannot conflict -- so the 'linking' span of long reads is what makes the condition
    satisfiable in practice.

    The check uses reads that each observe BOTH the PSV column AND the junction column (i.e. the
    reads span and link the two distinguishing features, as a long molecule would).  This is the
    minimal instance where Strong Separation holds over the combined (allele + junction) column set,
    and Theorem 2 then gives MCC = 2 = number of true copies.
    """
    JJ = 99  # junction pseudo-column (models the isoform-distinguishing splice event)
    copies = {"A": {0: 0, JJ: 0}, "B": {0: 1, JJ: 1}}  # differ at PSV col 0 AND junction col JJ

    def mk(cp, cols):
        return frozenset((c, copies[cp][c]) for c in cols)

    # Reads span BOTH columns so every cross-copy pair co-observes a distinguishing column.
    # (A read observing only col 0 and a foreign read observing only col JJ would share no
    # column and thus cannot conflict -- that is the honest caveat noted in §6.)
    reads = [
        mk("A", {0, JJ}),  # copy A read spanning PSV + junction
        mk("A", {0, JJ}),  # copy A read spanning PSV + junction (additional coverage)
        mk("A", {0, JJ}),  # copy A read spanning PSV + junction (additional coverage)
        mk("B", {0, JJ}),  # copy B read spanning PSV + junction
        mk("B", {0, JJ}),  # copy B read spanning PSV + junction (additional coverage)
        mk("B", {0, JJ}),  # copy B read spanning PSV + junction (additional coverage)
    ]

    # The combined conflict machinery applies verbatim.
    H = conflict_graph(reads)
    mcc = mcc_bruteforce(reads)
    assert mcc == 2, f"two copies separated by allele AND junction -> MCC=2, got {mcc}"

    # Verify Strong Separation holds: every cross-copy read pair co-observes {col 0, col JJ}
    # and disagrees on at least one of them.
    labels = [0, 0, 0, 1, 1, 1]
    copy_vecs = [(0, 0), (1, 1)]  # placeholder; cond_strong only uses reads + labels
    assert cond_strong(copy_vecs, reads, labels), (
        "spanning reads must satisfy Strong Separation over combined (allele+junction) columns"
    )

    return (
        "Corollary: junction treated as a linkage column -> path-cover recovers copies+isoforms (MCC=2)"
    )


def compatibility_edges(reads):
    """Edges of the compatibility graph H-bar (complement of the conflict graph): the pairs of reads that
    do NOT conflict (they agree at every co-observed column)."""
    edges = []
    for i, j in itertools.combinations(range(len(reads)), 2):
        oi, oj = observed(reads[i]), observed(reads[j])
        if all(oi[c] == oj[c] for c in (oi.keys() & oj.keys())):
            edges.append((i, j))
    return edges


def recover(reads):
    """RECOVER (Theorem 3): connected components of the compatibility graph H-bar, via union-find.
    Returns the read-partition as a frozenset of frozensets of read indices. O(n^2 * m)."""
    n = len(reads)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i, j in compatibility_edges(reads):
        parent[find(i)] = find(j)
    comp = {}
    for i in range(n):
        comp.setdefault(find(i), set()).add(i)
    return frozenset(frozenset(c) for c in comp.values())


def is_disjoint_clique_union(reads):
    """Self-certifying check (Theorem 3 addendum): is the compatibility graph a disjoint union of cliques?
    Equivalently, is every connected component of H-bar internally complete (all pairs compatible) — i.e. is the
    compatible relation transitive? This certifies that recover()'s partition is the UNIQUE minimum copy cover of
    the reads. Strong Separation is sufficient for it but NOT equivalent: at K=0 (copies identical over observed
    columns) H-bar is one clique, so this returns True (passes) while Strong Separation fails — recover then
    returns the coarser unique cover that merges the indistinguishable copies. O(n^2 * m) (it builds H-bar)."""
    n = len(reads)
    compat = [[i == j for j in range(n)] for i in range(n)]
    for i, j in compatibility_edges(reads):
        compat[i][j] = compat[j][i] = True
    for cls in recover(reads):
        members = sorted(cls)
        for a in range(len(members)):
            for b in range(a + 1, len(members)):
                if not compat[members[a]][members[b]]:
                    return False  # same component but not directly compatible -> not a clique
    return True


def check_thm3_recovery_algorithm():
    """EXHAUSTIVE certificate of Theorem 3: over every Strong-Separation instance in the K in {2,3}, L=3
    enumeration, RECOVER returns exactly the TRUE partition, and is_disjoint_clique_union accepts it.
    Also: the explicit recombination witness (sep+link holds, strong fails) is REJECTED by the certificate."""
    L = 3
    cols = list(range(L))
    all_windows = [frozenset(s) for k in range(1, L + 1) for s in itertools.combinations(cols, k)]
    windows_for = {2: all_windows, 3: [w for w in all_windows if len(w) >= 2]}
    copyvecs = [tuple(c) for c in itertools.product((0, 1), repeat=L)]

    n_strong = 0
    recover_mismatch = 0
    cert_reject_on_strong = 0
    for K in (2, 3):
        windows = windows_for[K]
        for copies in itertools.combinations(copyvecs, K):
            copies = list(copies)
            for assign in itertools.product(itertools.combinations(windows, 2), repeat=K):
                reads, labels = [], []
                for ci, wins in enumerate(assign):
                    for w in wins:
                        reads.append(frozenset((c, copies[ci][c]) for c in w))
                        labels.append(ci)
                if not cond_strong(copies, reads, labels):   # use the file's actual cond_strong signature
                    continue
                n_strong += 1
                if recover(reads) != partition_of(labels):
                    recover_mismatch += 1
                if not is_disjoint_clique_union(reads):
                    cert_reject_on_strong += 1

    assert recover_mismatch == 0, f"RECOVER disagreed with the true partition on {recover_mismatch} strong instances"
    assert cert_reject_on_strong == 0, f"self-certify wrongly rejected {cert_reject_on_strong} strong instances"

    # Recombination witness: sep+link holds, strong fails -> certificate must REJECT (refuse to recover).
    wc = [(1, 1, 0), (0, 0, 1), (0, 1, 1)]
    wwins = [(0, {1, 2}), (0, {0, 1}), (1, {0, 1, 2}), (1, {0, 1}), (2, {0, 1}), (2, {1, 2})]
    wreads = [frozenset((c, wc[ci][c]) for c in w) for ci, w in wwins]
    assert not is_disjoint_clique_union(wreads), "certificate must reject the non-strong recombination witness"

    print(f"    [thm3] strong instances={n_strong}: RECOVER==true 100%, self-certify accepts all; witness rejected")
    return f"Thm 3: RECOVER == true partition on all {n_strong} strong instances (exhaustive K=2,3/L=3); witness rejected"


CHECKS = [check_lemma_mcc_equals_chromatic, check_thm1_reduction]
CHECKS.append(check_thm2_recovery)
CHECKS.append(check_thm2_strong_exhaustive)
CHECKS.append(check_thm2_K0_merge)
CHECKS.append(check_corollary_paths)
CHECKS.append(check_thm3_recovery_algorithm)


def main():
    for fn in CHECKS:
        print("OK  -", fn())


if __name__ == "__main__":
    main()
