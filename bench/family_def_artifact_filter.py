#!/usr/bin/env python3
"""family_def_artifact_filter.py — VG-TOPOLOGICAL artifact filter for the multi-copy-family definition.

Drops ARTIFACT families by GRAPH STRUCTURE (no ORF / coding requirement — a real lncRNA family is KEPT):
  * repeat_cyclic   — members are internal repeats -> a CYCLE in the variation graph
                      (intrinsic signal: long-k-mer SELF-RECURRENCE; 4-mer complexity saturates and
                      misses these). Real family = acyclic bubble-chain.
  * overmerge_sparse— the component is NOT a clique (sparse bridge through an articulation point,
                      e.g. the comp-238 lncRNA-repeat blob). Real family = dense reciprocal clique.
  * te_highcopy     — dispersed TE/satellite array (too many copies).

DEFAULT-ON. Opt out with FAMILY_DEF_NO_ARTIFACT_FILTER=1 -> nothing is dropped, but every family is
still TAGGED with its class, so you can keep the repeat / over-merge / lncRNA families and VISUALISE
their graphs (e.g. to show an advisor what a repeat/lncRNA family graph looks like vs a clean family).

Thresholds are env-tunable: FAMILY_DEF_MAX_COPIES (80), FAMILY_DEF_MIN_DENSITY (0.30),
FAMILY_DEF_REPEAT_RECUR (0.15). No coding model, no arbitrary SEQUENCE threshold — pure topology.
"""
import os
from collections import Counter

DISABLE = os.environ.get("FAMILY_DEF_NO_ARTIFACT_FILTER") == "1"
MAX_COPIES = int(os.environ.get("FAMILY_DEF_MAX_COPIES", "80"))   # high backstop only (real cliques reach 40+)
MIN_DENSITY = float(os.environ.get("FAMILY_DEF_MIN_DENSITY", "0.30"))
REPEAT_RECUR = float(os.environ.get("FAMILY_DEF_REPEAT_RECUR", "0.15"))
KMER = 15


def kmer_self_recurrence(seq, k=KMER):
    """fraction of k-mers that RECUR within the sequence = internal-repeat = a cycle in the VG.
    median ~0 for real transcripts (unique k-mers, acyclic); high for tandem/internal repeats."""
    seq = seq.upper()
    if len(seq) < 2 * k:
        return 0.0
    c = Counter(seq[i:i + k] for i in range(len(seq) - k + 1))
    tot = sum(c.values())
    return sum(v for v in c.values() if v > 1) / tot if tot else 0.0


def classify_family(members, subgraph, seq_provider=None):
    """Return (class, metrics). `subgraph` = the family's induced graph (for cliqueness).
    `seq_provider(member)->seq` enables the repeat/cycle check (optional)."""
    n = len(members)
    e = subgraph.number_of_edges()
    density = (2 * e) / (n * (n - 1)) if n > 1 else 1.0   # clique = 1.0; over-merge bridge = low
    recur = None
    if seq_provider is not None:
        rs = [kmer_self_recurrence(s) for s in (seq_provider(m) for m in members) if s]
        if rs:
            rs.sort()
            recur = rs[len(rs) // 2]   # median over members
    metrics = dict(n=n, density=round(density, 3), recur=(round(recur, 3) if recur is not None else None))
    # DENSITY (cliqueness), not copy-count, is the artifact discriminator: a real family is a dense
    # reciprocal CLIQUE at any size (a 41-copy clique = a real expansion); an over-merge / dispersed
    # blob is SPARSE. Copy-count is only a high backstop for a dense-but-implausibly-huge component.
    if recur is not None and recur >= REPEAT_RECUR:
        klass = "repeat_cyclic"
    elif n >= 4 and density < MIN_DENSITY:
        klass = "overmerge_sparse"
    elif n > MAX_COPIES:
        klass = "te_highcopy"        # dense, non-repeat, yet implausibly large -> flag for review
    else:
        klass = "real"
    return klass, metrics


def partition_families(DG, seq_provider=None, min_size=2):
    """Classify every connected component (>= min_size). Returns a list of
    (members, class, metrics, kept) sorted by size desc. `kept` = (DISABLE or class == 'real')."""
    import networkx as nx
    out = []
    for comp in nx.connected_components(DG):
        if len(comp) < min_size:
            continue
        members = sorted(comp)
        klass, metrics = classify_family(members, DG.subgraph(comp), seq_provider)
        out.append((members, klass, metrics, DISABLE or klass == "real"))
    out.sort(key=lambda x: -len(x[0]))
    return out


if __name__ == "__main__":
    import networkx as nx, random
    rng = random.Random(0)
    def uniq():  # genuinely non-repetitive sequence (near-unique k-mers)
        return "".join(rng.choice("ACGT") for _ in range(2000))
    # self-test: a clique (real), a sparse bridge (over-merge), a high-copy (TE), a repeat-seq (cyclic)
    G = nx.Graph()
    G.add_edges_from([("a1", "a2"), ("a1", "a3"), ("a2", "a3")])                 # clique of 3 = real
    for i in range(1, 7):
        G.add_edge("bx", f"b{i}")                                                # 7-node star (density 0.29) = over-merge
    for i in range(20):
        G.add_edge("t0", f"t{i}")                                                # 21-node star = TE/high-copy
    G.add_edges_from([("r1", "r2"), ("r1", "r3"), ("r2", "r3")])                 # clique but repeat seqs
    seqs = {f"r{i}": "ACGTACGT" * 200 for i in (1, 2, 3)}                         # tandem repeat -> high recurrence
    for m in list(G.nodes()):
        seqs.setdefault(m, uniq())
    res = partition_families(G, seq_provider=lambda m: seqs.get(m))
    print(f"DISABLE={DISABLE}  (density-primary: sparse=over-merge, dense clique=real at any size)")
    for members, klass, metrics, kept in res:
        print(f"  {klass:16s} kept={kept} n={metrics['n']} density={metrics['density']} recur={metrics['recur']}  {members[:4]}")
    tally = Counter(k for _, k, _, _ in res)
    # the 21-node STAR is sparse -> over-merge (NOT te_highcopy); the repeat-clique -> repeat_cyclic; the
    # unique-seq clique -> real. (te_highcopy is now only a high-copy DENSE backstop, not exercised here.)
    assert tally["overmerge_sparse"] >= 1 and tally["repeat_cyclic"] >= 1 and tally["real"] >= 1, tally
    print("self-test OK:", dict(tally))
