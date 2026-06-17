#!/usr/bin/env python3
"""Upgrade #1 (OPTION, not a replacement): a per-family ALL-COPIES POA variation graph yields a clean,
formally-definable multi-copy-family criterion.

FORMAL DEFINITION (graph-based, level (theta, T)):
  Given member sequences S={s_1..s_N}, let G be their partial-order (POA) alignment graph with columns
  c_1..c_M in topological order. For column c let supp(c)=#{i : s_i is non-gap at c}. The CONSERVED CORE
  is the longest contiguous run R of columns with supp(c) >= theta*N for all c in R. S is a multi-copy
  family at level (theta,T) iff  |R| >= T * median_i |s_i|.
  -> family_score(S) = |R| / median|s_i|  in [0,1].

This is the N-way generalization of the validated pairwise contiguous-core: at N=2, theta*N=2 forces both
copies non-gap, so R is exactly the longest co-aligned block = the pairwise core (which we validated).
It is graph-native (one statistic from one graph, no O(N^2) pairwise / transitive closure), and it
exposes over-merge: a domain-hub chain has NO column shared by ~all members -> tiny core -> rejected.

Run with /home/juanfra/miniforge3/bin/python (pyabpoa). Deterministic.
"""
import csv
import math
import os
import statistics

import pyabpoa

BASE = os.path.dirname(__file__)
FA = "/tmp/gene_reps_gw.fa"
FAMS = os.path.join(BASE, "families.tsv")          # pairwise components (candidates)
OUT_TSV = os.path.join(BASE, "family_graph_scores.tsv")
THETA = 0.5       # MAJORITY core: a core column is shared by >= ceil(theta*N) copies (floor 2).
T = 0.13          # validated pairwise-core threshold (the graph def reduces to pairwise at N=2)
CAP = 30          # cap members per graph (longest, to bound POA memory)
LEN_CAP = 15000   # skip members longer than this (titin-class; POA matrix blows up)

# labeled validation sets
CURATED = {
    "RABL2": ["RABL2A", "RABL2B"], "DAZ": ["DAZ1", "DAZL"],
    "APOBEC3": ["APOBEC3B", "APOBEC3C", "APOBEC3D", "APOBEC3F", "APOBEC3G", "APOBEC3H"],
    "RFPL": ["RFPL1", "RFPL2", "RFPL3", "RFPL4A"],
}
DOMAIN_SHARERS = {
    "CDPF1/PPARA": ["CDPF1", "PPARA"], "CREB1/METTL21A": ["CREB1", "METTL21A"],
    "GCA/KCNH7": ["GCA", "KCNH7"], "CASP8/FLACC1": ["CASP8", "FLACC1"],
    "GPR39/LYPD1": ["GPR39", "LYPD1"], "ASDURF/ASNSD1": ["ASDURF", "ASNSD1"],
}


def load_fa(path):
    seqs = {}; name = None; buf = []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name = line[1:].split()[0]; buf = []
        else:
            buf.append(line)
    if name:
        seqs[name] = "".join(buf)
    return seqs


def graph_score(seqs, theta=THETA):
    """POA MSA -> column support profile -> longest contiguous core run / median length.
    Returns (None, 0) if the POA fails (e.g. memory on pathological inputs)."""
    if len(seqs) < 2:
        return 0.0, 0
    try:
        al = pyabpoa.msa_aligner()
        res = al.msa([s for s in seqs], out_cons=False, out_msa=True)
        rows = [r.decode() if isinstance(r, bytes) else r for r in res.msa_seq]
    except Exception:
        return None, 0
    if not rows:
        return 0.0, 0
    n = len(rows); M = len(rows[0])
    need = max(2, math.ceil(theta * n))   # majority core, >=2 copies (=> reduces to pairwise at N=2)
    best = run = 0
    for c in range(M):
        supp = sum(1 for r in rows if r[c] != "-")
        if supp >= need:
            run += 1; best = max(best, run)
        else:
            run = 0
    med = statistics.median(len(s) for s in seqs)
    return (best / med if med else 0.0), best


def main():
    fa = load_fa(FA)

    def score_set(members):
        seqs = [fa[g] for g in members if g in fa and len(fa[g]) <= LEN_CAP]
        if len(seqs) < 2:
            return None
        seqs = sorted(seqs, key=len)[:CAP]   # shortest first (bounds POA memory)
        sc, core = graph_score(seqs)
        if sc is None:
            return None
        return (round(sc, 3), core, len(seqs))

    # labeled validation
    print("=== labeled validation (graph core_score) ===")
    cur = {}
    for f, ms in CURATED.items():
        r = score_set(ms)
        if r:
            cur[f] = r; print(f"  CURATED {f}: score={r[0]} core={r[1]} N={r[2]}")
    dom = {}
    for f, ms in DOMAIN_SHARERS.items():
        r = score_set(ms)
        if r:
            dom[f] = r; print(f"  DOMAIN-SHARER {f}: score={r[0]} core={r[1]} N={r[2]}")

    # all pairwise-candidate families -> graph score (expose over-merge on big ones)
    rows = []
    with open(FAMS) as fh:
        fh.readline()
        for line in fh:
            fid, n, genes = line.rstrip("\n").split("\t")
            members = genes.split(",")
            r = score_set(members)
            if r:
                rows.append((fid, int(n), r[0], r[1], r[2]))
    rows.sort(key=lambda x: -x[1])   # by size

    import numpy as np
    sc_all = np.array([r[2] for r in rows])
    big = [r for r in rows if r[1] >= 25]
    small = [r for r in rows if r[1] <= 5]
    print(f"\n=== {len(rows)} pairwise-candidate families graph-scored ===")
    print(f"  graph_score: median={np.median(sc_all):.2f} ; >=T({T}): {int((sc_all>=T).sum())}")
    print(f"  big components (N>=25, over-merge suspects): {len(big)}, "
          f"median graph_score={np.median([r[2] for r in big]):.2f} (LOW => not one family)")
    print(f"  small components (N<=5): {len(small)}, "
          f"median graph_score={np.median([r[2] for r in small]):.2f}")

    with open(OUT_TSV, "w") as fh:
        fh.write("family_id\tn_members\tgraph_score\tcore_cols\tn_scored\n")
        for fid, n, sc, core, ns in rows:
            fh.write(f"{fid}\t{n}\t{sc}\t{core}\t{ns}\n")
    # stash labeled for the figure
    import json
    json.dump({"curated": cur, "domain_sharers": dom,
               "big_med": float(np.median([r[2] for r in big])) if big else 0,
               "small_med": float(np.median([r[2] for r in small])) if small else 0,
               "theta": THETA, "T": T},
              open(os.path.join(BASE, "family_graph_labeled.json"), "w"), indent=2)
    print(f"\n[wrote {OUT_TSV} + family_graph_labeled.json]")


if __name__ == "__main__":
    main()
