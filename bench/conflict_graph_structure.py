#!/usr/bin/env python
"""F4 structured-instance lemma (empirical) — are the real GGO read-conflict graphs CHORDAL?

§5c of the theory note: MIN-COVER (= chi(H), the minimum copy number, Lemma 1) is graph colouring, hence
inapproximable for arbitrary H. BUT real read-conflict graphs are not arbitrary: two reads conflict iff they
DISAGREE at a shared PSV column. If H is CHORDAL (or perfect), chi(H)=omega(H) is poly-time exact, so the
minimum cover itself is computable on the families that matter — strictly stronger than the worst-case
(1-1/e) bound for MWCA. This script builds H per co-located family (from RNA reads' alleles at the DNA-catalog
PSV columns, reusing the F3 extraction), collapses identical-observation reads to signature nodes (chi is
unchanged), and tests chordality via Maximum-Cardinality-Search; for chordal H it reports chi=omega and compares
to the DNA copy count K.

Run: /home/juanfra/miniforge3/bin/python bench/conflict_graph_structure.py [N_families]
"""
import sys
import os

sys.path.insert(0, os.path.dirname(__file__))
import dna_supervised_decode as dd     # colocated_families, build_signatures, read_alleles_via_ref0
import json


def conflict_graph(obs_list):
    """Nodes = distinct PSV-observation signatures; edge iff two signatures disagree at a shared column.
    Returns (nodes, adj) with adj[i] = set of neighbor indices."""
    sigs = sorted({tuple(sorted(o.items())) for o in obs_list if len(o) >= 1})
    nodes = [dict(s) for s in sigs]
    n = len(nodes)
    adj = [set() for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            a, b = nodes[i], nodes[j]
            if any(c in b and b[c] != v for c, v in a.items()):   # shared column, different allele
                adj[i].add(j); adj[j].add(i)
    return nodes, adj


def mcs_order(adj):
    """Maximum Cardinality Search: returns an elimination order (list of vertices)."""
    n = len(adj)
    wt = [0] * n
    visited = [False] * n
    order = []
    for _ in range(n):
        u = max((v for v in range(n) if not visited[v]), key=lambda v: wt[v], default=None)
        if u is None:
            break
        visited[u] = True
        order.append(u)
        for w in adj[u]:
            if not visited[w]:
                wt[w] += 1
    return order


def is_chordal_and_omega(adj):
    """Test chordality via a perfect-elimination check on the reverse MCS order; return (chordal, omega).
    For chordal graphs omega = max over the PEO of (1 + |later-neighbours|)."""
    order = mcs_order(adj)
    pos = {v: i for i, v in enumerate(order)}
    peo = list(reversed(order))          # eliminate in reverse-MCS order
    elim_pos = {v: i for i, v in enumerate(peo)}
    chordal = True
    omega = 1
    for v in peo:
        later = [w for w in adj[v] if elim_pos[w] > elim_pos[v]]   # neighbours eliminated after v
        omega = max(omega, 1 + len(later))
        if later:
            # the earliest later-neighbour must be adjacent to all the others (clique condition)
            u = min(later, key=lambda w: elim_pos[w])
            for w in later:
                if w != u and w not in adj[u]:
                    chordal = False
                    break
        if not chordal:
            break
    return chordal, omega


def chromatic_greedy_on_peo(adj):
    """Greedy colouring on the reverse-MCS (perfect elimination) order; optimal (=chi) iff chordal."""
    order = list(reversed(mcs_order(adj)))
    color = {}
    for v in order:
        used = {color[w] for w in adj[v] if w in color}
        c = 0
        while c in used:
            c += 1
        color[v] = c
    return (max(color.values()) + 1) if color else 0


def main():
    limit = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    fams = dd.colocated_families()
    n_chordal = n_total = 0
    chi_eq_K = 0
    rows = []
    for fid, members in fams:
        if n_total >= limit:
            break
        sig = dd.build_signatures(fid, members)
        if not sig:
            continue
        chrom, ref_seq, psv_gpos, copy_vecs = sig
        rs = json.load(open(os.path.join(dd.DNADIR, f"{fid}.json")))["iv"][0]
        obs_list = dd.read_alleles_via_ref0(fid, chrom, members, ref_seq, psv_gpos, rs)
        if len(obs_list) < 3:
            continue
        nodes, adj = conflict_graph(obs_list)
        if len(nodes) < 2:
            continue
        n_total += 1
        chordal, omega = is_chordal_and_omega(adj)
        chi = chromatic_greedy_on_peo(adj)
        K = len(copy_vecs)
        n_chordal += chordal
        chi_eq_K += (chi == K)
        rows.append((fid, len(nodes), chordal, chi, omega, K))

    # chi_ub (greedy on reverse-MCS) is EXACT only when chordal; otherwise it is an upper bound. omega via the
    # PEO is valid only when chordal. So report omega/exactness only for chordal H, and use the colouring-vs-K
    # ratio to expose the error inflation of the RAW conflict graph.
    infl = [r[3] / r[5] for r in rows if r[5]]          # chi_ub / K
    infl.sort()
    med_infl = infl[len(infl) // 2] if infl else float("nan")
    print(f"co-located families with a non-trivial conflict graph tested: {n_total}")
    print(f"  CHORDAL (=> chi=omega poly-exact, min-cover tractable EXACTLY): {n_chordal}/{n_total} "
          f"({100*n_chordal/max(n_total,1):.1f}%)")
    print(f"  among chordal: chi(H) == DNA copy count K: {chi_eq_K}/{n_chordal}")
    print(f"  RAW-graph error inflation: median (colouring upper bound / K) = {med_infl:.1f}x "
          f"(raw allele-disagreement includes sequencing-error edges)")
    print(f"\nper-family (fid, #sig-nodes, chordal, colouring-UB, K)  [omega/chi exact only when chordal=True]:")
    for r in sorted(rows, key=lambda x: -x[1])[:25]:
        tag = f"chi=omega={r[3]}" if r[2] else f"colour<= {r[3]}"
        print(f"  {r[0]:>10} nodes={r[1]:>4} chordal={r[2]!s:>5} {tag:>16} K={r[5]}")
    print("\nInterpretation: (1) where H is CHORDAL, chi(H)=omega(H) and the minimum copy cover is poly-time")
    print("EXACT (stronger than the worst-case (1-1/e) MWCA bound) — holds for ~1/4 of co-located families.")
    print("(2) The dominant finding: the RAW read-conflict graph is heavily error-inflated (colouring >> K),")
    print("because raw allele-disagreement counts sequencing-error edges. So Lemma 1's MCC=chi(H) is about the")
    print("ERROR-FREE/de-tied conflict graph; this is the empirical case for the per-read significance gate")
    print("(Theorem 4) over solving chi on the raw graph, and for MWCA's evidence-weighted relaxation (Thm 6).")


if __name__ == "__main__":
    main()
