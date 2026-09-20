#!/usr/bin/env python3
"""NPIP/TBC1D3 layer order — P-layer stability and variants (audit 2026-09-16). Read-only on the saved searches.

Inputs: light/work/P/{proteins.index.tsv, blastp.tsv, searched.txt, S.round0.txt} through
light/scripts/layer_protein_bounded.py load()/hops() (edges = bench/protein_families.pair_hsps + edges_from, the shipped
§6ko greedy-by-bitscore HSP order; weight identity x coverage; MCL = bench/mcl_port.mcl, inflation 2.8).

Outputs (integrate_slim/):
  P_stability_plain.tsv   plain §6ko rule, member-containing MCL clusters on the edges among N_k, k = 1..8
                          (k = 8 has 2 unsearched genes; they are dropped before MCL and reported)
  P_aa50.tsv              §6ko aa-identity >= 0.50 stratum (ledger §6ko qualification): N_k grown over aa>=0.50 edges,
                          member-containing MCL clusters, k = 1..5; and the exact connected components holding members
  P_N2_clusters.tsv       every gene of N_2 (plain rule) with its k = 2 MCL cluster (used by the whole-group JOIN)
  P_variant_labels.tsv    gene_id, name, pid, P_aa50_mcl (k=2 member cluster or ''), P_aa50_comp (member component or '')
"""
import collections
import csv
import sys
import os
# §6r9: repo root from THIS file, so the tool runs from any clone (it used to hardcode
# /mnt/c/Users/jfris/Desktop/Rustle, which only ever worked on one machine).
_RUSTLE_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

import time

sys.path.insert(0, "/mnt/linuxdisk/home/juanfraitu/layer_order/npip_tbc1d3/light/scripts")
sys.path.insert(0, os.path.join(_RUSTLE_REPO, 'bench'))
import layer_protein_bounded as lpb  # noqa: E402
import mcl_port  # noqa: E402

LIGHT = "/mnt/linuxdisk/home/juanfraitu/layer_order/npip_tbc1d3/light"
INT = "/mnt/linuxdisk/home/juanfraitu/layer_order/npip_tbc1d3/integrate_slim"


def tsv(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


def mcl_clusters(nodes, edges):
    Ek = {e: v[0] for e, v in edges.items() if e[0] in nodes and e[1] in nodes}
    cl = [list(c) for c in mcl_port.mcl(Ek)] if Ek else []
    seen = {x for c in cl for x in c}
    cl += [[x] for x in sorted(nodes - seen, key=int)]
    return cl, len(Ek)


def main():
    t0 = time.time()
    idx, E, raw, searched, adj = lpb.load()
    name = {p: idx[p]["name"] for p in idx}
    N0 = set(lpb.read_list(f"{lpb.W}/S.round0.txt"))
    memb = {r["gene_id"] for r in tsv(f"{LIGHT}/members.corrected.tsv")}
    g2p = {v["gene_id"]: p for p, v in idx.items()}
    assert N0 == {g2p[g] for g in memb if g in g2p}, "S.round0 != coding corrected members"
    print(f"[load] proteins {len(idx)}; searched {len(searched)}; §6ko edges from all searches {len(E)}; member proteins "
          f"{len(N0)}; {time.time() - t0:.1f}s", flush=True)

    # ---------------------------------------------------------------- plain rule, k = 1..8
    N, rep = lpb.hops(8, idx, E, raw, searched, adj)
    rows, ref = [], None
    for k in range(1, 9):
        Nk = N[k]
        uns = Nk - searched
        cl, ne = mcl_clusters(Nk & searched, E)
        mc = sorted((sorted(name[p] for p in c) for c in cl if set(c) & N0), key=lambda c: (-len(c), c))
        if k == 2:
            ref = mc
            with open(f"{INT}/P_N2_clusters.tsv", "w") as fh:
                fh.write("gene_id\tname\tpid\tk2_cluster\tcluster_size\tholds_member\n")
                for i, c in enumerate(sorted(cl, key=lambda c: (-len(c), sorted(name[p] for p in c)))):
                    for p in sorted(c, key=lambda p: name[p]):
                        fh.write(f"{idx[p]['gene_id']}\t{name[p]}\t{p}\tN2c{i}\t{len(c)}\t"
                                 f"{'yes' if set(c) & N0 else 'no'}\n")
        rows.append({"k": k, "N_k": len(Nk), "unsearched_dropped": len(uns), "edges": ne,
                     "member_clusters": " | ".join(f"{len(c)}:{','.join(c)}" for c in mc),
                     "identical_to_k2": "" if ref is None else ("yes" if mc == ref else "no")})
        print(f"[plain] k={k} |N_k| {len(Nk)} unsearched {len(uns)} edges {ne} member clusters "
              f"{[len(c) for c in mc]} identical_to_k2 {rows[-1]['identical_to_k2']} {time.time() - t0:.0f}s", flush=True)
    with open(f"{INT}/P_stability_plain.tsv", "w") as fh:
        cols = list(rows[0].keys())
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    # ---------------------------------------------------------------- aa identity >= 0.50
    E50 = {e: v for e, v in E.items() if v[1] >= 0.50}
    adj50 = collections.defaultdict(set)
    for u, v in E50:
        adj50[u].add(v)
        adj50[v].add(u)
    rows50 = []
    Nk = set(N0)
    lab_mcl = {}
    for k in range(1, 6):
        Nk = Nk | {v for u in Nk for v in adj50[u]}
        cl, ne = mcl_clusters(Nk, E50)
        mc = sorted((c for c in cl if set(c) & N0), key=lambda c: (-len(c), sorted(name[p] for p in c)))
        if k == 2:
            for i, c in enumerate(mc):
                for p in c:
                    lab_mcl[p] = f"P50mcl|{i}" if len(c) >= 2 else ""
        rows50.append({"kind": "mcl", "k": k, "N_k": len(Nk), "all_searched": "yes" if Nk <= searched else "no",
                       "edges": ne, "member_clusters": " | ".join(f"{len(c)}:{','.join(sorted(name[p] for p in c))}"
                                                                  for c in mc)})
        print(f"[aa50] k={k} |N_k| {len(Nk)} all searched {Nk <= searched} edges {ne} member clusters "
              f"{[len(c) for c in mc]}", flush=True)
    seen, comps = set(), []
    for s in sorted(N0, key=int):
        if s in seen:
            continue
        st, comp = [s], set()
        while st:
            x = st.pop()
            if x in comp:
                continue
            comp.add(x)
            st += list(adj50[x] - comp)
        seen |= comp
        comps.append(comp)
    lab_comp = {}
    for i, c in enumerate(sorted(comps, key=lambda c: (-len(c), sorted(name[p] for p in c)))):
        for p in c:
            lab_comp[p] = f"P50comp|{i}" if len(c) >= 2 else ""
        nb = {p: sorted((name[q], round(E50[(min(p, q, key=int), max(p, q, key=int))][1], 3)) for q in adj50[p])
              for p in c if p not in N0}
        rows50.append({"kind": "component", "k": "", "N_k": len(c), "all_searched": "yes" if c <= searched else "no",
                       "edges": sum(1 for e in E50 if e[0] in c and e[1] in c),
                       "member_clusters": f"{len(c)}:{','.join(sorted(name[p] for p in c))}; non-member aa>=0.50 "
                                          f"neighbours: {nb}"})
        print(f"[aa50] component {len(c)} all searched {c <= searched}: {sorted(name[p] for p in c)}; non-members {nb}")
    with open(f"{INT}/P_aa50.tsv", "w") as fh:
        cols = list(rows50[0].keys())
        fh.write("\t".join(cols) + "\n")
        for r in rows50:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    with open(f"{INT}/P_variant_labels.tsv", "w") as fh:
        fh.write("gene_id\tname\tpid\tP_aa50_mcl\tP_aa50_comp\n")
        for p in sorted(set(lab_mcl) | set(lab_comp), key=lambda p: name[p]):
            fh.write(f"{idx[p]['gene_id']}\t{name[p]}\t{p}\t{lab_mcl.get(p, '')}\t{lab_comp.get(p, '')}\n")
    # the aa of the plain-rule edges that aa>=0.50 removes from PC2
    pc2 = [r for r in tsv(f"{LIGHT}/P.edges.tsv") if "TBC1D26" in (r["u_name"], r["v_name"]) or
           "USP6NL" in (r["u_name"], r["v_name"])]
    for gname in ("TBC1D26", "USP6NL"):
        ids = [float(r["identity"]) for r in pc2 if gname in (r["u_name"], r["v_name"])]
        print(f"[aa50] plain-rule P.edges of {gname}: {len(ids)} edges, aa identity {min(ids):.3f}-{max(ids):.3f}")
    print(f"[done] {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
