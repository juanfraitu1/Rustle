#!/usr/bin/env python3
"""
family_def_junction_genomewide.py — angle A2, GENOME SCALE.

Question: does INTRON-ARCHITECTURE concordance reduce the cDNA-homology OVER-MERGE, and
specifically cut the DOMAIN-BRIDGE edges of the 389-member mega-family, while keeping the
airtight real families (RABL2 / APOBEC3 / RFPL) connected?

Pipeline:
  1. Build the cDNA homology 'truth' graph from rep_ava.tsv: an edge (gA,gB) iff
     id >= 0.90 AND max(cov_a,cov_b) >= 0.30  (reproduces the 1460-family / 389-mega-component
     baseline). Vertices = gene Names.
  2. Load the per-gene intron index (family_def_intron_index.py output).
  3. Junction-concordance edge filter — keep a homology edge only if the two genes share
     intron architecture. We test a LADDER of principled criteria (no magic numbers where
     avoidable):
        C0 baseline                      : keep all homology edges
        Cspliced : BOTH genes have >=1 intron   (kills retrocopy/intronless edges)
        Ccount.J : both spliced AND Jaccard of intron-count compatibility ...
        Carch    : both spliced AND ordered intron-length vectors concordant over the
                   homology-aligned span (DTW-free: compare the n_intron and the matched
                   intron-length profile with relative tolerance), principled tolerance =
                   15% rel / 20bp abs (same tol as the read-level concordance script,
                   reused for consistency — not a new tuned constant).
  4. MEASURE per criterion: edges kept/dropped, #connected components (families), MAX
     component size (does 389 break up?), edges removed INSIDE the mega-family vs OUTSIDE,
     and whether the airtight panel families stay connected.
  5. Compare to the disjoint-loci baseline (307 co-location edges removed, 389 -> 382).

Run: /home/juanfra/miniforge3/bin/python bench/family_def_junction_genomewide.py
"""
import collections
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from family_def_read_filters import dna_homology  # rep_ava.tsv based

INDEX = "/home/juanfra/winloci_scratch/gene_intron_index.json"
GENES_BED = "/home/juanfra/winloci_scratch/unmapped_poc/genes.bed"

ID_MIN = 0.90
COV_MIN = 0.30
ABS_TOL = 20
REL_TOL = 0.15

# airtight panel (from bench/family_def_airtight_panel.tsv)
PANEL_REAL = {
    "RABL2": ["RABL2A", "RABL2B"],
    "APOBEC3": ["APOBEC3D", "APOBEC3F"],
    "RFPL": ["RFPL1", "RFPL2", "RFPL3"],
}
# known retrocopy / processed-pseudogene counterexamples (cDNA id ~1.0, intronless copy)
RETRO_PROBE = [("EEF1A1", None), ("CNN2", None)]


# ---------- union-find over an edge set ----------
def components(edges):
    parent = {}
    def find(x):
        parent.setdefault(x, x)
        r = x
        while parent[r] != r:
            r = parent[r]
        while parent[x] != r:
            parent[x], x = r, parent[x]
        return r
    for a, b in edges:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    comp = collections.defaultdict(set)
    for g in list(parent):
        comp[find(g)].add(g)
    return [c for c in comp.values() if len(c) >= 2]


def arch_concordant(ra, rb):
    """ordered intron-length vectors concordant: same n_intron (>=1) and each matched
    intron length within rel/abs tolerance (5'->3')."""
    ia, ib = ra["intron_lengths"], rb["intron_lengths"]
    if len(ia) < 1 or len(ib) < 1:
        return False
    if len(ia) != len(ib):
        return False
    for x, y in zip(ia, ib):
        if abs(x - y) > max(ABS_TOL, REL_TOL * max(x, y)):
            return False
    return True


def count_compatible(ra, rb):
    """weaker: both spliced AND intron counts comparable (within +/-1 or 15%)."""
    na, nb = ra["n_intron"], rb["n_intron"]
    if na < 1 or nb < 1:
        return False
    return abs(na - nb) <= max(1, round(0.15 * max(na, nb)))


def main():
    H, dna = dna_homology()
    idx = json.load(open(INDEX))
    base_edges = sorted(dna)
    print(f"[baseline] cDNA homology edges (id>={ID_MIN}, maxcov>={COV_MIN}): {len(base_edges)}")
    base_comp = components(base_edges)
    base_sizes = sorted((len(c) for c in base_comp), reverse=True)
    print(f"  families (components |C|>=2): {len(base_comp)}   max component: {base_sizes[0]}   "
          f"top sizes: {base_sizes[:6]}")

    # identify the mega-component (largest)
    mega = max(base_comp, key=len)
    mega_edges = {(a, b) for (a, b) in base_edges if a in mega and b in mega}
    print(f"  MEGA component: {len(mega)} genes, {len(mega_edges)} internal edges")

    # index coverage
    have_idx = sum(1 for g in mega if g in idx)
    print(f"  mega genes with an intron-index entry: {have_idx}/{len(mega)}")

    # ------- panel connectivity check helper -------
    def panel_status(edge_set):
        comps = components(edge_set)
        g2c = {}
        for i, c in enumerate(comps):
            for g in c:
                g2c[g] = i
        out = {}
        for fam, members in PANEL_REAL.items():
            present = [m for m in members if m in g2c]
            if len(present) < 2:
                out[fam] = f"N/A(<2 in graph: present={present})"
            else:
                cids = {g2c[m] for m in present}
                out[fam] = "CONNECTED" if len(cids) == 1 else f"SPLIT into {len(cids)}"
        return out

    base_panel = panel_status(set(base_edges))

    # ------- criteria ladder -------
    def edge_keep(a, b, mode):
        ra, rb = idx.get(a), idx.get(b)
        if mode == "C0":
            return True
        if ra is None or rb is None:
            return False  # no architecture info -> cannot confirm shared architecture
        if mode == "Cspliced":
            return ra["n_intron"] >= 1 and rb["n_intron"] >= 1
        if mode == "Ccount":
            return count_compatible(ra, rb)
        if mode == "Carch":
            return arch_concordant(ra, rb)
        raise ValueError(mode)

    print("\n=== CRITERIA LADDER (genome-wide on the cDNA-homology graph) ===")
    hdr = f"  {'mode':10} {'edges':>6} {'dropped':>7} {'families':>9} {'maxcomp':>8} " \
          f"{'megaEdrop':>10} {'megaKept':>8} {'panel(RABL2/APOBEC3/RFPL)':>30}"
    print(hdr)
    results = {}
    for mode in ("C0", "Cspliced", "Ccount", "Carch"):
        kept = [(a, b) for (a, b) in base_edges if edge_keep(a, b, mode)]
        comp = components(kept)
        sizes = sorted((len(c) for c in comp), reverse=True)
        maxc = sizes[0] if sizes else 0
        mega_drop = sum(1 for (a, b) in mega_edges if not edge_keep(a, b, mode))
        mega_kept = len(mega_edges) - mega_drop
        ps = panel_status(set(kept))
        pstr = "/".join(ps[f][:4] for f in ("RABL2", "APOBEC3", "RFPL"))
        print(f"  {mode:10} {len(kept):6d} {len(base_edges)-len(kept):7d} {len(comp):9d} "
              f"{maxc:8d} {mega_drop:10d} {mega_kept:8d}   {pstr}")
        results[mode] = dict(edges=len(kept), dropped=len(base_edges) - len(kept),
                             families=len(comp), maxcomp=maxc, mega_edge_drop=mega_drop,
                             mega_kept=mega_kept, panel=ps,
                             top_sizes=sizes[:8])
    print("  (panel codes: CONN=connected, SPLI=split, N/A=<2 members in homology graph)")

    # ------- does it cut the mega-family into pieces? show the new components carved out -------
    print("\n=== MEGA-FAMILY BREAKUP under Carch (architecture-concordant) ===")
    mode = "Carch"
    mega_kept_edges = [(a, b) for (a, b) in mega_edges if edge_keep(a, b, mode)]
    sub = components(mega_kept_edges)
    sub_sizes = sorted((len(c) for c in sub), reverse=True)
    singletons = len(mega) - sum(len(c) for c in sub)
    print(f"  389-mega ({len(mega)} genes, {len(mega_edges)} edges) -> after Carch: "
          f"{len(sub)} sub-families, sizes {sub_sizes[:10]}, isolated singletons={singletons}")

    # name the biggest surviving sub-clusters (root-name composition)
    import re as _re
    for c in sorted(sub, key=len, reverse=True)[:6]:
        roots = collections.Counter(_re.sub(r"[0-9]+[A-Z]?$", "", m) for m in c)
        print(f"    sub n={len(c):3d}  roots={roots.most_common(3)}  e.g. {sorted(c)[:6]}")

    # ------- retrocopy probe: do intronless cDNA-homologous edges get cut? -------
    print("\n=== RETROCOPY / INTRONLESS PROBE (cDNA id high, but architecture absent) ===")
    intronless = [g for g in mega if g in idx and idx[g]["n_intron"] == 0]
    print(f"  intronless genes inside the mega-family: {len(intronless)} "
          f"(e.g. {sorted(intronless)[:8]})")
    # how many mega edges touch an intronless gene (these are the retrocopy/domain-bridge edges)?
    il_set = set(intronless)
    il_edges = [(a, b) for (a, b) in mega_edges if a in il_set or b in il_set]
    print(f"  mega edges incident to an intronless gene: {len(il_edges)} "
          f"({100*len(il_edges)/max(len(mega_edges),1):.0f}% of mega edges) — all dropped by Cspliced/Carch")

    json.dump(dict(baseline=dict(edges=len(base_edges), families=len(base_comp),
                                 maxcomp=base_sizes[0], mega=len(mega), mega_edges=len(mega_edges)),
                   criteria=results,
                   base_panel=base_panel,
                   mega_breakup_Carch=dict(subfamilies=len(sub), sizes=sub_sizes[:10],
                                           singletons=singletons),
                   intronless_in_mega=len(intronless),
                   il_incident_edges=len(il_edges)),
              open(os.path.join(HERE, "family_def_junction_genomewide.json"), "w"), indent=2)
    print(f"\n[+] wrote {os.path.join(HERE,'family_def_junction_genomewide.json')}")


if __name__ == "__main__":
    main()
