#!/usr/bin/env python3
"""family_def_retrocopy_filter.py — SHIPPED precision filter for cDNA-homology family edges.

A processed RETROCOPY (intronless copy of a spliced parent) is cDNA-homologous to its parent
(id ~0.95-1.0) but is NOT a transcribed tandem family member -- sequence identity cannot flag it,
but the splice structure can: the copy has ZERO introns while the parent is spliced. This is a
threshold-free, sequence-orthogonal, VG-native exclusion (a 0-intron path cannot thread a
multi-exon junction backbone).

CRITICAL (verified bench/family_def_vg_junction work, 2026-06-23): the valid pattern is
ONE-SIDE-intronless-vs-spliced ONLY. BOTH-intronless edges (3267 genome-wide) are NOT retrocopies
-- they contain real single-exon TANDEM families (IFNA interferon-alpha cluster, PCDHB
protocadherin-beta, ELOA, HNRNPCL). Dropping both-intronless would destroy real families at ~10:1
false:true. So we drop an edge iff exactly one endpoint is intronless and the other is spliced.

Genome-wide on the id>=0.90 & maxcov>=0.30 homology graph: 1931 one-side-intronless edges dropped
(retrocopies), 3267 both-intronless KEPT (real intronless tandems), 12212 both-spliced KEPT.

API:  is_retrocopy_edge(a, b, idx) ;  partition_edges(edges, idx) ;  build_family_graph(...)
Run (self-test + provenance report): python bench/family_def_retrocopy_filter.py
"""
import collections
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

INDEX = "/home/juanfra/winloci_scratch/gene_intron_index.json"
GENES_BED = "/home/juanfra/winloci_scratch/unmapped_poc/genes.bed"
ID_MIN = 0.90
COV_MIN = 0.30
# opt-out so a caller can reproduce the pre-filter (over-merged) truth byte-for-byte
DISABLE = os.environ.get("FAMILY_DEF_NO_RETROCOPY_FILTER") == "1"


def load_intron_index(path=INDEX):
    return json.load(open(path))


def n_intron(idx, g):
    """Intron count for gene g, or None if g has no annotated structure in the index."""
    v = idx.get(g)
    if v is None:
        return None
    return v.get("n_intron")


def is_retrocopy_edge(a, b, idx):
    """True iff (a,b) is a processed-retrocopy edge: exactly ONE endpoint is intronless and the
    other is spliced (>=1 intron). Both-intronless and both-spliced -> False. Unknown structure
    on either side -> False (cannot confirm the retrocopy pattern; keep the edge, stay conservative).
    Assumes the caller has already applied the homology gate (id>=ID_MIN & maxcov>=COV_MIN)."""
    na, nb = n_intron(idx, a), n_intron(idx, b)
    if na is None or nb is None:
        return False
    return (na == 0) ^ (nb == 0)        # XOR: exactly one side intronless


def partition_edges(edges, idx):
    """Split homology edges into kept vs retrocopy-dropped, with the diagnostic sub-classes.
    `edges` = iterable of (a,b) gene-name pairs. Returns a dict of edge lists."""
    out = collections.defaultdict(list)
    for a, b in edges:
        na, nb = n_intron(idx, a), n_intron(idx, b)
        if na is None or nb is None:
            out["unknown_kept"].append((a, b))
        elif na == 0 and nb == 0:
            out["both_intronless_kept"].append((a, b))   # real single-exon tandems -> KEEP
        elif (na == 0) ^ (nb == 0):
            out["retrocopy_dropped"].append((a, b))       # one-side-intronless -> DROP
        else:
            out["both_spliced_kept"].append((a, b))
    out["kept"] = (out["both_spliced_kept"] + out["both_intronless_kept"] + out["unknown_kept"])
    return out


def homology_edges(Hd, id_min=ID_MIN, cov_min=COV_MIN):
    return [(a, b) for (a, b), r in Hd.items()
            if r.get("id", 0) >= id_min and max(r["cov_a"], r["cov_b"]) >= cov_min]


def build_family_graph(Hd=None, idx=None, filter_retrocopy=None):
    """Build the cDNA-homology family graph (networkx) with the retrocopy filter applied by
    default. Returns (graph, provenance_dict). Set filter_retrocopy=False (or env
    FAMILY_DEF_NO_RETROCOPY_FILTER=1) to reproduce the pre-filter graph."""
    import networkx as nx
    if Hd is None:
        from family_def_read_filters import dna_homology
        Hd, _ = dna_homology()
    if idx is None:
        idx = load_intron_index()
    if filter_retrocopy is None:
        filter_retrocopy = not DISABLE
    edges = homology_edges(Hd)
    part = partition_edges(edges, idx)
    use = part["kept"] if filter_retrocopy else edges
    G = nx.Graph()
    G.add_edges_from(use)
    prov = dict(filter_retrocopy=filter_retrocopy, total_homology_edges=len(edges),
                retrocopy_dropped=len(part["retrocopy_dropped"]),
                both_intronless_kept=len(part["both_intronless_kept"]),
                both_spliced_kept=len(part["both_spliced_kept"]),
                unknown_kept=len(part["unknown_kept"]))
    return G, prov, part


def _selftest(Hd, idx):
    """Guard the load-bearing invariants with concrete known cases."""
    ok = True
    def check(name, cond):
        nonlocal ok
        print(f"    [{'PASS' if cond else 'FAIL'}] {name}")
        ok = ok and cond
    # one-side-intronless retrocopy probes (parent spliced, LOC copy intronless)
    for parent, copy in [("EEF1A1", "LOC109023808"), ("CNN2", "LOC129524764")]:
        if parent in idx and copy in idx:
            np_, nc = n_intron(idx, parent), n_intron(idx, copy)
            check(f"retrocopy {parent}({np_})~{copy}({nc}) flagged",
                  is_retrocopy_edge(parent, copy, idx) and np_ >= 1 and nc == 0)
    # real spliced tandem family must NOT be flagged
    if "RABL2A" in idx and "RABL2B" in idx:
        check("real family RABL2A~RABL2B NOT flagged", not is_retrocopy_edge("RABL2A", "RABL2B", idx))
    # XOR semantics: both-intronless NOT flagged (synthetic)
    fake = {"x0": {"n_intron": 0}, "y0": {"n_intron": 0}, "z3": {"n_intron": 3}}
    check("both-intronless x0~y0 NOT flagged", not is_retrocopy_edge("x0", "y0", fake))
    check("one-side x0~z3 flagged", is_retrocopy_edge("x0", "z3", fake))
    check("unknown side kept (None)", not is_retrocopy_edge("x0", "missing", fake))
    return ok


def main():
    from family_def_read_filters import dna_homology
    Hd, _ = dna_homology()
    idx = load_intron_index()
    print("=== retrocopy filter self-test ===")
    if not _selftest(Hd, idx):
        print("  SELF-TEST FAILED"); sys.exit(1)

    G_pre, prov_pre, _ = build_family_graph(Hd, idx, filter_retrocopy=False)
    G_post, prov_post, part = build_family_graph(Hd, idx, filter_retrocopy=True)
    import networkx as nx
    def fam_stats(G):
        comps = [c for c in nx.connected_components(G) if len(c) >= 2]
        sizes = sorted((len(c) for c in comps), reverse=True)
        return len(comps), (sizes[0] if sizes else 0), sizes[:6]

    print("\n=== edge partition (id>=0.90 & maxcov>=0.30 homology graph) ===")
    print(f"  total homology edges      : {prov_post['total_homology_edges']}")
    print(f"  DROPPED one-side-intronless: {prov_post['retrocopy_dropped']}  (processed retrocopies)")
    print(f"  KEPT  both-intronless      : {prov_post['both_intronless_kept']}  (real single-exon tandems: IFNA/PCDHB/...)")
    print(f"  KEPT  both-spliced         : {prov_post['both_spliced_kept']}")
    print(f"  KEPT  unknown-structure    : {prov_post['unknown_kept']}")

    npre, mpre, spre = fam_stats(G_pre)
    npost, mpost, spost = fam_stats(G_post)
    print("\n=== family-graph effect ===")
    print(f"  pre-filter : {npre} families, max comp {mpre}, top {spre}")
    print(f"  post-filter: {npost} families, max comp {mpost}, top {spost}")

    # panel real families must remain connected (the filter is precision-only, no real-family loss)
    panel = {"RABL2": ["RABL2A", "RABL2B"], "APOBEC3": ["APOBEC3D", "APOBEC3F"],
             "RFPL": ["RFPL1", "RFPL2", "RFPL3"]}
    print("\n=== panel real families still connected post-filter? ===")
    allok = True
    for fam, members in panel.items():
        present = [m for m in members if m in G_post]
        if len(present) < 2:
            print(f"  {fam:10}: N/A (<2 in graph)"); continue
        comp = nx.node_connected_component(G_post, present[0])
        conn = all(m in comp for m in present)
        allok = allok and conn
        print(f"  {fam:10}: {'CONNECTED' if conn else 'SPLIT (REGRESSION!)'}")

    # examples of dropped retrocopy edges (audit trail)
    print("\n  sample dropped retrocopy edges (parent_introns -> copy_introns):")
    for a, b in sorted(part["retrocopy_dropped"])[:6]:
        na, nb = n_intron(idx, a), n_intron(idx, b)
        spliced, intronless = (a, b) if nb == 0 else (b, a)
        print(f"    {spliced}({max(na,nb)} introns) -> {intronless}(0)  id={Hd[(a,b) if a<b else (b,a)]['id']:.3f}")

    json.dump(dict(provenance=prov_post,
                   families_pre=npre, families_post=npost, maxcomp_pre=mpre, maxcomp_post=mpost,
                   panel_all_connected=allok),
              open(os.path.join(HERE, "family_def_retrocopy_filter.json"), "w"), indent=2)
    print(f"\n[+] wrote {os.path.join(HERE,'family_def_retrocopy_filter.json')}")
    if not allok:
        print("  WARNING: a panel real family split -- filter is NOT precision-only here"); sys.exit(1)


if __name__ == "__main__":
    main()
