#!/usr/bin/env python
"""BYTE-PARITY golden-fixture generator for the Rust port of the O1 refinement
operator R (genome_family_def.refine_families / _refine_component / _induced_density
/ _split_once) as deployed by the SHIPPED RNA family definition (family_rna_refine.py).

It reconstructs the EXACT gated components + kept edges that
`family_rna_refine.build_catalog` feeds to `genome_family_def.refine_families`
(the default RNA operating point: core_recip>=0.19 AND aln_frac>=0.24 AND NOT
repeat-hub, then the gamma=0.20 / seed=0 gamma-quasi-clique refinement), then for
EVERY raw connected component emits:

  * its node set (chrom/start/end/name/biotype per node, LOCAL index space),
  * the induced edges (kept-edge pairs restricted to the component),
  * a `was_split` flag  (n>2 AND induced density < gamma  => the splitter fires),
  * the Python refine output for THAT component == the gamma-quasi-clique blocks
    that pass the >=2-distinct-loci multi-copy predicate (i.e. what
    refine_families would keep from this component).

MEASURED FACT (this build): 1933 components, exactly ONE (a 213-node blob,
density < gamma) invokes the splitter; the other 1932 are returned WHOLE by the
deterministic density-gate path.  So the Rust deterministic core must byte-match
(set of member-sets) on the 1932 clean components; on the 1 blob it need only
emit a VALID witness (every family a gamma-quasi-clique with >=2 loci) -- the
splitter partition is a non-unique certificate witness (seed-dependent in Python,
deterministic-community in the Rust option-B port), so the Rust blob family SET
may differ from Python's while both satisfy the certificate.

Run: PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/gen_refine_families_fixture.py
Writes: src/rustle/vg_family/testdata/refine_families_fixture.json
Deterministic (PYTHONHASHSEED=0; shipped gamma=0.20 / seed=0; networkx Louvain seed=0 for the
Python blob reference only -- the Rust port is networkx-free / RNG-free).
"""
import os
import sys

# determinism: pin the hash seed BEFORE the heavy imports (re-exec preserves argv)
if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

import json
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import family_er_pr as FP            # noqa: E402  shipped graph context / genes dict
import genome_family_def as G        # noqa: E402  refinement operator R (the port target)
import graph_def_refine_sweep as SW  # noqa: E402  components_from_edges + SEED
import rna_only_edge_oracle as RO    # noqa: E402  RNA feature loaders
import family_rna_refine as FR       # noqa: E402  the DEFAULT RNA operating point (gates/constants)
import networkx as nx                # noqa: E402  Python blob reference ONLY
from networkx.algorithms.community import louvain_communities  # noqa: E402

OUT = os.path.join(HERE, "..", "src", "rustle", "vg_family", "testdata",
                   "refine_families_fixture.json")


def reconstruct_comps_and_kept():
    """EXACT copy of the family_rna_refine.build_catalog prefix that produces the
    (comps, kept, genes) fed to refine_families -- default gates (repeat-hub ON,
    core>=0.19, aln>=0.24).  Mirrors build_catalog verbatim so the fixture is what
    the shipped catalog actually refines (no re-derivation)."""
    FR.rna_only_guard()
    gene_of_dn = RO.load_gene_of_dn()
    pair_core = RO.load_pair_core(gene_of_dn)
    univ_aln = RO.load_universal_aln()
    pair_repeat_mult = FR.load_repeat_mult()

    meta = FP.load_meta()
    annot = FP.load_annot()
    gene_of = FP.gene_of_factory(annot)
    raw_fams = FP.load_raw_families()
    edge_pairs = FP.load_edges()
    genes, gene_of_dn2, *_ = FP.build_genes_dict(raw_fams, meta, gene_of)
    all_nodes = set()
    for f in raw_fams:
        all_nodes.update(f)

    def core_aln_keep(k):
        c = pair_core.get(k)
        c = c if c is not None else 0.0
        a = univ_aln.get(k)
        a = a if a is not None else 0.0
        return (c >= FR.CORE_MIN) and (a >= FR.ALN_MIN)

    def repeat_hub(k):
        m = pair_repeat_mult.get(k)
        return (m is not None) and (m >= FR.REPEAT_MULT_MIN)

    Gr = nx.Graph()
    Gr.add_nodes_from(all_nodes)
    for a, b in edge_pairs:
        if a in all_nodes and b in all_nodes:
            Gr.add_edge(a, b)

    kept = set()
    for u, v in Gr.edges():
        ga, gb = gene_of_dn2.get(u), gene_of_dn2.get(v)
        if ga is None or gb is None or ga == gb:
            kept.add(frozenset((u, v)))
            continue
        k = frozenset((ga, gb))
        if core_aln_keep(k) and not repeat_hub(k):
            kept.add(frozenset((u, v)))

    comps = SW.components_from_edges(all_nodes, kept)
    return comps, kept, genes


def main():
    comps, kept, genes = reconstruct_comps_and_kept()
    all_edges = [tuple(e) for e in kept]

    # adjacency EXACTLY as refine_families builds it (restricted to same-component pairs)
    fam_of = {}
    for idx, m in enumerate(comps):
        for i in m:
            fam_of[i] = idx
    adj = defaultdict(set)
    for u, v in all_edges:
        if u in fam_of and v in fam_of and fam_of[u] == fam_of[v]:
            adj[u].add(v)
            adj[v].add(u)

    def induced_density(nodes):
        n = len(nodes)
        if n <= 1:
            return 1.0
        ns = set(nodes)
        e = sum(1 for u in nodes for w in adj.get(u, ()) if w in ns and u < w)
        return 2 * e / (n * (n - 1))

    # deterministic ordering of components (by their sorted node ids) for a stable fixture
    comps_sorted = sorted((sorted(m) for m in comps), key=lambda m: (len(m), m))

    out_components = []
    n_split = 0
    n_full_refined = 0
    for m in comps_sorted:
        member_ids = list(m)                       # DN ids, sorted -> LOCAL index = position
        local = {dn: i for i, dn in enumerate(member_ids)}
        n = len(member_ids)
        was_split = (n > 2 and induced_density(member_ids) < G.GAMMA)
        if was_split:
            n_split += 1

        nodes_field = [[genes[dn]["chrom"], genes[dn]["start"], genes[dn]["end"],
                        genes[dn]["name"], genes[dn]["biotype"]] for dn in member_ids]

        # induced edges (local index pairs, i<j), deterministic order
        edge_set = set()
        for dn in member_ids:
            for w in adj.get(dn, ()):
                if w in local:
                    a, b = local[dn], local[w]
                    if a < b:
                        edge_set.add((a, b))
        edges_field = sorted(edge_set)

        # Python refine output for THIS component (>=2-distinct-loci filtered blocks)
        py_blocks = [sorted(b) for b in
                     G._refine_component(set(member_ids), adj, G.GAMMA, louvain_communities, G.SEED)
                     if G.distinct_loci(b, genes) >= 2]
        py_fams_field = sorted([sorted(local[dn] for dn in b) for b in py_blocks])
        n_full_refined += len(py_fams_field)

        out_components.append(dict(
            was_split=was_split,
            nodes=nodes_field,
            edges=edges_field,
            python_families=py_fams_field,
        ))

    fixture = dict(
        gamma=G.GAMMA,
        purity=G.PURITY,
        locus_overlap=G.LOCUS_OVERLAP,
        n_components=len(out_components),
        n_split=n_split,
        n_full_refined=n_full_refined,
        components=out_components,
    )

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(fixture, fh, separators=(",", ":"), sort_keys=True)
        fh.write("\n")

    blob_sizes = [len(c["nodes"]) for c in out_components if c["was_split"]]
    print(f"[fixture] wrote {len(out_components)} components -> {os.path.normpath(OUT)}")
    print(f"          split (density<gamma & n>2): {n_split}   blob sizes: {sorted(blob_sizes, reverse=True)}")
    print(f"          total nodes: {sum(len(c['nodes']) for c in out_components)}   "
          f"python full-refined families (>=2 loci): {n_full_refined}")
    print(f"          gamma={G.GAMMA} seed={G.SEED} locus_overlap={G.LOCUS_OVERLAP} purity={G.PURITY}")


if __name__ == "__main__":
    main()
