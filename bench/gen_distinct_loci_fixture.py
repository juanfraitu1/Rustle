#!/usr/bin/env python
"""BYTE-PARITY golden-fixture generator for the Rust port of the O1 multi-copy
predicate `genome_family_def.distinct_loci`.

It reuses the SHIPPED build verbatim -- `genome_family_def.load_genes` for the
`genes` dict and the shipped edge-emission + single-linkage + `refine_families`
(gamma-quasi-clique) pipeline for real member-index families -- then, for every
family, emits the member coords and the Python `distinct_loci` golden count.
Synthetic edge-case families (singleton, duplicate coords, nested/reciprocal
boundaries, zero-length, cross-chrom, transitive chain) are appended so the Rust
port is exercised on the corners the real catalog may not hit.

The golden count is `genome_family_def.distinct_loci(members, genes)` computed on
the EXACT coords written to the row -- so the fixture is internally consistent by
construction and the Rust port has a single source of truth to match.

Run: PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/gen_distinct_loci_fixture.py
Writes: src/rustle/vg_family/testdata/distinct_loci_fixture.tsv
Deterministic (PYTHONHASHSEED=0; shipped seed=0 / gamma=0.20).
"""
import os
import sys
from collections import defaultdict
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import genome_family_def as G  # noqa: E402

OUT = os.path.join(HERE, "..", "src", "rustle", "vg_family", "testdata",
                   "distinct_loci_fixture.tsv")


def build_real_families():
    """Shipped pipeline -> (genes, [family as sorted member-idx list]).

    Mirrors genome_family_def.main(): load_genes -> load_sedef_pairs (RAW default
    oracle) -> parallel gene-gene edge emission -> single-linkage raw components
    (famA all-biotype, famPC protein-coding) -> refine_families (gamma=0.20,
    seed=0). Returns genes plus the deduped union of refined + raw families, which
    spans small/large, single-/cross-chrom, and same-locus-collapse cases."""
    genes, genes_by = G.load_genes(G.GFF)
    gstarts = {c: [x[0] for x in v] for c, v in genes_by.items()}
    pairs, _n_total, _n_kept = G.load_sedef_pairs(G.SEDEF, bailey=False)
    pairs_by_cA = defaultdict(list)
    for p in pairs:
        pairs_by_cA[p[0]].append(p)
    shards = sorted(pairs_by_cA.items(), key=lambda kv: -len(kv[1]))
    with Pool(8, initializer=G._init_worker, initargs=(genes_by, gstarts)) as pool:
        edge_lists = pool.map(G._edges_for_contig, shards)
    all_edges = [(u, v) for el in edge_lists for (u, v) in el]

    def build(edges, keep=None):
        uf = G.UF()
        for u, v in edges:
            if keep is None or (keep(u) and keep(v)):
                uf.union(u, v)
        fams = [sorted(set(m)) for m in uf.groups().values() if len(set(m)) >= 2]
        fams.sort(key=lambda m: (-len(m), genes[m[0]]["name"]))
        return fams

    famA = build(all_edges)
    famPC = build(all_edges, keep=lambda i: genes[i]["biotype"] == "protein_coding")
    refined = G.refine_families(famA, all_edges, genes, G.GAMMA, seed=G.SEED)
    return genes, famA, famPC, refined


def synthetic_families():
    """Hand-built edge cases the real catalog under-covers. Each is a list of
    (chrom, start, end) members; the golden count is computed by the SAME shipped
    distinct_loci below, so these assert the exact reciprocal-overlap / break /
    union-find / no-coord-singleton semantics."""
    return [
        # (label, [(chrom, start, end), ...])
        ("syn_singleton",             [("SYN1", 100, 200)]),                         # 1 member -> 1 locus (uf.find singleton)
        ("syn_two_distinct",          [("SYN1", 100, 200), ("SYN1", 500, 600)]),     # disjoint same chrom -> 2
        ("syn_dup_coords",            [("SYN1", 100, 200), ("SYN1", 100, 200)]),     # identical -> full overlap -> 1
        ("syn_dup_coords_x3",         [("SYN1", 100, 200)] * 3),                     # three identical -> 1
        ("syn_nested_ge50_both",      [("SYN1", 100, 200), ("SYN1", 110, 210)]),     # ov=90; 90>=50 both -> 1
        ("syn_partial_lt50",          [("SYN1", 100, 200), ("SYN1", 160, 260)]),     # ov=40; 40<50 -> 2
        ("syn_exact_50pct",           [("SYN1", 100, 200), ("SYN1", 150, 250)]),     # ov=50; 50>=50.0 both -> 1
        ("syn_just_under_50pct",      [("SYN1", 100, 200), ("SYN1", 151, 251)]),     # ov=49; 49<50 -> 2
        ("syn_nested_asym_fail",      [("SYN1", 100, 1100), ("SYN1", 100, 200)]),    # inner 100% but outer 10% -> reciprocal FAIL -> 2
        ("syn_nested_asym_pass",      [("SYN1", 100, 300), ("SYN1", 100, 200)]),     # inner len100 ov100>=50; outer len200 ov100>=100 -> 1
        ("syn_diff_chrom_overlap",    [("SYN1", 100, 200), ("SYN2", 100, 200)]),     # identical coords, diff chrom -> never collapse -> 2
        ("syn_zero_len_touch",        [("SYN1", 100, 100), ("SYN1", 90, 110)]),      # zero-width member: ov=0 (not >0) -> 2
        ("syn_zero_len_pair",         [("SYN1", 100, 100), ("SYN1", 100, 100)]),     # two zero-width same point: ov=0 -> 2
        ("syn_chain_transitive",      [("SYN1", 0, 100), ("SYN1", 30, 130), ("SYN1", 60, 160)]),  # A~B ov70, B~C ov70, A!~C ov40 -> union all -> 1
        ("syn_break_after_gap",       [("SYN1", 0, 100), ("SYN1", 30, 130), ("SYN1", 500, 600)]),  # {0,1} ov70 collapse, third apart -> 2
        ("syn_mixed_two_chrom",       [("SYN1", 0, 100), ("SYN1", 30, 130), ("SYN2", 0, 100)]),    # SYN1 pair (ov70) collapses, SYN2 alone -> 2
        ("syn_five_some_collapse",    [("SYN1", 0, 100), ("SYN1", 30, 130), ("SYN1", 400, 500),
                                       ("SYN1", 430, 530), ("SYN1", 900, 1000)]),    # {0,1}ov70->1 {2,3}ov70->1 {4}->1 -> 3
        ("syn_out_of_order_starts",   [("SYN1", 500, 600), ("SYN1", 100, 200), ("SYN1", 130, 230)]),  # unsorted input; {1,2} ov70 collapse -> 2
        ("syn_all_distinct_chroms",   [("SYN1", 0, 100), ("SYN2", 0, 100), ("SYN3", 0, 100)]),   # 3 chroms -> 3
    ]


def coords_field(members):
    return ";".join(f"{c}|{s}|{e}" for (c, s, e) in members)


def main():
    genes, famA, famPC, refined = build_real_families()

    # Dedupe families by their coord multiset (order-independent for distinct_loci).
    # `members_of` maps a family (member-idx list) -> its coord tuples.
    def coords_of(fam):
        return [(genes[i]["chrom"], genes[i]["start"], genes[i]["end"]) for i in fam]

    rows = []          # (label, members[list of (chrom,start,end)])
    seen = set()

    def add(label, members):
        key = tuple(sorted(members))
        if key in seen:
            return False
        seen.add(key)
        rows.append((label, members))
        return True

    # Canonical shipped refined catalog first (multi-copy families, >=2 loci).
    for k, fam in enumerate(refined):
        add(f"refined_{k}", coords_of(fam))
    # Raw single-linkage families (adds the large cross-chrom blob + extra collapse).
    for k, fam in enumerate(famA):
        add(f"rawA_{k}", coords_of(fam))
    # Protein-coding raw families (extra granularity / collapse variety).
    for k, fam in enumerate(famPC):
        add(f"rawPC_{k}", coords_of(fam))
    # Synthetic edge cases (corners the real catalog under-covers).
    for label, members in synthetic_families():
        add(label, members)

    # Emit. Golden = shipped distinct_loci over an index-space rebuilt from the row
    # coords (so the fixture is self-consistent regardless of family provenance).
    n_collapse = n_xchrom = n_single_chrom = n_singleton = 0
    max_n = 0
    with open(OUT, "w") as fh:
        fh.write("# byte-parity golden fixture for genome_family_def.distinct_loci "
                 "(LOCUS_OVERLAP=0.50). Generated by bench/gen_distinct_loci_fixture.py.\n")
        fh.write("# columns: label<TAB>n_members<TAB>distinct_loci<TAB>coords"
                 "  (coords = chrom|start|end;chrom|start|end;...)\n")
        for label, members in rows:
            local = [dict(chrom=c, start=s, end=e) for (c, s, e) in members]
            golden = G.distinct_loci(list(range(len(local))), local)
            fh.write(f"{label}\t{len(members)}\t{golden}\t{coords_field(members)}\n")
            # diversity accounting
            nch = len({c for (c, s, e) in members})
            if golden < len(members):
                n_collapse += 1
            if nch >= 2:
                n_xchrom += 1
            else:
                n_single_chrom += 1
            if len(members) == 1:
                n_singleton += 1
            max_n = max(max_n, len(members))

    print(f"[fixture] wrote {len(rows)} families -> {os.path.normpath(OUT)}")
    print(f"          collapse (distinct_loci < n_members): {n_collapse}")
    print(f"          cross-chrom (>=2 contigs): {n_xchrom}   single-chrom: {n_single_chrom}")
    print(f"          singletons (n=1): {n_singleton}   largest family: {max_n} members")


if __name__ == "__main__":
    main()
