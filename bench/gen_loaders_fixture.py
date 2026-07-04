#!/usr/bin/env python
"""BYTE-PARITY golden-fixture generator for the Rust port of the O1 family-definition
LOADERS (`vg_family/family_loaders.rs`), the migration increment that ports the
loaders feeding the already-ported `distinct_loci` / `refine_families`.

It reuses the SHIPPED Python loaders verbatim -- `family_er_pr.load_meta /
load_annot / gene_of_factory / load_raw_families / load_edges / build_genes_dict`
and `graph_def_refine_sweep.components_from_edges` -- and emits, self-contained:

  (a) a build_genes_dict / gene_of fixture: a REPRESENTATIVE, DIVERSE subset of
      >=300 DN loci (mapped-above-floor, mapped-below-floor -> NA, no-overlap,
      missing-meta-fallback, multi-underscore contigs, and hand-built gene_of
      TIES) with the input (raw_fams / meta / a REDUCED-but-answer-EXACT annot)
      and the Python golden output (genes dict, gene_of_dn, gene_of_dn_raw,
      miss_meta, ovstats).

  (b) a components_from_edges fixture: several (all_nodes, kept_edges) cases (real
      families + synthetic: disconnected singletons, chains, cross-chrom,
      endpoint-not-in-nodes) with the Python connected components.

  (0) small raw-TSV parser fixtures (load_meta / load_annot / load_raw_families /
      load_edges) with their golden parsed structures.

ANNOT REDUCTION (compact + provably answer-exact). The real per-contig `maxlen`
is ~1.5 Mb, so a full window would be huge. Instead the fixture ships, per contig,
only the union of the OVERLAPPING intervals (`start <= qe & end > qs`) of the
selected loci. gene_of's answer depends only on overlappers, and the break-prune
`start < qs - maxlen` never skips an overlapper (an overlapper's span <= the
slice's own maxlen, so its start > qs - maxlen). The generator ASSERTS
gene_of(reduced) == gene_of(full) for every selected locus, so the fixture is
guaranteed to reflect the real full-annotation behavior.

Run: PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/gen_loaders_fixture.py
Writes: src/rustle/vg_family/testdata/loaders_fixture.json
Deterministic (sorted selection; PYTHONHASHSEED=0).
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import family_er_pr as FP                     # noqa: E402
import graph_def_refine_sweep as SW           # noqa: E402  (components_from_edges)

OUT = os.path.join(HERE, "..", "src", "rustle", "vg_family", "testdata",
                   "loaders_fixture.json")

N_ABOVE = 200     # mapped-above-floor loci to sample (deterministic, sorted)
N_NOOV = 80       # no-overlap loci to sample


# ---------------------------------------------------------------- synthetic annot
def synthetic_annot():
    """Hand-built contigs exercising gene_of TIE-BREAKS + a mapped multi-underscore
    contig. Returned as a load_annot-style idx {contig: (starts, ivals, maxlen)}.
    Each interval is (start, end, biotype, gene)."""
    contigs = {
        # TIE by overlap: LOW[500,1500) ov500 span1000 vs HIGH[1500,3000) ov500 span1500
        # for query [1000,2000). Descending scan hits HIGH first; strict `>` keeps HIGH.
        "SYNTIE": [(500, 1500, "bt_low", "GENE_LOW"), (1500, 3000, "bt_high", "GENE_HIGH")],
        # TIE with EQUAL (start,end): AAA vs BBB for query [1000,2000); tuple sort puts
        # AAA before BBB, descending scan hits BBB first -> BBB wins.
        "SYNTIE2": [(900, 1600, "bt_a", "GENE_AAA"), (900, 1600, "bt_b", "GENE_BBB")],
        # a plain mapped multi-underscore contig (meta-backed) above the floor
        "SYN_multi_us": [(0, 4000, "protein_coding", "SYN_GENE1")],
    }
    idx = {}
    for c, rows in contigs.items():
        ivals = sorted((int(s), int(e), bt, g) for (s, e, bt, g) in rows)
        starts = [x[0] for x in ivals]
        maxlen = max((e - s) for s, e, _, _ in ivals)
        idx[c] = (starts, ivals, maxlen)
    return idx


def synthetic_loci():
    """(dn, in_meta, meta_or_None, category). Synthetic loci whose golden is computed
    by the SAME shipped build_genes_dict below."""
    return [
        # gene_of ties (meta-backed, mapped)
        ("DN_SYNTIE_1000_5", True, ("SYNTIE", 1000, 2000), "tie"),
        ("DN_SYNTIE2_1000_4", True, ("SYNTIE2", 1000, 2000), "tie"),
        # mapped multi-underscore contig (meta-backed)
        ("DN_SYN_multi_us_100_3", True, ("SYN_multi_us", 100, 3900), "above"),
        # missing-meta fallback with a MULTI-UNDERSCORE contig (parsed chrom = SYN_ctg_x)
        ("DN_SYN_ctg_x_5000_3", False, None, "missing_meta"),
        # missing-meta fallback, simpler contig
        ("DN_SYNZ_777_2", False, None, "missing_meta"),
    ]


# ---------------------------------------------------------------- main
def main():
    meta_full = FP.load_meta()
    annot_full = FP.load_annot()          # {contig: (starts, ivals, maxlen)}
    gene_of_full_real = FP.gene_of_factory(annot_full)
    raw_all = FP.load_raw_families()
    edges_all = FP.load_edges()

    # ---- classify every distinct real DN member into a category ----
    above, below, noov, missing = [], [], [], []
    seen = set()
    for fam in raw_all:
        for dn in fam:
            if dn in seen:
                continue
            seen.add(dn)
            mm = meta_full.get(dn)
            if mm is None:
                missing.append(dn)
                continue
            c, s, e = mm
            g, bt, ovf = gene_of_full_real(c, s, e)
            if g is None:
                noov.append(dn)
            elif ovf >= FP.OV_FLOOR:
                above.append(dn)
            else:
                below.append(dn)
    above.sort(); below.sort(); noov.sort(); missing.sort()

    sel_above = above[:N_ABOVE]
    sel_below = below[:]                # all below-floor (small)
    sel_noov = noov[:N_NOOV]
    sel_missing = missing[:]            # all missing-meta (small)
    real_selected = sel_above + sel_below + sel_noov + sel_missing

    # ---- combined FULL annot (real + synthetic) so gene_of works for all loci ----
    syn_annot = synthetic_annot()
    combined_full = dict(annot_full)
    combined_full.update(syn_annot)
    gene_of_full = FP.gene_of_factory(combined_full)

    # ---- curated meta (subset) + synthetic meta ----
    curated_meta = {}
    for dn in real_selected:
        if dn in meta_full:
            curated_meta[dn] = meta_full[dn]
    syn = synthetic_loci()
    for dn, in_meta, mm, _cat in syn:
        if in_meta:
            curated_meta[dn] = mm

    # ---- all query DN ids + their resolved coords (mirrors build_genes_dict) ----
    all_dns = real_selected + [dn for dn, *_ in syn]

    def resolve_coords(dn):
        mm = curated_meta.get(dn)
        if mm is not None:
            return mm
        parts = dn.split("_")
        start = int(parts[-2]); chrom = "_".join(parts[1:-2]); end = start
        return (chrom, start, end)

    # ---- REDUCED annot = per-contig union of overlapping intervals ----
    reduced_by_c = {}   # contig -> set of interval tuples
    for dn in all_dns:
        chrom, qs, qe = resolve_coords(dn)
        rec = combined_full.get(chrom)
        if rec is None:
            continue
        _starts, ivals, _maxlen = rec
        for iv in ivals:
            if iv[0] <= qe and iv[1] > qs:
                reduced_by_c.setdefault(chrom, set()).add(iv)
    # build reduced idx {contig:(starts,ivals,maxlen)} from the union
    reduced_idx = {}
    for c, ivset in reduced_by_c.items():
        ivals = sorted(ivset)
        starts = [x[0] for x in ivals]
        maxlen = max((e - s) for s, e, _, _ in ivals)
        reduced_idx[c] = (starts, ivals, maxlen)
    gene_of_reduced = FP.gene_of_factory(reduced_idx)

    # ---- VERIFY reduced == full for every selected locus (fixture faithfulness) ----
    n_verify = 0
    for dn in all_dns:
        chrom, qs, qe = resolve_coords(dn)
        assert gene_of_full(chrom, qs, qe) == gene_of_reduced(chrom, qs, qe), \
            f"reduced-annot gene_of diverges from full for {dn}"
        n_verify += 1

    # ---- curated raw_fams (partition + a DEDUP duplicate + synthetic families) ----
    # Chunk the real selection into families of ~40; append synthetic loci; then a
    # family that RE-LISTS an earlier dn (tests build_genes_dict first-occurrence dedup).
    raw_fams = []
    chunk = 40
    for i in range(0, len(real_selected), chunk):
        raw_fams.append(real_selected[i:i + chunk])
    raw_fams.append([dn for dn, *_ in syn])
    if real_selected:
        raw_fams.append([real_selected[0], real_selected[-1]])   # duplicates -> deduped

    # ---- GOLDEN via shipped build_genes_dict on the reduced annot ----
    genes, gene_of_dn, gene_of_dn_raw, miss_meta, ovstats = \
        FP.build_genes_dict(raw_fams, curated_meta, gene_of_reduced)
    # sanity: identical when run on the FULL annot (reduction is answer-exact)
    genes_full, god_full, godr_full, mm_full, ov_full = \
        FP.build_genes_dict(raw_fams, curated_meta, gene_of_full)
    assert genes == genes_full and gene_of_dn == god_full and gene_of_dn_raw == godr_full \
        and miss_meta == mm_full and ovstats == ov_full, "reduced vs full build_genes_dict diverge"

    # ---- diversity accounting (over the distinct genes dict) ----
    div = dict(mapped_above_floor=0, mapped_below_floor_na=0, no_overlap=0,
               missing_meta=0, ties=len([1 for _, _, _, cat in syn if cat == "tie"]))
    for dn, gm in genes.items():
        in_meta = dn in curated_meta
        if not in_meta:
            div["missing_meta"] += 1
            continue
        if dn in gene_of_dn:
            div["mapped_above_floor"] += 1
        elif dn in gene_of_dn_raw:
            div["mapped_below_floor_na"] += 1
        else:
            div["no_overlap"] += 1

    # ---- emit reduced annot as {contig: [[s,e,bt,g],...]} (Rust re-sorts) ----
    annot_json = {c: [[s, e, bt, g] for (s, e, bt, g) in ivals]
                  for c, (_starts, ivals, _ml) in reduced_idx.items()}

    build_genes = dict(
        raw_fams=raw_fams,
        meta={dn: [c, s, e] for dn, (c, s, e) in curated_meta.items()},
        annot=annot_json,
        golden=dict(
            genes={dn: dict(chrom=g["chrom"], start=g["start"], end=g["end"],
                            biotype=g["biotype"], name=g["name"], ov=g["ov"])
                   for dn, g in genes.items()},
            gene_of_dn=gene_of_dn,
            gene_of_dn_raw=gene_of_dn_raw,
            miss_meta=miss_meta,
            ovstats=dict(mapped_raw=ovstats["mapped_raw"],
                         mapped_floored=ovstats["mapped_floored"],
                         dropped_by_floor=ovstats["dropped_by_floor"],
                         lt20_overlap=ovstats["lt20_overlap"],
                         lt50_overlap=ovstats["lt50_overlap"]),
            diversity=div,
        ),
    )

    # ================= (b) components_from_edges cases =================
    # index the real DN edges for restriction to a node subset
    edge_pairs = [(a, b) for (a, b) in edges_all]

    def within(nodes):
        ns = set(nodes)
        kept = set()
        for a, b in edge_pairs:
            if a in ns and b in ns:
                kept.add(frozenset((a, b)))
        return kept

    def case(label, all_nodes, kept_edges):
        comps = SW.components_from_edges(all_nodes, kept_edges)
        return dict(label=label,
                    all_nodes=list(all_nodes),
                    kept_edges=[list(tuple(e)) for e in kept_edges],
                    golden=[sorted(c) for c in comps])

    components = []
    # real: a single family + its internal edges
    fam1 = sorted(set(raw_all[1]))                      # DNFAM1 members (28)
    components.append(case("real_family1", fam1, within(fam1)))
    # real: two families' union + internal edges (expect >=2 components)
    fam2 = sorted(set(raw_all[2]))
    union12 = sorted(set(fam1) | set(fam2))
    components.append(case("real_family1_2_union", union12, within(union12)))
    # synthetic: disconnected singletons (no edges)
    components.append(case("syn_singletons", ["n0", "n1", "n2"], set()))
    # synthetic: one chain
    chain = ["c0", "c1", "c2", "c3", "c4"]
    chain_e = {frozenset((chain[i], chain[i + 1])) for i in range(len(chain) - 1)}
    components.append(case("syn_chain", chain, chain_e))
    # synthetic: two separate chains -> 2 components
    twoc_nodes = ["a0", "a1", "a2", "b0", "b1"]
    twoc_e = {frozenset(("a0", "a1")), frozenset(("a1", "a2")), frozenset(("b0", "b1"))}
    components.append(case("syn_two_chains", twoc_nodes, twoc_e))
    # synthetic: cross-chrom-like labels, partial connectivity + isolated node
    xnodes = ["chrA_1", "chrA_2", "chrB_1", "chrB_2", "chrC_1"]
    xe = {frozenset(("chrA_1", "chrB_2")), frozenset(("chrA_2", "chrB_1"))}
    components.append(case("syn_cross_chrom", xnodes, xe))
    # synthetic: an edge endpoint NOT in all_nodes (networkx adds it)
    components.append(case("syn_edge_adds_endpoint", ["p"],
                           {frozenset(("p", "q")), frozenset(("r", "s"))}))

    # ================= (0) raw-TSV parser fixtures =================
    # small, hand-built TSV snippets + their golden parsed structures
    meta_tsv = ("id\tchrom\tstart\tend\tstrand\tn_exon\tn_reads\n"
                "DN_NC_1_100_5\tNC_1\t100\t2000\t+\t3\t50\n"
                "DN_ctg_a_b_500_2\tctg_a_b\t500\t900\t-\t1\t9\n")
    meta_golden = {"DN_NC_1_100_5": ["NC_1", 100, 2000],
                   "DN_ctg_a_b_500_2": ["ctg_a_b", 500, 900]}
    annot_tsv = ("chrom\tstart\tend\tbiotype\tgene\n"
                 "K\t300\t400\tlncRNA\tGB\n"
                 "K\t100\t250\tprotein_coding\tGA\n"
                 "K\t100\t250\tprotein_coding\tGA0\n"      # equal (start,end); tuple sort GA<GA0
                 "L\t0\t9000\tprotein_coding\tGBIG\n")
    # build annot golden by mirroring the SHIPPED load_annot (sort + starts + maxlen)
    def parse_annot_text(text):
        from collections import defaultdict
        by_c = defaultdict(list)
        for ln in text.splitlines()[1:]:
            f = ln.split("\t")
            by_c[f[0]].append((int(f[1]), int(f[2]), f[3], f[4]))
        idx = {}
        for c, ivals in by_c.items():
            ivals.sort()
            idx[c] = dict(starts=[x[0] for x in ivals],
                          maxlen=max((e - s) for s, e, _, _ in ivals),
                          ivals=[[s, e, bt, g] for (s, e, bt, g) in ivals])
        return idx
    annot_golden = parse_annot_text(annot_tsv)
    raw_fams_tsv = ("family_id\tn_copies\tmembers\n"
                    "F0\t3\tDN_a,DN_b,DN_c\n"
                    "F1\t2\tDN_d,DN_e\n")
    raw_fams_golden = [["DN_a", "DN_b", "DN_c"], ["DN_d", "DN_e"]]
    edges_tsv = ("a\tb\tcore_recip\n"
                 "DN_a\tDN_b\t0.5\n"
                 "DN_b\tDN_c\t0.31\n")
    edges_golden = [["DN_a", "DN_b"], ["DN_b", "DN_c"]]

    raw_parsers = dict(
        meta_tsv=meta_tsv, meta_golden=meta_golden,
        annot_tsv=annot_tsv, annot_golden=annot_golden,
        raw_fams_tsv=raw_fams_tsv, raw_fams_golden=raw_fams_golden,
        edges_tsv=edges_tsv, edges_golden=edges_golden,
    )

    fixture = dict(ov_floor=FP.OV_FLOOR,
                   build_genes=build_genes,
                   components=components,
                   raw_parsers=raw_parsers)

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(fixture, fh, indent=1, sort_keys=True)

    n_loci = len(genes)
    n_annot_ivals = sum(len(v) for v in annot_json.values())
    print(f"[fixture] wrote {os.path.normpath(OUT)}")
    print(f"          build_genes: {n_loci} distinct DN loci  "
          f"(above={div['mapped_above_floor']} below-NA={div['mapped_below_floor_na']} "
          f"no-ov={div['no_overlap']} miss-meta={div['missing_meta']} ties={div['ties']})")
    print(f"          reduced annot: {len(annot_json)} contigs / {n_annot_ivals} intervals "
          f"(verified gene_of(reduced)==gene_of(full) for {n_verify} loci)")
    print(f"          ovstats: {build_genes['golden']['ovstats']}  miss_meta={miss_meta}")
    print(f"          components: {len(components)} cases; raw-parser fixtures: 4")
    assert n_loci >= 300, f"only {n_loci} distinct loci (< 300)"


if __name__ == "__main__":
    main()
