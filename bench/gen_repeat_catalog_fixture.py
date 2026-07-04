#!/usr/bin/env python
"""BYTE-PARITY golden-fixture generator for the Rust port of the CATALOG-PRODUCTION
core of `bench/vg_repeat_catalog.py` (`vg_family/repeat_catalog.rs`) -- the O1
over-merge-gate migration (RETIREMENT_AND_MIGRATION.md Tier-2 #7, full-port path).

It imports the SHIPPED Python (`vg_repeat_catalog` for `minimizers` / `dn_exons` /
`load_skeletons` / K / W / M_OP, and `family_er_pr` for the meta + gene_of_dn
projection) and runs the REAL functions to emit a SELF-CONTAINED fixture to
`src/rustle/vg_family/testdata/repeat_catalog_fixture.json` with three parts:

  (a) "dn_exons":  >= 40 cases -- each ships the meta triple + skeleton introns
      (self-contained) and the REAL `dn_exons` output. Real multi-intron loci +
      hand-built edge cases (single-exon via blank introns, missing skeleton ->
      single exon, meta-miss -> None, leading/trailing/whole-span/unsorted/adjacent
      introns).
  (b) per-locus NODE SET: the first >= 30 "loci" carry `node_set` =
      sorted(set(v for v,_ in minimizers(exon_seq))) on the REAL soft-masked bytes.
  (c) "loci" + "node_section_text": a CLOSED-BY-CONSTRUCTION universe (the first
      CATALOG_N sorted gene-assigned loci). The node section (`# SECTION nodes
      (mult>=2)`, the exact block `vg_repeat_catalog.py:599-604` writes) is computed
      over EXACTLY these loci via the real build loop + writer logic, and shipped
      verbatim. Each locus also carries its (chrom, exons) so the Rust test can
      cross-check the case-carrying FASTA reader against pysam when the genome exists.

Run: PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/gen_repeat_catalog_fixture.py
Deterministic (sorted loci, order-independent accumulators, PYTHONHASHSEED=0).
"""
import os
import sys

# Deterministic hashing regardless of the caller's env (mirrors gen_minimizers_fixture.py).
if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

import json
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import vg_repeat_catalog as V  # SHIPPED minimizers / dn_exons / load_skeletons / K / W / M_OP
import family_er_pr as F       # SHIPPED meta + gene_of_dn projection

OUT = os.path.join(HERE, "..", "src", "rustle", "vg_family", "testdata",
                   "repeat_catalog_fixture.json")

CATALOG_N = 60   # closed universe size for part (c); yields both is_repeat paths
NODESET_N = 34   # part (b): loci carrying node_set (>= 30 required by the Rust test)
DNEXON_REAL = 42 # part (a): real loci covering multi-intron structure (>= 40 w/ synthetics)


def main():
    print("[load] meta / annot / gene projection / skeletons ...", flush=True)
    meta = F.load_meta()
    annot = F.load_annot()
    gene_of = F.gene_of_factory(annot)
    raw_fams = F.load_raw_families()
    _, gene_of_dn, _, _, _ = F.build_genes_dict(raw_fams, meta, gene_of)
    skel = V.load_skeletons()
    loci = sorted(gene_of_dn.keys())
    print(f"       {len(loci)} gene-assigned loci, {len(set(gene_of_dn.values()))} genes,"
          f" {len(skel)} skeletons", flush=True)

    import pysam
    fa = pysam.FastaFile(V.GENOME)

    # ---------------------------------------------------------------- part (a) dn_exons
    dn_cases = []

    def emit_dn(dn, meta_map, skel_map):
        """Run the REAL dn_exons; record its (chrom, exons) + the inputs that produced it."""
        chrom, exons = V.dn_exons(dn, meta_map, skel_map)
        mm = meta_map.get(dn)
        intr = None
        if mm is not None:
            intr = skel_map.get(mm)  # None => missing skeleton branch
        dn_cases.append(dict(
            dn=dn,
            meta=(list(mm) if mm is not None else None),
            introns=intr,
            chrom=chrom,
            exons=(None if exons is None else [list(e) for e in exons]),
        ))

    # (a.1) real loci with diverse exon counts (multi-intron structure).
    by_nexon = defaultdict(list)
    for dn in loci:
        _, ex = V.dn_exons(dn, meta, skel)
        if ex is None:
            continue
        by_nexon[min(len(ex), 6)].append(dn)
    real_dn = []
    for nex in sorted(by_nexon):          # spread across 2,3,4,5,6+ exon counts
        real_dn.extend(by_nexon[nex][:12])
    real_dn = real_dn[:DNEXON_REAL]
    for dn in real_dn:
        emit_dn(dn, meta, skel)

    # (a.2) synthetic edge cases -- pure functions of (meta triple, introns string),
    #       so the REAL dn_exons on synthetic maps is fully faithful + self-contained.
    syn_meta = {}
    syn_skel = {}
    synth = [
        ("DN_syn_single_blank", ("chrS", 100, 300), ""),              # blank introns -> 1 exon
        ("DN_syn_missing_skel", ("chrM", 100, 300), None),            # no skeleton key -> 1 exon
        ("DN_syn_meta_miss",    None,               None),           # not in meta -> None
        ("DN_syn_lead_eq_start",("chrL", 100, 300), "100-150"),       # a==start (no lead exon)
        ("DN_syn_trail_eq_end", ("chrT", 100, 300), "200-300"),       # b==end (no trail exon)
        ("DN_syn_whole_intron", ("chrW", 100, 300), "100-300"),       # whole span intron -> []
        ("DN_syn_unsorted",     ("chrU", 100, 300), "250-280;120-150"),
        ("DN_syn_adjacent",     ("chrJ", 100, 300), "120-150;150-200"),
        ("DN_syn_zero_len_intr",("chrZ", 100, 300), "150-150"),       # a==b
        ("DN_syn_three_introns",("chrE", 100, 400), "120-140;200-230;330-360"),
    ]
    for dn, mm, intr in synth:
        if mm is not None:
            syn_meta[dn] = mm
            if intr is not None:
                syn_skel[(mm[0], mm[1], mm[2])] = intr
    for dn, mm, intr in synth:
        emit_dn(dn, syn_meta, syn_skel)

    # ---------------------------------------------------------------- part (c) universe
    # First CATALOG_N sorted loci = the closed universe; fetch case-carrying spliced seq.
    universe = loci[:CATALOG_N]
    loci_out = []
    seqs = {}
    for dn in universe:
        chrom, exons = V.dn_exons(dn, meta, skel)
        s = "".join(fa.fetch(chrom, a, b) for a, b in exons)
        seqs[dn] = s
        loci_out.append(dict(
            dn=dn, gene=gene_of_dn[dn], seq=s,
            chrom=chrom, exons=[list(e) for e in exons],
        ))
    fa.close()

    # (b) node_set on the first NODESET_N loci (>= 30).
    for rec in loci_out[:NODESET_N]:
        node_set = sorted({v for v, _ in V.minimizers(rec["seq"])})
        rec["node_set"] = [int(v) for v in node_set]

    # Build the node catalog over EXACTLY the universe (verbatim from vg_repeat_catalog
    # :245-289) -- multiplicity counts DISTINCT GENES; cyclic = max within-tx run count.
    node_genes = defaultdict(set)
    cyc_max = defaultdict(int)
    lc_occ = defaultdict(int)
    tot_occ = defaultdict(int)
    n_tx = defaultdict(int)
    for dn in universe:
        g = gene_of_dn[dn]
        s = seqs[dn]
        if len(s) < V.K:
            continue
        runs = V.minimizers(s)
        loc_runcount = defaultdict(int)
        seen = set()
        for val, mk in runs:
            node_genes[val].add(g)
            seen.add(val)
            loc_runcount[val] += 1
            tot_occ[val] += 1
            lc_occ[val] += mk
        for val in seen:
            n_tx[val] += 1
            if loc_runcount[val] > cyc_max[val]:
                cyc_max[val] = loc_runcount[val]
    mult = {c: len(gs) for c, gs in node_genes.items()}
    node_sm = {c: (lc_occ[c] / tot_occ[c]) for c in mult if tot_occ[c] > 0}

    # Node section text -- byte-identical to vg_repeat_catalog.py:599-604.
    lines = []
    lines.append("# SECTION nodes (mult>=2)")
    lines.append("node_id\tmultiplicity\tn_genes\tcyclic\tis_repeat_m%d\tn_transcripts\tsoftmask_frac" % V.M_OP)
    for c in sorted((c for c in mult if mult[c] >= 2), key=lambda c: (-mult[c], c)):
        isrep = 1 if (mult[c] >= V.M_OP or cyc_max[c] >= 2) else 0
        lines.append(f"{c}\t{mult[c]}\t{mult[c]}\t{1 if cyc_max[c]>=2 else 0}\t{isrep}\t{n_tx[c]}\t"
                     f"{round(node_sm.get(c,0.0),4)}")
    node_section_text = "\n".join(lines) + "\n"

    n_rows = sum(1 for c in mult if mult[c] >= 2)
    n_cyc = sum(1 for c in mult if mult[c] >= 2 and cyc_max[c] >= 2)
    n_mult5 = sum(1 for c in mult if mult[c] >= V.M_OP)

    fixture = dict(
        K=V.K, W=V.W, M_OP=V.M_OP,
        universe_n=len(universe),
        n_node_rows=n_rows,
        dn_exons=dn_cases,
        loci=loci_out,
        node_section_text=node_section_text,
    )

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(fixture, fh, indent=1, sort_keys=True)

    print(f"[fixture] wrote {os.path.normpath(OUT)}")
    print(f"          dn_exons cases : {len(dn_cases)} ({len(real_dn)} real + {len(synth)} synthetic)")
    print(f"          node_set loci  : {NODESET_N}")
    print(f"          catalog universe: {len(universe)} loci -> {n_rows} node rows "
          f"({n_cyc} cyclic, {n_mult5} mult>={V.M_OP})")


if __name__ == "__main__":
    main()
