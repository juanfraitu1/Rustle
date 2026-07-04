#!/usr/bin/env python3
"""Golden byte-parity fixture generator for the O2 VG-materialization LOADERS + copy-dedup
Rust port (src/rustle/vg_family/o2_materialize.rs).

Ground truth = the SHIPPED Python functions, imported directly:
  * bench/psv_graph_genomewide.py :: dedup_copies + OVERLAP_MERGE + FAM_TSV
  * bench/o2_vg_visualization.py  :: load_families + load_gene_labels + FLAGSHIP + SUN_TSV

Run deterministically:
    PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/gen_o2_materialize_fixture.py

Emits src/rustle/vg_family/testdata/o2_materialize_fixture.json (self-contained: it
ships the raw input rows so the Rust #[test]s need no external files).

Fixture contents:
  overlap_merge   : the OVERLAP_MERGE constant (asserted == Rust const).
  flagship        : FLAGSHIP as [[id, gene, role], ...] in dict/insertion order.
  dedup_cases     : one case per REAL family (>=100) run through the REAL dedup_copies,
                    PLUS hand-built synthetic cases covering: single-locus, tandem-stays-
                    separate, isoform full-overlap merge, multi-chrom first-seen ORDER
                    fidelity, chain/3-way merge (max-end + nreads accumulation), and the
                    reciprocal-overlap boundary JUST-ABOVE / JUST-BELOW / at-equality of
                    OVERLAP_MERGE (both even and odd min-length, exercising the strict >).
                    Each case ships its input loci rows AND the exact returned regions
                    list (ORDER included).
  load_families   : the FULL validated-families TSV text + the full parsed dict.
  load_gene_labels: the FULL sun_identifiability.tsv text + the parsed (chrom,start)->gene map.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import psv_graph_genomewide as pg     # noqa: E402  dedup_copies / OVERLAP_MERGE / FAM_TSV
import o2_vg_visualization as ov      # noqa: E402  load_families / load_gene_labels / FLAGSHIP / SUN_TSV

OUT = os.path.join(HERE, "..", "src", "rustle", "vg_family", "testdata", "o2_materialize_fixture.json")


def region_list(regs):
    """dedup_copies returns [(c, start, end, nr), ...] -> JSON list-of-lists."""
    return [[c, s, e, nr] for (c, s, e, nr) in regs]


def loci_list(loci):
    """load_families gives [(lid, c, s, e, nr), ...] -> JSON list-of-lists."""
    return [[lid, c, s, e, nr] for (lid, c, s, e, nr) in loci]


def make_case(cid, loci):
    """loci = [(lid, c, s, e, nr), ...]; run the REAL dedup_copies and package."""
    regs = pg.dedup_copies(loci)
    return {"id": cid, "loci": loci_list(loci), "expect": region_list(regs)}


def main():
    fams = ov.load_families()              # REAL FAM_TSV -> {int fi: [(lid,c,s,e,nr),...]}
    dedup_cases = []

    # ---- (a1) every REAL family through the REAL dedup_copies (>=100) ----
    for fi in sorted(fams):
        dedup_cases.append(make_case(f"fam{fi}", fams[fi]))

    # ---- (a2) synthetic classes (lid is unused by dedup_copies; label them) ----
    syn = []

    # single locus -> one region
    syn.append(make_case("syn_single", [("L0", "cA", 100, 500, 3)]))

    # tandem, non-overlapping on same chrom -> STAY separate (2 regions)
    syn.append(make_case("syn_tandem_sep",
                         [("L0", "cA", 100, 200, 2), ("L1", "cA", 5000, 6000, 4)]))

    # overlapping but BELOW threshold on same chrom -> stay separate
    #   cur=(100,200) len100, next=(150,250) len100, ov=50, thresh=0.5*100=50, 50>50 False
    syn.append(make_case("syn_boundary_equal_even",
                         [("L0", "cA", 100, 200, 1), ("L1", "cA", 150, 250, 2)]))

    # JUST-ABOVE threshold (even min-length) -> MERGE (max end + nr accum)
    #   cur=(100,200) len100, next=(149,250) len101, ov=51, thresh=50, 51>50 True
    syn.append(make_case("syn_boundary_above_even",
                         [("L0", "cA", 100, 200, 1), ("L1", "cA", 149, 250, 2)]))

    # ODD min-length boundary, NO merge:
    #   cur=(100,201) len101, next=(151,252) len101, ov=50, thresh=0.5*101=50.5, 50>50.5 False
    syn.append(make_case("syn_boundary_odd_nomerge",
                         [("L0", "cA", 100, 201, 1), ("L1", "cA", 151, 252, 2)]))

    # ODD min-length boundary, MERGE:
    #   cur=(100,201) len101, next=(150,251) len101, ov=51, thresh=50.5, 51>50.5 True
    syn.append(make_case("syn_boundary_odd_merge",
                         [("L0", "cA", 100, 201, 1), ("L1", "cA", 150, 251, 2)]))

    # full reciprocal overlap (isoforms) -> merge to one region
    syn.append(make_case("syn_isoform_merge",
                         [("L0", "cA", 1000, 2000, 5), ("L1", "cA", 1010, 1990, 3)]))

    # merged end must keep the LARGER end (max), inner-contained next
    #   cur=(100,500,1), next=(150,300,2): ov=150, minlen=150, thresh=75 -> merge -> (100,500,3)
    syn.append(make_case("syn_merge_keep_maxend",
                         [("L0", "cA", 100, 500, 1), ("L1", "cA", 150, 300, 2)]))

    # 3-way chain merge -> nreads accumulate 1+2+4=7, end=max
    syn.append(make_case("syn_merge_3way",
                         [("L0", "cA", 100, 300, 1), ("L1", "cA", 120, 320, 2), ("L2", "cA", 140, 340, 4)]))

    # chain then break: (100,300,1)+(150,400,2) merge; (350,500,3) separates
    syn.append(make_case("syn_merge_then_break",
                         [("L0", "cA", 100, 300, 1), ("L1", "cA", 150, 400, 2), ("L2", "cA", 350, 500, 3)]))

    # MULTI-CHROM first-seen ORDER fidelity: chrB seen before chrA; return must be
    #   chrB region(s) FIRST (first-seen), NOT alphabetical. Within chrB two separate regions.
    syn.append(make_case("syn_multichrom_order",
                         [("L0", "cB", 1000, 2000, 3), ("L1", "cA", 1000, 2000, 4), ("L2", "cB", 5000, 6000, 2)]))

    # multi-chrom where each chrom has its own isoform merge, interleaved input order
    syn.append(make_case("syn_multichrom_merge",
                         [("L0", "cZ", 100, 400, 1), ("L1", "cY", 100, 400, 2),
                          ("L2", "cZ", 110, 390, 3), ("L3", "cY", 110, 390, 4)]))

    # input not pre-sorted within a chrom -> the internal items.sort() must reorder
    syn.append(make_case("syn_unsorted_input",
                         [("L0", "cA", 5000, 6000, 2), ("L1", "cA", 100, 200, 1), ("L2", "cA", 150, 250, 3)]))

    dedup_cases.extend(syn)

    # ---- (b) load_families: ship the FULL real TSV + the full parsed dict ----
    with open(pg.FAM_TSV) as f:
        fam_tsv_text = f.read()
    fams_expect = {str(fi): loci_list(loci) for fi, loci in fams.items()}

    # ---- (c) load_gene_labels: ship the FULL real SUN TSV + the parsed map ----
    gene_labels = ov.load_gene_labels()     # {(chrom, int start) -> gene}
    if os.path.exists(ov.SUN_TSV):
        with open(ov.SUN_TSV) as f:
            sun_tsv_text = f.read()
    else:
        sun_tsv_text = ""
    gene_expect = [[c, s, g] for (c, s), g in gene_labels.items()]

    fixture = {
        "overlap_merge": pg.OVERLAP_MERGE,
        "flagship": [[int(k), v[0], v[1]] for k, v in ov.FLAGSHIP.items()],
        "dedup_cases": dedup_cases,
        "load_families": {"tsv": fam_tsv_text, "expect": fams_expect},
        "load_gene_labels": {"tsv": sun_tsv_text, "expect": gene_expect},
    }

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(fixture, f)

    n_real = sum(1 for c in dedup_cases if c["id"].startswith("fam"))
    n_merge = sum(1 for c in dedup_cases if len(c["expect"]) < len(c["loci"]))
    print(f"wrote {OUT}")
    print(f"  dedup_cases: {len(dedup_cases)} ({n_real} real families + {len(syn)} synthetic); "
          f"{n_merge} exhibit >=1 merge")
    print(f"  load_families: {len(fams)} families, {len(fam_tsv_text.splitlines())} tsv lines")
    print(f"  load_gene_labels: {len(gene_labels)} labels, {len(sun_tsv_text.splitlines())} tsv lines")
    print(f"  OVERLAP_MERGE={pg.OVERLAP_MERGE}  FLAGSHIP={dict(ov.FLAGSHIP)}")


if __name__ == "__main__":
    main()
