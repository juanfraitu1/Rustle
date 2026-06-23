#!/usr/bin/env python3
"""family_def_completeness.py — are the RNA families COMPLETE? (DNA-backwards recall)

All prior validation tested PRECISION (are detected families real). This tests COMPLETENESS:
of the DNA-defined complete family set (genome-derived cDNA all-vs-all homology, NOT Compara),
what fraction of each family's members does the RNA pipeline recover -- and WHERE does the gap
come from? Funnel per DNA-family member gene:
  DNA family member -> EXPRESSED (has a de-novo locus) -> in an ~R∩~B RNA family.
The gap is categorized: SILENT (not expressed -> out of RNA scope by design), EXPRESSED-not-
recovered (resolvable / ~B-pruned), RECOVERED. Quantifies recall and shows most of the
incompleteness is by-design (RNA only sees expressed + read-confusable copies).
Run: python bench/family_def_completeness.py
"""
import collections
import json
import os
import sys

import networkx as nx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from family_def_genomewide import best_gene, GENES_BED
from family_def_read_filters import dna_homology

META = "/home/juanfra/winloci_scratch/denovo_transcripts.meta.tsv"
FAM_TSV = "/home/juanfra/winloci_scratch/validated_families.tsv"


def load_genes():
    by = collections.defaultdict(list)
    with open(GENES_BED) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 4:
                by[p[0]].append((int(p[1]), int(p[2]), p[3]))
    for c in by:
        by[c].sort()
    return by


def main():
    by = load_genes()
    # 1. DNA-defined complete family set (genome cDNA all-vs-all, id>=0.90 & recip-cov>=0.30)
    Hd, _ = dna_homology()
    DG = nx.Graph()
    for (ga, gb), r in Hd.items():
        if r.get("id", 0) >= 0.90 and max(r["cov_a"], r["cov_b"]) >= 0.30:
            DG.add_edge(ga, gb)
    dna_fams = [c for c in nx.connected_components(DG) if len(c) >= 2]
    dna_genes = set(DG.nodes())

    # 2. EXPRESSED genes = genes overlapped by a de-novo locus
    expressed = set()
    with open(META) as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            g = best_gene(by, p[1], int(p[2]), int(p[3]))
            if g:
                expressed.add(g)

    # 3. genes that landed in an ~R∩~B RNA family (the 196)
    rna_fam_genes = set()
    fams = collections.defaultdict(list)
    for line in open(FAM_TSV).read().splitlines()[1:]:
        fi, lid, c, s, e, nr = line.split("\t")
        fams[int(fi)].append((c, int(s), int(e)))
    for fi, loci in fams.items():
        gs = {best_gene(by, c, s, e) for c, s, e in loci} - {None}
        if len(gs) >= 2:
            rna_fam_genes |= gs

    # 4. completeness funnel over DNA-family member genes
    tot = len(dna_genes)
    expr = len(dna_genes & expressed)
    recov = len(dna_genes & rna_fam_genes)
    silent = len(dna_genes - expressed)
    expr_not_recov = len((dna_genes & expressed) - rna_fam_genes)
    print(f"=== COMPLETENESS (DNA-backwards): RNA recall vs the genome's cDNA paralog families ===")
    print(f"  DNA-defined paralog families (cDNA all-vs-all): {len(dna_fams)} families, {tot} member genes")
    print(f"\n  funnel over the {tot} DNA-family member genes:")
    print(f"    EXPRESSED (has a de-novo locus)     : {expr:5d} ({100*expr/tot:.0f}%)")
    print(f"    -> in an ~R∩~B RNA family (RECOVERED): {recov:5d} ({100*recov/tot:.0f}% of all; "
          f"{100*recov/max(expr,1):.0f}% of expressed)")
    print(f"\n  the gap, categorized:")
    print(f"    SILENT (not expressed -> out of RNA scope by design): {silent:5d} ({100*silent/tot:.0f}%)")
    print(f"    EXPRESSED but not in an RNA family (resolvable/pruned): {expr_not_recov:5d} ({100*expr_not_recov/tot:.0f}%)")
    print(f"    RECOVERED                                            : {recov:5d} ({100*recov/tot:.0f}%)")

    # per-family completeness (how many families recovered at all / fully)
    comp = []
    for fam in dna_fams:
        m = len(fam); e = len(fam & expressed); r = len(fam & rna_fam_genes)
        comp.append((m, e, r))
    any_recov = sum(1 for m, e, r in comp if r >= 2)
    expr2 = sum(1 for m, e, r in comp if e >= 2)        # families with >=2 expressed members (RNA-detectable)
    recov_of_detectable = sum(1 for m, e, r in comp if e >= 2 and r >= 2)
    print(f"\n  per-family: {len(dna_fams)} DNA families; {expr2} have >=2 EXPRESSED members "
          f"(RNA-detectable in principle); of those, {recov_of_detectable} recovered as an RNA family "
          f"({100*recov_of_detectable/max(expr2,1):.0f}%)")
    sizes = collections.Counter(m for m, e, r in comp)
    print(f"  DNA family sizes: " + ", ".join(f"{k}:{v}" for k, v in sorted(sizes.items())[:8]) + " ...")

    json.dump(dict(dna_families=len(dna_fams), dna_member_genes=tot,
                   expressed=expr, recovered=recov, silent=silent, expressed_not_recovered=expr_not_recov,
                   detectable_families=expr2, recovered_detectable=recov_of_detectable),
              open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "family_def_completeness.json"), "w"), indent=2)
    print("\n[+] wrote family_def_completeness.json")


if __name__ == "__main__":
    main()
