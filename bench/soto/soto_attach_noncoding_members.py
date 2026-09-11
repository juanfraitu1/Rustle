#!/usr/bin/env python3
"""Attach non-eligible-biotype genes (lncRNA, processed_pseudogene, and similar) to an ALREADY-FORMED
family via a shared exon, without letting them merge two different families together or found a family
on their own -- matching Soto's own text ("SD98 genes associated with other gene features, including
lncRNAs and processed pseudogenes, were also assigned a gene family ID [when they join one]") and the
spirit of their real algorithm's coding-gene-gated bridging (docs/o1_ledger.md §6if Finding 2: expansion
in their actual code only propagates through protein-coding/unprocessed-pseudogene genes, never through
a shared non-coding gene).

WHY THIS EXISTS. §6ie/§6if's family-clustering pipeline restricted its ENTIRE gene universe (both graph-
building input and scoring truth) to Soto's own 1,793-gene "family-eligible" set (In Table S1=Yes), which
excludes ~541 real genes (mostly lncRNA/processed_pseudogene, per docs/o1_ledger.md §6ih) that Soto's own
S1C table DOES place in a family. Scoring only over 1,793 genes when 2,185 (2,334 minus 149 with
ambiguous multi-family ground truth) have a real, checkable answer is the SAME shrinking-denominator
metric trap this project's own Soto replication already caught and retracted once before (2026-08-02,
"HEADLINE NUMBERS RETRACTED... calling a gene solitary IS a failure to reproduce, not an absence of
data"). This script is the fix's second half: given an eligible-only family assignment (the graph-
building/clustering step is unchanged), attach whichever extra genes can be, so the SAME complete,
honest universe can be scored both with and without this step.

NOTE ON WHY EXTRA GENES ARE NOT SIMPLY ADDED TO THE GRAPH-BUILDING STEP INSTEAD: tested (docs/o1_ledger.md
§6ih) and REJECTED -- letting lncRNA/processed_pseudogene genes act as full graph nodes (able to bridge
two otherwise-separate components together, not just join one) collapses precision (0.925->0.73 mean-MAD)
because these biotypes are more repetitive/promiscuous and create spurious merges, the same failure mode
as this project's own earlier mega-component bug (§6ie Bug 1). Attaching AFTER family formation, without
letting an attached gene ever cause two families to merge, avoids this.

Usage: soto_attach_noncoding_members.py --families replicated_families_1793_medianmad.tsv
       --shared shared_exons_2334.tsv --eligible-geneset soto_1793_geneset.tsv
       --full-geneset soto_2334_geneset.tsv --out replicated_families_medianmad_plus_attached.tsv
"""
import argparse, csv, sys
from collections import defaultdict


def attach(gene_family, extra_genes, shared_path):
    """For each extra (non-eligible) gene with >=1 shared-exon edge to a gene already in gene_family,
    attach it to whichever family it has the MOST such edges with (a disclosed tie-break for the rare
    case of edges to more than one family -- Soto's own text does not specify one). Returns
    {gene: family_id} for attached genes only.
    """
    edge_count = defaultdict(lambda: defaultdict(int))
    with open(shared_path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            a, b = r["gene_a"], r["gene_b"]
            if a in extra_genes and b in gene_family:
                edge_count[a][gene_family[b]] += 1
            if b in extra_genes and a in gene_family:
                edge_count[b][gene_family[a]] += 1
    return {g: max(fams.items(), key=lambda kv: kv[1])[0] for g, fams in edge_count.items() if fams}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--families", required=True, help="eligible-only family assignment TSV (gene_id, family_id columns)")
    ap.add_argument("--shared", required=True, help="shared-exon edge TSV over the FULL (eligible+extra) gene universe")
    ap.add_argument("--eligible-geneset", required=True, help="the eligible-only geneset used to build --families")
    ap.add_argument("--full-geneset", required=True, help="the expanded geneset (eligible + extra) --shared was built from")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    gene_family = {}
    with open(a.families) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("family_id"):
                gene_family[r["gene_id"]] = r["family_id"]

    with open(a.eligible_geneset) as fh:
        eligible = {r["gene_id"] for r in csv.DictReader(fh, delimiter="\t")}
    with open(a.full_geneset) as fh:
        full = {r["gene_id"] for r in csv.DictReader(fh, delimiter="\t")}
    extra_genes = full - eligible

    attached = attach(gene_family, extra_genes, a.shared)

    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene_id", "family_id", "status"])
        for g, f in sorted(gene_family.items()):
            w.writerow([g, f, "eligible_backbone"])
        for g, f in sorted(attached.items()):
            w.writerow([g, f, "attached_noncoding_member"])

    print(f"[done] {len(extra_genes)} extra genes considered, {len(attached)} attached "
          f"({len(extra_genes) - len(attached)} had no shared-exon edge to any formed family) -> {a.out}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
