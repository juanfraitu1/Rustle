# PREREG SHARED-EXON SPLIT (register 1496, 2026-11-02)

Why: alternative-promoter clusters (PCDHA 15, PCDHG 22, UGT1A 9 genes) are one locus in our assembly (A119b families.gtf: 502
shared-exon loci over >= 2 coding RefSeq genes), so their paralogous gene-specific exons can never form a family.
Change: `units_from_annotation.py --cut-shared-exon` cuts every intron between a private exon (one gene) and a shared exon (>= 2
genes); the list goes to `RUSTLE_BRIDGE_REGROUP=f1units RUSTLE_BRIDGE_UNITS_LIST`. Annotation = native CHM13 RefSeq (CAT collapses
these clusters into one gene, so it cannot be the guide), converted GTF -> GFF3 (gene/transcript/exon), `--readthrough both`
(structural precision 0.976 on RefSeq, r1492) so RefSeq readthrough models do not join the parts.
Arms (A119b, `--contig chr2,chr5`, assemble + families): A default (f1v2); C cut-shared.
Metrics: shared-exon cluster loci in PREFIX.families.gtf (merge_audit.py); pairscore coding precision/recall at simiiformes (default)
and, descriptive, at mammalia; families holding >= 2 PCDHA / PCDHG / UGT1A genes.
Decision: C is the recommended annotation-guided option if shared-exon loci fall >= 50% AND simiiformes precision and recall each
fall by <= 0.01. Else FAIL.
