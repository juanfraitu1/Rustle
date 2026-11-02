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

## RESULT r1496 (2026-11-02): ❌ FAIL
| arm | loci | shared-exon | fused | pred | hit | precision | recall |
|---|---|---|---|---|---|---|---|
| A f1v2 | 9,550 | 60 | 55 | 43 | 30 | 0.698 | 0.811 |
| C cut-shared | 9,986 | 101 (+68%) | 14 | 90 | 31 | 0.344 | 0.838 |
Shared-exon loci rose instead of falling (the cut units still overlap >= 2 genes through the constant exons), and precision halved
(47 more predicted pairs, 1 more hit). It did cut fused loci (55 -> 14), via the RefSeq guide + structural readthroughs. Stays opt-in;
not recommended. Next idea, not run: drop the shared constant-exon unit from the copy table, keep only the private units.
