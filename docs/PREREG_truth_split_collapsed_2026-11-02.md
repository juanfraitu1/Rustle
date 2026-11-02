# PREREG TRUTH SPLIT COLLAPSED (register 1494, 2026-11-02)

Why: CAT collapses genes that share exons into one gene (`collapsed_gene_ids`; docs/ANNOT_CAT_VS_REFSEQ_2026-11-02.md); famtruth
gives that gene ONE truth gene carrying every merged ENSG, so pairs between the merged genes can never be scored, and our assembly
makes the same merges (A119b families.gtf: 502 shared-exon / 543 fused loci over >= 2 coding RefSeq genes).

Change (opt-in env `DENOVO_SPLIT_REFSEQ=REFSEQ.gtf[.gz]`, famtruth.load_cat): a protein-coding CAT gene with collapsed ids is replaced
by one truth gene per coding RefSeq gene (NM_/XM_) on its strand that shares exonic bases with it AND whose HGNC symbol maps to an
ENSG among the CAT gene's (source + collapsed) ids, when there are >= 2 such RefSeq genes: id `CATID#SYMBOL`, exons = the RefSeq
gene's exon union, ens = {that ENSG}, biotype protein_coding. Otherwise the CAT gene is unchanged.

Run: pairscore.py --species human on A119b (~/runs/human_A119b/a1.families.gtf + a1.fam.copies.tsv) and OR6737 gorilla (truth via
symbols of the human CAT genes) with and without the env; also the counts of CAT genes split and truth pairs added.

Decision rule: if coding pair recall or precision moves by >= 0.005 absolute on A119b, the split truth becomes the default human
truth (env default on, headline numbers re-reported with both); otherwise it stays opt-in and the merges are recorded as negligible
for the pair score at the default level (simiiformes).

## RESULT r1494 (2026-11-02): below the bar -> stays opt-in
Split: 188 collapsed coding CAT genes -> 444 RefSeq genes. Truth pairs among expressed coding genes 2,183 -> 2,199 (+16), universe 16,394 -> 16,468.
| sample | arm | pred | hit | precision | recall |
|---|---|---|---|---|---|
| A119b | CAT truth | 2,228 | 1,617 | 0.7258 | 0.7407 |
| A119b | split truth | 2,229 | 1,618 | 0.7259 | 0.7358 |
| OR6737 | both (symbol truth) | 49 | 21 | 0.429 | 0.239 |
Recall change 0.00493 < 0.005 (precision +0.0001): the merged clusters are mostly ancient duplications, outside the simiiformes truth
level, so they hardly enter the pair score. Of the 16 added truth pairs our families find 1. Opt-in: `DENOVO_SPLIT_REFSEQ=REFSEQ.gtf`.
