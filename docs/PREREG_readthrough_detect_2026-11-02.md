# PREREG READTHROUGH DETECT (register 1492, 2026-11-02)

Question: can `bench/units_from_annotation.py` find an annotation's readthrough (fused) gene models on any annotation convention?
Today it drops a gene only when its RefSeq `description` contains "readthrough"; CHM13 CAT, gorilla RefSeq GFF/GTF and HG002 carry no
such label (0 / 0 / 7 incidental), so their readthrough models stay in as genes, overlap both parts, and the fused transcript is one
block = never cut (silent).

## Change
- `--readthrough labels|structural|both` (default `both`), `--readthrough-list FILE` (gene IDs or names, always honoured).
- LABEL: any attribute value of the gene or of a transcript under it contains `readthrough` (RefSeq `description`, GENCODE/Ensembl
  `tag=readthrough_transcript`, biotype, Note), case-insensitive.
- STRUCTURAL: gene R is readthrough when one of its spliced transcripts has an intron whose donor exon shares >= 1 base with exons of
  gene A and whose acceptor exon shares >= 1 base with exons of gene B, where A, B != R, same strand, A and B share no exonic base and
  do not overlap in span, each of A and B is not a small RNA (exon union >= 200 bp and biotype not in miRNA/snoRNA/snRNA/scaRNA/
  tRNA/rRNA/misc_RNA/Y_RNA/vault_RNA), and R's exons cover >= 1 exon of A and >= 1 exon of B that are not exons of the other.
  Thresholds are fixed here and not tuned.
- Side report gets a `reason` column (label / structural / both / list).

## Truth and arms (downloaded, public): GENCODE human (GRCh38) GFF3 and NCBI RefSeq GRCh38 GFF3; truth = LABEL genes.
Gene-level, protein-coding/lncRNA genes, chr1-22,X: structural calls vs label set.

## Decision rule (all must hold to adopt `both` as the default)
1. GENCODE: structural recall >= 0.80 and precision >= 0.80 vs GENCODE labels.
2. RefSeq: structural recall >= 0.80 and precision >= 0.80 vs RefSeq labels.
3. Identity: `--readthrough labels` on the RefSeq GFF gives a units list byte-identical to the current script (GTF = RefSeq transcripts).
4. Gorilla RefSeq (no labels): every structural call listed for hand review; reported, not a bar.
If 1 or 2 fails: `labels` becomes the default (it can only add labelled exclusions), structural stays opt-in, recorded FAIL.
