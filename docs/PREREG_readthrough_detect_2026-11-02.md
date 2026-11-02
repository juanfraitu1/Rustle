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

## RESULT r1492 (2026-11-02): ❌ FAIL (bar 1)
| annotation | labelled | structural | both | recall | precision |
|---|---|---|---|---|---|
| GENCODE v46 basic, chr1-22,X | 960 | 771 | 727 | **0.757** | 0.943 |
| RefSeq GRCh38.p14, NC chr1-22,X | 99 | 84 | 82 | 0.828 | 0.976 |
Bar 3 PASS: `--readthrough labels` units list byte-identical to the previous script on RefSeq chr16-22 (72 cuts each).
By the rule: default = `labels`; `structural` / `both` opt-in. Evaluator `bench/readthrough_detect/rt_eval.py`; outputs `/anvme/workspace/v138da13-rustle/annot/rt.*.txt`.
Found while reading the misses (post hoc, descriptive): GENCODE tags `readthrough_transcript` on single isoforms of normal genes; 32 of the
951 tagged genes are mixed (919 have every transcript tagged), and the any-record label would drop those 32 whole genes.

## PREREG r1493 (2026-11-02, before the run): gene-level label
LABEL = the gene record's own attributes contain `readthrough`, OR every transcript under it does (a mixed gene is kept).
Bars: (a) RefSeq chr16-22 `labels` list still byte-identical to the previous script; (b) GENCODE: the 32 mixed genes are kept and the
919 all-tagged genes still dropped. Pass both -> this is the `labels` definition (default). Fail -> `labels` stays gene-record-only
(the previous behaviour) and the transcript tag is ignored.

## RESULT r1493 (2026-11-02): ❌ FAIL as worded (bar b)
(a) PASS: RefSeq chr16-22 list byte-identical. (b) 919/919 all-tagged dropped, but 2 of the 32 mixed genes dropped too: ARHGAP11B and
CARMIL2 carry `tag=readthrough_gene` on the GENE record (the definition drops them by design; the bar's count missed gene tags).
By the rule: `labels` = gene record only (GENCODE: the 654 genes tagged `readthrough_gene`; ~265 transcript-only-tagged genes not caught).
Gorilla RefSeq GFF (no labels at all): structural calls 2 (LOC101124002, LOC115930417), for review; no gorilla readthrough is labelled.
Shipped: default `--readthrough labels` (gene record, any attribute: RefSeq description, GENCODE tag=readthrough_gene),
`--readthrough structural|both` and `--readthrough-list` opt-in, `--readthrough-report` (gene, reason).
