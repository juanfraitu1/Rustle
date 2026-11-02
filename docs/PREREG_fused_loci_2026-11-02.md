# PREREG FUSED LOCI (register 1495, 2026-11-02)

Why: A119b families input has 543 loci joining >= 2 coding RefSeq genes that share no exon (merge_audit.py; 684/691 in the assembly
carry a transcript that itself joins them). Today's default cutter is `--bridge-regroup f1v2`.
Arms (A119b, panel chr17 + chr19 = 105 of the 543, `rustle_pipeline.sh assemble --contig chr17,chr19` then `families`):
- A: default (f1v2).
- B: annotation-guided: `bench/units_from_annotation.py chm13v2.0_CAT_Liftoff.slim.gff3.gz A.gtf --readthrough labels
  --readthrough-list ~/data/denovo_eval/chm13_readthrough.ids` -> `RUSTLE_BRIDGE_REGROUP=f1units RUSTLE_BRIDGE_UNITS_LIST`.
Metrics: fused loci in PREFIX.families.gtf (merge_audit.py vs RefSeq coding genes); pairscore human coding precision/recall (panel).
Decision: B is the recommended annotation-guided mode (documented, opt-in: needs an annotation) if fused loci fall >= 25% AND panel
precision and recall each fall by <= 0.01. Else FAIL; the de novo route (read evidence, step 3 of the plan) is next either way.

## RESULT r1495 (2026-11-02): ✅ PASS (panel)
| arm | loci | fused | shared-exon | pred | hit | precision | recall |
|---|---|---|---|---|---|---|---|
| A f1v2 | 4,885 | 106 | 53 | 311 | 77 | 0.248 | 0.438 |
| B guided f1units | 5,086 | 55 (-48%) | 59 | 201 | 79 | 0.393 | 0.446 |
Units list: 594 transcripts cut at 609 introns (128 readthrough genes dropped by labels/list). Precision +0.145, recall +0.008.
B is the recommended annotation-guided mode (opt-in, needs an annotation): `RUSTLE_BRIDGE_REGROUP=f1units RUSTLE_BRIDGE_UNITS_LIST=`
units_from_annotation.py output. Panel only (chr17+chr19; panel truth = panel genes): a genome-wide confirmation is the next run.

## PREREG r1495 genome-wide confirmation (2026-11-02, before the run)
Same arms, A119b genome-wide (A default, B guided with the list built from A's PREFIX.gtf, no --contigs). Same bars: fused loci in
PREFIX.families.gtf fall >= 25% AND coding pair precision and recall each fall by <= 0.01. Pass -> documented as the recommended
annotation-guided mode; fail -> the panel result is recorded as not replicated.

## RESULT r1495 genome-wide (2026-11-02, job 4407908): ✅ PASS
| arm | loci | fused | shared-exon | pred | hit | truth | precision | recall |
|---|---|---|---|---|---|---|---|---|
| A f1v2 | 76,161 | 547 | 502 | 2,228 | 1,617 | 2,183 | 0.7258 | 0.7407 |
| B guided f1units | 77,705 | 282 (-48%) | 556 | 2,102 | 1,599 | 2,187 | 0.7607 | 0.7311 |
Units list: 5,044 transcripts cut at 5,234 introns (599 readthrough genes dropped). Precision +0.035; recall -0.0096 (bar 0.01, narrow).
B assemble+families 2,693 s vs A 4,717 s (A wall includes cold caches). Recommended annotation-guided mode (opt-in, needs an annotation).
