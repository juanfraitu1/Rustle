# PREREG row 1665: single-exon loci fed to `families` additively (2026-11-13)

Families-layer measurement (O1). Prototype outside the driver: `bench/denovo_eval/mono_feed.py`. Not adopted by this row.

## Why
Member funnel (A119b): 291 findable truth members have no model. 36 of them have >= 3 own unspliced primary alignments
and nothing in `PREFIX.families.gtf` (PPP1R2B 1,107, TTC30B 737). `PREFIX.final.gtf`'s single-exon cluster layer models
some, but `families` never reads it. Row 1490 replaced the single-exon models by that layer and lost copies
(3,299 -> 217: read-sized clusters do not become family copies). Here nothing is replaced.

## Rule (frozen before any families run)
Pile = overlapping same-strand cluster models (`mono_model "cluster"`) of `final.gtf`. A pile is added only when no
exon of `PREFIX.families.gtf` overlaps it on the same strand. One model per pile, two arms, both registered now:
- `most`: the cluster with the most reads.
- `span`: the pile's merged span, unless > 2x its longest cluster (row 1518 blobs): then the longest cluster.

## Step 0 (files on disk, done before the bars were run)
Stop if < 10 of the 291 get an added model on >= half their exons, own strand. Result: `most` 15, `span` 16; loci added
2,671 of 231,766 transcripts (A119b), 293 (KB3781). Go.

## Run
`families` only (no `assign`), fresh directories, default settings, baseline rescored through the same commands.
Funnel `--gtf` = the GTF the arm's families stage read.

Validity gate (not a result): every baseline locus still in `fam.loci.gff3`; added loci counted in `fam.loci.gff3` and
in `fam.copies.tsv`. If no added locus reaches `fam.loci.gff3`, the run says nothing about the rule.

## Bars, per arm, vs the baseline of the same job
A119b (PASS needs all):
1. members placed (member_funnel) >= baseline + 8;
2. pair precision >= baseline - 0.005 and recall >= baseline (pairscore, simiiformes, PAIRSCORE_SAME_SYMBOL=ignore);
3. true pairs lost <= 5.
KB3781 transfer (same arm, submitted together, no rule change in between): protein-hit precision >= baseline - 0.005,
recall >= baseline, members placed >= baseline. Gain there is reported, not required.
An arm that fails stays a prototype. A pass is proposed to the user as opt-in; no arm is chosen after seeing A119b
unless only one passes both.

## Result: FAIL as registered (recall bar), job 4468786

Validity gate passed: no baseline locus missing; added loci entered `fam.loci.gff3` (A119b 2,671, KB3781 293) and
109-111 / 5 of them became family copies.

| run | pred pairs | hit | false | precision | recall | truth pairs (expressed) | placed | no_model | true pairs gained / lost |
|---|---|---|---|---|---|---|---|---|---|
| A119b base | 1,960 | 1,768 | 192 | 0.902 | 0.808 | 2,188 | 482 | 291 | |
| A119b most | 1,974 | 1,779 | 195 | 0.901 | 0.804 | 2,213 | 495 | 274 | 11 / 0 |
| A119b span | 1,974 | 1,780 | 194 | 0.902 | 0.804 | 2,213 | 497 | 274 | 12 / 0 |
| KB3781 base | 493 | 108 | 385 | 0.219 | 0.450 | 240 | 89 | 64 | |
| KB3781 most = span | 495 | 108 | 387 | 0.218 | 0.446 | 242 | 89 | 62 | 0 / 0 |

Bars 1 and 3 pass on A119b (placed +13 / +15, 0 true pairs lost), precision holds. Bar 2 fails: recall 0.808 -> 0.804,
and 0.450 -> 0.446 on gorilla. The hits did not fall; the denominator grew, because pairscore counts truth pairs among
genes that have a model, and the added models make 25 more truth pairs (2 on gorilla) countable that are then not
linked. The bar was written without foreseeing that; it is not loosened here. Gorilla: no gain (2 members modelled,
none placed). Not adopted; prototype stays in `bench/denovo_eval/mono_feed.py`.
