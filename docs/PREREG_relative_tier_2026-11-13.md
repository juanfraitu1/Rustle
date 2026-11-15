# Row 1676 — a second tier of "relative" links between loci, families unchanged (2026-11-13)

User choice after row 1675: "Second tier of relatives". Families layer (O1). `fam.copies.*`, the families and assign
are NOT changed by this row; the tier is a separate table, opt-in if it passes. Not for `advisor_share`.

**Idea.** Rows 1674/1675: alignment features of gate-refused locus pairs find homologous gene pairs but cannot date the
duplication. So the pairs are reported as relatives (a link between loci or families), not merged into families.

**Seen before registration (disclosed):** row 1675's out-of-fold ranking, fitted on young-duplicate truth only, where
(true + related) / all picks was about 0.91 at 80 true pairs. Not seen: any fit on the label below, any gorilla score.

**Label (frozen).** A pair is *correct* when it is a simiiformes truth pair or related by pairscore's definition (a
Compara paralog pair at any level, or one HGNC gene group); *wrong* when both genes have a human Ensembl id and neither
holds; *undecidable* (set aside) when a gene has no human id (gorilla genes with no mapped human gene) or, on gorilla,
when both map to the same human gene. `bench/denovo_eval/edgerelated.py` adds the column.

**Model (frozen).** Row 1675's model and features (`edgescore.py`), positives = correct pairs, fit on A119b + HG002
pooled. Threshold = the lowest score whose 5-fold out-of-fold precision is >= 0.90 (`--minp 0.90`).
**Go/no-go:** at that threshold the out-of-fold picks hold >= 100 correct pairs of which >= 30 simiiformes truth pairs;
otherwise STOP.

**Test (frozen weights and threshold).** OR6737 + KB3781 pooled: precision >= 0.85 over decidable picks, >= 50 correct,
>= 20 simiiformes truth pairs among them. Void if more than half of the picks are undecidable.

**Reported, no bar:** recall of young-duplicate truth pairs carried by either tier (family or relative link) against
the family tier alone, per sample.

**If it passes:** an opt-in product proposal to the user (a `PREFIX.fam.relatives.tsv` table written from the loci and
their PAF with the frozen weights). Nothing is wired before the user says so. **If it fails:** the row closes.

## Result: PASS as registered (2026-11-13, `r1674_edgerule/{rel.out,rel.*.tsv,relpick.rel.{OR,KB}.tsv}`)

**Fit (A119b + HG002 pooled, 370,077 pairs, 6,315 correct).** Out-of-fold threshold at precision 0.90: score >= 0.0511,
1,208 correct / 134 wrong, 209 simiiformes truth pairs among them. Go/no-go (>= 100 correct, >= 30 truth): GO.
Per sample, out of fold: A119b 774 / 90, truth pairs 129 of 159 reachable; HG002 434 / 44, truth pairs 80 of 118.

**Gorilla test, frozen weights and threshold:**

| sample | population | picks correct | of which truth pairs | wrong | precision |
|---|---|---|---|---|---|
| OR6737 | 8,152 | 482 | 88 (of 111 reachable) | 36 | 0.931 |
| KB3781 | 226,839 | 382 | 71 (of 99 reachable) | 85 | 0.818 |
| pooled | | 864 | 159 | 121 | 0.877 |

Bars: precision >= 0.85 (0.877), >= 50 correct (864), >= 20 truth pairs (159), undecidable picks 0. **PASS.**
KB3781 alone is 0.818, under 0.85; the bar was registered on the pool.

**Reported, no bar — young-duplicate truth pairs carried by either tier (family or relative link):**
OR6737 768 -> 856 of 1,012 (0.759 -> 0.846); KB3781 108 -> 179 of 240 (0.450 -> 0.746); out of fold, A119b
1,765 -> 1,894 of 2,180 (0.810 -> 0.869), HG002 184 -> 264 of 414 (0.444 -> 0.638). Family precision is unchanged
because the families are unchanged; the relative tier's own precision against young-duplicate truth alone is low
(OR6737 88 of 518 picks) by design: most links are older relatives.

**Comparator (post hoc, descriptive):** the same number of picks ranked by the best single feature (shared exon
bases) is 0.680 / 0.537 correct on OR6737 / KB3781; by record count 0.386 / 0.002.

**Weights (standardized):** same_chrom +0.41, exons_max -0.37, aln_bp -0.35, n_records +0.28, sd_held +0.25,
exon_bp +0.16, identity -0.16, exon_len_ratio +0.12, span_ratio +0.11, rank_max -0.11; the rest within +-0.08.

**Limits.** "Related" = Compara any level or one HGNC group, mapped through human ids on gorilla (protein-hit truth);
pairs of non-coding or unexpressed genes are outside the judged population; the human figures are out-of-fold on the
samples the model was fitted on; one model, one threshold, two species. Nothing is wired; the product proposal
(`PREFIX.fam.relatives.tsv`, opt-in) is the user's call.

## Amendment 1 (2026-11-13, post hoc, reviewer check): the OR6737 table was not comparable; on a comparable one the bar is missed

The OR6737 loci used above (`~/runs/gorilla_OR6737/or.fam.loci`) were built under the old 0.6 coverage gate: 9,081
gene pairs with a record against 205k-227k on the other three samples, 9.4% positives against 1.5%. OR6737's families
were rebuilt with today's driver (`r1674_edgerule/or_all`, job 4470689; default score unchanged, 0.359 / 0.759 as
scored): 259,065 gene pairs with a record, 123 reachable missed truth pairs. Frozen model and threshold, no refit:

| sample | picks correct | truth pairs | wrong | precision |
|---|---|---|---|---|
| OR6737 rebuilt | 514 | 97 (of 123 reachable) | 98 | 0.840 |
| KB3781 | 382 | 71 (of 99) | 85 | 0.818 |
| pooled | 896 | 168 | 183 | **0.830** (bar 0.85) |

**Verdict to use: NOT passed.** The registered PASS (0.877) rests on the pre-filtered OR6737 table; on tables built the
same way as the fit, the pooled precision is 0.830, 0.02 under the bar. The count bars hold (896 correct, 168 truth
pairs). No threshold was moved.

**What "correct" is made of (picks per sample):**

| sample | young-duplicate truth | older Compara paralog | HGNC group only | wrong | ZNF-type pairs | precision without ZNF-type |
|---|---|---|---|---|---|---|
| OR6737 rebuilt | 97 | 239 | 178 | 98 | 119 | 0.801 |
| KB3781 | 71 | 173 | 138 | 85 | 121 | 0.754 |
| A119b (fitted on) | 128 | 271 | 372 | 84 | 389 | 0.820 |
| HG002 (fitted on) | 79 | 163 | 187 | 44 | 323 | 0.707 |

So the tier is: about 1 in 6 picks a young duplicate the families miss, about 2 in 5 an older Compara paralog, about
1 in 3 a pair that only shares an HGNC gene group (mostly KRAB zinc fingers, a domain-defined group), and 1 in 6-7
unrelated. Counting HGNC-only pairs as wrong, precision is 0.47-0.60. "Finds homologs" overstates it; it finds
relatives by sequence or by domain family.

**What stands:** the links carry most of the reachable missed young-duplicate pairs: OR6737 97 of 123, KB3781 71 of 99
(two-tier coverage of truth pairs 0.759 -> 0.855 and 0.450 -> 0.746). The family definition's own edge recall is
unchanged. Human figures are on the samples the model was fitted on and the two share human gene pairs across folds.
The picks files now carry named columns (`relpick.rel.*.tsv`; the first version wrote the wrong columns).

## Wired as the opt-in stage `relatives` (2026-11-14, register 1737; user decision, uncommitted)

No rule was changed. What was built:

- `tools/data/relative_links_model.json`: the model refit from `rel.A.tsv` + `rel.H.tsv` (means, SDs, weights, the log1p
  list, threshold 0.05106, and the training min / max of every feature). The refit reproduces the row's figures exactly
  (OR6737 514 / 98, KB3781 382 / 85).
- `bench/denovo_eval/edgefeat.py`: split into `pair_features` (locus pairs, no annotation) and the gene labelling that the
  row 1674-1676 tables use. `feat.OR2.tsv` regenerates byte-identical.
- `bench/denovo_eval/relative_links.py`: scores every refused locus pair with a record and writes the links. Standard
  library only.
- `tools/rustle_pipeline.sh relatives --bam B --fasta G --out PREFIX` (named only, never part of `all`) ->
  `PREFIX.fam.relatives.tsv`: `locus_a locus_b family_a family_b score stage share identity aln_bp sd_held in_fitted_range`.
  Needs `PREFIX.fam.loci.gff3`, `.loci.gated.paf`, `.copies.tsv` and a SEDEF BED (`RUSTLE_SEDEF_BED` or `FASTA.sedef.bed`;
  refused without one, since `sd_held` is a model feature). A119b: 273 s, 3.7 GB.

Population: a record, loci not overlapping on one chromosome, not both in one family (copy rows map to loci by
`gene_id`; an extra-representative piece is `<locus id>.pN`), stage not `edge`.

**Parity with this row's picks.** OR6737: 604 of 612 picks are links with scores within 0.0006 (the row scored features
rounded to 4 digits); the other 8 are 7 overlapping-locus pairs (refused on purpose; one is the truth pair APOBEC3C /
APOBEC3D) and 1 pair in one family at locus level. A119b: 819 of 855; of the 36 others 19 overlap, 11 are in one
family, 5 loci are absent from this run (the row's table came from the `port1593` loci) and 1 is unexplained
(ZNF729 / ZNF98).

**Satellite loci.** Without a set-aside, A119b gives 37,130 links and 27,843 of them have a locus lying >= 50% inside the
censat BED. `RUSTLE_EXCLUDE_BED` (the BED `assign` uses) sets those loci aside; rows named `ct*` (centromeric transition,
not satellite) are not counted, because the whole BED removed 364 of the 855 measured picks (NBPF, GOLGA, AGAP).
With it: 17,299 links (3,085 between two families, 7,554 family to non-family locus, 6,660 between non-family loci).

**What is measured and what is not.** The 0.830 of this row holds for protein-coding gene pairs, one locus pair per gene
pair: at most 819 of the 17,299 A119b links. All coding gene pairs among the links, scored at edge level with pairscore:

| sample | coding pairs | truth | related | unrelated | truth or related |
|---|---|---|---|---|---|
| OR6737 (held out for the model, no satellite BED) | 681 | 104 | 477 | 100 | 0.853 |
| A119b (fitted on) | 1,122 | 235 | 796 | 91 | 0.919 |

A further 621 OR6737 pairs involve a non-coding gene and are unjudged, as are links between loci with no annotated gene.
KB3781 was not rerun through the stage.

**Score order is not confidence at the top.** Training scores end at 7.2; links between long repeat-like loci reach
23-35 with exon share 0. `in_fitted_range` = 0 marks a link whose score or any feature lies outside the training range
(50 of 1,980 on OR6737, 207 of 17,299 on A119b).
