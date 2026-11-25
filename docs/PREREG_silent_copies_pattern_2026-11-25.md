# PREREG: do secondary alignments find silent copies, and can those be told from the rest? (row 1955, 2026-11-25)

The advisor (relayed by the user, 2026-11-25): multimappers are very useful for finding silent copies. The user: he has
a point, but they look indistinguishable from error. Check whether the alignments on silent copies are in the near-tied
pool, and if not, whether they have a pattern that could be modelled.

Known (row 1909, same planted set): with every secondary alignment a gene's reads reach 292 of 417 silent processed
pseudogenes, 141 of 150 silent tandem copies and 96 of 128 silent dispersed copies; within 0.98 of the best AS only 7, 49
and 10; and 1,061 places off any planted locus are reached too. Not measured then: the `ms` near-tie band, and whether the
places on silent relatives differ from the other places.

## Data
The planted set of row 1909 (CHM13; reads planted at 150 genes with 417 silent processed pseudogenes, 150 tandem pairs
and 128 dispersed pairs with one copy silent; `minimap2 -N 50 -p 0.1`; workspace `r1909_silent_members/G.p0.1.bam`).
For the near-identical regime, the six half-silent sets of row 1952's seeds 1967-1969 are not needed: row 1954 covers
them (a silent twin receives half of an expressed copy's primaries).

## Method (`bench/o2_truth/silent_reach.py`)
A PLACE = the secondary alignments of reads from one source locus whose aligned blocks overlap, off the source locus,
with at least 2 reads. Its class: `silent relative` (a planted silent member of the source's unit), `other planted locus`,
or `elsewhere`; `elsewhere` is split by the annotation (protein-coding gene, pseudogene, other gene, nothing annotated).
Per place: reads, best and median `ms` relative to the read's best, the share of the read covered, identity of the
aligned part, the share of spliced alignments, aligned bases.
1. **Near-tie pool:** of the silent relatives reached, how many are reached by a read within 10 `ms` of its best.
2. **Pattern:** for each feature, how well it separates `silent relative` from `elsewhere` (area under the ROC curve);
   a logistic model on all features, 5-fold cross-validated by source unit: AUC, and precision at recall 0.5 and 0.8.

## Claims
- **S1, not in the near-tied pool:** fewer than 10% of the reached silent pseudogenes and fewer than 20% of the reached
  silent dispersed copies have a read within 10 `ms`; the tandem copies are reported.
- **S2, a pattern exists:** the cross-validated model separates silent relatives from the places elsewhere with AUC >=
  0.85.
- **S3, it is usable:** at recall 0.5 the model's precision is >= 0.80.

## Predictions (mine)
S1 holds. S2 holds, carried by the share of the read that the alignment covers (a relative aligns over most of the
read; a chance hit covers a piece). S3: uncertain, about 0.7-0.9. The user expects no usable separation.
Caveat fixed beforehand: `elsewhere` is not "error". It contains real homologues the planted units do not list, so the
precision measured here is a LOWER bound; the split by annotation is reported to size that.

## Result (job 4505180; both scenarios of the planted set were run)

### Scenario G: genes and one copy of each pair expressed; pseudogenes and the other copy silent
Places with at least 2 reads: 4,547. Silent relative 511; other planted locus 2,148; elsewhere 1,888 (annotated
pseudogene 908, protein-coding gene 690, other gene 58, nothing annotated 232).

| silent relatives | reached | by a read within 10 `ms` of its best |
|---|---|---|
| processed pseudogenes | 291 | 0 |
| dispersed copies | 93 | 11 (12%) |
| tandem copies | 127 | 48 (38%) |

| feature of a place | median, silent relatives | median, elsewhere | AUC |
|---|---|---|---|
| best `ms` relative to the read's best | 0.60 | 0.86 | 0.38 |
| share of the read covered | 0.92 | 1.00 | 0.41 |
| identity of the aligned part | 0.92 | 0.97 | 0.36 |
| share of spliced alignments | 0.83 | 0.83 | 0.51 |
| reads / aligned bases | 60 / 1,108 | 60 / 1,031 | 0.53 / 0.56 |
Logistic model, cross-validated: AUC 0.68; precision 0.36 at recall 0.5 and 0.29 at recall 0.8. AUC for the pseudogenes
alone 0.78, for the silent copies alone 0.55.

### Scenario P: pseudogenes expressed, their gene silent
Places: 2,470. Silent gene 289 places (122 genes reached, 2 within 10 `ms`); elsewhere 1,359 (annotated pseudogene 1,146).
Share of spliced alignments: 1.00 at the silent genes against 0.00 elsewhere (AUC 0.90). Model AUC 0.91; precision 0.64 at
recall 0.5, 0.60 at recall 0.8.

### Claims
- **S1 HOLDS.** Silent processed pseudogenes are never reached through a near tie (0 of 291), dispersed copies in 12%;
  only tandem copies often are (38%). So the near-tied pool is not where silent relatives are found.
- **S2 FAILS for the question asked** (gene expressed, find its silent relatives): AUC 0.68, bar 0.85. It holds only in
  the reverse direction (0.91). My prediction was wrong: the share of the read covered does not separate them.
- **S3 FAILS in both** (precision 0.36 and 0.64 at recall 0.5; bar 0.80).

### Reading
- **The user's doubt is right about telling them apart, wrong about "error".** The places elsewhere are not noise:
  88% sit on annotated genes or pseudogenes, and they align BETTER than the listed silent relatives (median 0.86 of the
  best score against 0.60; identity 0.97 against 0.92). They are other homologous sequence. A secondary alignment marks
  homology; nothing in it says "family member of interest" as against any other homologue, and nothing says silent or
  expressed (every place here is silent by construction).
- **The advisor's point holds as reach, not as detection.** The reads of an expressed gene do reach most of its silent
  relatives (row 1909: 70% of pseudogenes, 94% of tandem copies). But for each place on a listed silent relative there are
  about eight other places, which score higher. A list made from secondaries is a list of sequence relatives, which a
  comparison of the genome with itself gives without any read.
- **One pattern can be modelled:** a read whose primary is intronless and whose secondary is SPLICED points at the parent
  gene of an expressed processed copy (spliced share 1.00 against 0.00, AUC 0.90). It is the mirror of the ribosomal-
  protein class of row 1953, and it finds a parent gene, not a silent copy in general.
- Limits: one planted set; the truth lists a unit's members only, so precision is a lower bound; real reads not run
  (row 1915 on A119b: places from secondaries gave gene pairs at 0.23 precision and 0.74 recall, the same picture).

### Addendum (description, after the result; user: "how to distinguish the useful ones from all others")
The planted lists are narrower than homology, so every place was judged again against the annotation of related genes
(workspace `annot_related/chm13.related.relations.tsv`: Compara paralogues, pseudogene parents, copies of one source gene).
**Scenario G (gene expressed):** 3,449 of the 4,547 places (76%) are on an annotated relative of the read's gene (recent
paralogue 2,433, pseudogene 486, copy of the same source 395, older paralogue 131); 879 are on a gene that is not a listed
relative and 219 (5%) on nothing annotated. By the place's best score relative to the read's best:

| score relative to the best | places | on an annotated relative |
|---|---|---|
| under 0.3 | 536 | 47% |
| 0.3-0.5 | 314 | 59% |
| 0.5-0.7 | 395 | 75% |
| 0.7-0.9 | 541 | 73% |
| 0.9 and above | 2,761 | 84% |

With this truth the model has AUC 0.70 and precision 0.84 at recall 0.5, 0.83 at recall 0.8.
**Scenario P (pseudogene expressed):** only 330 of 2,470 places are annotated relatives; most of the rest are other
pseudogenes, and the relation table does not list pseudogenes of one parent as relatives of each other, so this is an
undercount. The spliced-secondary signature keeps its AUC (0.84) and the model 0.88.
**So:** most secondary places ARE relatives; the earlier AUC of 0.68 measured "on my short list or not", not "useful or
not". What tells useful from not: (1) the score relative to the read's best is the only feature of the alignment that
helps, and it moves the share of relatives from about half to 84%; (2) no feature tells a family member from another
relative, nor silent from expressed; (3) the spliced-secondary signature names the parent of an expressed processed copy.

### Addendum 2 (description; user: "if we expanded the definition of gene family, could we absorb some of the useful information?")
An extended family = the gene plus every place its reads reach at or above a score threshold (best `ms` of the place
relative to the read's best). Planted set, scenario G:

| threshold | places added | on an annotated relative | silent pseudogenes brought in (of 417) | silent tandem copies (of 150) | silent dispersed copies (of 128) |
|---|---|---|---|---|---|
| 0.98 | 1,899 | 89% | 0 | 54 | 14 |
| 0.9 | 2,761 | 84% | 5 | 77 | 38 |
| 0.7 | 3,302 | 82% | 36 | 88 | 72 |
| 0.5 | 3,697 | 81% | 120 | 110 | 81 |
| 0.3 | 4,011 | 80% | 173 | 121 | 90 |
| none | 4,547 | 76% | 291 | 127 | 93 |

Real reads (A119b, the member places of row 1915, stage `members` at its thresholds: 5 reads, 0.3 of the best AS; the
source locus mapped to its annotated gene; 59,643 rows with an annotated source, 38,335 without):

| score relative to the best | rows | on an annotated relative | on another annotated gene | on nothing annotated |
|---|---|---|---|---|
| 0.3-0.5 | 15,289 | 6% | 66% | 28% |
| 0.5-0.7 | 16,969 | 10% | 72% | 19% |
| 0.7-0.9 | 15,867 | 16% | 74% | 11% |
| 0.9 and above | 11,518 | 35% | 62% | 3% |

On planted reads a wider family absorbs most silent copies at little cost; on real reads the same widening is mostly
places on genes that are not listed relatives. Not established: how many of those are relatives the table lacks (shared
domains, pseudogenes of one parent) and how many are other things (repeats in untranslated regions, overlapping genes).
The planted sources are 150 chosen family genes; the real sources are every expressed locus. No default changed.
