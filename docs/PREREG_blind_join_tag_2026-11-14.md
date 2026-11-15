# PREREG row 1712: a blind gene-join tag in place of F1's bundle-bridge tag (steps 1 and 2 only)

User, 2026-11-14: "maybe F1 could be instead gene joins?". Row 1707 found the `fusion_of` tag inside one annotated gene on
0.526 / 0.728 of tagged transcripts and absent from 0.915 / 0.917 of class a + b1 joins. Question: can a rule that reads
only `x.gtf` (no annotation) mark gene joins? This row scores such a rule on existing tables. No F1 code is edited, no
families run, no default changed. Written before any number of the rule was seen. HG002's class table is already known
(row 1707); HG002 is held out only for this rule's outcome.

## Rule (on x.gtf alone, same chromosome and strand)
A multi-exon transcript T is tagged when one of its introns (donor exon L, acceptor exon R, intron [s, e)) has
- a left model a: a spliced transcript other than T that shares an exon base with L and whose span ends before e;
- a right model b: a spliced transcript other than T that shares an exon base with R and whose span starts after s;
- a and b share no exon base;
- neither a's nor b's intron chain is a contiguous sub-chain of T's (truncated copies of T are not models of a side).
Sides cannot be merged by T or another joiner: a model that crosses the intron is on neither side.

Forms scored on A119b (all listed here; one is carried to HG002):
- F1: the rule as written.
- F2: a and b each have >= 5 reads (`reads` attribute).
- F3: a and b each have at least as many reads as T.
- F4: a's last exon or b's first exon lies wholly inside the intron.
Selection: among forms with inside-one-gene share <= 0.10 on A119b, the one leaving fewest class a + b1 joins untagged;
if no form reaches 0.10, the form with the lowest inside-one-gene share. That one form is frozen before HG002 is run.

## Step 0, ceiling (annotation, A119b)
Share of class a + b1 joining transcripts (row 1707 table) whose two flanking genes each have a transcript that joins
nothing. The rule needs a model on each side, so it cannot tag more than about this. A ceiling under 0.50 means bar 2 is
out of reach for any rule of this shape and the limit is missing separate models, not the rule; the forms are still scored.

## Bars (same kinds and denominators as row 1707, scored by the same code path; self-check: the current tag must
reproduce 662 of 1,259 and 1,469 of 1,605 on A119b)
1. Tagged transcripts inside one annotated gene <= 0.10 of all tagged.
2. Class a + b1 joining transcripts left untagged <= 0.50.
3. A119b only: of the crossing junctions of the 69 spliced fused tandem loci (`~/fused/tandem_fused.tsv`, rows 1516 /
   1686) that are present in this x.gtf (intron match within 1 bp), the share carried by a tagged transcript is at least
   that of the current `fusion_of` tag.
PASS = bars 1 and 2 on A119b and on HG002, and bar 3 on A119b. Otherwise FAIL; if the frozen form is not better than the
current tag on both bar 1 and bar 2, the wording is "adds nothing over the current tag".

Reported without verdict: tag rate per class (b1, b2, b0, w, c; b is never called false), number of tagged transcripts,
tags by content kind, genes left with no untagged model, tags per form.
Tool: `bench/annot_compare/blind_join.py`.

## RESULT, A119b (recorded before HG002 was run)

Self-check holds: the current tag reproduces 662 of 1,259 and 1,469 of 1,605. No form reaches 0.10 inside one gene; by the selection rule F3 (lowest share, 0.450) is frozen for HG002.

```
/anvme/workspace/v138da13-rustle/r1687_sdcut/D0/x.gtf
step 0 ceiling: a + b1 joining transcripts with a non-joining model on both flanking genes 933 of 1605 = 0.581
fused tandem junctions listed 60, present in this x.gtf 60
tag	tagged	inside_one_gene	share(bar<=0.10)	ab_untagged	share(bar<=0.50)	fused_junctions_tagged	kinds	tag_rate_by_class	genes_without_untagged_model
current fusion_of	1259	662	0.526	1469/1605	0.915	0/60	{'join': 285, 'exon outside any gene': 265, 'inside one gene': 662, 'touches two genes, no crossing intron': 47}	{'a': '73/1109', 'b1': '63/496', 'b2': '78/710', 'b0': '36/758', 'w': '34/506', 'c': '1/20'}	0
F1	41688	22713	0.545	745/1605	0.464	35/60	{'join': 1802, 'inside one gene': 22713, 'exon outside any gene': 16730, 'touches two genes, no crossing intron': 443}	{'a': '572/1109', 'b1': '288/496', 'b2': '392/710', 'b0': '329/758', 'w': '210/506', 'c': '11/20'}	20
F2	18138	8209	0.453	1045/1605	0.651	24/60	{'inside one gene': 8209, 'exon outside any gene': 8603, 'join': 1122, 'touches two genes, no crossing intron': 204}	{'a': '382/1109', 'b1': '178/496', 'b2': '289/710', 'b0': '143/758', 'w': '122/506', 'c': '8/20'}	3
F3	24740	11135	0.450	953/1605	0.594	30/60	{'join': 1444, 'exon outside any gene': 11813, 'inside one gene': 11135, 'touches two genes, no crossing intron': 348}	{'a': '430/1109', 'b1': '222/496', 'b2': '334/710', 'b0': '247/758', 'w': '200/506', 'c': '11/20'}	13
F4	29703	15250	0.513	802/1605	0.500	34/60	{'join': 1664, 'inside one gene': 15250, 'exon outside any gene': 12422, 'touches two genes, no crossing intron': 367}	{'a': '529/1109', 'b1': '274/496', 'b2': '358/710', 'b0': '290/758', 'w': '202/506', 'c': '11/20'}	13
```

## RESULT, HG002 (frozen form F3)

```
/anvme/workspace/v138da13-rustle/r1692_trunc5_hg/hg.gtf
step 0 ceiling: a + b1 joining transcripts with a non-joining model on both flanking genes 648 of 1072 = 0.604
tag	tagged	inside_one_gene	share(bar<=0.10)	ab_untagged	share(bar<=0.50)	fused_junctions_tagged	kinds	tag_rate_by_class	genes_without_untagged_model
current fusion_of	621	452	0.728	983/1072	0.917	NA	{'exon outside any gene': 59, 'inside one gene': 452, 'join': 108, 'touches two genes, no crossing intron': 2}	{'a': '48/671', 'b1': '41/401', 'b2': '14/293', 'b0': '0/89', 'w': '5/230', 'c': '0/15'}	0
F1	22527	16178	0.718	514/1072	0.479	NA	{'exon outside any gene': 5317, 'inside one gene': 16178, 'touches two genes, no crossing intron': 215, 'join': 817}	{'a': '349/671', 'b1': '209/401', 'b2': '129/293', 'b0': '45/89', 'w': '80/230', 'c': '5/15'}	16
F2	7318	5101	0.697	741/1072	0.691	NA	{'inside one gene': 5101, 'exon outside any gene': 1629, 'join': 469, 'touches two genes, no crossing intron': 119}	{'a': '179/671', 'b1': '152/401', 'b2': '84/293', 'b0': '16/89', 'w': '34/230', 'c': '4/15'}	6
F3	10703	6788	0.634	663/1072	0.618	NA	{'inside one gene': 6788, 'exon outside any gene': 3169, 'join': 598, 'touches two genes, no crossing intron': 148}	{'a': '243/671', 'b1': '166/401', 'b2': '91/293', 'b0': '25/89', 'w': '68/230', 'c': '5/15'}	9
F4	16734	11589	0.693	561/1072	0.523	NA	{'inside one gene': 11589, 'exon outside any gene': 4193, 'join': 761, 'touches two genes, no crossing intron': 191}	{'a': '330/671', 'b1': '181/401', 'b2': '121/293', 'b0': '45/89', 'w': '79/230', 'c': '5/15'}	15
```

## VERDICT: FAIL (bars 1 and 2, both samples)
| | A119b current tag | A119b F3 | HG002 current tag | HG002 F3 |
|---|---|---|---|---|
| tagged transcripts | 1,259 | 24,740 | 621 | 10,703 |
| bar 1: inside one annotated gene (<= 0.10) | 0.526 | 0.450 | 0.728 | 0.634 |
| bar 2: a + b1 joins untagged (<= 0.50) | 0.915 | 0.594 | 0.917 | 0.618 |
| bar 3: fused tandem junctions tagged | 0 of 60 | 30 of 60 | n/a | n/a |

The table of the 69 loci holds 60 distinct crossing junctions, all present in this x.gtf. Step 0 ceiling: 0.581 / 0.604 of
a + b1 joins have a separate model on both sides, so bar 2 was reachable. F3 is better than the current tag on every bar
but far from bar 1: half or more of its tags are inside one annotated gene, and it tags 20 times more transcripts (about
a tenth of the assembly; 13 genes on A119b would lose every untagged model). The plain form F1 reaches bar 2 on both
samples (0.464 / 0.479) at 0.545 / 0.718 inside one gene. Reading: from the models alone, a gene with alternative first and
last exons looks the same as two neighbouring genes joined by a transcript; "a model ends before, another starts after"
does not separate them, with or without read floors. Nothing was changed.
