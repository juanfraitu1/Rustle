# Rows 1707-1709 — are the readthrough rules too strict, and how well do we find real joins? (2026-11-14, user goal)

Builds on row 1706 (`bench/annot_compare/fusion_audit.py`, `join_support.py`). Development sample A119b (default run of
2026-11-13, `r1687_sdcut/D0`); held-out sample HG002 Kinnex on CHM13 (`r1692_trunc5_hg/hg.*`, default assemble of
2026-11-13, never looked at for joins); each is the other's recurrence witness (reads on CHM13). Gorilla KB3781 panel is
reported apart, one annotation, no witness, no verdict. Nothing pooled. No default changes.

**Not blind on A119b.** Row 1706 already showed the F1 tag counts per class, the read-support buckets and the secondary
fractions. So N and the bars below were chosen with those in view; only HG002 is a test. A verdict counts when HG002
gives the same side of the bar.

## Unit and classes (frozen)
Annotated gene = any RefSeq gene record with exons, all biotypes (so pseudogene joins count). A **join** is an assembled
transcript with a crossing intron: consecutive exons where the left one lies on a gene and the right one on another gene
that shares no exon base with it (same strand). Transcripts that touch two genes without such an intron are counted in
their own bin and enter no verdict. A transcript's class comes from its weakest crossing intron (fewest primary records).
- **a annotated**: the crossing intron is an intron of a transcript in the RefSeq GTF or in CAT/Liftoff
  (`chm13v2.0_CAT_Liftoff.slim.gff3.gz`), exact coordinates.
- Supported = >= 5 primary records (not secondary, not supplementary) carry exactly that intron in the sample's BAM,
  and its splice sites are canonical (GT-AG, GC-AG, AT-AC on the transcript strand).
- **b1** unannotated, supported, recurring: >= 2 primary records carry the intron in the witness BAM.
- **b2** unannotated, supported, not recurring, witness informative (>= 5 primary records overlap each flank).
- **b0** unannotated, supported, witness not informative.
- **w weak**: unannotated, canonical, 2-4 primary records. Own bin, no verdict.
- **c unsupported**: unannotated and (fewer than 2 primary records or non-canonical sites).

## Row 1707 — rules (A119b, then HG002)
| Rule | What it acts on | Too strict if | Too loose if |
|---|---|---|---|
| F1 f1v2 (tags, removes from the families input) | all tagged transcripts, by their own content | > 10% of tagged transcripts lie inside one annotated gene (every exon on one gene, no crossing intron) | it leaves untagged > 50% of the joins of classes a + b1 (a relation it should mark) |
| final.gtf layers (removes) | joins of x.gtf whose intron chain is absent from final.gtf | it removes > 5% of class a + b1 joins | > 50% of class c joins pass |
| readthrough_filter, emulated as secondary share >= 0.5 on the crossing intron (opt-in; labelled an emulation: the binary judges annotation records) | joins it would flag | > 10% of flagged joins are a or b1 | not judged (flags a minority by design) |
| catalog readthrough / mis-chain gates | not on the default path: 0 log lines in the default assemble; they drop single-exon spans engulfing >= 5 junctions in the legacy catalog | recorded, not tested | |
Sanity (reported, no bar): the crossing introns of RefSeq records PKD1P6-NPIPP1, PDXDC2P-NPIPB14P (expected kept by the
emulation) and PKD1P4-NPIPA8, PKD1P3-NPIPA1 (expected flagged). Also reported: joins added by final.gtf.

## Row 1708 — finding the real ones
- **Annotated readthrough genes.** Structural definition: a RefSeq gene with an intron whose flanking exons lie on two
  other genes that share no exon with each other. Expressed = >= 5 primary records on such an intron. Recall = expressed
  ones whose intron is carried by an assembled transcript, for x.gtf and final.gtf. Bar: recall >= 0.80 = found.
- **Precision of our joins** = the class table itself (share a; a + b1), with b2, b0, w, c listed, never called false.
- **Tandem paralog pairs.** Truth paralog pairs (famtruth, primates level), same chromosome and strand, no shared exon,
  <= 1 Mb apart, both spliced and both expressed (>= 5 primary records on an own annotated intron). Outcomes: two models
  (each gene has a transcript that joins nothing), fused only, one missing, both missing. Bar: two models >= 0.80.
- **Chimeric molecules.** No truth on healthy samples, so this is a report with no precision or recall: supplementary
  records genome-wide by gene pair (primary part, other part), molecules per pair, strand agreement, distance class,
  and whether the pair recurs in the witness. A pair is listed as a candidate when >= 5 molecules in the sample and
  >= 2 in the witness join the same two genes on consistent strands. Whether any deserves a model is the user's call.

## Row 1709 — the class b list
One table per sample of every b1 / b2 / b0 join (genes, intron, primary and secondary records, motif, both annotations,
witness records, F1 tag, in final.gtf), near-tie joins of row 1706 marked. Goes to the advisor share when final.

## Done
Each rule has a verdict on A119b confirmed or refuted on HG002; the three tables of row 1708 exist for both samples;
the b lists are on the advisor share; ledger, register and memory are updated.

## RESULT, A119b (development sample; recorded 2026-11-14 before the HG002 output was opened)

First job hit the 90-minute wall after the join table was written; the tail was rerun with the table reused (identical by cmp). Raw output:

```
/anvme/workspace/v138da13-rustle/r1687_sdcut/D0/x.gtf
transcripts 233025; joins (crossing intron) 3599; touching two genes without one 2007; CAT-only annotated introns 124320
row	a	b1	b2	b0	w	c	all
joins	1109	496	710	758	506	20	3599
F1-tagged	73	63	78	36	34	1	285
absent from final.gtf	39	14	31	19	2	0	105
secondary share >= 0.5	103	51	26	23	39	20	262
primaries mostly MAPQ 0	17	8	7	11	11	1	55
F1-tagged transcripts by content: {'exon outside any gene': 265, 'inside one gene': 662, 'join': 285, 'touches two genes, no crossing intron': 47} inside one gene 0.526 (too strict if > 0.10)
F1 leaves untagged a + b1 joins: 1469 of 1605 = 0.915 (too loose if > 0.50)
final.gtf removes a + b1 joins: 53 of 1605 = 0.033 (too strict if > 0.05); class c joins passing: 20 of 20 = 1.000 (too loose if > 0.50)
secondary-share emulation flags a or b1: 154 of 262 = 0.588 (too strict if > 0.10)
joins in final.gtf whose chain is not in x.gtf: 0
control PKD1P6-NPIPP1 no crossing intron
control PDXDC2P-NPIPB14P no crossing intron
control PKD1P4-NPIPA8 no crossing intron
control PKD1P3-NPIPA1 no crossing intron
annotated readthrough genes (structural) 177; expressed 87; recall x.gtf 66 of 87 = 0.759; final.gtf 66 of 87 = 0.759 (bar 0.80); not expressed but modelled 10
tandem truth pairs (same strand, <= 1 Mb, both spliced and expressed) 304: two models 272, fused only 1, one missing 30, both missing 1, two models and also a join 16; two models 0.895 (bar 0.80)
```

## RESULT, HG002 (held out; opened after the A119b block above was written)

```
/anvme/workspace/v138da13-rustle/r1692_trunc5_hg/hg.gtf
transcripts 155267; joins (crossing intron) 1699; touching two genes without one 1202; CAT-only annotated introns 124320
row	a	b1	b2	b0	w	c	all
joins	671	401	293	89	230	15	1699
F1-tagged	48	41	14	0	5	0	108
absent from final.gtf	16	13	4	0	0	0	33
secondary share >= 0.5	72	51	12	1	14	15	165
primaries mostly MAPQ 0	23	0	0	1	5	0	29
F1-tagged transcripts by content: {'join': 108, 'exon outside any gene': 59, 'inside one gene': 452, 'touches two genes, no crossing intron': 2} inside one gene 0.728 (too strict if > 0.10)
F1 leaves untagged a + b1 joins: 983 of 1072 = 0.917 (too loose if > 0.50)
final.gtf removes a + b1 joins: 29 of 1072 = 0.027 (too strict if > 0.05); class c joins passing: 15 of 15 = 1.000 (too loose if > 0.50)
secondary-share emulation flags a or b1: 123 of 165 = 0.745 (too strict if > 0.10)
joins in final.gtf whose chain is not in x.gtf: 0
control PKD1P6-NPIPP1 no crossing intron
control PDXDC2P-NPIPB14P no crossing intron
control PKD1P4-NPIPA8 no crossing intron
control PKD1P3-NPIPA1 no crossing intron
annotated readthrough genes (structural) 177; expressed 77; recall x.gtf 59 of 77 = 0.766; final.gtf 59 of 77 = 0.766 (bar 0.80); not expressed but modelled 3
tandem truth pairs (same strand, <= 1 Mb, both spliced and expressed) 110: two models 109, fused only 0, one missing 1, both missing 0, two models and also a join 6; two models 0.991 (bar 0.80)
```

## VERDICTS (2026-11-14)

Universe: all RefSeq genes (not only coding ones as in row 1706), so the join count is 3,599 on A119b against 1,294 there.
Script change after freezing: `join_truth.py` reuses an existing OUT.joins.tsv (wall-limit resume; the rule is unchanged, tables identical by cmp).

| Rule / bar | A119b | HG002 | Verdict |
|---|---|---|---|
| F1 f1v2: tagged transcripts inside one annotated gene (too strict if > 0.10) | 662 of 1,259 = 0.526 | 452 of 621 = 0.728 | too strict, confirmed |
| F1 f1v2: class a + b1 joins left untagged (too loose if > 0.50) | 1,469 of 1,605 = 0.915 | 983 of 1,072 = 0.917 | too loose, confirmed |
| final.gtf layers: class a + b1 joins removed (too strict if > 0.05) | 53 of 1,605 = 0.033 | 29 of 1,072 = 0.027 | fine, confirmed |
| final.gtf layers: class c joins passing (too loose if > 0.50) | 20 of 20 | 15 of 15 | too loose, confirmed (n = 20 and 15) |
| readthrough_filter emulation (secondary share >= 0.5): flagged joins that are a or b1 (too strict if > 0.10) | 154 of 262 = 0.588 | 123 of 165 = 0.745 | too strict, confirmed (an emulation, not the binary's filter) |
| catalog readthrough / mis-chain gates | not exercised: off the default path, 0 log lines | same | no verdict |
| assign readthrough certificate (`readthrough_into`) | not exercised: a column with no consumer (row 1706) | same | no verdict |
| Row 1708: recall of expressed annotated readthrough genes (bar 0.80) | 66 of 87 = 0.759 (x.gtf = final.gtf) | 59 of 77 = 0.766 | FAIL, confirmed |
| Row 1708: tandem truth pairs with a model each (bar 0.80) | 272 of 304 = 0.895 | 109 of 110 = 0.991 | PASS, confirmed |

Reading. The two F1 lines together say the `fusion_of` tag does not follow gene joins: most tagged transcripts sit inside
one gene and nine in ten annotated or recurring joins carry no tag. Cost of the tag inside one gene (post hoc): 198 (A119b)
and 83 (HG002) genes have a tagged transcript of their own, and none is left without an untagged one. The final.gtf layers
remove almost no join of any class (105 of 3,599; 33 of 1,699), so they neither harm the real ones nor remove the
unsupported ones. The secondary-share flag catches all class c joins (all have fewer than 2 primary reads) but most of what
it flags is annotated or recurring. The four NPIP controls are void, not passed: their component genes have no exon
records of their own in this RefSeq file, so no crossing intron exists (154 of the 177 structural readthrough genes carry
a hyphenated name, so the definition itself works).

Class table of our joins (the precision table; b is its own outcome):

| sample | a | b1 | b2 | b0 | w | c | all |
|---|---|---|---|---|---|---|---|
| A119b | 1,109 (221 in RefSeq, 888 CAT only) | 496 | 710 | 758 | 506 | 20 | 3,599 |
| HG002 | 671 (189 in RefSeq) | 401 | 293 | 89 | 230 | 15 | 1,699 |

Tandem pairs, post hoc (the registered outcome counts any non-joining transcript as a model): both genes with a model
on an own annotated intron 266 of 304 and 104 of 110; those with no locus in common 258 (0.849) and 101 (0.918). The class is
expressed truth pairs with both genes spliced; it is not reconciled with the 69 fused loci of rows 1516 / 1686. 16 and 6
pairs have two models and also a joining transcript.

Class (b) per crossing intron (row 1709, `r1707_joins/{a119b,hg002}.class_b.tsv`): A119b 298 b1 / 428 b2 / 432 b0;
HG002 222 / 227 / 55. Motifs on A119b: 1,121 GT-AG, 35 GC-AG, 2 AT-AC. On 514 of
1,158 (A119b) and 228 of 504 (HG002) introns the two sides lie in the span of at least one common CAT gene, and
249 / 87 have a side outside every CAT gene: for the user to judge, no claim on which annotation is right.
Readthrough recall by name: hyphen-named genes 54 of 74 and 49 of 67, others 12 of 13 and 10 of 10. The 21 A119b
misses have 5-104 primary reads on the best crossing intron (14 of them under 12); 6 genes are missed in both
samples (DUS4L-BCAP29, HNRNPUL2-BSCL2, SLMO2-ATP5E, ARHGAP19-SLIT1, CHKB-CPT1B, RAB4B-EGLN2).

Chimeric molecules (report only; supplementary records, both parts on two RefSeq genes): gene pairs with >= 5 molecules
1,203 (A119b) and 2,115 (HG002); also >= 2 in the other sample 528 and 762. By span relation (>= 5 molecules /
candidates): overlapping spans 639 / 333 and 838 / 439 (sense-antisense and nested genes, e.g. MALAT1-TALAM1: one
locus, not a fusion); adjacent same strand 137 / 48 and 304 / 97; adjacent opposite strands 171 / 82 and 364 / 134;
>= 1 Mb 103 / 42 and 175 / 56; other chromosome 153 / 23 and 434 / 36. A zero in the other sample is confounded by
expression (XIRP2 381 vs 0). Whether any class deserves a model is the user's decision.

Gorilla KB3781: not classified here (one annotation, no second sample); its row 1706 panel figures stand alone.
