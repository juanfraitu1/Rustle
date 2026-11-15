# PREREG row 1713: a read-end gene-join tag, then its effect on families

User goal 2026-11-14: "develop a rule using read ends (poly-A sites and transcript starts inside the crossing intron)
rather than models ... and the effect on families". Row 1712's model rule failed (0.450 / 0.634 of tags inside one gene).
Found while orienting (coverage count only, no proof column read): F1 already tests read ends, but only at introns that
are the sole link between two pieces of a gene_id; 230 of the 1,605 class a + b1 joins of A119b have their crossing intron
in F1's candidate table. This row applies a read-end test to every intron. Written before any number of the rule was
seen. HG002 has been scored against these bars twice (rows 1707, 1712); it is held out only for this rule's outcome.
No F1 code edit, no default change.

## Evidence (one pass over the BAM, primary records only, no MAPQ filter)
Per record: chromosome, transcript strand (read strand, flipped on ts:A:-), start, end. All positions below are in
transcript orientation. For an intron of a multi-exon transcript T (donor d, acceptor a):
- **end reads**: 3' end inside the intron, at least 10 bp from d and from a, and 5' start upstream of d;
- **start reads**: 5' start inside the intron, same margins, and 3' end downstream of a;
- an **end site** = a 50-bp window with >= N end-read 3' ends; a **start site** = a 50-bp window with >= N start-read 5' starts;
- **signature**: the first end site lies upstream of the last start site (a gene ends, then a gene starts, inside the intron).
T is tagged when one of its introns shows the signature.

Forms scored on A119b (all listed; one frozen for HG002):
- E1: N = 5.
- E2: N = 5, end sites dropped when >= 12 of the 20 genomic bases downstream of the window's first end are A (internal priming).
- E3: E2, and both sites hold at least as many reads as T's `reads`.
- E4: E2 with N = 10.
- M (comparator, models only): the same signature with the assembled transcripts of x.gtf as the "reads" and N = 1.
Selection among E1-E4: with inside-one-gene share <= 0.10 on A119b, the form leaving fewest a + b1 joins untagged; if none
reaches 0.10, the lowest inside-one-gene share.
Ceiling line: a + b1 joins showing the signature at N = 1 with reads.

## Bars (row 1707's kinds, denominators and code path; self-check: the current tag reproduces 662 of 1,259 and 1,469 of 1,605)
1. Tagged transcripts inside one annotated gene <= 0.10 of all tagged. Reported beside it: inside / (inside + join),
   because tags outside any gene dilute bar 1.
2. Class a + b1 joining transcripts left untagged <= 0.50.
3. A119b: fused tandem junctions (60 of `~/fused/tandem_fused.tsv`) carried by a tagged transcript >= the current tag's 0.
PASS = bars 1 and 2 on both samples and bar 3. STOP wording if the frozen form is not better than M on bars 1 and 2:
"adds nothing over the model rule of row 1712". Descriptive: tag rate per class, tagged count, genes left with no
untagged model, and F1's own proofs (up_proof and down_proof) on the joins its table holds.

## Step 3, families (run for the frozen form whatever the bars say; labelled descriptive if they fail)
Through the shipped path: `assemble` + `families` with RUSTLE_BRIDGE_REGROUP=f1units and RUSTLE_BRIDGE_UNITS_LIST naming
the cuts (transcript, crossing junctions), as arms E / S of row 1687 did. Arms on A119b: D0 (default f1v2, 0.902 / 0.808
pair precision / recall, pairscore default level, row 1687), E (F1's bridges as units, 0.910 / 0.809, row 1687), and
R = the frozen form's tags as the units list. Bars: R harms if precision or recall is more than 0.01 under D0; R helps
if one is >= 0.01 above D0 and the other not under it; else no effect. Reported: families, fused tandem loci
(classcount as row 1687). HG002 families effect is run only if R helps on A119b.
Tool: `bench/annot_compare/read_end_join.py`.

## RESULT, A119b, registered forms (recorded before HG002 and before the amendment below was run)

Self-check holds (662 of 1,259; 1,469 of 1,605). No form reaches 0.10; by the selection rule E3 has the lowest inside-one-gene share (0.553). E3 is worse than the models-only comparator M on bar 1 (0.553 vs 0.515) and bar 2 (0.552 vs 0.472).

```
/anvme/workspace/v138da13-rustle/r1687_sdcut/D0/x.gtf
step 0 ceiling: a + b1 joining transcripts with a non-joining model on both flanking genes 933 of 1605 = 0.581
fused tandem junctions listed 60, present in this x.gtf 60
tag	tagged	inside_one_gene	share(bar<=0.10)	inside/(inside+join)	ab_untagged	share(bar<=0.50)	fused_junctions_tagged	kinds	tag_rate_by_class	genes_without_untagged_model
current fusion_of	1259	662	0.526	0.699	1469/1605	0.915	0/60	{'exon outside any gene': 265, 'touches two genes, no crossing intron': 47, 'inside one gene': 662, 'join': 285}	{'a': '73/1109', 'b1': '63/496', 'b2': '78/710', 'b0': '36/758', 'w': '34/506', 'c': '1/20'}	0
ceiling N=1	178199	109648	0.615	0.970	108/1605	0.067	59/60	{'inside one gene': 109648, 'exon outside any gene': 64001, 'join': 3385, 'touches two genes, no crossing intron': 1165}	{'a': '1011/1109', 'b1': '486/496', 'b2': '696/710', 'b0': '711/758', 'w': '463/506', 'c': '18/20'}	10962
M models	36074	18593	0.515	0.910	758/1605	0.472	37/60	{'inside one gene': 18593, 'exon outside any gene': 15215, 'touches two genes, no crossing intron': 437, 'join': 1829}	{'a': '565/1109', 'b1': '282/496', 'b2': '412/710', 'b0': '340/758', 'w': '217/506', 'c': '13/20'}	20
E1	78605	50509	0.643	0.963	680/1605	0.424	33/60	{'inside one gene': 50509, 'exon outside any gene': 25661, 'join': 1967, 'touches two genes, no crossing intron': 468}	{'a': '595/1109', 'b1': '330/496', 'b2': '462/710', 'b0': '345/758', 'w': '229/506', 'c': '6/20'}	1013
E2	56577	35460	0.627	0.953	796/1605	0.496	32/60	{'inside one gene': 35460, 'exon outside any gene': 18978, 'join': 1752, 'touches two genes, no crossing intron': 387}	{'a': '519/1109', 'b1': '290/496', 'b2': '439/710', 'b0': '301/758', 'w': '197/506', 'c': '6/20'}	500
E3	38906	21532	0.553	0.930	886/1605	0.552	32/60	{'inside one gene': 21532, 'exon outside any gene': 15415, 'join': 1609, 'touches two genes, no crossing intron': 350}	{'a': '458/1109', 'b1': '261/496', 'b2': '408/710', 'b0': '279/758', 'w': '197/506', 'c': '6/20'}	37
E4	21920	12737	0.581	0.921	1088/1605	0.678	24/60	{'inside one gene': 12737, 'exon outside any gene': 7887, 'join': 1097, 'touches two genes, no crossing intron': 199}	{'a': '320/1109', 'b1': '197/496', 'b2': '294/710', 'b0': '151/758', 'w': '131/506', 'c': '4/20'}	67
```

## Amendment (post hoc on A119b, written before it was run)
E3 compares the two sites with T's own reads, but T can be a minor isoform sharing the intron with the gene's main
isoforms. Form **E5**: E2, and both sites hold at least as many reads as the summed `reads` of all x.gtf transcripts that
carry the intron (the join is weaker than the gene that ends and the gene that starts; inside one gene the intron's own
reads outnumber intronic ends and starts). Selection is extended to E1-E5 with the same rule; the frozen form is then
scored once on HG002. E5 was designed after seeing A119b's E1-E4, so A119b is not a test of it.

## RESULT, A119b with the amendment (recorded before HG002 was run). Frozen form: E5 (lowest inside-one-gene share, 0.380).

```
tag	tagged	inside_one_gene	share(bar<=0.10)	inside/(inside+join)	ab_untagged	share(bar<=0.50)	fused_junctions_tagged	kinds	tag_rate_by_class	genes_without_untagged_model
current fusion_of	1259	662	0.526	0.699	1469/1605	0.915	0/60	{'exon outside any gene': 265, 'inside one gene': 662, 'join': 285, 'touches two genes, no crossing intron': 47}	{'a': '73/1109', 'b1': '63/496', 'b2': '78/710', 'b0': '36/758', 'w': '34/506', 'c': '1/20'}	0
ceiling N=1	178199	109648	0.615	0.970	108/1605	0.067	59/60	{'exon outside any gene': 64001, 'inside one gene': 109648, 'join': 3385, 'touches two genes, no crossing intron': 1165}	{'a': '1011/1109', 'b1': '486/496', 'b2': '696/710', 'b0': '711/758', 'w': '463/506', 'c': '18/20'}	10962
M models	36074	18593	0.515	0.910	758/1605	0.472	37/60	{'inside one gene': 18593, 'exon outside any gene': 15215, 'join': 1829, 'touches two genes, no crossing intron': 437}	{'a': '565/1109', 'b1': '282/496', 'b2': '412/710', 'b0': '340/758', 'w': '217/506', 'c': '13/20'}	20
E1	78605	50509	0.643	0.963	680/1605	0.424	33/60	{'inside one gene': 50509, 'exon outside any gene': 25661, 'join': 1967, 'touches two genes, no crossing intron': 468}	{'a': '595/1109', 'b1': '330/496', 'b2': '462/710', 'b0': '345/758', 'w': '229/506', 'c': '6/20'}	1013
E2	56577	35460	0.627	0.953	796/1605	0.496	32/60	{'exon outside any gene': 18978, 'inside one gene': 35460, 'join': 1752, 'touches two genes, no crossing intron': 387}	{'a': '519/1109', 'b1': '290/496', 'b2': '439/710', 'b0': '301/758', 'w': '197/506', 'c': '6/20'}	500
E3	38906	21532	0.553	0.930	886/1605	0.552	32/60	{'inside one gene': 21532, 'exon outside any gene': 15415, 'join': 1609, 'touches two genes, no crossing intron': 350}	{'a': '458/1109', 'b1': '261/496', 'b2': '408/710', 'b0': '279/758', 'w': '197/506', 'c': '6/20'}	37
E4	21920	12737	0.581	0.921	1088/1605	0.678	24/60	{'exon outside any gene': 7887, 'inside one gene': 12737, 'join': 1097, 'touches two genes, no crossing intron': 199}	{'a': '320/1109', 'b1': '197/496', 'b2': '294/710', 'b0': '151/758', 'w': '131/506', 'c': '4/20'}	67
E5	13202	5015	0.380	0.821	1128/1605	0.703	30/60	{'exon outside any gene': 6933, 'inside one gene': 5015, 'join': 1096, 'touches two genes, no crossing intron': 158}	{'a': '287/1109', 'b1': '190/496', 'b2': '296/710', 'b0': '136/758', 'w': '181/506', 'c': '6/20'}	29
```

## RESULT, HG002 (frozen form E5; the other forms are printed for reference)

```
tag	tagged	inside_one_gene	share(bar<=0.10)	inside/(inside+join)	ab_untagged	share(bar<=0.50)
current fusion_of	621	452	0.728	0.807	983/1072	0.917
ceiling N=1	131260	95556	0.728	0.983	46/1072	0.043
M models	20222	14150	0.700	0.945	518/1072	0.483
E1	58100	44630	0.768	0.981	510/1072	0.476
E2	43443	33508	0.771	0.978	561/1072	0.523
E3	24956	17212	0.690	0.965	662/1072	0.618
E4	16199	12573	0.776	0.965	785/1072	0.732
E5	9777	5530	0.566	0.928	810/1072	0.756
```

## VERDICT, rule: FAIL (bars 1 and 2, both samples; bar 3 passes)
| | A119b current tag | A119b M (models) | A119b E5 | HG002 current tag | HG002 M | HG002 E5 |
|---|---|---|---|---|---|---|
| tagged transcripts | 1,259 | 36,074 | 13,202 | 621 | 20,222 | 9,777 |
| bar 1: inside one gene (<= 0.10) | 0.526 | 0.515 | 0.380 | 0.728 | 0.700 | 0.566 |
| inside / (inside + join) | 0.699 | 0.910 | 0.821 | 0.807 | 0.945 | 0.928 |
| bar 2: a + b1 untagged (<= 0.50) | 0.915 | 0.472 | 0.703 | 0.917 | 0.483 | 0.756 |
| bar 3: fused tandem junctions tagged | 0 of 60 | 37 of 60 | 30 of 60 | | | |

E5 is better than the models-only comparator on bar 1 and worse on bar 2, so it is not better on both: by the registered
wording it adds nothing over the model rule of row 1712. The signature itself does not discriminate: at one read per site
it covers 0.933 / 0.957 of the real joins and 109,648 / 95,556 transcripts inside one gene. Reads that end and reads that
start inside an intron are a general property of annotated genes (alternative ends and starts, unspliced fragments), at
any read floor tried. Among tags on annotated genes, the current F1 tag is the least often inside one gene (0.699 / 0.807),
because it only asks at introns that are the sole link of a locus.

## RESULT, step 3 families (A119b, descriptive since the rule failed; arm R = E5's 13,202 transcripts as the units list)

```
[job] R assemble 0 families 0 901 s
[audit] R loci 73414; spanning >= 2 coding RefSeq genes: 985 -> shared-exon cluster 527, fused 404, mixed 54
[class] R tandem-class pairs still in one fused locus: 39 of 69
[score] R /anvme/workspace/v138da13-rustle/r1713_ends/R/x.fam.copies.tsv	1963	50438	1779	184	116	0.906	0.813	214
[job] DONE
```

| arm | pair precision | pair recall | loci spanning >= 2 coding genes: fused | tandem-class pairs still in one fused locus |
|---|---|---|---|---|
| D0 default f1v2 (row 1687) | 0.902 | 0.808 | 544 | 69 of 69 |
| E, F1's bridges as units (row 1687) | 0.910 | 0.809 | 539 | 68 of 69 |
| R, read-end form E5 as units | 0.906 | 0.813 | 404 | 39 of 69 |

By the registered bars: no effect on pair precision / recall (+0.004 / +0.005, under 0.01), no harm. Not a bar, but the
largest movement any rule has made on the fused class: fused loci 544 -> 404 and tandem pairs in one fused locus 69 -> 39
(rows 1516 / 1517 reached 60-68, row 1687 65-68; their bar was <= 34). The cost side is not fully measured: 13,202
transcripts are cut into 26,812 units, 5,015 of them inside one annotated gene; the pair score did not fall, the effect on
members per family, on assign and on isoform models of cut genes was not measured. HG002 families effect not run (the
registered condition was a pair-score gain on A119b). No default changed.

## Checks after the first write-up (post hoc)
- **The 69 -> 39 count overstates the separation.** That count only asks whether a fused locus still holds both genes.
  Direct test on the families input (each gene of the pair has a model on an own annotated intron, and the two share no
  locus): D0 28 of 69, E 29, R 40. So arm R separates 11-12 more pairs than the path-matched baseline E, not 30. In all
  arms 18 pairs have a gene with no spliced model of its own. Row 1516's other bar (cuts inside one gene <= 10%) is not met:
  5,015 of R's 13,202 cut transcripts are inside one annotated gene (38%).
- **Loci.** R has 73,414 loci against 75,510 (E) and 74,701 (D0) although more transcripts are cut. Likely because R's
  list does not contain F1's own bridges, so their splits are not made; not verified. Coding genes with a non-joining
  model: 16,495 / 16,497 / 16,513, so no gene is lost. Families: D0 6,445 copies in 1,380, E 6,484 in 1,393, R 6,457 in 1,362.
- **F1's proofs on the joins its table holds (A119b).** 230 of the 1,605 a + b1 joins are in its candidate table; 183 have
  the poly-A proof, 183 the promoter proof, 147 both (= its bridges). So its reach is limited first by the candidate set
  (0.14 of real joins), then by the proofs (0.64 of those).
- On bar 1 as registered E5 beats the current tag (0.380 vs 0.526; 0.566 vs 0.728); among tags on annotated genes the
  current tag is less often inside one gene (0.699 vs 0.821; 0.807 vs 0.928). Both readings are given.
