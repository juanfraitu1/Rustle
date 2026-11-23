# PREREG: a member-discovery layer, the relatives of every expressed locus through all of its reads' secondaries (row 1910, 2026-11-22)

User (2026-11-22): "lets do that, the member discovery layer". The advisor's position: when one member of a family is
found the others should be available, a silent gene through its expressed pseudogene and back.
Written before any place is counted. Seen: row 1909 (silent loci of a fixed list reached by >= 2 reads: pseudogenes 70%,
genes 81%, mostly by records under half of the best score; 1,061 other places reached, unjudged).

## The layer (`tools/member_discovery.py ASSEMBLY.gtf BAM OUT_PREFIX`)
A read's SOURCE is the locus of its best-scoring record (AS). Every other record of the read, primary or secondary, that
does not lie on the source locus is a hit; hits of one strand are merged into PLACES (overlapping or within 1 kb). A
place is a member candidate of source S when >= `--min-reads` reads of S hit it with >= `--min-bases` aligned bases and
AS >= `--min-ratio` x the read's best. Output: `OUT.members.tsv` (place, coordinates, strand, source, reads, best and
median ratio, the assembled locus it coincides with or `.`) and `OUT.members.gtf` (one model per place that is not an
assembled locus, exons from its best record, `o2_status "member_candidate"`). Nothing is asserted about expression.

## Truth and control (planted reads, CHM13; `r1910_member_discovery`)
Per draw, 60 reads per gene, every relative silent, loci given to the layer = the expressed genes:
- 150 genes with named processed pseudogenes (`SYMBOL` + `P<number>`),
- 150 genes with a paralogue of last common ancestor Homo sapiens to Catarrhini,
- 150 CONTROL genes: protein-coding, no paralogue at any level of the Ensembl table, no named pseudogene.
A place of source S is: RELATIVE when it covers >= 50 exon bases of an annotated paralogue of S (any level) or of a
pseudogene named after S or after one of its paralogues; PSEUDOGENE, UNLISTED when it covers a gene of a pseudogene
biotype that is not on that list (most pseudogenes carry clone names; plausible relatives, not counted as right);
OTHER GENE (a protein-coding or RNA gene that is not a listed relative); UNANNOTATED.
Grid: min-ratio 0, 0.3, 0.5, 0.8, 0.98; min-reads 2, 5; min-bases 100.

## Procedure and claims
- Draw A (seed 1910) chooses the default: the loosest grid point at which at most 10% of the control genes have any
  place. If none qualifies, the layer has no usable default and that is the result.
- Draw B (seed 1911, drawn independently) tests it:
  - N1 (it finds relatives): at the default, >= 50% of the named processed pseudogenes and >= 50% of the recent
    paralogues of the expressed genes are reached.
  - N2 (it does not invent them): at the default, at most 10% of draw B's control genes have a place, and RELATIVE +
    PSEUDOGENE, UNLISTED are >= 80% of all places of the non-control genes.
  - N3 (near ties are not enough): at min-ratio 0.98 under half as many relatives are reached as at the default.

## Predictions (mine)
Control genes have places at ratio 0 for well over 10% (repeats and shared domains inside transcripts), so the default
lands at 0.3 or 0.5. N1: holds for recent paralogues, borderline for pseudogenes (row 1909: 119 of 417 at 0.5). N2: the
control half holds by construction of the default on A if B agrees; the 80% share is the risk. N3 holds.

## Result (2026-11-22, job 4495094, `r1910_member_discovery/{A,B}`): default = 5 reads at ratio >= 0.3; N1 FAILS narrowly for pseudogenes, N2 HOLDS, N3 HOLDS
Draw A, control genes with a place (of 150): ratio 0 / 2 reads 33, 0 / 5 30, 0.3 / 2 16, 0.3 / 5 14, 0.5 / 2 9, 0.8 / 2 5,
0.98 / 2 2. First setting at or under 15: ratio 0.3 with 5 reads. That is the default.
Draw B at the default (5 reads, ratio >= 0.3):
- places of non-control genes 1,350: listed relatives 969, annotated pseudogenes not on the list 246, other genes 65,
  unannotated 70. Relatives + pseudogenes 90.0% (bar 80%).
- control genes with a place 9 of 150 = 6% (bar 10%). N2 holds.
- reached: named processed pseudogenes 295 of 610 (48.4%, bar 50%), recent paralogues 770 of 1,066 (72%), older
  paralogues 33 of 3,475 (1%). N1 holds for paralogues and FAILS for pseudogenes by 1.6 points.
- near ties only (ratio 0.98, 5 reads): 2 pseudogenes, 334 recent paralogues, 0 older = 336 relatives against 1,098 at
  the default (31%). N3 holds.
- Without a score threshold (2 reads) draw B reaches 416 of 610 pseudogenes (68%) with 24 control genes hit (16%) and
  81.6% of places on relatives or pseudogenes.
Predictions: default at 0.3 right; paralogues pass, pseudogenes borderline right (they fail by 1.6 points); the 80%
share I called the main risk holds at 90%; N3 right.
Since then: the tool's defaults are 5 reads and ratio 0.3, and the driver has the opt-in stage `members`.
Limits: planted reads, sample = reference, every relative silent and every source at 60 reads; "unlisted pseudogene"
places are counted as plausible, not verified; the 1 kb merge can chain neighbouring places into one.

## Real sample, no truth (2026-11-22, jobs 4495096, 4495101; `r1910_member_discovery/or`)
Gorilla OR6737 through the driver's `members` stage: 4,404,440 reads, 5,611,233 hits off their source locus; 3,537 source
loci have a member place (706 have >= 5); 8,244 places, 1,562 on an assembled locus and 6,682 on none (the models of
`g.members.gtf`). First run 16 min and 195 GB in the collecting process (the scan kept every record's blocks); after
reading the blocks back only for the places that get a model: 11 min, 4.6 GB in the largest process, same counts, and
identical models on draw B. Not checked against gorilla's annotation.

# PREREG: the member candidates of a real sample against its annotation (row 1912, 2026-11-22)

Seen: row 1910 (planted: 6% of control genes get a place, 90% of places on listed relatives or annotated pseudogenes);
on gorilla OR6737 the counts only (3,537 source loci with a place, 8,244 places, 6,682 on no assembled locus). Not seen:
what those places are. Workspace `r1912_members_gorilla_annot`; `r1910_member_discovery/or/g.members.tsv` as written.

## Method
Gorilla RefSeq annotation (`GGO_genomic.gff`, loader of `famtruth.py`; human Compara paralogues through the gene's
symbol or its best human protein hit, as rows 1611 / 1855).
- A source locus is the annotated gene its transcripts map to (`map_copy`). WITH RELATIVES: a protein-coding gene with a
  paralogue at any level. CONTROL: a protein-coding gene with a human counterpart and no paralogue at any level.
  Other source loci (no annotated gene, no counterpart, non-coding) are counted and set aside.
- A place of a source is RELATIVE (>= 50 exon bases of an annotated paralogue of the source, any level), ANNOTATED
  PSEUDOGENE (RefSeq gives no parent, so plausible, not verified), OTHER GENE, or UNANNOTATED.

## Claims
- G1: among the places of sources WITH RELATIVES, RELATIVE + ANNOTATED PSEUDOGENE are >= 60%.
- G2: at most 20% of the expressed CONTROL genes have a place.
- G3: the share of source genes with a place is at least 3 x higher WITH RELATIVES than CONTROL.
Descriptive beside them: recent paralogues (Homo sapiens to Catarrhini) of expressed genes reached; places on an
assembled locus against new ones.

## Predictions (mine)
Weaker than on planted reads: G1 fails (about 40-50%: a real sample has repeats inside UTRs and unannotated copies),
G2 borderline (15-25%), G3 holds.

## Result of row 1912 (2026-11-22, `r1912_members_gorilla_annot/result.tsv`): G1 FAILS narrowly, G2 FAILS narrowly, G3 FAILS
| source class | source loci | with a place | place rows | RELATIVE | ANNOTATED PSEUDOGENE | OTHER GENE | UNANNOTATED |
|---|---|---|---|---|---|---|---|
| WITH RELATIVES | 12,362 | 2,088 (16.9%) | 7,025 | 712 | 3,347 | 1,720 | 1,246 |
| CONTROL (no paralogue) | 3,603 | 778 (21.6%) | 2,920 | 0 | 1,543 | 888 | 489 |
| no human counterpart | 462 | 163 | 977 | 0 | 406 | 391 | 180 |
| non-coding or pseudogene | 1,433 | 316 | 2,213 | 0 | 1,006 | 680 | 527 |
| no annotated gene | 1,815 | 192 | 1,349 | 0 | 428 | 304 | 617 |
- G1: relative + annotated pseudogene 57.8% of the places of sources with relatives (bar 60%). Listed relatives alone 10%.
- G2: 21.6% of the control genes have a place (bar 20%).
- G3: genes with relatives have a place LESS often than control genes (16.9% against 21.6%; bar 3 x more).
- Recent paralogues of expressed genes reached: 423 of 2,054 (21%; planted 72%).
Why G3 cannot be read as planned: this control only excludes PARALOGUES. RefSeq does not name a pseudogene's parent, so
genes with processed pseudogenes are in it, and 53% of the control's places are annotated pseudogenes (the planted
control excluded genes with named pseudogenes). A design fault of this row, found after.
What stands: on a real sample 53-58% of member places lie on an annotated relative or pseudogene; 42-47% lie on another
gene or on unannotated sequence, against 10% on planted reads. Those are unexplained here: shared domains and repeats
in transcripts, unannotated copies, or noise. The candidates of `PREFIX.members.gtf` are not to be read as members
without that check; the layer needs a tighter rule (aligned span, not only score ratio) and a control that excludes
pseudogene parents before any claim on real data.
Predictions: G1 fails right (I said 40-50%; 58%), G2 borderline right (21.6%), G3 holds wrong.

# PREREG: the member-discovery layer on a human real sample, what its places are (row 1913, 2026-11-22; written as 1912 before I saw that number taken)

Seen: row 1910 (planted reads: at 5 reads and ratio 0.3, 6% of control genes get a place, 90% of places lie on a listed
relative or an annotated pseudogene); the real gorilla run's counts (6,682 new places, 1,562 on assembled loci, no
classification). Not seen: any place of a human real sample.
Seen after writing the claims below and before the run started (the job was cancelled for the number clash, nothing ran):
row 1912, another session's classification of the gorilla run: relative + annotated pseudogene 57.8%, control genes with
a place 21.6%, recent paralogues reached 21%, all three of its bars failed. The claims and predictions below are
unchanged; read them knowing that.

## Run (workspace `r1913_members_real`)
`tools/member_discovery.py` at its defaults on A119b (human testis, CHM13, the usual `-p 0.1` BAM) with the assembly
`r1869_multicopy_stage/a1/x.gtf`. A source locus is given the annotated gene its transcripts map to (`famtruth.map_copy`,
most transcripts); sources on no annotated protein-coding gene are counted and set aside.
Each place x source is classified as in row 1910: RELATIVE (annotated paralogue of the source's gene at any level, or a
pseudogene named after it or after a paralogue), PSEUDOGENE, UNLISTED (any other gene of a pseudogene biotype), OTHER
GENE, UNANNOTATED. CONTROL sources = genes with no paralogue at any level and no named pseudogene.

## Claims
- P1: of the places of non-control protein-coding sources, RELATIVE + PSEUDOGENE, UNLISTED are at least 70% (planted: 90%;
  the bar is lower because real reads carry repeats and unannotated transcription).
- P2: at most 15% of the control sources have a place (planted: 6%).
- P3: at least half of the annotated recent paralogues (Homo sapiens to Catarrhini) of sources with >= 20 reads are reached.
Descriptive beside them: time, memory, places per source, the largest sources.

## Predictions (mine)
P1 fails or is close (I expect 50-75%: UTR repeats and readthrough into neighbours make OTHER GENE and UNANNOTATED
larger on real reads). P2 fails (I expect 20-40% of control genes hit, through repeats in long 3' UTRs). P3 holds.
If P1 or P2 fail, the next step is a repeat filter on the hit (soft-masked fraction) before the layer is used on real
samples, tested as its own row.

# PREREG: relatedness from the reads, cross-mapping links against the families stage and Compara (row 1914, 2026-11-22)

User (2026-11-22): "is there anything that the reads themselves can tell us? ... can the relationships be modelled
automatically from the reads?" Written before row 1913's member table of A119b exists. Seen: rows 1910, 1912 (place
classes); the families stage's pair scores on A119b (row 1858: 0.650 / 0.639 with near-tie secondaries).

## Measure (A119b, CHM13; workspace `r1914_read_relatedness`, `pairs.py`)
Universe = protein-coding genes with an assembled transcript that are in the Ensembl paralogue table; truth = paralogue
pairs of last common ancestor up to Simiiformes with both genes in the universe (as `pairscore.py`).
- READ LINKS: (gene of the source locus, gene under the place) from row 1913's `a1.members.tsv`, both in the universe;
  at the layer's default (>= 5 reads, ratio >= 0.3), and restricted to places whose best ratio is >= 0.8 and >= 0.98.
- FAMILIES: all pairs of genes inside one family of `r1858_seed_pool_families/a1/good/x.fam.copies.tsv`.
- Union and intersection of the default read links and the families pairs.
Pair precision and recall for each. Also, without a bar: read links whose target gene is NOT expressed (truth pairs with
one expressed gene), the part no expression-based family definition can reach.

## Claims
- X1: read links at the default have pair precision >= 0.5.
- X2: the union has recall at least 0.05 above the families stage alone.
- X3: read links restricted to ratio >= 0.8 have precision at least as high as the families stage.

## Predictions (mine)
X1 fails (0.2-0.4: links through repeats and shared domains among expressed genes). X2 fails (the families stage
compares the same sequences more sensitively; I expect under 0.03). X3 unsure, leaning fail.

## Amendment to rows 1913 and 1914 (2026-11-22, before any human result)
The run of row 1913 has been scanning one contig for over three hours (the other contigs were done within minutes; the
known satellite hot spot of A119b), with 2.5 h of its limit left. A second run is started beside it with
`member_discovery.py --exclude` on the satellite list assign uses (`chm13v2.0_censat_v2.1.satellite_only.bed`): records
starting inside a satellite interval are not read (new option; output unchanged with an empty list). Rows 1913 and 1914
are scored on the full run if it finishes, otherwise on the excluded run, said which. No number of either seen.

## Result of rows 1913 and 1914 (2026-11-22, job 4495226, the satellite-excluded run; the full run was cancelled at 4 h 50 min, still on one contig): the numbers below are NOT a measure of specificity
As registered. Row 1913: P1 holds (relative 26.7% + annotated pseudogene 44.3% = 71.0% of 23,164 places of non-control
sources; bar 70%), P2 holds (610 of 4,714 control sources with a place, 12.9%; bar 15%), P3 holds (2,030 of 3,584 recent
paralogues of sources with >= 20 reads reached, 56.6%). 68,901 of the 97,501 rows are on sources that are no annotated
protein-coding gene. 1 h 17 min, 25 GB.
Row 1914 (16,024 expressed genes, 2,112 truth pairs): read links at the default 41,688 pairs, precision 0.034, recall
0.679; at ratio >= 0.8 0.088 / 0.566; near ties 0.163 / 0.365; families stage 2,116 pairs, 0.816 / 0.818; union 0.046 /
0.931; intersection 0.818 / 0.565. X1 fails, X2 holds (recall + 0.113), X3 fails. Against paralogues of any age the
default links are 0.057 precise.
Why none of this can be read as intended (found after): a PLACE is built from the SPANS of the hits, merged within
1 kb. A spliced alignment's span covers its introns, so places chain neighbouring genes: 15,903 places, median span
21.5 kb, 5,834 over 100 kb, 258 over 1 Mb, the largest 6.0 Mb. A place that covers many genes is counted as RELATIVE
when any of them is one (inflating P1 and P3) and yields a pair with every gene under it (deflating row 1914's
precision). The planted rows are affected less (hits on intronless pseudogenes are short) but the same rule was used:
row 1910's 90% is an upper bound.
What stands: the families stage is 0.82 / 0.82 on this sample; cross-mapping reads reach relatives the families stage
does not (the union's recall), but with the present place rule the links cannot be told from neighbours.
Next (its own row): places from the aligned BLOCKS of the hits, not their spans, and the gene under a hit taken per hit
by exon overlap; then rows 1910, 1912, 1913 and 1914 scored again.
Predictions: P1 and P2 wrong (they hold), P3 right; X1 fails right, X2 wrong (it holds), X3 fails right; none of it
worth more than the measure.

# PREREG: the layer with places from aligned blocks, every earlier claim scored again (row 1915, 2026-11-22)

Rows 1910, 1913 and 1914 were scored with places merged from hit spans, which chains neighbouring genes. Since then
`member_discovery.py` builds a place from overlapping aligned BLOCKS (no slack), takes a read's source from its primary
record, counts secondary records only, and classification uses the place's block footprint, not its span. Truth for a
relative is now the relations table of `related_gtf.py` (Compara paralogues of any age, NCBI pseudogene parents, the
naming rule, copies of one source gene), not the name rule alone.
Seen: the old-rule numbers; one run of the new tool on planted draw B for its span statistics only (1,038 places,
median 1.6 kb, largest 435 kb), no class, no pair.
Thresholds stay 5 reads at ratio 0.3 (chosen under the old rule on draw A: said, not re-chosen).

## Runs (workspace `r1915_members_blocks`, `eval2.py`)
- Planted draws A and B of row 1910 (same BAMs).
- A119b (`--as-table`, `--exclude` satellites), assembly as row 1913.

## Claims, the earlier bars unchanged
- B1 (row 1910's N2 on draw B): <= 10% of control genes with a place; >= 80% of the places of the other genes on a
  listed relative or an annotated pseudogene.
- B2 (N1 on draw B): >= 50% of named processed pseudogenes and of recent paralogues reached.
- B3 (row 1913's P1, P2, P3 on A119b): >= 70%; <= 15% of control sources; >= 50% of recent paralogues of sources with
  >= 20 reads.
- B4 (row 1914's X1 on A119b): gene pairs from the places, pair precision >= 0.5 against paralogues up to Simiiformes;
  reported beside the families stage, their union and intersection.

## Predictions (mine)
B1 holds. B2 holds for paralogues, fails for pseudogenes (as before). B3: the 70% share fails (I expect 50-65% once
neighbours stop counting), the control and paralogue halves hold. B4 fails but far above 0.034 (I expect 0.15-0.35).

## Result of row 1915 (2026-11-22, job 4495547, `r1915_members_blocks`): B1 HOLDS, B2 HOLDS, B3 FAILS on the share, B4 FAILS
Planted, places from blocks, truth from the relations table:
| draw | control genes with a place | places on a relative + annotated pseudogene | recent paralogues reached | pseudogenes of the gene reached | place span median / largest |
|---|---|---|---|---|---|
| A | 14 of 150 (9.3%) | 81.9% + 12.8% = 94.6% | 823 of 1,099 (74.9%) | 267 of 699 (38.2%) | 2.4 kb / 618 kb |
| B | 9 of 150 (6.0%) | 79.4% + 13.2% = 92.7% | 824 of 1,066 (77.3%) | 348 of 691 (50.4%) | 1.8 kb / 435 kb |
B1 and B2 hold on draw B (pseudogenes by 0.4 points; draw A 38%).
A119b (53 min 38 s, 7.6 GB; 97,978 rows, 68,798 on sources that are no annotated protein-coding gene):
- places of non-control sources 24,199: relative 29.7%, annotated pseudogene 18.6%, other gene 5.4%, UNANNOTATED 46.4%.
  Relative + pseudogene 48.3% (bar 70%): B3 fails on the share. Control sources with a place 559 of 4,635 (12.1%,
  holds); recent paralogues of sources with >= 20 reads reached 2,292 of 3,584 (64.0%, holds); pseudogenes of the gene
  3,703 of 6,603 (56.1%).
- gene pairs (2,112 truth pairs among 16,024 expressed genes), precision / recall: read links 0.228 / 0.737; links at
  ratio >= 0.8 0.360 / 0.617; families stage 0.816 / 0.818; union 0.266 / 0.920; intersection 0.821 / 0.635; families +
  links at >= 0.8 0.427 / 0.874. B4 fails (bar 0.5); precision is 0.228, not row 1914's 0.034.
- Places still chain on real data: median span 122 kb, largest 5.4 Mb. A place is a connected set of hits, and spliced
  hits of different reads bridge neighbouring genes through shared exons. The footprint is exact, the unit is still too
  coarse.
Predictions: B1 right; B2 right for paralogues, wrong for pseudogenes on B (50.4%); B3 as predicted (48%); B4 fails
right, 0.228 inside my 0.15-0.35.
Reading: on planted reads the layer finds three quarters of recent paralogues and half the pseudogenes, 93-95% of its
places on known relatives, 6-9% of control genes hit. On a real human sample it reaches 64% of recent paralogues and 56%
of pseudogenes, but 46% of its places lie on no annotated exon and its gene pairs are 23% precise. The families stage,
which compares the assembled sequences, is 0.82 / 0.82 on the same sample: that is the working read-derived relatedness.
Cross-mapping adds recall (0.92 in union) at a precision no better than 0.27-0.43.
Open: links per (source, target locus) instead of connected places; what the unannotated half is (repeats in
transcripts, unannotated copies); the gorilla rescoring with the repaired tool.
