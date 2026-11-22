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
