# PREREG: a label for processed copies from high-score multimappers (row 1962, 2026-11-25)

The user's question: the near-tied group is useful for discovering loci; are there other useful rules in the high-score
multimappers? A census on OR6737 (workspace `multimapper_rules_2026-11-25/census.or.txt`, description, no bar) showed one
large pattern that no rule uses: the same read aligned SPLICED at one place and INTRONLESS at another (38,852 reads with
a secondary at 0.9 or more of the primary's `ms`; 334 intronless places with 2 or more reads, 181 on an annotated
single-exon gene or pseudogene, 65 unannotated). Row 1955 saw the same pattern on planted reads (AUC 0.90) and did not
turn it into a rule. Seeding from these alignments is NOT the question (rows 1953, 1959, 1961).

## The rule (fixed here; 0.9 and 2 reads are the census's values, not tuned)
- A read with exactly one primary alignment P. Every secondary S with `ms` >= 0.9 x P's `ms` makes a pair (the census
  looked at the best S only; the rule takes all of them, because one parent has many processed copies).
- A pair counts when one alignment has no intron, the other has at least one, and the two do not overlap on the reference.
- A PLACE = intronless alignments of such pairs that overlap each other, from at least 2 reads.
- Label: `processed_copy_of` = the gene under the spliced alignments of most of the place's reads, and `own_reads` = the
  reads whose primary is the intronless alignment.
- It labels; it removes nothing and seeds nothing.

## Data
- Real, not used for the census: human testis A119b on CHM13 (`reads/HSA/A119b.t2t.bam`). Truth: workspace
  `annot_related/chm13.related.gtf` (biotypes `processed_pseudogene`, `transcribed_processed_pseudogene`,
  `translated_processed_pseudogene`; attribute `pseudogene_parent` from NCBI Gene) and `chm13.related.relations.tsv`.
- Planted, read origin known: row 1909's scenario G (150 genes expressed, their 417 processed pseudogenes silent;
  `r1909_silent_members/G.p0.1.bam`). Seen before by row 1955, said here.
- Script `bench/o2_truth/processed_copy.py`, workspace `r1962_processed_copy_label`.

## Claims
- **L1, the places are processed copies (A119b):** at least 80% of the labelled places lie on an annotated processed
  pseudogene or on a gene all of whose transcripts have one exon (50 exon bases or more).
- **L2, the parent is right (A119b):** of the labelled places on a pseudogene with a listed parent, the label names that
  parent in at least 80%, and the parent or a listed relative of it in at least 90%.
- **L3, reach (A119b):** of the annotated processed pseudogenes whose listed parent has 50 or more spliced primary reads,
  at least 30% are labelled.
- **L4, planted:** of the labelled planted pseudogenes, the label names the planted source gene in at least 95%
  (the reach among the 417 is reported, no bar).
- Descriptive: the labelled places by `own_reads` share, and how many carry a transcript of today's assembly
  (`r1920_mono_from_primaries/real_a1/g.gtf`).

## Predictions (mine)
L1 fails narrowly (65-75%: the rest on exons of spliced genes and on nothing annotated). L2 holds (about 85% / 95%).
L3 is near the bar (20-35%): old processed copies have diverged below 0.9. L4 holds. Most places have under 5% own reads.

## Reading
L1 and L2 hold: the label is a candidate (opt-in; a default is the user's call) and the unannotated places are candidate
processed copies with a named parent. L1 fails, L2 holds: the pattern names a relative, not a processed copy. L3 only
sizes the rule. One real sample; a gorilla sample has no parent truth.

## Result (job 4508758, 73 min; tables in workspace `r1962_processed_copy_label/{planted,a119b}.{txt,places.tsv}`)
| claim | bar | result |
|---|---|---|
| L1 places are processed copies (A119b) | >= 80% | **FAILS: 13.8%** (948 of 6,886: processed pseudogene 505, single-exon gene 443, exon of a spliced gene 1,010, nothing annotated 4,928) |
| L2 parent right (A119b) | >= 80% / >= 90% | **HOLDS: 87.1% / 94.5%** (269 and 292 of 309) |
| L3 reach (A119b) | >= 30% | **FAILS: 5.4%** (283 of 5,226) |
| L4 planted source gene | >= 95% | HOLDS on 4 of 4 (4 of 417 planted pseudogenes labelled; planted places on a processed pseudogene or single-exon gene 32 of 83) |
- My predictions: L1 "fails narrowly at 65-75%" was far too high; L3 "20-35%" was too high; L2 and L4 right.
- Where L1 is lost: 4,466 of the 6,886 places have NO annotated gene at the spliced side (4,352 of them on nothing
  annotated at the intronless side either; chrY 1,088 places, the chr13 and chr21 short arms hold the largest by reads):
  repeat-rich regions where one of two alignments has a gap of intron length. Post hoc, places with a gene at the
  spliced side only: 926 of 2,420 (38%), still far under the bar.
- The census on OR6737 (54% on single-exon genes or pseudogenes) did not transfer to the human sample.
- Reading, as registered: the pattern names a relative correctly where a listed parent exists, and it is neither a
  processed-copy detector nor a way to reach them (old processed copies align under 0.9 of the best score). Not adopted.
  Do not re-offer an intronless-against-spliced label without a repeat filter and a stricter spliced side; even then
  the reach is about 5%.
