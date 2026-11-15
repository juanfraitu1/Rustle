# PREREG r1627: a sibling's structure as the second requirement for a low-support chain (2026-11-11)

User goal: "improve the family-wide assembly to share information across the family members". Written before any run.

## Idea
The assembler drops a spliced chain carried by fewer than 2 reads (r1614: 69% of the chains StringTie or FLAIR find
and we miss have one full-chain read). A single read alone is not enough (r1616, r1617). Here the second requirement
comes from the family: a chain at member B is accepted when
1. a default multi-exon model of another member A of the same family, aligned as a spliced template onto B's locus
   (r1519's frozen parameters: minimap2 splice:hq, cover >= 0.9, identity >= 0.85, same strand, within 5 kb of the
   locus), gives exactly this intron chain, and
2. at least one read with a primary alignment at B has exactly this chain and has no rival placement (the molecule has
   one record, or its second-best AS is below 0.98 of its best: `BAM.as_table.tsv`), and
3. the chain is not in the default GTF.
Call this set S. No parameter is tuned.

## How this differs from what is closed
- r1519 / r1522 emitted every projected model; here the projection is only a test on a chain a read at B already shows.
- FAMASM (r1479) pooled reads across copies and its certificate was not selective; here no read is moved and reads
  with a rival placement do not count (they belong to the assignment layer).
- r1617 asked for pair support and known ends at the same locus; here the support is a sibling's model.

## Controls (both registered)
- **Read origin:** chains whose only primary records have a rival placement (class T) are reported and never accepted.
- **Circular truth:** RefSeq models of paralogs are partly built by homology, so a projected chain can match the
  annotation whether or not B is expressed. Class N = projected chains with no record at B at all. If S is not
  clearly more often correct than N, the read adds nothing and this is guided assembly under another name.

## Measure
Correct = the intron chain equals a RefSeq chain on the same contig and strand. Baseline = the sample's default chains
carried by exactly 2 reads at family-member loci (same loci as the targets). Families excluded by the censat BED on
A119b (`a1.assign.excluded.txt`) are out. New chains only; nothing is dropped, so there are no lost chains.

## Bars
Step 0, gorilla KB3781 (`audit/g2/kb_r76` products; 941 copies):
- B1 (stop rule): correct chains in S >= 15. Below that: STOP, no A119b run.
- B2: share correct in S >= share correct among the baseline 2-read chains.
- B3: share correct in S >= 1.5 x share correct in N.
Held out, A119b (`audit/g2/r76e` products), same rule, run only if step 0 passes all three:
- B1: correct chains in S >= 100. B2 and B3 as above.
PASS both samples = worth an opt-in layer (a separate GTF, default products untouched); the default is the user's call.
FAIL = recorded, not ported.

## Tools
`bench/project/project_models.py` (adds the no-rival count per chain) and `bench/famcmp/sibling_chains.py` (classes
and scores). No Rust.

## Outcome, step 0 (2026-11-11, job 4465543, `r1627/kb`, 57 s): STOP as registered
323 families, 911 target copies, 4,275 templates; projected chains absent from the default GTF, by class:

| class | chains | correct | share |
|---|---|---|---|
| S (the rule) | 173 | 13 | 0.075 |
| S with exactly 1 clean read | 93 | 12 | 0.129 |
| S with 2 or more clean reads | 80 | 1 | 0.013 |
| T (primary records, all with a rival) | 90 | 3 | 0.033 |
| X (AS-tied secondary records only) | 130 | 7 | 0.054 |
| N (no record) | 2,556 | 28 | 0.011 |
| baseline: default 2-read chains at the same loci | 1,284 | 90 | 0.070 |

- **B1 FAIL: 13 correct < 15.** The stop rule fires; A119b was not run.
- B2 PASS (0.075 >= 0.070). B3 PASS (0.075 >= 1.5 x 0.011): the read does add something; chains with no read are
  almost never right, so the truth is not circular here.
- Seen after the numbers, not a result: the single-read part of S is right 0.129 of the time, nearly twice the
  2-read baseline, and carries 12 of the 13. The chains with 2 or more clean reads are ones the default polish
  already dropped, and it was right to (1 of 80). Reads with a rival (T) are worse than the baseline.
- The limit is yield, as in r1522: 13 chains on a sample whose default has about 26,000 correct chains.
