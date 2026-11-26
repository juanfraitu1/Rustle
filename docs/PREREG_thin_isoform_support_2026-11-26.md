# PREREG: relatives' reads as support for thin isoforms, measured over a whole gorilla sample (row 1965, 2026-11-26)

Row 1961: a secondary alignment at 0.9 or more of the read's best score, admitted only where the place has reads of its
own and only with introns those reads show (`pool_bam.py --rule janchor:S:0.9:K`), is safe on planted reads and gave 0-2
confirmed copies per real sample on the 49 family copies. The user: try that, adding support to thin isoforms, and
measure how beneficial it would be on a smaller sample (not human).

## Samples
OR6737 (gorilla, the sample of row 1961's clear case) and KB3781 (gorilla, second sample).

## Measure (`bench/o2_truth/thin_isoforms.py`), over the whole annotation, not the 49 family copies
An annotated ISOFORM = the intron chain of a spliced transcript of the NCBI gorilla annotation (`GGO_genomic.gff`).
It is ASSEMBLED when a transcript of the assembly has exactly that chain. Its OWN READS = primary alignments with
mapping quality >= 1 whose intron chain is exactly that chain (full-length reads placed there without a tie-break).
THIN = 1-4 own reads. Reported by own-read class (0, 1-2, 3-4, 5-9, 10 or more): annotated isoforms, assembled today,
assembled by each arm, gained, lost; and the share of the assembly's spliced transcripts that match an annotated isoform.

## Arms
Today's assembly; `janchor:as:0.9:3` and `janchor:ms:0.9:3` (built in row 1961); new: `janchor:as:0.9:1` and
`janchor:ms:0.9:1` (one own read at the place suffices: reaches thinner loci).

## Claims
- **T1, a net gain:** on both samples and both scores (K = 3) the annotated isoforms gained are at least twice those lost.
- **T2, no cost:** the share of spliced transcripts that match an annotated isoform falls by at most 0.002.
- **T3, the thin class is what gains:** at least half of the gained isoforms have 1-4 own reads.
- **T4, size:** at least 5% of the thin isoforms that today's assembly misses are assembled (K = 3, either score).
- **T5, a lower anchor helps:** K = 1 gains more annotated isoforms than K = 3 with a matching share within 0.002.
## Predictions (mine)
T1 and T2 hold. T3 holds. T4 fails (under 2%: most missed thin isoforms have no relative with the same introns).
T5: more gained, share slightly lower; uncertain.

## Result, the arms of row 1961 over the whole annotation (jobs 4509365 / 4509366, first part)
| | OR6737 | KB3781 |
|---|---|---|
| annotated isoforms (distinct spliced chains) | 95,810 | 95,810 |
| assembled today / share of today's spliced transcripts that match | 24,698 / 0.3400 | 23,582 / 0.3837 |
| strict rule, AS, K = 3: gained / lost / share | 37 / 20 / 0.3380 | 27 / 13 / 0.3814 |
| strict rule, ms, K = 3: gained / lost / share | 39 / 22 / 0.3379 | 26 / 12 / 0.3814 |
| anchored, AS: gained / lost / share | 66 / 21 / 0.3340 | 46 / 13 / 0.3784 |
| open AS 0.90: gained / lost / share | 141 / 22 / 0.3249 | 117 / 14 / 0.3678 |
| of the gained (strict, AS): 1-4 own reads / no own read | 19 / 16 | 20 / 6 |
| thin isoforms (1-4 own reads) not assembled today, and assembled by the strict rule | 10,855, 19 (0.18%) | 9,274, 20 (0.22%) |
Annotated isoforms by own full-length reads, OR6737 (assembled today): none 59,529 (277); 1-2 13,083 (3,105); 3-4
4,327 (3,450); 5-9 4,697 (4,096); 10 or more 14,174 (13,770).
- **T1** holds on KB3781 (27 v 13, 26 v 12), FAILS on OR6737 (37 v 20, 39 v 22). **T2** borderline: -0.0020 / -0.0021
  on OR6737, -0.0023 on KB3781 (bar 0.002: fails on three of four). **T3 HOLDS** (19 of 37, 22 of 39, 20 of 27, 19 of 26).
  **T4 FAILS** (0.2%).
- Reading: net about +15 annotated isoforms a sample (0.06% of those assembled), for about 400-470 more spliced
  transcripts. The pool of thin isoforms is large (about 10,000 annotated isoforms with 1-4 full-length reads of their
  own are not assembled), and relatives' reads reach 0.2% of it: most of those genes have no relative to borrow from.

## What speaks for the unannotated transcripts the arms add (user: how can we prove the extra ones are real; job 4509372,
`bench/o2_truth/extra_transcripts.py`; description)
EXTRA = an intron chain of the arm that today's assembly and the annotation lack. Evidence: a full-length read of the
sample's own with exactly the chain (primary, mapping quality >= 1); the same in the other gorilla; the chain assembled
in the other gorilla by today's pipeline; every intron annotated. (Every intron of every set is canonical, 100%: the
aligner and the assembler already enforce it, so that check says nothing.)

| set | chains | own full-length reads: 0 / 1-2 / 3+ | full-length read in the other gorilla | assembled in the other gorilla | any of these | every intron annotated |
|---|---|---|---|---|---|---|
| OR6737, today's unannotated | 47,941 | 2,445 / 16,597 / 28,899 | 51.4% | 28.6% | 95.4% | 41.3% |
| OR6737, extra of the strict rule (AS) | 887 | 584 / 290 / 13 | 16.0% | 7.9% | 42.5% | 42.1% |
| OR6737, extra of the strict rule (ms) | 921 | 610 / 298 / 13 | 16.4% | 7.9% | 42.2% | 41.6% |
| OR6737, extra of the anchored rule | 2,020 | 1,690 / 317 / 13 | 8.2% | 4.6% | 21.6% | 20.8% |
| OR6737, extra of the open 0.90 pool | 4,700 | 4,182 / 505 / 13 | 4.0% | 2.8% | 14.3% | 12.8% |
| KB3781, today's unannotated | 37,875 | 1,073 / 9,999 / 26,803 | 55.1% | 36.2% | 97.8% | 39.6% |
| KB3781, extra of the strict rule (AS) | 632 | 424 / 190 / 18 | 17.7% | 12.8% | 45.6% | 42.2% |
| KB3781, extra of the anchored rule | 1,247 | 1,012 / 217 / 18 | 10.3% | 7.9% | 27.3% | 23.6% |
- About 42-46% of the strict rule's extra transcripts are shown by a full-length read at their locus (this animal or
  the other) or assembled independently in the other animal; today's unannotated transcripts reach 95-98%. The rest
  (about 510 on OR6737, 344 on KB3781) have no full-length read anywhere: each of their introns is in a read of the
  locus (the rule requires it), the combination is not. Not shown real, not shown wrong.
- The evidence falls with the looseness of the rule: 42% strict, 22-27% anchored, 14% open.
- Not available here: junction counts from short-read RNA-seq of the same tissue, which is what would settle the rest.

## Result, K = 1 beside K = 3 (same jobs, second part)
| strict rule | OR6737: gained / lost / share / spliced transcripts | KB3781: gained / lost / share / spliced transcripts |
|---|---|---|
| AS, K = 3 | 37 / 20 / 0.3380 / 73,112 | 27 / 13 / 0.3814 / 61,859 |
| AS, K = 1 | 44 / 20 / 0.3380 / 73,137 | 30 / 13 / 0.3812 / 61,903 |
| ms, K = 3 | 39 / 22 / 0.3379 / 73,135 | 26 / 12 / 0.3814 / 61,873 |
| ms, K = 1 | 46 / 22 / 0.3379 / 73,163 | 29 / 12 / 0.3811 / 61,927 |
**T5 HOLDS** (more gained on both samples and both scores, share within 0.0003). The extra isoforms of K = 1 are in the
class with no full-length read of their own (OR6737: 22 against 16; KB3781: 8 against 6): one read at the place is enough
to let a relative's read in, and the isoform it then supports is one that no read of the place shows in full.
The 49 family copies are unchanged by K (38 and 31 modelled, junction precision 0.65 and 0.63).

## Summary of row 1965
- Benefit of the strict rule over a whole gorilla sample: net +14 to +24 annotated isoforms (0.06-0.1% of those
  assembled), for 400-530 more spliced transcripts, of which 42-46% are shown by a full-length read at their locus or by
  the other animal. T1 mixed, T2 borderline, T3 and T5 hold, T4 fails.
- The pool it was meant for is large and almost untouched: about 10,000 annotated isoforms with 1-4 full-length reads of
  their own are not assembled today; relatives' reads reach 0.2% of them. Thin isoforms are mostly a matter of the
  assembler's own support thresholds, not of borrowing from relatives.
