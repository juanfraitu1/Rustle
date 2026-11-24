# PREREG: the second assembly of the multi-copy stage only where the first left nothing (row 1935, 2026-11-23)

Speed. On HG002 the multi-copy stage assembles the whole genome a second time at a one-read floor (about 1.9 h) and
keeps 128 loci. Its `candidates` step only takes floor-1 transcripts at places the main assembly leaves empty on that
strand, so an alignment overlapping a transcript span of its strand cannot support a kept transcript. Proposed to the
user and accepted: give the second assembly only the other alignments (`tools/empty_place_bam.py`, driver opt-in
`RUSTLE_MULTICOPY_EMPTY=1`; the reads' genome-wide best scores still come from the full BAM's table).
(A first proposal, restricting to the families' own loci, was wrong: the rescued loci are new places that join a family.)

## Why it may not be identical (said before the run)
The assembler decides inside bundles: a chain's share of its locus, collapse into a longer chain, replacement of a weak
junction by a strong neighbour. Removing the reads of a neighbouring occupied place can let a chain through that the
full run collapsed or replaced. So extra candidates are expected, missing ones are not.

## Test (gorilla OR6737; full run on disk: `r1920_mono_from_primaries/real_or/g.multicopy.d/`, `g.multicopy.gtf`)
The stage rerun with the switch on a copy of that run's inputs; compared by intron chain (contig, strand, introns).
- P1: at least 0.99 of the full run's candidates (`cand.gtf`) are produced.
- P2: extra candidates are at most 0.05 of the full count.
- P3: the rescued transcripts of `PREFIX.multicopy.gtf` (tagged `multicopy_rescue`) differ in at most 0.02 of the full
  run's (missing + extra).
- P4: the floor-1 assembly takes under half the time.
All four = run KB3781 and HG002 the same way, then the default is the user's call. Otherwise recorded with the counts.

## Predictions (mine)
P1 and P4 hold. P2 is the risk (I expect 2-10% extras); P3 probably holds, since few extras should be family copies.

## Result (job 4499054, OR6737)
The filter keeps 4,699,246 of 10,705,377 records (0.44), not the small share expected. Stage 613 s against 1,163 s.

| by intron chain | full run | empty places only | shared | missing | extra |
|---|---|---|---|---|---|
| candidates | 11,798 | 18,771 | 11,317 (0.959) | 481 | 7,454 (0.63 of full) |
| rescued transcripts | 1,262 | 5,284 | 952 (0.754) | 310 | 4,332 (3.4 x) |

- **P1 FAILS** (0.959 against 0.99), **P2 FAILS** (0.63 against 0.05), **P3 FAILS** (missing + extra 3.7 x the full run's
  count against 0.02), P4 holds on stage time (0.53). My prediction of P1 and P3 was wrong; extras were far beyond the
  2-10% I expected.
- **Not equivalent, not offered.** The alignments on occupied places are the context that makes the assembler collapse,
  replace or reject chains beside them; without them a one-read floor lets thousands of fragments through, and many of
  those then join a family. The driver switch written for the test is removed; `tools/empty_place_bam.py` stays as a
  recorded arm with the result in its header, called by no stage.
- What is left for the cost of the second assembly: (a) the pipeline's other one-read assembly (the one-read tiers at
  the end of `all`) differs from this one only by a guard setting, and sharing one run would change one of the two
  products (a default, the user's call); (b) accept the cost and run samples side by side.
