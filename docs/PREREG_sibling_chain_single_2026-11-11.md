# PREREG r1628: the single-read form of the sibling-chain rule, tested on A119b and OR6737 (2026-11-11)

r1627 stopped on yield at step 0 (gorilla KB3781: 13 correct of 173). Seen there after the numbers: the chains
carried by one clean read were right 12 of 93 times (0.129, baseline 0.070), while those with two or more were chains
the polish had dropped, rightly (1 of 80). The user chose to register this form and test it on A119b
(AskUserQuestion, 2026-11-11). Written before any A119b or OR6737 run of it. KB3781 is the sample the rule was read
off, so it is not a test here.

## Frozen rule S1
A spliced chain absent from the default GTF is accepted at family member B when
1. a default multi-exon model of another member of the family projects onto B with exactly this intron chain
   (r1519 parameters, as r1627), and
2. exactly one primary record at B has this chain, and its molecule has no rival placement (one record, or second-best
   AS below 0.98 of the best).
So: a chain one read short of the 2-read floor, where the family supplies the second witness.

## Tests (same tools as r1627; correct = exact RefSeq intron chain; families in the censat exclusion are out)
- A119b (`audit/g2/r76e` products): B1 correct chains in S1 >= 100; B2 share correct in S1 >= share correct among
  default 2-read chains at the same loci; B3 share correct in S1 >= 1.5 x share correct in N (projected, no record).
- OR6737 (a fresh default `all` run with the same binaries): B1 >= 15; B2; B3.
- PASS = all three on both samples. PASS means: worth an opt-in layer that writes these chains to a separate GTF
  (default products untouched); making it default is the user's call. FAIL on either = recorded, not ported.
- Descriptive: the classes of r1627 (S, T, X, N), and S1's share among all default chains at the loci.

## Outcome (2026-11-11, jobs 4465550 and 4465551, `r1628/`): FAIL as registered (yield on A119b)

| sample | S1 chains | correct | share | 2-read baseline | no-read (N) | B1 | B2 | B3 |
|---|---|---|---|---|---|---|---|---|
| A119b (2,582 copies after the exclusion) | 306 | 20 | 0.065 | 0.045 (191 of 4,241) | 0.004 | **FAIL** (20 < 100) | PASS | PASS |
| OR6737 | 113 | 16 | 0.142 | 0.111 (223 of 2,013) | 0.016 | PASS (16 >= 15) | PASS | PASS |

- The quality holds on all three samples: a chain one read short of the floor, with a sibling's model as second
  witness, is right more often than the chains the assembler already accepts on two reads (1.3 to 2 times), and
  chains with no read are almost never right.
- The yield does not: 20 correct chains on A119b, 16 on OR6737, 12 on KB3781. The A119b bar of 100 was set from the
  r1617 bar and was far too high for what family members hold: after the exclusion only 2,582 copies are targets.
- Reads with a rival placement (T) are worse than the baseline on A119b and both gorillas, so leaving them to the
  assignment layer was right.
- Not ported. The rule is sound but adds a few dozen chains per sample.

## Ceiling of structure sharing at these parameters (after r1629, from the same files)
RefSeq chains at family members that the default GTF misses, and how many of them ANY sibling projection reaches
(classes S + T + X + N, whatever the read rule):

| sample | RefSeq chains at members | missed by default | reached by a sibling projection |
|---|---|---|---|
| KB3781 | 3,472 | 2,518 | 51 (2.0%) |
| A119b | 7,078 | 5,451 | 78 (1.4%) |
| OR6737 | 5,283 | 3,910 | 88 (2.3%) |

So a sibling's assembled structure can supply at most about 2% of what the assembler misses at family members. This
closes additive structure sharing at cover >= 0.9, not just the two rules. Not measured: the one-read chains at the
same members with NO agreeing sibling (the arm that isolates the sibling's lift); r1629 puts that lift at 0.057 ->
0.088 for 2-read chains, so "valid second witness" should be read as modest.

## Port (2026-11-11)
Ported 2026-11-11 (user decision, opt-in despite the failed yield bar): `tools/rustle_pipeline.sh siblings` -> PREFIX.sibling.gtf = PREFIX.gtf + r1628 chains (`sibling_chains.py emit`; never part of `all`; needs python3 + pysam). KB3781 through the driver: 82 chains added, 12 RefSeq-correct (0.146 vs 0.070 for 2-read default chains), default models byte-identical, 67 s. Uncommitted.
Port on the other two samples through the driver (2026-11-11): OR6737 113 added, 16 correct (0.142), 78 s, equal to the prototype; A119b (run r76e) 393 added, 24 correct (0.061), 1,263 s on 36 threads (prototype on the earlier run: 306 / 20 / 0.065); default models byte-identical in both.
