# Row 1856: are polyA tails in the reads evidence the assembler should use? (pre-registered 2026-11-21)

User (2026-11-21): there is no rule for polyA tails yet; see whether one can be made. State before this row: the driver's
read check warns when reads end in a polyA run and otherwise the tail is ignored (it is soft-clipped by the aligner).
Of the twelve datasets only the mouse ES library carries tails (60,299 of the first 100,000 primary reads; every other
library 0-121 of 100,000: trimmed before mapping). So a tail rule can be measured on ONE library and would act only on
untrimmed libraries. Measurement only. Workspace `r1856_polya_tails` (`tails.py`). Committed before the script was written
or run. (Registered first as row 1855; renumbered at once: another session had already claimed 1855 with its workspace.)

## Definitions, fixed
A read HAS A TAIL when, at its 3' end in the transcript's orientation, at least 16 of the last 20 sequenced bases are A and
at least 10 bases there are soft-clipped (not aligned: a templated A-run in the genome does not count). Reads = primary
records whose intron chain equals the transcript's. Tailed share of a chain = tailed reads / reads.

## Measured, mouse ES (today's `final.gtf`, mouse contigs, GENCODE vM27)
1. Read level: share of primary spliced reads with a tail; of reads ending in an A-run, the share where it is aligned
   instead of clipped.
2. Multi-exon chains by tailed share (0; under 0.2; 0.2-0.5; 0.5 and over): chains, correct, precision; the same within the
   row 1851 classes of wrong chains (piece, shifted site, combination, novel site, no annotated intron).
3. **Early-ending chains:** the chains `RUSTLE_POLISH_ISM_RATIO_3P=0` adds to today's default on mouse ES (row 1843's class:
   contained in a longer chain, cut at the 3' side only; one new run), by tailed share.
4. For chains with at least 5 reads: tailed share against whether the chain's 3' end is annotated (within 50 bp of an
   annotated transcript end on the strand).

## Signatures and reading, fixed now
- **NOTAIL** = at least 3 reads and none has a tail; **LOWTAIL** = at least 5 reads and tailed share under 0.2.
  A signature is *informative* when the chains holding it are at most half as precise as the rest, they are at least 10% of
  the wrong chains, and at most 2% of the correct chains hold it.
- **Early ends:** tails *separate* them when the added chains with tailed share >= 0.5 are at least twice as precise as the
  other added chains and are at least 20% of the correct added ones.
- An informative signature or a separation makes a rule worth building as an opt-in that switches itself on only when the
  read check finds tails; it could then be validated only on further untrimmed libraries, which we do not have. That limit
  is part of the verdict whatever the numbers.

## Expectation
Tails are added to every polyadenylated molecule, so a full-length read has one and a 3'-truncated or internally primed
read does not: I expect chains with no tailed read to be clearly less often annotated, and the early-ending chains to split
(true alternative ends tailed, truncations not). I am unsure about size: 40% of reads lack a tail, so absence is common.

## Limits
- One library. The simulated sets carry no tails by construction, so no complete truth is available for this question.
- GENCODE is incomplete (row 1851): an unannotated but real transcript with tails counts as wrong.
- The tail test reads the last 20 bases; adapter remnants after the tail would hide it.

## Result (2026-11-21): no rule; nearly every chain already has tailed reads, wrong or right

`r1856_polya_tails/tails.txt` (`tails.py`; one extra run `es/r3p0`, exit 0). Mouse ES, 6,721,797 primary spliced reads.

1. **Reads:** 59.9% carry a tail by the definition. Of the reads ending in an A-run, 0.7% have it aligned instead of clipped:
   the aligner clips tails almost always, so they do not enter the alignments.
2. **Chains by tailed share** (41,638 with reads; correct / all):
   | tailed share | 0 | under 0.2 | 0.2-0.5 | 0.5 and over |
   |---|---|---|---|---|
   | chains, precision | 216 / 1,623 = 0.133 | 52 / 140 = 0.371 | 1,453 / 4,501 = 0.323 | 15,352 / 35,374 = 0.434 |
   Share of each class with a tailed share of 0.5 or more: correct 0.899; combination of annotated introns 0.847; novel site
   0.819; piece 0.804; no annotated intron 0.716. Mean tailed share 0.61 (correct) against 0.55-0.62 (wrong classes).
   - **NOTAIL** (>= 3 reads, none tailed): 454 chains at 0.150 against 0.413 for the rest, but 1.6% of the wrong chains (bar
     10%): not informative. **LOWTAIL**: 225 chains at 0.289, 0.7% of the wrong: not informative.
3. **Early-ending chains** (the 1,207 chains the 3'-only exemption adds, 290 correct): tailed share >= 0.5: 224 / 949 =
   0.236; the others 66 / 258 = 0.256. **Tails do not separate them.** 79% of these chains have a majority of tailed reads.
4. **3' ends** (chains with >= 5 reads; end annotated within 50 bp): no tailed read 31 / 85 = 0.365; under 0.2: 0.643;
   0.2-0.5: 0.748; 0.5 and over: 19,701 / 24,642 = 0.799.

- Expectation: wrong on both counts. Chains without tailed reads are less often annotated, but there are almost none of
  them; and the early-ending chains do not split.
- **What it does say (reading, mine):** the transcripts GENCODE lacks are polyadenylated molecules just like the annotated
  ones: 80-85% of the "wrong" chains of every class have a majority of tailed reads. That fits row 1851 (canonical sites,
  seen by other tools) and speaks against truncation as the source of the early-ending chains on this library: a
  3'-truncated molecule would have no tail. They are real 3' ends the annotation does not list, or internal polyadenylation.
- A tail is good evidence for the END position (item 4: 0.80 against 0.37), not for the chain. The only use the table
  supports is the 3' end layer: on an untrimmed library, place a transcript's 3' end at the tailed reads' end. That is an
  end-accuracy rule, it cannot be tested on a second library here, and the 3' layers were just made opt-in (row 1846).

**Decision:** no filter rule from polyA tails. The read check keeps warning about untrimmed tails. Open, user's call: tails
as 3'-end evidence in the opt-in end layers, once a second untrimmed library exists to check it on.
