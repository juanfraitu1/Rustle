# PREREG: one read as enough to seed a locus (row 1987, 2026-11-27)

The user's goal: improve the assembler so that even one read can be enough to seed a locus, as StringTie, FLAIR,
IsoQuant and bambu allow; with long reads one read is sometimes all there is.

## What the other assemblers do with one read (their source, clones of 2026-11-25)
| tool | default | where |
|---|---|---|
| StringTie `-L` | ONE spliced read is enough: junction threshold 1, bundle coverage threshold 1; single-exon needs 4.75 | `stringtie.cpp:168-170` (`junctionthr=1`, `readthr=1`, `singlethr=4.75`) |
| FLAIR | 3 reads (`--support 3`); 1 only when asked | `flair_collapse.py:53` |
| IsoQuant | 2 reads for a novel isoform (`default_pacbio`: `min_novel_count` 2); strategies `pacbio_all` and `all` admit single-read novel isoforms, "at the cost of some novel precision" | `isoquant.py:1186-1206` |
| bambu | 2 reads (`min.readCount = 2`) | `bambu-processReads.R:64` |
So one tool of four does it by default and the other three offer it as a setting.

## What Rustle has, and what earlier rows found (not to be rediscovered)
- The floor is 2 distinct records (user's choice, row 1850). A one-read-floor assembly exists (`RUSTLE_MIN_CHAIN_READS=1`
  with `RUSTLE_GATE_MIN_READS=1`); its chains are the second tier of the ranked list, and the multi-copy mode adds its
  spliced transcripts at places the main assembly leaves empty WHEN the locus is a family copy (`multicopy_rescue.py`).
- Row 1805: a lone one-read chain in a locus with no other chain is 0.67-0.72 precise on simulated reads and 2-6%
  annotated on real libraries; junction and end proofs keep too few. Row 1806: no read-level score reaches the precision
  of two-read chains. Row 1867: genome-wide the empty-place candidates are 3-7 times less often annotated than the base.
- What those rows did not look at: (a) the candidate as a LOCUS and not as one more isoform; (b) whether the read is
  there by a clear lead or by a tie (this week: a place the sample puts no read of its own at is a relative's shadow);
  (c) evidence other than the annotation on real data (a replicate library).

## The rule under test: one read seeds a locus, never one more isoform
A spliced one-read transcript of the one-read-floor assembly is added only where the main assembly has NO transcript on
that strand (`multicopy_rescue.py candidates`, without the family condition). Variants, by what the lone read must show:
- **L0:** nothing more (StringTie's condition: one spliced read with canonical introns; the junction gate is already strict).
- **L1:** the read is placed without a rival: no other alignment within 10 `ms` of its best.
- **L2:** L1 and at least two introns.

## Step 0 and test in one run (`bench/o2_truth/one_read_loci.py`)
Non-ape first. Complete truth: fly, mouse, zebrafish simulated. Real: WTC11 replicate 1, judged by the annotation AND by
replication (a read with exactly the chain in replicate 2 or 3, the same cell line). Then mouse libraries and apes as checks.
Per variant: loci added, correct by exact intron chain (simulated) / annotated and replicated (real), and what it does to
the assembly's precision and to the count of correct chains. Beside it StringTie `-L` at the same places.

## Claims
- **O1, a usable class exists on complete truth:** for at least one variant, on at least 2 of the 3 simulated sets, the
  added loci are at least 0.85 correct and at least 100 in number.
- **O2, the tie matters:** on all 3 simulated sets L1's added loci are more often correct than the ones L1 leaves out,
  by at least 0.15.
- **O3, real data:** on WTC11 r1 the loci a variant adds are replicated in r2 or r3 at least half as often as the
  two-read chains of the main assembly are.
- **O4, like StringTie:** at the places where StringTie builds a transcript from one read, the variant that meets O1 is
  at least as often correct as StringTie on at least 2 of the 3 simulated sets.
## Predictions (mine)
O2 holds (the tied ones are the paralog shadows). O1 holds for L1 or L2. O3 fails or is borderline (row 1805's 2-6%
annotated says most real lone reads are something else). O4 holds.

## Result (jobs 4513311-4513314; built option 4513335; fix check 4513355)
| data | L0 loci, correct | L1 loci, correct | left out by L1, correct | main assembly precision | StringTie's chains the main assembly lacks, correct |
|---|---|---|---|---|---|
| fly simulated | 1,027, 0.915 | 999, 0.938 | 28, 0.107 | 0.943 | 380, 0.589 |
| mouse simulated | 2,354, 0.875 | 2,295, 0.885 | 59, 0.508 | 0.930 | 2,153, 0.693 |
| zebrafish simulated | 2,567, 0.804 | 2,433, 0.836 | 134, 0.209 | 0.842 | 1,201, 0.435 |
| WTC11 r1 (annotated) | 9,191, 0.125 | 9,016, 0.126 | 175, 0.074 | 0.565 | 17,630, 0.165 |
WTC11 r1, seen again in r2 or r3: L1 0.199, main assembly 0.895. L2 (>= 2 introns) adds nothing over L1 (0.945 / 0.884 / 0.843; real 0.142).
- **O1 HOLDS** for L1 on fly and mouse (0.938, 0.885; zebrafish 0.836, the level of its main assembly).
- **O2 HOLDS** (+0.83, +0.38, +0.63).
- **O3 FAILS** (0.199 against a bar of 0.448): on a real library four of five such loci are not seen again in two replicates.
- **O4 as written is NOT ANSWERED:** StringTie has only 21-126 of these chains (it builds other ones), too few to compare at the
  same places. Its own extra chains are less often correct on simulated reads (0.44-0.69) and about as often on the real library.
Simulated reads flatter the rule: every simulated read comes from a true transcript.

## Built (opt-in, user to decide on anything further)
`RUSTLE_ONE_READ_LOCI=1` on stage `multicopy` -> `PREFIX.oneloci.gtf` = `PREFIX.gtf` + the L1 loci, tagged `one_read_locus "1"`
(`tools/one_read_loci.py`, selftest). Through the pipeline (`all --to families`, then `multicopy --seed-pool good`):
fly 6,783 -> 7,720 correct at 0.939 -> 0.939; mouse 25,426 -> 27,451 at 0.925 -> 0.922; WTC11 r2 (not used above)
20,095 -> 21,442 annotated at 0.593 -> 0.495 (9,477 loci, 14% annotated). cargo test 1,591 passed.
Not run: mouse libraries, apes, single-exon one-read loci.

## Found on the way: a race in the default of row 1984 (commit 601efd7f)
Under `all` the layered assembly ends in the background; with own-evidence seeding its input `PREFIX.seedpool.d/pool.bam`
was removed first. Fly: the layered run failed (WARNING, no `PREFIX.final.gtf`); mouse and WTC11 finished by timing.
`PREFIX.gtf` is not affected. Fix: the inner call waits for the layered run (`RP_LAYERS_BG=0`); checked on fly and mouse
(no warning, `final.gtf` made, `PREFIX.gtf` identical), cargo test 1,591 passed.

## Addendum before looking (2026-11-27): is there a class of one-read loci that holds on REAL data?
O3 failed for L1 as a whole. Rows 1805 / 1806 looked for proofs and scores of one-read chains against the ANNOTATION; here
the judge is replication (a read with exactly the chain in another library of the same cell line) and the unit is a locus
at an empty place. Exploration on WTC11 r1 (replicates r2, r3) with features of the lone read only
(`bench/o2_truth/one_read_features.py`): introns, shortest and longest intron, shortest end exon in the read, soft clips,
divergence (`de`), MAPQ, splice motifs, distance to the nearest assembled transcript and whether one lies on the other strand.
A rule is then written down from r1 and judged ONCE on WTC11 r2 (replicates r1, r3) and on one mouse pair.
- **O5:** a rule from r1 keeps at least 500 loci on r2 that are seen again at >= 0.60 (two thirds of the main assembly's
  rate), and annotated at least as often as the main assembly's chains are precise minus 0.10.
Prediction: fails; the best class reaches about 0.4 (row 1806: read-level features rank but do not separate).
