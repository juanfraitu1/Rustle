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

## Result of the addendum (jobs 4513387, 4513406): no class holds on real data
WTC11 r1, 8,976 one-read loci (12.7% annotated, 20.0% seen again in r2 / r3). No feature of the lone read lifts a bin above
0.33 seen again: introns 0.16-0.22; shortest intron 0.13-0.21; longest intron 0.16-0.23; shortest end exon 0.17-0.21;
soft clips 0.12-0.22; divergence 0.16-0.26; MAPQ 60 0.20; a transcript on the other strand 0.185 against 0.205; distance to
the nearest transcript 0.18-0.24; aligned length 0.33 (< 1 kb, 919 loci) falling to 0.13 (> 4 kb). Every locus already has
canonical splice motifs. **O5 FAILS on the exploration set itself; no rule was written, and the r2 table was not opened.**
Calibration, main assembly of the same library by reads of the chain (annotated / seen again): 2 reads 0.300 / 0.681;
3 reads 0.388 / 0.860; 4-5 0.474 / 0.953; 6-10 0.576 / 0.994; > 10 0.844 / 1.000. One read at an empty place: 0.127 / 0.200.
The step from one read to two is the largest on the scale. A transcript present at about one read per library would be
seen again in two more libraries about 85% of the time; 20% says most of these are either artifacts or molecules far
rarer than one per library. The data here cannot tell those two apart.
(`bench/o2_truth/one_read_features.py`; the calibration script is `depth.py` in the workspace.)

## Second addendum before looking (2026-11-27): judge the LOCUS and its introns, and use the sample's other reads
The exact chain may be the wrong judge for one read: long reads were seen again LESS often (0.13 over 4 kb against 0.33
under 1 kb), which is what truncated ends would do to a real transcript. The goal is a real LOCUS.
Measures per one-read locus (`bench/o2_truth/one_read_introns.py`), WTC11 r1 explored, r2 judged once:
- introns seen again: every intron of the chain is in a primary read of r2 or r3 (any chain);
- introns annotated: every intron is in the annotation;
- in-sample support: how many OTHER primary reads of r1 itself share an intron with the chain, and whether every intron
  has one (reads that disagree with each other on the full chain, so the two-read floor never met them).
Calibration: the same measures on the main assembly's two-read chains.
- **O6:** the class "every intron also in another read of the same sample" holds >= 500 loci on r2 whose introns are all
  seen again at >= 0.90 of the rate of the main assembly's two-read chains, and all annotated at >= 0.90 of theirs.
Prediction: holds on replication, borderline on annotation. If it holds it is a rule from one sample's own reads.

## Result of the second addendum (job 4513408): O6 HOLDS
| WTC11 | loci | all introns seen again | all introns annotated | exact chain annotated | exact or a contiguous part of an annotated chain |
|---|---|---|---|---|---|
| r1 (explored): every intron also in another own read | 1,917 | 0.810 | 0.636 | 0.122 | 0.550 |
| r1: the other one-read loci | 7,059 | 0.25-0.32 | 0.21-0.27 | 0.128 | 0.206 |
| r1: two-read chains of the main assembly | 11,118 | 0.896 | 0.543 | | |
| **r2 (judged once)**: every intron also in another own read | 2,243 | 0.942 | 0.775 | 0.124 | 0.708 |
| r2: the other one-read loci | 7,234 | 0.33-0.48 | 0.26-0.32 | 0.148 | 0.249 |
| r2: two-read chains of the main assembly | 9,292 | 0.945 | 0.605 | | |
The class is as reproducible as the two-read chains and its introns are MORE often all annotated. Its chain is rarely an
annotated chain exactly (12%): six in ten are a contiguous part of one. These are real genes seen through reads that each
stop somewhere else, so no two agree on a full chain; the locus is real, the model is a fragment.
Built: tag `one_read_locus "supported"`, and `RUSTLE_ONE_READ_LOCI=supported` writes only those.

## Third addendum before running: the supported class on other species (rule fixed, nothing tuned)
`bench/o2_truth/one_read_score.py` on fly / mouse / zebrafish simulated, four mouse libraries, two gorillas.
- **O7:** on at least 3 of the 4 mouse libraries the supported loci have all introns annotated at >= 0.90 of the rate of
  the same library's two-read chains, and on each simulated set their chain is exact or a contiguous part of a true
  chain for >= 0.90. Gorillas are reported as a check, not part of the claim.
Prediction: holds.

## Result of the third addendum (jobs 4513436, 4513677): O7 HOLDS on 4 of 4 mouse libraries and all 3 simulated sets
Columns: loci; chain exact; exact or a contiguous part of an annotated chain; all introns annotated.
| data | supported one-read loci | two-read chains of the main assembly | other one-read loci |
|---|---|---|---|
| fly simulated | 25; 0.200; 1.000; 1.000 | 625; 0.915; 0.960; 0.973 | 973; 0.958; 0.985; 0.985 |
| mouse simulated | 144; 0.465; 1.000; 1.000 | 1,769; 0.881; 0.964; 0.970 | 2,143; 0.914; 0.970; 0.970 |
| zebrafish simulated | 170; 0.282; 0.965; 0.971 | 1,366; 0.795; 0.908; 0.919 | 2,256; 0.879; 0.949; 0.954 |
| mouse heart | 1,211; 0.173; 0.736; 0.842 | 4,863; 0.289; 0.425; 0.604 | 4,268; 0.257; 0.359; 0.414 |
| mouse adrenal | 1,694; 0.140; 0.773; 0.870 | 5,404; 0.273; 0.454; 0.639 | 4,691; 0.229; 0.369; 0.427 |
| mouse cortex | 1,461; 0.161; 0.636; 0.763 | 7,174; 0.206; 0.314; 0.544 | 5,499; 0.161; 0.229; 0.276 |
| mouse ES | 1,304; 0.100; 0.474; 0.574 | 8,027; 0.121; 0.227; 0.368 | 6,144; 0.101; 0.155; 0.174 |
| gorilla OR6737 (check) | 4,276; 0.078; 0.491; 0.531 | 18,219; 0.165; 0.393; 0.440 | 21,241; 0.039; 0.070; 0.074 |
| gorilla KB3781 (check) | 1,665; 0.094; 0.540; 0.593 | 11,423; 0.157; 0.360; 0.419 | 7,698; 0.061; 0.114; 0.121 |
On every real library the supported loci have all introns annotated more often than the two-read chains the assembler
already emits (ratio 1.2-1.6), and are 2-7 times better than the other one-read loci. Their EXACT chain is annotated
less often than a two-read chain's (0.08-0.17 against 0.12-0.29): they are fragments of real transcripts.
(The first ES submission failed on a flag the driver does not have; rerun without it.)

## Decision (user, 2026-11-27): the supported one-read loci go INTO `PREFIX.gtf` by default
`all`: after the main assembly, `stage_oneloci` adds them (tag `one_read_locus "supported"`); `RUSTLE_ONE_READ_LOCI=0` = off.
`assemble` named alone and `PREFIX.final.gtf` are unchanged. The floor-1 run it needs is the one the one-read tier reuses.
Check (job 4515330), exact intron chains of `PREFIX.gtf` before -> after:
| data | chains | correct | precision | loci added |
|---|---|---|---|---|
| fly simulated (whole `all`) | 7,221 -> 7,246 | 6,783 -> 6,788 | 0.9393 -> 0.9368 | 25 |
| mouse heart | 22,638 -> 23,849 | 13,086 -> 13,295 | 0.5781 -> 0.5575 | 1,211 |
| WTC11 r2 | 33,878 -> 36,121 | 20,095 -> 20,374 | 0.5932 -> 0.5640 | 2,243 |
The cost in exact-chain precision (0.02-0.03 on real libraries) is the fragments; filter on the tag to get the assembly
before. No warning, `PREFIX.final.gtf` made, the one-read tier reused the run, a second `all` redid nothing and kept the
count. cargo test 1,591 passed. NOT measured: the effect of the added loci on families and the later stages (they ran).
Every `PREFIX.gtf` made with `all` before this date lacks these loci: rerun or filter before comparing.
