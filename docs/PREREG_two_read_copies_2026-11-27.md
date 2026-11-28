# PREREG: family copies with two reads (row 1991, 2026-11-27)

The user's question: can a "two-read copy" be defined, a family copy accepted on two reads, if needed with rules that
keep precision from falling. (My reading of the question; the thin-copy benchmark so far used 3 reads per copy.)

## What exists
- The floor is two reads with the same chain; in a family copy with two reads the two can be one isoform each (the planting
  tool gives a copy's reads 2 : 1 to its two isoforms), so neither passes.
- The multi-copy mode adds one-read transcripts at places the assembly leaves empty when the locus is a family copy.
- The labels (`PREFIX.o2.gtf`) keep a locus only with >= 2 reads by the EM over tied reads (`RUSTLE_O2_MIN_EM`).
- Three-read copies today: isoform recall of the mode 0.65-0.69, of the main assembly 0.56-0.57 (row 1990).

## Step 0 = the measurement (no rule yet)
Planted families, fresh seeds 1994 / 1995, `--reads` 3, 2, 1; scenario E (every copy has the reads) and S (half the copies
silent, the others have the reads). Products of today's default `all`: main assembly, the mode (`multicopy.gtf`), the
labelled set (`o2.gtf`) at `RUSTLE_O2_MIN_EM` 2 (default) and 1.
Reported: copies found, false copies (silent copies given a transcript), isoforms found and isoform precision.
- **T1 (a two-read copy is found):** at 2 reads, scenario E, the mode finds >= 0.80 of the copies and >= 0.50 of the isoforms.
- **T2 (it is safe):** at 2 reads, scenario S, the labelled set at the default threshold gives a transcript to at most
  0.10 of the silent copies.
- **T3 (the threshold matters):** lowering `RUSTLE_O2_MIN_EM` to 1 raises the copies found at 2 reads by >= 0.10 and the
  false copies by less than that.
Prediction: T1 holds for copies and fails for isoforms; T2 holds; T3: copies up, false copies up as much (fails).

## Step 0 result (job 4516462, seeds 1994 / 1995) and a correction of my own scoring
My first table scored `PREFIX.o2.gtf` whole; that file LABELS transcripts (`o2_status`), it does not drop them. Rescored
on the transcripts labelled `expressed`.
| reads per copy | product | E: copies found | S: copies found | S: silent copies given a transcript (of 63) |
|---|---|---|---|---|
| 3 | mode, unlabelled | 0.962 / 0.970 | 0.971 / 0.943 | 42 / 37 |
| 3 | labelled expressed, threshold 2 | 0.571 / 0.602 | 0.571 / 0.571 | 0 / 1 |
| 3 | labelled expressed, threshold 1 | 0.594 / 0.617 | 0.586 / 0.571 | 0 / 1 |
| 2 | mode, unlabelled | 0.955 / 0.955 | 0.943 / 0.929 | 42 / 37 |
| 2 | labelled expressed, threshold 2 | 0.398 / 0.436 | 0.386 / 0.414 | 0 / 0 |
| 2 | labelled expressed, threshold 1 | 0.586 / 0.586 | 0.529 / 0.529 | 0 / 1 |
| 1 | mode, unlabelled | 0.895 / 0.895 | 0.929 / 0.871 | 39 / 37 |
| 1 | labelled expressed, threshold 2 | 0.015 / 0.030 | 0 / 0 | 0 / 1 |
| 1 | labelled expressed, threshold 1 | 0.421 / 0.406 | 0.457 / 0.357 | 0 / 1 |
Isoform recall at 2 reads (E): main 0.447 / 0.417, mode 0.614 / 0.662, expressed at threshold 2 0.241 / 0.303, at 1 0.377 / 0.399.
- **T1 HOLDS** for the unlabelled mode (copies 0.955, isoforms 0.61-0.66; my isoform prediction was wrong), but that set
  gives a transcript to 60-67% of the silent copies at every depth: it is not a safe definition.
- **T2 HOLDS** (0 of 63 in both seeds).
- **T3 HOLDS** (my prediction was wrong): threshold 1 raises the two-read copies found by 0.15-0.19 (S: 0.12-0.14) for
  0 / +1 silent copy. At one read per copy it finds 0.36-0.46 where threshold 2 finds none.
**Found on the way, rule A of row 1990:** at two reads per copy 2 of 63 silent copies got a clean one-read chain in each
seed (none at 3 reads or 1 read; none in row 1990's seeds). The repair named in row 1990, the rival test (no other
alignment of the read within 10 ms), removes all 4 and costs 16 of the 96 correct clean isoforms in the twelve sets.
Built into `one_read_loci.py clean` (the driver passes the score table).

## Confirmation before anything becomes a default (fresh seeds 1996 / 1997; reads 2 and 60)
- **C1:** threshold 1 against 2 at two reads: copies found (labelled expressed) up by >= 0.10 in E and S, silent copies
  with a labelled-expressed transcript up by at most 1.
- **C2:** at 60 reads threshold 1 changes neither the copies found nor the silent copies by more than 1.
- **C3:** with the rival test no clean chain sits on a silent copy, in all eight sets.
- **C4:** genome-wide the rival test costs under 5% of the clean chains (mouse heart, WTC11 r2) and no precision.

## Confirmation result (job 4516574, seeds 1996 / 1997)
Labelled expressed, threshold 2 -> 1 (copies found; silent copies with such a transcript, of 63):
| set | 1996 | 1997 |
|---|---|---|
| 2 reads, E | 0.421 -> 0.564; 0 -> 0 | 0.436 -> 0.564; 0 -> 0 |
| 2 reads, S | 0.429 -> 0.543; 0 -> 0 | 0.486 -> 0.543; 0 -> 0 |
| 60 reads, E | 0.609 -> 0.609; 0 -> 0 | 0.624 -> 0.624; 0 -> 0 |
| 60 reads, S | 0.600 -> 0.600; 0 -> 1 | 0.614 -> 0.614; 1 -> 2 |
- **C1 HOLDS in 3 of 4 sets** (+0.14, +0.13, +0.11; S of seed 1997 +0.057, under the 0.10 bar); no silent copy gained.
- **C2 HOLDS by its bound, with a cost:** at 60 reads nothing is gained and one more silent copy is labelled expressed in
  both seeds (2 more transcripts each).
- **C3 HOLDS:** with the rival test no clean chain is on a silent copy in any of the eight sets (6 chains left out).
- **C4 HOLDS:** mouse heart 644 -> 642 chains, WTC11 r2 558 -> 553, precision 0.5768 and 0.5908 (unchanged). cargo test 1,591 passed.
Even at 60 reads the labelled-expressed set holds only 0.60-0.62 of the expressed copies: the label is strict at every depth;
the threshold is not what limits it there.
State: the rival test is in rule A (working tree). `RUSTLE_O2_MIN_EM` stays 2; 1 is the user's decision. Planted data only.

## Real families, before running (2026-11-27, user: yes, test the threshold on the gorilla with haplotype truth)
KB3781 against its own haplotypes: `r1785_truth_repair/kb/reads.tsv` gives, for 23,710 decidable reads of the selected
families, the locus the read comes from among the loci it aligns to. A locus is TRULY EXPRESSED when at least one decidable
read originates there, and SILENT BY TRUTH when >= 2 decidable reads list it and none originates there.
Today's default `all` on KB3781, labels at `RUSTLE_O2_MIN_EM` 2 and 1; label loci matched to truth loci by overlap.
- **K1:** among the loci labelled expressed at threshold 1 and not at 2 that the truth can judge, at least 0.80 are truly
  expressed, and at most 0.10 are silent by truth.
- **K2:** that share is within 0.10 of the truly expressed share of the loci labelled expressed at threshold 2.
Prediction: K1 fails (real two-read loci are where the labels were brittle before, row 1921), K2 fails.
The truth covers the selected families only; a small number of judged loci is expected and will be stated.

## Real families result (jobs 4516721, 4516999; KB3781)
Labels: threshold 2 = 16,650 loci expressed, 1,619 homology_only, 42 shared; threshold 1 = 17,165 / 1,089 / 57. Gained 515
(505 of them with exactly ONE read placed there alone), lost 0. So on real data the threshold decides about loci holding one
read whose best placement is that locus alone; two-read loci are already expressed unless their reads are tied.
Two judges of mine were unsound and are NOT results: (1) "no read in the truth table originates here" (the table is a
sample of reads, a locus's own reads are often not in it); (2) every primary read overlapping the locus (labelled loci
overlap each other; even homology_only loci came out at 1.00).
Sound judge: the reads the label stage's own rule assigns to the locus (most exon bases on the record's strand), looked up
in the haplotype truth.
| class | loci on truth families | with a read the truth decides | of them with a read truly from the locus | decided reads from here / elsewhere |
|---|---|---|---|---|
| gained at threshold 1 | 49 | 9 | 9 | 9 / 0 |
| expressed at 2, placed 2 | 29 | 1 | 1 | 1 / 0 |
| expressed at 2, placed 3-5 | 75 | 10 | 10 | 21 / 0 |
| expressed at 2, placed >= 6 | 872 | 256 | 256 | 23,506 / 6 |
| homology_only at both | 108 | 0 | - | - |
- **K1 HOLDS, K2 HOLDS, on nine loci** (my prediction was wrong): every gained locus the truth can judge holds a read that
  truly comes from it. Nine is few; 466 of the 515 gained loci lie outside the truth families and are unjudged.
- This truth says a read placed alone is from that locus 99.97% of the time, so it tests the placement, not whether one
  read is enough to call a copy expressed; that part rests on the planted sets (0-1 silent copy of 63 gained).

## Decision (user, 2026-11-28): `RUSTLE_O2_MIN_EM` unset = 1
A family copy is labelled expressed when at least one read is left to it after the tied reads are shared out; 2 = the
labels before. Check (job 4517108): the new default gives the table of the tested threshold-1 run and
`RUSTLE_O2_MIN_EM=2` the table of the earlier default, both byte for byte (planted seed 1996, 2 reads, S: 38 / 20 / 65
against 30 / 6 / 87 expressed / shared / homology_only); cargo test 1,591 passed.
Known cost: at normal depth one more silent copy of 63 labelled expressed per planted set. Real support: nine loci.
Labels made before this date used 2.
