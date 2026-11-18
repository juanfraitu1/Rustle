# PREREG r1775: the adapted Eichler filter as the baseline for reads of multi-copy families, and Winnowmap as an arm (2026-11-17)

> Naming corrected 2026-11-18: "Eichler's rule" here always meant our adaptation to one reference of the Eichler lab's 10-point filter (Guitart 2024); the paper's own procedure is different and was not tested. Definition: `docs/copy_assignment_definition.md`. Rules, bars and predictions are unchanged.

Written before any count of this row. Workspace `r1775_eichler_baseline`, script `bench/o2_truth/neutral_truth.py`.

## Why
Advisor, 2026-11-17: work on multi-copy families needs a baseline; proposed baseline = the adapted Eichler filter (a read whose
best alignment leads every rival by more than 10 in alignment score is assigned to that copy); Winnowmap, VACmap and
DuploMap named as methods that may help. Earlier rows measured the rule on our own catalogue and panels (r1270, r1403,
r1404, r1477, r1509); this row measures it against read origin with a truth that depends on no catalogue and no arm.

## Step 0 (feasibility, done before the table)
Winnowmap 2.03 (commit 72403b3) on the KB3781 reads with a primary on chrX, `-ax splice:hq` with the flags of our BAMs
and a `meryl` k = 15 list (`distinct=0.9998`). GO when >= 0.95 of the reads map and the spliced share of primaries is
within 0.05 of minimap2's. Read from the source before the run: the `splice` presets set `SVaware = false`, so on
spliced reads Winnowmap runs its 2020 method (weighted minimizers) and not the Winnowmap2 method.
VACmap and DuploMap are not run: both take DNA reads only (no spliced mode; DuploMap also needs a human SD database).

## Sample and reads
Gorilla KB3781 (fibroblast; the T2T individual, haplotypes `o2contested/haps/KB3781.{mat,pat}.mmi`). Reads = molecules
with a PRIMARY record on the exons and strand of a member of a truth family (Compara, simiiformes, as r1772). They are
realigned to the reference with minimap2 (our flags) so that every record of each read is in hand.
- Placements of a read = its records merged into loci (same chromosome and strand, overlapping spans).
- lead = best AS minus the best AS at any other locus. Strata: S0 no other locus; S1 lead > 10; S2 lead <= 10.
- Hash-ordered caps: S0 5,000 reads, S1 30,000, S2 30,000.

## Truth (origin), independent of catalogue and arms
- The read is aligned to both haplotypes. It is usable when its best hit is near-perfect (NM <= max(2, 0.003 x aligned));
  its origin hits are the hits with that NM (slack 0) and >= 0.98 of the best hit's aligned bases.
- Each placement kept for a read (the primary, the best-AS locus, each arm's pick, and loci with AS >= 0.8 x best, at
  most 12) is turned into the reference sequence under its aligned blocks and aligned to both haplotypes; its haplotype
  loci are its hits within NM + 2 of its best on that haplotype, the best being near-perfect.
- A placement is the origin when one of its haplotype loci overlaps an origin hit. A read is decidable when exactly one
  placement is; "tie" when several are; "outside" when none is.

## Arms (each picks one placement or abstains)
- primary: minimap2's primary.
- eichler-keep: the primary when it is the best-AS locus and lead > 10; otherwise abstain.
- eichler-move: the best-AS locus when lead > 10, wherever the primary is (the advisor's wording).
- winnowmap: Winnowmap's primary (only if step 0 is GO).
Our own assignment (assign stage) is NOT an arm of this row: it needs a current-default families + assign run on
KB3781 and is scored by the same script in the next row.

## Readings, fixed now
1. Baseline quality: eichler-keep's wrong share among the decidable reads it decides. Called sound at <= 1%.
2. The advisor's wording: on decidable reads where eichler-move picks a locus other than the primary, right and wrong of
   each. eichler-move is refuted there if the primary is right more often.
3. What the baseline leaves: decidable S2 reads (the rule abstains on all of them): the primary's right share is the
   figure any method for multi-copy families has to beat.
4. Winnowmap gains when its right share on decidable S2 reads exceeds the primary's by >= 1 percentage point with
   >= 200 decidable reads; otherwise "no gain as a spliced aligner".

## Predictions
1: sound (< 1% wrong). 2: refuted, the primary right on > 90% of those reads. 3: primary right on 85-97% of decidable
S2 reads, and most S2 reads are ties or outside rather than decidable. 4: no gain.

## Limits stated beforehand
- One sample, one tissue; reads whose primary lies outside the truth families are not sampled.
- The truth compares mismatch counts; a read from an allele 1-2 edits from two copies is a tie by construction.
- minimap2 2.30 here against 2.31 for the stored BAM.
- A placement that is the origin under this truth can still be a copy the reference places elsewhere on the haplotype;
  no reciprocity check beyond the near-perfect rule (the known weakness of r1659).
