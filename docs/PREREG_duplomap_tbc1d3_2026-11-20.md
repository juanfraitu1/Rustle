# PREREG: DuploMap on Iso-Seq reads at TBC1D3 (row 1844, 2026-11-20)

User request: "lets do the demo of duplomap in TBC1D3, if I don't at least try it my advisor wont believe me". A demonstration
on one family, not a test with power: the truth has about 30 decidable near ties in TBC1D3 (row 1790).

Written before any DuploMap alignment was produced. Known at this point (setup only): DuploMap 0.9.5 (commit 32e2385) builds
with minimap2 v2.17; its `-x` preset takes only `pacbio` / `nanopore` (no spliced preset); `duplomap-prepare` on the 144 CHM13
SEDEF pairs with both ends on a TBC1D3 gene gives 4 clusters, and with its default (`--skip-copy-num 6 3000`) it REMOVES 11
regions of copy number >= 6, which are the TBC1D3 gene units themselves.

## Setup
- Sample: HG002 Kinnex on CHM13 (`reads/HSA/HG002_kinnex.t2t.bam`, minimap2 `-ax splice:hq -uf --eqx -Y -N 50 -p 0.1`).
  Used before (rows 1776-1791, 1829): reused, said.
- Input reads: every record of the full BAM that overlaps a database region or a removed high-copy region, primaries and
  supplementaries only (`-F 0x100`), so DuploMap learns its PSV reliabilities from all reads of the locus.
- Arms: A = database as built by default (high-copy regions removed); B = database with `--skip-copy-num 1000 3000` (TBC1D3
  units kept). Both with `duplomap -x pacbio`, other options default.
- Truth: `bench/o2_truth/neutral_truth.py` tables of the 25% read sample, repaired test first
  (`r1785_truth_repair/hg/reads.tsv`), old test beside it (`r1784_winnowmap2_spliced/hg/reads.tsv`). A DuploMap pick = the
  candidate locus of the read's table row that its realigned primary overlaps on the same strand; none = "elsewhere"
  (counted wrong when the read is decidable).
- Scored reads: reads of the table whose minimap2 primary overlaps a TBC1D3 gene (`TBC1D3`, `TBC1D3B`-`K`, `TBC1D3P*`).

## Step 0: does it run on spliced reads
Reported, per arm: reads in, reads in `realigned.bam`, reads DuploMap realigned (tag `re`) vs left as they were; among the
realigned, the share whose new alignment keeps >= 90% of the read's aligned bases of the input alignment, and the share
that still has an intron (`N`) or a deletion >= 50 bases where the input had an intron.
- APPLICABLE: >= 80% of the realigned reads keep >= 90% of their aligned bases.
- NOT APPLICABLE AS BUILT: fewer. Placement is still scored (a clipped alignment can sit at the right copy).

## Step 1: placement
On the decidable scored reads (all strata, and near ties = S2 alone), per arm: minimap2 primary right / wrong, DuploMap
right / wrong, fixes and breaks; DuploMap calls with MAPQ >= 30 and their right / wrong.
- GAIN: DuploMap right - primary right >= 3 reads on the decidable scored reads of the repaired table, and <= 1 read
  broken in the clear-lead stratum (S1).
- NO GAIN: otherwise. Fewer than 15 decidable scored reads = UNINFORMATIVE.
Also reported, no bar: PSVs in the clusters and how many DuploMap marks reliable (`rel` >= 0.95) from RNA reads.

## Predictions (mine)
Arm A leaves nearly all TBC1D3 reads untouched (their loci are removed from the database). Arm B realigns them without
introns: most alignments break at the first intron (NOT APPLICABLE AS BUILT), and placement does not beat the primary.

## What this cannot show
One family, about 30 decidable reads, one individual; the truth is the catalogue-free test with its known limits. A NO GAIN
here is evidence about DuploMap as released on spliced reads, not about PSV-based placement in general (rows 1783, 1829).
