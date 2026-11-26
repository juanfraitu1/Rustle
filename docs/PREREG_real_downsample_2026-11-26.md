# PREREG: do near ties recover real transcripts at low depth? Real reads, downsampled (row 1964, 2026-11-26)

The one positive result for the near-tie pool is on planted reads: at 3 reads per isoform it gives two to three times
the exact isoforms of a primaries-only assembly, at 60 reads nothing correct (rows 1952, 1963, `two_products_2026-11-25`).
Planted sets hold only multi-copy families (42% of primaries tie-broken; a real sample 0.26%). This row asks the same
question on real reads, with a truth that needs neither annotation of expression nor planting.

## Design (fixed here)
- Samples: gorilla OR6737 (testis) and KB3781 (fibroblast). Each read set is downsampled to 5%, 10% and 25% of its reads
  (`samtools view --subsample F --subsample-seed 1964`: a read keeps all of its alignments).
- Arms at each depth: `primaries` (`assemble --seed-pool primary`), `pool` (`--seed-pool good`: secondaries at AS >= 0.98
  of the best), `mode` (pool + families + stage `multicopy`). Settings of the full-depth runs.
- Truth: a spliced intron chain (contig, strand, introns) is CONFIRMED when the FULL-DEPTH PRIMARIES-ONLY assembly of the
  same sample has it (`r1920_mono_from_primaries/real_or/p.gtf`; KB3781's is made in this job). No secondary alignment
  takes part in the truth, so an arm with secondaries cannot confirm itself. Single-exon transcripts are not scored.
- Per depth: `added` = chains of the arm that `primaries` at that depth lacks; `lost` = the reverse.
- Scorer `bench/o2_families/depth_score.py`; workspace `r1964_real_downsample`.

## Claims (for `pool`; `mode` is reported with the same measures, no bar)
- **D1, the added chains are real:** at 5% and at 10%, at least 50% of the added chains are confirmed, in both samples.
- **D2, net gain:** at 5% and 10%, confirmed added chains are at least 3 times the confirmed lost chains, in both samples.
- **D3, a depth effect:** confirmed added chains, as a share of the confirmed chains of `primaries` at that depth, are
  larger at 5% than at 25%, in both samples.
- Descriptive: the size of the gain (that share), the precision of each arm, and the share of added chains that are
  exact RefSeq chains.

## Predictions (mine)
D1 fails narrowly (40-50% confirmed; the planted tied layers were 0.45-0.72 at low depth). D2 holds. D3 holds. The gain
is small genome-wide (under 2% of the confirmed chains at 5%), because near ties are rare outside families.

## Reading
D1-D3 hold: "at low coverage the near-tie pool recovers transcripts that primaries alone miss" may be said for real
reads, with the size. D3 holds and D1 fails: the pool adds more at low depth, and under half of what it adds is
confirmed. D3 fails: the planted depth effect does not show on real reads; near ties are then only "safe when labelled".
One seed per depth; the truth favours chains with primary support, so a chain that only tied reads can show at any
depth counts as unconfirmed.

## Result (job 4509341, 46 min; workspace `r1964_real_downsample/result.txt`)
Confirmed = in the full-depth primaries-only assembly. `primaries` itself is 0.90-0.92 confirmed at every depth.

| sample, depth | primaries: chains (confirmed) | pool adds (confirmed, share) | pool loses (confirmed) | gain, % of confirmed | mode adds (confirmed, share) | mode gain, % |
|---|---|---|---|---|---|---|
| OR6737 5% | 13,454 (12,401) | 108 (56, 0.52) | 2 (2) | 0.45 | 1,030 (344, 0.33) | 2.77 |
| OR6737 10% | 20,659 (18,630) | 209 (91, 0.44) | 3 (2) | 0.49 | 1,523 (385, 0.25) | 2.07 |
| OR6737 25% | 34,789 (31,397) | 411 (121, 0.29) | 14 (11) | 0.39 | 2,089 (323, 0.16) | 1.03 |
| KB3781 5% | 17,923 (16,429) | 89 (45, 0.51) | 2 (1) | 0.27 | 530 (191, 0.36) | 1.16 |
| KB3781 10% | 24,922 (22,413) | 135 (58, 0.43) | 7 (4) | 0.26 | 685 (174, 0.25) | 0.78 |
| KB3781 25% | 36,890 (33,211) | 271 (82, 0.30) | 16 (8) | 0.25 | 916 (156, 0.17) | 0.47 |
- **D1 FAILS:** 0.52 and 0.51 at 5%, but 0.44 and 0.43 at 10% (bar 0.50 at both depths).
- **D2 HOLDS:** 56 v 2, 91 v 2, 45 v 1, 58 v 4.
- **D3 HOLDS by the rule and is negligible in size for the pool:** 0.45% v 0.39% and 0.27% v 0.25%. For the mode the
  depth effect is clear: 2.8% -> 1.0% and 1.2% -> 0.5%.
- Predictions: D1 "fails narrowly" right, D2 and D3 right, the size (under 2%) right for the pool; the mode at 5% on
  OR6737 is above it (2.8%).
- The confirmed share of what is added falls with depth in every arm (pool 0.52 -> 0.29, mode 0.33 -> 0.16): the deeper
  the sample, the more of the addition is chains its own full-depth primaries never show.

## Reading
On real reads the near-tie pool adds little at any depth (0.3-0.5% more confirmed chains) and about half of what it
adds at 5% is confirmed, less at 10% and 25%. The planted result (two to three times the isoforms at 3 reads) does not
carry over as a genome-wide gain: it is a property of multi-copy families, which are a small part of a real sample. The
mode shows the depth effect (gain 2.8% at 5% against 1.0% at 25% on OR6737), with a third to a sixth of its additions
confirmed. What may be said: near ties are safe when labelled and layered; their gain is confined to families and
grows as depth falls. Not measured here: the same score restricted to family loci (the truth for those is weaker,
since full-depth primaries at identical copies are placed by the tie-break).
