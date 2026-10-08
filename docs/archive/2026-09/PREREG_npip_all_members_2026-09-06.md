# PREREG — every NPIP member accepted: two definition polishes and the read-through certificate (2026-09-06)

O1 evaluation from now on = `bench/o1_eval.py`: sensitivity (truth loci rediscovered), specificity (family members
that are truth), and the bipartite 1:1 coverage of each truth member by its rediscovered unit (§6ft). NPIP truth
= the 26 LCR16a loci (`adj/size/lcr16a.bed`); Soto truth = `bench/soto/80_fams.chr.bed`.
Baseline (`rna_units_v11`): NPIP sensitivity 24/26, specificity 24/24 members, coverage median 0.89 (units) /
1.00 (extents); Soto 272/362 = 0.751, specificity 0.703 (no member_status on that catalog), coverage 0.98.

## Polish 1 — an unexpressed member keeps its annotated model as its unit
The member at NC_073242.2:29,864,793 (full 22.5-kb core, `kept_full`, one primary read with no block in its
exons) gets no unit because the unit rule requires ≥ 1 read inside the emitted chain. Change: a member whose
chain is the GFF fallback is emitted with `n_reads 0` (`source gff_fallback`); `--no-units-keep-unexpressed`
restores the previous row set. O2 must accept a copy with 0 reads (checked before the run).
Prediction P1: NPIP sensitivity 25/26; genome-wide the count of new 0-read units is reported; O2 outputs on
the paired 35 unchanged except where a 0-read unit becomes a candidate (reported).

## Polish 2 — the majority counts the locus itself
The core rule keeps a member whose core is shared with ≥ half of the OTHER members ((n − 1)/2). The 7-kb
locus at 28,300,719 carries 22 % of the core shared with 15 of 31 others (15 < 15.5) and is dropped. Change:
"shared with at least half of the family, itself included" (depth + 1 ≥ n/2) — `--core-majority-inclusive`.
Predictions P2: NPIP sensitivity 26/26 with polish 1; Soto sensitivity ≥ 0.751 and specificity within 0.01 of
0.703; genome-wide the dropped → kept changes are counted and their read support reported (≥ 3 reads vs < 3).
Fail rule: Soto specificity falls by more than 0.01 ⟹ not adopted as default, reported.

## O2 — the read-through certificate (the 71 of §6fp)
A read whose alignment runs past the candidate's locus is rejected for its unaligned bases. Two classes:
(a) 48 real read-throughs: the unaligned tail is ANOTHER catalog unit's expressed sequence; (b) 15 artefacts:
a giant unsupported intron the aligner chained. Change: a read position the candidate leaves unaligned is
"explained" if any other candidate of the family aligns it (the read-star already knows every candidate's base
at every read position) or if it lies beyond a giant (> 50 kb) intron fewer than 3 molecules support (the O1
mis-chain rule applied to the molecule); explained positions leave the certificate's numerator and
denominator. The pairwise certificate then runs on the shared columns as before; a molecule explained by two
loci that covers different parts of the read is assigned to the candidate covering the larger part and marked
`readthrough_into = <other>` (a new column), never a K = 0 tie. Escape `--no-readthrough-certificate`.
Predictions P3: of NPIP's 71, ≥ 48 become assigned (the real read-throughs to unit 2 → MCL27:0 among them);
the 62 anchors stay 0 wrong; paired-35 agreement stays 1.0000; the simulation's precision stays ≥ 0.998.

## Status (2026-09-06 18:40)
P1 held (NPIP 26/26 with polish 1; O2 accepts 0-read copies). P2 failed (row 725; polish 2 OFF). P3 partial:
of the 72 reads past the extent 38 pass by a partner (the family's own candidates could not be partners for a
read-through into another family — partner rows added to the sweep), the mis-chain cut explained 0 of the 15
artefacts (open); anchors 0 wrong. Paired 35 under everything: running (`sweep_v17`).
