# Row 1679 — a dating feature: divergence on aligned, non-exon, non-repeat sequence (2026-11-13)

Goal: "find empirical rules to improve edge recall in the multi-copy gene families definition". Families layer (O1),
opt-in only if found, not for `advisor_share`. Rows 1674, 1675 and 1678 closed models on the 20 alignment features with
the same reading each time: they find relatives but cannot date the duplication. Row 1678's registration said I would
not propose another model on those features; this row adds features of a different kind, and says so.

**Mechanism tested.** A young duplicate keeps alignable intron and flank sequence; an older paralog keeps exons only.
`nonexon_frac` and `span_cov` already measure aligned non-exon sequence, but interspersed repeats (Alu, L1) align
between unrelated loci and pollute both. The genomes are soft-masked, so the repeat part can be removed.

**New features (frozen, `edgefeat.py --fasta GENOME`):** over the records of the pair's locus pair, positions aligned
gap-free that are outside the exons of both loci and upper-case (unmasked) in both: `neutral_bp` (their number, log1p
in the model), `neutral_id` (share of them that match; 0 when fewer than 50), `neutral_frac` (neutral_bp / aln_bp).

**Design = row 1675's, unchanged.** L2 logistic on the 20 + 3 features, positives = simiiformes truth pairs, fit on
A119b + HG002 pooled; threshold = lowest score with 5-fold out-of-fold precision >= 0.75 and >= 20 true pairs.
**Go/no-go:** no such threshold -> STOP, gorilla not scored.
**Test (frozen; comparable tables: rebuilt OR6737 + KB3781):** pooled >= 20 new true pairs at precision >= 0.60.
**If it passes:** row 1674's step 2 (family level through `--extra-edges`, same bars) decides whether a rule is found.
**Reported, no bar:** the three features' medians for truth pairs, older relatives and unrelated pairs.
**If it stops or fails:** the row closes; the dating idea is then closed for alignment-derived features.

## Result (2026-11-13, job 4470703, workspace `r1679_neutral`) — STOP

No out-of-fold threshold reaches 0.75 with 20 true. Out of fold: 20 true at 0.476 (22 false), 40 at 0.408, 80 at 0.357.
Row 1675 without the three features had 0.606 / 0.421 / 0.352, so the features add nothing. Gorilla not scored.

Reported medians (refused pairs; share of pairs with >= 50 usable positions):

| sample | truth pairs | older relatives | unrelated |
|---|---|---|---|
| A119b | 159: 712 bp, id 0.855, 0.673 | 3,311: 0 bp, 0.215 | 199,183: 0 bp, 0.012 |
| HG002 | 118: 1,526 bp, id 0.904, 0.797 | 2,727: 0 bp, 0.220 | 164,579: 0 bp, 0.082 |
| OR6737 rebuilt | 123: 793 bp, id 0.887, 0.756 | 3,708: 1 bp, 0.218 | 254,265: 0 bp, 0.127 |
| KB3781 | 99: 768 bp, id 0.885, 0.737 | 3,268: 0 bp, 0.193 | 223,472: 0 bp, 0.108 |

Reading: the signal is real in the medians (most truth pairs keep non-repeat intron or flank alignment, most others
keep none), but the classes are too unequal in size. About a fifth of ~3,000 older relatives and 1-13% of ~200,000
unrelated pairs also keep such sequence, which outnumbers the ~100-160 truth pairs. The row closes; the dating idea is
closed for alignment-derived features.
