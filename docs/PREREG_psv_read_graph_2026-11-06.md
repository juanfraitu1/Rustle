# PREREG: reads as nodes, PSV agreement as edges (row 1555, 2026-11-06 22:50, user's advisor's suggestion)

## Question
Row 1553's overlap graph fused because whole-sequence identity is dominated by PSV-free stretches. If edges count
only the informative sites, do components separate MCL1's copies?

## Set-up (MCL1; row 1544 ctl_genomic.paf, --eqx cigars; backbone = window g0)
Each read's best alignment on g0 gives ref/non-ref calls at every aligned base (X vs =; binary alleles).
Informative site = depth >= 20 and minor-allele fraction >= 5%. Edge between two reads iff they share >= 3
informative sites and agree at every shared site. Union-find components.
## Bars (as row 1553)
- Purity: no component holds baseline-assigned reads of two copies (<= 1 violation).
- No fusion: largest component <= 20% of the reads that reached g0.
- Coverage: components containing the baseline-assigned reads hold >= 50% of those reads.
- Time is not a bar (seconds; the PSV columns come from alignments already paid for): answer lever only.
## Prediction
Purity passes (strict all-agree rule). Fusion vs fragmentation is again the coin flip; the expectation is
fragmentation (1-2 kb reads cover few sites; binary alleles merge copies that differ from the backbone with
different bases), so coverage FAILS, and a component graph cannot replace per-read decisions. If instead the
components are both pure and large, that is new structure and the next row would be the KB3781 truth test.
