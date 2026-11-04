# Goal: family-wise assembly (drafted 2026-11-04 from registers 1469-1535)

## The goal
Emit one correct model per expressed copy of every multi-copy family on the r1521 panel (A119b, 324 families, 1,092 expressed
members under definition T), judged by family_panel.py as registered: complete families >= 209, fused <= 49, ordinary-member
hit rate >= 0.829, T-only hit rate >= 0.678, with the catalog pair score (r1530 scorer) not below precision 0.802 / recall 0.740.
A method counts only if it is registered before its numbers, runs the panel first under a wall limit, and never borrows
strength across copies (a copy's model is built from reads threaded to that copy alone).

## What we know about the misses (212 expressed members without a model)
1. Multimapper misses are at most 35 and the panel ceiling for any multimapper method is <= 10 (r1532): 18 are exon-identical to
   a sibling, so no transcript-level graph separates them; read-level rules (r1526, r1527), gate/transcript edges (r1528,
   r1529), pooled depth (r1531) and graph threading (r1532, r1533) are all closed with recorded verdicts.
2. 177 misses have uniquely placed primaries (r1534). They split: 51 single-exon truth genes whose own unspliced primaries (up to
   1,787 at CTAGE15, 4% internally primed) form no model in the default pass 1; 53 more whose reads sit unspliced inside the
   gene span (46 of 53); 52 whose own reads are spliced yet no spliced model appears (SPANXA2: 53 reads on one intron chain,
   MAPQ 13-22, zero models even in a clean primaries-only window); 20 unspliced-read members (the r1518/r1520 trade); 1 neighbour.
3. The single-exon block is reached today only by the r1489 `--mono-model cluster` layer (final.gtf), which the panel base and
   the families stage do not consume (SPLIT 2026-11-02). r1535 arm L1 measures it genome-wide; adopting it is a product decision.
4. The panel ruler is part of the picture: it counts any model, and one 113 kb TAF11L single-exon blob carries 30 of the 87
   T-only hits (r1535 D1: a span guard failed its own adoption bar at 0.872). Window-only gains are secondary-seeded blobs
   (r1534: vanish with secondaries stripped); never read a window gain as a lever without a primaries-only rerun.
5. Sequence does not help where structure is identical: HiFi duplicates carry no signal (r1525), IsoCon recovers 9.6% with 58%
   artifacts (r1515), junction snapping and compatible-chain admission fuse families (r1526, r1527).

## What the next rows are
- r1536: why pass 1 drops a 53-read one-chain locus (SPANXA2) and whether the same rule explains >= 10 of the 52 spliced misses.
- r1535 L1 result: does the mono layer reach >= 26 of the 51 single-exon members without raising fused.
- Decision (user): which GTF the panel and the families stage judge, x.gtf or final.gtf.
Everything else on the family layer is closed on this panel unless a new class of evidence appears.
