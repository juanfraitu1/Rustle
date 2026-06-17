# Integrated end-to-end: identify families → assign reads to copies (hard multimappers)

Runs the full pipeline and answers, honestly: can the method identify multi-copy families accurately
AND assign reads to copies — especially hard multimappers (MAPQ-0)? Validated WITH ground truth on the
synthetic 5-copy benchmark; censused on real GGO (where there is no per-read truth).

## End-to-end on the 5-copy benchmark (WITH truth) — bench/integrate_end_to_end.py
For each divergence level K: (1) POA variation graph of the 5 copies → graph core-score (family
detected?); (2) PSVs recovered from the graph's variant columns; (3) the MAPQ-0 *hard multimappers*
assigned by their PSV allele-vector → accuracy vs truth.

| K | family score | detected | PSVs from graph | hard MM (MAPQ-0) | PSV acc on hard MM |
|---|---|---|---|---|---|
| 0 (identical) | 1.00 | ✓ | 0 | 200 | 0 — UNASSIGNABLE (no info) |
| 1 | 1.00 | ✓ | 1 | 80 | 0.50 (1 column < 5 copies) |
| 2 | 1.00 | ✓ | 2 | 0 | — (minimap2 resolved) |
| 4 / 8 | 1.00 | ✓ | 2* | 0 | — (minimap2 resolved) |

- **Family detection works at every K** — the 5 near-identical copies are correctly called one family
  (score 1.0 ≫ T), and the PSVs are recovered from the graph (the true # of *varying* columns; the
  base-4 design saturates at ⌈log₄5⌉=2 distinguishing columns, so "K=4/8" add no new variation — the
  graph correctly reports 2).
- **Hard multimappers**: at K=0 (identical) they are information-theoretically unassignable (both the
  aligner and PSV); at K=1 PSV resolves 50% of the MAPQ-0 reads the aligner left tied; at K≥2 there are
  **no** MAPQ-0 reads because **minimap2 itself resolves copies once ≥2 PSVs are spanned by full-length
  reads**.

**The honest catch:** for full-length HiFi reads the aligner already resolves non-marginal cases, so
PSV-assignment's *distinctive* edge over minimap2 on hard multimappers is the **marginal (single-PSV)
regime** — plus it provides a principled, identifiability-bounded assignment and an explicit
"unassignable" verdict where the aligner just emits MAPQ-0.

## Real GGO census — how much hard-multimapper signal exists
Over the graph-defined multi-copy families (sampled 400; 174,459 primary reads at family loci):
- **hard multimappers (MAPQ-0) = 0.7%** of family reads (1,156); **6% of families** have ≥5; a few
  carry most (FAM156: 358).

So the hard-multimapper regime the method most distinctively addresses is **sparse in GGO** — consistent
with every prior finding (most paralogs are coordinate-separated; collapsed/co-located real copies are
rare). The 26 hard families are where it applies; there is no per-read truth on real data, so real-GGO
is a *demonstration* + census, and the accuracy proof comes from the synthetic benchmark.

## Bottom line (answering the question)
- **Identify multi-copy families accurately: yes** — the per-family POA-graph definition (validated,
  fixes over-merge) detects them, including on real data.
- **Assign reads to copies, incl. hard multimappers: yes, up to the identifiability boundary** —
  validated with truth on the synthetic hard case; it resolves the marginal hard multimappers PSVs can
  resolve and provably declares the rest unassignable.
- **Caveat:** real GGO under-exercises the hard regime (0.7%), and for full-length reads the aligner
  already handles ≥2-PSV cases — so the method's biggest distinctive payoff is in a **deep co-located
  dataset** (testis HiFi for DAZ/RBMY) where near-identical copies and hard multimappers are abundant.

## Reproduce
- `MINIFORGE python bench/integrate_end_to_end.py` ; `python3 bench/integrate_fig.py`
