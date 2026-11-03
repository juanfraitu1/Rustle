# PREREG r1522: projected models v2, judged set gated (hotspot excluded / transcribed target)

Written before any arm number. Parent: r1519 (docs/PREREG_project_models_2026-11-03.md), FAIL as registered because 99% of the
projections fell in the chr13 13-16 Mb tandem hotspot, where CAT annotates nothing and intron motifs are mutated. The mechanism,
templates, thresholds, tiers and wall bar are unchanged; what changes is which targets are judged.

## Arms (one run; two judged sets cut from the same projected.gtf and base.gtf)
- H: targets outside chr13:13,000,000-16,000,000 (the hotspot is reported, not judged).
- G: targets with T >= 2, where T = records overlapping the target locus (locus_start..locus_end, any strand, spliced or not)
  that are primary, or secondary with AS >= 0.98 x the molecule's genome-wide best (r1521's T). The hotspot stays in.
Base control per arm = the arm's own transcripts at the same judged targets.

## Bars (r1519's, per arm)
1. Reach: >= 50% of the judged low-read copies (<= 2 reads) gain >= 1 projected chain of any tier.
2. `own` tier: artifact share <= base's and FSM+ISM >= 0.9 x base's (SQANTI3 vs CAT chr13+chr20). Tier with < 20 models: reported, not judged.
3. `tied` tier: same two conditions, same n rule.
4. Motifs: canonical (GT-AG, GC-AG, AT-AC) at >= 95% of projected introns in the judged set, all tiers.
5. `none` tier: descriptive.
Wall <= 30 min for the projection (SQANTI excluded, as in r1519).
PASS = an arm passes bars 1, 2, 4 (and 3 when judged). PASS -> genome-wide run with that gate, then the Rust port as a layer.
FAIL -> prototype stays on the bench; the next family-layer attempt is the tied-seed targets + fused-model cuts, registered on its own.
Not in this row: tied-seed extra targets, fused-model cuts, Winnowmap.
