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

## RESULT (jobs 4435658 projection 613 s + 4435682 SQANTI; split_arms crashed once on exon lines, fixed before any class was read)
| arm | judged targets | bar 1 reach | bar 4 motifs | base FSM+ISM / artifact | own tier (n) | tied tier (n) |
|---|---|---|---|---|---|---|
| H hotspot out | 286 / 892 | 32/78 = 0.41 FAIL | 0.892 FAIL | 11.8% / 37.0% (799) | 41.3% / 17.5% (63) PASS bar 2 | 6.0% / 13.1% (84) FAIL bar 3 (FSM+ISM < 0.9 x base) |
| G T >= 2 | 892 / 892 | 0.726 PASS | 0.365 FAIL | 1.5% / 91.7% (6,532) | 1.0% / 98.0% (2,695) FAIL | 0.0% / 99.5% (15,647) FAIL |
None tier (descriptive): H 1,004 models, FSM+ISM 3.5%, artifact 35.3%; G 112,850, 99.4% artifact.
Verdict: FAIL as registered, both arms. The G gate is vacuous on this panel: every target has T >= 2 because unspliced reads and
AS-tied secondaries pile up on tandem copies by the thousand (hotspot copies with 2 assigned reads show T 2,600-4,200), so arm G
is r1519 again. Arm H is a real judged set and shows the honest shape of the method: the own tier is the only clean product
(41% FSM+ISM vs base 12%, artifact 17.5% vs 37%) and it is 63 models over 286 targets; the none tier (1,004 models, 35% artifact)
is annotation transfer that nothing in this sample can verify.
Decision: projection stays a bench tool for per-member own/tied support; the family layer moves to tied-seed targets + fused-model
cuts (read-backed), registered on its own. The G gate would need spliced placements only to mean anything; not pursued here.
