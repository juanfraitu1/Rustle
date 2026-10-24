# A locus model for overlapping genes on ideal reads: primary-first representatives and structural separators (2026-10-06)

**Status (arms produced 2026-10-07 00:04-00:08 PDT; independently verified, workflow wf_994ece76-8d2: every number reproduces): DEV, human CHM13 v2.0 / CAT-Liftoff, simulation, circular by construction, the same two windows (NPIP chr16, TBC1D3 chr17) x two read replicates as `docs/IDEAL_EXPRESSION_DEFAULT_2026-10-06.md`. Protocol: `docs/PREREG_locus_units_2026-10-06.md` (committed 4ce90b18; Amendment 1 bff7dbfe, written before any arm existed). Tools: `bench/entangled/{locus_units.py (12 unit tests), run.sh units2}`. Products: `/mnt/linuxdisk/tmp/entangled_2026-10-06/{NPIP,TBC1D3}/rep{1,2}/{q,sd,sq}.*`.**
User direction: option B (a locus-level model for junction-sharing overlapping genes; the advisor wants clean combinatorial structure, provable statements, no arbitrary thresholds).

## Answer

1. **Rule 2 (evidence-free structural separators) is refuted at the family level on the dev windows (its wording was chosen with a truth-chain oracle on the same windows).** The private-bridge rule removes the transcripts that are the only junction-link between groups: on NPIP replicate 1 its 30 separators include almost only a gene's own long or main isoform (NPIPB3, B4, B5, B12 and B13 with 47 reads, SNX29, MAPK3, MRTFB, SGF29, TUFM, ITGAL, ZNF785, MLKL, AATF, SRCIN1 ...), in TBC1D3 TBC1D3G, TBC1D3B and TBC1D3P5. A bridging isoform and a fusion have the same junction structure; (hypothesis: reads' PAS and promoter signals, F1v2's rule, could tell them apart; the evidence-based F1 arm of `docs/LOCUS_UNITS_LEVELS_2026-10-07.md` found no signal to act on in these windows). Pooled over the four runs S_D finds 38 of 54 reachable copies (D: 48) and 9 of 24 entangled (D: 12).
2. **Rule 1 (primary-first representative) repairs exactly the failures it targets and is not a candidate by the registered rule.** E3 (the representative has one of the copy's chains) goes from 50 to **54 of 54** reachable gene-runs: NPIPA2 (both replicates) and TBC1D3P4 / TBC1D3P3 are fixed. The family-level loss is one event: in NPIP replicate 2 the NPIP cluster splits and NPIPB6, NPIPB8 and NPIPB9 land in a sibling cluster (E4 fails for three copies; Compara F .733 -> .643). R found 49 against D's 48.
3. **Primaries-only seeding (arm P, not registered as a family-level candidate here) makes the registered ideal verdict YES for both families**: NPIP R found 14 / 13 (bar 13), TBC1D3 13 / 13 (bar 12); the price is the entangled copies (E found 7 of 24 against 12): their cluster membership came from the enlarged bodies that secondary-borne chains build.
4. **The family-level criterion E4 ('the holder is in the one cluster K*') is brittle to partition granularity**: each arm that changes a few representatives or bodies re-partitions NPIP and moves a subfamily to a sibling cluster (rep 2: MCL0 -> MCL2). This is (post hoc) the advisor's hierarchy question (nested levels, level-relative TP / FP) showing up in the measurement.

## Pooled results (four runs; R = 54 reachable gene-runs, E = 24 entangled)

| arm | R found | E found | R: E2 / E3 / E4 | E: E2 / E3 / E4 | NPIP R found r1 / r2 (bar 13) | TBC1D3 R found r1 / r2 (bar 12) | Compara F NPIP r1 / r2 | Compara F TBC1D3 r1 / r2 |
|---|---|---|---|---|---|---|---|---|
| D default | 48 | 12 | 50 / 50 / 54 | 16 / 18 / 22 | 12 / 12 | 12 / 12 | .733 / .733 | .753 / .737 |
| P primaries only | 53 | 7 | 53 / 53 / 54 | 17 / 17 / 14 | 14 / 13 | 13 / 13 | .746 / .724 | .753 / .737 |
| C sub-chain drop | 48 | 13 | 50 / 50 / 54 | 16 / 18 / 23 | 12 / 12 | 12 / 12 | | |
| PC | 53 | 7 | as P | as P | 14 / 13 | 13 / 13 | .746 / .724 | .753 / .737 |
| **Q** Rule 1 | 49 | 12 | 52 / **54** / 51 | 16 / 18 / 22 | 13 / **10** | 13 / 13 | .733 / **.643** | .753 / .737 |
| **S_D** Rule 2 | 38 | 9 | 51 / 47 / 46 | 17 / 15 / 19 | 11 / 9 | 9 / 9 | .690 / .690 | .753 / .737 |
| **S_Q** Rules 1+2 | 45 | 10 | 52 / 54 / 47 | 16 / 16 / 20 | 13 / 10 | 11 / 11 | .690 / .667 | .753 / .737 |

Registered rule (YES at >= ceil(0.9 R), lower of the replicates): D NPIP NO, TBC1D3 YES; P and PC YES, YES; Q NPIP NO (rep 2: 10), TBC1D3 YES; S_D NO, NO; S_Q NO, NO.

## Predictions and the decision rule (prereg sections 5 and 6)

Q1 (Q finds >= 52 of 54 reachable): **failed** (49). Q2 (Q finds >= 11 entangled): held (12). Q3 (Q gives YES for both families in both replicates): **failed** (NPIP replicate 2). Q4 (S_D finds >= 2 more entangled than D with no loss at R): **failed** (9 against 12; R -10). Q5 (S_Q at least Q at R and E, and one more at E): **failed**. Q6 (Compara F of Q and S_Q >= D minus .005): **failed** (NPIP replicate 2; S_Q NPIP both).
Decision rule: Q has R + E 61 > 60 and R >= D and E >= D - 1, but Compara F fails in NPIP replicate 2 and the direction does not hold in both replicates of both families: **no arm is a CANDIDATE**; no held-out or real-data check is earned.

## What the arms show (post hoc, mechanism)

- Rule 1 per copy against D: gains NPIPA2 (representative purity .454 -> 1.0, both replicates), TBC1D3P4 / P3 (representative displaced chain -> own chain); no entangled copy changes; losses only the three NPIP replicate-2 cluster moves (NPIPB6, NPIPB8, NPIPB9: MCL0 -> MCL2), locus purities unchanged.
- Rule 2 cuts 30-33 transcripts per NPIP replicate and 11 per TBC1D3 replicate (F1v2 cuts 1-2); the truth-chain oracle had predicted the false-cut cost (16 of 108 junction-sharing entangled genes resolved, 12 more non-entangled genes split than with components). At the family level the cuts also shrink bodies and split loci into several nodes: copies move out of the family cluster (TBC1D29P and TBC1D3B to a sibling cluster, NPIPB5 to MCL5).
- The entangled copies that fail in D split into (a) artifact-driven (NPIPA8, NPIPB3: locus purity .129 and .234, 5 and 24 artifact transcripts in the locus; NPIPA1: the representative is an artifact chain) and (b) truly fused (NPIPB14P with the readthrough PDXDC2P-NPIPB14P as representative, AC138894.1 with CLN3 and NPIPB7, PKD1P6 with NPIPP1); primaries-only seeding fixes NPIPB3's locus purity and NPIPA1's representative, but breaks NPIPA8 (representative purity .906 -> .295) and loses the cluster membership of the two AC138969.1 copies and NPIPA1.

## What this does not show, and the next registered step

Same two windows; the arms were motivated by what the dev windows showed; the primary-supported marks come from a second assembly. Nothing here is held out. What the data point to (not yet tested): (1) a family-level criterion that is hierarchy-aware (membership in the family at the coarser nested level, graph component beside the MCL cluster, with the precision of that level) so that a partition flip is not scored as a representation failure; (2) an evidence-based rule for the fusion cut (F1v2's) with the structural rule only as a candidate generator; (3) the seeding trade-off, which has opposite signs for reachable and entangled copies, on real data and on windows around other multi-copy regions.
