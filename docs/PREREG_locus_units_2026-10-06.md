# Pre-registration: a locus model for overlapping genes: primary-first representatives and structural separators (family level, ideal reads) (2026-10-06)

**Written before any arm of this study exists.** What was already seen (all on the same two dev windows, NPIP chr16 and TBC1D3 chr17, two read replicates each): the default arm D and the lever arms P (primaries-only seeding), C, PC at the FAMILY level (`bench/ideal_expression/score.py`, pooled over the four runs: reachable copies R, 54 gene-runs: D found 48, P 53, C 48, PC 53; entangled copies E, 24 gene-runs: D 12, P 7, C 13, PC 7), the per-copy differences D to P, and the anatomy of the entangled copies in D. DEV, human CHM13 / CAT-Liftoff, simulation, circular by construction. No held-out substrate is spent here; a held-out and a real-data check are the next registered step for any arm that passes.
User direction (2026-10-06): option B, a locus-level model for junction-sharing overlapping genes (units, multi-membership); the advisor wants clean combinatorial structure, provable statements, no arbitrary thresholds.

## 1. What the dev data say (the motivation, not a result of this study)

- A locus (the assembler's gene_id) is a connected component of the junction-sharing graph; replacing gene_ids by junction-sharing components reproduces the scores exactly. 216 of 256 entangled gene-runs share an exact junction with a partner gene (`docs/ENTANGLED_BASELINE_2026-10-06.md`): no component rule separates them.
- D loses reachable copies through the REPRESENTATIVE: NPIPA2 (representative is an artifact chain built from paralogous reads), TBC1D3P4 / P3 (the displaced-intron chain, carried only by the paralog's reads as secondary alignments, ties the copy's own chain at 10 reads), NPIPA7 and NPIPB3 (locus purity .112 and .234, bodies polluted with artifact chains). With primaries-only seeding all of them are found (R 53 / 54).
- Primaries-only seeding loses entangled copies from the family clusters (E4 22 -> 14: NPIPA1, AC138969.1 twice) and breaks others (NPIPA8 representative purity .906 -> .295): the secondary-borne chains enlarge the locus BODY toward its family, and the families stage aligns bodies. The representative and the body are different roles.
- Entangled copies in fused loci (NPIPB14P with the readthrough PDXDC2P-NPIPB14P as representative, AC138894.1 with CLN3 and NPIPB7, PKD1P6 with NPIPP1): the F1v2 read-share rule cannot cut a readthrough that carries as many reads as its parents (every simulated transcript has 10 reads).

## 2. The model (two rules, no free parameter)

**Rule 1, primary-first representative.** Within a locus a chain is PRIMARY-SUPPORTED iff the primaries-only assembly (the same run with `--no-seed-secondaries`) contains it exactly. The representative of a locus is the best primary-supported chain (most reads, then span, then index, the shipped key); only if the locus has none, the best chain overall. Bodies are untouched: every seeded transcript still belongs to its locus and enters the alignment of the locus body. Principle: secondary alignments may create a site and enlarge a body (O1 must not wait for an O2 coin toss), but structure and exactness come from the reads whose best placement is here.
**Rule 2, structural separators.** For a locus, the INCIDENCE GRAPH has the multi-exon transcripts and the junctions carried by >= 2 of them as vertices, an edge for each (transcript, junction) containment. A transcript is a SEPARATOR iff it is an articulation point of the graph (its removal splits the locus into >= 2 groups that each hold a transcript): the transcripts that are the only link between groups, the evidence-free generalisation of F1's bridge. Separators are removed from the families input (as F1's bridges are, listed in a side table), the remaining transcripts are regrouped by the components of the graph without them (gene_ids `<locus>.s<k>`), single-exon and junction-free transcripts attach to the component they overlap most, loci without a separator are unchanged. Facts used: the separator set and the components are canonical (order independent, linear time, Tarjan); the structure left is a tree of bridge-free groups; two groups linked by two vertex-disjoint paths cannot be separated by this rule (the boundary of a single-vertex rule, which the 216 junction-sharing genes of the dev windows largely sit behind).

## 3. Arms (the families input of each; every arm then runs `mcl_families` with the driver's families-stage flags, `--min-exonic-bp 1 --min-shared-exon-frac 0.60 --emit-units`)

| arm | families input |
|---|---|
| **D**, **P** | existing: the default `asm.families.gtf`, the primaries-only `lever_P.families.gtf` (same binaries) |
| **Q** | D's input with Rule 1: the `reads` attribute of every transcript is replaced by `reads + 1,000,000` when its chain is primary-supported (the representative key is (reads, span, -index); `reads` is used for nothing else in `fam_from_gtf`, checked in `src/family.rs`) |
| **S_D** | D's input with Rule 2 on all transcripts |
| **S_Q** | Q's input with Rule 2 where the graph holds the primary-supported transcripts only: the other multi-exon transcripts of a locus attach to the component sharing the most junctions with them (ties: the earliest), those sharing none form their own groups by junction sharing among themselves; separators are always primary-supported |

## 4. Measures

Family level, as the registered ideal test (`docs/PREREG_ideal_expression_2026-10-06.md`, `bench/ideal_expression/score.py`): IDEAL-FOUND (E2 and E3 and E4) per copy on stratum R (14 and 13 copies per run) and on the entangled copies E (10 and 2), E2-E4 separately, the registered rule (YES / PARTLY / NO, lower of the replicates), cluster precision, and `family_score` on the windows (Compara F, `bench/ideal_expression/fam_score.py`). Transcript level: `bench/entangled/arms_score.py` on each families input (chains recovered, artifacts, resolved genes). Pooled over the four runs and per family.

## 5. Predictions (fixed now)

- **Q1** Q finds >= 52 of the 54 reachable gene-runs (P: 53, D: 48). **Q2** Q finds >= 11 of the 24 entangled gene-runs (D: 12, P: 7): the body is kept. **Q3** Q gives YES for NPIP and for TBC1D3 in both replicates under the registered rule.
- **Q4** S_D finds >= 2 more entangled gene-runs than D with no loss at R. **Q5** S_Q finds at least as many gene-runs as Q at R and at E, and >= 1 more at E.
- **Q6** Compara F on the windows of Q and S_Q is >= D's minus .005 in both families.

## 6. Decision rule

An arm is a CANDIDATE iff, pooled over the four runs, IDEAL-FOUND (R plus E) is strictly above D's, R found is >= D's, E found is >= D's minus 1, and Compara F is >= D's minus .005 in both families; the direction must hold in both replicates of both families. A candidate earns a registered held-out check (ideal windows around other multi-copy regions, then A119b chr16 / chr17 and gorilla OR6737 with the copy-recovery instruments); nothing is claimed for real reads from this study. A non-candidate is reported with where it loses. Rule 2 is also reported by itself at the gene level against the oracle of the truth chains (`CUT-T`: 17 of 108 junction-sharing entangled genes resolved, 27 more non-entangled genes split than with components).

## 7. Limits declared in advance

Same two windows, same simulator, ideal reads, the arms were motivated by what the dev windows showed; the primary-supported marks come from a second assembly of the same reads (an approximation of 'the read's best placement is here'); Rule 2 uses no read evidence, so on real data it can cut a legitimate linking isoform (the truth oracle shows 5% of non-entangled genes split); the seeding gain for copies with no primaries is invisible in these windows.

## Amendment 1 (before any arm exists; the module and its unit tests were written after the section above)

Writing the unit tests showed that Rule 2 as worded ('a separator is an articulation point') also cuts an isoform that is the only link to a one-transcript attachment (A1 = j1 j2, A2 = j1 j2 j3, leaf = j3 j77: A2 would be removed from the families input). The truth-chain oracle of the dev windows (one replicate per family; connected components, three separator definitions; 108 junction-sharing entangled genes, 20 non-sharing, 513 non-entangled) chose the wording:

| rule | separators | junction-sharing entangled genes resolved / merged / split | non-entangled genes split (components: 128) |
|---|---|---|---|
| components (current loci) | 0 | 0 / 90 / 18 | 128 |
| articulation transcript (the first wording) | 92 | 17 / 64 / 27 | 155 |
| articulation transcript that contains the junction set of a transcript of two groups | 23 | 2 / 87 / 19 | 137 |
| **articulation transcript whose chain runs from a junction of one group, through a junction no group carries, to a junction of another group (private bridge)** | 38 | 16 / 68 / 24 | 140 |

**Rule 2 is the private bridge** (code: `bench/entangled/locus_units.py`, 12 unit tests, `python3 -m unittest test_locus_units`). It keeps nearly all of the articulation rule's gain with 12 false cuts instead of 27 and 38 separators instead of 92; it still resolves only 16 of 108 junction-sharing entangled genes (the other 68 are merged by shared junctions with no bridge: the boundary of a single-vertex rule, stated in section 2). The oracle is a dev analysis on the same windows and is disclosed as such; the family-level arms of section 3 are the test. Predictions Q1-Q6, the decision rule and Rule 1 are unchanged.
