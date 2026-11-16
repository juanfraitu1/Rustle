# PREREG row 1762: read-set Jaccard between de novo bundles as the family rule (2026-11-16)

> **Status 2026-11-16 (row 1762, A119b only, all uncommitted): every arm FAILS the registered bar (within 0.02 of the
> default 0.907 / 0.809 on both). Best per step, precision / recall:**
> components J >= 0.30 0.464 / 0.834 | bundlenode arms 0.000 / 0.93 | MCL I 2.0 0.479 / 0.795 (flat in t) |
> SEDEF gate + components 0.738 / 0.831 (t 0.10) | SEDEF gate + MCL 0.812 / 0.765 (t <= 0.10, best of 15 SD arms).
> Read sharing finds members (recall) but joins older relatives; the SD gate is the only lever that moved precision at
> kept recall. Gorilla never run. Open: SD link alone as a control, near-tie secondaries, Jaccard as a cut rule.

Written before any arm was scored. Advisor's question (via the user): are the family functionalities too complicated;
would "bundles like StringTie, then the Jaccard between bundles" define multi-copy gene families?

Earlier tests: sequence Jaccard + components (`archive/2026-09/JACCARD_SWEEP_RESULT_2026-09-20.md`), read-set Jaccard on
ANNOTATED genes, three contigs, bipartite F (`archive/2026-09/ADVISOR_JACCARD_ANCHORS_2026-09-25.md`; its caveat: the de
novo version was not run). This row runs the de novo version genome-wide with today's readout.

- **Nodes (bundles):** the 76,142 loci of the default run (`PREFIX.fam.loci.gff3`: overlapping assembled transcripts merged
  into one locus), unchanged. Only the edge rule changes.
- **Reads of a bundle:** molecules with a primary or secondary record (supplementary records ignored) placing >= 20 bp on
  the bundle's exons.
- **Shared molecule of A and B:** one record on A and not B, a different record on B and not A (multi-mapping; a single
  record spanning both never counts). Pairs whose spans overlap are never linked.
- **Edge:** J = shared / (n_A + n_B - shared) >= t, t in {0.01, 0.05, 0.10, 0.20, 0.30}, with shared >= 1 (literal) or
  shared >= 2 (control): 10 arms. Families = connected components. Tool `bench/denovo_eval/bundle_jaccard.py`.
- **Readout:** `pairscore.py` default level, pair precision and recall, against the default table of the same run
  (A119b `r1749_fam_darklist/base.fam.copies.tsv`, reproduced today 0.907 / 0.809).
- **Rule.** An arm MATCHES when precision >= default - 0.02 AND recall >= default - 0.02. No arm matches on A119b: STOP,
  the simple rule is not a replacement. Some arm matches: freeze the matching arm with the highest recall, run it once on
  gorilla OR6737 (`r1674_edgerule/or_all`, protein-hit truth), same bar: PASS / FAIL.
- **Descriptive:** false pairs that are `related`; arm edges held by one SEDEF pair (co-duplicated).
- Not controlled: the default clusters with MCL, the arms with components (09-20: operator worth +0.035 F).

## Amendment 1 (2026-11-16, written while the first job was still counting; no arm scored yet): bundlenode arms

User's proposal: cut each bundle into StringTie-style bundlenodes and compare those. The locus exons of
`PREFIX.fam.loci.gff3` are already merged, non-overlapping covered intervals (251,634 in 76,142 loci) = the bundlenodes.
- **Linked nodes:** node x (bundle A) and node y (bundle B) are linked when >= 2 molecules have one primary/secondary
  record with >= 20 bp on x and not on y and another record on y and not on x.
- **Edge A-B:** linked node bp of the bundle with the smaller exon length / that exon length >= T, T in {0.3, 0.5, 0.7}
  (the read-built form of the default exon-share gate, 0.5). Overlapping bundles never linked. Components.
- Same readout and the same MATCH rule; these three arms join the ten above for the A119b decision.

## Amendment 2 (2026-11-16, both jobs still running; no arm scored yet): strand

User: StringTie builds bundlenodes per strand. The first two jobs ignore strand: a read counts on any bundle it overlaps,
so an antisense bundle inherits the links of the gene it overlaps. Every locus is stranded (49,275 + / 26,867 -) and the
reads are mapped `-uf` (alignment strand = transcript strand). **Registered arms are now the STRANDED ones** (a record
counts only on bundles of its own strand; tags `a1s`, `a1ns`); the unstranded runs (`a1`, `a1n`) are reported beside
them as a control, not for the decision.

## Result (2026-11-16, A119b, stranded arms; workspace `r1762_bundle_jaccard`, `a1xs.score.tsv`, `a1xs.arms.tsv`)

**STOP as registered: no arm is within 0.02 of the default on precision.** Gorilla not run.

| arm | precision | recall | false pairs (related) | largest family (loci) |
|---|---|---|---|---|
| default table | 0.907 | 0.809 | 180 (110) | - |
| J >= 0.01 | 0.001 | 0.960 | 2,040,568 (40,626) | 19,667 |
| J >= 0.05 | 0.008 | 0.923 | 242,133 (9,214) | 7,474 |
| J >= 0.10 | 0.160 | 0.919 | 10,512 (3,203) | 5,000 |
| J >= 0.20 | 0.360 | 0.869 | 3,375 (2,370) | 4,995 |
| J >= 0.30 | 0.464 | 0.834 | 2,105 (1,783) | 4,863 |
| bundlenodes T 0.3 / 0.5 / 0.7 | 0.000 | 0.932 / 0.932 / 0.931 | 5.6-7.9 million | 22,843 / 21,980 / 20,799 |

- The shared >= 2 control equals shared >= 1 from t = 0.05 up. 13.9 million bundle pairs share a molecule; 53.4 million
  node pairs share >= 2. SEDEF-held edges 2.6% (J 0.01) to 9.0% (J 0.30) of ALL edges (dominated by the satellite block; says nothing about the false gene pairs).
- Recall is above the default at t <= 0.30; precision fails in every arm.
- Post hoc (descriptive, `ph.score.tsv`). The largest family at t 0.10-0.30 is the chr13:13.9-15.5 Mb satellite hotspot
  (4,398 loci, + 461 on chr21): it holds no coding pair, and the scores without it are identical. Outside the grid:
  t 0.5 = 0.503 / 0.766, t 0.7 = 0.581 / 0.594: precision keeps rising but recall falls below the default before
  precision reaches 0.6. At t 0.30 the false pairs sit in 270 families (default 76); one KRAB-ZNF family of 47 genes holds
  1,066 of them (0.64 precision without it); 1,783 of 2,105 are `related` (older paralogues / HGNC group).
- Not controlled: components vs MCL (default); no tie filter on the secondaries (tie-group recipe, rows 1603-1631).
- Cost: BAM pass 22 min at 36 processes; shared-molecule counting 2 h 45 min at 12 workers (64 at once ran out of memory);
  the default families stage takes ~7-14 min (this prototype, not the rule: most of the cost is the satellite block). Unstranded control not completed (timed out at 34 of 64 groups; not needed for
  the verdict). In the fast job a record counts on a bundle with >= 20 bp on one node (not summed over nodes).

## Amendment 3 (2026-11-16, AFTER the Result above; no MCL arm scored yet): MCL on the Jaccard graph

User request. The component scores above have been seen; this closes the "components vs MCL" control, it does not reopen
the STOP on the component arms.
- **Graph:** the same stranded bundle pairs (`a1xs.pairs.tsv`, shared >= 1, overlapping bundles never linked), edge
  weight = raw J, kept at J >= t, t in {0.01, 0.05, 0.10, 0.20, 0.30}.
- **Operator:** `bench/mcl_port.py` (`utilities mcl-port`: self-loops 1, prune 1e-9), run per connected component (exact:
  the matrix is block-diagonal). This is the bench MCL, NOT bit-identical to the shipped `mcl_families` MCL (no
  two-node rule, no fold). Clusters of one locus are dropped.
- **Decision arms: inflation 2.0** (the pipeline default). Inflation 4.0 is run beside it, descriptive only.
- **Feasibility rule (fixed now):** a component whose MCL has not finished in 50 min stays one family and is listed. At
  t 0.10-0.30 the largest component is the chr13 satellite block with no coding pair, so that cannot move those scores;
  at t 0.01 / 0.05 it can, and such an arm is reported as "largest component unclustered".
- **Bar unchanged:** MATCH = within 0.02 of the default on both (0.907 / 0.809). A MATCH triggers the gorilla run, not
  adoption. Also reported: hit pairs lost vs the component arm at the same t, and what happens to the 47-gene KRAB-ZNF
  family (1,066 of the 2,105 false pairs at t 0.30).

### Result of Amendment 3 (2026-11-16, job 4477110, 55 min; `a1xs.mclscore.tsv`, `a1xs.mclarms.tsv`)

**FAIL: no MCL arm is within 0.02 of the default on precision (best 0.479 vs 0.907).** Every component finished (none left
unclustered). Gorilla not run.

| t | components (seen) | MCL I 2.0 (decision) | false (related) | largest | MCL I 4.0 (descriptive) |
|---|---|---|---|---|---|
| 0.01 | 0.001 / 0.960 | 0.403 / 0.854 | 2,761 (2,381) | 4,974 | 0.430 / 0.797 |
| 0.05 | 0.008 / 0.923 | 0.420 / 0.841 | 2,540 (2,196) | 4,967 | 0.445 / 0.787 |
| 0.10 | 0.160 / 0.919 | 0.446 / 0.830 | 2,255 (1,951) | 4,425 | 0.440 / 0.727 |
| 0.20 | 0.360 / 0.869 | 0.464 / 0.807 | 2,033 (1,754) | 1,964 | 0.454 / 0.727 |
| 0.30 | 0.464 / 0.834 | 0.479 / 0.795 | 1,892 (1,654) | 1,470 | 0.464 / 0.710 |

- MCL removes the threshold dependence (precision 0.40-0.48 at every t, against 0.001-0.46 for components): it cuts the
  weak chaining links, not the false gene pairs. At t 0.30 it drops 213 false pairs and 84 hit pairs (1,820 -> 1,736).
- MCL only splits components, so at t 0.30 at least 853 of the KRAB-ZNF family's 1,066 false pairs survive (only 213
  false pairs fell in total). Illustration: the chr19 component of 129 loci (96 on chr19) stays one cluster of 118.
  MCL at I 2.0 does not separate these dense, high-J neighbourhoods; other partitioners were not run.
- **Inflation 4.0 is prune-confounded, not a result:** with the absolute prune 1e-9 a near-uniform block above ~178 nodes
  loses its columns at I 4.0. At t 0.30, 3,971 of the 4,863 satellite-component loci (6,657 of 15,361 overall) drop out
  as singletons, against 96 (866) at I 2.0. The decision arm is safe (limit ~31,600 nodes at I 2.0).
- Limit of the weights: raw J sits against self-loops of 1, so low-J edges are nearly inert (hence flat in t) and a lone
  weak pair is split (part of the 84 lost hits). That biases towards MORE splitting, so it cannot explain the precision
  failure.
- This replaces "Not controlled: components vs MCL" above. Operator = bench `mcl_port`, not the shipped MCL.

## Amendment 4 (2026-11-16, after the Amendment 3 result; no SD arm scored or counted yet): SDs on the Jaccard graph

User question: can SDs be included. The default ALREADY uses them (`--sd-damp 0.25`, row 1580: edges between loci no SEDEF
pair links are multiplied by 0.25 before MCL), so this gives the Jaccard graph the same information as the baseline.
- **SD-linked edge:** the default's span rule (as `bench/denovo_eval/sdbio.py` `holds`): bundle A sits in one copy of a
  SEDEF pair and bundle B in the OTHER copy (>= half of the exon bp and >= half of the span inside the copy). BED = the
  default run's, `~/data/chm13/chm13v2.0.fa.sedef.bed`. (The `edges_sedef_held` column of the first Result used a looser
  count that did not tell the two sides apart; it was descriptive only and is not reused.)
- **Arms, t in {0.01, 0.05, 0.10, 0.20, 0.30}, same graph as Amendment 3:**
  `SJ` = J >= t AND SD-linked, components; `SM` = the same edges, MCL I 2.0; `DM` = all J >= t edges, weight J x 0.25
  when not SD-linked, MCL I 2.0 (the default's own form; against self-loops of 1 it should behave close to the gate).
  No inflation 4.0 arm (prune-confounded, Amendment 3).
- **Decision:** all three arm families, bar unchanged (within 0.02 of 0.907 / 0.809 on both). MATCH triggers the gorilla
  run, not adoption. Recall of `SJ` at t 0.01 is read as the ceiling of an SD gate on this graph.
- Caveat fixed now: SEDEF is a self-alignment of the reference, not read evidence; with it the rule is no longer
  "bundles + Jaccard" alone. New outputs `a1xs.sdarms.tsv`, `a1xs.sdscore.tsv`; Amendment 3 files untouched.

### Result of Amendment 4 (2026-11-16, job 4477345, 48 min; `a1xs.sdscore.tsv`, `a1xs.sdarms.tsv`)

**FAIL as registered: no SD arm is within 0.02 of the default on both (bar 0.887 / 0.789). Best: SD gate + MCL
0.812 / 0.765 (t <= 0.10).** Gorilla not run. (A first job, 4477323, died in 34 s on a wrong boundary case in the
`holds` self-check; the rule itself was not changed.)

| t | SJ gate, components | false (related) | SM gate, MCL I 2.0 | false (related) | DM damp 0.25, MCL I 2.0 |
|---|---|---|---|---|---|
| 0.01 | 0.602 / 0.837 | 1,207 (563) | 0.812 / 0.765 | 387 (301) | 0.505 / 0.826 |
| 0.05 | 0.715 / 0.833 | 725 (536) | 0.812 / 0.765 | 387 (301) | 0.511 / 0.819 |
| 0.10 | 0.738 / 0.831 | 643 (510) | 0.812 / 0.764 | 386 (300) | 0.511 / 0.816 |
| 0.20 | 0.784 / 0.783 | 471 (362) | 0.814 / 0.729 | 364 (286) | 0.515 / 0.777 |
| 0.30 | 0.785 / 0.758 | 454 (358) | 0.811 / 0.673 | 343 (286) | 0.511 / 0.720 |

- The SD gate is the first lever that raises precision while recall stays at or above the default: false pairs
  2,105 -> 454 at t 0.30 (components), 1,892 -> 343 with MCL; the satellite block is gone (largest family 258-924 loci,
  against ~4,900).
- No arm passes both bars. Gate + components passes recall at t <= 0.10 (0.831-0.837; the gate's ceiling 0.837 is not
  the limit) but reaches only 0.738 precision; gate + MCL fails both (0.812 < 0.887, 0.765 < 0.789).
- Gate + MCL is flat in t up to 0.10 (the threshold is not needed). Against gate + components at t <= 0.10 it gives up
  0.07 recall for 0.07-0.21 precision; at t 0.30 0.085 recall for 0.026 precision.
- What separates the best arm from the default: 301 of its 387 false pairs are `related` (older paralogues / HGNC
  group) and SD-linked, against the default's 110 of 180; the unrelated remainder is similar (86 vs 70); and 96 fewer
  hits (1,669 vs 1,765). With equal SD information the remaining gap is what the sequence + exon-share test adds.
- Damp 0.25 on raw J does not cut these neighbourhoods (0.51 precision): where every edge of a neighbourhood is
  unlinked they are rescaled alike. Not evidence about the default's damping; the weight scale was not tuned.
- 0.812 / 0.765 is the best of 15 arms on a spent sample, each amendment written after the previous scores: not an
  expected performance. Not done: SD link alone without Jaccard (control), gorilla, which SD pairs hold the 301.
