# Pre-registration: attacking the named losses of the Soto 2025 replication, held out by Soto family

**Written 2026-09-29 (KEY=soto_losses) after a development pass on the DEV half only and before any held-out number.**
Human CHM13 v1.0 only (Soto's genome, annotation, gene universe and Table S1C famCN). This file binds once Amendment 1
records its sha1 and the frozen instruments' sha1s; until then nothing is scored on the held-out half. Every fix is
opt-in analysis code under `/mnt/linuxdisk/tmp/rustle_figures_dev/soto_losses/`; `bench/soto/soto_replication.py` and
`src/` are not edited. Whatever the verdict, this is concordance with Soto's own tables, not an independent
replication (register 858 / 1085: famCN, universe and the 71 curated genes are theirs).

## 0. The question

The best replication so far is Soto's literal recipe on native CHM13 v1.0 (`docs/archive/2026-09/SOTO_REPLICATION_STATUS_2026-09-28.md`
§1.1): ARI 0.7096 / 235 of 491 families exact (median-MAD), 0.7039 / 261 (mean-MAD). Its non-exact families were
attributed (median arm) to over-merge we cannot split (127), edges we lack (105: 62 fragments at 90-98% identity with no
map-back hit, **49 where SEDEF has a >= 98% row but minimap2 finds nothing**, 17 linked only via non-eligible members),
and copy-number split (24: 12 families that break Soto's own MAD < 1, **20 avoidable cuts by our greedy split**). Three
mechanisms learned since 09-28 were proposed against them: (a) readthrough / fusion bridges (the fusion-container
simulation and F1 / F1v2: one readthrough transcript merges two gene groups; Soto's STAR Methods drop
"readthrough-fusion transcripts" by hand), (b) map-back sensitivity (`-N 50`, `-p 0.5`, short exons, masking) for the
49, and (c) a threshold-free alternative to the greedy famCN split. Which of them explains a loss on DEV, and does the
resulting fix, stated in Soto's own vocabulary with no fitted constant, beat both the literal recipe and a matched
random null on the HELD-OUT half of Soto's families?

## 1. The split (fixed 2026-09-29 08:41, before any Phase-1 analysis)

`frozen/split.py` 7b39767d -> `frozen/split.tsv` 49bcbcfe; inputs S1C d008a179 and the BASE edge set `frozen/base_edges.tsv`
6cc14eeb (= `soto_v1_literal/arms/finalv1/edges.project.tsv`, 11,935 edges).

1. Every S1C multigene Family ID F: `h(F) = int(sha1("soto_losses_2026-09-29|" + F), 16) % 2`; 0 = DEV, 1 = HELD-OUT.
2. Units = connected components of a graph on the 2,334 S1C genes joining (i) S1C co-membership (a cover gene joins all
   its families) and (ii) every BASE shared-exon edge. Every Soto family, every BASE over-merge and every BASE predicted
   family therefore lies inside one unit (BASE never scores a predicted family across halves).
3. A unit whose families all hash to one half takes it. A cross-half unit takes the half of its LARGEST family (most
   S1C genes; tie -> smallest Family ID string). A unit with no multigene family takes h(smallest gene id).
4. DEV 184 units / 1,170 genes / 225 families (19 cross-half units); HELD-OUT 213 units / 1,164 genes / 266 families
   (24 cross-half units). Largest unit 220 genes (DEV).
5. A half is scored with `soto_replication.py`'s own `score()` and `bipartite_score()` on the clean (single-family) S1C
   genes of that half (DEV 1,089 scored genes). A fix that ADDS edges can put genes of both halves in one predicted
   family; those cross-half pairs are reported (`cross_pairs`) and counted as false positives in `P_guard`, the pair
   precision used by the no-harm clauses.

## 2. Development evidence (DEV half only; every number here was seen before this file)

Harness `frozen/lib.py` reproduces BASE (`cluster --full-geneset`) partition-for-partition in both MAD statistics
(asserted by `lib.check_base()` at the top of every run); DEV BASE: median ARI 0.5957 / 104 of 225 exact, mean 0.6118 / 119.

**(a) Readthrough / fusion bridges: NO mechanism.** Readthrough / bridge genes over the 2,334 = GENCODE
`readthrough_transcript` tag on any CAT v4 source transcript (63 genes, 246 transcripts; tags from
`gencode_chm13/chm13v2.0_gencode.gff3` joined on ENST), GENE1-GENE2 names (15), or an F1-style structural bridge (a
transcript whose exons overlap exons of two CAT genes that do not overlap each other, 66): 76 genes, `frozen/rt_genes.tsv`.
Most are Soto cover genes (PKD1P6-NPIPP1 sits in 6 families, AC243919.2 in 4) and non-eligible lncRNAs, which the
eligible-only backbone never lets bridge. Of 781 DEV cross-family edges, 9 touch one (1 of 726 backbone edges); deleting
all 76 disconnects 1 of 61 backbone-joined family pairs inside DEV over-merged predicted families (median arm; 1 of 20
mean). Demoting the 10 eligible ones to attach-only members changes DEV exact by 0 and ARI by 0.0000 (null
indistinguishable). The DEV over-merges are instead high-identity exon links between families at different copy number:
674 family-to-family backbone edges, median gap-excluded identity 0.9948 (569 >= 0.98), median famCN gap 15.1 (617 >= 2),
248 (37%) between families Soto itself links through a shared cover gene. Carried as a NEGATIVE-CONTROL arm (A), not a fix.

**(b) Map-back sensitivity: the cap is not it; two sub-mechanisms.** DEV holds 22 of the "SEDEF >= 98% row, no map-back"
fragments (12 families; 107 evidence gene pairs; 169 SEDEF rows; 48 units). No unit involved reaches the `-N 50` cap
(at most 6 secondaries). Rerunning 25 of the 48 units (those of ID_290-293 / ID_47 / ID_477, which share one unit, and of
ID_388, ID_131, ID_113; the >= 385 kb units of ID_321 / ID_328 / ID_363 were not rerun) against a freshly built full v1.0
index: Soto's settings reproduce the original 107 alignments (columns 1-12 identical); `-N 1000` gives the byte-identical
107; `-p 0` (with `-N 50` or `-N 1000`) gives 25,449 alignments.
- **(b-i) self-paired SD98 rows (7 families: ID_290-293, ID_47, ID_477, ID_388).** The SEDEF row pairs an interval with
  ITSELF (side1 = side2, strands +/-: an inverted duplication inside one SD98 unit). The unit's map-back reports only its
  self-hit: the internal copy aligns over at most part of the unit, scores below half of the self-hit, and `-p 0.5`
  drops it. SEDEF's own CIGAR walk links all 7. `-p 0` links all 7 (16 of the 89 distinct evidence pairs of the 12
  families, 7 of them AMY pairs of ID_131) but adds 124 edges beyond BASE from these 25 units alone, only 49 of them
  between genes sharing a Soto Family ID (0.40): removing `-p` is not a fix.
- **(b-ii) region co-containment only (5 families: ID_113 GOLGA8, ID_131 AMY, ID_321 / ID_328 CU6339xx, ID_363 FRG1BP).**
  A full-length map-back alignment onto the partner side EXISTS, but no exon of one gene projects onto an exon of the
  other: the "SEDEF >= 98% row contains both exons" evidence is region-level (one SD contains an exon of each gene, at
  different offsets). Neither the map-back nor SEDEF's own CIGAR has an exon-to-exon >= 98% correspondence. Not a
  sensitivity problem; no fix proposed. Short-exon chaining and minimizer masking were not separately tested: the loss is
  at the whole-unit alignment (b-i) or has no exon correspondence at all (b-ii).

**(c) The greedy famCN cut: a threshold-free split avoids it.** 11 DEV families (median arm; 2 mean) have all clean
eligible members in one backbone component, clean-member MAD < 1, and are cut by the ascending-famCN greedy walk. The
minimum-segmentation split (§3, fix C) keeps 10 of 11 whole (2 of 2 mean).

**DEV scores of the frozen arms** (`runs/dev.tsv`, `runs/dev.nullcmp.tsv`; 20 null seeds; "null>=" counts null draws at
least as good as the fix):

| arm | median ARI | median exact | median P / R / F1 | median MICRO P / R | undet | mean ARI | mean exact | mean MICRO R |
|---|---|---|---|---|---|---|---|---|
| BASE | 0.5957 | 104 | .812 / .474 / .599 | .767 / .672 | 54 | 0.6118 | 119 | .697 |
| A (control) | 0.5957 | 104 | .813 / .474 / .599 | .768 / .671 | 55 | 0.6118 | 119 | .695 |
| B1 | 0.5959 | 109 | .808 / .476 / .599 | .768 / .687 | 46 | 0.6118 | 124 | .710 |
| C | 0.6239 | 107 | .811 / .511 / .627 | .787 / .700 | 53 | 0.6129 | 121 | .699 |
| ALL (B1 + C) | 0.6233 | 112 | .805 / .513 / .626 | .786 / .714 | 46 | 0.6129 | 126 | .713 |
| nullB1 (median of 20) | 0.5962 | 107 | .804 / .478 / .599 | .766 / .682 | 49 | 0.6136 | 122 | .707 |
| nullC | 0.6177 | 105 | .802 / .507 / .621 | .779 / .693 | 55 | 0.6174 | 119 | .699 |

DEV null comparisons: B1 exact null>= 2/20 (both statistics), MICRO R 1/20 / 2/20, ARI 11/20 / 18/20. C median ARI 1/20,
exact 0/20, MICRO R 0/20, MICRO P 0/20; C mean ARI 17/20, exact 6/20. A: every metric within its null.
DEV exact gains / losses vs BASE (median): B1 +ID_290-293, ID_477 / none; C +ID_104, ID_280, ID_352, ID_485 / -ID_254.

## 3. The fixes (binding; code `frozen/lib.py` ff0d105a and `frozen/run_arms.py` 74bd75dc)

- **B1 — "a unit that overlaps its own paralog is read through its own alignment."** For every SEDEF row with
  fracMatch >= 0.98 (Soto's SD98 floor) whose two sides overlap on one chromosome (565 of 3,628 rows; 76 give edges), the
  shared-exon test (`-f 0.99` exon coverage) is made through the row's own CIGAR (`soto_replication.py`'s
  `build_blocks` / `find_shared_exons`, unchanged; per-row edges `frozen/cigar_rows.py` 7d08ed10 -> `cigar_rows.tsv`
  2d7175a9, whose union over all 3,628 rows equals the 09-28 CIGAR-walk arm's 4,384 edges, asserted). Their edges are
  added to BASE: **38 new edges**. Reason in Soto's words: the map-back sets each unit's bar at 0.5 x its own self-hit
  (`-p 0.5`), so a paralog lying inside or across the unit's own interval can never pass it.
- **C — "the fewest famCN-coherent groups."** A component whose famCN MAD >= 1 is split into the MINIMUM number of groups
  that are contiguous in sorted famCN and each satisfy Soto's own MAD < 1 (the component-level test and every other step
  are BASE's). Among minimal splits: the smallest total absolute deviation about each group's centre (median in the
  median arm, mean in the mean arm); remaining ties -> earliest cut positions. Exact dynamic programme (`split_dp`),
  checked against brute force on 300 random instances x 2 statistics (0 mismatches).
- **ALL** = B1 edges + C split.
- **A (negative control)** = the 10 eligible `rt_genes.tsv` genes that have a BASE backbone edge are demoted to
  attach-only members ("a readthrough is a member, not a bridge").

| element | value | source |
|---|---|---|
| SD98 floor | fracMatch >= 0.98 | Soto |
| exon coverage | 0.99 | Soto (`-f 0.99`) |
| "overlap" (B1) | >= 1 bp, same chromosome | definition, no constant |
| MAD < 1, statistic median / mean | 1.0 | Soto (code: median; prose: mean) |
| minimum groups; L1 tie-break; earliest cut (C) | — | design choices fixed in code BEFORE any DEV number was computed |
| readthrough set (A) | tag / name / structural | GENCODE tag, name, F1 structure |

## 4. Arms, nulls, metrics

Arms BASE, A, B1, C, ALL, each under median-MAD and mean-MAD. Nulls, seeds 0..19 (`run_arms.py`):
**nullB1** = BASE + as many edges as B1 adds (38) drawn at random from the SAME pool (every CIGAR-walk edge not in BASE,
66 edges); **nullC** = per split component, C's number of groups with cut positions drawn uniformly among all valid
(MAD < 1 per group) famCN-contiguous segmentations; **nullALL** = nullB1(seed) edges + nullC(seed) split; **nullA** = as
many random eligible non-R genes with a backbone edge demoted (10). Metrics per arm x statistic on the held-out half:
ARI, exact families (of 266), pair P / R / F1, `P_guard`, `cross_pairs`, bipartite MICRO and MACRO P / R (scipy tie
policy, register 1045), undetected families, predicted families. Also `all` (491 families) after the verdict, labelled
DEV + HELD-OUT.

## 5. Clauses (held-out half; median-MAD is primary — their released code's statistic and our headline)

"null>= k/20" counts null draws at least as good as the fix (ties count against the fix; for undetected, <=).
- **B1 EFFECTIVE** iff, median arm: (B1.1) exact Δ >= +2 vs BASE; (B1.2) exact null>= <= 2/20; (B1.3) MICRO R Δ > 0 and
  undetected Δ <= 0; and in BOTH arms (B1.4) `P_guard` Δ >= -0.010, ARI Δ >= -0.005, MICRO P Δ >= -0.010, exact Δ >= 0.
- **C EFFECTIVE** iff, median arm: (C.1) ARI Δ >= +0.005 and MICRO R Δ > 0; (C.2) ARI null>= <= 2/20; (C.3) exact Δ >= 0;
  and (C.4) in BOTH arms `P_guard` Δ >= -0.010; mean arm ARI Δ >= -0.005 and exact Δ >= 0.
- **ALL EFFECTIVE** iff B1 and C are EFFECTIVE, median ALL exact >= max(B1, C) exact, and nullALL exact null>= <= 2/20.
- **A (control) holds** iff |exact Δ| <= 1 and |ARI Δ| <= 0.002 in both arms. **Falsified** if exact Δ >= +3 or ARI Δ >=
  +0.005 in either arm (then readthrough bridges DO carry held-out over-merges and §2(a) is wrong).
- A fix that fails any clause is **NOT EFFECTIVE**; one that passes the Δ clauses but not its null is **NOT DISTINGUISHABLE
  FROM ITS NULL** (the curation-rule lesson, `SOTO_REPLICATION_STATUS` §2.1).

## 6. Predictions (this author's probabilities, before any held-out number)

B1.1 0.65 (the held-out half holds known union gains, §8.1); B1.2 0.40 (its pool null contains 58% of its own edges);
B1 EFFECTIVE 0.35. C.1 0.70; C.2 0.45; C EFFECTIVE 0.40. ALL EFFECTIVE 0.15. A holds 0.90. Mean-arm C: no effect.

## 7. Falsifiers of the mechanisms (reported whatever the verdict)

- **F-B1:** fewer than half of B1's held-out exact gains are families with a member incident to a B1 edge -> the gains are
  not the self-paired mechanism.
- **F-C:** of the held-out families cut by greedy although their clean MAD < 1 (median), C keeps fewer than half whole
  -> the avoidable-cut diagnosis does not generalise.
- **F-A:** more than 2% of held-out cross-family backbone edges touch a readthrough / bridge gene.
- **F-guard:** `cross_pairs` of B1 / ALL on the held-out half exceeds 1% of its predicted pairs -> the half-scoring is
  leaking and every B1 number is reported unjudged.

## 8. Hostile self-review (applied above)

1. **The held-out half is not blind for B1's pool.** On 09-28 the union of the literal recipe with ALL 66 CIGAR-walk-only
   edges (B1's pool) was scored on all 491 families and its per-family diff was read: union gains ID_118, 217, 317, 438,
   479 and losses ID_344/345 are ALL held-out families (checked from the split only). Whether B1's 38 edges are the ones
   behind them was NOT looked up. So B1.1 is partly foreseeable; B1.2 (beats random 38 of the same 66) is the decisive
   clause, and it is conservative (the null's draws share 58% of B1's edges on average).
2. **C's largest known cases are held-out.** The 09-28 diff named ID_356 (56 -> 53) and ID_154 (22 -> 17) as greedy-cut
   families; both are held-out. C.1 is partly foreseeable; C.2 (DP positions vs random valid positions with the SAME
   number of groups) is the decisive clause. C.2 tests the choice of cut positions, not the reduction in the number of
   groups: a gain that comes only from cutting fewer times is shared by nullC and fails C.2. DEV already shows most of C's
   ARI gain is shared by nullC (0.6177 vs 0.6239) — said here so a pass is not over-read.
3. **BASE's held-out numbers are derivable** (full 0.7096 minus DEV) and were known in aggregate; BASE is not a candidate.
4. **B1's definition was generalised after seeing DEV:** the 7 DEV families come from 10 identical-sides rows; the rule
   uses all 565 overlapping-sides rows (tandem shifts included) because the stated reason (`-p 0.5` against the unit's own
   self-hit) applies to both. Fixed here, not tuned on a score (the identical-sides-only variant was never scored).
5. **Multiple metrics.** Each fix has one primary Δ clause and one null clause, named in advance; the rest are no-harm.
6. **Nulls are matched in count, not in effect size of each element.** nullB1 draws from the same 66-edge pool
   (per-edge size varies); nullC keeps each component's group count. An edge-count-matched null proves nothing on its own
   (metric-trap list); here it is paired with the within-pool restriction.
7. **Cover genes are excluded by the scorer**, which caps the truth at 1.0 (not 409/491); famCN is Soto's; the 71
   curated genes are Soto's. The verdict is about concordance with S1C, not about biology.
8. **Tie policy of bipartite matching** is scipy's (register 1045); pair P / R / F1 and ARI are tie-free and carry the bars.
9. **Readthrough set from GENCODE on v2.0** (CHM13 v2.0 CAT/GENCODE GFF3) joined on the versioned ENST of CAT v4 (v1.0):
   a tag is a property of the source transcript, so the join is version-safe; transcripts without an ENST get no tag.

## 9. Order, stop rules, machine rules

Freeze (Amendment 1) -> `run_arms.py heldout runs/heldout` (light lock, ~15 s) -> clauses -> `run_arms.py all runs/all`
(reported, labelled DEV + HELD-OUT) -> Outcome. No arm, null, clause or bar changes after the freeze; an instrument fix
that changes no rule is recorded as an Amendment. Light jobs only (`tools/rlock.sh light`); the one heavy job of this
study (the 25-unit diagnostic map-back and its index) ran during Phase 1.

## 10. Not in this test

The over-merge bucket (127 families; 59 where the true family is a sub-half minority) has NO fix here: its mechanism
(§2a) is not readthrough; a split that sees a sub-half famCN minority needs a new rule (e.g. MAD < 1 under BOTH of Soto's
own statistics) and is left for a separate pre-registration. The 62 fragments at 90-98% identity and the (b-ii)
co-containment fragments are below the SD98 recipe by construction. The 12 families that break Soto's own MAD < 1 rule.
No default is changed; nothing is committed.

## Amendments

### Amendment 1 — the freeze (2026-09-29 09:05, written BEFORE any held-out command)

**This file's sha1 before this amendment:** `c10e0157ecb79b833f6c9cdc81a7d814c75c8501` (17,656 bytes; byte copy
`/mnt/linuxdisk/tmp/rustle_figures_dev/soto_losses/PREREG_soto_losses_2026-09-29.pre_amendment1.md`). The frozen text
equals this file with the Amendment and Outcome bodies removed.
**Frozen instruments:** `/mnt/linuxdisk/tmp/rustle_figures_dev/soto_losses/frozen/SHA1SUMS` sha1
`1842abbafff2efb07e6a4060162baff17cc5d902` (`sha1sum -c` OK at the freeze): `split.py` 7b39767d, `split.tsv` 49bcbcfe,
`base_edges.tsv` 6cc14eeb, `lib.py` ff0d105a, `cigar_rows.py` 7d08ed10, `cigar_rows.tsv` 2d7175a9, `rt_genes.tsv`
7581901e, `run_arms.py` 74bd75dc, `score_half.py` f0500ea0 (Phase-1 helper, not used by the run). Library
`bench/soto/soto_replication.py` 13ffb698 (unchanged), S1C d008a179, `final_v1_clean.bed` acc17cbe (the trailing header
row of `final_v1.bed` e4459a4d stripped). Python 3.14.4. Nothing had been scored on the held-out half or on `all` for any
arm other than BASE when this was written; `runs/` held only `dev.*`.
**Command:** `bash tools/rlock.sh light python3 frozen/run_arms.py heldout runs/heldout`, then `... all runs/all`.

### Amendment 2 — the run (2026-09-29 09:07-09:12; no rule, arm, null, clause or bar changed)

Both runs as written, `sha1sum -c` of the frozen set OK before each; BASE reproduced (asserted) in both statistics; no
instrument fix was needed. Held-out scored universe 1,096 clean genes. Cross-check: every arm's partition was written
out and re-scored by the separate `score_half.py` path (held-out, 8 arm x statistic cells: identical to 4 decimals) and
by `soto_replication.py score` on the full set (ALL / B1 / C median: identical).

## Outcome (2026-09-29)

**Verdict: no fix is EFFECTIVE.** B1 = NOT DISTINGUISHABLE FROM ITS NULL; C = NOT EFFECTIVE (fails its pre-registered
mean-arm no-harm clause, while passing every median-arm clause including its null); ALL = NOT EFFECTIVE (needs both);
the readthrough control A HOLDS (readthrough bridges do not carry the over-merges, held-out). No falsifier fired.

**Held-out half (266 families; `runs/heldout.tsv`, `runs/heldout.nullcmp.tsv`, 20 null seeds):**

| arm | median ARI | median exact | median P / R / F1 (`P_guard`) | median MICRO P / R | undet | mean ARI | mean exact | mean P / R | mean MICRO R | mean undet |
|---|---|---|---|---|---|---|---|---|---|---|
| BASE | 0.8183 | 131 | .845 / .796 / .820 | .800 / .763 | 59 | 0.7933 | 142 | .887 / .720 | .767 | 54 |
| A (control) | 0.8191 | 131 | .849 / .794 / .821 | .800 / .761 | 59 | 0.7931 | 142 | .890 / .718 | .765 | 54 |
| B1 | 0.8194 | **134** | .845 / .798 / .821 | .798 / **.774** | **53** | 0.7944 | 145 | .888 / .722 | .778 | 48 |
| C | **0.8520** | 132 | .859 / .848 / .853 | .811 / .776 | 60 | **0.7758** | 147 | .892 / .689 | .769 | 53 |
| ALL | 0.8529 | 135 | .859 / .850 / .854 | .808 / .788 | 54 | 0.7770 | 150 | .892 / .691 | .780 | 47 |
| nullB1 (median of 20) | 0.8186 | 132.5 | .843 / .798 / .820 | .799 / .772 | 55 | 0.7938 | 143.5 | .886 / .722 | .775 | 50 |
| nullC (median of 20) | 0.7524 | 130 | .777 / .733 / .755 | .789 / .756 | 61.5 | 0.7836 | 144 | .892 / .700 | .769 | 54 |
| nullALL (median of 20) | 0.7528 | 131 | .776 / .735 / .755 | .788 / .764 | 58 | 0.7846 | 146 | .890 / .702 | .777 | 50 |

`cross_pairs` = 0 for every fixed arm (nulls drawing CIGAR-walk edges: median 2), so `P_guard` = P.

**Clauses.**
- B1: B1.1 exact +3 (131 -> 134) PASS; **B1.2 exact null>= 4/20 FAIL** (null range 130-137); B1.3 MICRO R +.012, undetected
  -6 PASS; B1.4 PASS (both arms: `P_guard` +.0003, ARI +.0011, MICRO P -.0024 / -.0028, exact +3 / +3). The other null
  counts, not clauses: MICRO R 2/20, undetected 1/20, ARI 3/20 (median). -> **NOT DISTINGUISHABLE FROM ITS NULL.**
- C: C.1 ARI +.0337, MICRO R +.0136 PASS; C.2 ARI null>= 0/20 PASS (random valid segmentations with C's group counts score
  0.7524, BELOW BASE); C.3 exact +1 PASS; **C.4 FAIL: mean-arm ARI -0.0175** (0.7933 -> 0.7758; `P_guard` +.0138 / +.0047
  and mean exact +5 pass). -> **NOT EFFECTIVE.**
- ALL: NOT EFFECTIVE (B1, C not both EFFECTIVE). For the record: median exact 135, null>= 1/20; ARI 0/20; mean ARI -.0163.
- A: exact Δ 0 / 0, ARI Δ +.0008 / -.0002 -> **HOLDS** (gains ID_211, loses ID_83 in both arms).

**Falsifiers (none fired).** F-B1: 3 of 3 B1 gains (ID_317, ID_438, ID_479) have a member on a B1 edge. F-C: 5 held-out
families are cut by greedy although their clean MAD < 1 (median); C keeps 4 whole (mean arm, not a clause: 2 of 5).
F-A: 5 of 432 held-out cross-family backbone edges (1.2%) touch a readthrough / bridge gene. F-guard: 0 cross-half pairs.

**Predictions vs outcome.** B1.1 (0.65) yes; B1.2 (0.40) no; C.1 (0.70) yes; C.2 (0.45) yes; C EFFECTIVE (0.40) no — via
the mean arm, which I predicted would show "no effect" (it lost 0.0175 ARI and gained 5 exact families).

**Post hoc (after the verdict; not a clause).** The whole held-out mean-arm ARI loss is ONE family: ID_356 (56 genes,
famCN 31.1-42.1). Both the greedy walk and C need 3 groups; C's L1 tie-break moves the upper cut from 41.50|41.54 to
38.56|40.22, so ID_356's recovered pairs fall 951 -> 789 (median arm: C keeps ID_356 whole, +159 pairs). C's median-arm
gain is therefore real on both halves (DEV +.028, held-out +.034, both above nullC), but under the mean statistic the
choice among equally few groups is decided by the tie-break, and this tie-break is not safe. The B1 gains are three of
the five held-out gains the 09-28 union arm already showed (§8.1): B1 selects them, but a random 38 of the same 66 edges
does as well in 4 of 20 draws.

**Full set, labelled DEV + HELD-OUT (not a verdict; `runs/all.tsv`):** BASE 0.7096 / 235 (median), 0.7039 / 261 (mean).
ALL median ARI **0.7402**, exact **247/491**, pair P / R / F1 .836 / .666 / .741, MICRO P / R .797 / .751, undetected
100 (mean: 0.6952 / 276, .896 / .569 / .696, .812 / .747, 86). C median 0.7402 / 239; B1 median 0.7101 / 243.

**What this changes.** Nothing is adopted; the literal recipe (0.7096 / 235) stays the replication headline. What is
now measured: (1) readthrough / fusion bridges are NOT behind the over-merges (DEV and held-out); the over-merges are
>= 98% exon links between families at different famCN, 37% of them between families Soto itself links through a cover
gene; (2) the "SEDEF >= 98% but no minimap2" bucket is not the `-N 50` cap: part is `-p 0.5` against self-paired units,
part is region co-containment with no exon correspondence at all; (3) the greedy famCN cut is beatable by a
minimum-segmentation split in the median statistic (+0.03 ARI on both halves), but the tie-break among minimal splits
must be settled before it can be proposed — a new pre-registration, on a substrate other than these 491 families
(both halves are now spent for famCN-split rules).
