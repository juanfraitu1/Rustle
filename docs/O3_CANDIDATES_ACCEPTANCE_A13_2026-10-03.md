# `o3_candidates` re-run on the 53-family held-out (Amendment 13 + 13b-13e): A13-1 PASSES, A13-2 PASSES, A13-3 PASSES — default-on flip made in 1f49d0f0 and reverted in d04b6ae9: Amendment 14's no-deletion control FAILED (R22) — 2026-10-03

Prereg: `docs/PREREG_rna_allele_haplotype_count_2026-10-01.md`, Amendment 13 (`e3e9d4bf`) with 13b (`d69f0e02`), 13c (`065b2b46`),
13d (`f31f663e`) and 13e (`ff869c40`), each written before the A13 run; Amendment 12's rules (A12-1/2/3) unchanged except A13-1's
comparator (13b, ruling R17). Stage: `o3_candidates` rebuilt from HEAD `0f5824a7` (binary sha1 `0c97f623...`). Recipe:
`bench/rna_allele/accept_o3_candidates.sh` (`ACC=a13`, the default) + `bench/rna_allele/accept_o3_candidates.py`. Work dir
`/mnt/linuxdisk/tmp/rna_allele/a13/` (delta reruns in `half/`, `double/`). Scorer and helper outputs copied to
`docs/O3_CANDIDATES_ACCEPTANCE_A13_score.out.txt`. The previous acceptance: `docs/O3_CANDIDATES_ACCEPTANCE_2026-10-02.md` (A12).

## What changed since A12 (the stage under test)

- **Net (Amendment 13b/13c, rulings R16, R18, R19).** The 31-mer attribution (0 of 5,312 unmapped reads attributed in A12) is retired.
  BAM pass B collects the unmapped records >= 300 bp AND the reads in no net of this run whose primary is poorly placed (`de > 0.02` or
  MAPQ 0), >= 300 bp; they are aligned once (`minimap2 -x map-ont -c -N 5 -p 0.5`) to this run's net reads (`{family}|{read}`) plus every
  family's copy sequences, and a read joins the family of its best hit iff the hit covers >= 50% of the READ and `de <= 0.20`. Under
  `--families` the eligibility and the targets are the run's (R18), so batched runs can attribute one read in two batches.
- **Template (Amendments 13d/13e, rulings R20, R21).** Each cluster is polished on the medoid under a structural distance (indels >= 20 bp
  plus the partner's uncovered ends >= 20 bp, mean over aligned partners; eligible = aligned to >= min(0.5 x (n - 1), 50) others) instead
  of its longest read; votes with the splice preset; insertions >= 20 bp voted before the < 20 bp majority; a kept set whose template
  was split off is re-templated; an empty merged consensus undoes the merge (Amendment 13).
- Unchanged: delta 0.00958, the merge rule, `--min-support 6`, `--min-cluster 3`, the 0.98 tie ratio, the 1,000-read cap, R13.

## Substrate and scoring (Amendment 12's, unchanged)

Amendment 7's 53 families, one copy per family hard-masked (`linktest/masked.fa`); the 59,013 scored reads (17,286 D of the deleted
copies, 41,727 S of the 148 surviving copies) as `linktest/R.bam`; A12's copies table and FASTA (`A12.copies.*`, one row per surviving
copy, linked into `a13/`); the masked splice index; the five batches of A12 (`batches.txt`). Arm M = masked genome + one union per
flagged candidate (`iso_<family>_<k>`), all 59,013 reads realigned with the pipeline's flags (+ `-K 100M`), each candidate its own
component, scored by `merge_test.py score`; contig labels from the best hit in the UNMASKED genome (`GGO.splice.mmi`). No
`RUSTLE_CACHE_DIR` (every run computed). `--threads 4`.

## A13-3: wall time — PASSES

| batch | families | net (used) | `/usr/bin/time` Elapsed | stage's own clock (nets / clusters / genome) | peak RSS | A12 Elapsed |
|---|---|---|---|---|---|---|
| 0 | 6 | 5,674 (3,708) | 283.3 s | 254.5 s (12.7 / 73.2 / 168.0) | 14.5 GB | 149.3 s |
| 1 | 12 | 12,229 (9,464) | 306.6 s | 277.9 s (12.5 / 144.5 / 119.7) | 14.9 GB | 252.0 s |
| 2 | 12 | 9,934 (9,280) | 389.1 s | 349.5 s (16.2 / 260.9 / 72.0) | 15.7 GB | 392.2 s |
| 3 | 11 | 11,288 (9,138) | 305.2 s | 276.4 s (23.9 / 195.0 / 56.4) | 16.0 GB | 288.8 s |
| 4 | 12 | 11,588 (9,228) | 271.8 s | 246.6 s (11.2 / 187.2 / 48.0) | 15.7 GB | 278.7 s |
| **sum** | 53 | 50,713 (40,818) | **1,555.8 s = 25.9 min** | 1,404.9 s = 23.4 min | | 1,360.9 s = 22.7 min |

- **A13-3 PASSES: 25.9 min <= 40 min** (the registered sum of the batches' `/usr/bin/time` Elapsed).
- The nets phase grew from 4-7 s to 11-24 s per batch (the attribution alignment). The genome phase of batches 0 and 1 (168 / 120 s)
  is the 13.6 GB masked index read from a cold page cache; in the delta reruns the same phase took 31-55 s and the five batches summed to
  20.9 and 21.1 min. Clock note as in A12: `/usr/bin/time` (CLOCK_REALTIME) ran ~11% ahead of the stage's monotonic clock; the verdict
  holds on either.

## What the stage produced

**Attribution (pass B), per batch, from the stage logs** (the harness re-runs the same alignment from the BAM and reproduces every count
in every batch, and n_net / n_used of all 53 families: `accept_o3_candidates.py nets`):

| batch | unmapped >= 300 bp | poorly placed >= 300 bp (below the floor) | aligned: unmapped / poorly placed | attributed: unmapped / poorly placed | joined a net of the batch |
|---|---|---|---|---|---|
| 0 | 5,312 | 8,664 (4) | 556 / 6,856 | 21 / 1,144 | 91 |
| 1 | 5,312 | 7,530 (4) | 709 / 5,873 | 14 / 1,381 | 267 |
| 2 | 5,312 | 8,072 (3) | 1,291 / 6,034 | 54 / 1,337 | 112 |
| 3 | 5,312 | 8,759 (4) | 1,474 / 6,924 | 452 / 2,018 | 1,192 |
| 4 | 5,312 | 7,915 (3) | 1,388 / 5,868 | 27 / 775 | 66 |
| sum | 26,560 | 40,940 (18) | 5,418 / 31,555 | 568 / 6,655 | **1,728** |

- **Joined (a read entering a net of its batch): 1,728 joins of 1,667 distinct reads — 557 unmapped (555 = 99.6% the right family by
  `labels.tsv`) and 1,171 poorly placed (927 = 79.2% right).** All 1,482 right-family joiners are D reads of that family (18 families
  receive some); 976 of them fall under the 1,000-read cap. Attributed to a family of another batch (not joined): 5,484 poorly placed reads,
  97.8% to their own family, and 11 unmapped. The poorly placed set is mostly deleted-copy reads: of the 40,940 summed over the
  batches, 37,079 are D reads (36,061 selected by `de > 0.02`, 1,018 by MAPQ 0) and 3,861 S reads.
- From `reads.tsv` alone (as the brief asked): 1,113 of the 1,728 joins are in a reported cluster, 938 of them (84.3%) the right family.
  `reads.tsv` lists only the new-copy and linked clusters' reads, so it cannot show the joins left out by the cap or in in-reference /
  unreported clusters; the harness therefore takes net membership from its validated replication of pass B.
- **Ruling R18 in the counts:** 241 reads sit in nets of two batches, 193 of them attributed in at least one. Of the 246 wrong-family
  joins, 163 are reads that their own family's batch nets by pass A (129 D, 34 S): a single run over the 53 families would not have
  made them eligible.
- **Nets:** 50,713 reads (40,818 after the cap; A12 48,985 / 40,405). D reads in a net of their own family: **8,663** (7,181 with a
  record on a surviving copy + 1,482 attributed; A12 7,181 + 0).
- **Clusters and candidates:** 363 clusters (A12 349): 121 already in the reference, 126 linked, 116 new copy (A12 136 / 146 / 67) ->
  **109 candidates, 82 flagged (>= 6 reads) in 41 families** (A12 60 / 39 in 25). Flagged by label: **30 D-derived, 46 survivor-derived,
  6 elsewhere** (A12 26 / 12 / 1). **Deleted copies with >= 1 D-derived flagged candidate: 25/53** (A12 19/53; IsoCon 44/53 at any
  support, 41/53 at its >= 2-transcript floor); exactly one D-derived candidate in 20 of the 25 (2 in 5). Clusters per flagged candidate:
  1 in 75, 2 in 7; reads per flagged candidate median 158 (6-536).
- **Phase-1 counters (stage log, summed over the 53 families):** 268 read clusters >= `--min-cluster`; 0 empty consensus dropped;
  refinement split off 3,331 reads into 111 new clusters (0 fell under `--min-cluster`); **14 kept sets re-templated** (12 families:
  GWFAM23, 28, 98, 99, 149, 182, 272, 317, 331, 348, 402, 490); significance merge absorbed 16 clusters in 65 rounds, **0 absorptions
  undone**; **fallback templates (no eligible member): 15 the longest member with an aligned partner, 2 the longest member**; every
  other template choice was the structural medoid (or the kept template). 363 final clusters.

## A13-1: read placement — PASSES

| | R (masked) | **o3_candidates A13, arm M** | T (A13's D-derived candidates grouped WITH truth) | o3_candidates A12, arm M | IsoCon, arm M (Amendment 8) |
|---|---|---|---|---|---|
| D right | 0 | **7,898 (45.7%)** | 8,334 | 5,995 (34.7%) | 12,787 (74.0%) |
| D wrong | 10,218 | 4,190 | 4,190 | 4,880 | 3,599 |
| D unplaced | 7,068 | 5,198 | 4,762 | 6,411 | 900 |
| S stay | 41,051 | 37,872 | 37,872 | 40,194 | 40,421 |
| S unplaced | 670 | **3,527** | 3,527 | 1,182 | 1,278 |
| S false moves | 0 | **323 (0.77%)** | 323 | 346 (0.83%) | 25 (0.06%) |
| S elsewhere | 6 | 5 | 5 | 5 | 3 |

**The comparator C (Amendment 13b, ruling R17):** IsoCon's right D reads of Amendment 8 (arm M read by read, `merge_test.py score`
semantics on `linktest/RIL.bam`, `linktest/contigs.tsv`, `linktest/merge/components.tsv`; the replication equals `merge_test.py score`'s
totals) counted over the truth-free attainable D reads:

| D reads | n | IsoCon right | A13 stage right |
|---|---|---|---|
| (1) a record (primary or secondary) on a surviving copy of their family in `R.bam` | 7,181 | **C1 = 6,048** | 5,323 |
| (2) else attributed into their own family's net by this run (555 unmapped, 927 poorly placed) | 1,482 | **C2 = 1,255** | 1,417 |
| attainable, (1) + (2) | 8,663 | **C = 7,303** | 6,740 |
| not attainable | 8,623 | 5,484 | 1,158 |
| all | 17,286 | 12,787 | 7,898 |

- **A13-1 PASSES: D right 7,898 >= 0.80 x C = 5,842.4, and false moves 323 / 41,727 = 0.77% <= 5%.** On the attainable reads alone the
  stage places 6,740 right (92.3% of C), so the verdict does not rest on the registered rule's all-reads numerator. With part (2) read
  from `reads.tsv` alone (the attributed D reads that reached a reported cluster: 938 reads, C2 730, C 6,778, bar 5,422.4) it passes too.
  Counting supplementary records as well in part (1) (a record of any kind on a survivor: 7,286 reads; part (2) 1,416): C = 6,150 +
  1,189 = 7,339, bar 5,871.2 — verdict unchanged. The registered reading is the brief's "any record, primary or secondary" (pass A's
  scope: supplementary records name no read in the stage either).
- Beside, not decided on: **A12-1's original bar (10,230 = 80% of IsoCon's 12,787) is still not met** (7,898 = 77.2% of it).
  Amendment 5's overall rule, R vs arm M: wrong D 10,218 -> 4,190 (-59.0%), false moves 0.77% -> HELP.
- Where the right calls come from (post hoc, `decompose`; the A13 arm M read by read equals `merge_test.py score`): D reads with a record
  on a surviving copy 5,323 right (IsoCon 6,048; A12 5,247), mapped only elsewhere 1,293 (IsoCon 1,501; A12 248), unmapped in the masked
  genome 1,282 (IsoCon 5,238; A12 500). A handful of attributed reads is enough to seed a candidate that the realignment then fills:
  GWFAM21 has 7 attributed D reads in its net and 498 D reads right; GWFAM415 61 and 411.
- False moves: 295 GWFAM175 reads on a candidate labelled S:GWFAM175 (as in A12: probably the real reference-absent copy GWFAM175_B0,
  register row 1206, whose reads the held-out labels S), 19 GWFAM244 reads on its `elsewhere` candidate, 5 + 2 on D-derived candidates
  (GWFAM175, GWFAM269), 1 + 1 on survivor-derived ones (GWFAM331, GWFAM37).
- Per family against A12 (appendix): found now and not in A12 — GWFAM21 (0 -> 498), 23 (0 -> 500), 158 (0 -> 288), 192 (0 -> 110), 415
  (0 -> 411), 440 (0 -> 39), 317 (found, 0 right, IsoCon 0); higher — GWFAM247 158 -> 495, 331 1 -> 136, 490 88 -> 140; **lower —
  GWFAM100 477 -> 110** (its two D-derived candidates now tie for 368 reads, which abstain), GWFAM314 492 -> 433, 425 334 -> 308, 169
  10 -> 0 (its D clusters now link to a locus outside the family), 175 498 -> 495, 268 500 -> 498.

## A13-2: the union representative — PASSES

Measured exactly as A12-2 (`accept_o3_candidates.py keep`: each flagged candidate's reads, `reads.tsv` via `clusters.tsv`, aligned with
`minimap2 -c -x splice:hq -uf -N 10` — the arm-M preset as A12 used it — to the unions and to the cluster consensus sequences; kept = best
AS on its own union >= 0.98 x its best AS over its candidate's consensus sequences). The prereg does not choose between pooled targets
(all unions, all consensus sequences) and isolated targets (per candidate); A12's doc decided on pooled, and so does this one —
**ruling R24 (2026-10-03, after this run): A13-2's pooled reading is the registered one.**

| | pairs measured | kept | lost | kept fraction |
|---|---|---|---|---|
| **all 82 flagged candidates, pooled targets (decided)** | 15,364 (175 not measured) | 15,287 | 77 (1 without a record on the union) | **99.50%** |
| all 82, isolated targets | 15,503 (36 not measured) | 15,469 | 34 | 99.78% |
| 75 single-cluster candidates (union = the consensus), isolated | 13,590 | 13,590 | 0 | 100% |
| 7 two-cluster candidates, isolated | 1,913 | 1,879 | 34 | 98.2% (A12: 68.7%) |

- **A13-2 PASSES: 99.50% >= 95%** (A12: 90.9%). Of A12's four losers, three keep everything now: GWFAM37's deleted copy is a
  single-cluster D-derived candidate (`cand_GWFAM37_1`, 364/364; A12's intron-retaining two-cluster union kept 38/371), `cand_GWFAM269_0`
  (two clusters) keeps 374/374 (A12 338/369), and GWFAM244's D-derived candidate is a single cluster (288/288); `cand_GWFAM331_0` still
  fails (below). The 7 two-cluster unions keep 98.2% of their reads (A12 68.7%).
- **Per-component reading (not the registered verdict, reported because A12's doc reported it):** 2 of 82 candidates keep < 95% with
  isolated targets — `cand_GWFAM440_0` (D-derived, 2 clusters, 18 of 25 reads, 72.0%) and `cand_GWFAM331_0` (D-derived, 2 clusters, 127
  of 142, 89.4%); with pooled targets also `cand_GWFAM104_1` (survivor-derived, one cluster, 341 of 375 = 90.9%; 100% isolated: a single
  cluster's union IS its consensus, so the loss is the pooled targets competing for the `-N 10` slots). Under a strict "each component
  >= 95%" reading A13-2 would fail on those two.

## The deleted copies without a D-derived flagged candidate (28; A12 34), by cause

| cause (registered categories) | n | families |
|---|---|---|
| (i) no D read in the family's net | **23** (A12 25) | GWFAM54, 62, 99, 105, 112, 125, 144, 149, 161, 163, 173, 177, 181, 182, 185, 236, 246, 272, 282, 335, 348, 398, 439 |
| (ii) D reads in the net, in no reported cluster | 1 (A12 6) | GWFAM6 (2 attributed D reads) |
| (iii) the D reads' clusters linked as an allele | 3 (A12 3) | GWFAM28 (333 of 337 used D reads, linked to a survivor; D-to-survivor `de` 0.0057), GWFAM401 (273/276, survivor, 0.0035), GWFAM169 (222/228, linked to a locus outside the family's copies; 0.0095) |
| (iv) not registered: the D reads sit in a flagged candidate labelled survivor-derived | 1 | GWFAM4 (86 of 95; D-to-survivor `de` 0.0057) |

- From A12's 34: GWFAM21 (i) and GWFAM23, 158, 192, 415, 440, 317 (ii) are found now; GWFAM6 moved (i) -> (ii); GWFAM4 (iii) -> (iv);
  GWFAM169, found in A12 (10 right), is (iii) now.
- (i) is still most of the gap: the 23 families' deleted copies have 7,184 D reads (4,000 unmapped in the masked genome, 3,184 mapped
  only off the family's copies), and the attribution brings none of them into the family's net (Amendment 13b's measurement at the
  attribution step found only 135 unmapped + 302 poorly placed D reads attributable over all 53 families: no truth-free rule reaches
  IsoCon's label-scoped net).
- (iii)/(iv): the D reads lie within allele divergence of a survivor (median read `de` 0.35-0.95%, delta 0.958%): the designed
  abstention of Amendment 7's link rule.

## The >= 2-cluster floor (reported beside, ruling R1)

7 candidates (of the 82 flagged) have >= 2 clusters, all 7 D-derived; deleted copies with a D-derived >= 2-cluster candidate: **7/53**
(vs 25/53 at >= 6 reads; A12: 7 candidates, 5/53).

## delta/2 and 2 x delta (reported beside)

Same 53 families, same five batches, the stage re-run with `--delta 0.00479` and `--delta 0.01916`; the nets are identical in the three
runs (the replication checks n_net / n_used and the attribution counts in every batch). The flagged contig sets differ from the
registered run's (34 of 102 and 43 of 68 flagged unions are byte-identical to a registered one), so each run was scored with its own
arm M.

| | delta/2 = 0.00479 | **delta = 0.00958** | 2 x delta = 0.01916 |
|---|---|---|---|
| clusters: in reference / linked / new copy | 141 / 99 / 144 | 121 / 126 / 116 | 112 / 110 / 87 |
| candidates (flagged) | 130 (102) | 109 (82) | 85 (68) |
| flagged by label D / survivor / elsewhere | 38 / 57 / 7 | 30 / 46 / 6 | 27 / 36 / 5 |
| deleted copies with a D-derived flagged candidate | 28 | 25 | 23 |
| D right (arm M) | 8,546 | **7,898** | 7,392 |
| D wrong / unplaced | 3,891 / 4,849 | 4,190 / 5,198 | 4,407 / 5,487 |
| S false moves | 652 (1.56%) | **323 (0.77%)** | 228 (0.55%) |
| S unplaced | 4,061 | 3,527 | 2,658 |
| D right, D-derived candidates grouped with truth (T) | 9,098 | 8,334 | 7,761 |
| A13-1 against the same C (bar 5,842.4) | passes | **passes** | passes |
| >= 2-cluster candidates (D-derived); deleted copies with one | 14 (9); 9 | 7 (7); 7 | 2 (1); 1 |
| deleted copies without a D-derived flag: i / ii / iii / iv | 23 / 1 / 0 / 1 | 23 / 1 / 3 / 1 | 23 / 1 / 4 / 2 |
| stage wall time, 5 batches | 20.9 min | 25.9 min | 21.1 min |

- Cause (i) (23) does not move with delta. Halving delta finds 3 more deleted copies and 648 more D right, at twice the false moves and
  11 more survivor-derived flags; doubling it links more D clusters as alleles (iii: 4) and loses 2 detections. The registered delta sits
  between, and A13-1 holds at all three.

## The default flip: made in commit 1f49d0f0, reverted in d04b6ae9 (Amendment 14 failed, ruling R22)

Commit 1f49d0f0 (2026-10-03, after these three verdicts) flipped `tools/rustle_pipeline.sh` to `CANDIDATES` on by default — `all` runs
`candidates` between `families` and `assign`, and `assign` / `flag` use its products — with `--no-candidates` as the off switch;
`--candidates` (the opt-in switch of ruling R14) stays accepted; `--legacy-catalog` implies `--no-candidates` (an explicit
`--candidates` with it still exits 2); `cand_ready` unchanged. In that commit the driver e2e (`tests/fixtures/o3_candidates/driver/
run_e2e.sh`, on the 0f5824a7 binaries) passes with check (d) rewritten ("a plain assign uses the products": split made, MCL0 a 2-copy
family) and a new control (e) (`assign --no-candidates`: one line says the products are unused, no split, MCL0 1 copy); it also carries
the matching README, AGENTS, REPRODUCE, figures/README + samples.py comments, MODULE_STATUS row, module header and spec (§4, §7, §9,
§9b) lines.

**Ruling R22 (2026-10-03, after this acceptance): the flip does NOT ship until Amendment 14's no-deletion control holds for the stage**
(prereg, last section: Amendment 9's control with the stage in place of the IsoCon chain; C1' = families with >= 1 false flag <= 8 of
53, C2' = false moves without a deletion <= 5% of all reads). That control had not been run when this acceptance was written: if it
held, 1f49d0f0 would ship as committed; if it failed, 1f49d0f0 would be reverted and the stage stay opt-in (R14), this document's
verdicts above standing unchanged. The stage's cost on a full BAM was left to a separate task (ruling R23).

**Outcome (2026-10-03, `docs/O3_CANDIDATES_CONTROL_A14_2026-10-03.md`, register rows 1226-1230): C1' FAILED — 35 of the 53 families
carry a false flag with nothing deleted (bar 8; 54 of the 56 flags match neither of KB3781's haplotypes at 0.999) — and C2' held
(0.54% false moves), so 1f49d0f0 was reverted in d04b6ae9: the stage is opt-in again (R14). The verdicts above stand unchanged. R23:
one batch of 50 families on the full fibroblast BAM did not finish in a 10-minute call (spec §9b).**

**Revert recipe (2026-10-03), if R22 reverts the flip.** Everything in this document and in register rows 1221-1225 is worded to stay
true either way (the choice of this fix round: the true-either-way wording, plus this recipe because the gating commit itself — the one
after 1f49d0f0 that adds this paragraph — edits lines 1f49d0f0 introduced). Run `git revert 1f49d0f0`. Checked on 2026-10-03 in a
scratch worktree on top of the gating commit: it conflicts in exactly three files, five hunks, resolved as follows (and, so resolved,
the pre-flip e2e with its opt-in check (d) passes, 9/9):
- `tools/rustle_pipeline.sh`, 1 hunk — the comment above the `case "$STAGE" in *)` block (the default of `CANDIDATES`): take the
  pre-flip side (aad7aca3's text: the block only refuses `--candidates` with `--legacy-catalog`). The code (`CANDIDATES=0`, the block,
  the `all` dispatch) and the header revert by themselves; the header line "The stage's own cost on a full BAM is not yet measured;
  see Amendment 14 / R23." merges in and stays.
- `README.md`, 1 hunk — the pipeline paragraph: take the pre-flip side (`candidates` OPT-IN, `--candidates`) and re-add its clause
  "its cost on a full BAM is not yet measured, see Amendment 14 / R23".
- `docs/superpowers/specs/2026-10-02-o3-candidates-design.md`, 3 hunks — header and §9: take the pre-flip side and append "The default
  flip (1f49d0f0) was reverted on <date> by ruling R22 (Amendment 14 failed)."; §9b: keep the flip note and R22-R24 as written (they
  record what was decided) and append that dated line.
Everything else 1f49d0f0 changed reverts without conflict (the driver's default and `all` dispatch, the e2e checks, README / AGENTS /
REPRODUCE / figures / MODULE_STATUS / module-header wording). Then run `bash tools/rlock.sh heavy bash
tests/fixtures/o3_candidates/driver/run_e2e.sh --bin <release> --out <scratch>`.
**Applied on 2026-10-03 in d04b6ae9:** the same three files and five hunks, resolved as above except that, the cost having been measured
by then, the README's cost clause and the driver header state the R23 finding (citing "spec §9b, R23") instead of "not yet measured",
and the spec's header and §9 keep the record of the flip with the dated outcome line; the pre-flip e2e passed 9/9.

## Caveats

- **The flip raises the cost on S reads, which no registered rule bounds:** 46 of the 82 flagged candidates are survivor-derived (A12:
  12 of 39; median whole-length d 0.0513, 9 below 0.02), and S reads left unplaced in arm M rise from 1,182 (A12) to 3,527 (8.5% of S;
  IsoCon 1,278): the reads of a surviving copy tie between it and its own survivor-derived candidate and abstain (GWFAM47 491, GWFAM158
  428, GWFAM54 411, GWFAM173 280, GWFAM268 214, ...). They are not false moves (false moves fall to 323), but in O2 these candidates are
  extra copies. Amendment 14 then ran Amendment 9's no-deletion control for the A13 stage (ruling R22): 35 of the 53 families carry a
  false flag with nothing deleted (C1' fails), and 28 of these 46 survivor-derived flags recur there — the flip was reverted
  (d04b6ae9; `docs/O3_CANDIDATES_CONTROL_A14_2026-10-03.md`).
- **One deleted copy, several candidates:** exactly one D-derived candidate in 20 of the 25 found copies; in GWFAM100 the two D-derived
  candidates tie for 368 D reads (D right 477 in A12 -> 110).
- The registered A13-1 compares the stage's D right over ALL reads with C over the attainable reads; on the attainable reads alone the
  stage still clears the bar (6,740 >= 5,842.4).
- Batching (ruling R18): 241 reads sit in nets of two batches, and 163 of the 246 wrong-family joins exist only because a read's own
  family ran in another batch. A single run over all families would attribute differently (not measured here).
- A13-2's per-component reading fails on 2 of 82 candidates (above); the decided aggregate passes with margin.
- The harness change: `accept_o3_candidates.py nets` no longer re-implements the retired k-mer rule; it re-runs the A13 pass B from the BAM
  and is checked against the stage's own logs, families.tsv, reads.tsv and nets.fa (all equal in 5 batches x 3 runs). Part (2) of C
  uses those nets; the `reads.tsv`-only reading is reported beside and passes too.
- One individual (KB3781), fibroblast Iso-Seq, the 1,000-read cap; the registered run's genome phase started on a cold page cache (above).

- **A13 is not a held-out of its own rules** (final review, 2026-10-03): Amendment 13b's preset and thresholds were chosen from
  truth-labelled counts on these same reads; 13d/13e were fixed after smoke runs on these families (GWFAM105); R24, the reading under
  which A13-2 passes, was ruled after the run (the per-component reading fails on 2 of 82). The held-out of the corrected stage is
  Amendment 15's 30 families of Amendment 10 disjoint from these 53.

## Appendix: per family (registered run)

Flagged labels: D = D-derived, S = survivor-derived, e = elsewhere. "D reads in net (attributed)": D reads of the family in its net, and
how many came by attribution. D right from `decompose` (A13 / A12 / IsoCon's Amendment 8 arm M). Cause as in the table above; "found" =
a D-derived flagged candidate exists.

| family | batch | net (used) | D reads in net (attributed) | clusters: ref / linked / new | candidates (flagged) | flagged labels | D right: A13 / A12 / IsoCon | cause |
|---|---|---|---|---|---|---|---|---|
| GWFAM4 | 0 | 400 (400) | 95 (0) | 0 / 2 / 1 | 1 (1) | S | 0 / 0 / 9 | iv |
| GWFAM6 | 1 | 1023 (1000) | 2 (2) | 2 / 0 / 3 | 3 (2) | SS | 0 / 0 / 34 | ii |
| GWFAM21 | 2 | 1007 (1000) | 7 (7) | 2 / 1 / 1 | 1 (1) | D | 498 / 0 / 500 | found |
| GWFAM23 | 3 | 1995 (1000) | 495 (494) | 2 / 2 / 4 | 4 (3) | DSS | 500 / 0 / 478 | found |
| GWFAM28 | 0 | 1500 (1000) | 500 (0) | 0 / 4 / 0 | 0 (0) | - | 0 / 0 / 30 | iii |
| GWFAM37 | 4 | 1315 (1000) | 499 (0) | 0 / 6 / 2 | 2 (2) | DS | 496 / 496 / 496 | found |
| GWFAM47 | 1 | 2056 (1000) | 104 (11) | 5 / 7 / 6 | 6 (6) | DSSSSe | 107 / 107 / 107 | found |
| GWFAM54 | 2 | 1001 (1000) | 0 (0) | 3 / 0 / 1 | 1 (1) | S | 0 / 0 / 80 | i |
| GWFAM62 | 3 | 1000 (1000) | 0 (0) | 3 / 1 / 1 | 1 (0) | - | 0 / 0 / 433 | i |
| GWFAM98 | 4 | 1489 (1000) | 488 (1) | 2 / 1 / 4 | 3 (3) | DDS | 488 / 488 / 494 | found |
| GWFAM99 | 1 | 1000 (1000) | 0 (0) | 1 / 2 / 2 | 2 (2) | SS | 0 / 0 / 0 | i |
| GWFAM100 | 0 | 1311 (1000) | 469 (30) | 3 / 4 / 3 | 3 (2) | DD | 110 / 477 / 478 | found |
| GWFAM104 | 2 | 1298 (1000) | 486 (10) | 0 / 2 / 2 | 2 (2) | DS | 417 / 417 / 417 | found |
| GWFAM105 | 0 | 83 (83) | 0 (0) | 4 / 0 / 0 | 0 (0) | - | 0 / 0 / 496 | i |
| GWFAM112 | 4 | 521 (521) | 0 (0) | 3 / 1 / 0 | 0 (0) | - | 0 / 0 / 268 | i |
| GWFAM123 | 4 | 596 (596) | 46 (0) | 3 / 1 / 3 | 3 (1) | D | 39 / 39 / 30 | found |
| GWFAM125 | 2 | 680 (680) | 0 (0) | 2 / 1 / 1 | 1 (1) | S | 0 / 0 / 310 | i |
| GWFAM144 | 3 | 1000 (1000) | 0 (0) | 3 / 2 / 0 | 0 (0) | - | 0 / 0 / 116 | i |
| GWFAM149 | 1 | 234 (234) | 0 (0) | 4 / 1 / 0 | 0 (0) | - | 0 / 0 / 423 | i |
| GWFAM158 | 4 | 1474 (1000) | 26 (25) | 7 / 9 / 4 | 4 (4) | DSSe | 288 / 0 / 290 | found |
| GWFAM161 | 4 | 300 (300) | 0 (0) | 3 / 4 / 0 | 0 (0) | - | 0 / 0 / 0 | i |
| GWFAM163 | 3 | 687 (687) | 0 (0) | 0 / 6 / 3 | 3 (2) | SS | 0 / 0 / 476 | i |
| GWFAM164 | 0 | 2155 (1000) | 288 (0) | 1 / 4 / 4 | 4 (2) | DS | 287 / 287 / 287 | found |
| GWFAM169 | 1 | 679 (679) | 228 (1) | 3 / 3 / 0 | 0 (0) | - | 0 / 10 / 40 | iii |
| GWFAM173 | 3 | 815 (815) | 0 (0) | 6 / 1 / 2 | 2 (2) | SS | 0 / 0 / 240 | i |
| GWFAM175 | 1 | 2463 (1000) | 500 (0) | 4 / 1 / 8 | 8 (4) | DDSS | 495 / 498 / 499 | found |
| GWFAM177 | 1 | 610 (610) | 0 (0) | 2 / 0 / 1 | 1 (1) | S | 0 / 0 / 0 | i |
| GWFAM181 | 2 | 868 (868) | 0 (0) | 1 / 0 / 1 | 1 (1) | S | 0 / 0 / 71 | i |
| GWFAM182 | 4 | 287 (287) | 0 (0) | 3 / 0 / 0 | 0 (0) | - | 0 / 0 / 40 | i |
| GWFAM185 | 1 | 383 (383) | 0 (0) | 2 / 2 / 1 | 1 (1) | S | 0 / 0 / 0 | i |
| GWFAM192 | 1 | 642 (642) | 108 (107) | 4 / 3 / 1 | 1 (1) | D | 110 / 0 / 0 | found |
| GWFAM236 | 2 | 1012 (1000) | 0 (0) | 3 / 2 / 2 | 2 (0) | - | 0 / 0 / 329 | i |
| GWFAM244 | 3 | 1697 (1000) | 500 (436) | 3 / 2 / 5 | 5 (3) | DSe | 500 / 500 / 500 | found |
| GWFAM246 | 4 | 1000 (1000) | 0 (0) | 2 / 1 / 3 | 3 (1) | S | 0 / 0 / 498 | i |
| GWFAM247 | 1 | 1223 (1000) | 500 (2) | 0 / 3 / 2 | 2 (2) | DD | 495 / 158 / 497 | found |
| GWFAM268 | 4 | 795 (795) | 500 (0) | 0 / 1 / 5 | 4 (2) | DS | 498 / 500 / 498 | found |
| GWFAM269 | 2 | 1281 (1000) | 493 (19) | 2 / 4 / 3 | 2 (1) | D | 491 / 491 / 491 | found |
| GWFAM272 | 2 | 555 (555) | 0 (0) | 2 / 1 / 1 | 1 (1) | e | 0 / 0 / 0 | i |
| GWFAM282 | 3 | 1000 (1000) | 0 (0) | 1 / 2 / 1 | 1 (1) | S | 0 / 0 / 500 | i |
| GWFAM314 | 4 | 2000 (1000) | 500 (0) | 2 / 2 / 2 | 1 (1) | D | 433 / 492 / 492 | found |
| GWFAM317 | 3 | 548 (548) | 41 (38) | 5 / 4 / 4 | 4 (3) | DSS | 0 / 0 / 0 | found |
| GWFAM331 | 2 | 422 (422) | 152 (1) | 0 / 0 / 5 | 4 (3) | DSS | 136 / 1 / 103 | found |
| GWFAM335 | 3 | 606 (606) | 0 (0) | 4 / 3 / 4 | 4 (3) | SSe | 0 / 0 / 162 | i |
| GWFAM348 | 4 | 729 (729) | 0 (0) | 2 / 3 / 5 | 5 (4) | SSSe | 0 / 0 / 200 | i |
| GWFAM398 | 1 | 1000 (1000) | 0 (0) | 1 / 3 / 0 | 0 (0) | - | 0 / 0 / 0 | i |
| GWFAM401 | 2 | 1055 (1000) | 290 (0) | 2 / 2 / 0 | 0 (0) | - | 0 / 0 / 172 | iii |
| GWFAM402 | 2 | 526 (526) | 134 (0) | 3 / 2 / 2 | 2 (2) | DS | 118 / 118 / 111 | found |
| GWFAM407 | 3 | 1458 (1000) | 491 (212) | 3 / 6 / 2 | 1 (1) | D | 494 / 494 / 485 | found |
| GWFAM415 | 0 | 225 (225) | 119 (61) | 1 / 2 / 1 | 1 (1) | D | 411 / 0 / 53 | found |
| GWFAM425 | 4 | 1082 (1000) | 434 (0) | 1 / 2 / 3 | 3 (3) | DDS | 308 / 334 / 375 | found |
| GWFAM439 | 3 | 482 (482) | 0 (0) | 3 / 6 / 1 | 1 (1) | S | 0 / 0 / 0 | i |
| GWFAM440 | 2 | 229 (229) | 26 (25) | 2 / 0 / 3 | 2 (2) | DS | 39 / 0 / 39 | found |
| GWFAM490 | 1 | 916 (916) | 142 (0) | 1 / 4 / 3 | 3 (2) | DS | 140 / 88 / 140 | found |

Register rows 1221-1225.
