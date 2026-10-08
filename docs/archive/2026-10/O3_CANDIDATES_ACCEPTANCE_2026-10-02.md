# `o3_candidates` on the 53-family held-out (Amendment 12): A12-1 FAILS, A12-2 FAILS, A12-3 PASSES — 2026-10-02

Prereg: `docs/PREREG_rna_allele_haplotype_count_2026-10-01.md` Amendment 12 (commit 0458f928, written before the stage existed; flag floor
= ruling R1, `--min-support 6`). Stage: `o3_candidates` at HEAD a3564999 (rebuilt before the run). Recipe:
`bench/rna_allele/accept_o3_candidates.sh` (steps) + `bench/rna_allele/accept_o3_candidates.py` (helpers) +
`bench/rna_allele/panel_to_copies.py` (the copies table). Work dir `/mnt/linuxdisk/tmp/rna_allele/a12/` (delta reruns in `half/`,
`double/`). Scorer output copied to `docs/O3_CANDIDATES_ACCEPTANCE_score.out.txt`.

## Substrate and inputs (Amendments 7-8, unchanged)

- Amendment 7's 53 families (>= 3 copies), one copy per family hard-masked to N (`linktest/masked.fa`); the 59,013 scored reads
  (17,286 of the deleted copies = D, 41,727 of the 148 surviving copies = S) aligned to the masked genome with the pipeline's flags
  (`linktest/R.bam`); the masked splice index.
- Copies table `A12.copies.tsv` (`P.fam.copies.tsv` layout) from `panel.json`: one row per SURVIVING copy (148), `start-end` = the clean
  interval as one exon block, `n_reads` = primaries of `R.bam` on it (none is 0), `source panel`; `A12.copies.fa` = the interval from the
  masked genome. The deleted copies are not in the table: the stage must find them from the reads alone.
- Stage defaults: `--delta 0.00958 --max-reads 1000 --min-cluster 3 --min-support 6 --threads 4`; no `RUSTLE_CACHE_DIR` (every run
  computed).

## A12-3: wall time — PASSES

The stage ran in 5 foreground batches (`--families`, one call each, < 10 min; `batches.txt`). Batching changes no result: nets, clusters
and candidates are per family, the unmapped-read index holds every family's copies in every batch, the genome hits are per consensus. The
Python recomputation of the nets (`accept_o3_candidates.py nets`) equals the stage's `n_net` / `n_used` in 53/53 families.

| batch | families | used reads | `/usr/bin/time` Elapsed | stage's own clock (nets / clusters / genome) | peak RSS |
|---|---|---|---|---|---|
| 0 | 6 | 3,647 | 149.3 s | 134.8 s (3.9 / 69.2 / 61.8) | 15.7 GB |
| 1 | 12 | 9,228 | 252.0 s | 226.8 s (4.1 / 132.7 / 89.6) | 15.7 GB |
| 2 | 12 | 9,236 | 392.2 s | 348.5 s (6.9 / 253.8 / 83.6) | 14.9 GB |
| 3 | 11 | 9,088 | 288.8 s | 260.0 s (4.5 / 191.5 / 63.8) | 16.3 GB |
| 4 | 12 | 9,206 | 278.7 s | 249.6 s (5.4 / 200.3 / 42.9) | 15.7 GB |
| **sum** | 53 | 40,405 | **1,360.9 s = 22.7 min** | 1,219.7 s = 20.3 min | |

- **A12-3 PASSES: 22.7 min <= 40 min** (the registered sum of the batches; 2 x IsoCon's ~20 min).
- Clock note: `/usr/bin/time` reports CLOCK_REALTIME, which on this WSL2 VM advanced ~11% faster than the monotonic clock during these
  runs (batch 3: 288.8 s realtime vs 260.1 s `/proc/uptime`; the stage's own `Instant` timer agrees with the monotonic one). The verdict
  holds on either clock. Batching costs one genome-phase index load (13.6 GB `.mmi`, 43-90 s) and one BAM pass per batch; a single run
  would take ~16 min (monotonic sum minus 4 redundant genome and nets phases).

## What the stage produced

- Nets: 48,985 reads over the 53 families (40,405 after the 1,000-read cap). **Pass B attributed 0 of the 5,312 unmapped reads >= 300 bp**
  (all of them D reads) to any family: they share essentially no 31-mer with any surviving copy. 1,293 of the 5,312 share >= 1 indexed
  31-mer with some family (691 with their own); the fraction of a read's 31-mers found in its own family's copies is p99 0.008, max
  0.097 (counting the repeat k-mers the index drops; 0.006 / 0.044 without them), the best family's p90 0.002, max 0.076, against the
  0.30 floor (final review's recount, reproduced 2026-10-03; a sample of 300 mapped S reads: best-family median 0.85, 286/300 >= 0.30).
- 349 clusters: 136 already in the reference (identity x coverage >= 0.999), 146 linked (d <= delta), 67 new copy -> **60 candidates, 39
  flagged** (>= 6 reads) in 25 families. Flagged by label (best unmasked hit): **26 D-derived, 12 survivor-derived, 1 elsewhere**.
- **Deleted copies with >= 1 D-derived flagged candidate: 19/53** (IsoCon, Amendment 8: 44/53 — but at a different floor: IsoCon's 44
  counts a candidate of ANY support, >= 1 IsoCon transcript of >= 2 reads; at its >= 2-transcript floor, the floor R1 translated to
  >= 6 reads, IsoCon finds 41/53, `docs/RNA_ALLELE_CONTROL_2026-10-01.md`). Exactly one D-derived candidate in 13 of
  the 19 (2 in 5, 3 in 1). Candidates per family: 0 in 23 families, 1 in 15, 2 in 7, 3 in 4, 4 in 2, 5 and 6 in one each; clusters per
  flagged candidate: 1 in 32, 2 in 7; reads per flagged candidate median 107 (6-390).

## A12-1: read placement — FAILS

Arm M = masked genome + the 39 flagged unions (`iso_<family>_<k>` = `cand_<family>_<k>`), all 59,013 reads realigned with the R arm's
command (+ `-K 100M`), each candidate its own component (its union already is the merge), so M == RIL.

| | R (masked) | **o3_candidates, arm M** | T (D-derived candidates of a family grouped WITH truth) | IsoCon, arm M (Amendment 8) |
|---|---|---|---|---|
| D right | 0 | **5,995 (34.7%)** | 6,420 | 12,787 (74.0%) |
| D wrong | 10,218 | 4,880 | 4,880 | 3,599 |
| D unplaced | 7,068 | 6,411 | 5,986 | 900 |
| S stay | 41,051 | 40,194 | 40,194 | 40,421 |
| S unplaced | 670 | 1,182 | 1,182 | 1,278 |
| S false moves | 0 | **346 (0.83%)** | 346 | 25 (0.06%) |
| S elsewhere | 6 | 5 | 5 | 3 |

- **A12-1 FAILS: D right 5,995 against the bar of 10,230 (80% of IsoCon's 12,787)**; false moves 346 / 41,727 = 0.83% pass (bar 5%).
- Reported beside: Amendment 5's overall rule, R vs arm M: wrong D 10,218 -> 4,880 (-52.2%), false moves 0.83% -> HELP. (`merge_test.py`
  also prints Amendment 8's M1 / M2 lines; they are not Amendment 12's rules — with one component per candidate, M2's "D contigs in one
  component" is 0/6 by construction.)
- False moves: 295 are GWFAM175 reads placed on `iso_GWFAM175_1`, a candidate labelled as derived from another surviving copy
  (GWFAM175:2); 36 GWFAM331 reads likewise; 8 GWFAM244 reads on its `elsewhere` candidate; 7 on D-derived candidates.
  ⚠ The 295 are probably not false: GWFAM175 carries a real reference-absent copy, GWFAM175_B0 (Amendment 10, row 1206: 281 reads
  whose primaries sit on `_pri` copy GWFAM175:2), and the held-out labels a read by the copy of its primary alignment on the unmasked
  reference (Amendment 7), so B0's reads are labelled S of GWFAM175:2. `iso_GWFAM175_1` would then be B0 and these reads correctly
  placed — unverified (the candidate was not aligned to the haplotype assemblies here); the verdict does not depend on it (false moves
  pass either way).

## Why: the deleted copies without a D-derived flagged candidate (34), by cause

| cause (registered categories) | n | families |
|---|---|---|
| (i) no D read in the family's net | **25** | GWFAM6, 21, 54, 62, 99, 105, 112, 125, 144, 149, 161, 163, 173, 177, 181, 182, 185, 236, 246, 272, 282, 335, 348, 398, 439 |
| (ii) D reads in the net, but below the floors (< 3-read clusters / a 3-read component) | 5 | GWFAM23, 158, 192, 440 (1 D read each); GWFAM317 (3 D reads, one unflagged 3-read candidate) |
| (ii) D reads in no reported cluster (below `--min-cluster`, split off, or in an in-reference cluster) | 1 | GWFAM415 (58 D reads: 38 in no reported cluster, 20 linked; D-to-survivor `de` 0.0052) |
| (iii) D reads' clusters linked to a survivor | 3 | GWFAM4 (93 of 95 D reads), GWFAM28 (332/337), GWFAM401 (281/284); D-to-survivor `de` 0.0057 / 0.0057 / 0.0035 |

- The 34 = 25 (i) + 5 (ii, below the floors) + 1 (ii, GWFAM415: no reported cluster) + 3 (iii).
- **(i) is most of the gap.** The 25 deleted copies' reads are unmapped in the masked genome (4,498 reads) or map only to loci outside
  the family's copy list (3,225); no record of theirs touches a surviving copy, and the 31-mer attribution gives none of the unmapped
  ones to a family. Over all 53 families: 7,181 D reads have a record on a surviving copy (in a net), 4,792 map only elsewhere,
  5,313 are unmapped.
- (iii) and GWFAM415: the D reads sit within allele divergence of a survivor (median read `de` 0.35-0.57%, below delta 0.958%), so their
  consensus links as an allele — the designed abstention (Amendment 7's link rule). IsoCon's transcripts in these families still cleared
  delta (D right 9 / 30 / 172 / 53 for GWFAM4 / 28 / 401 / 415).
- (ii) 1-3 D reads in the net cannot make a 3-read cluster or a 6-read candidate.

**Post hoc, not a registered rule (`accept_o3_candidates.py decompose`; it reproduces `merge_test.py score`'s arm M read by read for
both runs, totals equal):**

| D reads, by where the R arm puts them | n | IsoCon right | o3_candidates right |
|---|---|---|---|
| a record on a surviving copy of the family (what a truth-free net sees) | 7,181 | 6,048 | **5,247 (86.8% of IsoCon's)** |
| mapped only elsewhere | 4,792 | 1,501 | 248 |
| unmapped in the masked genome | 5,313 | 5,238 | 500 |

IsoCon's input net (Amendment 7: "the reads with an R-arm record overlapping a surviving copy, plus the reads R leaves unmapped",
within the family's scored reads) gave each family its unmapped reads **by the truth label**; 5,238 of IsoCon's 12,787 right calls are
reads the R arm leaves unmapped (the stage: 500). The stage has to attribute unmapped reads by sequence, and these reads align nowhere in
the masked genome and share essentially no 31-mer with their family's surviving copies (own-family 31-mer fraction p99 0.008, max 0.097;
above); 4,498 of the 5,313 belong to the 25 cause-(i) families,
where no read of the deleted copy touches the family's copies at all. Where both runs find the deleted copy, the counts (IsoCon -> stage)
are near-identical (GWFAM37 496 -> 496, GWFAM100 478 -> 477, GWFAM104 417 -> 417, GWFAM164 287 -> 287, GWFAM175 499 -> 498, GWFAM244
500 -> 500, GWFAM268 498 -> 500, GWFAM269 491 -> 491, GWFAM314 492 -> 492, GWFAM407 485 -> 494); lower in GWFAM247 (497 -> 158), GWFAM331
(103 -> 1), GWFAM425 (375 -> 334), GWFAM490 (140 -> 88), GWFAM169 (40 -> 10). The bar was set on a comparator that held those unmapped reads by label; the registered verdict stands
as FAILS.

## A12-2: the union representative — FAILS

Each flagged candidate's reads (`reads.tsv` rows of its clusters; 5,416 read-candidate pairs) aligned with `minimap2 -c -x splice:hq -uf
-N 10` to the unions and to the cluster consensus sequences (`clusters.fa`); kept = best AS on the read's own union >= 0.98 x its best AS
over its candidate's consensus sequences; no record on the union = not kept; no record on any of its consensus sequences = not measured.
Pooled targets (all unions, all of `clusters.fa`: the registered verdict) and, as a check, isolated targets (per candidate its union
alone and its own consensus sequences alone, so that no other target competes for the `-N 10` slots).

| | pairs measured | kept | lost | kept fraction |
|---|---|---|---|---|
| **all flagged candidates, pooled** | 5,352 (64 not measured) | 4,866 | 486 | **90.9%** |
| all flagged candidates, isolated | 5,412 (4 not measured) | 4,929 | 483 | 91.1% |
| 32 single-cluster candidates (union = the consensus), isolated | 3,869 | 3,869 | 0 | 100% |
| 7 two-cluster candidates, isolated | 1,543 | 1,060 | 483 | 68.7% |

- **A12-2 FAILS: 90.9% < 95%** (91.1% with isolated targets; 4 of the 39 candidates are below 95% on their own, so the per-component
  reading fails too). A single-cluster candidate passes by construction: its union is its consensus.
- The four failing candidates, all two-cluster. Three fail because the union carries segments that most of its reads lack, which those
  reads then cross as introns (a splice penalty > 2% of their AS). **`cand_GWFAM37_0` (D-derived; 38/371 kept)**: the main cluster's
  consensus is seeded on the cluster's longest read, an intron-retaining one, and ruling R2 never applies a >= 20 bp deletion, so the
  consensus (6,762 bp) keeps four introns (655 / 630 / 1,402 / 1,095 bp) that the second cluster's consensus (3,023 bp, spliced) lacks;
  the union is that backbone, and 340 of the main cluster's own 347 reads align better to the spliced consensus, 314 of them by more than
  2%. **`cand_GWFAM269_0` (D; 338/369)**: a 593-bp segment of the backbone consensus that 27 of its own reads skip.
  **`cand_GWFAM244_0` (survivor-derived; 88/100)**: two segments (479 and 136 bp) inserted from the 3-read minority cluster, skipped by
  11 of the main cluster's 97 reads. **`cand_GWFAM331_0` (survivor-derived; 150/257 kept with isolated targets, 144 pooled)** fails the
  other way: two clusters of 7,413 / 7,410 bp joined into one component; the union is one of them and lacks the other's substitutions and
  a 3-bp junction shift, so 105 of the second cluster's 128 reads score lower on the union (10-13% of AS in the inspected cases).
- By label: D-derived candidates keep 3,964 / 4,328 (91.6%, isolated), survivor-derived 959 / 1,078.

## The >= 2-cluster floor (reported beside, ruling R1)

Every final cluster holds >= `--min-cluster` 3 reads, so a >= 2-cluster candidate always has >= 6 reads: the >= 2-cluster flags are a
subset of the registered ones. **7 candidates** (of the 39 flagged) have >= 2 clusters: 5 D-derived, 2 survivor-derived; deleted copies
with a D-derived >= 2-cluster candidate: **5/53** (vs 19/53 at >= 6 reads). On this substrate the stage's clusters already merge a copy's
isoforms, so the cluster-count floor would cut detection by three quarters.

## delta/2 and 2 x delta (reported beside)

Same 53 families, same nets (delta does not enter them), the stage re-run with `--delta 0.00479` and `--delta 0.01916` (9 batches each;
`half/`, `double/`), each with its own arm M (the flagged contig sets differ: 17 of 48 and 18 of 37 flagged unions are byte-identical to
a registered one), scored the same way.

| | delta/2 = 0.00479 | **delta = 0.00958** | 2 x delta = 0.01916 |
|---|---|---|---|
| clusters: in reference / linked / new copy | 158 / 112 / 89 | 136 / 146 / 67 | 126 / 128 / 50 |
| candidates (flagged) | 79 (48) | 60 (39) | 49 (37) |
| flagged by label D / survivor / elsewhere | 23 / 25 / 0 | 26 / 12 / 1 | 22 / 13 / 2 |
| deleted copies with a D-derived flagged candidate | 19 | 19 | 15 |
| D right (arm M) | 5,944 | **5,995** | 5,021 |
| D wrong / unplaced | 4,921 / 6,421 | 4,880 / 6,411 | 5,318 / 6,947 |
| S false moves | 312 (0.75%) | **346 (0.83%)** | 272 (0.65%) |
| D right, D-derived candidates grouped with truth (T) | 5,985 | 6,420 | 5,394 |
| >= 2-cluster candidates (D-derived); deleted copies with one | 10 (7); 7 | 7 (5); 5 | 1 (1); 1 |
| deleted copies without a D-derived flag: i / ii / iii / iv | 25 / 5 / 3 / 1 | 25 / 6 / 3 / 0 | 25 / 5 / 6 / 2 |

(iv), not a registered category: the plurality of the family's D reads under the 1,000-read cap sit in a FLAGGED candidate whose best
unmasked hit labels it survivor-derived or elsewhere (`accept_o3_candidates.py report` takes each family's plurality fate of those
reads; a tie resolves in the order (iv), (ii) below `--min-support`, (iii), (ii) in no reported cluster).

- No delta rescues A12-1 (bar 10,230): cause (i), 25 deleted copies with no read in the net, does not depend on delta. Halving delta
  doubles the survivor-derived flags (25 vs 12: alleles beyond delta/2 become "new copies") without finding more deleted copies (19:
  GWFAM28 gained, GWFAM104 lost); doubling it links three more deleted copies to a survivor (iii: 6) and loses 4 detections. The
  registered delta is at the top of the range for D right and found copies.

## Caveats

- One individual (the reference animal KB3781), fibroblast Iso-Seq, the 1,000-read cap per family (as IsoCon's chain).
- R11: the all-vs-all uses `--dual=no`, so each read keeps <= 100 hits to later-sorting reads; in a multi-copy net cross-copy hits can
  crowd out same-copy partners, and with `-p 0.1` relative to the query's self hit a pair > 10x apart in length gets no hit when the longer
  read sorts first. Not isolated here: in the 19 families where the stage finds the deleted copy, D right equals IsoCon's in most (above);
  the 25 cause-(i) families have no D read to cluster; its share in the other 9 is not measured.
- Pass A also admits spliced host-gene reads whose intron spans a copy; they dilute the 1,000-read cap (in-reference clusters, no flag).
- The comparator (IsoCon's 12,787) was computed with truth-scoped nets for unmapped reads (above); the registered bar was set against it.
- `contigs.tsv`'s `best_masked` is NA (the candidates were not re-aligned to the masked genome; `d` is the candidate's own, from the stage).

## Appendix: per family (registered run)

Flagged labels: D = D-derived, S = survivor-derived, e = elsewhere. D right from `accept_o3_candidates.py decompose` (IsoCon = Amendment 8's
arm M). Cause as in the table above; "found" = a D-derived flagged candidate exists.

| family | net (used) | D reads in net | clusters: ref / linked / new | candidates (flagged) | flagged labels | D right: stage / IsoCon | cause |
|---|---|---|---|---|---|---|---|
| GWFAM4 | 400 (400) | 95 | 0 / 3 / 0 | 0 (0) | - | 0 / 9 | iii |
| GWFAM6 | 1021 (1000) | 0 | 3 / 1 / 1 | 1 (0) | - | 0 / 34 | i |
| GWFAM21 | 1000 (1000) | 0 | 2 / 1 / 0 | 0 (0) | - | 0 / 500 | i |
| GWFAM23 | 1501 (1000) | 1 | 2 / 3 / 1 | 1 (0) | - | 0 / 478 | ii |
| GWFAM28 | 1500 (1000) | 500 | 0 / 3 / 0 | 0 (0) | - | 0 / 30 | iii |
| GWFAM37 | 1314 (1000) | 499 | 1 / 3 / 2 | 1 (1) | D | 496 / 496 | found |
| GWFAM47 | 2030 (1000) | 93 | 8 / 7 / 1 | 1 (1) | D | 107 / 107 | found |
| GWFAM54 | 1000 (1000) | 0 | 3 / 1 / 0 | 0 (0) | - | 0 / 80 | i |
| GWFAM62 | 1000 (1000) | 0 | 1 / 3 / 1 | 1 (0) | - | 0 / 433 | i |
| GWFAM98 | 1488 (1000) | 487 | 1 / 2 / 3 | 3 (2) | DS | 488 / 494 | found |
| GWFAM99 | 1000 (1000) | 0 | 1 / 3 / 1 | 1 (1) | S | 0 / 0 | i |
| GWFAM100 | 1281 (1000) | 439 | 3 / 5 / 3 | 3 (2) | DD | 477 / 478 | found |
| GWFAM104 | 1288 (1000) | 476 | 2 / 3 / 1 | 1 (1) | D | 417 / 417 | found |
| GWFAM105 | 83 (83) | 0 | 3 / 2 / 0 | 0 (0) | - | 0 / 496 | i |
| GWFAM112 | 521 (521) | 0 | 2 / 2 / 0 | 0 (0) | - | 0 / 268 | i |
| GWFAM123 | 596 (596) | 46 | 1 / 3 / 2 | 2 (1) | D | 39 / 30 | found |
| GWFAM125 | 679 (679) | 0 | 2 / 2 / 0 | 0 (0) | - | 0 / 310 | i |
| GWFAM144 | 1000 (1000) | 0 | 3 / 2 / 0 | 0 (0) | - | 0 / 116 | i |
| GWFAM149 | 234 (234) | 0 | 3 / 1 / 0 | 0 (0) | - | 0 / 423 | i |
| GWFAM158 | 1432 (1000) | 1 | 10 / 7 / 0 | 0 (0) | - | 0 / 290 | ii |
| GWFAM161 | 300 (300) | 0 | 3 / 3 / 0 | 0 (0) | - | 0 / 0 | i |
| GWFAM163 | 687 (687) | 0 | 0 / 5 / 3 | 3 (2) | SS | 0 / 476 | i |
| GWFAM164 | 2155 (1000) | 288 | 3 / 4 / 2 | 2 (1) | D | 287 / 287 | found |
| GWFAM169 | 550 (550) | 227 | 2 / 3 / 1 | 1 (1) | D | 10 / 40 | found |
| GWFAM173 | 815 (815) | 0 | 7 / 1 / 1 | 1 (1) | S | 0 / 240 | i |
| GWFAM175 | 2463 (1000) | 500 | 5 / 2 / 7 | 6 (4) | DDDS | 498 / 499 | found |
| GWFAM177 | 610 (610) | 0 | 2 / 1 / 0 | 0 (0) | - | 0 / 0 | i |
| GWFAM181 | 868 (868) | 0 | 4 / 0 / 0 | 0 (0) | - | 0 / 71 | i |
| GWFAM182 | 287 (287) | 0 | 3 / 1 / 0 | 0 (0) | - | 0 / 40 | i |
| GWFAM185 | 383 (383) | 0 | 2 / 4 / 0 | 0 (0) | - | 0 / 0 | i |
| GWFAM192 | 535 (535) | 1 | 4 / 2 / 0 | 0 (0) | - | 0 / 0 | ii |
| GWFAM236 | 1000 (1000) | 0 | 3 / 3 / 0 | 0 (0) | - | 0 / 329 | i |
| GWFAM244 | 1261 (1000) | 64 | 2 / 4 / 6 | 5 (3) | DSe | 500 / 500 | found |
| GWFAM246 | 1000 (1000) | 0 | 2 / 1 / 2 | 2 (0) | - | 0 / 498 | i |
| GWFAM247 | 1220 (1000) | 498 | 0 / 3 / 2 | 2 (2) | DD | 158 / 497 | found |
| GWFAM268 | 795 (795) | 500 | 1 / 2 / 4 | 4 (2) | DD | 500 / 498 | found |
| GWFAM269 | 1260 (1000) | 474 | 2 / 3 / 2 | 1 (1) | D | 491 / 491 | found |
| GWFAM272 | 540 (540) | 0 | 3 / 0 / 0 | 0 (0) | - | 0 / 0 | i |
| GWFAM282 | 1000 (1000) | 0 | 1 / 3 / 0 | 0 (0) | - | 0 / 500 | i |
| GWFAM314 | 2000 (1000) | 500 | 2 / 4 / 2 | 1 (1) | D | 492 / 492 | found |
| GWFAM317 | 506 (506) | 3 | 2 / 7 / 2 | 2 (1) | S | 0 / 0 | ii |
| GWFAM331 | 421 (421) | 151 | 0 / 0 / 5 | 4 (3) | DSS | 1 / 103 | found |
| GWFAM335 | 598 (598) | 0 | 4 / 3 / 2 | 2 (1) | S | 0 / 162 | i |
| GWFAM348 | 707 (707) | 0 | 2 / 6 / 1 | 1 (0) | - | 0 / 200 | i |
| GWFAM398 | 1000 (1000) | 0 | 2 / 2 / 0 | 0 (0) | - | 0 / 0 | i |
| GWFAM401 | 1038 (1000) | 290 | 0 / 4 / 0 | 0 (0) | - | 0 / 172 | iii |
| GWFAM402 | 524 (524) | 134 | 3 / 4 / 1 | 1 (1) | D | 118 / 111 | found |
| GWFAM407 | 1246 (1000) | 279 | 10 / 1 / 2 | 1 (1) | D | 494 / 485 | found |
| GWFAM415 | 164 (164) | 58 | 1 / 3 / 0 | 0 (0) | - | 0 / 53 | ii |
| GWFAM425 | 1082 (1000) | 434 | 1 / 2 / 2 | 2 (2) | DD | 334 / 375 | found |
| GWFAM439 | 482 (482) | 0 | 4 / 5 / 1 | 1 (1) | S | 0 / 0 | i |
| GWFAM440 | 204 (204) | 1 | 3 / 0 / 0 | 0 (0) | - | 0 / 39 | ii |
| GWFAM490 | 916 (916) | 142 | 2 / 3 / 3 | 3 (2) | DD | 88 / 140 | found |

Register rows 1216-1219.
