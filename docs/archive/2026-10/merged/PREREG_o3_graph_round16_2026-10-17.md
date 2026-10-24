# Round 16: integrated v12 (the O2 margin) on every unused family of round 15's pool, and the cost of full families in CPU seconds (2026-10-17)

Pre-registration, written before any world of round 16's families was built: no read of these families has been simulated.
Code: branch `o3-fast` (`bench/o3_graph/`, `bench/o3_unexplained/psvclust.py`), frozen at the commit that adds this file.

**Why.**
- **Round 15 failed 3 of 6 cells (register 1316).** Every failing bar was a few FAM90A errors where the current method makes
  none. Read by read, path O2's fixed margin of 2 edits:
  1. rewarded a copy for read-derived paths its siblings lacked;
  2. took a locus-fallback call on bases its own alignment clipped;
  3. took a call that no rival's paths could answer.
- **Round 12 A119b, the one development run v10 and v11 never passed,** adds two more mechanisms:
  4. a margin of two 1-bp homopolymer indels against 63 DUX4 rivals;
  5. a retained intron that only one copy's annotation absorbs.
- **Four rules in O2's genome step address them.** They were developed on spent runs only (section 16 of
  `docs/O1O2O3_FAST_REALISTIC_2026-10-16.md`).
- **Round 15's cost bar compared wall times at unequal threads:** the current method's minimap2 ran at 8 threads, each
  integrated step at 2. Round 16's cost bar is in CPU seconds, which `fastset` now records for every step (d6a8e6ba).

## The two systems

**Current:** copy_assign's O2 and `--psv-clusters` O3, unchanged since round 1.

**Integrated v12** = v11 + four O2 rules. Each acts in the genome step (`o2.integv`), so O3's input (`o2.integ`) is unchanged:
- **`--veto-abstain`** (1bd6693a, mechanism 1): a path call is abstained when a rival's worst-case genome cost (every clipped
  base paid) beats the call's best case (clips under 30 bases free) by the margin, and no confident winner takes the read.
- **`--locus-pairwise`** (10208f8c, mechanism 2): against each rival, a locus-fallback call pays its edits plus the read bases
  the rival aligns and it clips; it must still win by the margin.
- **`--o2-significance`** (1bd6693a, mechanism 4): each path edit is priced at -ln of its type's rate in the world's uniquely
  placed reads (mismatch, homopolymer indel, other indel). A call must beat every rival by ln(rivals / alpha), alpha = 1e-3,
  as copy_assign's certificate divides alpha by its rivals.
- **`--structure-check`** (10208f8c, mechanisms 3 and 5): a rival whose own locus explains the read better than its paths
  do, by the margin, is compared with the call on the genome; the call is abstained unless it wins there by the margin.

**On all 46 spent development runs (276 worlds), against v11's table** (register 1317):
- 36 -> 40 runs pass. The four that turned are round 15 A119b 502 and 503, round 15 public 503, and round 12 A119b. No world
  that passed fails.
- Correct -1,294 (-0.21%) and errors -139 (-2.1%), of v11's 603,857 and 6,617.
- The six still failing are the O2 proof's cells, on their OUTSIDE twins.
- The thinnest B1 margins left are 3-6 reads.

**Integrated flags (`$V12` below):**
```
--o2-distance gapcomp --member-genome-check --member-attribution --union-exact --o2-genome veto --per-read-detectability
--twins --error-floor --o1-hosts --branch-genome-novel --mosaic-shared --haplotype-check --haplotype-check-genome
--branch-linear --aln-cache --veto-abstain --locus-pairwise --o2-significance --structure-check
--sim-het-draw --heavy-current --step-timeout 1800 --het 0 --jobs 2 --threads 2
```

## Families (held out; no draw)

**Every family of round 15's pool with >= 3 copies that round 15 did not draw** (`pool16.py`, d6a8e6ba):
- **The pool** is round 15's (`r15/families15_all.json`): every unused CAT protein-coding stem with >= 3 genes, built by
  round 12's rule. 105 stems formed a family and 30 have >= 3 copies; round 15 drew 8 of those.
- **Round 16 takes the other 22.** Taking all of them leaves nothing to select.
- **22 families, 90 copies.** The largest are OR4A (7 copies), OR2L and VN1R (6). FAM90A- and DUX4-scale arrays are not in
  this set; they remain tested on spent runs only.
- **File:** `/mnt/linuxdisk/tmp/o3_graph/families16.json` (set `r16`, sha1 3698cae80929).

| family | copies | isoforms | chromosomes |
|---|---|---|---|
| CLEC | 5 | 20 | chr16 |
| CST | 3 | 4 | chr20 |
| GGTLC | 4 | 8 | chr20, chr22 |
| KCNJ | 3 | 4 | chr17 |
| LGALS | 4 | 16 | chr17 |
| MBD3L | 5 | 5 | chr19 |
| MKRN | 3 | 6 | chr7, chr12, chrX |
| MTRNR2L | 4 | 4 | chr3, chr5, chr11, chr17 |
| MYL | 3 | 9 | chr6, chr18 |
| NUTM | 4 | 10 | chr10 |
| OR13C | 3 | 3 | chr9 |
| OR2J | 4 | 5 | chr6 |
| OR2L | 6 | 9 | chr1 |
| OR2M | 5 | 6 | chr1 |
| OR4A | 7 | 7 | chr11 |
| OR4N | 3 | 3 | chr14, chr15 |
| SLC35G | 4 | 4 | chr8, chr17, chr18 |
| SLC6A | 3 | 10 | chr16, chrX |
| SPAG | 3 | 12 | chr8 |
| TAS2R | 4 | 4 | chr12 |
| TRPC | 4 | 4 | chrX, chrY |
| VN1R | 6 | 6 | chr13, chr21, chr22 |

## Worlds and runs

- **Worlds:** six world types (INTACT, CROSS, OUTSIDE, REMOVED, BRANCH, SILENT) under `profiles/A119b_testis_v2.json` and
  `profiles/public_testis_v2.json`.
- **Seeds:** 601, 602 and 603, fixed with no seed selection: 6 cells (seed x profile), 36 worlds.

```
F=/mnt/linuxdisk/tmp/o3_graph/families16.json; O=/mnt/linuxdisk/tmp/o3_graph/fast/heldout16
for S in 601 602 603; do for P in A119b public; do
  python3 bench/o3_graph/fastset.py run $O/${P}_s$S $F r16 $S --local $O/local $V12 --sim-profile bench/o3_graph/profiles/${P}_testis_v2.json
done; done
```

## Bars

**Correctness (as round 15).** Read outcomes are per-read detectability less each world's error floor, scored on the veto
table (`--o2 o2.integv.m2.tsv`). Bars per world:
- **B1:** integrated correct >= current correct.
- **B2:** integrated errors <= current errors.

**A cell PASSES when all six of its worlds pass both bars. Round 16 PASSES only when all six cells pass** (72 bar checks).

**Cost (C1), in CPU seconds.** Measured on the mainstay set (NPIP as one family, TBC1D3, GOLGA8, GOLGA6L, CT47A, XAGE,
SPATA31; seed 306; both profiles), as in round 15. These families are spent, but cost is not tuned on them.
- **Setup:** `fastset.py run` with v12's flags, the existing mainstay local reference, nothing else running:
  ```
  M=/mnt/linuxdisk/tmp/o3_graph/fast/r14dev/mainstay_npip306/local
  for P in A119b public; do
    python3 bench/o3_graph/fastset.py run /mnt/linuxdisk/tmp/o3_graph/fast/r16cost/v12/$P \
      /mnt/linuxdisk/tmp/o3_graph/families_mainstay_npip.json mainstay 306 --local $M $V12 \
      --sim-profile bench/o3_graph/profiles/${P}_testis_v2.json
  done
  ```
- **The measure:** each step's `cpu` in `timing.json` (fastset's `/usr/bin/time` %U + %S).
- **Bar:** v12's integrated CPU seconds (user + system, children included), summed over the six worlds (O1 + O2 +
  O2_genome + O3), at most the current method's summed CPU seconds, in each profile.
- **Reported, no bar:** wall times (at the flags' threads: current 8, each integrated step 2), peak RSS, and the CPU seconds
  of round 16's six cells.

**Rules:**
- Nothing is changed after the first run.
- A bug is recorded as an amendment, with its reason, before any rerun.
- A run that crashes is rerun once, and the crash is recorded.

## What a PASS would and would not claim

**A correctness PASS** would claim that v12 does at least as well as the current method in every world type, on 22
held-out families, under both profiles and three seeds.

**A cost PASS** would claim that, on NPIP-sized families, the integrated O1 -> O2 -> O3 uses no more CPU than the current
method.

**Neither would claim:**
- anything on real reads;
- that near-identical uncatalogued twins are handled (the O2 proof's end-to-end failure, register 1314-1315);
- anything about arrays of FAM90A's or DUX4's size beyond the spent runs.

## Results (2026-10-17, 20:33-21:16; nothing changed after the first run, no crash)

**Verdict: FAIL.** 4 of the 6 cells pass; the bar was all six. **C1 PASSES** in both profiles.

| cell | verdict | failing worlds (integrated vs current correct) | integrated / current CPU |
|---|---|---|---|
| A119b 601 | FAIL | INTACT 4,232 vs 4,238; SILENT 3,224 vs 3,231 | 148 / 240 s |
| public 601 | PASS | | 119 / 211 s |
| A119b 602 | FAIL | SILENT 2,638 vs 2,642 | 168 / 287 s |
| public 602 | PASS | | 132 / 270 s |
| A119b 603 | PASS | | 167 / 273 s |
| public 603 | PASS | | 135 / 245 s |

**B2 holds in all 36 worlds.** The integrated pipeline makes no error where the current method makes none, and far fewer
elsewhere: OUTSIDE 44-182 against 811-1,214, REMOVED 0-3 against 39-132, BRANCH 157-304 against 285-653. Over the 36 worlds it
makes 1,833 errors against 8,796, and 147,583 correct calls against 137,426 (+7.4%).

**Every failing bar is B1, by 4-7 correct reads, in INTACT or SILENT under A119b.** Post-hoc diagnosis (read-only, on the
round's tables):
- **The two methods nearly tie there.** In A119b 601 INTACT, 28 reads are correct under the current method only and 22 under
  the integrated pipeline only.
- **Most of the 28 are reads integrated O2 leaves `unexplained`:** more edits than the fixed divergence gate,
  (delta + floor) x length, allows. The rest are `thin` (margin 1) or locus ties. Round 10's set 10 failure was the same gate.
- **v12's new rules cost 5 of the 13 missing reads (`structure_only`), but not the verdict.** Rebuilt with v11's genome
  step, the cell still fails: INTACT 4,235, SILENT 3,226.
- **The integrated lead is thin in these worlds because the families are small** (3-7 copies). The current method's
  certificate assigns single-PSV reads among few rivals, and short noisy reads.

**C1, the cost on the mainstay set** (nothing else running): v12's integrated stages used 810 CPU seconds (A119b) and 660
(public), against the current method's 1,046 and 949.
- **That passes C1:** 23% and 30% less CPU.
- **Wall time, reported:** 449 s and 379 s against 312 s and 282 s, at the flags' threads (the current method 8, each
  integrated step 2).
- **O2 is two-thirds of integrated CPU:** 533 s and 491 s.
- On round 16's cells the integrated stages are cheaper and faster: 119-168 CPU seconds against 211-287, and 68-100 s wall
  against 118-138.
