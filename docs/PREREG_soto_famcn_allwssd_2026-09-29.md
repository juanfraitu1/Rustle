# Pre-registration: our own WSSD famCN from all local SGDP tracks in Soto's reconciled family recipe

**Written 2026-09-29 (KEY=soto_famcn269) BEFORE any new famCN value was computed or any new arm scored.** Human CHM13
v1.0 only. Follows `PREREG_soto_reconciliation_2026-09-29.md` (EXON × PAIR, register drafts 1156-1160) and its
Independent verification section. Nothing in `src/` or `bench/` is edited; nothing is committed. Scratch:
`/mnt/linuxdisk/tmp/rustle_figures_dev/soto_famcn269/`. New table: `winloci_data/soto_replication/famcn_ours_allwssd.tsv`
(new file; `famcn_ours_all.tsv` is never overwritten). This file binds once Amendment 1 records its sha1.

## 0. The question

The verification's ladder (ALL 491 families; held-out in brackets): sequence only 0.7307 / 345 exact; **our own famCN
(`famcn_ours_all.tsv`, 10 WSSD samples, r = 0.933 with S1C) 0.9197 / 373 (0.9315 / 209)**; S1C famCN 0.9698 / 479
(0.9681 / 263). Soto used SGDP n = 269. **Does computing our famCN from the whole local SGDP panel instead of 10
samples close the 0.9197 → 0.9698 gap, and how much of it?** A second, separate lever is Soto's own per-gene interval
definition (§2), reported as its own row, never folded into the sample-set answer.

## 1. How `famcn_ours_all.tsv` was built (established today from files, before this prereg)

- Code: `bench/soto/famcn_bulk.py` (removed in cleanup wave 2, commit 77df6c4c; text at `77df6c4c^`). Run 2026-08-01
  22:14 with defaults `--samples 10 --tool bigBedToBed`.
- **Samples:** the first 10 of the sorted `*_wssd.bb` names in `soto_wssd/`: LP6005441-DNA_A01, A03, A04, A05, A06,
  A08, A09, A10, A11, A12 (the uppercase `LP…` names sort before `chm13_` / `hg38_`).
- **Intervals:** `sd98_gene_exons.tsv` (21,768 rows, 5,154 genes, CRLF line ends) = **every CAT v4 exon of each of the
  5,154 autosomal SD98 genes, merged per gene** (checked today: identical set to the per-gene merge of `cat_v4.bed`
  blocks; 3,679 of the 21,768 intervals are not fully inside SD98).
- **Per sample:** length-weighted mean of the WSSD window CN (bigBed column 10, float) over all of the gene's
  intervals together; **famCN = median over samples**; `famCN_mad` = median |x − median| over samples.

## 2. What Soto did (their released code, `A_SD98_regions.md` §4 and §2.2; paper l.1041 / l.1047)

- **Samples:** "SGDP (n=269)" (paper). Code: all of `wssd-t2t/t2t_v1.0/sgdp_complete/*bed`, with the comment
  "**we removed LP6005442-DNA_A08.CN.bed because it had outlier copy numbers**". Non-human (t2tdp) tracks are copied
  separately and used only for the great-ape comparison; `/HG`, `/NA`, `/CHM` tracks are excluded.
- **Local panel:** `soto_wssd/` holds 271 bigBeds = exactly the UCSC hub listing fetched today (271): **269 SGDP
  (LP… / SS…)** + `chm13_wssd.bb` + `hg38_wssd.bb`. LP6005442-DNA_A08 is among the 269 (7.8 MB, the light track).
- **Intervals:** `CHM13.combined.v4.genes.SD-98.{protein_coding,unprocessed_pseudogenes}.bed` = the CAT v4 **gene
  body intersected with the merged SD98 regions**, one row per intersected piece (a gene can have several rows),
  genotyped per sample by `genotype_cn_parallel.py` (**not released**; its per-row aggregator is unknown).
- **Use:** notebook cell 6 `get_mad(elements)`: all WSSD rows of the pair's genes → per-row median over samples →
  `stats.median_abs_deviation` over those rows; cell 14 per-gene `median_CN` = median over the gene's rows of the
  per-row medians (presumably S1C "Median famCN").

## 3. Arms (the recipe is fixed: EXON edges `soto_reconcile/exon/edges_exon_v230.all5154.tsv` × PAIR closure,
## node universe = S1C's 2,334 genes; copy numbers given ONLY to the 1,793 S1C genes with a "Median famCN")

| arm | samples | intervals | per-gene value / gate | role |
|---|---|---|---|---|
| S1C | Soto's | Soto's | S1C Median famCN; pair \|Δ\|/2 < 1 | reference (0.9698 / 479) |
| OURS10 | 10 (§1) | all exons merged | median of 10; pair \|Δ\|/2 < 1 | control (0.9197 / 373) |
| I10 | same 10, recomputed by me (pyBigWig) | same | same | instrument check only (must reproduce OURS10) |
| **ALL268 (PRIMARY)** | **the 269 local SGDP minus LP6005442-DNA_A08** (Soto's code removes it) | all exons merged | median of 268 | **the answer: only the sample set changes** |
| ALL269 | all 269 local SGDP | same | median of 269 | sensitivity (A08 kept, = the paper's "n = 269" read literally) |
| SUB-n | random subsets of the 268, n ∈ {10, 25, 50, 100, 200}, 20 draws each (seeds 0-19) | same | median of n | sample-size curve + noise band for n = 10 |
| **SOTOIV-gene** (lever 2) | 268 | gene body ∩ merged autosomal SD98, one row per piece | per-row median over samples → per-gene median over rows; pair \|Δ\|/2 < 1 | interval lever alone |
| SOTOIV-rows (lever 2, literal) | 268 | same | pair kept iff MAD over ALL rows of both genes < 1 (median-MAD = their code; mean-MAD = their prose; both reported); a pair whose genes have no rows is dropped | their `get_mad` as written |

Instrument details (fixed now): per-row / per-gene per-sample value = length-weighted mean of window CN over the
interval set (same as §1; used for SOTOIV too because Soto's aggregator is unreleased — disclosed). Gene body = min
start / max end over the gene's CAT v4 transcripts in `cat_v4.bed` (approximates the GFF3 `gene` feature). SD98 =
`sd98_v1.bed` (merged UCSC v1.0 ≥ 0.98, 97,797,568 bp autosomal). No liftover (tracks and CAT v4 are both v1.0).
A gene/sample with no covered base gives no value for that sample (median over the samples that have one).

## 4. Scoring (reused, not modified)

`soto_reconcile_verify/v_common.py` + `v_nulls.py` (`fast_pair` / `run`, `pair_rule` for the mean-MAD check,
`collapse`, `score`, `cover_exact`); SOTOIV-rows uses my own gate function that differs from `fast_pair` only in the
MAD line. Frozen split `soto_losses/frozen/split.tsv` (49bcbcfe): DEV 225 / HELD-OUT 266 families; scored on the
clean S1C genes of each half. Per arm × {DEV, HELD-OUT, ALL}: **ARI (median-MAD and mean-MAD; identical by
construction for one-value-per-gene arms, verified numerically for ALL268), exact families, pair P / R / F1 (+ false
pairs), bipartite MICRO P / R / F and MACRO P / R, undetected**; on ALL also **cover-aware exact**, **flagship
families** (the 38 of the reconciliation §7, prefix rule of `v_fams.py`; status exact / split / merged / split+merged /
undetected as `recon_lib.family_status`) and **nesting** (share of Soto families with ≥ 2 clean members whose members
all sit in ONE of our clusters; share of our clusters that are unions of whole Soto families; `soto_nest.py` logic).
**CN agreement with S1C** over the 1,793 genes: Pearson, Spearman, median |Δ|, share |Δ| < 1 and < 2, and **gate
agreement** = share of the EXON pairs in the universe with both genes valued whose |Δ| < 2 decision equals S1C's.

## 5. Clauses

- **C1 (sample count is a lever).** ALL268 ALL ARI − OURS10 ALL ARI (0.9197) ≥ +0.010 AND ALL268 − OURS10 > 0 on DEV
  and on HELD-OUT. Pass → "MORE SAMPLES CLOSE x% of the ARI gap (y% of the exact gap)" with x = Δ / (0.9698 − 0.9197),
  y = (exact − 373) / (479 − 373). |Δ| < 0.010 → "SAMPLE COUNT IS NOT THE LEVER". Δ ≤ −0.010 → "MORE SAMPLES MOVE
  AWAY FROM S1C".
- **C1-noise.** C1's pass additionally requires ALL268's ALL ARI to exceed at least 19 of the 20 SUB-10 draws;
  otherwise "NOT DISTINGUISHABLE FROM THE CHOICE OF 10 SAMPLES". Where OURS10 sits among the SUB-10 draws is reported.
- **C2 (interval definition is a lever).** SOTOIV-gene ALL ARI − ALL268 ALL ARI ≥ +0.010 AND > 0 on both halves →
  "INTERVAL DEFINITION IS A LEVER (Soto's gene-body ∩ SD98 pieces)"; |Δ| < 0.010 → "not a lever"; ≤ −0.010 → "moves
  away". SOTOIV-rows is reported beside it, no clause.
- **C3 (reach).** Per arm, whether ALL ARI ≥ 0.9598 (within 0.010 of S1C) — descriptive, no selection among arms.
- No parameter is tuned; nothing is selected on HELD-OUT. The split is reused only so each number has a held-out copy.

## 6. Not blind (disclosed)

OURS10's ALL and HELD-OUT scores are known (0.9197 / 373; 0.9315 / 209) and its r = 0.933 with S1C; the 08-01 memory
reports that 5 vs 271 samples did not change famCN separations on 40 over-merged copies (a different question and
intervals). No famCN from more than 10 samples has been computed on these genes; no SOTOIV value has been computed.

## 7. Predictions (this author, before any new number)

- Pearson(ALL268, S1C) > 0.933: 0.85; ≥ 0.95: 0.50. Pearson(SOTOIV-gene, S1C) ≥ 0.98: 0.40.
- C1: pass 0.35; "not the lever" 0.55; "moves away" 0.10. ALL268 ALL ARI in [0.91, 0.94]: 0.70.
- ALL268 exact ≥ 400: 0.30.
- C2: pass 0.50.
- Any arm reaches ≥ 0.9598 (C3): 0.20.
- ALL269 vs ALL268 |Δ ARI| < 0.005: 0.80.

## 8. Order and machine rules

(1) freeze this file (sha1) → (2) build the per-sample matrices (exon intervals, SOTOIV rows) under `tools/rlock.sh
heavy` (foreground; one pass per sample), I10 reproduction asserted (max |Δ| vs `famcn_ours_all.tsv` ≤ 0.002 on
≥ 99.9% of genes and an identical OURS10 / I10 partition) → (3) Amendment 1 (instrument sha1s, counts only) → (4) one
scoring run under `tools/rlock.sh light` → (5) Outcome + draft register rows (after 1160; not appended). TMPDIR under
`/mnt/linuxdisk`; never `pkill -f`.

## Amendments

### Amendment 1 — the freeze (2026-09-29 13:58, written BEFORE any new arm was scored or compared with S1C)

**This file's sha1 before this amendment:** `169d3c988e43350e58ed854b47c990756a84239d` (frozen 13:48; byte copy
`soto_famcn269/frozen/PREREG_soto_famcn_allwssd_2026-09-29.pre_amendment1.md`). **Instruments** (scratch
`soto_famcn269/`): `build_intervals.py` 453117af, `compute_matrix.py` aa3501bf, `write_table.py` 092a20d0, `run.py`
1248af3a; data `exon_iv.bed` 5dcf61eb, `sotoiv_rows.bed` c6940790, `genebody.bed` 8c4b8695, `all269.npz` 5a58bbee,
`i10.npz` 303c24c4, `sotoiv_rows_famcn.tsv` 8e5ff7ed; **`famcn_ours_allwssd.tsv` 11daa3ce** (new; `famcn_ours_all.tsv`
089aad4b untouched). Scorer imported unmodified from `soto_reconcile_verify/` (`v_common.py`, `v_nulls.py`).

**Built and checked (counts only; no new famCN compared with S1C, no new arm scored):**
- Exon intervals: 21,768 / 5,154 genes / 8.41 Mb (= `sd98_gene_exons.tsv`). SOTOIV rows (gene body ∩ `sd98_v1.bed`,
  autosomal, already merged: 817 regions): 5,217 rows / 5,154 genes / 49.06 Mb; 5,095 genes have 1 row, 57 have 2,
  1 has 3, 1 has 5.
- Matrices: 269 SGDP tracks (LP… / SS…; `chm13_wssd.bb`, `hg38_wssd.bb` excluded), pyBigWig, 4 processes under the
  heavy lock, 133 s; 0 unparsable CN fields; 0 empty gene / row / sample cells (LP6005442-DNA_A08 included, column 86).
- **I10 reproduces `famcn_ours_all.tsv` exactly**: the same 10 samples (LP6005441-DNA_A01, A03-A06, A08-A12), all
  5,154 genes, max |Δ| = 0.000 at the table's 3 decimals, n_samples equal. The first 10 columns of the 269-sample
  matrix are identical to the I10 run. So the ONLY change from OURS10 to ALL268 is the sample set.
- The new table's `famCN` = ALL268 (median over the 268; `famCN_mad` over the same), plus `famCN_269`, `famCN_sotoiv`
  (per-row median over the 268, then median over the gene's rows) and `n_sotoiv_rows`.
- `run.py` asserts, before anything is printed: S1C arm 0.9698 / 479 (ALL) and 0.9681 / 263 (HELD-OUT); OURS10
  0.9197 / 373 (ALL); I10 partition = OURS10 partition; exactly 1,793 genes receive a copy number in every arm.

**Command:** `bash tools/rlock.sh light python3 run.py runs/main` (one run). A crash fixed without changing a rule is
recorded here as Amendment 2.

### Amendment 2 — the run (2026-09-29 14:02; no rule, arm, clause or bar changed)

The first launch crashed before computing anything: `v_nulls.py` reads `sys.argv[1]` as its null-seed count and got
`runs/main`. Fix: `run.py` clears `sys.argv` around the `exec` (sha1 1248af3a → **3410a221**; no other line changed).
One run, 18 s; every built-in assertion passed (S1C 0.9698 / 479 and held-out 0.9681 / 263; OURS10 0.9197 / 373;
I10 partition = OURS10; 1,793 genes valued, = the biotype-eligible set). Output `runs/main/results.json` c95fb5ca,
stdout `runs_main.stdout` 67024f3c. Post hoc diagnostics (after the numbers, not clauses): `posthoc/diag.py` 30f53a30,
`posthoc/diag2.py` 6cb81920.

## Outcome (2026-09-29)

**Verdict. C1 (sample count): "MORE SAMPLES MOVE AWAY FROM S1C" as registered (ALL ARI 0.9197 → 0.8853), but the
whole move is ONE knife-edge family, so read it as "sample count is not the lever". C2 (Soto's interval definition):
PASS. C3: no arm of ours reaches 0.9598.**

| arm (EXON × PAIR; CN on the 1,793 S1C-valued genes only) | ARI ALL (DEV / HO) | ARI ALL without FAM90A | exact ALL (DEV / HO) | pair P / R | MICRO P / R | undet. | cover-aware exact | flagship exact / 38 | nesting: Soto families inside one cluster / our clusters that are unions | Pearson / share \|Δ\|<1 vs S1C |
|---|---|---|---|---|---|---|---|---|---|---|
| S1C famCN (reference) | **0.9698** (.9708 / .9681) | 0.9650 | **479** (216 / 263) | .9995 / .9420 | .9995 / .9797 | 0 | 481 | 36 | 97.5% / 96.2% | — |
| OURS10 = `famcn_ours_all.tsv` (control) | 0.9197 (.9096 / .9315) | 0.9130 | 373 (165 / 209) | .9644 / .8797 | .9281 / .8972 | 37 | 373 | 21 | 85.1% / 81.6% | .932 / .55 |
| **ALL268** (primary: 269 SGDP − A08) | **0.8853** (.9039 / .8607) | **0.9088** | **375** (167 / 208) | .9630 / .8202 | .9228 / .8894 | 37 | 375 | 21 | 85.1% / 81.3% | .932 / .65 |
| ALL269 (A08 kept) | 0.9189 (.9091 / .9298) | 0.9120 | 375 (167 / 208) | .9653 / .8774 | .9304 / .8976 | 37 | 375 | 21 | 85.4% / 81.6% | .932 / .64 |
| **SOTOIV-gene** (gene body ∩ SD98, 268) | **0.9277** (.9227 / .9343) | **0.9251** | **411** (178 / 235) | .9764 / .8842 | .9580 / .9353 | 19 | 413 | 22 | 90.3% / 86.2% | **.977** / .81 |
| SOTOIV-rows (their `get_mad` over rows; mean-MAD) | 0.9280 (.9249 / .9325); mean 0.9259 | 0.9255 | 413 (182 / 233) | .9772 / .8841 | .9579 / .9348 | 20 | 413 | 23 | 90.8% / 86.6% | — |

ARI median-MAD = mean-MAD with identical partitions in every one-value-per-gene arm (checked). Macro P / R, F1s and
false-pair counts per half are in `results.json` / `runs_main.stdout` (e.g. false pairs ALL: 5 / 377 / 366 / 366 / 248 / 239).

**C1 in detail.** Δ ARI ALL268 − OURS10: DEV −0.0057, HELD-OUT −0.0708, ALL −0.0344 (−69% of the 0.9197 → 0.9698 gap);
exact +2 (+1.9% of the 373 → 479 gap). ALL268 beats only 7 of the 20 SUB-10 draws; OURS10 is beaten by only 2 of
them. **Post hoc cause:** FAM90A (ID_356, 56 clean genes, all HELD-OUT, famCN 31-42) holds together as 54 + 2 in
OURS10 / ALL269 but breaks 38 + 16 + 2 in ALL268; exactly **2 of its 1,533 internal pairs** (famCN ≈ 31-36, |Δ| 1.9-2.1)
flip the |Δ| < 2 gate when the single sample LP6005442-DNA_A08 is added or removed. Scored without FAM90A's 56 genes,
ALL268 − OURS10 = −0.0042 (ALL), −0.0057 (DEV), −0.0024 (HO), and ALL269 − OURS10 = −0.0010. The SUB-n curve without
FAM90A: median 0.903 (n = 10) → 0.906 (50) → 0.910 (200), range narrowing 0.892-0.910 → 0.907-0.913; FAM90A's largest
piece is 37-40 in 8-17 of 20 draws at every n from 10 to 200 (a coin flip that more samples do not settle). More
samples DO repair high-CN genes (share |Δ| < 1 at S1C famCN ≥ 50: 0.11 → 0.60; median |Δ| 0.84 → 0.53; the docstring's
"CN > ~100 unresolved at low sample counts" is confirmed) but those genes rarely decide a pair gate: Pearson 0.9325 →
0.9324, gate agreement 0.845 → 0.852. Consistent with 08-01 ("sample count irrelevant", 5 vs 271 samples).

**C2 in detail.** SOTOIV-gene − ALL268: ALL +0.0423, DEV +0.0189, HELD-OUT +0.0736 (the held-out figure includes
FAM90A partly re-joining, 53 + 2 + 1; without FAM90A +0.0163 ALL, +0.0188 DEV, +0.0143 HO); exact +36 (DEV +11,
HO +27); false pairs 366 → 248; undetected 37 → 19; Pearson 0.932 → 0.977; share |Δ| < 1 0.65 → 0.81 (CN 0-5: 0.96,
5-10: 0.89, 10-20: 0.78, 20-50: 0.68, ≥ 50: 0.60); gate agreement 0.852 → 0.885. It closes 36 / 104 = **35% of the
exact-family gap** and (without FAM90A) (0.9251 − 0.9088) / (0.9650 − 0.9088) = **29% of the ARI gap**. The literal
row-level `get_mad` adds 2 exact families (413); 98.9% of genes have one row, so it barely differs.

**What is left (S1C 0.9698 / 479 vs 0.928 / 411-413).** Our CN still disagrees with S1C by ≥ 1 on 19% of the 1,793 genes,
mostly at famCN ≥ 10; S1C-kept pairs we cut 402, S1C-cut pairs we keep 708 (of 9,663 valued pairs). Named, unmeasured
causes: Soto's per-row aggregator (`genotype_cn_parallel.py`, unreleased; we use a length-weighted window mean),
the gene body from the GFF3 `gene` feature (we use the transcript span in `cat_v4.bed`), and whether S1C's median
includes A08 (their code drops it; their paper says n = 269 and the hub holds 269 including it). None was searched.

**Predictions vs outcome.** Pearson(ALL268) > 0.933 (0.85): no (0.9324, unchanged); ≥ 0.95 (0.50): no. Pearson(SOTOIV-gene)
≥ 0.98 (0.40): no (0.977). C1 pass 0.35 / not-the-lever 0.55 / moves-away 0.10: registered "moves away" (one family);
without that family, "not the lever". ALL268 ARI in [0.91, 0.94] (0.70): no (0.8853; 0.9088 without FAM90A). ALL268
exact ≥ 400 (0.30): no (375). C2 (0.50): yes. C3 (0.20): no. |ALL269 − ALL268| < 0.005 (0.80): no (0.0336, the same
FAM90A flip).

**What this changes.** (1) The verification's rung "our own famCN 0.92 (held-out 0.93, 373 exact)" was a 10-sample
draw that happens to keep FAM90A whole (2 of 20 random 10-sample draws score higher). With the whole SGDP panel the
exon-interval famCN gives 0.885 (0.909 without FAM90A) and the same ~375 exact families. (2) With Soto's own interval
(gene body ∩ SD98) and all samples, our famCN reaches **ARI 0.928 (held-out 0.934), 411 of 491 exact (held-out 235 /
266), flagship 22 / 38, 90% of Soto families inside one cluster**; the remaining ~0.04 ARI / ~68 families still need
S1C's exact values. (3) Metric trap: on this recipe ALL-491 ARI swings ±0.035 (held-out ±0.07) on two FAM90A pairs, so
famCN rungs should be quoted with exact families and ARI-without-FAM90A beside ARI. Nothing adopted into `bench/` or
`src/`; nothing committed.

**Draft register rows (NOT appended; next numbers after 1160).**

| # | date | area | claim | verdict |
|---|---|---|---|---|
| 1161 | 2026-09-29 | Soto replication (famCN, sample set) | Computing our WSSD famCN from the whole local SGDP panel (268 = 269 − LP6005442-DNA_A08, Soto's code exclusion) instead of 10 samples moves the reconciled recipe (EXON × PAIR) toward S1C | ⛔ **Not the lever** (prereg `PREREG_soto_famcn_allwssd_2026-09-29.md`, 169d3c98; registered clause "moves away": ALL ARI 0.9197 → 0.8853, held-out 0.9315 → 0.8607). The whole drop is FAM90A (ID_356): 2 of 1,533 pairs at \|Δ\| ≈ 2 flip with one sample; without it −0.004. Exact 373 → 375; Pearson with S1C 0.932 → 0.932; high-CN genes improve (share \|Δ\|<1 at CN ≥ 50: 0.11 → 0.60) but rarely decide a gate. Table `famcn_ours_allwssd.tsv` (new) |
| 1162 | 2026-09-29 | Soto replication (famCN, intervals) | Measuring famCN over Soto's own interval (CAT v4 gene body ∩ merged SD98, one row per piece; `A_SD98_regions.md` §2.2/§4) instead of all merged exons | ✅ **A lever**: ALL ARI 0.8853 → 0.9277 (DEV +0.019, HO +0.074; without FAM90A +0.016), exact 375 → 411 (HO 208 → 235), false pairs 366 → 248, undetected 37 → 19, Pearson 0.932 → 0.977, nesting 85% → 90%; closes 35% of the exact gap to S1C (479). Literal row-level `get_mad`: 0.9280 / 413. Remaining gap (~68 families) = unreleased per-row aggregator / gene-feature bounds / A08 — unmeasured |
| 1163 | 2026-09-29 | Soto replication (metric) | "Our own famCN reaches ARI 0.92 (held-out 0.93)" (verification ladder) is a stable property of the 10-sample table | ⚠ **No — a lucky draw on one family.** ALL-491 ARI on EXON × PAIR swings ±0.035 (held-out ±0.07) on whether FAM90A (56 genes) stays whole; 20 random 10-sample draws span 0.870-0.923 with exact 370-382; only 2 of 20 beat the table. Quote famCN rungs as exact families + ARI-without-FAM90A beside ARI (metric trap). Current rung: 0.928 / 411 exact with Soto's intervals and 268 samples |
