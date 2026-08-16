# O1 — Can full-length read evidence repair E_r? (read-tiling vs. repeat bridges)

**Date:** 2026-08-15/16 · **Objective:** O1 (family definition, edge relation `E_r`)
**Companion to:** `docs/o1_coverage_repair.md` (the coverage-threshold impossibility sweep, wf_0c3783ec-0a8)
**Artifacts:** `/home/juanfra/winloci_scratch/o1_readfix/` — `routeA_edges_ALL.tsv` (731 rows), `ROUTEB_per_edge.tsv` (1,809 rows), `ROUTEB_sweep.tsv` (255 rows), `ROUTE_A_THRESHOLDS.txt` and `ROUTEB_PREDECLARATION.txt` (both written **before** scoring), `attack_circular/`, `attack_expr.py`, `attack_rarefy.py`, `attack_matched.py`.

---

## VERDICT

**It does not work.** Read evidence does not separate real homology from repeat bridges. Both routes are
refuted, in **both species**, on **three independent read arms**, with pre-registered falsifiers firing
and (after the attacks) a **non-circular** positive control passing. The recommendation is
**`reject-all`**: `E_r` **stays a pure sequence relation**. `impossibility_broken: false`.

This is worth knowing precisely because it closes the last obvious route. The previous sweep proved no
*threshold* on `c_long` works. This run proves the failure is **not** a consequence of our truncated
reps: repair the reps, or bypass them entirely with the reads, and **the ordering does not cross**.

The one result that points the right way — the genome-wide **secondary-record indicator S** — is a
**LEAD, not a result**: it breaks P1, fails the ceiling constraint, and inherits every read-dependence
defect below. **Do not promote it.**

---

## 1. ⚠ LIMITS FIRST

Read this section before any number below is quoted anywhere.

### 1.1 Evidence tier (T8) — everything here is PAF/BAM-level re-derivation
No result in this document is an end-to-end catalog run. Confirmation would require
`gw_family_catalog` run end-to-end with a read term inside `E_r`, through `refine_families` +
`distinct_loci >= 2`, on the full catalog, with panel-7/HERC2 inside the ceiling constraint.
**I do not recommend spending it.** The pair-level arm under-counts family splits ~2× and is blind to
whole-family deletion, so every "damage" figure here is a **lower bound** on the damage.

### 1.2 Arm sizes are small
GGO FP n=14 · GGO TP_LOST n=21 · HSA FP n=19 · HSA TP_LOST n=10 · GGO scored-core TP n=1,795 ·
GGO families n=375 · GGO catalog copies n=1,415 in 494 multi-copy families.
Every **FP-rejection rate** here carries a CI ~0.25 wide and **must not be quoted as a point estimate**.
The AUCs at p ≤ 0.001 survive the small n (permutation, wrong direction, replicated across two species
and three read arms) — the failure is an **ordering inversion, not a power problem**.

### 1.3 Species are never pooled; the honest arm is fibroblast
- **GGO testis (OR6737) is the arm that BUILT the catalog** (predates the fibroblast BAM by 4 weeks;
  ρ 0.829 vs 0.612). Its **0/115 unservable rate is circular by construction** and must never be
  quoted as the cost of read-dependence.
- **GGO fibroblast (KB3781)** is the honest arm — and it is the **donor-matched** arm, since
  mGorGor1 = KB3781.
- **HSA has exactly one RNA arm (A119b).** There is no second human library, and there is no gorilla
  read for a CHM13 locus. **HERC2's library-dependence is structurally unmeasurable.**

### 1.4 Structural coverage of the frozen scored core
Repair is gorilla IsoSeq, so **only 1,795/9,424 = 19.0%** of the scored-core TP arm and
**14/74 = 18.9%** of the FP arm are repairable at all. **81% of the frozen scored core is human and is
structurally unreachable by Route B.** Any Route-B statement about HERC2 would be **unfalsifiable, not
negative** — logged as a TRAP.

### 1.5 What read-dependence costs (the definitional price)
Honest arm (fibroblast), **COPY unit**, denominator = 1,415 GGO catalog copies fixed before scoring:

| stratum | count | rate [95% CI] |
|---|---|---|
| zero eligible reads | 205/1,415 | **0.1449 [0.1275, 0.1642]** |
| <3 reads | 259/1,415 | 0.1830 [0.1638, 0.2040] |
| no ≥3-read exon block | 263/1,415 | 0.1859 [0.1665, 0.2070] |

**FAMILY unit:** multi-copy families with ≥1 unservable member **114/494 = 0.2308 [0.1958, 0.2699]**.
**PAIR unit:** scored pairs with ≥1 unservable member **281/1,809 = 0.1553 [0.1394, 0.1728]**;
unservable edges in Route-A primary scoring **14/115 = 0.1217 [0.0739, 0.1940]** (fib), 0/115 (tes,
circular), 0/501 (HSA).

⚠ **Tissue and individual are perfectly confounded between the two GGO arms.** The 12.2% is a
*discordance/undefined* rate. **It must never be attributed to tissue.**

### 1.6 What the attacks broke (details in §6)
Five of this run's own claims do not survive:

1. **RETRACTED — "C4 positive control PASSES".** C4 is **circular by construction**: reads are
   selected *by* their own primary alignment inside A±5 kb, so self-`tile_frac` = 1.0 is inevitable
   (only 1.77% fib / 2.65% tes fail to align back at all; under T20 it excluded nothing).
   **Substitute the divergence-stratified cross-rep control** (§6.1) — which passes.
2. **RETRACTED — "the deficit shrinks 4× (6.14× → 1.53×)".** `c_long` and `tile_frac` are different
   statistics on different scales; the ratio is **not scale-invariant**. The order-**preserving** map
   `x → x^0.25` reproduces the entire "gain" (6.14× → 1.57×) while changing nothing about which edge
   dies first. The scale-free currency is **rank**: HERC2's quantile inside the 19-FP distribution is
   **0.0000 → 0.0000**. Zero movement.
3. **CORRECTED — MRPS17's mechanism.** `n_distinct = 1` is **not** one shared interval between the
   *pair*; it is **one transposable element** shared with **≥56/190 panel reps and ~70% of
   read-bearing loci**. MRPS17 is a **TE magnet**, not a two-locus repeat bridge.
4. **CORRECTED — Route B's loss mechanism.** All 1,415 repaired reps exist in all four arms; the
   60–69% TP loss is **over-extension**, not missing data, with a clean dose-response (§4.2).
5. **REFUTED (my own headline) — "unexpressed copies become unverifiable".** **0/1,415 copies and
   0/494 families are silent in BOTH GGO arms** at the 0-read threshold. The real cost is
   **library-conditionality**, not unverifiability (§6.2).

And one **decisive new defect the sweep never tested**: **no negative control was ever run.** There is
a **manufactured-tiling floor** (§6.1) and a **22.61% cross-library arm-dependence** at 4× the
sampling null (§6.2).

---

## 2. The idea, and whether it escapes the impossibility argument

### 2.1 The proposal
> *"These are IsoSeq FULL-LENGTH reads, and we have the genome. If a long region is covered by the same
> reads it might be real, not just a domain/repeat."*

The previous sweep proved no coverage **threshold** can work: HERC2, a real family, splits at
`c_long ≥ 0.034` while the first FP only dies at 0.05 — TP loss begins **before** FP rejection. The
stated cause was structural confounding: TP pairs lost have median `len_longer/len_shorter` **6.984**
vs **1.282** retained — which *is* the repeat-bridge signature.

The insight: **that asymmetry is an artefact of OUR rep construction, not of biology.** The reps are
truncated; the FLNC reads are not. Consulting reads should therefore split into two signatures:

- **Real family, truncated rep:** full-length reads from copy A **tile** copy B over its whole extent.
- **Repeat bridge:** reads from MRPS17 hit MDM2 only over the **same ~970 bp sub-interval**, however
  long the read is — all reads **pile**, they never tile.

### 2.2 Why it *would* escape R0–R5 (the argument is structurally sound)
- It is evidence from a **third source** (the reads), not another function of the two reference
  sequences — which is exactly why every rule in R0–R5 was unable to do it.
- It **stays local to the pair**: a function of (A, B, reads at A, reads at B). Unlike the promiscuity
  gate R5 it does **not** make `E_r` catalog-dependent, so **P1 seed-invariance is not threatened**.

Both halves of that argument are **correct in principle** and were verified here (§7.2, P1 audit).
The route does not fail on its logic. It fails on its **measurement**.

### 2.3 Two premises of the brief that are factually wrong on this substrate
- **The 6.984-vs-1.282 confound is a HUMAN-CURATED-TRANSCRIPT phenomenon.** MDM2's *gorilla node* is
  **3,232 bp**, not 12 kb — the 12,238 bp is `ref_MDM2`, a human curated transcript. Curated stratum:
  TP_LOST 7.131 / FP_LOST 6.055. **GGO stratum: 3.474 / 2.755.** Gorilla FP pairs have median rep
  ratio **1.725** — i.e. on gorilla **the FPs were never the asymmetric pairs**. TRAP.
- **"The reads are a third source" is false as stated for the GGO scored pairs.** Median
  `tile_frac_long − c_long_rep` = **+0.0003** (tes TP_LOST), **+0.0016** (tes FP);
  **0/14 tes FP and 0/13 fib FP gain >0.10**. A read aligned to a rep re-derives that rep's alignment.
  ⚠ **Caveat that must be stated, not smoothed:** F3 does **not** hold uniformly — fib TP_LOST median
  **+0.0319** with 5/18 gaining >0.10, and on **HSA the FP arm gains >0.10 in 6/19 = 0.3158**, median
  **+0.0407**. And on a fresh all-vs-all panel the reads reach *further* than rep-vs-rep by
  **0.07–0.10** at identity 0.60–0.90. **The reads DO add coverage — they just add it on TE magnets as
  readily as on real homology.** "Not a third source" is a **GGO-scored-pairs statement only**.

---

## 3. Is the confound ours? — the truncation measurement

**Yes, reps are truncated — substantially. And it does not matter.**

`rep_len / read-supported extent`, GGO catalog nodes (n=1,415):

| denominator | fibroblast median | testis median |
|---|---|---|
| union of read exon blocks, ±5 kb | **0.264** | **0.478** |
| union, copy's own bounds | 0.333 | 0.621 |
| **longest single FLNC read** (path-comparable) | 0.634 | 0.776 |
| longest single read, own bounds | **0.727** | **0.920** |

⚠ **TRAP — the isoform-union denominator inflates apparent truncation 0.920 → 0.478 on identical
data.** A rep is **one linear path**; `rustlib.single_path_ceiling`'s own docstring calls the union
denominator unfair. Against the fair path-comparable denominator the median rep already carries
**0.727–0.920** of the longest single read.

**Truncation is not concentrated in the lost pairs.** rep/longest-read median (testis):
TP_LOST **0.715** · TP_RETAINED **0.775** · FP **0.880**. The pairs the rule loses are *not* the
pairs whose reps are worst.

Node substrate confirmed spliced, not genomic: median `rep_len/genomic_span` = **0.1292** (n=1,415).

**The decisive split fails.** Repair *does* collapse the length ratio, exactly as predicted —
testis, best-read repair: TP_LOST **3.474 → 1.872**, TP_RETAINED 1.305 → 1.305 unmoved;
AUC(TP_LOST vs TP_RETAINED) **0.9972 → 0.7905**. **The brief's mechanism is partly right.**
But it yields no discriminator:

| substrate | ratios | AUC(TP_LOST > FP) |
|---|---|---|
| original reps | TP_LOST 3.474 / FP 1.725 | **0.9898**, p=0.0005 |
| after repair | TP_LOST 1.872 / FP 1.518 | **0.6190**, p=0.2509 |

The rule needs AUC ≪ 0.5. It is **0.47–0.66 across all four repair denominators and both arms, never
significant.** The original asymmetry runs the **wrong way** (the lost TPs are *more* asymmetric than
the FPs).

**Why:** the lost TPs are textbook paralogue pairs — GSTM3×GSTM1, TPM3×TPM2, PPP2CA×PPP2CB,
PPP1CA×PPP1CB, KIF5A×KIF5C, KMT2A×KMT2B, NAA11×NAA10, ING4×ING5, SERPINB9×SERPINB6 (19 distinct
families / 21 pairs — independent). Their reps *are* truncated and repair *does* collapse the ratio
(KMT2A×KMT2B 4.14→1.19, GSTM3×GSTM1 4.21→2.13) — but read breadth stays **0.125** and **0.149**.
**The homology genuinely stops.** These pairs share only the CDS, and the full-length reads
**confirm** that rather than overturn it. **Ratio collapse and partner coverage are decoupled.**

---

## 4. Both routes, scored on the frozen arms

### 4.1 ROUTE A — read tiling as a kill-only term on R0 edges

`ROUTE_A_PRIMARY` = `tile_frac_long ≥ 0.50 OR (n_distinct_long ≥ 2 AND pile_conc_long < 0.80)`,
pre-registered in `ROUTE_A_THRESHOLDS.txt` before any score.

**Discrimination — the statistic that would have to work.** UNIT = pair, species never pooled:

| statistic | GGO tes AUC(TP_LOST>FP) | GGO fib | HSA node | needs |
|---|---|---|---|---|
| `tile_frac` | **0.1259** p=0.0005 | 0.2308 p=0.0125 | **0.0947** p=0.0010 | ≫0.5 |
| `pile_conc` | **0.7109** p=0.0260 | 0.5662 p=0.5432 | **0.7816** p=0.0135 | ≪0.5 |
| `n_distinct` | 0.5340 p=0.5722 | 0.5491 p=0.6277 | 0.3947 p=0.2864 | ≫0.5 |

**Pre-registered falsifier F1 fires for `tile_frac` AND `pile_conc`, in BOTH species.**
P-A1/P-A2 are **inverted**: the pairs the coverage rule loses tile the partner **less**
(GGO tes MED 0.188 vs FP 0.358; HSA 0.191 vs 0.409) **and pile harder**
(0.878 vs 0.744; 0.984 vs 0.583) than the false positives do.

`n_distinct` carries no signal at all — **median 1 in EVERY group including TP_RETAINED**; real
families do not give multiple disjoint intervals. It is a **length statistic in disguise**:
Pearson **r = 0.9246** with `len(longer rep)` (n=717).

**FP rejection is bought, not earned.** UNIT = pair, kill-only on R0 edges:

| arm | FP rejected | TP_LOST retained | TP_RETAINED retained |
|---|---|---|---|
| GGO tes | 9/14 = 0.6429 [0.3876, 0.8366] | **2/21 = 0.0952 [0.0265, 0.2891]** | 60/80 = 0.7500 — **deletes 25%** |
| GGO fib | 7/14 = 0.5000 [0.2680, 0.7320] (1 FP unservable) | **5/21 = 0.2381 [0.1063, 0.4509]** | 56/80 = 0.7000 — **deletes 30%** |

**The declared 240-cell sweep (t × q × d, emitted in full, no cherry-picking):** 189/240 cells reject
≥1 FP. **Only 5 cells reject ≥1 FP while retaining ALL 21 TP_LOST — and in EVERY one of those 5 the
single FP rejected is BLMH × SLC38A6, which the shipped R2@0.20 already rejects.**

> ⭐ **ROUTE A REJECTS 0 NOVEL FPs AT ANY TP-SAFE CELL.** In the testis arm it does not even reject
> that one (BLMH × SLC38A6 `tile_frac` 0.1905 → kept), i.e. **a net regression on the only gorilla FP
> rejection the shipped rule has.**

**HSA-node FP floor:** **0/19 = 0.0000 [0.0000, 0.1682]** HSA FPs sit below HERC2's `tile_frac`
bottleneck — 0 rejectable while HERC2 survives.

**FAMILY unit vs. constraint C1** (≥95% of R0 in *both* species: GGO ≥352/375, HSA ≥257/274). Only the
HSA panel-7 arm has families computable here:

| t | families intact | FP rejected |
|---|---|---|
| 0.00 / 0.05 / 0.10 | 6/7 | **0** |
| 0.20 | **5/7 = 0.833 — below C1, and the family that drops IS HERC2** | — |
| ≥0.60 | 3/7 | — |
| 0.80 | 1/7 | — |

**No Route-A operating point with any FP rejection satisfies C1.** Deleting 25–30% of retained TP
edges at the pair level cannot leave the GGO family ceiling intact.

**Scale-free damage accounting** (attack 3's correction — TP edges that die *before* the first FP dies,
kill-only, unit = pair, never pooled):

| arm | `c_long` | `tile_frac` (Route A) | `R6` (§6.3) |
|---|---|---|---|
| HSA node (TP 90 / FP 19) | 10/90 = 0.111 | **9/90 = 0.100** | 7/90 = 0.078 |
| GGO tes (TP 101 / FP 14) | 17/101 = 0.168 | **11/101 = 0.109** | 10/101 = 0.099 |

**The honest movement is 10→9 and 17→11. Never to 0.**

⚠ Named counterexample: **VIM × PRPH is TP_RETAINED with `tile_frac_long = 0.0000`** (fib) — a real
edge that any threshold > 0 deletes.

### 4.2 ROUTE B — repair the reps, then apply the shipped rule

Two repairs, both with the minus-strand guard: **B1** = union of read exon blocks (the brief's
prescription), **B2** = single longest FLNC read (path-comparable, the route's strongest form).

| substrate | TP kept @t=0 (n=1,795) | FAM intact @t=0 (n=375) | first TP loss | first FP death | ordering |
|---|---|---|---|---|---|
| ORIGINAL | 1795 = 1.0000 | **370 = 0.9867 [0.9692,0.9943] — PASS** | 0.13 | 0.20 | FAILS |
| B1_union_fib | 563 = 0.3136 [0.2926,0.3355] | **70 = 0.1867 — FAIL** | 0.03 | 0.05 | FAILS |
| B1_union_tes | 717 = 0.3994 [0.3770,0.4223] | **71 = 0.1893 — FAIL** | 0.09 | 0.09 | FAILS |
| B2_path_fib | 653 = 0.3638 [0.3418,0.3863] | **75 = 0.2000 — FAIL** | 0.09 | 0.20 | FAILS |
| B2_path_tes | 707 = 0.3939 [0.3715,0.4167] | **84 = 0.2240 — FAIL** | 0.11 | 0.13 | FAILS |

**Repair deletes 60–69% of true GGO edges and 78–81% of GGO families BEFORE `c_long` is consulted.**
C1 needs ≥352/375; repair gives 70–84/375 — **FAIL by ~5×, at t=0, before any threshold.**

⚠ **TRAP — apparent FP rejection at t=0 is the SUBSTRATE, not the rule:** B1_fib 12/14 = 0.8571,
B1_tes 8/14 = 0.5714, B2_fib 12/14 = 0.8571, B2_tes 10/14 = 0.7143 — with **no `c_long` term
whatsoever**. **FP rejection is only interpretable CONDITIONAL on the edge surviving the substrate.**
Conditional on surviving repair the FP arm is **n = 2, 6, 2, 4** — below the pre-declared floor of 10.
**The ordering question is not scoreable on repaired reps.** Stopped there, as pre-declared; no
threshold quoted.

**The mechanism runs BACKWARDS.** Paired, same pairs, median Δ`c_long`: **−0.2376** (B1_fib),
**−0.2306** (B1_tes), **−0.0407** (B2_fib), **−0.0648** (B2_tes); moved **up** in only 6.9%–25.0% of
TP pairs. Lengthening a rep adds sequence the partner does not cover, so coverage-**of-the-longer**
falls. **Truncation was suppressing the DENOMINATOR.**

**The loss is OVER-EXTENSION, not missing data** (attack 3's correction — `rb_final.py:88-90` writes
`NA` when no edge survives, *not* when a rep is missing; all 1,415 repaired reps exist in all four
arms, 1070/1070 TP nodes present). Dose-response by quintile of `max(len_repaired/len_shipped)`:

| substrate | loss rate by inflation quintile | quintile median inflation | r(inflation, LOST) | r(inflation, Δc_long) among survivors |
|---|---|---|---|---|
| B1_union_fib | 0.253 0.688 0.777 0.799 **0.914** | 1.15 → 7.73 | +0.3300 | **−0.5089** |
| B1_union_tes | 0.242 0.613 0.694 0.666 **0.788** | 1.27 → 4.41 | +0.2697 | **−0.5597** |
| B2_path_fib | 0.189 0.699 0.602 0.827 **0.864** | 1.03 → 2.85 | +0.2983 | **−0.5442** |
| B2_path_tes | 0.228 0.643 0.641 0.747 **0.772** | 1.09 → 2.23 | +0.2287 | **−0.6250** |

**Route B replaced truncation with over-extension and destroyed edges in proportion to how much it
lengthened the rep** — exactly the boundary-work failure mode (recall 0.698→0.788, precision
0.709→0.269).

**Two further disqualifiers:**
- **B1 rebuilds an object the frozen SPEC already rules out.** Union reps are median **3.795×** (fib)
  / **2.094×** (tes) longer than shipped exon-sum reps — the `q915_exon.fa` class explicitly marked
  *"NOT A VALID E_r SUBSTRATE"*.
- **Repair raises BLMH × SLC38A6 from 2.75 → 1.64**, destroying the only gorilla FP rejection the
  shipped rule has.

**On the ORIGINAL substrate the ordering is measurable and fails:** AUC(first-lost TP decile > FP) =
**0.2709, perm p = 0.0015**, where the corrected rule inequality needs ≫0.5. Same direction as Route
A's 0.1259. *(Disclosure: pre-declaration §5 P2 carries a sign typo, corrected in §9 **before**
interpretation.)*

**Controls that pass.**
- **C1 instrument reproduces the frozen sweep exactly**: R0 FAM 370/375; R2@0.20 TP 1774/1795,
  FP_rej 1/14, FAM 355/375. **The movement is the repair, not the code.**
- **C2 minus-strand guard**, correct vs. deliberate per-record-revcomp bug, best **colinear** span:
  B2 tes 0.9876 vs 0.5081, B2 fib 0.7966 vs 0.4444, **40/40 wins, 0 losses**, identity 1.0000.
  ⚠ The first version of this self-test used union-over-all-records coverage and was **blind** to the
  bug (1.00 vs 1.00) — **any reuse must use the colinear form.**

---

## 5. HERC2 and MRPS17 by name — has the ordering crossed?

### 5.1 HERC2 — **NO. The threshold MOVES; the ORDERING DOES NOT.**
Family unit; connectivity bottleneck = largest threshold keeping the family one component. Verified
from the per-edge table, not asserted.

| statistic | HERC2 bottleneck | HSA FP floor (min over 19) | ordering |
|---|---|---|---|
| shipped `c_long_rep` | **0.0325** (reproduces the stored 0.034) | 0.1998 (PTPN22 × LRR1) | **FAILS** |
| Route A `tile_frac_long` | **0.1593** | 0.2437 (RPLP1 × MAP7) | **FAILS** |
| `R6` length-normalised (§6.3) | **0.5838** | 0.5079 (PTPN22 × LRR1) | **FAILS** (7/19 FPs still below) |

**Every threshold that keeps HERC2 intact still rejects 0/19 = 0.0000 [0.0000, 0.1682] of the HSA
false positives.** HERC2's bottleneck edge is `L~chr15_25853560_26064943 × L~chr15_28620350_28638850`
(114,109 bp rep, 120 reads, `tile_frac` 0.1593, `pile_conc` 0.0723, `n_distinct` 44).

⚠ **RETRACTED: "the split threshold moves 4.90× and the deficit shrinks 4× (6.14× → 1.53×)".** The
ratio is not scale-invariant (§1.6 item 2). **Rank quantile inside the FP distribution: 0.0000 →
0.0000.** Zero movement. `gap/IQR(FP)` is no defence either (only affine-invariant): +1.319 → +0.304
as reported, but the same order-preserving relabelling sends it to +3.903.

⚠ **The comparison is n=1 vs n=1.** HERC2's bottleneck node `L~chr15_28620350_28638850` has
**degree 1** — a pendant leaf; its single edge is a bridge. Leave-one-edge-out over the 15 HERC2
edges: removing it sends the bottleneck to **0.0000** (2/15 removals change it at all). The HSA FP
floor is likewise a minimum over n=19 (`c_long` 2nd-lowest 0.2119, bootstrap [0.1998, 0.2990];
`tile` 2nd-lowest 0.2704, bootstrap [0.2437, 0.3010]). **The magnitude should never have been quoted.**

⚠ **NEW TRAP — HERC2 and the FP floor are not on the same substrate, and HERC2 is at its arithmetic
ceiling.** `rep_len/genomic_span`: panel-7 nodes MED **0.6148** (HERC2's ten nodes 0.3083–0.6718, one
rep is **114,109 bp**) vs the GGO catalog's **0.1292** and the HSA-node FP arm's reps of 511–4,488 bp.
Since `c_long ≤ len_short/len_long`: HERC2's bottleneck edge (5,704 vs 114,109 bp) **cannot exceed
0.0500** — measured 0.0325 = **65.0% of the attainable maximum**. The FP floor pair PTPN22 × LRR1
(1,335 vs 3,394 bp) has ceiling **0.3933**, i.e. **7.9× more headroom**, and uses only 50.8% of it.
**"HERC2 splits at 0.034 while the first FP dies at 0.20" is substantially a statement about a 20× rep
-length ratio in a non-shipped rep class.** Re-state the impossibility argument's flagship on a
length-comparable pair, or the advisor's first question kills it.

⚠ **TRAP — Route B cannot test HERC2 at all.** Its 15 scored-core rows are `T4_PANEL7`, species HSA,
CHM13 loci. There is no gorilla read for a human locus, so its 0.034 **cannot move in either
direction** under a gorilla-read repair. Declared before scoring. **"The threshold didn't move" under
Route B would be unfalsifiable, not negative.**

⚠ **"HERC2 remains the FIRST family to die" is statistic-specific, not a property of the catalog** —
under R6 the first family to die is **NPIP (0.5608)**, not HERC2 (0.5838).

### 5.2 MRPS17 — **the flagship FP IS rejected, but only by a rule that also kills 25–30% of retained true edges and drops HERC2.**

All four GGO MRPS17 hub edges (GWFAM210) have `tile_frac_long < 0.50` **and** `n_distinct_long = 1`,
so `ROUTE_A_PRIMARY`'s OR-clause fails and **all four are killed**:

| arm | pair | `c_long_rep` | `tile_frac_long` | `pile_conc` | `n_distinct` |
|---|---|---|---|---|---|
| tes (11 reads) | ×MDM2 | 0.2970 | 0.2970 | 0.727 | 1 |
| tes | ×PIGX | 0.3268 | 0.3267 | 0.778 | 1 |
| tes | ×EIF3J | 0.2068 | 0.2089 | 0.429 | 1 |
| tes | ×TRAPPC2 | 0.3621 | 0.4382 | 0.389 | 1 |
| fib (100 reads) | ×MDM2 | 0.2970 | 0.3013 | 0.9295 | 1 |
| fib | ×PIGX | 0.3268 | 0.3306 | 0.9307 | 1 |
| fib | ×EIF3J | 0.2068 | 0.2133 | 0.4924 | 1 |
| fib | ×TRAPPC2 | 0.3621 | 0.4423 | 0.4761 | 1 |
| HSA node (62 reads) | IVD×MRPS17 | 0.2119 | 0.2839 | 0.340 | 1 |
| HSA node | KCTD2×MRPS17 | 0.2990 | 0.3010 | 0.432 | 1 |
| HSA node | MRPS17×TRAPPC2 | 0.3621 | 0.4456 | 0.525 | 1 |

**THE COST:** that same rule retains **2/21 = 0.0952** of TP_LOST and 60/80 = 0.7500 of TP_RETAINED
(tes), 5/21 and 56/80 (fib), and takes the HSA family ceiling to **5/7 = 0.833 by killing HERC2**.
**At every TP-safe sweep cell MRPS17 SURVIVES.**

⚠ **Statistically, the "4/4" is n=1.** Length-matched TP control (TP edges with `len_long` within
0.8–1.25× of each MRPS17 row) is killed at **tes 32/87 = 0.368 · fib 21/79 = 0.266 · hsanode
22/81 = 0.272** (all-TP rates 0.386/0.307/0.278). Treating the 4 hub edges as independent gives
one-sided Fisher p = 0.0220/0.0069/0.0241 — **but they are one hub = one of the frozen FP set's 24
independent mechanism components**, so n=1 and p collapses to the baseline **0.27–0.37, n.s.**
**The flagship rejection is one draw from a rule that kills a quarter to a third of length-matched
true edges.**

⚠ **CORRECTED MECHANISM — MRPS17 is a TE magnet, not a two-locus repeat bridge.** Its rep is 1,833 bp.
Reads from **116/164 (fib)** and **89/190 (tes)** *other-family* nodes tile it ≥0.20. The four shipped
GWFAM210 partners score MDM2 0.5406 / EIF3J 0.5445 / TRAPPC2 0.5412 / PIGX 0.5412 — **and unrelated
loci score the same**: ZKSCAN3 0.5406, CCNC 0.5423, FYN 0.5423, RPL39 0.5406, KIF5A 0.5406, ING4
0.5423, VIM 0.5308, GSTM1 0.5325, SLC38A6 0.5336. They hit the **identical interval**: MDM2
(756,1747), TRAPPC2 (750,1742), FYN (750,1744), RPL39 (751,1742), ZKSCAN3 (752,1743), VIM (767,1740),
GSTM1 (762,1738), SLC38A6 (757,1739). **One ~990 bp element = 0.54 of the rep; excised and realigned
it hits 56/190 panel reps.** This *is* the brief's "~970 bp sub-interval" — but it is **not a bridge
between MRPS17 and MDM2**, it is a **locus-independent element shared with ~70% of read-bearing
loci**. Magnet census (foreign-family nodes tiling ≥0.20, fib): TRAPPC2 129/164 = 0.787, MRPS17 0.707,
MDM2 0.646, PIGX 0.348, EIF3J 0.201 — **five of the top five magnets are GWFAM210**.
**`n_distinct = 1` is measuring one transposable element, not one relationship.**
⭐ **Honest counterweight:** the #2 magnet in both arms, **GWFAM449_1, is a TP_LOST node** (0.726 fib /
0.505 tes). **Magnetism does not separate FP from TP either.**

**Under Route B** the hub loses its edges to the **substrate**, not to `c_long` — an uninterpretable
"rejection": ×TRAPPC2 0.3621 → no edge (B1_fib 0.0877); ×MDM2 0.2970 → no edge (B1_fib 0.1006);
×PIGX 0.3268 → B2_tes 0.3034 (**survives**); ×EIF3J 0.2068 → B2_tes 0.2095 (**survives**). All shipped
values sit **above** the first TP loss at 0.13. Highest-`c_long` GGO FP overall remains
**DNAJA1 × LAGE3 at 0.7398** — whose *repaired* ratio is 8.19 with read breadth **0.755, the highest
breadth in the entire panel, and it is a false positive.**

---

## 6. The attacks and how each fared

### 6.1 Attack — circularity of the read placement (**PARTIALLY lands; strengthens `reject-all`**)
- **Limb 1 — "reads were placed by the same aligner": CONFIRMED as a mechanism, but SMALL.** Over the
  190 Route-A node windows, best-scoring rep is **not** the read's own node in 0.0517 (fib) / 0.0351
  (tes); is in a **different family** in 0.0306 / 0.0233; truncation-robust criterion **0.0116 /
  0.0107**. The naive form is **refuted** — MAPQ==0 is 12/199,463 = 0.0001 (fib) and 108/28,284 =
  0.0038 (tes); placement is not a coin flip. ⭐ But **96.4% (fib) / 96.1% (tes)** of reads whose best
  rep is in another family carry **MAPQ ≥ 60** — an **independent replication on a fresh substrate of
  the excision finding that minimap2 is confidently wrong and MAPQ is uninformative.**
- **Limb 2 — "circularity produced the inverted AUC": THE ATTACK FAILS.** Re-derived from scratch
  (tes 0.1293 vs 0.1259; fib 0.2393 vs 0.2308). Deleting every read whose best rep is not its own
  node: tes 0.1259 → **0.1224** (p=0.0005); fib 0.2308 → **0.2094** (p=0.0050). **Removing the
  circularity makes the inversion slightly worse. The negative result is not an artefact.**
- **Limb 3 — ⚠⚠ NO NEGATIVE CONTROL WAS EVER RUN.** `grep` for null/random/negative control in both
  pre-declarations returns nothing. The manufactured-tiling floor (ordered pair, different family,
  **zero rep-vs-rep alignment record** under shipped `E_r` params):

  | arm | n pairs | tile≥0.10 | tile≥0.20 | tile≥0.50 | max |
  |---|---|---|---|---|---|
  | fib | 26,911 | 0.0800 | **0.0219** | 0.0017 (47) | **0.7026** |
  | tes | 31,646 | 0.0482 | **0.0112** | 0.0006 (20) | 0.6623 |

  Cross-family pairs that *do* have a rep alignment but fail R0: tile≥0.20 in **0.2583** (fib) /
  **0.2281** (tes). Home-filtering the reads cuts the null ≥0.20 rate to 0.0144/0.0069 — **foreign
  reads account for ~35–55% of the top manufactured signal; the rest is the read carrying repeat
  sequence its own rep does not have.**
- **Limb 5 — C4 is circular; the control that isn't PASSES.** Sensitivity at the operating divergence
  (R0 rep pairs on the panel, restricted to `c_long ≥ 0.85` so full tiling is expected):

  | identity | fib median `tile_frac` | frac <0.50 | tes median | frac <0.50 |
  |---|---|---|---|---|
  | 0.60–0.80 (n=7) | 0.9494 | **0/7** | 0.9410 | **0/7** |
  | 0.80–0.90 (n=60/62) | 0.9630 | **0/60** | 0.9617 | **0/62** |

  **The instrument is sensitive in the regime. The sweep's conclusion is right; its stated reason was
  not.**
- ⚠⚠ **HARD PROHIBITION:** the insight naturally invites promoting tiling from kill-only to an **ADD**
  rule. At tile≥0.50 that adds 47 (null) + 137 (sub-R0) = **184 cross-family edges in fib** and
  20 + 97 = **117 in tes**, against 1,021/1,086 same-family R0 pairs on a 190-node panel — an
  **11–18% cross-family inflation** before scaling to 1,415 copies. **Never as an add rule.**

### 6.2 Attack — expression / library dependence (**PARTIALLY lands; kills the surviving LEAD**)
- **CLAIM "unexpressed copies become unverifiable": REFUTED (the attacker's own headline dies).**
  COPY unit, both arms: **silent in BOTH arms 0/1,415 = 0.0000 [0.0000, 0.0027]** at 0 reads;
  0/1,415 at <3 reads; 102/1,415 = 0.0721 at <10 reads. **FAMILY unit, ≥1 both-silent member:
  0/494 = 0.0000 [0.0000, 0.0077]**. Spearman(fib depth, tes depth) = 0.6672. **Two libraries jointly
  serve the entire GGO catalog.** The 14.49% is a **library-conditionality** rate, not an
  unverifiability rate.
- ⭐ **CLAIM "the same pair gains/loses an edge between substrates": CONFIRMED — and it is the real
  defect.** PAIR unit, 115 GGO pairs, `ROUTE_A_PRIMARY`, fib vs tes:
  **ARM-DEPENDENT 26/115 = 0.2261 [0.1592, 0.3107]** (12 FLIP + 14 ONE-UNDEF); discordance among
  both-decidable **12/101 = 0.1188**; agreement only **89/115 = 0.7739**. Not a knife-edge — FLIP rate
  rises monotonically 1/101 (t=0.10) → 16/101 (t=0.80).
  **It is not read count.** Depth-MATCHED (B=200, seed fixed, 96 pairs), vs the same-library
  sampling-noise null:

  | group | D_cross (different library, matched depth) | D_within (null) |
  |---|---|---|
  | TP_LOST | 449/3600 = 0.1247 | 269/7200 = 0.0374 |
  | TP_RETAINED | 782/13200 = 0.0592 | 491/26400 = 0.0186 |
  | **FP** | **443/2400 = 0.1846** | **74/4800 = 0.0154** |
  | ALL | 1674/19200 = **0.0872** | 834/38400 = **0.0217** |

  Ratio ALL = **4.01×**; per-pair means give **FP 10.68×**. Paired sign test: cross > within in
  **17/18 non-tied pairs, two-sided p = 1.45e-04**.
  ⭐ **The instability is worst exactly where the route has to earn its keep.**
  **NAMED: BLMH × SLC38A6, the only gorilla FP the shipped rule rejects, FLIPS** (fib KEEP: tile
  0.0913, pile 0.5000, nd 2, 3 reads; tes KILL: tile 0.1905, pile 0.8750, nd 1, 10 reads).
  ⚠ **Concessions:** concentrated, not pervasive — 79/96 = 0.8229 of pairs have D_cross exactly 0 and
  the top 10 carry 87.0% of all disagreement; **the MRPS17 hub is stable** (KILL on all 4 edges in
  both arms), so this attack does **not** reach the flagship FP. **Attribution is unresolvable** —
  tissue, individual and date are confounded, and fib is the donor-matched arm. The **definitional**
  consequence holds either way; **the word "tissue" must never be attached to it.**
- **CLAIM "the rule becomes a function of sequencing effort": CONFIRMED but WEAK.** Rarefaction
  P(decision ≠ full-depth): 10% depth 0.0951 (fib) / 0.0888 (tes); 50% 0.0288 / 0.0291; 75% 0.0131 /
  0.0188. Pairs unstable at half depth **11/100 = 0.1100** (fib). **Direction matters more than
  magnitude:** `tile_frac` is a union over reads, hence **non-decreasing in depth** — the deeper arm
  gives the higher `tile_frac` in **57/64 = 0.8906**. Route A being kill-only, **a shallow library
  deletes edges a deep one keeps**: at 4× shallower, 5.4% (fib) / 4.1% (tes) fewer edges survive.
  **Undersequencing a real family splits it.** ⚠ `getreads.py` caps at **CAP=150 reads/node,
  LONGEST FIRST** — an engineering constant that is currently a **hidden parameter of `E_r`**.
- ⚠ **NEW DEFECT — "unservable ⟹ keep" is NOT a neutral convention.** The exempt stratum is not
  random: fib 0-read copies have median rep_len **2,370 bp (n=205)** vs expressed **3,061 bp
  (n=1,210)**, **perm p < 0.0001**; enriched in families of ≥4 copies **0.4439 vs 0.3182, p = 0.0005**
  (⚠ **non-monotone**: 0.1106 / 0.1523 / **0.2526** / 0.1489 by 2 / 3 / 4–6 / ≥7 copies — driven by
  the 4–6 band, not a clean gradient). Consequence: **15.53% of E_r edges are structurally immune to
  the read term in the honest arm**, and immunity is concentrated on short reps in mid-sized families
  — the population where a coverage call is hardest and where O1's false-merge risk lives.
  **Route A's scrutiny is conditioned on expression; its charity is biased.**

### 6.3 Attack — scale-invariance and a genuinely untested fix (**PARTIALLY lands; adds 1 REFUTED**)
Covered in §5.1 (rank invariance, pendant node, arithmetic ceiling) and §4.2 (over-extension
dose-response). The attack also **built and tested the fix its own defect implies**:

**R6 = `c_long × len_long/len_short`** (length-normalised coverage), kill-only on R0 edges. This is
genuinely outside R0–R5: R0 gates coverage-of-the-**shorter** at 0.50, R1/R2 raise coverage-of-the-
**longer**; R6 raises neither.

- HERC2 bottleneck **0.5838**; HSA FP floor **0.5079**. **7/19 FPs still sit below HERC2**, and
  **7/90 HSA TP edges sit below the FP floor.** **Ordering still fails.**
- Largest HERC2-safe point t=0.55 → HSA FP rejected 4/19 = 0.211 at the cost of **13/90 = 0.144** of
  HSA-node TP edges; t=0.60 rejects 11/19 but kills 22/90 = 0.244 **and drops HERC2** (panel-7 family
  ceiling 6/6 at t=0.55, **4/6 at t=0.60**). GGO: t=0.55 rejects 6/14 FP at the cost of
  17/101 = 0.168 TP. **C1 is violated at every point with any FP rejection.**
- Side finding: **under R6 the first family to die is NPIP (0.5608), not HERC2.**

### 6.4 Attack scoreboard

| attack | verdict | what it changed |
|---|---|---|
| circularity of read placement | **PARTIALLY** — headline limb FAILS | retracts C4; adds the missing negative control; corrects MRPS17 to a TE magnet; hard prohibition on an ADD rule |
| expression / library dependence | **PARTIALLY** — own headline REFUTED | refutes "unexpressed ⟹ unverifiable"; establishes 22.61% arm-dependence at 4.01× the null; kills the LEAD |
| scale-invariance | **PARTIALLY** | retracts the "4× deficit shrink"; adds the panel-7 arithmetic-ceiling TRAP; restates Route B's loss as over-extension; adds R6 as REFUTED |

**All three concur with `reject-all`. None breaks `impossibility_broken: false`.**

---

## 7. ⭐ RECOMMENDATION — and what it does to the DEFINITION

### 7.1 Recommendation: `reject-all`. `E_r` stays a pure sequence relation.

The user has said the definition needs to **WORK** and that theorems are expendable. That is the right
standard, and it is the standard under which this route dies: **a read-dependent `E_r` does not work.**
It is not rejected to protect P1. It is rejected because, measured on the frozen arms:

1. **It does not discriminate.** `tile_frac` AUC 0.1259 / 0.2308 / 0.0947 — significantly the **wrong
   direction** in both species and all three read arms, with a **non-circular** control passing.
2. **It rejects 0 novel FPs at any TP-safe operating point** (0/240 sweep cells; the only FP ever
   rejected safely is one R2@0.20 already rejects — and in testis it isn't even rejected).
3. **It violates C1 everywhere it does anything**: 25–30% of retained TP edges deleted at the pair
   level (Route A), 60–69% of TP edges and 78–81% of families deleted before any threshold (Route B).
4. **It buys instability**: 22.61% arm-dependence at 4.01× the sampling null, worst on the FP arm at
   ~11×; monotone in depth with kill-only semantics; a biased 15.53% exempt stratum; a hidden
   `CAP=150` parameter.
5. **It has a manufactured-signal floor** (2.19% / 1.12% at tile≥0.20, reaching 0.70) that no control
   in the pre-declaration would have caught.
6. **Its mechanism, where it appears to work, is misidentified** — MRPS17 is a TE magnet, and the #2
   magnet is a TP.

### 7.2 The three options, priced honestly

| option | what survives | what it costs | verdict |
|---|---|---|---|
| **(A) `E_r` stays a sequence relation** (recommended) | **P1, P2, P3, P4 and the λ-certificate all survive unchanged.** The shipped catalog, `lambda`/`cut_certified` in `families.tsv`, and the 788 passing tests are untouched. | Nothing is gained: the impossibility argument stands, HERC2 still splits before the first FP dies, MRPS17 remains a shipped FP. **The known false-merge rate stays 2/150 = 1.33% [0.37, 4.73].** | **ADOPT** |
| **(B) `E_r` becomes a sequence+read relation** | Nothing usable. **P1 survives only in a weakened form**: `E_r(A,B) → E_r(A,B \| L)` for a named library L, and seed-invariance must be restated as *invariance given fixed L* — it is no longer a property of the catalog alone. λ, P2, P3, P4 are all computed on a graph that is now library-conditional, so **the λ-certificate no longer certifies a property of the family, only of the family-under-L**. | Everything in §7.1 items 1–6. Plus: depth normalisation is **necessary but not sufficient** (4.01× survives exact depth matching); `CAP=150 / longest-first` must be removed or swept; the unservable convention must be scored **both ways**; a both-arms-agree rule abstains on **22.61%** of pairs and can only ever run on the **19.0%** of the scored core that has two libraries. | **REJECT** |
| **(C) unchanged, with read evidence as a non-binding FLAG** | Everything in (A). The flag never enters an edge decision, so P1–P4 and λ are untouched by construction. | Honest but near-useless: the flag's discrimination is AUC 0.13–0.23 **the wrong way**, so a flag built on `tile_frac` would mark real edges more often than repeat bridges. **Only defensible as a diagnostic annotation on the FP audit, never as a confidence score.** | **Optional, low value** |

### 7.3 The surviving LEAD, and why it must NOT be promoted
The secondary-record indicator **S** = "≥1 primary read at one node carries a secondary record inside
the other node's interval", already present in the `-N 50 --secondary=yes` BAMs with **no realignment**:

| arm | S+ TP_RETAINED | S+ TP_LOST | S+ FP |
|---|---|---|---|
| fib | 52/80 = 0.6500 | 7/21 = 0.3333 | **1/14 = 0.0714** |
| tes | 59/80 = 0.7375 | 7/21 = 0.3333 | **1/14 = 0.0714** |

Correct direction (Fisher TP_LOST vs FP p = 0.1078, **n.s.**), and **not** a function of the two reps
— minimap2's own genome-wide multimapping separates repeat bridges from paralogues better than any
forced pairwise realignment does; **the aligner declines to place FP-partner reads at all**, which is
the opposite of what the repeat-bridge story predicts.

**Three independent reasons not to promote it:**
1. **It BREAKS P1.** Whether a read at A carries a secondary record inside B is decided by
   genome-wide `-N 50` competition — a **third locus elsewhere in the genome** can take a record slot
   and flip S for the pair (A,B). S is a function of (A, B, reads, **THE REST OF THE GENOME**). It is
   not catalog-dependent in the R5 sense (it never reads the partition), but it is **not** a function
   of (A, B, reads at A and B), and any promotion must declare that.
2. **It fails C1 anyway** — as a keep-rule it deletes **26% (tes) / 35% (fib)** of retained TPs.
3. **It is read-derived**, so it inherits all of §6.2: ≥22.6% arm-dependence, monotonicity in depth
   under kill-only semantics, and the biased 15.5% exempt stratum.

### 7.4 P1 audit of the statistics actually tested
**Clean under P1** (function of (A, B, reads at A and B) only): `tile_frac`, `pile_conc`, `n_distinct`
(reads are selected by A's genomic interval — a read-**selection** device only — then aligned to rep B;
no other family enters); Route B's B1/B2 repairs (built per-copy from that copy's own reads).
`n_distinct` is degenerate (r = 0.9246 with rep length) but still local.
The 240-cell sweep, the HSA FP floor, and the HERC2 bottleneck are **calibration over the frozen arms,
not part of any rule** — they never enter an edge decision. *(Stated explicitly because a catalog-wide
floor used as a **live** threshold WOULD be the R5 failure mode.)*
**Breaks P1:** the indicator **S** (§7.3). **Minor:** `ROUTE_A_PRIMARY`'s "unservable ⟹ keep" is a
convention, and the HSA family-ceiling table scores an UNDEF edge as kept — **charitable to Route A,
and still does not save it.**

### 7.5 What this closes, and what it does not
**Closes:** read evidence as a source of separation for `E_r`. Both the direct route (tile the
partner) and the indirect route (repair the rep, re-apply coverage) are dead, and the third statistic
an attacker constructed (R6 length-normalised coverage) is dead too. **The impossibility argument was
NOT conditional on truncation** — this run tested that directly and the ordering did not cross.

**Does not close:** the impossibility argument's own **flagship framing**. HERC2 vs the HSA FP floor is
an n=1-vs-n=1 comparison of two extreme order statistics, across a **20× rep-length gap in a
non-shipped rep class**, on a **pendant edge** whose removal sends the bottleneck to 0.0000. The
conclusion is right; **the example must be re-stated on a length-comparable pair before it is shown to
the advisor**, or the first question kills it.

---

## 8. Register entries to draft

**REFUTED (5)**
1. Read tiling (`tile_frac`) as an `E_r` repeat-bridge discriminator — AUC 0.1259 (GGO tes) / 0.2308
   (GGO fib) / 0.0947 (HSA), p ≤ 0.0010, **wrong direction**, both species, three read arms,
   non-circular control passing. **High redo-risk.**
2. `pile_conc` (200 bp start concentration) — AUC 0.7109 / 0.7816, needs ≪0.5; the lost TPs pile
   **harder** than the repeat bridges.
3. Repaired-rep `c_long` (Route B) as a constraint-satisfying `E_r` rule — orderings fail to separate
   on **4/4** repaired substrates; repair moves `c_long` **down** (median −0.24 to −0.04).
4. **R6 = `c_long × len_long/len_short`** (length-normalised coverage) — genuinely untested in R0–R5;
   fails C1 at every FP-rejecting point; 7/90 and 10/101 TP edges still below the FP floor; NPIP dies
   before HERC2.
5. "A read-dependent `E_r` leaves unexpressed copies unverifiable" — **0/1,415 copies and 0/494
   families are silent in BOTH GGO arms** at 0 reads; the real cost is library-conditionality
   (205/1,415 = 0.1449 one-arm-silent).

**TRAP (8)**
6. `n_distinct` is a rep-length statistic in disguise — r = 0.9246 with `len(longer rep)`; median 1 in
   **every** group including TP_RETAINED.
7. "The reads are a third source" — median(`tile_frac_long` − `c_long_rep`) = +0.0003, 0/14 FP gain
   >0.10 **on the GGO scored pairs**; ⚠ but **not general** — HSA FP gains >0.10 in 6/19, and on a
   fresh panel the reads reach 0.07–0.10 **further** than rep-vs-rep. Any read statistic aligned to a
   rep largely re-derives that rep's alignment.
8. The **isoform-union denominator** inflates apparent rep truncation **0.920 → 0.478** on identical
   data. A rep is one linear path; use the longest-single-read denominator.
9. The **6.984 length-ratio confound is a human-curated-transcript artefact** — gorilla FPs are the
   **symmetric** pairs (median ratio 1.725).
10. A repaired/expanded rep substrate "rejects FPs" at t=0 by **deleting 60–69% of true edges** —
    **FP rejection must always be read CONDITIONAL on the edge surviving the substrate.** The loss
    mechanism is **over-extension**, with dose-response r(inflation, LOST) = +0.23 to +0.33 and
    r(inflation, Δc_long) = −0.51 to −0.63.
11. **HERC2 and 81% of the frozen scored core are HUMAN**; a gorilla-read repair is structurally
    unable to move them — "the threshold didn't move" would be **unfalsifiable, not negative**.
12. Any read-based `E_r` term requires a **cross-family / zero-rep-alignment negative control** — the
    manufactured-tiling floor is **2.19% (fib) / 1.12% (tes) at ≥0.20 and reaches 0.7026**. Related:
    **self-alignment positive controls for read statistics are circular**, because the read set is
    defined by its own primary alignment.
13. "No reads ⟹ keep the edge" is **not a neutral convention** — the exempt stratum is biased (shorter
    reps p<0.0001; enriched in 4–6-copy families p=0.0005, **non-monotone**) and covers
    **281/1,809 = 0.1553** of scored GGO pairs.
14. Comparing thresholds across **different statistics by ratio** is meaningless — `x → x^0.25`
    reproduces the reported "4× deficit shrink" while changing no ordering. **Use rank/quantile, or
    "TP edges that die before the first FP".** Related: the panel-7 substrate puts a **114,109 bp**
    quasi-genomic rep against ~2 kb spliced reps, so HERC2's 0.0325 is **65.0% of its arithmetic
    ceiling of 0.0500** while the FP floor pair has 7.9× more headroom.

**REFUTED — locality/stability (1)**
15. "The read statistic is local to (A,B) and therefore stable" — depth-matched cross-library
    disagreement **0.0872 vs 0.0217 same-library null = 4.01×** (paired sign test 17/18,
    p = 1.45e-04), worst on the FP arm at ~11×; **BLMH × SLC38A6, the only gorilla FP the shipped rule
    rejects, flips between arms.**

**DEAD-END with a LEAD (1)**
16. HERC2's split threshold **does** move 0.0325 → 0.1593 under reads — but its **rank inside the FP
    distribution is 0.0000 → 0.0000** and 0/19 HSA FPs are rejectable while HERC2 survives, under
    `c_long`, `tile_frac`, **or** R6. **LEAD:** the genome-wide secondary-record indicator S (already
    in the `-N 50` BAMs) is the only statistic pointing the right way — S+ FP 1/14 vs TP_LOST 7/21 in
    both arms — but it **breaks P1**, fails C1 as a keep-rule (deletes 26–35% of retained TPs), and
    inherits every read-dependence defect. **Do not promote.**

**HARD PROHIBITION (1)**
17. **Never promote read tiling from kill-only to an ADD rule.** At tile≥0.50 it adds **184
    cross-family edges (fib) / 117 (tes)** against 1,021–1,086 same-family R0 pairs on a 190-node
    panel — an **11–18% cross-family inflation** before scaling to the full 1,415 copies.

**Independent replication worth recording**
18. Reads whose best-scoring rep is in a **different family** carry **MAPQ ≥ 60 in 96.4% (fib) /
    96.1% (tes)** of cases, while MAPQ==0 is 0.0001/0.0038 over the same windows — an independent
    replication, on a fresh substrate, of the excision finding that **minimap2 is confidently wrong
    and MAPQ is uninformative**.
