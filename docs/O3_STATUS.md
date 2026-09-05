# O3 — reference-absent copies: the evidence, and what it does and does not support

**2026-09-03.** Built from a 3-agent inventory of the full ledger, the negative-results register, memory,
and code reachability. Ledger §6dk–§6dq and the sections cited below.

---

## ⭐⭐ THE HEADLINE — AND IT IS NOT THE OBVIOUS ONE

> **O3 is NOT refuted, and "well tested but lacking support" is NOT defensible.**
> What has been refuted is a set of NAMED DETECTORS. The phenomenon itself has **two DNA-level positive
> instances**, and **no shipped binary can emit a reference-absent copy at defaults** — every admission
> route is behind a default-false flag.
>
> ⭐**The defensible sentence: "O3's RNA detectors are refuted, vacuous or unvalidated; O3 itself is
> UNTESTED IN THE REGIME WHERE IT WOULD OCCUR."**

**Absence of output is a consequence of the gating, not a measurement.**

---

## 1. THE THREE-WAY DECOMPOSITION — "O3" NAMES THREE SEPARABLE THINGS

| | status |
|---|---|
| **specific RNA statistics** | **genuinely REFUTED** — real measurements went against them |
| **the headline zeros** | **VACUOUS / BLOCKED / OFF-TARGET BY CONSTRUCTION** — not refutations |
| **the phenomenon** | **DEMONSTRATED at DNA level, never tested on RNA** |

### 1a. Genuinely refuted — named statistics
- **k+1 MEC** — four independent grounds: objective degenerate under missing data (GWFAM115 k=6 gives
  MEC = 0 at **15.3% of cells scored** vs a genuine k=5 fit of 1,076 at 100%); the ONLY clean negative had
  the LARGEST k→k+1 drop; the k=n control failed (ARI **0.7621**, **0.2718** against a pre-registered 0.80);
  the panel holds **≈0.2 expected true positives** and *"zero instances of the object under test"* (§6ai).
- **Between-copy RNA depth** — signal **1.27x** against between-copy expression variance of **109.96x**
  (GWFAM225) and 55.75x (GWFAM163). Buried two orders of magnitude down (§6s).
- **O2 near-tie ambiguity as a proxy** — median depth skew **1.49 vs 1.56** (identical) while max identity
  differs **0.9949 vs 0.7720**. It tracks copy SIMILARITY, not a missing copy (§6s).
- **WGS k-mer multiplicity** — paralogues at median identity 0.8234 share **~1.7%** of their 21-mers, so
  copy number reads as ~1 regardless of family size (§4u).
- **`collapse_gate`** — fires on `EEF1A1`, reporting chi(H) = 7 for a **one-copy** locus
  (`collapse_gate.rs:24-26`).

### 1b. Vacuous or blocked — the instrument cannot move
- ⭐**The shipped transcript-side detector (FARDIV ∧ FARCLIP) fires 0 of 915 loci — and CANNOT fire.**
  FARCLIP's median is **0.0006** against a gate of **0.05**. The source's own words: the mini-reference
  *"did not merely flatter the rate, it MANUFACTURED the signal"*; the 0/915 *"says nothing about the
  phenomenon"* (§6t).
- **The unmapped-read route** — detection power in the collapsible class is **1/35 = 0.0286
  [0.0007, 0.1492]**, and **0/26 [0, 0.1323]** at coverage ≥ 0.8. *"Vacuous there, not merely loose."*
  ⚠**And the collapsible class IS O3's target by definition** (§4y).
- **`absent_copy` gate 5** — compares a HOST-DERIVED reconstruction against the HOST: **0.9961** where the
  true copy-to-copy identity is **0.9754**, i.e. it reproduces **~18%** of the real divergence, and only
  **2 of 19** families land in the recoverable window. It rejects **172/320 = 53.8%** of candidates on
  "≥98% remap identity". **Near-unpassable by construction.**
- **`em_copy_assign`** — its `logl` is an immutable reference-built input, so it **structurally cannot form
  a novel consensus**.
- **S2 divergence-mixture** — power **0.4500** above 0.01 divergence but **0.0588** below, while **45.78%**
  of true positives lie below it. Near-blind in exactly the near-identical class that collapses (§4y).

### 1c. The phenomenon — demonstrated, at DNA level, in one animal
- **GOLGA6L7** (§6u, reproduced §6dk): MAT 2 / PAT 3 units in one array; extra unit **6,035 bp, cov 1.000,
  id 0.9673** against siblings 0.9670/0.9677, regular **40,819 / 40,808 bp** period; **9/9 autosomal
  controls concordant**. ⚠**It is the POSITIVE CONTROL, not a finding** — the primary is PAT-derived there,
  so the reference already contains the extra unit.
- ⭐**`LOC134757045`** (§6do) — uncharacterized lncRNA, chr14, contig pair 17/17 CLEAN. mat **5** / pat **3**;
  three units match one-to-one (0.9957↔0.9957, 0.9780↔0.9780, 0.9735↔0.9739); **two maternal units at
  2,298 bp, cov 1.000, id 0.9961 have NO paternal counterpart**, and pat carries **3 passing / 3 UNFILTERED**
  records — nothing hidden. The primary is PAT-derived ⟹ **those two copies are genuinely ABSENT from the
  reference.** ⚠**It has never been put in front of an RNA detector.**

---

## 2. ⛔ NO SHIPPED BINARY CAN EMIT A REFERENCE-ABSENT COPY AT DEFAULTS

All four admission routes are gated default-false: `--absent-copies`, `--vg-realign`, `--collapse-gate`,
`--tied-seed` (`denovo_pipeline.rs:2235/:2466/:2482`, `gw_family_catalog.rs:885/896/907`). `absent_copy.rs`
is 856 lines behind a flag that has never been on in a production run; its certificate leg (`linearize.rs`)
is behind another.
⟹ ⭐**The pipeline has never attempted the detection. There is no null result to interpret.**

---

## 3. ⚠⚠ THE DO-NOT-QUOTE LIST — every one of these has been retracted or is being mis-cited

| ⛔ do not say | why |
|---|---|
| *"289 candidates, 0 admitted"* | **RETRACTED as over-generalised** (§6dc). Correct: *"0 of 289 **in the 12-family MEC panel**"*. On 73 real gorilla families it is **8/328 = 2.4%**, 3 of the 8 in a family the register voids. ⚠It was still wrong at `O3_STATUS.md:50` and `MEMORY_DIGEST.md:39`. |
| *"TPR 0.7965 / FPR 0.0000 / AUC 0.8924"* | **RETRACTED (§6t)** — an ORACLE statistic requiring BOTH references, i.e. already knowing where the copy was removed. Only the mechanism survives (1.2744 vs control 1.0000). |
| *"8/101 = 7.92% haplotype CNV rate"* | §6w says verbatim **"Do not quote 7.92%"** — half were fragment artefacts. Defensible: **≤ 4/101 = 0.0396 [0.0155, 0.0977]**, tandem-confirmed **1/101**. |
| *"8/9/9 vs 5/6/8 haplotype deficit"* | Does not reproduce at either `-p`; deficits shrink 3,3,1 → ~1,2,0. **Direction survives, magnitude does not.** |
| *"2.78% mat-vs-pat difference rate"* | The 08-19 run is UNINFORMATIVE by its own pre-registration: span-matched null **0.1512** vs signal **0.0278** — a floor **5.4x the signal**. |
| *"the DNA-side 0/817 is an independent second zero"* | Explicitly forbidden as a repeat of the denominator-error pattern. ⚠`MEMORY.md` and `MEMORY_DIGEST.md:129` still title it *"AN INDEPENDENT ZERO"* — **the 08-15 caution is the later record and it governs.** |
| *"the NPIP chr16:80,195,301 novel gene"* | **It is NPIPB15**, to the base. It scored as unannotated only because the comparison used a wrong coordinate (§6cc → §6cd). |
| *"multi-copy families do not differ between individuals"* | The 0/226,615 zero is on **99 housekeeping loci** selected by the very filter that makes the test valid; **NPIP, GOLGA6 and TBC1D3 can never pass it** (§6dq). |
| *"GOLGA6L7 is a reference-absent finding"* | Its extra unit is on the haplotype the reference came FROM. It is the positive control. |

---

## 4. THE ONE WELL-POWERED ZERO, AND WHY IT DOES NOT GENERALISE

§6dq: **0 positions with ≥3 alleles across 226,615 positions**, 99 loci, two animals — no collapsed extra
copy in either, and no reference collapse. Genuinely well powered (~2,290 positions per locus).
⛔**But the 99 loci are `NCOR1`, `RPL10`, `SUMO2`, `HMGN2`, `XPO6`, `DAP3`, `GSPT1`…** — housekeeping genes,
because **requiring expression in BOTH tissues (necessary to control the expression confound) selects for
constitutively expressed genes**, which are the least likely to be recently duplicated. **49 of 3,310
families were testable, and they are the wrong 49.**
⭐**A negative on a population the hypothesis excludes is not evidence against the hypothesis.**

---

## 5. WHAT WOULD ACTUALLY SETTLE IT

1. ⭐**MATCHED-TISSUE IsoSeq FROM A SECOND GORILLA.** The blocker is comparing two *tissues* from two
   animals: controlling expression costs exactly the gene class the hypothesis concerns. With tissue held
   constant, "expressed in both" stops selecting for housekeeping genes and NPIP/GOLGA6/TBC1D3 become
   testable. **This is an acquisition, not an analysis.**
2. ⭐**PUT `LOC134757045` IN FRONT OF THE RNA DETECTOR.** Two copies known absent from the reference, in the
   animal the fibroblast RNA came from. A detector that misses them is falsified; one that finds them has
   passed a real test. ⚠It is a lncRNA — **if silent in fibroblast the test is UNINFORMATIVE, not negative.**
3. **Fix gate 5 before quoting any admission rate.** It measures a reconstruction against its own source.
   Until then no number from `absent_copy` means anything.
