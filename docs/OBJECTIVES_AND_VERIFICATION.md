# Thesis objectives and their verification checklist

> ⚠ **CATALOG PROVENANCE.** Figures quoted per **/494** (families) or **/1415** (copies) describe the
> **superseded** 2026-07-17 catalog, which **no invocation of the current binary reproduces** — it was
> built with `refine` on, a default removed on 2026-08-20.
> ⚠⚠ **SECOND BREAK, 2026-08-30: the 627-family / 2,019-copy figure is ALSO superseded.** The O1 node
> floor moved 3 → 2 (`NODE_MIN_READS`, ledger §6ac), so **no pre-08-30 catalog count is reproducible
> by the current binary** — set `RUSTLE_GATE_MIN_READS=3` to reproduce one. Re-measure before quoting:
> see [`NUMBERS.md`](NUMBERS.md) and [`o1_catalog_provenance.md`](o1_catalog_provenance.md).

⚠ **THE "CLOSED" STATUS BELOW IS STALE (it reads 2026-08-11).** Substantial work has landed since —
see **§ Answers to the advisor's standing questions** immediately below, and ledger §6ac–§6ba.

Status as of **2026-08-11 — FINAL PASS, CLOSED. The conclusions are stated below and are the thing to
quote; everything under them is the audit trail that earns them.** Three objectives, agreed with the
advisor on 2026-08-07.

> ⭐⭐ **THE ONE-LINE STATE OF THE THESIS.** **O1 is CONCLUDED with stated limits** — the definition is
> seed-free, P1 is a **theorem**, its properties are now measured at the **tier the binary ships**, and
> it has a falsifiable false-merge rate of **2/150 = 1.33% [0.37, 4.73]** with **measured power**; its
> binding constraint (blind node construction) is **named and not delivered**. **O2 is CONCLUDED with
> stated limits** — the objective **provably decomposes**, so the shipped per-read gate **is** the
> optimum and no joint estimator can beat it; K = 0 abstention is **entailed**, not chosen; accuracy is
> **synthetic-only**, and on real data what ships is a **commitment rate + an independently recomputed
> certificate** on a denominator outside the method's control. **O3 is PARKED** awaiting
> **same-individual** data (assembly + reads from one animal), which does not exist for this project and
> which no computation substitutes for; what O3 holds without it is stated in its section.
> ⚠⚠ **THIS PARKING REASON IS FALSE AS OF 2026-08-13 AND SUPERSEDED AGAIN ON 2026-08-14 — see the O3
> section.** The same-individual data **does** exist (mGorGor1 = KB3781 "Jim", and the matched IsoSeq is
> its own BioSample) and has been run: it returned **no reference-absent copy**. The follow-up question —
> *does that zero mean "nothing to find" or "our instrument found nothing"?* — was then put to the
> field's own **S1** standard (assembly-vs-assembly over three assemblies of that one animal) and the
> answer is **STILL AMBIGUOUS, but now measured and bounded**: no collapse-shaped deficit at **816
> evaluable loci** at a **0/817 [0.00, 0.47]%** false-positive floor, against a literature prior of only
> **0.47–0.94** collapsed loci in this stratum ⟹ **a perfect screen returns zero 39–63% of the time**.
> **O3's honest claim is a bounded negative about a compartment**, not a null result about a method —
> `/home/juanfra/winloci_scratch/o3_collapse/O3_COLLAPSE_VERDICT.md`, rows 3.15–3.19.

Full conclusions documents: `/home/juanfra/winloci_scratch/final/CONCLUSIONS.md` and
`/home/juanfra/winloci_scratch/close_o1o2/O1_O2_CONCLUDED.md` (O1 + O2 final, 2026-08-11 evening).
Test baseline behind every claim in this file: **771 passed / 0 failed / 11 ignored, exit 0**
(`CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo test --release --all-targets`,
captured to a file, `/home/juanfra/winloci_scratch/x4/test_final.log`).

> **2026-08-11 — the four unpinned claims (P5), settled.** Three of the four were **one `ls` away**, and
> the fourth is gone for good. **O2.4 restored to DONE** — the demand for a genome was a *category error*;
> 2752 is an exhaustive combinatorial enumeration with no substrate by construction, re-run today,
> exact. **O2.5 restored to DONE** — GORILLA, `GGO_mm.bam` (the hunt failed because `GGO.bam` is a
> *symlink* to it); 3599/3599 and 4030/4030 re-derived exactly today, and the 386 denominator recovered
> to the read as 311 + 75 (⚠ it is definition-sensitive — 1,692 under a looser rule — so the definition
> is now written into the row). **O2.7 restored to DONE** — the artifact was never missing;
> `o2_definitive.assignments.tsv` recounts to 104,147 rows / 43,239 reads / 74 families and max assigned
> p = **0.0006869**, with **0 of 72,212 assigned rows reaching α**. **O1.13's 13.6% / 0.237 WITHDRAWN** —
> no artifact survives, the one candidate file is 0 bytes, and **three** different purities were
> circulating for one experiment. Two disclosure debts also closed: the O2 corpus is now **split by
> species** (the 99.9990% is unmoved but **all 16 residuals are human**; the α plateau **does** move,
> pooled −0.384% = gorilla −0.303% / human −0.570%), and the **24/6 vs 27/7 panel gap is resolved with an
> exact bridge** — different sets, different dates, no residual.

> ⭐ **Added 2026-08-11 (X.4 / X.2 / X.8) — ONE TIER, and it is now a `diff` rather than a claim.**
> `refine_families_exon_sum` called `primary_seed_args()` unconditionally at two sites, so the
> 2026-08-07 sensitive-only default never reached O2: refine ran `-x asm20` @ **0.80** while O1 ran
> `-k 11 -w 5` @ **0.60**. Both sites now read one predicate (`er_sensitive_only`) and one selector
> (`er_primary_tier`), and both emit a **data-free `rule.tsv`** from one `er_rule_rows()`.
> **`diff O1.rule.tsv O2.refine.rule.tsv` is EMPTY on 5/5 refine calls**, and prints **5 lines** when the
> O2 arm is re-run at `RUSTLE_ER_SENSITIVE_ONLY=0` — so the instrument is shown able to fail.
> **88/88** refine calls across the 25-region `copy_assign` panel now log the sensitive tier; **0** asm20
> invocations remain anywhere. Regression: O1 control panel **75/75 byte-identical** (⚠ *by construction* —
> D1 made refine opt-in there, so it could not have moved) and the `copy_assign` panel **0/100 files
> changed** (⚠ *weak* — only **5 copysets** in 88 refine calls ever reach the edge rule). Tests **771
> passed / 0 failed / 11 ignored, exit 0** (was 768/0; +3 new X.4 tests).

> **What moved in this pass, in one line each.** A **full provenance audit** pinned genome / annotation /
> BAM / species behind every headline — verdict: **no number is produced by a cross-species or
> cross-annotation computation**, but four *disclosure* defects and one withdrawn sentence came out of it
> (new section below) · `O3.14` retired the denominator the detector controlled and **every O3 rate is
> restated on a fixed denominator** · `O3.8`'s reported "cost" is **retracted — the sign reverses**, the
> envelope fix was an improvement · `O2.10` was **measured for the first time**: the EM changes **zero**
> of the gate's decisions on reads that carry evidence, so the shipped gate is **empirically adequate**
> and that is a result · `O1.10` moves from "n = 2 with no power" to **n = 0 — precision is not
> measurable with annotation-derived negatives**, and the honest proxy is stated instead.
>
> ⭐ **Added 2026-08-11.** `O1.10` moves again, and **upward**: the *pair* route stays dead, but changing
> the **unit** from a labelled pair to a gene-tight **window** and the **truth** from gene groups to
> **coordinates** yields a negative stratum with **demonstrated power** (1,630 eligible windows, 150 run,
> 108 with a real merge opportunity, the rule firing 28 times) and O1's first false-merge rate —
> **2/150 = 1.33% [0.37, 4.73]**, both failures certified without an aligner or a gene label. ⚠ It
> **contradicts** the earlier c-sweep reading: the two failures sit at coverage **0.5571** and **0.9660**,
> so `c = 0.51` fixes neither (rows 1.7 / 1.10 / 1.15). · A sixth provenance defect, **P6**, records that
> `bench/soto/known_family_bench.py` builds a **gorilla** benchmark whose roster *and* copy-number
> denominator are **human** — no number in this document rests on it, but the memory record does, and
> those completeness figures are **withdrawn as gorilla results**. ⚠ The "no cross-species computation"
> verdict above is therefore true **of this document**, not of the project as a whole. · **O3.13 closes
> the gap between the fixed statistic and the published tables**: `detector.rate()` and
> `detector.famboot()` now **raise**, one shared `outcome()` backs every report, all 20 call sites in
> `report_det.py` plus `final_table.py`'s inline reimplementation are repointed, and every O3 table is
> regenerated with a fixed-denominator header, the four outcomes per cell, and **ARM and LOCUS split
> into separate blocks**. No O3 number changed — they were already computed on the fixed denominator in
> `O3_FIXED.txt`; what changed is that the condemned convention can no longer be produced by accident.
>
> ⚠ **A correction to the previous edition of this file.** It presented "O3 TPR 41.1 → 37.1, O3 locus
> 40.4 → 37.5" as a deliverable *fall*. **Both movements were artefacts of a denominator the detector
> itself controlled** (row 3.14): on the fixed denominator the same fix, same data, same resamples, moves
> the arm yield **29.48% → 31.09%, +1.62 pp [+0.27, +3.20]** — strictly more true detections at a
> **bit-identical** matched FPR. What *is* genuinely lower is the level of every O3 rate, because
> abstentions now sit in their own denominator (arm yield **31.09%**, not the called-rate 37.14%). The
> retired O2 headline (99.28%) remains withdrawn.

Legend: **DONE** = measured and survived adversarial review · **PARTIAL** = established with a stated
limit · **TODO** = not attempted or not finished · **REFUTED** = tried and found false.

> ⚠ A subtask is only DONE if its evidence is (i) reproducible today and (ii) not derived from the
> pipeline being evaluated. Six metrics have died here for failing (ii) — either truth containing the
> prediction by construction, or a denominator conditioned on the prediction.

---


## INDEX

> **Index.** 9 sections; this is the map. **The titles carry the verdicts** — no tag is derived
> here. ⚠ In `o1_ledger.md` an earlier auto-derived verdict tag scored **11/22 = 50%** against
> sections whose outcome was known first-hand, so tags were removed rather than shipped. Search a
> heading to jump.


- Answers to the advisor's standing questions (current as of 2026-09-01)
- The three objectives
- THE THREE CONCLUSIONS — final, 2026-08-11
- Provenance and substrate hygiene — read this before quoting any number below
- O1 — Define
- O2 — Assign
- O3 — Detect
- Cross-cutting
- What is blocking a conclusion — FINAL, revised 2026-08-11
- Final note — how to use this file
