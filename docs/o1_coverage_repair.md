# O1 — the coverage-denominator repair: can we fix the bugs, and is O1 defensible?

**Status 2026-08-15.** Supersedes the "strongest objection still open" paragraph of
[`o1_error_case_census.md`](o1_error_case_census.md) §5.

> **2026-08-16 follow-up:** a new threshold-free axis, transcript alignment
> orientation, rejects 29/74 frozen false edges while preserving the named families
> after repairing the historical panel's strand bug. It is an RNA-only hard-rule
> candidate now implemented behind `--rna-forward-only`, pending an end-to-end
> GGO/HSA catalog run before any default change. See
> [`o1_false_positive_rules.md`](o1_false_positive_rules.md). The negative conclusions
> below about coverage, length, soft-mask, and promiscuity thresholds are unchanged.

**Short answer, in one line:**

> **The bugs are fixed and guarded. The repeat-bridge failure class is real, named, and reproducible —
> but it is NOT fixable by any coverage-denominator threshold that O1 can afford: every threshold that
> kills a single false edge already breaks HERC2, and the best constraint-satisfying candidate buys
> 2 false catalog merges for 21 destroyed families. The deep-paralog case is OUT OF SCOPE for a
> sequence-based rule, and the census's headline objection built on it is withdrawn. O1 stands as a
> definition; what changes is that the repeat-bridge class must be EMITTED as a per-edge flag, not
> excluded by the membership rule.**

---

## 1. ⚠ LIMITS FIRST

Read this section before quoting any number below it.

| limit | statement |
|---|---|
| **T8 — nothing here is pipeline-confirmed** | Every sweep number is **offline PAF re-derivation** on frozen arms through the `rustlib.paf_pairs` / `er_edges` mirror. It reproduces R0 exactly (74/74 FP, 9,424/9,424 TP core; independent re-derivation matched the sweep's panel-7 figure 325/392 exactly). It is still a **hypothesis generator**. **Nothing has been run through `gw_family_catalog`.** See §5 for the run that would confirm. |
| **FP arm is small, and smaller than it looks** | 78 rows / **62 distinct gene pairs** / **27 independent mechanism components** (24 in the scored core). 30/74 scored rows = 40.5% are ONE component; **14/24 components are singletons**. The pair-level binomial CI is an independence error: cluster-bootstrapped over components, 24/74 = 0.3243 has CI **[0.1429, 0.4578]**, not [0.2286, 0.4373]. Quote the mechanism unit alongside the pair unit, always. |
| **Power** | Paired McNemar on 74 rows separates candidates differing by ≳ 6–8 discordant pairs ≈ **10 pp of FP-kill rate**. The sweep can answer *"does changing the denominator kill the repeat-bridge class without destroying short-gene families"*. It **cannot** answer *"is c = 0.55 better than c = 0.60"*. |
| **Species never pooled** | FP arm: HSA 52 rows, GGO-involving 26. The headline 24/74 is earned almost entirely on the **human curated-transcript** stratum, which ships nothing. |
| **What the three attacks broke** | (a) the binomial CI (§4); (b) the "changes 2 real catalog merges" headline — it is **0 at c=0.19** (§4); (c) the "≥95% of family ceiling in BOTH species" constraint — computed on census catalogs only; **with panel-7 included, HERC2 fails at every c_long ≥ 0.034, below the first FP kill at 0.05** (§4, §5); (d) the sweep's family metric itself — the TP-pair graph is not the shipped family, and the true structural cost is ~2× what was reported (§5). |
| **What the attacks FAILED to break** | The "threshold was tuned on the cases it is scored against" objection is **measurably false** — leave-one-mechanism-out over 24 folds picks the same argmax **24/24**, held-out FP rejection 0.3243 = in-sample; half-split cluster CV over 2,000 reps gives optimism **−0.0022**. Do not run that objection. The P1/R5 finding (§5) is independent of the FP arm and stands. |
| **`fp_risk_pool.tsv` is not the FP arm** | 61,587 unrelated pairs with a record but no shipped edge. Relevant only to candidates that ADD edges. |
| **`q915_exon.fa` is NOT a valid E_r substrate** | Same intervals as the shipped copies (915/915) but whole-locus exon-*block* sequence, **median 4.618× longer** (0/915 identical lengths). It recovers 159/2,165 = 7.3% of within-family pairs where the shipped reps give 1,334/2,165 = 61.6%. Its pairs are parked as `G6_Q915_WRONG_SUBSTRATE`. **Do not use it as a bulk arm.** |

**Frozen arms:** `/home/juanfra/winloci_scratch/o1_fix/arms/`, `SHA256SUMS.txt` sha256
`a0aa73955c39e33e36718a45d14a26d071935d75ddddf07dc8506b214136efb7`
(`fp_set.tsv` 78 · `tp_set.tsv` 20,528 · `fp_risk_pool.tsv` 61,587 · `grey_not_scored.tsv` 1,522 ·
`excluded_deep_paralog.tsv` 2 · `records/*.paf` 11 files / 64,080 records · `SPEC.txt`).
Builders: `o1_fix/build_arms.py`, `truth.py`, `finalize.py`. Sweep: `o1_fix/sweep2.py`.
Attacks: `o1_fix/attack_overfit.py`, `attack_tp.py`, `struct_attack{,2,3,4}.py`, `p2margin.py`.

---

## 2. The `mkreps` bug — fixed, guarded, blast radius mapped

### 2.1 The bug and the fix

`bedtools getfasta -s` reverse-complements **each BED record individually**. Concatenating the records
back in file order therefore gave every **minus-strand** gene sense-correct exons in **reversed order**.

The fix in `/home/juanfra/winloci_scratch/o1_errorcensus/mkreps.py`: extract every exon on the **plus**
strand (no `-s`), concatenate ascending-genomic, then reverse-complement the **whole** sequence once.

```python
rep = revcomp(plus) if st == '-' else plus     # ONE flip, on the WHOLE sequence
check_rep(g, st, parts[g], rep)                # guard, unconditional
```

BED names carry the exon start (`{gene}|{start}`) so records stay addressable regardless of emission order.

### 2.2 The guard — two layers, both verified to FAIL on the old code

* **Layer 1 — `check_rep()`, on every gene on every invocation.** Asserts
  `revcomp(rep) == plus-strand genomic-order concatenation`, plus an explicit assert that `rep` is not
  the per-exon-revcomp string. Negative control: `ATP4A (−, 22 exons)` → **GUARD FIRES**
  (`"exons are in the WRONG ORDER"`); `ATP1A1 (+, 23 exons)` → guard silent, correctly (buggy == good
  on the plus strand).
* **Layer 2 — `mkreps.py --selftest`.** Rebuilds ATP4A(−), ATP1A3(−) and ATP1A1(+, control) and
  requires each to align back to its **own** locus as one colinear spliced chain at cov ≥ 0.95 /
  id ≥ 0.99. Fixed ATP4A **cov 1.0000**; buggy ATP4A **cov 0.1647** → hard fail.
  **Re-run 2026-08-15 for this document: all three PASS (cov 1.0000 / id 1.0000).**

### 2.3 Blast radius — THREE buggy generators, not one

Exhaustive sweep: only 4 `getfasta -s` sites exist in the tree, all in `o1_errorcensus`. A separate
search for the Python-side disguise (per-exon revcomp inside a loop) returned **0 hits**.

| script | status |
|---|---|
| `o1_errorcensus/mkreps.py` | **BUGGY → FIXED** |
| `o1_errorcensus/verdict/mkclean.py` | **BUGGY → FIXED**, 41/77 genes corrupted |
| `o1_errorcensus/verdict/mkggo.py` | **BUGGY → FIXED**, 8/18 genes corrupted |
| `o1_errorcensus/verdict/addspan.py` | uses `-s` but max 1 BED line per gene ⟹ no concatenation ⟹ **not affected** |
| `mkreps_fixed.py`, `cls3/mkreps_ggo.py` | already correct |
| `o1_blind/nodegap.py`, `gw/gw_fasta.py`, `o1_bridge/build_edges.py`, `seedfam/.../rna_cmp.sh` | no `-s`, all-plus-strand, never revcomp ⟹ **not affected** |

⚠ **The fix was written but never propagated**: `mkreps_fixed.py` at 20:09, `mkclean.py` ran **21:00 —
51 minutes later, still buggy** — and `mkclean.py` is the most consequential of the three.

### 2.4 Corrected numbers

Rebuild diff was surgical: **41 genes changed, all 41 minus-strand, 0 plus-strand.**

| number | reported | **corrected** |
|---|---|---|
| ATP1A1 × ATP4A | "no edge, cov 0.0443" | **EDGE — cov 0.5689, id 0.7163** (clears the shipped floors) |
| E_r edges, `clean.fa` (77 genes) | 81 | **127** (+48 falsely absent, −2 spurious) |
| E_r edges, `all_clean.fa` (88 seqs) | 104 | **152** (+50, −2) |
| `pairs.json` case-pair verdicts | — | **18/111 flip, ALL `no → EDGE`, zero the other way** |
| CASE3 ANKRD18B ~ ANKRD20A1 | no edge, cov 0.2302 | **EDGE, cov 0.6954** (id 0.7571 → 0.7574) |
| CASE4 FAR2P1 ~ CYP4F30P | no edge, cov 0.3560 | **EDGE, cov 0.7747** (id 0.9913 → 0.9862) |
| `verdict/rates.txt` — all four rates | 4/8, 1/8, 3/8, 7/126 | **INVALID, do not quote** |

Other flips: ANKRD18A~ANKRD18B 0.2206→**0.9870**; LINC01297~PSLNR 0.4693→**1.0000**;
TUBA1A~TUBA3C 0.4425→**0.8902**; TUBA4A~TUBA3C 0.4425→**0.8856**; ZNF670~ZNF669 0.4211→**0.5188**.
Verified **unchanged**: CASE1 TUBA3C~MZT2A (no record either way), CASE2 CXADRP2~POTEB2 (EDGE,
byte-identical), CASE5 BCRP2~POM121L8P (no edge, byte-identical).

⭐ **Diagnostic signature: identity essentially unchanged, coverage roughly doubles.** Each exon aligns
fine; only one exon-run can chain. **Any past "low coverage, high identity, no edge" verdict on a
minus-strand gene is suspect.**

⚠ **The correction moves the census in the WRONG direction for O1.** All 18 flips are `no → EDGE`.
The census's 30 definitional failures were a floor before; the corrected reps **raise the edge count
by ~57%** on the case set. `rates.txt` must be redone: 4/8, 1/8, 3/8 are human judgement and must be
re-classified **by hand** (CASE3 and CASE4 both have their central claim reversed). For 7/126, the
author's exact rule could not be reproduced, but under the natural "top-hit gene ≠ label gene" rule
**12 of 126 nodes were falsely refuted** — since 12 > 7, **7/126 cannot survive**. Falsely-refuted
nodes concentrate in the census-dominant families GWFAM230, GWFAM164, GWFAM149.

⚠ **Quoting hazard:** `0.0556` also appears in memory as the O1 **false-omission** headline (9/162).
That is a different number from a different pipeline (ARM 3 excision), does not use `mkreps`, and is
**unaffected**.

### 2.5 Shipped Rust and `bench/` — CLEAN (verified, not assumed)

* **No `getfasta` anywhere in `Rustle/src/`.**
* All Rust multi-segment builders use concat-then-revcomp-whole:
  `src/rustle/vg_family/denovo_assemble.rs:1002-1013`, `denovo_pipeline.rs:4604-4629`
  (`cothread_locus_geometry`), `:4667-4692` (`union_locus_reps`).
* `family_graph.rs:399`, `from_genome.rs:55/71/95` — single interval, plus strand, no revcomp. Correct.
* `bench/extract_gene_reps.py:80-88` — concat-then-revcomp-whole. Correct.
* `bench/copy_specific_junctions.py:77`, `bench/extract_intron_chains.py:66` — explicit
  `ordered = exons if strand == "+" else exons[::-1]`. Correct.
* `bench/crossspecies/{seed_cds.sh:31, isoseq_probe.sh:36, graph_vs_graph.sh:43}`,
  `seedfam/dnapr/meeting/rna_cmp.sh:16` — `getfasta` **without `-s`**, never revcomp; an
  all-plus-strand convention that is self-consistent (minimap2 aligns both strands). Not the bug.

### 2.6 Files changed / regenerated / still stale

* **Changed:** `o1_errorcensus/mkreps.py`, `o1_errorcensus/verdict/mkclean.py`,
  `o1_errorcensus/verdict/mkggo.py`.
* **Regenerated correctly:** `verdict/{clean.fa, clean.bed, all_clean.fa, all_clean.er.paf,
  node_vs_clean.er.paf, ggo_ann.fa, x.bed}`.
* ⚠ **Still stale — do not read:** `verdict/rates.txt`, `verdict/fam_nodes.er.paf`, and the pre-20:09
  root artifacts `atp.fa`, `case_all.fa`, `golga_ex.fa` (superseded by the `*_fx` versions; the
  un-suffixed originals are still on disk).

---

## 3. The two failure modes, kept apart

The census's §5 verdict conflated them. They have different answers.

### (i) REPEAT-BRIDGE — real, named, reproducible, and the coverage denominator IS the mechanism

A ~970 bp dispersed Alu-like element in `MRPS17[726-1695]` bridges unrelated genes. Under the shipped
rule (coverage of the **shorter**) it passes; under coverage of the **longer** it collapses by an order
of magnitude, not marginally. Reproduced against the frozen mirror on curated transcripts:

| edge | cov(shorter) | cov(longer) |
|---|---|---|
| MRPS17 × MDM2 | 0.561 PASS | **0.0792** |
| MRPS17 × GREB1 | 0.582 PASS | **0.0689** |
| MRPS17 × TRAPPC2 | 0.561 PASS | **0.3442** |

⚠ The census's ATP row (qlen 3938 / tlen 3608, cov 0.551) is the **NODE** pair, not the curated one;
curated `ref_ATP1A1 × ref_ATP4A` is 3764/3721, covmin 0.5689, **covmax 0.5624**.

**But the class is much narrower than the row count suggests.** Of 78 FP rows, 52 are repeat-bridge —
yet they are only **8 of 27 mechanism components**: MRPS17 (17 gene pairs), TTC6/DNAH14 (6),
MCFD2/ENTPD1 (5), SDHAF4 (4), ZNF480 (4), TDRD5 (3), DNAJA1 (2), SLC38A6 (2). Two hubs are new and
were not in the census — **SDHAF4** (aln 567 bp, covmin 0.51, covmax 0.026–0.19) and **ZNF480**
(aln ~2,600 bp, covmax 0.13–0.40). ZNF480 proves the mechanism is **not confined to ~1 kb elements**,
so "short repeat" is not a sufficient description of the class.

`len_shorter ≤ 2 kb`: 55/78 = 0.7051 [0.5962, 0.7948]. **Only 42/78 are merges the shipped catalog
actually makes** — the rest exist on curated transcripts only.

### (ii) DEEP-PARALOG — **OUT OF SCOPE, and the census's headline objection dies with it**

ATP1A1 and ATP4A **share an mmseqs protein cluster** (`protclust_cluster.tsv`: both → `ATP1A1`).
HGNC splits them on **function** ("Na+/K+ transporting" vs "H+/K+ transporting"); by sequence they are
homologs. Both rows are moved to `excluded_deep_paralog.tsv`.

⭐ ATP1A1 × ATP4A was the **only** case behind the census's
*"no value of τ or c orders these two correctly"*. **That objection is not supported by a sequence-level
truth and is withdrawn.** A sequence-based membership rule cannot be faulted for merging two sequences
that a protein-clustering oracle also merges; the disagreement is with a functional nomenclature, not
with homology.

GFPT1 × GFPT2 sits in TP on three substrates (curated tx covmin 0.9475 / covmax 0.929; GGO node covmax
0.4312; HSA node covmax 0.4292) — **the node-level numbers are the hard case, not the curated ones**,
which is a node-construction observation, not a definitional one.

**⟹ The honest position: the repeat-bridge class is the only definitional failure mode on the table,
and the deep-paralog case is out of scope for a sequence-based rule.**

---

## 4. The candidate sweep — FP rejection NEVER without its TP cost

**R0 (incumbent, shipped):** one PAF record with `identity(1−de) ≥ 0.60` **AND**
`aligned-span / min(qlen,tlen) ≥ 0.50`, floors per record.
Verified on the frozen arms: **FP kept 74/74** (24 components), **TP core kept 9,424/9,424**
(T1 HSA curated tx 5,806 · T3 GGO node 1,795 · T3 HSA node 1,416 · T4 panel-7 392 · T2 panel POS 15),
risk pool **0/61,587**. Family ceiling (TP-pair graph one component): **GGO 370/375 = 0.9867,
HSA 270/274 = 0.9854** — 5 GGO / 4 HSA families are already multi-component under R0.
**Every candidate must beat R0, not merely differ from it.**

### Candidates

| id | rule |
|---|---|
| **R1** | coverage-of-the-**longer** alone, `aligned / max(qlen,tlen) ≥ c_long` (drops the cov_min floor) |
| **R2** | R0 **plus** a third per-record condition `aligned / max(qlen,tlen) ≥ c_long`, same record — a strict restriction of R0 |
| **R3** | absolute aligned-base floor `B` |
| **R4** | soft-mask-fraction of the aligned blocks |
| **R5** | block promiscuity (partner count of the aligned block over the catalog) |

**R1 is DOMINATED and must not be adopted in any form.** It is not a restriction: it **adds**
335/61,587 = 0.0054 unrelated risk-pool edges at c_long 0.15 (94 at 0.20, 8 at 0.30) plus 4,278/11,104
T6 edges at 0.15, while attaining FP rejection and TP retention identical to R2 to within 1–2 pairs at
every threshold ≥ 0.15.

### The trade-off table (unit = pair; FP denominator 74 scored-core, TP denominator 9,424 core)

| candidate | FP rejected /74 | (repeat-bridge /52) | TP core retained | real catalog merges changed /39 |
|---|---|---|---|---|
| R0 (incumbent) | 0 | 0 | 9,424 = 1.0000 | 0 |
| R2 @ 0.05 | 2 | 2 | 9,180 = 0.9741 | — |
| R2 @ 0.10 | 10 | 9 | 8,806 = 0.9344 | — |
| R2 @ 0.15 | 15 | 14 | 8,595 = 0.9120 | — |
| R2 @ 0.19 | 19 | — | 8,349 = 0.8859 | **0** |
| **R2 @ 0.20** | **24 = 0.3243** | **17** | **8,300 = 0.8807** | **2 = 0.0513** |
| R2 @ 0.21 | 28 | — | 8,226 = 0.8729 | 3 |
| R2 @ 0.25 | 32 | 24 | 7,916 = 0.8400 | — |
| R2 @ 0.30 | 41 | 30 | 7,514 = 0.7973 | — |
| R2 @ 0.40 | 58 | 43 | 6,374 = 0.6763 | — |
| R3 @ B=600 | 8 | 4 | 7,741 = 0.8214 | — |
| R3 @ B=1000 | 44 | 34 | 6,457 = 0.6852 | — |
| R3 @ B=2000 | 67 | 46 | 2,398 = 0.2545 | — |
| **null: `len_longer/len_shorter ≥ 3`** | **20 = 0.2703** | — | (worse on TP) | — |

R2 @ 0.20 is the best FP rejection reachable while keeping ≥ 95% of R0's **arm** family ceiling in both
species, and R2 adds **0/61,587** risk-pool and **0/11,104** T6 edges at every threshold.
McNemar (paired, n=74): R2@0.20 vs R2@0.15 b=9/c=0 p=0.0039; vs R3@600 b=17/c=1 p=1.5e-4;
vs R2@0.25 b=0/c=8 p=0.0078.

### ⚠ Five things that disqualify R2 @ 0.20 as an ADOPTED rule

1. **Its one real-world consequence evaporates at the 4th decimal.** 17/74 = 23% of FP rows sit in
   `[0.15, 0.25]`; four rejections are decided in `[0.1990, 0.1999]`. The two "real catalog merges"
   are `BLMH × SLC38A6` (GGO, stat **0.1903**) and `PTPN22 × LRR1` (HSA, stat **0.1998** — clears by
   0.0002), and **neither is repeat-bridge**. **Retract "changes 2 real catalog merges" as a stable
   result; quote "0–3/39 across c ∈ [0.19, 0.21], not distinguishable from zero."**
2. **It does essentially nothing on the thesis substrate.** HSA 21/49 = 0.4286, **GGO 1/14 = 0.0714
   [0.0127, 0.3147]**, GGO×HSA 2/11. Delete the framing "the repair works on human" — it is untested
   on gorilla.
3. **It is not separable from a rule that reads no alignment at all.** `len_longer/len_shorter ≥ 3`
   rejects 20/74 vs R2's 24/74; McNemar **b=5/c=1, p = 0.2188**, 91.9% agreement — below the arm's own
   power threshold. Lost-pair length asymmetry: FP 6.180 vs 1.686 kept; TP 7.120 vs 1.285. **The FP and
   TP losses sit on the same axis at the same magnitude.** R2 is a *length-asymmetry guard*, not a
   *repeat-bridge guard*. (R2 does beat the null on the TP arm, 49 vs 252 one-sided losses among
   disagreements — so it is the better rule overall; the FP-side specificity claim is what fails.)
4. **It does not generalise off the flagship.** It kills only 12/30 of the MRPS17/ORAI2 component —
   **the flagship survives**. Repeat-bridge rows excluding the flagship: **5/22 rejected, 1/5
   mechanisms fully killed**; TTC6 0/6, MCFD2 0/5, TDRD5 0/3 untouched. Provenance of the 24
   rejections: 12 flagship + 4 SRGAP1 + 6 singleton one-offs + 1 ZNF480 + 1 SLC38A6 — replicated,
   non-flagship evidence is **6 rejections across 3 mechanisms**. Mechanism unit: fully killed
   **7/24 = 0.2917 [0.1491, 0.4917]**, any rejection 10/24, **median per-mechanism rejection = 0**.
5. **The constraint it satisfies was computed on the wrong universe.** See §5.

**Out-of-sample FP rejection, corrected for clustering (there is no selection optimism to shrink):**
pair-level all species **~0.31, CI [0.14, 0.46]**; **gorilla ~0.07, CI [0.00, 0.25]** — indistinguishable
from zero; for a fresh mechanism, P(any rejection) 0.42 [0.24, 0.61], P(full kill) 0.29 [0.15, 0.49],
**modal outcome 0**.

---

## 5. The recommendation, and what it costs

> ### `recommendation: adopt-as-flag-only`
> ### `repeat_bridge_fixed: false` · `deep_paralog_fixed: false` (out of scope, §3)

**No E_r membership change is adopted.** R2 is the least-bad numeric candidate and it still fails.

### 5.1 The exchange rate on the substrate the catalog actually ships

Restricted to the catalog-node strata (the only substrate `gw_family_catalog` runs on):

| species | FP node rows killed | TP node pairs killed | families destroyed | copies orphaned |
|---|---|---|---|---|
| GGO | **1/14 = 0.0714** | 21/1,795 | **15** | 24/1,070 = 0.0224 |
| HSA | **1/19 = 0.0526** | 10/1,416 | **6** | 14/793 = 0.0177 |

**2 false node edges removed for 21 real families destroyed — 10.5 real families per false edge.**
Named GGO casualties include PPP2CA×PPP2CB, PPP1CA×PPP1CB, TPM2×TPM3, **KIF5A×KIF5C (cov_min 0.952,
id 0.786, rejected at cov_max 0.194)**, SERPINB6×SERPINB9, GMPR×GMPR2, ING4×ING5, KMT2A×KMT2B,
SCN7A×SCN9A, MYL9×MYL12B, BTF3×BTF3L4, NAA10×NAA11, FABP4×PMP2, TMSB4X×LOC129530270, and ZKSCAN3/
ZSCAN21/ZSCAN30 (5 nodes → 3 components). HSA: ACVR1B×TGFBR1, MBNL2×MBNL3, IQCF2×IQCF3, PHKA1×PHKA2,
RPS27×RPS27L, PIP4K2A×PIP4K2B. λ damage inside surviving families: GGO GWFAM34 (GSTM3/GSTM1/GSTM4)
λ 2→1, GWFAM454 λ 2→1; HSA GWFAM230 λ 2→1.

### 5.2 ⚠⚠ There is no constraint-satisfying threshold — HERC2

The sweep's `famintact` iterated only the two census source_keys, so **the ≥95%-of-ceiling constraint
never saw the seven hand-curated panel families.** Include them:

* HERC2 components **1 → 4**; the 114 kb full-length parent `L~chr15_25853560_26064943` goes
  **degree 4 → 0** — **the parent gene is expelled from its own family of partial duplicates**, which
  is exactly the "full-length copy vs partial duplicate" case the repair was supposed to respect.
* HERC2 splits at **any c_long ≥ 0.034**; the parent is orphaned at **c_long > 0.143**.
  **The first FP row dies at c_long 0.05.**

⟹ **No threshold kills a single false edge and leaves HERC2 intact.** Retract
*"R2@0.20 keeps ≥95% of R0's family ceiling in BOTH species."*

### 5.3 It rejects perfect containment — the two canonical great-ape expansions are in the loss set

* **NPIP edges retained 197/260 = 0.7577 [0.7021, 0.8058]** (the family metric hides this — the
  partition survives as one component). **14 of the 63 lost NPIP edges have `cov_min = 1.0000`** — the
  shorter copy 100% contained in the longer — at identity 0.9754–0.9805.
  **Never quote "NPIP edge-Jaccard 1.0000" next to a c_long rule; under R2 it is 0.7577.**
* Lost pairs at `cov_min ≥ 0.95 and identity ≥ 0.95`: **29** (T1 15, panel-7 14) — NBPF10×NBPF19
  (1.0000 / 0.9887), NBPF10×NBPF14, NBPF10×NBPF20, NBPF20×NBPF9, OR4F16×OR4F3, GAGE1×GAGE5, + 14 NPIP.
  **NBPF/Olduvai and NPIP are the two canonical great-ape expansions.** A rule that rejects a perfectly
  contained 98%-identical copy is not a homology rule.
* MAGEA 60/60, TBC1D3 55/55, GSTM 1/1, RABL2 1/1 panel pairs are fully retained — the damage is
  specific, not diffuse.

### 5.4 The short-gene population

Real, and the binding cost — but **smaller at c_long 0.20 than the absolute-floor scenario feared**,
because R2 is a ratio, not a length floor.

* Unit = pair, TP core: shorter ≤ 2 kb **3,931/4,877 = 0.8060 [0.7947, 0.8169]** vs shorter > 2 kb
  **4,369/4,547 = 0.9609** — a **15.5 pp penalty** concentrated on short reps.
* On the shipped node substrate the short-gene damage is modest: GGO ≤2 kb 417/436 = 0.9564,
  HSA 389/399 = 0.9749.
* Unit = family: **133/376 = 35.4% of gorilla TP families (116/274 = 42.3% human) are held together
  ONLY by edges whose shorter rep is ≤ 2 kb** (LDHA/LDHB/LDHC GWFAM272; GFPT1×GFPT2 at node level).
  The feared wipe-out does **not** happen at 0.20 (15 GGO / 6 HSA broken); it starts at **0.25**
  (GGO 336/375).
* ⚠ **R3, the literal per-base floor, is where the fear is fully realised:** at B=1000 short-rep pairs
  retain 1,910/4,877 = 0.3916; at B ≥ 2000, **exactly 0/4,877 survive**, families collapse to
  GGO 109/375. R3's one virtue: panel-7 reps are long, so it retains 392/392 up to B=1100, where
  R2@0.20 has already lost 67.
* ⚠ **The TP loss is not "our truncation".** Lost pairs have median `len_longer/len_shorter` **6.984**
  vs 1.282 retained (n_lost 1,124 / n_kept 8,300) — a 5.4× asymmetry. That is *also* the repeat-bridge
  signature, which is precisely why the two cannot be separated. And 39 lost pairs sit at
  cov_min ≥ 0.90 & id ≥ 0.90 (§5.3), which truncation does not explain.

### 5.5 Threat to P1 and P2 — the structural verdict

**P1 (seed-invariance theorem).** R0–R4 are functions of the two sequences alone, so the theorem
survives them. **R5 breaks P1 outright and must never be a membership condition.** The discriminator is
genuinely strong (human curated tx, ≥50% block overlap: MRPS17[726-1695] hits **50** partners,
SDHAF4[438-643] **73**, ZNF480[875-1319] **105**; TP controls GFPT1 1, GFPT2 1, TRAPPC2 1, LDHA 3,
LDHB 3) — **but the count is a function of the universe**: over the full 7,730-sequence catalog MRPS17's
block scores 50; over 20 random 50% subsets 18–29 (median 25); at 20% median 11; at 10% median 5; at 5%
median 2; **over the 4-node seed {MRPS17, TRAPPC2, MDM2, GREB1} it is 1.** ⟹ under any fixed
promiscuity threshold **the same pair (MRPS17, TRAPPC2) is REJECTED when the catalog is run whole and
ACCEPTED when it is run from a seed.** E_r would stop being a relation between two sequences and become
a function of the node set — the exact negation of the locality P1 rests on. Second defect: MDM2's
densest block hits **291** partners while forming **0** shipped edges ⟹ promiscuity is a property of a
**BLOCK**, not of a gene or an edge.

⚠ **R2 has a weaker version of the same disease.** R0's coverage clause is provably invariant to the
extent of the **longer** node; R2's new clause is **strictly decreasing** in it. P1's own refined
statement is *"membership is seed-invariant; locus **extent** inherits the seed's extent"*
(`seeded_family_definition.md:981`) — **R2 couples membership to the one quantity P1 explicitly excludes
from invariance.** Measured exposure at the observed scale of extent variation (NPIPB8: same 27 members,
median span 9,894 vs 20,314 bp = 2.05×): **GGO 651/2,700 = 0.2411** and **HSA 443/2,574 = 0.1721** of
shipped edges die under a 2× inflation of the longer rep; under the SPEC's q915 extraction (4.618×),
**88.07% GGO / 88.66% HSA** fall. And R2 nearly **doubles substrate discordance**: on the 1,332 related
gene pairs present on both the node-rep and curated-transcript substrates, discordant verdicts go
**R0 298/1,332 = 0.2237 → R2 549/1,332 = 0.4122**, McNemar b=255/c=4, **p = 4.0e-70**; 12/20 of the gene
pairs whose gorilla families R2 destroys are still ACCEPTED by R2 on curated transcripts.
**R5 makes E_r a function of the node set; R2 makes E_r a function of which representative you extract.
Same genre of defect.**

**P2 (γ-quasi-clique refinement, `family_definition.rs` GAMMA=0.20).** R2@0.20 pushes
**GGO 19/494 = 0.0385** and **HSA 8/394 = 0.0203** families below the γ gate, so their emitted partition
stops being canonical and becomes splitter-witness-dependent (the code certifies only "a VALID witness";
exact max-γ-quasi-clique is NP-hard). Flagship margins collapse: **NPIP d 0.688 → 0.521**;
**HERC2 d 0.333 → 0.244, one component → four already at c_long 0.15, and at 0.30 d = 0.178 < γ ⟹ P2
fails outright.**

**λ certificates.** λ decreases for **32/494 GGO, 16/394 HSA**; `cut_certified` flips TRUE→FALSE for
**5 GGO, 3 HSA** (GWFAM64 λ 7→0, GWFAM230 λ 3→0). Every shipped `lambda`/`density`/`cut_certified`
value is a function of E_r and would have to be re-emitted.

**T7 — an edge-only restriction silently changes the NODE SET.** `refine_families` keeps a block only if
`distinct_loci(members) ≥ 2`. At c_long 0.20: **GGO 47/1,415 = 0.0332 copies and 16/494 families vanish
entirely**; HSA 20/1,220 = 0.0164 and 8/394. Of the 24 deleted families, **0 contain an FP edge; 13
contain a TP edge**. The population hit is two-copy families — **16/348 = 0.0460 of GGO two-copy
families cease to exist**, i.e. exactly the population of the excision / false-omission arms (162
two-copy families, 9/162 = 5.6%). **Every node-level headline (915 multi-copy loci, MAPQ-0 0.0004 inside
them, 5.6% false omission, 0.55 reach) is quoted against a node set R2 would change.**

**Retire the sweep's family metric.** "The TP-pair graph stays one component" is not the shipped family.
Scored through the full catalog rather than the arm, the true cost is **GGO 29/494 = 0.0587 splits +
16 families deleted + 47 copies deleted; HSA 14/394 = 0.0355 + 8 + 20** — roughly double the reported
15/6. Monotone: at c_long 0.30, GGO 92 splits / 51 families gone / 156 copies deleted.

One claim survives honestly: R2 is a strict restriction and **empirically adds no new FPs** —
searched for the merge-under-restriction witness, found **0 in both species**. But it is **not a theorem
at the partition level**: the γ gate is non-monotone (a smaller node set can have *higher* density; two
30-node stars joined via u–a, u–v, v–b have d=0.032 < γ and split u from v, but deleting u–a and v–b
leaves the block {u,v} with n ≤ 2, kept whole ⟹ u,v merged). **Prove it or bound it; do not assert it.**

### 5.6 What IS adopted: the flag

Emit the repeat-bridge diagnostic as a **per-edge flag on the emitted family record. Emission is not
definition** — a flag changes no edge, no density, no λ, no node, so P1/P2/λ/T7 are all safe by
construction. Requirements:

1. **Per-edge and per-block**, attributed to the specific aligned block of the specific edge
   (MDM2: 291 partners, 0 edges — a gene-level or edge-level attribution is wrong).
2. **Suppressed by containment.** An edge with `cov_min ≥ 0.95 and identity ≥ 0.90` must **never** be
   flagged (29 such pairs, incl. all 14 NPIP and 4 NBPF, would be flagged at c_long 0.20).
   **Emitting `cov_max` alone reproduces the error in a column.**
3. If the R5-style promiscuity value is emitted, it must carry **its universe (N and a node-set hash)**,
   and **no downstream filter may consume it** — otherwise non-locality re-enters through emission.

### 5.7 The pipeline run that would confirm anything here

Everything above is **T8**. To promote any of it:

1. **For the flag (the recommendation):** run `gw_family_catalog` end-to-end on both catalogs with flag
   emission on, and assert **byte-identical `families.tsv` membership, `density`, `lambda` and
   `cut_certified`** against the current shipped catalog — the flag must be provably inert.
2. **If an E_r change is ever revisited:** score the candidate through `refine_families` +
   `distinct_loci ≥ 2` on the **full** catalog, not the TP arm; re-emit `lambda`/`density`/
   `cut_certified`; and re-derive every node-level headline against the new node set. Add three hard
   structural gates, all of which R2@0.20 fails: **zero families deleted by `distinct_loci ≥ 2`; zero
   families newly below γ = 0.20; zero `cut_certified` TRUE→FALSE flips.** Include the seven
   hand-curated panel families (HERC2 above all) in the family-ceiling constraint.
3. **For the census itself:** re-run `rates.txt` by hand on the corrected reps (§2.4) and re-run the
   three attacks the session limit killed — **30 definitional failures is still a floor, and the
   corrected reps push it up.**

### 5.8 Register entries to add (`NEGATIVE_RESULTS_REGISTER.md`)

* **TRAP** — *coverage-of-the-longer is a length-ratio statistic, not a homology statistic*: it rejects
  perfect containment (cov_min = 1.0000 at id 0.98, 14 NPIP + 4 NBPF pairs), it is not separable from
  `len_longer/len_shorter ≥ 3` on the FP arm (McNemar p = 0.2188), and it flips verdicts between
  substrates (discordance 0.2237 → 0.4122, p = 4e-70).
* **DEAD-END** — *block promiscuity as a membership condition (R5)*: separates cleanly (50/73/105 vs 1/1/1)
  but is a function of the node set (MRPS17 block: 50 whole-catalog, median 25 at 50%, 1 from a 4-node
  seed) ⟹ negates P1. Safe only as a flag.
* **TRAP** — *scoring a rule change on a TP-pair-graph "family" instead of through `refine_families`*:
  under-counts splits ~2× and is blind to whole-family deletion by `distinct_loci ≥ 2` (24 families,
  67 copies).
* **REFUTED** — *"no value of τ or c orders ATP1A1×ATP4A correctly, therefore O1's definition fails"*:
  the two share an mmseqs protein cluster; the truth label is functional (HGNC), not sequence-level.
* **REFUTED** — *"`bedtools getfasta -s` gives transcript-order exons"*: it reverse-complements each
  record individually; every minus-strand multi-exon rep came out with reversed exon order.

---

## 6. ⭐ Is O1 defensible? — the paragraph for the advisor

**Yes, with one honest amendment, and the amendment is smaller and more precise than the census implied.**
The hole is real and it is singular: `E_r`'s coverage clause is normalised on the **shorter** sequence,
which makes it scale-free, so a dispersed repeat that occupies most of a short transcript clears the
coverage floor against essentially anything — an `MRPS17` element bridges unrelated genes at identity
0.75–0.80, reproducibly, in both human and gorilla, with the locus builder removed from the circuit.
That is one mechanism with many instances, not a diffuse defect, and identity has never once been the
culprit. We built a frozen 78-pair false-positive / 20,528-pair true-positive evaluation and swept every
denominator repair we could name. **The finding is that the repair costs more than the defect.** The best
constraint-satisfying candidate — keep the shipped floors, add coverage-of-the-longer ≥ 0.20 — rejects
about a third of the false edges on human curated transcripts, **one of fourteen on gorilla**, is
statistically indistinguishable from a rule that reads no alignment at all and merely compares the two
lengths, and does not even kill its own flagship. Against that it destroys 21 real gene families on the
shipped substrate, deletes 67 copies from the catalog outright, breaks 27 γ-certificates and 8 cut
certificates, rejects perfectly-contained 98%-identical copies in **NPIP and NBPF** — the two canonical
great-ape expansions this thesis exists to describe — and expels the full-length **HERC2** parent from
its own family at any threshold above 0.034, while the first false edge only dies at 0.05: **there is no
threshold that removes one error and leaves HERC2 intact.** What the repair *cannot* fix is worse than
what it can: the true-positive losses and the false-positive kills lie on the *same* length-asymmetry
axis at the *same* magnitude, so a coverage denominator cannot tell "short repeat inside a long gene"
from "short real paralog of a long gene". Separately, the census's strongest-sounding objection —
that no τ or c can order ATP1A1×ATP4A against GFPT1×GFPT2 — **is withdrawn**: those two genes share an
mmseqs protein cluster, and HGNC splits them on *function*. A sequence-based membership relation is not
defective for agreeing with sequence. **So the definition stands, and it stands on its structure, not on
its error rate**: P1 seed-invariance is a theorem precisely because `E_r` is a relation between two
sequences and nothing else, and every candidate repair we tested bought a few percent of false-positive
rejection by making membership depend on something outside that pair — on which representative you
extracted (R2), or on which nodes happen to be in the catalog (R5). We are not willing to trade the one
formal result O1 has for that. **The repeat-bridge class is therefore handled where it belongs — as an
emitted, per-edge, containment-suppressed FLAG on the family record, changing no partition** — and the
honest statement of O1's limit is: *the definition admits a named and now-measurable class of edges in
which a dispersed repeat is a majority of a short representative; we report them rather than silently
excluding them, because every exclusion rule we could construct removed more true copies than false
ones.* Two caveats stated plainly: all sweep numbers are offline PAF re-derivation and none has been run
through `gw_family_catalog` — the flag's inertness must be demonstrated by a byte-identical catalog
re-run — and the census's count of 30 definitional failures remains a **floor**, one that the `mkreps`
correction (§2, 18/111 case-pair verdicts flipping `no → EDGE`) pushes **up**, not down.
