# The method, one way

One canonical statement per objective. If the code or a slide says it differently, the code/slide is
wrong — fix it to this. (Status of the code vs this spec is in §"Consistency status" below.)

## The one-paragraph version

A **multi-copy gene family** is a set of ≥2 genomic loci that are *homologous* — a γ-quasi-clique of the
transcribed-homology graph. Within a family, the number of **copies** is `χ_H`, the chromatic number of the
read-conflict graph (a lower bound). Each read is **assigned** to a copy only when a significance test says
it can be — otherwise it abstains. That is the whole method: **homology defines the family; conflict counts
the copies; a significance test assigns the reads.**

## One sentence per objective

⚠ **Renumbered 2026-08-19.** The thesis is **three** objectives. **ASJ (allele-specific junctions) is
DROPPED** — it was the only objective ever rated ATTAINED, so say what was given up. What older docs
call "O3" means three different things; the numbering below is the current one.

- **O1 — family definition.** A family is a **γ-quasi-clique of the transcribed-homology graph `E_r`**
  (mutual exon-level homology, ≥2 physically-distinct loci). *Homology is the only membership criterion.*
  `E_r` = identity ≥ 0.60 **and** coverage ≥ 0.50 of the shorter, on **one single record**, over
  exon-sum spliced representatives. **P1 (seed-invariance) is a theorem.**
- **O1 — copy number.** `χ_H` = the **chromatic number of the read-conflict graph** = the Minimum Copy
  Cover (Lemma 1). A read-conflict edge joins two reads that disagree at a shared PSV; a copy is a colour
  class; a *lower bound* (identical copies collapse — the K=0 floor). *The read-conflict graph counts
  copies; it never defines the family.*
- **O2 — copy assignment under ambiguity.** ⚠ **The target population is ALIGNMENT-SCORE NEAR-TIES, not
  MAPQ-0** — MAPQ-0 is 0.0004 inside the multi-copy loci, while **21.75%** of reads there are within 5% of
  the primary's alignment score. Assign a read to a copy iff the significance certificate passes;
  otherwise **abstain**. **Never split a read 1/k.** ⚠ **Defend O2 on ABSTENTION, not reassignment** —
  held-out abstention gives TPR 0.5066 / FPR 0.0280, AUC 0.7995, where **MAPQ is at chance (0.4944)**.
  **Never claim "assigns better than minimap2"**: net headroom is ~0.1%.
- **O3 — reference-absent / unannotated copies.** A copy the reads *demand* but the assembly lacks.
  **Detect-and-flag with a measured FPR, STRATIFIED by whether the orphaned reads have anywhere to go**;
  copy-vs-allele needs DNA. The stratification is forced by the excision control: a deleted copy is
  **ORPHANED 33.3%** of the time (median 92.7% of its reads unmapped) or **ABSORBED 64.2%** (reads land
  on the best paralogue at **1.75× depth**). *Expression is not the constraint — where the reads go is.*
  | stratum × route | bound |
  |---|---|
  | unique sequence, unmapped-read | **M ≤ 6.4** missing expressed copies |
  | paralogous, unmapped-read | ⚠⚠ **vacuous** (π = 1/35, **formally unbounded**) — and O3's target class lives here |
  | paralogous, depth (S2) | **TPR 0.2703 / FPR 0.0200** held out; set by **divergence not abundance**, ⚠ 45.78% of positives below the 0.01 divergence where it works |

  ⚠ The signature is **UNMAPPED READS, not clipping**. **Never claim completeness.**

## The one theorem and the one test

- **Theorem.** ⚠ **Corrected 2026-08-19.** Copy number = `χ_H` = MCC (Lemma 1). **Assignment is NOT a
  facility-location / max-weight cover** — that framing is retracted: O2 **decomposes**, so **per-read
  argmax is optimal**, and loci must **never** be built with facility location or bipartite matching.
  The load-bearing theorem is **P1, seed-invariance of `E_r`**: membership depends only on the two
  sequences, so the partition cannot be a function of the node set or the seed order.
- **Test (why believe it).** A **planted-genome simulation**: known copy number → exact recovery, and
  K=0 identical copies → *certified Tied, not guessed*. Non-circular (truth planted in read names,
  de-novo discovery, multimapper-only scoring). Corroborated on real data: the known families are
  recovered exactly (GSTM 3, PCDHB 5, MAGEA 2, DAZ 2), and the same method tracks *species-specific*
  copy number on held-out human (MAGEA 2→11, TSPY 5→33) — so it is not overfit.

## What is measured, and the one named hole

> ⚠ **Every quotable figure, with its substrate, lives in [`NUMBERS.md`](NUMBERS.md). Look a number up
> there before quoting it.** The two headline rates below are on **different species** — false-merge is
> HUMAN, false-omission is GORILLA — and must never be pooled.

⚠⚠ **QUOTE THE PROVENANCE WITH THE NUMBER.** The rates below were measured on the **shipped 494-family
catalog**, which was built with `refine` on and which **no invocation of the current binary can
reproduce** (`o1_catalog_provenance.md`). The current default emits **627 families**. Rates needing
truth labels have NOT been re-measured; structural properties HAVE.

**Rates — each on its OWN substrate; read the row, do not assume the catalog:**

| | |
|---|---|
| false-merge rate | **2/150 = 1.33%** [0.37, 4.73], power measured. ⚠ **HUMAN CHM13 v2.0 / A119b, 150 gene-tight single-locus windows — NOT the GGO catalog.** A **specificity and a LOWER bound**, not a precision (no positive stratum ⟹ no prevalence). ⭐**RE-MEASURED 2026-08-20 under the new defaults: unchanged at 2/150, same two windows** — while spurious E_r edges on the same panel fell **28 → 3** (`o1_false_merge_remeasured.md`) |
| false-omission rate | **9/162 = 5.6%** [0.0295, 0.1022] |
| identity-clause failures | **0/728** — the failure mode is localised to the coverage clause |
| DNA vs RNA partition, same loci | **identical, 7/7** |
| reach | ~0.55 of families genome-wide (chr1 22/40, representative at Fisher p = 0.6090) |

**Structure — RE-MEASURED 2026-08-20 on the 627-family catalog, from the pipeline's own certificates:**

| | 494 (offline) | **627 (pipeline)** |
|---|---:|---:|
| 2-copy share — no split possible | 0.7045 | **0.7018** |
| n ≥ 3, the hierarchy ceiling | 0.2955 | **0.2982** |
| complete graphs — γ provably inert | 0.1012 | **0.0893** |
| real reach (density < 1) | 0.1923 | **0.2089** |
| λ = 1 with n ≥ 3 | 0.1599 | **0.1786** |

⭐ **The offline re-derivation was faithful.** Every quantity agrees to within **±0.019** across a
different catalog (+27% families), a different transcript set (94,257 vs 79,569 skeletons), a different
code path (refine on vs off) and a different method (offline reconstruction vs the shipped Rust). The
structural claims survive their substrate change — including the one that matters most:

> **γ is provably inert on 0.7018 + 0.0893 = 79.11% of the catalog** (two-copy families plus complete
> graphs), so the graph-theoretic content of the definition applies to ~21% of its output.

⭐ **Only 75/627 = 11.96% of families are `cut_certified` (λ ≥ 2).** Among n ≥ 3 families,
**112/187 = 59.9% have λ = 1** — a single alignment record holds them together. (For 2-copy families
λ = 1 is arithmetic, not a defect.)

**The hole.** ~30 of 105 classified bad-family cases are **definitional**, and they are **one
mechanism**: the min-length coverage denominator is **scale-free**, so a ~1 kb dispersed repeat is
≥ 0.50 of any node under 2 kb, and 24.88% of gorilla copies are ≤ 2 kb. Exposure ceiling
**41/494 = 8.30%** of gorilla families, and **30 is a floor**.

**Five repair routes are closed**, all recorded in `NEGATIVE_RESULTS_REGISTER.md`: thresholds
(impossibility — HERC2 splits at c_long ≥ 0.034 before the first FP dies at 0.05), coverage-of-longer
(dominated), rep repair (no discriminator), read tiling (AUC 0.1259, significant in the **wrong**
direction), and rare-anchor replacement (no operating point — see below).

**Two guards exist, both pair-local by construction so P1 is untouched. Neither has run through the
shipped binary (T8), and the definition is UNCHANGED.**

1. **Transcript-orientation guard** — 29/74 FPs rejected for 4 lost edges of 9,032.
   `o1_false_positive_rules.md`. Blocked on a whole-genome GGO/HSA comparison.
2. **Genome-anchored repeat veto** — `min_shared_gmult`, the minimum genome-wide occurrence count of a
   shared canonical 21-mer. 10/12 scored FPs at **0/135** TP cost; seed-invariance demonstrated
   (**0/147** vs the catalog-counted analogue's **94/147**). `o1_genome_anchored_repeat_gate.md`.
   ⚠ It is a **veto, never an admission criterion** — replacing coverage with it has no operating point.
   ⚠ Ship as a **flag** first: as a gate it changes `E_r`, so every `density`/`λ`/`cut_certified` value
   must be re-emitted and γ margins move (the analogous R2 pushed 19/494 GGO families below γ).

## Consistency status (code vs this spec, 2026-07-21) — HISTORICAL

⚠ This table is a **2026-07-21 snapshot** and is not maintained. Its "O3 published number"
row refers to the **dropped ASJ objective**. Read it as a record of that pass, not as current state.

Honest map of where the shipped code already matches the spec and where two approaches still coexist.
Deliverable B attempted to consolidate each coexisting pair under a byte-identity gate
(`bench/mechanism/byte_identity_gate.sh`): **the items that consolidated are byte-identical (gate-verified);
the three that were attempted (EM, ≥2-loci, refiner) each DIVERGE — they were never equivalent — so they are
kept-both and deferred to deliverable C, where the correct choice gets a before/after. Deliverable B changed
no number.** Divergence detail and reproduction steps: `bench/mechanism/consolidation_divergences.md`.

| Item | Canonical (this spec) | Still coexisting → action |
|---|---|---|
| Family definition | `gw_family_catalog` + `copy_assign`: E_r homology gate (`refine_families_exon_sum`, default ON) | ✅ **relabeled (B2.1, doc-only, byte-identical, gate passed):** `family_define` binary is now documented as a **legacy parity fixture** reproducing a frozen Python catalog from precomputed TSVs; `gw_family_catalog` is the sole live catalog. |
| "family = ?" wording | homology defines; conflict counts | ✅ fixed: `read_conflict.rs` + `gw_family_catalog.rs` docstrings now say this (were calling conflict "the definition"). |
| γ-quasi-clique refiner | one refiner | ❌ **DIVERGES (B2.4) → NOT consolidated, deferred to deliverable C** (see `consolidation_divergences.md`). Still **two**: `family_definition::refine_component` (CNM, parity-tested) vs `family_split::gamma_quasi_clique_partition` (Louvain, live path). Fixture-proven to partition differently (star-graph case: CNM merges hub+6 leaves into one block, Louvain shatters into 15 singletons); the gate corpus never exercises this divergently (its only reachable call site short-circuits on a K3 triangle, density ≥ γ, before either splitter runs) — a vacuous PASS was confirmed and reverted rather than landed. |
| "≥2 distinct loci" | one predicate | ❌ **DIVERGES (B2.3) → NOT consolidated, deferred to deliverable C** (see `consolidation_divergences.md`). Still **two**: `distinct_locus_reps` (any-overlap, live) vs `distinct_loci` (reciprocal-50%). Fixture-proven to collapse loci differently (nested-fragment case: any-overlap collapses to 2 loci, reciprocal-50% keeps 3, a spurious extra locus). The gate's raw PASS on this item is vacuous — the corpus reaches the collapse function exactly once, on 5 mutually non-overlapping copies, so it is a no-op under both predicates by construction. |
| Copy number name | `χ_H` = chromatic number = MCC | ✅ fixed (was "min path cover" / "facility-location count" in the glossary/figures). |
| Abundance (soft leg) | one EM (`em_assign_family`, convergent, junction-aware) | ❌ **DIVERGES (B2.2) → NOT consolidated, deferred to deliverable C** (see `consolidation_divergences.md`). Still **three** paths (`soft_quantify_em` default — QUANT_ERROR=0.01/100-iter, `em_assign_family` under `--vg-realign` — error_rate=0.003/eps=1e-6/200-iter, a third file under `--em`); the first two give different quant on the gate corpus itself (gstm.quant.tsv abundance 0.0564 vs 0.0565 — different error model, entangled with the 5-epsilon inconsistency). |
| O3 published number | shipped `asj`/`asj_verify` binaries | the 54-call genetic core (SOR filter) lives in **unwired** modules; the binaries emit a *different* TSV — **wire the genetic core into the binary** so the tool reproduces the published number. |
| Assembler (41 modules) | not in the thesis build | unreachable from the thesis binaries but still compiled; pinned by 12 tendril symbols — carve per `RETIREMENT_AND_MIGRATION.md`. |
| Hardcoded scratch paths (B1a) | explicit `--meta/--annot/--skeletons/--genome` | ✅ **removed (byte-identical, gate passed):** `family_define` no longer silently defaults to `/home/juanfra/winloci_scratch` — those flags are now required. |
| copy_assign caps (B1b) | opt-in `--poa-cap`/`--read-cap`, defaults unchanged | ✅ **added (byte-identical, gate passed).** FINDING: neither is on copy_assign's main path today — `--read-cap` is a no-op warn (its `o2_materialize` target is a python-parity port imported by no binary); `--poa-cap` only bites under the opt-in `RUSTLE_INTRON_PSV=1` path (wired via `RUSTLE_POA_CAP`). |

**Gate-corpus coverage gap:** the byte-identity gate corpus (`copy_assign` on the GSTM region + `gw_family_catalog`
on 2 regions) does **not** exercise the code paths where B2.3/B2.4 diverge — GSTM's homology graph is a dense
triangle that short-circuits the refiner before either implementation runs, and no overlapping same-chrom loci
in the corpus ever reach the ≥2-loci predicate. Both divergences were proven with targeted fixtures, not the
gate. Deliverable C must extend the gate corpus (add an overlapping-loci case and a low-density/star-graph case)
before either consolidation can be decided and gate-verified for real.

**Nothing is safe to blindly delete right now** — every remaining "dead" module is compile-pinned or a tested
deliverable. The remaining simplification is *consolidation* (pick one of each pair above, in deliverable C)
and the assembler carve, both refactors that must keep the numbers identical.
