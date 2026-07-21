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

- **O1 — family definition.** A family is a **γ-quasi-clique of the transcribed-homology graph `E_r`**
  (mutual exon-level homology, ≥2 physically-distinct loci). *Homology is the only membership criterion.*
- **O1 — copy number.** `χ_H` = the **chromatic number of the read-conflict graph** = the Minimum Copy Cover
  (Lemma 1). A read-conflict edge joins two reads that disagree at a shared PSV; a copy is a colour class; a
  *lower bound* (identical copies collapse — the K=0 floor). *The read-conflict graph counts copies; it never
  defines the family.*
- **O2 — copy assignment.** Assign a read to a copy iff the **IsoCon significance certificate** passes
  (`min_p < α/(n−1)`, Bonferroni over the PSV + junction likelihood); otherwise **abstain** (certify Tied).
  **Never split a read 1/k.** One calibrated `α` (1e-3), one per-base error `ε`; no hand-set thresholds.
- **O3 — allele-specific junctions.** On a single molecule, link the **allele** to the splice **junction** it
  co-occurs with (Fisher exact + BH-FDR + transversion/SOR filters). No phasing.
- **O4 — reference-absent copies.** A copy the reads *demand* but the assembly lacks: `χ_H` exceeds the
  annotated loci, **or** depth (`E_fam/λ`) exceeds one copy. Detect-and-flag; copy-vs-allele needs DNA.

## The one theorem and the one test

- **Theorem.** Copy number = `χ_H` = MCC (Lemma 1); assignment = a facility-location / max-weight cover with
  a per-read identifiability certificate. One clean combinatorial object; everything else is engineering.
- **Test (why believe it).** A **planted-genome simulation**: known copy number → exact recovery, and K=0
  identical copies → *certified Tied, not guessed*. Non-circular (truth planted in read names, de-novo
  discovery, multimapper-only scoring). Corroborated on real data: the known families are recovered exactly
  (GSTM 3, PCDHB 5, MAGEA 2, DAZ 2), and the same method tracks *species-specific* copy number on held-out
  human (MAGEA 2→11, TSPY 5→33) — so it is not overfit.

## Consistency status (code vs this spec, 2026-07-21)

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
