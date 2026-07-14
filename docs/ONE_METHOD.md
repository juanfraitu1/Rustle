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

## Consistency status (code vs this spec, 2026-07-13)

Honest map of where the shipped code already matches the spec and where two approaches still coexist and are
being **consolidated to one** (none of these change the numbers — see the verified regression checks).

| Item | Canonical (this spec) | Still coexisting → action |
|---|---|---|
| Family definition | `gw_family_catalog` + `copy_assign`: E_r homology gate (`refine_families_exon_sum`, default ON) | `family_define` binary reproduces a **frozen Python catalog** from precomputed TSVs — demote to `--legacy` parity fixture; `gw_family_catalog` is the sole catalog. |
| "family = ?" wording | homology defines; conflict counts | ✅ fixed: `read_conflict.rs` + `gw_family_catalog.rs` docstrings now say this (were calling conflict "the definition"). |
| γ-quasi-clique refiner | one refiner | **two** (`family_definition::refine_families` CNM, parity-tested, vs `family_split::gamma_quasi_clique_partition` Louvain, on the live path) — consolidate onto the parity-tested one. |
| "≥2 distinct loci" | one predicate | **two** (`distinct_loci` reciprocal-50% vs `distinct_locus_reps` any-overlap) — pick the reciprocal-50% rule. |
| Copy number name | `χ_H` = chromatic number = MCC | ✅ fixed (was "min path cover" / "facility-location count" in the glossary/figures). |
| Abundance (soft leg) | one EM (`em_assign_family`, convergent, junction-aware) | **three** paths (`soft_quantify_em` default, `em_assign_family` under `--vg-realign`, a third file under `--em`) — consolidate onto `em_assign_family`. |
| O3 published number | shipped `asj`/`asj_verify` binaries | the 54-call genetic core (SOR filter) lives in **unwired** modules; the binaries emit a *different* TSV — **wire the genetic core into the binary** so the tool reproduces the published number. |
| Assembler (41 modules) | not in the thesis build | unreachable from the thesis binaries but still compiled; pinned by 12 tendril symbols — carve per `RETIREMENT_AND_MIGRATION.md`. |

**Nothing is safe to blindly delete right now** — every remaining "dead" module is compile-pinned or a tested
deliverable. The remaining simplification is *consolidation* (pick one of each pair above) and the assembler
carve, both refactors that keep the numbers identical.
