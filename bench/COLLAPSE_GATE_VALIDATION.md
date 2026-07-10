# Collapse gate: built, tested, and refuted by a control

**Date:** 2026-07-09. Implements `docs/superpowers/specs/2026-07-09-collapse-gate-design.md` (Tasks 1–4).
**Shipped OFF by default.** Tasks 5–6 are cancelled — the reason is below, and it is not a tuning problem.

## What was built

`src/rustle/vg_family/collapse_gate.rs`: an ambiguity statistic (`k` MAPQ-0 primaries out of `n`), a binomial
upper-tail p-value against a genome-wide background, and a two-leg verdict that consults haplotypes **only** if
the collapse leg fires. Wired into `detect_and_assign` for reps that no family claimed. `--collapse-gate` enables
it; `--eps-amb` overrides the background. 866 tests.

## It works on DAZ

```
[detect_and_assign] collapse gate: NC_073248.2:42783128-42859657 ambiguous 1096/1274 p=0.00e0
                    -> chi(H)=2 copies (reads certified TIED; chi(H) is a LOWER BOUND)
[copy_assign] 1 families, 1274 read assignments   (status: tied=1274)
```

**χ(H) = 2 — exactly the annotated DAZ copy number**, from a window that previously produced zero families.

## And it fires on a single-copy control

```
[detect_and_assign] collapse gate: NC_073229.2:97606864-97612687 ambiguous 67/3588 p=1.95e-52
                    -> chi(H)=7 copies
```

**EEF1A1.** One copy at that locus. `χ(H) = 7`.

Tracing its MAPQ-0 reads: they also align to **NC_073224.2** and **NC_073227.2** — EEF1A1's processed pseudogenes,
on other chromosomes. The reads are ambiguous because a paralog exists *elsewhere*, not because this locus is
collapsed.

## Why this is a design error, not a threshold error

**MAPQ 0 means "this read maps equally well somewhere else". It does not mean "this locus is collapsed."** The two
coincide only when the somewhere-else is the same locus.

The logic actually inverts. If a copy is genuinely **absent** from the reference — the case a collapse gate exists
to catch — its reads pile onto the copy that *is* present, at **high** mapping quality. A true reference collapse
therefore produces **depth excess and no ambiguity at all**.

That is exactly why SDA detects collapses by read depth:

> "We defined collapsed regions as those with a mean sequence coverage > 3 s.d. beyond the mean sequence coverage
> of the de novo assembly…" — Vollger et al. 2019, Methods

The spec argued: *in DNA a collapse shows as excess depth; in RNA it shows as reads the aligner cannot place
uniquely.* The first clause is right. **The second is wrong**, and the EEF1A1 control is the proof. What ambiguity
measures is **unresolvable paralogy** — which is precisely what `read_conflict`'s `E_c` oracle already measures.
The gate reinvented E_c and mislabelled it.

DAZ does not rescue the design. DAZ2 **is** in the reference, 20 kb from DAZ1; its reads are ambiguous between two
*present* loci. The gate got χ(H) = 2 for a locus whose real defect is that **DAZ2 never assembled** (20 primary
reads). Right answer, wrong reason.

## Consequence for χ(H)

χ(H) was to be reported as a **lower bound** on copy number. At EEF1A1 it returns **7** for a one-copy locus, so it
is not a lower bound on copies — it is a count of read **haplotypes**, inflated by het alleles, RNA editing and
isoform noise. This is the same copy-vs-allele wall as `bench/COLLAPSED_COPY_GATE.md`, reached from the other side.

## Status

| | |
|---|---|
| Tasks 1–4 | **done**, 866 tests, gate **OFF by default** |
| Task 5 (planted 2 copies × 2 alleles vs 4 copies) | **cancelled** — the instrument fails before the sim can adjudicate |
| Task 6 (`--admit-collapsed-copies`) | **cancelled** — materialising copies from these haplotypes would be indefensible |

The code is kept, off, because the statistic and its tests are correct for the question it *actually* answers
(is there an unresolvable paralog?). It must not be used for copy number.

## What would actually recover DAZ2

Not a gate. DAZ2's defect is **assembly**: 20 primary reads, 19 at MAPQ 0, no rep. Two candidate routes, both
recorded elsewhere:

1. **Depth excess** (SDA's real instrument, Bailey's WSSD) — needs an RNA expression model, because RNA depth is
   copy number × expression. Recorded as DNA-bound in `bench/COLLAPSED_COPY_GATE.md`.
2. **Reference-free candidate generation** (IsoCon step 1: nearest-neighbour read clustering in edit-distance
   space) — recovers a copy the aligner refused to place. This was the original item 3 in
   `bench/PRIOR_ART_LEVERS.md`, and it is still the right one.

Related: `docs/superpowers/specs/2026-07-09-collapse-gate-design.md`, `bench/COLLAPSED_COPY_GATE.md`,
`bench/YAG_CHECK.md`, `reference_sda_vollger`.
