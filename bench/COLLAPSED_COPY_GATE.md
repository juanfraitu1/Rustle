# Can collapsed copies be wired into the family gate? Not as stated.

**Date:** 2026-07-09. Follow-up to `YAG_CHECK.md` / `PRIOR_ART_LEVERS.md`. **No code changed.**

## The idea

DAZ2 and BPY2 exist in the data but produce no family. After the readthrough filter, the DAZ window leaves a
single assembled rep (DAZ1), and `colocated_families` requires `>= min_copies` **reps**. Yet DAZ2's reads are
present — 20 primaries (19 at MAPQ 0) and 1119 secondary placements over its span — piled onto DAZ1 because the
two copies are near-identical.

`recover_collapsed_copies` already reports **3 PSV-distinguishable collapsed copies** at that rep (it emits every
identifiable haplotype except the host, so 3 emitted = 4 identifiable). Its return value is **`eprintln!`'d and
discarded**. The proposal was therefore to count copies rather than reps:

> a locus whose reads carry `>= min_copies` PSV-distinguishable haplotypes is a multi-copy family, χ(H) ≥ 2

This is the advisor's own "even unassignable reads indicate copies", made operational. It needs no new machinery.

## The measurement that refutes it

Gorilla is diploid. Two het alleles at an ordinary single-copy gene are two haplotypes, and the detector cannot
tell them from two copies. Measured on single-copy control genes, with `--recover-copies`:

| locus | collapsed candidates | AS-tied secondaries |
|---|---|---|
| **DAZ** — true 2-copy, collapsed | **3** | **977** |
| **TSPYL1** — single-copy intronless gene | **12** | 4 |
| **DERPC** — single-copy gene | **5** | 1 |
| GSPT2 | 0 | 78 |
| ATXN7L3B | 0 | 0 |

**TSPYL1, a single-copy gene, reports four times more "collapsed copies" than DAZ.** The haplotype count alone is
not merely noisy — it is *anti-correlated* with truth on this sample. Wiring it into the family gate would admit
het alleles, RNA-editing sites and isoform noise as copies genome-wide. This is the same
`copy-vs-allele` wall recorded in `project_reference_absent_wiring` ("MECHANISM-ONLY: 0 real copies").

## What does separate them

The **collapse evidence**, not the haplotype count. DAZ has **977 AS-tied secondary placements**; the single-copy
controls have 0–78. A genuinely collapsed copy is one whose reads map equally well somewhere else — that is what
makes it collapsed. The haplotype split is only meaningful *conditional on* that.

This is exactly SDA's design (`reference_sda_vollger`): step 1 detects collapsed segdups by **read-depth excess**,
and only then are PSVs defined and correlation-clustered. We skipped step 1. Bailey's WSSD (depth → copy number,
R² = 0.96) is the same instrument, and we already compute its RNA analogue as `depth_cn` in
`readonly_copy_number.rs`.

## Revised proposal (not implemented)

Admit a single-rep locus as a multi-copy family only when **both** hold:

1. **Collapse evidence** — the locus shows depth/multimapping excess (`depth_cn >= 2`, or a dominant fraction of
   its reads carrying AS-tied secondary placements). This is the SDA/WSSD step we never ran.
2. **Haplotype evidence** — its reads then split into `>= min_copies` PSV-distinguishable, identifiable
   haplotypes.

On the five loci above this separates cleanly (DAZ passes both; TSPYL1 and DERPC fail (1); GSPT2 fails (2)). That
is **n = 5**, which is not a validation — it is a hypothesis worth a spec.

## Open risks

- Depth excess in RNA is confounded by **expression**, not just copy number. `depth_cn` normalises by a
  single-copy anchor (λ_global), but a highly expressed single-copy gene still has high depth. The correct
  statistic is depth *relative to comparable loci*, and that needs care — this is the weakest link.
- AS-tied secondary count is a proxy, not a quantity with a null distribution. `≥ 2 reads`-style floors will not
  do; the gate should be a test, not a threshold.
- Two het alleles of a *collapsed* copy pair give 4 haplotypes, which is what DAZ shows. Distinguishing 2 copies ×
  2 alleles from 4 copies is the same problem, one level up, and DNA parCN is the only clean oracle.

## Verdict

The cheap wiring is wrong. The guarded version is plausible and paper-backed, but it turns "wire up an existing
count" into "run the collapse-detection step we skipped", which deserves a design pass rather than an inline edit.

Related: `bench/YAG_CHECK.md`, `bench/PRIOR_ART_LEVERS.md`, `reference_sda_vollger`, `reference_bailey_2002_segdups`.
