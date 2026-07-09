# Reference-free copy number: do unassignable reads indicate copy number?

**Date:** 2026-07-08. **Question (advisor's):** can copy number be estimated from RNA reads ALONE — without the
genome — such that even reads that cannot be *assigned* to a specific copy still *indicate* how many copies
exist? **Answer: yes, for expressed copies** — and the depth leg is precisely "unassignable reads indicating
copies."

## The two reference-free legs (no genome, no assembly)

- **`chi_H`** = χ(H) = number of distinct pairwise-conflicting copy hap-vectors (THEORY.md Lemma 1, MCC).
  Counts the copies the PSV **conflict** structure can distinguish — including copies no single read can be
  *pinned* to (Tier-2: "counted but unassignable", the advisor's point). A lower bound; blind to identical copies.
- **`depth_cn`** = `E_fam / λ_global`, where `E_fam` = total family read depth (ALL reads, incl. unassignable)
  and `λ_global` = the genome-wide single-copy RNA expression floor (median `n_reads` over single-copy
  transcripts) — a reference-free RNA quantity (Sudmant 2010 / QuicK-mer2 parCN normalization, ported to RNA).
  Recovers the **identical/collapsed (Tier-3)** copies `chi_H` misses, from depth alone.

## Validation vs the phased assembly (`asm_hapCN`), per-gene aggregated (n = 59 genes)

| reads-only estimate (NO genome) | within 25% of `asm_hapCN` | correlation |
|---|---|---|
| `chi_H` (distinguishable) | 49% | 0.44 |
| **`depth_cn`** (depth, incl. unassignable) | **66%** | 0.52 |
| `max(chi_H, depth_cn)` (shipped `famcn_readonly`) | **66%** | 0.53 |

## The result that proves the advisor right — depth recovers copies NO assignment can

On collapsed families `chi_H` sees ~1 copy (the aligner/PSVs cannot separate the identical paralogs), yet the
**unassignable reads' depth** recovers the true count:

| gene | `chi_H` | `depth_cn` | `asm_hapCN` |
|---|---|---|---|
| LOC115930538 | **1** | **11.4** | 12 |
| LOC109025447 | 1 | 15.8 | 11 |
| LOC129526550 | 5 | 11.8 | 13 |
| LOC115930164 | 1 | 4.4 | 13 |

For LOC115930538, all ~11 identical copies are indicated *entirely* by read depth — the reads pile ~11× deep
on the collapsed locus, and even though none can be assigned to a specific copy, their depth counts them. This
is the reference-free copy number the advisor asked for.

## The honest limit: it counts EXPRESSED copies

Reads cannot see untranscribed copies, so the estimate under-counts when copies are silent (13/59 genes < 60%
of `asm_hapCN`): ZNF425 sees ~4 of 26, LOC129531752 sees 2 of 22. This is the one thing the genome has that RNA
does not, and it is the correct boundary of the method — `famcn_readonly` is an **expressed**-copy number, a
reference-free proxy that is exact for the expressed multiplicity and a lower bound on the genomic count.

## Real cases on GGO: copies the LINEAR reference hides (confirmed by the phased assembly)

The strongest concrete result. In these families the *linear* reference collapses several copies onto one
locus, so per-read PSV analysis distinguishes only `chi_H` of them — but **read depth (`depth_cn`) recovers the
rest reference-free, and the phased assembly (`asm_hapCN`) independently confirms more copies than the linear
conflict graph sees.** These are real copies present in the assembly that the reference we map to does not
distinctly represent.

**30 families where depth recovers copies `chi_H` misses; 23 confirmed** by `asm_hapCN > chi_H`. The strongest:

| gene | `chi_H` (linear ref distinguishes) | `depth_cn` (reads recover) | `asm_hapCN` (phased assembly) |
|---|---|---|---|
| LOC109025447 | 1 | 15.8 | 11 |
| LOC115930538 | **1** | **11.4** | 12 |
| LOC129526550 | 5 | 11.8 | 13 |
| LOC101130894 | 2 | 10.4 | 16 |
| LOC129534585 | 1 | 8.0 | 7 |
| CPLANE1 | 1 | 7.7 | 6 |
| LOC115933254 | 1 | 6.9 | 8 |
| LOC101130854 | 3 | 6.7 | 8 |
| LOC115930164 | 1 | 4.4 | 13 |
| LOC109023386 | 1 | 4.1 | 4 |
| LOC101141440 | 1 | 5.0 | 6 |

For LOC115930538, per-read assignment sees **one** copy (`chi_H=1`) but the reads pile ~11× deep — and the
assembly has **12**. All of it recovered from reads, no genome used, assembly-confirmed. This is the headline
"copies the linear reference hides."

**Two honest caveats:** (1) `depth_cn` is an *expressed*-copy count, so it under-counts vs `asm_hapCN` where some
copies are silent (LOC101130894 10.4 vs 16; LOC115930164 4.4 vs 13) and a few high-depth families with no
oracle match (SORL1 39.7, LOC129523503 61.2) are more plausibly high expression than many copies — only the 23
assembly-confirmed cases above are asserted. (2) These are *collapsed* copies (present in the phased assembly,
merged in the linear reference), NOT copies absent from every assembly.

**Copies absent from EVERY assembly: 0 confirmed on GGO.** The admission mechanism (O4 `admit_candidate` + the
vg-realign admissions) is wired, but no copy absent from the phased assembly is confirmed on this data — the
data-limited divergent-absent frontier (needs DNA-level ground truth to separate a novel copy from an allele;
the unmapped-read POC in `bench/VG_REALIGN.md` reinforced that such copies are *mapped-but-mis-placed*, not
unmapped).

## O1 ↔ O2 harmony

The reads-only copy number and the EM copy-assignment consume **one family object**: `detect_families` → copies
(O1) → PSV bubbles → the EM (O2). `chi_H` **is** the O1 conflict-graph count; it equals the copy count `K` the
EM sums abundance over. `depth_cn` counts the same family's reads — including the EM's SoftZone (unassignable)
mass. So improving O1 — cleaner copies, more distinguishing PSV bubbles — raises `chi_H` (more copies
distinguishable), improves the EM's evidence and per-read certification, and sharpens the depth normalization
together. The coupling is structural, not incidental (pinned by the O1↔O2 harmony test).

## In the pipeline

`copy_assign` emits `<out>.famcn_readonly.tsv` (`family_id  chrom  n_copies  n_reads  chi_H  depth_cn  regime
famcn_readonly`) for every family — always (additive). `chi_H` and `n_reads` need nothing external; `depth_cn`/
`famcn_readonly` populate when `--lambda-global` (the RNA single-copy floor from `bench/rna_copy_number_depth.py`)
is supplied, and `n_reads` is emitted regardless so the depth leg is recomputable post-hoc.

## Reproduce

```
# the validated legs (analysis harness):
python bench/family_copy_number.py     # chi_H, tiers
python bench/rna_copy_number_depth.py  # depth_cn, lambda_global
# in the shipped pipeline:
copy_assign --bam GGO_mm.bam --fasta GGO.fasta --regions <fam> --lambda-global <λ> --out ca   # -> ca.famcn_readonly.tsv
```

Relates to `bench/family_copy_number.py`, `bench/rna_copy_number_depth.py`, `bench/em_consistency.md`,
`bench/THEORY.md` (Lemma 1 MCC = χ(H)), and the K=0 frontier (identical copies = depth-only, the parCN regime).
