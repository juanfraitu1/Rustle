# Non-Exon-Signal Rescue POC — Design

**Date:** 2026-07-21
**Status:** approved design, pre-implementation
**Type:** analysis / proof-of-concept (NOT a pipeline change).

**Motivation.** The largest genuine-DNA-only floor group is "exon-homogenized K=0" — copies
whose *exons* are homogenized by gene conversion so RNA reads tie in exon space. On PPIAL4 we
showed that when the family's full-length reads are aligned to the **full-genomic** copy
references (UTR + intron + flank), **32% become uniquely assignable** (290/302 via soft-clip
flank) — signal the exon-space PSV search never sees. This POC quantifies that rescue across the
families that actually need it, to give a concrete number for the advisor meeting.

**Scope (per user).** Only the families **below 100% member sensitivity** — no point testing the
42 already-perfect ones. Of the 40 incomplete families, **24 carry exon-homogenized / young-
identical (K=0) members** — those are the rescue targets. The POC runs on those 24.

---

## 1. What it measures

For each of the 24 K=0-bearing families: how many of its reads (and its K=0 members) become
**distinguishable** when the read is scored against the **full genomic sequence of every copy in
the family** — i.e., using UTR/intron/flank, not just the homogenized exon.

**Distinguishable** = the read aligns uniquely to one copy (minimap2 primary MAPQ > 0 against a
multi-copy reference of all the family's copies). This is the honest, all-copies test — stricter
than the 2-copy PPIAL4 probe.

---

## 2. Method (per family)

1. Get the family's copies from `bench/soto/80_fams.chr.bed` (all members of that `ID_*`).
2. For each copy, extract its **genomic sequence + flank** from `chm13v2.0.fa`
   (`samtools faidx chrom:start-FLANK..end+FLANK`, FLANK = 1500 bp) → a multi-copy reference FASTA
   with each copy named `copy<idx>`.
3. Extract the family's reads (all copy loci, deduped by read name) from `A119b.t2t.bam` → fastq.
4. Align reads to the multi-copy reference: `minimap2 -ax map-hifi --eqx`.
5. Count primary alignments by MAPQ: **distinguishable (MAPQ > 0)** vs **tie (MAPQ 0)**; record
   which copy each distinguishable read hit, and whether it carries soft-clip flank.
6. A K=0 member is **rescued** if ≥ `min_reads` (3, the existing floor) distinguishable reads hit
   its copy.

**Output:** `bench/soto/nonexon_rescue.tsv` — per family: `family_id, gene, n_copies, n_reads,
n_distinguishable, pct_distinguishable, n_softclip_flank, k0_members, k0_members_rescued`. Plus a
one-paragraph aggregate in `bench/soto/nonexon_rescue.md`: *"M of the K=0 members across the 24
families become distinguishable with non-exon signal; X% of family reads are rescued on average."*

---

## 3. Honesty rails (must be in the output)

- **This shows the information exists, not that the pipeline uses it.** MAPQ > 0 from genomic
  alignment ≠ a copy-assignment; wiring non-exon columns into `copy_assign`'s PSV search is the
  separate production follow-up.
- **All-copies test** (not 2-copy) so "distinguishable" means resolvable among all siblings.
- **A rescued member must map to a real annotated Soto member** (it does — copies come from the
  Soto BED); report per-member, no fabricated rescues.
- Families where copies are too few (<2) or reads too sparse are reported as `n/a`, not forced.

---

## 4. Success criteria

- `nonexon_rescue.tsv` covers the 24 K0-bearing families (or reports why any is skipped).
- The aggregate rescue number (K=0 members rescued / total K=0 in these families) is computed
  from real alignments, with the soft-clip-flank fraction (the mechanism).
- The output states the honesty rails (info-exists ≠ assigned; all-copies; follow-up is the
  production wiring).

---

## 5. Out of scope

- Any change to `copy_assign_pipeline.rs` / the PSV search / the significance certificate — that
  is the production feature, a separate spec.
- The genome-wide seeding-fix member-sensitivity number (the interrupted catalog run) — can be
  re-run later; this POC is about the non-exon rescue signal specifically.
- The 16 incomplete families whose misses are seeding-bug (not K=0) — the seeding fix addresses
  those; they are not non-exon-rescue targets.
