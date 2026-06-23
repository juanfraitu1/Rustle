# Genome-wide PSV resolution: where great-ape gene-family copies are read-resolvable

The PSV-aware variation graph (`psv_graph_demo.py`, the two-family demo) scaled to **every**
validated multi-copy gene family. This is the empirical, genome-wide instantiation of the
identifiability arc in `family_to_copy_bridge.md`: detection hands each family a backbone, the
bubbles on that backbone are the columns, and threading the reads through the bubbles either
assigns them to a copy or hits the **K-frontier**.

## Method

Input: the 196 backbone-reinforced families from `family_def_vg_reinforce.py` (de-novo loci,
`~R ∩ ~B`). For each family (`psv_graph_genomewide.py`):

1. **Copy-regions.** Collapse de-novo loci that are the same genomic copy (reciprocal overlap
   > 0.5 = isoforms / read-through nests) into one region. Families whose members all collapse
   to one region are **not genomically multi-copy** and are skipped (42 of 196 — these are
   single loci with isoform-multiplicity, e.g. a 335 kb transcript nested over a 178 kb one).
2. **Backbone = longest copy.** Align the other copies (`minimap2 asm20`) and the reads
   (`minimap2 splice:hq`) onto it.
3. **PSV bubbles** = backbone columns where the copies carry ≥2 alleles **and** the column is
   read-supported (≥3 reads — the recurrence test that separates a paralog difference from a
   one-read HiFi error).
4. **Copy paths** = each copy's allele string; **K** = number of distinct paths (resolvable
   haplotypes). **Thread** each read by best PSV-allele match (HiFi-tolerant).

## Verification (`psv_graph_verify.py`)

Two corrections and one load-bearing check, before any number is reported:

- **Cross-family dedup.** Nine families were the *same* genomic copies re-detected as several
  de-novo isoform-loci (e.g. 168/169/170 = one 2-copy locus counted three times). Collapsed by
  region-set match: **154 → 145 unique families**, and read totals recomputed over uniques so
  no locus is counted twice.
- **`no_psv` triage.** A family with 0 PSVs is only a *real* collapse if its copies actually
  align to the backbone. Re-aligning each collapsed family's copies and grading aligned-fraction
  + identity splits the 24 zero-PSV families into **18 genuine K=0** (every copy aligns ≥0.9 of
  its length at ≥0.99 identity, yet no read-supported divergent column exists → truly exonically
  identical) and **6 align_fail** (a copy does not map to the backbone — a genomic-span method
  limit, *not* a biological collapse).

## Result (verified, 145 unique families)

| class | families | % | meaning |
|---|---|---|---|
| **fully resolvable** (K = #copies) | 118 | 81.4 | every copy read-distinguishable |
| **partial** (2 ≤ K < #copies) | 3 | 2.1 | some copies separate, some collapse |
| **genuine K = 0** | 18 | 12.4 | copies exonically identical → provably unresolvable from RNA |
| **indeterminate** (align_fail) | 6 | 4.1 | genomic-span alignment failed; needs exon-union models |
| **read-resolvable total** | **121** | **83.4** | |

**Reads:** 64,066 threaded across the unique families — 34,020 (53 %) assigned to **one copy**,
660 to a collapsed group, 22,095 (34.5 %) cover **no PSV** (the *coverage* face of the
K-frontier — a full-length read that happens to span only shared exons), the rest
unexplained/ambiguous. Single-copy assignments **agree with the independent best-mapping copy
95.3 %** of the time.

## Why this matters

- **Independent convergence.** 12.4 % genuine K=0 from RNA read-threading lands almost exactly on
  the **12 % K=0** from the combinatorial copy-assignment census (`copy_assignment_theory.md`,
  reached with no read-threading at all). Two unrelated methods agree on the size of the
  unresolvable core.
- **The dichotomy is structural, not anecdotal.** 81 % of families lie on the `K = #copies`
  diagonal (Fig. B): the bubbles separate every path. The K=1 row is the frontier. The demo's
  RABL2 (resolvable) and RFPL4A (collapsed) were not cherry-picked — they are the two ends of a
  genome-wide distribution.
- **Detection ⊇ resolution, measured.** Every one of these 145 families was *detected* (it is a
  `~R ∩ ~B` component); 83 % are also *resolvable*. The 12 % genuine-K0 gap is exactly the
  bridge note's prediction: detectable as a family, provably unresolvable into copies.

## Honest caveats

- **Genomic-span, not exon-union.** The scan aligns copies by their genomic span (introns
  included), so copies with fast-diverging introns can fail to align (the 6 indeterminate). The
  family/PSV definition is exonic; an exon-union re-run would reclassify these. They are reported
  as a separate bucket, not folded into either side.
- **Validation is a proxy, not ground truth.** "Agreement with the best-mapping copy" treats
  minimap2's primary as truth, which is exactly what is unreliable for paralogs. 95.3 % global
  agreement is reassuring; the ~5 % disagreement concentrates in size-heterogeneous families
  (family 32: a 124 kb backbone bundled with ~30 kb paralogs, where mapping itself is arbitrary)
  — there, PSV disagreement may be a *correction*, not an error.
- **De-novo loci are imported** from detection; the copy-region dedup (overlap > 0.5) is the one
  free parameter added here.

Artifacts: `psv_graph_genomewide.py` · `psv_graph_verify.py` · `psv_graph_genomewide.png` ·
`psv_graph_genomewide_perfamily.tsv` · `psv_graph_{genomewide,verify}.json`.
