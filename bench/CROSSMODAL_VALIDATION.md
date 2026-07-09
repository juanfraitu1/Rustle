# Cross-modal validation of the RNA E_r family catalog

**Date:** 2026-07-08  **Substrate:** gorilla (GGO) HiFi Iso-Seq → `GGO_mm.bam`; genome `GGO.fasta`.
**Catalog:** `gw_family_catalog --homology-primary --min-identity 0.98 --enumerate-copies` (genome-wide) = **211
families / 520 copies** (39 cross-chromosome). Companion to the homology-primary design spec.

Three orthogonal, published witnesses corroborate the RNA de-novo family catalog, and each shows where the
RNA/long-read method uniquely contributes.

---

## Axis A — vs Liftoff `-copies` (DNA + annotation), head-to-head

**Method** (`bench/sim_liftoff_headtohead.py`): plant 7 copies of ONE gene at 0–30 % divergence, give BOTH
tools the same annotated gene model, threshold 0.85.

| per-copy divergence | copy identity | Liftoff `-sc 0.85` | exon-sum + minimap2 (`-x splice -p 0.01`) |
|---|---|---|---|
| 0 % | ~1.00 | ✅ | ✅ |
| 5 % | ~0.96 | ✅ | ✅ |
| 10 % | ~0.90 | ❌ | ✅ |
| 15 % | ~0.85 | ❌ | ✅ |
| 20–30 % | ≤0.80 | ❌ | ❌ |
| **total found** | | **2 / 7** | **4 / 7** |

**Result:** exon-sum + minimap2 finds **more** copies and reaches further into divergence (down to ~0.84
identity), recovering the 90 % and 85 % copies Liftoff missed *even though they clear Liftoff's own 0.85
threshold*. Liftoff's default `-sc 1.0` reports only 100 %-identical copies; on real data `-sc 0.90` found **0**
copies of the ~19-copy GSTM cluster (its paralogs are annotated separately and sit below 0.90). `asm20`/`asm10`
presets found 0 — confirming `-x splice` is the correct preset for a spliced query vs genomic copies.

**Conclusion:** the exon-sum method matches the proven Liftoff on recent copies and **extends past it** into the
divergent regime — where Liftoff, by design, cannot go. (Requisite: `-p 0.01`; the default `-x splice -p`
suppresses divergent secondaries — a bug we fixed in `project_family_copies`.)

## Axis B — vs SEDEF SD98 (DNA segmental duplications), the Soto regime

**Method:** SEDEF `final.bed` filtered to gorilla-vs-gorilla pairs ≥98 % identity = **5,774 SD98 segdups**;
a RNA-98 family is DNA-confirmed if ≥1 copy overlaps a SD98 interval (`bench/crossmodal_validate.py`).

**Result:** **131 / 211 RNA-98 families (62.1 %) are DNA-confirmed** by SEDEF SD98. The **80 RNA-only** families
are either more divergent than the DNA ≥98 % regime, or expressed copies SEDEF's ≥1 kb/≥98 % criterion did not
pair. This is a circularity-free RNA-vs-DNA cross-check (independent modalities, same 98 % threshold).

## Axis C — vs the phased assembly (Soto famCN / parCN)

**Method:** project each family's consensus onto the genome (one minimap2 pass, `-p 0.01`), and count genomic
copy loci in two identity/coverage bands, then join to the phased-assembly oracle
(`bench/diploid_cn_oracle.tsv`, `asm_hapCN` = assembly copy number):
- **famCN@0.98** (identity ≥ 0.98, coverage ≥ 0.90) = the recent-duplication count (Soto SD98 famCN).
- **totalCN@0.80** (identity ≥ 0.80, coverage ≥ 0.50) = the full copy number including divergent paralogs.

**Result** (51 oracle-matched families):

| metric | within 25 % of assembly `asm_hapCN` |
|---|---|
| famCN@0.98 | **37 / 51 (72 %)** |
| totalCN@0.80 | **36 / 51 (70 %)** |

25 families had divergent genomic copies recovered (`totalCN > n_rna_expressed`). Examples:

| family | expressed (RNA) | famCN@0.98 | totalCN@0.80 | assembly `asm_hapCN` | gene |
|---|---|---|---|---|---|
| GWFAM21 | 3 | **22** | 22 | **22** | LOC129531752 (perfect) |
| GWFAM9  | 8 | 13 | **21** | **19** | GSTM2 (totalCN recovers the divergent copies) |
| GWFAM18 | 11 | 12 | 16 | 11 | LOC115932756 |
| GWFAM186 | 9 | 21 | 22 | 12 | LOC129529768 |

**Conclusion:** the RNA famCN tracks the assembly copy number (~70–72 % within 25 %), and **totalCN recovers the
divergent copies** that RNA expression and famCN@0.98 alone miss (GSTM: 8 expressed / 13 near-identical / **21
total ≈ 19 assembly**) — Soto's famCN/parCN, produced reference-free from long RNA + genome projection.

---

## What RNA uniquely adds (the contribution)

- **Divergent families** (down to ~0.60 nucleotide / protein-level) that Liftoff (`-sc 1.0` = identical) and
  SEDEF SD98 (≥98 %) cannot reach — e.g. the KRAB-ZNF cluster (0.62–0.74 identity, 30-copy GWFAM8).
- **Expression** — which of the genomic copies are transcribed (the RNA copies), not just present.
- **Reference-free** — families and famCN are built from de-novo RNA consensus + genome projection, without
  reference annotation, so no annotation-centric bias (unlike Liftoff, which is annotation-projected).

## Caveats (honest)

- The oracle→gene join is imperfect (RefSeq naming; 51 of 169 oracle multi-copy families matched to an RNA
  family). `asm_hapCN` is a diploid/assembly count, so ±25 % tolerance reflects real haplotype/assembly slop.
- 80 RNA-only families (Axis B) are not independently confirmed here — some are genuine divergent/expressed
  families, some may be over-merges; the DNA/annotation witnesses cannot adjudicate the divergent tail.
- `--enumerate-copies` in the binary now uses a single-pass batch projection (the earlier per-family genome
  re-indexing was a perf bug, fixed). famCN@0.98 requires full-length (cov ≥ 0.90) to avoid fragment inflation.

## Reproduce

```
# RNA-98 catalog + famCN
gw_family_catalog --bam GGO_mm.bam --fasta GGO.fasta --homology-primary --min-identity 0.98 \
  --enumerate-copies --out rna98
# Axis A head-to-head
python bench/sim_liftoff_headtohead.py     # + liftoff -copies -sc 0.85 + minimap2 -x splice -p 0.01
# Axes B + C
python bench/crossmodal_validate.py rna98   # SEDEF SD98 = final.bed col21>=0.98 (gorilla), oracle = diploid_cn_oracle.tsv
```
