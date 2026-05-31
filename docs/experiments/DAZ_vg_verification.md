# DAZ multi-copy VG verification (2026-05-31)

Concrete end-to-end test of `rustle --vg` on a multi-copy gene family whose second copy is dominated by
multi-mapping reads — the thesis's core "do what StringTie can't" capability. Data: subset of GGO.bam
(gorilla, IsoSeq) on chrY `NC_073248.2:42,700,000–43,000,000`, with `--genome-fasta ../GGO.fasta`.

## Setup (matches the advisor design doc 2026-05-28)
- DAZ1: 42,783,133–42,859,657 (− strand) — 167 unique reads + 42 multi-mappers.
- DAZ3 (LOC129530216): 42,879,918–42,945,552 (+ strand, **inverted** vs DAZ1) — **0 unique reads, 216
  multi-mappers**. Confirmed in the subset: DAZ3 region = 216 MAPQ=0 + 1 MAPQ=1.
- DAZ1 and DAZ3 are ~20kb apart on opposite strands → SEPARATE bundles (unlike the same-bundle GOLGA6L7
  tandem), so VG mode can link them across the bundle boundary.

## Result — VG recovers the starved copy; StringTie gets a fragment

| At DAZ3 (multi-mapper-dominated copy) | transcripts |
|---|---|
| StringTie | 1 (truncated: 42,899,569–42,945,549, starts ~20kb into the gene) |
| Rustle baseline (no --vg) | 1 |
| **Rustle --vg** | **5 full-length isoforms** (span 42,879,743–42,945,549) |

The 5 VG isoforms have 25–27 exons (consistent with the DAZ repeat-array / tandem exon-7 structure),
dominant cov 112.8. VG gave DAZ1 *fewer* isoforms than baseline (9 vs 10) — it redistributed
multi-mappers toward the starved copy (correct direction). **Positive example: VG recovers the full
DAZ3 gene + isoforms where StringTie produces one truncated fragment.**

## Honest mechanism caveat (the part to disclose, not hide)
The recovery did NOT use the intended EM reweighting. VG diagnostics:
```
[VG] Discovered 1 gene family group (4 linked bundles, 195 multi-map reads)
[VG] graph-k-mer-Jaccard filter: 1 → 0 families (dropped, low_kmer_jaccard 0.000 < 0.05)
[VG] preserved secondaries on 4 bundles whose raw family was filter-dropped
```
The DAZ family was DISCOVERED but DROPPED by the k-mer-Jaccard family filter (jaccard 0.000) because
DAZ1/DAZ3 are **inverted** — forward k-mers don't overlap. DAZ3's reads survived only via a
secondary-preservation fallback, so:
- **Open #1 (fixable):** the k-mer-Jaccard family filter (`--vg-family-min-kmer-jaccard`, default 0.05;
  `vg.rs`) must compute reverse-complement k-mers for inverted paralogs, else EM never runs on them.
- **Open #2:** because EM didn't run, the shared multi-mappers aren't weight-split between DAZ1/DAZ3
  (the design doc's 53/47) — DAZ3's cov may be inflated. Running EM after the filter fix calibrates it.

## Counterexamples / limitations (from prior experiments, docs/experiments/GGO_family_vg.md)
- **AMY (amylase, chr1):** copies diverged → reads map uniquely (secondaries rare). EM not needed
  (Zone 1). Rustle ≈ StringTie.
- **GOLGA6L7 (chr19, tightest tandem, ~85kb):** L7_2 starved (6 primary / 70 secondary) → 0 predictions.
  VG mode does NOT help — the paralogs are close enough that the bundle-builder MERGES them into one
  bundle, and VG operates at bundle boundaries (same-bundle paralogs invisible). Needs an intra-bundle
  sub-cluster step.
- **chr19 whole:** VG finds 23 families, reweights 416 reads, but is −4 matches vs baseline; EM
  converges trivially (delta=0) — most multi-mappers have identical compatibility across copies (Zone 3).

## Advisor framing (examples + counterexamples)
- **It works (DAZ3):** a copy with 0 unique reads, recovered as 5 full-length isoforms vs StringTie's 1
  truncated fragment — separate-bundle, multi-mapper-linked, inverted paralog.
- **It doesn't yet (GOLGA6L7):** same-bundle tandem paralogs — needs an intra-bundle splitter.
- **It's not needed (AMY, X-Y pairs):** diverged copies map uniquely.
- **Honest uncertainty (Zone 3):** truly identical equally-expressed copies → EM gives 50/50.
- **Two concrete fixable gaps surfaced:** inverted-paralog k-mer Jaccard; EM weight-split validation.

## Reproduce
```
samtools view -b ../GGO.bam NC_073248.2:42700000-43000000 > /tmp/daz.bam && samtools index /tmp/daz.bam
./target/release/rustle --vg --genome-fasta ../GGO.fasta -L /tmp/daz.bam -o /tmp/daz_vg.gtf
./tools/stringtie/stringtie -L /tmp/daz.bam -o /tmp/daz_st.gtf
# transcripts at DAZ3 42879918-42945552 +strand: vg=5, stringtie=1
```
