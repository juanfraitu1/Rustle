# The family definition under proper multimapper sampling (new -N 50 -p 0.1 BAM)

The definition's single biggest empirical caveat was that `~R` (read-confusability) was built on
`GGO.bam`, aligned with minimap2's **default secondary cap** — so most cross-mapping was never
emitted, and `~R` was measured on a small fraction of the real read confusion. A re-alignment
with `-N 50 -p 0.1 --secondary=yes --eqx -Y` (`GGO_mm.bam`) surfaces the suppressed secondaries
(de tags intact). This re-validates the definition on the data it should always have had.

Method (`bench/family_def_newbam_validate.py`): per chromosome, re-run the **exact** `~R` scan
(`family_def_genomewide.scan` logic, region-restricted) on OLD vs NEW, then apply `~B` (cached /
on-demand-built exon-union copy models, reciprocal coverage ≥ 0.30). `~B` copy models are
BAM-sampling-independent (built from primary alignments), so the comparison isolates the `~R`
effect.

## Result (three paralog-dense chromosomes)

| chrom | secondaries OLD→NEW | multimapper reads | `~R` candidate edges | **`~B` real-copy edges** | `~B`-pruned bridges |
|---|---|---|---|---|---|
| NC_073244.2 (RFPL) | 4,818 → 313,883 (65×) | 751 → 5,113 (6.8×) | 18 → 33 | **2 → 2** | 16 → 31 |
| NC_073248.2 (chrY-like) | 5,732 → 92,889 (16×) | 182 → 353 (1.9×) | 28 → 49 | **21 → 35** | 7 → 14 |
| NC_073240.2 | 22,309 → 258,319 (12×) | 2,487 → 10,821 (4.4×) | 68 → 128 | **58 → 69** | 10 → 59 |

## What this shows

1. **The undercount caveat was real and large.** The old BAM emitted 12–65× fewer secondaries and
   1.9–6.8× fewer multimapper reads than a proper `-N 50 -p 0.1` alignment. `~R` was genuinely
   undersampled.

2. **`~R` gains recall — it recovers hidden real families.** On the *dispersed*-paralog
   chromosomes, `~B`-validated real-copy edges **grow** (chrY 21→35, +67%; NC_073240 58→69, +19%):
   genuine paralogs whose cross-mapping the default cap had hidden, now detected. On RFPL the count
   is stable (2→2) because its copies are co-located/near-identical and already cross-mapped even
   under the old cap — exactly where the secondary cap did *not* bite.

3. **`~B` remains the precision filter, and it scales.** Every extra candidate edge the richer
   sampling surfaces beyond the real copies is a bridge, and `~B` prunes it (bridges 10→59 on
   NC_073240). Without `~B`, "candidate families" would balloon with the sampling depth (5→14 on
   RFPL); with `~B`, the validated set tracks real copies, not read counts.

**Conclusion.** The two-relation design is vindicated by the data it was missing: `~R` supplies
recall (now properly sampled, recovering hidden dispersed paralogs), `~B` supplies precision
(load-bearing, absorbing the extra bridges that richer sampling inevitably surfaces). The
definition is not merely robust to multimapper sampling — it improves with it, because `~B`
anchors family identity on copy-model homology rather than on read counts.

## Honest scope
- Three chromosomes, not genome-wide (the full `~R` scan on the 11.7 GB BAM is memory-heavy; a
  per-chromosome serial genome-wide pass is the natural next step).
- "Real copy" = `~B` pass (reciprocal coverage ≥ 0.30 of the exon-union models), the note's copy-
  level metric (validated as the correct copy-level criterion in `family_def_metric_compare.py`;
  contiguous-core is too strict on exon-union models — it fragments real copies).
- `~B` copy models reused from the old-BAM primary alignments (sampling-independent by design);
  on-demand models for new candidate loci built identically from primary reads.
