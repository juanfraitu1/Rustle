# The K=0 wall is a property of the READ, not the COPY

**Date:** 2026-07-09. **Sim:** `bench/sim_k0_flank.py` (planted, deterministic). **Question (advisor):** for
exonically-identical copies (K=0), is there *absolutely nothing* in the BAM that determines a read's provenance?

## Design — perfect coverage isolates the wall

One chromosome, four co-located copies of one gene, **all reads full-length (perfect coverage)** so no read can
fail for lack of coverage:

| copies | exons | introns / 3′ flank | expectation |
|---|---|---|---|
| A, B | diverged (d=0.01) → **exonic PSVs** | diverged | assignable |
| **C, D** | **IDENTICAL (d=0)** → **no exonic PSV** | diverged (introns d=0.02, flank d=0.05) | **K=0** |

Two read classes per copy (60 each): **`_exon`** = the spliced transcript only; **`_flank`** = the transcript
plus that copy's own 3′ flank (readthrough). Aligned `minimap2 -ax splice -N 50` (secondaries, as in GGO_mm.bam).

## Result

Per **unique read** (placements deduplicated):

| true copy / class | reads | assigned | tied |
|---|---|---|---|
| A_exon | 60 | 60 (100%) | 0 |
| **C_exon** (identical exons, no flank) | 60 | **0** | **60 (100% tied)** |
| **D_exon** | 60 | **0** | **60 (100% tied)** |
| **C_flank** | 60 | **60 (100%)** | 0 |
| **D_flank** | 60 | **60 (100%)** | 0 |

**Separation is exact:** all 60 assigned `C_flank` → copy 1; all 60 `D_flank` → copy 2. Zero cross-assignment.

**Direct information-content test** (pipeline-independent: align only the read's flank portion to flank_C vs
flank_D): **40/40 = 100%** recover the true copy.

## Conclusions

1. **With perfect coverage, K=0 is the ONLY wall.** No read abstained for coverage reasons; the tied mass is
   exactly the exon-only reads of the exonically-identical copies. This isolates the two causes of abstention
   observed on real data — (i) K=0 (no distinguishing column exists) vs (ii) coverage (columns exist, reads
   don't reach them) — and shows (ii) vanishes under full-length reads.
2. **The wall belongs to the READ, not the COPY.** The *same* K=0 copies are perfectly separable from reads that
   carry flanking sequence. For a read confined to the identical exonic span, provenance is
   information-theoretically absent (100% tied, as it must be). For a read that reads through into the divergent
   flank, provenance is fully recoverable (100% assigned, 100% correct).
3. **The pipeline already exploits it.** C and D exist as *separate copies at all* only because the flank-bearing
   reads extended the de-novo transcripts into the divergent flank, yielding 56 PSV columns there. Nothing new
   had to be built — the flank became exonic in the de-novo model.

## Bearing on real GGO data

47.7% of real reads in a sampled family region carry **≥20 bp soft-clipped sequence — median 162 bp, 0% polyA,
0% low-complexity, 100% real mixed sequence** (`samtools view` audit). That is exactly this flank material,
currently discarded by exon-column PSV calling. The sim says: wherever a read carries it, a K=0 copy is
recoverable. **Estimating how much of the real ~25% tied mass is flank-bearing (rescuable) vs exon-confined
(provably lost) is the natural follow-up** — it converts a single "abstained" number into "X% provably
impossible, Y% recoverable with the sequence already in the BAM."

## Honest caveats

- **B was excluded from the family** (its reads map at MAPQ 60 → uniquely placed → no de-tie conflict edge).
  This is the known `E_c`-no-edge limitation of the conflict-graph family definition — the reason the `E_r`
  homology-primary definition exists. A family-*definition* artifact, not an assignment failure.
- Planted flank divergence is 5% over 300 bp — generous. Recent duplicates may have near-identical flanks, in
  which case the flank carries less (or no) signal. The real-data soft-clip audit shows the sequence is *there*;
  whether it *discriminates* on GGO is the open measurement.
- Introns diverge too (d=0.02) but spliced reads carry no intron sequence (`N` skips it); only intron-**retained**
  reads would expose intron PSVs. Junction *positions* are already used (O3).

## Reproduce

```
python bench/sim_k0_flank.py          # builds genome + reads + BAM, prints the direct flank test
copy_assign --bam .../k0.bam --fasta .../k0.fasta --region <printed> --min-copies 2 --skip-poa-diagnostic --out .../k0
```
Related: `bench/READONLY_COPY_NUMBER.md`, `bench/em_consistency.md` (the K=0/SoftZone floor), `project_k0_frontier_unresolvable`.
