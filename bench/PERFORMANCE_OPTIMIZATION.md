# Copy-assignment pipeline speedup (2026-06-27)

Goal: make the genome-wide copy-assignment sweep faster **without changing a single assignment** (the O2
headline is byte-frozen). Every change below was verified output-identical by re-running real regions and
diffing the assignments + abundance (`diff` line count 0), plus the 659-test suite stays green.

## Where the time actually went (measured, not guessed)

Profiled the heaviest real family (GAGE, `NC_073228.2:144928569-145085403`, 20,067 mapped reads → 8,461
assignments, 5 copies, 3,238 PSV columns) with env-gated stage timers (`RUSTLE_TIMING=1`):

| stage | time | note |
|-------|------|------|
| `reads_in_region` (BAM I/O, 20k reads) | 0.6 s | cheap |
| pass1 + gate + collapse | 0.0 s | cheap |
| **`build_read_placements` (20k × 10 reps)** | 0.0 s | the audit feared O(B·S); it's negligible |
| **POA homology *diagnostic*** | **~458 s** | dominant — and *purely diagnostic* |
| **`discover_psvs` (poasta, 4 pairwise DP)** | **69.7 s** | the real compute after the diagnostic is gone |
| per-read assignment (8,461 reads) | 1.1 s | already parallel; was never the bottleneck |

The original audit's asymptotic guess (per-read assignment is the O(R·C·N) hot path) was **wrong in
practice**: assignment is ~1 s. The wall-clock is the two alignment stages.

## The changes (all output-identical)

1. **Skip the POA homology *diagnostic* on sweeps — `--skip-poa-diagnostic`.** Families come from the de-tie
   conflict graph; the POA pairwise pass only writes `.fallback.tsv` and the diagnostic edge counts. It was
   left ON in the definitive run and cost ~85% of wall-clock. Skipping it: **537 s → 79 s (6.8×)**,
   assignments + abundance byte-identical. (Wraps the existing `RUSTLE_SKIP_POA_DIAGNOSTIC` env as a real,
   documented flag.)
2. **Parallelize `discover_psvs`** (`copy_assign_pipeline.rs`). The (n−1) ref-vs-copy poasta alignments are
   independent; run them concurrently and merge (per-index amaps + set-union of diffs = identical result).
   poasta is a pure function (thread-safe); the opt-in minimap2 path stays serial (its temp-file nonce can
   collide between equal-length copies under threads). **79 s → 43 s**, peak RSS 1.66 GB (safe on the 19 GB
   box), diff 0.
3. **Parallelize the per-read assignment loop** + **iterate only spanned PSV columns** in the significance
   gate and the editing pre-pass (`copy_assign.rs`, `copy_assign_pipeline.rs`). Reads are independent;
   rayon's indexed collect preserves order. The gate previously re-scanned all 3,238 columns per competitor
   per read — now it walks only the ~tens a read spans. Marginal on the heavy region (assignment is ~1 s
   there) but real on small/dense families and at cluster scale.
4. **`BamIndexCache`** (`denovo_assemble.rs`): parse the multi-MB `.bai` + header ONCE and reuse across all
   regions instead of per-region. ~free on 79 regions; compounds over thousands at genome scale.

**Heaviest region: 537 s → 43 s = 12.5×, byte-identical.** End-to-end 79-region sweep (`o2_verify`):
**3,854 s vs the old 8,130 s = 2.1×**, and the full result is **byte-identical** — all 74 families / 104,147
assignments diff = 0, abundance diff = 0, silver 24,660/24,682 = 99.9% unchanged. (The 2.1× is a *conservative
lower bound*: that sweep ran contended with the A\*/old-binary benchmark runs competing for the 5 cores; the
clean per-region 12.5× is the reliable figure. A clean re-run would land between the two.)

## Also tested, gated OFF (changes output → revalidate-grade, not a free speedup)

**Exact-A\* (`RUSTLE_POA_ASTAR=1`)** swaps poasta's `AffineDijkstra` for `AffineMinGapCost` (same optimal
*cost*, admissible heuristic, fewer explored states). Faster, but A\* returns a *co-optimal-but-different*
traceback: on the heavy family it changed **1,254 / 8,461 assignments** and found 3,239 vs 3,238 columns. So
it is **not** byte-identical and is gated **off** by default (the O2 result stays exactly poasta-Dijkstra). The
same co-optimal-traceback hazard applies to banded/WFA aligners — so the only *truly* output-identical
`discover_psvs` speedup is the parallelization above; going faster means accepting a small, revalidate-grade
output change (benchmark vs the poasta 3,238-col / silver-99.9 headline before adopting).

## Rejected (changes output)

**minimap2 PSV discovery** (`RUSTLE_PSV_MINIMAP2=1`) is ~270× faster (43 s → 2 s) but is a heuristic,
coarser aligner: on this family it found **109 PSV columns vs poasta's 3,238**, and **all 8,461 assignments
changed**. It trades accuracy for speed and would move the O2 headline, so it is NOT used for the frozen
result. (It remains a valid opt-in for a fast/approximate pass.)

## Recommendation

Genome-wide / cluster sweeps: run `copy_assign … --skip-poa-diagnostic`. The parallel `discover_psvs` and
per-read assignment are automatic. Keep poasta (default) for the PSV columns — it is the accurate engine the
O2 numbers are computed with. Use `RUSTLE_TIMING=1` to re-profile if a region is unexpectedly slow.
