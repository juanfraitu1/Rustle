# PREREG row 1568: one best-AS table per BAM, shared by every run on it

## Why
The table depends on the BAM alone, but its default home was PREFIX.molecules.tsv, so every new --out rebuilt it:
~200 s of I/O on the 96 GB A119b BAM (row 1567: the pass is bound by the /anvme link, not CPU).

## Change (tools/rustle_pipeline.sh stage_assemble, good seeding pool)
- Default table BAM.as_table.tsv next to the BAM when its directory is writable (else PREFIX.molecules.tsv as before);
  --as-table FILE still names it.
- Reused when its header names this BAM and it is not older than the BAM; otherwise built under a temporary name and
  renamed into place with its .asbin (rename keeps the inode/mtime the sidecar identity uses; a concurrent run never
  reads a half-written table). A failed as_table now stops the stage (exit 3) instead of continuing.
- PREFIX.molecules.tsv (+ .asbin) becomes a symlink to the shared table (tools/rustle_reassemble.sh reads that name).
No binary change.

## Run (one job, singlenode, current main build)
A. no shared table; `assemble --contig chr22 --out A` (builds it); B. `assemble --contig chr22 --out B` (reuses);
C. genome-wide `assemble --out C` (reuses).

## Bars (all must hold)
1. The shared table byte-identical to r1566's a1.molecules.tsv.
2. B logs `reusing`, writes no B.as_table.log, and B's six products equal A's.
3. C's six products byte-identical to r1566 (a1.gtf, families.gtf, final.gtf, bridges.tsv, bridge_junctions.tsv,
   quant.tsv).
4. C stage wall <= 470 s (r1566: 632 s with the pass); B.molecules.tsv resolves to the shared table.

## Prediction
PASS; C ~430 s (632 - ~200 s of table pass).
