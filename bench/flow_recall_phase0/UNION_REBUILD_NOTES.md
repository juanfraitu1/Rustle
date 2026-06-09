# In-flow union-baseline re-bundle — notes

## Task 4 integration (NC_073224.2 rna-XM_063708549.1) — slice can't reproduce
Both OFF and ON byte-identical on the slice; flag didn't trigger (the over-collapse needs the
full-chromosome family/secondary context, absent in a narrow slice — consistent with the 94%
`context_only` attribution). Recovery validation must be genome-wide.

## Task 5 genome-wide validation — RESULT: NET-NEGATIVE

5 high-regression chromosomes, full VG config (`--vg --vg-snp` + TANDEM + DECISIVE_GATE),
opt-in `RUSTLE_VG_UNION_BASELINE=1`, gate `MIN_LONGCOV=2`:

| chrom | regressions | recovered | FP added | lost |
|---|---:|---:|---:|---:|
| NC_073247.2 | 17 | 3 | 2 | 0 |
| NC_073229.2 | 15 | 0 | 2 | 0 |
| NC_073242.2 | 14 | 0 | 81 | 0 |
| NC_073224.2 | 13 | 0 | 91 | 0 |
| NC_073228.2 | 10 | 0 | 39 | 0 |
| **TOTAL** | 69 | **3** | **215** | 0 |

**recall:FP ≈ 0.01 — catastrophically net-negative.** NC_073247 (3 real / 2 FP) was an
unrepresentative outlier. On the other four chromosomes the re-bundle recovers ZERO regressions
but unions 39–91 read-supported-but-NON-ANNOTATED chains that pass the `longcov>=2` gate.

### Diagnosis
The in-flow re-bundle (re-split secondary-bearing bundles' primary reads, assemble) produces many
read-supported-but-spurious chains (wrong-strand / fragmented sub-bundle assemblies) the longcov
gate can't separate from real recoveries — the same indistinguishability rock, on the GENERATION
side. build-fresh made the assembly *clean* (strictly additive, 0 lost) but not *correct* (re-split
fragments ≠ real baseline bundles).

### Verdict
- The in-flow fix is committed (default-OFF, default byte-identical, suite 277/0) but is **NOT a
  net-positive win** — stays opt-in/off, do NOT promote. "Ship in-flow now" is RETRACTED.
- Architecture is wrong: re-split fragments ≠ real baseline assembly. The validated path is the
  **post-process union** (separate clean full `baseline -L` run; union novel chains scoped to
  family/secondary-bearing spans; primary-gated) — measured recall:FP **1.27**, 17/17 on
  NC_073247. Uses the real baseline pipeline. Build that next.
