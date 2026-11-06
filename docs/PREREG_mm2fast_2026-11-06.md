# PREREG: mm2-fast as the read-star aligner (row 1549, 2026-11-06 18:45, user-requested)

## Question
Read-star minimap2 chaining is 92% of chr13's CPU. mm2-fast (Intel, github bwa-mem2/mm2-fast, 14fe36c, reports
2.24-r1122) vectorises chaining with AVX-512; Fritz singlenode is Ice Lake (avx512bw). We run minimap2 2.30-r1287.
Does swapping the binary (`RUSTLE_MINIMAP2`, no code change) cut the panel wall with identical tables?

## Arms (hot9 panel, singlenode, hot10 recipe, bin10 flag unset, staged BAM, wall limit 1 h)
- base: h9out/auto products (1,654-1,664 s on this recipe, 104 Assigned).
- mm2fast: `RUSTLE_MINIMAP2=/anvme/workspace/v138da13-rustle/r1549/mm2-fast/minimap2.avx512bw`.
## Bars
- Five tables identical to h9out/auto (the version gap 2.24 vs 2.30 is the expected failure mode, not the SIMD).
- Panel wall <= 1,100 s (chaining-bound: biggest family 1,372 s at 16 threads).
## Prediction
Coin flip on identity: minimap2 2.24 -> 2.30 changed chaining defaults for some presets; if `-x splice` output
differs the row FAILS on identity regardless of speed, and the lever is closed unless a 2.30-based mm2-fast
exists. If identical: wall 900-1,200 s. Tie-break if tables differ: report the assignment diff size, no adoption.
