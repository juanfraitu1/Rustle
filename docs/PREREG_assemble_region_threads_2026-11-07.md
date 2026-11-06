# PREREG row 1564: run the assemble stage's region sweep in parallel (`--region-threads`)

Follows row 1560 (assemble 2,588 s on A119b, node ~7 of 72 cores: the stage passes no `--region-threads`, so the
25 per-contig regions are fetched and assembled one at a time in each of the two concurrent copy_assign runs).

## Change (tested on a copy of the pipeline; the repo script is edited only after a PASS and after the
## r1550 e2e jobs, which read the live script, have finished)
Both assemble copy_assign calls (main + final_run) get `--region-threads 8`, with
`RUSTLE_REGION_READ_BUDGET=20000000` for the assemble calls only (the pipeline's 4M default is sized for assign,
OOM 250 GB at 4.4M records there; at 4M, chr1/chr2/chr13 (5-6M reads each) would each run alone and the sweep
would stay near serial). The per-contig polish loop (copy_assign_bin.rs:10532) is NOT touched: it runs after the
sweep, outside `--region-threads`, and stays serial (~1,170 s) — that is the next row.

## Run
Same binary as r1560 (HEAD 1b699e7b build, /anvme/workspace/v138da13-rustle/r1560/target/release), same BAM,
genome, defaults, RUSTLE_TIMING=1, singlenode, same line stamps + 10 s CPU sampler, wall limit 50 min.
Only the flag and the budget differ from r1560.

## Bars (all must hold)
1. Identity: a1.gtf, a1.families.gtf, a1.final.gtf, a1.bridges.tsv, a1.bridge_junctions.tsv, a1.quant.tsv
   byte-identical to r1560.
2. Sweep: as_table end -> last `reads_in_region` line <= 400 s (r1560: 1,118 s).
3. Stage wall <= 1,900 s (r1560: 2,588 s).
4. Peak RSS (sum of user processes) <= 120 GB.

## Prediction
Identical products; sweep 250-400 s; stage wall ~1,700 s (as_table ~190 + sweep ~350 + polish loop ~1,170 + ~90);
RSS 30-60 GB. The serial polish loop then becomes ~70% of the stage.
