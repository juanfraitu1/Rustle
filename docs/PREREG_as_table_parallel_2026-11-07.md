# PREREG row 1567: parallel `as_table` scan (ms columns kept)

## Why (r1566)
Assemble is 632 s; the best-AS table pass is 206 s of it (scan ~157 s at one record decoder, 18 of 64 cores; sort +
write ~38 s; binary sidecar 9 s).

## Change (src/as_table.rs)
- With a `.bai` and no `--records`: every reference cut into 10 Mb chunks, taken by threads/2 workers (own reader,
  one inflate thread each); a record is folded only by the chunk holding its alignment start (queries also return
  overlapping records that start earlier). Folds kept per worker in 256 name-hash shards, merged shard by shard in
  parallel. Rows sorted with a parallel sort and formatted in parallel blocks, written in name order.
- The fold's two order-dependent rules (LAST primary wins; FIRST best-AS secondary is the rival) now use a file-order
  key per record (chunk << 32 | ordinal; file order in a coordinate-sorted BAM, which a .bai implies), so partial
  folds merge to the sequential answer in any order. Sequential path (`--records`, or no index) keeps file order and
  reduces to the old fold exactly. Unit test: 2,000 random molecules, chunked + shuffled merge == sequential.
- Output format unchanged: all 17 columns, including best_ms / second_ms / primary_ms / rival_ms. The .asbin
  sidecar still holds best AS only (its one consumer, the seeding pool, reads nothing else).

## Run
`as_table --bam stage/A119b.t2t.bam --threads 64` on singlenode (the pipeline passes --threads 64), wall limit 15 min.

## Bars (all must hold)
1. The TSV byte-identical to r1566's a1.molecules.tsv (header records=67476966 molecules=21391966 included; every
   column including the four ms columns).
2. as_table process wall (scan + write + sidecar) <= 75 s (r1566: 206 s).
3. RSS <= 120 GB.

## Prediction
PASS; scan ~40 s (32 workers, I/O-bound share unknown), sort/write ~10 s, sidecar 9 s: ~60 s; RSS ~15-25 GB.
