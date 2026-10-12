# Running the pipeline genome-wide fast on a small machine (2026-10-12)

Goal (user, 2026-10-12): integrate and simplify the codebase so the pipeline runs genome-wide in a short time with limited
resources. Machine of record: WSL2, 5 cores, 25 GB RAM (Claude Code's background reaper kills jobs well before that).

## Baseline (recorded stage runs, `/mnt/linuxdisk/tmp/rustle_figures/runs/*/<stage>.time`)

| stage | gorilla OR6737 / KB3781 | human A119b | what dominates |
|---|---|---|---|
| assemble (streaming) | 1.7-3 min, 1.7-2.2 GB | 7.2 min, 2.9 GB | BAM decoding |
| families | 2 min, 6.3 GB | **45 min, 12.6 GB** | all-vs-all `minimap2 -x asm20` over 2 GB of locus spans (introns included; 4.5 M of the pairs it finds are later dropped for no exonic evidence) |
| assign | 9.4 min, 5.1 GB | never completed (reaped 2026-10-11); measuring | per-region read loading; cross-chromosome families (MCL3: 260 copies on 23 chromosomes) |
| flag | not run | not run | `minimap2 -x splice:hq` of a few thousand consensus sequences against the WHOLE splice index (12-14 GB `.mmi` loaded at once: the 16-19 GB peak) |
| splice index build | — | 7 min, 19-20 GB | one-part index of the whole genome |

## Targets

One command, resumable, genome-wide on the human sample: assemble -> families -> assign -> flag, every process <= 12 GB peak RSS
(a low-memory profile <= 8 GB), total well under 2 h on 4 threads; gorilla well under 30 min. Products byte-identical to the stages
run one by one, except where a documented, measured approximation is opted into.

## Work items

1. flag: always pass `--split-prefix` to minimap2 (identical on a one-part index — checked), so a MULTI-PART splice index
   (`minimap2 -I`) can be used; the driver builds one on request. A multi-part index is NOT byte-identical (toy test with paralogs:
   2% of PAF lines differ — MAPQ, a few chain ends), so it is an opt-in low-memory profile, measured against the one-part verdicts.
2. families: left unsharded by default — its all-vs-all PAF is already cached under `PREFIX.cache` (a repeated run replays it);
   `RUSTLE_MINIMAP2=tools/mm2_shard.sh` adds mid-run resume (cmp-identical) at the cost of a shard index on disk. Its 45 min on human
   is profiled below.
3. assign: measure on human; bound the cross-chromosome families' memory if they are the peak.
4. One entry point: `all` gets `merged`'s per-phase resume, per-stage time / peak-memory records and a resource summary.
5. Measure the whole run on human and gorilla and fill in the table above.

## How to run it (implemented 2026-10-12)

    tools/rustle_pipeline.sh all --bam READS.bam --fasta GENOME.fa --out OUT/sample --threads 4 [--low-mem] [--no-flag]

- `all` runs assemble -> families -> assign -> index (when no `--index` is given) -> flag, each stage as its own driver call under
  `/usr/bin/time`: `OUT/sample.<stage>.time` holds the wall seconds, the peak RSS of the stage's whole process tree and the exit
  code, and the run ends with a stage / wall / peak-memory table. Re-running the same command resumes: a stage whose product exists,
  whose record says exit 0 and whose settings match `OUT/sample.all.env` is skipped, until one stage runs (every later stage then
  runs too). `--no-resume` reruns everything.
- `assign` skips the POA homology diagnostic by default (assignments, families and abundance byte-identical — checked on gorilla
  KB3781: 7 and 20 families, cmp-identical; `RUSTLE_POA_DIAGNOSTIC=1` restores it) and takes `RUSTLE_REGION_THREADS` regions at
  once (default 1; byte-identical; memory grows with it).
- `index` builds `OUT/sample.splice.mmi` (`minimap2 -x splice -I <N>G`, `--index-gb N`, default 8 = one part, as the figures' index).
  `flag` always passes `--split-prefix`, so a multi-part index works.
- `--low-mem`: `--index-gb 1` (flag holds ~1 Gbp of index at a time) and one assign region at a time. The multi-part index is
  near-identical, not byte-identical, to one part (toy genome with paralogs: 2% of PAF lines differ in MAPQ or a chain end).
- `--no-flag` stops after assign (no index needed).

## assign on human: where the memory went and the fixes (2026-10-12)

- The default assign stage (POA diagnostic skipped) still reached 20 GB on human A119b. `RUSTLE_TIMING=1` showed why: each region
  loads every record in its copy windows (copy +- 50 kb) with full sequence and CIGAR, e.g. 1,244,342 records for chr1 (8.7 GB), and
  the acrocentric short arms hold millions more (chr13:5.2-16 Mb: 4,769,917 records, 76,512 of them primary).
- **`copy_assign --skip-single-placement MOLECULES.tsv`** (byte-identical): a record whose molecule has one placement genome-wide
  (`as_table` n_records == 1, ~70% of molecules) feeds the primary pool but is never held as a full record; the AS-tied gate drops it
  anyway. cmp-identical assignments / families / abundance on gorilla KB3781 (20 families). The driver passes it whenever the
  assemble stage's molecules table (or `--as-table`) matches the BAM.
- **`--max-window-records N`** (driver; `--low-mem` = 500,000): `utilities window-depth` estimates each family's records per merged
  copy window from the BAM index (compressed bytes / mean bytes per record; 0.8 s on human) and the families over N are not assigned
  (`--skip-families`) and their regions are not swept. On human A119b the 62 such families all lie entirely on the acrocentric short
  arms (chr13/14/15/21/22 < 17 Mb; e.g. MCL1, a 3.5 M-record hull). Only reads those families share with others can change.
