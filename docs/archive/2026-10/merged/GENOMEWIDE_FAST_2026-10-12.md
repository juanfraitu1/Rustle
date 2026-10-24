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
- A first fix, dropping single-placement molecules' records at read time, was WITHDRAWN: inside multi-copy families almost every read
  has several placements (chr1 region: 1,174,977 of 1,244,342 records kept), and it hid copies supported only by unique reads from
  the per-copy read-presence check (whole-gorilla run: MCL188 failed).
- **`copy_assign --tied-prefilter`** (driver default; byte-identical): each region is read twice — first LIGHT (name, scores, span;
  no sequence, qualities or full CIGAR, ~100 bytes a record), the AS-tied molecules decided exactly as the AS-tied gate does
  (`as_evidence_per_read` + `as_tied`), then FULL records loaded for those molecules only; the per-copy read-presence check reads the
  light pass's spans. The gate keeps exactly those records anyway, so families, assignments and abundance are unchanged (the gate's
  log counts then cover tied molecules only). Memory follows the tied reads, not every read in the copy windows. **Checked on
  gorilla KB3781 genome-wide (326 families): all 7 output tables cmp-identical (assignments, families, quant, family_join,
  famcn_readonly, copy_conversions, params); peak 5.2 -> 3.6 GB, wall 572 -> 639 s.**
- **`--max-window-records N`** (driver; `--low-mem` = 500,000): `utilities window-depth` estimates each family's records per merged
  copy window from the BAM index (compressed bytes / mean bytes per record; 0.8 s on human) and the families over N are not assigned
  (`--skip-families`) and their regions are not swept. On human A119b the 62 such families all lie entirely on the acrocentric short
  arms (chr13/14/15/21/22 < 17 Mb; e.g. MCL1, a 3.5 M-record hull). Only reads those families share with others can change.

## Status and what is left (2026-10-12, PAUSED while the PSV engine is being fixed in another session)

**Done and committed** (branch `machine2/soto-evidence`): `03184537` timed, resumable `all` (per-stage `PREFIX.<stage>.time`,
`PREFIX.all.env` resume, final stage / wall / peak-memory table), `index` stage, `--low-mem`; `a36b9251` `flag` against the genome
FASTA indexed on the fly in 1-Gbp parts (`missing_copy_flag --index-gb`, `--split-prefix` always); `98119eeb` + `a30ca851`
`utilities window-depth` (index estimate, 0.8 s on human) + driver `--max-window-records`, and `copy_assign --tied-prefilter`
(driver default). `--skip-single-placement` was withdrawn in `a30ca851`.

**Measured so far**

| what | result |
|---|---|
| assign, POA diagnostic skipped vs not (gorilla KB3781, 7 + 20 families) | cmp-identical tables |
| assign `--tied-prefilter` vs not (gorilla KB3781 genome-wide, 326 families) | 7 tables cmp-identical; peak 5.2 -> 3.6 GB; 572 -> 639 s |
| human assign, old loading | killed at 16-20 GB (chr1 region 1.24 M records = 8.7 GB; chr13 short arm another ~7 GB) |
| human assign, `--low-mem` + prefilter | stopped by hand at 1.5 min (1.0 GB; the old loading was at 8.7 GB by then); never finished |
| `all` on one gorilla contig (NC_073244.2, `--no-flag`) | works; second run skips every stage |

**Left to do, in order, once the PSV engine is settled**

1. Rebuild and test (the full lib suite has NOT been run since `98119eeb`):
   `CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo build --release` and `cargo test --release --lib`
   plus `cargo test --release --bin missing_copy_flag` (output to a file). If the PSV engine default changed, re-run the gorilla
   identity check of `--tied-prefilter` (the prefilter is upstream of the PSV code, so it should still hold):
   copy_assign on `/mnt/linuxdisk/tmp/gw_fast/gkb/regions.txt` with and without `--tied-prefilter`, then `cmp` every table.
2. Human assign with the new loading (scratch `/mnt/linuxdisk/tmp/gw_fast/human/`, the copy table symlinked from the run cache):
   `tools/rustle_pipeline.sh assign --bam .../A119b.t2t.bam --fasta .../chm13v2.0.fa --out /mnt/linuxdisk/tmp/gw_fast/human/human
   --bin .../rustle_target/release --threads 4 --low-mem` under a watchdog that stops `copy_assign` above 16 GB (sample
   `ps -C copy_assign -o rss=`; do NOT use `pgrep -f` with a pattern that also matches the watcher's own command line). Record wall
   and peak; the 62 skipped families are listed in `human.assign.depth_skipped.tsv`.
3. `flag` on gorilla, both index paths: the `.mmi` (`--index /mnt/linuxdisk/home/juanfraitu/winloci_data/GGO.splice.mmi`, ~16-19 GB
   expected) and `--low-mem` (FASTA, 1-Gbp parts); compare the two verdict tables (near-identical expected, not byte-identical) and
   record wall / peak.
4. The whole pipeline in one command: `tools/rustle_pipeline.sh all --bam ... --fasta ... --out OUT --threads 4 --low-mem` on gorilla
   OR6737, then human A119b; copy the final table into this document. Disk: ~5 GB free on `/mnt/linuxdisk` (2026-10-12); the human
   run writes ~5-8 GB (loci FASTA + PAF + GTFs) — free space first.
5. Human `families` (45 min, 12.6 GB, the largest remaining cost): profile whether the minimap2 all-vs-all of 2 GB of locus spans or
   `mcl_families`' PAF handling holds the memory; anything that changes the PAF (intron masking, exon-only sequences) changes the
   validated family definition and needs a pre-registration.
6. Then: a register row, the memory note, and the stage table at the top of this document.

**Caveats to carry**: `--max-window-records` changes results for the skipped (acrocentric-array) families and for reads they share with
other families; under `--tied-prefilter` the AS-tied gate's log counts cover tied molecules only (tables unchanged); a multi-part splice
index is near-identical, not byte-identical, to one part.
