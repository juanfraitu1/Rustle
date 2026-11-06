# PREREG row 1560: where does the assemble stage's wall go (A119b, singlenode 72 cores)

Written before any number of this run was read. One number from an older log was seen while learning the log
vocabulary (gw_off a1.assemble.log tail: BRIDGE REGROUP 12.5 s); no other assemble timing was read.

## Run
- Binary: clean `git archive` of HEAD 1b699e7b (the repo working tree carries another session's uncommitted
  `src/` edits, so `target/release` is not used), built to /anvme/workspace/v138da13-rustle/r1560/target.
- `tools/rustle_pipeline.sh assemble --threads 64 --no-cache` on /anvme/workspace/v138da13-rustle/stage/A119b.t2t.bam,
  chm13v2.0.fa, out /anvme/workspace/v138da13-rustle/r1560/a1, RUSTLE_TIMING=1, shipped defaults (no RUSTLE_*
  speed variables, as the e2e run r1543/e2e). The seeding pool is `good`, so the stage first builds its own
  best-AS table (`as_table`, one BAM pass; the gw_off table names a different BAM path and must not be reused,
  the stage would overwrite it).
- Instrumentation outside the pipeline (no script edits): every line of a1.as_table.log, a1.assemble.log and
  a1.final.log stamped with seconds since stage start (`tail -F`); a 10 s sampler of node-wide busy cores from
  /proc/stat and, per `copy_assign` process (main vs `final_run`), cores used (utime+stime delta) and threads.

## The four candidate sinks and the markers that bound them
1. BAM pass for the best-AS table: as_table start -> exit (serial before both assembly runs).
2. Region sweep (BAM read + per-region assembly, interleaved): main copy_assign start -> first `ASSEMBLY POLISH`
   line. The BAM-read part = sum of `[timing] reads_in_region` seconds (if the assemble-only path prints them);
   per-region assembly has no built-in timer under --assemble-only (`detect_and_assign` timing is gated off), so
   it is the remainder, attributed from the CPU trace, not from a timer.
3. Polish + TPM + bridge regroup + write: `ASSEMBLY POLISH` line -> main copy_assign exit.
4. final.gtf layer run: runs CONCURRENTLY with the main run (same inputs + 4 layers); its cost to the stage wall =
   final_run exit - main exit (if > 0), plus whatever it slows the main run by sharing the node.

## Scheduling-vs-CPU discriminator (fixed now)
A sink is a SCHEDULING problem (not CPU) if, over its interval, node-wide busy cores average < 36 of 72 (50%)
while the interval is >= 20% of the stage wall; or if one process sits at <= 2 cores for >= 50% of an interval
that is >= 10% of the wall (a serial tail). It is a CPU problem if busy cores average >= 54 (75%). In between =
mixed, reported as such.

## Prediction (structural, from the stage code, before the run)
`stage_assemble` passes neither `--region-threads` (binary default 1 = the serial region path) nor `--threads`
(default 4 BAM-reading threads) to either copy_assign run. Predicted: the region sweep is >= 60% of the stage wall
and is a SCHEDULING problem (busy cores well under 36: two serial processes, each ~1-5 cores); the as_table pass is
the second sink (CPU-bound, it gets --threads 64); polish third. If this holds, the lever is `--region-threads`
on the assemble call (the assign stage already passes it and it is documented byte-identical), tested as its own
row with a cmp of a1.gtf / a1.final.gtf.

## Report
Top three sinks in seconds and % of stage wall, each labelled scheduling / CPU / mixed by the rule above.
No bar to pass or fail: a profile. The prediction is scored right or wrong.
