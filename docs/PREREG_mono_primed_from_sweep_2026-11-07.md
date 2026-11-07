# PREREG row 1566: `--polish-mono-primed` evidence from the sweep, not a second BAM pass

## Why (r1565)
Assemble is 1,024 s; the mono-primed BAM pass is still ~317 s of each run's serial polish loop. The sweep already
decodes every record of each contig (streaming pass 1), then the polish loop re-reads the same records from the BAM.

## Change
- TssEvidence gains `mono` (chromosome -> every primary, non-QC-fail, unspliced record as [a, a + M/D/=/X span), NOT
  deduplicated: the BAM pass counts each record) and `mono_window` (chromosome -> its one fetch window; a second window
  or region marks it unusable). Collected in `push_alignment` (called on every record of the streaming reader before
  the pool's rules) when `--polish-mono-primed` is on; the buffered path does not collect.
- Polish: a contig fetched in ONE window that covers [min start, max end] of its single-exon models is counted from
  memory; any other contig takes the old BAM pass. Both use one counting function (`mono_primed_count`), so the
  per-record rule is shared. The MONO PRIMED log line names the source (sweep | bam).
- Exact by construction under the guard: one window covering the models' span returns every record the BAM query
  returns, once.

## Run
r1565 harness (shipped pipeline, singlenode, wall limit 40 min), build of main + this change. Unit test
`mono_evidence_from_sweep_counts_every_record_and_guards_windows`; copy_assign tests 102/102.

## Bars (all must hold)
1. a1.gtf, a1.families.gtf, a1.final.gtf, a1.bridges.tsv, a1.bridge_junctions.tsv, a1.quant.tsv byte-identical to r1565.
2. Every MONO PRIMED line says `evidence sweep` (genome-wide: one whole-contig window per contig).
3. Stage wall <= 850 s (r1565: 1,024 s).
4. RSS <= 120 GB.

## Prediction
PASS; mono-primed ~317 s -> < 10 s per run; stage wall ~720 s; RSS +1-3 GB.
