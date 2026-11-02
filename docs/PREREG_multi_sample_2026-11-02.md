# PREREG MULTI-SAMPLE POOL (register 1503, 2026-11-02, before the run)

Feature: `rustle_pipeline.sh all --bam a,b` (tools/rustle_samples.sh; spec docs/superpowers/specs/2026-11-02-multi-sample-design.md).
Identity (gate, before the panel): chr20, A119b alone, driver at a29ac838 vs 9cd6313b (same binaries): every PREFIX.* product of
`all --to families` and of `assign` byte-identical.
Panel: CHM13, `--contig chr17,chr19`, arms A = A119b alone, H = HG002 Kinnex alone (`all --to families`), P = both pooled
(`all`, per-sample assign + counts); each arm under `timeout 7200`. gffcompare of each PREFIX.final.gtf vs chains/ref.gtf
restricted to chr17+chr19.
Bars: (1) for each sample, the PREFIX.counts.tsv column sum equals the log's compatible count (to 0.01 per 10^4 reads, the 2-decimal
rounding); (2) P `=` >= max(A `=`, H `=`); (3) P chain precision >= min(A, H). Report: per-sample assign agreement (same family and
copy, read in both) between P's shared-catalog assign and an assign on that sample's own panel catalog; walls.
Pass -> genome-wide run next. Fail on (2)/(3) -> recorded; the pooled mode stays available but documented as not better.
