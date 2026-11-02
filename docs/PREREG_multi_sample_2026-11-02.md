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

## Amendment 1 (2026-11-02, before any panel output; user: time)
All arms read slim BAMs: A119b chr20 (identity gate), A119b and HG002 chr17+chr19 (panel), sliced with samtools from the full BAMs.
The best-AS molecule table then sees only alignments on those contigs, identically for every arm. Bars unchanged.

## Identity gate result (job 4412353) and amendment 2 (before any panel output)
35 products compared, 32 identical. The 3 differing: `x.assign.params.tsv` and `x.fam.params.tsv` differ only in the lines that
record their own output paths (`id.old/` vs `id.new/`); `x.molecules.tsv.asbin` is the cache sidecar, which stores the TSV's mtime
and inode (AsTsvIdentity), while `x.molecules.tsv` itself is identical. Read as PASS on products; amendment 2: the gate compares
products, excluding `*.params.tsv` path rows and `*.asbin`. Panel run as job below.

## RESULT r1503 (2026-11-02, jobs 4413201 + 4416521): ❌ FAIL (bar 3); bars 1 and 2 PASS
| arm (chr17+chr19 vs RefSeq) | transcripts | = | chain sens / prec |
|---|---|---|---|
| A119b alone | 20,483 | 4,438 | 24.9 / 22.0 |
| HG002 alone | 17,945 | 4,618 | 26.0 / 25.9 |
| pooled | 29,807 | 5,957 | 33.4 / 20.2 |
(1) counts: A119b column sum 878,130.81 vs 878,131 compatible; HG002 2,352,151.44 vs 2,352,152 (2-decimal rounding) PASS.
(2) pooled `=` 5,957 >= 4,618 PASS (+29%). (3) pooled chain precision 20.2 < 22.0 FAIL: the pool is the union of both samples'
isoforms, noise included. Report: reads assigned on the shared catalog vs the sample's own: A119b 332 vs 60 (13 in both), HG002
1,391 vs 41 (23 in both), of the reads in both runs. First pooled run (4413201) stopped in HG002's assign: copy_assign refused
copies with no HG002 read (catalog built from both samples); fixed by RUSTLE_SHARED_CATALOG=1, set by rustle_samples.sh.
Verdict as pre-declared: the pooled mode ships, documented as more sensitive, not more precise.
