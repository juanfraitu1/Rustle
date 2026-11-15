# Transcript ends, read checks and the fraction rule: state on 2026-11-14

Summary of the session of 2026-11-14 (register rows 1717-1724, 1731, 1733). Figures are copied from the pre-registrations
named in each section; those files hold the bars and the full tables. Nothing here is committed yet.

## What changed in the pipeline

| change | status | switch | record |
|---|---|---|---|
| `PREFIX.earlyends.gtf`: `PREFIX.gtf` plus nested early polyA ends | default last step of `all` (side file) | `RUSTLE_EARLYENDS=0` | rows 1721, 1722; families check row 1740 |
| `PREFIX.capstarts.gtf`: `PREFIX.gtf` plus nested starts with a capped read cluster | default last step of `all` (side file), user 2026-11-14 | `RUSTLE_CAPSTARTS=0` | rows 1717, 1720 |
| `tools/check_reads.py`: read pre-flight at the start of `assemble` | default, warnings only | none | this file |

`PREFIX.gtf`, `PREFIX.final.gtf` and everything families and assign read are unchanged. Both side files are human-tested for
starts and four-sample-tested for ends; merged into one output they cost 0.007 chain precision on HG002 (bar 0.005), which is
why they stay separate files and are not folded into `PREFIX.gtf`.

## Reads: polyA tails and orientation

- Rustle does no polyA trimming and the four read sets need none: at most 0.12% of primary reads still end in a polyA run
  (first 300,000 primary reads of each BAM). The reads are FLNC; the BAM headers record only minimap2.
- Every read carries `ts:A:+` because the BAMs were made with `minimap2 -uf` (reads declared to be in transcript
  orientation). The reads are in fact oriented: reads with a leading polyT run are 0.01-0.14%, spliced reads whose introns
  are mostly antisense 0.02-0.09%, and no read has both.
- Weak spot: input that was not oriented or trimmed upstream. Aligned with `-uf`, reversed reads land on the wrong strand
  unnoticed; aligned without it, unspliced reads have no `ts` tag and are taken at read orientation.
- `tools/check_reads.py BAM FASTA [N]` reads the first 100,000 primary reads (about 5 s) and warns, never fails, when
  >= 1% of reads end in a polyA run, >= 1% start with a polyT run, >= 2% of spliced reads have mostly antisense introns, or
  >= 2% carry `ts:A:-`. Thresholds are about ten times what the four samples show. `--selftest` checks the logic; a chr22
  slice with 25 A appended to each read fires the tail warning. The two orientation warnings have not been seen on real
  data (no unoriented read set on disk).
- What already handles polyA artifacts after trimming: `RUSTLE_POLISH_MONO_PRIMED` (default on, internally primed
  single-exon models), `RUSTLE_POLISH_TES=off|tag|pas-end` (opt-in, A-rich 3' ends, `src/end_artifacts.rs`),
  `--polish-tes-extend on` in `PREFIX.final.gtf` (row 1488), and the `altends` / `earlyends` stages.

## Cap G starts (`capstarts`)

Rule: a nested start is added when >= 3 reads of the chain carry a non-templated 5' G within 10 bp of one position and they
are at least half of the informative reads there. Rows 1717 / 1720 (`PREREG_cap_g_starts_2026-11-14.md`,
`PREREG_alt_rules_2026-11-14.md`): nested-start pairs with both chains emitted 0.572 -> 0.702 (A119b) and 0.521 -> 0.700
(HG002), chain precision unchanged within 0.003.

The stage checks the library first (capped share >= 0.5 of reads starting in a first exon, <= 0.15 inside internal exons)
and writes nothing when the check fails. Both gorilla libraries fail, for different reasons (exploratory table in
`PREREG_cap_g_starts_2026-11-14.md`):

- OR6737: the extra G is not in the reads (89% of reads at annotated starts have none); it was removed before alignment.
- KB3781: an extra G is on about 77% of reads even when they start inside an internal exon; a "two or more G" rule gives
  0.455 vs 0.194, below both bars.

Neither is fixable by a threshold on these files. `all` reports a refusal as a note and does not retry it while `PREFIX.gtf`
is unchanged.

## Early ends (`earlyends`)

Rule: a nested early end is added when the chain's reads end at one position (>= half within 10 bp), a polyA signal lies
upstream and no A-rich run follows. Rows 1721 / 1722: early-end pairs with both chains emitted 0.458 -> 0.669 (A119b),
0.434 -> 0.630 (HG002), 0.671 -> 0.865 (OR6737), 0.561 -> 0.745 (KB3781); precision -0.001 to -0.006.

Chains added by the stage (row 1722 runs), which is why it became a default step:

| sample | chains before | added | of the added, in RefSeq | stage wall, 72 threads |
|---|---|---|---|---|
| A119b | 219,517 | 2,930 | 281 | 168 s |
| HG002 | 153,033 | 4,071 | 294 | about 380 s |
| OR6737 | 73,348 | 842 | 68 | about 95 s |
| KB3781 | 72,102 | 1,643 | 102 | about 100 s |

In `all` it runs after the last stage (`all_earlyends` in the driver), rewrites the file only when it is missing or older
than `PREFIX.gtf`, and a failure (for example no `pysam`) is reported without failing the run. Not run end to end on a real
sample through `all` yet; the off switch has no test.

Families non-regression (row 1740, `PREREG_earlyends_families_2026-11-14.md`): with the chains added to the families
input, A119b stays 0.902 / 0.808 and OR6737 0.514 / 0.200 (35 pairs, weak check); clusters tables byte-identical; in the copy
tables 2 representatives are replaced by an early-end chain (1 per sample) and A119b gains 2 and loses 1 extra-representative row. PASS, so merging the chains into `PREFIX.gtf` is eligible. Not done: the
merge (user decision) and an `assign` run on the changed copy table.

## What other tools do (bambu)

bambu groups reads by corrected intron chain, scores read classes with a classifier trained on the sample, and assigns by
compatibility with an EM; it does not model alternative starts or ends, and its junction correction removes small real
shifts of internal exon ends. Its one transferable idea, counting full-length matches of a short isoform separately from
fragments of a longer one, is what rows 1724 / 1731 / 1733 tested in our terms. Written from memory of the paper and
package, not re-read. Earlier benchmarks: `PREREG_bambu_benchmark_2026-11-12.md` (row 1641) and `PREREG_isoquant_benchmark_2026-11-12.md` (row 1646).

## The fraction rule and full-length reads (rows 1724, 1731, 1733): closed

Question: are chains dropped by `--polish-isoform-fraction 0.02` real when their own reads are full-length at both ends
(capped 5' cluster and peaked polyA 3' end)? HG002 held out; `PREREG_fraction_exempt_2026-11-14.md`.

| measure, HG002 | chains with the signature | control | annotated chains with the signature |
|---|---|---|---|
| annotated in RefSeq (row 1724) | 0.071 of 10,016 | 0.090 (top-read, no signature) | default precision 0.234 |
| 5' in a CAGE peak (row 1731, unannotated) | 0.907 | 0.799 | 0.879 |
| 3' at a polyA site (row 1731, unannotated) | 0.838 | 0.521 | 0.854 |
| all introns annotated (row 1733, unannotated) | 0.068 | 0.095 (read-count matched) | 1.000 |
| chain read in the other human sample (row 1733) | 0.374 | 0.394 | 0.599 |

All three rows FAIL as registered. The ends of these chains are real, their intron chains are no better than the rest of
what the rule drops, so there is no exemption. Left open, not tested: about 37% of all chains the rule drops recur exactly
in an independent individual and tissue; whether those are minor isoforms or systematic artifacts is a question about the
rule itself. No short-read junction set is on disk.

## Also closed this day

- Row 1719: sibling GT-AG rule for shifted internal-exon ends, 0.071 annotated against a 0.352 bar.
- Row 1723: `RUSTLE_JUNCTION_MAJORITY=1` recovers the lost non-canonical minor chains but adds chains at 0.014-0.021
  annotated; `--assembly-junctions strict` stays.
- Row 1718: `PREFIX.final.gtf` already carries the better 5' end (in a CAGE peak 0.864 HG002 against 0.417 for
  `PREFIX.gtf`); quote `PREFIX.final.gtf` for ends.

## Yardsticks on disk (human, CHM13)

refTSS CAGE peaks `r1716_sqanti/refTSS.hs1.bed`; GENCODE V44 polyA sites `r1719_polya/sites.hs1.bed` (both in the
workspace, lifted with UCSC `hg38ToHs1`). No gorilla lift of either. RefSeq chain identity understates real unannotated
ends about tenfold, so score ends by these two files and chains by RefSeq.

## Known driver oddities seen, not investigated

Repeating `all --to assemble --contig chr22` re-ran `assemble` on three of four repeats, and on the repeat where
`assemble` was skipped `--to assemble` was not honoured (families, assign, index and flag ran).
