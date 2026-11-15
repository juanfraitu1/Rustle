# Shifted ends of internal exons: rules tried and the `shiftdrop` stage, state on 2026-11-15

Summary of the sessions of 2026-11-14 / 15 under the user goal "derive empirical rules to improve assembly focusing on
shifted internal exons" (register rows 1725-1730, 1732, 1734, 1738, 1739). Figures are copied from the two
pre-registrations, which hold the bars, the full tables and the corrections made after review:
`PREREG_shifted_exon_ends_2026-11-14.md` (row 1725) and `PREREG_shifted_exon_rule_2026-11-14.md` (all other rows).
Nothing here is committed. The goal was cleared by the user on 2026-11-15.

## Outcome in three lines

- No rule on this class improves the assembly at the level row 1677 asks for (+0.02 precision over a read cut of the
  same size). On human the data bound it: an oracle that drops every class chain not in RefSeq gains only +0.007 /
  +0.010 (A119b / HG002). On gorilla the oracle clears the bar (+0.022 / +0.023, OR6737 / KB3781), but the rules
  found reach +0.0006 to +0.003 over a read cut.
- One rule survives as a small, safe clean-up of the isoform catalogue on all four samples: the indel-at-junction
  rule, shipped as the opt-in stage `shiftdrop` (row 1739).
- The recall side was closed earlier (rows 1714 / 1719 / 1723: the default already emits both isoforms of 0.90-0.98
  of annotated shifted-site pairs).

## What changed in the pipeline

| change | status | switch | record |
|---|---|---|---|
| `PREFIX.shiftdrop.gtf`: `PREFIX.gtf` without weaker shifted-site chains whose junction reads carry an indel | opt-in stage, named only, never part of `all` | `tools/rustle_pipeline.sh shiftdrop --bam B --fasta G --out PREFIX [--threads N]` | rows 1734, 1739 |

`PREFIX.gtf`, `PREFIX.final.gtf` and every default are unchanged; no other stage reads `PREFIX.shiftdrop.gtf`. The
stage needs `PREFIX.gtf` (run `assemble` first), the BAM, and python3 with pysam (`RUSTLE_PYTHON`); it uses no
annotation. Log: `PREFIX.shiftdrop.log` (last line = summary). Wall on 72 threads: 42 s OR6737, 95 s KB3781, 172 s
A119b. To remove it: delete `stage_shiftdrop`, its usage lines and dispatch entry in `tools/rustle_pipeline.sh`, the
test `shiftdrop_needs_the_assembly_and_then_writes_the_shiftdrop_gtf`, and `bench/famcmp/shiftdrop.py`.

## Definitions

- **S1**: an emitted chain that is the weaker copy (fewer reads; ties: the later chain) of a stronger emitted chain
  on the same contig and strand, equal except one intron that shares exactly one boundary: one internal exon end
  moved.
- **S2**: equal except two adjacent introns, the exon between them moved at both ends and still overlapping.
- **S1t**: the moved boundary belongs to a terminal exon. Reported apart, never touched by any rule.
- **Rule A, frame rule (script only)**: S1 / S2, shift <= 10 bp, not a multiple of 3, weaker chain under 0.05 of the
  partner's reads.
- **Rule B, indel rule (the stage)**: S1, and at least half of the primary reads carrying the moved intron have an
  insertion or deletion within 10 bp of either end of that intron (at most 300 reads per intron).
- **Comparator**: the RefSeq share of genome-wide 2-read chains (0.170 OR6737, 0.154 KB3781, 0.066 A119b, 0.066
  HG002). A rule is compared with a plain read cut that drops the same number of chains.
- **Yardsticks**: chain in RefSeq; moved intron annotated (RefSeq, plus CAT on CHM13, plus the T2T primates CAT
  consensus on gorilla since row 1738); row 1714's annotated shifted-site pairs (both isoforms emitted).

## The class (rows 1725-1727; gorilla intron share from row 1738)

| | OR6737 | KB3781 | A119b | HG002 on CHM13 |
|---|---|---|---|---|
| S1 chains | 11,838 | 13,511 | 26,591 | 25,880 |
| S1 in RefSeq (whole assembly) | 0.276 (0.352) | 0.292 (0.358) | 0.188 (0.176) | 0.238 (0.234) |
| moved intron annotated | 0.612 | 0.638 | 0.606 | 0.633 |
| oracle ceiling over the read-count curve | +0.022 | +0.023 | +0.007 | +0.010 |

Gorilla "moved intron annotated" is with the second annotation of row 1738 (RefSeq alone: 0.473 / 0.495). The weaker
chain of a shifted-site pair is annotated about as often as the average emitted chain on human and somewhat less
often on gorilla; about 0.6 of the moved introns are annotated. Shifts that keep the frame are in RefSeq about twice
as often as the others (0.386 vs 0.157 OR6737, 0.242 vs 0.137 A119b): a graded signal, not a cut. No GTF-level
stratum (shift size, frame, reads ratio, splice motif) holds a class that is rarely annotated on both discovery
samples.

## Verdicts, as registered

| row | question | verdict | figures |
|---|---|---|---|
| 1725 | class and ceiling; existing `RUSTLE_POLISH_JUNCTION_SNAP=reads` genome-wide | FAIL | ceiling STOP on A119b (+0.007 vs +0.02); snap changes 637 / 385 chains, near harmless and inert (reach 0.22 / 0.085 vs 0.25) |
| 1726 | rule A derived on OR6737, confirmed on A119b, frozen, HG002 held out | FAIL held out | dropped chains in RefSeq 0.060 vs the read cut's 0.062 on HG002: adds nothing over a read cut. With the corrected comparator A119b would have missed the derivation bar (0.56 x vs 0.5 x) |
| 1727 | rule A on KB3781 | PASS, gorilla libraries only | 549 chains at 0.020 (0.13 x comparator); any shift 1,611 at 0.038, precision 0.3584 -> 0.3658 (read cut 0.3631). Not a blind test (row 1677 had shown the direction) |
| 1728 | read-level features at the junction, A119b step 0 (10 cells) | GO, one cell | rule B: 578 chains, 0.033 in RefSeq, moved intron annotated 0.100 (class 0.606). Qualified only on the rounded bar (19 of 578 = 0.03287 vs 0.03285). Homopolymer and intron read ratio: no signal |
| 1729 | rule B frozen, HG002 held out | FAIL by two chains | 473 chains, 0.036 in RefSeq (17; the bar allows 15), moved intron annotated 0.099 (class 0.633) |
| 1730 | rule B frozen on gorilla | PASS both | OR6737 655 chains, 0.026 in RefSeq, 0.040 intron annotated (RefSeq); KB3781 308, 0.023, 0.042 |
| 1732 | rules A + B applied to the gorilla assemblies | FAIL on the pair bar by two pairs per sample | precision 0.3517 -> 0.3561 (read cut 0.3542), 0.3584 -> 0.3619 (0.3606); pairs int3 <= 10 bp 101 -> 99 of 103, 174 -> 172 of 187 |
| 1734 | port of rules A + B as `shiftdrop`; families / assign with the filtered input | PASS (no harm, no gain) | parity exact; assignments identical to base on OR6737 and KB3781; copies tables identical on KB3781, one field of one row on OR6737 |
| 1738 | second gorilla annotation (T2T primates CAT consensus, same assembly) | FAIL for rule A on KB3781; rule B PASS both | moved intron annotated: rule B 0.072 / 0.097; rule A 0.120 / 0.202 (bar 0.15). CAT takes RefSeq as an input: wider, not independent |
| 1739 | `shiftdrop` cut to rule B, four samples | PASS all five bars | table below |

### Row 1739: the shipped stage

| | OR6737 | KB3781 | A119b | HG002 on CHM13 |
|---|---|---|---|---|
| chains dropped | 655 | 308 | 578 | 473 |
| dropped chains in RefSeq; moved intron annotated | 0.026; 0.072 | 0.023; 0.097 | 0.033; 0.100 | 0.036; 0.099 |
| correct chains lost | 0.07% | 0.03% | 0.05% | 0.05% |
| row 1714 internal pair strata, largest fall (bar 0.01) | 0.0097 | 0.001 | 0.004 | 0.006 |
| precision before -> after | 0.3517 -> 0.3546 | 0.3584 -> 0.3599 | 0.1764 -> 0.1767 | 0.2336 -> 0.2342 |
| gain over a read cut of the same size | +0.0013 | +0.0006 | +0.0001 | +0.0001 |

Row 1739 is a no-harm confirmation, not a held-out validation: the RefSeq / intron bar and the human precision
figures were seen before it was registered. New were the annotated pairs on OR6737, KB3781 and HG002 and the gorilla
precision against the read cut. OR6737 passes the pair bar by one pair. Rows 1728 / 1729 stay FAIL against their own,
stricter comparator bar. On the human samples the stage is equal to a read cut of the same size to three decimals.

## What the rows establish

- A weaker chain that differs from a stronger emitted chain by one moved internal exon end is not a rarely
  annotated class: 0.19-0.29 of it is in RefSeq and about 0.6 of the moved introns are annotated, so dropping the
  class removes annotated chains.
- Two features mark a rarely annotated part of it: a frame-breaking shift of <= 10 bp at under 0.05 of the partner's
  reads (gorilla libraries only), and indels beside the moved intron in most of its reads (all four samples; holds
  within read strata: at 5+ junction reads the moved intron is annotated in 0.09 vs 0.66-0.67 of the rest on human).
- "Rarely annotated" is all that is shown. Nothing shows the removed chains are alignment errors.
- The gorilla / human difference of rule A sits mostly in the comparator: 2-read chains are in RefSeq 0.170 / 0.154
  on gorilla and 0.066 on human, while the rule A cell is 0.021 / 0.020 on gorilla and 0.037 / 0.060 on human.
- The existing switch `--polish-junction-snap reads` was measured genome-wide for the first time: close to harmless,
  changes 0.9% / 0.2% of chains, and already removes 206 of the 578 A119b chains of rule B.

## Not done, not known

- Families and assign were not re-run for the rule B stage. Row 1734 measured the two-rule version on gorilla only;
  for the shipped stage the fact is that no other stage reads its output.
- Rule A is a script (`shift_rule.py`), gorilla-only, and fails the wider annotation on KB3781. The union of rules A
  and B was never run on human.
- The A119b assembly is the one of rows 1714-1723 (`r1617/default/a.gtf`), not today's default. HG002 and KB3781 are
  spent for this class; a further rule needs a new sample or a new yardstick.
- The 0.36-0.39 of moved introns that no annotation holds have no yardstick here.
- Whether `shiftdrop` ever joins `all` is the user's decision; the measured gain does not argue for it.
- Two `all` dispatch tests in `tests/pipeline_driver_dispatch.rs` failed when last run because `capstarts` now runs
  inside `all` (another session's change); the `shiftdrop` test and 18 others pass.

## Scripts and files

| script (`bench/famcmp/`) | use |
|---|---|
| `shiftdrop.py GTF BAM OUT_GTF [THREADS]` | the stage (rule B, annotation-free); `--selftest` |
| `shift_class.py REF_GTF GTF FASTA [EXTRA_GFF3 \| -] [OUT_TSV]` | class table per emitted chain, strata, oracle ceiling; `--selftest` |
| `shift_rule.py` | rule A grid (24 cells) and scoring against the read cut; can write a filtered GTF |
| `shift_reads.py CLASS_TSV BAM FASTA OUT_TSV [THREADS]` | read-level features of the moved intron (10 cells) |
| `shift_filter.py CLASS_TSV READS_TSV IN_GTF OUT_GTF union\|indel` | scored filter the stage is checked against |

Workspace `/anvme/workspace/v138da13-rustle`: `r1725_shift/` (class tables, snap arm), `r1726_shiftrule/` (rule A),
`r1728_shiftreads/` (read features, filtered GTFs, pair tables), `r1733_shiftdrop/` (two-rule stage, families /
assign arms), `r1735_ggo_cat/` (second gorilla annotation: `cat.gff3`, built from
`mGorGor1_v2.0_CatRefConsensus_v3.bb`, GenomeArk), `r1739_shiftdrop_indel/` (shipped stage on four samples).

Inputs: OR6737 `r1628/orrun/or.gtf`, KB3781 `r1622/kb/kb.gtf`, A119b `r1617/default/a.gtf`, HG002 on CHM13
`r1692_trunc5_hg/hg.gtf`. Gorilla and human samples were never pooled.
