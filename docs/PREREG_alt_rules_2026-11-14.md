# PREREG rows 1719-1721 (2026-11-14): rules for lost alternative isoforms, applied and measured through the output

User goal 2026-11-14 (as row 1717). The stop check on that goal: row 1717's rule was not applied, early ends have no
yardstick, and nothing was found for shifted ends of internal exons. Three rows, each written before its tables.

## Seen before this was written
- Rows 1714-1718 and their tables.
- Post hoc on row 1714's pairs (minor >= 5 reads, only the major chain emitted; scratch script, A119b / OR6737): the
  lost chain of non-nested pairs is under 0.02 of its locus's top emitted chain in 140 of 365 / 26 of 62 (shifted
  site), 253 of 448 / 39 of 58 (cassette), 87 of 220 / 32 of 65 (alternative first exon); a sub-chain of an emitted
  chain in 45 / 13, 22 / 4, 65 / 25; within 10 bp of an emitted chain 69 / 7 (shifted site); unexplained 111 / 16,
  160 / 15, 61 / 8. So the largest single cause is the isoform-fraction rule (`--polish-isoform-fraction 0.02`).
- Row 1616: exempting chains from that rule by read count alone (3 / 5 / 10 reads) gains chains at 0.07 annotated on
  OR6737 (FAIL). Row 1614: the rule off gains 947 correct chains at precision 0.352 -> 0.301 on OR6737.
- GENCODE V44 polyA track downloaded (50,520 `polyA_site` features, hg38); nothing scored against it.

## Row 1719: single-event sibling exemption from the isoform-fraction rule (shifted sites, cassette exons)
CANDIDATE = a chain in the fraction-off assembly (`--polish-isoform-fraction 0`, all else default) that is not in
the default `PREFIX.gtf`, with `reads` >= 5. SIBLING = it differs from a chain of the default `PREFIX.gtf` (same
strand) by exactly one event: one intron replaced by one that shares exactly one boundary with it (shifted donor or
acceptor, or an alternative first / last exon), or one intron replaced by two that keep its outer boundaries (a
cassette exon in or out). CANONICAL = every intron of the candidate that is not in the sibling is GT-AG on its strand.
Rule: keep candidates that are SIBLING and CANONICAL (`bench/famcmp/sibling_frac.py`).
Samples: OR6737 (`r1614/nofrac`, discovery; its totals were seen), A119b and HG002 on CHM13 (new fraction-off runs,
held out). Yardstick: RefSeq chain; reported also "every intron annotated".
PASS on a sample: >= 100 annotated chains kept, and the kept set's annotated share >= the default assembly's own
chain precision on that sample (the bar of row 1616) ; verdict = PASS on both held-out samples. Reported: non-sibling
candidates' share (the rule must beat it >= 2 x to count as a rule and not a read floor); the share of row 1714's
lost shifted-site / cassette pairs it restores (A119b, OR6737); sanity: default chains absent from the fraction-off run.

## Row 1720: the cap-G rescue as an opt-in stage, measured on its output
Port of row 1717 A1 exactly as registered there: driver stage `capstarts` (`capg.py stage`) writes
`PREFIX.capstarts.gtf` = `PREFIX.gtf` plus every not-emitted nested start (>= 5 reads, 3'-terminal sub-run of an
emitted chain) that has a capped cluster, with its 5' end at the cluster and its 3' end at the container's. It refuses
to write on a library that fails a self-check needing no annotation: capped share of informative spliced reads
starting in the first exon of an emitted model >= 0.5 and of those starting inside an internal exon <= 0.15.
Measured on HG002 and A119b (all chromosomes; the rule was frozen in row 1717): row 1714's nested-start pairs
(minor >= 5) with both chains emitted, default vs capstarts; exact-chain precision / recall against RefSeq
(`curve.py`); the self-check on all four samples.
PASS: the nested-start both-share rises by >= 0.10 on both samples, RefSeq chain precision falls by <= 0.005, and the
self-check accepts both human samples and refuses both gorilla ones. Opt-in whatever the verdict; no default changes.

## Row 1721: early ends against annotated polyA sites
Yardstick: GENCODE V44 `polyA_site` features lifted to CHM13 (`liftOver`, row 1716's chain); a model end is
SUPPORTED when a site on its strand lies within 25 bp; background = the same test 200 bp upstream in transcript
sense. Models: row 1716's early-end (NEST3) models on A119b (`r1716_sqanti/a.models.gtf`, 3' end = commonest read
end) and the same built for HG002. Seen: their RefSeq shares (row 1715) and SQANTI3 motif shares (row 1716).
Gate: annotated emitted early ends supported >= 3 x their background, else the yardstick is void and the row stops.
Rule (keep side): a not-emitted early end is kept when >= half of its reads end within 10 bp of the commonest end,
a polyA signal lies upstream and no A-rich run follows (`alt_ends.flags`), the signature row 1715 found most enriched.
PASS (HG002; A119b reported): >= 100 kept, supported share >= half of the annotated emitted early ends' share and
>= 3 x background. Reported: the complement (not-emitted without the signature) and the unannotated emitted ones.
A PASS is a candidate opt-in stage; a FAIL closes early ends with a yardstick this time.

## Outcome of row 1719 (2026-11-14, login node; `r1719_sibling/or.txt`): FAIL on the discovery sample, held-out runs not made
OR6737, candidates = 1,486 chains with >= 5 reads that only the fraction-off run emits (11 default chains are absent
from that run):

| set | candidates | annotated | every intron annotated |
|---|---|---|---|
| sibling, GT-AG (the rule) | 1,290 | 92 (0.071) | 0.102 |
| sibling, not GT-AG | 45 | 3 (0.067) | 0.067 |
| not a sibling | 151 | 9 (0.060) | 0.073 |

Bar: annotated share >= 0.352 (the default's chain precision) and >= 100 annotated. The rule keeps 0.071 and 92, and
does not separate (0.071 vs 0.060 for non-siblings): nearly every chain the fraction rule drops is a one-event sibling
of an emitted chain. It would restore 31 of 72 lost shifted-site pairs and 48 of 65 cassette pairs of row 1714, at
1,198 unannotated chains. **Stopped here:** a rule that fails by 5 x on the sample it was designed on is not run on
A119b / HG002 (no fraction-off assemblies made). Deviation from the registration, stated.
Post hoc breakdown of the siblings (same sample, not a test): shifted end of an internal exon 6 of 355 annotated
(0.017), and 0 of 202 at shifts <= 10 bp; shifted end of a terminal exon 49 of 419 (0.117); cassette exon skipped
18 of 564 (0.032); cassette exon included 62 of 558 (0.111), 0.170 when its length is a multiple of 3 against 0.072.
No stratum reaches 0.2. Reading: for shifted ends of internal exons the fraction rule is doing its job (what it
removes is 98% unannotated), and the default already emits both isoforms at the cassette-control rate (row 1714).
There is no rule to add for this class; the measured state is the result.

## Outcome of row 1720 (2026-11-14, job 4471506; `r1720_capstarts/`): PASS, stage `capstarts` built (opt-in)
| | A119b | HG002 |
|---|---|---|
| nested starts added | 1,298 | 2,568 |
| row 1714 nested-start pairs (minor >= 5), both emitted: default -> capstarts | 0.572 -> 0.702 (788 -> 967 of 1,377) | 0.521 -> 0.700 (802 -> 1,076 of 1,538) |
| same, all expressed pairs | 0.484 -> 0.571 | 0.441 -> 0.567 |
| RefSeq chain precision (chains / correct) | 0.176 -> 0.176 (220,815 / 38,833, +118 correct) | 0.234 -> 0.231 (155,601 / 35,922, +177 correct) |
| self-check, first exon / internal exon | 0.734 / 0.020 accepted | 0.807 / 0.021 accepted |

Self-check on the gorilla samples: OR6737 0.058 / 0.005 refused, KB3781 0.826 / 0.669 refused. All three bars met
(+0.130 and +0.179 vs +0.10; precision -0.000 and -0.003 vs -0.005; accept / refuse as required). Other classes
unchanged (cassette 0.912 -> 0.913, 0.889 -> 0.891). Limits: human libraries only; RefSeq counts the added chains as
0.07-0.09 correct while CAGE supports 0.45-0.71 of them (row 1717); families / assign never read this file.

## Outcome of row 1721 (2026-11-14, same job; `r1719_polya/{a,hg}.polya.txt`): PASS
Share of early-end models with a GENCODE polyA site within 25 bp (background 200 bp upstream), HG002 (A119b):

| set | with the signature | without it |
|---|---|---|
| emitted, annotated | 0.851 of 1,037 (0.802 of 1,326) | 0.358 (0.361) |
| not emitted, annotated | 0.772 of 294 (0.765 of 281) | 0.238 (0.227) |
| emitted, unannotated | 0.563 of 1,180 (0.460 of 2,018) | 0.044 (0.055) |
| not emitted, unannotated | 0.586 of 3,777 (0.493 of 2,649) | 0.038 of 27,200 (0.035 of 18,620) |

Gate: annotated emitted early ends 0.591 vs background 0.011 (0.594 vs 0.008): the yardstick works. Rule: 4,071
(2,930) not-emitted early ends carry the signature, 0.599 (0.519) of them at an annotated polyA site, against a bar of
0.296 (0.297) and a background of 0.011: PASS. The signature separates 15 x among the unannotated not-emitted ones.
This reverses row 1715's reading for early ends: the same signature scored 0.08-0.10 "annotated" by RefSeq chain,
because RefSeq lacks the chain, not because the end is false. Limits: GENCODE sites are annotation too (human
curated, partly from the same cDNA evidence as RefSeq), and its curators use the polyA signal, so signature and
yardstick are not fully independent; human only. A first run of the scorer on A119b printed zeros (a counting bug in
the script, fixed before any figure was read).

## Row 1722 (pre-registered 2026-11-14, after rows 1720 / 1721): the early-end rule as an opt-in stage, on its output
Port of row 1721's rule: driver stage `earlyends` (`capg.py stage --ends`) writes `PREFIX.earlyends.gtf` =
`PREFIX.gtf` plus every not-emitted early end (>= 5 reads, 5'-terminal sub-run of an emitted chain) with the
signature; 3' end at the commonest read end, 5' end at the container's. No library check (it reads no cap).
Measured on HG002, A119b, OR6737, KB3781: row 1714's early-end pairs (minor >= 5) with both chains emitted, default
vs earlyends; RefSeq chain precision. Seen: the stage on A119b chr21 + chr22 adds 94 models (count only).
PASS: the early-end both-share rises by >= 0.10 and precision falls by <= 0.005 on HG002 and A119b; the gorilla
samples are reported as transfer (no polyA-site yardstick there). Opt-in whatever the verdict.

## Outcome of row 1722 (2026-11-14, job 4471530; `r1722_earlyends/`): PASS, stage `earlyends` built (opt-in)
| | A119b | HG002 | OR6737 (transfer) | KB3781 (transfer) |
|---|---|---|---|---|
| early ends added | 2,930 | 4,071 | 842 | 1,643 |
| row 1714 early-end pairs (minor >= 5), both emitted | 0.458 -> 0.669 (690 -> 1,008 of 1,506) | 0.434 -> 0.630 (917 -> 1,331 of 2,112) | 0.671 -> 0.865 (194 -> 250 of 289) | 0.561 -> 0.745 (264 -> 351 of 471) |
| RefSeq chain precision | 0.1764 -> 0.1753 (+281 correct) | 0.2336 -> 0.2294 (+294 correct) | 0.352 -> 0.349 (+68) | 0.358 -> 0.352 (+102) |

Bars met on both human samples (+0.211 and +0.196 vs +0.10; precision -0.001 and -0.004 vs -0.005). The gorilla
samples move the same way (+0.19, +0.18) with no polyA-site yardstick there; KB3781 loses 0.006 of RefSeq precision.
Other classes unchanged. Limits: RefSeq counts the added chains as 0.07-0.10 correct while 0.52-0.60 sit at an
annotated polyA site (row 1721, human); `capstarts` and `earlyends` write separate files and were not combined;
families / assign never read either file; nothing was run on a fifth sample.

**Both stages together** (`r1722_earlyends/{a,hg}.both.gtf`, `curve.py`, measured after the two verdicts): A119b
223,730 chains, 39,112 correct, 0.175 (-0.002 vs default); HG002 159,660 chains, 36,216 correct, 0.227 (-0.007).
Each stage met the -0.005 precision bar alone; on HG002 the two together do not. RefSeq chain precision is the
pessimistic reading (rows 1717 / 1721), but that is the bar that was registered.
**Row 1719, limits added after review:** it tested only the fraction-rule class and was scored on RefSeq chain
identity alone, the yardstick that hid real ends in row 1715. `bench/famcmp/lost_why.py` (the breakdown cited under
"Seen", now with membership in the filter-off runs of rows 1614 / 1616) on OR6737's 62 lost shifted-site chains:
26 fraction class (all back with the fraction rule off); 13 sub-chains (7 back with the ISM rule off); 16 other
(6 back with the fraction rule off, 10 with the polish off, 6 in no run); 7 within 10 bp of an emitted chain, of
which 1 is back with the polish off, so that class is lost before the polish.
Attribution of the pre-polish class (code reading by a delegated agent, then checked on the pairs; no new row):
`assemble` runs `--assembly-junctions strict` (`tools/rustle_pipeline.sh:752`), which drops a chain when any of its
introns is not GT-AG / GC-AG / AT-AC on its strand (`src/family.rs:51622-51631`), whatever its read count. Check on
row 1714's non-nested minor chains with >= 5 reads: emitted ones with a non-canonical intron 0 of 13,132 (A119b) and
0 of 7,095 (OR6737); lost ones 62 of 904 and 12 of 162 are non-canonical. So about 7% of the lost annotated minor
isoforms are RefSeq chains with a non-canonical intron; the existing switch for them is `RUSTLE_JUNCTION_MAJORITY=1`
(strict under `--assemble-only` was a decision of row 1074; not re-tested here). The rest of the lost canonical
chains fall to the fraction rule, the ISM rule and other polish rules as counted above.
**Scope of the "no rule" statement for shifted ends of internal exons:** it covers the fraction class (tested, FAIL)
and the non-canonical class (existing switch, a past decision); the ISM-class and the 10 "other" chains per 62 were
not given a rule. The default emits both isoforms of 0.90-0.98 of such pairs, at the cassette control.

## User rulings 2026-11-14 (question dialog)
- `capstarts` and `earlyends`: "Keep opt-in". Neither is folded into a default output.
- Median 5' end for `PREFIX.gtf`: "Leave as is". `PREFIX.final.gtf` is the file to quote for start accuracy.
- Shifted ends of internal exons, non-canonical class: "Measure the switch" -> row 1723.

## Row 1723 (pre-registered 2026-11-14, user go): majority junctions on the non-canonical class
Seen: 62 of 904 (A119b) and 12 of 162 (OR6737) lost non-nested minor chains with >= 5 reads have a non-canonical
intron; no emitted chain has one; rows 959 / 960 / 1074 (strict under `--assemble-only` was decided there).
Run: `assemble` with `RUSTLE_JUNCTION_MAJORITY=1`, all else default, on OR6737 and A119b (current binary).
Measured: chains only in the majority run (new) and their annotated share; default chains missing from it; how many
of the lost non-canonical minor chains it emits (`lost_why.py`); row 1714's shifted-site pairs with both emitted.
PASS on a sample: >= half of the lost non-canonical minor chains emitted, and the new chains' annotated share >= the
default's chain precision (0.352 OR6737, 0.176 A119b; the bar of row 1616). Verdict needs both samples. A PASS is a
proposal to the user (the switch exists); a FAIL leaves the strict setting.

## Outcome of row 1723 (2026-11-14, job 4471745; `r1723_majority/`): FAIL on the precision bar, strict stays
| | OR6737 | A119b |
|---|---|---|
| lost non-canonical minor chains emitted by the majority run | 19 of 19 | 98 of 106 |
| lost canonical ones emitted | 0 of 173 | 0 of 955 |
| new chains, annotated share (bar: default precision) | 3,011 at 0.021 (bar 0.352) | 14,342 at 0.014 (bar 0.176) |
| default chains missing from the majority run | 263 | 1,444 |
| shifted internal-exon ends, both emitted (4 strata, minor >= 5) | 0.954-0.981 -> 0.961-0.982 | 0.902-0.913 -> 0.916-0.937 |

Bar 1 is met (the switch does recover that class), bar 2 fails by more than 10 x: 63 and 194 annotated chains come
with 2,948 and 14,148 unannotated ones, and 263 / 1,444 default chains are lost. The counts of lost non-canonical
chains (19 / 106) are larger than the 12 / 62 under "Seen" because this count takes every lost minor chain, not only
pairs whose major chain is emitted. Reading: the strict junction setting costs 1-2 points of shifted-site pairs and
keeps out about 3,000-14,000 chains that are 98% unannotated. No rule for shifted ends of internal exons survives a
test; the class stays at the cassette-control level. Limit: scored on RefSeq chain identity only; a splice site has
no CAGE- or polyA-like independent yardstick here.
