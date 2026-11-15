# PREREG row 1714 (2026-11-14): are alternative TSS / TES and shifted ends of internal exons modelled?

User question 2026-11-14, then "ok lets do that". Measurement, no code change. Written before any pair-level count.
Seen before writing: the class aggregates of rows 1618 / 1677 (junction "wobble", 5' truncation class) and the
end-accuracy table of the 2026-11-12 ends study. No default is proposed from this row.

## Truth and events (`bench/famcmp/alt_events.py`)
RefSeq (`GGO.GCF_029281585.2_RefSeq.gtf.gz`, `A119b.chm13v2.0_RefSeq.gtf.gz`), unique multi-exon intron chains per
contig and strand. One row per unordered pair of annotated chains on one contig and strand:
- **SITE**: equal chains except one intron, which shares one boundary. Sub-classes by the exon end that moves:
  `int3` (donor, the exon is internal), `int5` (acceptor, the exon is internal), `first3` (donor of the first exon,
  first exons overlap), `last5` (acceptor of the last exon, last exons overlap). Shift strata <= 10 bp and > 10 bp.
- **AFE / ALE**: as `first3` / `last5` but the two first (last) exons do not overlap (most extended annotated ends).
- **NEST5 / NEST3**: one chain is the other's contiguous run reaching its 3' (5') end: an internal alternative start /
  an early alternative end. This is where a real alternative start meets the 5' truncation class.
- **CASSETTE (control)**: equal chains except two consecutive introns of one are one intron of the other.
- **SAME3 / SAME5**: one annotated chain with two annotated 3' (5') ends >= 100 bp apart.

## Expressed
Full-chain reads = primary spliced records whose intron chain equals the annotated chain. A pair is EXPRESSED when both
chains have >= 2 (the floor of row 1616's default; `RUSTLE_MIN_CHAIN_READS=1` is the opt-in below it). Strata on the
minor chain: 2-4 and >= 5 reads; minor share of the pair's reads < 0.1 and >= 0.1. SAME3 / SAME5: EXPRESSED when each
annotated end has >= 3 full-chain read ends within 25 bp.

## Scored
Per sample, class and stratum: both chains emitted / major only / minor only / neither, for `PREFIX.gtf` (default),
`PREFIX.final.gtf`, `PREFIX.trunc5.gtf`, StringTie and FLAIR (same BAMs; contig + chain match, tool strand ignored).
Samples OR6737 and A119b; KB3781 ours only. SAME3: the default emits one model per chain (two is impossible by
construction); scored = the emitted 3' end is within 25 bp of one annotated end and the opt-in `altends` stage reports
`alt_tes` within 25 bp of the other. SAME5: counted only, nothing reports alternative 5' ends.

## Reading rule (relative; fixed here)
At minor >= 5 reads, on OR6737 and A119b both: a class is "modelled like ordinary alternative splicing" when its
both-emitted share in `PREFIX.gtf` is within 0.05 of CASSETTE's, and "behind the tools" when more than 0.05 below
StringTie's or FLAIR's on the same pairs. Classes with < 50 expressed pairs in a sample are reported, not read.
Limits stated up front: RefSeq ends are weak truth (gorilla models are predicted, A119b 3' ends); a pair needs both
isoforms expressed in one tissue, so the counts are small; large genes contribute many pairs.

## Outcome (2026-11-14, job 4471351; tables `r1714_altevents/{or,a,kb}.txt`, pairs `*.pairs.tsv`)
Both-emitted share in `PREFIX.gtf` at minor >= 5 reads (expressed pairs), OR6737 / A119b / KB3781; StringTie and FLAIR
on the same pairs for OR6737 / A119b:

| class | OR6737 | A119b | KB3781 | StringTie | FLAIR |
|---|---|---|---|---|---|
| CASSETTE (control) | 0.975 (3,531) | 0.912 (7,382) | 0.942 | 0.799 / 0.695 | 0.916 / 0.876 |
| internal exon 3' end, <= 10 bp | 0.981 (103) | 0.902 (225) | 0.930 | 0.049 / 0.044 | 0.893 / 0.889 |
| internal exon 3' end, > 10 bp | 0.966 (469) | 0.913 (908) | 0.938 | 0.471 / 0.423 | 0.893 / 0.874 |
| internal exon 5' end, <= 10 bp | 0.974 (783) | 0.912 (1,248) | 0.956 | 0.033 / 0.031 | 0.931 / 0.897 |
| internal exon 5' end, > 10 bp | 0.954 (690) | 0.909 (1,333) | 0.934 | 0.443 / 0.421 | 0.916 / 0.897 |
| first exon 3' end (both strata) | 0.983 / 0.979 | 0.959 / 0.931 | 0.966 / 0.956 | 0.051-0.670 | 0.818-0.888 |
| last exon 5' end (both strata) | 1.000 / 0.978 | 0.938 / 0.907 | 0.980 / 0.930 | 0.040-0.688 | 0.780-0.828 |
| AFE | 0.942 (1,874) | 0.902 (4,294) | 0.901 | 0.764 / 0.586 | 0.853 / 0.796 |
| ALE | 0.949 (138) | 0.903 (401) | 0.892 | 0.833 / 0.589 | 0.848 / 0.736 |
| NEST5 (internal alternative start) | 0.433 (215) | 0.572 (1,377) | 0.518 | 0.381 / 0.333 | 0.242 / 0.259 |
| NEST3 (early alternative end) | 0.671 (289) | 0.458 (1,506) | 0.561 | 0.567 / 0.440 | 0.412 / 0.268 |

Reading by the fixed rule:
- **Shifted ends of internal exons, shifted ends of terminal exons, AFE, ALE: modelled like ordinary alternative
  splicing** (every class within 0.05 of CASSETTE on both samples) and behind neither tool. StringTie emits both
  chains of a <= 10 bp shift in 3-5% of pairs.
- **NEST5 and NEST3: not modelled like ordinary splicing** (0.43-0.67 against 0.91-0.98); not behind the tools, which
  are lower. `PREFIX.final.gtf` costs NEST5 0.02-0.04; the opt-in `trunc5` stage halves it (0.274 / 0.248): the 5'
  truncation drop removes annotated internal starts.
- **Same chain, two ends: not measurable with RefSeq.** Annotated chains with two 3' ends >= 100 bp apart: 4 (gorilla),
  9 (A119b), 1 expressed. The default cannot emit two by construction; `altends` was run on both samples (252 s, 620 s)
  but has nothing to be scored against. Nothing reports alternative 5' ends.
At all expressed pairs (minor >= 2) the default shares are lower (CASSETTE 0.896 / 0.780; SITE 0.78-0.97; NEST 0.38-0.55):
most single-emitted pairs keep the major chain. Not separated: in NEST pairs, whether the surviving chain is the long
or the short one; pairs are not independent (large genes give many).
