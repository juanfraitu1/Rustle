# PREREG row 1717 (2026-11-14): the cap G of Iso-Seq reads as evidence for a transcript start

User goal 2026-11-14: "find empirical rules to find lost nested isoforms and overall improve alternative TES/TSS as
well as alternative 5'/3' acceptors for internal exons". Rows 1715 / 1716: read-end positions cannot tell a nested
start from a 5' truncation, and the not-emitted nested chains are mostly not in a CAGE peak. New evidence tried here:
a full-length cDNA of a capped mRNA carries a non-templated G at its 5' end, visible as a 1-3 base all-G soft clip.

## Seen before this was written (disclosed)
- A119b chr19, read level, against the refTSS peaks of row 1716: all-G clip on 0.72 of reads starting in a peak, 0.28
  elsewhere (spliced), 0.05 (unspliced).
- Gate `capg.py gate` (RefSeq only; reads starting within 10 bp of a RefSeq 5' end vs inside an internal RefSeq exon
  >= 500 bp from any; reads whose genome base just upstream is G are uninformative and left out), capped share:
  A119b chr19+chr22 0.814 vs 0.044; HG002 chr19+chr22 0.819 vs 0.050; OR6737 NC_073224.2 0.059 vs 0.013 (read as: Gs trimmed); KB3781 NC_073224.2 0.808 vs 0.704 (read as: a G on nearly every read).
- No table of the signal by candidate class (annotated / emitted) on any sample. Row 1716's A-rich and polyA-motif
  shares of early ends on A119b (arm C uses them).
**Consequence fixed now:** the signal exists in the two human libraries only. The rule is not testable on the gorilla
samples; that is a stated limit, and any port must carry the gate as a self-check.

## Definitions (`bench/famcmp/capg.py rule`)
- Read states: `cap` (5' soft clip of 1-3 bases, all G), `uninf` (genome base just upstream is G), `no`.
- CAPPED CLUSTER of a chain: >= 3 capped reads of that exact chain starting within 10 bp of the commonest capped
  start, and capped >= 0.5 of the informative reads of the chain starting in that window. Found greedily, largest first.
- Candidates as row 1715: read chains with >= 5 primary reads that are a terminal sub-run of a chain emitted in
  `PREFIX.gtf`; NEST5 / NEST3 by the container's strand.
- Yardsticks: RefSeq chain (annotated), and the refTSS CAGE peaks lifted to CHM13 (row 1716) with the same test 100 bp
  downstream as background. Thresholds are fixed here, a priori (not fitted).
- Samples: A119b without chr19, chr21 and chr22 (chr19 + chr22 seen in the gate, chr21 + chr22 in the script trial, whose
  table was seen after this was written); HG002 on CHM13
  (`r1692_trunc5_hg/hg.gtf`) all chromosomes, as the held-out sample. Verdicts are read on HG002; A119b is reported.

## Arms and bars
**A1 rescue (lost nested starts):** not-emitted NEST5 candidates that have a capped cluster. PASS: >= 100 rescued,
their cluster position in a CAGE peak >= 0.40 and >= 5 x background. Reported: annotated-lost recovered, annotated
share of the rescued.
**A2 drop (unsupported nested starts):** emitted NEST5 candidates with >= 5 informative reads and no capped cluster.
PASS: >= 0.90 of the annotated emitted ones kept, >= 30% of the unannotated emitted ones dropped, and the dropped
unannotated ones in a CAGE peak (commonest start) <= 2 x background. Reported next to `trunc5` on the same table
(annotated kept, unannotated dropped); rescued and dropped are scored separately, never net.
**B same-chain starts:** emitted chains with >= 5 reads. B1: CAGE-peak share of the emitted 5' end, of the commonest
read start, and of the largest capped cluster (chains that have one); PASS for the cluster as the model's start when
it is >= 0.05 above the emitted end on annotated chains. B2: further capped clusters >= 50 bp from the largest; PASS
as reportable alternative starts when >= 100 of them, in a CAGE peak >= 0.40 and >= 5 x background.
**C early ends, drop side only** (the keep side failed in row 1715; the cap G says nothing about 3' ends): emitted
NEST3 candidates whose commonest end is A-rich with no polyA signal (`alt_ends.flags`). PASS: >= 0.90 of the
annotated emitted ones kept and annotated <= 0.05 of the dropped.
**D shifted ends of internal exons:** no arm. Row 1714 has them at the cassette control in every stratum (minor 2-4
reads 0.79-0.92 vs 0.80 on OR6737); what is lost there is the general low-read loss that rows 1616-1618 closed.

Reading: a PASS is evidence for an opt-in rule on libraries that pass the gate, to be ported only on the user's go.
A FAIL closes the cap G for that use. Nothing changes a default.

## Outcome (2026-11-14, job 4471387; `r1717_capg/{a,hg}.{gate.txt,txt,tsv}`)
Gate, whole genome: capped share at RefSeq 5' ends vs inside internal exons, spliced reads: A119b 0.823 vs 0.043;
HG002 0.845 vs 0.050. Verdicts are read on HG002 (held out); A119b (without chr19 / 21 / 22) in brackets.

| arm | result on HG002 (A119b) | verdict |
|---|---|---|
| A1 rescue of not-emitted nested starts | 2,568 of 72,129 have a capped cluster (1,185 of 72,857); cluster in a CAGE peak 0.710 vs background 0.043 (0.450 vs 0.044); 177 of 432 annotated-lost recovered (106 of 392); the rest, commonest start in a peak 0.026 (0.023) | **PASS** |
| A2 drop of emitted nested starts without a cluster | annotated kept 0.871 (0.818), bar 0.90; unannotated dropped 0.716 (0.792); dropped ones in a peak 0.032 vs background 0.014 = 2.3 x (0.021 vs 0.009), bar 2 x | **FAIL** on both bars |
| B1 start of a chain | annotated chains with a cluster (0.901 of them): emitted 5' end in a peak 0.412 (0.288), commonest read start 0.919 (0.788), largest capped cluster 0.926 (0.802) | **PASS** as registered; see reading |
| B2 further starts of a chain | 17,016 clusters on 10,317 annotated chains, in a peak 0.503 vs background 0.217 = 2.3 x (0.402 vs 0.224); unannotated chains 0.556 vs 0.180 | **FAIL** on the 5 x bar |
| C drop of A-rich early ends without a polyA signal | annotated kept 0.748 (0.830), bar 0.90; annotated among the dropped 0.095 (0.067), bar 0.05 | **FAIL** |

Reading:
- **A1 is the rule the goal asked for.** A not-emitted nested start with >= 3 capped reads at one position is CAGE-
  supported 16 x over background on the held-out sample: 0.710 (0.450), against 0.856 (0.729) for annotated nested
  starts that are emitted and kept by the same rule; without the cluster it is at 0.026. It recovers 41% (27%) of the annotated nested starts the ISM rule loses
  and adds about 2,400 (1,100) unannotated ones, 70% (43%) of them CAGE-supported. RefSeq alone would have scored
  this rule as 0.07-0.09 "correct".
- **A2:** against `trunc5` on the same table the cap G keeps more annotated nested starts (0.871 vs 0.729; 0.818 vs
  0.719) and drops fewer unannotated ones (0.716 vs 0.856); what `trunc5` drops is 0.101 in a peak, what the cap G
  drops 0.032. Better targeted, still short of the bar. The unannotated emitted ones that have a cluster are 0.635
  (0.358) in a peak.
- **B1:** the gain is over the emitted end, which is the most extended read end and so mostly lies upstream of the
  narrow peak. The commonest read start does as well as the capped cluster (0.919 vs 0.926): the lever is "commonest
  start, not most extended", and the cap G adds 0.007. This is a peak-membership test, not a distance to the
  annotated start; `--polish-tss-end median` (row 1487) is the existing switch in that direction and was not
  compared here.
- **B2:** further capped clusters are enriched (0.50 vs 0.22) but the 100 bp background is high next to a real
  start, so the 5 x bar was not reachable as written. Not a reportable set as registered.
- **C:** the drop side fails. The keep side of early ends (row 1715) was scored on RefSeq share, which this row
  shows understates real ends about tenfold for starts (rescued 0.07 annotated, 0.71 CAGE-supported): it has no
  independent yardstick here (no polyA-site atlas), so it is unmeasured, not closed.
- **Limits:** human libraries only (gate on the gorilla samples 0.059 vs 0.013 and 0.808 vs 0.704; why is not established);
  cap G and CAGE both mark a capped 5' end, so recapped cleavage products would pass both, and A1 was not split by
  whether the cluster lies inside the container's exon; overlap with the default `--polish-tss rescue` unmeasured; CAGE peaks are
  a multi-tissue aggregate; A1 was not run through the assembler, so its effect on the precision curve is unmeasured.
Later the same day: A1 was ported as the opt-in stage `capstarts` and measured on its output (row 1720, PASS,
`docs/PREREG_alt_rules_2026-11-14.md`); the early-end keep side got a yardstick and a rule (rows 1721 / 1722, PASS).

# PREREG row 1718 (2026-11-14): does the existing median 5' end reach the commonest-start level?
Written after row 1717's tables (B1: commonest read start in a CAGE peak 0.919 HG002 / 0.788 A119b on annotated chains
with a capped cluster, emitted most-extended end 0.412 / 0.288). `PREFIX.final.gtf` already carries
`--polish-tss-end median` (row 1487; `tools/rustle_pipeline.sh:713`), so no new assembly is run.
Measured (`capg.py ends`): annotated chains with `reads` >= 5 present in both `PREFIX.gtf` and `PREFIX.final.gtf`
(HG002 `r1692_trunc5_hg/hg`, A119b `r1617/default/a` without chr19 / 21 / 22): 5' end in a CAGE peak and within 50 bp
of one, both files, with the 100 bp downstream background; next to the commonest read start on all annotated chains
>= 5 reads from row 1717's table (in peak only).
Reading (HG002; A119b reported): median end within 0.05 of the commonest start's in-peak share = "the existing switch
is the fix, nothing to port"; further below = a commonest-start (mode) rule is a candidate port, user's call.
Measurement only; no default changes.

## Outcome of row 1718 (2026-11-14, login node; `r1717_capg/ends1718.txt`): the existing median end is the fix
Annotated chains with >= 5 reads, HG002 29,394 (A119b 25,731):

| 5' end | in a CAGE peak | within 50 bp | background in peak |
|---|---|---|---|
| most extended, `PREFIX.gtf` | 0.417 (0.288) | 0.781 (0.620) | 0.245 (0.218) |
| median, `PREFIX.final.gtf` | 0.864 (0.711) | 0.948 (0.863) | 0.073 (0.082) |
| commonest read start (row 1717 table) | 0.891 (0.737) | not computed | |

The median end is 0.027 (0.026) below the commonest start, inside the 0.05 bar: nothing to port. `PREFIX.final.gtf`
already carries it; in `PREFIX.gtf` it is the opt-in `RUSTLE_POLISH_TSS_END=median`, proposed as default in row 1487
and still the user's decision. Limit: the commonest-start figure comes from a slightly different chain set (all
annotated chains with >= 5 primary reads) and has no within-50 bp column.

## Why the gorilla libraries fail the cap gate (exploratory, 2026-11-14, user question; no decision rule, not a register row)

Count of non-templated 5' G (all-G soft clip) on reads whose upstream genome base is not G, first three annotated contigs:

| library | reads at a RefSeq 5' end: 0 / 1 / 2 G | reads starting inside an internal exon: 0 / 1 / 2 G |
|---|---|---|
| HG002 | 0.136 / 0.818 / 0.020 | 0.903 / 0.064 / 0.002 |
| OR6737 (`GGO_mm.bam`) | 0.886 / 0.055 / 0.005 | 0.897 / 0.015 / 0.002 |
| KB3781 (`GCA_029281585.2_flnc_mm.bam`, SRR27178663) | 0.148 / 0.364 / 0.423 | 0.235 / 0.534 / 0.181 |

- OR6737: the G is not in the reads at all (89% of reads at annotated starts have no extra G; 44% of all reads start with G,
  the base rate). Whatever produced that FASTQ removed it; nothing in the BAM can restore it.
- KB3781: nearly every read has at least one extra G whether or not it starts at an annotated 5' end (77% of
  internal-exon starts), and capped starts carry one more only part of the time. A shifted rule (cap = 2 or more G) gives
  0.455 at annotated starts vs 0.194 inside internal exons, below both gate bars (0.5 / 0.15).

So the two failures have different causes and neither is fixable by a threshold on these files. The wet-lab or trimming
reason behind each is not established (it is not recorded in the BAM headers).

## User ruling 2026-11-14: `capstarts` is a default side product

"ok, lets make it default then": `all` runs `capstarts` as a last step and writes `PREFIX.capstarts.gtf` (`RUSTLE_CAPSTARTS=0` = off), like `earlyends`. `PREFIX.gtf` and every stage that reads it are unchanged. On a library without usable cap G the stage refuses, which `all` reports and does not retry while `PREFIX.gtf` is unchanged. Evidence is two human samples (rows 1717, 1720). Trial on A119b chr22: 25 starts added, +16 s at 8 threads; switch, refusal and no-retry paths exercised.
