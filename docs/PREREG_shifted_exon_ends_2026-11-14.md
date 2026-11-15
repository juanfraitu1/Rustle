# PREREG row 1725 (2026-11-14): shifted ends of internal exons, the emitted side

User goal: "derive empirical rules to improve assembly focusing on shifted internal exons". Written before any table
of this row.

Seen before writing: rows 1714 / 1719 / 1723 (the default emits both isoforms of 0.90-0.98 of annotated shifted-site
pairs; no rule restores the lost ones); row 1618 (read-level 10 bp snap loses 1,598 correct chains); row 1677's
`wobble` class (an intron within 10 bp of one with >= 5 x its reads: share correct 0.70 / 0.99 / 0.77 of the rest on
OR6737 / A119b / KB3781, not a class on A119b); the help text of `--polish-junction-snap reads` (existing opt-in,
designed on three dev contigs, never run genome-wide, no register row). Not seen: how many emitted chains are the
weaker side of a shifted-site pair, and how many of those are annotated.

Question: the recall side is closed. Is there a precision-side rule: emitted chains that are the weaker copy of a
stronger emitted chain with one exon end moved?

Samples: discovery OR6737 and A119b (never pooled), the assemblies of rows 1714-1723 (`r1628/orrun/or.gtf`,
`r1617/default/a.gtf`). Held out: HG002 on CHM13 (`r1692_trunc5_hg/hg.gtf`), opened only if step 1 passes on both
discovery samples and after A119b is recorded. KB3781 is spent for chain rules.

## Step 0: the class and its ceiling (GTF + annotation only, `bench/famcmp/shift_class.py`)
- **Class S1**: an emitted chain c with an emitted chain d of the same contig and strand, reads(d) >= reads(c), equal
  to c except one intron that shares exactly one boundary with c's (row 1719's shifted-site sibling, both chains
  emitted). First / last intron of the chain excluded when the moved boundary is the terminal exon's (that is a
  terminal-exon end, row 1714's `first3` / `last5`); reported apart as S1t.
- **Class S2**: equal except two adjacent introns, where the exon between them has both ends moved and still
  overlaps (one internal exon shifted at both ends).
- Yardsticks, both reported: chain in RefSeq; c's differing intron(s) annotated (RefSeq, plus CAT on CHM13).
- Strata: shift <= 10 bp / 11-50 / > 50; shift a multiple of 3; reads(c) / reads(d) < 0.1 / 0.1-0.5 / >= 0.5; splice
  motif of c's differing intron (GT-AG / GC-AG / AT-AC; strict assembly emits no other).
- **Ceiling**: assembly precision (chain in RefSeq) when every S1 + S2 chain not in RefSeq is dropped, against the
  read-count curve at the same number of chains dropped (lowest `reads` first, ties by file order).
- **Bar (assembly level, row 1677's)**: ceiling - curve >= +0.02 on both samples. Below it: no rule on this class can
  be an assembly-precision rule, and the verdict for that use is STOP.
- **Bar (class level)**: step 1 runs if the S1 + S2 chains with an unannotated differing intron are >= 1% of the
  emitted chains on both samples. Below it the row ends at step 0.
- A stratum is called a rule candidate only if it holds >= 500 chains, <= 0.02 of them in RefSeq and <= 0.05 with the
  differing intron annotated, on both samples. Strata are read on OR6737 first; a candidate must then hold on A119b
  as read there, with no cut moved.

## Step 1: the existing read-equivalence snap, genome-wide
One arm per sample: `RUSTLE_POLISH_JUNCTION_SNAP=reads` under the shipped `assemble` (no code change). It moves a
junction onto a stronger one within 10 bp when the weaker one's reads do not align better at their own site.
- changed set = default chains absent from the arm. Bars, each on both samples:
  1. <= 0.02 of the changed set is in RefSeq (class base rate otherwise);
  2. correct chains (RefSeq) lost <= 0.2% of the default's;
  3. row 1714's shifted-site pairs (minor >= 5, four strata) both emitted: no stratum falls by more than 0.01
     (`alt_events.py` on the arm);
  4. the changed set is >= 0.25 of step 0's S1 + S2 chains at <= 10 bp with an unannotated differing intron (reach).
- PASS = a safe opt-in clean-up, already shipped as a switch; it is not a default candidate without the user. FAIL on
  1-3 = the switch harms; FAIL on 4 only = safe but inert.
- A stratum candidate from step 0 (if any) is scored the same way as a plain drop rule on the GTF, bars 1-3, and
  bar 4 replaced by "drops >= 0.25 of the class's unannotated chains".

No default changes, no driver stage, nothing committed.

## RESULT step 0 (2026-11-14; `r1725_shift/{or,a}.txt`, per-chain tables `{or,a}.class.tsv`)
| | OR6737 | A119b |
|---|---|---|
| emitted chains, precision | 73,348, 0.352 | 219,517, 0.176 |
| S1 chains (weaker side, one internal exon end moved) | 11,838 | 26,591 |
| S1 in RefSeq / differing intron annotated | 0.276 / 0.473 | 0.188 / 0.606 (RefSeq + CAT) |
| S1 <= 10 bp | 4,451: 0.317 / 0.516 | 8,274: 0.224 / 0.621 |
| S1 shift multiple of 3 vs not (in RefSeq) | 0.386 vs 0.157 | 0.242 vs 0.137 |
| S1 reads ratio < 0.1 / 0.1-0.5 / >= 0.5 (in RefSeq) | 0.094 / 0.382 / 0.391 | 0.118 / 0.226 / 0.202 |
| S1 <= 10 bp and ratio < 0.1 | 1,646: 0.075 / 0.119 | 2,329: 0.109 / 0.301 |
| S1 GC-AG (in RefSeq) | 281: 0.178 | 814: 0.107 |
| S2 (one exon moved at both ends) | 80: 0.300 / 0.525 | 193: 0.218 / 0.508 |
| S1t (terminal exon end, reported apart) | 13,263: 0.361 / 0.491 | 40,715: 0.199 / 0.467 |
| S1 + S2 with an unannotated differing intron, share of emitted | 6,273 = 0.086 | 10,567 = 0.048 |
| ceiling vs read-count curve | 0.3985 vs 0.3761 = +0.022 | 0.1958 vs 0.1886 = +0.007 |

- **Assembly-level bar: STOP** (+0.02 needed on both; A119b +0.007). Even an oracle that drops every unannotated
  weaker shifted-site chain adds nothing over a plain read cut on A119b.
- **Class-level bar: met** (0.086 / 0.048 >= 0.01), step 1 runs.
- **No stratum is a rule candidate**: the lowest registered stratum (<= 10 bp and reads ratio < 0.1) is 0.075 / 0.109
  in RefSeq and 0.119 / 0.301 intron-annotated, against bars 0.02 / 0.05.
- Reading: no registered stratum isolates an artifact class. On A119b the weaker chain is in RefSeq as often as the
  average emitted chain (0.188 vs 0.176); on OR6737 less often (0.276 vs 0.352), and OR6737 clears the ceiling bar
  (the STOP is A119b's). The 39-53% of moved introns that are unannotated have no yardstick in this row. On A119b the moved intron is annotated in 0.61; at <= 10 bp with comparable reads
  the moved intron is annotated in 0.82 / 0.79. Shifts that keep the frame (multiple of 3) are in RefSeq 1.8-2.5 x as
  often as the others, a graded signal, not a cut. OR6737's intron yardstick is RefSeq only (no CAT), so its 0.47 is
  not comparable with A119b's 0.61.

## RESULT step 1 (2026-11-14; job 4471946, `r1725_shift/{or,a}.snap.txt`, `{or,a}.events.txt`)
`RUSTLE_POLISH_JUNCTION_SNAP=reads`, layers off, against the assemblies above.
| bar | OR6737 | A119b |
|---|---|---|
| chains, precision: default -> arm | 73,348, 0.3517 -> 72,720, 0.3545 | 219,517, 0.1764 -> 219,157, 0.1766 |
| 1. changed set in RefSeq <= 0.02 | 14 of 637 = 0.022 FAIL | 10 of 385 = 0.026 FAIL |
| 2. correct chains lost <= 0.2% | 0.05% PASS | 0.03% PASS |
| 3. row 1714 shifted-site pairs (minor >= 5), largest fall <= 0.01 | 0.003 (int5 > 10 bp) PASS | 0.002 PASS |
| 4. reach of the <= 10 bp unannotated class >= 0.25 | 419 of 1,905 = 0.220 FAIL | 215 of 2,515 = 0.085 FAIL |

**VERDICT row 1725: FAIL.** The existing read-based snap is close to harmless (bar 1 missed by a hair, on 14 and
10 chains) and nearly inert: it changes 0.9% / 0.2% of emitted chains. No stratum of step 0 qualified as a drop
rule, and the assembly-level ceiling stops on A119b (+0.007). HG002 not opened. No default changed, no code changed.
Reach is counted as changed chains carrying a class intron (the class table stores the differing intron, not the
chain), so it can only overstate.

What the row does establish, for reading the output rather than filtering it (no rule was derived): a weaker chain
that differs from a stronger emitted chain by one moved internal exon end is annotated as often as the average
emitted chain on A119b and somewhat less often on OR6737; the moved intron is itself annotated in about half (gorilla, RefSeq) to 0.6 (A119b, RefSeq + CAT);
frame-keeping shifts are about twice as often annotated as frame-breaking ones; and only below 0.1 of the partner's
reads does the class fall to 0.09-0.12, which the read count alone already expresses.

Corrections after review: bar 3 is scored on the four internal strata (int3 / int5), where the largest fall is 0.003
/ 0.002; the terminal stratum first3 <= 10 bp fell 115 -> 113 of 117 on OR6737, an observation outside the bar. The
changed set is the switch, not binary drift: `params.tsv` counts 1,564 / 1,063 collapsed transcripts before the
polish (net 628 / 360 fewer chains after it), and only 9 / 25 chains are new. The default assemblies are those of
rows 1714-1723, not re-run today. Scope of the close: GTF-level strata (shift size, frame, reads ratio, motif) and
the existing edit-distance snap; other read-level features and the frame signal as a score were not tested.
