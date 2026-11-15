# PREREG r1626: single-exon models in the de novo and guided outputs (2026-11-10)

User (2026-11-10): "can we focus on fixing the single exon representation in both de novo and guided assemblies?"
Bars and the guard below were written after the gorilla OR6737 diagnosis and before any A119b single-exon number was read.

## Measure (`bench/famcmp/mono_score.py`, `mono_reads.py`)
- A reference (RefSeq) single-exon transcript is EXPRESSED at k when >= k unspliced primary alignments on its strand overlap
  it by >= 50% of both lengths. It is FOUND by a single-exon model overlapping it by >= 50% of both lengths on its strand.
- E3 = expressed-at-3 transcripts found. P = single-exon models that match a reference single-exon transcript / all
  single-exon models. RefSeq is incomplete for single-exon genes, so P is a relative figure only.
- Guard set S (the multi-copy class): reference single-exon transcripts the DEFAULT `PREFIX.gtf` finds that have < 2
  unspliced primary alignments of their own (found through secondary alignments or spliced primaries).

## Diagnosis on gorilla OR6737 (seen)
RefSeq single-exon transcripts 10,675; expressed at 2 / 3 / 10 reads: 510 / 402 / 227.
| output | single-exon models | E3 (of 402) | found at >= 10 reads (of 227) | guard set S found |
|---|---|---|---|---|
| ours `PREFIX.gtf` (default) | 563 | 75 | 60 | 86 (definition) |
| ours `PREFIX.final.gtf` (layers, `--mono-model cluster`) | 1,546 | 134 | 107 | 0 |
| StringTie | 758 | 97 | 93 | - |
| FLAIR | 14,850 | 283 | 165 | - |
The 120 missed by `final.gtf` at >= 10 reads: inside an exon of a same-strand spliced model 44, inside the span of one 22,
nothing at the locus 35 + 10 mitochondrial, opposite-strand spliced exon over it 6, other 3. The default's 167 misses:
41 / 27 same-strand exon / span, 30 an over-long single-exon model, 30 an opposite-strand spliced exon over it (the shadow
rule's either-strand clause, not a strand error), 36 nothing, 3 other. Of S, 55 have no primary evidence at all and 30 are
covered by spliced primaries: `final.gtf` loses all 86 because cluster mode builds from primary unspliced records only and
replaces every single-exon model (the E5 loss of r1490).

## Scope
`PREFIX.gtf` and `PREFIX.families.gtf` (the families input) are NOT changed: cluster mode in the families input failed on
A119b (r1490). The fix targets `PREFIX.final.gtf`, and the guided output is then built from it (decided with the result).

## Step 0 (attribution, descriptive): one gorilla run with `RUSTLE_MONO_DUMP`
Every cluster of >= 3 reads with each filter's verdict (same-strand spliced exon overlap, internally primed, canonical PAS).
Offline: the expressed reference transcripts each filter kills, and the unannotated clusters each filter currently blocks.
Regression: with nothing set, the rebuilt binary's `or.final.gtf` must equal the existing one byte for byte.

## Arms (chosen after step 0, from: keep the default's single-exon models that no cluster replaces; PAS not required;
3'-end-only clustering; never 3'-end-only together with a relaxed overlap rule)
## Bars for an arm, gorilla (dev) then A119b (held out, once)
- B1 E3 >= 1.2 x `final.gtf`'s.
- B2 P >= `final.gtf`'s - 0.05.
- B3 guard: S found >= 0.9 x the default's (gorilla 86; A119b counted from its own default).
- B4 multi-exon chains identical to `final.gtf`'s.
- B5 single-exon models >= 20 kb: none added.
PASS on both samples: proposed as `final.gtf`'s single-exon rule and the guided input. Otherwise reported; defaults stay.

## Step 0 outcome (2026-11-10, job 4459480, `mono/g0.dump.tsv`; regression: both GTFs byte-identical)
Clusters of >= 3 reads: 22,460; 1,546 pass every filter. The 402 expressed-at-3 reference transcripts:
| | expressed at 3 (402) | at 10 (227) | clusters matching no expressed reference transcript (21,070) |
|---|---|---|---|
| a cluster passes | 134 | 107 | 1,207 |
| no cluster (ends disagree by > 50 bp) | 119 | 9 | - |
| blocked by same-strand spliced exon overlap only | 48 | 44 | 13,871 |
| blocked by no PAS only | 45 | 34 | 376 |
| blocked by primed only | 20 | 12 | 76 |
| primed and no PAS | 15 | 6 | 423 |
| exon overlap with another filter | 21 | 15 | 5,117 |
The exon-overlap rule protects against 13,871 clusters for 48 transcripts: not touched. PAS and priming cost 80 transcripts
for 875 clusters. Ragged ends cost 119, almost all at low depth.

## Arms (frozen before any arm was run; `RUSTLE_MONO_RULE`, unset = today's rule)
- `keep`: an existing single-exon model with < 3 primary unspliced records on its strand and a span < 20 kb is kept (guard S).
- `keep,pas`: also, a cluster needs no canonical PAS (still dropped when primed).
- `keep,either`: also, a cluster is dropped only when it is both primed and PAS-less.
- `keep,ends3`: clusters on the 3' end alone; 5' end = the 5'-most start that max(2, 10%) of the reads reach.
- `keep,either,ends3`.
The arm proposed = the passing arm with the highest E3; ties by P.

## Gorilla arm outcome (2026-11-10, job 4459540): every arm FAILS as registered, on B3 alone for two arms
Base `or.gtf` byte-identical in every arm. Bars: E3 >= 161, P >= 0.169, S >= 77 of 86.
| arm | single-exon models | E3 (of 402) | P | S found | chains = final | models >= 20 kb | verdict |
|---|---|---|---|---|---|---|---|
| default (`or.gtf`) | 563 | 75 | 0.298 | 86 | - | 6 | - |
| final today | 1,546 | 134 | 0.219 | 0 | yes | 0 | - |
| keep | 1,743 | 134 | 0.239 | 75 | yes | 0 | FAIL B1, B3 |
| keep,pas | 2,536 | 180 | 0.329 | 75 | yes | 0 | FAIL B3 (75 < 77) |
| keep,either | 2,653 | 199 | 0.330 | 75 | yes | 0 | FAIL B3 (75 < 77) |
| keep,ends3 | 1,814 | 167 | 0.158 | 75 | yes | 0 | FAIL B2, B3 |
| keep,either,ends3 | 2,835 | 243 | 0.190 | 75 | yes | 0 | FAIL B3 |
No A119b run (no arm passed). Post hoc, the 11 lost guard transcripts: none is a >= 20 kb model; all have default models whose
span holds >= 3 primary unspliced records by plain overlap (so `keep` hands them to the cluster rule, which emits nothing)
while the measure counts < 2 by 50% reciprocal overlap; 6 are one tandem array (NC_073246.2 3.75-4.06 Mb).
`ends3` halves precision (0.330 -> 0.190 with `either`): 3'-only clustering is closed.

## Amendment 1 (2026-11-10, post hoc: written after the gorilla arm table; user approved)
The 11 lost guard transcripts are a counting mismatch between `keep` (plain overlap) and the measure (50% reciprocal).
New switch `fallback` (replaces `keep`): an existing single-exon model with a span < 20 kb is removed only when a cluster
model is emitted overlapping it on its strand; otherwise it stays. Arms: `fallback`, `fallback,pas`, `fallback,either`
(no `ends3`). Bars B1-B5 unchanged, same gorilla numbers (161 / 0.169 / 77). Gorilla is now a twice-used dev set, so the
verdict that counts is A119b, run once, on the single arm chosen by the frozen order (highest E3 among passing, ties by P).

## Amendment 1 gorilla outcome (job 4459600): all three arms PASS; base `or.gtf` byte-identical
| arm | single-exon models | E3 (of 402) | P | S found (of 86) | chains = final | models >= 20 kb | verdict |
|---|---|---|---|---|---|---|---|
| fallback | 1,902 | 162 | 0.240 | 83 | yes | 0 | PASS |
| fallback,pas | 2,655 | 196 | 0.324 | 83 | yes | 0 | PASS |
| fallback,either | 2,768 | 212 | 0.325 | 83 | yes | 0 | PASS |
Chosen by the frozen order: `fallback,either`. A119b held out, run once (baseline `r1617/default/a.{gtf,final.gtf}`).

## A119b held-out outcome (2026-11-10, job 4459747): FAIL on B2 (precision); B1, B3, B4, B5 pass
RefSeq single-exon expressed at 3: 1,561. Bars: E3 >= 370, P >= 0.116, S >= 80 of 89. Base `a.gtf` byte-identical.
| | single-exon models | E3 | P | S found | chains = final | models >= 20 kb |
|---|---|---|---|---|---|---|
| default (`a.gtf`) | 13,508 | 265 | 0.027 | 89 | - | 517 |
| final today | 9,214 | 308 | 0.166 | 0 | yes | 0 |
| fallback,either | 35,784 | 609 | 0.080 | 87 | yes | 0 |
Found doubles (609 vs 308) but models quadruple and precision halves. Not adopted; `RUSTLE_MONO_RULE` stays an
experimental opt-in, defaults unchanged. Not split here: how many of the 26,570 added models come from `fallback`
(A119b's default single-exon set is at P 0.027, gorilla's at 0.298) and how many from `either`.

## Post hoc split of the `fallback,either` models by source (2026-11-10; no new run, not a registered test)
| source | A119b models | P | E3 (of 1,561) | S (of 89) | gorilla models | P | E3 (of 402) | S (of 86) |
|---|---|---|---|---|---|---|---|---|
| clusters already in final today | 9,214 | 0.166 | 308 | 0 | 1,546 | 0.219 | 134 | 0 |
| new clusters (`either`) | 14,712 | 0.079 | 307 | 0 | 910 | 0.503 | 93 | 0 |
| kept default models (`fallback`) | 11,858 | 0.016 | 83 | 87 | 312 | 0.327 | 13 | 83 |
| `either` alone (all clusters) | 23,926 | 0.112 | 526 | 0 | 2,456 | 0.325 | 199 | 0 |
| `fallback` alone | 21,072 | 0.081 | 391 | 87 | 1,858 | 0.237 | 147 | 83 |
`either` alone on A119b: found 526 (bar 370 passes), precision 0.112 (bar 0.116 misses by 0.004), guard 0. The kept default
models are the precision cost on A119b (11,858 at 0.016, for 83 expressed and the 87 guard transcripts).

## Post hoc context, same day (not registered tests; existing outputs only)
**Cluster models of `fallback,either` by read count** (precision vs RefSeq single-exon):
| reads | A119b models | P | gorilla models | P |
|---|---|---|---|---|
| 3-4 | 16,146 | 0.066 | 1,430 | 0.197 |
| 5-9 | 5,671 | 0.154 | 622 | 0.359 |
| 10-49 | 1,958 | 0.329 | 371 | 0.717 |
| 50+ | 151 | 0.722 | 33 | 0.818 |
About 60% of the unmatched low-read clusters lie inside a same-strand annotated transcript span in both samples. Reading:
a fixed 3-read floor means less in the deeper library (depth inferred from cluster counts, not measured); a Kinnex
artefact is not ruled out.

**Trusted baselines on the same measure** (existing benchmark GTFs, not rerun):
| | A119b models | E3 (of 1,561) | P | gorilla models | E3 (of 402) | P |
|---|---|---|---|---|---|---|
| StringTie | 38,635 | 281 | 0.008 | 758 | 97 | 0.135 |
| FLAIR | 468,993 | 990 | 0.004 | 14,850 | 283 | 0.028 |
Every tool loses precision on A119b; final today (0.166) is the most precise there.

**The 29 single-exon IsoCon-only transcripts of r1515 (chr20, A119b):** found by default 0, final 0, `fallback,either` 0,
StringTie 1, FLAIR 6. Own unspliced primary reads (50% reciprocal): none 17 (miRNA/snRNA/pseudogenes of 50-500 bp; credited
to IsoCon by SQANTI's loose single-exon match), 1-2 reads 3, >= 3 reads 9 (GTSF1L 150, FNDC11 45, AL121832.3 30, MOCS3 17,
SNX5 15, AL031666.3 7, RPS3P2 6, ID1 6, C20orf144 5). At most of the 9 a same-strand spliced model covers the locus
(GTSF1L: spliced models with 5 and 2 reads); the blocking filter was not confirmed per transcript. Possible lever, not
sized: let a cluster through when its reads far exceed the overlapping spliced models'.
