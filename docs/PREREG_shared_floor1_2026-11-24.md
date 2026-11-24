# PREREG: one one-read assembly for both the one-read tiers and the multi-copy stage (row 1937, 2026-11-24)

Speed. `all` assembles at a one-read floor twice: for the one-read tiers (`oneread`, with the pooling-guard bypass
`RUSTLE_GATE_MIN_READS=1`, default since row 1868) and for the multi-copy stage (`multicopy`, without the bypass; its
rows 1865-1869 judged it that way). On HG002 each costs about 1.9 h. Sharing one run means the multi-copy stage takes
its candidates from the run WITH the bypass (`RUSTLE_MULTICOPY_GATE=1`, an existing switch, "2.4 x the candidates on
OR6737, not judged for this mode"). The user chose to try this (the other direction would undo their row 1868 decision).
Row 1935 showed how sensitive this step is to its input, so the change is judged on what the mode is for, not assumed.

## Test
- Planted human named families, fresh seeds 1948, 1949, 1950, scenarios E and S (sample = reference, ragged ends), the
  mode built both ways on the same BAM; scored by copies (`score.py`), units (`units_score.py`) and labels
  (`em_filter.py --label` at today's default + `label_score.py`).
- Real: OR6737 and KB3781, stage `multicopy` both ways on the runs of rows 1920 / 1928; rescued transcripts, and the
  real-sample arm of rows 1867-1869 (`toolcompare.py --species gorilla`: annotated copies with a spliced transcript,
  junction precision).
- That the shared file is the same: `multicopy`'s own run with the bypass against the `oneread` stage's floor-1 GTF on
  OR6737, by intron chain.

## Claims (bypass against no bypass)
- S1 planted: expressed units found not lower in any set; silent units labelled `expressed` or `shared` higher by at
  most 1 per set and not higher in the sum by more than 2; isoform precision of the mode not lower by more than 0.01.
- S2 real, both samples: annotated copies with a spliced transcript not lower; junction precision not lower by more
  than 0.01.
- S3: the two floor-1 runs give the same transcripts (>= 0.999 shared by chain), so one can stand for the other.
S1 + S2 + S3 = propose that `all` makes one floor-1 assembly and both stages use it; the default is the user's call.

## Predictions (mine)
S3 holds (same command but for the single-exon switch, which is off). S1: units found hold; false labels uncertain,
since more candidates mean more loci for tied reads to land on, though the labels exist to catch exactly those.
S2: copies hold; junction precision is the risk (the bypass admits single-read chains the guard was written to stop).
After row 1935 I give my own predictions here little weight.

## Result (job 4499074; a first submission died on a syntax error of the job script before anything ran)
**S3 HOLDS exactly:** the `oneread` stage's floor-1 run and `multicopy`'s run with the bypass give the same 204,426
transcripts on OR6737 (1.00000 shared by chain).

**S1 HOLDS** (planted, seeds 1948-1950; no bypass / bypass):

| set | expressed units found | silent units with a transcript | of them labelled `expressed` or `shared` | isoform precision of the mode |
|---|---|---|---|---|
| 1948 E, 1949 E, 1950 E | 104 / 104 each | none | none | 0.624, 0.644, 0.629: unchanged |
| 1948 S | 58 / 58 | 27 / 31 | 5 / 5 | 0.424 / 0.418 |
| 1949 S | 58 / 58 | 29 / 31 | 3 / 3 | 0.412 / 0.409 |
| 1950 S | 57 / 57 | 33 / 35 | 3 / 3 | 0.431 / 0.427 |

The bypass rescues two to four more silent units per half-silent set; every one of them is labelled `homology_only`.
Without the labels the mode's false units rise (27 -> 31, 29 -> 31, 33 -> 35): the labels are what makes this safe.

**S2 HOLDS** (real, `toolcompare.py`, the named families' annotated copies with a spliced transcript / junction precision):

| sample | no bypass | bypass | rescued loci (transcripts) |
|---|---|---|---|
| OR6737 | 34 / 0.66 | 37 / 0.67 | 445 (1,262) -> 1,688 (2,992) |
| KB3781 | 27 / 0.63 | 28 / 0.63 | 161 (364) -> 514 (807) |

- Predictions: S3 right; S1 and S2 better than I feared (junction precision did not fall).
- **By the registered rule: proposed.** `all` would make ONE floor-1 assembly (the `oneread` stage's) and the multi-copy
  stage would take its candidates from it. On HG002 that removes about 1.9 h of 5. What changes for the user: three to
  four times as many rescued loci on real samples (1,688 against 445 on OR6737), all carrying a label. Not done; a
  default, the user's call. The stage with the bypass alone takes as long as without (OR6737 1,120 s against 1,163 s):
  the saving exists only when the run is shared.

**Default since 2026-11-24 (user: "yes lets make it default").** `RUSTLE_MULTICOPY_GATE` unset = 1 (0 = the mode without
the bypass, as rows 1865-1869 judged it). The multi-copy and the one-read stage share one floor-1 assembly: each run
records its settings (`f1.settings`: bypass, single-exon switch, duplicate counting, seeding pool and ratio, contigs,
score table, BAM, polishing options) and the stage that runs second takes the other's GTF when the settings are equal
and the file is newer than `PREFIX.gtf`; `RUSTLE_F1_SHARE=0` = each assembles for itself. End to end on OR6737 (job
4499093), both orders: `PREFIX.multicopy.gtf`, `.multicopy.tsv`, `.oneread.gtf`, `.oneread.ranked.gtf`, `.ranked.gtf`
byte-identical to the separate runs. Times: multi-copy then one-read 943 + 308 s; one-read then multi-copy 859 + 21 s;
separate runs about 1,980 s. The reuse rule itself checked on four cases (same settings; other bypass; switch off; GTF
older than the assembly): only the first reuses. Planted results of the mode quoted before this date had no bypass.
