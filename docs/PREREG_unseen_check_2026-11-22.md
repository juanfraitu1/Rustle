# Row 1887: the two defaults of row 1886 on libraries neither rule has seen (pre-registered 2026-11-22)

Row 1886 made the short-exon repair a default of `assemble` and writes `PREFIX.ranked.gtf`. Both were judged on datasets that
were in hand; the user's standing rule (2026-11-10 / reuse ruling) asks for one more unseen dataset per new rule. ENCODE
ENCSR507JOF WTC11 replicates 1 and 3 (`reads/HSA/WTC11_r1_mm.bam`, `_r3_mm.bam`) were mapped as check data in row 1862.
Replicate 3 was assembled once (row 1866, one-read promotion) and replicate 1 was a partner library there; neither has been run
with the repair, and no other assembler has been run on either. Human, so this checks transfer to unseen libraries, not to
another species (the non-primate side of both rules is rows 1878 and 1880: 7 of 7 datasets each).

## Runs (driver of commit 1886, nothing set unless said; `RUSTLE_PTC=0`)
Per replicate: `assemble` with `RUSTLE_SHORTEXONS=0` (before), `assemble` (default), `oneread` on the default; StringTie 3.0.3
`-L` and IsoQuant 4.1.0 with the commands of `apetools/job.sbatch`. Truth: RefSeq on CHM13 (`ptr_ppa_rules/gffcmp/hsa.ref.gtf`).
Workspace `r1887_unseen_check`.

## Reading, fixed now
- **Repair CONFIRMED** on a replicate when the default `final.gtf` has at least as many correct chains and at least the
  precision of the run without it, one strictly better (row 1878's "pays"). Confirmed overall = both replicates.
- **List CONFIRMED** on a replicate when `PREFIX.ranked.gtf`, cut at StringTie's and at IsoQuant's number of multi-exon chains,
  holds at least as many correct chains as each, and its `list_sensitive` set has at least StringTie's correct chains and at
  least its precision. Confirmed overall = both replicates.
- A miss is reported as a miss; it does not by itself undo a default (the user's call).

## Result (2026-11-22): both rules CONFIRMED on both unseen replicates

Jobs 4494803 (r1) and 4494804 (r3), every step exit 0; driver of the row 1886 commit, `r1887_unseen_check/check.txt` (`check.py`).
`final.gtf` = chains / correct / precision.
| replicate | without the repair | default (repair on) | new (correct) | lost (correct) | repair |
|---|---|---|---|---|---|
| WTC11 r1 | 41,830 / 25,000 / 0.5977 | 41,847 / 25,048 / 0.5986 | 140 (63) | 123 (15) | CONFIRMED |
| WTC11 r3 | 43,556 / 25,423 / 0.5837 | 43,571 / 25,468 / 0.5845 | 128 (64) | 113 (19) | CONFIRMED |

| replicate | StringTie | list at its size | IsoQuant | list at its size | sensitive cut | list |
|---|---|---|---|---|---|---|
| WTC11 r1 | 49,903 / 24,209 / 0.4851 | 26,607 (+2,398) | 31,406 / 18,765 / 0.5975 | 21,598 (+2,833) | 48,114 / 26,284 / 0.5463: ahead of StringTie on both axes | CONFIRMED |
| WTC11 r3 | 50,541 / 24,648 / 0.4877 | 26,734 (+2,086) | 31,034 / 18,538 / 0.5973 | 21,391 (+2,853) | 49,968 / 26,652 / 0.5334: ahead of StringTie on both axes | CONFIRMED |

- The repair's new chains are annotated 3 to 4 times as often as the ones it removes (45-50% against 12-17%), as on WTC11 r2.
- The default `final.gtf` alone has more correct chains and higher precision than StringTie on both replicates (25,048 / 0.599
  and 25,468 / 0.585 against 24,209 / 0.485 and 24,648 / 0.488). Against IsoQuant (18,765 / 0.598 and 18,538 / 0.597) it has
  a third more correct chains on both; it is also more precise on r1 (0.5986 against 0.5975) and 0.013 less precise on r3
  (0.5845 against 0.5973). The registered comparison is the list at IsoQuant's size, which is ahead on both.
- Limits: both replicates are the library type of WTC11 r2, on which the repair and the list were already seen to work; this
  is a check on unseen reads, not on an unseen species. WTC11 r1 and r3 are now used for both rules.
