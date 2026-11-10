# PREREG r1608: tie groups — unspliced loci carry no strand (2026-11-10)

Follows r1607 (`docs/PREREG_tie_groups_repeat_sd_2026-11-10.md`, FAIL on B5 only: largest group 86 loci = the TAF11L array,
46 loci on `+` and 40 on `-`). An unspliced alignment has no splice site to give it a strand, so the two strand copies of one
unit are one locus. Written and committed before the arm runs.

## Rule (`tie_groups.py build PREFIX as98 0 1 RMSK_BED:0.5 unstranded`)
A tied record with a single aligned block goes to strand `.`; loci are formed per (chrom, strand) as before, so unspliced
records of both orientations fall in one locus. Spliced records keep their strand. A `.` model is emitted with strand `.`
and, in the scorer, can reach or overlap genes on either strand. Everything else as the r1607 filter-on arm.

## Bars (as r1607; M = the 72 indistinguishable members)
- B1g grouped reach >= 58 of 72; B2s purity (co-duplicated-only groups counted as pure) >= 0.90; B3 TAF11L models <= 5 kb
  reach >= 20 of 31; B5 largest group <= 60 loci; B4 job <= 1 h.
- PASS = all five. Control = r1607 filter-on (70 / 0.982 / 31 / 86), not rerun.
- If it passes, the rules are frozen for the gorilla OR6737 transfer row.

## Outcome (2026-11-10, job 4456417, `/anvme/workspace/v138da13-rustle/r1603/uns.*`): FAIL as registered (B5)
| arm | loci kept / set aside as repeat | groups | B1g (>= 58) | B2s (>= 0.90) | B2f strict | B3 (>= 20) | B5 (<= 60) |
|---|---|---|---|---|---|---|---|
| r1607 filter-on (control) | 2,509 / 3,032 | 1,141 | 70 | 220/224 = 0.982 | 0.884 | 31 | 86 |
| unstranded | 5,116 / 4,349 | 2,311 | 71 | 902/907 = 0.994 | 0.917 | 31 | **151** |

- The rule is defective as written: taking the unspliced records out to `.` splits every locus that has spliced and
  unspliced records into two (9,465 loci before the repeat filter against 5,541), so loci and groups grow (largest: 151
  loci on chr16, 127 on chr1). Models by strand: `.` 3,804, `+` 679, `-` 633.
- It does what it was for on TAF11L: 86 loci -> 56 (48 on `.`).
- Judged groups rise 224 -> 907 because a `.` model reaches panel genes on either strand; purity figures are therefore not
  comparable with the control.
- Next row: form loci without strand and give each locus the strand of its spliced records.
