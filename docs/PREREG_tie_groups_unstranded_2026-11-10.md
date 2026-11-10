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

## Outcome
(empty until the run)
