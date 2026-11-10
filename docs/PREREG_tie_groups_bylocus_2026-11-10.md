# PREREG r1609: tie groups — loci formed without strand, strand from the spliced records (2026-11-10)

Follows r1608 (`docs/PREREG_tie_groups_unstranded_2026-11-10.md`, FAIL: the rule split loci holding spliced and unspliced
records). Written and committed before the arm runs.

## Rule (`tie_groups.py build PREFIX as98 0 1 RMSK_BED:0.5 locus`)
Loci are formed per chromosome with the strand ignored (same seed-overlap rule). A locus then takes the strand most of its
spliced records have; with no spliced record, or an even split, its strand is `.`. A `.` model can reach or overlap genes on
either strand in the scorer. Everything else as the r1607 filter-on arm.

## Bars (as r1607; M = the 72 indistinguishable members)
- B1g grouped reach >= 58 of 72; B2s purity (co-duplicated-only groups counted as pure) >= 0.90; B3 TAF11L models <= 5 kb
  reach >= 20 of 31; B5 largest group <= 60 loci; B4 job <= 1 h.
- PASS = all five. Controls on file: r1607 filter-on (70 / 0.982 / 31 / 86), r1608 (71 / 0.994 / 31 / 151).
- Known cost, reported: two genes on opposite strands that overlap become one locus.
- If it passes, the rules are frozen for the gorilla OR6737 transfer row.

## Outcome
(empty until the run)
