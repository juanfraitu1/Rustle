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

## Outcome (2026-11-10, job 4456425, `/anvme/workspace/v138da13-rustle/r1603/byl.*`): FAIL as registered (B5)
| arm | loci kept / repeat | groups | B1g (>= 58) | B2s (>= 0.90) | B2f strict | B3 (>= 20) | B5 (<= 60) |
|---|---|---|---|---|---|---|---|
| r1607 stranded (control) | 2,509 / 3,032 | 1,141 | 70 | 0.982 | 0.884 | 31 | 86 |
| by locus | 1,583 / 2,327 | 673 | 61 | 200/207 = 0.966 | 0.826 | 28 | **93** |

- B4: 1,365 s - pass. Models by strand: `+` 531, `-` 482, `.` 570.
- The TAF11L group is now 46 loci (86 stranded), close to the array's 39 annotated units. B5 is missed by another group:
  93 loci on chrY (the stranded run had it as three groups of 46 / 40 / 32), then 73 on chr1.
- Costs against the stranded run: grouped reach 70 -> 61 (lost DEFA1, DEFA1B, DEFB105A, H2AB3, MBD3L4, SPDYE13, three
  TAF11L), strict purity 0.884 -> 0.826.
- Reading: merging strands fixes the doubled TAF11L loci and loses reach elsewhere; the largest groups in both runs are
  single-chromosome arrays (chr5, chrY, chr1), not chaining. A fixed limit of 60 loci does not separate chaining from a
  large real family.
