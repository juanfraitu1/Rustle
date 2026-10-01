# IsoCon (reference-free) on KB3781: simulation and real fibroblast reads (Amendments 1-3), 2026-10-01

Prereg `docs/PREREG_rna_allele_haplotype_count_2026-10-01.md` (Amendments 1-3); truth frozen in `docs/RNA_ALLELE_TRUTH_FROZEN_2026-10-01.md`.
IsoCon 0.3.3 (PyPI; env `isocon`: python 3.8, networkx 2.3, numpy 1.19.5), `IsoCon pipeline --nr_cores 4`, defaults. Scripts
`bench/rna_allele/isocon_sim.py`, `isocon_score_sim.py`, `isocon_score_real.py`; work dir `/mnt/linuxdisk/tmp/rna_allele/isocon/`.
Descriptive arm: no pass/fail.

## Simulation (truth known exactly; 20 reads per haplotype-copy transcript, 1% errors, shortened ends)

| | NPIP (3' 4 kb window, Amendment 3) | TBC1D3 |
|---|---|---|
| truth transcripts (haplotype copies) | 44 (25 copies on `_pri`, 18 B alleles, 1 B-only) | 28 (14 copies, 14 B alleles) |
| IsoCon outputs | 77 | 54 |
| recovered (identity >= 0.999) | **44 / 44** (43 exact) | 23 / 28 (22 exact) |
| copies whose two alleles differ: separated / merged | **10 / 0** | 8 / 5 |
| copies with identical alleles (cannot be told apart) | 8 | 1 |
| outputs tied between two different copies (paralogs merged) | 0 | 5 |
| reference-absent (B-only) copy recovered | 1 / 1 | - |
| outputs per recovered transcript | 1-6 (same sequence, different 5' ends) | 1-10 |

Closest pair of different copies: NPIP 749 edits, TBC1D3 6 edits; allele differences: NPIP 0-1,017, TBC1D3 0-21. Where paralogs are as
close as alleles (TBC1D3), IsoCon merges both. The output count is not a copy count (77 for 44, 54 for 28).

The NPIP simulation at full annotated length (6-27 kb transcripts) did not finish in 10 minutes (Amendment 3).

## Real fibroblast reads (KB3781)

| | NPIP | TBC1D3 |
|---|---|---|
| reads in the net (any alignment record on a copy) | 2,951 (all used; 93 s) | 1,264 |
| copies expressed (>= 2 primaries) | 24 of 25 | 0 of 14 (not expressed in fibroblasts) |
| IsoCon outputs | 108 | 86 |
| outputs on family copies | 31 (9 copies; 18 on NPIPA2) | 0 |
| outputs that are other genes caught by the net (>= 0.99 elsewhere) | 76 (mostly PDXDC1, which shares sequence with the NPIP blocks) | 86 |
| outputs matching nothing | 1 | 0 |
| haplotype copies recovered (>= 0.999) | **13 of 43** (haplotype copies of the 24 expressed copies) | - |
| expressed copies whose alleles differ (13): separated / merged / not recovered | **1 / 1 / 11** | - |
| outputs tied between copies | 5 | - |
| reference-absent maternal copy recovered | 0 / 1 | - |

## Reading

- With ideal reads, reference-free clustering does what the advisor expects on NPIP: every haplotype copy, every differing allele pair
  and the reference-absent copy come out. It fails where paralogs are as close as alleles (TBC1D3): 5 of 13 allele pairs merged and 5
  outputs merge different copies.
- With real reads it recovers a third of the expressed NPIP haplotype copies and separates alleles for 1 of 13 copies. Real reads are
  5'-truncated and spread over isoforms at ~30-180 primaries per copy, so few reads converge on one sequence (IsoCon needs >= 2).
- IsoCon needs a read set to start from, and any net wide enough to catch reference-absent copies also catches the genes that share
  sequence with the family (PDXDC1): 76 of its 108 NPIP outputs are other genes. Its output count is not a copy count in either arm.
