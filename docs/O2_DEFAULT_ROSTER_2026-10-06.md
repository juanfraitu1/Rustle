# O2 with read-level truth on the default families roster: NPIP, TBC1D3, Y (2026-10-06)

**Status: chr16 (NPIP), chr17 (TBC1D3) and chrY (approximation: the 173-copy family left out, Amendment 2) scored. The full chrY run (that family included) is NOT done.**

Protocol: `docs/PREREG_o2_default_roster_2026-10-06.md` (committed c8753249, amended 9e2deefa for chrY before any chrY assignment product existed). Runner and report: `bench/o2_default_roster/`.
Products: `/mnt/linuxdisk/tmp/o2_default_roster_2026-10-06/{chr16,chr17,chrY}/` (`report.txt`, `report.json`, `score/`). HEAD build in `rustle_target_m2/release` (`copy_assign` 87824d91, `mcl_families` a6308244).
Human A119b, CHM13 v2.0, seed 20260925. **Simulation, DEV, first measurement; no rule was changed after a number was seen.**

## Roster provenance and checks

The default roster (HEAD `assemble` + `families` on the contig) is **byte-identical** to the 2026-09-30 e163d955 products wherever those exist: chr16 and chr17 (assembled GTF, `families.gtf`, clusters, loci, copy table), chrY (clusters, loci, copy table).
Rosters: chr16 372 copies in 111 families, chr17 226 copies in 78 families, chrY 458 copies in 64 families (one family of 173). Simulated reads: 8,225 / 5,843 / 10,224.
Mapping parts (read-disjoint, they only bound one call's wall time): chr16 4, chr17 3, chrY 8.
GS1 holds (one primary-or-unmapped record per simulated read: 8,225 and 5,843). GS2 holds (the per-read tallies equal the ALL rows of `score.py reads` for OWN, ANY and UNION).
**GS4 (Amendment 1):** `copy_assign --skip-poa-diagnostic --region-threads 4` gives `assignments.tsv`, `families.tsv` and `family_join.tsv` byte-identical to the default `assign` stage's on chr16 and on chr17.

## NPIP (chr16 roster; target = roster copies overlapping a CAT/Liftoff NPIP copy: 41 copies)

The aligner decides almost every read of these copies: of 1,099 simulated reads, 371 map at MAPQ 60, 695 at MAPQ 1-59 and **33 at MAPQ 0**. The 33 contested reads come from 5 copies: 3 reads from the NPIP family itself (MCL1: copies 11 and 15), 30 from three other families
that overlap NPIP exons (10 reads each, the depth floor). One of the 41 target copies (MCL1|25) enters through the territory of NPIPB13, the one CAT copy with no exons in the truth GTF; it has no MAPQ-0 read.

| reading (33 MAPQ-0 reads) | correct | wrong | conflict | abstain | no row | assigned | accuracy |
|---|---|---|---|---|---|---|---|
| OWN (the true family's row) | 1 | 0 | 0 | 30 | 2 | 1 | 1.00 |
| PRIMARY (rows with `primary_local = 1`) | 1 | 11 | 0 | 19 | 2 | 12 | 0.08 |
| ANY (any assigned row) | 0 | 28 | 1 | 2 | 2 | 29 | 0.00 |
| UNION (arm U2, beside) | 0 | 0 | 0 | 31 | 2 | 0 | n/a |

**Bars: B1 UNDERPOWERED (1 assigned); B2 PRIMARY UNDERPOWERED (12 assigned); B2 ANY UNDERPOWERED (29 assigned, one under the registered 30), with 0 of 29 correct; B3: nothing assigned.** Aligner-primary on the same 33 reads: 0.36.

All 931 MAPQ-0 reads of the chr16 roster (every family; context with power):

| reading | correct | wrong | conflict | abstain | assigned | accuracy | coverage |
|---|---|---|---|---|---|---|---|
| OWN | 27 | 0 | 0 | 902 | 27 | 1.000 | 0.029 |
| PRIMARY | 48 | 78 | 4 | 799 | 130 | 0.369 | 0.140 |
| ANY | 11 | 393 | 106 | 419 | 510 | 0.022 | 0.548 |
| UNION | 0 | 0 | 0 | 929 | 0 | n/a | 0.000 |

OWN: all 27 assigned reads are in the 99.5-100% identity band, none wrong; 188 of the 931 reads come from copies with an identical sibling and none of them is assigned. The recorded legacy-roster run (2026-09-25) shows the same pattern: OWN 163/0, PRIMARY 112/205, ANY 60/432/156, union 0.

## TBC1D3 (chr17 roster; target = roster copies overlapping a CAT/Liftoff TBC1D3 copy: 22 copies)

Of 887 simulated reads from the 22 target copies, 105 map at MAPQ 60, 782 at MAPQ 1-59 and **none at MAPQ 0**: there is no contested read to assign, so B1, B2 and B3 are UNDERPOWERED with 0 assigned.
All 134 MAPQ-0 reads of the chr17 roster (other families): OWN 6 correct / 0 wrong (coverage .045); PRIMARY 26 / 0; ANY 26 / 0; union 0 assigned. On this contig the readings a consumer sees are right wherever they assign (26 of 26).

## Y (chrY roster; target = roster copies overlapping the body of a RefSeq DAZ, RBMY, TSPY, BPY2, HSFY, VCY, CDY, PRY or XKRY gene: 40 copies in 12 families)

**What the 173-copy family is (user question, 2026-10-06): not a gene family but the Yq12 satellite array.** MCL0's 173 members are 161 single-exon loci (median 9.2 kb, up to 57 kb) tiling chrY 31.4-60.8 Mb; their sequence is 33% GGAAT/ATTCC pentamers (other families: 1%) and compresses 3x better (zlib ratio 0.107 against 0.32): the (GGAAT)n DYZ1 / HSat3 heterochromatin. The default `assemble` turns the reads that align over the array into single-exon pseudo-genes and `mcl_families` joins them (96-99.9% identical; a single exon meets the shared-exon test trivially). The simulation is faithful to that roster and is dominated by it: 5,194 of the 10,224 simulated reads (51%) come from MCL0. The other large Y families are real arrays (MCL1 46, MCL2 36, MCL3 35 copies, all at 5.9-10.2 Mb, the TSPY / FAM197Y region; none single-exon-satellite-like).
**Approximation (Amendment 2): the 173-copy family MCL0 is left out of the assignment run** (it holds none of the 40 target copies, and it is not a gene family). The full assignment stalled in one minimap2 child (`-x splice -N 173`, 614 sequences against 173 copies, five cores, 15 minutes without finishing) and was stopped; skipping a family is byte-identical for the other families on chr17
and changes 0.6% of the assignment rows on chr16, so the numbers below carry that uncertainty.

Of 10,224 simulated reads, 609 MAPQ-0 reads come from MCL0 (excluded) and 2,425 from the other 63 families; **452 MAPQ-0 reads come from the 40 target copies** (134 from copies with an identical sibling, 201 from the 99.5-100% band, 117 from the 98-99% band).
**O2 assigns none of them under any reading** (OWN, PRIMARY, ANY and the union arm all give 0 assigned; 438 abstain, 14 have no row), and none of the 2,425 MAPQ-0 reads of the 63 families either. B1, B2 and B3 are UNDERPOWERED (0 assigned). The aligner-primary baseline puts 41% of the 452 reads on their source copy.
Even the 117 reads of the 98-99% band (copies 1-2% diverged from their closest sibling, where chr16 coverage was 0.80) are not assigned: on this roster the certificate has no decisive column for them.

## Reading, with its limits

- Under the registered read model the contested reads of these families are rare (NPIP 33 of 1,099, TBC1D3 0 of 887): reads simulated from a copy's own sequence with 0.13% error align best to that copy, so the aligner decides them. The registered bars therefore cannot be tested at NPIP and TBC1D3 on this
  design, and "UNDERPOWERED" is the outcome. Real A119b reads at NPIP are contested far more often (24% of primaries at MAPQ 0, `docs/PSV_CEILING_2026-10-03.md`): that excess comes from sequence variation that the simulation does not have.
- Where the simulation does have contested reads (all 931 on chr16), the certificate is exact (27 of 27) but assigns 2.9% of them, and the consumer readings are poor (PRIMARY .37, ANY .02): the cross-family claims of `docs/PREREG_o2_read_truth_2026-09-23.md` persist on the default roster. At NPIP, 0 of 29 ANY-assigned reads are right.
- P1 (OWN 0 wrong, >= 30 assigned at NPIP and TBC1D3): 0 wrong holds, the >= 30 does not (1 and 0). P2 (consumer readings below .95 at NPIP and TBC1D3): holds at NPIP (0.08 and 0.00 on 12 and 29 reads), untestable at TBC1D3. P3 (Y UNDERPOWERED on B1): holds (0 assigned, approximation). P4 (union assigns almost nothing): holds (0 of 931, 0 of 134, 0 of 2,425).
- Upper bound on real-data accuracy: no allelic variation, no readthrough, no uncatalogued copies, one error model, CHM13 sequence as the source.
