# Row 1876: what is left on complete truth at today's defaults? (step 0, pre-registered 2026-11-21)

User goal (2026-11-21): more empirical rules that raise sensitivity and precision, ahead of the other assemblers, not fitted to
primates. Every earlier class table (rows 1811, 1818, 1838) predates four defaults (junction ratio 20, spliced share 0.05, every
read counted, one-read gate). This row only measures; it builds nothing and changes no default. Workspace `r1876_whats_left`
(`left.py`). Registered before the script ran on any dataset.

## Data: the four datasets whose truth is complete
Today's default assemblies of row 1849, arm `dc` (`r1849_dup_count_rerun/NAME/dc/x.gtf`, `x.final.gtf`), no new assembly.
| | dataset | truth |
|---|---|---|
| fly simulated | `reads/DME/DME_sim_mm.bam` | `r1779_simulated_reads/dme/truth.gtf` |
| mouse simulated | `reads/MMU/MMU_sim_mm.bam` | `r1779_simulated_reads/mmu/truth.gtf` |
| zebrafish simulated | the BAM named in the arm's log | `r1800_snap_ism/dre/truth.gtf` |
| SIRV (real reads, in the mouse ES library) | `reads/MMU/MMU_ES_mm.bam`, contigs `SIRV*` | `genomes/MMU/MMU.ref.gtf`, contigs `SIRV*` |
All four are development data for every rule that follows from this row. Unit: multi-exon intron chains (`bench/famcmp/curve.py load`).

## Measures, fixed now
**Reads per chain.** Primary records (flag & 0x904 = 0) whose whole intron list equals the chain: records, and distinct
(start, end) pairs.

**A. Sensitivity by support.** Truth chains by distinct full-chain records (0, 1, 2, 3-5, 6-20, > 20): how many, and the share in
`final.gtf`.

**B. Missed truth chains with >= 2 distinct records** (the floor is met, a rule removed them). First class that fits:
1. FILTERED: in `x.gtf`, not in `final.gtf`.
2. SUB-CHAIN: a contiguous run of a chain in `x.gtf`; reported with records / the container's `reads`.
3. LOW SHARE: records < 0.02 x the largest `reads` among `x.gtf` chains sharing an intron with it (the isoform fraction).
4. NEAR JUNCTION: one of its introns is absent from `x.gtf` and lies within 10 bp of an `x.gtf` intron.
5. INTRON ABSENT: another of its introns is in no `x.gtf` chain.
6. OTHER: every intron reported, chain not.

**C. Wrong chains of `final.gtf`.** First class that fits:
1. PIECE of a true chain: cut at the 5' side, the 3' side, or both.
2. COMBINATION: every intron is a true intron, the combination is not; sub-class SKIP when it equals a true chain with one or more
   inner exons left out and the spanning intron true elsewhere.
3. NEAR: a new intron within 10 bp of a true one (1-2 bp, 3 bp, 4-10 bp, by the smallest total shift).
4. NEW INTRON, some true: by the new intron's two sites (donor true / acceptor true / both true, pairing new / neither).
5. NO TRUE INTRON: overlaps a truth transcript on the same strand, on the other strand only, or none.
Each with the chain's `reads` (2, 3-5, 6-20, > 20) and the share whose reads are all MAPQ < 20.

## Reading, fixed now
- A class of B or C is a LEAD when it holds >= 15% of its table on at least 3 of the 4 datasets, or >= 30% on the real one (SIRV)
  and one simulated set. Leads get their own registered row; nothing is proposed from this row alone.
- Table A is the ceiling statement: the share of missed truth chains that have fewer than 2 distinct full-chain records is what no
  rule on full-chain reads can recover.

## Limits, stated in advance
- Simulated reads carry the simulator's errors, not a library's (no intronic or antisense signal, uniform truncation model).
- SIRV has 60 spliced truth chains: its shares are small numbers and are printed as counts.
- Classes are assigned to the first that fits; a chain with two faults is counted once.
- Table B's class 3 uses intron sharing as the locus, an approximation of the assembler's locus.

## Result (2026-11-21): the rules that drop are spent; two leads, one of them tested in row 1878

Jobs 4494561-4494564, exit 0; tables `r1876_whats_left/left.NAME.txt`, per-chain lists `NAME.missed.tsv`, `NAME.wrong.tsv`.
`final.gtf`, correct / precision / sensitivity: fly 6,753 / 0.973 / 0.720; mouse 25,077 / 0.970 / 0.578; zebrafish 12,129 /
0.912 / 0.665; SIRV 48 / 1.000 / 0.800.

**A. The ceiling.** Truth chains by distinct full-chain records, and the share found:
| records | fly | mouse | zebrafish | SIRV |
|---|---|---|---|---|
| 0 | 221, 0 | 2,123, 0 | 1,495, 0 | 0 |
| 1 | 1,770, 0 | 8,035, 0 | 3,539, 0 | 1, 0 |
| 2 | 861, 0.787 | 3,890, 0.500 | 1,633, 0.833 | 0 |
| 3-5 | 1,115, 0.838 | 5,253, 0.556 | 2,127, 0.871 | 0 |
| 6-20 | 2,060, 0.900 | 9,020, 0.686 | 3,671, 0.911 | 4, 0.250 |
| > 20 | 3,352, 0.980 | 15,102, 0.928 | 5,786, 0.962 | 55, 0.855 |
76% (fly), 55% (mouse), 82% (zebrafish) of the missed truth chains have fewer than two distinct full-chain records: no rule on
full-chain reads reaches them. Chains with exactly one record are 19% of the truth on all three simulated sets, the largest
reservoir, and the one rows 1805-1808 and 1866 already worked (one-read tiers).

**B. Missed truth chains with >= 2 distinct records** (635 / 8,188 / 1,088 / 11), share per class:
| class | fly | mouse | zebrafish | SIRV | lead |
|---|---|---|---|---|---|
| LOW SHARE (isoform fraction) | 0.630 | 0.343 | 0.313 | 0 of 11 | LEAD by the reading; closed, see below |
| SUB-CHAIN of a reported chain | 0.137 | 0.468 | 0.298 | 3 of 11 | LEAD by the reading; closed by rows 1839-1843 |
| INTRON ABSENT from `x.gtf` | 0.113 | 0.061 | 0.179 | 2 of 11 | no |
| NEAR JUNCTION | 0.068 | 0.032 | 0.067 | 1 of 11 | no |
| FILTERED | 0.046 | 0.039 | 0.128 | 2 of 11 | no |
| OTHER | 0.006 | 0.056 | 0.015 | 3 of 11 | no |

**C. Wrong chains of `final.gtf`** (187 / 787 / 1,170 / 0), share per class:
| class | fly | mouse | zebrafish | lead |
|---|---|---|---|---|
| PIECE (5' cut only) | 0.230 (0.214) | 0.390 (0.355) | 0.313 (0.218) | LEAD by the reading; closed by rows 1064, 1486, 1715, 1840, 1841 |
| NEW INTRON, some true | 0.193 | 0.219 | 0.215 | LEAD; its largest sub-class is "both sites true, pairing new" (0.102 / 0.141 / 0.066) |
| NO TRUE INTRON | 0.305 | 0.086 | 0.213 | no (2 of 3) |
| NEAR a true intron | 0.128 | 0.154 | 0.152 | no (2 of 3; taken up in row 1879) |
| COMBINATION of true introns | 0.144 | 0.151 | 0.108 | no |

**Post hoc looks (after the tables; `probe.py`, `frac_look.py`, `attrib.py`, `combo.py`, `shortint.py`).**
- **Skipped short exons.** Wrong chains with a new intron that skips exactly one true exon: 19 / 121 / 77, the exon 15 bp or
  shorter in 19 / 111 / 64. This is the aligner dropping a short exon: 10%, 14%, 5% of the wrong chains. Taken to row 1878.
- **LOW SHARE on real libraries** (fraction off against today's default, existing arms of row 1864): the chains the fraction
  removes are 0.055 correct on mouse ES (library 0.410) and 0.258 on WTC11 (0.639). Exempting chains made only of introns that
  reported chains carry does not help (0.049 / 0.282), nor does a read count (> 20 reads: 0.19-0.27 on ES). On simulated reads
  every low-share chain is true, so the simulated sets show only the cost. The fraction stays a trade (rows 1837, 1849).
- **Which default removed the rest** (class INTRON ABSENT, mouse / zebrafish / fly): the spliced-share rule 114 / 35 / 27 chains
  (its known cost, row 1847), the strict junction gate most of the others: 190 / 134 / 28 of these chains carry an intron that is
  not GT-AG, GC-AG or AT-AC (0.24% / 1.31% / 0.16% of annotated introns).
- **NO TRUE INTRON** is reads placed on an unannotated copy: all reads MAPQ < 20 in 43 of 57 (fly), 25 of 68, 155 of 249
  (zebrafish). A MAPQ rule would remove them and is NOT proposed: on real genomes those are duplicated genes (row 1854).
- **Zebrafish truth holds 67 "introns" of 30 bp or less** (2-18 bp, annotation frame fixes) that reads cross as deletions: 72
  wrong chains and as many missed ones are this artifact of the annotation, not of the assembler.

**Reading.** Every class that passes the lead bar on the drop side (fraction, sub-chain, 5' pieces) is a registered trade or a
closed lever. The one new lead with both axes in it is the skipped short exon: a wrong chain that becomes its true one.
