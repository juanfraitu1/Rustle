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
