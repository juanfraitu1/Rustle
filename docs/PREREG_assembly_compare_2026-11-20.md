# Row 1838: where do other assemblers find true chains that Rustle misses, primate against non-primate? (pre-registered 2026-11-20)

User (2026-11-20): a quick comparison of assemblies, to see whether other levers need work to remove overfitting to primates.
Rows 1793, 1811 and 1831-1837 switched levers one at a time. This row asks the other way round: which true chains another
de novo assembler has and Rustle lacks, and where Rustle lost them. Diagnostic, no new run, no default changes. Workspace
`r1838_assembly_compare` (`compare.py`). Committed before the script was run on any dataset.

## Data (existing files only)
| | dataset | Rustle default (`base/x.gtf`, `x.final.gtf`) | StringTie, IsoQuant de novo | truth |
|---|---|---|---|---|
| N | fly simulated | `r1837_real_count_grid/dme` | `apetools/DME_sim_mm` | `r1779_simulated_reads/dme/truth.gtf` |
| N | mouse simulated | `.../mmu` | `apetools/MMU_sim_mm` | `r1779_simulated_reads/mmu/truth.gtf` |
| N | mouse ES, mouse contigs | `.../es` | `apetools/MMU_ES_mm` | GENCODE vM27 |
| P | chimpanzee testis | `.../ptr` | `apetools/PTR_mm` | `ptr_ppa_rules/PTR.refseq.gtf` |
| P | bonobo testis | `.../ppa` | `apetools/PPA_mm` | `ptr_ppa_rules/PPA.refseq.gtf` |
| P | gorilla KB3781 fibroblast | `.../kb` | `apetools/GGO_KB3781_fib_mm` | `apetools/GGO.refseq.gtf` |
The Rustle files are today's defaults (repo `78daaa6b` or later, arm `base` of row 1837). The other tools' files are the de novo
runs of rows 1779 / the ape tool benchmark on the same BAMs. All six datasets were used before (reuse, user ruling).

## Counts (multi-exon, exact contig + strand + intron chain)
1. Per assembly: chains, correct, precision, sensitivity.
2. **Missed:** true chains that StringTie (and, apart, IsoQuant) has and `x.final.gtf` lacks, each put in the first class that fits:
   - FILTERED: in `x.gtf`, absent from `x.final.gtf` (chain score filter, 5' truncation rule, final sub-chain drop).
   - PIECE: a contiguous part of a chain in `x.gtf` (sub-chain collapse, isoform fraction on a contained chain).
   - LONGER: a chain of `x.gtf` is a contiguous part of it (Rustle stopped short: no extension, ends).
   - NEAR: a chain of `x.gtf` has the same number of introns, every boundary within 10 bp (junction snap / shifted junction).
   - OTHER FORM: none of the above, and every intron of it is in some `x.gtf` chain (an alternative combination was dropped:
     isoform fraction, retained-intron rule, read floor).
   - NEW INTRON: at least one of its introns is in no `x.gtf` chain (junction gate, read floor, seeding).
   Beside the classes: the StringTie `longcov` of the missed chains (< 2, 2-5, > 5).
3. The mirror: true chains Rustle `final.gtf` has and StringTie lacks (count only).
4. Wrong chains of `x.final.gtf` and of StringTie, same classes against the truth: PIECE of a true chain, LONGER than one, NEAR
   one, OTHER FORM (all introns true), NEW INTRON.

## Reading, fixed now
- A class points at a lever to look at when its share of the missed chains is at least twice as large on every N dataset as on
  every P dataset (non-primate loss), and when it holds at least 20% of the missed chains on at least two N datasets.
- The reverse pattern (larger on every P) is reported as a primate-side loss.
- Anything else: no class separates the groups; the table is the result.
- A finding here is a pointer for a pre-registered lever test, not a result about a lever.

## Limits, stated in advance
- Three and three datasets; N has two simulated sets (no degraded reads), P has two testis libraries. Tissue and data type are
  confounded with the group.
- Real-library truths are incomplete; the "missed" sets use annotated chains only, so they are true but not all that is true.
- The classes name the stage where a chain is absent, not the rule that removed it; OTHER FORM and NEW INTRON each cover
  several rules.
- StringTie and IsoQuant ran with their own defaults; a chain they keep at 1 read is below Rustle's floor by design.
