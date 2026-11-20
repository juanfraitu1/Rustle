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

## Result (2026-11-20): no class separates the groups by the registered reading; the only clean difference is on the primate side

Tables `r1838_assembly_compare/cmp.NAME.txt` (`compare.py`), Rustle files of row 1837's arm `base`.

**Assemblies** (`final.gtf`: correct chains / precision; StringTie; IsoQuant):
| | fly sim | mouse sim | mouse ES | chimp testis | bonobo testis | KB3781 |
|---|---|---|---|---|---|---|
| Rustle `final.gtf` | 7,041 / 0.964 | 27,159 / 0.953 | 17,731 / 0.347 | 20,847 / 0.416 | 17,884 / 0.337 | 25,745 / 0.373 |
| StringTie | 6,549 / 0.963 | 23,705 / 0.959 | 16,402 / 0.367 | 21,532 / 0.350 | 17,288 / 0.281 | 22,087 / 0.380 |
| IsoQuant | 6,592 / 0.972 | 24,098 / 0.967 | 13,443 / 0.471 | 15,761 / 0.441 | 13,109 / 0.369 | 15,785 / 0.427 |
| true in StringTie, not in `final.gtf` | 160 | 1,101 | 899 | 3,989 | 2,272 | 1,688 |
| true in `final.gtf`, not in StringTie | 652 | 4,555 | 2,228 | 3,304 | 2,868 | 5,346 |

**Where the chains StringTie has are lost** (share of the missed chains):
| class | fly sim | mouse sim | mouse ES | chimp | bonobo | KB3781 | reading |
|---|---|---|---|---|---|---|---|
| FILTERED | 5.0% | 18.5% | 10.0% | 0.5% | 1.1% | 3.3% | larger on every N, not twice (5.0 against 3.3) |
| PIECE | 13.1% | 23.8% | 26.5% | 2.4% | 4.0% | 7.2% | larger on every N, not twice (13.1 against 7.2) |
| LONGER | 6.9% | 11.4% | 16.0% | 47.8% | 38.9% | 46.9% | **primate-side loss** (twice on every pair) |
| NEAR | 5.6% | 2.0% | 1.6% | 1.2% | 1.3% | 1.8% | none |
| OTHER FORM | 6.2% | 5.6% | 9.8% | 7.2% | 6.6% | 7.3% | none |
| NEW INTRON | 63.1% | 38.6% | 36.2% | 40.9% | 48.2% | 33.6% | none |
StringTie `longcov` < 2 of the missed chains: 9% / 12% / 42% on N, 77% / 73% / 64% on P.

- **Registered verdict:** no class meets the non-primate bar. PIECE and FILTERED lean that way on all nine pairs and miss the
  factor 2 on fly. Both name levers already known and toned down: the sub-chain collapse ratio (row 1810: 0.3 still costs 7.6%
  on mouse simulated; 0.1 is the non-primate choice) and the `final.gtf` filter. With IsoQuant as the reference PIECE is larger
  still (22% / 41% / 60% against 8% / 14% / 24%), also short of the factor 2.
- **Primate side:** on the three ape libraries 39-48% of the missed chains are ones Rustle built shorter, and two thirds to
  three quarters of all missed chains have under 2 reads of StringTie coverage. That is the transfrag extension (off since row
  1812) and the two-read floor: deliberate recall trades, not overfit rules.
- **Size:** the missed set is 2.3% / 4.1% / 5.1% of Rustle's correct chains on N and 19% / 13% / 6.6% on P. Rustle has more
  correct chains than StringTie on five of six (not chimpanzee) and the higher precision on fly, chimpanzee and bonobo.
  By this comparison the defaults do not favour the primates; the remaining recall gap is on the ape testis libraries.

**Wrong chains** (not registered as a reading; share of each tool's wrong chains that are NEAR a true chain, within 10 bp):
Rustle 25.7% / 18.7% / 4.0% / 3.1% / 2.6% / 5.1%; StringTie 15.9% / 17.6% / 0.7% / 0.6% / 0.5% / 1.0%. On the four real
libraries Rustle carries 4-7 times StringTie's number of near-junction wrong chains (mouse ES 1,321 against 197, KB3781 2,196
against 375). Not primate-specific. The opt-in junction ratio (rows 1819 / 1835) acts on this class.
On mouse ES 61% of Rustle's wrong chains (StringTie 55%) contain an annotated chain as a part: GENCODE lists many partial
transcripts, so much of this is the annotation.

Limits as stated; in addition the mouse ES truth is GENCODE and the ape truths RefSeq, which differ in how many partial
transcripts they list, so PIECE / LONGER shares are not comparable between those two without that caveat.

## Follow-up, registered 2026-11-20 before it was computed: is IsoQuant's precision an operating point on these six datasets?
User (2026-11-20): IsoQuant had very good metrics; is there anything to port? Row 1765 answered "operating point, nothing to
model" at junction level on OR6737 and A119b only. Same question at chain level on the six datasets of this row (`iq_point.py`).
- Rustle `x.final.gtf` multi-exon chains ranked by `reads` (ties: all taken); cut at the smallest read threshold whose kept set
  still has at least IsoQuant's number of correct chains. **OPERATING POINT** on a dataset when Rustle's precision at that cut
  is >= IsoQuant's - 0.02 (row 1765's bar). The same against StringTie, beside it.
- On a dataset that fails: Rustle's wrong chains at the cut by the classes of this row (PIECE, LONGER, NEAR, OTHER FORM, NEW
  INTRON) beside IsoQuant's, to name the class in which IsoQuant is cleaner. A class with at least half of the excess wrong
  chains is the pointer for a port; a port is a new pre-registered row.
- IsoQuant's rules read from its source (4.1.0, `isoquant_lib/model_construction`) are listed with their Rustle counterparts
  whatever the outcome.
