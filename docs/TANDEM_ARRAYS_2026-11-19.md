# Tandem arrays: how StringTie and Rustle assemble them, read artifacts, and the size of a member (2026-11-19, PSU, register 1826)

User's goal: explore how StringTie assembles transcripts in tandem arrays; find whether Iso-Seq reads that map to tandem arrays
can be told from mapping artifacts; derive a rule that gives the size of each member. Descriptive, no bar. Data: A119b testis,
sliced at the 227 arrays of row 1820 (`r1820_undecidable_copies/arrays`: P.bam = primaries, asm_P / asm_D = Rustle primaries /
default pool, StringTie `-L` on the primaries), and the famsim arrays of `MULTICOPY_FAMILIES_2026-11-19.md` (truth known).
Code `bench/famwise_eval/tandem.py` (`tandem_arrays`, `fold_member`, `classify_model`); workspace `goal_multicopy_dev/tandem`
(linked from `r1826_tandem_arrays`).

## Regular arrays
An array of row 1820 is kept as TANDEM when it has >= 4 units, >= 80% on one strand and an even spacing of unit starts
(coefficient of variation <= 0.35); period P = median spacing. 25 arrays pass in CHM13 (TSPY fails: other genes sit inside its
span; its reads appear inside FAM197Y's). The filter still admits gene clusters with long spacing (ABCA, FMO, RARRES2P 70 kb)
and the protocadherins (overlapping genes: the control for real multi-unit splicing).

## 1. Read artifacts
- **Period jumps are rare.** A read of one unit aligned across two (an "intron" of k x P with donor and acceptor at the same
  offset inside their units) occurs in 5 of 15,115 spliced primaries.
- **What looks like an artifact is mostly real transcription the annotation does not describe:**
  - USP17L: all 120 unspliced reads are ANTISENSE to the genes, 1.8-4.4 kb long against a 1.6 kb gene, within the 4.7 kb
    repeat unit (antisense transcription of the repeat). This is why no assembler gives USP17L a model.
  - RARRES2P / AGGF1P (one region on chr16), GAPDHP, OR5H: ~2 kb reads over 0.5-1 kb pseudogene pieces: the transcribed unit
    is bigger than the annotated gene.
  - FAM197Y's array carries TSPY reads on the other strand, folding to exactly TSPY's 2.8 kb.
- The annotated gene is therefore the wrong yardstick for a member; the reads have to give it.

## 2. A rule for the size of a member: fold the reads onto the period
Each primary read is placed on the unit whose start is nearest at or before it; its start and end become offsets inside that
unit. Per strand, the member is the contiguous run of folded positions covered by >= 20% of the peak folded coverage; reads
whose folded span exceeds one period are MULTI-UNIT and set apart.
- Simulations (10 TSPY3 units, identical to 0.5% diverged, 1 to 10 expressed): folded size 2,800 bp in all six (true 2,797);
  0 multi-unit reads.
- A119b, short-period arrays: CT47A 3,325 bp (annotated 3,323), SPATA31D 6,475 (6,455), ELOA 3,125 (annotated 1,641: a
  longer transcribed unit), USP17L 3,650 (antisense), TSPY inside FAM197Y 2,800.
- Not usable where one period holds several genes (long-period clusters fold into one merged interval), and multi-unit reads
  dominate where splicing is really multi-unit (protocadherins: 940-1,140 of ~1,700 reads).

## 3. How the assemblers build models in arrays (folded member as the yardstick)
| array | StringTie | Rustle primaries | Rustle default pool |
| CT47A (member 3.3 kb) | 9 models, 6 multi-unit, median 1.6 x member | 7, 0 multi-unit, 1.00 x | 7, 0, 1.11 x |
| DEFA | 6, 1.54 x | 2, 1.04 x | 4, 1.34 x |
| USP17L (3.7 kb) | 2 multi-unit, ~30 kb | none | one 265 kb blob (row 1518's artifact) |
| FAM197Y / TSPY | 1, 2.1 x | 1, 2.25 x | 104 models at 3.0 x |
| SPATA31D | 35, 17 multi-unit, 6 x | 36, 23, 7 x | 36, 23, 7 x |
StringTie joins neighbouring units more readily than Rustle's primaries-only assembly; Rustle's default pool (secondaries
seeding) makes array-wide blobs and many over-long models.

## 4. Artifact signatures of a model (`classify_model`)
- **jump:** an intron of ~k x P (2%), donor and acceptor at the same unit offset, non-canonical motif. Example: StringTie's
  USP17L model, one 42.7 kb intron (~9 P), AG-CG.
- **repeat:** the same intron (folded) in >= 2 units: members concatenated. Example: StringTie's CT47A model, introns
  1,138 / 887 / 1,138 / 887 (canonical), two members joined.
- **bridge:** an exon covering a whole inter-unit gap.
Counts (multi-unit models flagged by any signature): StringTie 7 of 50 in regular arrays, Rustle primaries 0 of 33, Rustle
default pool 1 of 70; protocadherin control 0-1 of 111-128 (real multi-unit splicing is left alone). Most multi-unit models
carry no signature (e.g. SPATA31D's canonical, varied introns) and cannot be judged without truth.

## Proposed rule (not built into the assembler, not judged)
1. Per regular array, fold the primaries onto the period per strand -> member interval and size.
2. Split a model at the member boundary when it carries a jump, repeat or bridge signature, or spans > 1.5 x the member
   without a canonical junction crossing the boundary.
3. Do not seed arrays from secondaries (the default pool's 265 kb USP17L blob and FAM197Y's 104 over-long models); use
   `RUSTLE_SEED_SPLICED=support` (row 1821) and leave identical units to the copy-set call.
4. Report antisense and annotation-free transcribed units as such: they are not artifacts.

## What would judge it
No origin truth exists for expressed testis arrays (A119b has no assembled haplotypes; HG002 is a cell line). Options: a testis
Iso-Seq sample of an individual with assembled haplotypes; or a famsim array with the observed read classes added (antisense
reads, reads longer than the gene, period jumps at the observed rate) so that the split rule can be scored against truth.
