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

## 5. Why StringTie joins units: a causal test on A119b
Removing only the 1,488 primaries longer than one period inside the regular arrays (0.3% of the slice; protocadherins left
untouched) and rerunning StringTie: multi-unit models in regular arrays 50 -> 4 (CT47A 6 -> 1, longest 13.5 -> 4.6 kb;
USP17L 2 -> 0; SPATA31D 17 -> 1); protocadherin control 114 -> 114. StringTie chains neighbouring units through the few reads
that physically span from one unit into the next (read-through and run-on molecules); Rustle's primaries-only assembly does not.

## 6. The member-size rule, validated on simulated arrays with the real read classes
`tandem.member_sizes`: reads on the array's strand (the strand of its spliced reads), assigned to the unit nearest at or before
their start; antisense reads and reads longer than one period set apart; per unit with >= 5 reads, the member = from the mode
of the reads' 5' ends to the mode of their 3' ends (25 bp bins, inside the unit's window); fewer reads: the folded interval.
Simulation: 10 TSPY3 units 20 kb apart, 5 expressed with different members (full, 5' exon truncated, middle exon deleted, 3'
exon truncated), plus 100 antisense, 50 run-on and 40 read-through molecules (`arrays_sim/inject.py`, `score_rich.py`):
| estimator | 0.3% diverged units: members within 5% | silent units given a member | identical units: within 5% | silent |
| rule (read-end modes) | **5 / 5** (median error 0.4%) | 0 / 5 | 5 / 5 | 5 / 5 |
| rule v1 (coverage run at 20% of peak) | 4 / 5 (run-on reads inflate one member) | 1 / 5 | 1 / 5 | 5 / 5 |
| StringTie | 0 / 5 (models 24-47 kb) | 3 / 5 | 0 / 5 | 5 / 5 |
| Rustle primaries | 5 / 5 | 0 / 5 | 5 / 5 | 5 / 5 |
Identical units scatter their reads, so every method puts a member on silent units: that case stays with the copy-set call.
Real data (no truth): CT47A, rule 3,325 bp per member (5 members) vs annotation 3,323 and Rustle 3,348; StringTie 13,541.

## 7. The rule on REAL arrays, against the curated annotation (protein-coding members)
No haplotype truth exists for expressed testis arrays; the closest independent yardstick is the annotated full-length span of
PROTEIN-CODING members (real gene models, unlike the pseudogene pieces of section 1). A119b, every array of row 1820 with >= 2
protein-coding units, members with >= 5 sense reads; StringTie and Rustle take the same-strand model with the largest overlap
with the member (`goal_multicopy_dev/tandem/eval_real*.py`).
First version (end modes alone): tandem families 48 / 77 within 20%; its errors were members with 5-10 reads sized 9-22 kb (an
end mode taken from one run-on read) or shorter partial fragments. Two safeguards were then added, chosen AFTER seeing these
real errors (so the numbers below are development, not an independent test): each end must be shared by >= 2 reads, and a member
cannot exceed one period; otherwise the unit takes the folded array consensus. Simulations unchanged (5 / 5).
| | tandem families (CT45A, CT47A, GAGE, MAGEA, TSPY, AMY, SSX, PRAMEF, SPATA31D, DEFA, TBC1D3): within 20% of annotation | all 919 protein-coding cluster members |
| rule | **73 / 100** (17 longer, 10 shorter) | 285 / 919 |
| StringTie | 44 / 100 | 240 / 919 |
| Rustle primaries | 44 / 100 | 237 / 919 |
Examples: CT47A 3,325 per member (annotated 3,323; StringTie up to 13,541); TSPY members 2,775-2,800 (2,797); MAGEA3 3,600 (3,570);
GAGE12 members 7,325 (7,290). Disagreements left: PRAMEF (annotated 3-9 kb, reads give 2-25 kb), members whose reads cover only a
fragment, MAGEA11 (annotated 29 kb, reads 5.2 kb: a long annotated isoform). On diverged gene clusters (ZNF, KRT, olfactory
receptors) every method agrees with the annotation for about a third of members: member size is not the problem there.

## 8. Held-out tests of the frozen rule (OR6737 gorilla testis; predictions in `PREDICTION_tandem_heldout_2026-11-19.md`)
Both written before running, both FAILED:
| test | members | rule | StringTie | Rustle primaries |
| 1: arrays by name stems, >= 3 protein-coding units, length CV <= 0.2 | 20 (only large diverged clusters: ABCA, ADH, NLRP, CCL, ZNF) | 0.60 | 0.80 | 0.75 |
| 2: arrays by NCBI description, same + median unit <= 15 kb | 44 in 20 arrays | 0.55 | 0.55 | 0.66 |
| A119b, test 2's generic definition (development) | 131 in 30 arrays | 0.47 | 0.43 (0.45 without reads > 1 period) | 0.34 |
(share of members within 20% of the longest annotated transcript span)
Reading: the read-end rule is exact on simulated arrays and right on near-identical testis arrays (CT47A, TSPY in both species,
GAGE12, MAGEA3), but it is not a general member-size rule: on long genes its read-end modes size members too short (reads cover
fragments), and where reads run past a gene it doubles the member. Across real samples no method agrees with the annotation for
more than about two thirds of members, and which method is best changes with the sample; the annotation (longest isoform) is
itself an imperfect yardstick. What holds on every data set: StringTie joins neighbouring units through reads longer than one
period (section 5), and Rustle's primaries-only assembly does not.

## The rule as it stands (exact on simulation; not a general rule on real arrays, section 8)
1. Per array, per member: the read-end modes of `member_sizes` (section 6), each end shared by >= 2 reads and the member no longer
   than one period; otherwise the folded array consensus.
2. Split a model at the member boundary when it carries a jump, repeat or bridge signature, or spans > 1.5 x the member
   without a canonical junction crossing the boundary.
3. Do not seed arrays from secondaries (the default pool's 265 kb USP17L blob and FAM197Y's 104 over-long models); use
   `RUSTLE_SEED_SPLICED=support` (row 1821) and leave identical units to the copy-set call.
4. Report antisense and annotation-free transcribed units as such: they are not artifacts.

## What would judge it on real data
No origin truth exists for expressed testis arrays (A119b has no assembled haplotypes; HG002 is a cell line). Options: a testis
Iso-Seq sample of an individual with assembled haplotypes; or a famsim array with the observed read classes added (antisense
reads, reads longer than the gene, period jumps at the observed rate) so that the split rule can be scored against truth.
