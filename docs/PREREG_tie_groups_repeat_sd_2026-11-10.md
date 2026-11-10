# PREREG r1607: tie groups — interspersed repeats set aside, co-duplicated loci as their own category (2026-11-10)

Follows r1606 (`docs/PREREG_tie_groups_floor_2026-11-10.md`, FAIL). Post hoc there: the largest group is a transcribed
interspersed repeat; other impure groups hold different genes of one segmental duplication. User (2026-11-10): download the
T2T RepeatMasker track; segmental duplications are already in use (SEDEF). Written and committed before any arm runs.

## Inputs
- RepeatMasker: `chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed` (T2T public bucket, downloaded 2026-11-10 to
  `/anvme/workspace/v138da13-rustle/annot/dl_rmsk/`; 5,590,282 rows).
- SEDEF pairs: `~/data/chm13/chm13v2.0.fa.sedef.bed` (the pipeline's `RUSTLE_SEDEF_BED` sidecar).

## Rules
- **Repeat locus** (`tie_groups.py build PREFIX as98 0 1 RMSK_BED:0.5`): a locus whose model has >= 50% of its exonic bases
  inside LINE, SINE, LTR, DNA, Retroposon or RC elements (merged). It is written to `PREFIX.repeat.tsv`, is not emitted as
  a model and takes no part in linking. Simple repeats, low complexity and satellites do not count.
- **Co-duplicated locus** (scorer only): a foreign locus that overlaps one side of a SEDEF pair when a family locus of the
  same group overlaps that side or the other side. A group whose foreign loci are all co-duplicated is counted with the
  pure groups in B2s. The strict B2f figure is still printed.

## Arms
Repeat filter off (control = r1605/r1606 baseline) and on. Ties as98, link fraction 0, floor 1, scorer `full`.

## Bars (the filter-on arm; M = the 72 indistinguishable members)
- B1g grouped reach >= 58 of 72; B3 TAF11L models <= 5 kb reach >= 20 of 31; B5 largest group <= 60 loci; B4 job <= 1 h.
- B2s purity, co-duplicated-only groups counted as pure, >= 0.90.
- PASS = all five. The control must reproduce 70 / 0.859 (B2f) / 31 / 831.
- Descriptive: loci set aside as repeat (count, molecules), B2f strict, the remaining impure groups by name.

## Outcome
(empty until the run)
