# PREREG r1829: how many near-tied reads exist, and what could separate them (2026-11-19)

Written on PSU (Roar) while Fritz is in maintenance; workspace `~/ws/r1829_near_tie_census` (claimed on PSU after the
register tail on origin/main was 1828; the Fritz handoff guide's support-filter task now takes the next free number).
Design accepted by the user on 2026-11-19: layer 2 on HG002 only; the span class in scope; row claimed on PSU. Descriptive
census, no rule is changed by it. Nothing of this row was counted before this document was committed; changes after the
commit are amendments, dated, in a section of their own.

## Why
User, 2026-11-19: "can we make a study about how many near tied multimapping reads actually exist? ... we know what
happens if we use most of the secondaries for assembly and it's not all positive, but can we also have an idea of what is
the pool of undecidable, separable via PSV or via any other means?"

What bears on it (not re-run):
- Near ties are 5.7% of FAMILY reads on HG002 / CHM13, 1.3% on KB3781 on its own assembly (1786), 0.6% on the Sumatran
  orangutan (1787); half of HG002's sit in 4 families, 14 families hold 81% (1790).
- Of the true ties in those 14 families 74.9% fit several HG002 copies equally (SEVERAL COPIES), 69.4% over all (1791);
  of all near ties there: 44% several copies, 15% ties the truth cannot judge, 34% no home, 6% decidable.
- Own-read copy differences (r1783) certify 427 of 710 decidable near ties; reference lead of 1-5 is wrong 41.6% of the
  time on the best-AS locus (1786).
- Seeding with the good pool (AS >= 0.98 x best) helps on held-out gorilla (1060 / 1100) but adds shadow loci and
  wrong-copy chains (1767, 1820, 1821, 1825).
Missing: a census over ALL reads (not only family reads) on every sample, with the three tiers side by side, and for the
near ties a split by WHAT could separate them.

## Inputs
- Layer 1, one library per sample, the `as_table` of each BAM (`BAM.as_table.tsv`, checked against the BAM's size in its
  header; rebuilt if it does not match): HSA A119b testis and HG002 Kinnex (both on CHM13 v2.0); GGO OR6737 testis
  (`GGO_mm.bam`) and KB3781 fibroblast (gorilla reference); PTR testis / AG18354 LCL; PAB testis / AG06213 fibroblast;
  PPA testis / PR00251 fibroblast (`PPA_PR00251_fib_mm.bam`); PPY testis / AG05252 fibroblast; each ape on its own
  species reference under `~/genomes`. All BAMs: `minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`
  (checked in the `@PG` lines of the human and gorilla BAMs; the rest come from the same `map_reads.sbatch`). So at most
  50 secondaries per read: arrays with more copies are cut short, and a rival with ms < 0.1 x best is never recorded.
- CHM13 location tracks: SEDEF `chm13v2.0.fa.sedef.bed`, the exclusion list `chm13v2.0.fa.exclude.bed` (satellite /
  rDNA), the truth families of `bench/denovo_eval/famtruth.py` (genes CAT `$DENOVO_CAT_TSV`), CAT biotypes. No SEDEF or
  exclusion track exists for the ape references on PSU: their location split is not made.
- Layer 2: HG002, hash residues 0-24 (`neutral_truth.selected`, offset 0), the SEEN quarter of rows 1776-1791. All four
  HG002 quarters are spent; reuse is allowed (user ruling) and said here. This quarter is chosen because rows 1786 / 1790
  / 1791 counted it, which gives the calibration checks below. Truth: `bench/o2_truth/neutral_truth.py score` with
  `NT_REPAIR=1` against HG002 v1.1 mat / pat (`~/ws/r1776_hg002_room/{mat,pat}.mmi`), then `tie_split.py`.

## Definitions, fixed now
- **Molecule:** one `as_table` row (mapped, non-supplementary records of one read name).
- **Places:** the primary and the rival (`as_table`'s best-AS non-primary record). The rival is the SAME place when it is
  on the primary's contig with `|rival_pos - primary_pos| < 5,000` (counted and reported, not a near tie).
- **Tiers** (a molecule with a rival at another place):
  - `exact`: rival ms = primary ms (and the rival is not worse on AS by more than 0);
  - `ms10`: |best ms - rival ms| <= 10 (rows 1786 / 1790's near tie; includes `exact`);
  - `good`: rival AS >= 0.98 x best AS (the default seeding pool; reported whether or not within ms10).
  Shares are of all molecules with a primary (one-place molecules included in the denominator).
- **Location (CHM13 samples only), first that applies, both places considered:** `excluded` (either place in the
  exclusion list) / `family` (the primary on a truth-family member, as `neutral_truth.select` defines it) / `SD` (both
  places in SEDEF intervals) / `SD one side` / `other`.
- **Structure tag (all samples, truth-free):** `primary_spliced != rival_spliced` (one place spliced, the other not).
- **Layer 2 selection:** HG002 residues 0-24, molecules in tier `ms10`; at most 30,000 by hash order (`NT_CAP`); their
  records extracted from `HG002_kinnex.t2t.bam` (all records of the read names, nothing re-aligned) and scored by
  `neutral_truth.py score` (truth status per read and its own S0 / S1 / S2 tier; the agreement of its S2 with `ms10` is
  printed).
- **Layer 2 classes, a cascade, each read in the first class that applies:**
  1. `truth limit`: status `outside`, `no near-perfect hit`, `competitor unlocatable`, or a tie split ONE PLACE / UNSPLIT
     (the truth cannot say where the read came from). Reported FIRST in every table.
  2. `undecidable`: status `tie`, split SEVERAL COPIES (the read fits two or more places of one HG002 haplotype equally).
     Split again by span (in scope by the user's choice): `undecidable, copies differ elsewhere` when the CHM13 genes at
     the read's two best places (CAT, same strand, largest exonic overlap) have longest transcripts whose spliced
     sequences differ (any base or length; a longer molecule could decide); `undecidable, identical transcripts` when
     they are equal; `undecidable, no gene` when a place has no gene. A CHM13 proxy: HG002's copies may differ from
     CHM13's, said in the result.
  3. `reference lead`: decidable, and the place with the higher ms (lead >= 1) is the origin.
  4. `own variant`: decidable, not 3, and the r1783 rule (`own_psv.py` as committed) calls the origin. If `own_psv.py`
     cannot run on this selection with input plumbing alone, this class is reported as NOT MEASURED and its reads fall
     to 5 / 6.
  5. `structure`: decidable, not 3-4, the origin's record and the competitor's differ in spliced / unspliced or in the
     number of introns, and preferring the record with more introns picks the origin.
  6. `decidable, no signal`: decidable and none of 3-5.
  Beside (not in the cascade): decidable reads whose reference lead points to the WRONG place (the reference misleads),
  and structure differences that would pick the wrong place.

## Calibration checks (must hold before layer-2 numbers are read as this row's)
- K1: the `ms10` share of HG002 FAMILY molecules (residues 0-24) from `as_table` lies within 4.2-7.2% (row 1786: 5.7%).
  If not, the as_table definition is reported as not comparable to 1786 and the layer-1 numbers are given with that note.
- K2: on the family reads of the layer-2 selection, SEVERAL COPIES / (SEVERAL COPIES + ONE PLACE) lies within 60-80%
  (row 1791: 69.4% over all). If not, the cascade is reported but marked as not reproducing 1791.

## Predictions (mine, before any count)
- P1: `ms10` share of all molecules: A119b and HG002 between 1% and 4%; every ape sample on its own reference below 2%.
- P2: testis above the same species' cell line in at least 4 of the 5 species pairs.
- P3: `good` holds at least twice as many molecules as `ms10` in every sample.
- P4: on CHM13, at least 60% of `ms10` molecules have both places in SDs (`SD` + `excluded` with both in SD counted
  apart); `family` holds at most half of them (families are a subset of the duplicated genome).
- P5 (HG002 layer 2): `truth limit` at least 30% of the selection; among the rest, `undecidable` at least 50%.
- P6: of `undecidable`, at least 30% `copies differ elsewhere` (part of the pool is a read-length limit).
- P7: of the decidable near ties, classes 3-5 together place at least 70%; `reference lead` alone at least 40%.
Each is reported HOLDS / FAILS with its number.

## Output
`~/ws/r1829_near_tie_census/{census.tsv, location.tsv, layer2/classes.tsv, summary.txt}`, scripts
`bench/o2_truth/near_tie_census.py` and `bench/o2_truth/near_tie_classes.py`, a result section below, one register row.

## Time to answer (estimate)
Layer 1: minutes per sample (tables exist). Layer 2: extraction + two haplotype alignments of <= 30,000 reads + tie split,
under an hour on one node.
