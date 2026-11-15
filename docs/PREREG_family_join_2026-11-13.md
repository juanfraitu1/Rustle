# PREREG row 1662: join MCL families through strong exon-sharing edges (2026-11-13)

Goal: more family members linked in the O1 families stage (edge recall), without paying more precision than lower
MCL inflation would.

## Step 0 (done before this registration; `bench/denovo_eval/edgemiss.py`, workspace `r1662_edgemiss/`)

First stage of the shipped edge rule that drops each expressed coding gene pair (Compara simiiformes truth, different
symbols), on the shipped baselines `port1593/{Ab,Hb}` (SD-damped, SHARE >= 0.5, I 2.0):

| stage | A119b truth missed | A119b false apart | HG002 truth missed | HG002 false apart |
|---|---|---|---|---|
| no record | 224 | - | 92 | - |
| short / low identity | 3 | 953 | 3 | 978 |
| no exon-to-exon base | 50 | 199,045 | 31 | 164,545 |
| SHARE < 0.5 | 202 | 2,473 | 107 | 1,734 |
| ... of which SHARE over the union of records >= 0.5 | 31 | 52 | 12 | 29 |
| coverage < 0.3 | 1 | 60 | 0 | 67 |
| passes every gate, not in one family | 198 | 346 | 84 | 266 |

Code facts: SHARE is the best single record (coverage and identity are pooled over a pair's records); the exon gate
does not align pairs with no chain holding exon bases of both loci, so "no record" includes "not attempted".
No gate has a workable price: SHARE < 0.5 holds 12 to 16 false pairs per true one, the union variant 2 per true one.
The pairs that pass every gate split by where their genes sit:

| passes every gate | A119b truth / false | HG002 truth / false |
|---|---|---|
| both genes in no family | 10 / 170 | 13 / 139 |
| one gene in no family | 12 / 90 | 11 / 85 |
| genes in two different families | 176 / 86 | 60 / 42 |
| ... and SHARE >= 0.75 | 87 / 9 | 20 / 3 |

Family pairs joined by such edges are nearly pure (A119b 49 all-true, 55 all-false, 6 mixed of 110).
Threshold 0.75 was chosen on this table: A119b and HG002 are development samples for this row.

## Rule (prototype `bench/denovo_eval/famjoin.py`, no binary change)

After MCL, two families become one when at least one edge of the shipped graph with SHARE >= 0.75 joins a locus of
one to a locus of the other (edges from `mcl_families --min-exon-share 0.75 --dump-graph`, every other flag as
shipped; joins are transitive). Loci in no family are not attached.

Not the closed levers: components over all edges, inflation, Louvain, two-level MCL, SD cluster merging (all on the
0.5 or coverage graph). This one uses only the strong-share subset between existing families.

## Scoring

`pairscore.py`, precision and recall only, PAIRSCORE_SAME_SYMBOL=ignore, default level. Comparator: the inflation
curve of the shipped rule, I in {1.4, 1.6, 1.8, 2.0}, linear interpolation in recall; when the arm's recall is above
the curve's range, the comparator is the curve's highest-recall point and the arm must have precision >= it.

1. Development go/no-go (A119b `a1.families.gtf`, HG002, SD-damped): on both, recall >= baseline + 0.01 and
   precision >= the curve at the arm's recall. Otherwise STOP, gorilla not run.
2. Held out: gorilla OR6737 and KB3781, truth with DENOVO_GGO_ENS_MAP (r1611 map), scored as row 1655's
   `score_ensmap`. That truth counts two gorilla genes mapping to one human gene as false; registered now: bars are
   judged on the figures as scored, and the same figures with those pairs set aside are reported beside them.
   PASS needs, on both samples: (a) recall >= baseline + 0.01; (b) precision >= the curve at the arm's recall;
   (c) largest family <= 2x the baseline's largest.
   Known risk: gorilla misses are mostly genes in no family (row 1655), which this rule does not touch; (a) may fail.
3. PASS -> report, binary port as an opt-in flag proposed to the user. FAIL -> register, nothing ported, no second
   threshold tried.

Descriptive readout: expressed truth genes in a family, before and after.

## Result (2026-11-13): FAIL at the development go/no-go; step 0 table corrected

| sample | arm | pairs | precision | recall |
|---|---|---|---|---|
| A119b | I 1.4 / 1.6 / 1.8 | 2,027 / 1,905 / 1,879 | 0.869 / 0.903 / 0.910 | 0.807 / 0.788 / 0.783 |
| A119b | shipped I 2.0 | 1,860 | 0.918 | 0.782 |
| A119b | join | 1,951 | 0.875 | 0.782 |
| HG002 | I 1.4 / 1.6 / 1.8 | 275 / 235 / 220 | 0.575 / 0.668 / 0.691 | 0.382 / 0.379 / 0.367 |
| HG002 | shipped I 2.0 | 209 | 0.703 | 0.355 |
| HG002 | join | 316 | 0.465 | 0.355 |

The join added no true pair on either sample (hits 1,707 and 147 unchanged) and 91 / 107 false ones. Families
1,395 -> 1,237 and 537 -> 482. Gorilla not run. Nothing ported; `famjoin.py` kept as the record.

**Why step 0 misled.** `edgemiss.py` called a truth pair "linked" from its own gene-to-family mapping
(`famtruth.load_catalog`), which is stricter than `pairscore.py`: about 200 A119b pairs that pairscore counts as hits
read as missed, and most of the "two different families" pool was those. With pairscore's hit list as the definition
(`--hits`, added), the missed truth pairs are:

| first stage that drops the pair | A119b (473 missed) | HG002 (267 missed) |
|---|---|---|
| no record | 224 | 92 |
| SHARE < 0.5 | 166 (28 with union SHARE >= 0.5) | 103 (12) |
| no exon-to-exon base | 36 | 31 |
| passes every gate (clustering) | 43 (8 with SHARE >= 0.75) | 38 (6) |
| short / low identity / coverage | 4 | 3 |

So clustering holds under 10% of the misses on A119b; the loss is before the graph: pairs with no alignment record
(about half) and pairs under the exon-share gate (about a third). The false-pair counts of the step 0 table stand.

## Step 0b, registered before running: are the no-record pairs lost to the all-vs-all settings?

The all-vs-all is `minimap2 -x asm20 -c -X -N 50 -p 0.1` behind the exon gate. A paralog whose chain scores under
0.1 of the query locus's best hit is never reported, and the gate skips pairs with no chain on exons of both loci.
Test on A119b (development): query = the loci of the genes in no-record missed truth pairs where both genes have a
locus; target = all loci (`Ab.fam.loci.fa`). Arms: S0 stock flags without the gate; S1 `-N 500 -p 0.01`;
S2 S1 plus `-k 15 -w 10`. Each PAF is read by `edgemiss.py --hits` (shipped rule: stage `edge`).
Stop rule: an arm goes on to a family-level test only if it brings >= 30% of those pairs to stage `edge` and the
unlinked false pairs it brings to `edge` (same query loci) are at most 1 per true pair. Otherwise aligner
sensitivity stays closed and the no-record class is recorded as below the DNA alignment floor.

### Step 0b result: STOP

153 of the 224 no-record pairs have a locus for both genes. Brought to stage `edge`: S0 0, S1 0, S2 1 (bar: 46).
Still without any record: 146 / 146 / 141. S1 gives the same 41,082 records as S0, so `-N 50 -p 0.1` and the exon
gate drop nothing here; shorter seeds (S2) add records but not for these pairs. The no-record class is below the DNA
alignment floor of locus spans; aligner sensitivity stays closed.
At least 35 of the 153 are GAGE x XAGE cross pairs (truth-set paralogs with no alignable DNA).

### Step 0c: the exon-share class re-priced with pairscore's hit list (post hoc, no bar)

Unlinked pairs under SHARE < 0.5, by SHARE over the union of a pair's records and pooled identity (true / false):

| condition | A119b | HG002 |
|---|---|---|
| all | 166 / 2,473 | 103 / 1,734 |
| union >= 0.5 | 28 / 52 | 12 / 29 |
| union >= 0.5, identity >= 0.85 | 8 / 10 | 0 / 6 |
| union >= 0.75 | 9 / 5 | 1 / 9 |

For comparison, lowering inflation from 2.0 to 1.6 adds 13 true / 32 false pairs on A119b and 10 / 16 on HG002.
No condition is clearly above that curve on both samples. No lever registered; stop and ask the user.

## Step 0d, registered before running: protein-level pairs genome-wide (user choice 2026-11-13)

Row 1590 kept protein pairs only between SD-linked loci. Here: every locus pair of row 1590's A119b DIAMOND
all-vs-all (`prot/allvall.tsv`; e <= 1e-10, both covers >= 50%, identity >= floor), mapped to expressed coding gene
pairs as `edgemiss.py` does, crossed with the shipped run's hits and stages. Development sample, existing output.
Stop rule: a family-level test is registered only if some identity floor in {50, 70, 90} holds, among gene pairs
not linked today, >= 40 true pairs at >= 1 true per false pair. Otherwise protein edges are closed genome-wide too.

### Step 0d result: STOP (and row 1589 had already answered it)

Row 1589 (2026-11-08) ran protein edges genome-wide on A119b and stopped at new-pair precision 0.035; I recommended
this step without finding that row. Run here on row 1590's locus-keyed output (`r1590dev/pe.allvall.tsv`) against
the exon-sharing default, gene pairs not linked today (true / false):

| identity floor | all | with a DNA record under SHARE < 0.5 | passing every DNA gate |
|---|---|---|---|
| 70 | 67 / 2,420 | 35 / 164 | 18 / 74 |
| 90 | 44 / 408 | 21 / 74 | 14 / 56 |

Bar was >= 40 true at >= 1 true per false. Not met at any floor, nor by the two conjunctions with the DNA stages
(post hoc). Protein edges are closed genome-wide as well as between SD-linked loci.

## Step 0e, registered before counting: absorbed genes (user choice 2026-11-13), is row 1654's A119b precision bar reachable?

Row 1654's F2 arm adds +82 true / +55 false pairs on A119b and +59 / +19 on HG002; it fails only A119b's
"precision >= baseline" bar (allows about 7 new false at 1,789 hits). `piece_attach.py` counts a sibling piece's
parent locus as a family mate; pieces attached only through that path (SMN x SERF1, NBPF x NOTCH2NL) could be
removed by an attach rule, pieces sitting with unsplit family members cannot. Check on row 1654's files: attach a
piece only when its family in the multi-representative run holds an unsplit member of the target family; score with
pairscore. Stop rule: go on to a registered arm only if A119b keeps <= 7 new false pairs and >= 60 of the 82 new
hits. Otherwise no further attach variant; the opt-in decision on F2 goes back to the user.

### Step 0e result: STOP, the A119b precision bar is not reachable by an attach rule

| arm (damped default) | A119b pairs added, true / false | A119b precision / recall | HG002 added | HG002 precision / recall |
|---|---|---|---|---|
| shipped | - | 0.918 / 0.782 | - | 0.703 / 0.355 |
| F2 (row 1654) | +82 / +55 | 0.896 / 0.820 | +59 / +19 | 0.718 / 0.498 |
| direct pieces only | +58 / +27 | 0.907 / 0.809 | +37 / +7 | 0.727 / 0.444 |

Bar: <= 7 new false and >= 60 new hits on A119b; got 27 and 58. The parent-ID path does carry 28 of the 55 false
pairs (SMN x SERF1, NXF2 x TCP11X, AGAP x PTPN20, most NBPF x NOTCH2NL) but also 24 hits. The 27 that remain sit
with unsplit family members: USP6 x TBC1D3 (9), GCC2 x RGPD (6), NBPF x NOTCH2NLA (6), PKD1 x NPIPA (2), CBWD (2),
TCP10L x TCP10L2, LIMS1. USP6 and the RGPD genes are chimeric genes that carry TBC1D3 and GCC2 sequence, and CBWD1/2/6
are paralogs: these are pairs the truth set scores false, not attach errors. No further attach variant.
For the user's decision: F2 passes row 1654's four bars on HG002, fails A119b on precision against baseline only,
and on both samples sits above the inflation curve (A119b curve tops out at 0.869 / 0.807; HG002 at 0.575 / 0.382).

## Truth audit of the false pairs (user choice 2026-11-13), rule fixed before counting

For every pair pairscore scores false in the shipped, F2 and direct arms (A119b, HG002): spliced exon sequences of
the two annotated genes, `minimap2 -x asm20 -c` all-vs-all (`bench/denovo_eval/falseaudit.py`). Classes, first that
applies: `related` (pairscore's own flag: Compara at any level, HGNC group); `shared` = records covering >= 0.5 of
the shorter gene at pooled identity >= 0.85 (seqscore.py's truth rule); `partial` = at least one record >= 300 bp at
>= 0.70 identity (the families stage's floor); `none`. Reported: precision as scored, with `related` + `shared` set
aside, and with `partial` also set aside. A report, not a bar; the F2 decision stays with the user.

### Truth audit result

False pairs by class (related / shared / partial / none) and precision as scored, with related + shared set aside,
and with partial also set aside:

| sample | arm | hits | false | related | shared | partial | none | as scored | w/o related, shared | w/o partial too |
|---|---|---|---|---|---|---|---|---|---|---|
| A119b | shipped | 1,707 | 153 | 98 | 1 | 15 | 39 | 0.918 | 0.969 | 0.978 |
| A119b | F2 | 1,789 | 208 | 114 | 3 | 26 | 65 | 0.896 | 0.952 | 0.965 |
| A119b | direct | 1,765 | 180 | 110 | 2 | 24 | 44 | 0.907 | 0.963 | 0.976 |
| HG002 | shipped | 147 | 62 | 44 | 0 | 5 | 13 | 0.703 | 0.891 | 0.919 |
| HG002 | F2 | 206 | 81 | 51 | 0 | 13 | 17 | 0.718 | 0.873 | 0.924 |
| HG002 | direct | 184 | 69 | 50 | 0 | 5 | 14 | 0.727 | 0.906 | 0.929 |

Two thirds of every arm's false pairs are older paralogs or same HGNC group (USP6 x TBC1D3 and CBWD fall here).
Pairs with no shared exon sequence added over the shipped run: F2 +26 on A119b (NBPF x NOTCH2NL, PTPN20 x AGAP,
NXF2 x TCP11X, ANXA8 x FAM25G: neighbours in one locus) and +4 on HG002; direct +5 (three NBPF x NOTCH2NLA,
PKD1 x NPIPA8, TCP10L x TCP10L2) and +1. The audit does not lift F2 over the shipped run on A119b under any reading;
the direct variant is within 0.002 to 0.006 of it on A119b and above it on HG002, for +58 and +37 true pairs.
`pairscore.py`'s false dump gained a fifth column (related flag) for this.

## Wired as opt-in (user decision 2026-11-13): driver stage `famreps`

`tools/rustle_pipeline.sh famreps --bam B --fasta G --out PREFIX`, after `families`: `locus_reps.py --disjoint`, the
families stage again on that input (`PREFIX.famreps.fam.*`), `piece_merge.py`, `piece_attach.py --keep-all --direct`
-> `PREFIX.fam.reps.copies.tsv`. `PREFIX.fam.*` untouched; never part of `all`; fails row 1654's A119b bar as scored.
End-to-end check on A119b (`r1662_edgemiss/e2e`, 465 s at 72 threads): output identical to the scored direct table,
default copy table byte-identical, pairscore 0.918 / 0.782 -> 0.907 / 0.809. Not run through the driver on HG002 or
gorilla; downstream stages (assign) do not read the new table.

### `famreps` through the driver on the other samples (2026-11-13, `r1662_edgemiss/e2e_{hg,or,kb}`)

| sample, truth | shipped precision / recall (hits) | famreps |
|---|---|---|
| HG002 | 0.703 / 0.355 (147) | 0.727 / 0.444 (184); identical to the scored direct table |
| OR6737, symbol truth | 0.514 / 0.205 (18) | unchanged |
| OR6737, protein-hit truth | 0.359 / 0.744 (753) | 0.359 / 0.759 (768; +15 true, +28 false) |
| KB3781, symbol truth | 0.632 / 0.176 (12) | 0.619 / 0.191 (13; +1 true, +1 false) |
| KB3781, protein-hit truth | 0.209 / 0.408 (97) | 0.223 / 0.450 (107; +10 true, +6 false) |

Gorilla equals row 1655's F2 figures: the direct filter removes nothing there. Wall 222 to 318 s at 24 threads.

## Assign regression on the table with extra representatives (A119b, 2026-11-13; opt-in measurement, read-assignment layer)

User asked for `famreps` as the default with `assign` reading it, after a regression run. `reps_table.py` writes the
complete copy table (tsv + fa + regions). Both runs start from the reference run `audit/g2/r76e` (`r1662_edgemiss/asg_*`).

| | shipped table (control) | table with extra representatives |
|---|---|---|
| copies | 6,445 | 6,639 |
| read x family rows | 200,817 | 208,413 |
| assigned rows | 4,815 (identical to the reference run) | 5,198 |
| assigned reads | 3,662 | 3,599 |
| table checks | all pass | 1 fails: 3 transcripts were copies of two families (fixed in `reps_table.py` after this run) |

392 of the 3,662 assigned reads lose their assignment and 329 reads gain one. Of the 475 lost rows, 254 have no row
at all in the new run and 221 became tied or ambiguous. One cause found: MCL215 grows from 8 to 26 copies and then
falls under the satellite exclusion rule (629 vs 628 families excluded), dropping its 129 assigned rows. MCL564 (SMN
region) goes from 4 to 6 copies, loses 141 rows and gains 591. A119b has no read-level truth, so the run cannot say
whether the changed assignments are better. Default NOT changed.

## Default since 2026-11-13 (user decision after the regression): `families` runs famreps, `assign` reads its table

`stage_families` runs `stage_famreps` unless `RUSTLE_FAMREPS=0`; `PREFIX.fam.copies.{tsv,fa,regions}` hold the extra
representatives, `PREFIX.fam.base.copies.*` the first run's table (byte-identical to the reference run's). The
satellite exclusion in `assign` is judged on the base table. Regression through that path (`r1662_edgemiss/asg_def`
vs the control `asg_base`, A119b, reference `audit/g2/r76e`):

| | before | default now |
|---|---|---|
| copies | 6,445 | 6,636 |
| families wall (72 threads) | about 410 s | 837 s |
| read x family rows | 200,817 | 207,377 |
| assigned rows | 4,815 | 5,399 |
| assigned reads | 3,662 | 3,797 (365 lose their assignment, 500 gain one) |
| families excluded as satellite | 628 | 628 |
| `check_run_tables.py` | all pass | all pass |

Of the 445 rows assigned only before, 164 became tied, 119 ambiguous and 162 have no row; MCL564 (SMN region,
4 -> 6 copies) loses 141 rows and gains 591, MCL215 (8 -> 26 copies) loses 103 and gains 184. A119b has no read
truth: the changed assignments are not graded. This is the read-assignment layer; not for `advisor_share`.
Driver tests (stub binary, 31 tests) pass but do not exercise famreps.

### Checks after the default flip (same day)

- Score of the tables the default path wrote (`asg_def`, pairscore, PAIRSCORE_SAME_SYMBOL=ignore): base table
  0.912 / 0.782 (1,710 true, 165 false), default table **0.902 / 0.808** (1,768 true, 192 false): +58 true, +27 false.
  The 0.918 -> 0.907 quoted above came from the older `port1593` table, which still scores 0.918 / 0.780 with this
  GTF: the pipeline's base table has drifted since (15 more pairs, 12 of them false; cause not looked for).
- The 162 assigned rows with no row in the new run are two families. The evidence below (positions, AS from the
  as_table, two reads viewed, the 50 kb pad) is consistent with one mechanism each; the gate in `copy_assign` that
  drops the row was not traced, and the 40 MCL215 rows had `in_copy=true` before, which is not explained. MCL215 (40 rows): the
  read's best-scoring alignment (its primary, unique best AS) now lies on a newly added copy, so the read is no longer
  tied. MCL564 (122 rows): a new copy at chr5:71,177,003-71,211,430 brings chr5:71.23 Mb, an off-catalog locus, into
  the 50 kb read window (`COPY_READ_PAD`); 106 reads have their primary there with a unique best AS (e.g. 1,578 vs
  1,558 on the copy they were assigned to), 10 have a better record on a base copy. Before, these reads were assigned
  on alignments other than their best; now they are not assigned. Skip lists are identical (628 families).
- Helpers moved to `tools/` (`locus_reps.py`, `piece_merge.py`, `piece_attach.py`, `reps_table.py`; symlinks left in
  `bench/denovo_eval/`), so the driver finds them beside itself.
- `families` and `all` now need python3 by default: both refuse at the start without it and name `RUSTLE_FAMREPS=0`.
- Driver tests: a stub that writes a whole copy table checks that the default leaves `fam.base.copies.*` and
  `RUSTLE_FAMREPS=0` does not; the python refusal is tested in the dispatch file.
- Other sessions: byte-for-byte comparisons of `PREFIX.fam.copies.*` and `PREFIX.assign.*` against runs made before
  this date will differ. Compare `PREFIX.fam.base.copies.*`, or run with `RUSTLE_FAMREPS=0`.
- KB3781 `all` from scratch with the driver run from a copy (`kb_all`): exit 0 in 642 s (families 294 s), table checks
  pass, a `families` rerun reproduces both tables. Scores base -> default: symbol truth 0.600 / 0.176 -> 0.591 / 0.191,
  protein-hit truth 0.205 / 0.408 -> 0.219 / 0.450. The HG002 and OR6737 figures above predate the duplicate fix in
  `reps_table.py`. The stage leaves about 1.2 GB of `PREFIX.famreps.*` files per gorilla run.
- `cargo test --release`: all 40 test binaries pass.
- Members, not pairs: genes placed in a family beside at least one true partner, base -> default: A119b 500 -> 506,
  KB3781 (protein-hit truth) 97 -> 102. Most of the +58 A119b pairs join genes that already had a partner.

## Step 0f (user choice 2026-11-13, registered before the count): no-record pairs by protein hit AND proximity

Step 0d failed on protein hits alone (67 true / 2,420 false at identity >= 70). Not tried: the conjunction with
position. Many no-record truth pairs are old tandem duplicates (GAGE x XAGE), close on one chromosome, where chance
protein hits between unrelated genes are rare. Rule counted: an unlinked gene pair gets an edge when a protein hit
(row 1590's DIAMOND output, identity >= 50) joins loci on the same chromosome within D, D in {0.2, 1, 5} Mb, read at
identity floors 50 / 70. A119b only, `port1593` tables (as step 0d). No pipeline change.
**Stop rule (frozen):** go on only if some (D, floor) cell gives >= 30 new true pairs at new-pair precision >= 0.75
(false counted as pairscore counts them). Otherwise STOP and the no-record class stays closed to sequence evidence.

**Result: STOP.** Only 73 of the missed truth pairs have a same-chromosome protein hit at all. Best cells: 51 new true /
689 false (5 Mb, identity >= 50; precision 0.069), 39 / 299 (5 Mb, >= 70; 0.115), 24 / 441 (0.2 Mb, >= 50). No cell
reaches precision 0.75. The false pairs are mostly real but older paralogs (ATAD3A x ATAD3B x ATAD3C, FBXO2 x FBXO6 x
FBXO44, AKR7A2 x AKR7A3), which the simiiformes-level truth scores false: proximity plus protein similarity finds
old tandem families, not the young ones the truth asks for. No-record class closed to sequence evidence.
