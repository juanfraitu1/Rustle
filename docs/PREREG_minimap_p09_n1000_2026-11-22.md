# PREREG: minimap2 with -p 0.9 -N 1000 instead of -p 0.1 -N 50, does it recover more family members (row 1885, 2026-11-22)

User (2026-11-22): "can we run minimap again for a couple of small samples in which we can use the same parameters as
usual except -p 0.9 to keep only the near tied multimappers and a very high number of possible multimappers in that stage
like -N 1000 just to see if that would recover all members?"
Written before any realignment. Seen: row 1871 (2.4% of molecules at the cap of 51 records are all near-tied, none on
a named copy); row 1881 (joining secondaries score below 0.98 x best). Workspace `r1885_minimap_p09_n1000`.

## Arms (everything else as usual: `minimap2 -ax splice:hq -uf --eqx -Y --secondary=yes`)
- USUAL: `-N 50 -p 0.1`.  NEW: `-N 1000 -p 0.9`.
Both arms are realigned here from the same reads, so the comparison is like for like.

## Small samples
- Planted family reads of row 1874, scenarios E (all 133 copies expressed) and S (63 silent): the truth is known.
- A119b (human testis): the 46,825 molecules with a spliced record on a named-family copy (`mode.a1.tsv.reads`), realigned
  to the whole genome.
- OR6737 (gorilla testis): the 5,318 such molecules (`mode.or.tsv.reads`).
Limit of the two real subsets: a molecule with no record on any named copy under the usual parameters is not in them.

## Measures
- REACHED: copies with >= 2 molecules that have a near-tie record there (AS >= 0.98 x the molecule's best in that BAM;
  spliced records for the real samples). Planted S: silent copies reached (a cost, not a gain).
- Records: secondaries written, molecules with > 51 records, largest number of records of one molecule.
- MODELED (real subsets): `assemble --seed-pool good` + stage `multicopy` on each arm's BAM, named copies with a spliced
  transcript / junction precision (`toolcompare.py`).
- "Recovers more members" = NEW reaches or models at least 3 more copies than USUAL in a real subset, or at least 3 more
  expressed copies in planted E. Fewer by 3 or more = NEW loses members. Otherwise NO DIFFERENCE.

## Predictions (mine)
NO DIFFERENCE in copies reached and modeled on all samples; planted E: all 133 reached in both arms. NEW writes far
fewer secondaries (under half). Risk I see: -p compares chaining scores, not AS, so 0.9 may drop a record whose AS is
within 2%; if so NEW loses a few copies.

## Result (2026-11-22, job 4494750): NO DIFFERENCE in members; NEW writes fewer records
| sample | arm | secondaries | molecules at 51 / above 51 | largest | copies reached | modeled, mode (copies / junction precision) |
|---|---|---|---|---|---|---|
| planted E | usual | 139,393 | 65 / 0 | 51 | 133 of 133 | - |
| planted E | new | 111,896 | 0 / 0 | 46 | 133 of 133 | - |
| planted S | usual | 69,967 | 5 / 0 | 51 | 70 of 70 expressed, 41 of 63 silent | - |
| planted S | new | 57,197 | 0 / 0 | 46 | 70 of 70, 41 of 63 | - |
| A119b family reads | usual | 733,741 | 648 / 0 | 51 | 127 of 133 | 126 / 0.48 (secondaries only: 123 / 0.49) |
| A119b family reads | new | 120,004 | 0 / 15 | 225 | 127 of 133 | 126 / 0.48 (123 / 0.49) |
| OR6737 family reads | usual | 54,214 | 56 / 0 | 51 | 38 of 49 | 34 / 0.66 (27 / 0.58) |
| OR6737 family reads | new | 29,962 | 0 / 14 | 301 | 38 of 49 | 35 / 0.66 (27 / 0.58) |
Copies reached are the same copies per family in both arms. One more copy modeled on OR6737 (under the bar of 3).
A119b: the BAM of these reads is 1.26 GB against 0.27 GB, alignment 297 s against 141 s. The 15 and 14 molecules
above 51 records show the cap of 50 does cut real near ties for a few molecules, without changing a copy.
Deviation: the gorilla measures were computed after the job (the job lacked `DENOVO_GGO_GFF`); same scripts, same BAMs.
Predictions: no difference right; "under half the secondaries" right on A119b (16%), wrong on OR6737 (55%) and the
planted reads (80%); the feared loss from -p on chaining scores did not happen.
Not tested: a whole sample (reads with no record on a named copy under the usual parameters), genome-wide family
scores, other tissues. -p 0.9 removes the records below the near-tie band, so nothing that needs them (the row 1881
kind of census, tools reading all secondaries) can run on such a BAM.

# PREREG: the limit of useful multimappers, a whole sample realigned at four values of -p (row 1888, 2026-11-22)

User (2026-11-22): "yes lets do that, maybe with this we know what is the limit of useful multimappers, I originally
used -p 0.1 which was way too much". Written before any whole-sample realignment. Seen: row 1885 (family-read subsets:
no difference at -p 0.9). Workspace `r1888_multimapper_limit`.

## Arms, gorilla OR6737 testis, all reads (`minimap2 -ax splice:hq -uf --eqx -Y --secondary=yes`)
- USUAL: the existing BAM (`-N 50 -p 0.1`), reassembled here with today's code.
- `-N 1000` with `-p 0.5`, `-p 0.8`, `-p 0.9`, `-p 0.95`: realigned from the primaries' sequences of the existing BAM.
Each arm: `assemble --seed-pool good`, stage `families` (RUSTLE_FAMREPS=0, as row 1858), stage `multicopy`.

## Measures
- Records: secondaries, near-tie secondaries (AS >= 0.98 x the molecule's best: the only ones the mode uses), molecules
  above 51 records, BAM size, alignment time. In the USUAL BAM, the spread of secondaries by AS / best.
- RETAINED: share of USUAL's near-tie secondaries (molecule, contig, start, strand) present in the arm; and the arm's
  near-tie secondaries that USUAL lacks (what -N 50 cut).
- Named copies modeled / junction precision (`GGO_MAINSTAYS`), assembly with secondaries and the mode.
- Family pair precision / recall, level simiiformes (as row 1858).
- Transcripts, and intron chains shared with USUAL's assembly.

## Reading
An arm is SUFFICIENT when it retains >= 99% of USUAL's near-tie secondaries, models within 1 named copy of USUAL (mode),
and has pair recall and pair precision within 0.01 of USUAL (not lower by more). The limit = the highest -p that is
sufficient. Nothing here changes the mapping recipe: that is the user's call, and one sample.

## Predictions (mine)
0.5, 0.8 and 0.9 sufficient. 0.95 not (it drops more than 1% of the near-tie secondaries: -p acts on chaining scores).
Under 20% of USUAL's secondaries are near-tie. -N 1000 adds near-tie secondaries for under 1% of molecules.

## Result of row 1888 (2026-11-22, jobs 4494838-4494842, 4494856): the limit by the registered rule is -p 0.5
OR6737, 4,404,440 molecules. USUAL = `-N 50 -p 0.1`; the others `-N 1000`.
| arm | BAM GB | secondaries | near-tie secondaries | USUAL's near-tie retained | molecules > 51 records (largest) | named copies, mode | pair precision / recall | transcripts | intron chains shared with USUAL |
|---|---|---|---|---|---|---|---|---|---|
| USUAL | 11.7 | 6,266,238 | 110,012 | - | 0 (51) | 34 / 0.66 | 0.090 / 0.394 | 70,181 | 69,647 |
| -p 0.5 | 7.0 | 1,698,759 | 116,339 | 99.7% | 1,050 (942) | 34 / 0.66 | 0.090 / 0.394 | 70,242 | 69,643 of 69,646 |
| -p 0.8 | 5.9 | 520,581 | 112,385 | 98.3% | 428 (889) | 35 / 0.66 | 0.090 / 0.394 | 70,196 | 69,641 of 69,642 |
| -p 0.9 | 5.7 | 292,547 | 105,211 | 92.7% | 286 (880) | 35 / 0.66 | 0.090 / 0.394 | 70,114 | 69,629 of 69,632 |
| -p 0.95 | 5.6 | 180,167 | 87,663 | 76.7% | 227 (857) | 35 / 0.66 | 0.090 / 0.390 | 70,048 | 69,601 of 69,611 |
- Of USUAL's secondaries 1.8% are near-tie (110,012); 70% score under half of the best.
- SUFFICIENT by the rule (>= 99% retained + copies + pairs): -p 0.5 only. -p 0.8 misses the retention bar by 0.7 points
  with every downstream measure equal to USUAL's; -p 0.9 loses 7.3% of the near-tie secondaries, still with equal
  downstream measures; -p 0.95 loses 23% and the first family pairs (190 hits against 192).
- -N 1000: 6,625 near-tie secondaries at -p 0.5 that USUAL lacks; 1,050 molecules (0.02%) exceed 51 records.
- Alignment 21-24 min per arm on 56 threads.
Predictions: 0.5 sufficient right; 0.8 and 0.9 sufficient wrong (retention); 0.95 not sufficient right; "under 20%
near-tie" right (1.8%); "-N adds near ties for under 1% of molecules" right.
Deviations: the assemblies ran on the driver's own derived BAM (`g.shortexons.d/reads.bam`, a default from another
session, the same in all arms); the record counts were taken on each arm's `x.bam` in a second job (4494856) because
the score table is no longer written next to the input BAM.
One sample, one tissue; the mapping recipe is unchanged.

# PREREG: minimap2 two-pass alignment (--write-junc, then --pass1) and the recipe -N 50 -p 0.8 (row 1894, 2026-11-22)

User (2026-11-22): "what if we tried the --write-junc / --pass1". Written before any two-pass alignment of these
reads. Seen: a toy run on 100 planted reads showing the two options run in 2.30 (`--write-junc` writes one line per
read junction; `--pass1 FILE` reads it); row 1888. Workspace `r1894_minimap_two_pass`.
Also settled here: the user is right that our BAMs carry the sequence in secondary records (`-Y`; 2,000 of 2,000
secondaries checked in `GGO_mm.bam`), so `--secondary-seq` adds nothing.

## Arms (`minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.8 --secondary=yes`)
- ONE: one pass.  TWO: a first pass with `--write-junc` (same options, no SAM), then the same command with `--pass1`.
Samples: planted family reads of row 1874 (E, S: truth known); gorilla OR6737, all reads.

## Measures
- Planted: reads whose primary is on the true locus, placed alone on a wrong locus, tied (ms, exact); share on the
  true locus EM / EQUAL; expressed and silent copies reached (row 1885's measure).
- OR6737: primaries whose position or CIGAR differs between the arms; molecules with a near-tie secondary; distinct
  introns of the primaries and the share on RefSeq introns; named copies / junction precision; family pair precision /
  recall; transcripts. ONE is also read against row 1888's USUAL (`-N 50 -p 0.1`) for the recipe.
- TWO HELPS if, on the planted reads, wrong placements + tied reads fall by >= 10% (relative) in E and S with no
  expressed copy lost, or on OR6737 the named copies rise by >= 1 or their junction precision by >= 0.01 with pair
  precision and recall not lower by 0.01. TWO HURTS on the mirror image. Otherwise NO DIFFERENCE.
- Recipe: ONE is as good as USUAL when named copies are within 1 and pair precision / recall within 0.01.

## Predictions (mine)
NO DIFFERENCE: under 1% of the primaries change; the share of read introns on RefSeq introns rises by under 0.5
points. ONE as good as USUAL.

## Amendment to row 1894 (2026-11-22, before any two-pass result)
minimap2 2.30 aborts when `--pass1` meets `--eqx` (`jump.c:198: mm_jump_split: Assertion (opt->flag & 0x4000000) == 0`;
reproduced on 100 reads, gone without `--eqx`). The two-pass arm cannot be run with our usual options. Both arms are
therefore run again WITHOUT `--eqx` (`one_m`, `two_m`) and compared with each other; Rustle reads =/X operators in some
stages, so the first `one` (with `--eqx`, already aligned) stays as the recipe check against USUAL and shows what the
missing =/X costs. Seen before this amendment: `one` on the planted reads only (identical to row 1874's numbers).

## Result of row 1894 (2026-11-22, jobs 4494859, 4494891, 4494892, 4494893): two-pass NO DIFFERENCE; -N 50 -p 0.8 as good as USUAL
Two-pass against one pass, both without `--eqx`:
| | one pass | two pass |
|---|---|---|
| planted E: on true locus / wrong locus / tied | 4,651 / 26 / 3,285 | 4,653 / 29 / 3,280 |
| planted S: on true locus / wrong locus / tied | 2,364 / 6 / 1,814 | 2,366 / 6 / 1,812 |
| planted copies reached, E; S expressed, silent | 133; 70, 41 | 133; 70, 41 |
| OR6737 distinct introns of primaries, on RefSeq | 585,667, 36.64% | 585,669, 36.64% |
| intron uses on RefSeq | 97.70% | 97.70% |
| introns seen once | 290,063 | 289,554 |
| near-tie secondaries | 108,310 | 108,310 |
| named copies / junction precision (mode) | 34 / 0.66 | 34 / 0.66 |
| family pairs predicted / hit | 2,135 / 192 | 2,138 / 192 |
| transcripts | 70,167 | 70,160 |
71,173 of 4,404,421 primaries (1.6%) differ in position or CIGAR. Alignment 39 min against 22. Two-pass cannot be
combined with `--eqx` in minimap2 2.30 or 2.31.
Recipe: `-N 50 -p 0.8` one pass (with `--eqx`) against USUAL `-N 50 -p 0.1`: named copies 34 / 0.66 both, family pairs
2,135 / 192 both, transcripts 70,167 against 70,181, near-tie secondaries 108,310 against 110,012 (98.5%), secondaries
473,973 against 6,266,238, BAM 5.8 GB against 11.7. AS GOOD AS USUAL by the rule.
Dropping `--eqx` changed none of these measures (same transcripts, copies, pairs); stages that read =/X (copy
assignment, phasing) were not run here.
Predictions: no difference right; "under 1% of primaries change" wrong (1.6%); RefSeq intron share rise under 0.5
points right (0.00); recipe right.
