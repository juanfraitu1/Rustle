# PREREG EXON-SPACE PLACEMENT (register 1642, 2026-11-12)

Why: user question (2026-11-12): is there an alignment form that is faster and keeps the whole read in one alignment, so
that multi-mapping regions are placed better than by minimap2's spliced genome alignment? Candidate: align the read
UNSPLICED to the family copies' spliced sequences (`kb.fam.copies.fa`, 812 copies, 2.5 Mb): no intron cost, every copy
scored end to end on one scale, a 2.5 Mb index instead of a genome.

How it differs from closed rows: 1410/1411 scored the GENOME alignments' tags or a likelihood over them (agreed with ms
on 220/222); O2 (`copy_assign --copies-fa`) uses the same copy sequences but reads only PSV columns through a
certificate. This arm is a fresh whole-read alignment score in exon space, nothing else.

Expected outcome, stated before the run: NULL (the arm equals the primary). The test set is ms near ties, ms carries no
intron cost already, and 1410/1411 found the wrong primaries are other-allele reads whose sequence favours the wrong copy.

Sample and set: KB3781 (Jim) on mGorGor1, the r1514 rows: `o2mstie/pof_panel` + `pof_heldout` contested rows, truth =
`strict_truth.py` at slack 0 on `truth_all10/hits.tsv` (dumps regenerated with `STRICT_DUMP_OUTSIDE=1`, decidable part
byte-identical to the r1514 dumps). Counted per MOLECULE (distinct read name), not per row.
- decidable: 1,207 panel + 476 held-out rows (primary right 1,123/84 and 351/125).
- negative set: rows with a near-perfect haplotype hit on no catalog copy (panel 1,062, held-out 792 rows).
Read sequences: `truth_all10/query.fa`.

Arm (one, `bench/o2_contested/exon_place.py`): `minimap2 -cx map-hifi --eqx -N 100 -p 0.1` reads vs all 812 copy
sequences. Per read, a copy's score = the best AS of its records. RULE (the only registered one): assign the read to the
top copy iff (a) its AS exceeds by >= 1 the best AS of any copy at another locus (copies whose catalog spans overlap are
one locus) and (b) the top record covers >= 0.90 of the read; else abstain. Right = the copy is in the truth set.

Bars (all pre-declared):
1. Reach: the arm assigns >= 30% of decidable molecules.
2. Paired vs the primary, decidable molecules the arm assigns: (arm right, primary wrong) >= 2 x (arm wrong, primary
   right), and the first count >= 10.
3. Bar 2's ratio holds on held-out alone and with the largest family removed (first count >= 5 each).
4. Forced placement: the arm assigns <= 10% of negative-set molecules.
5. Cost: alignment + scoring wall <= 10 min on 8 threads.
Readout, no bar: vs O2 on molecules both assign; per family; rows where the truth copy got no record; margins >= 5 and
>= 10; cover gate off; discordance by the primary's splice shape.

Verdicts: bars 1-4 pass = PASS, which licenses ONLY a row on a sample that is not the reference individual (reads,
reference and truth haplotypes are all Jim here), never adoption. Bar 2 fails with discordant counts both small
(sum < 15) = NULL as expected: record, close exon-space re-alignment as a placement lever. Bar 2 fails otherwise or bar 4
fails = FAIL. No second threshold is promoted after the fact.

Outputs: `/anvme/workspace/v138da13-rustle/r1642_exonplace/`.

## RESULT r1642 (2026-11-12, login node, 1 s; `r1642_exonplace/result.txt`)
Molecules: 1,570 decidable (panel 1,130, held-out 440), 1,419 negative.
- Bar 1 PASS: assigns 1,097 of 1,570 (0.699); abstains 196 tie, 157 cover, 120 no record.
- Bar 2 FAIL: arm right / primary wrong 15, arm wrong / primary right 44 (needed 15 >= 88). On the assigned molecules
  the arm is 1,029 / 68, the primary 1,058 / 39.
- Bar 3 FAIL: held-out 11 vs 21; minus MCL115 15 vs 44.
- Bar 4 FAIL: assigns 647 of 1,419 negative molecules (0.456, bar 0.10): forced placement, the r1510 class.
- Bar 5 PASS: 1 s on 8 threads.
**FAIL as registered** (not the expected null: 59 discordant molecules, three to one against the arm).
Readout: margin >= 10 gives 12 vs 32, cover gate off 32 vs 78: no threshold turns it. Against O2 on the 429 molecules
both assign: arm 399 right, O2 411. By the primary's shape: unspliced primary 3 vs 19, spliced primary 12 vs 25. MCL5
alone is 1 vs 16. The truth copy got no record on 172 of 1,570 molecules (one model per copy: other isoforms and
unspliced reads do not align end to end).
Checks (`r1642_exonplace/split.txt`): scoring the arm by locus overlap, as the primary is scored, changes nothing (15 vs
44). The 44 split: truth copy has no record 15, a record under 0.90 cover 10, a full record with a lower AS 19. O2
assigns 230 of the same 1,419 negative molecules (0.162).
Reading: 25 of the 44 are the catalog's one model per copy not holding the read, 19 are reads whose exon sequence
prefers the wrong copy. Exon-space re-alignment against a one-model-per-copy catalog is closed as a placement lever; a
multi-isoform reference was not tested. This set excludes the spliced-vs-unspliced class (1403/1404).
