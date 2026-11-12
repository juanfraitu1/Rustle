# PREREG PSEUDOALIGNMENT PLACEMENT, lr-kallisto (register 1647, 2026-11-12)

Why: user question (2026-11-12): can a method that is not minimap2 place multi-mapping long reads better? Rows 1642-1644
closed exon-space re-alignment WITH minimap2. kallisto's long-read mode is the one installed method that shares no code
with it: k-mer pseudoalignment (k 63) against transcript sequences, a read gets the set of sequences it is compatible
with.

Samples and sets (existing files, nothing rerun): KB3781 = r1642's reads, dumps and truth (1,570 decidable molecules,
1,419 negative); HG002 = r1644's (2,412 decidable, 8,098 negative). Per molecule, all primary shapes; shape classes are
readout.

Scorer: `exon_place.py` with `EXON_PLACE_RIGHT=locus` (a pick whose span overlaps a truth copy is right, as the primary
is scored). The literal rule mis-scored 171 of 214 on HG002 (r1644); the locus rule is the registered one here.

Reference (new, `multi_iso_ref.py ... other`): r1643's sequences per catalog copy (cat, tx, pre) PLUS every other
transcript of the sample's assembly (not on a catalog copy, or over two loci) under the name OTHER. Reason: a
catalog-only reference failed the forced-placement bar three times; a quantifier is always run on a whole
transcriptome.

Arm: `kallisto index -k 63`; `kallisto bus -x bulk --long --threshold 0.8 -n` (kallisto 0.52.0; the documented long-read
settings, no sweep). RULE: assign the read to a copy iff every sequence in its compatibility class belongs to copies of
one locus and none is OTHER; else abstain (`bench/o2_contested/kallisto_place.py`).
Control (reported, no bar): minimap2 `-cx map-hifi --eqx -N 100 -p 0.1` on the SAME reference, r1642's rule, top = OTHER
abstains. It separates the reference's effect from the tool's.

Validity: kallisto must pseudoalign >= 80% of the reads in each sample, else the run is INVALID (not FAIL).

Bars, per sample, both must pass:
1. Reach: assigns >= 30% of decidable molecules.
2. Tool right / primary wrong >= 2 x (tool wrong / primary right), first count >= 10.
3. Negative molecules assigned <= 10%.
4. Wall <= 30 min per sample on 8 threads.

Expected, stated before the run: bar 2 fails. r1642 found 19 of 44 losses are reads whose sequence prefers the wrong
copy, r1645 (FLAIR-like realignment) also lost to the primary; a sequence-only method inherits both. Exact k-mer
compatibility may also abstain on most near ties (bar 1).

Verdicts: all bars on both samples = PASS (licenses a design row, not adoption). Otherwise FAIL; with the control it is
recorded whether the reference with OTHER sequences fixes forced placement.

Not in this row: oarfish (its read mode maps with rammap and then runs an abundance EM = the r1412 class; needs every
read of the panel regions, not the contested reads alone) and graph alignment (step 0, separate row).

Outputs: `/anvme/workspace/v138da13-rustle/r1647_kallisto/`.

## RESULT r1647 (2026-11-12, login node, 8 threads; `r1647_kallisto/{kb,hg}.{result,control}.txt`)
Reference: KB3781 74,328 sequences (68,810 OTHER), HG002 162,255 (135,831 OTHER).
Validity FAILS on both samples: kallisto pseudoaligned 1,234 of 2,989 reads (0.413) on KB3781 and 7,341 of 10,510
(0.698) on HG002, bar 0.80. **INVALID as registered.** The gate was meant to catch a broken run; what it caught is the
tool's own behaviour (a read without enough exact 63-mers gets no class), so the bars are reported below as read.

| | KB3781 | HG002 |
|---|---|---|
| Bar 1 reach (>= 0.30) | 937 of 1,570 (0.597) pass | 684 of 2,412 (0.284) FAIL |
| Tool right / wrong on its molecules | 919 / 18 | 678 / 6 |
| Primary right / wrong on the same | 898 / 39 | 612 / 72 |
| Bar 2: tool right, primary wrong vs reverse | 23 vs 2 pass | 68 vs 2 pass |
| Bar 3 negatives assigned (<= 0.10) | 125 of 1,419 (0.088) pass | 4,681 of 8,098 (0.578) FAIL |
| Bar 4 wall | 97 s | 261 s |
| vs O2 on shared molecules | 247 vs 247 of 259 | 346 vs 345 of 346 |
| minimap2 control, same reference: right/primary-wrong vs reverse | 32 vs 22 | 25 vs 42 |
| control negatives assigned | 0.574 | 0.291 |

Reading:
- Bar 2 holds on both samples and on every cut made (KB3781 held-out 14 vs 0, panel 9 vs 2; HG002 minus MCL338 57 vs 1),
  against the stated expectation. Where kallisto commits, it is as accurate as O2 and better than the primary.
- It is not the reference: minimap2 on the same reference stays level with the primary and still force-places.
- It is a selection: 41-70% of reads get a class, and the rule keeps one-locus classes only.
- Caveat that may explain all of it: reads, reference and truth are one individual, so a read matches its true copy
  exactly, and exact 63-mer compatibility is nearly the truth's own criterion (unique best edit distance on the same
  individual's haplotypes). On a sample that is not the reference individual, alleles break exact matches.
- HG002's negative rate (0.578) is unexplained; KB3781 passes that bar.
Verdict: INVALID as registered; no PASS under any reading (HG002 fails bars 1 and 3). Open, not closed: a row with the
validity gate replaced by a reach bar, on a non-self sample, plus a look at the HG002 negatives.

## POST HOC CHECKS (2026-11-12, flagged; `bench/o2_contested/place_checks.py`, `r1647_kallisto/{kb,hg}.checks.txt`)
1. Scoring asymmetry: the primary is right only if its leftmost position lies in the truth span. Rechecked with the
   primary RECORD's overlap: 8 of KB3781's 23 and 9 of HG002's 68 were right after all. Corrected: **15 vs 2 and 59 vs 2**.
2. Self-reference, HG002 (reference = paternal haplotype), by the haplotype of the read's best hit: reach paternal 224 of
   522 (0.43), both 444 of 1,740 (0.26), maternal 16 of 150 (0.11); the 16 are all right. Reach falls four-fold on reads
   from the haplotype that is not the reference; accuracy shows no loss on a tiny count. The caveat is real for reach.
3. HG002 negatives: for 4,474 of the 4,681 assigned, the pick's span holds the primary. Aligner and tool agree and the
   haplotype truth says outside the catalog: a truth or catalog question, not forced placement by the tool (KB3781: 48 of 125).
4. Read-number parsing: reads in reverse order give identical picks.
5. Beyond O2: tool assigns where O2 abstains 672 right / 6 wrong (KB3781) and 332 / 6 (HG002); O2 assigns where the
   tool abstains 270 / 36 and 135 / 44.
6. 10-AS rule on the same decidable molecules: decides 179 of 1,570 (165 / 14) and 179 of 2,412 (149 / 30).
The verdict stays INVALID as registered. What the checks leave standing: exact long k-mer compatibility with a one-locus
rule is a precise, cheap certificate (wrong 2% and 1%), complementary to O2, with reach tied to how close the reads are
to the reference.
