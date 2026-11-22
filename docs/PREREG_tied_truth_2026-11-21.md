# PREREG: which rule shares tied reads correctly, tested where the origin is known (rows 1873, 1874, 2026-11-21)

User (2026-11-21): "What can we try that gives us an objective answer" to the open decision of row 1871 (keep the EM
share of tied reads, or split them equally inside a tie group); then "yes lets do that".
Written before any fragment or simulated read exists. Seen: row 1871's tables on A119b (679 tie groups; inside groups
517,601 placed molecules on the exact band, 430,008 on the wide band; the EM gives 0.0 to identical copies without a
placed read).

## The two rules
For a read tied between a set C of loci (minimap2 `ms`, exact band, `tools/eq_classes.py` as the stage runs it):
- EM (current): share of locus L = theta_L / sum of theta over C, theta from `eq_classes.em`.
- EQUAL: share = 1 / |C|.
Both leave placed reads where they are. `outside` (no locus) is a member like any other.

## Row 1873: real reads with the answer hidden (workspace `r1873_tied_truth_clip`)
- Parents: A119b molecules whose primary is their best record, with no second record within 10 of the best `ms` (the
  wide band), lying on a locus of row 1871's `a1/x.multicopy.gtf` that is in a tie group. Truth = that locus. This is
  the aligner's placement with a margin, not a proof of origin.
- Fragments: each parent's sequence cut in consecutive 500-base pieces (a remainder under 500 is dropped), aligned with
  the sample's command (`minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`) to CHM13, `as_table`, then
  classes over the same GTF. The fragment set stands alone (no other reads): abundances come from the fragments.
- Test set: fragments that are tied (class of >= 2 members).
- Primary measure: mean share a rule gives the true locus, over the test set (0 when the true locus is not in the class).
- Secondary: misallocation = sum over loci of |predicted tied - true tied| / (2 x test set); fragments placed alone on a
  locus other than the parent's (how wrong `placed` can be); the named families separately (descriptive).
- A rule WINS row 1873 if its primary measure is higher by >= 0.05 pooled and it is higher in the majority of tie groups
  with >= 50 test fragments. Otherwise NO WINNER.
- Cannot answer: copies with no distinguishing stretch have no parent (the 14 identical TSPY copies).

## Row 1874: planted reads (workspace `r1874_tied_truth_sim`)
- Copies: the annotated protein-coding members of `MAINSTAYS` + `YAGS` on CHM13 (`famtruth.py`), one transcript per copy
  (most exons, then longest). Loci for counting = these transcripts (gene_id = copy): the counting rule is isolated from
  the assembly.
- Reads: the transcript's sequence; 70% full length, 30% cut at the 5' end at a uniform point leaving >= 500 bases;
  errors per base 0.1% substitution, 0.2% insertion, 0.2% deletion; seed 1874. Same aligner command and tables.
- Scenarios, per family, fixed here: E every copy 60 reads; D one copy (drawn with the seed) 600, the others 20;
  S a drawn half of the copies 60, the other half 0 (families of >= 2 copies; odd sizes round the silent half down).
- Measures per scenario: the primary measure and misallocation of row 1873 on the tied reads; silent copies given >= 2
  reads (false expressed), expressed copies given < 2 (false silent), by rule; whether the true count lies in
  [placed, max] for every copy.
- A rule wins a scenario by >= 0.05 on the primary measure.
- Simplifications, disclosed: the sample equals the reference (no private variants), error and truncation models are
  simple, three expression scenarios of my choosing. A scenario, not a mimic of real data (rows 1797, 1809).

## Decision
- EQUAL replaces EM in `reads_tieadj` only if it wins row 1873 and at least 2 of the 3 scenarios.
- EM stays if it wins row 1873 and at least 2 of 3.
- Anything else: NO RULE WINS; EM stays as the reported share, and for copies without a placed read the documentation
  says to quote the range. A hybrid suggested by the results would be post hoc and need its own test.
The default is the user's call either way.

## Predictions (mine)
1873: EM wins (real expression is uneven). 1874: E no winner or EQUAL; D EM; S EM on the primary measure, while EQUAL
calls most silent copies expressed and EM calls some expressed identical copies silent. True counts inside
[placed, max] everywhere. Overall: EM stays.

## Result (2026-11-21, jobs 4494121, 4494124): NO RULE WINS; EM stays, quote the range
Row 1873 (430k parents -> 1,829,391 aligned pieces; 297,936 tied = test set): EM WINS. Mean share on the true locus EM
0.507, EQUAL 0.396; misallocation 0.228 / 0.321; EM higher in 184 of 244 tie groups with >= 50 test pieces (medians
0.508 / 0.440). In 43,838 test pieces (15%) the true locus is not in the tied class at all. Of 1,531,455 pieces placed
alone, 1,216,086 are on the parent's locus, 143,166 (9%) on another locus, 172,203 on none: a 500-base piece can fall
on an overlapping locus of the GTF, so this is an upper bound on wrong placement, not a rate of paralogue errors.

Row 1874 (133 copies, 15 families; primary measure EM / EQUAL; per-copy count error EM / EQUAL / placed only):
| scenario | reads | tied | primary | winner | count error | silent copies called expressed | expressed called silent |
|---|---|---|---|---|---|---|---|
| E all 60 | 7,980 | 3,285 | 0.279 / 0.337 | EQUAL | 0.084 / 0.021 / 0.208 | - | EM 15, EQUAL 0 of 133 |
| D 600 + 20s | 11,360 | 4,733 | 0.313 / 0.320 | none | 0.204 / 0.217 / 0.210 | - | EM 15, EQUAL 0 of 133 |
| S half silent | 4,200 | 1,814 | 0.305 / 0.338 | none | 0.158 / 0.143 / 0.218 | EM 10, EQUAL 25 of 63 | EM 6, EQUAL 0 of 70 |
The planted count is outside [placed, max] for 8 / 7 / 4 copies (NPIP, GOLGA6L, TSPY3; mostly by 1-10 reads; reads
aligned outside the locus or placed on a neighbour): the range is not a guarantee. 51 of 133 copies have no placed read
in E (tied only).

Decision rule: EM wins 1873 but no scenario; EQUAL wins one scenario. NO RULE WINS. EM stays as `reads_tieadj`.
Predictions: 1873 right; E right; D wrong (EM does not win: when the dominant copy is identical to others it has no
placed read to pull the tied reads); S wrong on the primary measure, right that EQUAL calls more silent copies expressed
(25 against 10) and that EM calls expressed identical copies silent (6); "inside the range everywhere" wrong.
Reading (mine): the two tests differ because on real data the copies with more placed reads do hold more of the tied
ones, while in my scenarios expression was drawn independently of which copies are distinguishable.

# PREREG: other near-tie cases, gene + processed pseudogene and tandem pairs (row 1875, 2026-11-21)

User (2026-11-21): "can we try to model some other cases of nearly tied secondary alignments like a gene and a
pseudogene or tandem copies?" Written before any locus is drawn or read written. Same planted-read method as row 1874
(`r1874_tied_truth_sim/sim.py` generalised to a plan file; same error and truncation model, aligner command, scoring).
Workspace `r1875_tied_truth_cases`, seed 1875.

## Sets (CHM13 CAT/Liftoff annotation; one transcript per locus, most exons then longest)
- PSEUDO: protein-coding genes with one locus whose symbol S has processed pseudogenes named `S P<number>` (biotype
  processed_pseudogene or transcribed_processed_pseudogene). 150 parents drawn; all their pseudogenes, at most 10 each.
- TANDEM: Ensembl paralogue pairs with last common ancestor Homo sapiens to Catarrhini, both protein-coding, same
  chromosome and strand, not overlapping, gap <= 200 kb. 150 pairs drawn, no gene twice.
- DISPERSED (the control for TANDEM): same ancestor levels, different chromosomes or gap >= 5 Mb. 150 pairs, no gene twice.
No gene is in two sets. Fewer than 150 available = all taken, said in the result.

## Scenarios (one alignment each, all sets together)
- E: PSEUDO gene 60 reads, each pseudogene 20; TANDEM / DISPERSED both copies 60.
- S: PSEUDO gene 60, pseudogenes 0; TANDEM / DISPERSED one copy of the pair (drawn) 60, the other 0.

## Measures, per set and scenario
- Share of planted reads that are tied, and placed on a wrong locus, under three readings of a tie: `ms` equal (what the
  stage uses), `ms` within 10, AS equal.
- On tied reads: mean share on the true locus, EM and EQUAL (the measure of rows 1873-1874), a win = 0.05.
- Silent loci given >= 2 reads, expressed loci given < 2, by EM / EQUAL / placed only.
- TANDEM against DISPERSED only: planted reads whose primary record covers >= 50 exon bases of both copies of the pair
  (an alignment running across two copies).
Descriptive: no default hangs on it. If a rule wins both scenarios in a set where it lost in row 1874, said so.

## Predictions (mine)
- PSEUDO, S: under `ms` equal under 5% of the gene's reads tie with a pseudogene and EM calls fewer than 5% of silent
  pseudogenes expressed; under AS equal the tied share is at least twice the `ms` share (AS charges for introns).
  EQUAL calls more silent pseudogenes expressed than EM under every reading.
- TANDEM / DISPERSED: E no winner or EQUAL; S EM wins (a pair's expressed copy mostly has placed reads of its own).
  Tied shares similar in the two sets; alignments across two copies above 0 in TANDEM, 0 in DISPERSED.

# PREREG: the no-locus member of a tie, pooled or per class (row 1877, 2026-11-21; 1876 was taken by another session)

Found while reading row 1875's output (seen: its tables; a rescoring of row 1875 with the change, `ms` equal: TANDEM
EM 0.170 -> 0.306 / 0.326, DISPERSED 0.211 -> 0.492, 0.240 -> 0.560). `eq_classes.py` gives every tied record that lies on
no locus the single member `outside`, so the EM treats all such places in the genome as one locus: it gathers weight
from every class and takes the tied reads away from real loci (row 1875, S: `outside` theta 1,928). Candidate: one
`outside` member per class (the places off the loci are different places).
Not seen: the candidate on row 1873's real pieces or on row 1874's scenarios. Test = exactly those, same measure.
- The candidate is worth proposing to the user as a change to `eq_classes.py` if, against the current EM, it is higher
  by >= 0.02 on row 1873's pieces and not lower by more than 0.01 in any of row 1874's three scenarios.
- No code of `tools/` changes in this row; a variant of the scoring script only (`UNPOOL=1`).
Prediction (mine): small gain on row 1873 (under 0.02: there the GTF is a whole assembly, few ties reach a place with
no locus), gains in row 1874. So: not enough to propose on the registered bar.

## Result of row 1875 (2026-11-21, job 4494558, `r1875_tied_truth_cases`)
Drawn: 150 tandem pairs (of 2,332), 128 dispersed pairs (all that are disjoint, of 982), 150 genes with 417 processed
pseudogenes. Reads E 50,700, S 25,680. Tie reading `ms` equal unless said.

| set | scenario | tied share of reads | placed on a wrong locus | EM / EQUAL | silent called expressed EM / EQUAL | expressed called silent EM / EQUAL |
|---|---|---|---|---|---|---|
| gene + pseudogene | E | 0.005 (ms within 10: 0.023; AS: 0.020) | 0 of 17,320 | 0.228 / 0.500 (88 reads) | - | 3 / 1 of 567 |
| gene + pseudogene | S | 0.000 (ms within 10: 0.000; AS: 0.031) | 0 of 9,000 | - | 0 / 0 of 417 (AS: 1 / 5) | 0 / 0 |
| tandem | E | 0.305 | 9 of 18,000 | 0.170 / 0.280 | - | 48 / 0 of 300 |
| tandem | S | 0.288 | 6 of 9,000 | 0.171 / 0.281 | 16 / 42 of 150 | 22 / 0 of 150 |
| dispersed | E | 0.103 | 62 of 15,360 | 0.211 / 0.445 | - | 11 / 0 of 256 |
| dispersed | S | 0.101 | 54 of 7,680 | 0.240 / 0.455 | 6 / 10 of 128 | 6 / 0 of 128 |
No primary record ran across both copies of a pair, in TANDEM or DISPERSED.
- Gene + pseudogene is not a tie problem under `ms`: no read of the gene tied with a silent pseudogene, none of 417
  silent pseudogenes got a read. Under AS 3.1% of the gene's reads tie and 1 (EM) / 5 (EQUAL) pseudogenes are called.
- Tandem pairs tie three times as often as dispersed ones (0.29-0.31 against 0.10).
- EQUAL beats EM in every pair row, by 0.11-0.23. Cause, found by reading the classes: the pairs come from larger
  families whose other members are not loci of the test GTF, so their ties hold the member `outside`; the tool pools
  every `outside` in the genome into one EM member, which takes the tied reads (S: theta 1,928). POST HOC with one
  `outside` member per class: tandem 0.306 / 0.326, dispersed 0.492 / 0.560 (EQUAL unchanged), expressed called silent
  5 and 1 (from 48 and 22). Tested as row 1877.
Predictions: pseudogene under 5% and AS at least twice ms right; EQUAL calls more silent pseudogenes right only under
AS (0 = 0 under ms). Tandem / dispersed: "S EM wins" wrong, "tied shares similar" wrong (3 x), "alignments across two
copies above 0 in tandem" wrong (0).
Limit added: a pair drawn from a larger family leaves the family's other members off the GTF; the pair rows measure the
tool on an incomplete set of loci, which is what an assembly with missing copies is too.

## Result of row 1877 (2026-11-21, job 4494582, `r1877_outside_unpooled`; first job crashed on my variant's group lookup, fixed, no number seen before): BAR MET, proposed to the user
Mean share on the true locus, EM pooled (current) -> EM with one `outside` per class; EQUAL for reference.
| test | EM current | EM per-class outside | EQUAL |
|---|---|---|---|
| row 1873 real pieces | 0.507 | 0.542 (ahead of EQUAL in 218 of 244 groups, was 184) | 0.396 |
| row 1874 E | 0.279 | 0.340 | 0.337 |
| row 1874 D | 0.313 | 0.365 | 0.320 |
| row 1874 S | 0.305 | 0.376 | 0.338 |
Gain 0.035 on the real pieces (bar 0.02), no scenario lower (bar -0.01). Expressed copies called silent 15 -> 4 (E),
15 -> 5 (D), 6 -> 1 (S); silent copies called expressed 10 -> 16 of 63 (EQUAL 25). Misallocation on the real pieces
0.228 -> 0.191. My prediction (under 0.02 on the real pieces) was wrong.
Against EQUAL the changed EM still wins only row 1873 by the 0.05 margin (D +0.046, S +0.038, E +0.003): it no longer
loses anywhere. `tools/eq_classes.py` is unchanged; it also feeds stage `eqclasses`, so the change is the user's call.
