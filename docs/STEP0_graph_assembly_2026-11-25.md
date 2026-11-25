# STEP 0: ASSEMBLY IN THE FAMILY GRAPH'S COORDINATES (register 1944, 2026-11-25)

Why: user, 2026-11-25, after row 1943 (the graph does not replace the secondary alignments): "Test assembly on a graph".
Pooling reads or structure across copies was closed (FAMASM, rows 1628, 1822 B1, multi-copy ledger step 11: borrowed chains
15 of 90 expressed); this row reopens one form on the user's word. It cannot rebuild the assembler on a graph; it tests the
decision a graph assembly would make differently, on chains the assembler already built.

## The graph form
In the family's exon map (`exonmap`, row 1595: exon columns shared by the members = the nodes of the family graph) every
transcript is a chain of columns. A graph assembly counts support per column chain over the family, not per copy. Its one
new decision: a chain with ONE read at a copy is kept when the same column chain is a two-read transcript at another
member. Today the same job is done by the multi-copy stage: a second (one-read) assembly, candidates at EMPTY places, a
second families run, keep the candidates that join a family.

## Measurement (OR6737 cached run `r1920_mono_from_primaries/real_or`, development; multi-copy stage without the bypass)
- Families and members: `g.multicopy.d/u.fam` (assembly + candidates), families of >= 2 members; `exonmap` on it.
- Universe: multi-exon transcripts of the one-read assembly `g.multicopy.d/f1.gtf` whose intron chain is not in `g.gtf`,
  inside a member of such a family.
- Rule G: every exon maps to a column of its member (most shared bases), and the column chain equals that of a `g.gtf`
  transcript at ANOTHER member of the family with every internal exon within 3 bp in length. Loose arm (reported): columns only.
- Today's set D: the multi-exon transcripts tagged `multicopy_rescue` in `g.multicopy.gtf`.
- Annotated = intron chain equal to a gorilla RefSeq chain (`GGO_genomic.gff`).

## Stop rule (fixed before the run)
PROCEED (rule inside the stage, one unseen sample) only if both hold:
- (1) lost: of D's annotated chains, G keeps >= 0.80;
- (2) additions: the annotated share of G's chains >= the annotated share of D's chains.
Else STOP. Reported, not judged: G at empty and at occupied places, the loose arm, chains per family.

Prediction (mine): STOP on (1). Rescued copies rarely carry a chain a sibling has in the same columns (8% of isoforms are
shared between members, row 1822); G keeps tens to a few hundred chains at a higher annotated share than D (as row 1628:
0.14 against 0.11), and under 0.3 of D's annotated chains.

Outputs: `/anvme/workspace/v138da13-rustle/r1944_graph_assembly_step0/`. Script `bench/o2_families/graph_assembly_step0.py`.

## RESULT (2026-11-25, job 4504559; `r1944_graph_assembly_step0/step0.summary.txt`, `step0.tsv`)
Exon map of the families with the candidates: PASS by its own check (collision 0.0037, closure 0.979).
Universe: 5,824 new one-read chains inside a family member (59,036 one-read chains are outside any member).
| bar | value | verdict |
|---|---|---|
| (1) of D's annotated chains, G keeps | 3 of 54 (0.056); loose arm 4 of 54 | FAILS (bar 0.80) |
| (2) annotated share of G against D | 7 of 94 (0.074) against 54 of 1,262 (0.043); loose arm 12 of 165 (0.073) | holds |

**STOP, on (1).** Prediction held (under 0.3 of D's annotated chains; G keeps 94 chains at a higher annotated share).
- G and D barely meet: 22 chains are in both. 72 of G's 94 chains are at occupied places or outside D (4 annotated).
- 3,719 of the 5,824 chains have an exon in no column of the family's exon map: the one-read chains at a copy mostly use
  exons the family does not share, so there is no sibling chain to compare them with.
- Reading: the copies of a family do not carry the same chains in shared coordinates often enough for support counted over
  the family to stand in for the second assembly. Same finding as rows 1628 and 1822 (8% of isoforms shared), now for the
  rescue step itself. One sample, development; the multi-copy stage of the cached run is the one without the bypass.

## SECOND SAMPLE, A119b on CHM13 (2026-11-25, job 4504560, same script and bars; truth = `annot_related/chm13.related.gtf` chains)
(1) 2 of 28 (0.071), FAILS; (2) 17 of 524 (0.032) against 28 of 1,166 (0.024), holds. 14,835 of 18,662 one-read chains inside a
member have an exon in no shared column. Same verdict as OR6737: STOP on (1).
