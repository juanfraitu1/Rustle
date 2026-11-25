# PREREG: IsoQuant with `--use_secondary`, and bambu's supplementary alignments (row 1947, 2026-11-25)

Closes the two open items of `docs/MULTIMAPPERS_ACROSS_ASSEMBLERS_2026-11-25.md`.

## A. IsoQuant `--use_secondary` as a registered arm
IsoQuant has a switch to read secondary alignments. In its source a read with several placements is flagged
`multimapper` and then builds nothing (intron graph, paths). On the one set looked at for the study (seed 1951, S; NOT part
of this test) the switch took IsoQuant from 123 to 30 transcripts.
Data: the planted sets of row 1946 (seeds 1964-1966; E and S at 60 reads, E at 3 reads), on which IsoQuant has not been run.
Arms: `isoquant -d pacbio_ccs` (as in the benchmark) and the same with `--use_secondary`.
- **Q1:** with `--use_secondary` IsoQuant finds FEWER expressed units than its default in each of the six sets at 60
  reads, and its isoform recall is lower in each.
- Reported: units, false units, isoform recall and precision of both arms on all nine sets.
Prediction (mine): Q1 holds, by a wide margin.

## B. bambu and supplementary alignments (description, no claim)
bambu loads supplementary alignments (the second part of a split read); the other tools drop them. On gorilla OR6737:
how many supplementary alignments there are, how many touch a family copy, and whether bambu's models at family loci
change when they are removed. Reads overlapping the family regions of the run of row 1920 (`g.fam.copies.regions`), bambu
de novo at NDR 1 with and without the supplementary alignments; transcripts compared by exon structure.

## Result (jobs 4504644 bambu, 4504681 IsoQuant; the first IsoQuant attempt wrote nothing because of a bug in my job script, fixed and rerun)

### A. IsoQuant
| | default | `--use_secondary` |
|---|---|---|
| expressed units found, all expressed, 60 reads (of 108) | 102, 99, 100 | 35, 33, 37 |
| expressed units found, half silent (of 60, 61, 61) | 56, 56, 56 | 19, 22, 24 |
| expressed units found at 3 reads (of 108) | 66, 59, 69 | 12, 11, 11 |
| false units, half silent (of 48, 47, 47) | 13, 12, 11 | 2, 3, 3 |
| isoform recall, all expressed | 0.611, 0.598, 0.594 | 0.134, 0.142, 0.138 |
| isoform recall, half silent | 0.615, 0.595, 0.571 | 0.148, 0.190, 0.143 |
| isoform precision, all expressed | 0.802, 0.822, 0.821 | 0.667, 0.739, 0.647 |
| transcripts (E / S / 3 reads, seed 1964) | 193 / 124 / 79 | 48 / 28 / 12 |

**Q1 HOLDS** in all six sets, by a wide margin: with its secondary switch on, IsoQuant finds about a third of the units
and a quarter of the isoforms it finds by default. The mechanism is in its source: a read with several placements is
flagged `multimapper` and builds nothing, so reading the secondary alignments REMOVES the family reads from model
construction. Its false units fall too (13 / 12 / 11 -> 2 / 3 / 3), because it models far fewer loci at all.

### B. bambu and supplementary alignments (gorilla OR6737)
34,699 supplementary alignments beside 4,409,946 primary ones (0.8%); at the 1,345 family loci (5 kb each side) 1,473 beside
117,693 (1.3%). bambu at those loci: 1,873 transcripts with the supplementary alignments, 1,862 without, 1,843 with the same
exon structure; 30 only with them, 19 only without. So bambu's habit of loading supplementary alignments changes about
2.6% of its models at family loci on this sample. It is a difference from the other tools, a small one, and not a use of
multi-mapping evidence (a supplementary alignment is a second piece of one split read, not a second placement).
