# PREREG: named families one by one, loci and isoforms per tool (row 1859, 2026-11-21)

User (2026-11-21, after row 1857): "we need to look at it on a family by family way, for example the mainstay gene families
or YAGs. Which tool finds more loci and isoforms?" Written before running; descriptive, small counts.

## Measurement
`bench/denovo_eval/toolcompare.py` (2026-10-31) as it is. Families: human `MAINSTAYS` (NPIP, TBC1D3, GOLGA8, GOLGA6L, CT47A,
XAGE, SPATA31) + `YAGS` (the nine Y-ampliconic families) on A119b; `GGO_MAINSTAYS` on OR6737. GTFs = row 1855's inputs
(`r1847_family_benchmark/S/ARM/x.gtf`: rustle, stringtie, flair, isoquant, bambu), and row 1858's `primary` / `good`
Rustle assemblies when they exist. Per family and tool: copies with a spliced transcript (loci), spliced transcripts
(isoforms), share of their introns on annotated exon boundaries (junction precision).

## How it is read
- Loci: per family the tool(s) with most copies modeled; summary = families where Rustle is first or tied, and total copies.
- Isoforms: counts are reported WITH junction precision; a higher count at lower junction precision is not called better.
- "Rustle finds more loci" is said only if its total of modeled copies exceeds every other tool's on A119b.

## Prediction (mine)
Copies: Rustle and StringTie within a few of each other, FLAIR highest. Isoforms: FLAIR most, at the lowest junction precision.

## Limits
Annotated copies only (CAT on CHM13); a copy is "found" by any spliced transcript mapped to it, complete or not; one
testis sample per species.

## Result (2026-11-21): Rustle does not find more loci; StringTie does
Job 4491266, tables `r1859_family_by_family/{a1,or}.tsv`. Copies modeled / junction precision / spliced transcripts.

A119b, 110 copies expressed in any tool: stringtie 107 / 0.47 / 958; flair 98 / 0.35 / 1,291; rustle 93 / 0.56 / 699;
isoquant 85 / 0.65 / 258; bambu 61 / 0.78 / 94. Rustle is first or tied in 10 of 15 families. The 14 copies it lacks against
StringTie: TSPY 10 against 19, BPY2 1 / 3, PRY 0 / 1, TBC1D3 8 / 9, XAGE 3 / 4. In the seven mainstays Rustle has 67 of 69
expressed copies, StringTie 69. Rustle's junction precision is above StringTie's in 9 of 15 families.

OR6737, 30 copies: stringtie 28 / 0.59 / 155; rustle 25 / 0.61 / 97; flair 23 / 0.55 / 116; isoquant 23 / 0.56 / 83; bambu
18 / 0.68 / 36. Row 1858's assemblies: `primary` 25 / 0.61 / 96, `good` 27 / 0.58 / 132 (+1 LRPAP1, +1 JMJD7-PLA2G4B).

Predictions: FLAIR most isoforms at the lowest junction precision held; FLAIR most copies and Rustle within a few of
StringTie did not (StringTie most copies, 14 ahead of Rustle on A119b). A119b `primary` / `good` arms were not ready.
