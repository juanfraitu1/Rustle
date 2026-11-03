# PREREG SECONDARY USEFULNESS + MULTIMAPPER FAMILY GRAPH (registers 1511-1512, 2026-11-03)

Why (user, 2026-11-03): the BAMs were aligned with `-p 0.1 -N 50` on purpose, to learn up to what secondary/primary score
ratio secondaries carry information; minimap2's default `-p 0.8` is a candidate lower bound, but with 50 secondaries per
read some below it may still be useful. Second question (the advisor's): can the secondaries alone give a family's full
membership, silent copies included (register 1421 saw PGM5's silent pseudogenes through secondaries), i.e. is the
multimapper graph a family-discovery method, not only a near-tie detector?

Data: `as_table --records` dump of every primary and secondary record of KB3781 (name, secondary, chrom, pos, end, AS,
ms, spliced), `/anvme/workspace/v138da13-rustle/o2contested/KB3781.records.tsv`. Ratio = secondary ms / primary ms of the
same molecule (ms is minimap2's ranking score, the quantity `-p` thresholds approximate; AS ratio reported beside it).

## r1511 Secondary usefulness curve (report, pre-declared reading)
Per secondary record, by ratio bin (< 0.5, 0.5-0.8, 0.8-0.9, 0.9-0.95, 0.95-0.99, >= 0.99): where it lands:
(a) an annotated gene (DENOVO_GGO_GFF exon overlap) that is a paralog of the primary's gene by the pair truth O1 uses
(pairscore universe, primates tier); (b) an annotated gene that is NOT a known paralog; (c) no annotated gene; and
whether its locus is a catalog copy. "Useful" = (a) or lands on a catalog copy of the primary's family. Reading fixed
now: the bin below which useful secondaries are < 5% of the bin's records is the lower bound to recommend; if that bin is
at or above 0.8, minimap2's default would have lost nothing; if useful secondaries persist below 0.8 at >= 5%, `-p 0.1`
was right and the recommended bound is reported.

## r1512 Multimapper family graph (pre-registered bars)
Nodes = annotated genes (exon-overlapping loci); edge A-B when >= K molecules have their primary in A and a secondary in B
with ratio >= R (K = 3; R swept over the r1511 bins). Families = connected components. Scored with
`bench/denovo_eval/pairscore.py --species gorilla` on gene pairs exactly as O1's families are (pair precision / recall),
and membership of SILENT genes (annotated, zero primaries): how many silent genes join a component (precision of those
joins against the pair truth).
Bars: (1) at some R the multimapper families reach pair precision >= O1's current families (`kb.fam` pairscore) with recall
>= 0.5 x O1's: then the graph is a usable discovery method and is reported as such; (2) silent-gene joins precision
>= 0.8 at that R: the secondaries find silent members. Fail either: report the curve, no claim. IsoCon-style read
partitioning applies only where the reference LACKS the copy (reads' secondary mismatch patterns cluster into unseen
loci); that is `candidate_copies`' job and is out of this row.
