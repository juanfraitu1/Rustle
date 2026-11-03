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

## RESULT (2026-11-03, job 4421542, KB3781: 13,516,443 molecules, 34,762,361 records; workspace o2contested/r1511/)
r1511 useful fraction by ms-ratio bin (useful = known paralog by the primates pair truth, or a catalog copy of the
primary's family): <0.5 4.8% (14.45M records) | 0.5-0.8 6.9% (4.57M) | 0.8-0.9 16.2% (1.27M) | 0.9-0.95 14.9% (469k) |
0.95-0.99 33.4% (325k) | >=0.99 71.4% (159k). By the pre-declared reading the lowest bin at >= 5% is 0.5-0.8, so the
recommended lower bound is ms ratio 0.5: minimap2's default `-p 0.8` would have dropped 4.57M records of which 316k are
useful (23% of all useful secondaries); `-p 0.1` was right to keep them, and below 0.5 the 14.45M records hold 689k useful
(4.8%, under the bar). In every bin most secondaries land on an annotated gene that is neither a known paralog nor a catalog
copy ("other_gene": 68-96%); the AS-ratio median exceeds 1 above ratio 0.9 (spliced primary vs unspliced rival: the AS-gap
class of r1508). Standing rule for the records: ms ratio >= 0.5 is the informative set (verdict: REPORT, bound = 0.5).

r1512 multimapper gene graph (K = 3 molecules; nodes = any annotated gene with an exon under the record), scored at the
primates tier against O1's families (`kb.fam`: pair precision 0.341, recall 0.136 on 103 truth pairs):
R >= 0.99: precision 0.032 recall 0.019 | 0.95: 0.022 / 0.039 | 0.9: 0.019 / 0.097 | 0.8: 0.016 / 0.155 | 0.5: 0.001 /
0.447 | all: 0.000 / 0.981. Silent protein-coding genes joined with a truth pair in their component: 0 of 27 (R 0.99), 0 of
151 (0.95), 1 of 295 (0.9), 7 of 537 (0.8), 14 of 1,041 (0.5), 29 of 2,235 (all). Both bars FAIL at every R: the graph's
edges are dominated by secondaries on unrelated genes (repeats and long genes under the record), so connected components
are not families and silent members are not found this way. ❌ FAIL. The multimapper signal that holds is per locus (r1513,
docs/PREREG_secondary_pool_2026-11-03.md), not a gene-level graph.
