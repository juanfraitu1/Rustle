# PREREG: do secondary alignments reveal silent family members, a silent gene from its expressed pseudogene and back (row 1909, 2026-11-22)

User (2026-11-22), the advisor's view: more secondary alignments are useful than we acknowledge; "if there were a
supposedly silent gene but we find a very similar pseudogene, we should be able to find the gene even if it is not
expressed and vice versa".
Written before any count. Seen: row 1875 (150 genes + 417 named processed pseudogenes planted; with the gene expressed
and the pseudogenes silent no read TIES with a pseudogene on ms, 3.1% on AS; nothing was counted below the near-tie
band); row 1888 (only 1.8% of a real sample's secondaries are near ties; the slim mapping mode `-p 0.8` drops most of
the others).

## Question
Not "is the member expressed" but "does the read set show that the member exists": a silent locus REACHED by secondary
records of reads that come from its expressed relative. How many are reached, at which score relative to the read's
best, and what the aligner settings keep.

## Runs (workspace `r1909_silent_members`; loci of row 1875: 150 genes, 417 processed pseudogenes, 150 tandem and 128
dispersed paralogue pairs; minimap2 as usual, `-N 50 -p 0.1`, and the slim mode `-p 0.8`)
- G: genes 60 reads, pseudogenes silent; one copy of every pair 60, the other silent (row 1875's scenario S, replanted
  with seed 1909).
- P: each pseudogene 20 reads, genes silent (the reverse).
A silent locus is REACHED when >= 2 reads have a record (any class) with >= 50 exon bases on it. Reported by the best
AS of such records relative to the read's best: >= 0.98 (the near-tie band the mode uses), 0.9-0.98, 0.8-0.9, 0.5-0.8,
< 0.5, cumulatively; and under `-p 0.8`.
Also: places reached by >= 2 reads' secondaries that are on no planted locus, merged within 1 kb (candidates for members
the annotation list lacks, or noise; no truth for them, a count only).

## Claims
- M1 (the advisor's): with every secondary of the usual mapping, at least half of the silent processed pseudogenes are
  reached from their gene's reads (G), and at least half of the silent genes from their pseudogenes' reads (P).
- M2 (they are not near ties): under 10% of those silent loci are reached by a record within 2% of the read's best.
- M3 (what the slim mode costs): `-p 0.8` reaches at most half as many silent pseudogenes / genes as `-p 0.1`.
Silent mates of tandem and dispersed pairs are reported the same way, without a claim.

## Predictions (mine)
M1: holds for G (a spliced gene read aligns to its intronless pseudogene end to end), unsure for P (a pseudogene read
on the gene needs a spliced alignment across the gene's introns; I expect fewer, 30-60%). M2 holds. M3 holds: most of
these records score 0.5-0.9 of the best. If M1 and M3 hold, the slim mode is not "equivalent" for member discovery and
row 1902's wording needs that exception.

## Result (2026-11-22, job 4495091, `r1909_silent_members`): M1 HOLDS both ways, M2 HOLDS, M3 HOLDS
Reads: G 25,680 (317,719 records, 292,039 secondaries), P 8,340.
| silent loci, reached from the reads of their expressed relative | n | within 2% of the best | >= 0.9 | >= 0.8 | >= 0.5 | any score | slim mode `-p 0.8`, any score |
|---|---|---|---|---|---|---|---|
| processed pseudogenes, from their gene's reads | 417 | 7 | 14 | 31 | 119 | 292 (70%) | 20 |
| genes, from their pseudogenes' reads | 150 | 1 | 2 | 7 | 64 | 122 (81%) | 13 |
| silent copy of a tandem pair | 150 | 49 | 74 | 82 | 118 | 141 (94%) | 84 |
| silent copy of a dispersed pair | 128 | 10 | 35 | 54 | 82 | 96 (75%) | 53 |
- M1 holds: 70% of the silent pseudogenes and 81% of the silent genes are reached by >= 2 reads of their relative.
- M2 holds: near-tie records reach 1.7% of the pseudogenes and 1 of 150 genes.
- M3 holds: the slim mode reaches 20 of the 292 pseudogenes and 13 of the 122 genes.
- 173 of the 292 pseudogenes are reached only by records scoring under half of the read's best.
- Places on no planted locus reached by >= 2 reads through secondaries: 1,061 (G) and 327 (P) with `-p 0.1`, 312 and 33
  with `-p 0.8`. No truth for them: unlisted relatives and noise are not separated here.
Predictions: G right; P wrong (I said 30-60%, it is 81%); M2, M3 right.
Consequences: the advisor's statement holds on planted reads. The low-scoring secondaries carry the existence of silent
relatives, and the slim mapping mode (row 1902) removes it: "equivalent" there is for assembling expressed copies only.
Not shown: how many of the reached places are real members on a real sample (the specificity), which needs a member
list as truth.
