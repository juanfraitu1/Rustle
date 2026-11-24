# Correlation clustering of the reads as the candidate step of the sequence test: stopped at development (row 1932, 2026-11-23)

User: "how can we apply correlation clustering to O2; isn't it what IsoCon does or an extension of it?" then "lets try the
registered experiment". IsoCon does not do correlation clustering (nearest-neighbour graph, consensus, a statistical test
per pair of candidates, sequential merging); correlation clustering is the global form of the same pairwise idea. The
measured weakness it was meant to address: `recurrence_test.py` folds a shorter read into the one longer pattern that
contains it, and about 30% of the sequences it keeps hold reads of more than one true sequence (rows 1926, 1927).

## What was built (`recurrence_test.py --cluster cc`, opt-in, self-tested)
Nodes = reads with the same covered sites and alleles. Pair score = log-likelihood ratio of one sequence against two
over the sites both cover (agree: the chance two different sequences agree there, from the allele frequency; disagree:
the chance of an error in either read, IsoCon's per-variant rates). Greedy merging of the pair of clusters with the
largest positive total, then single-node moves; clusters whose consensus agrees at every shared site are merged;
IsoCon's test then runs on the clusters as before.

## Development (seeds already spent: 1930-1932 homozygous + control, 1939-1941 heterozygous; jobs 4498739, 4498752)
| exact cells | folding (the tool as it is) | cc, first version | cc, with the consensus merge |
|---|---|---|---|
| homozygous copies (128) | 98 | 87 | 97 |
| heterozygous, against haplotypes (199) | 106 | 104 | 108 |
| heterozygous, against copies (199) | 57 | 56 | 60 |
| control, no variants (42) | 40 (2 over) | 31 (11 over) | 40 (2 over) |
| purity, homozygous (kept sequences of one true sequence) | 375 of 533 (0.70) | 390 of 566 (0.69) | 334 of 504 (0.66) |

First version: an error shared by two reads splits them off, reads not covering the site join them, the consensus
equals the neighbour's, no site differs, the test never sees the pair: 11 control cells over-counted. Fixed by the
consensus merge. Run time about 1.5 x the folding's.

## Decision: no fresh-seed test
After the fix correlation clustering equals the folding on every count (97 / 98, 108 / 106, 40 / 40) and is slightly
below it on purity, the one thing it was built to improve. A registered test on fresh seeds would be a test of "no
difference" at the cost of three seeds; not run. **STOP at development.** The option stays in the tool, off, as a
recorded arm; it is not a candidate default.
What this says: the impurity is not caused by the order of merging. A read that covers none of the sites separating two
sequences is ambiguous for any clustering; the pairwise evidence is simply absent. Better use of partial reads would
need other evidence (the read's ends, its isoform, linkage through a third read), not a different objective.
