# PREREG row 1693 (2026-11-13): soft mask at the edge level, step 0

User go 2026-11-13 ("Edge-level mask"). Row 1690 tested the soft mask as a family-level skip list only and left the edge
level as a separate row: a family link should rest on aligned sequence that is not soft-masked. This row is the step 0
count on existing products. No new alignment, no re-clustering, no default changes.

## Question
Of the false gene pairs the default links through an edge, how many rest on masked sequence, and how many true pairs
would the same rule cut?

## Data and measure
A119b default run (`r1662_edgemiss/asg_def/a1`, CHM13 soft mask). Gene pairs, truth and stages exactly as
`bench/denovo_eval/edgemiss.py` (expressed protein-coding genes, simiiformes truth, `--hits` from pairscore). For every
gene pair whose best locus pair reaches stage `edge`: the exon-to-exon aligned bases of that locus pair's records, both
sides, read from the genome FASTA; LOW = share of those bases that are lowercase. An edge is MASKED when LOW >= 0.9
(the threshold row 1690 froze for exons). Script `bench/denovo_eval/edgemask.py`.

## Bars (written before the count)
GO to a re-clustering row only when both hold:
- R: MASKED edges are at least 20% of the false linked pairs at stage `edge`;
- T: MASKED edges are at most 1% of the true linked pairs at stage `edge`.
Otherwise STOP: the edge-level mask is closed as a precision lever on protein-coding pairs. Reported without a bar: the
same two shares at LOW >= 0.5, and the LOW distribution of both classes.
Limits stated up front: only edges of the base locus graph (extra representatives' edges are not in `a1.fam.loci`);
non-coding and unannotated loci are outside the truth universe, which is where row 1690's satellite families sit.

## Outcome (2026-11-13, job 4470868): STOP
| class (stage `edge`, linked) | pairs | LOW >= 0.9 | LOW >= 0.5 |
|---|---|---|---|
| false | 103 | 2 (0.019) | 12 (0.117) |
| true | 1,529 | 10 (0.007) | 184 (0.120) |

- R fails: 2 of 103 false linked pairs rest on masked sequence (bar 20%). T holds (0.7%, bar 1%).
- At LOW >= 0.5 the two classes are masked equally often (0.117 vs 0.120): the mask does not tell false links from true
  ones among protein-coding pairs. Median LOW is higher for true pairs (0.38) than for false ones (0.09).
- Closed: the edge-level mask is not a precision lever on protein-coding pairs. The soft mask's use stays the family-level
  skip list of row 1690 (satellite / LINE families outside the truth universe). No re-clustering row is opened.
Table: workspace `r1693_edgemask/a.edges.tsv`.
