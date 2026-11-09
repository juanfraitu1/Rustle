# PREREG row 1595: a family exon map (which exon of each copy matches which), from the shipped families' own alignments

Written 2026-11-09, before any column is built or counted. Advisor context: Soto 2025 defines duplicated genes by SD98 +
shared exons. The families stage now keeps an edge when the two loci share exons (row 1593, `--min-exon-share 0.5`, default
since e08aab2b), but it only answers yes or no per pair. This row turns the same exon matches into a per-family map:
exon "columns" shared across copies, and which copies carry each column.

**Not a family change.** Families are read, never changed. Not a rerun of `bench/family_container.py` either: the
container marks each member's exon blocks core or accessory (does it align to any exon of any other member). New here is
the cross-member correspondence (exon k of copy A = exon j of copy B = ...) and the classes of those columns.

## Inputs (A119b, the shipped default)

- Families: `/anvme/workspace/v138da13-rustle/port1593/Ab.*` = damped (`--sd-damp 0.25`, CHM13 SEDEF) share >= 0.5 at
  I 2.0, the pipeline default. 1,395 families, 6,495 member loci, all families >= 2 members.
- **Member** = one row of `Ab.fam.copies.tsv` (the node's `locus_rep`); its exons = that row's `exons`. Its PAF name is
  the node key `CONTIG:START-END` (`Ab.fam.clusters.tsv`; `Ab.fam.loci.tsv` maps the 8,237 folded annotations onto it as
  `representative`). The PAF aligns node sequences only, so a folded record never appears as a separate member, and two
  exons of one member in one column are a collision, never false sharing.
- Alignments: `Ab.fam.loci.gated.paf` (CIGAR `cg:Z`).

## Definitions (fixed now; `bench/denovo_eval/exonmap.py`)

- **Exon match**: exonshare.py's rule, reused (`blocks_from_cigar`, `seg_map`, `project`): exon x of locus A matches exon
  y of locus B when x's bases projected through the record overlap y by >= 50% of the length of x and of y.
- **Which record**: for each member pair (A, B) of a family, the one PAF record with the largest SHARE (ties: earlier in
  the PAF). This is the record that admits the pair under the shipped gate. Pairs with no record contribute no matches.
  Self-records and records between two loci of different families are not used. Known effect: where minimap2 split a
  pair into several records (an inversion or a large indel), exons outside the chosen record go unmatched. That lowers
  triangle closure and inflates private columns.
- **Columns**: connected components of the graph whose nodes are (locus, exon) of one family and whose edges are exon
  matches. Every exon is a node, so an unmatched exon is a column of one.
- Column classes: **core** = holds an exon of every member; **accessory** = exons of >= 2 members, not all; **private**
  = exons of 1 member only.
- **Collision**: a non-private column holding >= 2 exons of the same locus.
- **Triangle**: loci A, B, C of one family with exon matches a-b (A-B record) and b-c (B-C record), and an A-C record
  exists. **Closed** iff the A-C record matches a with c. Counted over (a, b, c) exon triples, every triangle once.

## Bars

- **STOP, no verdict**, if fewer than 1,000 triangles exist (no test of transitivity).
- **PASS** iff collision rate <= 0.10 of non-private columns AND triangle closure >= 0.80.
- FAIL either way: the matches are pairwise facts but do not chain into coherent columns. Only pairwise exon matches are
  reported, and no column product ships.

No dev or held-out split: nothing is fitted (shipped families, fixed 50% rule, no threshold chosen from these data).
Tandem exon duplications are legitimate collisions; the bar allows for them, and colliding columns are listed.

No external truth: CAT slim GFF has no CDS, and blastn (for an independent exon-to-exon truth on CAT paralogs) is not on
this host. The row makes no claim about whether the columns are biologically right, only that they are coherent.

## Reported, no bar

- Per family (`r1595/map.tsv`): members, columns per class, core exons per member, and per column the member exons it
  holds, with their lengths (length ratio max/min) and coordinates.
- Genome-wide: core / accessory / private column and exon counts; families whose members share no core column.
- Agreement with family_container.py's per-member call on the same families: share of container-core blocks in a
  non-private column, and of container-accessory blocks in a private one (the container uses every record, so some
  disagreement is expected).
- The same counts by family size (2, 3-5, >= 6 members).

## What follows

- PASS: `exonmap.py` becomes an opt-in families-stage product (PREFIX.fam.exonmap.tsv). Next is the per-family graph:
  core columns as the shared path, accessory columns as bubbles, feeding `family_gfa` (thesis variation-graph model).
- FAIL: register; only pairwise exon matches are reported.
- STOP: register the triangle count; no verdict.

Outputs: `/anvme/workspace/v138da13-rustle/r1595/`.

Implementation note (written before the run): members are the clusters.tsv node keys. Their loci.gff3 exons equal the
`copies.tsv` exons for 6,463 of 6,495 copies; the other 32 are copies whose span lies in two nodes by containment, so
they could not be matched.

## Outcome (2026-11-09, job 4453983, 31 s): PASS

171,206 same-family records; 108,936 member pairs with a record.

| class | columns | exons |
|---|---|---|
| core (every member) | 2,512 | 9,219 |
| accessory (>= 2 members, not all) | 1,227 | 4,826 |
| private (1 member) | 2,810 | 2,810 |

- Collision: 18 / 3,739 non-private columns = **0.0048** (bar <= 0.10).
- Triangle closure: 14,504,715 / 14,632,737 = **0.9913** (bar >= 0.80). **PASS** as registered.
- Families with no core column: 136 of 1,395.
- By family size (columns):

  | members | families | core | accessory | private |
  |---|---|---|---|---|
  | 2 | 735 | 1,623 | 0 | 1,283 |
  | 3-5 | 458 | 684 | 630 | 933 |
  | >= 6 | 202 | 205 | 597 | 594 |

- Container-style agreement: 14,045 / 15,201 container-core exons are in a non-private column. 1,654 / 1,654
  container-accessory exons are private. The second figure is close to tautological: an exon with no aligned base on
  another member cannot be matched.

Caveat: triangles are counted over exon triples, so large families of near-identical copies dominate the 14.6 M.
Closure was not registered by family size and is not reported by size here.

Per column: `r1595/em.columns.tsv`. Per family: `r1595/em.map.tsv`.
