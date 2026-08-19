# O1 single-outgroup rooting proof of concept

This experiment asks whether the two proposed human source contexts satisfy the
prerequisites for later provenance rooting against one ape species. It uses both phased
KB3781 gorilla haplotypes independently. Gorilla annotations are not used.

| source | linked human GOLGA-family loci | maternal synteny | paternal synteny | result |
|---|---:|---|---|---|
| GOLGA2 | 8 | TWO_SIDED_UNIQUE | TWO_SIDED_UNIQUE | ROOT_CANDIDATE_SINGLE_OUTGROUP |
| ITSN2 | 6 | TWO_SIDED_UNIQUE | TWO_SIDED_UNIQUE | ROOT_CANDIDATE_SINGLE_OUTGROUP |

Both source contexts pass the proof-of-concept prerequisites: recurrent human family
homology and a unique, two-sided syntenic gorilla locus in each haplotype. They remain
`UNROOTED`, because the queries are whole-locus proxies and the shared intervals are
pairwise witnesses rather than stable multi-locus block classes.

Of 18 audited family probes, 3 have unique two-sided synteny in both gorilla haplotypes:

- FAM_010 (GOLGA8A)
- FAM_012 (chr15:73646286-73657758)
- FAM_017 (chr5:7237265-7244593)

The other 15 probes lack the strict two-flank certificate. This is compatible with
duplication/rearrangement but is not treated as proof of gorilla absence: repetitive-region
mapping failure and assembly structure remain alternatives.

## Fixed proof-of-concept rules

- human source-to-family witness: aligned span >=1,000 bp and identity >=0.60;
- flank length: 25,000 bp on each side;
- qualifying flank anchor: coverage >=0.80, identity >=0.90, MAPQ >=20;
- two anchors must occur on one ape sequence, in compatible order and orientation; and
- collinear split anchors may be chained with query gap <=1,000 bp and target gap <=5,000 bp.

The split-anchor rule was added after an audit showed that minimap2 represented the ITSN2
right flank as two adjacent records. It is substrate-neutral, applied to every probe and
haplotype, and does not change any alignment threshold. Prebuilt asm20 indexes and FASTA
asm5 runs use different contig labels; target length plus coordinates reconcile those
labels without changing placements.

## Claim boundary

This is a positive feasibility result for single-outgroup rooting, not evidence that
GOLGA2 or ITSN2 has already been proven ancestral by Rustle. Production direction requires
stable human block classes, block-specific source/derived paths, assembly-gap checks and a
rooting certificate. Until those exist, no `DERIVED_FROM` edge is emitted.

Normative evidence: `probes.tsv`, `human_source_family_links.tsv`, `gorilla_synteny.tsv`,
`gorilla_locus_hits.tsv`, and `rooting_candidates.tsv`. Raw PAF files are retained under
`work/`. `rooting_candidate.gfa` is a typed visual projection and deliberately contains no
`DERIVED_FROM` edge.
