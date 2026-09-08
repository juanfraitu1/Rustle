# PREREG — the annotation-degradation ablation for O1 (2026-09-07)

**Written before any arm was run.** md5 of this file recorded in `mcl_ann/adj/ablation/PREREG.md5`.

## The objection this answers
"The SD-core definition is over-reliant on the annotation." Two claims are already measured and are NOT the
subject of this pre-registration: (i) the annotation carries no family information — the 25 NPIP members carry
25 distinct gene symbols, 24 of them `LOC*`, 8 distinct descriptions whose largest group is "titin-like" (9),
and genome-wide 95.2 % of members of 1,021 product-defined families are `LOC`-named; (ii) without any
annotation the seed-free mode finds every locus and delineates none (nodes 2.45–16.62× too long, chr1 18/40
vs an all-singletons 17/40). What is NOT yet measured, and is the subject here: **how good does the
annotation have to be** for the definition to recover the family.

## Design
One substrate (gorilla 3-contig: NC_073241.2, NC_073242.2, NC_073244.2), one truth (the 26 LCR16a loci,
`mcl_ann/adj/size/lcr16a.bed`), one scorer (`bench/o1_ablation_score.py`, written before the arms run).
**Only the GFF varies.** Fixed in every arm: the genome FASTA, SEDEF (`GGO_sedef_final.bed` — an
annotation-independent instrument), the reads (`npip3.bam`), the repeat library, and every `mcl_families`
threshold (0.70 / 0.30 / 300 bp / size 2 / "half" / 3 reads / 50 kb).

Each arm regenerates the whole pipeline from its degraded GFF: gene-span FASTA (`bedtools getfasta`) →
all-vs-all `minimap2 -x asm20 -c -X -N 50 -p 0.1 -t 4` → `mcl_families --min-exonic-bp 1 --sedef --bam
--fasta`. Degradations are deterministic (`numpy.random.default_rng(1337)`), applied to every gene and
pseudogene record on the three contigs, exons clipped to the new span and records shorter than 500 bp dropped.

| arm | degradation |
|---|---|
| A0 CONTROL | GFF unchanged, full regeneration |
| A1 JIT5 | each record's start and end shifted independently by U[-5 kb, +5 kb] |
| A2 JIT20 | the same at U[-20 kb, +20 kb] |
| A3 DROP50 | 50 % of records deleted at random |
| A4 SPANONLY | every record's exon structure replaced by one exon spanning the gene |
| A5 TRUNC50 | every record replaced by the middle 50 % of its span |

## Metrics
- **M1 fragmentation** (selection-free): the number of clusters holding ≥1 truth locus.
- **M2 sensitivity**: truth loci recovered by the dominant cluster, and by any cluster.
- **M3 precision**: members of the dominant cluster that are truth / its members.
- **M4 bipartite**: Hungarian 1:1 core-hull matching — median truth coverage, in-band 0.5–2× fraction.
The dominant cluster is the one with the most truth-overlapping members (ties broken by member count). It is
truth-selected, which is why M1 is reported alongside it and is not.
For A3, M2 is reported twice: over all 26 loci, and over the truth loci whose annotation record survived.

## Predictions
| # | arm | prediction |
|---|---|---|
| **P0** | A0 | the dominant cluster's member set is **exactly** the shipped 25 members of `rna_units_v12` MCL1. **If P0 fails the ablation is void** and nothing below may be quoted. |
| **P1** | A1 | M1 = 1; M2 (dominant) ≥ 24/26; M3 ≥ 0.95 |
| **P2** | A2 | M1 ≤ 2; M2 ≥ 20/26; M3 ≥ 0.90 |
| **P3** | A3 | M3 ≥ 0.95, and M2 over surviving records ≥ 0.90 — the family must degrade by losing deleted members, not by gaining false ones |
| **P4** | A4 | the member set differs from A0 by ≤ 2 members (the core rule reads genomic spans, not exons) |
| **P5** | A5 | M3 ≥ 0.90. M2 is deliberately NOT predicted: this is the arm expected to lose members. |

## Interpretation fixed in advance
- P1 and P3 holding ⟹ the annotation supplies **approximate intervals**, not family information, and MCL plus
  the core rule earn their keep. That is the sentence the advisor gets.
- P5 failing ⟹ **gene-model extent matters more than gene-model presence**; this is a limit of the definition
  and must be reported as one, not omitted.
- Any arm that IMPROVES on A0 is a warning that the shipped GFF is a poor node, not a success.
