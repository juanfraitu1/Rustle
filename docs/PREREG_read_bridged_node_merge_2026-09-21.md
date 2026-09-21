# PREREG — read-bridged node merge (RBM)

**Written 2026-09-21, before any scoring of the rule.** Goal: make de novo nodes more faithful to
actual family members by repairing the split-locus pathology with direct read evidence.

## Diagnosis this rule responds to (chr16 de novo, A119b vs CHM13)

| observation | value |
|---|---|
| annotated genes owning >=1 de novo locus | 967 |
| owning >=2 (SPLIT) | 404 = 41.8% |
| adjacent owned-locus pairs | 1,193 |
| pairs whose spans OVERLAP (gap <= 0) | 641 = 53.7%, median gap -1,575 bp |
| junctions used by >1 de novo locus | **0 of 13,685** -> junction-disjointness is BY CONSTRUCTION |
| absorbable mono fragments (mono + spliced, same strand) | 12 = 0.9% (confirms register 487: class is empty) |
| intra-gene same-strand spliced-spliced exon-overlapping pairs | 728 |
| ... overlapping only PARTIALLY in exon (no exact shared exon) | 685 = 94.1% -> the break is MID-EXON |
| ... bridged by >=1 primary read (`-F 2308`, >=25 bp on each side's UNIQUE exons) | 508 = 69.8% |
| ... bridged by >=1 **MAPQ-60** read | 473 = 75.1% of 630 testable |
| MAPQ-0 share of the 65,240 bridging reads | 16.2% |

So the split pieces are one transcriptional unit under uniquely-mapping read evidence, not paralogs.

## Relation to refuted prior art (checked before writing)

- **Register 331** refuted linking split halves via **supplementary alignments**: 37 spanning reads,
  "every one at a single read, indistinguishable from a chimeric artifact". RBM uses **primary**
  alignments only (`-F 2308` excludes supplementary) and requires a read COUNT; median bridge support
  here is 21 reads. Different evidence, and 331's stated cause of death does not apply.
- **Register 487** (locus-stitching of orphan de-novo fragments) is a no-op because the class is empty.
  Re-measured here: 12 pairs = 0.9%. RBM does not target orphans; it targets exon-overlapping spliced pairs.
- **Register 372** (span-containment fix for split loci) fixed 1 of 4 — RBM is not a span rule.
- **§6p1's 5 kb same-strand merge** is the incumbent and is a DISTANCE rule; it consolidated 42.7% of
  loci but dropped dominant-cluster coverage 20/21 -> 17/21 (register 879). 53.7% of the split pairs
  here overlap in span, so distance was never the discriminator. RBM replaces distance with evidence.

## The rule

For each primary read r (`-F 2308`, MAPQ >= Q), let H(r) = the set of de novo loci for which r covers
>= 25 bp of exonic sequence UNIQUE to that locus (not shared with the other locus of the pair).
Every unordered pair in H(r) receives one bridge vote. Merge two loci (union-find, transitive) iff:

    same strand  AND  bridge votes >= N

No distance term, no span-overlap term, no threshold on sequence identity.
Parameters swept: N in {1, 3, 5, 10}, Q in {0, 1, 60}. Incumbent = no merge.

## Scorer (fixed-universe, one-to-one, collisions are misses)

Three traps bind this and are designed out:
- *"never judge a node-definition change on node-level metrics"* -> loci-per-gene is DIAGNOSTIC ONLY
  and is not an endpoint.
- *"a copy->node mapping that lets several truth copies share one node rewards node-merging"* ->
  the assignment is ONE-TO-ONE and collisions count as misses for BOTH genes.
- *node-shrink mirror (register, 09-18)*: merging GROWS nodes, so a grown node clears overlap floors
  against more targets. Guarded by the collision rule and by the false-merge bar below.

- **Universe (FIXED at baseline, never recomputed):** every annotated chr16 gene owning >= 1 de novo
  locus in the UNMERGED output. Fixed so the merge cannot move its own denominator (register 770).
- **Assignment:** each truth gene takes the locus with maximum exonic overlap. If two truth genes take
  the SAME locus, that is a collision and BOTH are scored incorrect.
- **PRIMARY ENDPOINT — per-copy correctness:** a gene is correct iff (a) no collision AND (b) >= 50% of
  its assigned locus's exonic bp fall inside that gene's span. Reported as correct / |universe|.
- **Reported alongside:** false-merge rate (loci whose exonic bp hit >= 2 distinct annotated genes at
  >= 100 bp each), junk loci (overlapping no gene at >= 100 bp), merges applied, median loci per gene.

## Bars — stated before looking

1. **PRIMARY (chr16, development):** per-copy correctness must rise by **>= 2.0 percentage points**
   over the unmerged baseline. Below that the rule is not worth a pipeline change.
2. **GUARD A:** false-merge rate must not rise by more than **1.0 percentage point**.
3. **GUARD B:** NPIP dominant-cluster coverage must stay **>= 20/21** — the exact number the §6p1
   distance merge broke (17/21). A merge rule that breaks it again is refuted regardless of endpoint 1.
4. **HELD-OUT — chr19, run LAST, never inspected during development.** Per-copy correctness must also
   rise there, by any positive amount, at the SAME (N, Q) chosen on chr16. **If chr19 fails, RBM is not
   adopted and this document records it as refuted** ([[feedback_hold_a_substrate_back]]).

Parameter choice is made on chr16 ONLY, then frozen before chr19 is touched.

## What a negative result means

If bridged merging does not clear bar 1, the split-locus pathology is not repairable by read evidence at
the node-construction stage, and the 508 bridged pairs are recording something other than a construction
defect. That is a publishable negative and goes to the register either way.

---

## OUTCOME (appended 2026-09-21, after scoring)

**REFUTED.** Bar 1 required >= +2.0pp per-copy correctness on chr16; every arm was negative and monotone
in merge aggressiveness (Q>=60: N>=1 **-26.03pp**, N>=3 -15.43pp, N>=10 -7.75pp, N>=25 -3.91pp), with
collisions rising 214 -> 283 and the false-merge RATE rising 26.6% -> 33.6%. The post-hoc
exon-overlap-restricted variant was also negative throughout.

Bars 2-4 were never reached: **chr19 was not run**, because the rule failed on development. chr19
remains unexamined and available as a held-out substrate for the next rule.

Cause: the bridges are readthrough transcription, not split loci — at N>=1 only 49.3% of bridged pairs
share a home gene, at a median 30 kb separation. Full result and the diagnosis reversal it forced:
`docs/NODE_CONSTRUCTION_OVERMERGE_2026-09-21.md`. Register rows 937-940.
