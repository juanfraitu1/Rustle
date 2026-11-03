# PREREG PROJECTED MODELS, MULTI-COPY LAYER (register 1519, 2026-11-03)

**Written before any number of the arm.** User (2026-11-03): "multi-copy assembly as a different layer is a non-negotiable for my
thesis ... the idea is the 'borrow info across members' so things with low coverage or not expressed are still found, probably the
multimappers are key" -> "ok lets do it".

## Scope (first, because of r1441/1455/1457/1462/1465/1467/1468/1469)
The family catalog is arm A's `x.fam.copies.tsv` (r1495 protocol, `/anvme/workspace/v138da13-rustle/fusedgw/A`), UNCHANGED. Nothing here
discovers families, edges or copies. Projection adds TRANSCRIPT MODELS at existing members; no pair-level score can move. Unexpressed
copies absent from the catalog (`--tied-seed` TSFAM loci, candidate_copies output) are v2, not this row. Every catalog copy is a
`locus_rep` with >= 1 model, so "copies gaining a first model" is 0 by construction; reach is measured on isoforms at LOW-READ copies.

## Mechanism (`bench/project/project_models.py`; structure is borrowed, reads are never pooled; r1479 FAMASM pooled reads and FAILED)
1. Templates: every unique MULTI-EXON intron chain among a family's member models (any member), as spliced sequence in transcript
   orientation (`liftmode.templates`' recipe). Single-exon templates are out (r1518 showed what those are made of).
2. Targets: every member on the panel contigs, locus span padded by the family's longest template span.
3. `minimap2 -ax splice:hq -uf --eqx -N <members> -p 0.1 --secondary=yes` templates vs targets (the isocon_targets flags); best AS
   per (template, member) within the template's own family. Accept: same strand as the member, template coverage >= 0.9, identity
   >= 0.85, >= 1 intron at the target. A chain already among the member's own models is skipped (`already`).
4. Tier by the reads at the target, exact intron-chain match: `own` = >= 1 PRIMARY record; `tied` = no primary, >= 1 SECONDARY record
   whose molecule's primary lies on another member of the same family (r1514's --primary-on-family definition; a secondary from
   outside the family is foreign evidence and does not count); `none` = neither.
5. Output: a separate layer (`projected.gtf`, attrs tier/template/from_copy/cov/ident/motif), one GTF per tier, `base.gtf` = the arm's
   own transcripts at the same target copies (the SQANTI control).

## Panel
A119b, families with >= 2 members and a member on chr13 or chr20 (285 families, 1,197 copies on those contigs, 201 with <= 2 reads).
SQANTI3 QC of base.gtf and each tier GTF vs the CAT Liftoff annotation on chr13+chr20 (built from the slim gff3 like chr20_ref_cat).
Wall <= 30 min on one node.

## Bars (fixed now)
1. Reach: >= 50% of the 201 low-read target copies gain >= 1 projected chain of any tier.
2. Precision, `own` tier: artifact share (genic, antisense, intergenic, genic_intron, fusion) <= base's at the same loci, and
   FSM+ISM share >= 0.9 x base's.
3. Precision, `tied` tier: the same two conditions as bar 2 (a tier with < 20 transcripts is reported, not judged).
4. Mechanism: canonical splice motifs (GT-AG, GC-AG, AT-AC) at >= 95% of projected introns, all tiers.
5. `none` tier: reported descriptively only (FSM/ISM/artifact vs CAT). CAT is itself a Liftoff projection from GRCh38, so an FSM there
   means two projectors agree, not that the isoform exists; its truth is expression in another sample (KB3781/HG002 are SPENT for
   families; a new sample is the user's call).
PASS = bars 1, 2, 4 and (3 when judged). PASS -> genome-wide run, then the Rust port as a layer. FAIL -> stays a bench prototype.

## Wall note (before any arm number; three runs hit the 30-min timeout inside minimap2)
The chr13 mega-families (1,200 templates x 23 targets of tandem repeats) took 12-14 min each: the repeat minimizers (average
occurrence 105) make the splice chaining quadratic. `-f 2000` (ignore minimizers occurring > 2,000 times in the family's targets)
cuts 100 templates of MCL0 from 28.7 s to 2.5 s with the same records and 6% MORE (template, target) pairs at coverage >= 0.9.
Adopted; the bar stays 30 min.
