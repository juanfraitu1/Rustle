# PREREG: the three tool-comparison pages restricted to multi-copy genes (row 1941)

Date 2026-11-24. Workspace `r1941_multicopy_pages`. Asked by the user: the assembler comparison, the SQANTI3 category page and
the gffcompare class-code page, "but only comparing multi-copy gene families; ideally that is where I expect Rustle to do
better than all others".

## What is measured
Nothing is assembled again. The assemblies, gffcompare runs and SQANTI3 classifications behind the three pages are reused
(`artifact_refresh_2026-11-21/runs/ID/score/g{1,2,3}`, `artifact_refresh_2026-11-16/sqanti`); only the set of genes changes.
Consequence, to be stated on the pages: Rustle guided is the version of 2026-11-21 (before row 1934's change), and the SQANTI3
page uses the assemblies of 2026-11-17, as the pages they mirror do.

## Multi-copy gene, defined from the annotation alone, the same way in every species
Per RefSeq annotation: the longest transcript of every gene (spliced sequence, `gffread -w`), all against all with
`minimap2 -c -x asm20 -X -N 50 -p 0.1`. Gene A is multi-copy when its transcript aligns to the transcript of a gene B whose
locus does not overlap A's, over >= 50% of A's transcript, at identity >= T (1 - `de`). Two tiers, both reported:
T = 0.90 ("related copy") and T = 0.98 ("near-identical copy"). Pseudogenes and non-coding genes are in (user rule: families
include them). Fixed before any tool is scored on the set; no tuning of T or the coverage on the result.

## Scores
- A tool transcript belongs to the set when gffcompare assigns it to a multi-copy gene (`ref_gene_id` of the tmap). Transcripts
  with no reference gene (code `u`) are therefore outside; said on the page.
- Page 1: distinct multi-exon intron chains; correct = equal to a chain of a multi-copy gene's transcript; sensitivity over the
  union of every tool's correct chains (as on the genome-wide page); precision = correct / chains at multi-copy genes.
  Blocks kept separate: assembler alone (primary alignments), with tied secondaries, guided (hidden chains found and precision
  outside the guide, hidden as in `pages/gffcompare_classes_2026-11-23/agg.py`).
- Page 2: SQANTI3 structural categories of the transcripts whose `associated_gene` is a multi-copy gene.
- Page 3: gffcompare class codes of the transcripts at multi-copy genes, with and without a guide.

## Claims, tier 0.90, 12 samples (tier 0.98 reported beside, not a second chance)
- C1 Assembler alone: Rustle more precise than StringTie and FLAIR on >= 10 samples and more correct chains than StringTie on >= 8.
- C2 (the user's expectation) Rustle has the most correct chains of the five tools on >= 10 samples, in the block with tied
  secondaries.
- C3 Tied secondaries matter more here than genome-wide: the relative gain in correct chains (with secondaries / primary only)
  is larger at multi-copy genes than genome-wide on >= 10 samples.
- C4 Guided: Rustle finds the most hidden chains at multi-copy genes on >= 8 samples.

## My prediction
C1 holds. C2 FAILS: FLAIR will have the most correct chains on most samples here as genome-wide (register row 1670: parity with
StringTie / FLAIR in family members on the default assembly). C3 holds. C4: about as genome-wide; not sure.
The pages show whatever comes out; a failed C2 is written in the headline, not left out.

## Current state
Registered, not run.
