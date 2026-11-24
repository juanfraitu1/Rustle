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

## Amendment before any tool was scored (2026-11-24)
The definition as registered also called genes copies when they only share a transposon (seen in the gene list itself, e.g.
ABHD12B / LRRK1, MIR4713HG / TRIM56; human 3,008 genes). Added: the aligned part of A's transcript must hold >= 50% of the
transcript's non-repeat bases (upper-case in the soft-masked genome) and at least 100 of them. A first form (repeat bases
turned to N before aligning) was dropped because it lowered the identity of true near-identical copies and lost TBC1D3
(kept in `sets_maskN/`; the registered form in `sets_unmasked/`). Final sets: human 2,004 genes (1,200 near-identical),
gorilla 2,621 (1,452), chimpanzee 2,330 (1,257), bonobo 2,363 (1,421), Bornean 2,096 (919), Sumatran 2,187 (992).
Chosen by looking at gene names only; no tool output had been read against any of the three forms.

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

## Outcome (2026-11-24, `r1941_multicopy_pages/{data.json,claims.py}`)
Tier 0.90, 12 samples (tier 0.98 in brackets):
- C1 HOLDS: more precise than StringTie and FLAIR on 12 (12); more correct chains than StringTie on 8 (6).
- C2 FAILS BY ONE at the registered tier: with tied secondaries Rustle has the most correct chains on 9 (10); FLAIR on 3 (2).
  With primary alignments only Rustle has the most on 1 (2): FLAIR 7 (5), StringTie 4 (5).
- C3 HOLDS: tied secondaries add 3.8% to 29.5% correct chains at multi-copy genes (8.9% to 51.2%) against 0.1% to 1.0%
  genome-wide, larger on 12 of 12. Cost in precision at these genes 0.006 to 0.049.
- C4 HOLDS: most hidden chains found, guided, on 9 (8); highest precision outside the guide on 6 (8), bambu on the others.
My prediction: verdicts right (C1, C3 hold, C2 fails), reason wrong: I expected FLAIR to lead on most samples in the block
with secondaries; Rustle leads on 9.
Limits: same assemblies as the genome-wide pages (Rustle guided of 2026-11-21, SQANTI3 on the assemblies of 2026-11-17);
transcripts with no reference gene are outside; RefSeq is the truth, weaker for the apes; "found" is an annotated chain
assembled, not proof that this copy is the expressed one.
Pages (private): comparison ba43d097-d3e0-4cfd-84bf-86e5fdc16d65, SQANTI3 b35489d4-f6d1-42dd-b97b-9123c7df8d0b, gffcompare
codes d75b8770-8bb6-4886-b875-b74f3b43eb13; build `pages/multicopy_2026-11-24/build.py` in the workspace.

## Current state
Run; C1, C3, C4 hold, C2 fails by one sample at the registered tier.
