# PREREG FAMILY PANEL ACROSS TOOLS (register 1521, 2026-11-03; descriptive, predictions stated)

**Written before any number.** User (2026-11-03): "focus on checking in stringtie, isocon, flair what did they emit for each of the
transcripts and isoforms for multi-copy gene families (at the single, base mode), did it produce all the members, did it only find
some, did it find none, did it find merged artifacts? Since the secondary alignments seem necessary for this and stringtie and flair
don't use them I expect to see poor performance there."

## Inputs
Truth: famtruth primates-level components over CAT (420 families, 1,487 genes). Tools, A119b genome-wide, each in its single-sample
mode: ours = `fusedgw/A/x.gtf` (r1495 arm A, RUSTLE_ASSEMBLY_LAYERS=off, no family stage), StringTie `-L` and FLAIR collapse from the
lab recipes (`benchmark_collapse/A119b.{stringtie,flair.isoforms}.gtf.gz`). IsoCon is NOT in this row: it was only run on the chr20
single-copy pools (r1482), there is no genome-wide IsoCon output.

## Definitions (fixed; `bench/famcmp/family_panel.py`)
Expressed member = >= 2 primary records on the gene's strand overlapping its exons. Hit = >= 1 spliced model (>= 2 exons) on the
strand overlapping the exons. Families judged = >= 2 expressed members, on the expressed members only: complete / partial / none.
Fused = a tool model overlapping exons of >= 2 members. Member hit rate reported for members whose primaries are >= 50% MAPQ 0
(the multimapper-dominated copies) vs < 50%.

## Predictions (the user's, written here so the result can contradict them)
1. ours complete-family count > StringTie > FLAIR; none-family count the reverse.
2. The member hit-rate gap between ours and the others is concentrated in the MAPQ0-dominated members.
3. Fused families: ours <= StringTie and FLAIR (the span-overlap definition catches readthrough and tandem merges alike).
No decision bar: nothing ships or changes on this row.

## Amendment 1 (before any number; the first run was cancelled unseen)
User: "we cannot trust primary alignments at multi-copy gene families ... it could have been one of the secondary alignments that
were not used for assembly with stringtie, flair". A primary-only expression definition judges primary-only tools on exactly the
copies their input can see. Two definitions are therefore reported side by side: P = >= 2 molecules with a primary record on the
member; T = >= 2 molecules with any placement there (primary, or secondary at AS >= 0.98 x the molecule's genome-wide best, the
good pool's tie set). T-only members (reachable only through a tied secondary) get their own hit-rate column. Prediction 4: the
tools' hit rate on T-only members is near 0 for StringTie and FLAIR; ours is above 0 only through the good seeding pool. No PSV
certificate is used (arm A has no assign stage), so neither definition is the truth of which copy is expressed; T is the
upper bound of what any tool could see, P the lower.
