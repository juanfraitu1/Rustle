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
