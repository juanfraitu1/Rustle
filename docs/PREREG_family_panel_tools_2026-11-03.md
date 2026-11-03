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

## Result under the registered definition (spliced models only), and Amendment 2 (POST HOC, after seeing it)
Registered run (`~/famcmp/run.4435467.out`): complete / partial / none / fused, P: ours 201/91/25/48, StringTie 227/69/21/57,
FLAIR 207/87/23/68; T: ours 199/95/30/48, StringTie 217/82/25/57, FLAIR 202/94/28/68. Member hit rate on T-only members (87):
ours 0.322, StringTie 0.000, FLAIR 0.011 (prediction 4 confirmed); on MAPQ0-dominated members ours 0.925 vs 0.827 / 0.759
(prediction 2 confirmed); fused ours < StringTie < FLAIR (prediction 3 confirmed). Prediction 1 CONTRADICTED: StringTie completes
more families. Two causes, found post hoc: (a) the registered "hit" required a SPLICED model, and many truth members are
single-exon genes (TAS2R, OR, TSPY, SPANX, TUBB8) where ours emits a single-exon model that did not count while StringTie's
spliced readthrough did; (b) at 3-9 primary molecules ours hits 0.34-0.47 vs StringTie 0.53-0.64 (the pass-1 floor of 3 and the
polish), at 10+ the two are within 0.08. Amendment 2: hit = any model (single-exon included); rerun reported below, labelled post hoc.
