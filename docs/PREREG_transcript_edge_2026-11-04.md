# PREREG r1529: a transcript-level edge for sibling loci the shared-exon gate rejects (user goal, 2026-11-04)

User's goal text, verbatim in the register row; mechanism: pairscore.py's non-coding transcript rule (spliced sequences aligned
with minimap2 -c -x asm20, shorter covered >= 50% at >= 85% identity) applied as a families-stage edge on top of the 0.60 DNA gate.
Step 0 stop rules, fixed in the goal before measuring: fewer than 60 singleton expressed members with a passing sibling pair, or
more than 20% of the base catalog's false coding pairs passing -> stop, record the ceiling, no arm.

## Step 0 RESULT (~/r1529/step0_tx.py; spliced locus sequences = loci.gff3 exon union on chm13; one minimap2 run)
| measure | value | stop rule |
|---|---|---|
| singleton expressed members with a sibling locus | 238 (378 distinct locus pairs) | |
| sibling pairs passing the transcript rule | 17 | |
| members with >= 1 passing pair | **30** | STOP (< 60) |
| base false coding pairs (611 gene pairs -> 584 locus pairs) passing | 235 = **0.404** | STOP (> 0.20) |
Both stop rules fire. **No arm; row closed at step 0. Ceiling: 30 members, mostly recent duplicates (EIF3C/EIF3CL, GLUD2,
MICA/MICB, RBMX/RBMXL1, UTP14A/C, ZXDA/B).** The singleton pairs the DNA gate rejects are rejected by the transcript rule too:
their spliced sequences are below 85% identity or cover less than half of the shorter one, i.e. the members are diverged at the
transcript level, not just across introns.

### After the numbers (information for the scorer, no verdict changed)
Of the 235 passing "false" pairs, 205 are SAME-SYMBOL gene pairs (AMY1A~AMY1A, GAGE12C~GAGE12C, CDY1~CDY1: CAT duplicate gene ids
for identical copies, which Compara cannot name as a pair and pairscore therefore counts as false). Excluding them, 30 of 379 =
0.079 pass. The transcript rule's false rate as registered is 0.404; its rate on scorable pairs is 0.079. That changes nothing
here (the member ceiling fails alone) but says the 611 "false" pairs overstate the catalog's error: ~1/3 are identical-copy pairs
the truth cannot score (see memory cat-annotation-merges).
