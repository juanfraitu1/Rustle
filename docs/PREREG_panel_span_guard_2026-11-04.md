# PREREG r1535: panel span guard (scorer correction) + the mono layer measured genome-wide (2026-11-04)

From r1534: the base emits no model at 51 single-exon truth members (class A) and the panel counts any model (Amendment 2);
46 of the 53 class-B4 members have >= 50% of their panel reads' aligned bases inside the gene span (unspliced, own), 7 are
neighbours' reads (~/r1534/check_b4_span.out). No own-read rule is registered: B4 is class A's sibling, not a ruler artifact.

## Step 0 (done before this file: ~/r1534/check_blob_hits.out)
Expressed members (T >= 2) of the r1521 panel whose only base hits are single-exon models with span > 2x the gene span:
**39** (TAF11L3/4/5/8/9/10 x 25 in one 113 kb blob, gene 597 bp); 885 members keep a hit under the guard, 212 have none.

## Definition change D1 (scorer, applied to every tool, base included)
`PANEL_MONO_SPAN=2`: a single-exon model counts as a hit (and towards fused) only if its span <= 2 x the member's exon span.
Adoption bar (as r1530: a correction adopted only if it removes what it claims): >= 0.9 of the hits it removes come from models
with span > 10 x the gene span. Base before/after is reported; the after figure replaces the r1521 base numbers if adopted.
Numbers move DOWN for the base (up to 39 members, fused families may drop too); no assembler change.

## Arm L1: the layered assembly, genome-wide (measurement only; defaults untouched)
One A119b genome-wide `assemble` with the 05:01 binary and RUSTLE_ASSEMBLY_LAYERS=on (the default): writes x.gtf (unlayered)
and x.final.gtf (--polish-tss-end median, --polish-tes-extend on, --mono-model cluster). Output /anvme/.../fusedgw/L1.
- Gate (before the layer is judged): x.gtf spliced chains kept >= 0.99 vs fusedgw/R1/x.gtf (chain_cat.py). FAIL -> the
  baseline moved with the binary; re-pin before judging, no bar is read.
- Judged with D1 on, full BAM + a1.molecules.tsv: base = L1/x.gtf, arm = L1/x.final.gtf.
- Bars: (1) class A members (the 51 of ~/r1534/step0_members.tsv) with a hit >= 26; (2) fused families <= base; (3) complete
  >= base; (4) ordinary-member hit rate (mapq0 share < 0.5) >= base; (5) wall <= 4 h.
- Whether the panel or the families stage should consume final.gtf is the user's decision, not decided by this row.

## D1 RESULT (~/r1535/panel_R1_guard.out; base = fusedgw/R1/x.gtf, same 1,136 members)
| definition T | complete | partial | none | fused | member hit | ordinary | T-only (n=87) |
|---|---|---|---|---|---|---|---|
| r1521 (guard off) | 209 | 94 | 21 | 49 | 0.832 | 0.829 | 0.678 |
| D1 guard on | 206 | 94 | 24 | 48 | 0.797 | 0.822 | 0.333 |
Hits removed 39; 34 of them from single-exon models > 10x the gene span -> share 0.872 < 0.9: **adoption bar FAIL as
registered**. D1 is not adopted; the r1521 numbers stay the figure. The arm L1 is still judged with the guard on, as
registered, and reported with the guard off after the numbers (labelled). Note: the T-only hit rate was carried by the
113 kb TAF11L blob (0.678 -> 0.333 under the guard); 30 of the 87 T-only members are TAF11L copies hit by one model.
Incident: `assemble --as-table` rebuilt ~/runs/human_A119b/a1.molecules.tsv in place (its header path did not match
`readlink -f` of the BAM; 6 min, same BAM, same molecule count 21,391,966); a panel run that read the half-written table
was discarded.
