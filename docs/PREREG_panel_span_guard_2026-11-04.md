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

## Arm L1 RESULT (job 4441808, 57 min; ~/r1535/run.4441808.out, panel.*, panel_noguard.*, classA_hits.out)
Gate: L1 x.gtf = R1 x.gtf chain-for-chain (219,503 kept = 1.000). final.gtf 211,392 transcripts vs x.gtf 233,011.
| bar (guard on, definition T) | base L1/x.gtf | arm L1/x.final.gtf | verdict |
|---|---|---|---|
| (1) class A members with a hit >= 26 of 51 | 0 | **15** | FAIL |
| (2) fused <= base | 48 | 48 | PASS |
| (3) complete >= base | 206 | 213 | PASS |
| (4) ordinary hit rate >= base | 0.822 | 0.829 | PASS |
| (5) wall <= 4 h | | 57 min | PASS |
**FAIL as registered (bar 1).** Hits: CTAGE15 (135-read model), PPP1R2B, GK3P, TTC30B, CTAGE6, TAF1L, RTL8B, RAB6D, HCAR2,
PCDHB10, PRR23C, TAS2R31, DDTL, H3C1, KRTAP9-2. Misses with many own reads: FAM47A 394, EIF5AL1 265, NANOGP8 66 (the mono
cluster rule needs a PAS, no internal priming and no multi-exon exon on the strand; which of these holds per miss is unchecked).
Classes B2 (20) and B4 (53) gain nothing. Guard off (after the numbers): arm T-only 0.678 -> 0.322 and member hit 0.832 ->
0.805, because the mono layer removes the 113 kb TAF11L single-exon transcript and rebuilds none at those 30 T-only copies;
class A hits 16. The layered GTF trades the blob hits for 15 real single-exon members.
After the numbers (~/r1535/classA_miss_why.out): the mono cluster rule's PAS / no-priming tests explain the big misses. FAM47A
381 same-strand unspliced reads, 365 ending in one 50 bp window, 0 primed, but only 2 reads with a canonical PAS -> no model
(the rule requires a PAS). EIF5AL1 311 reads, 215 internally primed. NANOGP8 61 / 60 primed, PPIAL4G 19 / 19, H4C3 14 / 14
(a histone; no poly(A) by biology). TAS2R30 49 / 35 primed. None of the misses overlaps a spliced model's exon. Relaxing the
PAS test would add FAM47A-type members; relaxing priming would re-admit the r1518 pseudogene fill-in. Not registered here.
