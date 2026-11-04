# PREREG r1528: the shared-exon fraction gate and catalog coverage of family members

Registered 2026-11-04 after step 0, before any arm number (user: "ok, lets try that now").

## Step 0 (measured before this file; ~/r1528/step0_cov.py, step0_why.py, step0_pairs.py, step0_gate.py)
Provenance route (assemble a member from reads O2 certifies to it) has ceiling 0 on the 28 missed T-only members: 146 of their
201 tied reads never enter the assign table, 15 of the 28 members overlap no catalog copy. So the catalog was checked first.
Of the panel's 1,136 expressed members (T >= 2), 831 overlap a catalog copy (a1.fam.copies.tsv, 2026-10-29); 305 do not, and 299
of those HAVE a base model: 265 are loci that stayed singletons after the sketch (257 spliced), 34 sit on a copy or locus of the
other strand only. Of 332 singleton expressed members, 55 have no sibling locus; 277 do. Of those 277: 137 pairs are in the gated
PAF but the best record's shared-exon fraction is < 0.6 for 120 (65 at 0.3-0.6, 55 below; CEACAM5/6, PSG3, TAS2R, LILR, PCDHA);
140 pairs are not in the gated PAF, of which 68 do not align at asm20 (PRAMEF, GAGE/PAGE/XAGE: protein-level families), 48 align
below 0.6 and 24 at or above (the alignment-time gate missed them). The gate `--min-shared-exon-frac 0.60` (r903: held-out
against Soto's families, kept) is therefore the first limit on member coverage, ahead of O2 and the assembler.

## Fixed before any number
- Arms: the families stage alone on the A119b assembly a1 (copied as x.gtf + x.families.gtf per arm), `rustle_pipeline.sh
  families --no-cache` with RUSTLE_MIN_SHARED_EXON_FRAC = 0.60 (rerun base), 0.45, 0.30 (new pipeline switch, unset = 0.60).
  The 0.60 rerun must be byte-identical to a1.fam.copies.tsv; if not, it is the base and the difference is reported.
- Scores: `bench/denovo_eval/pairscore.py --species human` (pair precision / recall, primates level), `bench/famcmp/catalog_cov.py`
  on ~/r1527/panel.genes.tsv (covered expressed members, largest family), families-stage wall.
- Bars for an arm to PASS (all): pair precision >= base - 0.01; pair recall >= base + 0.03; covered expressed members >= base + 80
  (base 831 of 1,136); largest family <= 1.5 x base (base 220); wall <= 15 min.
- PASS -> the lowest fraction that passes becomes the pipeline default and the provenance route (certified-read assembly of
  T-only members) is re-tested on the new catalog as its own row. No arm passes -> the gate stays at 0.60; row closed, and the
  coverage limit is recorded as the family layer's ceiling on this panel.

## RESULT (job 4436270; ~/r1528/run.4436270.out; arms in /anvme/workspace/v138da13-rustle/fam1528/)
Correction noted AFTER the numbers: "base 831" above was a transcription error; catalog_cov.py on the base catalog gives 631
covered expressed members (625 hit + 6 miss, the step-0 table's own numbers). The bar is base + 80 as written, i.e. >= 711.
| bar | 0.60 (base, byte-identical rerun) | 0.45 | 0.30 |
|---|---|---|---|
| pair precision (>= 0.716) | 0.726 | 0.670 FAIL | 0.662 FAIL |
| pair recall (>= 0.771) | 0.741 | 0.753 FAIL | 0.779 PASS |
| covered expressed members (>= 711 of 1,136) | 631 | 676 FAIL | 713 PASS |
| largest family (<= 330) | 220 | 220 PASS | 220 PASS |
| stage wall (<= 900 s) | 309 | 323 PASS | 326 PASS |
Predicted pairs 2,228 -> 2,451 -> 2,571; false pairs 611 -> 808 -> 870 (related-but-not-truth 524 -> 712 -> 760), so the new
edges are two thirds domain-level relatives, not Compara paralogs; families 1,565 -> 1,655, copies 7,987 -> 8,554.
**Verdict: FAIL as registered (both arms fail precision). The gate stays at 0.60; row closed.** Coverage limit of the family
layer on this panel: 631 of 1,136 expressed members are catalog copies at 0.60; raising it to 713 costs 6.4 points of pair
precision. The provenance route is not re-tested.
