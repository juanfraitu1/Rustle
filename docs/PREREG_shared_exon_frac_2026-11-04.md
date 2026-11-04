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
