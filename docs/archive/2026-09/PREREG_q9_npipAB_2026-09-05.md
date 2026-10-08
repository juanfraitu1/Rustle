# PREREG — Q9: are NPIPA and NPIPB distinct subfamilies in the gorilla catalog? (2026-09-05)
Family = rna_units_v3 MCL3 (29 loci). Labels = CHM13 landing of each locus's records (all46_chm13_labels.json):
17 land on NPIPB stems only ("B"), 3 land on NPIPA2 AND NPIPB13 ("AB"), 9 land on ABCC1/SORL1 (chimeric models,
5 kept_trimmed to a 24-kb LCR16a core, 4 dropped = core 0). No locus lands on NPIPA only.
Edges = allgenes.asm20.paf records between the 29, identity = matches/block, cov_longer = aligned/exonic length
of the longer record (mcl_families' rule). MCL (bench/mcl_annotation.py) on the 29-node weight matrix
(identity x cov_longer), inflation 1.4, 2.0, 2.8, 4.0, 6.0, 8.0, prune 1e-9.
Predictions:
- P1 identity does NOT separate AB from B: mean(AB–B) within 0.01 of mean(B–B).
- P2 coverage separates weakly at best: mean cov(AB–B) < mean cov(B–B), but the AB loci are not a graph cut.
- P3 at inflation <= 4.0 all 20 NPIP-proper loci (B + AB) stay in ONE cluster; the first split at higher
  inflation separates the ABCC1/SORL1 chimeras, not A from B.
Reading if P1–P3 hold: under the core rule NPIPA-landing and NPIPB-landing loci share the LCR16a core and are
one family; the human A/B subfamily labels do not partition the gorilla loci (they land on B stems, three on
both), so the subfamily question is a human 16p-position question, not a gorilla-catalog cut. If P3 fails
(AB separate at <= 4.0): the catalog does carry the A/B split and must report it.
