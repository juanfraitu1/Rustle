# PREREG — G3: sensitivity of the catalog to the edge thresholds identity 0.70 / cov_longer 0.30 / min_bp 300 (2026-09-05)
Substrates: (1) gorilla 3-contig, new defaults (fold within clusters, exon-less span, exon-to-exon, min_size 2), no
reads/SEDEF (clusters.tsv is independent of them); (2) Soto slice, same defaults, scored with soto_adjudicate.py.
Grid, one-at-a-time from the defaults: identity {0.60,0.65,0.70,0.75,0.80,0.85,0.90}; cov_longer
{0.10,0.20,0.30,0.40,0.50,0.60}; min_bp {100,200,300,500,1000,2000}. 19 runs per substrate.
Anchor metrics (gorilla; anchors from rna_bp1_p9, the audited catalog): cohesion = share of the anchor's records
in its modal cluster, for NPIP (MCL2, 44), MCL1 (48), MCL3 (36), MCL4 (29); NPIP∥LCR16u = NPIP's modal cluster
contains none of LCR16u's (MCL7, 14) records; L1 blob (MCL0, 104) modal-cluster share (dissolved = ≤ 0.10).
Soto metrics (50 % floor): detection, band-[0.90,1) precision, recall|both, family exact.
Predictions:
- P1 invariance band: every anchor's cohesion ≥ 0.95, NPIP∥LCR16u holds and the blob stays dissolved for identity
  0.60–0.80, cov_longer 0.10–0.40, min_bp 100–500 (the thresholds sit inside a plateau, like inflation §6ec).
- P2 walls: identity ≥ 0.85 or cov_longer ≥ 0.50 fragments NPIP (cohesion < 0.95) — NPIP records pair at
  0.94–0.98 identity but chimeric models have low cov_longer; min_bp ≥ 1000 removes short-exon pairs (Soto
  detection drops ≥ 5 %).
- P3 Soto band precision varies by < 0.02 over identity 0.60–0.80 (the band is above every tested floor; only
  transitive membership changes).
Reading: the defaults are reported as a point inside a measured plateau with its walls; if P1 fails at the
defaults' own neighbours (0.65/0.75, 0.20/0.40, 200/500) the threshold is a decision and is stated as one.
