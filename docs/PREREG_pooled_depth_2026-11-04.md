# PREREG r1531: depth, not rules: does pooling HG002 Kinnex reads reach the family panel's misses? (2026-11-04)

Every rule-level route on the r1521 panel is closed (r1526-r1530). The remaining lever is read depth: the multi-sample pool
mode (r1503, ships as more sensitive, not more precise). HG002 Kinnex (lymphoblastoid line, 206 M mapped records on CHM13) vs
A119b (testis, 68 M). The misses are largely testis/placenta genes (PRAMEF, TSPY, SPANX, CGB, PSG, GAGE), so step 0 asks whether
HG002 has reads there at all before any assembly is run.

## Step 0 (decides whether an arm runs; fixed before measuring)
- Member set M = expressed panel members with no base model (~/r1526/panel.genes.tsv rows base = 0 and tied_molecules >= 2;
  212 rows) UNION the 97 one-chain members (~/famcmp/assemblable.tsv, ours = 0, maxP_chain < 2, maxT_chain < 2,
  n_local_spliced_P >= 1).
- Measure: for each member, HG002 primary, on-strand, exon-overlapping reads (family_panel.py's fetch, primaries only: no HG002
  as-table), from ~/work/jfim/rustle_data/reads/HG002_kinnex.t2t.bam. Report members with >= 2 such reads, and their
  chromosome distribution.
- **Stop rule: fewer than 25% of M reachable (>= 2 HG002 reads) -> stop, record "HG002 depth does not reach the misses" as the
  ceiling, no arm.** Otherwise the arm prereg is written as an amendment below, with the contigs chosen from the reachable
  members' chromosomes, BEFORE any assembly.

## Arm (only if step 0 passes; details fixed then)
Slices of both BAMs on the chosen contigs -> pool stage (name-collision check + merge) -> two single-BAM `assemble --contig`
runs with the same binary, layers off, strand vote default, each building its own as-table: base = A119b slice alone, arm =
pooled slice. Panel judged on the FULL A119b BAM + a1.molecules.tsv (expression truth stays A119b), restricted to families
whose every gene lies on the chosen contigs (one-line filter in family_panel.py, env PANEL_CONTIGS). Bars relative to the in-job
base: complete >= base + N (N fixed from the panel size), fused <= base, ordinary and T-only hit rates >= base, base chains kept
>= 0.95, new chains' exact CAT match >= 0.05, job wall <= 2 h. PASS -> one genome-wide pooled row vs fusedgw/R1; ship = pooled
input documented as the recommended family-layer input (not a default). FAIL -> depth closed on this panel.
