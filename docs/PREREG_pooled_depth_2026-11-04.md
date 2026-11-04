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

## Step 0 RESULT (~/r1531/step0_hg002.py -> step0.out, step0_hg002.tsv)
Members 215 (212 no-model + 3 one-chain not already among them); reachable with >= 2 HG002 primary on-strand exon reads: **59 =
0.274 (stop rule < 0.25 not fired, narrowly)**; 0 of the one-chain members. Reachable by chromosome: chr6 9, chr19 8, chr12 7,
chrX 5, chr1 4, chr17 4, chr4 4, then <= 3 each. The reachable genes are immune/blood genes (LILR, SIGLEC, CCL3, APOBEC3B,
HLA-DQ), histones, TAS2R and OR: HG002's tissue, not A119b's misses in PRAMEF/TSPY/SPANX.

## Amendment 1 (arm fixed BEFORE any assembly; written after step 0, before any run)
- Contigs: chr1,chr4,chr6,chr12,chr17,chr19,chrX (41 of the 59 reachable members). Panel: 165 families fully on these contigs
  (genome-wide base complete 102, expressed members 471).
- Reads on the contigs: A119b 22.1 M, HG002 70.9 M mapped records.
- Products: slices via `samtools view -b` of both BAMs, `pool` stage on the slices (name-collision check + merge), then
  `assemble --contig <list>` twice with ~/rustle_target_test/release, layers off, strand vote default, no --as-table (each
  builds its own from its slice): base = A119b slice (~/.../depth/A/x.gtf), arm = pooled slice (~/.../depth/AH/x.gtf).
- Panel judged on the FULL A119b BAM + a1.molecules.tsv with PANEL_CONTIGS set (expression truth stays A119b).
- Bars (all must hold), base = in-job A119b slice: complete >= base + 5 (definition T); fused <= base; ordinary hit rate
  >= base; T-only hit rate >= base; base spliced chains kept >= 0.95 and new chains' exact CAT match >= 0.05 (chain_cat.py);
  pooled assemble wall <= 7200 s. Info only: step-0 reachable members gaining a same-strand model.
- PASS -> genome-wide pooled row next vs fusedgw/R1; ship = documented recommended input, not a default. FAIL -> depth closed.
