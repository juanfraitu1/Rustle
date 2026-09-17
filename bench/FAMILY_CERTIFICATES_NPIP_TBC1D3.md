# Family certificates on NPIP and TBC1D3: where each family is an exact connected component, and with how much margin

Written by the orchestrator from `family_cert/cert/certificates.tsv` (the building agent was blocked from writing this
file; every number here is copied from that table, and the independent verifier recomputed all 210 of its rows).
Development families, one assembly, one annotation: this is a certificate about THIS evidence, not a validation of the
lattice. Nothing in `src/` was changed. Data: `/mnt/linuxdisk/home/juanfraitu/family_cert/` (`DECLARATIONS.txt` written
09:17:35, before any evidence; addenda A1-A3 dated and disclosed).

## 0. What is certified

For a set S in a graph whose edges carry a weight w, with G_x = the edges with w >= x (§0★★★ tests supply w):

- **h_join(S)** = the largest weight of an edge with exactly one endpoint in S (0 if none);
- **h_split(S)** = the largest x at which S is still connected using only its internal edges (the bottleneck of S's
  maximum spanning tree);
- **S is an exact connected component of G_x for every x in (h_join, h_split]** — and for no other x.

The interval is the margin. A family with a wide interval is a cluster of the single-linkage filtration, not an
artefact of one cut; a family whose interval is empty cannot be isolated by any cut on that axis.

## 1. Read first

| set | level | axis (shipped cut) | members present | h_join | h_split | exact interval | shipped cut inside | component at the shipped cut | outside |
|---|---|---|---|---|---|---|---|---|---|
| **NPIP** | **P protein** | union cover, e <= 1e-5 (0.30) | 21/27 | **0.128250** | **0.840580** | **(0.128250, 0.840580]** | **yes** | **21** | **0** |
| TBC1D3 | P protein | union cover, e <= 1e-5 (0.30) | 9/19 | 0.688525 | 0.900000 | (0.688525, 0.900000] | no | 5,275 | 5,266 |
| NPIP | L1 family | identity, cov >= 0.50 (0.80) | 26/27 | 1.000000 | 0.971665 | empty | no | >= 122 | 96 |
| TBC1D3 | L1 family | identity, cov >= 0.50 (0.80) | 12/19 | 1.000000 | 0.840488 | empty | no | >= 31 | 19 |
| NPIP | L1 family | coverage, id >= 0.80 (0.50) | 26/27 | 1.000000 | 1.000000 | empty | no | >= 122 | 96 |
| TBC1D3 | L1 family | coverage, id >= 0.80 (0.50) | 12/19 | 1.000000 | 0.841369 | empty | no | >= 31 | 19 |
| NPIP | L2 shared-exon unit | f_ex (0.30) | 26/27 | 1.000000 | 0.989040 | empty | no | >= 90 | 64 |
| TBC1D3 | L2 shared-exon unit | f_ex (0.30) | 12/19 | 1.000000 | 0.832037 | empty | no | >= 22 | 10 |
| NPIP | L3 >= 0.98 identity unit | w_98 (0.98) | 26/27 | 1.000000 | 0.980559 | empty | no | >= 83 | 57 |
| TBC1D3 | L3 >= 0.98 identity unit | w_98 (0.98) | 12/19 | 0.961397 | 0.867946 | empty | no | 12 | **0** |
| NPIP | L0 = P or D | binary, shipped cuts | 26/27 | 1.000000 | 1.000000 | empty | no | >= 358 | 332 |
| TBC1D3 | L0 = P or D | binary, shipped cuts | 12/19 | 1.000000 | 1.000000 | empty | no | >= 5,450 | 5,438 |

**One sentence per family.**

- **NPIP, protein level: CERTIFIED.** The 21 NPIP proteins are an exact connected component of the protein graph for
  every coverage cut in (0.128250, 0.840580]; the shipped 0.30 lies inside, and the component holds those 21 proteins
  and nothing else. The nearest outside protein is ACSM2A/ACSM2B at cover 0.128250 (through NPIPA2); NPIP's own weakest
  internal link is 0.840580.
- **TBC1D3, protein level: NOT certified at the shipped cut.** Its 9 proteins are exact only for cover in
  (0.688525, 0.900000]; TBC1D26 joins at 0.688525. At the shipped 0.30 the component has 5,275 proteins. This is the
  control: the protein test does not call every family clean.
- **No DNA level certifies either family.** On every DNA axis the strongest outside link reaches weight 1.000000 (a
  neighbouring record as close as the members are to each other), which is above every internal bottleneck.
  TBC1D3 at L3 is the one place where no outsider is inside the component — but its 12 nodes sit in 3 components at
  w_98 = 0.98 ({TBC1D3, B, D, E, F, G, H, I, K, P2}, {TBC1D3P5}, {LOC124905656}), so h_split (0.867946) is below the
  cut and the interval is still empty.

**Who joins, and at what weight** (top boundary edges, primary variant):

| set / level | joining nodes |
|---|---|
| NPIP P | ACSM2A 0.128250 (via NPIPA2), ACSM2B 0.128250, ACSM1 0.102253, ACSM5 0.068729 |
| TBC1D3 P | TBC1D26 0.688525 (via TBC1D3), TBC1D10C 0.510949, GRTP1 0.433515 |
| NPIP L1 identity | CLN3 1.000000 (via NPIPB10P), EIF3CL 1.000000 (via NPIPB9), PKD1P2 0.999370, PKD1P5-LOC105376752 0.999363, LOC131696449 0.999356 |
| TBC1D3 L1 identity | NPEPPSP1 1.000000 (via TBC1D3G), TBC1D3P1-DHX40P1 0.988213, TBC1D29P 0.961397, USP6 0.885467 |
| NPIP L2 | CLN3, LOC100190986, LOC124907830, LOC124907845, LOC128966632, all 1.000000 |
| TBC1D3 L2 | USP6 1.000000 (via TBC1D3P2), TBC1D29P 0.811940 |
| NPIP L3 | CLN3 1.000000, LOC100190986 1.000000, PKD1P5-LOC105376752 0.999363, LOC131696449 0.999356 |
| TBC1D3 L3 | TBC1D29P 0.961397, USP6 0.927995 |

CLN3 here is a 253-bp record of biotype "other", and LOC100190986 a 2.5-kb lncRNA: small records embedded in an NPIP
copy, aligning at identity 1.000. No cut on any declared axis can remove them; only a node rule could.

## 2. E-value axis (protein, cover cut 0.30)

| set | e <= 1e-5 | 1e-10 | 1e-20 | 1e-50 | 1e-100 |
|---|---|---|---|---|---|
| NPIP: component size at cover 0.30 | 21 | 21 | 21 | 21 | 21 |
| NPIP: exact interval | (0.128, 0.841] | (0.071, 0.841] | (0.000, 0.773] | (0.000, 0.372] | (0.000, 0.372] |
| TBC1D3: component size at cover 0.30 | 5,275 | 3,695 | 107 | 18 | 13 |
| TBC1D3: exact interval | (0.689, 0.900] | (0.689, 0.900] | (0.689, 0.900] | (0.689, 0.813] | (0.689, 0.813] |

NPIP stays exactly the 21 proteins at every E-value cut. TBC1D3's component shrinks as the cut tightens but never
becomes its 9 proteins at cover 0.30.

## 3. Subgroups (are the literature subfamilies exact anywhere?)

| set | protein | L1 identity | L2 f_ex | L3 w_98 |
|---|---|---|---|---|
| NPIPA (8) | empty | empty (h_join 0.999370) | empty | empty (h_join 0.999363) |
| NPIPB (19) | empty | empty | empty | empty |
| A6-9 (4) | empty | empty | empty | empty |
| B3-5 (4) | empty | empty | empty | empty |
| B6-9 (4) | empty | empty | empty | empty |
| B12/13 (3) | empty | empty | empty | empty by 0.000066 (joined 0.999349, splits 0.999283) |
| **B15 (3)** | empty | **(0.986296, 0.999517]** | empty | **(0.992435, 0.999517]** |
| NPIP u TBC1D3 | never connected | never connected | never connected | never connected |

Only B15 is an exact component of its own on a DNA axis, and its interval does not contain the shipped cut. NPIPA
fails only because readthrough and PKD1P records join it at 0.999363. The 7 NPIPA proteins are exact on the E-value
axis for cover in (0.845261, 0.859296] at e <= 1e-10, never at cover 0.30.

Internal order: NPIP first splits into Dishuck's NPIPA (with PKD1P6-NPIPP1) and NPIPB, at L1 identity 0.971665 and at
L3 w_98 0.980559.

## 4. Sensitivities (verdict changes only where stated)

| variant | what changes |
|---|---|
| tx records only (drop greedy gene-body chains) | B12/13 at L3 becomes exact in (0.997129, 0.999283]. Components shrink (NPIP L1 122 -> 74, L3 83 -> 59); TBC1D3 L2 h_split 0.832 -> 0.553. No other verdict moves |
| exon-less members as span nodes (adds 8 records: NPIPB14P, TBC1D3P1/P3/P4/P6/P7, LOC100420289, LOC100420311) | components grow (NPIP L1 122 -> 143, L3 83 -> 114; TBC1D3 L1 31 -> 67); no verdict changes |
| shipped chain denominator (min(query body, extrapolated span)) | no verdict changes (NPIP L1 component 122 -> 117) |

## 5. Method, in one block

- **Nodes.** Protein: the 20,088-protein database of §6ko (longest annotated CDS, r2 biotype rule, >= 10 aa). DNA:
  58,563 RefSeq gene and pseudogene records of `chm13v2.0_RefSeq_full.gff.gz` (RS_2025_08); 47,965 have >= 1 exon
  feature and are primary nodes. Members present: NPIP 21 proteins / 26 DNA nodes of 27; TBC1D3 9 / 12 of 19.
- **Protein evidence (H2/H3 clean).** blastp 2.17.0+, all 20,088 proteins as queries against the same database,
  `-evalue 1e-5 -max_target_seqs 100000 -dbsize 11710993 -num_threads 5`, max_hsps unlimited; 33 min, 752 MB peak.
  3,769,713 HSPs, 575,198 pairs. The fixed `-dbsize` makes E-values independent of the database (H2); the weight is the
  **union** of all HSP intervals on the longer protein, max over the two query directions (monotone, H3) — not the
  shipped greedy cover, so this graph is denser than §6ko's (397,031 pairs at cover >= 0.30).
  No query came near `max_target_seqs` (max 733 subjects).
- **DNA evidence.** Each node's representative transcript (`-x splice -uf -c -N 50 -p 0.1`) and gene body
  (`-x asm20 -c -N 50 -p 0.1`) aligned to the whole CHM13 genome with the prebuilt indexes; records on the query's own
  locus skipped; exon-record and greedy gene-body-chain witnesses, sense condition when both copies are spliced, exactly
  as `shared_definition.rs` (cross-checked: the shipped edge loop reproduces the witness pair set, symmetric difference 0).
  Query set: the members, then two hops of every node touched by any record (265 query nodes; a third hop would add 364).
- **Axes.** L1 identity with coverage 0.50 and coverage with identity 0.80; L2 f_ex over records; L3 single-record
  gap-excluded w_98 with shared exons >= 0.30; L0 as a binary join.

## 6. Caveats

- **Descriptive, not validation.** NPIP and TBC1D3 are development families; one haplotype (CHM13), one annotation.
  Adding proteins or genomes can only raise h_join (monotone tests, H1-H2), so a certificate can be lost, never gained,
  by more evidence — and that monotonicity argument covers the protein rows, the exon-record disjunct and the f_ex /
  w_98 maxima only. Every L1/L2/L3/L0 row depends on the greedy gene-body chains, which are not monotone (H3); §4 gives
  the tx-only sensitivity.
- **DNA component sizes are lower bounds.** Evidence from a node that no member query ever touched is assumed absent
  (declared). The tables record how many component nodes have no mapped query: 16 of 122 for NPIP L1; 5,395 of 5,450
  for TBC1D3 L0.
- **The protein weight is not the shipped test.** The shipped §6ko edge uses a greedy non-overlapping HSP cover; this
  report uses the monotone union cover, which admits a superset of the shipped edges.
- **CLN3 / lncRNA / readthrough boundary edges at weight 1.000** are what blocks every DNA certificate. Whether such
  records should be nodes is the open node-rule question of §0★★★.1, not a threshold question.

## Verification (independent recompute)

Agent 4 (independent verifier), 2026-09-17, 10:23-11:03 -07:00. No script under `family_cert/` outside
`family_cert/verify/` was read or imported; the verifier's code is in `family_cert/verify/` (`protein/cover_check.py`,
`dna/vnodes.py`, `dna/vqueries.py`, `dna/vwitness.py`, `dna/vpairs.py`, `dna/vcompare.py`, `cert/vcert.py`,
`cert/vcmp.py`, `cert/vsweep.py`). The report `bench/FAMILY_CERTIFICATES_NPIP_TBC1D3.md` did not exist at 11:03 -07:00
(40 minutes after the last write under `family_cert/cert/`), so no report sentence or number could be checked or
corrected; everything below was recomputed from the declarations and the evidence tables.

1. **Order.** `DECLARATIONS.txt` was written 09:17:35, before the first evidence file (`protein/queries/b000.faa`
   09:17:50; `protein/batches/b000.tsv` later). Addendum A2 (10:01:05) precedes the reuse index (10:03:46) and the first
   minimap2 run (10:03:50); A3 (10:16:04) precedes the first cert file (10:16:37). A1, A2R and A3 are dated and
   disclosed as post-evidence. `dna/nodes.tsv` (09:58:50) predates A2, as A2 itself states.

2. **Protein.** 109 queries re-run with the recorded flags (all 30 NPIP/TBC1D3 proteins, all 29 further proteins with
   any HSP touching them, 50 random others; blastp 2.17.0+, `-evalue 1e-5 -max_target_seqs 100000 -dbsize 11710993
   -num_threads 5`): 8,285 non-self HSPs at e <= 1e-5, **byte-identical to the 8,285 rows of `protein/hsps.tsv`** for
   those queries. `blastdbcmd -info` and an independent FASTA recount both give 20,088 sequences / 11,710,993 residues =
   the `-dbsize` constant. Max distinct subjects for any query in `hsps.tsv` is 733 (+ self), far below
   `max_target_seqs`. Union covers recomputed from the re-run HSPs for **all 495 pairs touching the two families**:
   identical key set and 0 mismatches in `w_P`, the four eps covers and `best_hsp_pident`. Global recount of
   `hsps.tsv`/`pairs.tsv`: 3,769,713 HSPs, 1,131,015 directions, 575,198 pairs, 19,381 one-direction-only, 397,031
   pairs at cover >= 0.30, 7,188 components, largest 5,275 — all as reported.

3. **DNA.** Nodes rebuilt independently from the same GFF: 58,563 records, 47,965 with >= 1 exon, and every node's
   chromosome, span, strand, body, exon union, representative transcript and representative exons **identical** to
   `dna/nodes.tsv` (0 differing fields). All 530 query sequences regenerated from `idx/target.fa`: 0 md5 or length
   mismatches. 24 queries re-mapped (12 tx incl. 3 NPIP + 2 TBC1D3, 12 body incl. 3 NPIP + 2 TBC1D3, from both reused
   union PAFs and this run's batches): 96/96 tx and 145/145 body records identical, in the same order. Witnesses
   re-derived with the verifier's own chain, exon-block, sense and sx code: 5,027 member-endpoint rows, kind x variant
   counts identical, 960 primary pairs / 137 exonless_span pairs, and every pair weight (`id_w`, `cov_w`, `t1`, the
   tx-only variants, `t1_shipden`, `f_ex`, `w_98`, row counts) agrees with `cert/dna_pairs.members.tsv` (1,097 pairs)
   and `cert/dna_pairs.extended.tsv` (2,303 pairs) to 1.5e-6; the only differences are 7 + 44 pairs where the tables
   write `f_ex = NA` and the verifier writes 0.0 (pairs with chain-only evidence; no certificate changes). Two-hop
   expansion reproduced exactly: touched(hop0) 119 -> hop1 81, touched(hop1) 250 -> hop2 138, touched(hop2) 516 ->
   would-be hop3 364; hop record counts 1,529/105, 1,572/142, 1,051/255.

4. **Certificates.** All 210 rows of `cert/certificates.tsv` recomputed from the verifier's own evidence with an
   independent implementation of D1 (h_join = max boundary weight; h_split = maximum-spanning-forest bottleneck;
   interval, exactness, shipped-cut membership, component at the shipped cut, component nodes outside S, number of
   components touching S, component nodes without a mapped query): **0 differences** in `n_present`, `h_join`,
   `h_split`, `exact_nonempty`, `shipped_inside`, `comp_size_at_shipped`, `comp_outside_n`,
   `comp_nodes_without_query` and `comp_n_components_touching_S` across every set (NPIP, TBC1D3, NPIPA, NPIPB, A6-9,
   B3-5, B6-9, B12/13, B15, NPIP+TBC1D3), level (P with its four eps axes and the report-only identity axis, L1
   identity and coverage, the shipped-denominator sensitivity, L2, L3, L0 and the L0 sweep) and variant (primary,
   exonless_span, txonly_sensitivity). In all 196 rows with a joining edge, the top joining edge's weight equals
   h_join.

5. **Scope.** The verification covers the arithmetic and the evidence, not the declarations' reach: X -> S evidence
   from nodes never touched by an S query is assumed absent (D3), so every DNA and L0 component size remains a lower
   bound, and D5 stands unchanged.
