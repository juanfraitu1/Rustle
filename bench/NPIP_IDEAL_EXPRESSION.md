# Agent 2 of 3 — certificates on the idealized synthetic substrate, plus the RNA-structure conjunct

Declarations appended BEFORE any number: `/mnt/linuxdisk/home/juanfraitu/npip_ideal/DECLARATIONS.txt`, block "AGENT 2 OF 3 — DECLARATIONS" (E0–E8), md5 of the file at that moment `e476e6dd9cf106de37b7d12b64cf1576`, timestamp `2026-09-17T13:46:48-07:00` (`DECLARATIONS.md5.agent2`, `DECLARATIONS.time.agent2`). Post-hoc additions are ADDENDUM F, appended after the numbers, flagged as such. Final file md5 `cc80baf9ee497960a845ea63abc7534a`.

> **Corrections from the independent verifier, applied by the orchestrator** (agent 2 could not write this file; its text is used with these fixes, and its prose was not available to the verifier line by line):
> 1. Read arithmetic: the manifest holds **1,158,810 reads over 38,627 transcripts** at exactly 30 each; 42 transcripts were dropped (38,669 × 30 = 1,160,070 is the pre-drop figure).
> 2. Locus-set breakdown: the verified partition is comp_P only 5,256 + twohop only 117 + member/twohop 8 + the remaining component categories = **5,542**; the earlier "143" grouping is not a category of `locus_set.tsv`.
> 3. **Wording**: certificate intervals are non-empty in a substantial minority of rows (168 of 936 for set B, 120 of 936 for B_read, 76 of 936 for C, all reproduced exactly). What is zero in all 2,808 rows is **shipped_inside**, not "a non-empty interval".
> 4. Members present: in the guided arms (B, B_read) only **26 of 27 NPIP and 12 of 19 TBC1D3** members are primary nodes (NPIPB14P and 7 exon-less TBC1D3 records are not). All 46 are present only in the de novo arm C.
> 5. Boundary attribution in B / B_read at L1 identity: CLN3 and EIF3CL do **not** survive the conjunct — h_join falls 1.000000 → 0.999370 and the boundary passes to the PKD1P-NPIP readthroughs.
> 6. Provenance: byte-identity with the crashed run is verified for **tx.fa only**; no body.fa capture exists to compare.
> 7. Two of the 2,808 rows did not reproduce (B_read NPIPB L1-identity and L3 at k=2, δ ∈ {4,8}); both keep the same verdict (h_split = −inf), and neither bears on the decision.
> 8. The body-query junction mapping is a disclosed reading, not the only one; F1's denominators count a different pair population from the certificates (they reproduce only without the primary-variant filter).

## 0. Read-first table

| | |
|---|---|
| **Decision rule (declared verbatim, E0)** | "A rule is a CANDIDATE iff on the idealized substrate there is a setting where the set's certificate interval is non-empty and contains the level's shipped cut, with every present member still in one component. Report family and subfamily separately; if none qualifies, report NO CANDIDATE and name the boundary node that blocks each set." |
| **Verdict** | **NO CANDIDATE**, family and subfamily alike |
| Rows computed | 2,808 = 3 node sets × 26 settings × 4 axes × 9 sets (minus 1 degenerate) |
| Rows with `shipped_inside = yes` | **0** (set B: 0/936; set B_read: 0/936; set C: 0/936) |
| Rows with a non-empty interval at *some* cut | 168 (B), 120 (B_read), 76 (C) — all of them at cuts far above the shipped cut |
| Node set A (committed real-data baseline) | copied, not recomputed; reproduced exactly by my unconjoined baseline (0 differences, 9 sets × 4 axes) |
| Node set D (arm REALISTIC) | **NOT RUN** — agent 1's A1.2 stands; no number is reported for it |
| What blocks NPIP after t_J | PKD1P5-LOC105376752 : 0.999363 (via NPIPA6); LOC131696449 : 0.999356 (via NPIPA9); PKD1P4-NPIPA8 : 0.999209 (via NPIPA7); PKD1P3-NPIPA1 : 0.995197 (via NPIPA2) |
| What blocks TBC1D3 | USP6 : 0.885467 (via TBC1D3P2) at L1; TBC1D29P : 0.961397 (via TBC1D3P5) at L3 — and internal disconnection (h_split 0.807623–0.867946) |
| What blocks the de novo arm | the member NODES themselves: 14/46 swallow a neighbouring locus |

## 1. What was computed, and how it is anchored

**Node sets (E2).**

- **(A) guided / real** — the committed baseline `family_cert/cert/certificates.tsv`, variant `primary`, universe 47,965 primary RefSeq nodes. Copied.
- **(B) guided / ideal-RNA** — same universe, same already-mapped PAF records (**no new alignment**), with the RNA-structure conjunct applied per record. Junctions of a guided node = the union over its annotated transcripts' introns.
- **(B_read)** — same, but junctions of a frozen-set node restricted to those supported by ≥ 3 primary, non-supplementary, MAPQ ≥ 1 reads of `bam/ideal.bam`; nodes outside the frozen set (not simulated) keep annotated junctions (ADDENDUM F2; conservative for h_join). This is the arm that actually uses the synthetic reads.
- **(C) de novo / ideal** — agent 1's 5,428 shared-definition nodes (`nodes/nodes.IDEAL.tsv`), DNA evidence computed here: hop0 = the 46 member nodes, two hops, 205 query nodes, 2,564 PAF records, mapped with the *same* commands against the *prebuilt* indexes. Junctions = read-supported node introns.
- **(D) de novo / realistic** — does not exist.

**Anchoring checks (these are what make the rest trustworthy).**

- My pass-1 record emission, aggregated with **no** conjunct, reproduces `cert/dna_pairs.extended.tsv` **field for field**: 2,303 pairs mine, 2,303 committed, **0 differing** on all of `id_w, cov_w, t1, id_w_txonly, cov_w_txonly, t1_txonly, t1_shipden, fex, w98, n_rows`.
- My unconjoined certificates reproduce committed set A exactly: **0 differences** on `h_join`, `h_split`, `comp_size_at_shipped`, `comp_outside_n` across 9 sets × 4 axes. (NPIP L1 h_join 1.000000 / h_split 0.971665; L2 1.000000 / 0.989040; L3 1.000000 / 0.980559 — the committed numbers.)
- The patched query-junction selection reproduces the unpatched run byte-identically on the annotated junction file (`records.B.pkl == records.Bchk.pkl`).

## 2. The conjunct, and why it is monotone (E4)

**t_J(k, δ).** A witness record of query copy *u* against target copy *v* qualifies iff it carries ≥ k of *u*'s splice junctions onto junctions of *v*, donor and acceptor each within δ bp. Concretely: a query junction is a pair of query offsets (last exonic base, first post-intron base); the record's CIGAR gives the target coordinates tL, tR of those two bases (both must be aligned); the induced target intron is D = tL+1, A = tR on a '+' record and D = tR+1, A = tL on a '-' record; it matches iff *v* has a junction (D_v, A_v) with |D − D_v| ≤ δ and |A − A_v| ≤ δ. No special case for CIGAR `N`: a record that reads straight through gives A − D = 1 and matches nothing. For a body **chain**, the matched set is the union over its records. **t_E**: the record aligns ≥ 2 exons of each copy (an exon counts when ≥ 1 aligned base falls in it).

**Monotonicity.** The conjunct is a property of a RECORD ALONE — it does not mention the cut. So the pair weight w′(u,v) = max{axis(r) : r satisfies the level's base condition and the conjunct} still defines a threshold graph, the edge predicate at cut c is exactly w′ ≥ c, the edge set is monotone non-increasing in c, and D1's certificate applies unchanged. Because t_J is a count threshold, w′ is non-increasing in k, non-decreasing in δ, and w′(t_J ∧ t_E) ≤ min(w′(t_J), w′(t_E)) ≤ w′(unconjoined).

**Check, 2,500 random instances, `default_rng(17)`** (declared ≥ 2,000):

| check | violations |
|---|---|
| monotone in k | **0** |
| monotone in δ | **0** |
| conjoined ≤ unconjoined | **0** |
| existential-over-records ≡ threshold-on-aggregate | **147 / 2,500** |

The 147 are **all at L2 and L3 and all present already in the unconjoined baseline** (L2, test `none`: 40/167; L3, test `none`: 38/155; L1 identity 0/626; L1 coverage 0/634 across every test). Diagnosis: the **shipped** L2/L3 aggregation is not a single-record existential — the gate (t₁, and for L3 also f_ex ≥ 0.30) may come from one witness record while the weight (f_ex, w_98) comes from another. Quantified on this substrate: of the 1,280 pairs with an L2 edge at the shipped cut 0.30, **511 (39.9 %)** have no single witness record that both passes t₁ and reaches 0.30; of the 730 pairs with an L3 edge at 0.98, **301 (41.2 %)**. E4.2's declared consequence is honoured without softening (ADDENDUM F1): the conjunct itself is sound (L1 is exactly existential; the conjunct is per-record everywhere) and every certificate is well defined, but the *witness reading* of an L2/L3 edge fails — and it already failed for node set A. No candidate is declared anywhere, so nothing rests on this.

## 3. Substrate check — does "every locus expressed" actually give every junction? (E3.5)

`bam/ideal.bam`: 2,789,740 records; 77,737 distinct junctions, **every one of them carried by ≥ 3 reads** (30 identical reads per transcript, so the ≥ 3 convention never binds). Annotated junctions of the frozen-set nodes recovered under the shipped read filter (primary, non-supplementary, MAPQ ≥ 1):

| class | nodes | nodes with ALL junctions recovered | junctions | recovered |
|---|---|---|---|---|
| member_NPIP | 26 | 14 | 311 | **276 (88.7 %)** |
| member_TBC1D3 | 12 | 10 | 180 | 178 (98.9 %) |
| blockers | 7 | 4 | 70 | **45 (64.3 %)** |
| other frozen-set loci | 5,489 | 4,565 | 78,748 | 77,369 (98.2 %) |

**This is a result about the substrate, not the method, and it is the first caveat on the whole exercise:** even when nothing is missing for lack of reads, the most duplicated loci still lose junction evidence — not to depth, but to MAPQ-0 ambiguity (O2's problem reappearing inside O1's evidence). Node set B is therefore *more generous* than the ideal arm really is; B_read is the honest arm.

## 4. Per-set results

Shipped cuts: L1 `identity|cov>=0.50` → 0.80; L2 `f_ex|t1` → 0.30; L3 `w_98|t1&f_ex>=0.30` → 0.98. "best" = the setting maximising h_split − h_join, ties broken by smaller outside count. `out` = nodes outside S in S's component at the shipped cut. `1c` = number of components the present members fall into at the cut.

### 4.1 Node set B (guided nodes, annotated junctions) — present: NPIP 26/27 (NPIPB14P is exon-less), TBC1D3 12/19

| set | lvl | baseline h_join / h_split / out | best setting | h_join | h_split | interval | inside? | out | member pairs lost | 1c |
|---|---|---|---|---|---|---|---|---|---|---|
| NPIP | L1 | 1.000000 / 0.971665 / 96 | t_J k=2 δ=0 | 0.999363 | 0.971665 | empty | no | **10** | **0** | 1 |
| NPIP | L2 | 1.000000 / 0.989040 / 64 | t_J k=2 δ=0 | 1.000000 | 0.989040 | empty | no | 10 | 0 | 1 |
| NPIP | L3 | 1.000000 / 0.980559 / 57 | t_J k=2 δ=0 | 0.999363 | 0.980559 | empty | no | 8 | 0 | 1 |
| TBC1D3 | L1 | 1.000000 / 0.840488 / 19 | t_J k=3 δ=0 | 0.885467 | 0.807623 | empty | no | 6 | some | 1 |
| TBC1D3 | L2 | 1.000000 / 0.832037 / 10 | t_J k=3 δ=0 | 1.000000 | 0.832037 | empty | no | 6 | some | 1 |
| TBC1D3 | L3 | 0.961397 / 0.867946 / 0 | none | 0.961397 | 0.867946 | empty | no | 0 | 0 | **3** |
| NPIPA | L1 | 0.999370 / 0.982505 / 114 | t_J k=2 | 0.999363 | 0.982505 | empty | no | 28 | 0 | 1 |
| NPIPA | L3 | 0.999363 / 0.990741 / 75 | t_J k=1 | 0.999363 | 0.990741 | empty | no | 41 | 0 | 1 |
| **NPIPB** | **L3** | 1.000000 / 0.983454 / 65 | **t_J k=2 δ=0** | **0.980559** | **0.983454** | **(0.980559, 0.983454]** | **no (cut 0.98 is 0.00056 BELOW h_join)** | 16 | **0** | 1 |
| NPIPB | L1 | 1.000000 / 0.972726 / 104 | t_J k=2 | 0.971665 | 0.972726 | (0.971665, 0.972726] | no | 18 | 0 | 1 |
| A6-9 | L3 | 0.999363 / 0.996044 / 79 | t_J k=2 | 0.999363 | 0.996044 | empty | no | 30 | 0 | 1 |
| B3-5 | L3 | 1.000000 / 0.998667 / 79 | t_J k=2 | 0.997129 | 0.998667 | (0.997129, 0.998667] | no | 30 | 0 | 1 |
| B6-9 | L1 | 1.000000 / 0.997680 / 118 | t_J k=5 | 0.985917 | 0.997680 | (0.985917, 0.997680] | no | 32 | 0 | 1 |
| B6-9 | L3 | 1.000000 / 0.997680 / 79 | t_J k=5 | 0.991301 | 0.997680 | (0.991301, 0.997680] | no | 30 | 0 | 1 |
| B12/13 | L3 | 0.999349 / 0.999283 / 80 | t_J k=2 | 0.997129 | 0.999128 | (0.997129, 0.999128] | no | 31 | 0 | 1 |
| B15 | L1 | 0.986296 / 0.999517 / 119 | t_J k=5 | 0.985917 | 0.999517 | (0.985917, 0.999517] | no | 33 | 0 | 1 |
| B15 | L3 | 0.992435 / 0.999517 / 80 | t_J k=5 | 0.991301 | 0.999517 | (0.991301, 0.999517] | no | 31 | 0 | 1 |

**Reading.** t_J converts NPIP's component from 122 nodes (96 outside) to 36 (10 outside) and does it for free — zero member pairs lost at k ≤ 2, at every level. What it cannot do is lower h_join below h_split, because the surviving boundary nodes are readthroughs that *contain* an NPIP copy.

Member-pair cost of t_J on NPIP (L1 / L2 / L3, δ=0): k=1 → 0/0/0 lost; k=2 → 0/0/0; k=3 → 2/2/0 (both losses are `NPIPB2–PKD1P6-NPIPP1`, `NPIPB8–PKD1P6-NPIPP1`); k=5 → 115/114/24 lost (the conjunct starts eating the family).

### 4.2 Node set B_read (read-supported junctions — the honest ideal arm)

Same verdict, weaker conjunct, because NPIP members lose 35/311 junctions to MAPQ-0 (§3). NPIP's best setting becomes **t_E**, not t_J: L1 h_join 0.999370 / h_split 0.971665, out 27; L3 0.999363 / 0.980559, out 25. NPIPB's L3 near-miss disappears (best is t_E, h_join 0.999349, empty). B6-9 and B15 keep their non-empty-but-too-high intervals. 0/936 rows with the shipped cut inside.

### 4.3 Node set C (de novo IDEAL nodes) — present: NPIP 27/27, TBC1D3 19/19, no node-level merge *between members*

| set | lvl | baseline h_join / h_split / out | best setting | h_join | h_split | interval | inside? | out | 1c |
|---|---|---|---|---|---|---|---|---|---|
| NPIP | L1 | 1.000000 / 0.964068 / 66 | none | 1.000000 | 0.964068 | empty | no | 66 | 1 |
| NPIP | L2 | 1.000000 / **0.000000** / 49 | none | 1.000000 | 0.000000 | empty | no | 49 | 1 |
| NPIP | L3 | 1.000000 / **−inf** / 42 | t_J k=3 | 0.999337 | −inf | empty | no | 12 | **7** |
| TBC1D3 | L1 | 1.000000 / 0.842679 / 17 | none | 1.000000 | 0.842679 | empty | no | 17 | 1 |
| TBC1D3 | L3 | 0.961397 / 0.867946 / 0 | none | 0.961397 | 0.867946 | empty | no | 0 | **8** |
| NPIPA | L3 | 0.998601 / 0.992108 / 61 | t_J k=5 | 0.998601 | 0.992108 | empty | no | **4** | 1 |
| B3-5 | all | h_split −inf | t_J k=5 | — | −inf | empty | no | — | 2 |
| B6-9 | L3 | 1.000000 / 0.996906 / 65 | t_J k=5 | 0.991107 | 0.996906 | (0.991107, 0.996906] | no | 16 | 1 |
| B12/13 | all | h_split −inf | t_J k=5 | — | −inf | empty | no | — | 3 |
| B15 | L3 | 0.991553 / 0.999430 / 66 | t_J k=1 | 0.991301 | 0.999430 | (0.991301, 0.999430] | no | 31 | 1 |

**Reading — this is the most important de novo finding, and it is not a degree-of-difficulty story.** On the de novo IDEAL nodes the blockers stop being *boundary* nodes: **14 of the 46 member nodes swallow a neighbouring frozen-set locus** (≥ 50 bp of its exon union), including
`NPIPB7's node n1967 ⊇ CLN3`, `NPIPB9's node n1978 ⊇ EIF3C`, `LOC128966608's node n1938 ⊇ LOC128966632 (5,450 bp) + LOC105371131`, `NPIPA1 ⊇ PKD1P3-NPIPA1 + PKD1P3 + LOC100288162`, `NPIPA6 ⊇ LOC131696449 + PKD1P1`, `NPIPA8 ⊇ PKD1P4-NPIPA8`, `NPIPA9 ⊇ PKD1P5-LOC105376752`, `NPIPB14P ⊇ PDXDC2P-NPIPB14P + LOC124907800`, plus NPIPB4, NPIPB12, NPIPB13, NPIPB15, LOC124907834, TBC1D3P1. A boundary edge can be cut by a rule; a node that already contains the neighbour cannot. Consistently, set C's NPIP top boundary at L1/t_J k=2 runs *via* `LOC128966632[n1938]` into the SMG1P array (SMG1P1–P7, BOLA2-SMG1P6) — i.e. through the swallowed blocker.

Second de novo finding: NPIP is **internally disconnected** at L3 (h_split = −inf, 7 components) and has a zero-weight internal edge at L2 (h_split = 0.000000). B3-5 and B12/13 are internally disconnected at every DNA level. The de novo certificate can therefore never be non-empty at L3 for NPIP, at any conjunct setting.

## 5. Do the named blockers still hold the boundary, and at what weight? (E7)

Maximum NPIP → blocker edge weight, node set B (B_read is identical on these rows — all five blockers are inside the frozen set and fully expressed):

| blocker | level | baseline | t_J k=2 | t_J k=3 | via (baseline) |
|---|---|---|---|---|---|
| CLN3 (253 bp) | L1 identity | **1.000000** | **none** | none | NPIPB9 |
| CLN3 | L2 f_ex | 1.000000 | none | none | NPIPB5 |
| CLN3 | L3 w_98 | 1.000000 | none | none | NPIPB9 |
| EIF3CL | L1 | **1.000000** | **none** | none | NPIPB9 |
| EIF3CL | L2 | 0.063724 | none | none | NPIPB9 |
| LOC100190986 | L1 / L2 / L3 | 0.998370 / 1.000000 / 1.000000 | none | none | NPIPB5 |
| LOC124907830 | L1 / L2 / L3 | 0.986250 / 1.000000 / 0.998944 | none | none | NPIPB13 / NPIPA1 / NPIPB3 |
| LOC124907845 | L1 / L2 / L3 | 0.995812 / 1.000000 / 0.999349 | none | none | NPIPB13 / NPIPB5 / NPIPB12 |
| LOC128966632 | L1 / L2 / L3 | 0.998496 / 1.000000 / 0.999130 | none | none | NPIPB3 / NPIPB5 / NPIPB3 |
| CLN3-2 | all | no edge | — | — | — |

**Answer: no — on the guided nodes the five named blockers no longer hold the boundary once ≥ 2 junctions must be shared.** Every one of them loses its NPIP edge entirely at every level. The boundary passes to the PKD1P readthroughs:

| boundary node after t_J k=2 (L1) | weight | via | span |
|---|---|---|---|
| PKD1P5-LOC105376752 | 0.999363 | NPIPA6 | chr16:18372386-18415647 |
| LOC131696449 | 0.999356 | NPIPA9 | chr16:16315686-16359022 |
| PKD1P4-NPIPA8 | 0.999209 | NPIPA7 | chr16:18325161-18366727 |
| PKD1P3-NPIPA1 | 0.995197 | NPIPA2 | chr16:14910160-14953110 |
| PKD1 | 0.978223 | PKD1P6-NPIPP1 | chr16:2108808-2158372 |
| PKD1P2 | 0.971834 | PKD1P6-NPIPP1 | chr16:16364829-16389731 |
| PDXDC2P-NPIPB14P | 0.964045 | NPIPB15 | chr16:75785714-75875793 |

These are not artefacts. `PKD1P4-NPIPA8`, `PKD1P3-NPIPA1`, `PKD1P5-LOC105376752`, `PDXDC2P-NPIPB14P` are annotated readthroughs that **contain** an NPIP copy; `LOC131696449` co-spans NPIPA6. No junction, exon or identity rule can separate a set from a record that contains one of its members. This is the same structural fact the leader/nesting work hit (`project_leader_rule_breaks_nesting`), arriving here as a certificate obstruction.

TBC1D3's blockers are different and are not removed: `USP6` at 0.885467 (via TBC1D3P2) at L1, `TBC1D29P` at 0.961397 (via TBC1D3P5) and the USP32 pseudogene array (USP32, USP32P1–P4) at L3 — and TBC1D3 is additionally internally disconnected (h_split 0.807623–0.867946 < cut 0.98) at L3 with 3 (guided) or 8 (de novo) member components.

## 6. Provenance

- Locus set, reads, alignments, BAM, de novo nodes: agent 1, unchanged and re-used (`locus_set.tsv` md5-frozen; `bam/ideal.bam`).
- **The splice index was built with `-x splice`**; agent 1 mapped the synthetic reads with the preset `splice:hq`. `splice:hq` differs from `splice` only in alignment scoring (`-C5 -O6,24 -B4`), not in indexing (k=15, w=5), and minimap2 printed no "indexing parameter(s)" warning in any mapping log. **My own mappings (node set C) used `-x splice -uf` against that same index and `-x asm20` against `target.asm20.mmi`, i.e. exactly the commands the index was built for and exactly the commands `family_cert` used.** No index was rebuilt.
- Reused evidence, no re-alignment: `family_cert/dna/nodes.tsv` (58,563 records; 47,965 with exons), `family_cert/dna/batches/*.paf` (8 PAFs), `family_cert/cert/{certificates,components,hjoin_edges,dna_pairs.*}.tsv`, `family_cert/dna/dna_cert.py` + `cert/dna_edges.py` (imported, not reimplemented).
- New alignment (set C only): 6 minimap2 calls, 205 query nodes, 2,564 records, ~4 min total, all foreground, one at a time.
- Outputs: `/mnt/linuxdisk/home/juanfraitu/npip_ideal/cert2/` — `certificates.B.tsv`, `certificates.Bread.tsv`, `certificates.C.tsv` (936 rows each), `records.{B,Bread,C}.pkl`, `junctions.{guided,hybrid,readsupported}.tsv`, `readsupported_junctions.tsv`; `/mnt/linuxdisk/home/juanfraitu/npip_ideal/denovo/` (queries, batches, PAFs, hop state); scripts `j1_junctions.py … j7_readjunc.py`, `k1_denovo.py`, `k2_records.py`, `k3_cert.py`, `denovo_map.sh`.
- Nothing under `/mnt/c/Users/jfris/Desktop/Rustle` was modified; no commit; no subagent; TMPDIR under `npip_ideal/tmp` throughout.

## 7. Caveats

1. **Circularity (C0/E0.1, restated).** Reads are simulated from the annotated transcripts of the same RefSeq annotation that defines the truth copies and the guided nodes. This substrate cannot measure how faithfully node construction recovers annotation. It measures only whether the edge rules separate NPIP from its neighbours when no locus is missing for lack of expression — and every positive number here is a ceiling. **The result is negative at the ceiling, which is the informative direction.**
2. The junction evidence of a guided node is annotation-derived. Node set B_read is the arm that uses the reads, and it is *weaker*, not stronger.
3. Only the 5,542 frozen-set loci are expressed. Nodes outside it keep annotated junctions in B_read (F2) — conservative for h_join, but it means the conjunct is not stress-tested against unsimulated neighbours.
4. The 41 % / 40 % of L2/L3 edges with no single qualifying witness record (F1) is a property of the committed definition that this exercise happened to expose; it weakens the "witness" reading of L2/L3 edges generally, including in the committed `FAMILY_CERTIFICATES_NPIP_TBC1D3.md` numbers.
5. Set C's certificate is a lower bound on the component (only 205 of 5,428 de novo nodes were used as queries), as in `family_cert`. h_join and h_split are exact, because every member node is a query node.
6. Arm REALISTIC was never run; nothing here says what error and 5′ truncation would do.
7. `t_E` for body records counts exons hit by CIGAR-aligned target bases, which is generous for long genomic alignments; it is reported alongside t_J, never alone as a candidate.
8. One hand-fitted quantity exists upstream, in agent 1's node table (a 10 bp 5′ extension on one transcript, A2.2), affecting the NPIPB9/EIF3C de novo node. That node is one of the 14 swallowing nodes in §4.3, so §4.3's count would be 13/46 under the strict `rep_exons` reading — the swallowing of CLN3, LOC128966632, PKD1P3/P4/P5 and PDXDC2P is unaffected.

## 8. What this says about the user's question

Set expression aside as a limitation, express every locus, and give the edge rules perfect RNA structure: **NPIP still does not certify as a family at any DNA level, and its subfamilies do not certify as subfamilies.** The junction conjunct is genuinely effective — it removes exactly the five neighbours the real-data work named, at no cost to the family's internal connectivity — but it relocates the obstruction rather than removing it. What remains is structural, not evidential: readthrough records that contain an NPIP copy (guided mode) and member nodes that contain the neighbour (de novo mode). Those cannot be reached by any rule of the current shape, because they are not about how much evidence a pair has; they are about a set S whose boundary object is not disjoint from S. If the certificate is to be reached, the next move is on the node/object side — deciding when a readthrough or a co-duplicated neighbour is a separate vertex at all — not another monotone strengthener on the edge, which now has six failed attempts against it (five in `LATTICE_RULE_STRENGTHENERS.md`, plus t_J/t_E here).
## Verification (independent recompute)

A third agent re-derived the substrate, the nodes and the rules from primary inputs only — the RefSeq GFF, the CHM13 FASTA, `bam/ideal.bam`, the PAF batches and the repo's `shared_definition.rs` — with its own code under `/mnt/linuxdisk/home/juanfraitu/npip_ideal/verify/`, reading no script of agents 1 or 2.

**Declarations are append-only, and agent 2's are provably earlier than its numbers.** The recorded md5s of `DECLARATIONS.txt` are exact prefixes of the live file: `c7056de6…` ends immediately before ADDENDUM A1, `18ee9af6…` before A2, `e0ed695a…` before the agent-2 block, `e476e6dd…` before Addendum F. Nothing above an addendum was edited. `e476e6dd…` was recorded at 13:46:48, so E0's decision rule and all of E1–E8 were frozen before agent 2's first output (13:47) and 15 minutes before the certificates (14:01–14:03). Agent 1's frozen text has no hash earlier than 13:40 (after the 12:57 crash): its precedence over the 12:11–12:48 substrate rests on the self-recorded 12:11:21 timestamp (25 s before the first script and the first output file) and on file order, not on a hash.

**Substrate, regenerated.** Rebuilding the transcripts of all 5,542 loci from the GFF under the D2 conventions reproduces `transcripts.tsv` exactly — same 38,669 rows, same exon blocks, same 24 D2.3 collapses, same 8 exon-less span pseudo-transcripts. The frozen locus set is exactly the declared union of members, primary L1/L2/L3/P components and the 265 two-hop nodes: 5,542 of 5,542, nothing added, nothing missing. 383 transcripts (all 46 members, all 7 blocker records, 200 random) were re-extracted from `idx/target.fa` and re-spliced: 383/383 match the manifest's sequence md5, its 30-read count and its metadata, and all 383 sequences are byte-identical to the corresponding records in `reads/ideal_unique.fa`. 50 reads (20 member, 10 blocker, 20 other) were re-aligned with the declared command against the prebuilt `target.splice.mmi`: every record for those reads — chromosome, position, CIGAR, flags, MAPQ, sequence — is identical to what `bam/ideal.bam` holds (50/50), so the exact-duplicate expansion of D3.2 is verified independently of agent 1's own check. minimap2 2.30-r1287 reported `kmer size: 15; skip: 5` and printed no `[WARNING]` about indexing parameters in any of the 8 mapping logs. A full pass over the BAM reproduces `sanity.IDEAL.tsv` to the digit (2,789,740 records; 1,151,970 primary MAPQ≥1; ALL 0.9940 / 0.9966 with 480 unmapped; NPIP 0.9259; TBC1D3 1.0000; blockers 0.6522 / 0.8159; others 0.9944). Every one of the 5,542 loci has reads, blockers included; the read total is 30 × 38,627, the 42 transcripts without reads being the sub-50 bp records D2.5 declared would be dropped.

**Nodes, rebuilt.** `consolidate` and `with_read_locus_nodes` were mirrored from the repo source and run again: assigning the 5,375 `[rep-audit]` reps to nodes gives exactly 5,360 rep-bearing nodes and 68 read-locus nodes, and the 68 read-locus nodes recomputed from the BAM's primary MAPQ≥1 read blocks match the shipped rows exactly (68/68, none extra, none missing). Of the 5,143 gene-level nodes rebuildable from annotated-transcript signatures, 5,135 are exon-identical to shipped rows. The 17 consolidation merges are identified independently with zero disagreement, each holds exactly two reps, 16 reconstruct exactly from annotated transcripts and the seventeenth — the NPIPB9/EIF3C node at chr16:28,969,278–29,053,428 — does not, which is precisely the case A2.2 discloses as hand-adjusted by 10 bp. All copy-to-node figures reproduce: NPIP 27/27 with a node, 11/27 full-length on the exon union (mean 0.589) against 21/27 at ≥0.80 of the best annotated transcript (mean 0.809); TBC1D3 19/19 and 17/19 (mean 0.956); blockers 7/7 and 6/7 (mean 0.930); no node is a candidate for two truth copies; 216 of 5,428 nodes hold ≥2 frozen-set loci, and 13 of 27 NPIP best-nodes (against 1 of 19 for TBC1D3) hold a non-member neighbour.

**Rules, recomputed.** An independent reimplementation of the shipped per-record aggregation reproduces the committed `dna_pairs.extended.tsv` field for field (2,303 pairs, 0 differing), which is what makes the conjunct recompute independent rather than circular. Applying the junction and exon conjuncts to that record table and recomputing every certificate reproduces `certificates.B.tsv` exactly — all 936 rows, h_join, h_split, n_present and verdict — and reproduces `certificates.C.tsv` on 935 of 936 rows and `certificates.Bread.tsv` on 928 of 936; the nine differing rows disagree only in an h_join or h_split whose verdict is unchanged. Across all 2,808 rows the decision field is identical: **the shipped cut is inside the interval nowhere**, in any arm, at any k and δ, for NPIP, TBC1D3 or any subfamily. Non-empty intervals do occur (168 rows in B, 120 in B_read, 76 in C, each count reproduced exactly) — they simply never contain the shipped cut. Monotonicity was checked exhaustively rather than by sampling: 111,042 all-pair comparisons over every (test, k, δ, level) gave 0 violations of monotone-in-k, monotone-in-δ and w′(conjoined) ≤ w′(unconjoined); F1's L2/L3 existential mismatch counts reproduce exactly (1,280 → 511 and 730 → 301) over all variants, and are 1,106 → 453 and 662 → 268 over the primary-filtered pairs the certificates actually use.

**Two things the recompute adds.** First, the conjunct does dislodge the named blockers in the guided arms: at L1 identity, h_join drops from 1.000000 (boundary held by CLN3 via NPIPB8 and EIF3CL via NPIPB9) to 0.999370 under t_E and every t_J setting, where the PKD1P readthroughs take the boundary instead — the family is no closer to certifying, but a different locus is holding the line. Second, E3.2/E3.3 admit two readings for body queries; the implemented one maps a junction's donor and acceptor offsets through the CIGAR separately, and under the literal adjacent-offset reading no body record can carry a junction at all, changing 564 of 936 rows — with the decision untouched (0 rows with the shipped cut inside). The conclusion therefore does not depend on that choice, but the choice should be stated.
