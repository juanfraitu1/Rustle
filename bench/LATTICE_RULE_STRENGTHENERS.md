# Strengthening the DNA levels of the nested edge-test lattice: five strengtheners, swept, on NPIP and TBC1D3

Agent 1 of 2. Evidence reused unchanged from `/mnt/linuxdisk/home/juanfraitu/family_cert/`; no new alignment, no change to `src/`, nothing committed. Pre-registration: `/mnt/linuxdisk/home/juanfraitu/lattice_rules/DECLARATIONS.txt` (written 11:14:05 -07:00, before any strengthener was computed; addenda A1-A5 dated and disclosed, A4/A5 written after the first run and affecting S5's node definition only). Tables: `lattice_rules/sweep.tsv` (320 rows), `margins.tsv`, `nodedrops.tsv`, `S5_readthroughs.tsv`, `records.tsv` (8,673 per-record witnesses), `alt_top.tsv`, `tables.md`.

**NPIP and TBC1D3 are development families, one assembly (CHM13v2.0), one annotation (RefSeq RS_2025_08). Everything below is descriptive. Nothing here is a validation of anything.**

> **Corrections applied by the orchestrator after the independent verification** (the building agent could not write this file; its draft text is used, with these fixes): (1) eight non-empty intervals exist, not one — all TBC1D3 at L2, all above the shipped cut; (2) S2-nochain's 0.970232 is the lowest among single-strengthener rows only, C3 reaches 0.540550 at L2; (3) S1 also breaks NPIP into 3 components at L3 for θ ≥ 0.40; (4) the S4 bipartite-F claim ("rises to L = 500 then falls back") holds only for ΔF_soto at L2: at L1a ΔF_soto still rises at L = 1000 and ΔF_lit is largest at L = 2000, with members already dropped; (5) the surviving 1.000000 boundary edge under S5 is EIF3CL–NPIPB9 at L1a and MIR6511A4–NPIPB2 at L1b, not LOC100190986/LOC128966632 (those hold L3, and L2 under dropRT); (6) **undisclosed deviation**: every `bipF_lit` was computed on the Dishuck `subfamily` field (NPIPA 8 / NPIPB 18 present), not the declared `project_level1` field (NPIPA 7 / NPIPB 14); recomputed as declared the NPIP baseline is 0.195804, not 0.243243, and every lit value falls similarly — no ΔF sign and no condition changes; (7) the addendum clock times A3 11:20, A4 11:26, A5 11:30 are impossible (DECLARATIONS.txt was last written 11:19:53); precedence itself is fine, every result file postdates the frozen declaration.

## Read first

**NO CANDIDATE.** Across 5 strengtheners, 9 threshold values each where a grid applies, 4 DNA levels and 3 post-hoc combinations — 320 certificate rows — **not one produces a non-empty NPIP certificate interval containing the level's shipped cut**, so every strengthener fails condition (i) of the decision rule first, and conditions (ii)-(iv) are never reached as the binding constraint. The stronger form: **NPIP's margin `h_split − h_join` is never positive in any of the 160 (strengthener, θ, level) cells** — its best value is exactly `0.000000`, on the coverage axis L1b, where h_join and h_split both saturate at the axis ceiling 1.000000. So the failure is not a badly chosen cut: on this evidence no cut on any of these axes, under any of these strengtheners, makes NPIP an exact component.

Three readings that the sweep — not the single baseline — is what shows:

1. **The blocker is a queue, not a record.** Removing the 253-bp `gene-CLN3` fragment (S4 at L ≥ 300 bp) does not lower h_join by a single digit: EIF3CL, LOC100190986, LOC124907830/845, LOC128966632 and MIR6511A4 all sit at exactly 1.000000 behind it. Under the strongest edge-only rule the last outsider standing is **LOC128966632, a 5,598-bp protein-coding "SMG1-like" gene** — a full-size co-duplicated SD neighbour, not a small embedded record. The node question ("should CLN3 be a node?") is therefore not the question; the axis is. This is the same over-merge of co-duplicated SD neighbours that the Soto benchmark named (§6kr).
2. **Strengthening lowers h_split as fast as it lowers h_join.** Every conjunct removes evidence, and it removes the family's own weakest internal links as readily as its boundary. NPIP's margin at L1a moves from −0.0283 (baseline) to −0.0287 (S2 θ = 0.20) to −0.1655 (S2 θ ≥ 0.25) — it gets *worse* — and past θ ≈ 0.25 on S2-nochain and on every aggressive combination the family shatters (h_split = −inf) instead of separating. TBC1D3 shows the same flip: under the everything-at-once combination C3 its h_join reaches **0.000000** — the family touches nothing outside at all — but h_split is −inf, two internal parts. Over-merge becomes fragmentation; it never passes through exactness.
3. **The strengtheners cost nothing on the evaluation metrics, which is itself the warning.** Bipartite F against the Soto families and against the NPIPA|NPIPB literature groups **rises in every row and never falls** (max ΔF_soto +0.358, min +0.000). Condition (iv) is vacuous here, and an F that only improves while the certificate never appears is a reminder that F measures agreement with a truth partition, not isolation in the filtration. Prior art holds: at θ = 0.50, S1 costs 38 of NPIP's 310 L1a member edges and S2 costs 224, with the certificate no closer — the regime in which a symmetric 0.50 coverage-of-longer floor already destroyed NPIP in the de novo E_r graph (register 295/304/330).

Non-empty intervals appear in **eight rows, all TBC1D3 at L2** (verifier): S1 θ = 0.25 `(0.811940, 0.832037]`; S2 θ = 0.20 and 0.25 `(0.508554, 0.832037]`; S2-nochain θ = 0.10 `(0.513433, 0.553243]` and θ = 0.15 / 0.20 / 0.25 `(0.508554, 0.553243]`; C2 `(0.508554, 0.832037]`. All lie above the shipped L2 cut of 0.30. That interval is real and wide (0.32), but it lies **above** the shipped L2 cut of 0.30, so condition (i) still fails. If the L2 cut were a free parameter, TBC1D3's 12 DNA nodes would be an exact component of the f_ex filtration anywhere in 0.51-0.83 once S2 at θ = 0.20 is applied. That is a cut question, not a strengthener question, and it is the single most actionable thing the sweep found. NPIP has no such window at any level.

### What each strengthener fails first

| strengthener | first failing condition | detail |
|---|---|---|
| S1 two-sided coverage (θ 0.00-0.50) | **(i)** at every θ and level | NPIP h_join stays 1.000000 on L1a/L1b/L2 at every θ; only L3 drops it to 0.999349 (θ ≥ 0.30) against h_split 0.980559. At L3 with θ ≥ 0.40 the 26 members also fall into 3 components, so (ii) fails there too |
| S2 shared exon of the longer copy (θ 0.00-0.50) | **(i)** at every θ and level | best NPIP h_join 0.996905 (L1a, θ ≥ 0.20) against h_split 0.968187; L1b and L2 h_join pinned at 1.000000 throughout |
| S2-nochain (the literal D2 reading) | **(i)**, then also **(ii)** from θ = 0.25 | NPIP h_join reaches 0.970232, the lowest among single-strengthener rows (C3 reaches 0.540550 at L2), but h_split = −inf: the family is in 2-9 pieces |
| S3 reciprocity (binary) | **(i)** | NPIP h_join 1.000000 unchanged at all four levels; 0 member edges lost; component 122 → 61 at L1a |
| S4 minimum exonic length (L 0-2000 bp) | **(i)**, then also **(ii)** from L = 1000 | NPIP h_join 1.000000 at every L including 2000, where 8 of 26 NPIP members are themselves dropped |
| S5 contained records (both variants) | **(i)**, then also **(ii)** | NPIP h_join 1.000000; keepRT drops 3 NPIP members, dropRT 4 |
| C1, C2, C3 (post hoc combinations) | **(i)** | C3 shatters both families (h_split −inf) while still leaving LOC128966632 on NPIP's boundary at 0.980662 |

Condition (iii) (TBC1D3 gains no outside members) is satisfied by every strengthener at every value — TBC1D3's outside count only ever falls, by up to 19 at L1a. Condition (iv) is satisfied everywhere, for the reason given above.

## 1. What was swept, and why each form is monotone

Each strengthener is an extra **conjunct on the witness record**, so at every cut the strengthened edge set is a subset of the baseline's: nesting (T1) and H1 hold by construction, and a certificate can only be lost, never gained, by more evidence. For a pair (u, v) each level's weight is

    w_level(u, v; θ) = max { g_level(r) : r a witness record for (u, v) with P_θ(r) }

— an existential over records / a maximum over records. `P_θ(r)` depends only on r and on the two node records, so adding evidence can only enlarge the qualifying set, hence only raise the max, hence only raise h_join.

| id | conjunct on a record r | grid |
|---|---|---|
| S1 | `cov_longer(r) ≥ θ`; tx records: aligned transcript bp / **longer** exonic length; body records and chains: aligned query bp / **longer** body | 0.00 … 0.50 |
| S2 | `fex_longer(r) = sx(r) / max(u_exon_len, v_exon_len) ≥ θ`, `sx = min(u exon bp in r, v exon bp in r)` | 0.00 … 0.50 |
| S2-nochain | as S2, and chains never qualify for θ > 0 (the literal D2 reading; kept because it separates "the exon conjunct" from "drop the chains") | 0.00 … 0.50 |
| S3 | a qualifying record must exist with u as query **and** with v as query | binary |
| S4 | *node rule*: drop nodes with exon-union length < L | 0, 200, 300, 500, 1000, 2000 bp |
| S5 | *node rule*: drop a node whose span lies strictly inside another node's span on the same chromosome | keepRT / dropRT |

S3 is a conjunction of two monotone existentials, so it is monotone too. **S4 and S5 are not monotone in the same sense**: they delete nodes, and a future annotation can lengthen a node's exon union or create a container, so the node set is not a function that only grows with evidence. They are reported as node-rule tables.

**Which rows depend on the non-monotone greedy gene-body chains.** L1a and L1b are built from `tx_record` and `body_chain` witnesses, so **every L1a / L1b row inherits the chains' non-monotonicity (H3)**. L2 and L3 take their f_ex and w_98 maxima from `tx_record` and `body_record` only — but their gate is `t_1`, computed over chains, so **L2 and L3 rows depend on the chains through the gate, not through the weight**. The S2-nochain rows are the closest chain-free reading of L1 available here.

**Reciprocity testability (declared limitation).** 1,670 of 1,986 pairs (84.1%) have both endpoints among the 265 mapped query nodes and are testable; the 316 untestable edges were kept (the conservative choice, against the strengthener). Of the testable pairs, 356 fail reciprocity at L1. **0 of NPIP's 325 member-member pairs are non-reciprocal; 16 of TBC1D3's 62 are** — which is why S3 costs NPIP nothing and costs TBC1D3 its connectivity.

## 2. The curves

`h_join / h_split` per level; `comp` = component size at the shipped cut with (outsiders) and `p` = number of components the members fall into; `lost` = member pairs that had an edge at the shipped cut under the current rule and no longer do; ΔF = change in bipartite F against Soto / against NPIPA|NPIPB, relative to the same level's θ = 0 row. **θ = 0 (and L = 0) is the current rule, and reproduces `family_cert/cert/certificates.tsv` exactly** — h_join, h_split, component size and members present match on all 8 family × level rows, and all 1,986 pair weights match. The `inside?` column is L1a/L1b/L2/L3.

**S1_cov_longer — NPIP**

| θ | n | L1a h_join/h_split | L1b | L2 | L3 | comp L1a | comp L2 | comp L3 | lost L1a/L2/L3 | ΔF_soto L1a/L2/L3 | ΔF_lit L1a/L2/L3 | inside? |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.00 | 26 | 1.000000 / 0.971665 | 1.000000 / 1.000000 | 1.000000 / 0.989040 | 1.000000 / 0.980559 | 122 (96) p1 | 90 (64) p1 | 83 (57) p1 | 0/0/0 | +0.000/+0.000/+0.000 | +0.000/+0.000/+0.000 | no/no/no/no |
| 0.05 | 26 | 1.000000 / 0.971665 | 1.000000 / 1.000000 | 1.000000 / 0.989040 | 1.000000 / 0.980559 | 122 (96) p1 | 90 (64) p1 | 83 (57) p1 | 0/0/0 | +0.000/+0.000/+0.000 | +0.000/+0.000/+0.000 | no/no/no/no |
| 0.10 | 26 | 1.000000 / 0.971665 | 1.000000 / 1.000000 | 1.000000 / 0.989040 | 1.000000 / 0.980559 | 122 (96) p1 | 90 (64) p1 | 83 (57) p1 | 0/0/0 | +0.000/+0.000/+0.000 | +0.000/+0.000/+0.000 | no/no/no/no |
| 0.15 | 26 | 1.000000 / 0.971665 | 1.000000 / 1.000000 | 1.000000 / 0.989040 | 1.000000 / 0.980559 | 121 (95) p1 | 89 (63) p1 | 83 (57) p1 | 1/1/1 | +0.001/+0.002/+0.000 | +0.002/+0.003/+0.000 | no/no/no/no |
| 0.20 | 26 | 1.000000 / 0.971665 | 1.000000 / 1.000000 | 1.000000 / 0.989040 | 1.000000 / 0.980559 | 121 (95) p1 | 89 (63) p1 | 83 (57) p1 | 4/4/2 | +0.001/+0.002/+0.000 | +0.002/+0.003/+0.000 | no/no/no/no |
| 0.25 | 26 | 1.000000 / 0.971665 | 1.000000 / 1.000000 | 1.000000 / 0.989040 | 1.000000 / 0.980559 | 118 (92) p1 | 87 (61) p1 | 82 (56) p1 | 5/5/3 | +0.005/+0.006/+0.006 | +0.007/+0.008/+0.003 | no/no/no/no |
| 0.30 | 26 | 1.000000 / 0.971665 | 1.000000 / 1.000000 | 1.000000 / 0.989040 | 0.999349 / 0.980559 | 118 (92) p1 | 87 (61) p1 | 82 (56) p1 | 7/7/5 | +0.005/+0.006/+0.006 | +0.007/+0.008/+0.003 | no/no/no/no |
| 0.40 | 26 | 1.000000 / 0.942387 | 1.000000 / 1.000000 | 1.000000 / 0.981818 | 0.999349 / 0.949617 | 100 (74) p1 | 87 (61) p1 | 82 (56) p3 | 26/27/15 | +0.034/+0.006/+0.006 | +0.042/+0.008/+0.133 | no/no/no/no |
| 0.50 | 26 | 1.000000 / 0.942387 | 1.000000 / 1.000000 | 1.000000 / 0.981818 | 0.999349 / 0.949617 | 99 (73) p1 | 86 (60) p1 | 82 (56) p3 | 38/38/25 | +0.035/+0.008/+0.009 | +0.045/+0.011/+0.133 | no/no/no/no |

**S1_cov_longer — TBC1D3**

| θ | n | L1a h_join/h_split | L1b | L2 | L3 | comp L1a | comp L2 | comp L3 | lost L1a/L2/L3 | ΔF_soto L1a/L2/L3 | ΔF_lit | inside? |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.00 | 12 | 1.000000 / 0.840488 | 1.000000 / 0.841369 | 1.000000 / 0.832037 | 0.961397 / 0.867946 | 31 (19) p1 | 22 (10) p1 | 12 (0) p3 | 0/0/0 | +0.000/+0.000/+0.000 | NA | no/no/no/no |
| 0.05 | 12 | 1.000000 / 0.840488 | 1.000000 / 0.841369 | 1.000000 / 0.832037 | 0.961397 / 0.867946 | 31 (19) p1 | 22 (10) p1 | 12 (0) p3 | 0/0/0 | +0.000/+0.000/+0.000 | NA | no/no/no/no |
| 0.10 | 12 | 1.000000 / 0.840488 | 1.000000 / 0.841369 | 1.000000 / 0.832037 | 0.961397 / 0.867946 | 31 (19) p1 | 22 (10) p1 | 12 (0) p3 | 1/1/0 | +0.000/+0.000/+0.000 | NA | no/no/no/no |
| 0.15 | 12 | 1.000000 / 0.840488 | 1.000000 / 0.841369 | 1.000000 / 0.832037 | 0.961397 / 0.867946 | 31 (19) p1 | 22 (10) p1 | 12 (0) p3 | 2/2/0 | +0.000/+0.000/+0.000 | NA | no/no/no/no |
| 0.20 | 12 | 1.000000 / 0.840488 | 1.000000 / 0.841369 | 1.000000 / 0.832037 | 0.950106 / 0.867946 | 25 (13) p1 | 16 (4) p1 | 12 (0) p3 | 3/3/0 | +0.087/+0.140/+0.000 | NA | no/no/no/no |
| 0.25 | 12 | 1.000000 / 0.807623 | 1.000000 / 0.822709 | **0.811940 / 0.832037** | 0.950106 / 0.867946 | 22 (10) p1 | 13 (1) p1 | 12 (0) p3 | 4/4/0 | +0.143/+0.214/+0.000 | NA | no/no/no/no |
| 0.30 | 12 | 1.000000 / 0.800429 | 1.000000 / 0.567305 | 0.811940 / 0.405910 | 0.921421 / 0.867946 | 22 (10) p1 | 13 (1) p1 | 12 (0) p3 | 6/6/0 | +0.143/+0.214/+0.000 | NA | no/no/no/no |
| 0.40 | 12 | 1.000000 / 0.800429 | 1.000000 / 0.567305 | 0.811940 / 0.405910 | 0.921421 / 0.867946 | 22 (10) p1 | 13 (1) p1 | 12 (0) p3 | 6/6/0 | +0.143/+0.214/+0.000 | NA | no/no/no/no |
| 0.50 | 12 | 1.000000 / 0.800429 | 1.000000 / 0.567305 | 0.811940 / 0.405910 | 0.921421 / 0.867946 | 20 (8) p1 | 13 (1) p1 | 12 (0) p3 | 6/6/0 | +0.186/+0.274/+0.076 | NA | no/no/no/no |

**S2_fex_longer — NPIP**

| θ | n | L1a h_join/h_split | L1b | L2 | L3 | comp L1a | comp L2 | comp L3 | lost L1a/L2/L3 | ΔF_soto L1a/L2/L3 | ΔF_lit L1a/L2/L3 | inside? |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.00 | 26 | 1.000000 / 0.971665 | 1.000000 / 1.000000 | 1.000000 / 0.989040 | 1.000000 / 0.980559 | 122 (96) p1 | 90 (64) p1 | 83 (57) p1 | 0/0/0 | +0.000/+0.000/+0.000 | +0.000/+0.000/+0.000 | no/no/no/no |
| 0.05 | 26 | 1.000000 / 0.971665 | 1.000000 / 1.000000 | 1.000000 / 0.989040 | 1.000000 / 0.980559 | 100 (74) p1 | 68 (42) p1 | 63 (37) p1 | 0/0/0 | +0.034/+0.053/+0.056 | +0.042/+0.073/+0.074 | no/no/no/no |
| 0.10 | 26 | 1.000000 / 0.971665 | 1.000000 / 1.000000 | 1.000000 / 0.989040 | 1.000000 / 0.980559 | 76 (50) p1 | 68 (42) p1 | 63 (37) p1 | 7/4/4 | +0.084/+0.053/+0.056 | +0.110/+0.073/+0.074 | no/no/no/no |
| 0.15 | 26 | 0.999363 / 0.971665 | 1.000000 / 1.000000 | 1.000000 / 0.989040 | 1.000000 / 0.980559 | 73 (47) p1 | 65 (39) p1 | 63 (37) p1 | 40/24/13 | +0.092/+0.062/+0.056 | +0.120/+0.085/+0.074 | no/no/no/no |
| 0.20 | 26 | 0.996905 / 0.968187 | 1.000000 / 1.000000 | 1.000000 / 0.989040 | 0.999349 / 0.980559 | 73 (47) p1 | 62 (36) p1 | 60 (34) p1 | 86/56/24 | +0.092/+0.071/+0.065 | +0.120/+0.099/+0.088 | no/no/no/no |
| 0.25 | 26 | 0.996905 / 0.831362 | 1.000000 / 1.000000 | 1.000000 / 0.985305 | 0.999349 / 0.980559 | 67 (41) p1 | 58 (32) p1 | 57 (31) p1 | 129/92/45 | +0.108/+0.085/+0.078 | +0.144/+0.118/+0.103 | no/no/no/no |
| 0.30 | 26 | 0.996905 / 0.831362 | 1.000000 / 1.000000 | 1.000000 / 0.985305 | 0.999130 / 0.980168 | 56 (30) p1 | 50 (24) p1 | 49 (23) p1 | 179/140/59 | +0.144/+0.114/+0.108 | +0.196/+0.163/+0.150 | no/no/no/no |
| 0.40 | 26 | 0.996905 / 0.831362 | 1.000000 / 1.000000 | 1.000000 / 0.931572 | 0.999130 / 0.980168 | 53 (27) p1 | 49 (23) p1 | 48 (22) p1 | 208/171/72 | +0.155/+0.118/+0.112 | +0.212/+0.170/+0.156 | no/no/no/no |
| 0.50 | 26 | 0.996905 / 0.831362 | 1.000000 / 1.000000 | 1.000000 / 0.855654 | 0.998966 / 0.976856 | 46 (20) p1 | 46 (20) p1 | 44 (18) p3 | 224/186/77 | +0.184/+0.131/+0.064 | +0.257/+0.190/+0.127 | no/no/no/no |

**S2_fex_longer — TBC1D3**

| θ | n | L1a h_join/h_split | L1b | L2 | L3 | comp L1a | comp L2 | comp L3 | lost L1a/L2/L3 | ΔF_soto L1a/L2/L3 | ΔF_lit | inside? |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.00 | 12 | 1.000000 / 0.840488 | 1.000000 / 0.841369 | 1.000000 / 0.832037 | 0.961397 / 0.867946 | 31 (19) p1 | 22 (10) p1 | 12 (0) p3 | 0/0/0 | +0.000/+0.000/+0.000 | NA | no/no/no/no |
| 0.05 | 12 | 0.996612 / 0.840488 | 1.000000 / 0.841369 | 1.000000 / 0.832037 | 0.961397 / 0.867946 | 31 (19) p1 | 22 (10) p1 | 12 (0) p3 | 0/0/0 | +0.000/+0.000/+0.000 | NA | no/no/no/no |
| 0.10 | 12 | 0.996612 / 0.807623 | 1.000000 / 0.822709 | 1.000000 / 0.832037 | 0.927995 / 0.830224 | 29 (17) p1 | 21 (9) p1 | 12 (0) p3 | 2/2/0 | +0.026/+0.014/+0.000 | NA | no/no/no/no |
| 0.15 | 12 | 0.885467 / 0.807623 | 1.000000 / 0.822709 | 1.000000 / 0.832037 | 0.927995 / 0.830224 | 21 (9) p1 | 21 (9) p1 | 12 (0) p3 | 3/3/0 | +0.017/+0.014/+0.000 | NA | no/no/no/no |
| **0.20** | 12 | 0.845438 / 0.807623 | 0.838258 / 0.822709 | **0.508554 / 0.832037 → interval (0.508554, 0.832037]** | 0.881057 / 0.830224 | 13 (1) p1 | 13 (1) p1 | 12 (0) p3 | 4/4/0 | +0.166/+0.163/+0.000 | NA | no/no/no/no |
| **0.25** | 12 | 0.845438 / 0.807623 | 0.838258 / 0.822709 | **0.508554 / 0.832037 → interval (0.508554, 0.832037]** | 0.881057 / 0.830224 | 13 (1) p1 | 13 (1) p1 | 12 (0) p3 | 5/5/0 | +0.166/+0.163/+0.000 | NA | no/no/no/no |
| 0.30 | 12 | 0.845438 / -inf | 0.838258 / -inf | 0.508554 / -inf | 0.881057 / -inf | 13 (1) p2 | 13 (1) p2 | 12 (0) p3 | 9/9/0 | +0.276/+0.274/+0.076 | NA | no/no/no/no |
| 0.40 | 12 | 0.845438 / -inf | 0.838258 / -inf | 0.508554 / -inf | 0.881057 / -inf | 13 (1) p2 | 13 (1) p2 | 12 (0) p3 | 23/23/10 | +0.346/+0.343/+0.076 | NA | no/no/no/no |
| 0.50 | 12 | 0.000000 / -inf | 0.000000 / -inf | 0.000000 / -inf | 0.000000 / -inf | 12 (0) p2 | 12 (0) p2 | 12 (0) p3 | 29/29/14 | +0.346/+0.343/+0.076 | NA | no/no/no/no |

**S2nc_fex_longer_nochain — NPIP**

| θ | n | L1a h_join/h_split | L1b | L2 | L3 | comp L1a | comp L2 | comp L3 | lost L1a/L2/L3 | ΔF_soto L1a/L2/L3 | ΔF_lit L1a/L2/L3 | inside? |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.00 | 26 | 1.000000 / 0.971665 | 1.000000 / 1.000000 | 1.000000 / 0.989040 | 1.000000 / 0.980559 | 122 (96) p1 | 90 (64) p1 | 83 (57) p1 | 0/0/0 | +0.000/+0.000/+0.000 | +0.000/+0.000/+0.000 | no/no/no/no |
| 0.05 | 26 | 1.000000 / 0.934476 | 1.000000 / 1.000000 | 1.000000 / 0.862385 | 1.000000 / 0.980168 | 37 (11) p1 | 37 (11) p1 | 36 (10) p1 | 161/123/26 | +0.171/+0.119/+0.112 | +0.328/+0.261/+0.250 | no/no/no/no |
| 0.10 | 26 | 1.000000 / 0.934476 | 1.000000 / 1.000000 | 1.000000 / 0.862385 | 1.000000 / 0.980168 | 37 (11) p1 | 37 (11) p1 | 36 (10) p1 | 162/124/28 | +0.171/+0.119/+0.112 | +0.328/+0.261/+0.250 | no/no/no/no |
| 0.15 | 26 | 0.999363 / 0.934476 | 1.000000 / 1.000000 | 1.000000 / 0.862385 | 1.000000 / 0.980168 | 36 (10) p1 | 36 (10) p1 | 35 (9) p1 | 171/133/33 | +0.188/+0.135/+0.129 | +0.337/+0.270/+0.260 | no/no/no/no |
| 0.20 | 26 | 0.978223 / 0.934476 | 1.000000 / 1.000000 | 1.000000 / 0.862385 | 0.988804 / 0.980168 | 35 (9) p1 | 35 (9) p1 | 34 (8) p1 | 189/151/41 | +0.192/+0.139/+0.133 | +0.347/+0.280/+0.270 | no/no/no/no |
| 0.25 | 26 | 0.971834 / -inf | 1.000000 / -inf | 1.000000 / -inf | 0.988804 / -inf | 35 (9) p2 | 35 (9) p2 | 34 (8) p2 | 214/175/58 | +0.313/+0.260/+0.254 | +0.380/+0.313/+0.303 | no/no/no/no |
| 0.30 | 26 | 0.971834 / -inf | 1.000000 / -inf | 1.000000 / -inf | 0.988804 / -inf | 35 (9) p2 | 35 (9) p2 | 34 (8) p5 | 241/202/74 | +0.313/+0.260/+0.168 | +0.380/+0.313/+0.103 | no/no/no/no |
| 0.40 | 26 | 0.970232 / -inf | 1.000000 / -inf | 1.000000 / -inf | 0.988804 / -inf | 35 (9) p4 | 35 (9) p4 | 34 (8) p6 | 256/217/85 | +0.253/+0.201/+0.151 | +0.249/+0.181/+0.103 | no/no/no/no |
| 0.50 | 26 | 0.970232 / -inf | 1.000000 / -inf | 1.000000 / -inf | 0.987970 / -inf | 35 (9) p6 | 35 (9) p6 | 34 (8) p9 | 262/223/92 | +0.157/+0.105/+0.076 | +0.183/+0.116/+0.103 | no/no/no/no |

**S2nc_fex_longer_nochain — TBC1D3**

| θ | n | L1a h_join/h_split | L1b | L2 | L3 | comp L1a | comp L2 | comp L3 | lost L1a/L2/L3 | ΔF_soto L1a/L2/L3 | ΔF_lit | inside? |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.00 | 12 | 1.000000 / 0.840488 | 1.000000 / 0.841369 | 1.000000 / 0.832037 | 0.961397 / 0.867946 | 31 (19) p1 | 22 (10) p1 | 12 (0) p3 | 0/0/0 | +0.000/+0.000/+0.000 | NA | no/no/no/no |
| 0.05 | 12 | 0.988213 / 0.823558 | 1.000000 / 0.841369 | 0.811940 / 0.553243 | 0.961397 / 0.867946 | 17 (5) p1 | 13 (1) p1 | 12 (0) p3 | 8/8/0 | +0.262/+0.214/+0.076 | NA | no/no/no/no |
| **0.10** | 12 | 0.988213 / 0.807623 | 1.000000 / 0.647746 | **0.513433 / 0.553243 → (0.513433, 0.553243]** | 0.896263 / 0.830224 | 16 (4) p1 | 13 (1) p1 | 12 (0) p3 | 10/10/0 | +0.291/+0.214/+0.076 | NA | no/no/no/no |
| **0.15** | 12 | 0.845438 / 0.807623 | 0.838258 / 0.647746 | **0.508554 / 0.553243** | 0.881057 / 0.830224 | 13 (1) p1 | 13 (1) p1 | 12 (0) p3 | 11/11/0 | +0.217/+0.214/+0.076 | NA | no/no/no/no |
| **0.20** | 12 | 0.845438 / 0.807623 | 0.838258 / 0.647746 | **0.508554 / 0.553243** | 0.881057 / 0.830224 | 13 (1) p1 | 13 (1) p1 | 12 (0) p3 | 12/12/0 | +0.217/+0.214/+0.076 | NA | no/no/no/no |
| **0.25** | 12 | 0.845438 / 0.807623 | 0.838258 / 0.647746 | **0.508554 / 0.553243** | 0.881057 / 0.830224 | 13 (1) p1 | 13 (1) p1 | 12 (0) p3 | 13/13/0 | +0.217/+0.214/+0.076 | NA | no/no/no/no |
| 0.30 | 12 | 0.845438 / -inf | 0.838258 / -inf | 0.508554 / -inf | 0.881057 / -inf | 13 (1) p2 | 13 (1) p2 | 12 (0) p3 | 16/16/0 | +0.346/+0.343/+0.076 | NA | no/no/no/no |
| 0.40 | 12 | 0.845438 / -inf | 0.838258 / -inf | 0.508554 / -inf | 0.881057 / -inf | 13 (1) p4 | 13 (1) p4 | 12 (0) p5 | 33/33/17 | +0.172/+0.169/−0.106 | NA | no/no/no/no |
| 0.50 | 12 | 0.000000 / -inf | 0.000000 / -inf | 0.000000 / -inf | 0.000000 / -inf | 12 (0) p5 | 12 (0) p5 | 12 (0) p5 | 34/34/17 | +0.203/+0.201/−0.106 | NA | no/no/no/no |

**S3_reciprocal (binary)**

| family | n | L1a h_join/h_split | L1b | L2 | L3 | comp L1a | comp L2 | comp L3 | lost L1a/L2/L3 | ΔF_soto L1a/L2/L3 | ΔF_lit L1a/L2/L3 | inside? |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| NPIP | 26 | 1.000000 / 0.971665 | 1.000000 / 1.000000 | 1.000000 / 0.989040 | 1.000000 / 0.980559 | 61 (35) p1 | 59 (33) p1 | 57 (31) p1 | 0/0/0 | +0.127/+0.081/+0.075 | +0.171/+0.113/+0.103 | no/no/no/no |
| TBC1D3 | 12 | 1.000000 / -inf | 1.000000 / -inf | 0.811940 / -inf | 0.961397 / -inf | 21 (9) p2 | 13 (1) p2 | 12 (0) p3 | 16/16/0 | +0.235/+0.243/+0.000 | NA | no/no/no/no |

**S4_min_exonlen — NPIP**

| L | n | L1a h_join/h_split | L1b | L2 | L3 | comp L1a | comp L2 | comp L3 | lost L1a/L2/L3 | ΔF_soto L1a/L2/L3 | ΔF_lit L1a/L2/L3 | inside? |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 26 | 1.000000 / 0.971665 | 1.000000 / 1.000000 | 1.000000 / 0.989040 | 1.000000 / 0.980559 | 122 (96) p1 | 90 (64) p1 | 83 (57) p1 | 0/0/0 | +0.000/+0.000/+0.000 | +0.000/+0.000/+0.000 | no/no/no/no |
| 200 | 26 | 1.000000 / 0.971665 | 1.000000 / 1.000000 | 1.000000 / 0.989040 | 1.000000 / 0.980559 | 86 (60) p1 | 68 (42) p1 | 63 (37) p1 | 0/0/0 | +0.061/+0.053/+0.056 | +0.078/+0.073/+0.074 | no/no/no/no |
| 300 | 26 | 1.000000 / 0.971665 | 1.000000 / 1.000000 | 1.000000 / 0.989040 | 1.000000 / 0.980559 | 85 (59) p1 | 57 (31) p1 | 55 (29) p1 | 0/0/0 | +0.063/+0.088/+0.082 | +0.081/+0.123/+0.114 | no/no/no/no |
| 500 | 26 | 1.000000 / 0.971665 | 1.000000 / 1.000000 | 1.000000 / 0.989040 | 1.000000 / 0.980559 | 70 (44) p1 | 54 (28) p1 | 52 (26) p1 | 0/0/0 | +0.100/+0.099/+0.093 | +0.132/+0.140/+0.131 | no/no/no/no |
| 1000 | 24 | 1.000000 / 0.971665 | 1.000000 / 1.000000 | 1.000000 / 0.936113 | 1.000000 / 0.980559 | 64 (40) p1 | 51 (27) p1 | 49 (25) p1 | 48/48/25 | +0.105/+0.098/+0.091 | +0.120/+0.116/+0.108 | no/no/no/no |
| 2000 | 18 | 1.000000 / 0.958290 | 1.000000 / 1.000000 | 1.000000 / 0.880795 | 1.000000 / 0.979280 | 49 (31) p1 | 42 (24) p1 | 40 (22) p2 | 171/171/71 | +0.086/+0.060/+0.126 | +0.175/+0.156/+0.290 | no/no/no/no |

**S4_min_exonlen — TBC1D3**

| L | n | L1a h_join/h_split | L1b | L2 | L3 | comp L1a | comp L2 | comp L3 | lost | ΔF_soto L1a/L2/L3 | ΔF_lit | inside? |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 12 | 1.000000 / 0.840488 | 1.000000 / 0.841369 | 1.000000 / 0.832037 | 0.961397 / 0.867946 | 31 (19) p1 | 22 (10) p1 | 12 (0) p3 | 0/0/0 | +0.000/+0.000/+0.000 | NA | no/no/no/no |
| 200 | 12 | 1.000000 / 0.840488 | 1.000000 / 0.841369 | 1.000000 / 0.832037 | 0.961397 / 0.867946 | 31 (19) p1 | 22 (10) p1 | 12 (0) p3 | 0/0/0 | +0.000/+0.000/+0.000 | NA | no/no/no/no |
| 300 | 12 | 1.000000 / 0.840488 | 1.000000 / 0.841369 | 1.000000 / 0.832037 | 0.961397 / 0.867946 | 31 (19) p1 | 22 (10) p1 | 12 (0) p3 | 0/0/0 | +0.000/+0.000/+0.000 | NA | no/no/no/no |
| 500 | 12 | 1.000000 / 0.840488 | 1.000000 / 0.841369 | 1.000000 / 0.832037 | 0.961397 / 0.867946 | 31 (19) p1 | 22 (10) p1 | 12 (0) p3 | 0/0/0 | +0.000/+0.000/+0.000 | NA | no/no/no/no |
| 1000 | 11 | 1.000000 / 0.840488 | 1.000000 / 0.822709 | 1.000000 / 0.861255 | 0.927995 / 0.886144 | 28 (17) p1 | 20 (9) p1 | 11 (0) p2 | 7/7/0 | +0.040/+0.029/+0.000 | NA | no/no/no/no |
| 2000 | 11 | 0.957267 / 0.840488 | 1.000000 / 0.822709 | 1.000000 / 0.861255 | 0.927995 / 0.886144 | 18 (7) p1 | 17 (6) p1 | 11 (0) p2 | 7/7/0 | +0.190/+0.214/+0.167 | NA | no/no/no/no |

**S5 contained records (binary)**

| variant | family | n | L1a h_join/h_split | L1b | L2 | L3 | comp L1a | comp L2 | comp L3 | lost L1a/L2/L3 | ΔF_soto L1a/L2/L3 | ΔF_lit L1a/L2/L3 | inside? |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| keepRT | NPIP | 23 | 1.000000 / 0.969630 | 1.000000 / 1.000000 | 1.000000 / 0.989040 | 1.000000 / 0.980559 | 79 (56) p1 | 58 (35) p1 | 55 (32) p1 | 69/50/31 | +0.084/+0.095/+0.089 | +0.051/+0.060/+0.054 | no/no/no/no |
| keepRT | TBC1D3 | 12 | 1.000000 / 0.840488 | 1.000000 / 0.841369 | 1.000000 / 0.832037 | 0.961397 / 0.867946 | 29 (17) p1 | 20 (8) p1 | 12 (0) p3 | 0/0/0 | +0.026/+0.029/+0.000 | NA | no/no/no/no |
| dropRT | NPIP | 22 | 1.000000 / 0.969630 | 1.000000 / 1.000000 | 1.000000 / 0.989040 | 1.000000 / 0.980559 | 70 (48) p1 | 46 (24) p1 | 46 (24) p1 | 79/60/38 | +0.135/+0.215/+0.209 | +0.083/+0.131/+0.111 | no/no/no/no |
| dropRT | TBC1D3 | 12 | 1.000000 / 0.840488 | 1.000000 / 0.841369 | 1.000000 / 0.832037 | 0.961397 / 0.867946 | 24 (12) p1 | 20 (8) p1 | 12 (0) p3 | 0/0/0 | +0.064/+0.140/+0.167 | NA | no/no/no/no |

**Combinations (POST HOC — see §4 for how they were chosen)**

| combo | family | n | L1a h_join/h_split | L1b | L2 | L3 | comp L1a | comp L2 | comp L3 | lost L1a/L2/L3 | ΔF_soto L1a/L2/L3 | ΔF_lit L1a/L2/L3 | inside? |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| C1 = S1@0.30 ∧ S2@0.20 | NPIP | 26 | 0.996905 / 0.968187 | 1.000000 / 1.000000 | 1.000000 / 0.989040 | 0.999349 / 0.980559 | 73 (47) p1 | 61 (35) p1 | 59 (33) p1 | 88/58/26 | +0.092/+0.074/+0.072 | +0.120/+0.103/+0.093 | no/no/no/no |
| C1 | TBC1D3 | 12 | 0.845438 / -inf | 0.838258 / -inf | 0.482115 / -inf | 0.858223 / -inf | 13 (1) p2 | 13 (1) p2 | 12 (0) p3 | 7/7/0 | +0.217/+0.274/+0.000 | NA | no/no/no/no |
| C2 = S2@0.20 ∧ S4@500 | NPIP | 26 | 0.996905 / 0.968187 | 1.000000 / 1.000000 | 1.000000 / 0.989040 | 0.999349 / 0.980559 | 61 (35) p1 | 54 (28) p1 | 52 (26) p1 | 86/56/24 | +0.127/+0.099/+0.093 | +0.171/+0.140/+0.131 | no/no/no/no |
| C2 | TBC1D3 | 12 | 0.845438 / 0.807623 | 0.838258 / 0.822709 | **0.508554 / 0.832037 → (0.508554, 0.832037]** | 0.881057 / 0.830224 | 13 (1) p1 | 13 (1) p1 | 12 (0) p3 | 4/4/0 | +0.166/+0.163/+0.000 | NA | no/no/no/no |
| C3 = everything | NPIP | 20 | 0.980662 / -inf | 1.000000 / -inf | 0.540550 / -inf | 0.996395 / -inf | 31 (11) p2 | 31 (11) p2 | 31 (11) p4 | 249/210/91 | +0.358/+0.306/+0.091 | +0.227/+0.160/+0.376 | no/no/no/no |
| C3 | TBC1D3 | 11 | **0.000000** / -inf | 0.000000 / -inf | 0.000000 / -inf | 0.000000 / -inf | 11 (0) p2 | 11 (0) p2 | 11 (0) p2 | 17/17/0 | +0.476/+0.474/+0.167 | NA | no/no/no/no |

## 3. The node rules, and what they drop

`nodedrops.tsv`. "Blockers" are the outside nodes that held weight 1.000000 (or the top boundary weight) on a family in the baseline certificate.

| rule | param | nodes dropped (of 47,965) | NPIP members dropped | TBC1D3 members dropped | blockers dropped |
|---|---|---|---|---|---|
| S4 min exon-union length | 0 | 0 | — | — | — |
| | 200 | 4,855 | — | — | — |
| | 300 | 5,896 | — | — | CLN3 |
| | 500 | 9,127 | — | — | CLN3 |
| | 1000 | 16,457 | NPIPB10P, NPIPB1P | LOC124905656 | CLN3, LOC124907830, TBC1D29P |
| | 2000 | 25,515 | NPIPA1, NPIPA6, NPIPA7, NPIPA8, NPIPB10P, NPIPB1P, NPIPB2, NPIPB8 | LOC124905656 | CLN3, LOC124907830, LOC124907845, NPEPPSP1, TBC1D29P, TBC1D3P1-DHX40P1 |
| S5 contained, keepRT | — | 10,618 | LOC124907834, LOC128966608, NPIPB4 | — | CLN3, LOC124907830, LOC124907845 |
| S5 contained, dropRT | — | 10,827 | LOC124907834, LOC128966608, NPIPB4, PKD1P6-NPIPP1 | — | CLN3, LOC124907830, LOC124907845, PKD1P5-LOC105376752, LOC131696449, TBC1D3P1-DHX40P1 |

**The residual blocker, by name and size** (this is the point of the section):

| node | biotype | exon union | span | what it is |
|---|---|---|---|---|
| `gene-CLN3` | other | **253 bp** | 253 bp, chr16:28,747,319 | the single-exon fragment; *not* the real CLN3, which is `gene-CLN3-2`, protein_coding, 4,079 bp, chr16:28,903,666 |
| `gene-EIF3CL` | protein_coding | 3,625 bp | 47.1 kb | survives every node rule tried |
| `gene-LOC100190986` | lncRNA | 2,453 bp | 2.5 kb | survives S4 up to L = 2000 |
| `gene-LOC124907830` | lncRNA | 768 bp | 1.5 kb | dropped at L = 1000 and by S5 |
| `gene-LOC124907845` | lncRNA | 1,106 bp | 1.9 kb | dropped at L = 2000 and by S5 |
| `gene-LOC128966632` | **protein_coding** | **5,598 bp** | 44.5 kb | "SMG1-like"; **survives every rule tested**, and is the last node on NPIP's boundary under C3 (0.980662) |
| `gene-LOC131696449` | ncRNA_pseudogene | 7,747 bp | 43.3 kb | PKD1P1-NPIPA5L readthrough |
| `gene-MIR6511A4` | miRNA | 67 bp | 67 bp | dropped at L = 200 |

Readthrough set for S5 = the **209** nodes whose RefSeq record carries a readthrough description; this reproduces `family_cert/dna/nodes.tsv`'s own `readthrough` column exactly. The hyphenated-name fallback declared in D2 was dropped (addendum A4) because it selected 5,212 nodes, almost all antisense/divergent transcripts (`A1BG-AS1`, `ACAD9-DT`).

## 4. The two combinations (POST HOC, and labelled)

The declared selection procedure was "the two (strengthener, θ) pairs with the smallest NPIP h_join at any level, ties broken by larger h_split, intersected". Applied literally that picks S2-nochain at θ = 0.40 (h_join 0.970232) and θ = 0.25 (0.971834), both rows where NPIP is *already shattered* (h_split = −inf) and so trivially cannot satisfy (ii) — a degenerate answer. I therefore applied the same procedure **restricted to rows where NPIP keeps all 26 members in one component**, and say so; this is a disclosed deviation from the letter of the declaration.

- **C1 = S1 @ 0.30 ∧ S2 @ 0.20** — the two edge conjuncts that individually moved NPIP's h_join most while leaving it connected.
- **C2 = S2 @ 0.20 ∧ S4 @ 500** — the same exon conjunct plus the node rule that removes CLN3 without losing a member. This is the best all-round row in the study: NPIP keeps 26 members in 1 component, its component falls 122 → 61, bipF rises +0.127 / +0.171, and TBC1D3 carries the L2 interval.
- **C3 = S1 @ 0.50 ∧ S2 @ 0.30 ∧ S3 ∧ S4 @ 1000 ∧ S5-dropRT** — an extra, descriptive ceiling row, also post hoc. It is the strongest statement of the negative result: with every rule applied at once, NPIP is down to 20 members in 2 parts and *still* has LOC128966632 on its boundary at 0.980662, while TBC1D3 achieves perfect isolation (h_join 0.000000, 0 outsiders) and is in 2 pieces.

## 5. The alternative top level (descriptive, separate)

`alt_top.tsv`. Protein edges at the shipped cover cut 0.30: 397,031 gene-level pairs; DNA `t_1` edges: 1,337; **meet: 316**; join: 398,052.

| top level | NPIP: seeds / component / outside / parts | TBC1D3: seeds / component / outside / parts | protein certificate |
|---|---|---|---|
| (a) current L0 = t_P **OR** t_1 | 26 / **358** / 332 / 1 | 12 / **5,450** / 5,438 / 1 | **destroyed** — the top-level component is nowhere near the 21 NPIP proteins |
| (b) the meet t_P **AND** t_1 | 26 / **26** / **0** / **6** | 12 / **12** / **0** / **4** | not intact at the top — h_join = 0 (the only place in the study where NPIP touches nothing outside) but h_split = −inf: perfectly precise, shattered |
| (c) two parallel stacks (protein chain; DNA chain L1-L3), no join | protein 21 / 21 / **0** / 1, **exact**; DNA 26 / 122 / 96 / 1 | protein 9 / 5,275 / 5,266 / 1; DNA 12 / 31 / 19 / 1 | **intact** — the protein stack's top *is* the protein graph, so the certified interval (0.128250, 0.840580] with the shipped 0.30 inside survives verbatim |

**Only (c) keeps the protein certificate intact.** (a) destroys it outright. (b) is the interesting failure: the meet is the only construction here that gives NPIP h_join = 0 — it admits no outsider at all — but it fragments the family into 6 pieces, the same over-merge → fragmentation flip the aggressive strengtheners produce. (c) does not fix TBC1D3 either; it simply declines to let the DNA side contaminate the protein side.

## 6. Readings

- **The DNA axes cannot isolate NPIP, and it is not a threshold problem.** h_split ≤ h_join in all 160 cells. Whatever is chosen as the cut, and whatever conjunct is bolted onto the witness, some outside record is at least as close to an NPIP copy as NPIP's weakest internal link is. The next step, if the lattice is to gain a certifying DNA level, has to change *what the axis measures* — not tighten it.
- **The coverage axis L1b is structurally incapable of certifying anything here.** Its weight is `min(coverage, 1.0)`, and both the boundary and the internal bottleneck saturate at 1.000000 in every single row, giving the empty interval (1.0, 1.0]. A bounded axis whose values pile up at its cap cannot separate; that is worth recording independently of NPIP.
- **The one real opening is TBC1D3 at L2 with S2 at θ ≈ 0.20, if the L2 cut may move.** `(0.508554, 0.832037]` is a 0.32-wide window in which TBC1D3's 12 DNA nodes are an exact component with exactly one outsider excluded (TBC1D29P, via LOC124905656). The shipped cut 0.30 sits below it. Nothing comparable exists for NPIP.
- **S3 is free on NPIP and expensive on TBC1D3, and that asymmetry is informative.** 0 of NPIP's 325 member-member pairs are non-reciprocal; 16 of TBC1D3's 62 are. NPIP's copies mutually recover each other; TBC1D3's do not. That is a property of the families, visible at no cost, and it is the kind of statement the certificate machinery was built to make.
- **Do not read the bipartite-F gains as progress.** F rises monotonically almost everywhere, peaking at +0.358 (C3, Soto, L1a) — while the certificate gets no closer and the family falls apart. On this substrate F and the certificate are not measuring the same failure.

## 7. Provenance

- **Evidence.** `family_cert/dna/nodes.tsv` (47,965 primary nodes + 10,598 exon-less records), `family_cert/dna/batches/*.paf` (265 query nodes, tx `-x splice -uf -c -N 50 -p 0.1` and body `-x asm20 -c -N 50 -p 0.1`, prebuilt indexes untouched), `family_cert/protein/pairs.tsv` (575,198 pairs, blastp 2.17.0+, fixed `-dbsize 11710993`). **No new alignment was run.**
- **Recompute.** `lattice_rules/records.py` re-derived the witness records from the existing PAFs, mirroring `cert/dna_edges.py cmd_extended` and adding the per-record two-sided quantities (`cov_longer`, `ux`, `vx`, `sx`, `fex_longer`): **8,673 records, 1.0 s**, 0 same-locus rows skipped.
- **Baseline check** (`engine.py check`): the θ = 0 aggregate reproduces `cert/dna_pairs.extended.tsv` on **1,986 / 1,986 primary pairs with 0 mismatches** in `id_w`, `cov_w`, `t1`, `fex`, `w98`; and reproduces `cert/certificates.tsv` on all 8 NPIP/TBC1D3 × level rows — h_join, h_split, component size and members present all identical.
- **Evaluation truth.** Soto: `layer_order/npip_tbc1d3/light/truth_soto_families.tsv`, 12 family_ids over 89 RefSeq genes (ID_154 is the NPIP family, 21 genes). Literature: `docs/lit_subclusters_npip_dishuck_check.tsv`, `project_level1` = NPIPA | NPIPB; Iso-Seq `level2` groups reported as `bipF_iso` in `sweep.tsv`. Bipartite F = one-to-one Jaccard matching (`scipy.optimize.linear_sum_assignment`), per family, universe = truth members present ∪ all nodes of every component containing one, precision over predicted sizes so that over-merging is penalised. **Evaluation only — never part of any rule.**
- **Pre-registration.** `DECLARATIONS.txt` written 11:14:05 -07:00 (md5 `db38fd04…`), before `records.tsv` existed. Addenda: A1/A2 11:16 (chain `sx`, before any result); A3 11:20 (per-family bipartite F, before any result); **A4 11:26 and A5 11:30, after the first sweep run** — both disclosed, both affecting the S5 readthrough/containment definition only; the first run's S5 rows were discarded and recomputed. `DECLARATIONS.md5` / `.time` record each version.

## 8. Caveats

- **Development families, one assembly, one annotation. Descriptive, not validation.** A strengthener that "qualified" here would still be unvalidated; none did.
- **DNA component sizes are lower bounds.** Evidence from a node no query ever touched is assumed absent (family_cert D3), so a strengthener that appears to shrink a component may only be shrinking the visible part.
- **Chain dependence.** Every L1a/L1b row depends on the non-monotone greedy gene-body chains; every L2/L3 row depends on them through the `t_1` gate. The S2-nochain rows are the chain-free reading; they are the rows that shatter NPIP soonest.
- **S4/S5 are node rules and are not monotone in evidence.** They do not carry the H1 guarantee the record conjuncts do.
- **The grid ends at 0.50 and at L = 2000 bp.** Both edges are already destructive (S2 at 0.50 costs 224 of NPIP's 310 L1a member edges; L = 2000 removes 8 of 26 NPIP members), so the curves say nothing about what lies beyond them — and there is no reason from these curves to look.
- **Condition (iv) never binds**, so this study gives no information about the cost side of the decision rule. It should not be read as evidence that the strengtheners are cheap in general.
- **Correction to `bench/FAMILY_CERTIFICATES_NPIP_TBC1D3.md`:** there are two nodes named CLN3. The blocker is `gene-CLN3` (biotype `other`, 253 bp, chr16:28,747,319), as that report said; `gene-CLN3-2` is the real 4,079-bp protein-coding CLN3 at chr16:28,903,666 and is not involved.
## Verification (independent recompute)

Agent 2 of 2 (independent verifier), 2026-09-17. No builder script under `/mnt/linuxdisk/home/juanfraitu/lattice_rules/` was read (only `DECLARATIONS.txt` and the data outputs `sweep.tsv`, `margins.tsv`, `tables.md`, `alt_top.tsv`, `nodedrops.tsv`, `S5_readthroughs.tsv`, `records.tsv`). The verifier's own code is in `/mnt/linuxdisk/home/juanfraitu/lattice_rules/verify/` (`wit.py`, `sweep.py`, `xcheck_wit.py`, `xcheck_pairs.py`, `diff.py`, `mono.py`, `alt.py`, `spot.py`), building on `family_cert/verify/` (independently verified earlier). No new alignment was run; everything is re-derived from the PAFs, `dna/nodes.tsv` and `dna/witnesses.tsv` already on disk. The report `bench/LATTICE_RULE_STRENGTHENERS.md` did not exist when this was written, so nothing could be corrected in place; the corrections below must be applied to the report text before it is published.

**1. Order of declarations.** `DECLARATIONS.txt` was hash-stamped four times (11:14:05, 11:14:45, 11:19:08, 11:19:53 -07:00); its current md5 `68497dce…` equals the last stamp, so it has not been touched since 11:19:53, which precedes every result file (`sweep.tsv` 11:21:39, `alt_top.tsv` 11:22:09, `margins.tsv` 11:24:21, `tables.md` 11:24:41). A1/A2 (11:16) sit inside the stamp trail. **A3 (11:20), A4 (11:26) and A5 (11:30) carry clock times that postdate the file's own last write (11:19:53) and are therefore wrong**; the addenda were in fact in place before any surviving result was produced. A4's disclosure that it followed a *discarded* first sweep run cannot be checked, because that run's outputs were overwritten; it is taken on the record as disclosed.

**2. Full recompute of the sweep.** Records were re-derived from the PAFs with the verifier's own chain, exon-block, sense, sx, coverage-of-longer and shared-exon-of-longer code: 7,822 primary witness rows over 1,986 primary pairs. Cross-checks: all 4,575 primary member-endpoint rows of `dna/witnesses.tsv` reproduced exactly for identity, gap-excluded identity, coverage, target exon overlap, `sx_frac`, chain shipped-pass, and for the two strengthener quantities recomputed from the witness columns (`cov_longer`, `fex_longer`) — the only difference is `sense_ok` on `body_record` rows, a field the aggregation never reads. The baseline aggregation reproduces the shipped pair tables `cert/dna_pairs.members.tsv` / `.extended.tsv` (`id_w`, `cov_w`, `t1`, `fex`, `w98`, `n_rows`) with no substantive mismatch (47 pairs print `NA` where the verifier prints `fex = 0.0`; identical downstream).

Recomputing all 320 certificate rows (5 strengtheners × 9 grid values × 4 levels, plus S3, the two node rules and the three combinations) gives **0 mismatches against `sweep.tsv` in `n_present`, `h_join`, `h_split`, interval emptiness, `shipped_inside`, `comp_size`, `comp_outside`, `comp_parts`, `member_edges`, `member_edges_lost`, `bipF_lit`, `bipF_iso`**, and 0 mismatches in the 160 NPIP margins of `margins.tsv`. `nodedrops.tsv` is reproduced exactly (S4: 0 / 4,855 / 5,896 / 9,127 / 16,457 / 25,515 nodes; CLN3 (253 bp) leaves at L = 300; NPIPB10P, NPIPB1P at 1,000; eight NPIP members at 2,000. S5: 10,618 keepRT / 10,827 dropRT, with exactly the listed member and blocker casualties), and the 209-node readthrough set of `S5_readthroughs.tsv` is set-identical to an independent "GFF description contains readthrough" selection. S3's counts reproduce: 1,986 pairs, 1,670 testable (84.1%), 316 untestable edges kept, NPIP 0 member edges lost, TBC1D3 16 of 62.

Two readings had to be fixed to reproduce the builder's numbers, and both should be stated in the report: (i) in S3, "a qualifying record exists in BOTH query directions" must mean a **level-candidate** record (sense-ok `tx_record`, or a chain with ≥ 1 exon base on the target); under the looser "any witness record" reading NPIP's L1a/L1b component is 65 (39) rather than 61 (35), with identical `h_join`, `h_split` and member-edge losses. (ii) in C3, reciprocity must be evaluated on the records that survive S1/S2/S4/S5; evaluated on the baseline record set C3's NPIP row reads 0.996905 / 0.942387, one part, 32 (12) instead of 0.980662 / −inf, two parts, 31 (11). NPIP is uncertified under both.

**3. Monotonicity.** By construction each of S1, S2, S2-nochain is a per-record conjunct under a maximum over records, so the qualifying set, every level weight and `h_join` can only rise when records are added; S3 is a conjunction of two such existentials. Empirically: 25 (strengthener, θ) settings × 2,500 random instances = 62,500 record-addition tests, plus 2,500 for S3 — **0 failures** (no edge ever disappeared and no weight ever fell). Two caveats the declaration gets wrong: **S1 and S2 are monotone in records but not under annotation growth**, because `cov_longer` and `fex_longer` divide by `max(u, v)` exon-union / body length and a node's exon union only grows as annotation is added — a 20 % longer exon union on the longer side disqualifies 334 of 3,621 records that qualify at S1 θ = 0.30 and 303 of 2,379 at S2 θ = 0.20; conversely **S4 is monotone under annotation growth** (exon_bp only rises, so a retained node stays retained), so declaring it non-monotone alongside S5 is conservative but inaccurate. S5's containment rule is genuinely non-monotone (a newly annotated longer gene can contain an existing node and delete all its edges).

**4. Alternative tops and the protein certificate.** Recomputed from `protein/pairs.tsv` and the DNA `t_1` predicate: the join `t_P` ∨ `t_1` has 398,052 gene edges, the meet `t_P` ∧ `t_1` has 316 — both exact. (a) NPIP 358 nodes, 332 outside, 1 part; TBC1D3 5,450 / 5,438 / 1 part. (b) meet: NPIP exactly 26 with 0 outside in 6 parts; TBC1D3 exactly 12 with 0 outside in 4 parts. (c) protein-only stack: NPIP = 21 proteins, 0 outside, 1 part, certificate interval (0.128250, 0.840580] containing the shipped 0.30; TBC1D3's protein component 5,275 (5,266 outside); DNA-only NPIP 122 (96), TBC1D3 31 (19). Every number in `alt_top.tsv` and the protein-certificate claim is confirmed, including that the meet leaves NPIP with `h_split = −inf` (6 internal parts) and so is not exact either.

**5. Boundary-edge spot checks.** Forty boundary edges (the top five per level for each family) were traced back to the raw PAF lines. Examples: NPIPB8–CLN3 at L1a, tx record `chr16:28928821-28929074`, `nm = bl = 253` → identity exactly 1.000000; NPIPA8–PKD1P2, `nm/bl = 1587/1588 = 0.999370`; NPIPB5–LOC100190986 at L3, `nm = 2450`, M-bases 2450 → `w_98 = 1.000000` with `nm/bl = 0.998370`. Every stored identity, gap-excluded identity and coverage matched the arithmetic on the raw record.

**6. Verdict.** The headline result stands on an independent recompute: **no strengthener produces a non-empty NPIP certificate interval containing the shipped cut at any level or grid value**, condition (i) fails first in all 160 cells, and NPIP's margin `h_split − h_join` is never positive — its maximum over the whole study is exactly 0.000000, on the L1b coverage axis where both quantities saturate at 1.000000. The sentences listed below are nevertheless wrong or imprecise and must be fixed.

**Corrections to the report text**

1. *"The only non-empty interval anywhere is TBC1D3 at L2 under S2 θ = 0.20/0.25 (and C2)"* — wrong. `sweep.tsv` itself holds **eight** non-empty intervals, all TBC1D3 at L2: S1 θ = 0.25 (0.811940, 0.832037]; S2 θ = 0.20 and 0.25 (0.508554, 0.832037]; S2-nochain θ = 0.10 (0.513433, 0.553243] and θ = 0.15/0.20/0.25 (0.508554, 0.553243]; C2 (0.508554, 0.832037]. All lie above the shipped 0.30 cut, so the decision is unaffected.
2. *"S2-nochain … gets NPIP h_join down to 0.970232 — the lowest anywhere"* — only the lowest among single-strengthener rows (L1a). Across the study C3 reaches NPIP `h_join` 0.540550 at L2.
3. *"S1 … NPIP stays in one component (26 members) throughout"* — false at L3 for θ ≥ 0.40, where the 26 members fall into three components at the 0.98 cut (`h_split` 0.949617); the S1 table prints "82 (56) p3" there. S1 fails (ii) as well as (i) in those two cells.
4. *"S4 … bipF rises monotonically to L = 500 then falls back once members start dropping"* — false. At L1a ΔF_soto is still rising at L = 1000 (+0.105 vs +0.100) although two members have gone, and ΔF_lit peaks at L = 2000 (+0.175) with eight members gone. Only ΔF_soto at L2 peaks at L = 500.
5. *"[S5] … still leave h_join at 1.000000 via LOC100190986 and LOC128966632"* — level-specific. Those two hold at L3 (LOC100190986, both variants) and at L2 for dropRT; at L1a the surviving 1.000000 boundary edge is EIF3CL–NPIPB9 in both variants, and at L1b it is MIR6511A4–NPIPB2.
6. **bipF_lit is not the declared metric.** D4(b) declares the literature truth as the `project_level1` field (NPIPA 7 / NPIPB 14 present genes, five genes unlabelled); every reported `bipF_lit` was in fact computed on the `dishuck_group` **subfamily** field (NPIPA 8 / NPIPB 18, all 26 present genes). Recomputed as declared, baseline NPIP `bipF_lit` = 0.195804, not 0.243243, and every other lit value falls similarly. No ΔF sign and no condition changes (iv never binds), but the column must be relabelled or recomputed, and the deviation disclosed.
7. Addenda A3/A4/A5 carry impossible clock times (11:20 / 11:26 / 11:30) that postdate the last write of `DECLARATIONS.txt` (11:19:53); fix the stamps.

**Not reproduced (evaluation only).** `bipF_soto` differs in the eight rows where a node rule deletes Soto truth genes — S5_dropRT (builder 0.359551 / 0.492308 at L1a/L2, verifier 0.344086 / 0.463768) and C3 (0.583333 / 0.375000 vs 0.571429 / 0.392857). The difference lies entirely in how truth genes removed by the node rule are counted in the A3 universe; every other bipF value in all 320 rows agrees to 1e-6, the sign of each ΔF is unchanged, and condition (iv) is never the binding constraint. The report should not quote those four ΔF_soto figures to three decimals without saying which convention produced them.