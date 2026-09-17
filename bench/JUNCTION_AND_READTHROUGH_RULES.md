# Adopting the two proposed definition changes: what each form costs, and why none of them is adoptable

Agent 1 of 2. Declarations frozen **before any number**: `/mnt/linuxdisk/home/juanfraitu/rule_adopt/DECLARATIONS.txt`, written 2026-09-17T14:47:18-07:00, md5 at that moment `06897978bdcbd5601a122eccbdad5d31` (recorded in `DECLARATIONS.md5`/`DECLARATIONS.time`). Addenda A1-A3 are appended, dated, and each says whether it precedes or follows the numbers it concerns. Final file md5 `089cde5f12a3b675bd9606294ca95424`. No new alignment, nothing under `/mnt/c/Users/jfris/Desktop/Rustle` touched, nothing committed, no subagent.

**NPIP and TBC1D3 are development families. One assembly (CHM13v2.0), one annotation (RefSeq RS_2025_08). Everything below is descriptive. Nothing here validates the lattice or the definition.**

---

> **Verifier corrections, applied by the orchestrator** (agent 1 could not write this file; 1,332 of 1,332 comparable arm-REAL rows reproduced field for field and the verdict survived the recompute):
> 1. **De novo member-loss list was mislabelled (substantive).** Under STRICT on the de novo arm the members destroyed are **NPIPB4, NPIPB12 and NPIPB13**; LOC124907844, LOC101929894 and LOC112268174 are not members at all. This strengthens the report's own argument against STRICT.
> 2. **The four certifying rows have h_join = 0.000000** — after deleting PKD1P6-NPIPP1 the remaining 25 members have *no* edge to any outside node at any level, so the "certificate" is (0, h_split] and contains every cut trivially.
> 3. **"Flat in delta" holds only for k ≤ 2**: at k = 3, six cells per arm change between δ = 4 and δ = 8 (member pairs lost 2 → 1).
> 4. **The GUARDED explanation is partly wrong on the facts**: LOC124907830 has 2 annotated junctions and LOC124907845 has 1, so they are not junction-less; the guard still re-admits them for other reasons.
> 5. **Addendum A2.1's guarded outside-count range is 216-238**, not 220-235.
> 6. **`bipF_soto` is not reproducible** under the verifier's reading of the Soto truth; `bipF_lit` and `bipF_iso` reproduce exactly (40/40 checked). Do not quote the Soto column.
> 7. **Arm IDEAL-READ's record table has ≥ 68 body rows with junction counts wrongly zeroed** (an error, not a file-matching question), and **arm IDEAL-DENOVO is not independently reproducible** (5,274 vs 5,781 records). Neither arm should carry argumentative weight until those are closed.
> 8. **NPIPB at L3 confirmed**: h_join 0.980559, h_split 0.983454, shipped cut 0.980000, gap **0.000559**; no form moves it strictly below 0.98.

## 0. Read first — adoptability verdict per form

The decision rule was frozen verbatim before any number (R2): a form is **ADOPTABLE iff on the real evidence some level certifies NPIP or a named subfamily with the level's shipped cut inside, and no present member is lost or disconnected.**

| # | form | parameters | certifies on REAL? | member lost? | member disconnected? | **verdict** |
|---|---|---|---|---|---|---|
| 1 | **RULE 1 (a) STRICT** | k=1,2,3 x delta=0,4,8 | **no**, 0 of 108 rows | no (0 isolated on REAL) | no (1 component) | **NOT ADOPTABLE** — fails (i) everywhere |
| 2 | **RULE 1 (b) GUARDED** | k=1,2,3 x delta=0,4,8 | **no**, 0 of 108 rows | no | no | **NOT ADOPTABLE** — and strictly dominated; h_join never leaves 1.000000 |
| 3 | **RULE 2 (a) REDUNDANT** | f=0.50,0.80,0.90,1.00 | **no**, 0 of 48 rows | **no, at any f** | no | **NOT ADOPTABLE** — safe but invisible alone |
| 3P | RULE 2 (a) permissive | same | no | no | no | identical drop sets to 3 at every f; the reading does not matter |
| 4 | **RULE 2 (b) BLUNT** | drop all readthroughs | **no**, 0 of 12 rows | **YES — PKD1P6-NPIPP1** | no | **NOT ADOPTABLE** — R2.1(ii), pre-committed at R8.3 |
| 5 | **RULE 1 STRICT k2d0 + RULE 2 (a)** | f=0.50,0.80,0.90,1.00 | **no** | no | no | **NOT ADOPTABLE** — closest miss: 0.0066 at L1 identity |
| 5P | same, post-hoc f=0.45, 0.40 | flagged A2.2 | no | no | no | **NOT ADOPTABLE**; and any such f is fitted to one record |
| 6 | **RULE 1 STRICT k2d0 + RULE 2 (b)** | — | **YES, all four DNA levels** | **YES — PKD1P6-NPIPP1** | REAL no / **IDEAL-READ yes** | **NOT ADOPTABLE** — certifies only by deleting a member, and does not replicate |
| 7 | RULE 1 GUARDED k2d0 + RULE 2 (a) f=0.50 | — | no | no | no | **NOT ADOPTABLE** |
| D | DIAG: STRICT k2d0 + drop non-member readthroughs | circular | no | no | no | **not a candidate** — uses the answer; reported for attribution only |

**3,420 certificate rows** were computed: 3 evidence arms x 33 rule forms x 4 levels x 9 sets. **Exactly 4 rows have `shipped_inside = yes`**, all of them form 6, arm REAL, set NPIP, one per level. Every one of those 4 rows also has `members_dropped = PKD1P6-NPIPP1`.

**The answer to the user's "lets add those then" is therefore: on this evidence, adding either rule does not buy a certificate, and the one cell that does buy one pays for it with a member.** The run is not empty, though — it produces two cut windows and one structural result that are worth more than the rules were.

---

## 1. What was measured, and what anchors it

**Three evidence arms (R3), all reusing already-mapped records.**

| arm | nodes | records | junction sets | is it "real" or "ideal"? |
|---|---|---|---|---|
| **REAL** | 47,965 primary RefSeq nodes (`family_cert/dna/nodes.tsv`) | `family_cert` PAF batches, 8,673 witness tuples | union of the node's annotated transcripts' introns (`cert2/junctions.guided.tsv`) | **real-data guided evidence** — nothing synthetic enters it |
| **IDEAL-READ** | same | same | junctions restricted to those carried by >= 3 primary MAPQ>=1 reads of the idealized synthetic BAM (`cert2/junctions.hybrid.tsv`) | ideal-substrate |
| **IDEAL-DENOVO** | 5,428 de novo shared-definition nodes | `npip_ideal/denovo`, 5,781 tuples | the node's own exon structure, >= 3-read filter | ideal-substrate |

Two disclosures. First, arm IDEAL-READ's junctions are read-derived, which is what RULE 1's own clause "never from the evidence" forbids; it is a sensitivity on the ideal substrate, not an instance of the rule as written (A1.4). The rule as written is instantiated by arms REAL and IDEAL-DENOVO. Second, the circularity statement C0 of `npip_ideal/DECLARATIONS.txt` applies unchanged to both ideal arms and is inherited, not re-argued.

**Anchoring (R3.5), recorded before the sweep in A1.3.** With no rule applied:

- arm REAL reproduces the committed `family_cert/cert/certificates.tsv` (variant `primary`) with **0 differences** in `h_join`, `h_split`, `n_present`, `comp_size_at_shipped`, `comp_outside_n` across all 9 sets x 4 axes (NPIP L1 1.000000/0.971665, L2 1.000000/0.989040, L3 1.000000/0.980559 — the committed numbers);
- arm IDEAL-READ reproduces `cert2/certificates.Bread.tsv` (test `none`) with **0 differences**;
- arm IDEAL-DENOVO reproduces `cert2/certificates.C.tsv` (test `none`) with **0 differences**.

**Monotonicity (R4.4), recorded before the sweep in A1.2.** The declared check was ">= 2,000 random pair-level instances". Fewer than 2,000 pairs exist (1,986 in REAL and IDEAL-READ, 1,312 in IDEAL-DENOVO), so the check was run **exhaustively over every pair of every arm** — 23,832 + 23,832 + 15,744 = 63,408 (pair, level, delta) cells. **0 violations** of monotone-in-k, **0** of monotone-in-delta, **0** of the chain `w(STRICT) <= w(GUARDED) <= w(unconjoined)`. R8.2 is not triggered. The argument is the one declared: the conjunct is a predicate on the record and the two nodes' fixed definitions, it never mentions the cut, so the edge set at cut c is still exactly `{w >= c}`, the graph is still a threshold graph, and D1 applies unchanged; the guard adds a pair-level, cut-independent predicate, which preserves that. The inherited caveat R4.5 stands: the shipped L2/L3 aggregation is not a single-record existential, that was already true of the unconjoined baseline, and every L2/L3 row here inherits it.

---

## 2. RULE 1 — the shared-junction conjunct

### 2.1 STRICT works, and is free, on annotated junctions — and still never certifies

Arm REAL, set NPIP, all 26 present members, delta makes no difference at all (k=1/2/3 x delta=0/4/8 give bit-identical certificates within each k):

| form | L1 identity h_join / h_split | L1 coverage | L2 | L3 | outside (L1/L2/L3) | member pairs lost | one component? |
|---|---|---|---|---|---|---|---|
| baseline | 1.000000 / 0.971665 | 1.000000 / 1.000000 | 1.000000 / 0.989040 | 1.000000 / 0.980559 | 96 / 64 / 57 | 0 | yes |
| STRICT k=1 | **0.999370** / 0.971665 | 1.000000 / 1.000000 | 1.000000 / 0.989040 | **0.999363** / 0.980559 | 25 / 25 / 23 | **0** | yes |
| **STRICT k=2** | **0.999363** / 0.971665 | 1.000000 / 1.000000 | 1.000000 / 0.989040 | **0.999363** / 0.980559 | **10 / 10 / 8** | **0** | yes |
| STRICT k=3 | 0.999363 / 0.971665 | 1.000000 / 1.000000 | 1.000000 / 0.989040 | 0.999363 / 0.980559 | 10 / 10 / 8 | **2** | yes |
| GUARDED, any k | **1.000000** / 0.971665 | 1.000000 / 1.000000 | **1.000000** / 0.989040 | **1.000000** / 0.980559 | 63-58 / 59-53 / 53-47 | 0 (2 at k=3) | yes |

`shipped_inside` is `no` in all 108 rule-1 rows of arm REAL, and in all 216 rule-1 rows of the two ideal arms.

**Who holds the boundary, and what STRICT does to them** (arm REAL, L1 identity):

| | baseline | STRICT k=1 | STRICT k=2 | GUARDED k=2 |
|---|---|---|---|---|
| CLN3 (253 bp, 1 block, biotype "other") | 1.000000 via NPIPB10P | **gone** | **gone** | **1.000000, still there** |
| EIF3CL | 1.000000 via NPIPB9 | **gone** | **gone** | gone |
| LOC100190986 (2,453 bp lncRNA, 1 block) | 1.000000 (L2/L3) | gone | gone | **1.000000, still there** |
| LOC124907830 / LOC124907845 (2 blocks each) | 1.000000 (L2) | gone | gone | **still there** |
| LOC128966632 (5,598 bp "SMG1-like") | 0.998496 | 0.998496 | gone | still there |
| PKD1P5-LOC105376752 | 0.999363 | 0.999363 | **0.999363 — holds the boundary** | 0.999363 |
| LOC131696449 | 0.999356 | 0.999356 | 0.999356 | 0.999356 |
| PKD1P4-NPIPA8 | — | 0.999209 | 0.999209 | 0.999209 |

This independently reproduces the `npip_ideal` result that the five named blockers all lose their NPIP edge at k >= 2, and that the boundary then passes to PKD1P readthroughs.

### 2.2 GUARDED is not weaker — it is refuted

GUARDED leaves NPIP's `h_join` at exactly **1.000000** at every level, every k, every delta, in all three arms. The reason is structural and reads directly off the table above: the guard exempts a pair when either copy has no junction in its fixed definition, and the records pinning NPIP at 1.000000 are exactly the junction-less ones — CLN3 has one exon block, LOC100190986 one, LOC124907830/845 two, and under GUARDED at L2 the top boundary is `CLN3:1.000000; LOC100190986:1.000000; MIR6511A1/A2/A3:1.000000(via PKD1P6-NPIPP1)`. **The guard re-admits precisely what the conjunct exists to remove.** It buys safety that STRICT did not need on this arm (STRICT isolates no member on REAL either) at the price of the entire effect.

### 2.3 Which members lose every edge under STRICT — the task's explicit question

| arm | members with 0 junctions in the fixed definition | members with **no edge at all** under STRICT |
|---|---|---|
| **REAL** | **0** | **none**, at any k or level |
| IDEAL-READ | 1: **NPIPB12** | **NPIPB12**, at every k and every level |
| IDEAL-DENOVO | 10 | **11**: LOC100420289, LOC100420311, NPIPB10P, LOC124907844, LOC101929894, **LOC112268174**, TBC1D3P1, TBC1D3P3, TBC1D3P4, TBC1D3P6, TBC1D3P7 |

LOC112268174 is the informative one: it *has* junctions in its de novo node, but none of them match, so it is isolated by the matching requirement rather than by junction-lessness. GUARDED isolates nobody in any arm.

### 2.4 The cost that decides it: STRICT does not survive the held-back substrate

| arm | NPIP under STRICT k=2 | member pairs lost | components |
|---|---|---|---|
| REAL | h_join 0.999363, h_split 0.971665 | **0** | **1** |
| IDEAL-READ | h_join 0.999363, **h_split -inf** | **24** (already 24 at k=1) | **2** |
| IDEAL-DENOVO | h_join 0.998065, **h_split -inf** | **34** (29 at k=1, 69 at k=3) | **5** |

TBC1D3 on IDEAL-DENOVO under STRICT k=2: h_split -inf, 8 components, 48 member pairs lost. The pattern is unambiguous: **the conjunct is free exactly where the junctions are annotation-derived, and destructive wherever they come from evidence.** On the ideal arms GUARDED keeps the family intact at zero cost and achieves nothing; STRICT achieves something and shatters the family. There is no setting of (k, delta) at which it does both.

### 2.5 Delta is dead

k = 1, 2, 3 differ. delta = 0, 4, 8 produce bit-identical certificates at every k, every level, every set, in all three arms. Junctions here either transfer exactly or not at all. If the conjunct is ever written down, **delta = 0** is the only defensible value and the parameter should not exist.

---

## 3. RULE 2 — the readthrough node rule

209 of the 47,965 primary nodes carry `readthrough = yes` (the GFF description contains "readthrough"; the same selector as `LATTICE_RULE_STRENGTHENERS` S5 after its A4). Exactly one is a member: **PKD1P6-NPIPP1** (NPIP, and a member of Dishuck's NPIPA).

### 3.1 What each form drops

| variant | f | nodes dropped genome-wide | members dropped | PKD1P-NPIP records that disappear |
|---|---|---|---|---|
| (a) redundant | 0.50 | **178** | **none** | PKD1P3-NPIPA1, PKD1P4-NPIPA8 |
| (a) redundant | 0.80 | **128** | **none** | PKD1P3-NPIPA1 |
| (a) redundant | 0.90 | **92** | **none** | PKD1P3-NPIPA1 |
| (a) redundant | 1.00 | **20** | **none** | none |
| (a) permissive | all f | 178 / 128 / 92 / 20 | none | identical to (a) at every f |
| **(b) blunt** | — | **209** | **PKD1P6-NPIPP1** | PKD1P3-NPIPA1, PKD1P4-NPIPA8, **PKD1P6-NPIPP1** |

The declared non-recursive reading (a readthrough may only be justified by a non-readthrough node) and the permissive reading give **identical drop sets at every f**, so that modelling choice carries no weight.

### 3.2 Why (a) spares the member, and (b) should not be used

Largest fraction of another node's exonic bases contained by each relevant readthrough:

| readthrough | best containment | of |
|---|---|---|
| PKD1P3-NPIPA1 | **0.9657** | PKD1P3 (also 0.9419 of NPIPA1) |
| LOC131696449 | **0.9268** | PKD1P1 (also 0.9159 of NPIPA6) |
| PKD1P4-NPIPA8 | **0.7550** | **NPIPA8 — a family member** |
| PKD1P5-LOC105376752 | **0.4897** | NPIPA9 — **misses f = 0.50 by 0.0103** |
| PDXDC2P-NPIPB14P | **0.1227** | LOC124907800 |
| **PKD1P6-NPIPP1** | **0.0038** | LOC100505915 |

PKD1P6-NPIPP1 contains essentially nothing, so variant (a) leaves it alone at any f above 0.004 — the member is safe by a wide margin, not by luck. Variant (b) deletes it by construction. **Variant (a) is the right shape of the rule; variant (b) is not.**

### 3.3 RULE 2 alone is invisible

NPIP's `h_join` is **1.000000 at every level under every f and both variants, including blunt.** The records holding the 1.000000 layer are not readthroughs (CLN3, EIF3CL, LOC100190986, LOC124907830/845, LOC128966632); the readthroughs sit at 0.999xxx underneath. The best RULE 2 does on its own is move NPIP's L1 outside count from 96 to 92. **RULE 2 only becomes measurable after RULE 1 has removed the layer above it** — which is exactly why the combination is the interesting run.

### 3.4 The rule fails where it is needed: PDXDC2P-NPIPB14P

Under RULE 1 STRICT k=2 + RULE 2 (a) at the post-hoc f = 0.40 (189 nodes dropped, all 26 members kept), NPIP's L2 boundary is held at exactly **1.000000 by PDXDC2P-NPIPB14P, via NPIPB1P**. That record **contains NPIPB14P — an NPIP member** — but NPIPB14P has no exon feature in RefSeq, so it is not a node, so the readthrough is not "redundant" under the rule, so it survives every f.

**A containment rule that can only test against nodes cannot remove a readthrough whose contained copy failed to become a node.** This is the same class of defect as the de novo arm's swallowing nodes, arriving from the opposite direction, and it is the single most transferable thing this run found about RULE 2's shape.

---

## 4. The combination, and the one cell that certifies

Best RULE 1 form under R6.1, recorded in A2.1 before the combination was run: no form satisfies (i), NPIP's best margin `h_split - h_join` is exactly **0.000000** for every form (attained on L1 coverage where both saturate at the axis ceiling 1.000000), so the tie-breaks decide — 0 member pairs lost (k <= 2), then NPIP's summed outside count (STRICT k=2: 38; STRICT k=1: 98; GUARDED: 220-235), then smaller k, then smaller delta, then STRICT. **BEST FORM = RULE 1 STRICT, k = 2, delta = 0.** This agrees with the independent `npip_ideal` selection.

**Arm REAL, NPIP:**

| form | L1 identity | L1 coverage | L2 | L3 | members | outside (L1/L2/L3) | pairs lost | inside? |
|---|---|---|---|---|---|---|---|---|
| STRICT k2 + (a) f=1.00 | 0.999363 / 0.971665 | 1.000000 / 1.000000 | 1.000000 / 0.989040 | 0.999363 / 0.980559 | 26 | 10 / 10 / 8 | 0 | no |
| STRICT k2 + (a) f=0.90 | 0.999363 / 0.971665 | 1.000000 / 1.000000 | 1.000000 / 0.989040 | 0.999363 / 0.980559 | 26 | 8 / 8 / 6 | 0 | no |
| STRICT k2 + (a) f=0.80 | 0.999363 / 0.971665 | 1.000000 / 1.000000 | 1.000000 / 0.989040 | 0.999363 / 0.980559 | 26 | 8 / 8 / 6 | 0 | no |
| STRICT k2 + (a) f=0.50 | 0.999363 / 0.971665 | 1.000000 / 1.000000 | 1.000000 / 0.989040 | 0.999363 / 0.980559 | 26 | 7 / 6 / 5 | 0 | no |
| *STRICT k2 + (a) f=0.45* (post-hoc) | **0.978223** / 0.971665 | 1.000000 / 1.000000 | 1.000000 / 0.989040 | 0.988804 / 0.980559 | 26 | 6 / 6 / 4 | 0 | no |
| *STRICT k2 + (a) f=0.40* (post-hoc) | **0.978223** / 0.971665 | 1.000000 / 1.000000 | 1.000000 / 0.989040 | 0.988804 / 0.980559 | 26 | 6 / 6 / 4 | 0 | no |
| **STRICT k2 + (b) blunt** | **0.000000** / 0.971665 | **0.000000** / 1.000000 | **0.000000** / 0.989040 | **0.000000** / 0.980559 | **25** | **0 / 0 / 0** | **10 / 10 / 7** | **YES x4** |

Form 6 gives NPIP a complete certificate at all four DNA levels — the component at each shipped cut is **exactly the 25 surviving members and nothing else**, with intervals (0, 0.971665], (0, 1.000000], (0, 0.989040], (0, 0.980559]. This is the first DNA-level certificate anywhere in this project. It is **not adoptable**, for two independent reasons, the first of which was pre-committed at R8.3:

1. It deletes **PKD1P6-NPIPP1**, a present NPIP member (and an NPIPA member — NPIPA goes 8 -> 7), and loses 10 member pairs at L1/L2 and 7 at L3. R2.1(ii).
2. **It does not replicate.** On arm IDEAL-READ the same form gives NPIP `h_split = -inf`, **2 components, 34 member pairs lost**. R2.1(iii).

### 4.1 Attribution: the last obstruction is a member

The diagnostic form D (drop every readthrough **except** family members; circular as a definition, labelled DIAG everywhere, A3.2) isolates the cause exactly. Arm REAL, STRICT k=2, all 26 members present, one component, **0 member pairs lost**:

| level | h_join | h_split | interval | outside | who is outside |
|---|---|---|---|---|---|
| L1 identity | 0.978223 | 0.971665 | empty **by 0.0066** | **4** | PKD1, PKD1P1, PKD1P2, PKD1P3 |
| L1 coverage | 1.000000 | 1.000000 | empty | 4 | same four |
| **L2 f_ex** | **0.746986** | **0.989040** | **(0.746986, 0.989040]** | 4 | same four |
| L3 w_98 | 0.988804 | 0.980559 | empty by 0.0082 | 4 | same four |

Top boundary, L1 identity: `PKD1:0.978223(via PKD1P6-NPIPP1); PKD1P2:0.971834(via PKD1P6-NPIPP1); PKD1P3:0.869490(via PKD1P6-NPIPP1); PKD1P1:0.868219(via PKD1P6-NPIPP1)`. **Every surviving boundary edge of NPIP runs through its own member PKD1P6-NPIPP1's PKD1 half.** The family's last link to the outside is not a neighbour that can be deleted; it is a member that is itself half of something else. That is why form 6 works and why form 6 is not a rule: the only way to cut the last edge is to cut a member.

This is the certificate-side statement of the same fact recorded in `project_leader_rule_breaks_nesting` and in `npip_ideal` §8 — a set S whose boundary object is not disjoint from S cannot be isolated by any predicate on pairs.

### 4.2 A real NPIP window at L2

The DIAG row also produces something new: **NPIP is an exact connected component of the f_ex filtration for any cut in (0.746986, 0.989040]** — width 0.242, with all 26 members. The shipped L2 cut is 0.30, **below** the window. This is the NPIP analogue of the TBC1D3 L2 window `LATTICE_RULE_STRENGTHENERS` called the most actionable thing its sweep found. It is reached only through a circular node rule, so it is a reading about where NPIP lives in the filtration, not a route to a certificate.

---

## 5. Sets, subfamilies, and the NPIPB near miss

### 5.1 The declared R7.2 question, answered

**Does any form move NPIPB's L3 h_join strictly below 0.98? No.** Under every RULE 2 form and every GUARDED form it is **1.000000**. Under STRICT at k >= 2 and every delta it is **exactly 0.980559** — the shipped 0.98 misses by **0.000559**, reproducing the `npip_ideal` near miss on the real arm.

And the boundary holder is now named: `NPIPA1:0.980559(via NPIPB2); PKD1P3-NPIPA1:0.980559(via NPIPB2); LOC131696449:0.980077; NPIPA6:0.980077; NPIPA2:0.979989`. **NPIPB's obstruction is no longer a foreign record — it is NPIPA**, i.e. the NPIPA/NPIPB split edge itself, at exactly the weight the committed report gives for that split. No edge rule can lower it without splitting NPIP.

### 5.2 Subfamilies: many non-empty intervals, none containing a shipped cut

Arm REAL, L3 w_98 (shipped cut 0.98), STRICT k=2 delta=0, 0 member pairs lost:

| set | baseline | STRICT k=2 | STRICT k=2 + (a) f=0.40 | miss vs 0.98 |
|---|---|---|---|---|
| NPIPA (8) | empty | empty | **(0.988804, 0.989838]** | +0.0088 |
| NPIPB (18) | empty | **(0.980559, 0.983454]** | (0.980559, 0.983454] | **+0.00056** |
| A6-9 (4) | empty | empty | **(0.995987, 0.996044]** | +0.0160 |
| B3-5 (4) | empty | **(0.997129, 0.998667]** | (0.997129, 0.998667] | +0.0171 |
| B6-9 (4) | empty | **(0.994729, 0.997680]** | (0.994729, 0.997680] | +0.0147 |
| B12/13 (3) | empty | **(0.997129, 0.999128]** | (0.997129, 0.999128] | +0.0171 |
| B15 (3) | (0.992435, 0.999517] | (0.992435, 0.999517] | (0.992435, 0.999517] | +0.0124 |

**The new positive.** Under **RULE 1 STRICT k=2 delta=0 alone** — 26/26 members, 0 member pairs lost, one component, no node deleted — **four of the five Dishuck Iso-Seq groups (B3-5, B6-9, B12/13, B15) are simultaneously exact connected components of the w_98 filtration, for any single cut in (0.997129, 0.997680]**. In the unconjoined baseline **no two of the five are ever simultaneously exact**. Only A6-9 is excluded (its upper end, 0.996044, falls below B3-5's lower end). The shipped 0.98 is not in the window, so this is a cut question, not a rule question, and it holds on arm REAL only — on IDEAL-READ B12/13 is internally disconnected.

### 5.3 TBC1D3

Under RULE 1 STRICT k=2 + RULE 2 at any f, TBC1D3's L1-identity h_join falls 1.000000 -> **0.961397** (NPEPPSP1 and TBC1D3P1-DHX40P1 removed; TBC1D29P via TBC1D3P5 and USP6 via TBC1D3P2 remain) against h_split 0.807623 — still empty. At L3 TBC1D3 has **0 nodes outside** under every form, and still does not certify, because its 12 nodes sit in **3 components** at the cut. Its problem was never over-merge at L3; it is internal fragmentation, and neither rule addresses that.

---

## 6. Bipartite F — reported as cost, and a live demonstration of the metric trap

Arm REAL, NPIP, against the Soto families (`truth_soto_families.tsv`), the Dishuck **subfamily** field and the Dishuck **Iso-Seq group** field. Following the correction flagged in `LATTICE_RULE_STRENGTHENERS` item (6), the field used is disclosed in advance (R7.1): `bipF_lit` is the **subfamily** field, `bipF_iso` the **isoseq_group** field; `project_level1` is not used.

| form | L1 F_soto | L1 F_lit | L1 F_iso | L3 F_soto | L3 F_lit | L3 F_iso |
|---|---|---|---|---|---|---|
| baseline | 0.228188 | 0.243243 | 0.054054 | 0.309091 | 0.330275 | 0.073394 |
| STRICT k=1 | 0.435897 | 0.467532 | 0.103896 | 0.447368 | 0.480000 | 0.106667 |
| **STRICT k=2** | 0.539683 | 0.580645 | 0.129032 | 0.557377 | 0.600000 | 0.133333 |
| GUARDED k=2 | 0.295652 | 0.315789 | 0.070175 | 0.326923 | 0.349515 | 0.077670 |
| (a) f=0.50 alone | 0.234483 | 0.250000 | 0.055556 | 0.330097 | 0.352941 | 0.078431 |
| (b) blunt alone | 0.244275 | 0.266667 | 0.059259 | 0.340426 | 0.367347 | 0.081633 |
| STRICT k2 + (a) f=0.50 | 0.566667 | 0.610169 | 0.135593 | 0.586207 | 0.631579 | 0.140351 |
| STRICT k2 + (a) f=0.40 | 0.576271 | 0.620690 | 0.137931 | 0.596491 | 0.642857 | 0.142857 |
| **STRICT k2 + (b) blunt** | **0.695652** | **0.720000** | **0.160000** | **0.695652** | **0.720000** | **0.160000** |

TBC1D3 F_soto: 0.487805 baseline -> 0.689655 under the combination at L1; **1.000000 at L3 under every form including the baseline**, so L3 F carries no signal for TBC1D3 at all.

**F rises in every row and never falls** — +0.467 on F_soto for NPIP at L1 — exactly as `LATTICE_RULE_STRENGTHENERS` reading 3 warned. R2.1(iv) and R8.4 were declared in advance precisely so this could not be read as support, and it is not: **the only form that certifies is form 6, and it is disqualified for a reason F cannot see.**

The Iso-Seq column makes the trap concrete rather than abstract. Under STRICT k=2, `bipF_iso` reads **0.129** — the Iso-Seq groups look hopeless as components at the shipped cut. At that same setting, **four of those five groups are exact connected components at a cut of 0.9975** (§5.2). F at a fixed cut and exactness in the filtration are measuring different things, and here they point in opposite directions.

---

## 7. Readings

1. **Neither rule is adoptable, and the decision rule caught the one case that would have looked like success.** Form 6 produces a clean, complete DNA certificate for NPIP — and does it by deleting a member. R8.3 was written before any number precisely to stop that being reframed later. It also fails to replicate on the ideal arm, so even the pre-commitment was not load-bearing.

2. **NPIP's final obstruction is a member, not a neighbour.** With every non-member readthrough removed and t_J applied, NPIP's whole boundary is four nodes — PKD1, PKD1P1, PKD1P2, PKD1P3 — every one reached through PKD1P6-NPIPP1. No monotone edge conjunct and no node-deletion rule can reach that without deleting a member. This is the seventh failed strengthener (five in `LATTICE_RULE_STRENGTHENERS`, t_J/t_E in `npip_ideal`, these two here) and the first one that names the obstruction as internal.

3. **The guarded variant is refuted, not merely weak.** Exempting junction-less copies exempts exactly the blockers. Any future guard of this shape will do the same.

4. **The redundancy test inherits the node set's blind spots.** PDXDC2P-NPIPB14P contains an NPIP member and cannot be removed, because that member (NPIPB14P) is exon-less and therefore not a node. The node question and the rule question are not separable.

5. **STRICT is free on annotated junctions and destructive on real ones.** 0 / 24 / 34 member pairs lost across REAL / IDEAL-READ / IDEAL-DENOVO at k=2. A rule validated on arm REAL alone would have looked costless. This is `feedback_hold_a_substrate_back` firing on a rule that was about to be adopted.

6. **delta is a parameter with no content here** — 0, 4 and 8 give bit-identical certificates everywhere.

7. **Two cut windows are the run's actual yield**, and both are reporting facts, not rule changes: NPIP exact at L2 in (0.746986, 0.989040] (via a circular node rule, so a reading only); and four of five Iso-Seq groups simultaneously exact at L3 in (0.997129, 0.997680] under STRICT k=2 with zero cost — which the baseline cannot do for any two of them.

8. **NPIPB misses by 0.000559 and the miss is now explained.** Its boundary at L3 under STRICT is NPIPA1 at 0.980559 via NPIPB2 — the NPIPA/NPIPB split. Chasing that 0.00056 means moving the L3 cut to about 0.981, which would be fitting the cut to one subfamily of one development family.

---

## 8. Provenance

- **Declarations** `/mnt/linuxdisk/home/juanfraitu/rule_adopt/DECLARATIONS.txt`, frozen 14:47:18 -07:00 (md5 `06897978bdcbd5601a122eccbdad5d31`, recorded in `DECLARATIONS.md5`/`.time`) before the first number; final md5 `089cde5f12a3b675bd9606294ca95424`. R0.2 discloses the two annotation counts (209 readthrough nodes; PKD1P6-NPIPP1 the one member readthrough) that were observed while surveying the node table before the file was written. A1 records anchoring, monotonicity and the junction-file choice before the sweep; A2 records the best-form selection and flags the post-hoc f = 0.45 / 0.40 probe; A3 closes the record and flags the DIAG form as circular.
- **Evidence, reused unchanged, no new alignment:** `family_cert/dna/nodes.tsv` (47,965 primary nodes), `family_cert/dna/witnesses.tsv`, `family_cert/dna/batches/*.paf`, `family_cert/cert/{certificates,components,dna_pairs.extended,dna_query_nodes}.tsv`; `npip_ideal/cert2/records.{B,Bread,C}.pkl`, `cert2/junctions.{guided,hybrid}.tsv`, `npip_ideal/nodes/*`, `npip_ideal/denovo/*`. `family_cert/cert/dna_edges.py` (its `Agg`) and `family_cert/dna/dna_cert.py` and `npip_ideal/k1_denovo.py` were imported, not reimplemented. Truth sets: `layer_order/npip_tbc1d3/light/members.corrected.tsv`, `truth_soto_families.tsv`, `docs/lit_subclusters_npip_dishuck_check.tsv`. Bipartite F follows `lattice_rules/engine.py:bip_f`.
- **Outputs**, all under `/mnt/linuxdisk/home/juanfraitu/rule_adopt/`: `sweep.tsv` (3,420 rows, md5 `a3c8f2da83ccfa67c46ddc7fb187e7ed`; 1,368 REAL + 1,368 IDEAL-READ + 684 IDEAL-DENOVO), `sweep_phase1.tsv`, `sweep_combo.tsv`, `sweep_diag.tsv`, `rule2_drops.tsv`, `rule2_why.tsv`; scripts `engine.py` (md5 `8ca20d6d2cfa21eab0feb94232b47dcf`), `anchor.py`, `anchor2.py`, `mono.py`, `sweep.py`, `combo.py`, `diag.py`.
- **Compute**: all foreground, one at a time; total wall time about 1 minute; no background job, no `pkill`. `TMPDIR=/mnt/linuxdisk/home/juanfraitu/rule_adopt/tmp` throughout.

## 9. Caveats

1. **Descriptive.** Two development families, one haplotype, one annotation. Adding evidence can only raise h_join under the monotone tests, so a certificate can be lost and never gained by more data — which cuts against form 6's certificate as much as against anything else here.
2. **Component sizes are lower bounds.** Evidence from a node no member query ever touched is assumed absent (`family_cert` D3). Every `comp_size` and `comp_outside_n` in this run inherits that.
3. **L2 and L3 rows inherit the aggregation caveat (R4.5).** The shipped L2/L3 gate may come from one record while the weight comes from another, so the "witness" reading of an L2/L3 edge is not sound. This was already true of the unconjoined baseline and of the committed `FAMILY_CERTIFICATES_NPIP_TBC1D3.md` numbers; nothing here fixes it, and §4.2's L2 window rests on it.
4. **L1a/L1b rows depend on the greedy gene-body chains**, which are not monotone (H3). No tx-only sensitivity was run here.
5. **Arm IDEAL-READ is not an instance of RULE 1 as written** — its junction sets are read-derived, which the rule's own text forbids (A1.4). It is a sensitivity.
6. **RULE 2 was not run on the de novo arm** (R5.6), so the combination has only two arms.
7. **f = 0.45 and 0.40 are post-hoc and fitted to one record** (A2.2). No conclusion rests on f < 0.50.
8. **The DIAG form is circular** and is reported only to attribute form 6's certificate; it is not a candidate definition (A3.2).
9. **bipF_lit uses the Dishuck subfamily field**, disclosed in advance to correct the undisclosed deviation flagged in `LATTICE_RULE_STRENGTHENERS` item (6). Values are not comparable with that document's pre-correction numbers.
10. **The four-Iso-Seq-group window (§5.2) is arm-REAL only.** It does not survive on the read-supported arm, where B12/13 is internally disconnected. Treat it as a lead, not a result.
## Verification (independent recompute)

Agent 2 (independent verifier), 2026-09-17, 15:00-15:25 -07:00. No script under `rule_adopt/` outside `rule_adopt/verify/` was read or imported; the builder's `engine.py`, `anchor*.py`, `mono.py`, `sweep.py`, `combo.py` and `diag.py` were never opened. The verifier's code is in `/mnt/linuxdisk/home/juanfraitu/rule_adopt/verify/` (`w2_ann.py` annotation re-derivation, `w4_engine.py` rule + D1 engine, `w5_rule2.py` readthrough containment, `w6_sweep.py`/`w8_c.py` sweeps, `w7_diff.py` diff, `w9_checks.py` claim checks, `wA_bipf.py` bipartite F), plus re-executions of the two previously verified record derivations (`w1b_real_rows.py`, `w1c_bread_rows.py`, `w1d_c_rows.py`, from `npip_ideal/verify/v6_rules.py` and `v9_bread.py`/`v8_setC.py`). Outputs: `w_real.tsv` (1,368 rows), `w_bread.tsv` (1,368), `w_c.tsv` (684).

1. **Order.** `DECLARATIONS.md5` was written 14:47:18, before the first result file (`engine.py` 14:49:14, `rule2_drops.tsv` 14:50:41, `sweep_phase1.tsv` 14:50:45). The current `DECLARATIONS.txt` (mtime 14:55:48) is longer than the pinned version, but its first 18,029 bytes - everything through R8.5, i.e. the whole pre-run body - hash to exactly the recorded `06897978bdcbd5601a122eccbdad5d31`, so the append-only claim is confirmed: A1, A2 and A3 are additions and nothing above them was edited. The individual write times of A1/A2/A3 are not independently checkable; only the pre-addendum body is hash-pinned. Addenda A1.1 (exhaustive instead of 2,000 pairs), A1.4 (read-derived junctions on arm IDEAL-READ are a sensitivity, not the rule as written), A2.2 (f = 0.45/0.40 post-hoc) and A3.2 (the DIAG form is circular) are properly dated and disclosed.

2. **Annotation re-derived.** Nodes, junction sets and the readthrough flag were rebuilt from `chm13v2.0_RefSeq_full.gff.gz` with the verifier's own parser. Junction sets (union over all annotated transcripts' introns, per gene/pseudogene): **0 differences on all 47,965 rows** of `npip_ideal/cert2/junctions.guided.tsv`. Readthrough flag (`'readthrough'` in the GFF description): **0 mismatches** against `family_cert/dna/nodes.tsv`, 209 of 47,965 primary nodes, and **exactly one member node carries it - gene-PKD1P6-NPIPP1** (R0.2 confirmed). All 38 present member nodes carry at least one annotated intron, so no member is junction-less on arm REAL.

3. **Records.** The arm-REAL record table was regenerated from the `family_cert` PAFs with the corrected body-junction offsets: 8,673 records, **byte-identical** to `npip_ideal/verify/v6b_rows.pkl`, which reproduces `cert2/certificates.B.tsv` with 0 differences. (The saved `v6_rules.py` still carries the uncorrected body path - it maps donor and acceptor to the same query offset - and gives 564 differing certificate rows; that script's output must not be used.) Arm IDEAL-READ regenerated byte-identical to `v9_rows.pkl`, arm IDEAL-DENOVO byte-identical to `v8_rows.pkl`.

4. **RULE 2 recomputed from the exon unions.** Genome-wide drops **178 / 128 / 92 / 20** at f = 0.50 / 0.80 / 0.90 / 1.00 and 186 / 189 at the post-hoc f = 0.45 / 0.40; the strict reading (containers must be non-readthrough) and the permissive one give **identical drop sets at every f**; blunt drops 209 including the single member PKD1P6-NPIPP1. Every containment fraction in the report is confirmed: PKD1P6-NPIPP1's largest containment of any node is **0.003849** (10 of LOC100505915's 2,598 bp), so variant (a) never drops it; PKD1P3-NPIPA1 contains 0.9657 of PKD1P3 and 0.9419 of NPIPA1; LOC131696449 0.9268 of PKD1P1 and 0.9159 of NPIPA6; PKD1P4-NPIPA8 0.7550 of NPIPA8; PKD1P5-LOC105376752 0.489718 of NPIPA9; PDXDC2P-NPIPB14P at most 0.1227 (of LOC124907800), and **NPIPB14P has no exon feature and is therefore not a node**, so the readthrough that contains an NPIP member survives every f - the sharpest structural finding is confirmed exactly.

5. **Certificates.** All 1,332 comparable arm-REAL rows were recomputed with an independent D1 implementation: **0 differences** in `n_present`, `members_dropped`, `n_nodes_dropped`, `h_join`, `h_split`, `interval`, `exact_nonempty`, `shipped_inside`, `comp_size_at_cut`, `comp_outside_n`, `n_member_pairs_at_cut` and `member_pairs_lost_vs_base`, across the baseline, both RULE 1 variants over k in {1,2,3} x delta in {0,4,8}, both RULE 2 variants over f, the combinations and the DIAG rows. The decisive count is confirmed: **exactly 4 rows have the shipped cut inside a non-empty interval**, all four NPIP under STRICT k=2 delta=0 + RULE 2 blunt on arm REAL, all four with PKD1P6-NPIPP1 dropped - and all four with **h_join = 0.000000**, i.e. the set has no boundary edge whatsoever once that member and the other readthroughs are deleted. On arm IDEAL-READ the same recompute agrees on h_join/h_split/exact/inside/member-pairs except 4 NPIPB h_join cells (see corrections); on arm IDEAL-DENOVO the builder reproduces `certificates.C.tsv` exactly but an independent record derivation does not, so that arm is marked.

6. **Costs, confirmed number by number (arm REAL).** STRICT k=2 delta=0: NPIP L1-identity h_join 1.000000 -> **0.999363**, h_split unchanged at 0.971665, outside count summed over the four levels **96 -> 10 at L1** (313 -> 38 over all levels), **0 member pairs lost** at k <= 2 and **2** at k = 3, all 26 present members in one component. The post-STRICT boundary is exactly as reported: PKD1P5-LOC105376752 0.999363, LOC131696449 0.999356, PKD1P4-NPIPA8 0.999209, PKD1P3-NPIPA1 0.995197, PKD1 0.978223 (via PKD1P6-NPIPP1). RULE 2 alone leaves NPIP's h_join at 1.000000 at every level and every f and moves the outside count 96 -> 92 at best. Blunt costs NPIP 26 -> 25 present, NPIPA 8 -> 7, and 10/10/10/7 member pairs at L1-id/L1-cov/L2/L3. The f = 0.40 probe: 189 nodes dropped, all 26 members kept, L1-identity h_join 0.978223 against h_split 0.971665 (empty by 0.006558), L3 0.988804 against 0.980559, and **PDXDC2P-NPIPB14P alone holds the L2 boundary at exactly 1.000000 via NPIPB1P**.

7. **Monotonicity, checked with the verifier's own instances.** Exhaustive over every weight-bearing pair of every arm (1,469 / 1,469 / 942 pairs x 4 levels x 3 deltas = 17,628 / 17,628 / 11,304 cells): **0 violations** of monotone-in-k, **0** of monotone-in-delta, and **0** of the chain w(STRICT) <= w(GUARDED) <= w(unconjoined), on all three arms. R8.2 is not triggered. The STRICT/GUARDED distinction behaves as declared on junction-less copies: on arm REAL no present member has |J| = 0, so STRICT isolates none; on arm IDEAL-READ exactly one does - **NPIPB12, |J| = 0 under the >= 3-read hybrid set** - and STRICT removes all its edges; on arm IDEAL-DENOVO ten member nodes have 0 read-supported introns and an eleventh has introns that never match.

8. **The two real positives.** The GUARDED variant's futility is confirmed structurally: NPIP's h_join stays at exactly 1.000000 at every level, every k and every delta, held at L1-identity by **gene-CLN3 alone** (|J| = 0, 253 bp, one exon block, biotype 'other') - though at L1-coverage and L2 junction-rich readthroughs also sit at 1.000000, so the guard is not the only reason (see corrections). The new positive reproduces exactly: under STRICT k=2 delta=0 on arm REAL the four Iso-Seq groups B3-5 (0.997129, 0.998667], B6-9 (0.994729, 0.997680], B12/13 (0.997129, 0.999128] and B15 (0.992435, 0.999517] are simultaneously exact for any cut in **(0.997129, 0.997680]**, a window that does not contain the shipped 0.98; at the unconjoined baseline the same intersection is empty (B3-5 and B6-9 both have h_join = 1.000000).

9. **Scope.** This verification covers the arithmetic, the annotation and the evidence tables, not the declarations' reach: D3's lower-bound component assumption and R8.5 stand unchanged, `bipF_soto` is unreproduced, and arm IDEAL-DENOVO carries an inherited record-set discrepancy. Nothing here changes the run's verdict - neither rule is adoptable under the pre-registered decision rule, and the only form that certifies does so by deleting an NPIP member and leaving the set with no boundary at all.