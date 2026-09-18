# A world with no readthroughs: it does not rescue the family definition, and on real data it cannot be built without deleting a member

Agent 1 of 2. Declarations frozen **before any number of this run**: `/mnt/linuxdisk/home/juanfraitu/no_readthrough/DECLARATIONS.txt`, written 2026-09-18T12:39:25-07:00, md5 at that moment `a6f544c392bede6d5797d4b73147940d` (recorded in `DECLARATIONS.md5`/`DECLARATIONS.time`). Three dated addenda are appended; A1 corrects a transcription error before the first certificate existed, A2 discloses a discarded run and two scope reductions, A3 records what was not run. Final md5 `5464335ea9d50a9113e3a3a16893e200`. All outputs under `/mnt/linuxdisk/home/juanfraitu/no_readthrough/` (`out/`, `nodes/`, `logs/`, `code/`), `TMPDIR` under `no_readthrough/tmp` throughout. Nothing under `/mnt/c/Users/jfris/Desktop/Rustle` was modified (no file under `src/` has a 2026-09-18 mtime), nothing was committed, no subagent was spawned, no MAPQ gate was lowered, no index was rebuilt. One heavy job at a time, foreground, killed by PID once.

---

> **Verifier corrections, applied by the orchestrator** (verifier ok = true; the answer stands, several supporting numbers do not):
> 1. **"No read spans two annotated genes by construction" is FALSE.** 783 of 38,594 surviving transcripts (2.03%), from 273 loci, still overlap ≥ 50 bp of another gene's exons — annotated genes overlap each other independently of readthroughs, so a world without readthrough *structure* cannot be built by deleting readthrough records.
> 2. **The "new boundary class" sentence is level-dependent** and is not uniformly the small lncRNA/"other" class at L1.
> 3. **The node-key comparison mixed two exon conventions** (control table with reconstructed merges vs treatment from the query FASTA), which inflates the control-only / treatment-only counts.
> 4. **Two headline deltas shrink under matched construction**: NPIP mean copy coverage is 0.541 → 0.518 (not 0.591 → 0.518) and mean node purity 0.672 → 0.768.
> 5. **The per-copy sentence is wrong on four of six copies**: LOC128966608 is 0.034 in both arms, and the claimed collapses/perfections do not all hold under matched construction.
> 6. **Arm A's family numbers are quoted from a column a previous verifier retired** (`bipF_soto`, not reproducible) — they carry no weight.
> 7. Minor: the transcript denominator is 38,627, not 38,669.

## 0. The answer, first

**No.** A counterfactual world without readthroughs does not HELP the family definition under the reading rule frozen at R2, and the reasons are more interesting than the verdict.

| | condition | result |
|---|---|---|
| (i) | NPIP's certificate interval becomes non-empty **and** contains the shipped cut | **only in the arm that deletes the member PKD1P6-NPIPP1** — R2.4 fires, recorded as NOT A HELP |
| (ii) | failing that, h_join falls strictly below NPIP's internal bottleneck at some level | **yes, at L2 only**, and only in the member-preserving variant with t_J: (0.746986, 0.989040]. The shipped L2 cut 0.30 is **below** the window |
| (iii) | the family-level numbers do not fall | bipartite F **rises everywhere** — which R5.6 named in advance as the mechanical effect of deleting nodes; on the simulated arm the node panel is **two-sided**, not an improvement |

And the required statement of where the boundary goes: **removing readthroughs hands NPIP's boundary to the small junction-poor record class** — CLN3 (253 bp, one exon block, biotype "other"), EIF3CL, LOC100190986 (2,453 bp lncRNA, one block), LOC124907830 (768 bp, two blocks), LOC124907845 (1,106 bp, two blocks), LOC128966632 ("SMG1-like", 5,598 bp). Every one of them sits at exactly 1.000000 and **none of them is a readthrough.** When the shared-junction conjunct removes that class too, the boundary passes to **PKD1, PKD1P1, PKD1P2 and PKD1P3 — all four reached through NPIP's own member PKD1P6-NPIPP1's PKD1 half.** The family's last link to the outside is not a neighbour that a counterfactual can delete; it is a member that is itself half of a readthrough.

---

## 1. What "readthrough" means here, and how many are members

Frozen at R1.1 before any number: a RefSeq record is a readthrough iff its gene/pseudogene line's `description` contains "readthrough" — exactly the selector of `bench/JUNCTION_AND_READTHROUGH_RULES.md` §3.

I re-derived it independently from `chm13v2.0_RefSeq_full.gff.gz` rather than trusting the flag: **209 of the 58,563 gene/pseudogene records**, agreeing with `family_cert/dna/nodes.tsv`'s `readthrough` column on **209/209**, with **0 node-only and 0 GFF-only** discrepancies.

**Exactly one of the 209 is an NPIP or TBC1D3 member: `gene-PKD1P6-NPIPP1`** (NPIP; Dishuck subfamily NPIPA). TBC1D3 has none. The other readthroughs that matter to these families are neighbours, not members: PKD1P3-NPIPA1, PKD1P4-NPIPA8, PKD1P5-LOC105376752, LOC131696449 (described "PKD1P1-NPIPA5L readthrough"), PDXDC2P-NPIPB14P, and TBC1D3P1-DHX40P1.

The shipped Rust catalog has its **own, different** readthrough filter (single-exon transcripts engulfing ≥ 5 distinct junctions). It is not R1.1, it was not modified, and it stayed enabled in both arm-B runs; where it fires is reported separately (§3.1) so the two are never confused.

---

## 2. ARM A — real guided evidence, readthrough nodes deleted

Node set: the 47,965 primary RefSeq nodes of `family_cert/dna/nodes.tsv`. Machinery reused, not reimplemented (A1.2): `rule_adopt/engine.py`, whose arm-REAL rows were recomputed field-for-field by that run's verifier with 0 differences. **Anchoring: my control reproduces the committed `family_cert/cert/certificates.tsv` (variant primary) on all 36 comparable (set, level, axis) cells with 0 differences in h_join/h_split.**

### 2.1 NPIP, all four DNA levels

| arm | level | h_join | h_split | interval | shipped inside | members | member dropped | outside | pairs lost |
|---|---|---|---|---|---|---|---|---|---|
| CONTROL | L1-id | 1.000000 | 0.971665 | empty | no | 26 | — | 96 | 0 |
| CONTROL | L1-cov | 1.000000 | 1.000000 | empty | no | 26 | — | 96 | 0 |
| CONTROL | L2 | 1.000000 | 0.989040 | empty | no | 26 | — | 64 | 0 |
| CONTROL | L3 | 1.000000 | 0.980559 | empty | no | 26 | — | 57 | 0 |
| **NO_RT** | L1-id | **1.000000** | 0.971665 | empty | no | **25** | **PKD1P6-NPIPP1** | 85 | 10 |
| **NO_RT** | L1-cov | **1.000000** | 1.000000 | empty | no | 25 | PKD1P6-NPIPP1 | 85 | 10 |
| **NO_RT** | L2 | **1.000000** | 0.989040 | empty | no | 25 | PKD1P6-NPIPP1 | 51 | 10 |
| **NO_RT** | L3 | **1.000000** | 0.980559 | empty | no | 25 | PKD1P6-NPIPP1 | 48 | 7 |
| CONTROL + t_J | L1-id | 0.999363 | 0.971665 | empty | no | 26 | — | 10 | 0 |
| CONTROL + t_J | L3 | 0.999363 | 0.980559 | empty | no | 26 | — | 8 | 0 |
| **NO_RT + t_J** | L1-id | **0.000000** | 0.971665 | (0.000000, 0.971665] | **yes** | **25** | **PKD1P6-NPIPP1** | 0 | 10 |
| **NO_RT + t_J** | L1-cov | 0.000000 | 1.000000 | (0.000000, 1.000000] | yes | 25 | PKD1P6-NPIPP1 | 0 | 10 |
| **NO_RT + t_J** | L2 | 0.000000 | 0.989040 | (0.000000, 0.989040] | yes | 25 | PKD1P6-NPIPP1 | 0 | 10 |
| **NO_RT + t_J** | L3 | 0.000000 | 0.980559 | (0.000000, 0.980559] | yes | 25 | PKD1P6-NPIPP1 | 7 pairs lost | 0 outside |
| **NO_RT keep member** | L1-id | 1.000000 | 0.971665 | empty | no | **26** | — | 85 | **0** |
| **NO_RT keep member + t_J** | L1-id | **0.978223** | 0.971665 | empty by 0.006558 | no | 26 | — | **4** | 0 |
| **NO_RT keep member + t_J** | L1-cov | 1.000000 | 1.000000 | empty | no | 26 | — | 4 | 0 |
| **NO_RT keep member + t_J** | **L2** | **0.746986** | **0.989040** | **(0.746986, 0.989040]** | **no** (cut 0.30 below window) | 26 | — | 4 | 0 |
| **NO_RT keep member + t_J** | L3 | 0.988804 | 0.980559 | empty by 0.008245 | no | 26 | — | 4 | 0 |

Three readings.

1. **Removing readthroughs alone is invisible.** NPIP's `h_join` does not move off 1.000000 at any level. The only thing that changes is the outside count (96 → 85 at L1, 64 → 51 at L2, 57 → 48 at L3) — and it costs a member and 10 member pairs to buy that.
2. **The certificate that appears is vacuous and illegal.** NO_RT + t_J certifies at all four levels, but `h_join = 0.000000` means the 25 survivors have **no boundary edge at all** after the member was deleted; the interval contains every cut trivially. R2.4 was pre-committed for exactly this case.
3. **The member-preserving variant is the only real positive, and it is one level wide.** Drop the 208 non-member readthroughs, keep PKD1P6-NPIPP1, apply t_J: all 26 members present, one component, 0 member pairs lost, and L2 becomes an exact window of width 0.242. But the shipped cut is at 0.30, below the window, so condition (i) still fails; L1-id misses by 0.0066, L3 by 0.0082.

### 2.2 What the boundary becomes, named

| configuration | NPIP's top boundary at L1 identity |
|---|---|
| control | `CLN3:1.000000(via NPIPB10P); EIF3CL:1.000000(via NPIPB9); PKD1P2:0.999370(via NPIPA8); PKD1P5-LOC105376752:0.999363(via NPIPA6); LOC131696449:0.999356(via NPIPA9)` |
| **NO_RT** | `CLN3:1.000000(via NPIPB10P); EIF3CL:1.000000(via NPIPB9); PKD1P2:0.999370(via NPIPA8); PKD1P1:0.999356(via NPIPA9); LOC128966632:0.998496(via NPIPB3)` |
| control + t_J | `PKD1P5-LOC105376752:0.999363; LOC131696449:0.999356; PKD1P4-NPIPA8:0.999209; PKD1P3-NPIPA1:0.995197; PKD1:0.978223(via PKD1P6-NPIPP1)` |
| **NO_RT keep member + t_J** | `PKD1:0.978223(via PKD1P6-NPIPP1); PKD1P2:0.971834(via PKD1P6-NPIPP1); PKD1P3:0.869490(via PKD1P6-NPIPP1); PKD1P1:0.868219(via PKD1P6-NPIPP1)` |

The class the boundary is handed to is precise: **small, junction-poor lncRNA and fragment records**, none of them a readthrough. CLN3 at the boundary is the 253 bp single-block "other" record, not the 15-exon protein-coding CLN3. And once the conjunct removes that class, what is left is the **PKD1 gene family reached through NPIP's own readthrough member**.

### 2.3 TBC1D3 does not move at all

| level | control h_join / h_split | NO_RT h_join / h_split | outside |
|---|---|---|---|
| L1-id | 1.000000 / 0.840488 | **1.000000** / 0.840488 | 19 → 14 |
| L1-cov | 1.000000 / 0.841369 | **1.000000** / 0.841369 | 19 → 14 |
| L2 | 1.000000 / 0.832037 | **1.000000** / 0.832037 | 10 → 10 |
| L3 | 0.961397 / 0.867946 | **0.961397** / 0.867946 | 0 → 0 |

Bit-identical `h_join` at every level. TBC1D3's obstruction is NPEPPSP1 (via TBC1D3G), TBC1D29P (via TBC1D3P5) and USP6 (via TBC1D3P2), and at L3 it has **zero nodes outside** and still does not certify, because its 12 nodes sit in several components at the cut. **TBC1D3's problem was never readthroughs and never over-merge; it is internal fragmentation, and the counterfactual does not touch it.**

### 2.4 Family-level numbers

Disclosed substitution (R2.1(iii)): arm A builds no nodes, so the family-level numbers reported are bipartite F of the certificate components against the three truth partitions.

| configuration | NPIP L1 F_soto / F_lit / F_iso | NPIP L3 F_soto | TBC1D3 L1 F_soto |
|---|---|---|---|
| control | 0.228188 / 0.243243 / 0.054054 | 0.309091 | 0.487805 |
| NO_RT | 0.244275 / 0.266667 / 0.059259 | 0.340426 | 0.555556 |
| control + t_J | 0.539683 / 0.580645 / 0.129032 | 0.557377 | 0.689655 |
| NO_RT + t_J | 0.695652 / 0.720000 / 0.160000 | 0.695652 | 0.689655 |
| NO_RT keep member + t_J | 0.596491 / 0.642857 / 0.142857 | 0.596491 | 0.689655 |

**F rises in every row and never falls.** R5.6 was written before any number precisely so this could not be read as support: deleting nodes shrinks the precision denominator. The one configuration with the highest F is the one disqualified for deleting a member.

---

## 3. ARM B — the simulated world rebuilt with no readthroughs

### 3.1 Building the world

40 of the 5,542 ideal loci are readthroughs under R1.1 (91 of 38,669 transcript rows, 2,730 reads, **11,100 of 2,789,740 alignment records = 0.40%**). Exactly one is a member (PKD1P6-NPIPP1).

The declared compute equivalence (R3.3, stated before use): error-free per-transcript reads aligned independently against the same index mean the treatment BAM is exactly the control BAM minus every record whose read name begins with a readthrough `gene_id`. **Verification passed in full, not on a sample:** 11,100 records removed, **0 residual readthrough-locus records** in the filtered BAM (whole-file scan), header byte-identical excluding `@PG`, 5,502 distinct source genes remaining. **No read in the treatment BAM spans two annotated genes by construction.**

| | control | treatment (no readthroughs) |
|---|---|---|
| shipped catalog's own RT filter | dropped **2** single-exon transcripts → 37,488 | dropped **1** → 37,401 |
| — the one it stops firing on | `chr16:75785805-75805631` (the PDXDC2P-NPIPB14P locus) | locus no longer exists |
| skeletons → reps | 38,052 → **5,375** over 24 contigs | 37,961 → **5,389** over 24 contigs |
| shared-definition nodes | 5,375 reps → 5,360 gene-level + 68 read-locus = **5,428** | 5,389 reps → 5,374 gene-level + 70 read-locus = **5,444** |

**Removing readthroughs produces more reps and more nodes, not fewer.** 5,375 node keys are identical in both sets; 53 are control-only and 69 treatment-only — so the change is small but not purely local: the homology/rep stage shifts globally when 40 loci disappear.

### 3.2 Per-copy node recovery, matched 5,502-locus universe

| | control | **treatment** | cheap control (control nodes minus the 38 readthrough-derived nodes) |
|---|---|---|---|
| nodes | 5,428 | **5,444** | 5,390 |
| NPIP copies with a node | 26/26 | **26/26** | **22/26** |
| NPIP full-length (≥0.95 cover, ≤0.05 outside) | 7 | **9** | 7 |
| NPIP full at the looser ≥0.80 | 11 | **10** | 9 |
| NPIP mean fraction of the copy covered | 0.5912 | **0.5175** ↓ | 0.4678 |
| NPIP mean fraction of the node that is the copy | 0.6979 | **0.7682** ↑ | 0.6779 |
| NPIP strand correct | 23 | **22** | 18 |
| copies sharing a node | 4 | **4** (same pairs) | 4 |
| nodes shared by ≥2 copies | 3 | **3** | 3 |
| TBC1D3 (19 copies) | 19 with node, 15 full, 0.9558 / 1.0000 | **identical in every cell** | identical |
| **nodes spanning ≥2 annotated loci** | **196** | **179** | **177** |

### 3.3 How much of it is just deleting readthrough nodes? (R3.6, answered plainly)

**For the multi-locus node count: all of it, and more.** The cheap control alone takes 196 → 177; the genuine rebuild lands at 179, i.e. it puts two multi-locus nodes *back*. Of the control's 216 multi-locus nodes (full universe), 43 involve a readthrough record and only **20** exist *because* one is annotated.

**For copy recovery: none of it — and the cheap control is actively wrong.** Deleting the readthrough-derived nodes costs **four NPIPA copies their node entirely**: NPIPA1, NPIPA6, NPIPA8, NPIPA9 all lose theirs (26 → 22), because in the control the only node covering them *is* the readthrough node. The real rebuild keeps all 26. **The difference between "delete the readthrough nodes" and "a world with no readthroughs" is exactly those four NPIPA copies** — and that difference is the whole reason the counterfactual had to be built rather than simulated by deletion.

### 3.4 Who wins, who loses (fraction of the copy covered / fraction of the node that is the copy)

| member | control | treatment | |
|---|---|---|---|
| NPIPA1 | 0.942 / **0.139** | **1.000 / 1.000** | escapes the PKD1P3-NPIPA1 node |
| NPIPA8 | 0.755 / **0.169** | **1.000 / 1.000** | escapes the PKD1P4-NPIPA8 node |
| NPIPA9 | 0.490 / **0.174** | 0.541 / **1.000** | escapes the PKD1P5-LOC105376752 node |
| NPIPA6 | **0.916** / 0.188 | **0.066 / 0.015** | loses the LOC131696449 node, gets a fragment |
| LOC128966608 | **1.000** / 0.400 | **0.034 / 0.022** | collapses |
| NPIPB9 | 0.369 / 0.314 | **0.020 / 0.024** | collapses |
| NPIPB14P | 0.122 / 0.577 | 0.015 / 0.725 | marginal |
| other 19 NPIP members | — | **unchanged node, to the base pair** | — |

Three copies are rescued from being a minority passenger inside a readthrough node and become exact; three collapse to fragments. That is the honest shape of the result: **not an improvement, a redistribution.**

### 3.5 Readthrough loci do not leave the node set

58 treatment nodes still overlap a (now non-existent) readthrough record's exons by ≥ 50 bp, against 44 in the control. The readthrough footprint does not disappear — it **fragments into the nodes of its constituent genes**. This is the node-side statement of §2.2: the sequence that makes the boundary edge is still there, it is just distributed differently.

### 3.6 What arm B did not answer

**NOT RUN, declared at A2.3 and invoked at A3.1:** the de novo NPIP/TBC1D3 certificates at every DNA level on the treatment node set, and the FAMILY R / P strict / F strict panel. Both need the DNA edge graph over the de novo nodes, which the shipped binary does not emit here (its edge stage is stubbed, `npip_ideal` A1.3), so they would require re-running the hop0–hop2 query/mapping/record/certificate chain against a new node set — more than the session's budget could adapt and validate, and a half-checked certificate is worse than none. **This is a gap, not a null result, and must not be quoted as "no change".** For context only, the control de novo certificates (`npip_ideal/cert2/certificates.C.tsv`, test `none`) already give NPIP `h_join = 1.000000` and an empty interval at every level, with `h_split = -inf` and 24+3 internal parts at L3.

---

## 4. ARM C — are the V4 family splits readthrough-related?

**Gorilla: the question is not even askable, which is itself the finding.** `GGO_genomic.gff` contains **zero** records described as readthrough — 0 occurrences of the string in 693 MB, and 0 among the 4,477 gene/pseudogene records on the three analysed contigs. R1.3's declared fallback fires and the structural surrogate (a node whose exons overlap the exons of ≥ 2 distinct same-strand annotated gene records) is used, labelled as a surrogate.

| V0 family | size → V4 destinations | fragment nodes | readthrough record? | surrogate multi-locus nodes |
|---|---|---|---|---|
| **2** | 21 → fam 2 (19) + fam 55 (2) | `NC_073242.2:99634383-99636626 + 1ex` (LOC129527693), `NC_073242.2:103962940-103964751 + 1ex` (LOC101134557) | **NO** — none exists | **0 of 21** |
| **4** | 11 → fam 5 (8) + fam 3 (3) | `:18446771-18453348`, `:18512682-18512869`, `:66803281-66808181`, all `+` 1ex | **NO** | **0 of 11** (every node overlaps *no* annotated record at all) |
| **11** | 4 → fam 45 (2) + fam 7 (2) | `:16302400-16303132 + 1ex` (LOC129527569), `:31483934-31485966 + 1ex` (LOC115933349) | **NO** | **0 of 4** |

**Chimp: testable, and also no.** `PTR_genomic.gff` carries exactly **one** readthrough record genome-wide — **CORO7-PAM16** — and it lies on NC_072416.2, the contig the chimp node sets cover (1,582 gene/pseudogene records there). Two V0 families split under V4 (43 → 42 families): family 0 (42 → 41 + 1; breakaway `NC_072416.2:36607651-36630313 − 16ex`, locus LOC112207901) and family 37 (2 → 1 + 1; breakaway `:44223894-44226204 + 1ex`, no annotated locus). **Neither fragment nor any holding node overlaps CORO7-PAM16.**

One nuance worth keeping: in chimp family 0, **five holding nodes do span two same-strand annotated loci** (`:21129062-21205261` over LOC129137282+LOC112207762; `:30970490-30985206` over LOC749096+LOC454016; `:33348609-33363597` over LOC129137438+LOC129137434; `:33859267-33894887` over LOC129137443+LOC129137444; `:35598715-35633908` over LOC112207360+LOC129137455). None is a readthrough record, and every one of them is *holding the family together*, not breaking it apart. On gorilla not even that happens.

**Verdict: the V4 splits are not readthrough-related on either species.** They are what the V4 report already described — wrong-strand single-exon placeholders superseded by correct-strand spliced nodes — and this diagnostic adds that on gorilla the readthrough hypothesis is not merely unsupported but *unavailable*.

---

## 5. Readings

1. **Readthroughs are not what breaks the family definition; they are the last layer of something that has several.** Strip them and NPIP's `h_join` does not move off 1.000000 — the small junction-poor lncRNAs underneath were always there. Strip those too (t_J) and the boundary lands on PKD1 through NPIP's own member. Each removal reveals the next holder; none of them is the cause.
2. **The counterfactual is not constructible on real data without deleting a member.** PKD1P6-NPIPP1 *is* an NPIP copy and *is* a readthrough. A world with no readthroughs is a world with 25 NPIP copies, not 26. That is the sharpest thing this run found, and it is a fact about the family, not about the method.
3. **"Delete readthrough nodes" and "a world with no readthroughs" are different operations, and the difference is four NPIPA copies.** Post-hoc deletion strips NPIPA1/A6/A8/A9 of any node; the genuine rebuild gives all 26 copies a node. Any future rule that removes readthrough nodes from a *built* catalog will pay that price; a pipeline that never saw them will not.
4. **The rebuild makes nodes cleaner and copies less covered.** Node purity 0.698 → 0.768, full-length nodes 7 → 9, multi-locus nodes 196 → 179 — against mean copy coverage 0.591 → 0.518 and three copies collapsing to fragments. Readthrough transcripts were supplying the long full-length evidence for part of NPIPA; without them those copies are recovered by their own, shorter reads.
5. **TBC1D3 is untouched at every level by every form of the counterfactual.** If readthroughs were the general obstruction, TBC1D3 would have moved. It did not; its failure is internal fragmentation.
6. **The one genuine positive is a cut question, not a rule question.** The member-preserving world plus t_J makes NPIP an exact component of the L2 filtration for any cut in (0.746986, 0.989040] — width 0.242 — with all 26 members and zero member pairs lost. The shipped 0.30 is below it. Chasing that would mean fitting a cut to one development family.
7. **The readthrough footprint does not leave the node set when the annotation does.** 58 treatment nodes still sit on it, against 44 in the control. Removing the record removes the label, not the sequence.

---

## 6. Provenance

- **Declarations** `no_readthrough/DECLARATIONS.txt`, frozen 12:39:25 −07:00, md5 `a6f544c392bede6d5797d4b73147940d` before the first number; final `5464335ea9d50a9113e3a3a16893e200`. Addendum A1 (L1-identity shipped cut is 0.80 not 0.90; machinery reuse) precedes the first certificate. A2 discloses a discarded run launched with `RUSTLE_SD_READ_LOCUS_SPLIT=1` (killed at the rep stage before any node table; log at `logs/rust_nort.discarded.err`), the dropped control re-derivation, the family-panel operationalisation and the "readthrough-derived node" definition. A3 records what was not run.
- **Evidence reused unchanged:** `family_cert/dna/{nodes.tsv,witnesses.tsv,batches/*.paf}`, `family_cert/cert/{certificates,components,dna_pairs.extended}.tsv`, `npip_ideal/{locus_set.tsv,transcripts.tsv,manifest.IDEAL.tsv.gz,bam/ideal.bam,nodes/_nodes.json,cert2/*}`, `v4_gorilla/out/{graph,nodes,panel_detail}.pkl`, `v5_retire/out/ptr_{graph,nodes}.pkl`, `winloci_data/{GGO_genomic.gff,PTR_genomic.gff,Reference/chm13v2.0_RefSeq_full.gff.gz}`, `npip_ladder/idx/target.{fa,splice.mmi,asm20.mmi}`. `rule_adopt/engine.py` was imported, not reimplemented.
- **New outputs:** `out/armA_certificates.tsv` (md5 `3f71e067c6ec8a7857721e18f939bf57`, 144 rows), `out/armA_diag.tsv` (`e6b54d23b709ace9010259f729c42dbc`), `out/armB_nodes.json` (`d34799c1956bf259a9303e2ab7e894e7`), `out/armC_gorilla.json` (`ceffdab273496ee58bc2c3dd87051491`), `out/armC_chimp.json` (`cbe9e83099a719a37db08fddfd9e03b7`), `out/rt_loci_ideal.txt`, `nodes/sd_nodes_tx.fa` (`ab7d6bbc7ce9a6ef2d862f49e1426529`, 5,444 records) and `nodes/sd_nodes_body.fa`, `bam_nort.bam` + `.bai`. Scripts in `code/`.
- **Compute:** one minimap2-free catalog run of the shipped `gw_family_catalog.bin` (≈ 6 min wall, peak RSS 21.9 GB, rc 1 at the stubbed edge stage exactly as the control), one BAM filter (48 s), everything else Python in seconds. No `minimap2` was run at all — every query was served from captured PAFs or not needed. The whole-genome `target.fa` (3.1 GB) was written to `no_readthrough/tmp` and deleted after the node FASTAs were copied out.

## 7. Caveats

1. **Arm B's certificates and family panel were NOT RUN** (§3.6). The certificate question for the counterfactual is answered on arm A alone, and condition (iii) is reported on the node panel instead of FAMILY R / P strict / F strict.
2. **Arm B's control was quoted, not re-derived** (A2.2). Same binary, flags, wrapper and indexes; only the stderr banners are compared directly.
3. **Arm A reuses `rule_adopt/engine.py`** (A1.2) — an independently verified code path, re-anchored here against the committed certificates with 0/36 differences, but not an independent implementation.
4. **Arm B inherits `npip_ideal` C0 verbatim:** reads are simulated from the annotation that also defines the truth copies, so the arm cannot measure how faithfully node construction recovers annotation. Every arm-B number is a ceiling.
5. **Bipartite F and FAMILY R rise mechanically when nodes are deleted** (R5.6). Every F rise here is subject to that.
6. **Two development families, one haplotype, one annotation.** Component sizes inherit `family_cert` D3 and are lower bounds; L2/L3 rows inherit the aggregation caveat.
7. **Arm C's chimp tables were read from another agent's concurrently running study** (`v5_retire/out/`) as read-only inputs; the splits and fragments were re-derived here, the pickles were not verified.
8. **Gorilla's arm-C answer rests on a structural surrogate**, because gorilla's annotation has no readthrough class. A surrogate negative is weaker than an R1.1 negative.
9. **Nothing here recommends dropping readthrough records from real analyses.** The user's framing — readthroughs are real biology and cannot be ignored — is accepted, and this run's own finding is that deleting them on real data deletes an NPIP copy.
## Verification (independent recompute)

Agent 2 of 2, 2026-09-18, 20:02-20:31 UTC. My declarations were written before my first number (`/mnt/linuxdisk/home/juanfraitu/no_readthrough/verify/DECLARATIONS.verifier.txt`, md5 `ffafb35d75d43d4e5c89fd46e476ef55`, 20:04:38Z). No script under `no_readthrough/code/` was read or imported; my code is in `/mnt/linuxdisk/home/juanfraitu/no_readthrough/verify/code/` (`v1_rt.py`, `v2_armA.py`, `v3_boundary.py`, `v5_reads.py`, `v6_span.py`, `v6b_span.py`, `v7_nodes.py`, `v7b_nodes.py`, `v8_percopy.py`, `v9_support.py`, `v10_mirror.py`), outputs in `verify/out/`. Where I reuse earlier runs' verifier code (`rule_adopt/verify/w1b_real_rows.pkl` record rows, `npip_ideal/verify/v4_nodes.py` node mirror) it is disclosed below: those recomputes are independent of agent 1, not of every previous agent.

**1. The readthrough definition, re-derived (agrees).** Parsing `chm13v2.0_RefSeq_full.gff.gz` myself: 58,563 gene+pseudogene records, **209** with "readthrough" in the case-folded description. All 209 are primary nodes, and the set agrees with `family_cert/dna/nodes.tsv`'s `readthrough` column **209/209, with 0 node-only and 0 GFF-only disagreements**. Exactly one is a member: **gene-PKD1P6-NPIPP1** ("PKD1P6-NPIPP1 readthrough", NPIP); TBC1D3 has none. A leak test the run did not do: all **177** fusion-named records (`A-B` where both halves are gene names, e.g. MROH7-TTC4) are inside the 209, so the description selector misses none of them.

**2. Arm A, every certificate row recomputed (agrees).** From the record table `rule_adopt/verify/w1b_real_rows.pkl` with my own drop sets and my own implementation of D1 (h_join = max boundary weight, h_split = max-spanning-forest bottleneck, component at the shipped cut), 4 levels x 9 sets x 6 arms, written to `verify/out/v2b_armA.tsv`. Control anchor reproduced exactly: NPIP 1.000000/0.971665 (L1-id), 1.000000/1.000000 (L1-cov), 1.000000/0.989040 (L2), 1.000000/0.980559 (L3); TBC1D3 1.000000/0.840488, 1.000000/0.841369, 1.000000/0.832037, 0.961397/0.867946. **NO_RT alone: h_join stays at exactly 1.000000 at all four levels**, members 26 -> 25, outside 96 -> 85 (L1), 64 -> 51 (L2), 57 -> 48 (L3) — every figure as reported. **NO_RT + t_J (STRICT k=2, delta=0): h_join = 0.000000 at all four levels**, 0 outside, shipped cut inside at all four — the certificate is vacuous and bought with the member, as the report says. **Member-preserving + t_J**: L1-id 0.978223/0.971665 (empty by 0.006558), L1-cov 1.000000/1.000000 (empty), **L2 0.746986/0.989040 — non-empty, shipped cut 0.30 below the window**, L3 0.988804/0.980559 (empty by 0.008245); the outside is exactly {PKD1 0.978223, PKD1P2 0.971834, PKD1P3 0.869490, PKD1P1 0.868219}, every edge via PKD1P6-NPIPP1. TBC1D3's h_join is bit-identical control vs NO_RT at all four levels, outside 19 -> 14 / 10 -> 10 / 0 -> 0. One caveat on method: two independent record tables exist for the junction conjunct, and they disagree (`w1_real`/`v6` puts NPIP's t_J h_split at 0.934476, `w1b`/`v6b` at 0.971665). I used `w1b`, the corrected one that reproduces `engine.py` and the committed tables; on the earlier table the L2 window would read (0.674086, 0.853211] instead of (0.746986, 0.989040]. The verdict is identical either way.

**3. Arm B substrate and BAM, checked on the whole file (agrees).** My readthrough set intersected with the 5,542-locus set is the same 40 loci, name for name. Manifest: 38,627 transcripts, 91 readthrough, 2,730 reads. BAM: control 2,789,740 records, treatment 2,778,640, **difference exactly 11,100**; **0** records of any readthrough locus survive anywhere in the treatment file; distinct source genes 5,542 -> 5,502, and the removed set is exactly my 40 loci; header byte-identical excluding @PG (md5 `2c985bdf90807980d21adc7bf04f3e97` both). The @PG chain confirms the declared shortcut (a `samtools view` filter of `ideal.bam`), as R3.3 states. I regenerated 60 random surviving reads from `idx/target.fa` + GFF exons: **60/60 md5-identical to the manifest**.

**4. Arm B nodes, rebuilt with my own mirror and rescored (agrees, with the construction caveat of corrections 3-5).** Node exon structures taken from the FASTA headers and the control node table, locus exon unions derived by me from the GFF (validated by reproducing the control table's per-locus overlap sets on **5,425 of 5,428** nodes). Reproduced: nodes spanning >= 2 loci 216 (control, 5,542 universe), 196 (control, 5,502), 177 (cheap control, 38 nodes dropped), 179 (treatment); NPIP with_node 26/26 control and treatment vs 22/26 cheap; full-length 7 -> 9; means 0.5912/0.6979 -> 0.5175/0.7682. Node construction mirrored independently from the binary's own `[rep-audit]` stderr (5,389 treatment reps, 5,375 control) through a reconstruct/consolidate mirror: **5,149 of agent 1's treatment nodes reproduced exactly by rep-exon key**, with the same residual profile as the same mirror run on the published control (5,135 matched, 8 mine-only, 293 unreconstructable reps + read-locus nodes) — the treatment node set is no less faithful to the binary than the control is. The mechanism reproduces directly: the control's single 42,950 bp / 39-exon node at chr16:14910160 (the PKD1P3-NPIPA1 readthrough) is replaced in the treatment by a 14,597 bp / 8-exon node that is NPIPA1 alone, and likewise for NPIPA8 and NPIPA9. 12/12 sampled treatment-only nodes carry 30-300 primary reads in the treatment BAM and each is the shortened form of a longer control node.

**5. Is the improvement real, or the removal of the nodes being scored?** Arm A: **the removal, entirely.** `NO_RT` *is* the control restricted to the surviving node set, and on that restricted set NPIP's h_join does not move off 1.000000 at any level; what moves is the precision denominator — the component at the shipped cut loses 11 (L1), 13 (L2) and 9 (L3) non-member nodes while the member subgraph is unchanged apart from the deleted member. That is exactly the rise R5.6 pre-named, and it carries nothing. The one genuine graph change in arm A is the member-preserving L2 window (h_join 1.000000 -> 0.746986), and it is unusable: the shipped cut sits below it and the surviving boundary is internal, through PKD1P6-NPIPP1's PKD1 half. Arm B: **mixed, and the report's own attribution holds.** On the metric being scored, restricting the control to the surviving node set already gives 177 nodes spanning two loci against a matched-construction control of 194 and a rebuilt treatment of 179 — i.e. deletion alone explains more than the rebuild achieves. What deletion cannot do, and the rebuild can, is give NPIPA1 and NPIPA8 clean full-length nodes: the cheap control destroys four NPIP copies' nodes (26 -> 22) where the rebuild keeps all 26. So the only real improvement in the whole run is two nodes, bought against one copy (NPIPA6, 0.916 -> 0.066) that is lost.

**6. Traps.** (i) A certificate non-empty only because a member was deleted: present and correctly flagged — the four certifying rows have h_join = 0.000000 over 25 members with no boundary edge at all. (ii) A boundary that moves to another class: present, and the class is misnamed at L1 (correction 1). (iii) A family metric improving by denominator shrinkage: present in every arm-A F number; quantified above. (iv) Quoting a retired metric: present (bipF_soto, correction 6).

**7. Scope of this verification.** I did not re-run the Rust catalog binary; the arm-B node sets are verified against the binary's own rep audit and against the treatment BAM, not by a fresh binary run. I did not reproduce any bipF number. Arm C is verified only at the level that carries its answer: gorilla `GGO_genomic.gff` has **0** readthrough-described records in 41,193 gene+pseudogene lines, chimp `PTR_genomic.gff` has exactly **1**, gene-CORO7-PAM16, in 41,815; the per-node surrogate rows of `armC_gorilla.json` / `armC_chimp.json` were not independently recomputed. Everything remains descriptive on one assembly, one annotation and two development families.