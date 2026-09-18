# Agent 3 of 4 — the certificates

Declarations appended to `/mnt/linuxdisk/home/juanfraitu/ggo_npip/DECLARATIONS.txt` as sections **A3.0–A3.7** at 2026-09-17T17:08:37-07:00 (file md5 at that moment `877e853e187a8de3872670bc8fac84f8`), **before any certificate, weight, h_join, h_split, component count or boundary listing existed**. They contain the A3.1 decision rule verbatim, the truth sets, the member→node map, the four axes, the conjunct definitions, the D1 certificate, and three declared deviations (D-a, D-b, D-c) from the human `family_cert` engine. Outputs: `/mnt/linuxdisk/home/juanfraitu/ggo_npip/cert/{cert.py,run.py,GGO.certificates.tsv,PTR.certificates.tsv,MANIFEST_agent3.txt}` — 324 + 288 = 612 certificate rows. No new alignment was run; nothing under `/mnt/c/Users/jfris/Desktop/Rustle` was written; no commit; no subagent.

> **Verifier corrections, applied by the orchestrator** (agent 3 could not write this file; the verifier returned ok = false with the following, none of which overturns the verdict):
> 1. **Chimp labelled count is 3, not 4.** NC_072416.2:30,970,428-31,260,911 is claimed by three liftoff records (NPIPB10P, NPIPB7, NPIPB9) and is UNRESOLVED under the declared rule.
> 2. **`native_introns` is not the declared quantity** (max over the record's transcripts of exons − 1); recomputed from the GFF it differs at 9 of 26 gorilla loci.
> 3. **Three chimp rows violate the declared h_split definition**: body-chain rows were given f_ex = 0.0 instead of being excluded, so h_split prints 0.000000 where the definition gives −inf (S internally disconnected).
> 4. **The junction-match count undercounts**: recomputed from the captured PAFs, 16 of 2,500 de novo and 35 of 8,595 guided transcript rows disagree, nearly all by exactly 1; 8 de novo and 9 guided rows flip the k = 2 admission. Treat every "0 of N pairs carry 2 junctions" statement as ±1 row.
> 5. **`comp_size_at_cut` is the total size of all components holding a member**, not one component's size, whenever the set is split.
> 6. **Declaration ordering is verifiable for agent 3 only** (DECLARATIONS.txt has a single mtime); agents 1-2 sections cannot be timestamp-checked.
> 7. **Human-side framing**: the repo doc says five named blockers lose their NPIP edge at k ≥ 2; six records hold the 1.000000 layer. Do not write "all six blockers removed".

## 0. Substrate honesty (declared in A3.0, repeated here as required)

Every gorilla NPIP locus in agent 1's table sits on **NC_073242.2 (24), NC_073244.2 (1), NC_073241.2 (1)** — the **same contig trio the shared definition and its node-construction rules were developed on** (ledger §6jy–§6kb, prereg Addenda AB–AD). Gorilla therefore holds out the **species** and the **conjunct** t_J (developed on human NPIP only); it does **not** hold out node construction. **Chimp (PTR, GCF_028858775.2) is the only substrate here never used for any of this work.** Gorilla de novo edges can only be found inside the contig trio; chimp de novo edges only inside NC_072416.2 (agent 2's A2.7). Both are declared search-space restrictions.

Per the user's instruction, **liftoff names are carried but EXPLORATORY**: `T_strict` = the 2 gorilla loci both arms name identically (NPIPA2, NPIPB2); `T_expl` = the 25 gorilla loci where both arms land, membership only. `T_strict` has |S| = 2 in gorilla, so its h_split is the weight of a single internal edge — every T_strict row is structurally weak and is never presented as equivalent to the human 26-member certificate.

## 1. PRIMARY RESULT — the gorilla de novo arm

25 T_expl loci; **12 have a de novo node, 13 have none** (agent 2's node-recovery result, inherited). 0 loci share a node. Levels use the shipped cuts 0.80 / 0.50 / 0.30 / 0.98.

| level / axis | BASELINE h_join / h_split | STRICT k=2 δ=0 | interval | member pairs at cut | pairs lost | internal parts | outside | labelled isolated |
|---|---|---|---|---|---|---|---|---|
| L1 identity (0.80) | 0.997868 / −inf | **0.000000 / −inf** | empty | 21 → 0 | **21** | 11,1 → **1×12** | 14 → 0 | 0 → **2** |
| L1 coverage (0.50) | 1.000000 / −inf | **0.000000 / −inf** | empty | 21 → 0 | **21** | 11,1 → 1×12 | 14 → 0 | 0 → **2** |
| L2 f_ex (0.30) | 1.000000 / −inf | **0.000000 / −inf** | empty | 17 → 0 | **17** | 11,1 → 1×12 | 12 → 0 | 0 → **2** |
| L3 w_98 (0.98) | 0.997868 / −inf | **0.000000 / −inf** | empty | 5 → 0 | **5** | 8,2,1,1 → 1×12 | 5 → 0 | 0 → **2** |

On the 2-locus declared truth `T_strict`: BASELINE L1-identity 0.997039 / 0.986948 (empty); STRICT 0.000000 / −inf, the one member pair lost, **both labelled copies isolated**.

**Verdict against A3.1.** h_join is strictly lower — but that is the only clause satisfied. A labelled copy **is** isolated (both are), the set is **not** in one component (12 singletons), and the boundary is emptied only because the member subgraph is emptied. **NOT CONFIRMED.**

### Why — and it is node quality, not rule strictness

Junction counts of the 12 gorilla de novo member nodes (from the node's own exon structure, the conjunct's only legal source):

| node | junctions | exploratory name |
|---|---|---|
| NC_073242.2:29391554-29415846 | 20 | NPIPA5 |
| NC_073242.2:32231716-32439027 | 9 | NPIPA2 |
| NC_073242.2:21074187-21095824 | 5 | NPIPB4 |
| NC_073242.2:32148199-32150876 | 2 | LOC124907808 |
| NC_073242.2:21765397-21769412 | 1 | NPIPB11 |
| 7 further nodes | **0** | NPIPB8, NPIPB2, NPIPB3, NPIPB15, NPIPA9, LOC124907807, LOC128966608 |

- Only **4 of 12** member nodes can satisfy k = 2 even in principle.
- Those 4 share **zero** witness records with one another.
- All **56** witness rows between member nodes have `jmatch_d0 = 0`; 39 have u_junc = 0, 39 have v_junc = 0, 24 have both = 0.
- The node with the most junctions (NPIPA5, 20) has **degree 0 at baseline** — it is already isolated before the conjunct is applied, so the de novo T_expl set was never one component to begin with.

So on gorilla de novo nodes the junction evidence and the homology evidence are **disjoint across the NPIP node set**: the junction-bearing nodes have no edges and the edge-bearing nodes have no junctions.

**The conjunct is not globally vacuous.** Over the whole gorilla de novo graph it keeps **286 of 1,248 pairs (22.9 %)**. It empties the NPIP part specifically.

GUARDED on this arm behaves as in human — it achieves almost nothing (h_join stays 0.997868 at L1 identity, 2 pairs lost, 0 isolated).

## 2. SECONDARY — gorilla guided arms

### 2a. Guided-annotated (GA), 25/25 loci present — the human result replicates

| level / axis | BASELINE | STRICT k=2 | GUARDED | interval (STRICT) | pairs lost | outside | isolated | components |
|---|---|---|---|---|---|---|---|---|
| L1 identity | 1.000000 / 0.953568 | **0.988399** / 0.953568 | 0.988399 | empty | **0** | 48 → **3** | 0 | 1 |
| L1 coverage | 1.000000 / 0.989399 | 1.000000 / 0.989399 | 1.000000 | empty | 0 | 48 → 3 | 0 | 1 |
| L2 f_ex | 0.981969 / 0.683746 | 0.981969 / 0.683746 | 0.981969 | empty | **10** | 25 → 3 | 0 | 1 |
| L3 w_98 | 1.000000 / 0.897445 | **0.989975** / 0.897445 | 0.989975 | empty | **2** | 20 → 2 | 0 | 1 |

Free at L1 (0 member pairs lost, 0 isolated, one component), h_join falls, and it **never certifies** — h_join stays above h_split at every level. This is the same shape as the human arm REAL (1.000000 → 0.999363, outside 96 → 10, 0 members lost, never certifies), with the gorilla drop being larger in absolute terms (0.0116 vs 0.000637) but reaching a still-far-above-h_split value.

The L2 cost, named: 10 pairs lost, all involving LOC109023568 (7 of them), plus LOC101134557–LOC129527693, LOC115933350–NPIPB11, LOC129523555–LOC129527720. L3 loses LOC109023568–LOC115933079 and LOC109023568–NPIPB11.

### 2b. Boundary composition, compared with human

**Removed by STRICT (45 of the 48 baseline outside nodes):** 20 are 1-or-2-exon records; median exon_bp of the 45 is 1,310; 2 are fully embedded inside a member's copy interval. The smallest:

| node | exons | exon bp | native description |
|---|---|---|---|
| LOC115931156, LOC115932719, LOC115933357, LOC129527842-6 (8 records) | 1 | 104 | small nucleolar RNA U13 |
| **LOC115932701** | 2 | **258** | uncharacterized; sits at NC_073242.2:32451238-32456187, **inside the NPIPA2 copy interval**; held h_join at 1.000000 |
| LOC115932780, LOC129527572, LOC129527571, LOC134757426 | 1 | 477–700 | large ribosomal subunit protein uL22-like |
| LOC101131206 | 4 | 766 | lncRNA, at NC_073242.2:28301414-28305185, **inside the NPIPB8 copy interval** |
| LOC101134912 / LOC129527726 / LOC129527593 (+8 more) | 21 / 6 / 12 | 3,096 / 1,164 / 3,695 | serine/threonine-protein kinase SMG1 and SMG1-like |

**Survivors (the new boundary):** LOC101130854 (MRP1, 33 exons, 6,765 exon bp, 190 kb span), LOC115932992 (uncharacterized, 41 exons, 5,679 bp, 180 kb), LOC129528916 (MRP1-like, 5 exons, 721 bp, 54 kb).

**Comparison with the human side.** Human removed CLN3 (253 bp, 1 block), EIF3CL, LOC100190986 (2,453 bp lncRNA, 1 block), LOC124907830/845 (2 blocks each), LOC128966632 (5,598 bp "SMG1-like"); the boundary then passed to **PKD1P–NPIP readthroughs** (0.999363). Gorilla reproduces the *removal* class exactly — 104 bp snoRNAs, 258 bp/2-exon embedded record (the size twin of CLN3's 253 bp), 477–700 bp single-exon ribosomal pseudogenes, an embedded lncRNA, and the same SMG1-like family — and reproduces the *replacement* clause: what remains is **full-size multi-exon neighbours, not small embedded records**. The one structural difference is that **no readthrough queue takes over in gorilla**, because the gorilla PKD1P6-NPIPP1 readthrough locus is the liftoff-only locus and is in neither truth set.

**GUARDED is NOT refuted in gorilla the way it is in human.** It gives the *same* h_join as STRICT at L1-identity (0.988399) and L3 (0.989975). The reason is a concrete cross-species difference: the human 1.000000 blockers were junction-**less** (CLN3 1 block, LOC100190986 1 block), so the guard re-admitted them; gorilla's top blocker LOC115932701 has **exactly one** junction, so the guard does not exempt it and k = 2 still fails it. GUARDED still re-admits the junction-less crowd (outside 40 vs STRICT's 3) — it just does not re-admit the one that holds the boundary.

### 2c. Guided read-supported (GR) — destructive, worse than human

| level | BASELINE | STRICT | pairs lost | parts | labelled isolated |
|---|---|---|---|---|---|
| L1 identity | 1.000000 / 0.953568 | 0.988395 / **−inf** | **269 of 270** | 2 + 23 singletons | **2** |
| L2 | 0.981969 / 0.683746 | 0.981969 / −inf | 97 of 98 | idem | 2 |
| L3 | 1.000000 / 0.897445 | 0.989975 / −inf | 22 of 23 | idem | 2 |

Cause: **14 of the 25 member native records carry ZERO junctions supported by ≥ 3 primary (-F 2308) Iso-Seq reads** (GGO_ds.bam is 37 %-downsampled). Human's IDEAL-READ arm had exactly one such member (NPIPB12). GUARDED on GR restores h_join to 1.000000 with 39 pairs lost and 0 isolated — **the human "the guard re-admits precisely what the conjunct exists to remove" pattern, replicated on a new species.**

## 3. CHIMP — the never-used substrate

16 both-arms loci on NC_072416.2 (4 of them labelled), all 16 present in all three arms.

| arm | BASELINE h_join / h_split | STRICT h_join / h_split | pairs lost | parts | outside | labelled isolated |
|---|---|---|---|---|---|---|
| de novo | 1.000000 / −inf | **1.000000** / −inf | 57 of 85 | 14,1,1 → 8 + 8 singletons | 50 → 22 | 0 → **1** |
| guided-annotated | 1.000000 / 0.960879 | **1.000000** / −inf | 29 of 98 | 16 → 13,3 | 111 → 63 | 0 |
| guided-read-supported | 1.000000 / 0.960879 | **1.000000** / −inf | 60 of 98 | 16 → 9,3,1,1,1,1 | 111 → 49 | 0 |

**h_join does not move on any chimp arm, at any level.** The reason is in the native annotation and is not a rule failure: the nodes holding 1.000000 under STRICT — **LOC112204066, LOC129137441, LOC112205827, LOC750787, LOC129135381** — are all described by chimp RefSeq itself as *"nuclear pore complex-interacting protein family member B12 / B12-like / B15-like"*. **29 of the 63 outside nodes that survive STRICT on the guided arm are natively annotated NPIP records.** They are outside S only because agent 1's two-arm label rule has a ceiling of 27 human source records and never proposed them as candidates. The conjunct cannot cut an edge to a genuine family member, and should not. This is **truth incompleteness, declared, not evidence about the rule** — and it means chimp gives no usable read on the "h_join strictly lower" clause. What chimp *does* show is the cost side unambiguously: the conjunct shatters the set (2–6 components) and isolates a labelled copy on the de novo arm even where the chimp BAM is undownsampled (249,539 primary reads vs gorilla's 140,223 — declared in A2.7 as not like-for-like).

The gorilla contrast is instructive: the gorilla guided baseline outside set contains **0** records the gorilla annotation names as NPIP-family (gorilla RefSeq calls NPIP copies "titin-like"/"SRRM2-like"/"NACA-like"), so gorilla's 48 → 3 is a genuine reduction in foreign records, while chimp's 111 → 63 is dominated by real members the truth could not reach.

## 4. Certificates — the count, and the degeneracy

Of **612 rows** (2 species × 3 arms × 3 variants × 4 axes × 8–9 sets):

- **7 rows** have a non-empty interval. **All 7** are on the gorilla de novo `T_strict` set, |S| = 2.
- **4 rows** have the shipped cut inside: gorilla de novo `T_strict`, variants **BASELINE** and **GUARDED**, levels **L2** (0.000000, 0.798879] and **L3** (0.000000, 0.988792].
- **All four have h_join = 0.000000**, because that 2-node set has no boundary edge at all at those levels — the interval (0, h_split] contains every cut trivially. This is exactly the degeneracy the human verifier flagged as correction #2 of `bench/JUNCTION_AND_READTHROUGH_RULES.md`. **They are not certificates of anything and are reported only so the count is complete.**
- **No STRICT row certifies, in either species, on any arm, at any level.** The human finding "it never certifies NPIP" replicates.

**Exploratory subfamily windows (A3.1 marks these exploratory only).** On the gorilla guided-annotated arm STRICT lowers h_join for every window at zero pair cost — NPIPA 0.997033 → 0.988399 (h_split 0.967514), NPIPB 0.999105 → 0.995090 (0.942293), iso A6-9 0.997033 → 0.979408 (0.967514), iso B12/13 0.999105 → 0.977901 (0.963376), iso B15 1.000000 → 0.988395 (0.881233), iso B3-5 0.991606 → 0.979408 (0.968655), iso B6-9 0.995943 → 0.972573 (0.935022). **Every one is still empty**: h_join stays above h_split in all 7. The two human §0★★★.7d cut windows (the four Dishuck groups at L3, NPIPB missing 0.98 by 0.000559) **do not replicate in gorilla** — 0 windows found. On the de novo arm every window collapses to singletons.

## 5. What the advisor should be told about the primary arm

The de novo arm's answer is not "the rule is too strict". It is that **the de novo nodes for gorilla NPIP copies do not carry junctions where the homology is**. Seven of twelve are single-exon fragments; the four junction-rich nodes have no homology records between them; the richest one (NPIPA5, 20 junctions) is already degree-0 at baseline. Any predicate that reads a node's splice structure — t_J or otherwise — is reading a structure the de novo node construction did not produce. That is the same conclusion §6j1/§6je reached by a different route (node presence and node construction are 58 % of the de novo↔guided gap), now reached at the level of a specific proposed rule, on a new species, with the rule's own certificate as the instrument.

The guided-annotated arm is where the rule works, and gorilla confirms that independently of human: 48 foreign records at the boundary become 3, at zero cost in members, pairs or components, and the records it removes are the same *kind* of record the human run named — 104 bp snoRNAs, a 258 bp two-exon record embedded in a copy, single-exon ribosomal pseudogenes, an embedded lncRNA, SMG1-like. It still does not buy a certificate, in either species.

## 6. Declared deviations, restated (A3.4 D-a/D-b/D-c)

- **D-a**: agent 2's witness table holds only records that already pass the shipped filters, so h_join on the identity axis is a **lower bound**; the declared audit (re-scan the captured PAF for sub-threshold boundary records) was conditional on a certificate appearing on the identity axis, and **none did**, so it was not triggered.
- **D-b**: agent 2's body rows are gene-body **chains**, so L2/L3 weights here are **tx-only** and are a lower bound on the human-style f_ex/w_98. This applies identically to BASELINE, STRICT and GUARDED, so every on/off comparison — which is what A3.1 decides on — is internally exact.
- **D-c**: the sense condition and own-locus exclusion were applied by agent 2 when the table was built and were not re-applied.

One refinement of D-a, made while implementing and consistent with it: rather than assuming every row is a t1 row, the per-row checks `identity ≥ 0.80` and `cov_query ≥ 0.50` are applied explicitly on every axis, exactly as the human engine does. No threshold was changed.
## Verification (independent recompute) - agent 4 of 4

Scope: I re-derived the labels, the junction/witness arithmetic and the certificates with my own code under `/mnt/linuxdisk/home/juanfraitu/ggo_npip/verify/` (`v_labels.py`, `v_native.py`, `v_junc*.py`, `v_cert.py`). I did not read any other agent's script; inputs were the data outputs, `DECLARATIONS.txt`, the original `GGO_genomic.gff`, the captured PAFs/FASTAs and `GGO.3ctg.bam`. No src/ change, no commit, nothing written under /mnt/c.

### What reproduces exactly

**Labels, arm (ii) re-derived from the PAF (§4.4/§4.5).** Independent merge-within-100 kb, identity = sum(nmatch)/sum(alnlen), query coverage = union of query blocks / qlen, admit at id >= 0.90 and qcov >= 0.50, then overlap-clustering with the liftoff intervals: **694 raw landings, 522 admitted, 25 identity loci, 26 loci, 2 LABELLED / 23 disagree / 1 liftoff-only**. The interval set is identical to `labels/GGO.truth.tsv` and **all 26 rows agree on status, liftoff name, identity winner and admitted-landing count (0 differences)**. Chimp reproduces 946 / 715 / 31 / 32 with one rule-application correction (see corrections).
**No liftoff-only locus entered the truth**: the PKD1P6-NPIPP1 landing inside gorilla PKD1 carries `both_arms=False` and is in neither T_strict nor T_expl. `PTR.truth.1ctg.tsv` is exactly `PTR.truth.tsv` restricted to NC_072416.2 (28 rows, identical).
**Native annotation spot check, all 26 gorilla loci (not just 5+3), recomputed from the original `GGO_genomic.gff`: 26/26 native names match** (LOC101141990, LOC115933071, LOC129527568, ... , LOC115933039), 26/26 have an overlapping gene record, 26/26 carry >= 2 introns under every reading. Intron *counts* deviate from the declared rule at 9 loci (correction 2).

**Node sets.** The de novo node table is faithful to the run that produced it: its 2,664 node ids are exactly the `body.fa` headers the catalog's own minimap2 call consumed, every `rep_exons` chain appears as a `tx.fa` header, and the `consolidate` invariants hold on the emitted table - 0 intra-node gaps above MAX_INTRON (271,359), 0 nodes below MIN_PIECE (100 bp), **0 same-strand node pairs with overlapping exons**. I could not mirror `shared_definition::build` end-to-end from the BAM, because step 1 consumes the catalog's upstream read-supported reps, not the BAM directly; I verified steps 1-2's output invariants plus a BAM-side probe of step 2 instead (see "what the de novo result actually is", point 3).

**Conjunct arithmetic.** Recomputing jmatch_d0 from the captured PAFs (map each node junction to its spliced-query offset, walk the record's CIGAR, require an N at that offset whose (donor, acceptor) is exactly a junction of v): **2,485/2,500 de novo tx rows and 8,560/8,595 guided tx rows agree exactly**; the residue is correction 4 and changes no certificate row.

**Certificates.** My own D1 engine (h_join = max boundary weight, h_split = maximum-spanning-tree bottleneck via descending Kruskal, +inf at |S|=1, -inf when S is internally disconnected) reproduces **all 72 gorilla T_strict/T_expl rows exactly** - h_join, h_split, n_present, member pairs at the cut, pairs lost vs baseline, one-component - and 69/72 chimp rows (the 3 exceptions are correction 3). Specifically the primary numbers hold: gorilla de novo BASELINE h_join 0.997039 (T_strict) / 0.997868 (T_expl), STRICT h_join 0.000000 with h_split -inf, 0 member pairs, 1 and 21 pairs lost, both labelled copies isolated.

### What the de novo result actually is (the priority arm)

The verdict stands - **t_J at k=2 does not confirm on the de novo arm** - but the gorilla de novo arm does not test the conjunct at all, and the report should say so:

1. **h_join falls to 0 by emptiness, not by blocker removal.** Of the 12 gorilla de novo member nodes, **8 have zero junctions in their fixed definition**, one has 1, one has 2, one has 5, one has 9. Zero of the 185 records incident to member nodes reaches jmatch >= 2; for the two T_strict members, 0 of 5 incident records. At k=2 at most 3 of 12 members could ever be admitted. "h_join strictly lower" here is the prediction-conditioned-denominator trap: every member simply goes to degree 0.
2. **Node recovery, not the conjunct, is the bottleneck.** 12/26 loci have any de novo node, **0/26 are full length**, best frac 0.40. The NPIPA2 member node is NC_073242.2:32231716-32439027 - a 207 kb span, 10,956 exonic bp of which only 2,085 (19%) lie inside the labelled copy, natively annotated LOC115932699; the T_strict baseline boundary runs to **PARN** through that node. The NPIPB2 member node is 1,606 bp = 10.7% of that copy's native exonic bp and has 0 junctions.
3. **A code-level cause, verified from the BAM.** At NC_073242.2:29,415,572-29,453,211 (NPIPA5; 1,443 primary reads, 1,442 spliced, all MAPQ >= 10 - not a MAPQ-0 loss) the depth>=2 read locus is 12,443 exonic bp spanning 29.39-29.59 Mb, but its first segments overlap the pre-existing node 29,391,554-29,415,846, and `with_read_locus_nodes` discards the **whole** read locus when **any** segment touches an existing node. The locus ends up contributing 274 exonic bp. That is Addendum AC node construction - precisely the part gorilla does not hold out.
4. **88% of de novo nodes use the fallback chain.** 2,341 of 2,664 gorilla nodes (and 1,414 of 1,589 chimp nodes) have no `copies.tsv` row, so their `exons` column - the source of t_J's junctions - is the *representative* transcript chain, with `n_reads` written as 0. A2.1 allowed this and said the count would be reported; the count is in neither manifest. Two of the 26 loci, including one T_strict member (NPIPB2), sit on fallback nodes.
5. **GUARDED is not a gorilla positive.** On T_strict it lowers L1 identity 0.997039 -> 0.974347 and coverage 0.835112 -> 0.677312 at 0 pairs lost, but only because the guard re-admits junction-less nodes and one of the two members has 0 junctions; on T_expl h_join is unchanged at 0.997868 with 2 pairs lost. On chimp, GUARDED is byte-identical to BASELINE on every de novo row - reproducing the human "GUARDED achieves nothing" finding on the held-out substrate.

### Chimp is the arm that actually tests t_J

The chimp de novo node set is not degenerate: 26/28 loci carry a node, member nodes carry 6-21 junctions, and 410 of 1,126 incident records reach jmatch >= 2. On that substrate - the only one never used for this work - **STRICT fails the A3.1 confirmation criteria**: h_join moves only once (L1 identity, T_strict, 1.000000 -> 0.996527) and stays at 1.000000 on all other axes and on every T_expl row, while 1 labelled copy is isolated, S is never in one component, and 57 of 85 T_expl member pairs (24 of 39 at L2) are lost. This is the clean cross-species negative; the gorilla de novo arm should be reported as uninformative about the conjunct, with this chimp result carrying the de novo claim.

### Substrate-usage caveat, where it bites

All 26 gorilla loci are on NC_073242.2 (24) / NC_073244.2 (1) / NC_073241.2 (1), the trio the shared definition and its node rules were developed on. Because the gorilla de novo outcome is decided entirely by node construction (points 1-4 above), the arm holds out neither the rules that produced the outcome nor, in effect, the conjunct - t_J is vacuous on junction-less fragment nodes. The gorilla/chimp node-quality gap is additionally confounded by depth (GGO_ds.bam is 37% downsampled, PTR_mm.bam is not; declared in A2.5/A2.7), so "chimp de novo nodes are better" is not a like-for-like statement either. The guided gorilla arms are unaffected by this particular caveat - their nodes are native GFF records - but they are secondary by the user's own ordering.