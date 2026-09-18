# Strand fix for the shared-definition node set: four variants on two substrates (agent 1 of 2)

Declarations written BEFORE any number: `/mnt/linuxdisk/home/juanfraitu/strand_fix/DECLARATIONS.txt`, md5 of the frozen text `e8f354f1573aeaf8f4485c17ec421a85` at 2026-09-17T23:12:28-07:00; ADDENDUM A1/A2 appended at 2026-09-17T23:38 (before any variant metric; md5 now `0b3eafc27c8c11029e4a6a9b7f1405d3`). All outputs under `/mnt/linuxdisk/home/juanfraitu/strand_fix/` (`out/`, `logs/`, `map/`, `reads/`, `aln/`, `rust/`, `scripts/`), `TMPDIR` under `strand_fix/tmp`. Binary `/mnt/linuxdisk/home/juanfraitu/rustle_target/release/gw_family_catalog`, sha256 `fde7cf2b6eaa5ec371eb6fb164901f03f91793d28c7e2220c83c367e8a69f59d` (commit 94ed994e). minimap2 2.30-r1287, samtools 1.22.1. **Nothing in `src/` was modified, nothing was committed, no subagents were spawned, no MAPQ gate was lowered anywhere.** No index was rebuilt; `npip_ladder/idx/target.{splice,asm20}.mmi` were reused.

---

> **Verifier corrections, applied by the orchestrator** (verifier ok = true; all 16 arm × variant node sets reproduced key-for-key, every published number identical):
> 1. **Clause (iii) cannot credit a strand fix — it is a metric trap.** Locus precision is matched/nodes against a fixed expressed-record set, so a correct new node at a locus whose record is already matched by the wrong-strand node lowers P by construction. Every "FAIL (iii)" in the table below is partly this instrument, not harm.
> 2. **V4 (strand-aware suppression restricted to blockers whose strand was measured) was built and scored by the verifier, not the builder**: it reproduces V1/V3's copy gains exactly and still fails the declared rule — so the proposed next step does not pass it either.
> 3. **Clause (i) is not discriminating on the realistic arms**: NPIPB13 and NPIPB4 are already correct there in the baseline, and NPIPB14P is an extent problem no strand rule touches. The copy that is actually strand-blocked on realistic input is NPIPB12.
> 4. **The stated root cause is half right**: the nodes blocking NPIPB13 and NPIPB4 on the idealized substrate are SPLICED models whose strand came from junction motifs, not single-exon placeholders.
> 5. **The V2 "no ts evidence" test was vacuous** — minimap2 under `-uf` tags every primary alignment `ts:A:+` — so the operational test became "no read with an N in the CIGAR"; this was an addendum after the declarations.
> 6. **MAPQ-0 counts in the substrate table are MAPQ < 10 counts**; the exact MAPQ-0 primaries are 2,299 / 2,614 / 2,925.
> 7. **Undeclared mirror approximation**: for nodes Rust formed by merging ≥ 2 pieces (31 / 35 / 24 per arm), the mirror uses the representative chain rather than the merged union.
> 8. The family blast radius covers chr16/17/18 only (48-57 families), not the 121-family panel the old unstranded-collapse kill was measured on.

## 0. The answer to the question that triggered this run

**Yes — put it in the simulation, and it changes the answer.** The idealized substrate says V3 is a clean win. The realistic substrate, built here for the first time, says it is not. The reason is mechanical: on the idealized substrate only **451 of 5,360** base nodes (8.4%) carry the `denovo_assemble.rs` '+' placeholder, so any strand rule only touches a thin edge of the graph. On the realistic substrate **43% / 61% / 71%** of base nodes carry it (U00 / U20 / U40), because 5'-truncated reads land inside a single exon and lose their junctions. That is the regime in which the prior kills (`RUSTLE_READ_STRAND`, `RUSTLE_COLLAPSE_UNSTRANDED`) were made, and it is the regime the idealized substrate cannot reach.

---

## 1. What was verified before anything was measured (parity)

Everything below is a **mirror** number (Python mirror of `shared_definition.rs`: `overlap_groups` / `consolidate` / `ExonIndex` / `with_read_locus_nodes`, plus the `iw_lib` read-isoform widening port). V0 parity was asserted first on every substrate:

- **S-IDEAL, P1**: 5,375 rep-audit reps; 68 read-locus nodes identified (matching the Rust's own `5375 reps -> 5360 gene-level loci + 68 read-locus nodes = 5428 nodes`); mirror V0 = 5,360 + 68 = **5,428 nodes, node-set identity PASS** against the frozen table. The mirror's widening count is **3,896 of 5,428** — bit-identical to the shipped stderr line.
- **S-IDEAL, P2**: `consolidate` is idempotent on the recovered base (5,360 in, 5,360 out, 0 symmetric difference), which is what licenses the declared two-stage consolidate used by V2.
- **S-IDEAL, P3**: V0 restricted to chr16/17/18 = 710 nodes vs the frozen substrate-(a) 710 nodes, **symmetric difference 6 (3 nodes differ)** — exactly the previously disclosed merged-node reconstruction residue.
- **S-REAL**: mirror V0 reproduces the Rust node count and node set **exactly on all three arms** (7,207 / 8,292 / 8,565, node-set identity PASS).
- Downstream, V0 reproduces the published K5 family row on substrate (a): FAMILY R 0.8519, P strict 0.8214, F strict 0.8364, 48 families, 23/27 full-length.

One disclosed deviation was needed on S-REAL: the Rust arms were run with `RUSTLE_SD_READ_ISOFORM=0` so that `tx.fa` carries **one rep chain per node** (with widening on it carries one entry per admitted chain, which is not a node table); the mirror then applies the shipped k=5 widening itself. A first U40 run made without this flag failed parity (16,682 tx entries for 8,565 nodes) and was discarded, not reported.

## 2. ADDENDUM A1 — the "ts evidence" test in the task spec is vacuous, and had to be replaced

The task defined V2's unknown-strand condition partly as "whose reads carry no ts evidence". **Under `-uf`, minimap2 tags every primary alignment `ts:A:+`, including single-block alignments with no junction.** On the ideal BAM all 1,151,970 primary MAPQ>=1 reads carry a `ts` tag, 15,990 of them with no N in the CIGAR; on all three realistic arms `frac_ts` is 1.0000 for both read classes. A `ts` tag is therefore not a measurement and cannot be the test.

The operational test actually used, declared in ADDENDUM A1 before any variant metric: **a node has junction/strand evidence iff at least one of its reads has an N operation in its CIGAR.** This is also the condition the Rust itself branches on (`build_spliced_seq`: `if introns.is_empty()` -> `read_strand.or(strand).unwrap_or('+')`, with `read_strand` inert because `RUSTLE_READ_STRAND` is off by default). A second consequence, reported not acted on: the mirror's and the Rust's "flip the read strand when `ts == '-'`" branch is **dead code** on this data — read strand is exactly FLAG 0x10.

Direct confirmation of the placeholder, measured not assumed: **every single-exon-representative base node is '+' and none is '-'** — 451/451 on S-IDEAL, 3,076/3,076, 5,031/5,031 and 6,072/6,072 on U00/U20/U40.

## 3. The premise about the three blocked copies is only half right

Measured on the frozen ideal node table, the nodes that block NPIPB13, NPIPB4 and NPIPB14P are **not single-exon placeholders**:

- NPIPB13 (copy '-'): blocker is node #2027, `chr16:30595499-30631816` **'+' with 3 exons** — a spliced model whose strand came from junction motifs.
- NPIPB4 (copy '+'): blocker is node #1944, `chr16:22359779-22383315` **'-' with 2 exons** — also spliced.
- NPIPB14P (copy '-'): its best node is **already same-strand** (`chr16:75785714-75875793`, '-', 26 exons) and covers only 0.122 of the copy. Its defect is node EXTENT, not strand, and no strand variant moves it on any substrate (0.122 in every arm, every variant).
- NPIPB12: 0 same-strand primary MAPQ>=1 reads on S-IDEAL, so it stays at 0.042 in every ideal variant. **It is not an O2-only problem on real-shaped data** — see section 6.

So on S-IDEAL the strand fix has to work against *measured* opposite strands, which is exactly what makes V3 (rather than V1) the right shape there.

## 4. The mechanism, in one number pair

Of the 5,338 read-locus candidates on S-IDEAL, 68 are admitted by V0, 5,096 are blocked by a **same-strand** node (correct), and **174 are blocked only by an opposite-strand node**. Of those 174:

- **171** have at least one blocker that is an unmeasured '+' placeholder;
- **3** have blockers that are all spliced.

V1 admits all 174. V3 admits exactly the **3** (because under V2 an 'U' node is compatible with either strand and therefore keeps blocking). The 174 V1 adds are **171 single-exon nodes + 3 multi-exon nodes**: V1 is, in bulk, creating a mirror-image node beside every unstranded locus — the same double-counting shape that killed `RUSTLE_COLLAPSE_UNSTRANDED`'s sibling. The 3 nodes V3 adds are `chr16:22359821-22382362+` (8 exons, 3,640 bp — NPIPB4, exon_bp 3,634), `chr16:30609373-30634670-` (9 exons, 5,956 bp of NPIPB13's 6,144) and `chr16:15095461-15110248+`; V2's one merge replaces a '+' single-exon node and a '-' two-exon node at `chr11:19031498` with one '-' node.

On the realistic arms the same split is measurable and shifts: opposite-strand-only blocks rise to 232 / 292 / 440, of which 56 / 116 / 246 have all-spliced blockers. So V3's admission set grows from 3 (ideal) to 59 / 121 / 251 — that growth is where its realistic blast radius comes from.

## 5. The realistic substrate (S-REAL)

Same 5,542-locus set and same 38,669 transcripts as `npip_ideal`. Per locus: 20 spliced-eligible + 8 single-exon reads emitted into one pool (154,729 reads, 386 Mbp), seed 20260917. 5' truncation: p=0.45 none, else floor(L*Uniform(0,0.7)), clamped at 200 bp, 3' intact. An "unspliced" read is a read that lies inside ONE exon (3'-anchored, length Uniform{300..1500}) — a realistic Iso-Seq single-block read, not a pre-mRNA genomic span. 0.5% independent substitutions, no indels. 2% antisense. Aligned once with `-ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes -t 4` against the prebuilt splice index; the three arms are exact BAM subsets (arm membership is a pure function of read class and index), so no read was aligned twice. One foreground `gw_family_catalog` run per arm with the stubbed minimap2 (8.5 min, 4.9 GB peak RSS each).

Sanity numbers are in `sim_table`. The two that matter: the **realized** single-block primary fraction is **0.28 / 0.43 / 0.57**, well above the designed unspliced fraction, because 28% of "spliced" reads are truncated into a single exon; and 98.6-98.7% of primary alignments map back to their source locus, so the substrate is not dominated by mismapping.

The substrate is also visibly harder than the ideal one and that is the point: locus-level P falls from 0.9465 (ideal) to 0.747 / 0.617 / 0.600, and node counts rise from 5,428 to 7,207 / 8,292 / 8,565 against 701 expressed records — fragmentation, driven by truncation, is the dominant defect there.

## 6. What the realistic arm changes about the earlier conclusions

1. **The three "blocked" copies are an artefact of the idealized substrate.** On all three realistic arms NPIPB13 (1.000) and NPIPB4 (0.953-1.000) are already correct in V0. The copy that is strand-blocked on real-shaped data is **NPIPB12**: V0 gives it a wrong-strand node at frac_copy 0.042 (U00) / 0.430 (U20, U40), and V1/V3 give it a same-strand node at **1.000** on every arm. The strand fix is therefore real and it generalises — but the declared clause (i), which names only B13/B4/B14P, cannot see it. Under the strict "gain relative to this arm's own V0" reading declared in section 5 of DECLARATIONS, clause (i) FAILS for every variant on S-REAL (nothing can be gained where nothing is lost); under the state reading used in the table it PASSES trivially, including for V0. Both readings are reported; neither rescues any variant, because (ii) and (iii) decide.
2. **NPIP wrong-strand nodes drop to 0 under V1 on every substrate** (3->1 ideal, 1->0, 3->0, 4->0 realistic) and NPIP full-length rises 22->24, 19->20, 17->20, 16->19. That is a genuine, replicated gain.
3. **But the cost is now visible.** On S-REAL, V1 adds 232 / 292 / 440 nodes into a graph that is already over-fragmented, and locus-level F falls on every arm (0.8451->0.8322, 0.7574->0.7451, 0.7440->0.7271). V2's unknown-strand consolidate **merges** (removes 8 / 21 / 20 nodes, adds 3 / 5 / 5) and F rises on every arm (0.8462, 0.7610, 0.7481) — a positive that S-IDEAL, where V2 is a near no-op (one merge), completely hides.
4. **V2's merging is what loses families.** On S-REAL, V2 and V3 lose 1 / 3 / 1-2 families at best-match Jaccard < 0.5, while V1 loses 0 on U00 and U20 (and 1 on U40). No variant on any substrate loses a family's every copy span, and the largest family's best-match Jaccard never drops below 0.9032 (ideal) / 0.9851 (realistic) — so the catastrophic shape of the old `RUSTLE_COLLAPSE_UNSTRANDED` kill (34 of 121 families lost, 25 losing every copy span) does **not** reproduce here. The declared V2 hazard, '+' and '-' pieces chained through a 'U' piece, never fired: `strand_bridged_groups` is 0 on all four substrates.
5. **FAMILY R is uninformative on S-REAL**: it is 1.0000 for V0 on all three arms (all 27 NPIP copies already in one family), so clause (iv) passes by tie everywhere. FAMILY P strict falls for V1/V3 on every realistic arm (0.6000->0.5870, 0.3375->0.3253, 0.4030->0.3803). FAMILY P NPIP-only is 1.000 in every single row and was used in no decision, per the metric-traps register.
6. **k could not be measured here either**, and was not claimed: the widening counts are 3,896 (ideal) and 3,409 / 4,015 / 4,243 (realistic) at the shipped k=5, with no sweep run.

## 7. Clause-by-clause verdicts

See `strand_table`. Summary: **V3 is a CANDIDATE on S-IDEAL and on no realistic arm. V1 and V2 are candidates nowhere.** Since the decision rule requires both substrates, **no variant qualifies.**

## 8. Anti-trap columns (reported as traps, not targets)

A1 (nodes overlapping no expressed record) and A2 (nodes overlapping >=2 records by >=50 bp) are reported without being scored as progress. On S-IDEAL A1 is 31 for every variant and A2 moves 63->64 under V1/V3; the annotated ceiling's own A2 is 89, so moving A2 toward 89 is not progress. On S-REAL A1 is the one place V2 looks clearly good — 49->40 (U40), 45->35 (U20), 34->32 (U00) — i.e. the unknown-strand merge is deleting nodes that overlap no expressed record at all. That is suggestive, not a pass: it is a trap column and the same merge is what loses families.

## 9. Scope limits of this run

- Every variant number is a mirror number; the Rust was used only to produce the rep/node set, with widening off, and parity was asserted first.
- The family-level blast radius is measured on substrate (a) (chr16/17/18), not genome-wide; the node-level panel is genome-wide.
- The expressed-record set used for the locus metric is held fixed at the frozen ideal substrate-(a) set (701 records) on both substrates, so the two are comparable; it is not re-derived from S-REAL expression.
- Merged nodes' exon unions on S-REAL are approximated by the representative chain (24 / 31 / 35 nodes per arm, <0.5%), as declared.
- S-REAL shares its annotation with its truth exactly as S-IDEAL does, so it is not an independent test of node-vs-annotation fidelity; what it adds is read SHAPE (truncation, single-block reads, errors, antisense), which is the variable under test.
## Verification (independent recompute) — agent 2 of 2

My declarations were written before any verifier number: `/mnt/linuxdisk/home/juanfraitu/strand_fix/verify/vdecl.txt` (md5 `5307c86ff3edca945a1b67d75f07270d`, 2026-09-18T00:52:08-07:00). I read no script under `strand_fix/scripts/`; my mirror of `shared_definition.rs` (`overlap_groups`, `consolidate`, `read_blocks`, `depth2_exons`, `with_read_locus_nodes`, `widen_with_read_isoforms`, `edges`, `gene_body_chains`, `triangle_leaders`) was written from `/mnt/c/Users/jfris/Desktop/Rustle/src/rustle/vg_family/shared_definition.rs` at commit 94ed994e. Code and outputs: `/mnt/linuxdisk/home/juanfraitu/strand_fix/verify/` (`vmirror.py`, `vmirror_w.py`, `v1_rederive.py`, `v3_v2v3.py`, `v5_ideal.py`, `v7_copies.py`, `v8_locus.py`, `v9_queries.py`, `v10_family.py`, `v11b_blast.py`, `v4_bamsanity.py`). Nothing in `src/` was touched, nothing committed, no MAPQ gate changed, no index rebuilt (`npip_ladder/idx/target.{splice,asm20}.mmi` reused), no subagents.

**1. The code claims are true.** `overlap_groups` is keyed `(chrom, strand)` and `consolidate` groups pieces on that key (shared_definition.rs:130, 177); `ExonIndex` carries only `chrom` (:333) and `with_read_locus_nodes` suppresses a candidate whenever `bidx.hits(...)` is non-empty on any strand (:366); widening is strand-strict (`nodes[v].strand == r.strand`, :~500). The placeholder is `denovo_assemble.rs` `build_spliced_seq_with`: `let strand = if introns.is_empty() { read_strand.or(strand).unwrap_or('+') }` — with the in-source note that all 5,928 single-exon reps in the shipped dump are '+' and none '-'.

**2. V0 parity with Rust, and node-for-node reproduction of all four variants.** The Rust arms ran with widening off, so `rust/<arm>/tx.fa` headers are the exact rep chains of the Rust node set. All 7,207 / 8,292 / 8,565 tx keys match the mirror's base nodes exactly (0 only-Rust, 0 only-mirror), every rep exon-sum matches the tx sequence length, and the read-locus counts match the Rust log lines (68 IDEAL, 56 U00, 25 U20, 22 U40). I then re-derived every variant from the BAMs myself:

| arm | V0 | V1 | V2 (base) | V3 | V4 (mine, unrun by builder) |
|---|---|---|---|---|---|
| S-IDEAL | 5428 / +68 / 3896 widened | 5602 / +242 / 3897 | 5427 / +68 (base 5359) | 5430 / +71 / 3897 | 5431 / +71 / 3897 |
| U00 | 7207 / +56 / 4243 | 7439 / +288 / 4263 | 7202 / +56 (base 7146) | 7258 / +112 / 4263 | 7263 / +112 / 4263 |
| U20 | 8292 / +25 / 4015 | 8584 / +317 / 4080 | 8276 / +25 (base 8251) | 8392 / +141 / 4080 | 8408 / +141 / 4080 |
| U40 | 8565 / +22 / 3409 | 9005 / +462 / 3530 | 8550 / +22 (base 8528) | 8796 / +268 / 3529 | 8811 / +268 / 3529 |

Every cell equals `out/variants_*.json`. Node sets are identical key-for-key including widened `tx_chains` (V0/V1: exact; V2/V3: exact except that the builder emits the leftover unknown-strand nodes as '+' — 424 / 480 / 503 / 568, matching its own `n_U_left_placeholder_plus`). Evidence sets reproduce: single-exon-rep base nodes 451 / 3,076 / 5,031 / 6,072, **100 % of them claiming '+' in every arm**; U sets 425 / 485 / 519 / 583; strand-bridged groups 0 everywhere.

**3. "V1/V2 change only what they claim" — confirmed.** V1 removes nothing and alters no V0 node (V0-nodes-absent = 0 on all four arms); it only adds 174 / 232 / 292 / 440 nodes. V2's only structural change is the U-merge (removes 2 / 8 / 21 / 20, adds 1 / 3 / 5 / 5) plus strand labels. These are exactly the builder's blast-radius node rows. One fact worth putting on the table: on S-IDEAL only 3 of V1's 174 additions are exon-sets that did not already exist in V0 — 171 are the *same locus a second time on the opposite strand* (U00 56 new of 232; U20 116 of 292; U40 244 of 440). That is why V1 costs precision.

**4. Metrics recomputed from scratch.** 46-copy table (truth = NPIP/TBC1D3 exon unions from `npip_ideal/locus_set.tsv` + `transcripts.tsv`) reproduces every published cell: NPIP full-length 22→24 (IDEAL), 19→20 (U00), 17→20 (U20), 16→19 (U40); wrong-strand nodes 3→1, 1→0, 3→0, 4→0; all-46 41→43, 37→38, 34→37, 31→34. Per copy: NPIPB13 0.028→0.944, NPIPB4 0.033→1.000, NPIPB14P 0.122 in all variants and all arms, NPIPB12 0.0423 (IDEAL, U00 V0) / 0.4303 (U20, U40 V0) → 1.000 under V1/V3/V4 on the realistic arms and unchanged on S-IDEAL. Locus-level bipartite on substrate (a) (frozen `locus_width/out/sub_a.pkl`, 701 records, max-weight matching) reproduces R, P, F, A1 and A2 for all sixteen arm×variant cells to 4 decimals. Family panel: I rebuilt the chr16-18 query fastas myself (sequences byte-identical to Rust's for all 913 shared keys), ran my own minimap2 (`-c -N 50 -p 0.1 -x splice -uf`, 94.9 s, 15.2 GB; `-x asm20`, 70.7 s, 12.9 GB) against the prebuilt indexes, and ported `edges`/`gene_body_chains`/`triangle_leaders`: U00 gives 50 / 50 / 49 / 49 families for V0 / V1 / V2 / V3, families lost 0 / 0 / 1 / 1, mean best Jaccard 1.0000 / 0.9960 / 0.9733 / 0.9733, largest family 45 at Jaccard 1.0000 — an exact match to `out/blast_U00.json`. The lost family is concrete: a 2-node family of single-exon '+' nodes, chr16:21,605,134-21,606,680 and chr16:22,888,332-22,889,842; under V2 the second is U-merged with an overlapping '-' node into chr16:22,888,332-22,891,320 ('-', 2 exons) and the pair edge dies. That is the prior kill (`project_engulfment_strand_placeholder`) reappearing, and it comes from V2's merging half, exactly as the builder says.

**5. The realistic reads hold their declared parameters,** measured from the BAM and manifest, not from the generator: substitution rate X/(=+X) = 0.00503 (declared 0.005), indel rate 0.00065 (none simulated — alignment artefact), antisense 0.0198 (0.02), 45.50 % of spliced reads untruncated (0.45), 5' truncation fraction f ≈ Uniform(0, 0.7) (median 0.353, max 0.700), single-exon reads 3'-anchored and capped at 1,500 bp, arm composition exactly as declared (U00 = SP indices 0-19; U20 = SP 0-15 + SE 0-3; U40 = SP 0-11 + SE 0-7, i.e. subsets of one alignment), primary/unmapped 110,365/150, 110,359/162, 110,348/178, maps-back-to-source 0.9864 / 0.9865 / 0.9782 (exact), realized single-block primary fraction 0.2811. I could not re-run the generator (its script is builder code, off-limits), so seed reproducibility was replaced by these conformance checks plus the subset-structure check; in that sense the seed claim is *unverified*, not *verified*.

**6. Does the realistic substrate reproduce the placeholder mechanism? Yes.** 3,076 of 7,151 base nodes (43.0 %), 5,031 of 8,267 (60.9 %) and 6,072 of 8,543 (71.1 %) have a single-exon representative and every one of them claims '+', against 451 of 5,360 (8.4 %) on the idealized substrate. The sim design is right on this axis. Two caveats: (a) the "0 % unspliced" arm is not a control — 43 % of its nodes are already placeholders, because 28 % of the *spliced* reads align as a single block after 5' truncation, so the sweep only spans 43 → 71 % and never a placeholder-free state; (b) the copies the fix was designed for (NPIPB13, NPIPB4) are not broken on S-REAL at all — the realistic arms exercise a different instance (NPIPB12). The two substrates therefore test the same mechanism at different loci, which is a strength for generalisation but means the pre-registered clause-(i) copy list was substrate-specific.

**7. Does the recommendation follow?** Yes, with one amendment. "Adopt none of V1/V2/V3 now" follows from numbers I reproduced independently: V1 fails (iii) on every arm, V2 fails (ii) on every realistic arm and fixes nothing, V3 fails both. "Put the no-reads-on-the-correct-strand issue into the simulations" is supported and is the single change that flips the verdict — I confirm the flip is real and not an artefact of the mirror. What does **not** follow is the closing step as stated: I built V4 and it fails clause (iii) on all three realistic arms (0.8430 / 0.7537 / 0.7357 vs V0 0.8451 / 0.7574 / 0.7440) while passing (i), (ii) and (iv) on U00. So V4 is the right *mechanism* to pre-register — it keeps every per-copy gain of V1/V3 with a V0-identical family panel — but it will be rejected by the rule as registered. Either clause (iii) is re-registered so that a correct second node at an already-matched locus is not counted purely as a precision loss, or V4 should be pre-registered explicitly as a change that trades locus precision for copy recovery, with the size of that trade (-0.002 to -0.008 F) declared in advance.

Not independently recomputed, and therefore not endorsed by me: M-C FAMILY R / FAMILY P strict on the frozen NPIP scorer (I verified family *construction*, not the NPIP-truth scorer), the family-level blast radius on S-IDEAL, U20 and U40 (I did U00 only), and the read generator's seed reproducibility.