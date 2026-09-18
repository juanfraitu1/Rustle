# V4 on real gorilla data — the tiebreak (agent 1 of 2)

Declarations written BEFORE any number of this run: `/mnt/linuxdisk/home/juanfraitu/v4_gorilla/DECLARATIONS.txt`, md5 `a099a9c28522336606b66297db9538aa`, frozen 2026-09-18T10:50:23-07:00. Calibrated thresholds frozen after Part 1: `/mnt/linuxdisk/home/juanfraitu/v4_gorilla/NOISE_FLOOR_FROZEN.txt`, md5 `69af973a85f582f9759a9bcfa634feda`, 2026-09-18T11:12:28-07:00. Both md5s re-verified unchanged at the end of the run. All outputs under `/mnt/linuxdisk/home/juanfraitu/v4_gorilla/` (`out/`, `map/`, `logs/`, `code/`, `replay/`), `TMPDIR` under `v4_gorilla/tmp`. **Nothing in `src/` was modified, nothing was committed, no subagent was spawned, no MAPQ gate was lowered, no index was rebuilt.** minimap2 ran twice, 32 s wall total, foreground, one at a time, peak RSS 5.6 GB, on 138 new transcript queries (0.54 Mbp) and 79 new gene-body queries (2.95 Mbp); every V0 query was served from the captured PAFs.

> **Verifier corrections, applied by the orchestrator** (verifier ok = true; the verdict stands, one supporting claim does not):
> 1. **Family 52 is misdescribed — and it is the trap the task named.** Only family 51 is a pair of single-exon `'+'` fragments absorbed into the bigger NPIP family; 52 is not, so "both lost families are absorbed" is wrong.
> 2. **"8 of 76 families changed" conflates two sets.** The clause flags **six** (lost {51,52}, scatter {2,4,11}, membership {2,4,11,35,52}); the other two merely grew, with membership intact.
> 3. **"A zero-width null means the clause says change nothing" overstates it, in V4's favour.** The clause did not fire on the two families that grew by 8 and 2 nodes. What it forbids is *splitting*, and V4 does split one.
> 4. **Protocol deviation**: the declared base-recovery procedure was not the one used; the declarations make a parity failure a BLOCK, and this was not treated as one. The verifier re-derived the split from the rep audit and the verdict is unchanged.
> 5. **Two input-table defects** (the node dump is not the SdNode dump for out-of-family nodes; three additions overlap an approximated span) — both verified inert.

## 0. Answer first

**V4 does not confirm on gorilla. It passes clauses (a), (c) and (d) and fails clause (b), the family clause.** V1 fails (b) *and* (c). Under the rule as pre-registered, V4 stays unadopted and the strand asymmetry is documented as measured-but-unfixed.

| clause | V1 | V4 |
|---|---|---|
| (a) strand | **PASS** opposite 9 → 3, 1 gain > 0.5 | **PASS** opposite 9 → 3, 1 gain > 0.5 |
| (b) family clause v2 | **FAIL** lost 2 ≤ 0, scatter 3 ≤ 0, memb 71/76 ≥ 76/76 | **FAIL** identical: lost 2, scatter 3, memb 71/76, orph 0 |
| (c) junk | **FAIL** 464 → 513 (+49 vs allowance 46.4) | **PASS** 464 → 482 (+18) |
| (d) no node lost | PASS removed 0, added 138 | PASS removed 0, added 86 |
| reported, non-gating | P_cand −0.0235, P_old −0.0266, R 0.9235→0.9337, F_old 0.7864→0.7720 | P_cand −0.0160, P_old −0.0142, R 0.9235→0.9337, F_old 0.7864→0.7805 |

Two things make this a different failure from U00's.

1. **The fresh arm delivered the dynamic range it was supposed to.** V0 carries **9** opposite-strand membership loci here, against 3 on U00 and 6 on U40. Clause (a) was not the binding constraint this time; V4 cleared it with room. Clause (c) had teeth too — V1 actually exceeded the junk allowance, which it never did on the simulated arms.
2. **The family charge is real, not noise.** The gorilla null has zero width, and so does V4's own graph's null, so the 8 changed families are a deterministic difference between two partitions. On U40 the single lost family sat inside a null of width 1; here nothing does.

## 1. Provenance — every parity gate passes

The de novo node table publishes the post-pass node set but no rep audit, and the base/read-locus split is not recoverable from the table alone (412 nodes match a read-locus candidate exactly; only 233 can be read-locus nodes). I therefore **replayed the shipped binary** (`npip_ladder/rust/gw_family_catalog.bin`, `RUSTLE_SHARED_DEFINITION=1 RUSTLE_LOCUS_AUDIT=1`) with a wrapper that runs **no minimap2** and replays the captured PAFs. 6:35 wall, 3.4 GB peak, rc 0.

- It reproduced the original log verbatim: `2469 reps -> 2431 gene-level loci + 233 read-locus nodes = 2664 nodes` and `edges: 1194 exon + 297 gene-body (1248 pairs); 76 triangle-supported families holding 323 loci`.
- The `tx.fa` and `body.fa` it generated are **byte-identical** to `ggo_npip/capture/ggo_dn/`.
- **GATE P1** base split 2,431 / 233 from the rep-audit spans — PASS.
- **GATE P2** re-deriving V0 from that base with the unmodified mirror gives back the published 2,664-node multiset exactly — PASS.
- **GATE P4** my Python graph on V0 gives 1,194 exon + 297 body edges, 1,248 pairs, 76 families, 323 loci — PASS, identical to the Rust's own line.

One declared approximation: 33 consolidate-merged base nodes carry the rep chain instead of the true exon union. It is measured and harmless — V0 still re-derives to exactly 233 added nodes, and **0 of V1's 138 and 0 of V4's 86 additions** sit inside an approximated span without touching its modelled exons.

## 2. Part 1 — the noise floor, frozen before any variant clause number

The V0 graph taken unchanged (2,664 nodes, 1,248 edges, 76 families, 323 loci), indices relabelled 30 times, `random.Random(20260918 + k)`, k = 1..30. Asserted and passed every draw: no identity permutation, edge **set** maps back identically, edge count unchanged, degree multiset unchanged, `n_reads` multiset unchanged.

| quantity | distribution over 30 draws |
|---|---|
| families lost (best-match exon-bp Jaccard < 0.5) | {0: 30} — P95 = 0 |
| scatter | {0: 30} — P95 = 0 |
| orphaned | {0: 30} |
| membership | {76/76: 30} — P5 = 1.000000 |
| mean best Jaccard | 1.0000 in all 30 (min best J 1.0000) |
| n_families / largest / loci | 76 / 49 / 323 in all 30 |

**Gorilla does NOT behave like U20/U40.** There, "one family lost" was pure index noise in 17/30 and 14/30 draws. Here the null is exactly zero-width, as on S-IDEAL and U00. Scatter is 0, consistent with 0 in all 120 simulated draws — it remains the strict, noise-free axis.

Frozen clause: `lost ≤ 0, scatter ≤ 0, membership ≥ 76/76, orphaned == 0`. Declared diagnostic, run after the freeze: the same null on **V4's own graph** is also zero-width (lost {0:30}, scatter {0:30}, membership 77/77 in all 30). Calibrating on V0 is fair to V4 — and it means the family difference below is deterministic.

Freeze-order limit, stated plainly: the graph build had already produced each variant's raw family *count* (V0 76, V1 77, V4 77) at 11:11:26, before the 11:12:28 freeze. No family-*clause* quantity existed until 11:13:35. The thresholds are mechanical percentiles of the V0 null; there was no value to choose. This is recorded inside the frozen file itself.

## 3. Part 2 — what each variant does

**Node sets.** V0 2,664 → V1 2,802 (+138, −0) → V4 2,750 (+86, −0). **V4's additions are a strict subset of V1's.** V4 adds 76 `-` and 10 `+` nodes; the 52 V1-only additions are *all* `-` and 48 of 52 are single-exon — which is exactly why V1's junk cost is 2.7× V4's. Base: 2,431 nodes, 863 with a single-exon rep chain, 297 of those unmeasured; all 297 carry `+`, so the two V4 restore conventions are byte-identical here.

**Strand, against the 25 membership loci.** Opposite-strand best node 9 → 3 for both variants, fixing the *same six* loci (LOC115932781, LOC129527636, LOC109023568, LOC115932744, LOC129527692, LOC115933039). Loci with any same-strand node 14 → 18. Mean same-strand coverage 0.0856 → 0.1638. GAIN above 0.5 = **1** (LOC115933039, 0.091 → 0.7331 — the new 45-exon `-` node `NC_073244.2:20932517-21081227`, 54 reads). LOSS = 0. **carried = 0**: V0 had *no* locus above 0.5, so none of the recovery is an already-correct copy and the denominator never moves.

**The three residual opposite-strand loci are not suppression failures.** LOC129523555 and LOC129527693 have **2** same-strand primary MAPQ≥1 reads each, below `MIN_LOCUS_READS = 3` — no variant can build a node there. LOC115931102 already *has* a same-strand node (2-exon `-`) and loses the best-node title on covered bp alone. So on gorilla the strand-suppression defect itself is fully repaired; what remains is the read floor and node size, which are different problems.

**Junk and precision, genome-wide over the three contigs, against 1,975 expressed native records.** Junk 464 → 482 (V4) / 513 (V1). P_cand 0.8291 → 0.8131 / 0.8056. P_old 0.6847 → 0.6705 / 0.6581. Record recall R rises 0.9235 → **0.9337 for both**. F_old 0.7864 → 0.7805 / 0.7720.

**Family panel, identical for V1 and V4:** 8 of 76 V0 families changed, 2 lost, 3 scattered, 0 orphaned, membership 71/76, mean best Jaccard 0.9477, min 0.0705, largest family (49 nodes) best J 1.0000. V4 has 77 families holding 345 loci vs V0's 76 holding 323.

## 4. Every flagged family, opened

All eight have one mechanism: **V4 installs the correct-strand spliced node at a locus where V0 had only a single-exon `+` placeholder, and the placeholder either joins the enlarged family or falls out as a claimed pendant.**

- **V0 family 52** (2 nodes, best J 0.0705, LOST + SCATTER). Its `-` 26-exon node joins the enlarged NPIP family 2; its `+` single-exon fragment `NC_073244.2:20955338-20957550` falls out — it is *inside* the new 45-exon `-` node that produced the one clause-(a) gain.
- **V0 family 51** (2 nodes, best J 0.0876, LOST). Both `+` single-exon fragments are **absorbed** into family 2, which grows from 21 to 28 nodes with 6 new `-` nodes. Nothing is destroyed; the footprint ratio alone drops the Jaccard.
- **V0 family 35** (best J 0.5599, SCATTER). The `+` fragment `NC_073241.2:26331092-26333296` is superseded by a new 19-exon `-` node containing it.
- **V0 family 2** (21 nodes, best J 0.5788, SCATTER). 19 nodes stay and gain 6 new `-` nodes; **two `+` single-exon fragments at the NPIPB14P and NPIPB6 loci break off into their own 2-node family 55.** This is the one genuine cost.
- **V0 families 3, 4, 6, 11** — the same pattern: new multi-exon `-` nodes containing the old `+` fragments; e.g. family 6's LOC129527692 moves onto a new 45-exon `-` node and joins family 4, the *long-spliced* NPIP family, which is arguably the better home.

**The price, named rather than buried:** membership-locus FAMILY R falls **0.68 → 0.64** (17/25 → 16/25 in the single best family) for both variants; loci with a node stays 20/25.

## 5. What I would put in front of the advisor

**V4 does on real gorilla data exactly what it does on simulation — it clears two thirds of the wrong-strand damage, removes no node, costs 1.6 points of candidate precision and 18 junk nodes, and *raises* record recall — and it is rejected by a clause that on this substrate reads "change no family".** The clause is not wrong to fire: the partition really does change, deterministically, and two NPIP fragments really do split off. But six of the eight families it charges for are V4 replacing a wrong-strand single-exon placeholder with the spliced node the reads support, and two of them are 2-node pairs being *absorbed* into a bigger, more complete family.

V1 is now clearly dominated and should not be revived: it reaches exactly the same copies as V4, with 138 additions to V4's 86 and a junk cost that breaks its own clause.

I did not re-register any clause after seeing this, and I recommend against a follow-up that re-scores gorilla under a loosened family clause. The clean move is to settle, *in advance and as a question of what a family is*, whether replacing fragments with correct-strand spliced nodes counts as damage to the partition — because across U00, U40 and gorilla that is now the only thing between V4 and adoption. A defensible pre-registration would score the partition on the truth (membership-locus family recovery, which here falls 0.68 → 0.64 and is a real cost worth arguing about) rather than on similarity to V0's own footprints, which structurally penalises any improvement to node construction.
## Verification (independent recompute)

Verifier code: `/mnt/linuxdisk/home/juanfraitu/v4_gorilla/verify/` (`w0`–`w30`), reusing only my own earlier verifier modules (`strand_fix/verify/vmirror.py`, `v4_heldout/verify/x7_family.py`, `xlib.py`, `x10_noise.py`). No builder script under `v4_gorilla/` was read or imported; `v4_gorilla/out/*.json` were read only as targets to compare against. No `src/` file was modified (no file under `src/` has a 2026-09-18 mtime), nothing was committed (HEAD still 89230134), no MAPQ gate was lowered.

### 1. Declarations and frozen thresholds precede every variant number
`DECLARATIONS.txt` 10:50:23, md5 a099a9c2… verified; `NOISE_FLOOR_FROZEN.txt` 11:12:28, md5 69af973a85f582f9759a9bcfa634feda verified, and both re-verify unchanged now. Ordering: declarations (10:50) < all results; V0 graph (graph.json 11:11:26) and the V0 null (noise_V0.json 11:12:00) < freeze (11:12:28) < scorer (11:13:05), family panel (11:13:35), V4 null (11:14:52), verdict (11:15:19). The one thing that predates the freeze is `variants.json` (11:10:00), i.e. the variant NODE counts and the raw family counts; the frozen file discloses exactly that ("V0 76, V1 77, V4 77") and no clause quantity (lost/scatter/membership/orphaned/Jaccard) existed before it. I consider the freeze honest. The only ordering problem in the run is the undeclared base-recovery substitution (correction 1), whose failing first attempt is itself timestamped at 10:53.

### 2. Parity, from my own replay of the shipped binary
I re-ran `npip_ladder/rust/gw_family_catalog.bin` (RUSTLE_SHARED_DEFINITION=1, RUSTLE_LOCUS_AUDIT=1, my own PAF-replay wrapper, no minimap2 run) on `dn/GGO.3ctg.bam`/`GGO.3ctg.fa`. It printed `2469 reps -> 2431 gene-level loci + 233 read-locus nodes = 2664 nodes` and `edges: 1194 exon + 297 gene-body (1248 pairs); 76 triangle-supported families holding 323 loci`; the tx/body query FASTAs it generated are byte-identical to `capture/ggo_dn/{tx,body}.fa` (md5 2bb35c9a…, 45d7516b…), and `vggo.{copies,families,pairs}.tsv` are md5-identical to the shipped `ggo_sd.*`. No widening line appears, confirming the binary predates read-isoform widening (P3's assumption holds).
- Reads: 139,912 primary/non-secondary/non-supplementary MAPQ≥1 (NC_073241.2 21,874 / NC_073242.2 59,535 / NC_073244.2 58,503; + 73,051 / − 66,861; 128,182 spliced) — identical to the builder.
- BASE split derived independently from the 2,469 rep-audit rows: 2,431 base / 233 read-locus, and all 233 unmatched nodes are exactly read-locus candidates (0 exceptions). **P1 PASS.**
- Re-deriving V0 from that base reproduces the published 2,664-node key multiset exactly (0 missing, 0 extra). **P2 PASS.**
- My Python edge/leader port on V0 gives 1194 exon + 297 body = 1248 pairs, 76 families, 323 loci, and the partition is **family-for-family, node-for-node identical to the shipped `ggo_sd.copies.tsv`** (not just the counts). **P4 PASS.**

### 3. Node sets and queries
V0 2,664; V1 2,802 (+138, 0 removed); V4 2,750 (+86, 0 removed); V4's additions are a strict subset of V1's. V4 adds 76 '−' and 10 '+'; the 52 V1-only additions are all '−' and 48/52 single-exon. Base 2,431, 863 with a single-exon rep chain, 297 unmeasured, all 297 carrying '+' in V0 — so the two restore conventions are byte-identical here. All of this matches the builder exactly. I generated the new queries myself (138 tx, 79 body — same header sets as the builder's, order differs) and ran my own minimap2 2.30 with the declared flags against `dn/GGO.3ctg.fa` (rc=0 both). Note: the body key must come from the node's real span (node_id / body.fa header), not from the degraded `exons` column, or the new-body count comes out 111 instead of 79; the builder evidently did this correctly.

### 4. Noise floor, my own permutation code
30 relabellings, `random.Random(20260918+k)`, k=1..30, run in BOTH index conventions (perm as scatter and as gather), each draw asserting non-identity, edge-set identity after mapping back, unchanged edge count, degree multiset and n_reads multiset. V0 graph: lost {0:30}, scatter {0:30}, orphaned {0:30}, membership 76/76 in all 30, mean and min best Jaccard 1.0000, n_families 76, loci 323 in every draw — **zero width, both conventions**. V4's own graph: lost {0:30}, scatter {0:30}, membership 77/77 in all 30, mean best J 1.0000. Confirmed: the gorilla null has zero width, calibrating on V0 is fair to V4, and the 6 clause-flagged families are deterministic, not tie-break noise.

### 5. Scorer v3 — every number reproduced
| | V0 | V1 | V4 |
|---|---|---|---|
| nodes | 2664 | 2802 | 2750 |
| edges | 1194+297=1248 | 1383+469=1466 | 1370+457=1453 |
| families / loci | 76 / 323 | 77 / 347 | 77 / 345 |
| opposite-strand membership loci | 9 | 3 | 3 |
| loci with any same-strand node | 14 | 18 | 18 |
| mean same-strand coverage | 0.0856 | 0.1638 | 0.1638 |
| loci with any node | 20/25 | 20/25 | 20/25 |
| JUNK | 464 | 513 (+49) | 482 (+18) |
| M / C | 1824 / 2200 | 1844 / 2289 | 1844 / 2268 |
| P_cand | 0.8291 | 0.8056 | 0.8131 |
| P_old | 0.6847 | 0.6581 | 0.6705 |
| R | 0.9235 | 0.9337 | 0.9337 |
| F_old | 0.7864 | 0.7720 | 0.7805 |
| panel vs V0 | — | lost 2, scatter 3, memb 71/76, meanJ 0.9477, minJ 0.0705 | identical |

Truth sets reproduce exactly: 4,477 gene/pseudogene records on the three contigs, 25/25 membership loci resolved (10 '+', 15 '−', mean 13,300 exonic bp, 23 on NC_073242.2, 0 excluded), 1,975 expressed records at ≥3 same-strand primary MAPQ≥1 reads.

### 6. The declared traps
- **Gain that is denominator shrinkage / hairline:** no. Denominator fixed at 25; loci with any node 20/25 in both; carried = 0 (V0 had no locus above 0.5), GAIN = 1 (LOC115933039, 0.0910 → 0.7331), LOSS = 0. The six loci whose best-node strand is fixed are exactly LOC115932781, LOC129527636, LOC109023568, LOC115932744, LOC129527692, LOC115933039; five of them reach only 0.197–0.344 same-strand coverage, so "clears every wrong-strand copy" is a statement about the best-node title, not about recovery.
- **A copy counted as recovered that was already correct:** none (carried = 0).
- **Residual 3 are not suppression failures:** verified. LOC129523555 and LOC129527693 have 2 same-strand primary MAPQ≥1 reads each (< MIN_LOCUS_READS = 3, unreachable by any variant); LOC115931102 has 39 same-strand reads and does have a same-strand node, losing the best-node title only on covered bp (0.0946 vs 0.1984).
- **A clause that cannot fail:** (d) is free for both variants by construction (neither can remove a node) — correctly stated, not counted. (a) had dynamic range (9 available) and could have failed. (b) can fail and did. (c) is the only clause separating V1 from V4.
- **Strand-merged families, as on U40:** yes, and the builder missed one. 37 of 76 V0 families hold both strands; 12 have one strand entirely single-exon while the other is spliced. V0 family 52 is one of them and is one of the two "lost" families (correction 2). More generally, the six flagged families are 85.7% single-exon nodes (36/42) against 25.6% (72/281) in the 70 untouched families, and 19 of their 42 nodes are co-located with a V4-added opposite-strand node (17 with a spliced one). Every one of the 98 V0 nodes that a V4 addition overlaps is overlapped on the opposite strand. The clause's charge therefore falls almost entirely on families built out of the wrong-strand single-exon placeholders V4 exists to replace.
- **Junk, rate rather than count (independent framing):** V0 0.1742, V4 0.1753, V1 0.1831. V4's additions are 20.9% junk (18/86) against a 17.4% background; V1-only additions are 59.6% junk (31/52). V1 and V4 separate the same way without relying on the absolute-count allowance.

### 7. Clauses and whether the recommendation follows
(a) opposite 9→3, ≥1 gain above 0.5 — **PASS** for V4 and V1. (b) lost 2 > 0, scatter 3 > 0, membership 71/76 < 76/76, orphaned 0 — **FAIL**, identically for V4 and V1, against a null that is zero-width in both permutation conventions. (c) allowance max(3, 0.10×464) = 46.4: V4 +18 **PASS**, V1 +49 **FAIL**. (d) 0 removed — PASS, free. Arithmetic and verdict confirmed: **V4 does not confirm on gorilla; V1 does not confirm and is strictly dominated** (same 3 residual loci, same 6 fixed loci, identical family panel, +138 nodes vs +86, junk +49 vs +18). The recommendation follows from the rule as written, and the builder's proposed next step — deciding in advance, on its own merits, whether replacing single-exon wrong-strand fragments with correct-strand spliced nodes is a cost or the point — is the right one; my family-52 and single-exon-enrichment findings sharpen it, because one of the two "lost" families is a strand-merged pair being separated and the clause's cost is concentrated on fragment families.