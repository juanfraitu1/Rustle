# SCORER v3 + genome-wide family panel: V1 / V3 / V4 on S-IDEAL, U00, U20, U40 (agent 1 of 2)

Declarations written BEFORE any number of this run: `/mnt/linuxdisk/home/juanfraitu/strand_v4/DECLARATIONS.txt`, md5 `07148b973ae1100db1fe69e9385c13c2`, frozen at 2026-09-18T08:34:32-07:00 (md5 re-verified unchanged at the end of the run). All outputs under `/mnt/linuxdisk/home/juanfraitu/strand_v4/` (`out/`, `map/`, `logs/`), `TMPDIR` under `strand_v4/tmp`. minimap2 2.30-r1287, samtools 1.22.1. **Nothing in `src/` was modified, nothing was committed, no subagent was spawned, no MAPQ gate was lowered, no index was rebuilt** (`npip_ladder/idx/target.{splice,asm20}.mmi` reused). 20 minimap2 runs, 23.7 min wall total, peak RSS 21.6 GB, one at a time, foreground.

> **Verifier corrections, applied by the orchestrator** (verifier ok = true; the headline changes materially):
> 1. **The S-IDEAL "lost family" is an improvement scored as a loss.** The four nodes carry 5 of the 6 possible edges; V0 splits that near-clique across three families and V4 unites it. Jaccard falls to 0.4933 purely because the family grew.
> 2. **The U20 failure is below the leader rule's own noise floor.** Re-running triangle leaders on the UNCHANGED V0 graph with node indices relabelled — identical nodes, edges, read counts and degrees — produces family changes of the same magnitude. Clause (ii) cannot distinguish V4's effect from index-order noise on that substrate.
> 3. **The new precision does not repair the trap it was written for** — it slightly aggravates it, because nodes added at real loci enter the candidate denominator unmatched. Junk exclusion itself is exact (asserted in code, all 8 cells).
> 4. **Clause (v) is vacuous on U20** (FAMILY R = 1.0000 for every arm) and **clause (iii) is a 1-node test on S-IDEAL** (5 junk nodes, so a 10% bar means "add exactly zero").
> 5. **Two V3 membership cells do not reproduce** under any natural key convention (250/252 and 249/264 become 251/252 and 258/264, or 229/252 and 232/264, depending on the 'U' key).
> 6. **The proposed membership clause needs a threshold before it is pre-registered**, calibrated against the measured noise floor rather than set at 100%.
> 7. Node tables were reused from the previously verified mirror pickles rather than re-derived from the BAMs; V0 parity against the Rust capture was re-asserted.

## 0. Answer first

**V4 fails the pre-registered rule, and it fails only clause (ii), by one family out of 252 on S-IDEAL and one out of 264 on S-REAL U20.** It passes (i), (iii), (iv) and (v) on both substrates. Both clause-(ii) failures were opened up and neither is a destroyed family: on S-IDEAL the "lost" family is intact and merely absorbed two newly admitted same-strand copies (Jaccard falls because the family GREW); on U20 one degree-1 node keeps its only edge but is left unassigned by the greedy triangle-leader rule. The rule as registered still says: do not adopt.

## 1. Instrument defects this run was built to fix

**(1) Locus precision could not credit a strand fix.** SCORER v3 adds `P_cand` = matched / (nodes overlapping >= 1 expressed record), reported beside the old `P_old` = matched / all nodes and labelled. The honest result is that the trap fix does **not** rescue the variant that the old instrument punished hardest: V1's new nodes are mirror-image copies at real loci, so they all enter the candidate denominator and cannot be matched — `P_cand` charges V1 *more* than `P_old` did on S-IDEAL (-0.0282 vs -0.0257). For V4 the two conventions agree to 0.0005 (-0.0014 IDEAL, -0.0048 U20). So the previous round's "clause (iii) is a trap" objection was correct in principle, and removing the trap changes no verdict.

**(2) The family blast radius was 48-57 families of chr16/17/18.** It is now **genome-wide: 252 families (S-IDEAL) and 264 families (U20)**, built with the shipped two-run procedure over every node on every chromosome — larger than the 121-family panel on which the old `RUSTLE_COLLAPSE_UNSTRANDED` kill was measured. Query volume: 37,969 tx / 5,428 body keys (177 + 359 Mbp) for S-IDEAL; 21,427 tx / 8,396 body keys (62 + 260 Mbp) for U20. 617 + 196 (IDEAL) and 2,496 + 803 (U20) md5 query keys were served from the captured PAFs under `strand_fix/map/`, `npip_ladder/union/` and `npip_ideal/`; everything else was mapped fresh with the shipped flags (`-c -N 50 -p 0.1 -x splice -uf` / `-x asm20`).

**External check on the new panel:** genome-wide V0 on S-IDEAL reproduces the published FAMILY R of the chr16-18 panel exactly (0.8519, 23 of 27 NPIP copies in the best family), so the wider panel is not a different instrument on the NPIP axis.

## 2. Parity (asserted before any variant metric)

V0's rep-chain key set equals the Rust capture **exactly on all four arms** (0 only-Rust, 0 only-mirror): 5,428 / 7,207 / 8,292 / 8,565 keys, matching the frozen node counts. Node tables are the previously verified mirror pickles (`strand_fix/verify/vnodes_*`), not re-derived from the BAMs here.

## 3. Per-copy strand correctness (the quantity the old scorer could not express)

The strand fix is real, replicated and located at different copies on the two substrates: on S-IDEAL it rescues **NPIPB13** (same-strand frac_copy 0.000 -> 0.944) and **NPIPB4** (0.000 -> 1.000); on the realistic arms it rescues **NPIPB12**, **NPIPB1P** and **LOC124907834** (all 0.000 -> 1.000). N_opposite falls on every substrate and for every variant (IDEAL 5 -> 3/2/3 for V1/V3/V4; U20 5 -> 2/1/2; U00 3 -> 2/1/2; U40 6 -> 2/1/2), and NPIP wrong-strand best nodes go to 0 on all three realistic arms.

One finding the old scorer hid: **V3 loses the same-strand node of 4 (IDEAL) / 5 (U20) TBC1D3 copies, and in 4/4 and 5/5 of those cases the node is still there at frac_copy 1.000 but carries the 'U' (unknown-strand) label** from its consolidate. Under the declared reading ('U' is not the copy's strand) that is a cost; if 'U' were read as compatible the cost vanishes. V4 has no such cases, because it changes no strand label.

## 4. Junk, counted separately

Genome-wide junk (nodes overlapping no expressed record of the 5,542-locus set) is tiny and essentially unmoved: 5 -> 5/5/5 on S-IDEAL and 26 -> 26/25/26 on U20; on the chr16-18 panel 31 -> 31 (IDEAL) and 45 -> 46/36/46 (U20). Clause (iii) passes everywhere for everything. This clause is weak by construction here (see concerns) — V3's -20% panel junk is the U-merge deleting placeholder nodes, and is reported, not scored as progress.

## 5. Genome-wide family panel — the decisive clause

No variant on either substrate loses a family's every copy span (0 everywhere), and the largest family's best-match exon-bp Jaccard never drops below 0.9605. But families "lost" at Jaccard < 0.5 are: IDEAL V1 2, V3 1, V4 1; U20 V1 1, V3 4, V4 1. Clause (ii) therefore FAILS for all three variants on both substrates.

Opening the two V4 losses:

- **S-IDEAL family #157** (2 nodes) — under V4 both nodes are still in one family, which now also holds two newly admitted same-strand copies of the same paralog (`chr16:29746756-29749368 '-'`, `chr16:30590837-30593449 '-'`). Exon-bp Jaccard 0.4933 is produced entirely by the family growing. Calling this a "lost family" is the instrument, not the fix.
- **U20 family #114** (3 nodes) — two nodes stay together (and gain a 27-exon chr9 relative); the third, `chr15:30285739-30291421 '-'`, ends up in no family. I re-derived its edges in both node sets: **degree 1 with the same single neighbour in V0 and in V4**. A degree-1 node can only be placed with that neighbour, so it is dropped by the greedy triangle-leader assignment after the neighbour is claimed by an earlier leader. No homology evidence was lost. This is the known leader-rule instability, surfacing as a family loss.

The membership-preserving diagnostic (post-hoc, reported as a diagnostic, not a clause): V0 families whose every node survives and stays in ONE variant family — **V4: 251/252 (IDEAL) and 262/264 (U20)**; V1: 250/252 and 262/264; **V3: 250/252 and 249/264** (V3 is much worse on U20 because its U-merge deletes 21 V0 nodes outright). V4 is the only variant that removes no V0 node at all on any arm (containment: 0 V0 nodes absent, +3 / +56 / +116 / +246 new).

## 6. FAMILY R on the frozen NPIP truth

S-IDEAL: 0.8519 -> **0.9259** for V1, V3 and V4 (23 -> 25 of 27 copies in one family). U20: 1.0000 in every arm, so clause (v) passes by tie and carries no information there. FAMILY P strict is reported for continuity only (IDEAL 0.8214 -> 0.7812 / 0.8065 / 0.8065); FAMILY P NPIP-only is not reported at all.

## 7. Verdicts, as registered

V1: fails (ii) on both and (iv) on S-IDEAL, U00 and U40. V3: fails (ii) on both, and fails (i) on U00 (only 1 gaining copy). **V4: passes (i), (iii), (iv), (v) on both substrates; fails (ii) on both.** No variant is adoptable. Every clause is reported for V1, V3 and V4 on every substrate in `panel_table`; U00 and U40 have no family panel, so (ii) and (v) are NOT EVALUATED there and are never counted as passes.

## 8. What I would pre-register next (not done here, not scored here)

Clause (ii) as written cannot tell "the family was destroyed" from "the family grew" or "the greedy leader rule reshuffled one degree-1 node". A replacement that survives both failure modes: (a) every surviving node of each V0 family must stay in ONE variant family; (b) no V0 family may have its nodes scattered across >= 2 variant families; (c) the largest family's best-match Jaccard >= 0.80. V4 satisfies (a) at 251/252 and 262/264 and (c) at 1.0000 / 0.9752 — but those numbers have now been seen, so per the hold-a-substrate-back rule the replacement clause must be registered first and then tested on a substrate held back (U00/U40 or a fresh gorilla arm), not re-scored on S-IDEAL and U20.

## 9. Scope limits of this run

Mirror numbers throughout; the Rust was used only for the node/rep capture that parity was asserted against. Merged Rust nodes (n_pieces > 1, < 0.5%) are approximated by their representative chain, as inherited. V2 is not a candidate and was not given a family panel. Both substrates share their annotation with their truth, so nothing here tests node-vs-annotation fidelity. Nodes labelled 'U' were read as '+' for the panel (matching what the builder emits) and as not-same-strand for the copy clause (as declared).

## Verification (independent recompute)

Declarations first: `/mnt/linuxdisk/home/juanfraitu/strand_v4/verify/VDECL.txt` (md5 `49ad0071364889cb1c3dd937a519ad8b`, written 2026-09-18T16:13:15Z, before any number below). Full result file: `/mnt/linuxdisk/home/juanfraitu/strand_v4/verify/VERIFICATION.md`. No builder script under `strand_v4/` was read or imported; no `src/` change, no commit, no MAPQ gate change, no index rebuild, no subagent. My mirror was grounded line-by-line in `/mnt/c/Users/jfris/Desktop/Rustle/src/rustle/vg_family/shared_definition.rs` (consolidate, read_blocks, overlap_groups, depth2_exons, ExonIndex, with_read_locus_nodes, widen_with_read_isoforms, tx_key/body_key, transcript_hits, gene_body_chains, edges, triangle_leaders, build).

### 2. Node rebuild and diff

Rebuilt from the BAMs myself: `npip_ideal/bam/ideal.bam` → 1,151,970 primary MAPQ≥1 reads / 37,959 distinct; `strand_fix/aln/U20.bam` → 107,745 / 84,802.

n_nodes: S-IDEAL V0 5428, V1 5602, V3 5430, V4 5431; U20 V0 8292, V1 8584, V3 8392, V4 8408 — all equal to the builder's panel table. Node **key sets are identical** to `strand_fix/verify/vnodes_{IDEAL,U20}_{V0,V1,V3,V4}.pkl` in all 8 cells (only-mine 0, only-stored 0).

Worth stating: those stored pickles are simultaneously the previous run's output *and* the strand_v4 builder's input — the builder did not rebuild nodes. This is the first independent node rebuild in the chain.

### 3. Scorer-v3 recompute

Every copy-level, junk and locus-level cell reproduces exactly (S-IDEAL V0/V1/V3/V4 N_opposite 5/3/2/3, ss>0.5 40/43/38/42, mean ss_frac 0.8628/0.9268/0.8181/0.9051, JUNK 5/5/5/5, P_cand 0.9897/0.9615/0.9883/0.9883, R 0.9586/0.9615/0.9615/0.9615; U20 N_opposite 5/2/1/2, ss>0.5 34/38/32/37, JUNK 26/26/25/26, P_cand 0.6433/0.6263/0.6381/0.6385, R 0.9800/0.9800/0.9786/0.9800). GAIN sets and containment (0/174, 2/4, 0/3; 0/292, 21/121, 0/116) reproduce exactly.

Junk exclusion: **confirmed** — `C_panel == N_panel − junk_panel` in all 8 cells (asserted in code). Conditioning: **not fixed, slightly worse** — see Correction 4.

### 4. Genome-wide family panel

Mapped independently: 54,017 distinct tx md5 + 10,771 body md5 from my own node tables; 2,891 tx + 888 body reused from captured PAFs under `strand_fix/map`, `npip_ladder/union`, `npip_ideal`; the remaining 51,126 tx (220 Mbp) and 9,883 body (483 Mbp) mapped here against the prebuilt `target.{splice,asm20}.mmi` with `-c -N 50 -p 0.1` (`-x splice -uf` / `-x asm20`). All 19 chunks exit 0, 0 unmapped queries. The builder's `strand_v4/map/` was **not** reused.

Every panel cell reproduces: pairs, exon/body split, family count, largest, loci-in-families, families lost, losing-every-copy-span, largest-family J, families changed, mean J, FAMILY R and FAMILY P strict. The declared max-over-families matching and a one-to-one maximum-weight bipartite matching give the **same** families_lost in all 8 cells.

It really is genome-wide: families touch all 24 contigs; chr16/17/18 hold only 61 of 252 (S-IDEAL) and 68 of 264 (U20); 63–69 families are multi-chromosome; node coverage 972/5428 = 17.9% and 1351/8292 = 16.3%. That is >2× the 121-family panel the `RUSTLE_COLLAPSE_UNSTRANDED` caution was measured on, so the blast-radius complaint is genuinely answered.

### 5. Clause-by-clause, from my numbers

| substrate | variant | (i) | (ii) | (iii) | (iv) | (v) |
|---|---|---|---|---|---|---|
| S-IDEAL | V1 | PASS (5→3, 3 gains) | **FAIL** (lost 2) | PASS (5→5, 0%) | **FAIL** (−0.0282) | PASS (0.852→0.926) |
| S-IDEAL | V3 | PASS (5→2, 2 gains) | **FAIL** (lost 1) | PASS | PASS (−0.0014) | PASS |
| S-IDEAL | V4 | PASS (5→3, 2 gains) | **FAIL** (lost 1) | PASS | PASS (−0.0014) | PASS |
| U20 | V1 | PASS (5→2, 4 gains) | **FAIL** (lost 1) | PASS (26→26) | PASS (−0.0170) | PASS (vacuous, 1.000) |
| U20 | V3 | PASS (5→1, 3 gains) | **FAIL** (lost 4) | PASS (26→25) | PASS (−0.0051) | PASS (vacuous) |
| U20 | V4 | PASS (5→2, 3 gains) | **FAIL** (lost 1) | PASS (26→26) | PASS (−0.0048) | PASS (vacuous) |

**Does the recommendation follow? Yes.** V4 fails clause (ii) as pre-registered on both substrates, by exactly one family of 252 and one of 264, so "NOT ADOPTABLE under the rule as pre-registered, do not adopt now" is the correct verdict, and re-registering a clause after seeing the number would be the wrong move. The builder's stated plan — pre-register a membership clause plus a no-scatter test and evaluate on a held-back substrate (U00/U40 or a fresh gorilla arm) — is the right next step.

**But the blocking evidence is much weaker than the report implies.** On U20 it is inside the leader rule's own tie-break noise (Correction 3). On S-IDEAL it is an edge-free re-partition of four pre-existing nodes that V4 arguably groups better (Corrections 1–2). And the proposed replacement clause fails V4 too, on different families (Correction 6). The honest summary is not "V4 breaks a family" but "clause (ii) cannot currently tell a strand fix from the greedy leader rule's own instability", and that should be fixed before the held-out run, or the held-out run will inherit the same ambiguity.

### 6. Flags — things that would only hold on the idealized substrate

1. **The leader rule is stable on S-IDEAL and unstable on U20.** Noise floor 0 lost families (5/5 seeds) vs 1 lost family (2/5 seeds). So "V4 loses a family" is only *measurable* on the idealized arm; on the realistic arm it is noise.
2. **Clause (iii) is a zero-tolerance test only because S-IDEAL has 5 junk nodes.** At any realistic junk level the 10% rule is far looser than what was actually tested.
3. **Clause (v) is informative only on S-IDEAL.** On U20 FAMILY R is pinned at 1.000 by an 80–87-node family with FAMILY P strict ≈ 0.33.
4. **The 0.02 absolute tolerance in (iv) means different things at P=0.99 and P=0.64.** On U20, 381 of 1068 candidate nodes are already unmatched at V0 — fragmentation dominates and the strand signal is a ~1–2% perturbation on top of it.
5. **V4's reach differs by an order of magnitude between arms.** It recovers 3 of V1's 174 new nodes on S-IDEAL (1.7%, all on chr16, 2 '+' / 1 '−') but 116 of V1's 292 on U20 (40%, 23 chromosomes, 115 of 116 on '−'). The "blockers whose strand was measured" restriction is nearly inert on the error-free substrate, because almost every base node there has a junction read.
6. **The whole mechanism is defined against one catalog convention** — every unspliced-only locus is given strand '+' (S-IDEAL 451/451 single-exon rep base nodes, U20 5031/5031). V4 is, in effect, "do not let an unspliced-only '+' node block a '−' read locus". Nothing here tests what V4 does if a real catalog assigns '−' to some unspliced loci.
7. **Both substrates share their annotation with their truth**, so junk means "off the frozen 5,542-locus set", R ≈ 0.96–0.98 by construction, and none of this tests node-vs-annotation fidelity. (The builder declares this; I confirm it bounds every number in the table.)
8. **One cross-chromosome family edge to watch**: U20 V4's new `chr9:140462942-140483121 '−'` (27 exons) joins family #114 with degree 1, through a 265 bp single-exon `chr15:30285739-30286004 '+'` node. Single-exon nodes acting as cross-chromosome bridges is exactly the failure mode the strand work is trying to remove, and V4 creates one here.