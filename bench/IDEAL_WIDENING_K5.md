# Read-isoform widening at k = 5 on the idealized simulated substrate (agent 1 of 2)

Declarations written BEFORE any number: `/mnt/linuxdisk/home/juanfraitu/ideal_widened/DECLARATIONS.txt`, md5 of the frozen text `826914f0bd3ca66ff4f97244ce37efa0`, written 2026-09-17T21:51:19-07:00; ADDENDUM A5 appended after the arms ran and before this report (file md5 now `42c9b07d204508be784e42ec2a65b893`). Outputs: `out/GAP_CLOSED.tsv`, `out/{locus,copies,family}.json`, `out/nodes_*.pkl`, `logs/{OFF,K5,K3,K8}.err`, `logs/parity_nodes.txt`, `cap/<arm>/{tx.fa,body.fa}`, `map/*.paf`. Binary `/mnt/linuxdisk/home/juanfraitu/rustle_target/release/gw_family_catalog`, sha256 `fde7cf2b6eaa5ec371eb6fb164901f03f91793d28c7e2220c83c367e8a69f59d`. Nothing in `src/` was modified, nothing committed, no subagents.

> **Verifier corrections, applied by the orchestrator** (verifier ok = true; the widening implementation was re-derived genome-wide over all 5,428 nodes with diff 0 in both directions):
> 1. **k was not measured and cannot be on this substrate.** The minimum junction support over all 37,260 (node, chain) instances is exactly 30, because the simulation emits 30 identical reads per transcript — so every k ≤ 30 admits the same chains. k = 3/5/8 outputs are byte-identical. The k = 5 default rests on the real-data measurement, not on this run.
> 2. **The circularity is tighter than the generic disclaimer**: 4,062 of the 4,241 distinct spliced node-query chains (95.8%) are byte-exact annotated transcript intron chains, because the reads were simulated from those transcripts. The containment gains are partly "recovering what was simulated".
> 3. **The account of the 5 short NPIP copies is corrected.** Per copy: NPIPB13 0.028 → 0.028, NPIPB4 0.033 → 0.033, NPIPB12 0.042 → 0.042, NPIPB14P 0.122 → 0.122, NPIPB15 0.455 → 0.704. NPIPB12 genuinely has 0 same-strand primary MAPQ ≥ 1 reads; the others are blocked by an opposite-strand node at the locus, which strand-strict widening cannot touch.
> 4. **Anti-trap columns were scored with the wrong sign**: moving toward the annotated arm's count of nodes overlapping ≥ 2 records is not progress — for a trap column the annotation's value is not a target.
> 5. **The summary omitted the regressed rows**: component-level FAMILY P strict 0.3418 → 0.3375 and F strict 0.5094 → 0.5047.
> 6. **The one disagreeing copy-table cell** (NPIPB8) is attributed to the wrong cause; the verifier's recompute agrees with the published table, not with the builder.
> 7. Unresolved and bounded: the frozen substrate node set and the published dump disagree on one node's interior exon union by 16 bp.

## 0. CIRCULARITY (verbatim, as required)

"The reads of this substrate were simulated FROM the same annotation that defines the truth, so this substrate cannot measure node-vs-annotation fidelity in general — it measures how much of the annotated structure the construction can recover when nothing is missing for lack of expression, i.e. a CEILING."

Consequences honoured below: a metric that reaches the annotated arm here is an upper bound, not validation; a metric that FAILS to reach it here is the informative direction, because the ceiling already fails. This run does not validate k = 5; the default was set in `src/` before the run, and this run only measures it.

## 1. Reading rule (verbatim, declared before any number)

"PERFECT means the arm reaches the annotated arm's own value on that metric (locus bipartite F 1.000, isoform containment 1.000, full-length copies 27/27, FAMILY R / P strict / F strict equal to the annotated arm's). Report the gap-closed fraction (variant - baseline) / (annotated - baseline) for every metric, and state for each whether it is PERFECT, CLOSED (>= 80% of the gap), PARTIAL, or FLAT."

Labels as declared: PERFECT = equal to the annotated arm; CLOSED = gc >= 0.80; PARTIAL = 0.05 <= gc < 0.80; FLAT = |gc| < 0.05; REGRESSED = gc < -0.05.

⚠ Two of the rule's nominal targets are NOT attainable as stated, and this was found, not assumed: the annotated arm's own full-length count on this substrate is **26/27, not 27/27** (NPIPB12 has **0** same-strand primary MAPQ>=1 reads on the ideal BAM — 50 MAPQ-0 reads — so it is not an expressed record and the ANN arm has no node for it), and the annotated arm's FAMILY P strict / F strict (0.619 / 0.754) are **below** both de novo arms (0.821 / 0.836), so "reaching the annotated arm" there would be a regression. Both are reported as they are.

## 2. Arms and how they were run

| arm | setting | Rust `[shared-definition]` stderr | wall | peak RSS |
|---|---|---|---|---|
| OFF | `RUSTLE_SD_READ_ISOFORM=0` | `5375 reps -> 5360 gene-level loci + 68 read-locus nodes = 5428 nodes` · `read-isoform widening OFF` | 6:24.69 | 21.70 GiB |
| K5 | shipped default (k = 5) | same 5428 nodes · `read-isoform widening k=5: 3896 of 5428 nodes widened, 37796 spliced queries` | 7:08.98 | 21.70 GiB |
| K3 | `RUSTLE_SD_ISOFORM_K=3` | identical line with `k=3`, same 3896 / 37796 | 7:05.98 | 21.70 GiB |
| K8 | `RUSTLE_SD_ISOFORM_K=8` | identical line with `k=8`, same 3896 / 37796 | 6:28.32 | 21.70 GiB |
| ANN | expressed annotated records as nodes (ceiling) | n/a | n/a | n/a |

All four Rust arms ran `RUSTLE_SHARED_DEFINITION=1` on `npip_ideal/bam/ideal.bam` against `npip_ladder/idx/target.fa`, `TMPDIR` under `ideal_widened/tmp` (so the 3.1 GB `target.fa` the build writes never touched the C: vhdx; it was deleted after each arm).

**DISCLOSED DEVIATION (declared in advance, A1.4).** `build()` computes and widens the node set and writes `target.fa`/`tx.fa`/`body.fa` BEFORE either minimap2 call; only `edges()` and `triangle_leaders()` consume the mappings. The whole-genome edge stage does not fit the 10-minute call cap (it is what crashed the original run), so the Rust arms ran with a stub minimap2 that captures the query files and exits non-zero. Edges were then computed by the frozen Python mirror of `shared_definition.rs` (`npip_ladder/scripts/ladder.py`, the mirror that produced the committed A0–A4 ladder rows) on the SAME node sets, with minimap2 2.30 run separately — `-c -N 50 -p 0.1 -x splice -uf -t 4` and `-c -N 50 -p 0.1 -x asm20 -t 4` — against the prebuilt `target.{splice,asm20}.mmi`. No index was rebuilt. minimap2 cost: tx 4 batches (17.4 Mbp new) 58.2 + 66.0 + 45.4 + 32.6 s at 14.6–16.6 GB RSS; body 3 batches (42.1 Mbp new) 47.4 + 40.7 + 35.5 s at 12.7–13.1 GB RSS.

**PARITY (the check that the new code path is inert when disabled).**
- P1: OFF `tx.fa` is **byte-identical** to the published capture — md5 `7d8de17f6a4b2d2f7e84b1266a672ea6` both sides; the node line reproduces `5428 nodes` exactly. PASS.
- P2: the mirror node set equals the Rust node set **exactly** in both arms: chr16/17/18 (710 nodes) tx-key-set diff 0 and body-key-set diff 0 (OFF 710 tx keys; K5 4296 tx keys); chr1+chr4 (747 nodes) diff 0 in both arms (K5 5216 tx keys). One disclosed correction was needed (ADDENDUM A5.2): 5 of the 710 frozen substrate-(a) nodes are `consolidate`-merged nodes whose `rep_exons` the earlier substrate build had replaced with the reconstructed exon union; the true rep chain was recovered from this run's own OFF `tx.fa` by unique containment. Node exon unions unchanged. PASS.
- P3: the OFF copy table reproduces the published `nodes/copy_to_node.IDEAL.tsv`: 46/46 identical `full_length` labels, 45/46 identical `frac_copy`; the single difference (NPIPB8, 1.000 vs 0.9876) comes from the merged-node exon reconstruction.
- P4: **K3 == K5 == K8 byte-for-byte** (tx.fa md5 `f9a6027…`, body.fa md5 `ddc10c7…`).

## 3. ADDENDUM A5.1 — the k sensitivity is VACUOUS on this substrate

Arm IDEAL emits 30 identical error-free reads per transcript, so every junction observed at a node carries >= 30 reads and every k <= 30 admits exactly the same chains. K3, K5 and K8 are therefore the same run. **This substrate cannot choose k, and no claim that "k = 5 is robust" can be made from it.** The k question belongs on real reads (the committed read-isoform run's dev/held-out substrates), not here. This is disclosed, not designed.

## 4. Results

### 4.1 Copy recovery, 46 truth copies (27 NPIP + 19 TBC1D3), npip_ideal D6 conventions

| | OFF | K5 | ANN |
|---|---|---|---|
| NPIP copies with a node | 27/27 | 27/27 | 26/27 |
| **NPIP full-length (>= 0.8 of copy exon union)** | **11/27** | **22/27** | **26/27** |
| NPIP mean copy coverage | 0.5897 | 0.8442 | 0.9630 |
| TBC1D3 full-length | 17/19 | **19/19** | 19/19 |
| TBC1D3 mean copy coverage | 0.9558 | 0.9993 | 1.0000 |
| all 46 full-length | 28/46 | 41/46 | 45/46 |
| all 46 mean copy coverage | 0.7409 | 0.9082 | 0.9783 |
| copies sharing a node | 0 | 0 | 0 |
| nodes swallowing a neighbouring locus (>= 50 bp) | 14 | 14 | 13 |
| copies whose node is on the opposite strand | 5 | 5 | 0 |

Widening never merged two truth copies into one node and never increased swallowing. Note the anti-trap reading: the ANN ceiling itself has 13 nodes overlapping another locus-set locus by >= 50 bp, so 14 is 1 above the ceiling, not 14 above zero.

### 4.2 Isoform containment (substrate (a), chr16/17/18, 701 expressed records)

C1 = fraction of the record's transcripts with >= 0.999 of their exonic bases inside the node exon union; C2 = fraction whose every junction is in the node's query junction set (rep chain + admitted chains). Declared convention: computed over ALL expressed records (an unmatched record contributes 0); matched-only means reported alongside.

| | OFF | K5 | ANN |
|---|---|---|---|
| C1, all records | 0.5791 | 0.9338 | 1.0000 |
| C2, all records | 0.4601 | 0.9193 | 1.0000 |
| C1, matched only | 0.6032 | 0.9741 | 1.0000 |
| C2, matched only | 0.4792 | 0.9590 | 1.0000 |
| records with ALL transcripts contained (C1) | 265 | 633 | 701 |
| records with ALL chains contained (C2) | 202 | 612 | 701 |
| annotated junctions of the scope present in some node query (of 9,677) | 6,858 (0.7087) | 9,528 (0.9846) | 9,677 (1.000) |

The OFF C1-matched 0.6032 reproduces the published 0.603 for this substrate.

### 4.3 Locus level, bipartite, with both anti-trap columns

| | OFF | K5 | ANN |
|---|---|---|---|
| nodes / records / matched | 710 / 701 / 673 | 710 / 701 / 672 | 701 / 701 / 701 |
| R / P / **F** | 0.9601 / 0.9479 / **0.9539** | 0.9586 / 0.9465 / **0.9525** | 1.000 / 1.000 / **1.000** |
| matched-pair exon recall / precision | 0.8583 / 0.9907 | **0.9929** / 0.9904 | 1.000 / 1.000 |
| A1 nodes overlapping NO expressed record | 31 | 31 | 0 |
| A2 nodes overlapping >= 2 records (>= 50 bp) | 61 | 63 | 89 |
| expressed records with no node | 28 | 29 | 0 |

This is the sharpest single statement in the run: **widening moves the exon content of a matched node from 86% to 99% of its record, and moves the locus-level F by −0.0014.** The bipartite F counts matches, not content, so it is blind to exactly the thing widening fixes. Conversely the 31 spurious nodes and the 28 records with no node — the whole locus-level gap — are untouched by construction, because widening cannot create, delete or merge a node.

### 4.4 FAMILY level, frozen NPIP scorer (Dishuck-checked truth), triangle-supported leaders

| | OFF | K5 | ANN |
|---|---|---|---|
| exon / body edges, pairs | 587 / 623, 795 | 637 / 640, 829 | 675 / 803, 945 |
| families (triangle) | 47 | 48 | 46 |
| **FAMILY R** | **0.8519** | **0.8519** | **0.9630** |
| FAMILY P NPIP-only | 1.000 | 1.000 | 1.000 |
| **FAMILY P strict** | **0.8214** | **0.8214** | **0.6190** |
| **FAMILY F strict** | **0.8364** | **0.8364** | **0.7536** |
| matched family = copies + others | 23 + 5 | 23 + 5 | 26 + 16 |
| components: FAMILY R / P strict | 1.000 / 0.3418 | 1.000 / 0.3375 | 0.9630 / 0.3714 |
| full-length copies (scorer's own M1) | 12/27 | 23/27 | 26/27 |

Every triangle-level FAMILY number is **bit-identical** between OFF and K5 despite 34 extra pairs: the copies that widening lifts were already inside the family, and the copies outside it stay outside. Under components both de novo arms already put **all 27 copies in one component** (R 1.000), above the ANN arm's 0.963 (which loses NPIPB12). ⚠ FAMILY P NPIP-only is 1.000 in every row and is **vacuous** here, exactly as the metric-traps register says; only R, P strict and F strict discriminate. Component-level P strict moves 0.3418 → 0.3375 (one more node pulled in): the only cell that regresses, and it is 0.004.

### 4.5 Gap-closed table (the declared reading)

Full table in `out/GAP_CLOSED.tsv`. Summary by label:

- **PERFECT** (equals the annotated arm): TBC1D3 full-length 19/19; copies sharing a node 0.
- **CLOSED (>= 80%)**: C1 all records 0.843 · C2 all records 0.851 · C1 matched 0.935 · C2 matched 0.921 · junction presence 0.947 · records with all transcripts contained 0.844 · records with all chains contained 0.822 · matched-pair exon recall 0.950 · TBC1D3 mean coverage 0.983.
- **PARTIAL**: NPIP full-length 0.733 (11 → 22 of a possible 26) · NPIP mean coverage 0.682 · all-46 full-length 0.765 · all-46 mean coverage 0.705 · scorer full-length 0.786 · scorer mean coverage 0.747 · A2 0.071.
- **FLAT**: locus R (−0.036), P (−0.027), **F (−0.031)** · matched-pair precision (−0.027) · A1 (0.000) · records with no node (−0.036) · FAMILY R, P strict, F strict, and components R (all exactly 0.000).
- **REGRESSED**: components P strict, gc −0.144 in a metric whose absolute move is 0.3418 → 0.3375.

### 4.6 What is left — and it is not isoform structure

The 5 NPIP copies that stay short under K5, with their read evidence on `ideal.bam` (primary, same-strand, MAPQ >= 1 / MAPQ 0 over the copy's exon union) and the nodes that exist at their locus:

| copy | txs | exon bp | mq>=1 same-strand | mq0 | frac OFF → K5 | node(s) at the locus |
|---|---|---|---|---|---|---|
| NPIPB12 | 3 | 3,737 | **0** | 50 | 0.042 → 0.042 | one `+` node only (copy is `-`) |
| NPIPB13 | 3 | 6,144 | 120 | 56 | 0.028 → 0.028 | one `+` node only (copy is `-`) |
| NPIPB4 | 1 | 3,634 | 30 | 11 | 0.033 → 0.033 | one `-` node only (copy is `+`) |
| NPIPB14P | 1 (exon-less span) | 19,826 | 60 | 0 | 0.122 → 0.122 | `-` node covering 0.122 of the 19.8 kb span |
| NPIPB15 | 7 | 4,340 | 90 | 66 | 0.455 → **0.704** | `+` node = 1.000 of itself, 0.704 of the copy |

Mechanism, stated as a finding and not a fix: `with_read_locus_nodes` suppresses a new read-locus node whenever the read's blocks touch ANY existing node's exons — `ExonIndex` is **strand-blind** — while `widen_with_read_isoforms` only admits chains from **same-strand** reads. An antisense neighbour node therefore both blocks the copy's own node from being created and is immune to widening by the copy's reads. That asymmetry, plus NPIPB12's total MAPQ-0 starvation (an O2 problem appearing inside O1's evidence, on a substrate where every locus is fully expressed), is the entire residual on the NPIP side. NPIPB14P is the exon-less pseudo-transcript record and is a substrate convention, not a defect.

### 4.7 Cost of the default

Widening multiplies the whole-genome spliced query set: 5,428 → 37,796 queries (37,790 distinct after `tx_seen` dedup), and the tx query FASTA 23.8 MB → 187 MB (7.9×). The gene-body query set is unchanged: the 710 chr16-18 body keys and the 747 chr1/chr4 body keys are **identical** between arms, i.e. widening never extended a node's span on this substrate. The Rust arm's own wall time rose 6:24.7 → 7:09.0 with peak RSS unchanged at 21.7 GiB; the extra cost lands on the downstream edge stage, which is ~8× more spliced query sequence to map genome-wide.

## 5. Answer to the question asked

**With k = 5 on simulated data, which metrics are perfect?** TBC1D3 full-length copies (19/19, equal to the annotated arm) and copies-sharing-a-node (0). Nothing else is exactly perfect.

**Which are close?** Everything that measures what is inside a node: isoform containment C1 0.934 and C2 0.919 over all 701 expressed records (ANN 1.000), 0.974 / 0.959 over matched records, 98.5% of annotated junctions now present in some node query, matched-pair exon recall 0.993. All CLOSED, 84–95% of the gap.

**Which are not?** (i) Per-copy completeness on NPIP: 22/27 full-length against a 26/27 ceiling, mean coverage 0.844 against 0.963 — PARTIAL, ~70% of the gap. (ii) Everything at locus and family level is FLAT: locus F 0.9539 → 0.9525 (ANN 1.000) and FAMILY R / P strict / F strict identical to the un-widened arm (0.852 / 0.821 / 0.836; ANN 0.963 / 0.619 / 0.754).

**What is left?** Three things, none of them isoform structure. (1) Node identity, not content: widening cannot create, delete or merge a node, so the 31 nodes matching no record and the 28–29 records with no node — the whole locus-level F gap — are untouched by design. (2) The strand-blind suppression / strand-strict widening asymmetry described in §4.6, which costs NPIPB13 and NPIPB4 their own node. (3) MAPQ-0 read loss: NPIPB12 has zero usable reads even though every locus is expressed at depth 30 — the ceiling substrate already fails there, which is the informative direction. And (4) as a measurement fact rather than a defect: this substrate cannot choose k at all (§3), so the k = 5 default has to keep resting on the real-read runs.

## Verification (independent recompute) — agent 2 of 2

Code and outputs: `/mnt/linuxdisk/home/juanfraitu/ideal_widened/verify/` (`w1_prenodes.py` … `w11_leak.py`, `vlocus.json`, `vfamily.json`, `vwide_K5.pkl`, `vper_node.pkl`). I did not read or import anything under `ideal_widened/scripts/`; I reused my own prior verifier code (`npip_ideal/verify/v3_reads.npz` read store, `npip_ladder/verify/v_edges.py` mirror). Nothing in `src/` was touched, nothing committed, no index rebuilt, no subagent spawned.

### 1. Declarations precede results — PASS
`DECLARATIONS.time` = 2026-09-17T21:51:19-07:00; earliest result file `cap/OFF/tx.fa` 21:58. The text the addendum freezes is intact byte-for-byte: md5 of the first 11,948 bytes (everything through A4.4) is `826914f0bd3ca66ff4f97244ce37efa0`, exactly the value recorded at 21:51 in `DECLARATIONS.md5` and quoted at the A5 line; the current whole-file md5 is `42c9b07d…`, i.e. only appended. Binary sha256 `fde7cf2b…f59d` as stated (built 21:45); `src/rustle/vg_family/shared_definition.rs` mtime 21:40:35, before the declarations; last commit 20:51, none after; `npip_ladder/idx/*.mmi` untouched (Sep 17 00:29–00:39).

### 2. OFF-arm parity — PASS, and verified more strongly than claimed
- `cap/OFF/tx.fa` md5 `7d8de17f6a4b2d2f7e84b1266a672ea6` = `npip_ideal/rust/capture/rust_tx.fa`. The `[shared-definition]` lines reproduce `5375 reps -> 5360 gene-level loci + 68 read-locus nodes = 5428 nodes` plus `read-isoform widening OFF`.
- Independently of the builder's scope: I rebuilt the pre-widening node set from the published `npip_ideal/nodes/_nodes.json` exon unions, recovering each node's rep chain from the Rust OFF `tx.fa` itself (exact exon-union match for 5,411 nodes; unique containment for the 17 `reconstructed_merge` nodes, 0 ambiguous). The result reproduces the Rust OFF capture **genome-wide**: tx keys 5,428/5,428 diff 0, body keys 5,427/5,427 diff 0 (the builder asserted parity only on 710 + 747 nodes).
- Copy table: my OFF arm reproduces `copy_to_node.IDEAL.tsv` on **46/46** `full_length` labels and **46/46** `frac_copy` at published precision (see correction 6 about NPIPB8).

### 3. Widening re-derived from the BAM with my own implementation — EXACT MATCH, invariants hold
From `v3_reads.npz` (1,151,970 primary / non-supplementary / MAPQ>=1 reads, `ts:A:-` flipped; 37,959 distinct read block-tuples) I implemented the declared rule from scratch (per-node same-strand read assignment by >=1 bp exon overlap, chain support and min-start/max-end, junction support summed over the node's chains, `exon_blocks` reconstruction, merge into the exon union) and ran it at k = 3, 5, 8, genome-wide:

| | mine | Rust K5 capture |
|---|---|---|
| nodes widened | 3,896 / 5,428 | 3,896 / 5,428 |
| spliced queries | 37,796 | 37,796 |
| unique tx keys | 37,790 | 37,790 (diff 0) |
| body keys | 5,425 | 5,425 (diff 0) |

k = 3, 5 and 8 give byte-identical output in my implementation too, independently confirming A5.1. Invariants checked on all 5,428 nodes: node count unchanged (5,428 → 5,428), **0** nodes whose exon union shrank or lost a pre-existing exon, **0** nodes whose rep chain was dropped, no node created/removed/merged. Every admitted chain satisfies the >= k conjunct by construction — and, as the junction audit shows, trivially so (min junction support = 30 everywhere). Scope cross-check of the builder's parity line: on chr16/17/18, 472 nodes gain a chain of which 409 also gain exon bases, total exon bp added 798,259 — the builder's `exon bp added: total 798259` exactly.

### 4. Metrics recomputed from my own node sets

M-A (46 copies, D6 conventions), M-B/M-C (701 expressed records on chr16/17/18 — I re-derived expressedness from the BAM and got 701 independently), M-D (edges → triangle leaders / components → FAMILY):

| metric | OFF | K5 | ANN | gc | label | builder |
|---|---|---|---|---|---|---|
| NPIP full-length /27 | 11 | 22 | 26 | 0.733 | PARTIAL | same |
| NPIP mean copy coverage | 0.5893 | 0.8437 | 0.9630 | 0.681 | PARTIAL | 0.5897 / 0.8442 |
| TBC1D3 full-length /19 | 17 | 19 | 19 | 1.000 | PERFECT | same |
| all 46 full-length | 28 | 41 | 45 | 0.765 | PARTIAL | same |
| C1 all / C2 all records | 0.5777 / 0.4615 | 0.9324 / 0.9200 | 1.0 | 0.840 / 0.851 | CLOSED | 0.5791 / 0.4601 → 0.9338 / 0.9193 |
| C1 / C2 matched only | 0.6017 / 0.4807 | 0.9727 / 0.9597 | 1.0 | 0.931 / 0.922 | CLOSED | ±0.002 |
| annotated junctions in a node query | 0.7087 | 0.9846 | 1.0 | 0.947 | CLOSED | identical |
| locus bipartite R / P / F | 0.9601 / 0.9479 / 0.9539 | 0.9586 / 0.9465 / 0.9525 | 1.0 | −0.03 | FLAT | identical to 7 dp |
| matched-pair exon recall / precision | 0.8583 / 0.9908 | 0.9929 / 0.9905 | 1.0 | 0.950 / −0.027 | CLOSED / FLAT | identical |
| A1 / A2 / records with no node | 31 / 61 / 28 | 31 / 63 / 29 | 0 / 89 / 0 | — | FLAT / (see correction 4) / FLAT | identical |
| FAMILY R / Ps / Fs (triangle) | 0.8519 / 0.8214 / 0.8364 | identical | 0.9630 / 0.6190 / 0.7536 | 0.000 | FLAT | identical |
| FAMILY R / Ps / Fs (components) | 1.0 / 0.3418 / 0.5094 | 1.0 / 0.3375 / 0.5047 | 0.9630 / 0.3714 / 0.5361 | −0.144 / −0.179 | REGRESSED | Ps identical |

Every PERFECT / CLOSED / PARTIAL / FLAT / REGRESSED label in `out/GAP_CLOSED.tsv` is reproduced. Containment cells differ by ≤ 0.002 and two counts by 1 (`C1_records_all_tx` 264 vs 265, `C2` 203 vs 202), attributable to the reconstructed exon unions of the 5 merged nodes on the scope. The junction fraction matches exactly (0.7087 / 0.9846 / 1.0000) under the declared "present in SOME node query" convention; under the stricter matched-node convention it is 0.9756 for K5 — still CLOSED.

M-D caveats: my edge stage reuses the existing md5-keyed PAF stores rather than re-running minimap2 (I verified every one of my 4,542 tx and 799 body queries had been mapped with identical flags against the same prebuilt index, and that queries present in two stores carry identical records — 0 disagreements in 200 sampled). My exon-edge entry counts are 569 (OFF) / 619 (K5) vs the builder's 587 / 637, and distinct pairs 790 / 824 vs 795 / 829 (0.6%); body edges (623 / 640 / 803) and **all** family-level outputs agree exactly once the leader ordering uses the node's MAPQ>=1 read count (`n_reads_mq1`) — with that convention I get OFF 47 families / 214 loci / 23 matched / strict 28 and K5 48 / 219 / 23 / 28 and ANN 46 / 205 / 26 / 42, i.e. the builder's numbers to the digit. Note that the triangle FAMILY numbers are sensitive to that ordering (with a different read-count proxy ANN's strict size moves 42 → 47, Ps 0.619 → 0.553), while the OFF-vs-K5 comparison is not: they stay identical under both conventions.

### 5. New anti-trap I added (cross-locus leakage)
Widening assigns a read to every same-strand node any of its blocks touches, so a chain can leak into a neighbour. Measured on the scope: of the 798,259 exon bases widening adds, 7,986 (1.0%) fall inside a *different* expressed record's exon union, and 17 widened nodes add >= 50 bp of another record. That is the mechanism behind A2 61 → 63 and it is small — but it is a real, one-directional cost that grows with read noise.

### 6. Does the answer follow from the numbers?
Yes, with the qualifications above. Widening buys node CONTENT and nothing else: copy coverage, containment, junction recovery and matched-pair exon recall close 68–95% of the gap to the annotated ceiling; node COUNT/IDENTITY is untouched (locus F 0.9539 → 0.9525, A1 31 → 31, records with no node 28 → 29) and every FAMILY-level number is bit-identical between OFF and K5 apart from two component-level precision/F cells that fall slightly. "Perfect" is reached on exactly two cells (TBC1D3 19/19 full-length, 0 copies sharing a node), and two of the reading rule's nominal targets are unreachable as written: the annotated arm's own full-length count is 26/27 (NPIPB12 is not expressed), and the annotated arm's FAMILY P strict / F strict sit *below* both de novo arms. Simulation artefacts to keep front and centre: (i) 30 identical error-free reads per annotated transcript make the k gate inert and make the chain set the annotation's own chain set (95.8% byte-exact), so the CLOSED cells are closer to a self-consistency check than to a fidelity measurement; (ii) the residual failures that survive are the ones the simulation cannot flatter — missing expression (NPIPB12), the strand-blind node-suppression rule (NPIPB13, NPIPB4), and whatever keeps NPIPB14P at 12% — and those are the informative direction, because they fail even at the ceiling.