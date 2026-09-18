# Read-isoform locus and node-admission floor, pre-registered — agent 1 of 2

Declarations written before any number of this run: `/mnt/linuxdisk/home/juanfraitu/readisoform/DECLARATIONS.txt`, mtime 2026-09-17 19:56:33 local (`DECLARATIONS.time` = 2026-09-18T02:56:33Z), md5 `b3564799977a57bd68c135d35ac01c7b`. Earliest result file 19:57:04; `out/R3_CHOICE.txt` 20:07:10 precedes every R3 number (locus 20:07:14 / 20:07:33, family 20:08:24). Nothing in `src/` was modified, nothing was committed, no subagents were spawned.

Outputs: `/mnt/linuxdisk/home/juanfraitu/readisoform/out/SUMMARY_readisoform.tsv` (every substrate × variant row with the anti-trap columns), `out/locus_{dev,held,burn}.json`, `out/family_dev.json`, `out/R3_CHOICE.txt`, scripts in `scripts/` (`ri_build_h.py`, `ri_lib.py`, `ri_locus.py`, `ri_family.py`), new alignments in `map/` (`tx_r0{0,1,2}.paf`, `body_r00.paf`).

> **Verifier corrections, applied by the orchestrator** (verifier ok = true; all 570 locus cells and 9 family rows reproduced):
> 1. **The pre-registered decision rule is INFEASIBLE — the annotation itself fails it.** Running the expressed annotated records as the node set (the exact target the goal names) through the five clauses on the held-out substrate fails clauses (ii) and (v). Every arm that moves locus F fails "no record loses its node"; every arm that respects it moves F by ≤ 0.005. The next rule must be stated RELATIVE TO THE ANNOTATED ARM, not in absolutes.
> 2. **R1 widens but does not preserve footprints**: R1 exons are a superset of R0's for every node on all three substrates (0 shrink, node count unchanged), but at k = 3 on the held-out substrate 85 R1 nodes fully contain another shipped node's exon set, against 39 for R0.
> 3. **"The W9 gain came from anchoring on the rep's splice structure" is an untested post-hoc explanation** — no arm here isolates that factor, and it is not what the W9 code did.
> 4. **Clause (v) fails more widely than the verdict says**: all four admission arms fail on FAMILY R (0.889-0.963), R3 at 0.963, and R1k2/k3 fail on P strict (0.290 → 0.153 / 0.175).
> 5. **R2's locus-F gain is 87% denominator shrinkage** (restricting R0's own matching to the surviving node set already gives F 0.795 of the 0.823), and the deletions are not all noise: 276 of 848 deleted nodes overlap an expressed record.

## 0. What was pre-registered, and parity

The arms below were POST HOC in the locus-width run (verifier correction 8) and therefore unquotable. This run fixes them in advance: the exact R1/R2/R3 definitions, the assignment rule, the chain-block convention, two written isoform-containment conventions (C1 exon-bp, C2 intron chain), the two anti-trap columns, and the decision rule verbatim.

- **dev parity asserted in code**: the R0 node set is byte-for-byte equal, in order, to the frozen A4 node set in `npip_ladder/verify/nodes.pkl` (560 nodes; `assert`, not prose). Its family row reproduces the frozen figures exactly — 1349 exon edges / 1885 body edges / 1919 pairs / 68 families / FAMILY R 1.000, Pm 1.000, P strict 0.290, F strict 0.450, 5/27 full-length, mean copy coverage 0.517.
- **held-out = chimpanzee NC_072416.2**, never used for any read-isoform or admission arm: 1,589 shipped de novo nodes (`ggo_npip/nodes/PTR.dn.nodes.tsv`), 64,255 primary MAPQ ≥ 1 reads from `ggo_npip/dn/PTR.1ctg.bam` (the NC_072416.2 slice of `winloci_data/PTR_mm.bam`), 1,571 annotated records with exons of which **841 expressed** (≥ 3 same-strand primary MAPQ ≥ 1 reads), 5.54 transcripts per expressed record.
- **gorilla is reported and marked BURNED** (post-hoc W7/W9 were already run there); it carries no weight in the decision rule.
- minimap2 was re-run only for query md5 keys absent from the captured PAFs (`npip_ladder/union/`, `locus_width/map/`): 5,418 new tx md5 (14.1 Mbp, 3 batches, 106+113+86 s, `-c -N 50 -p 0.1 -x splice -uf -t 4`) and 266 new body md5 (18.7 Mbp, 64 s, `-c -N 50 -p 0.1 -x asm20 -t 4`) against the prebuilt `npip_ladder/idx/target.{splice,asm20}.mmi`.

## 1. HELD-OUT substrate (chimp), locus level — the table the rule is read on

Anti-trap columns: **A1** = nodes overlapping NO expressed annotated record; **A2** = nodes overlapping ≥ 2 expressed records by ≥ 50 exonic bp each.

| variant | nodes | matched | R | P | **F** | ΔF | pair R | pair P (median) | **A1** | **A2** (rel) | C1 | C2 | recs w/o node | lost vs R0 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| R0 | 1589 | 743 | 0.883 | 0.468 | **0.612** | — | 0.649 | 0.841 (0.998) | 604 | 150 | 0.079 | 0.229 | 98 | 0 |
| R1k2 | 1589 | 749 | 0.891 | 0.471 | 0.616 | +0.005 | 0.744 | 0.622 (0.615) | 592 | 234 (+56.0%) | 0.161 | 0.411 | 92 | 8 |
| R1k3 | 1589 | 746 | 0.887 | 0.469 | 0.614 | +0.002 | 0.724 | 0.645 (0.656) | 599 | 199 (+32.7%) | 0.145 | 0.370 | 95 | 6 |
| R1k5 | 1589 | 744 | 0.885 | 0.468 | 0.612 | +0.001 | 0.707 | 0.667 (0.683) | 602 | 173 (+15.3%) | 0.131 | 0.317 | 97 | 6 |
| R2_U3 | 1059 | 713 | 0.848 | 0.673 | 0.751 | +0.139 | 0.655 | 0.851 (0.998) | 259 | 143 (−4.7%) | 0.079 | 0.242 | 128 | 30 |
| R2_U5 | 861 | 669 | 0.795 | 0.777 | 0.786 | +0.175 | 0.669 | 0.862 | 121 | 139 (−7.3%) | 0.082 | 0.249 | 172 | 74 |
| R2_U10 | 764 | 656 | 0.780 | 0.859 | 0.817 | +0.206 | 0.672 | 0.868 | 48 | 136 (−9.3%) | 0.082 | 0.254 | 185 | 87 |
| R2_Jonly | 741 | 651 | 0.774 | 0.879 | **0.823** | +0.211 | 0.671 | 0.868 | 32 | 133 (−11.3%) | 0.082 | 0.256 | 190 | 92 |
| **R3 = R1k3+R2_Jonly** | 741 | 657 | 0.781 | 0.887 | **0.831** | +0.219 | 0.763 | 0.632 | 27 | 182 (+21.3%) | 0.160 | 0.430 | 184 | 94 |
| ANN (reference, self-scoring) | 841 | 841 | 1.000 | 1.000 | 1.000 | +0.388 | 1.000 | 1.000 | 0 | 185 (+23.3%) | 1.000 | 1.000 | 0 | 0 |

C1 = fraction of a matched record's transcripts with ≥ 0.999 of their exonic bases inside the node exon union; C2 = fraction whose every junction is in the node's query (rep + kept chain) junction set. Records with all transcripts contained, chimp: C1 34 → 50 (R1k3) → 47 (R3) of 841; C2 94 → 154 → 157.

## 2. DEVELOPMENT substrate (human real testis, 18 windows, 560 nodes, 166 expressed records)

| variant | nodes | matched | R | P | **F** | pair R | pair P | A1 | A2 | C1 | C2 | recs w/o node | lost vs R0 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| R0 | 560 | 139 | 0.837 | 0.248 | 0.383 | 0.569 | 0.566 | 282 | 42 | 0.132 | 0.149 | 27 | 0 |
| R1k2 | 560 | 143 | 0.861 | 0.255 | 0.394 | 0.729 | 0.333 | 265 | 103 | 0.247 | 0.322 | 23 | 4 |
| R1k3 | 560 | 143 | 0.861 | 0.255 | 0.394 | 0.714 | 0.350 | 271 | 90 | 0.233 | 0.308 | 23 | 4 |
| R1k5 | 560 | 141 | 0.849 | 0.252 | 0.388 | 0.712 | 0.370 | 275 | 72 | 0.228 | 0.288 | 25 | 4 |
| R2_U3 | 336 | 130 | 0.783 | 0.387 | 0.518 | 0.557 | 0.578 | 132 | 37 | 0.133 | 0.153 | 36 | 9 |
| R2_U5 | 259 | 123 | 0.741 | 0.475 | 0.579 | 0.563 | 0.591 | 85 | 32 | 0.132 | 0.153 | 43 | 16 |
| R2_U10 | 198 | 118 | 0.711 | 0.596 | 0.648 | 0.558 | 0.596 | 40 | 30 | 0.138 | 0.160 | 48 | 21 |
| R2_Jonly | 170 | 114 | 0.687 | 0.671 | 0.679 | 0.551 | 0.609 | 25 | 29 | 0.125 | 0.149 | 52 | 25 |
| **R3** | 170 | 117 | 0.705 | 0.688 | **0.696** | 0.754 | 0.321 | 14 | 77 | 0.273 | 0.360 | 49 | 30 |
| ANN | 166 | 166 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 0 | 47 | 1.000 | 1.000 | 0 | 0 |

### FAMILY level, development substrate only, frozen NPIP scorer, triangle leaders

FAMILY P NPIP-only (Pm) is **1.000 in every row including the rows whose family structure collapses** — vacuous by construction, never used to discriminate. Only FAMILY R, P strict and F strict do.

| arm | nodes | exon / body edges | pairs | fams | **FAMILY R** | Pm | **P strict** | **F strict** | full-length | mean copy cov | copies with a node |
|---|---|---|---|---|---|---|---|---|---|---|---|
| R0 | 560 | 1349 / 1885 | 1919 | 68 | **1.000** | 1.000 | 0.290 | 0.450 | 5/27 | 0.517 | 27 |
| R1k2 | 560 | 3664 / 4373 | 4966 | 30 | 1.000 | 1.000 | 0.153 | 0.265 | **19/27** | **0.843** | 27 |
| R1k3 | 560 | 3006 / 3889 | 4163 | 38 | 1.000 | 1.000 | 0.175 | 0.298 | 18/27 | 0.822 | 27 |
| R1k5 | 560 | 2488 / 3222 | 3348 | 47 | 0.963 | 1.000 | 0.351 | 0.515 | 17/27 | 0.767 | 27 |
| R2_U3 | 336 | 734 / 1088 | 1104 | 35 | 0.926 | 1.000 | 0.342 | 0.500 | 5/27 | 0.497 | 26 |
| R2_U5 | 259 | 533 / 817 | 830 | 29 | 0.963 | 1.000 | 0.473 | 0.634 | 5/27 | 0.494 | 26 |
| R2_U10 | 198 | 440 / 648 | 657 | 22 | 0.926 | 1.000 | 0.543 | 0.685 | 5/27 | 0.486 | 26 |
| R2_Jonly | 170 | 382 / 546 | 554 | 20 | 0.889 | 1.000 | **0.585** | **0.706** | 5/27 | 0.486 | 25 |
| **R3** | 170 | 1236 / 1288 | 1476 | 10 | 0.963 | 1.000 | 0.371 | 0.536 | **18/27** | 0.804 | 26 |
| ANN | 166 | 455 / 678 | 732 | 23 | 0.963 | 1.000 | 0.667 | 0.788 | 27/27 | 1.000 | 27 |

## 3. R3 parameter choice (fixed on dev before any held-out R3 number)

`out/R3_CHOICE.txt`, 20:07:10. Dev locus F: R1k2 0.39394, R1k3 0.39394 (exact tie at 5 dp), R1k5 0.38843 → tie broken by the declared rule (larger dev FAMILY F strict: 0.265 vs 0.298) → **k\* = 3**. Dev locus F for the floor: R2_U3 0.51793 < R2_U5 0.57882 < R2_U10 0.64835 < R2_Jonly 0.67857 → **U\* = Jonly** (junction clause alone). **R3 = R1(k=3) + R2(Jonly)**.

## 4. Decision rule, clause by clause (held-out chimp; clause v on dev)

| variant | (i) ΔF ≥ 0.05 | (ii) A2 ≤ +10% | (iii) ΔpairP ≤ 0.05 | (iv) 0 records lost | (v) FAM R and P strict not down | **closes?** |
|---|---|---|---|---|---|---|
| R1k2 | +0.005 FAIL | +56.0% FAIL | −0.219 FAIL | 8 FAIL | R 1.000 ok / Ps 0.153 FAIL | no |
| R1k3 | +0.002 FAIL | +32.7% FAIL | −0.196 FAIL | 6 FAIL | R 1.000 ok / Ps 0.175 FAIL | no |
| R1k5 | +0.001 FAIL | +15.3% FAIL | −0.174 FAIL | 6 FAIL | R 0.963 FAIL / Ps 0.351 ok | no |
| R2_U3 | +0.139 PASS | −4.7% PASS | −0.010 PASS | 30 FAIL | R 0.926 FAIL / Ps 0.342 ok | no |
| R2_U5 | +0.175 PASS | −7.3% PASS | −0.021 PASS | 74 FAIL | R 0.963 FAIL / Ps 0.473 ok | no |
| R2_U10 | +0.206 PASS | −9.3% PASS | −0.027 PASS | 87 FAIL | R 0.926 FAIL / Ps 0.543 ok | no |
| R2_Jonly | +0.211 PASS | −11.3% PASS | −0.027 PASS | 92 FAIL | R 0.889 FAIL / Ps 0.585 ok | no |
| **R3** | +0.219 PASS | +21.3% FAIL | −0.209 FAIL | 94 FAIL | R 0.963 FAIL / Ps 0.371 ok | no |

**No variant closes the gap.** Every clause is reported above for every variant on both substrates.

## 5. Where the movement actually is, and why the rule cannot be met by these arms

1. **R2's F gain is denominator shrinkage, and the rule's own clause (iv) forbids it.** On chimp the junction-only floor deletes 848 of 1,589 nodes; **572 of the deleted nodes overlap no expressed annotated record** (A1 604 → 32) — that is the 50.4%-style junk the report flagged as the cap on bipartite precision — but **114 of them were matched to a record**, and 92 expressed records end up with no node. On dev: 390 of 560 deleted, 257 overlapping no record, 36 matched. The floor is a real and cheap precision instrument, and it is the *only* thing that moves locus F; it is also, by construction, the thing clause (iv) rejects. Clauses (i) and (iv) are close to mutually exclusive on a node set where 38% (chimp) / 50% (dev) of nodes overlap nothing expressed.
2. **R1 widens the node, and about two thirds of what it adds is intronic.** Of the exonic bp R1k3 adds: chimp 1,749 kb — 29.1% annotated exon, **67.4% intronic inside a gene span**, 3.6% intergenic; dev 1,009 kb — 28.0% / 60.7% / 11.3%. Against the *expressed-record* set, chimp: 15.7% falls in the node's own matched record, 13.1% in a different expressed record (swallowing), **71.2% outside every expressed record**. That is why per-pair precision collapses (chimp median 0.998 → 0.656 — the typical node degrades, not just the tail) and why A2 rises 150 → 199. The junction-read-support anchor does not prevent it: the kept chains use junctions that reads carry but the annotation does not.
3. **R1 is nevertheless the only arm that does what the user asked.** Isoform containment roughly doubles on both substrates (C1 chimp 0.079 → 0.145/0.161, dev 0.132 → 0.233/0.247; C2 chimp 0.229 → 0.370/0.411, dev 0.149 → 0.308/0.322), per-pair recall rises 0.649 → 0.724 (chimp) and 0.569 → 0.714 (dev), and at family level NPIP copies with a full-length node go **5/27 → 18-19/27** with mean copy coverage 0.517 → 0.82-0.84 and **no copy left without a node** (FAMILY R stays 1.000 at k = 2, 3). R2 alone never moves full-length (5/27 everywhere) — node admission and node width are orthogonal defects, as the previous round said.
4. **The earlier post-hoc family headline does not survive pre-registration in this form.** The post-hoc W9k3 arm reported P strict 0.480 / F strict 0.623 with FAMILY R 0.889. The pre-registered R1 (a chain is kept iff *every junction of it* is carried by ≥ k reads) gives P strict 0.175 / F strict 0.298 at k = 3 with R 1.000. Those are different rules: W9 required a chain to **share a junction with the shipped rep chain**, i.e. to corroborate the node's own splice structure, whereas R1 only requires read support. The family-level gain therefore belongs to the *anchoring to the rep*, not to read-isoform reconstruction generically, and "read-isoform reconstruction raises F strict" must not be quoted. Monotone evidence for the same conclusion inside this run: F strict recovers as the chain filter tightens (k = 2 → 3 → 5: 0.265 → 0.298 → 0.515) and P strict rises with it, while full-length falls only 19 → 18 → 17.
5. **R3 is the best compromise measured and still fails.** Held-out F 0.831 (the highest of any arm), A1 604 → 27, containment C1 0.160 / C2 0.430, dev FAMILY F strict 0.536 vs R0's 0.450 and P strict 0.371 vs 0.290 with 18/27 full-length — but FAMILY R 0.963 (NPIPB12 loses its node entirely: copy coverage 0.477 → 0.000, deleted by the junction floor), A2 +21.3%, pair precision −0.209, 94 records lost. It is a genuine two-sided trade, not a free win.
6. **Per-copy NPIP detail (dev, R0 → R1k3 → R3)**: the copies R1 rescues are the ones the width defect was starving — NPIPB4 0.051 → 1.000, LOC128966608 0.036 → 0.982, NPIPA9 0.245 → 1.000, PKD1P6-NPIPP1 0.352 → 0.944, NPIPB6 0.375 → 0.838, NPIPB7 0.432 → 0.871, NPIPB10P 0.501 → 0.928. Four copies do not move at all (NPIPB13 0.506, NPIPB11 0.529, NPIPA8 0.680, LOC124907807 0.395) — those are not width-limited and need a different diagnosis.

## 6. Gorilla (BURNED — post-hoc W7/W9 were run there; no weight in the rule)

2,664 nodes, 139,912 reads, 1,975 expressed records. R0 F 0.786 (pair R 0.683, pair P 0.906, A1 464, A2 248, C1 0.076, C2 0.314). R1k3 F 0.788, pair P 0.759, A2 324, C1 0.134, C2 0.453, 7 records lost. R2_Jonly F 0.896, A1 34, A2 221, 157 records lost. R3 F 0.896, pair R 0.765, pair P 0.743, A2 297, C1 0.146, C2 0.506, 163 records lost. Qualitatively identical to chimp in every direction, which is the one piece of cross-substrate agreement this run can offer.

## 7. What I would do next (not run, not measured here)

- The evidence points at a **rep-anchored** chain filter (W9's clause) *plus* the junction-support floor, i.e. keep a chain only if it both has ≥ k read support on every junction **and** shares a junction with the node's shipped rep chain. That is a new variant and must be pre-registered before it is scored; it is not in this run's declarations.
- The admission floor is worth separating from this question entirely: it is the only lever on bipartite F, it is cheap, and its cost is measurable per record (114 matched nodes deleted on chimp). Whether those 114 are real loci or annotation noise is a truth question, not a metric question.
- Clause (i) as written is not reachable by any width change on either substrate; a future pre-registration should read the width objective on the per-pair and containment quantities (where R1 moves 0.079 → 0.16 and 0.229 → 0.43) and keep bipartite F as a guard, not as the primary.
## Verification (independent recompute) — agent 2 of 2

I rebuilt every substrate, every variant node set and every metric from primary inputs only (BAM, shipped node TSV, GFF, frozen `nodes.pkl`), writing my own code from the text of `DECLARATIONS.txt`. I did not read or import any script under `readisoform/scripts/`. My code: `/mnt/linuxdisk/home/juanfraitu/readisoform/verify/{vr_ri.py, vr_run.py, vr_fam.py, vr_keys.py, vr_prep_held.py}`, reusing my own prior verifier modules `vr_reads.py`, `vr_gff.py`, `vr_metrics.py` (from `locus_width/verify/`) and `npip_ladder/verify/v_edges.py` for the shared-definition edge/leader mirror. Results: `verify/vr_rows_{dev,held,burn}.json`, `verify/vr_fam_dev.json`.

### 1. Declarations precede every result; the R3 choice precedes the held-out numbers

Full mtime ordering of the workspace: `DECLARATIONS.txt` 19:56:33.343 and `DECLARATIONS.time` 19:56:33.347 (md5 `b3564799977a57bd68c135d35ac01c7b`, matches the report) come before the earliest script (`ri_build_h.py` 19:57:01) and the earliest data file (`out/sub_h.pkl` 19:57:04). `out/R3_CHOICE.txt` 20:07:10 precedes `out/locus_dev.json` 20:07:14, `out/locus_held.json` 20:07:33, `out/family_dev.json` 20:08:24, `out/locus_burn.json` 20:09:05 and the summary 20:10:53. **No held-out R3 number exists in any file written before the choice file.**

Two honest limits on that evidence: (a) the dev locus table and dev family table the tie-break is read on were *serialised* after `R3_CHOICE.txt` (20:07:14 and 20:08:24), so file mtimes alone do not prove the dev numbers existed first — only that they were computed in the same process. (b) I closed that gap by content instead: the choice reproduces exactly from my independent dev numbers. Dev locus F for R1k2 and R1k3 are *identically* equal (both 143 matched / 560 nodes / 166 records → F = 0.393939…), a genuine tie; the declared tie-break is dev FAMILY F strict, which I recompute as 0.2647 (k2) vs 0.2983 (k3) → k* = 3. Dev locus F for the floors: 0.5179 / 0.5788 / 0.6484 / 0.6786 → U* = Jonly. Choice verified.

### 2. Node sets rebuilt, diffed node by node

- **dev**: `npip_ladder/verify/nodes.pkl` arms['A4'] is byte-for-byte equal, in order, to my previously verified base (560 nodes; 550/560 have `exons == rep_exons`). R0 parity holds; importing my frozen edge mirror re-derived A4 = 1,349 exon / 1,885 body edges / 1,919 pairs, exactly the frozen figures.
- **held**: rebuilt from `PTR.dn.nodes.tsv` → 1,589 nodes (1,588/1,589 `exons == rep_exons`); reads from `PTR.1ctg.bam` with `-F 2308 -q 1` → **64,255** reads; `PTR.1ctg.gff` → **1,571** records with exons, **841** expressed at ≥3 same-strand reads, **5.5398** transcripts per expressed record. Every substrate figure in the report reproduces.
- **burn**: 2,664 nodes from `GGO.dn.nodes.tsv`, identical to my prior verified gorilla node set; 1,975 expressed records.
- **R1 monotonicity**: for every node on all three substrates, `ov_bp(R1_exons, R0_exons) == bp(R0_exons)` → R1 exons are a strict superset; 0 nodes shrink; node count never changes; nodes that grow: dev 182/151/129, held 709/603/509, burn 1525/1311/1109 at k = 2/3/5. **R1 never merges two shipped nodes into one node** — but see correction 2: R1 node footprints do fully engulf other shipped loci (held k3: 85 cases, 18 same-strand, vs 39/0 at R0).
- **R2** deletes only; surviving nodes' exon sets are bit-identical to R0 (true by construction in my reimplementation and confirmed by pair-level metrics moving ≤0.03).

### 3. Every number recomputed — 570/570 locus cells and 9/9 family rows reproduce

I recomputed all 20 reported columns for all 10 variants on all 3 substrates and diffed against `out/SUMMARY_readisoform.tsv`: **570 cells compared, 0 mismatches** (tolerance 5e-5). This includes both anti-trap columns A1 and A2, per-pair means and medians, records without a node, records with no overlap, and records lost vs R0.

Isoform-containment conventions used, as declared: **C1** = a transcript is contained iff ≥0.999 of its exonic bases lie in the node's exon union, reported as the mean over matched records of the contained fraction. **C2** = every junction of the transcript, both coordinates equal, appears in the union of the junction sets of the node's queries (shipped `rep_exons` + kept chain models); single-exon transcripts fall back to C1. Both reproduce exactly (e.g. held R0 C1 0.0794 / C2 0.2289; held R1k3 0.1445 / 0.3704; dev R0 0.1315 / 0.1494). The A3b chain-block convention (median-low 5′ start and 3′ end of the reads carrying that exact chain, cut at the chain's introns, zero-length blocks dropped) reproduces the builder's exon sets exactly — the edge counts below would not match otherwise.

Family level on dev, recomputed through the frozen `npip_ladder` edge mirror + frozen NPIP scorer, using content-addressed (md5-keyed) query alignments only: all 7,641 tx keys and 847/854 body keys I needed were already present in the captured PAFs; the 7 absent body keys are present in a submitted query FASTA with zero PAF lines (attempted, no hits), so no minimap2 run was needed. Recomputed rows (all identical to the builder's):

| variant | exon | body | pairs | fams | FAM R | Pm | **Ps** | **Fs** | full | cov |
|---|---|---|---|---|---|---|---|---|---|---|
| R0 | 1349 | 1885 | 1919 | 68 | 1.0000 | 1.000 | 0.2903 | 0.4500 | 5/27 | 0.5172 |
| R1k2 | 3664 | 4373 | 4966 | 30 | 1.0000 | 1.000 | 0.1525 | 0.2647 | 19/27 | 0.8434 |
| R1k3 | 3006 | 3889 | 4163 | 38 | 1.0000 | 1.000 | 0.1753 | 0.2983 | 18/27 | 0.8217 |
| R1k5 | 2488 | 3222 | 3348 | 47 | 0.9630 | 1.000 | 0.3514 | 0.5149 | 17/27 | 0.7674 |
| R2_U3 | 734 | 1088 | 1104 | 35 | 0.9259 | 1.000 | 0.3425 | 0.5000 | 5/27 | 0.4970 |
| R2_U5 | 533 | 817 | 830 | 29 | 0.9630 | 1.000 | 0.4727 | 0.6341 | 5/27 | 0.4937 |
| R2_U10 | 440 | 648 | 657 | 22 | 0.9259 | 1.000 | 0.5435 | 0.6849 | 5/27 | 0.4860 |
| R2_Jonly | 382 | 546 | 554 | 20 | 0.8889 | 1.000 | 0.5854 | 0.7059 | 5/27 | 0.4858 |
| R3 | 1236 | 1288 | 1476 | 10 | 0.9630 | 1.000 | 0.3714 | 0.5361 | 18/27 | 0.8041 |
| ANN_reference | 455 | 678 | 732 | 23 | 0.9630 | 1.000 | 0.6667 | 0.7879 | 27/27 | 1.0000 |

FAMILY Pm is 1.000 in every row, as the declaration says it must be; it is reported and not used.

### 4. Decision rule, clause by clause, on the held-out substrate

Baseline R0(held): F 0.6115, A2 150, per-pair precision 0.8411, 0 records lost. Limits: (i) ΔF ≥ +0.05, (ii) A2 ≤ 1.10×, (iii) per-pair precision drop ≤ 0.05, (iv) 0 records lost, (v) dev FAMILY R and P strict do not fall.

| variant | ΔF | i | A2 (rel, % of nodes) | ii | Δ pairP | iii | lost | iv | dev R / Ps | v | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|
| R1k2 | +0.0050 | N | 234 (1.560×, 14.7%) | N | −0.2187 | N | 8 | N | 1.000 / 0.153 | N | FAIL |
| R1k3 | +0.0025 | N | 199 (1.327×, 12.5%) | N | −0.1959 | N | 6 | N | 1.000 / 0.175 | N | FAIL |
| R1k5 | +0.0008 | N | 173 (1.153×, 10.9%) | N | −0.1739 | N | 6 | N | 0.963 / 0.351 | N | FAIL |
| R2_U3 | +0.1390 | Y | 143 (0.953×, 13.5%) | Y | +0.0104 | Y | 30 | N | 0.926 / 0.342 | N | FAIL |
| R2_U5 | +0.1746 | Y | 139 (0.927×, 16.1%) | Y | +0.0205 | Y | 74 | N | 0.963 / 0.473 | N | FAIL |
| R2_U10 | +0.2059 | Y | 136 (0.907×, 17.8%) | Y | +0.0269 | Y | 87 | N | 0.926 / 0.543 | N | FAIL |
| R2_Jonly | +0.2115 | Y | 133 (0.887×, 17.9%) | Y | +0.0266 | Y | 92 | N | 0.889 / 0.585 | N | FAIL |
| R3 | +0.2191 | Y | 182 (1.213×, 24.6%) | N | −0.2095 | N | 94 | N | 0.963 / 0.371 | N | FAIL |
| **ANN_reference** | +0.3885 | Y | 185 (1.233×, 22.0%) | **N** | +0.1589 | Y | 0 | Y | 0.963 / 0.667 | **N** | **FAIL** |

**The builder's verdict — no variant closes the gap — is confirmed clause by clause.** The last row is the correction above: the annotation itself fails the same rule on clauses (ii) and (v).

Note on clause (ii) as written: it compares raw A2 counts, so a variant that deletes half the nodes passes it while nearly doubling the *rate* of neighbour-swallowing nodes (R2_Jonly 9.4% → 17.9% of nodes overlap ≥2 expressed records). Both forms are in the table.

### 5. Is the R2/R3 gain denominator shrinkage? Yes, ~87% of it

R2 cannot improve a node — it only removes nodes. Direct test: take R0's own matching and restrict it to the 741 nodes that survive R2_Jonly on held. That alone gives F 0.7952 (P 0.8489, R 0.7479), i.e. **+0.184 of the +0.211 total F gain (87%) is the denominator**; the residual +0.028 comes from records freed by deleted nodes re-matching to survivors in the bipartite step, not from better exon sets. Dev: restricted-R0 F 0.6131 of R2_Jonly's 0.6786, against an R0 baseline of 0.3829 (78% of the gain is shrinkage). Corroborating: per-pair recall 0.649→0.672, per-pair precision 0.841→0.868, C1 0.079→0.082 — the surviving nodes are not measurably better, there are simply fewer bad ones.

What is deleted is not all noise: R2_Jonly removes 848 held nodes, of which **572 overlapped no expressed record but 276 did**, and the deleted set includes a node with 1,247 assigned reads (median 2, mean 3.9). On dev it removes 390 nodes, 133 of which overlap an expressed record.

By contrast R1's containment gain is not a denominator effect — node count is fixed and widening is monotone: held C1 0.079→0.145, C2 0.229→0.370 (k=3); dev NPIP full-length 5/27→18/27 and mean copy coverage 0.517→0.822. The cost is equally real and is not denominator-driven either: dev FAMILY P strict 0.290→0.175 with exon edges 1,349→3,006 and families 68→38.

### 6. No annotation enters node construction

My reimplementation of R0/R1/R2/R3 consumes only: the shipped node TSV (or the frozen A4 node set), the BAM reads, and the constants k and U. It reproduces all 570 builder locus cells and all 9 builder family rows exactly — which is strong evidence that the builder's construction is likewise annotation-free, since any annotation leak would have to be reproduced by an annotation-free reimplementation. The annotation enters only in the evaluation denominators (expressed records, identical across all rows of a substrate: 166 / 841 / 1,975) and in the ANN_reference row, which is labelled as a reference and not a construction. The read→node assignment is computed once from the R0 exons and shared by every variant, so no variant gains reads by growing; I verified the assigned-read totals are the same across variants (held 60,095; dev 49,927).

### 7. Spot-checks of the report's specific claims

All verified against my own recompute: held F 0.612 → 0.823 for the junction-only floor (+0.211); 92 expressed records lose their node under it; dev FAMILY R 1.000 → 0.889; R1 moves held F by only +0.005 / +0.002 / +0.001 at k = 2/3/5 while A2 rises 56.0% / 32.7% / 15.3% and per-pair precision falls 0.219 / 0.196 / 0.174; R3 has the best held F (0.8306) and fails (ii), (iii), (iv), (v); NPIP full-length 5/27 → 18–19/27 with mean copy coverage 0.517 → 0.82–0.84 and dev P strict 0.290 → 0.175 at k = 3. The gorilla rows are reproduced too and are correctly marked BURNED.