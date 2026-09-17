# TBC1D3 subfamily truth corrected from Guitart 2024 Fig 6B/6C

2026-09-16. This is a **post-hoc truth amendment** and must be disclosed as one. Nothing was committed and no existing file was edited.

**New truth table:** `docs/lit_subclusters_tbc1d3_guitart_fig6c.tsv` (9 rows, 1-based GFF coordinates; the old truth's `start0` equals `start` − 1).
- Generator and assertions: `WC/make_fig6c_tsv.py`.
- A drop-in copy in the old column format: `R00/truth_guitart_fig.tsv` (md5 e3411e2e).

**Path shorthands**
- Directories: `TG` = `/mnt/linuxdisk/home/juanfraitu/tbc1d3_guitart_truth`; `WC` = `TG/writeup_checks`; `R0n` = `TG/rescore/0n`; `RA` = `TG/rescore_A`; `RB` = `TG/rescore_B`; `AU` = `TG/audit_lean`; `LIT` = `/mnt/linuxdisk/home/juanfraitu/o1_falsemerge/lit`.
- Output files: `glo` = `RA/guided_lo/rescore_guided_lo.out`; `gt` = `RA/guided_t/run_FS.out`; `d2` = `RA/denovo_d2/rescore_d2.out`; `ph` = `RA/phap/rescore_phap.out`; `lo` = `RB/inv22_layer_order/rescore_layer_order.out`; `lat` = `RB/inv23_lattice/rescore_lattice_truth.out`; `snap` = `WC/correction_state_snapshot.json`.

**Copy letters:** B, I, G, H, F, E, K and D are the RefSeq TBC1D3B…TBC1D3D. **T** is unsuffixed TBC1D3.

---

## §0 Summary

### Old vs corrected mapping

| CHM13 RefSeq copy | cluster | old level 2 (by name) | Guitart group (6C label) | sequence | S variant |
|---|---|---|---|---|---|
| TBC1D3B | 1 | B | **M** (M) | supported | M |
| TBC1D3I | 1 | I | I (I) | supported | I |
| TBC1D3G | 1 | G | B (B) | unresolved | B |
| TBC1D3H | 1 | H | **M** (M2) | supported | M |
| TBC1D3F | 1 | F | G (G) | unresolved | G |
| TBC1D3E | 2 | **AE** | O (O) | unresolved | O |
| TBC1D3K | 2 | **CDKL** | **CDKL** (CDKL) | supported, moderate | CDKL |
| TBC1D3D | 2 | **CDKL** | AE (AE, by colour) | **unresolved; leans CDKL** | excluded |
| TBC1D3 (T) | 2 | **AE** | **CDKL** (CDKL2) | supported | CDKL |

- **Multi-copy groups.** Old: AE={T,E} and CDKL={D,K}. Corrected (F): **M={B,H}** and **CDKL={K,T}**; the other 5 copies are singletons. The old and new pair sets share no pair.
- **Scale of the change.** 7/9 level-2 letters change, and ARI(old, F) = −0.059. Level 1 is unchanged (ARI 1.000) (`R00/rescore_truth.out:61,64`).
- **Letters are not RefSeq symbols.** RefSeq G is group B and RefSeq F is group G. Key the truth by coordinates only.

### Verification

Both verifiers returned **CONFIRMED_WITH_CAVEATS** (figure: `snap:5`; sequence: `snap:75`).

**Figure** (vector geometry, PDF p.45)
- The CHM13 track has 9 arrows for the 9 RefSeq copies, and all 9 strands match.
- For every arrow, the nearest legend colour is at dE 0.3-5.1 and the second nearest at ≥ 21.5 (`R00/fig6b_check.out:22-30`). An independent re-read gives the same map (`:40`).
- **Caveats**
  - Three cluster-2 arrows look **hand-placed**: O on E, CDKL on K, and CDKL2 on T. They miss the gene-end convention by +16.6, −3.1 and −11.1 kb (`snap:66`). The E→O assignment therefore rests on strand and order, not on position.
  - D's AE arrow is a **data-driven** one: its left edge sits 0.3 kb past the gene end, as in cluster 1.
  - 6C is a schematic with staggered labels, so colour (not label position) is what assigns AE.

**Sequence** (36 CHM13 pairs with minimap2 asm20, 11 GRCh38 anchors, population trees of 369 and 300 copies)
- **M={B,H} supported.** B and H are mutual nearest copies at p 0.00239 (`TG/verify_sequence/pairs_asm20.tsv:4`), and they form a population clade with no GRCh38 copy.
- **T in CDKL supported.** T is 0.00092 from GRCh38 TBC1D3L (`pairs21_asm20.tsv:141`), the closest CHM13-GRCh38 match of all.
- **Old pairings refuted.** Old AE={T,E} is refuted (0.00248, `pairs_asm20.tsv:34`), and so are the old B and H singletons.
- **Old {D,K} not refuted.** D and K are mutual nearest at 0.00135 (`:35`). D's nearest GRCh38 copies are the CDKL members L (0.00147) and D (0.00151) (`pairs21_asm20.tsv:129-130`).
- **G, F and E unresolved.** Their letters do not change any CHM13 partition, because each is a singleton either way.

**Text anchor** (neither verifier used it): "for 67 of the 69 assembled haplotypes, this expressed [CDKL] paralog is the last copy in cluster 2" (`WC/guitart_2024_biorxiv_pages.txt:484`). T is the last cluster-2 copy.

### Robustness tiers and the three truth variants

- **Robust:** M={B,H}, and T in CDKL. Figure, sequence and text agree.
- **Moderate:** K in CDKL. K's arrow is hand-placed. The smallest population-tree clade holding K and T has 115/100 tips (support 38/40), against 42/90 tips for D+K (`snap:133`).
- **Unresolved:** D (figure AE, sequence CDKL).

Every claim is scored under three truths:
- **F:** figure, all 9 copies.
- **S:** F with D excluded. S is not "sequence-safe", because it keeps K in CDKL.
- **V3:** F with D moved into CDKL.

### Net effect over the 33 inventory items

| action | items | meaning |
|---|---|---|
| **retract** | 16 | the claim, or part of it, fails under F, S and V3 |
| **suspend** | 5 (only), 10 (incl. partial) | the claim fails under F but holds under V3; do not quote it either way until D is resolved |
| **relabel** | 11 | labels or numbers change; the conclusion stands |
| **unaffected** | 1 | |
| **heavy** | 0 | no item needed a new tree or alignment |

**Verdicts that stand**
- Addendum AA / register row 822 is **NOT SUPPORTED** under every variant. IQ-TREE vs DT "either"-class recoveries:

  | truth | IQ-TREE | DT | source |
  |---|---|---|---|
  | OLD | 85 | 73 | `gt:40` |
  | F | 100 | 84 | `gt:81` |
  | S | 97 | 84 | `gt:122` |
  | V3 | 94 | 84 | `WC/aa_v3.out:5` |
  | M only (CDKL dropped) | 94 | 84 | `WC/aa_v3.out:11` |
  | TBC1D3 groups removed | 84 | 73 | `R01/aa_rescore.out:290` |

- Every positional cluster1|cluster2 call.
- TBC1D3 family-level 1.000/1.000.
- The R7 and F4 selections.

**The one robust new positive: M={B,H}**
- CHM13 reference trees: intron 100/100 (the MAFFT and projection trees share one input) and exon 77.7/71.
- Leave-out runs recover it in 9-10 of 10.
- The HG002 and gorilla panel trees recover it.
- It is exactly the root split of the guided identity UPGMA.
- Its calls are the same under F and V3, because D carries a label in both.

**Figure-dependent only: CDKL={K,T}**
- On the CHM13 reference it has SH-aLRT 83.5-84.2 and UFBoot 76-77, all from one intronic dataset.
- "RECOVERED" checks SH only (`bench/guided_pipeline.py:535`).
- Under S the registered rule calls it "unsupported". Under V3 it fails (0.0/39-43).

---

## §1 Derivation

1. **Old truth.** `docs/lit_subclusters_npip_tbc1d3_truth.tsv:24-32` (md5 a79fabc1) has three byte-identical copies: `LIT/lit_truth.tsv`, `LIT/guided_lo/truth.tsv` and `LIT/guided_t/truth.tsv`. The name map was pre-registered as a "declared assumption" (`docs/PREREG_core_definition_2026-09-12.md:223-225`).
2. **Figure reading** (Guitart et al., bioRxiv 2024.03.12.584650, posted 2024-03-13, PDF p.45).
   - Words and vector drawings were extracted with pymupdf.
   - Legend colours come from the Panel A key. The 6B axis was calibrated from the tick strokes (91.7 pt/Mb).
   - Each arrow's direction comes from its apex vertex, and its group from a CIE Lab colour match.
   - Arrows were matched to copies 1:1 in order. minimap2 asm20 of RefSeq T and B finds exactly 9 copies of about 10.9 kb in chr17:36.9-39.7 Mb (`TG/verify_figure/hits.paf`).
   - Two independent readings agree (verifier 1; `R00/fig6b_check.py`).
3. **Sequence.**
   - CHM13 pair distances: `TG/verify_sequence/pairs_asm20.tsv`. CHM13-vs-GRCh38 distances: `pairs21_asm20.tsv`.
   - Population trees `tree_pool` (369 copies) and `tree_c12` (300 copies) were built with IQ-TREE GTR+F+G4 and rooted on chimp.
   - Record fix: verifier 2's report swaps the §6js reference-tree site counts (`snap:134`). The files say exon 1,833 and intron 8,645 sites (`LIT/guided_t/tree_t/ref_TBC1D3_{exon,intron}.iqtree`, "Input data" line). No call is affected.
4. **Rescoring rule.**
   - Every item reused the original script or its printed outputs.
   - Each rescore first reproduced the original numbers under the OLD truth (the gate passed for every gated item) before scoring F and S.
   - V3 numbers come from the audit's independent recomputation (`AU/r1…r8`) and from `WC/aa_v3.py`.
   - No tree or alignment was rebuilt.
   - `clade_calls` lets the last matching split win and checks SH only (`bench/guided_pipeline.py:521-543`). "Any-match" is reported as a sensitivity.
5. **How independent the reference-tree results are.** They are one dataset, not several confirmations.
   - All four CHM13 reference trees read one input FASTA (md5 0ce19224): `LIT/guided_lo/tree/ref_TBC1D3.fa`, `tree_proj/ref_TBC1D3.fa`, `tree_union/ref_TBC1D3_intron.fa` and `guided_t/tree_t/ref_TBC1D3_intron.fa`.
   - Projection, §6jr intron and §6js intron share one alignment (`*.proj.fa` md5 72577a4c, 8,645 sites).
   - `tree_proj` and `tree_union` intron treefiles are byte-identical (md5 92abc5bb).
   - So "§6jp projection 83.5/77", "§6jr intron 84/77" and "§6js intron 83.5/77" are **one result**, and MAFFT (8,644 sites, 84.2/76) is the same data re-aligned.
   - Only the exon trees are separate data: §6js exon (1,833 sites) and §6jr exon (2,152 sites).
   - CDKL UFBoot is below 80 on every CHM13 reference tree, 56-82 in leave-out runs (`AU/r2_trees/r2.out:90,138,154,218,234`) and 62-63 in de novo trees (`AU/r7_d2/r7.out:38,48`).

---

## §2 Impact of the correction on all 33 inventory items

Each row gives the item's location, the old claim, the result under F, S and V3, and the action. Unless a line says otherwise, line numbers are current `docs/o1_ledger.md` lines.

### inv01 — truth TSV
**Where:** `docs/lit_subclusters_npip_tbc1d3_truth.tsv:24-32`, plus 3 copies.
**Old claim:** AE={T,E}, CDKL={D,K}, other copies singletons.

| truth | result |
|---|---|
| F | M={B,H}, CDKL={K,T}; ARI(old, F) −0.059 (`R00/rescore_truth.out:64`) |
| S | F without D |
| V3 | adds pairs {D,K} and {D,T}; ARI(old, V3) 0.280 (`:92`) |

**RELABEL.** New TSV written. {T,E} and the B/H singletons are refuted by sequence; {D,K} is not.

### inv02 — core prereg
**Where:** `PREREG_core…:223-225, 579, 671, 772`; AA at 960-975.
**Old claim:** the name map, and AA totals that include the TBC1D3 groups.

| truth | result |
|---|---|
| F | AA IQ-TREE 100 vs DT 84 (`gt:81`) |
| S | 97 vs 84 (`gt:122`) |
| V3 | 94 vs 84 (`WC/aa_v3.out:5`) |

**RELABEL** through a new addendum. The AA verdict is unaffected.

### inv03 — known-subclusters prereg
**Where:** `PREREG_known_subclusters…:22, 29, 36` (P2-P4).
**Old claim:** "within human, the copies are ONE clade"; P3's {D,K} cut is a sister pair.

| truth | result |
|---|---|
| F | P2/P3/P4 literal verdicts unchanged. Every group nests inside a cluster (`RA/inv04/run.stdout:45`; paper l.405-406). {D,K} crosses AE/CDKL. |
| S | groups nest; {D,K} untestable (`:47`, `:61`) |
| V3 | groups nest (`AU/r6_6gw/r6.out:22`); {D,K} is within CDKL |

**RETRACT** "ONE clade". **SUSPEND** P3's reading. The verdicts are unaffected.

### inv04 — §6gw and register row 778
**Where:** ledger 14928-14934; `NEGATIVE_RESULTS_REGISTER.md:1387`.
**Old claim:** no bimodality; D/K is the top pair (0.9986), "a recent tandem sister pair"; the method "found the strongest real sequence signal"; "no sequence signal to find".

| truth | result |
|---|---|
| F | Table and detectors reproduce (`RA/inv04/run.stdout:2-4,16`). Within-group identity p 0.0688 (`:10`). The I=5.0 cut scores sensitivity 0.500, precision 0.045, ARI −0.021 (`:31`). |
| S | p 0.0238 (`:12`); I=5.0 ARI −0.050 (`:37`) |
| V3 | p 0.0040; I=5.0 ARI −0.042 (`AU/r6_6gw/r6.out:18,20`) |

**RETRACT** "no signal to find". **SUSPEND** "{D,K} is the real signal". The I=5.0 cut breaks CDKL in every variant, but that rests on K being in CDKL, which is only moderate.

### inv05 — §6jg (done)
**Where:** 18891, 18909-18911, 18917.
**Old claim:** "CDKL={D,K} recovered exactly", ARI 0.654; de novo 0.364; "the tightest published group appears".

| truth | result |
|---|---|
| F | Guided cut ARI −0.038, 0/2 groups exact (`R04/variants.out:104`); de novo −0.061 (`:108`) |
| S | 0.000 with no predicted pairs (`:152`); de novo −0.068 (`:157`) |
| V3 | 0.372, sensitivity 0.25, 0/2 exact (`AU/r1_upgma/r1.out:56`) |

**RETRACT.** New observation: the guided root split {B,H}|rest equals M exactly under F, S and V3 (`r1.out:53,55,57`).

### inv06 — §6jh level 2 (done)
**Where:** 18971-18977, 18985, 18993-18995.
**Old claim:** de novo identity 7/7, "TBC1D3 fully by de novo identity"; guided coverage × identity 5/7; "is D-K the closest pair"; gene spans 7/7.

| truth | result |
|---|---|
| F | De novo identity 4/7; guided coverage × identity 4/7; guided identity 3/7; the best is de novo coverage at 5/7. Every pairwise sensitivity is 0 (`R05/summary_L2_tbc1d3.tsv:21-27`). On the shared 8 copies every unit type gives 4/7 (`:35-37`). |
| S | De novo identity 3/6, guided 2/6 (`RB/inv29_denovo_gap/cov_S.out:12,29`); all unit types 3/6 (`diag_S.out:28,34,40`) |
| V3 | Fails. V3 needs a B-H merge, and B-H ranks 6/28 (`R05/diag_fig6.out:28`; `AU/r8…out:15`) |

**RETRACT.** "Not a mode effect" still holds.

### inv07 — §6ji level B
**Where:** 19042.
**Old claim:** 5/5 on the shared 6 copies for R0, R2, R3, R7 and gene spans.

| truth | result |
|---|---|
| F | All 6 copies are singletons: trivial 6/6 (`RA/inv07/rescore_levelB_FS.out:12,61,75`) |
| S | trivial 5/5 (`:13,62,76`) |
| V3 | The restricted truth equals OLD ({D,K} plus singletons), so 5/5 stands (`AU/r8_restricted_partitions.out:1-5`) |

**SUSPEND.** Selection R7 is unaffected (`:82,84,86`).

### inv08 — lit_modes_j level B
**Where:** `LIT/lit_modes_j.out:108-116` (never quoted).
**Old claim:** 6/6 on the shared 7 copies for R0, F4 and gene spans.

| truth | result |
|---|---|
| F | trivial 7/7 (`RA/inv08/rescore_modes_j_FS.out:5,17,29`) |
| S | trivial 6/6 (`:6,18,30`) |
| V3 | equals OLD (`r8…out:6-10`) |

**SUSPEND.** The J1 candidate F4 is unaffected (`:36,41,46`).

### inv09 — §6jo G2
**Where:** 19397-19399.
**Old claim:** G2 no split is "correct"; one leave-out run splits {H,K}.

| truth | result |
|---|---|
| F | No-split runs: pair 1.000/0.056, bipartite 0.222 (`glo:354`). half_3 {H,K}: 0.000/0.000 (`:366`). |
| S | 1.000/0.071, bipartite 0.250; half_3 0.000 (`:355,367`) |
| V3 | {H,K} still crosses M and CDKL |

**RELABEL.** Positional "correct" is unaffected. G2 recovers no Guitart group in 11/11 runs.

### inv10 — §6jp
**Where:** 19438-19440, 19453-19456, 19476.
**Old claim:** AE 1/10 (MAFFT) and 0/10 (projection); CDKL 2/10 and 3/10; no split on the reference.

| truth | result |
|---|---|
| F | M: reference 100/100 in both; leave-out 10/10 MAFFT, 9/10 projection (`glo:56,111,228,254`). CDKL: reference 84.2/76 and 83.5/77; leave-out 2/10 and 4/10 (`glo:227,253`). |
| S | M as F. CDKL: leave-out 3/10 and 4/10; reference "unsupported" by the last-match rule, RECOVERED by any-match (`glo:58-59,113-114,243,269`). |
| V3 | M as F. CDKL reference unsupported, 0.0/43 and 0.0/39 (`AU/r2_trees/r2.out:10,23`). |

**RETRACT** the AE/CDKL cells. M is a robust positive; CDKL holds under F only. The positional calls are unaffected (`glo:229,232,255,258`). Both methods use one input FASTA (§1.5).

### inv11 — §6jq table cell
**Where:** 19522.
**Old claim:** "none stable (positional)".

| truth | result |
|---|---|
| F | M stable; CDKL on the reference only |
| S | M stable; CDKL by any-match only |
| V3 | M stable; CDKL not recovered |

**RETRACT** the cell.

### inv12 — §6jr
**Where:** 19580, 19591-19594.
**Old claim:** AE/CDKL no/no; intron leave-out 0/3; "coding shows none".

| truth | result |
|---|---|
| F | S2 exon alignment (2,152 sites): no split on the reference, 0/10 leave-out (`glo:299-303`). Intron: M 100/100 and CDKL 84/77 on the reference (`:221`); leave-out M 9/10, CDKL 4/10 (`:325-326`). |
| S | CDKL reference intron: 0/39 by the last-match rule, 84/77 by any-match (`glo:223-224`); leave-out 4/10 (`:341`) |
| V3 | M as F; CDKL reference unsupported (`r2.out:49`) |

**RETRACT** the cells. The positional reading is unaffected. Do **not** read the groups as intron-only: the §6js exon alignment recovers M at 77.7/71 (`r2.out:57`).

### inv13 — §6js B2 cell
**Where:** 19626.
**Old claim:** "AE/CDKL not literature clades (0 before too)".

| truth | result |
|---|---|
| F | M: reference exon 77.7/71 and intron 100/100. Leave-out either class: IQ-TREE 9/10 (exon 4, intron 8); DT 10/10 plus the reference. CDKL: IQ-TREE intron reference 83.5/77, leave-out 5/10; DT 0 (`gt:48-75`; `RA/guided_t/ref_supports.out:2,5`). |
| S | M as F. CDKL: 3/10 registered, 5/10 any-match; reference unsupported (`gt:90,110,175`). |
| V3 | M 9/10 strict; CDKL 0/10 (`r2.out:266-273`) |

**RETRACT** the cell text. The bar verdicts are unaffected: B2 excludes TBC1D3, and the B3 exon tree is never positional (`gt:60-61,76-77`).

### inv14 — §6jt D2
**Where:** 19843-19844.
**Old claim:** the DN0 and DN1 cells include CDKL.

| truth | result |
|---|---|
| F | CDKL absent (no T copy); DN1 intron M unsupported (`d2:9,15`) |
| S | same as F (`d2:10,16`) |
| V3 | CDKL RECOVERED as the partial pair {D,K}: DN0 83.6/73, DN1 82.0/75 (`AU/r7_d2/r7.out:10,20`) |

**SUSPEND.** The positional statement is unaffected (`d2:8,14`).

### inv15 — §6jx and register row 822
**Where:** 19899-19907; register `:1433`.
**Old claim:** TBC1D3 AE/CDKL IQ-TREE 0 / 1 of 10, DT 0/0; totals 85 vs 73, conflicts 72 vs 49; NOT SUPPORTED.

| truth | result |
|---|---|
| F | M: IQ-TREE 9/10 plus reference, DT 10/10 plus reference. CDKL: IQ-TREE 5/10 plus reference, DT 0. Totals 100 vs 84; conflicts 67 vs 47 (`gt:48-82`). |
| S | CDKL: IQ-TREE 3/10, reference unsupported. Totals 97 vs 84; conflicts 60 vs 44 (`gt:89-123`). |
| V3 | Totals 94 vs 84; conflicts 62 vs 44 (`WC/aa_v3.out:4-6`) |

**RELABEL** the TBC1D3 row. The verdict is unaffected: M-only gives 94/84 (`WC/aa_v3.out:11`), and any-match under F and S gives 101/84 (`gt:166,207`).

### inv16 — §6kp
**Where:** 20630-20634, 20641-20643.
**Old claim:** the bridge cut makes AE/CDKL evaluable, and they are "not recovered".

| truth | result |
|---|---|
| F | CDKL RECOVERED in the intron trees of DN1r+B (85.2/62) and DN0+B (85.0/63). M unsupported (59.4/58) (`d2:39,45`; `r7.out:38,48`). CDKL absent in DN0, DN1 and DN1r. |
| S | same as F (`d2:40,46`) |
| V3 | DN0+B and DN1r+B: CDKL no split. DN0, DN1, DN1r: recovered (`r7.out:10,20,30,40,50`). The original reading comes back. |

**SUSPEND.** "No de novo catalog recovers M" holds under F, S and V3.

### inv17 — §6kq
**Where:** 20862 (heading), 20880, 20882, 20885-20886, 20891-20892, 20913, 20919-20921.
**Old claims:**
- HG002 panel: "CDKL recovered (not on CHM13 alone)", AE no split.
- Gorilla panel: AE and CDKL no split.
- Single linkage: {D,K} plus 12 HG002 copies is "CDKL as a population group".
- 23 haplotypes: CDKL recovered.
- AE is never recovered.

| truth | result |
|---|---|
| F | Human panel: M RECOVERED (97.7/100, 91.8/98, 100/100), CDKL no split (`ph:2`). Gorilla panel: M 100/100, CDKL 76.2/81 (`ph:5`). HPRC: M unsupported by the last-match rule, RECOVERED by any-match (99.8/100); CDKL no split (`ph:8`). Human single linkage forms neither group (`ph:15`). CHM13 alone recovers CDKL (`glo:56,111`). |
| S | Tree calls as F (`ph:3,6,9`). CHM13-alone CDKL by any-match only. |
| V3 | CDKL: human no split, gorilla unsupported 0.0/68, HPRC no split (`AU/r3_phap/r3.out:9,21,33`). CHM13 reference unsupported, leave-out 0/10 (`r2.out:10,23,266`). |

**RETRACT:**
- "The HG002 and 23-haplotype panels recover CDKL".
- "{D,K} plus 12 copies is the CDKL population group". T sits outside that component (`ph:13`), so the group fails under every variant.
- "AE is never recovered". This is untestable: AE has at most 1 CHM13 copy.

**SUSPEND** "not recoverable from CHM13 alone" and the gorilla CDKL call.

**New positive:** M, in the human and gorilla panels.

### inv18 — §6ky Q9 summary
**Where:** 21450-21453.
**Old claim:** "TBC1D3-CDKL recovered, TBC1D3-AE never, even at 23 haplotypes".
**F, S, V3:** as inv17.
**RETRACT** the clause. Family-level 1.000/1.000 is unaffected.

### inv19 — §6kz framing
**Where:** 21494-21497.
**Old claim:** the Guitart groups are "literature facts".
**RELABEL** the wording. This item is not rescorable.

### inv20 — ADVISOR_QUESTIONS
**Where:** `ADVISOR_QUESTIONS.md:775-776`.
**Old claim:** "passes every bar … TBC1D3-CDKL is recovered exactly (0.9944 median divergence)".

| truth | result |
|---|---|
| F | The bars stand. {D,K} crosses groups (ARI −0.038). 0.9944 is the median gene-span identity (`LIT/lit_analysis.out:48`), not a divergence. |
| S | ARI 0.000 |
| V3 | ARI 0.372, sensitivity 0.25 |

**RETRACT** the CDKL sentence. The bars are unaffected.

### inv21 — ADVISOR_QUESTIONS
**Where:** `ADVISOR_QUESTIONS.md:810-812, 933`.
**Old claim:** a second haplotype recovers CDKL, which is "not recoverable from CHM13 alone"; "AE never recovered".

| truth | result |
|---|---|
| F | HG002 and HPRC: CDKL no split. CHM13 alone recovers it (`glo:56,111`). |
| S | CHM13 alone: any-match only |
| V3 | CHM13 alone: not recovered. Panels: no split. |

**RETRACT** "the second haplotype recovers CDKL" and "AE never". **SUSPEND** "not from CHM13 alone".

### inv22 — LAYER_ORDER report
**Where:** `bench/LAYER_ORDER_NPIP_TBC1D3.md:18, 19, 75, 131-135, 155, 251, 253-254, 263, 388`.
**Old claims:**
- Literature C has 0 TBC1D3 pairs.
- C_L1 ⊇ C_fine holds by construction.
- "AE and CDKL have no supported split".
- C_tree_top 0.444 and C_tree_min 0.667 for TBC1D3.

| truth | result |
|---|---|
| F | Literature C = M (exon 77.7/71, intron 100/100; compatible) plus CDKL (intron 83.5/77; conflicts with {E,K}) (`lo:8-10,79-80`). 52 tournament rows and 18 verdicts change, all TBC1D3; no pooled verdict changes (`lo:71`). C_L1 ⊇ C_fine is 0/2 for TBC1D3. C_tree_top 0.556, C_tree_min 0.778 (`lo:93,99`). |
| S | Literature C as F. C_tree_top 0.625, C_tree_min 0.875 (`lo:94,100`). |
| V3 | Literature C = M only (`r2.out:61-62,74-75`). C_tree_top 0.667, C_tree_min 0.889 (`AU/r5_lattice/r5.out:11-12`). |

**RETRACT** "0 pairs", "no supported split" and "by construction" for TBC1D3. **RELABEL** the numbers. "P and D split 0 clade pairs" is unaffected (`lo:114-115`).

### inv23 — NESTED_LATTICE report
**Where:** `bench/NESTED_LATTICE_NPIP_TBC1D3.md:442, 684-685`.
**Old claim:** lattice F 0.222 at every level; report-only C_tree_min 0.667 and C_tree_top 0.444; L2 groups AE (T,E) and CDKL (D,K).

| truth | result |
|---|---|
| F | Lattice 0.222 (`lat:38`); C_tree_min 0.778, C_tree_top 0.556 |
| S | Lattice 0.250 (`lat:63`); 0.875 / 0.625 (`lat:68-69`) |
| V3 | Lattice 0.333; 0.889 / 0.667 (`r5.out:10-12`) |

**RELABEL.** "The lattice does not resolve TBC1D3 subgroups" holds: all 9 copies form one group at every level (`lat:80-86`).

### inv24 — code: guided_pipeline
**Where:** `bench/guided_pipeline.py:509`.
**Old behaviour:** `literature_groups()` hard-codes ("AE","CDKL").

| truth | result |
|---|---|
| F | Returns AE={D} and CDKL={K,T}; **M is silently dropped** (`RB/code_checks/code_checks.out:5`) |
| S | AE={}; M dropped (`:6`) |
| V3 | same defect |

**RELABEL** the code (fix in §6).

### inv25 — code: guided_tree_width
**Where:** `bench/guided_tree_width.py:33-34`.
**Old behaviour:** GROUPS hard-codes AE/CDKL.
**F:** M dropped (`code_checks.out:8`). **S:** M dropped (`:9`). **V3:** same defect.
**RELABEL** the code.

### inv26 — code: lit scripts
**Where:** `LIT/lit_tbc_diag.py:33`; lit_analysis, lit_coverage, lit_modes and lit_modes_j read `LIT/lit_truth.tsv`.
**Old behaviour:** hard-coded truth pairs (D,K) and (T,E).
**F:** neither is a truth pair (`code_checks.out:18`). **S:** same (`:19`). **V3:** (D,K) is a truth pair.
**RELABEL** the code.

### inv27 — code: layer-order
**Where:** `truths_universe.py:211-212`; `lattice_truth.py:13`; `lo_*`; `lattice_edges.py:738-754`.
**Old behaviour:** annotations say "mapped by name (AE = TBC1D3+E, CDKL = D+K)" and "positional AE/CDKL".
**F, S, V3:** the propagation code is truth-agnostic; inv22 and inv23 reproduce the originals through it. Only `layer_c.py`, through `literature_groups`, drops M.
**RELABEL** the annotation text.

### inv28 — memory: population unit
**Where:** `project_tbc1d3_subclusters_population_unit.md:14, 21-24`; `MEMORY.md:83`.
**Old claim:** "AE and CDKL each merge 2 genes"; §6kq CDKL recovered (not from CHM13 alone), AE never; 23 haplotypes still give CDKL.
**F, S, V3:** as inv17. CDKL merges 4 GRCh38 genes (C, D, K, L); AE merges 2 (A, E) (`snap:108`).
**RETRACT** the panel-CDKL claims and "each merge 2". **SUSPEND** "not from CHM13 alone".

### inv29 — memory: de novo gap
**Where:** `project_denovo_vs_annotated_gap.md:391-393, 399-400, 403-404`.
**Old claim:** as inv05 and inv06. **F, S, V3:** as inv05 and inv06.
**RETRACT.**

### inv30 — memory: two modes
**Where:** `project_two_modes_scope.md:198`.
**Old claim:** the bridge cut "loses TBC1D3 AE/CDKL". **F, S, V3:** as inv16.
**SUSPEND.**

### inv31 — memory: leader nesting
**Where:** `project_leader_rule_breaks_nesting.md:26-27`.
**Old claim:** ADVISOR misstates CDKL; "tree 3/10".

| truth | result |
|---|---|
| F | The conclusion holds and is stronger. Projection CDKL 4/10 plus reference (`glo:253,256`). |
| S | 4/10; reference by any-match only (`glo:269,272,277`) |
| V3 | reference unsupported (`r2.out:23`) |

**RELABEL.**

### inv32 — memory: layer order
**Where:** `project_npip_tbc1d3_layer_order.md:19-20`.
**Old claim:** "TBC1D3 has no literature clades".
**F:** literature C = M + CDKL (`lo:8-10`). **S:** as F. **V3:** M (`r2.out:61,74`).
**RETRACT** that clause. The clause-5 groups are unaffected.

### inv33 — positional cluster1|cluster2 calls
**Where:** ledger §6jg-§6jx; `lattice_filtration.py:211`; `LAYER_ORDER:75`.
**Old claim:** the boundary is recovered exactly 0/2; the reference intron tree supports it ("wrong"); leave-out runs are mostly "correct".
**F:** identical to OLD in every rescore (`glo:229,232,255,258,301,304,327,330`; `gt:63,79`; `d2:8,14`). **S:** identical (`glo:245,248,271,274`; `gt:104,120`). **V3:** the positional scorer uses level 1 only.
**UNAFFECTED.** The interpretation is weaker: every Guitart group nests inside one cluster (`RA/inv04/run.stdout:45,47`).

**Tally:** retract (any part) = inv03, 04, 05, 06, 10, 11, 12, 13, 17, 18, 20, 21, 22, 28, 29, 32 (16). Suspend only = inv07, 08, 14, 16, 30 (5). Relabel = inv01, 02, 09, 15, 19, 23, 24, 25, 26, 27, 31 (11). Unaffected = inv33 (1). Heavy = none.

---

## §3 Audit findings and how each was handled

1. **Important: retractions that rest only on D=AE.**
   - inv07, 08, 14, 16 and 30 are reclassified from retract to **suspend**. So are the matching parts of inv03, 04, 17, 21 and 28.
   - V3 numbers were added from `AU/r1…r8`.
   - The AA totals were recomputed under V3 (`WC/aa_v3.out:5-6`: 94 vs 84, NOT SUPPORTED, gate PASS).
   - The audit's list of safe retractions is kept as it was.
   - **Residual:** D itself is unresolved (§4).
2. **Important: S over-trusted K.**
   - S is renamed "D excluded", with an explicit note that it is not sequence-safe.
   - Robustness tiers were added to §0, and every CDKL positive is worded as figure-dependent.
   - AA was re-checked with CDKL dropped (M only: 94 vs 84, `WC/aa_v3.out:11`).
   - The text anchor for T in CDKL is now cited.
   - **Residual:** K's membership is open (§4).
3. **Important: one dataset counted as several confirmations.** Fixed in §1.5 using the md5 and site-count evidence verified here. The CDKL UFBoot range is stated, and §6 wording says "one intronic dataset".
4. **Minor: "exon shows neither" is alignment-dependent.** The inv12 row and the §6 wording no longer say "intron-only".
5. **Minor: inv13 was labelled relabel.** Now **retract** (cell text); the bar verdicts are unaffected.
6. **Minor: other consumers of the old map.** Added to §6 (items 8, 9, 32). **Residual:** nothing is applied yet (§4).
7. **Minor: inv19 wording still owed.** §6 item 35.
8. **Minor: independent recomputation coverage.** The verifier-2 site-count swap is corrected in §1.3. The analyses that were not independently recomputed are listed in §4.

---

## §4 Open issues

- **D's group is unresolved** (remainder of audit finding 1). It needs Guitart's per-copy supplementary tables S4-S10 (not on disk) or the authors' tree. Until then the 10 suspended parts cannot be settled either way.
- **V3 not computed for the §6jp leave-out counts** (guided_lo MAFFT and projection). Only the reference trees (`AU/r2`) and the §6js leave-outs have V3.
- **Audit finding 2, verbatim:** "The mapping is over-trusted in one specific way: S is called 'sequence-safe' but keeps K in CDKL. K's 6B arrow is hand-placed and sequence trees do not support {K,TBC1D3} as a clade, while the D=AE arrow is the data-driven one.
  - Robust (figure and sequence agree): M={B,H}, and TBC1D3 in CDKL.
  - TBC1D3 in CDKL also has a textual anchor no batch item cites: Guitart p.15-16 says the expressed CDKL paralog is the last, telomeric-facing copy of cluster 2 in 67/69 haplotypes.
  - K in CDKL is moderate; D is unresolved.
  All CDKL-based positives (CHM13-alone CDKL, bridge-cut CDKL, C_fine CDKL) should be worded as figure-dependent. Only M is a robust new positive."
- **Audit finding 6, verbatim:** "Several places still carry the old mapping and no item covers them.
  - Two audit scripts hard-code AE/CDKL group names and would silently drop M on a rerun.
  - Lattice per-gene tables still print the old labels.
  - The §6kq section heading still claims CDKL is recovered.
  - No truth file has been corrected on disk yet; this is expected with no edits allowed, but every action is still pending."
- **Audit finding 8, verbatim:** "I rebuilt 7 analyses from the original inputs with my own code (Newick parser, UPGMA, ARI, exhaustive Jaccard matching, permutation test). None of the F or S numbers disagrees with the batch or done results. Not recomputed:
  - guided_lo MAFFT/projection leave-out counts (inv10/inv12/inv25), which need guided_tree_width's candidate relabelling
  - the DT split-dominance totals (inv15)
  - the §6jh coverage/de novo arms (inv26/inv29)
  - the pooled NPIP+TBC1D3 lattice columns
  A separate record error: verifier 2 swapped the guided_t reference-tree site counts (it lists exon 8,645 / intron 1,833; the files say exon 1,833 / intron 8,645). No call is affected."

---

## §5 What to tell the advisor now

**Say:**
1. The per-copy TBC1D3 subfamily truth I used was a name-based guess, and it was wrong for 7 of the 9 CHM13 copies. I re-derived it from Guitart 2024 Fig 6B, whose CHM13 track is a coordinate plot, and checked it against sequence.
2. One group is robust: **M = {TBC1D3B, TBC1D3H}**, which Guitart reports as absent from GRCh38. Our CHM13 trees recover it (intron 100/100, exon 77.7/71; 9-10 of 10 leave-outs). So do the HG002 and gorilla panel trees, and it is exactly the first split of the guided identity tree.
3. Guitart's groups sit inside the positional clusters (their text: "specific to either cluster 1 or 2"). I withdraw "within human the copies are one clade"; the positional calls themselves do not change.
4. The split-dominance-tree verdict (NOT SUPPORTED) holds under every version of the truth.
5. I withdraw three claims:
   - CDKL = {D,K} "recovered exactly". The "0.9944 median divergence" quoted with it was an identity.
   - "Extra haplotypes recover CDKL".
   - "AE is never recovered".

**Do not say:**
- That CDKL {K, TBC1D3} is recovered. That rests on the figure alone: K's arrow is hand-placed, UFBoot is below 80 on one intronic dataset, and the call fails if D is CDKL.
- Anything about {D,K}, either as a "sister pair" or as a "cross-group merge", or about CHM13-alone vs panels for CDKL. All of these depend on D.
- "AE never recovered". AE has at most one CHM13 copy.
- "The groups come only from introns". That depends on the alignment.
- "Literature fact" for the per-copy labels. They were read post hoc off a bioRxiv figure, and they are population-level groups (≥ 10 paralogs, 1.5× allelic cut), not an identity partition of CHM13.

**Ask:** can Guitart's per-copy group tables (supplementary S4-S10) be obtained? That alone resolves D.

---

## §6 Edits still needed

Apply only after user approval. Line numbers are as of HEAD b378435a.

Rules:
- Pre-registrations and the negative-results register are append-only: add an addendum or row and mark the old text.
- The ledger gets inline markers plus a new §6l7 that carries §2 of this report.
- Code edits change only how groups are derived.

### A. Truth files

1. **`docs/lit_subclusters_npip_tbc1d3_truth.tsv:24-32`**, and its byte-identical copies `LIT/lit_truth.tsv`, `LIT/guided_lo/truth.tsv` and `LIT/guided_t/truth.tsv`.
   - **Recommended:** keep these files, because they reproduce the registered numbers, and point new runs at `R00/truth_guitart_fig.tsv` (F) or an S copy with TBC1D3D's level2 blank.
   - **If overwritten,** level2 becomes: TBC1D3→CDKL, B→M, D→AE, E→O, F→G, G→B, H→M, I→I, K→CDKL.
2. **`/mnt/linuxdisk/home/juanfraitu/layer_order/npip_tbc1d3/light/truth_literature_subfamilies{,.corrected}.tsv` and `lattice/nodes.tsv` (`lit_level2`).** Regenerate from the new truth; the scripts are truth-agnostic (inv27).

### B. Code

3. **`bench/guided_pipeline.py:509`** → `vals = sorted({t["level2"] for t in T if t["level2"]}); g = {v: {t["name"] for t in T if t["level2"] == v} for v in vals}; g = {k: s for k, s in g.items() if len(s) >= 2}`
4. **`bench/guided_tree_width.py:33-34`** → the same construct for `GROUPS["TBC1D3"]`.
5. **`LIT/lit_tbc_diag.py:33`** → `for a, b in (("TBC1D3B", "TBC1D3H"), ("TBC1D3K", "TBC1D3")):`
6. **`…/light/scripts/truths_universe.py:211-212`** → "Guitart 2024 Fig 6B/6C per-copy group read from the CHM13 track (M = B+H, CDKL = K+TBC1D3; I, B (=RefSeq G), G (=RefSeq F), O (=E), AE (=D) singletons; post hoc; D figure-only, sequence leans CDKL)"
7. **`bench/layer_order/lattice_truth.py:13`** → "TBC1D3 Guitart Fig 6B/6C groups M, CDKL (phylogenetic, read post hoc; not positional)"
8. **`…/npip_tbc1d3/verify_slim_dna/verify_c.py:68-69` and `…/audit_slim_recompute/a4_ctree.py:109`** → replace the literal "AE"/"CDKL" names with level2-derived multi-member groups.
9. **`…/npip_tbc1d3/lattice/report_tables.md:829-837`** → regenerate; it still prints "lit L2 AE/CDKL".

### C. Pre-registrations and register (append only)

10. **`docs/PREREG_core_definition_2026-09-12.md`** → add a new ADDENDUM (next free letter after AP): "post-hoc TBC1D3 level-2 truth amendment".
    - It gives F, S and V3, and states that no registered reading changes: AA is NOT SUPPORTED in all variants, B2 excludes TBC1D3, and O-4 and P-3 are positional.
    - Tag lines 223-225, 579, 671 and 772 with "[superseded by the TBC1D3 truth addendum]".
11. **`docs/PREREG_known_subclusters_2026-09-09.md`**
    - `:22` → append "[2026-09-16: withdrawn — Guitart's phylogenetic groups nest inside clusters 1/2 ('specific to either cluster 1 or 2'); see bench/TBC1D3_GUITART_TRUTH_CORRECTION.md]"
    - `:36` → append "[suspended: TBC1D3D's group is unresolved (figure AE, sequence CDKL)]"
12. **`docs/NEGATIVE_RESULTS_REGISTER.md:1387`** (row 778) → replace "and there is no sequence signal to find" and the "⭐ Our only cut, {D,K} … sister pair" sentence with: "The clusters show no bimodality (between-cluster median 0.9944 inside the within-cluster ranges). Guitart's phylogenetic groups nest inside them (M={B,H} in cluster 1, CDKL={K,TBC1D3} in cluster 2), so 'no signal to find' is withdrawn. Whether the {D,K} cut is a sister pair or a cross-group merge depends on TBC1D3D's unresolved group."
13. **`docs/NEGATIVE_RESULTS_REGISTER.md:1433`** (row 822) → append "(TBC1D3 truth corrected 2026-09-16: 100 vs 84 figure truth, 97 vs 84 D excluded, 94 vs 84 D in CDKL — still not supported)".

### D. `docs/o1_ledger.md` (inline markers; add §6l7 = §2 of this report)

14. **`:14933-14934`** → after "a recent tandem sister pair." insert "[§6l7: 'no sequence signal' withdrawn — Guitart groups nest inside the clusters; '{D,K} = strongest real signal' suspended pending TBC1D3D's group]".
15. **`:18891`** → "Guitart/Eichler 2024 phylogenetic groups read per copy from Fig 6B/6C (M = B+H, CDKL = K+TBC1D3, rest singletons; D's AE label figure-only) [§6l7: the original name map was wrong for 7/9 copies]"
16. **`:18909`** → "at 0.0023 the cut is {D, K}, which recovers no Fig 6C group (ARI −0.038; 0.000 with D excluded; 0.372 with D in CDKL, still 0/2 groups exact); the 2-way root split {B, H} vs rest is exactly group M."
17. **`:18911`** → "cut gives {D, E, K} and {G, H}; ARI −0.061 vs Fig 6C (0.364 under the superseded name map)."
18. **`:18917`** → "…the tightest published NPIP groups (A6-9, B3-5 pairs) appear; for TBC1D3 only the root split (= group M) matches a published group; the major boundaries…"
19. **`:18971-18973`** → rows under Fig 6C:
    - guided coverage × identity 4/7 (sens 0, micro 0.778/0.778)
    - guided identity 3/7 (micro 0.667/0.750)
    - de novo identity 4/7 (micro 0.750/0.857)
20. **`:18975-18977`** → "Under Fig 6C every TBC1D3 pairwise sensitivity is 0 (no cut recovers M or CDKL); de novo coverage (5/7) beats de novo identity (4/7)."
21. **`:18985`** "TBC1D3 fully by de novo identity" → "TBC1D3 by no arm (best 5/7, all exact matches singletons)". **`:18986`** "name-mapped TBC1D3 groups" → "figure-derived TBC1D3 groups (§6l7)".
22. **`:18993-18995`** → "level 2 is a single decision: is B-H the closest pair — it never is (rank 16/28 spliced, 6/28 gene span and de novo); all three unit types give 4/7 on the shared 8 copies, so the 'not a mode effect' reading holds."
23. **`:19042`** → append "[§6l7: under Fig 6C the 6 shared records are all singletons (trivial); 5/5 holds only if D is CDKL — suspended]".
24. **`:19399`** → append "; the {H, K} split crosses Fig 6C groups M and CDKL".
25. **`:19438-19439`** → replace with:
    - "| TBC1D3 | M {B,H} | 10/10 | recovered (100/100) | 9/10 | recovered (100/100) |"
    - "| TBC1D3 | CDKL {K,TBC1D3} (figure-dependent) | 2/10 | recovered (84.2/76) | 4/10 | recovered (83.5/77) |"

    Add the note "same intronic input for both methods; CDKL UFBoot < 80; CDKL fails if D is CDKL".
26. **`:19522`** → "M {B,H} stable (ref 100/100, leave-out 9-10/10); CDKL {K,TBC1D3} reference only, figure-dependent"
27. **`:19580`** → "| TBC1D3 | M / CDKL | no / no | recovered 100/100 / recovered 84/77 | 0 / 0 | 9 / 4 | 9 / 4 |". **`:19593`** → keep the positional sentence, then append "(in this exon alignment neither Fig 6C group is recovered; the §6js exon alignment recovers M at 77.7/71 — do not read the groups as intron-only)".
28. **`:19626`** B2 cell → "Fig 6C M {B,H} recovered (ref exon 77.7/71, intron 100/100; leave-out 9/10); CDKL {K,TBC1D3} intron ref 83.5/77, leave-out 5/10, figure-dependent — B2 still excludes TBC1D3"
29. **`:19843-19844`** → tag the CDKL cell "[§6l7 suspended: absent under Fig 6C; recovered only if D is CDKL]".
30. **`:19899`** → "| TBC1D3 | M / CDKL | 9/10 + ref / 5/10 + ref | 10/10 + ref / 0/10 |". **`:19903-19904`** → append "(Fig 6C truth 100 vs 84, conflicts 67 vs 47; D excluded 97/84; D in CDKL 94/84 — NOT SUPPORTED throughout)".
31. **`:20633` and `:20643`** → "where the bridge cut adds the ninth locus, Fig 6C CDKL becomes evaluable and is recovered in the intron tree (DN1r+B 85.2/62, DN0+B 85.0/63), M unsupported; if D is CDKL this reverses — suspended (§6l7)". Tag the CDKL cells at `:20630-20632` the same way.
32. **`:20862`** heading → "…extra haplotypes make allelic divergence measurable; TBC1D3 subgroup claims corrected in §6l7"
33. **§6kq body lines**
    - `:20880` → "M recovered; CDKL {K,TBC1D3} no split"
    - `:20882` → "M recovered (100/100); CDKL 76.2/81 (figure-dependent)"
    - `:20885-20886` → "CHM13 {D,K} join 12 HG002 copies, but TBC1D3 stays apart, so Fig 6C CDKL does not form in any truth variant; B and H also fall in separate components"
    - `:20891-20892` → delete the CDKL-recovered and AE sentences
    - `:20913` → "M unsupported (registered; any-match 99.8/100); CDKL no split"
    - `:20919-20921` → delete "the tree clade recovers CDKL" and "AE stays unrecovered in every analysis"
34. **`:21452`** → "(TBC1D3 group M recovered by CHM13 trees and the HG002/gorilla panels; CDKL figure-dependent; AE untestable, §6l7)"
35. **`:21496`** → "where Dishuck/Bolognini-Yilmaz groups are literature facts and Guitart's per-copy TBC1D3 labels were read post hoc from a figure (one copy unresolved)"

### E. `docs/ADVISOR_QUESTIONS.md`

36. **`:776`** → "Guitart Fig 6C group M {TBC1D3B, TBC1D3H} is a supported clade on the CHM13 reference (intron 100/100, exon 77.7/71; leave-out 9/10). [Withdrawn: 'CDKL recovered exactly (0.9944 median divergence)' — that cut {D,K} is not a Fig 6C group, and 0.9944 is an identity.]"
37. **`:810-812`** → "…makes allelic divergence directly measurable. With the Fig 6C truth the HG002 and gorilla panels recover group M; CDKL {K,TBC1D3} is not recovered with HG002 or 23 haplotypes. [Withdrawn: 'recovers CDKL, not from CHM13 alone' (wrong name map) and 'AE never recovered' (AE has one CHM13 copy).]"
38. **`:933`** "TBC1D3-AE never recovered" → "TBC1D3 group M recovered; CDKL figure-dependent".

### F. Reports in `bench/`

39. **`LAYER_ORDER_NPIP_TBC1D3.md:18`** → "Literature C has 2 TBC1D3 pairs from Guitart Fig 6C: M={B,H} (supported, compatible) and CDKL={K,TBC1D3} (intron tree 83.5/77, conflicts with {E,K}); both lie inside P's and D's TBC1D3 group. With D in CDKL only M remains."
40. **`:19` and `:155`** → add "fails for TBC1D3: C_L1 ⊇ C_fine 0/2 (positional L1 cleared); pooled 0.889 [18] 4/6".
41. **`:75`** last cell → "Literature C: 2 TBC1D3 pairs (M, CDKL). C_tree: 12 top pairs".
42. **`:131-135`** → regenerate from `RB/inv22_layer_order/containment_old_vs_new.tsv`: P, D, C_tree_top and C_tree_min each contain C_fine/C_mid at 1.000 [2] 2/2.
43. **Truth-agreement rows**
    - `:251` → TBC1D3 "1.000 / 1.000 / 1.000 (9, circular)"; pooled "1.000 / 1.000 / 1.000 (31)"
    - `:253` → "0.167 / 1.000 / 0.556 (9)"; pooled "0.310 / 1.000 / 0.633 (30)"
    - `:254` → "0.286 / 1.000 / 0.778 (9)"; pooled "0.545 / 0.333 / 0.815 (30)"
44. **`:263`** → "TBC1D3 M is supported (exon 77.7/71, intron 100/100) and compatible; CDKL is supported only in the intron tree (83.5/77, UFBoot 77) and conflicts with {E,K} (exon 78.0/70)."
45. **`:388`** → "TBC1D3 L2 groups read per copy from Guitart Fig 6B/6C (post hoc); TBC1D3D unresolved (figure AE, sequence CDKL)."
46. **`NESTED_LATTICE_NPIP_TBC1D3.md`**
    - `:442` last cell → "C_tree_min 0.778; C_tree_top 0.556"
    - `:684-685` → "its literature L2 has two non-trivial groups, M (B, H) and CDKL (K, TBC1D3), read post hoc from Guitart Fig 6B/6C; TBC1D3D's group is unresolved."

### G. Memory (`~/.claude/projects/-mnt-c-Users-jfris-Desktop/memory/`)

47. **`project_tbc1d3_subclusters_population_unit.md`**
    - `:14` → "AE merges 2 GRCh38 genes (A, E); CDKL merges 4 (C, D, K, L)"
    - `:21-22` → "The name map is REFUTED for 7/9 copies by Guitart Fig 6B/6C: M={B,H} (robust), CDKL={K,TBC1D3} (figure-dependent), D=AE figure-only (sequence leans CDKL). See bench/TBC1D3_GUITART_TRUTH_CORRECTION.md."
    - `:23-24` → "§6kq under Fig 6C: HG002 panel recovers M, not CDKL; gorilla panel both; 23 haplotypes CDKL no split, M any-match only. 'CDKL not from CHM13 alone' suspended (depends on D). AE untestable."
48. **`MEMORY.md:83`** → "…name-mapped AE/CDKL truth REFUTED (7/9 copies) by Guitart Fig 6B/6C: M={B,H} robust, CDKL={K,TBC1D3} figure-dependent, D unresolved."
49. **`project_denovo_vs_annotated_gap.md`**
    - `:391` → "Guitart groups read from Fig 6C (M={B,H}, CDKL={K,TBC1D3})"
    - `:392` "TBC1D3(AE copy) no rep" → "unsuffixed TBC1D3 (a CDKL copy) no rep"
    - `:393` → "guided cut {D,K} recovers no Fig 6C group (ARI −0.038); its root split {B,H} = group M"
    - `:399-400` → "de novo identity TBC1D3 4/7 under Fig 6C (no arm recovers M or CDKL; best de novo coverage 5/7); guided cov x id 4/7"
    - `:403-404` → "L2 hinges on B-H, never the closest pair; all unit types 4/7 on the shared 8 copies; 'not a mode effect' holds"
50. **`project_two_modes_scope.md:198`** → "…but loses NPIP B3-5; TBC1D3 subgroup cells suspended (under Fig 6C the bridge cut makes CDKL evaluable and recovers it; if D is CDKL the original 'not recovered' holds)"
51. **`project_leader_rule_breaks_nesting.md:27`** → "misstates TBC1D3-CDKL: the 'exact' {D,K} cut is not a Fig 6C group (ARI −0.038); 0.9944 is identity; Fig 6C CDKL projection tree 4/10 + reference (figure-dependent)"
52. **`project_npip_tbc1d3_layer_order.md:19-20`** → "TBC1D3 literature clades (Fig 6C, post hoc): M={B,H} supported and compatible; CDKL={K,TBC1D3} supported (intron) but conflicts with {E,K}; both nest in the clause-5 groups {B,F,G,H} and {TBC1D3,D,E,K}; C_L1 ⊇ C_fine fails for TBC1D3."
