# NPIP membership re-run of the shared-junction conjunct t_J (k = 2, δ = 0)

Agent 1 of 2. Workspace `/mnt/linuxdisk/home/juanfraitu/npip_membership/`. No new alignment, no new catalog run, no `src/` change, no commit, no subagent, nothing written under `/mnt/c/Users/jfris/Desktop/Rustle`.

**STEP 0.** `/mnt/linuxdisk/home/juanfraitu/npip_membership/DECLARATIONS.txt` was written at **2026-09-17T17:35:41-07:00**, md5 **`1089aba18f9259cc85e140399bb46eda`**, **before any number in this run existed**. It fixes the three truth sets, the four fixes, the decision rule verbatim, the operationalization of its four clauses (S1.1, including the anti-trap that an h_join fall to 0.0 on an emptied subgraph is recorded as VACUOUS, not as clause 1 satisfied), the arms/levels/weights/certificate, and the node-quality definition. Nothing in it was changed afterwards.

Outputs and manifest: `out/*.certificates2.tsv`, `out/verdicts.tsv`, `out/boundary.tsv`, `out/*.nodequality.tsv`, `out/*.witness2.tsv`, `out/jmatch_compare.json`, `out/PTR.T_native.tsv`, `out/member_induced_components.tsv`, `MANIFEST.txt`.

---

> **Verifier corrections, applied by the orchestrator** (verifier returned ok = false; the de novo "confirmation" is WITHDRAWN):
> 1. **T_native is 47 records (42 in scope), not 46 (41)** — LOC112205920 matches the declared rule through a child feature. The confirming row restated: members 42, **5 of 758 member pairs lost** (not 3 of 719), boundary outside nodes **30 → 9**.
> 2. **WITHDRAW "CONFIRMED on the primary de novo arm (chimp T_strict3)".** That row's h_join fall is produced by cutting edges to records chimp RefSeq itself annotates as NPIP family — i.e. by cutting true family members that the 3-copy truth set does not contain.
> 3. **Clause C4(a) is vacuous on T_native by construction** (the NPIP-named records *are* the truth set, so "no NPIP-family survivor" is identically true). The one genuine NPIP-family survivor of the confirming row is LOC749627, "putative NPIP-like protein", at weight 0.993277.
> 4. **FIX 1 is disputed**: the verifier's independent map-through implementation reproduces the *previous* run's junction counts, not the "fixed" ones. The 358 moved rows are gene-body chain rows where donor and acceptor can come from different records of one chain; which reading is correct is unresolved. No h_join or member-pair count changes either way.
> 5. **`n_boundary_edges` is a node count, not an edge count.** True boundary-edge counts for the confirming row: 444 → 184 (46-record set) / 418 → 152 (47-record set).
> 6. **"6 T_member loci swallow 3-5 T_native records" is a miscount**: 5 loci (3, 3, 3, 4, 5 records), widths 76.0-290.5 kb.
> 7. **"29 of 49 surviving boundary nodes are NPIP-family" becomes 30 of 49** on the corrected set; the artefact is unchanged in direction and force.
> 8. **Structural point for the next pre-registration**: on a truth set that holds a small minority of the real family, most genuine copies sit at the *boundary*, so any rule that cuts boundary edges will look good. Size the truth set to the family before testing a boundary-cutting rule.

## 1. The three truth sets, verified from the label tables before use

| set | gorilla | chimp | in arm scope |
|---|---|---|---|
| **T_strict** (both projections name the SAME human copy) | **2** (NPIPA2 `NC_073242.2:32,426,774-32,456,482`; NPIPB2 `35,235,861-35,267,752`) | **4 as written**, 3 after removing the triple-claimed locus | gorilla 2, chimp 4 |
| **T_member** (both projections land, no name) — **new primary** | **25** of 26 candidates (23 + 1 + 1 across the contig trio) | **19** | gorilla 25, chimp **16** (3 out of scope) |
| **T_native** (chimp only; PTR RefSeq product/description names the NPIP family) | n/a | **46** records | **41** on NC_072416.2 |

The task text's "2 gorilla, 3 chimp" for T_strict is right for gorilla; the chimp label table carries **4** `LABELLED` rows, the fourth being `NC_072416.2:30,970,428-31,260,911`, claimed by three liftoff records (NPIPB10P, NPIPB7, NPIPB9). Both readings are scored side by side as **T_strict4** and **T_strict3**; every chimp claim below says which.

T_member counts match the task exactly (25 / 19). T_native is **46**, not the 45 found earlier — the declared string matches both the hyphenated and unhyphenated spellings.

**T_native ∩ T_member (the overlap the task asked for).** 18 of the 19 T_member loci contain at least one T_native record; the exception is `NC_072416.2:36,607,600-36,626,323`, whose native record LOC112207901 is *polycystin-1-like*. **32 of the 46** T_native records fall inside some T_member locus; **14 fall outside every one of them.** Six T_member loci each swallow 3–5 T_native records:

- `17,521,703-17,700,578` → LOC129135381, LOC112205827, LOC750787, LOC112204066, LOC112205833
- `30,970,428-31,260,911` → LOC749096, LOC129135420, LOC129137416, LOC129137415
- `31,799,523-31,977,765` → LOC129137423, LOC112205839, LOC129137424
- `33,030,791-33,106,820` → LOC129137437, LOC129137441, LOC129137436
- `33,308,180-33,386,432` → LOC129137435, LOC750420, LOC129137438
- `35,599,590-35,649,707` → LOC112207360, LOC129137455

That is the whole chimp artefact in one table: the projection truth merges tandem arrays into single intervals up to 179 kb wide and leaves 14 genuine members outside the set.

**Gorilla's own annotation does name NPIP records** — 5 of them (`out/GGO.native_npip_named.tsv`): LOC115934567, LOC109028586 (off the contig trio), LOC101141990, NPIPB11, LOC101134557. All three in-scope ones are natives of T_member loci, so gorilla's T_member boundary contains **0** NPIP-named records — the previous run's claim holds, and now with the reason.

---

## 2. What each of the four fixes actually moved

**Fix 1 — jmatch recomputed from the captured PAFs, independently (S3.1).** All 42,622 witness rows re-derived; the row set is identical (0 key-column mismatches).

| table | tx rows | tx differing | body rows | body differing | k=2 flips (body) |
|---|---|---|---|---|---|
| GGO.dn | 2,500 | **0** | 441 | 4 | 0 |
| GGO.ga | 8,595 | **0** | 6,127 | 29 | 15 (10 reject, 5 admit) |
| GGO.gr | 8,595 | **0** | 6,127 | 0 | 0 |
| PTR.dn | 890 | **0** | 1,497 | 22 | 10 (7 admit, 3 reject) |
| PTR.ga | 1,458 | **0** | 2,467 | 175 | 145 (72 reject, 73 admit) |
| PTR.gr | 1,458 | **0** | 2,467 | 128 | 92 (40 reject, 52 admit) |

**The previous verifier's defect (1) is not reproduced as stated.** On every one of the 23,496 transcript rows my independent map-through recompute is *identical* to the shipped `jmatch_d0`. The "16 of 2,500 / 35 of 8,595 transcript rows" it reported is the gap between the literal definition (map the donor and acceptor bases through the alignment) and an N-strict reading (require the record's own CIGAR to have an intron there): that contrast differs on **15** GGO.ga and **6** PTR.ga tx rows (2 admission flips), and **0** rows elsewhere. It is a definitional choice, not an arithmetic error.

The real movement is entirely in **gene-body chain** rows: a chain holds several PAF records, and the shipped code can take a junction's donor from one member record and its acceptor from another. 358 rows and 262 k=2 admissions move. A post-hoc sensitivity that requires one single member record to carry both ends moves a further 5 / 48 / 5 / 22 / 154 / 110 body rows (158 flips); it is reported, never used for a verdict.

**Net effect on every previously reported quantity: none.** Across the 144 comparable certificate rows (2 species × 2 sets × 3 arms × 3 variants × 4 axes) **h_join and the member-pair count at the cut are byte-identical** to the previous run, and the old `comp_size_at_cut` is reproduced by the new `comp_outside_total` (gorilla guided T_member 48 → 3 replicates exactly).

**Fix 2 — body chains excluded from L2/L3.** Exactly **3** rows change, all chimp, all L2: `PTR ga BASELINE T_strict`, `PTR gr BASELINE T_strict`, `PTR gr GUARDED T_member`, h_split `0.000000 → -inf`. Precisely the 3 rows the verifier named.

**Fix 3 — native_introns as the declared quantity.** Differs at **9 of 26** gorilla loci (PKD1 25→45, LOC115933071 6→28, LOC129527568 6→8, LOC115931102 8→9, LOC115932781 6→25, LOC115932779 6→70, LOC115933079 6→26, NPIPB11 6→25, LOC115933039 6→61) and at **9 of 32** chimp loci. No T_member locus in either species has fewer than 2 introns under either reading, so the k=2 precondition is untouched. The declared quantity is itself wrong for pseudogenes whose exons hang off the gene record (3 in-scope T_native records score 0 while carrying 6–7 exons); the guided node's junction count, which is what t_J reads, is unaffected.

**Fix 4 — components.** Per-component member sizes, per-component outside counts and the component count are reported separately, plus a stricter count over the member-induced subgraph. Both readings agree on all four confirming rows.

---

## 3. PRIMARY ARM — de novo

### Gorilla (still vacuous; the arm does not test the conjunct)

12 of the 25 T_member copies have a de novo node, 0 collapse. Under STRICT, every level goes to h_join **0.000000 with zero boundary edges and zero member edges** — the anti-trap fires: **VACUOUS**, not clause 1. 11 of 12 members are newly isolated; the set goes to 12 singletons at every level (member-induced count identical). Baseline for reference: L1 identity 0.997868 / h_split −inf, 21 member pairs; L2 1.000000, 17 pairs; L3 0.997868, 5 pairs. On T_strict (|S| = 2) both copies are isolated. This reproduces the previous run exactly and for the same reason (node quality, §4).

### Chimp (the arm that tests it)

| set | \|S\| | BASELINE h_join | STRICT h_join | newly isolated | components (graph / member-induced) | member pairs | C1–C4 |
|---|---|---|---|---|---|---|---|
| **T_strict3** | 3 | 1.000000 | **0.996527** | **0** | 1 → 1 / 1 → 1 | 3 → 3 | **all four pass** |
| T_strict4 | 4 | 1.000000 | 0.996527 | 1 | 1 → 2 | 6 → 3 | C2, C3 fail |
| T_member | 16 | 1.000000 | **1.000000** | 7 | 2 → 9 | 85 → 28 | C1, C2, C3, C4 fail |
| T_native | 32 | 1.000000 | 0.993468 | 13 | 1 → 14 | 398 → 154 | C2, C3 fail |

So on the primary arm the conjunct **confirms only on the 3-copy T_strict3 set**, and only at L1 identity (L1 coverage, L2, L3 all stay at 1.000000 there). On the larger chimp sets the de novo arm pays the cost the human run predicted: it shatters the set (9–14 components) and isolates 7–13 copies. The T_strict3 positive is real but small, and it still never certifies (h_split 0.938667 < h_join).

---

## 4. SECONDARY — guided arms

### Chimp guided-annotated (GA), T_native — the confirming row

| level / axis | BASELINE h_join / h_split | STRICT | GUARDED | members isolated | components | member pairs | outside at cut |
|---|---|---|---|---|---|---|---|
| **L1 identity (0.80)** | 1.000000 / 0.919079 | **0.998233** / 0.884983 | 0.998233 | **0 of 41** | 1 → **1** | 719 → **716** | 86 → **33** |
| L1 coverage (0.50) | 1.000000 / 0.973191 | 1.000000 | 1.000000 | 0 | 1 → 1 | 719 → 716 | 86 → 33 |
| L2 f_ex (0.30) | 1.000000 / −inf | 1.000000 | 1.000000 | 2 | 3 → 4 | 311 → 307 | 45 → 32 |
| L3 w_98 (0.98) | 0.999111 / −inf | 0.999111 | 0.999111 | 2 | 10 → 10 | 147 → 145 | 44 → 32 |

The three lost member pairs are LOC100608310–LOC129137444, LOC112207414–LOC749096, LOC112207668–LOC129137437.

**Boundary composition at L1 identity** — direct boundary nodes 28 → 9, 19 removed, **0 NPIP-named in either group**:

*Removed* (top by baseline weight): LOC101057162 (3 exons, 989 bp, uncharacterized, **held h_join at 1.000000**), LOC134808749 (3 exons, 1,141 bp, uncharacterized, also 1.000000), LOC112207762 (28 exons), LOC134808806 (**2-exon lncRNA**, 923 bp), LOC129135256 (EIF3C-like), LOC129137434, LOC453981, LOC454016, LOC112205917, LOC112207897, LOC134808741 (1-exon pseudogene). Median exonic bp 2,473; 3 of 19 have ≤ 2 exons.

*Survivors*: LOC129137282 (15 exons, uncharacterized) and **LOC739285 (23 exons, polycystin-1)** at 0.998233, LOC129137403 (SMG1-like), LOC749627 (**"putative NPIP-like protein"**), LOC112205920 (ERVK env), LOC107973150 (SMG1-like), LOC112207901 (polycystin-1-like), LOC744474. Median exonic bp 3,066; **0 with ≤ 2 exons, 0 embedded**.

This is the human pattern replicated on a species that never saw the rule: the conjunct removes small uncharacterized records, a 2-exon lncRNA and a single-exon pseudogene, and the boundary passes to **full-size multi-exon neighbours, headed by PKD1 (polycystin-1) records** — exactly the "PKD1P–NPIP readthroughs take the boundary" queue the human run reported at 0.999363. It still does not certify (0.998233 > h_split 0.884983).

### Chimp guided arms, T_member — the artefact, measured

h_join stays at **exactly 1.000000** on GA and GR at every level. **29 of the 49 surviving boundary nodes are records chimp RefSeq itself calls NPIP-family** (LOC112204066, LOC129137441, LOC112205827, LOC750787, LOC129135381, LOC112205839, LOC112205845, LOC112207360, …), against 29 of 111 at baseline. The conjunct cannot cut an edge to a genuine family member and should not. Section 1 quantifies why they are outside the set. Clause C4 also flags 2 non-NPIP survivors (LOC129135256, 24 exons; LOC112207897, 10 exons) as "embedded" — but they are embedded only in the sense of falling inside a 179 kb truth interval; they are not small embedded records.

### Gorilla guided-annotated (GA), T_member — all four clauses pass, but gorilla cannot confirm

| level | BASELINE | STRICT | GUARDED | pairs lost | outside at cut | isolated | components |
|---|---|---|---|---|---|---|---|
| L1 identity | 1.000000 / 0.953568 | **0.988399** | 0.988399 | **0** of 270 | 48 → **3** | 0 | 1 → 1 |
| L1 coverage | 1.000000 / 0.989399 | 1.000000 | 1.000000 | 0 | 48 → 3 | 0 | 1 → 1 |
| L2 f_ex | 0.981969 / 0.683746 | 0.981969 | 0.981969 | 10 of 98 | 25 → 3 | 0 | 1 → 1 |
| L3 w_98 | 1.000000 / 0.897445 | **0.989975** | 0.989975 | 2 of 23 | 20 → 2 | 0 | 9 → 10 (fails C3) |

Removed at L1 identity: LOC115932701 (2 exons, 258 bp, embedded in the NPIPA2 copy, held h_join at 1.000000), LOC101131206 (4 exons, 766 bp, embedded lncRNA), LOC115932780 / LOC129527571 (1-exon ribosomal pseudogenes), LOC101134912 (SMG1), LOC129527593 (SMG1-like); 7 of 16 removed have ≤ 2 exons, median exonic bp 3,096. Survivors: LOC129528916 (MRP1-like, 5 exons, 0.988399), LOC115932992 (41 exons), LOC101130854 (MRP1, 33 exons) — **0 with ≤ 2 exons, 0 embedded**. Per S1.1, gorilla is not eligible to CONFIRM; it is reported as a four-clause pass on a substrate that holds out the species and the conjunct.

### Guided read-supported (GR)

Gorilla: destructive, as before — 269 of 270 member pairs lost at L1, 22 of 25 copies newly isolated, 23 components; 14 of the 25 member records carry zero ≥3-read junctions in the 37 %-downsampled BAM. Chimp: T_native loses 291 of 719 pairs and isolates 9 of 41; T_member loses 60 of 98 and isolates 3. GUARDED on GR restores h_join to 1.000000 in gorilla at 39 pairs lost — the human "the guard re-admits precisely what the conjunct exists to remove" pattern, replicated.

---

## 5. Node quality (S4.7)

| species | arm | set | copies (in scope) | with a node | node ≥ 2 junctions | node 0 junctions | median exonic fraction |
|---|---|---|---|---|---|---|---|
| GGO | dn | T_strict | 2 | 2 | 1 | 1 | 0.186 |
| GGO | dn | **T_member** | 25 | **12** | **4** | **7** | **0.146** |
| GGO | ga | T_member | 25 | 25 | 25 | 0 | 0.471* |
| GGO | gr | T_member | 25 | 25 | 7 | **14** | 0.471* |
| PTR | dn | T_strict4 / T_strict3 | 4 / 3 | 4 / 3 | 3 / 3 | 1 / 0 | 0.739 / 0.500 |
| PTR | dn | **T_member** | 16 | **16** | **9** | 7 | **0.801** |
| PTR | dn | **T_native** | 41 | **32** | **22** | 8 | **0.706** |
| PTR | ga | T_member / T_native | 16 / 41 | 16 / 41 | 16 / 41 | 0 / 0 | 0.998* / 1.000* |
| PTR | gr | T_member / T_native | 16 / 41 | 16 / 41 | 13 / 33 | 2 / 7 | 0.998* / 1.000* |

\* on the guided arms the node **is** the native record, so the fraction is ≈ 1 by construction and carries no information; only the de novo rows are informative. The gorilla/chimp de novo gap (0.146 vs 0.706–0.801; 12/25 vs 32/41 with a node) is confounded by depth (GGO_ds.bam 37 %-downsampled, PTR_mm.bam not) and is not a like-for-like comparison.

---

## 6. Certificates

216 certificate rows (GGO 72 + PTR 144). **7 rows have a non-empty interval, all on the gorilla de novo `T_strict` set with |S| = 2; 4 have the shipped cut inside** — BASELINE and GUARDED at L2 `(0.000000, 0.798879]` and L3 `(0.000000, 0.988792]`. All four have h_join = 0.000000 because that 2-node set has **no boundary edge at all** at those levels, so the interval contains every cut trivially. They are not certificates of anything and are reported only so the count is complete. **No STRICT row certifies, in either species, on any arm, at any level, on any of the three truth sets.** The human finding "it never certifies NPIP" replicates on both species and under all three truth definitions.

---

## 7. Clause-by-clause verdict (72 rows in `out/verdicts.tsv`)

13 rows have a **non-vacuous** h_join fall. Four rows pass all four clauses:

| species | set | arm | level | h_join | C1 | C2 | C3 | C4 | eligible to confirm? |
|---|---|---|---|---|---|---|---|---|---|
| **PTR** | **T_native** | **ga** | **L1 identity** | 1.000000 → 0.998233 | ✓ | ✓ (0 of 41) | ✓ (1 → 1) | ✓ (0 ≤2-exon, 0 embedded) | **yes — CONFIRMED** |
| **PTR** | **T_strict3** | **dn** | **L1 identity** | 1.000000 → 0.996527 | ✓ | ✓ (0 of 3) | ✓ (1 → 1) | ✓ | **yes — CONFIRMED (\|S\| = 3)** |
| PTR | T_strict3 | gr | L2 | 1.000000 → 0.996327 | ✓ | ✓ | ✓ trivially | ✓ | degenerate — 0 member pairs at the cut in both arms; do not quote |
| GGO | T_member | ga | L1 identity | 1.000000 → 0.988399 | ✓ | ✓ (0 of 25) | ✓ (1 → 1) | ✓ | no — gorilla developed the node rules |

Everything else fails at least one clause. The de novo arm fails on gorilla by vacuity and on chimp T_member/T_native by isolation and shattering; L1 coverage, L2 and L3 essentially never move h_join on chimp.

---

## 8. What this changes relative to the previous run

1. The previous run's chimp conclusion — "h_join does not move on any chimp arm, at any level; chimp gives no usable read" — was **an artefact of the truth set, and it is now fixed**. With the species' own annotation as the truth set, chimp moves at L1 identity on all three arms (ga 1.000000 → 0.998233, gr → 0.997015, dn → 0.993468) and the guided arm satisfies all four clauses.
2. The previous run's reported numbers survive the four fixes intact: the only changes are 3 h_split values (0.000000 → −inf) and the component reporting. 358 witness rows and 262 k=2 admissions moved without moving a single h_join.
3. The "29 of 63 outside nodes are natively annotated NPIP records" statement becomes **29 of 49** under the corrected jmatch, and its cause is now measured: 14 T_native records outside every T_member locus, 6 T_member loci swallowing 3–5 records each.
4. The gorilla claim "the gorilla annotation names no NPIP records" needs restating: it names **5**; the 3 in scope are all natives of T_member loci, which is why the gorilla boundary is free of them.

## Verification (independent recompute) — agent 2 of 2

All verifier code and outputs are under `/mnt/linuxdisk/home/juanfraitu/npip_membership/verify/`
(`v_native2.py`, `v_jm2.py`, `v_chain.py`, `v_jm3.py`, `v_cert2.py`, `v_run2.py`; outputs
`PTR.T_native_verify.tsv`, `{GGO,PTR}.{dn,ga,gr}.jm3.tsv`, `cert_*.tsv`). No builder script under
`npip_membership/scripts/` was read or imported; the only inputs used were `DECLARATIONS.txt`, the
data files in `out/`, the inherited evidence in `ggo_npip/` (labels, nodes, captured PAFs, previous
witness and certificate tables) and the original `GGO_genomic.gff` / `PTR_genomic.gff`.

### 1. Declarations precede the results
`DECLARATIONS.txt` md5 `1089aba18f9259cc85e140399bb46eda`, matching MANIFEST.txt, mtime
2026-09-17 17:35:41; the first result file is `out/GGO.native_introns.tsv` at 17:36:19 and the last
`out/member_induced_components.tsv` at 17:43:49. Every script mtime post-dates the declarations.
The file contains no addendum, i.e. nothing was added after a number existed. **PASS.** The
per-record "sensitivity" the builder reports is not in S3.1 (only PRIMARY and the N-STRICT
CONTRAST are); the builder labels it post-hoc in its concerns, which is the correct disclosure.

### 2. Truth sets re-derived from the label tables and PTR_genomic.gff
- GGO `GGO.truth.tsv`: 26 rows, `status == LABELLED` → **2** (T_strict), `both_arms == True` → **25**
  (T_member), all on the declared contig trio. **Matches.**
- PTR `PTR.truth.tsv`: 32 rows, `LABELLED` → **4** (T_strict4; the fourth,
  NC_072416.2:30,970,428-31,260,911, is claimed by three liftoff records NPIPB10P/NPIPB7/NPIPB9, so
  T_strict3 is justified), `both_arms == True` → **19**, of which 16 on NC_072416.2 and 3 out of
  scope (2 × NC_072404.2, 1 × NC_072417.2). **Matches.**
- T_native re-derived independently from `PTR_genomic.gff` (41,815 gene/pseudogene records): I get
  **47** records, 42 on NC_072416.2 — one more than the builder (see corrections). Every other
  record, and every `introns` value, is identical. Gorilla NPIP-named records: **5** (3 in scope),
  identical to `GGO.native_npip_named.tsv`.
- T_native/T_member overlap (my 47-record set): **15** T_native records lie outside every T_member
  locus; **5** T_member loci swallow ≥ 3 T_native records (3, 3, 3, 4, 5; widths 76.0–290.5 kb).
  On the builder's own 46-record set the "14 outside" figure reproduces exactly; the "6 loci"
  figure does not (see corrections).

### 3. The four declared fixes
**FIX 1 (jmatch).** I re-derived jmatch for all 42,622 witness rows from the captured PAFs with my
own code, including exact identification of each body chain (subset of PAF records whose Σbl, Σnm,
min ts and max te reproduce the row; 1,475 direct / 19 solved / 3 unresolved for PTR dn, similar
elsewhere). Row keys are identical (0 mismatches; note 2 duplicate keys in `GGO.dn.witness2.tsv`,
both jmatch 0, which is why the file has 2,941 data rows for 2,939 distinct keys).
- **Transcript rows (23,496): 0 disagreements** with either the shipped `jmatch_d0` or the builder's
  recompute. The previous verifier's defect (1) as stated is **not** reproduced — confirming the
  builder.
- **Body-chain rows: my values equal the OLD `jmatch_d0` on 42,619 of 42,622 rows overall and
  disagree with the builder on 360 rows.** The builder's 358-row / 262-admission movement versus the
  old table reproduces exactly, per table (4 / 29 / 0 / 22 / 175 / 128 moved; admission flips
  GGO ga 10 reject + 5 admit, PTR dn 7 admit + 3 reject, PTR ga 72 + 73, PTR gr 40 + 52). The
  disagreement is entirely the undeclared conflict rule inside a merged column map (see
  corrections), not chain membership.
- **Materiality: none.** Recomputing every certificate row under my convention, under last-writer-
  wins, and under the previous run's values changes **no h_join, no h_split, no component count and
  no isolation anywhere**; the only movement is ±1 member pair on 8 chimp T_native rows (confirming
  row: 2 lost under my convention, 3 under the builder's, 4 under last-wins).

**FIX 2 (body chains excluded from L2/L3).** Comparing all 144 rows comparable to the previous run
(GGO T_strict↔T_strict, GGO T_expl↔T_member, PTR T_strict↔T_strict4, PTR T_expl↔T_member, axis
labels normalised): **exactly 3 rows change, all chimp, all L2 — PTR ga BASELINE T_strict4,
PTR gr BASELINE T_strict4, PTR gr GUARDED T_member — h_split 0.000000 → −inf.** No h_join and no
member-pair count changes on any of the 144. **Exactly as claimed.**

**FIX 3 (native_introns).** Recomputed from the original GFFs as max over transcript children of
(exons − 1): **identical to the builder at all 26 gorilla and all 32 chimp loci**; differs from the
old column at **9/26** and **9/32**. The three zero-intron pseudogenes (LOC112207411, LOC129137441,
LOC112207360) are confirmed: their exons hang off the pseudogene record (7, 7 and 6 exons), so the
declared quantity returns 0. Disclosed correctly.

**FIX 4 (per-component reporting).** Present and correct: `n_components_at_cut`,
`comp_member_sizes`, `comp_outside_sizes`, `comp_outside_total`, plus the stricter member-induced
count in `member_induced_components.tsv`. `comp_outside_total` reproduces the old
`comp_size_at_cut` (gorilla guided T_member 48 → 3 replicates exactly). One display defect:
`comp_member_sizes` appears truncated to 8 entries on rows with more parts (e.g. PTR dn T_native
STRICT prints 8 sizes for 14 parts).

### 4. Certificates and node quality, recomputed from scratch
My own D1 engine (weights rebuilt from the witness tables per S4.3, body rows excluded from L2/L3,
h_join = max boundary weight, h_split = maximum-spanning-tree bottleneck, +inf for |S| = 1, −inf if
internally disconnected) reproduces **all 216 rows** of `{GGO,PTR}.certificates2.tsv` on h_join,
h_split, n_components_at_cut, comp_member_sizes, comp_outside_total, n_member_pairs_at_cut,
n_present / n_inscope / n_absent / n_outofscope, `vacuous_fall` and `shipped_inside`. The only
systematic mismatch is `n_boundary_edges`, which is a node count (correction 5).
- **Shipped-cut-inside rows: 4**, all gorilla de novo T_strict with |S| = 2 and no boundary edge at
  all (L2 and L3, BASELINE and GUARDED). **No STRICT row certifies in either species at any level.**
  Confirmed.
- **Chimp T_member h_join = 1.000000 on all 36 rows** (3 arms × 3 variants × 4 axes). Confirmed.
- **Gorilla de novo:** 12 of 25 T_member copies have a node, 7 of them carry zero junctions; STRICT
  collapses the member subgraph to 12 singletons, 21 → 0 member pairs, 11 newly isolated, h_join
  0.997868 → 0.000000 flagged vacuous. Confirmed.
- **Gorilla guided T_member:** h_join 1.000000 → 0.988399, 0 isolated, 1 component, 270 → 270 member
  pairs, comp_outside_total 48 → 3. Confirmed (gorilla is not eligible to confirm).
- **Node-quality table: all 18 rows reproduce exactly**, including the median exonic fractions.
- **Member-induced components** for both confirming rows: 1 → 1, and 1 → 1 with the corrected
  47-record T_native.

### 5. The verdict sentence against the numbers
- **Chimp, guided-annotated, T_native, L1 identity — CONFIRMED stands**, with restated numbers
  (correction 2): h_join 1.000000 → 0.998233, h_split 0.919079 → 0.884983, 0 of 42 copies isolated,
  1 component under both readings, 5 of 758 member pairs lost, boundary outside nodes 29 → 9, none
  with ≤ 2 exons and none embedded in a truth interval (verified node by node). The removal side is
  clean: **none of the 20 removed boundary nodes is NPIP-named or NPIP-like** (the two 1.000000
  edges removed go to LOC134808749 and LOC101057162, both "uncharacterized", 3 exons each). The
  caveats that must travel with the claim: the set is annotation-derived, node set and truth set
  come from the same GFF, C4(a) is vacuous on T_native by construction, and one survivor
  (LOC749627, "putative NPIP-like protein") is NPIP-family annotated.
- **Chimp, de novo, T_strict3 — the positive does not survive scrutiny** (correction 3): the h_join
  fall is achieved by cutting an edge to a node overlapping LOC112205918, and the new h_join is set
  by an edge to a node overlapping LOC750787, both records chimp RefSeq names NPIP-family; 15 of 21
  removed and 17 of 20 surviving boundary nodes overlap NPIP-named records.
- **Chimp T_member does not confirm on any arm or level** — confirmed independently, with the
  artefact slightly larger than reported (30 of 49 surviving boundary nodes NPIP-named).
- **PTR gr T_strict3 L2** passes all four clauses but has 0 member pairs in both arms; the builder
  already flags it as not quotable. Agreed.
- The summary sentence "they move 358 witness rows and 262 k=2 admissions but change no h_join and
  no member-pair count anywhere" is right about h_join and about the 144 comparable rows, but is
  **not** exactly right about member pairs on the new T_native rows under a different (equally
  declared-compatible) merge convention: 8 rows move by 1.