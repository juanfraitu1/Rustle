# NPIP subfamily truth checked against Dishuck et al. 2025

2026-09-16. Read-only check. No truth file was edited and nothing was committed. Nothing below was pre-registered.

**Per-copy table:** `docs/lit_subclusters_npip_dishuck_check.tsv`. It has 27 rows (the 22 truth records plus 5 CHM13 NPIP copies the truth omits) and uses 1-based GFF coordinates. The truth file's `start0` equals `start` − 1.

**Paper:** Dishuck et al., bioRxiv 2025.02.04.636496v1. Page numbers below are PDF pages.

**Checked truth rows:** the NPIP rows of `docs/lit_subclusters_npip_tbc1d3_truth.tsv`. They are byte-identical to `LIT/lit_truth.tsv` and `LIT/guided_lo/truth.tsv`.

**Path shorthands**
- `LIT` = `/mnt/linuxdisk/home/juanfraitu/o1_falsemerge/lit`
- `CK` = `/mnt/linuxdisk/home/juanfraitu/npip_dishuck_check` (scripts, outputs, figure crops)

---

## §0 Verdict

**Group membership is faithful; the label, the names and the copy list are not.** Unlike TBC1D3, where 7/9 labels changed, no NPIP record moves group. Matched by coordinates, every one of the 22 records sits in the same level-1 subfamily and the same level-2 group that the paper gives that copy.

There are three defects.

### 1. The level-2 groups are described wrongly

Addendum G calls A2/3, A6-9, B3-5, B6-9 and B12/13 "Dishuck paralog groups". They are not the paper's 28 phylogenetic paralog groups, which are defined at SH-aLRT > 75 (p.3; Fig 2B row labels, p.13).

- **What they actually are:** the paper's Iso-Seq read-assignment and expression groups. The text says "Grouping highly similar paralogs (A2/3, A6-9, B3-5, and B12/13)…" (p.7). Fig 6A (p.19) has 14 rows, and the project's 13 groups are exactly those rows minus A4.
- **Only one is a phylogenetic group.** B12/B13 is one of the 28 (p.3).
- **All five are clades in the published tree.** Each multi-copy group is monophyletic in Fig 2A (p.13), but the figure draws no support values.

### 2. RefSeq symbols name the wrong paralog for 8 of 22 records

In every case the level-2 group is unchanged.

| RefSeq record | Dishuck paralog |
|---|---|
| NPIPB3 | **B5** |
| NPIPB5 | **B35L1** (labelled "B3" in Fig 1A/4B) |
| NPIPB6 | **B7** |
| NPIPB7 | **B9L1** (labelled "B9" in Fig 1A) |
| NPIPB8 | **B7** |
| NPIPB9 | **B7** |
| NPIPA7, NPIPA8 | joint paralog **A7/A8** |

Any paralog-level use of these RefSeq names is wrong. In CHM13 there is no phylogenetic B3, B6, B8 or B9 copy.

### 3. The truth has 22 of CHM13's 27 NPIP copies

The paper counts 27 (p.3; Fig 1A, p.12). The truth kept only records whose symbol starts with NPIP. The 5 omitted copies are:

| omitted record | Dishuck paralog | group it belongs to |
|---|---|---|
| PKD1P6-NPIPP1 | A4 | its own group, in NPIPA |
| LOC128966608 | B5 | B3-5 (RefSeq calls it "B13-like") |
| LOC124907834 | B13 | B12/13 |
| LOC124907808 | B15 | B15 |
| LOC124907807 | B15 | B15 |

A complete truth has:
- **Level 1:** NPIPA 8, NPIPB 19.
- **Multi-copy groups:** B3-5 has 4 copies, B12/13 has 3, and B15 has 3 (B15 is a singleton in the current truth).
- **A new singleton:** A4.

### Level 1 (NPIPA|NPIPB)

This is the paper's subfamily split: nomenclature plus "IGC occurs within but not between" (p.4).
- NPIPA is a clade in Fig 2A (p.13), and reads as one in the timetree (Fig S1, p.26).
- NPIPB is **paraphyletic** in both rooted trees. In Fig 2A, the B6-B10/B14/B15 lineage lies outside the clade that joins NPIPA with B1/B2/B3-5/B11/B12-13. Fig S1 was read at lower confidence.
- So quote NPIPA|NPIPB as an unrooted split. Do not call it "two clades".

---

## §1 How it was checked

### Text and figures

- **Text:** pymupdf text layer.
- **Figures are raster images.** They carry no vector drawings, so the vector-geometry method used for TBC1D3 does not apply. Figures were rendered at 4-40× and read by eye; crops are in `CK/figs/`.

### Fig 3A (p.14) positions

- The centres of the 21 red bars were measured in pixels (the A6/A7 pair renders as one bar).
- Calibration used the 1 Mbp ticks: 12 Mbp at px 487, 207.4 px per Mbp.
- Every bar lies within −7 to +24 kb of one CHM13 RefSeq start (`CK/fig3a_positions.tsv`). The closest pair of copies is 47 kb apart.
- Labels were assigned through the two-row stagger.
- The label order was confirmed on:
  - the CHM13 tracks of Fig 3B-E, which have their own axes;
  - Fig 4B-D (p.16), which has CHM13 coordinate axes with a grey highlight at every copy, including the three B15 copies on 16q.
- Copy spacing in the Fig 1A blocks (100 kbp scale bar) matches RefSeq:
  - B5–B5: 341 kb vs 343 kb;
  - B13–B13: 104 kb vs 102 kb;
  - B15: 119 and 115 kb vs 114 and 113 kb;
  - B10–B11–B12: 348 and 100 kb vs 344 and 102 kb.

### Strands

Strands were not used. The DupMasker arrows are 3-5 px wide and cannot be read reliably. Positions alone are unambiguous here; TBC1D3 needed strands because its copies are closer together.

### Sequence check (no new tree)

- Every leaf of the existing NPIP trees was mapped back to CHM13 (minimap2 asm5 on the chr16/chr18 NPIP regions):
  - 44 treefiles in 4 sets (`LIT/guided_lo/{tree,tree_proj,tree_union}`, `LIT/guided_t/tree_t`);
  - each set holds 10 leave-out trees and 1 reference tree.
- Each leaf hits exactly one of the 27 copies at a substitution rate ≤ 0.0005. All 22 named leaves hit their own record (`CK/located*.tsv`).
- The trees were then rescored twice:
  - with the original groups — this reproduces the ledger counts exactly (§6jp P1 and Q, §6jr intron);
  - with complete, coordinate-keyed Dishuck groups.
- Scripts and outputs: `CK/rescore*.py` and `CK/rescore*.out`.
- Rescoring took under 1 minute. No alignment or tree was rebuilt.

---

## §2 Per-copy summary

Columns: RefSeq record (CHM13) | project level 1 / level 2 | Dishuck paralog / Iso-Seq group | agreement

| record | project | Dishuck | agreement |
|---|---|---|---|
| NPIPA1 | A / A1 | A1 / A1 | same |
| NPIPA2 | A / A2/3 | A2 / A2/3 | same (CHM13 has no A3) |
| **PKD1P6-NPIPP1** | — | **A4 / A4** | **not in truth** |
| NPIPA5 | A / A5 | A5 / A5 | same |
| NPIPA6 | A / A6-9 | A6 / A6-9 | same |
| NPIPA7 | A / A6-9 | **A7/A8** / A6-9 | renamed |
| NPIPA8 | A / A6-9 | **A7/A8** / A6-9 | renamed |
| NPIPA9 | A / A6-9 | A9 / A6-9 | same |
| NPIPB3 | B / B3-5 | **B5** / B3-5 | renamed |
| **LOC128966608** | — | **B5 / B3-5** | **not in truth** |
| NPIPB4 | B / B3-5 | B4 / B3-5 | same |
| NPIPB5 | B / B3-5 | **B35L1** ("B3" in Fig 1A/4B) / B3-5 | renamed |
| NPIPB6 | B / B6-9 | **B7** / B6-9 | renamed |
| NPIPB7 | B / B6-9 | **B9L1** ("B9" in Fig 1A) / B6-9 | renamed |
| NPIPB8 | B / B6-9 | **B7** / B6-9 | renamed |
| NPIPB9 | B / B6-9 | **B7** / B6-9 | renamed |
| NPIPB10P | B / B10 | B10 / B10 | same |
| NPIPB11 | B / B11 | B11 / B11 | same |
| NPIPB12 | B / B12/13 | B12/B13 / B12/13 | same |
| **LOC124907834** | — | **B12/B13 (label B13) / B12/13** | **not in truth** |
| NPIPB13 | B / B12/13 | B12/B13 / B12/13 | same |
| NPIPB14P | B / B14 | B14 / B14 | same (evidence: Fig 1A and text only) |
| NPIPB15 | B / B15 | B15 / B15 | same |
| **LOC124907808** | — | **B15 / B15** | **not in truth** |
| **LOC124907807** | — | **B15 / B15** | **not in truth** |
| NPIPB1P | B / B1 | B1 / B1 | same |

**Totals:** 14 same, 8 renamed (group unchanged), 0 in a different group, 5 not in the truth.

---

## §3 Impact on the NPIP claims to be quoted

### (a) §6jp — "intron tree recovers NPIPA|NPIPB, B3-5, B6-9, B12/13 in 10/10 leave-out runs"

**The claim stands; relabel it.** Every leave-out tree already contained the 4 LOC copies as unlabelled candidate leaves, and 3 of the 10 contained A4 (projection set). Under complete groups the counts do not move:

| intron leave-out set | L1 A\|B | B3-5 | B6-9 | B12/13 | B15 (3 copies) | A6-9 |
|---|---|---|---|---|---|---|
| §6jp Q projection (quoted) | 10/10 → **10/10** | 10 → **10** | 10 → **10** | 10 → **10** | — → **10** | 5 → 5 |
| §6jr intron | 10 → 10 | 10 → 10 | 10 → 10 | 10 → 10 | — → 10 | 6 → 6 |
| §6js intron | 10 → 10 | 10 → 10 | 10 → 10 | 10 → 10 | — → 10 | 5 → 5 |
| §6jp P1 MAFFT | 8 → 8 | 10 → 10 | 5 → 5 | 10 → 10 | — → 10 | 5 → 5 |

The named β-helix subfamily, rescored with the 2 LOC copies added, is recovered 10/10 in every set.

**Reference trees are not re-tested.** They contain only the 22 records, so complete groups cannot be scored on them without a 27-leaf rebuild (not done).

**Wording to use:**
- "the NPIPA|NPIPB split and the Iso-Seq paralog groups B3-5, B6-9 and B12/13 of Dishuck 2025 (each a clade in their Fig 2A)";
- not "Dishuck's phylogenetic groups".
- Never name a copy's paralog by its RefSeq symbol.

**Extra check (post hoc) at the paper's 28-paralog level.** The same intron trees agree with the paper's per-copy labels:

| paralog group | members (RefSeq records) | Q + §6jr + §6js | MAFFT |
|---|---|---|---|
| A7/A8 | NPIPA7, NPIPA8 | 30/30 | 9/10 |
| B5 | NPIPB3, LOC128966608 | 30/30 | 8/10 |
| B7 | NPIPB6, NPIPB8, NPIPB9 (excluding NPIPB7 = B9L1) | 24/30 | 8/10 |

B12 vs B13 is not resolved (a 2-copy B13 clade is supported in 10/30 trees). The paper cannot separate them either.

### (b) §6jo — "one identity cut recovers NPIPA/NPIPB in 0/10"

**Stands.** It is a negative on exact recovery. Adding group members cannot turn a failed exact partition into a success.

Relabel the two keep-1 three-way cuts:
- **RefSeq {B3, B4, B5, B11, B12, B13}** = Dishuck {B5, B4, B35L1, B11, B12/B13}. This is the β-helix clade (Fig 2A, 3.1 mya) minus LOC128966608 and LOC124907834.
- **RefSeq {B6-B9, B10P, B15}** = Dishuck {B7 ×3, B9L1, B10, B15}. This is the signal-peptide clade (2.6 mya) minus the two B15 LOC copies.

### (c) §6jg / §6ji — "de novo: NPIP in 5 families", "all 22 present"

**The membership statements stand**, because they are keyed by record coordinates.

**Relabel:**
- de novo "B5+B9" = B35L1 + B7;
- union "B4+B5" = B4 + B35L1;
- §6jg level-B "A7+A8" is exactly the paper's A7/A8 paralog (worth quoting);
- "A6+A9" are sister tips in Fig 2A.

**Scope "all 22 present":** it means 22 of 27 CHM13 NPIP copies.
- The other 5 lie outside the input windows (truth records ±50 kb).
- No guided node and no de novo copy overlaps them (`CK` check against `LIT/dn_default` and `LIT/dn_union` `*.copies.tsv` and `LIT/guided.loci.tsv`).
- Testing them needs new windows and reruns (> 10 min, not done).

### Other users of these truth rows

Artifacts that read the NPIP truth rows (§6jh, §6kl DN rows, §6kq panels) have correct membership for the 22 records, but they carry the same 5-copy gap. Separately, §6kq already reports B6-9 as *unsupported* once HG002 tips are added. That result does not depend on this check, but it bears on how strongly to quote B6-9.

---

## §4 Caveats

- **One haplotype vs 169.**
  - The paper defines its groups over 4,665 copies. CHM13 carries 19 of the 28 paralog groups; A1L1, A3, A6/A9, B3, B6, B8, B9, B16 and B345L1 are absent.
  - A clade over CHM13 tips is necessary evidence for a group, not sufficient.
  - The Iso-Seq groups are pragmatic read-assignment aggregates. A6-9, B3-5 and B6-9 are not SH-aLRT-defined units, and Fig 2A shows no support values for them.
- **The trees are not independent confirmations.**
  - The 4 CHM13 reference trees read one 22-sequence input FASTA (md5 07074ea9): §6jp P1 MAFFT, §6jp Q, §6jr intron and §6js intron.
  - Three of them (Q, §6jr, §6js) also share one projected alignment (md5 3efe6ed6).
  - The 40 leave-out trees are 4 method variants over the same 10 seed splits of the same genome; half_2 has identical input in Q and §6jr.
  - Report this as one experiment with robustness variants, not as 4 or 40 confirmations.
- **Limits of reading the figures.**
  - Figures are raster, labels were read by eye, and strands were not read.
  - The identities of B1, B14 and the three B15 copies rest on Fig 1A, Fig 4D and the text. Fig 3A covers only 16p.
  - The paper uses two naming layers: phylogenetic names (B35L1, B9L1) in Figs 2, 3 and 4A, and GRCh38-best-match names ("B3", "B9") in Fig 1A and Fig 4B. The selection text on p.5 mixes both layers. Quote the phylogenetic names.
- **Post hoc.** The rescoring and the 28-paralog checks were done after the results were known. If a 27-copy truth is adopted, disclose it as a truth amendment, as was done for TBC1D3.
- **Preprint.** The paper is a bioRxiv v1 preprint. Per-copy tables (Supplementary Table S2 and others) are not in the PDF.
