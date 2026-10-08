# Pre-registration: families as a partition of duplicon units, genes as paths — does it predict Soto's multi-family genes? (2026-09-30, KEY=unitcover)

Written before any statistic below was computed. Human CHM13 v1.0, Soto's 2,334-gene table (S1C). Nothing in `src/` or `tools/`
changes. Script: `bench/soto_m2/soto_m2_unit_cover.py`, written after this file is committed.

## 1. Question

Our clustering keeps families as a strict partition of genes; Soto's table is a cover (149 genes in two or more families). Every attempt
to put genes in several families inside the clustering failed (r845, r846, r1165-1169). The proposed structure moves the partition down
one level: **each duplicon (ancestral duplication unit) belongs to exactly one family, and a gene is a path through duplicons, so its
families are the owners of the duplicons its exons touch.** If this structure is right, it should predict which of Soto's genes sit in
several families, and which families, from the single-family genes alone.

## 2. Inputs (frozen)

- Genes, exons, Soto families: the meeting page's `families_cn.json` (`genes`: S1C's 2,334 genes, CAT v4 exons, CHM13 v1.0; `sf` = the
  gene's S1C Family IDs). `Unassigned_*` labels are not families.
- Units: DupMasker duplicons, `winloci_data/duplicons/chm13.draft_v1.0_plus38Y_dupmasker_colors.bed` (CHM13 v1.0, same coordinates as
  the exons); every overlapping record counts.
- Split: `bench/soto/soto_split_2026-09-29.tsv` (gene -> dev / held-out).

## 3. The structure (fixed now)

- **Clean gene:** in exactly one Soto family. **Multi gene:** in two or more (expected 149).
- **Units of a gene `U(g)`:** the duplicon IDs overlapping at least 1 bp of its exons.
- **Owner of a duplicon `D`:** the Soto family whose clean genes have the most exonic bases on `D` (ties: the lower family number);
  a duplicon on no clean gene's exons has no owner. When the gene being predicted is itself clean, it is left out of the ownership
  count (leave-one-out). Owners form a partition of the owned duplicons.
- **Prediction:** `P(g)` = the owners of the owned duplicons in `U(g)`. Truth `S(g)` = the gene's Soto families.

## 4. Statistic, baselines, decision rule

- **Statistic:** mean Jaccard |S ∩ P| / |S ∪ P| over the multi genes (empty union counts 0), per split half.
- **Location baseline:** `P_loc(g)` = the Soto families of clean genes whose gene span overlaps g's; if none, the family of the nearest
  clean gene on the same chromosome. Same statistic.
- **Permutation null:** owners shuffled among the owned duplicons (each family keeps its number of owned duplicons), 1,000 permutations,
  seed 20260930; p = (1 + #{permuted mean >= observed}) / 1,001.
- **Decision (held-out half decides; the dev half reported beside it; if they disagree the verdict is SPLIT):**
  - **HOLDS:** p < 0.01 and the mean Jaccard is above the location baseline's.
  - **PARTIAL:** p < 0.01, not above the location baseline.
  - **FAILS:** p >= 0.01.

## 5. Secondary (descriptive)

- Exact-set rate, recall |S ∩ P| / |S| and precision |S ∩ P| / |P| for multi genes; the same for the location baseline.
- Clean genes (leave-one-out): share predicted exactly {own family}, share predicted two or more families (false multi), share empty.
- The NPIP side (families ID_149-ID_155): each multi gene's S and P listed.

## 6. Seen before (disclosed)

- ID_151's five shared genes and their Soto families (PKD1P6-NPIPP1, AC126755.6, PDXDC2P-NPIPB14P, AP001120.2 in ID_149 + ID_151-155;
  MSTRG.2119 in ID_41, ID_151-153, ID_169); their SD98 regions and duplicon tracks on the page.
- KEY=cnduplicon (FOLLOWS): within sequence families, Soto's copy-number cut follows duplicon composition. KEY=npipfusion (dev only):
  NPIP fusion partners sit on co-duplicated duplicons. The PKD1P-NPIP fusion duplicon strings (SD9474 ... SD9613 | SD9449, SD9443).
- No ownership, prediction or Jaccard has been computed.
- Expectation written now: the NPIP-core duplicons (e.g. SD9443) are shared by every NPIP subfamily, which Soto separates by copy number,
  so a strict partition of duplicons may under-predict NPIP's multi genes.

## 7. Result

(Filled in after the run, below this line, without editing anything above.)

Run on 2026-09-30 after this file was committed (`cedce819`, sha1 `39ae6af4942b045551d5dbfbcb03a18d71b1f919`);
`bench/soto_m2/soto_m2_unit_cover.py`, 1 s, light lock. Output `docs/archive/2026-09/SOTO_UNIT_COVER_2026-09-30.md`, per gene
`docs/SOTO_UNIT_COVER_2026-09-30.tsv`.

**VERDICT (section 4): HOLDS** — held-out and dev alike.

| half | multi genes | mean Jaccard, structure | mean Jaccard, location baseline | permutation p |
|---|---|---|---|---|
| dev | 81 | 0.439 | 0.221 | 0.001 (the floor at 1,000 permutations) |
| held-out | 68 | 0.473 | 0.277 | 0.001 |

2,071 clean and 149 multi genes; 1,546 duplicons owned by 392 families.

Secondary: over all 149 multi genes the structure gets the exact family set for 10.7% (location 5.4%), recall 0.555 (0.286), precision
0.718 (0.426). Clean genes, leave-one-out: 54.6% predicted exactly their own family, 30.3% two or more families, 2.0% none. NPIP side: the
shared genes get 1-3 of their 5-6 Soto families (PKD1P6-NPIPP1 and AC126755.6: ID_149 + ID_154 of six; AP001120.2: ID_154 only), as
expected in section 6: the NPIP subfamilies share their core duplicons and Soto separates them by copy number, which a partition of
duplicons cannot express.

**Reading.** The structure carries real signal: it predicts Soto's multi-family memberships from single-family genes alone far better
than location or chance, on the held-out half. It is not a complete account: exact sets are rare, a third of single-family genes touch a
second family's duplicon (the 1 bp touch rule counts slivers), and where families differ by copy number rather than by duplicon (NPIP)
the duplicon level is too coarse. A finer unit (duplicon x copy-number class) or a usage rule that ignores slivers would be the next
pre-registration, not a change to this one.
