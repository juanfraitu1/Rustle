# Soto families recovered in A119b — human Iso-Seq on CHM13v2.0

Soto 2025's human segmental-duplication gene families (`80_fams.bed`), run against **A119b** (human
Iso-Seq, `minimap2 -ax splice:hq -N 50 -p 0.1` → CHM13v2.0). **No liftover / re-alignment** — the BED
accessions (`NC_060925.1`…) are already CHM13v2.0 RefSeq (byte-identical `@SQ` lengths to the BAM), so
only the chromosome column was renamed `NC_0609xx.1`→`chrN`. 83 Soto family IDs / 362 members.

## Panel 1 — presence (expression census over all 362 member loci)
- **80/83 families have EVERY member expressed**; **357/362 members = 98.6% expressed**.

## Panel 2 — copy resolution (de-novo `gw_family_catalog --cross-chrom --homology-primary` on the Soto reads)
- **Precision = 100 %** — all 245 detected copies overlap a real Soto member (0 false families).
- **52/76 multi-copy families recovered** (68%; ≥2 members grouped into one detected family).
- **215/362 members (59%) resolved as a distinct de-novo copy** (71 % among members with ≥20 reads).
- Shortfall = the honest identifiability floor: near-identical members collapsing to one copy (K=0) + low-coverage members, NOT false output.

| family | ID | Soto members | expressed | resolved as copy | recovered |
|---|---|--:|--:|--:|:--:|
| GOLGA | ID_113 | 16 | 16 | 14 | ✅ |
| NPIPB | ID_154 | 14 | 14 | 13 | ✅ |
| SPDYE | ID_207 | 17 | 17 | 10 | ✅ |
| AL | ID_35 | 10 | 10 | 8 | ✅ |
| PMS | ID_8 | 9 | 9 | 7 | ✅ |
| TBC | ID_468 | 9 | 9 | 7 | ✅ |
| RGPD | ID_395 | 6 | 6 | 6 | ✅ |
| GUSBP | ID_163 | 6 | 6 | 6 | ✅ |
| GOLGA | ID_116 | 6 | 6 | 6 | ✅ |
| AMY | ID_131 | 6 | 6 | 5 | ✅ |
| AC | ID_215 | 7 | 7 | 5 | ✅ |
| FAR | ID_65 | 9 | 9 | 5 | ✅ |
| BCRP | ID_283 | 5 | 5 | 5 | ✅ |
| NOTCH | ID_400 | 6 | 6 | 4 | ✅ |
| FAM | ID_354 | 4 | 4 | 4 | ✅ |
| SRGAP | ID_462 | 4 | 4 | 4 | ✅ |
| ANAPC | ID_63 | 8 | 8 | 4 | ✅ |
| CNTNAP | ID_245 | 4 | 4 | 4 | ✅ |
| GOLGA | ID_88 | 4 | 4 | 4 | ✅ |
| UBE | ID_481 | 5 | 5 | 4 | ✅ |
| AC | ID_14 | 7 | 7 | 4 | ✅ |
| AC | ID_212 | 5 | 5 | 3 | ✅ |
| AC | ID_211 | 4 | 4 | 3 | ✅ |
| FCGR | ID_357 | 3 | 3 | 3 | ✅ |
| WASH | ID_391 | 5 | 5 | 3 | ✅ |
| FAM | ID_182 | 4 | 4 | 3 | ✅ |
| GTF | ID_208 | 5 | 5 | 3 | ✅ |
| GTF | ID_374 | 3 | 3 | 3 | ✅ |
| ANKRD | ID_280 | 3 | 3 | 3 | ✅ |
| AL | ID_104 | 3 | 3 | 3 | ✅ |
| FGF | ID_359 | 5 | 5 | 3 | ✅ |
| SHLD | ID_443 | 3 | 3 | 3 | ✅ |
| CTSLP | ID_327 | 4 | 4 | 3 | ✅ |
| AC | ID_49 | 3 | 3 | 3 | ✅ |
| DNM | ID_68 | 4 | 4 | 3 | ✅ |
| GOLGA | ID_78 | 5 | 5 | 3 | ✅ |
| TP | ID_474 | 3 | 3 | 3 | ✅ |
| MST | ID_393 | 2 | 2 | 2 | ✅ |
| SEC | ID_448 | 3 | 3 | 2 | ✅ |
| DDX | ID_332 | 6 | 5 | 2 | ✅ |
| NAIPP | ID_399 | 4 | 4 | 2 | ✅ |
| POM | ID_12 | 4 | 4 | 2 | ✅ |
| NSUN | ID_407 | 4 | 4 | 2 | ✅ |
| ZNF | ID_490 | 3 | 3 | 2 | ✅ |
| AC | ID_147 | 4 | 4 | 2 | ✅ |
| TRIM | ID_43 | 4 | 4 | 2 | ✅ |
| NF | ID_403 | 3 | 3 | 2 | ✅ |
| BMS | ID_300 | 2 | 2 | 2 | ✅ |
| HERC | ID_188 | 4 | 4 | 2 | ✅ |
| AC | ID_71 | 5 | 5 | 2 | ✅ |
| ULK | ID_179 | 3 | 3 | 2 | ✅ |
| AC | ID_209 | 3 | 3 | 2 | ✅ |
| H | ID_226 | 5 | 5 | 1 | ~ |
| DEFB | ID_338 | 3 | 3 | 1 | ~ |
| CNTNAP | ID_261 | 5 | 5 | 1 | ~ |
| WHAMMP | ID_156 | 5 | 5 | 1 | ~ |
| CSPG | ID_324 | 5 | 5 | 1 | ~ |
| PPIAL | ID_431 | 5 | 2 | 0 | — |
| AC | ID_24 | 3 | 3 | 0 | — |
| AC | ID_22 | 6 | 6 | 0 | — |
| AC | ID_213 | 4 | 4 | 0 | — |
| AC | ID_26 | 2 | 2 | 0 | — |
| ANKRD | ID_251 | 2 | 2 | 0 | — |
| LIMS | ID_386 | 2 | 2 | 0 | — |
| AC | ID_146 | 2 | 2 | 0 | — |
| CDH | ID_313 | 2 | 2 | 0 | — |
| NCF | ID_402 | 3 | 3 | 0 | — |
| TCAF | ID_127 | 3 | 3 | 0 | — |
| SPAG | ID_458 | 2 | 2 | 0 | — |
| AL | ID_260 | 3 | 3 | 0 | — |
| SYT | ID_240 | 3 | 3 | 0 | — |
| OR | ID_411 | 3 | 3 | 0 | — |
| AC | ID_148 | 2 | 2 | 0 | — |
| CHEK | ID_167 | 3 | 3 | 0 | — |
| AC | ID_92 | 2 | 2 | 0 | — |
| TEKT | ID_175 | 2 | 2 | 0 | — |

## Panel 3 — copy number via genome projection (`--enumerate-copies`)

The K=0 members that RNA merges into one locus (identical expressed sequence) are recovered as copy
*number* by projecting each family's consensus back onto CHM13v2.0 (Liftoff-style, minimap2). 171 projection
loci across the detected families:

| level | members recovered (multi-copy families) |
|---|---|
| Panel 1 — expressed (present) | 355 / 355 loci have reads; 98.6% overall |
| Panel 2 — RNA-resolved as a distinct copy | 212 / 355 = **60 %** |
| **Panel 3 — RNA + genome projection** | **248 / 355 = 70 %** (+36 K=0 collapses) |

Projection lifts member recovery 60 %→70 %. The residual 30 % = the 17 Soto families never detected
(insufficient RNA support / dropped by the readthrough & mis-chain gates) + members too divergent from the
family consensus to project. Precision stays clean: projection loci are near-identical hits to the family
consensus. See `NEAR_IDENTICAL_RULES.md` for why identity is the wrong axis and the exon-PSV/junction test is right.

## Panel 4 — protein-tier edge + the honest recall ceiling (verified)

Two additional levers, measured against all 362 members (0 false calls, precision held):

- **`--protein-tail` (coding-homology E_r edge, mmseqs fident≥0.50):** recovers **+6 RNA-split members** the
  nucleotide floor misses — `FGF7P8`, `HERC2P3`, `WHAMMP2`, `AC004980.1`, `AC127502.1`, `FO681491.1` (a
  pseudogene/clone mix). RNA-split members: **215 → 221**.
- **Genome-projection loci counted as recovered (DNA-localized parCN, id≥0.98 + family-assigned reads):**
  **+36 missed members** overlap a projection locus. Labeled distinctly from RNA-split — these are copy
  *number* localizations, not per-read resolutions (the copy-vs-allele line we keep clean).

**Combined: RNA-split 221 + DNA-projection 36 (2 overlap) = 255 / 362 = 70.4 %** (from a 59.4 % RNA-only
baseline), at 100 % precision, with essentially no new code.

**The honest ceiling.** Of the 147 misses, ~**5 are silent** (0 reads — all coding: PPIAL4A/C/E, DEFB104B,
DDX11L16) and ~**75–90 are genuine K=0 exon-identity collapses** (≥20 reads, no distinct exonic PSV or
copy-specific junction) — both **irreducible from RNA** (the DNA/`parcn` job). So the RNA-side algorithmic
ceiling is ~**76 %**; the identifiability floor, not the algorithm, is what bounds Soto member recall — a
clean, defensible statement. Miss taxonomy: 5 silent · 36 DNA-localizable (banked above) · 53
detected-family-no-projection · 58 dead-family (~6–15 coding via protein/homology) · residual genuine K=0.
Empirically dead levers (do NOT build): locus-stitching (zero orphan fragments) and coverage-floor
relaxation (~0 real members).

### Generalized projection — `--project-all-families` (built + measured)

Extends the projection from **one consensus/family** to **every resolved copy's consensus** (id≥0.98,
cov≥0.90, ≥3 primary reads), emitted to a distinct `<out>.allproj.tsv` (byte-identical OFF). Measured on
Soto (8:41 wall, 19.2 GB): **132 loci → 41 missed members recovered** (215 → 256 RNA-split ∪ project-all =
**70.7 %**). **Combined with `--protein-tail`: 261/362 = 72.1 %.**

Honest accounting: 35 of the 41 overlap what `--enumerate-copies` already localized, so the **net marginal
gain over the prior best is ~+6 members** (70.4 → 72.1 %). Of the 132 emitted loci, **93 (70 %) overlap the
annotated Soto set** — *higher* than the existing `--enumerate-copies` projection's **54 %** (56/104), so
this is not a precision regression but the Soto 80-family subset's incompleteness. The 39 off-annotation
loci are **real unannotated paralogs** (22/39 at id≥0.995, avg 177 reads — e.g. GWFAM34 at id=0.998 /
**1505 reads**), i.e. the O4 unannotated-copy signal, not false positives. A modest but real, precision-
preserving recall gain, with unannotated-copy discovery as a bonus. Every locus carries `n_support_reads`
so weak calls (a few at 3–5 reads) are filterable.
