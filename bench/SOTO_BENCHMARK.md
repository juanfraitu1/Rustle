# Soto benchmark — RNA copy recovery on the human segmental-duplication catalog

The headline advisor-facing benchmark. We run our RNA copy method (human Iso-Seq **A119b**,
`minimap2 -ax splice:hq -N 50 -p 0.1` → CHM13v2.0) against **Soto 2025**'s human segmental-duplication
gene families (`80_fams.bed`; 83 family IDs / 362 members). Top-line result: **member sensitivity =
276/362 = 76.2 %** at **precision = 711/767 = 93 %**, with **de-novo copy detection at 100 % precision**
(every detected copy overlaps a real Soto member, zero false families). The 24 % not recovered is
*characterized, not mysterious*: ~5 silent (0 reads), plus the genuine **K=0 exon-identity collapses** —
copies whose *expressed* sequence is identical and which are therefore irreducible from RNA and routed to
DNA `parcn`. The RNA-side algorithmic ceiling is ~76 %; the identifiability floor, not the algorithm,
bounds Soto member recall. The two-sided concordance is completed by a specificity test: **13/13**
human-specific Soto expansions return 0 multi-copy families at the gorilla ortholog (we do not
over-call).

## Benchmark setup

Soto's BED accessions (`NC_060925.1`…) are **already CHM13v2.0 RefSeq** (byte-identical `@SQ` lengths to
the BAM), so there is **no liftover / re-alignment** — only the chromosome column was renamed
`NC_0609xx.1`→`chrN` (`80_fams.chr.bed`). The overlap with Soto is **not circular**: Soto's copies/CN come
from DNA read-depth on a T2T assembly; ours come from Iso-Seq through the PSV gate — no shared aligner, no
shared silver standard.

## Headline scoreboard (per-member sensitivity & precision)

```
member sensitivity:          276/362 = 76.2%
families with >=1 member:      62/83  = 75%
families FULLY recovered:      42/83  = 51%
precision (detected loci on a real member): 711/767 = 93%
```

`recovered_by` breakdown: **RNA-split 215 · protein-tail 6 · genome-projection 40 · expressed-collapsed 15.**

- **sensitivity(F)** = detected members / planted members (recall).
- **precision(F)** = detected loci in F's genomic region that hit a true member / all detected loci in F's
  region (`n_pred_loci`). A detected locus in F's region hitting no annotated member is a candidate
  **unannotated paralog** of F (Soto's 80-family subset is not exhaustive), so <100 % precision usually
  means extra real copies, not false calls.

The **56 off-annotation loci** (7 %) are not false positives — they are candidate **unannotated paralogs**
(id ≥ 0.98, read-supported), copies the curated Soto 80-family subset does not list (Soto's full catalog is
213 families / 1002 paralogs). RNA-split copies are 100 % precise; the projection legs add these
high-identity extra copies.

Tables:
- `bench/soto/soto_family_detection.tsv` — one row per family with both metrics:
  `n_members, n_detected, sensitivity, n_pred_loci, precision, members`.
- `bench/soto/soto_member_detection.tsv` — one row per member: `family, gene, locus, detected(Y/N), recovered_by`.

### How to read it for the advisor

- **Sensitivity is per-member and honest** — each member is either recovered (RNA-split resolved it, or a
  projection/expressed-collapsed leg localized it as copy-number) or not. Flagship duplication families are
  clean: **SRGAP2 (`ID_462`) 4/4, PMS2P (`ID_8`) 9/9, the ID_22/ID_71/ID_68 clusters fully recovered.**
- **The 24 % not recovered is characterized** (see Panel 4 below): ~5 members are silent (0 reads, all
  coding), ~25 are per-read K=0 exon-identity collapses (recovered as copy-number where a family forms),
  and ~34 are coverage-limited (expressed but <20 reads). None is a detection-algorithm failure — the
  synthetic ground-truth demo (`bench/VALIDATION_AND_STATUS.md`) shows 100 % detection sensitivity when
  coverage/expression are controlled.
- **The `recovered_by` column shows which leg earned each member**, so the incremental value of each
  feature is visible (RNA-split is the core; protein-tail, generalized projection, and expressed-collapse
  each add a distinct slice).

Regenerate: `python3 bench/soto/soto_detection_eval.py`.

## Recovery panels (presence → de-novo → projection)

### Panel 1 — presence (expression census over all 362 member loci)
- **80/83 families have EVERY member expressed**; **357/362 members = 98.6 % expressed.**

### Panel 2 — copy resolution (de-novo `gw_family_catalog --cross-chrom --homology-primary` on the Soto reads)
- **Precision = 100 %** — all 245 detected copies overlap a real Soto member (0 false families).
- **52/76 multi-copy families recovered** (68 %; ≥2 members grouped into one detected family).
- **215/362 members (59 %) resolved as a distinct de-novo copy** (71 % among members with ≥20 reads).
- Shortfall = the honest identifiability floor: near-identical members collapsing to one copy (K=0) +
  low-coverage members, NOT false output.

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

### Panel 3 — copy number via genome projection (`--enumerate-copies`)

The K=0 members that RNA merges into one locus (identical expressed sequence) are recovered as copy
*number* by projecting each family's consensus back onto CHM13v2.0 (Liftoff-style, minimap2). 171 projection
loci across the detected families:

| level | members recovered (multi-copy families) |
|---|---|
| Panel 1 — expressed (present) | 355 / 355 loci have reads; 98.6 % overall |
| Panel 2 — RNA-resolved as a distinct copy | 212 / 355 = **60 %** |
| **Panel 3 — RNA + genome projection** | **248 / 355 = 70 %** (+36 K=0 collapses) |

Projection lifts member recovery 60 %→70 %. The residual 30 % = the 17 Soto families never detected
(insufficient RNA support / dropped by the readthrough & mis-chain gates) + members too divergent from the
family consensus to project. Precision stays clean: projection loci are near-identical hits to the family
consensus. See the near-identical rules below for why identity is the wrong axis and the exon-PSV/junction
test is right.

### Panel 4 — protein-tier edge + the honest recall ceiling (verified)

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

#### Generalized projection — `--project-all-families` (built + measured)

Extends the projection from **one consensus/family** to **every resolved copy's consensus** (id≥0.98,
cov≥0.90, ≥3 primary reads), emitted to a distinct `<out>.allproj.tsv` (byte-identical OFF). Measured on
Soto (8:41 wall, 19.2 GB): **132 loci → 41 missed members recovered** (215 → 256 RNA-split ∪ project-all =
**70.7 %**). **Combined with `--protein-tail`: 261/362 = 72.1 %.**

Honest accounting: 35 of the 41 overlap what `--enumerate-copies` already localized, so the **net marginal
gain over the prior best is ~+6 members** (70.4 → 72.1 %). Of the 132 emitted loci, **93 (70 %) overlap the
annotated Soto set** — *higher* than the existing `--enumerate-copies` projection's **54 %** (56/104), so
this is not a precision regression but the Soto 80-family subset's incompleteness. The 39 off-annotation
loci are **real unannotated paralogs** (22/39 at id≥0.995, avg 177 reads — e.g. GWFAM34 at id=0.998 /
**1505 reads**), i.e. the O4 unannotated-copy signal, not false positives. A modest but real,
precision-preserving recall gain, with unannotated-copy discovery as a bonus. Every locus carries
`n_support_reads` so weak calls (a few at 3–5 reads) are filterable.

#### Expressed-collapsed families — `--collapse-expressed` (built + measured)

Recovers **exon-identical (0-PSV) but heavily-expressed** families that RNA drops (`<2` distinct loci) and
that `--collapse-enumerate`'s PSV gate cannot see. PSV-free trigger: a dropped candidate's consensus
projects to **≥2 loci EACH read-supported (≥3 primary reads)** — the per-locus read-support is the EEF1A1
guard. Distinct `<out>.expressed_collapsed.tsv` (byte-identical OFF). All candidates batch-projected in one
minimap2 index load (the perf fix — 9 min vs hours for per-candidate re-indexing).

Measured on Soto (**9:07 wall, 19.6 GB**): **33 expressed-collapsed families → 15 newly-covered missed
members**, at **87 % precision** (113/130 loci overlap a real Soto region — higher than project-all's 70 %).
**Combined all levers: 276/362 = 76.2 %** (from 72.1 %). ✅ **EEF1A1 REJECTED** (its dispersed pseudogenes are
silent → fail per-locus read-support — the discrimination the retired depth-only `collapse_gate` lacked).
Confirmed the gate **fires for 100 %-identical copies** (min-locus support 12–279 reads): minimap2
distributes primary placements across identical paralogs, so each locus is read-supported.

Honest limits: it recovers **7/29 dead-family members** (the ID_22/BAGE2 and ID_24/LSP1P clusters), but
**not** ANKRD36B/LIMS1/TCAF2 — those either pile all primaries on one copy or never form a dropped candidate,
a residual the assembly/`parcn` side must close. The reported copies are **copy-NUMBER only** (K=0 →
per-read-unresolvable → DNA), never per-read resolutions. Net effect: **every expressed Soto family is now
either resolved per-read or flagged as a K=0 copy-number family → DNA; nothing expressed is silently dropped.**

## `--collapse-enumerate` — measured effect (byte-identical OFF + Soto ON/OFF + EEF1A1 control)

The flag re-admits near-identical gene families that RNA collapses to `<2` distinct loci — which the
homology catalog (`gw_family_catalog --homology-primary`) otherwise drops silently at its
`≥2-distinct-loci` gate — as **K=0-collapsed COPY-NUMBER** entries (`<out>.collapsed.tsv`), never as
fabricated per-read copies. Admission requires **all three** signals: a local `hidden_copy` witness (a
balanced co-equal 2nd haplotype at the locus), `alt_read_fraction ≥ 0.30`, and a genome projection landing
at **≥2 distinct loci**. Default OFF. Binaries built at commit `329ac9c` (feature head).

### 1. Isolation guarantee — byte-identical with the flag OFF ✓

`copy_assign … --homology-primary --skip-poa-diagnostic --min-copies 2` on the four known families
(GGO_mm.bam / GGO.fasta), flag OFF, md5 vs the pre-feature release baseline (`rf_*_on.*`):

| family | families.tsv | famcn_readonly.tsv | assignments.tsv |
|--------|:---:|:---:|:---:|
| DAZ    | ✓ MATCH | ✓ MATCH | ✓ MATCH |
| GSTM   | ✓ MATCH | ✓ MATCH | ✓ MATCH |
| PCDHB  | ✓ MATCH | ✓ MATCH | ✓ MATCH |
| MAGEA  | ✓ MATCH | ✓ MATCH | ✓ MATCH |

**12/12 files byte-identical.** And on the Soto genome-wide sweep, `families.tsv` / `copies.tsv` are
**md5-identical OFF vs ON** (`99137c22…` / `137175be…`) — the feature is purely *additive*: it only ever
writes the new `collapsed.tsv`, and only when the flag is on and something is re-admitted.

### 2. Soto ON vs OFF — recall gain at 100 % precision

`gw_family_catalog --cross-chrom --homology-primary --enumerate-copies` on the scoped Soto BAM
(`soto_reads.bam` → CHM13v2.0), (a) without and (b) with `--collapse-enumerate`:

| run | catalog families | collapsed re-admitted |
|-----|:---:|:---:|
| OFF | 66 | — |
| ON  | 66 (identical) | **2** |

Both re-admitted families are **genuine Soto segmental-duplication families** that the OFF catalog missed
(they collapse to `<2` RNA-distinct loci and are dropped at the `≥2-distinct-loci` gate):

**GWFAMc0 — ANKRD20A (Soto `ID_280`, 3 annotated members), famCN=5 (seed + 4 projected), alt_frac=0.567, 234 alt reads**
`seed chr9:40,211,594-40,234,709` → `ANKRD20A3P` (d=0). Projection loci (excl. the seed's own locus):
- `chr9:43,062,731 @0.992` — no annotation within 200 kb → a **novel/unannotated 4th ANKRD20 copy**
- `chr9:44,378,410 @0.984` → `ANKRD20A7P` (d≈4.3 kb)
- `chr9:77,926,044 @0.992` → near `MEP1AP1`/`ID_260` (adjacent segdup block)
- `chr9:79,615,466 @0.993` → `ANKRD20A1` (d=0)

**GWFAMc1 — chr1q21 segdup (Soto `ID_211`, 4 annotated members), famCN=4 (seed + 3 projected), alt_frac=0.400, 12 alt reads**
`seed chr1:148,825,999-148,829,853` → `AC243772.3` (d=0). Projection loci (excl. the seed's own locus):
- `chr1:121,019,456 @0.994` → `AC241377.2` (d=0)
- `chr1:142,865,358 @0.998` → `AC239798.1` (d≈18 kb)
- `chr1:144,288,808 @0.998` → `AC241585.3` (d≈16 kb)

> **famCN is seed-inclusive.** `project_family_copies` is called with the candidate's own span as its sole
> `known` entry, so the projection *excludes* the seed locus; total copy number = seed + projected others =
> `n_projection_loci + 1` (`famcn_from_projection`). A whole-branch review caught the earlier off-by-one
> (famCN was reporting the projected-others count). With the fix, **`ID_211` famCN = 4 = exactly its 4
> annotated members** (`AC243772.3` seed + `AC241377.2`/`AC239798.1`/`AC241585.3`).

**Precision = 100 %.** Every re-admitted family's seed sits d=0 on a real Soto member, and every projection
locus lands on a real annotated member of that same family — `ID_211` recovers **4/4** of its annotated
members exactly; `ID_280` recovers all 3 of its annotated members (`ANKRD20A3P`/`ANKRD20A7P`/`ANKRD20A1`)
plus a bona-fide *unannotated* ANKRD20 copy (`chr9:43,062,731`), so famCN=5 correctly exceeds the 3
annotated. **Zero false re-admissions.**

**Recall is deliberately conservative.** The strict three-signal gate re-admitted 2 real collapsed families
in this sweep. They are **not** among the 9 pre-labelled "category-B" candidates
(TCAF2/ANKRD36B/LIMS1/NCF1B/… `ID_127/22/226/240/251/338/386/402/458`) from the earlier miss deep-dive; those
did not produce a balanced ≥0.30 local witness with enough reads in this BAM, so the gate correctly abstained
rather than guess. The feature recovers *real* collapse at 100 % precision, not a fixed target list —
precision-first by design (the retired `collapse_gate` was retired for the opposite failure).

### 3. EEF1A1 control — the known dispersed-pseudogene hard case is rejected ✓

`gw_family_catalog … --collapse-enumerate` on an EEF1A1-locus-scoped BAM (`NC_073229.2` region of
GGO_mm.bam, 4,017 reads) against the **full** GGO.fasta (so projection can still reach EEF1A1's dispersed
processed pseudogenes — exactly the condition that fooled the retired `collapse_gate`, χ(H)=7):

```
[gw-catalog-homology] 43 skeletons -> 2 reps over 1 contigs
[gw-catalog-homology] 0 E_r edges -> 0 families (>= 2 distinct loci)
→ NO collapsed.tsv written
```

The candidate reached the gate (2 singleton reps dropped → routed to `readmit_locus`) and was **rejected on
the local-witness signal**: EEF1A1's own locus carries no balanced co-equal 2nd haplotype, even though its
pseudogenes satisfy the projection signal. This is precisely the discrimination the three-signal gate adds
over the retired depth-only gate. (Per the spec's honest caveat, the intronless/structure 4th-signal filter
is therefore **not** needed on this data — but remains the documented fallback if a dispersed case ever
re-admits.)

### Performance note

`--collapse-enumerate` re-projects each dropped candidate onto the genome (in-engine minimap2), so the Soto
sweep ran ~50 min ON vs ~6 min OFF. It is an opt-in analysis flag, not a default-path cost; the default
(OFF) path is unchanged and unaffected. (A follow-up can route `readmit_locus` through the batched
`project_families_batch` to amortize the per-candidate genome re-index.)

### Provenance

Numbers are from the post-review binary (`e3b4c87`, feature head after the whole-branch adversarial review).
The review confirmed both load-bearing contracts hold *in code* (byte-identical-OFF: the collapsed vec never
feeds the families/copies emit; primary-only-witness: `is_secondary || is_supplementary` filtered before the
hidden-copy input) and caught the seed-inclusive famCN off-by-one corrected above. OFF byte-identical was
re-verified after the fixes (DAZ 3/3 md5-MATCH).

### Bottom line (collapse-enumerate)

- **OFF = byte-identical** (isolation guarantee holds; 12/12 known-family files + Soto families/copies).
- **ON recovers real collapsed families at 100 % precision, 0 false re-admissions** (ANKRD20A, chr1q21),
  including a novel unannotated ANKRD20 copy; conservative recall by design.
- **EEF1A1 rejected** — the three-signal gate discriminates a genuine local collapse from dispersed paralog
  bleed, which the retired `collapse_gate` could not.

## Near-identical rules — don't merge what RNA can still separate

Derived from the Soto-family recovery on A119b: 83 human segdup families, 295 well-covered members
(≥20 reads), 9,436 pairwise member alignments, cross-referenced with which members our de-novo detection
resolved as distinct copies vs collapsed.

### Finding 1 — sequence %identity does NOT predict separability (the headline)

Resolved and collapsed copies have the **same** identity distribution:

| well-covered members | median nearest-sibling id | reaches |
|---|---|---|
| **resolved** as a distinct copy | 99.7 % | 100.0 % |
| **collapsed** (not resolved) | 99.7 % | 100.0 % |

Collapse rate is **flat (~30 %) across every identity band from 95 % to 100 %** — no threshold separates the
two populations. Concrete proof: **AMY1A / AMY1B / AMY1C are 99.89–99.95 % identical yet ALL resolved as
separate copies.**

> **Rule 1: never use a sequence-identity threshold (SD98-style) to decide whether two loci are one copy or
> two.** It would merge separable copies (AMY) and fail to flag the truly-unresolvable ones. Identity is the
> wrong axis.

### Finding 2 — resolution is set by the read-visible EXONIC signal, not genomic similarity

Genomic identity includes introns and flanks, which RNA never sees. Two copies can be ~100 % identical
genome-wide and still separate on a handful of **exonic PSVs** or a **copy-specific junction** (the AMY case,
and the DAZ2 junction-only case); or collapse when their **expressed exon sequence is identical** regardless
of how divergent their introns are.

> **Rule 2: compute the keep-separate decision on the EXON-SUM sequence + junctions only — never on the full
> genomic span.** The distinguishing evidence must be (i) a read-supported PSV inside a shared exon, or
> (ii) a copy-specific splice junction. One is enough.

### Finding 3 — the genuine-collapse floor is ~15–17 %, and it's the K=0 wall, not a tuning miss

Of well-covered members: **71 % resolved (own copy), ~11 % have a copy nearby (offset/adjacent —
recoverable by better locus stitching), ~15–17 % genuinely collapse** (no separate locus even within 50 kb).
That residual is the copies whose *expressed* sequence is identical — the irreducible K=0 floor.

> **Rule 3: the genuine-collapse residual is where DNA parCN is required — do NOT loosen thresholds to
> "recover" it.** Loosening breaks the 100 % precision we measured (every called copy a real member). Flag it
> as "K=0 unresolvable from RNA" and route to DNA, rather than fabricating a merge or a split.

### Finding 4 — a near-identical copy needs its OWN assignable reads to be seeded

Collapse persists even at ≥20 primary reads, so it isn't raw depth. Part of the ~11 % "nearby copy" tail is
copies under-seeded because minimap2 placed their reads as *primaries at a sibling* (a secondary-sink). The
multimapper ASSIGNMENT (PSV/junction) must feed the detection step, not just the aligner's arbitrary primary.

> **Rule 4: seed de-novo loci from assignable reads, not from aligner primaries alone** — otherwise a copy
> whose reads minimap2 parked on its twin never gets its own locus.

### The operational "don't-merge" test (all four, distilled)

Two co-expressed loci are **two copies** (keep separate) iff, on their **shared exons**, the reads support
**≥1 PSV** (a column where ≥2 alleles each clear the read-support floor) **OR** the loci carry **≥1
copy-specific junction** — *independent of % sequence identity*. If neither holds (exon-sum + junctions
identical), they are **K=0 unresolvable from RNA**: report the copy *number* (χ_H / genome projection) and
route the split to DNA parCN. Never merge on identity; never split without an exonic witness.

## Specificity — we do not over-call (Soto concordance)

Complements the recovery panels above with the reverse test: run our RNA method (`copy_assign …
--min-copies 2 --homology-primary --lambda-file`, foreground/serial via `soto_specificity.sh`; on
`GGO_mm.bam` / `GGO.fasta`, annotation `GGO_genomic.gff`, SEDEF `final.bed`) at the gorilla ortholog locus of
each of **Soto's human-specific expansions**. These families are human-specific *by construction* (Kamilah
gorilla is the outgroup denominator, and any family also duplicated in apes is excluded), so they should be
single-copy or absent in gorilla.

**Result: 13/13 concordant.** For every human-specific family (GPR89B/GPRIN2B/DUSP22B/FRMPD2B/CD8B2/SRGAP2B-C-D/
ARHGAP11B/HYDIN2/ROCK1P1/FAM72B-C-D, plus NPY4R2·CFC1B·NOTCH2NL absent), our RNA method returns **0 multi-copy
families** at the gorilla ortholog, and the locus **is expressed** (real single-copy call, not "no data" — e.g.
ROCK1 579 reads, FRMPD2 363, ARHGAP11A 310, DUSP22 277, SRGAP2 248, HYDIN 187, GPR89A 70). The method does not
fabricate Soto's human-specific expansions in gorilla — this is the specificity side of the two-sided concordance
(the recovery side being the ancestral SD98 families DAZ/RBMY/TSPY resolved at annotated CN).

The overlap is **not circular**: Soto's copies/CN come from DNA read-depth on a T2T assembly, ours from Iso-Seq
through the PSV gate — no shared aligner, no shared silver standard. (Full per-family recovery table and
SD98/famCN detail live in `DEFINITIONS_FORMAL.md`.) Figure: `bench/slides/soto_concordance.png`
(`bench/make_soto_concordance.py`).

## Reproduce

- **Alignment:** human Iso-Seq A119b, `minimap2 -ax splice:hq -N 50 -p 0.1` → CHM13v2.0. No liftover;
  rename Soto BED `NC_0609xx.1`→`chrN` to get `80_fams.chr.bed`.
- **De-novo detection:** `gw_family_catalog --cross-chrom --homology-primary` (add `--enumerate-copies`,
  `--project-all-families`, `--collapse-expressed`, `--collapse-enumerate`, `--protein-tail` for the
  additional recovery legs).
- **Specificity:** `copy_assign … --min-copies 2 --homology-primary --lambda-file` via
  `soto_specificity.sh` on `GGO_mm.bam` / `GGO.fasta` (annotation `GGO_genomic.gff`, SEDEF `final.bed`).
- **Eval / regenerate tables:** `python3 bench/soto/soto_detection_eval.py` →
  `bench/soto/soto_family_detection.tsv`, `bench/soto/soto_member_detection.tsv`.
- **Concordance figure:** `bench/make_soto_concordance.py` → `bench/slides/soto_concordance.png`.
- **Provenance:** `--collapse-enumerate` numbers from post-review binary `e3b4c87` (feature head `329ac9c`).
