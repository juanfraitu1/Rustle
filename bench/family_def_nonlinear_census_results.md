# Non-linear (tandem-repeat CNV) family census — adversarial verdicts

**Claim under test:** In 10 "K-frontier" families the copies are near-identical in sequence
(cDNA id ≥ 0.98) yet differ in **tandem-repeat unit count**, so the SNP/indel-PSV model
(a *linear* variation graph) collapses them but the tandem-repeat copy-number **structure**
(a *cyclic* vg graph) distinguishes them → structure moves the identifiability frontier.

**Method (per family, real code):** minimap2 self-alignment to count internal tandem-repeat
units independently of vg; pairwise alignment of the two most-similar members with the repeat
region masked to test whether the *unique* (non-repeat) sequence is really ≥98% identical;
CDS-vs-UTR via proteins.fa (CDS-only) length tracking; genomic arrangement from genes.bed;
repeat-array length vs a 1–10 kb FLNC read. Default to **refuted** unless numbers support it.

## Per-family table

| Family | Members / chrom | Repeat-count diff real? | Unit counts (key members) | Most-similar pair: unique-region id | Collapsed by SNP? | CDS/UTR | Arrangement | Read-spanable? | Verdict |
|---|---|---|---|---|---|---|---|---|---|
| RCF_102 | 6 / mixed | yes | 6/4/2/2 (~2.3 kb unit) | **47.8%** (CNV-differing pair has 4.5 kb unique flank) | no | UTR | mixed (4 co-located ≤99 kb + 2 dispersed) | yes | **REFUTED** |
| RCF_432 | 12 / 5 chrom | no | only high-id pair 496/567 = equal (369 bp/13.67 u) | 99.7% but **+2 unique-region SNPs**, counts equal | no | CDS | mixed | yes | **REFUTED** |
| RCF_6 | 23 / 10 chrom | yes | 5/4/3/14 (~482 bp unit) | 99.6% core but **704 bp member-specific indel flank** | no | mixed | dispersed (near-id pair 2.98 Mb apart) | yes | **REFUTED** |
| RCF_277 | 6 / 1 chrom | yes | ZF 11/11/9/7/6/5 | 99.95% but pair has **equal ZF count (11=11)**; differs by 294 bp unique N-term | no | CDS | mixed | yes | **REFUTED** |
| RCF_53 | 13 / 1 chrom | yes | 2..7 (~800 bp unit) | 99.65% shared body but **400–620 bp member-private 5′ CDS** | no | CDS | mixed | yes | **REFUTED** |
| RCF_718 | 7 / 1 chrom | yes | 2.7–5.5 (~105–120 bp unit) | **92.1%** (CNV-differing pair); only ≥0.98 pair is byte-identical | no | CDS | co-located | yes | **REFUTED** |
| RCF_611 | 4 / 1 chrom | yes | ~9/~2 units (~30 bp unit) | whole-cDNA **≤91.8%**; unique region still has **4 callable SNPs** | no | CDS | co-located | yes | **REFUTED** |
| RCF_135 | 3 / dispersed | yes | 1/2/3 (743 bp unit) | **0%** homologous unique region (unique = fully member-specific) | no | mixed | dispersed (136.6 Mb apart) | yes | **REFUTED** |
| RCF_13 | 8 / mixed | yes | 4/3/3/2/1.. (~291 bp unit) | co-located pairs cover 31–53%; 3′UTRs **don't align**; high-id pair on diff chrom + 189 bp indel | no | UTR | mixed | yes | **REFUTED** |
| RCF_30 | 9 / 8 chrom | **no** (terminal dup, not array) | 2/2 (terminal direct dup) + 7×1 | 99.34% but **20 scattered SNPs**, repeat-CNV not real | no | mixed | dispersed (≥7.6 Mb) | yes (moot) | **REFUTED** |

## Synthesis

### Genuine K-frontier cases (genuine_kfrontier = true)
**0 of 10.** None. Every family fails the conjunction "unique region ≥98% id **AND** distinguished
only by tandem-repeat count **AND** co-located **AND** read-spanable."

### How the 10 were refuted (the recurring failure modes)
- **Repeat dominates the alignment → inflated id (the central artifact): 6 families**
  (RCF_102, RCF_6, RCF_53, RCF_277, RCF_718, RCF_13). The census ≥0.98 came from best-*local*
  chain identity or the gap-compressed `de` tag, where repeat units tile against each other.
  When the repeat is masked, every **CNV-differing** pair drops to **48–92% global id** and/or
  reveals a large **member-specific unique flank** (e.g. RCF_102's 4.5 kb 5′ flank, RCF_53's
  400–620 bp private 5′ CDS, RCF_277's 294 bp unique N-terminus, RCF_6's 704 bp indel flank).
  Those unique flanks already carry SNP/indel PSVs → **linear PSVs already resolve them**.
- **Where sequence IS collapsed, the repeat counts are EQUAL (dichotomy is fatal): 4 families**
  (RCF_432 496/567 both 13.67 u; RCF_277 both 11 ZF; RCF_718 byte-identical pair; RCF_102
  134~136 both 2 u). The one near-identical pair per family has the *same* unit count, so
  tandem-CNV structure distinguishes **nothing**; they are separated (if at all) by a handful
  of unique-region SNPs, not by structure.
- **Unique regions genuinely differ → SNPs already resolve: all 10** (the union of the two above;
  RCF_611 had 4 callable SNPs even after masking; RCF_30 had 20 scattered SNPs).
- **Dispersed → not a multimapping/copy-assignment problem: 3 families** (RCF_6 near-id pair
  2.98 Mb apart; RCF_135 136.6 Mb apart; RCF_30 9 members on 8 chromosomes, only same-chrom pair
  7.6 Mb apart). No co-location ⇒ reads don't multimap ⇒ there is no copy-assignment ambiguity to
  resolve regardless of structure.
- **Repeat-diff was not real (closest to a vg artifact): 1 family** (RCF_30 — a single head-to-tail
  terminal direct duplication in 2/9 members, periodicity ≈ 0.28 = random-DNA baseline; not a
  variable-count tandem array). RCF_432's apparent 13-vs-0 was also a self-align detection-threshold
  artifact (the high-id pair is actually equal at 13.67 u).
- **Repeat too long to span: 0 families decisively** (one array in RCF_102, ~13 kb, exceeds a 10 kb
  read, but its 4.5 kb unique flank is read-spannable anyway, so spanning is never the binding
  constraint).

### CDS vs UTR breakdown of the genuine ones
**N/A — there are no genuine cases.** For completeness, of the (refuted) families the variable
repeat was: CDS in 5 (RCF_432, RCF_277, RCF_53, RCF_718, RCF_611), UTR/VNTR in 2 (RCF_102, RCF_13),
mixed (constant copy in CDS, variable copies in UTR) in 2 (RCF_6, RCF_135), and not-a-real-array
in 1 (RCF_30). Both CDS-repeat and UTR-VNTR are read-detectable; this breakdown would matter only
had any case survived.

### Bottom line — concrete identifiability-frontier gain
**Zero families' worth.** A structure-aware ("structural PSV") VG that adds tandem-unit count as a
distinguishing variant would newly resolve **0** of these 10 currently-SNP-collapsed candidates,
because in every one either (a) the copies that differ in unit count are already ≥8–52% divergent in
their unique sequence (SNP/indel-PSVs in the unique flanks already separate them), or (b) the copies
that are truly sequence-collapsed have **identical** unit counts (structure adds no axis), or (c) they
are dispersed across chromosomes (no multimapping, so nothing to resolve). The tandem-repeat CNV is
frequently real, but it is **redundant** with — never the **sole** — distinguishing axis. The
identifiability frontier does not move.

### What a "structural PSV" is, and how copy-assignment would use it
A **structural PSV** is a copy-distinguishing variant whose state is a **tandem-unit count** (a cyclic
bubble's traversal multiplicity in the variation graph) rather than a base substitution/indel — i.e.
the allele a read carries at this site is "how many times it traversed the repeat unit." Copy-assignment
(#2) would read the count off a single FLNC molecule that spans the array and thread the read to the
copy whose unit count matches, exactly as it threads SNP/indel alleles to copies — adding one extra
column (a count, not a base) to the per-read allele profile used in the path-cover/threading.

### Recommendation
**Do NOT wire structural-PSV (tandem-unit count) into copy-assignment (#2) now.** Across all 10
hand-picked best candidates it resolves nothing the linear SNP/indel-PSV model doesn't already
resolve (the unit-count signal is redundant where present and absent where it would be needed), and
the near-identical-yet-CNV-differing cell of the dichotomy is empirically empty. Park it: revisit only
if a *co-located* family with a **read-spannable** array, **equal unique-region sequence (≥0.98 after
masking)**, **and** a genuine unit-count difference is ever found — none exists in this census.
