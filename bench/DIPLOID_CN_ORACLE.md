# Phased-diploid copy-number oracle (RNA-independent corroboration)

Use the **phased diploid mGorGor1 assembly** — maternal (chrX) + paternal (chrY)
haplotypes, built from DNA WGS (HiFi/ONT), **fully independent of the IsoSeq RNA
reads** — as a DNA-derived copy-number oracle to (a) give the non-circular
corroboration A1/SEDEF could not, (b) resolve copy-vs-allele, and (c) sharpen the
RNA multi-copy catalog.

- Deliverables: `bench/diploid_cn_oracle.py` (counter, TSV md5 `9300e681…`),
  `bench/diploid_cn_oracle.tsv` (196 families), `bench/diploid_cn_corroborate.py`
  → `bench/diploid_cn_corroborate.json` (md5 `6ee643d4…`),
  `bench/diploid_cn_oracle_sensitivity.py` (threshold robustness).
- Inputs: `validated_families.tsv` (196 co-located families = per-copy reference
  loci), `GGO.fasta` (GCF_029281585.2), `a1_read_consensus_o1.tsv` (K_read/χ(H)/
  sedef), PAFs `cn_oracle_paf/{mat,pat}.paf`.
- Determinism: `minimap2 -t` fixed, `PYTHONHASHSEED=0`, set-based counts, sorted
  output; two independent runs are **byte-identical**.

## Method

For each family: extract every reference copy locus from `GGO.fasta`, align all
loci to **each** haplotype (`minimap2 -cx asm20 -N 50 -p 0.1 --secondary=yes`),
keep alignments at **identity ≥ 0.90**, then count **distinct full-length gene
copies** per haplotype.

The counter is the corrected core (the previous single-`merge_gap` /
absolute-1 kb-floor counter was biased in two directions; both are fixed):

1. **Length-proportional copy criterion (kills domain inflation).** A target
   region counts as a copy only if the **aggregate collinear query coverage** of
   that region reaches **≥ 0.5 · L**, where `L = median reference-locus length` of
   the family. A shared sub-gene domain/repeat (e.g. a 2 kb module) no longer
   counts as a full copy of a multi-kb gene. Sub-0.5·L regions are tallied
   separately as `frag_mat/frag_pat` — the inflation flag is folded **into** the
   shipped count, not left in a side file.
2. **Tandem-aware copy segmentation (kills tandem collapse).** Alignments of one
   query that **overlap in query space** are kept as **separate** copies (adjacent
   tandem-array units each re-align the whole gene → distinct copies), while
   query-**disjoint** collinear blocks are stitched into one copy (a single
   alignment split by an indel). No single `merge_gap` is used, so a ~1 kb
   inter-unit spacing no longer collapses a tandem array to CN = 1.
3. **Sex-aware.** mGorGor1 is **male**: maternal carries chrX, paternal chrY. An
   X-gene reads `hap_CN_pat = 0`, a Y-gene `hap_CN_mat = 0` — resolved as male
   hemizygosity, **not** an inter-haplotype CNV.

`asm_hapCN` = sex-aware haploid CN (X→mat, Y→pat, autosome→max(mat,pat)).
Robustness: sweeping the 0.5·L threshold over {0.4, 0.5, 0.6}·L leaves every
headline family stable (tandem fam42 = 12, over-splits, het_risk, MAGEA9,
LOC115935025 all invariant); only heavily domain-laden families shift, and those
carry large `frag` flags regardless.

## Assembly haploid copy-number distribution (196 families)

| asm_hapCN | 1 | 2 | 3–4 | 5–9 | 10–19 | 20+ |
|-----------|---|---|-----|-----|-------|-----|
| families  | 26 | 73 | 49 | 23 | 18 | 7 |

## (1) Non-circular corroboration (A1 subset, n = 62; Spearman ρ + 10 k-permutation null)

| DNA asm_hapCN vs | ρ | perm_p | reading |
|---|---|---|---|
| **K_read** (only truly RNA-independent axis) | **+0.313** | **0.0124** | positive, above null |
| n_loci_ref (A1 subset) | +0.456 | 0.0005 | clearly above null |
| **χ(H)** | **+0.278** | **0.0296** | above null |
| sedef_partners | +0.615 | 0.0001 | strong (duplication-content proxy) |
| n_loci_ref (all 196) | +0.279 | 0.0003 | above null |

`asm_hapCN ≥ K_read` in **56/62** families — K_read is mostly a **censored lower
bound** (the K = 0 read-identifiability floor), so the DNA oracle **bounds the RNA
resolution ceiling** rather than tracking the read count tightly. The 6 exceptions
are the read method *over*-counting relative to DNA, and 2 of them (DHRSX fam58,
LOC129530050 fam91) are exactly the het_risk allele families where a heterozygous
site was read as a 2nd copy — which the DNA correctly calls single-locus. This is
the corroboration **SEDEF structurally cannot give**: SEDEF is reference
self-alignment (co-varies with n_loci by construction), whereas the phased diploid
is an independent DNA observation. (Fixing the tandem-collapse bug *strengthened*
every correlation — e.g. K_read rose from ρ ≈ 0.24/p ≈ 0.06 to ρ = 0.313/p = 0.012
— because fam42's true CN = 12 now matches its higher K_read.)

## (2) Copy-vs-allele resolution

Class distribution (sex-aware): **169 multi_copy / 26 single_locus_allele / 1 cnv**.
Sex-linked: 14 X + 2 Y families resolved as male hemizygous (present-haplotype CN),
**not** mislabelled CNV. Autosomal CNV is now surfaced explicitly: **1** family with
one haplotype < 2 copies (fam40, 1/2 = `cnv`), plus **27** families flagged
`cnv_autosomal=1` (both haplotypes ≥ 2 but the counts differ, e.g. ZNF425 21/26) —
these were previously hidden inside `multi_copy` by a precedence bug.

- **A1 candidate families DNA-confirmed real multi_copy: 60/62.** The **only 2 not
  confirmed** are **exactly the two RNA `het_risk = 1` flags** — **DHRSX (fam58)**
  and **LOC129530050 (fam91)** — which the DNA independently calls single-locus
  (1/1). Their RNA variation is allele/heterozygosity, not copies.
- **26** families the RNA/reference method listed as multi-locus collapse to **1
  full copy per haplotype** in DNA ⇒ reclassified **copy → allele/het** (includes
  reference over-splits such as the 32-locus fam2 → 1).
- **A1's two diploid-het-risk survivors both resolve as REAL multi_copy:**
  - **MAGEA9 (fam94):** X-linked, mat = 2 / pat = 0. Two copies on a single male X
    cannot be alleles of one locus ⇒ necessarily **paralogous** ⇒ multi_copy.
  - **LOC115935025 (fam145):** autosomal, mat = 2 / pat = 4 (both ≥ 2, also
    `cnv_autosomal`) ⇒ **multi_copy**.

## (3) Reference-collapsed copies the phased genome recovers, and over-splits

- **16 genuine reference-collapsed full-gene copies** (`asm_hapCN ≥ 2·n_loci_ref`,
  Δ ≥ 3), all full-length (fragments already excluded): e.g. **fam49 7 → 30**,
  fam188 2 → 25, **ZNF425 5 → 26**, LOC129531752 (fam139) 2 → 22, fam44 7 → 21, and
  the recovered tandem array **LOC129529768 (fam42) 1 → 12** that a single merge-gap
  had collapsed. Non-circular support for the O4 "more copies than reference" thread.
- **46 reference over-splits** (`n_loci_ref ≥ 2·asm_hapCN`): the reference lists many
  co-located loci that map to a few true DNA copies — e.g. fam2 32 → 1, GSTM2 47 → 19,
  fam3 31 → 3, TUBGCP5 (fam5) 23 → 2. The diploid gives the lower, true CN.
- **25 domain-heavy families** (≥ 5 distinct sub-gene domain/repeat regions in
  `frag_*`) whose shared modules are now **rejected** from the copy count rather than
  inflating it — e.g. fam76 (22 raw → 1 copy + 21 frag), fam80 (37 raw → 4 + 33 frag),
  LOC115930164 (fam14, 37 raw → 13 + 32 frag).

## Verdict

The phased diploid mGorGor1 genome delivers a genuinely RNA-independent
copy-number oracle. **(a) Non-circular corroboration — YES, honest/partial:** it
independently confirms the A1 families are real paralogy and its copy count is
positively, above-null rank-correlated with the truly read-derived K_read
(ρ = 0.313, perm_p = 0.012) and with χ(H) (ρ = 0.278, perm_p = 0.030), while
sitting ≥ K_read in 56/62 families as a censored-ceiling bound — the check SEDEF
(reference self-alignment) structurally could not provide. **(b) Copy-vs-allele —
YES:** 60/62 A1 families DNA-confirmed multi_copy, the only 2 unconfirmed are
exactly the RNA het_risk flags (DHRSX, LOC129530050) called single-locus by DNA,
both het-risk survivors (MAGEA9, LOC115935025) resolve as real copies, 26 families
reclassified copy→allele, and autosomal CNV is now surfaced (1 class-cnv +
27 flags) instead of hidden. **(c) Catalog sharpened — YES, both directions:**
16 reference-collapsed real copies recovered (incl. the fam42 tandem array a
buggy merge-gap had erased), 46 reference over-splits corrected, and 25
domain-heavy families de-inflated with shared modules rejected from the count.

**Caveat (species-fixed vs donor-exact):** copy number / paralogy is largely
**species-fixed**, so the diploid assembly is a valid RNA-independent copy-number
oracle **regardless** of whether the IsoSeq donor is mGorGor1. Individual-exact
het/allele genotyping of a *specific* variant is donor-exact only if RNA ==
mGorGor1 (unknown). Lead with the robust copy-number / copy-vs-allele use, not
variant-level genotyping.
