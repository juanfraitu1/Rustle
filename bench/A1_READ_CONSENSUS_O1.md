# A1 — Genuinely read-derived (non-reference-shared) corroboration of O1

**Question.** O1 is the multi-copy family / copy structure. Two prior corroborations both
collapsed onto a **shared reference substrate**:

- **L1** (`GENOME_RNA_OVERLAY.md`): the de-novo "transcripts" are reference substrings, so grouping
  them against SEDEF was reference-vs-reference.
- **A2** (`genome_rna_overlay_readcontent.py`): the PSV columns came from a minimap2 **asm20**
  self-alignment of the reference copies. The break-test (`genome_rna_overlay_breaktest.py`, 92 min)
  showed the **read gate was INERT** (ungated ρ −0.420 ≈ gated −0.443) → the PSV axis was
  asm20-vs-SEDEF concordance carrying **zero read information**.

A1 asks whether the copy/PSV structure can be derived from **read heterogeneity itself**
(SDA / Vollger 2019 correlation-clustering) — not from a reference-copy alignment — and whether that
read-derived structure genuinely escapes the reference substrate that sank L1 and A2.

Artifacts: `bench/a1_read_consensus_o1.py` (+ `a1_read_sda_smoketest.py` core), outputs
`bench/a1_read_consensus_o1.tsv` / `.json` / `.rows.json`. Run:
`PYTHONHASHSEED=0 python bench/a1_read_consensus_o1.py`.

---

## 1. Method — read-heterogeneity PSV calling (SDA, no asm20)

At a **collapsed** reference locus (one reference locus, MAPQ-0 reads from *K* paralogs piling up),
columns are **discovered from read disagreement**, never seeded by an asm20 alignment:

1. **Allele matrix** — `pysam` pileup of `GGO_mm.bam` (−N50 −p0.1), `read × pos → base`,
   base-qual ≥ 20, ACGT only, **secondary alignments kept** (the MAPQ-0 collapse mixture).
2. **Het columns** — depth ≥ 12, minor allele ≥ 4 reads **and** VAF ≥ 0.20. The reference base is
   **never compared**; a "het" column requires the *reads to disagree among themselves*.
3. **SDA linkage** — binarize (major = 0 / minor = 1), pairwise **φ** over co-covering reads, keep
   columns co-segregating |φ| ≥ 0.60 with ≥ 2 others (errors correlate with nothing → dropped).
   The kept columns are the **read-discovered PSV set**.
4. **Read partition → `K_read`** — union-find on reads agreeing ≥ 80 % over the linked columns;
   `K_read` = number of clusters with ≥ 4 reads.

No FASTA, no minimap2/mappy, no asm20 `copy_cols` enter column discovery (verified: `asm20`/`mappy`
strings in these two files occur only in docstrings; the imported SEDEF module has no module-level
aligner). asm20 enters **only** as the baseline A1 must beat (Control 1).

**A→I RNA-editing veto** (implemented here, previously absent). Canonical A>G (sense) / T>C
(antisense) transition columns are dropped and `K_read` is re-derived on the remainder
(`K_read_noedit`). This is a **conservative over-veto** — it also removes genuine A/G, T/C transition
SNPs (the commonest SNP class) — so a survivor that keeps `K_read ≥ 2` with **zero** transition
columns has a copy structure that cannot be an editing artifact.

---

## 2. Decisive independence control — does A1 escape the reference (unlike A2)?

**Control 1 (read-only vs asm20-seeded).** The *airtight* subset = the **27** `no_psv` families where
the asm20 self-alignment supplies **0 PSV columns** (`K_asm20 == 1` by construction, verified
0 exceptions in `psv_graph_genomewide.json`). Any `K_read ≥ 2` there is copy structure the reference
self-alignment **cannot** produce — het requires reads to disagree among themselves, which a single
reference locus cannot supply.

| funnel on the 27 airtight `no_psv` families | n |
|---|---|
| asm20 supplies 0 PSVs (`K_asm20 == 1`) | **27** |
| `K_read ≥ 2` | 13 |
| … survives the cross-read shuffle null | **11** |
| … + external SEDEF segdup partner (genome, never sees a read) | 10 |
| … + copy-like MAF (`balanced_frac < 0.35`, not diploid) | **8** |
| … + editing-robust (`K_read ≥ 2` after dropping all A>G/T>C columns) | **8** |
| het-risk survivors **excluded** (MAGEA9 bal 0.74; LOC115935025 bal 0.50) | 2 |

**Control 2 (shuffle null) — the A1 analog of the A2 break-test.** Permute alleles across reads within
each column (preserves per-column MAF, destroys cross-read co-segregation) and re-run steps 3–4.
Frontier linked columns collapse **1829 → 262 (14.3 % retained)** and `K_read → 1` in the null. This
is the **opposite** of A2: there the read gate was inert (ρ −0.420 ≈ −0.443); here the read linkage is
**load-bearing** — remove it and the copy structure vanishes. The 2 frontier families whose null
*itself* yields `K_read ≥ 2` (CPLANE1 fam75 `null_K = 4`; fam11 `null_K = 2`) are **per-family
excluded** by `survives_null`, so the surviving set is null-clean (each `null_K ≤ 1`).

**Editing veto result.** All **8/8** clean copy-like survivors remain `K_read ≥ 2` after every
canonical A>G/T>C column is removed (`frac_edit_psv` ≈ 0.60–0.74, i.e. most PSVs *are* transitions,
yet enough transversion structure remains to reconstruct ≥ 2 copies). The one survivor whose structure
collapses under the veto (LOC115935025 fam145, `K_read_noedit = 0`) is *also* the diploid-balanced,
no-external-partner case — triple-excluded, consistently.

**Verdict on independence: YES, A1 escapes the reference substrate** on this bounded subset — the A2
failure mode (asm20-seeded columns, inert read gate) is absent, and the shuffle null proves the read
linkage is the load-bearing signal.

---

## 3. Cross-modal corroboration vs SEDEF (with null) — MARGINAL, not established

Leg 2 asks the *separate, weaker* question: does the read-derived copy **count** `K_read` agree with
the genome copy count above a distance-matched permutation null, on the 21-family collapsed core
(`collapsed_excess > 0`)? Read-heterogeneity (SDA) vs genome self-alignment (SEDEF) are different data
+ method, so agreement would be genuinely cross-modal — **but it does not robustly clear the null**:

| K_read vs genome axis | Spearman ρ | perm-p (B = 5000) |
|---|---|---|
| SEDEF distinct partner regions | 0.443 | **0.052** (boundary) |
| n_loci | 0.378 | 0.092 (fail) |
| collapsed_excess + 1 | 0.259 | 0.250 (fail) |

Worse, `K_read` is largely a **read-depth proxy**: ρ(K_read, n_read_psv) = 0.842,
ρ(K_read, n_reads) = 0.773. Controlling `med_depth`, the partial correlation
ρ(K_read, SEDEF | med_depth) = **0.40, p = 0.08** — not significant. Half the core (7/21 with
`K_read = 0`) contributes no read structure, so the effect rests on ≈ 14 informative points at n = 21.

**Leg 2 is reported for transparency only. The A1 headline does NOT rest on it** — it clears the null
on at most 1 of 3 axes and only at the boundary, and it is depth-confounded. Language downgraded from
"corroborates SEDEF" to **"consistent with SEDEF partner *presence* at the null boundary, not
established as above-null copy-count agreement."** The strong, independent corroboration is the
**qualitative presence of an external SEDEF partner** on 10/11 survivors (Control 1), not the
quantitative `K_read`↔count correlation.

---

## 4. Honest scope — where A1 is genuinely read-derived

| regime | families | A1 status |
|---|---|---|
| collapsed **frontier** (`no_psv` + `partial`) | 30 (**19.5 %** of 154) | genuine read-derived independence lives here |
| collapsed **core** (`collapsed_excess > 0`) | 21 | Leg-2 comparison set (~14 informative) |
| airtight `no_psv` (asm20 = 0 PSVs) | 27 | Control-1 foothold |
| **defensible clean copy set** | **8–10** (**~5–6.5 %** of 154) | read-derived, null-surviving, editing-robust, SEDEF-partnered |
| **resolved majority** (`fully_resolvable`) | 124 | **reference-shared** — method-validation only |

On a reference-**resolved** T2T, uniquely-mappable reads match the reference, so read-derived ≈
reference-derived there — read consensus *is* the reference. A1's genuine independence is therefore
**concentrated in the collapsed subset**. This is disclosed self-critically: `K_read ≥ 2` actually
fires **more** on the resolved set (frac 0.625) than on the frontier (0.533), so `K_read ≥ 2` alone is
**non-diagnostic** — independence rests specifically on **asm20 = 0**, not on `K_read ≥ 2` being rare.

### Residual reference dependence (disclosed; categorically weaker than A2)

A1 is **not** reference-free, but its residual dependence is at the coordinate/recruitment level, not
the signal level:

- **(a) Read recruitment** is reference-determined — which reads pile up at a collapsed locus is set by
  minimap2 mapping to the T2T reference (MAPQ-0 multimappers). The **set** of reads is
  reference-recruited; the **partition** into copies is read-derived.
- **(b) Column coordinates** are `col.reference_pos` — positions are reference-aligned; the alleles,
  het test and cross-read linkage are read-only.

Both fix **which** reads and **where**, **not** the copy-distinguishing signal (het + linkage), which
the single reference locus (`K_asm20 = 1`) cannot produce. This is a coordinate-frame + recruitment
dependence — **categorically weaker** than A2's signal-level asm20 seeding of the PSV columns.

- **(c) Copy-vs-allele call** still borrows the genome: SDA read-heterogeneity alone cannot separate
  paralog *copies* from diploid *haplotypes*. The copy interpretation of the survivors uses the
  external SEDEF partner / `collapsed_excess`, and the two diploid-balanced survivors (MAGEA9,
  LOC115935025) are excluded via the `balanced_frac`/`het_risk` guard.

---

## 5. Bottom-line verdict

**A1 PARTIALLY CLOSES the O1 circularity risk — for the *independence* claim, on the collapsed
frontier only.**

- **Closed (solid):** read-derived copy structure that is genuinely **independent of the reference**
  exists on **8–10 families** (of the 27 airtight `no_psv`, where asm20 provably supplies 0 PSVs).
  It survives the cross-read shuffle null (frontier linkage collapses to 14.3 %; the read gate is
  **load-bearing**, the exact inverse of A2's inert gate), is editing-robust (8/8 hold after removing
  all transition columns), and 10/11 carry an external SEDEF segdup partner. This is a real escape from
  the shared-reference substrate that sank L1 and A2.
- **Not closed:** the *quantitative* corroboration that `K_read` equals the genome copy **count** above
  a distance-matched null is **not established** (best axis perm-p = 0.052; fails 2/3 axes;
  depth-confounded, partial ρ = 0.40, p = 0.08).
- **Scope:** genuine independence is confined to ≈ 19.5 % of families (the collapsed frontier) and the
  defensible clean set is ≈ 5–6.5 % (8–10 / 154). On the reference-resolved majority (124 families,
  81 %) read ≈ reference and A1 provides only method-validation, no new independence.
- **Residual:** read recruitment (minimap2 → T2T) and column coordinates (`reference_pos`) remain
  reference-tied, and the copy-vs-allele call borrows the genome — all disclosed, none of which seeds
  the read-only het + linkage signal.

**In one line:** A1 is *not* reference-seeded and its read gate is *not* inert (unlike A2); it
delivers genuinely read-derived, reference-independent copy structure on a bounded ~6 % collapsed
subset, closing the O1 circularity risk **there** — but it does **not** establish above-null agreement
with the genome copy count, so the closure is **partial and scope-limited, not whole-catalog**.

---

### Determinism & reproducibility
`PYTHONHASHSEED=0` asserted; per-family `default_rng(1000 + family)`, perm-null `default_rng(20260701)`;
rows sorted before writing. Two independent full runs produced **byte-identical** per-family pileup
values; the `A1_ROWS_CACHE=1` aggregation path reproduces the TSV/JSON byte-identically. Fixed
thresholds: VAF ≥ 0.20, φ ≥ 0.60, `balanced_frac` copy/het cutoffs 0.35 / 0.40, null-retain 0.30.
