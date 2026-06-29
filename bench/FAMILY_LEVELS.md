# Multi-copy gene families at three SEPARATE levels (RNA / DNA / protein)

The thesis defines families at the **RNA level** (the contribution). DNA and protein are **independent
orthogonal verifiers**, deliberately kept SEPARATE — each level stands on its own and can confirm/extend the
others without contaminating the RNA definition. (Origin: the project began as an assembler improving multi-copy
families → the RNA-level family *definition* was requested → DNA/protein levels are the natural follow-ups.)

| Level | What it measures | Tool (standalone) | Reach on real GGO |
|---|---|---|---|
| **RNA** (the definition) | identifiability — *which loci a read cannot tell apart* (read-conflict graph + POA core) | `gw_family_catalog` (in-pipeline) | ~**81–87%** identity (links divergent-flank paralogs via conserved-exon read-ties) |
| **DNA** | segmental-duplication / shared-exon + copy-number | `bench/soto_family_validate.py` (Soto-2025) · SEDEF segdup map (cluster) | shared ≥200 bp / ≥90% block + famCN |
| **PROTEIN** | coding paralogy — shared ORF homology | `bench/protein_family_verify.py` (this) | down to ~**57–59%** (mmseqs fident, far past nucleotide) |

All three are **annotation-free** (de-novo loci / ORFs), so each is a true orthogonal check, not a re-statement
of the same evidence.

## The protein harness (`bench/protein_family_verify.py`)

Standalone, mirrors the in-pipeline protein tier (`batch_protein_edges`): each copy's **longest 6-frame ORF**
(≥40 aa, from the catalog's emitted SPLICED `copies.fa`) → one **mmseqs easy-search** all-vs-all (`-s 7.5`) →
protein-homology edges (`fident≥0.50`, `min(qcov,tcov)≥0.50`). Reports three things about a catalog:

1. **PRECISION** — what fraction of RNA families are a single connected component in the protein graph (their
   copies are mutual protein-paralogs) → orthogonal confirmation the family is a real coding-paralog group.
2. **REACH** — within-family protein-identity distribution: how divergent protein still confirms paralogy.
3. **EXTENSION** — cross-family protein homology: copies in different RNA families that are protein-paralogs =
   loosely-related members the nucleotide RNA definition split apart (the protein level's unique contribution).

### Result on the same-chrom catalog (`gw_off`, 81 families, 205 spliced ORFs)
- **PRECISION: 24/81** families are a single protein-paralog component. *Not "70% wrong"* — it is the orthogonal
  finding that read-conflict families (shared *sequence* = identifiability) and protein-paralog families (shared
  *coding* = paralogy) are **different objects**: ~30% are both; the rest are linked by non-coding / shared-exon
  multimapping, not coding paralogy. This independently confirms the project's standing claim that read-conflict
  measures *identifiability-relevance, not evolutionary paralogy* ([[project_family_detection_validation]]).
- **REACH: median 94%, min 59%** protein identity (10 within-family edges < 70%) — protein confirms paralogy
  **far below** the nucleotide ~81–87% ceiling.
- **EXTENSION: 15 family-pairs** are protein-paralogs (`fident≥0.50`), down to **57–59%** (GWFAM62~65 59%,
  GWFAM62~77 57%) — genuinely loosely-related paralogs that ONLY the protein level connects; plus several at
  100% (e.g. GWFAM36~37, 57~58) = near-identical paralogs read-conflict split because reads don't cross-map
  between their loci. These are the candidates a protein-level family definition would merge.

## ⭐ Authoritative cross-tab — every read-conflict family × {protein, DNA, SEDEF} (`family_levels_crosstab.py`)

The single FP/FN artifact: each `gw_off` family scored at all three levels at once — **protein** (ORF→mmseqs),
**DNA** (Soto pairwise asm20, ≥200 bp/≥90% block), **SEDEF** (its two copies are the two ends of a segmental-
duplication pair in the independent `final.bed` segdup map, which never saw the RNA catalog).

| protein | DNA | SEDEF | families | % |
|:---:|:---:|:---:|---:|---:|
| ✓ | ✓ | ✓ | 20 | 25% |
| ✓ | ✓ | ✗ | 4 | 5% |
| ✗ | ✓ | ✓ | 26 | 32% |
| ✗ | ✓ | ✗ | 6 | 7% |
| ✗ | ✗ | ✓ | 13 | 16% |
| ✗ | ✗ | ✗ | **12** | **15%** |

Marginals: **protein 30% · DNA 69% · SEDEF-pair 73% · SEDEF-≥2-copies-in-segdup 90%.**
**Confirmed by ANY level (a real duplication family): 69/81 = 85%. By ALL three: 25%. By NONE: 12/81 = 15%.**

**What this resolves about FP/FN:**
- **The read-conflict catalog is ~85% confirmed as real duplications** by at least one orthogonal level; **SEDEF
  alone — the most independent check (DNA segdup map, built without the RNA catalog) — confirms 73%** as segdup
  pairs and 90% as sitting in duplicated regions. So the **false-positive ceiling is ~15%** (the 12 families no
  level confirms = shared non-coding / one exon / repeat — identifiability-only, not a duplication).
- **The protein-only 30% was MISLEADING.** The cross-tab shows why: 26/81 (32%) are DNA✓ **and** SEDEF✓ but
  protein✗ — they are **real segmental duplications that don't encode a conserved protein** (non-coding /
  pseudogene / UTR duplications). That is a **protein-level FN (protein is blind to non-coding copies)**, NOT a
  read-conflict FP. SEDEF/DNA, not protein, are the right arbiters of "is this a real duplication."
- **SEDEF > DNA-pairwise (73% vs 69%)** because SEDEF's Jaccard+chaining is more sensitive than end-to-end asm20
  — it recovers 13 families the pairwise check missed.

**Verdict:** as a *duplication-family* detector the read-conflict definition is ~85% precise (≤15% FP), not the
~70%-wrong the protein-only number implied. As an *identifiability* detector (its actual purpose) it is sound
(0 FP on labeled domain-sharers, FN-by-design on resolvable copies). The "errors" are mostly the LEVELS
measuring different things (protein blind to non-coding; nucleotide RNA blind to >87%-divergent), which is
exactly why the levels are kept separate.

`final.bed` = SEDEF run on the GGO genome (cluster, 253,030 segdup pairs), at `~/Desktop/final.bed`.

## DNA copy-number axis (famCN / parCN): scope and why deferred

SEDEF, famCN, and parCN answer **different** questions — they are not substitutes:

| | answers | status |
|---|---|---|
| **SEDEF** | *where* are the duplications ("are these two loci a segmental duplication of each other?") — structural map from the assembly | ✅ used — the family-reality check (73% of RNA families confirmed; closed the defense P0) |
| **famCN** | *how many* copies of a family exist genome-wide (read-depth over the family) | not run — but covered: SEDEF gives a structural copy count, and the cross-level famCN observation (RNA captured 2 of a 9–11-copy DNA family) already gives the expression-restriction story |
| **parCN** | *how many of each specific paralog, per individual* (paralog-distinguishing k-mers across a cohort) | not run — it is the **copy-vs-allele** resolver |

**Was SEDEF enough? For this thesis's claim, yes.** The thesis claims an RNA-level family *definition* +
copy-*assignment*; the orthogonal validation that needs is "are the families real paralog groups?", which SEDEF
answers independently. famCN would add a count, not a confirmation, and is already covered qualitatively.

**parCN is the one genuinely-additional thing, deferred for a principled reason:** it is the DNA paralog-CN axis
(Soto / QuicK-mer2), not the RNA contribution, and it is precisely what would *resolve* the copy-vs-allele
boundary this work flags rather than guesses. It needs **DNA WGS + a population/cohort** (Soto used SGDP + 4
archaics), which is a separate data acquisition from one gorilla's IsoSeq + the T2T assembly. Clean positioning:
*RNA assignment answers which copies are transcribed and how they splice; parCN answers how many exist in DNA per
individual — complementary, cited to Soto/QuicK-mer2, required only to close copy-vs-allele.* So parCN is named
**future work / DNA-side collaboration**, not a gap in the RNA result.

## Why separate
Keeping the levels independent means: (a) the RNA definition stays the clean thesis object (no protein/DNA
dependency baked in); (b) each level is a falsifiable orthogonal check of the others; (c) when asked for the DNA
or protein definition, each already exists as its own harness. Run any of them on any catalog prefix.

Reproduce: `MINIFORGE python bench/protein_family_verify.py <catalog_prefix>` (needs `mmseqs`, biopython).
