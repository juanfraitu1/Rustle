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

## Why separate
Keeping the levels independent means: (a) the RNA definition stays the clean thesis object (no protein/DNA
dependency baked in); (b) each level is a falsifiable orthogonal check of the others; (c) when asked for the DNA
or protein definition, each already exists as its own harness. Run any of them on any catalog prefix.

Reproduce: `MINIFORGE python bench/protein_family_verify.py <catalog_prefix>` (needs `mmseqs`, biopython).
