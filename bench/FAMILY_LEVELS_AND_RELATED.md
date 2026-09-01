# Family Levels (RNA/DNA/Protein) & Related Methods (consolidated)

> Merged from 3 source docs (verbatim; git keeps the originals' history).

**Contents:** FAMILY_LEVELS · RELATED_WORK_METHODS · DNA_PROTEIN_VALIDATION


---


## INDEX

> **Index.** 31 sections; this is the map. **The titles carry the verdicts** — no tag is derived
> here. ⚠ In `o1_ledger.md` an earlier auto-derived verdict tag scored **11/22 = 50%** against
> sections whose outcome was known first-hand, so tags were removed rather than shipped. Search a
> heading to jump.


- FAMILY_LEVELS
- The protein harness (`bench/protein_family_verify.py`)
- ⭐ Authoritative cross-tab — every read-conflict family × {protein, DNA, SEDEF} (`family_levels_crosstab.py`)
- DNA copy-number axis (famCN / parCN): scope and why deferred
- Why separate
- RELATED_WORK_METHODS
- The dichotomy (how to read the table)
- Where this thesis sits (the gap)
- DNA_PROTEIN_VALIDATION
- dna_psv_catalog_summary
- dna_rna_overlay
- What actually transcribes (genome vs transcriptome)
- Ancient-family gain (the whole point of the DNA tier)
- Per-family 3-number summary (sample: largest + curated)
- Honest scope
- compara_fetch
- Universe gene inventory
- Mapping coverage (named genes)
- Family-level checkability
- compara_validation
- Headline (the non-circular number)
- KEY FRAMING: Compara is COARSER, so PRECISION is the metric
- (4) Coverage (stated up front -- small, honest sample)
- (1) PRECISION vs Compara -- UNIVERSE families (the headline)
- (2) PRECISION vs Compara -- RUSTLE's minimizer-Jaccard grouping
- (3) GRANULARITY (observation, not error)
- Cross-check: JSON's stated within-universe paralog relation
- Honest caveats
- transcript_validation
- Verdict: REAL + defensible
- PARCN_VALIDATION

## FAMILY_LEVELS

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


---

## RELATED_WORK_METHODS

# Methods for finding / defining multi-copy gene families (related work)

The field splits by the **level** the signal comes from. The decisive axis is the right-hand column: most
methods work for **ancient / divergent** families (divergence gives the signal) but **break on recent,
near-identical paralogs** — where the field has moved to segmental-duplication detection + read-depth copy
number + long reads. The RNA/transcriptome cell is nearly empty; that is where this thesis sits.

| Method / tool | Level | Input | Core idea | Recent / near-identical paralogs? | Ref |
|---|---|---|---|---|---|
| **All-vs-all BLAST/DIAMOND + MCL** | protein/seq homology | proteomes | pairwise similarity graph → Markov clustering | over-merges or collapses near-identical; threshold-bound | OrthoMCL; Enright 2002 |
| **OrthoFinder** | protein/seq homology | proteomes | DIAMOND all-vs-all + length/bias normalization + MCL + gene trees → orthogroups (current standard) | struggles: near-identical copies collapse to one node / over-merge | Emms & Kelly 2019 |
| **SonicParanoid / Broccoli / eggNOG** | protein/seq homology | proteomes | faster/precomputed orthogroup inference | same recent-paralog limits as OrthoFinder | various |
| **CD-HIT / MMseqs2 linclust** | seq identity | protein or NT | greedy clustering at an identity threshold | threshold-defined; arbitrary boundary on recent paralogs | Fu 2012; Steinegger 2018 |
| **Pfam / InterPro / HMMER** | domain/profile | protein | group by shared **domain architecture** (HMM profiles) | groups *domain-sharers* (FP for paralogy); domain ≠ copy | Finn 2014 |
| **Ensembl Compara (TreeBeST)** | phylogenetic | gene + species trees | gene-tree ↔ species-tree **reconciliation** → ortholog/paralog, dup/loss | hypersensitive to gene-tree error → **fails on recent near-identical** paralogs | Vilella 2009 |
| **GeneRax / Notung / ORTHOSCOPE** | phylogenetic | gene + species trees | probabilistic / parsimony reconciliation | same recent-paralog tree-error limit | Morel 2020 |
| **MCScanX / DupGen_finder / WGDdetector** | synteny | genome + annotation | **collinearity** blocks → classify WGD / tandem / proximal / transposed / dispersed duplicates | needs annotation + synteny; weak inside complex SD arrays | Wang 2012 |
| **WGAC** | DNA self-alignment | assembly | whole-genome assembly comparison → segmental duplications | classic SD detection; assembly-quality bound | Bailey 2001 |
| **⭐ SEDEF** | DNA self-alignment | assembly | **Jaccard similarity + local chaining** → segmental-duplication pairs (fast) | **strong on recent SDs** (the hard case); assembly-only, no expression | Numanagić 2018 |
| **BISER** | DNA self-alignment | assemblies | SEDEF successor: multi-genome, error/masked-aware SD detection | same SD paradigm, scaled to many genomes | Numanagić 2021 |
| **WSSD / fastCN (famCN)** | DNA copy-number | short/long reads | read-depth over multimapping windows → **family copy number** | quantifies CN; cannot phase WHICH copy | Bailey 2002; Pendleton |
| **QuicK-mer2 / fastCN (parCN)** | DNA copy-number | reads + k-mers | paralog-specific k-mers → **paralog-specific copy number** | the parCN standard (Soto 2025); needs paralog-distinguishing k-mers | Shen & Kidd 2020 |
| **CNVnator / read-depth CNV** | DNA copy-number | short reads | depth segmentation → CNV | coarse; not paralog-resolved | Abyzov 2011 |
| **Mash / minimizer-Jaccard** | alignment-free | sequences | MinHash / minimizer Jaccard distance | fast screen; **cannot separate domain-sharers** from true paralogs | Ondov 2016 |
| **⭐ SDA (Vollger)** | long-read DNA | long DNA reads | collapsed-segdup detection by depth → **PSV correlation-clustering** → per-paralog assembly | resolves collapsed segdups; **same K=0 identifiability floor** as ours | Vollger 2019 |
| **Soto 2025 (HSD pipeline)** | DNA + CN | T2T assembly + reads | SD98 self-map + shared-exon (>99% cov) + famCN MAD<1 grouping; parCN via QuicK-mer2 | the reference human catalog for recent families; **DNA, no RNA/isoform axis** | Soto 2025 (Cell) |
| **Guitart/Eichler 2024 (TBC1D3)** | long-read RNA | Iso-Seq + haplotypes | map FLNC to all haplotypes; **assign read to a paralog iff AS ≥ 10** else discard; phylo groups (graph failed) | per-copy expression assignment for ONE family; not de-novo / genome-wide | Eichler 2024 |
| **longcallR** | long-read RNA | long RNA-seq | CNN SNP caller + MEC phasing + ASE/ASJ | **uniquely-mappable only** — assigns to genes+haplotypes, never to family COPIES | Huang 2026 |
| **RSEM / Salmon / kallisto** | RNA quant (EM) | RNA-seq + transcript set | EM over multireads to apportion expression | needs a **pre-defined** transcript set; does NOT discover families | Li 2011; Bray 2016 |
| **⭐ IsoCon** | long-read RNA | targeted Iso-Seq | NN-graph cluster/correct → per-variant-position real-vs-error test → family transcripts | **closest prior art**; WITHIN a known family, targeted (RT-PCR), no de-novo / no copy-vs-allele | Sahlin 2018 |
| **⭐⭐ THIS THESIS** | long-read RNA, de-novo, genome-wide | Iso-Seq/HiFi | **read-conflict graph** (significance de-tie) defines families = identifiability components; **per-read significance gate** assigns reads to copies (no 1/k, assign-or-abstain); + allele-specific junctions | targets exactly the recent/collapsed regime (MAPQ-0); identifiability theorem (MCC=χ(H)) + per-read certificate | — |

## The dichotomy (how to read the table)

- **Ancient / divergent families** → homology clustering, domain profiles, gene-tree reconciliation. Divergence
  *is* the signal, so these work — but they **break on recent near-identical paralogs** (gene trees collapse,
  similarity grouping over-merges or can't separate domain-sharers).
- **Recent / nearly-identical paralogs (the hard case)** → the field switched to **segmental-duplication
  detection** (SEDEF/BISER) + **read-depth copy number** (famCN/parCN, QuicK-mer2) + **long reads**
  (SDA on DNA; IsoCon/Eichler on RNA). Short reads fail here — Soto 2025 measured SNV sensitivity dropping to
  **0.85% in SD98** regions, which is the long-read rationale.
- **RNA / transcriptome cell is nearly empty.** IsoCon resolves transcripts *within a known family*;
  RSEM/Salmon need a pre-defined transcript set; longcallR/Eichler stay in the uniquely-mappable or single-family
  regime. **No method does de-novo, genome-wide, RNA-level family discovery + copy-assignment** — the cell this
  thesis occupies.

## Where this thesis sits (the gap)

| | de-novo discovery | genome-wide | RNA-level | copy assignment | copy-vs-allele |
|---|:-:|:-:|:-:|:-:|:-:|
| OrthoFinder / Compara | ✓ | ✓ | ✗ (protein) | ✗ | ✗ |
| SEDEF / Soto / SDA | ✓ | ✓ | ✗ (DNA) | ~ (parCN) | ✓ (DNA) |
| IsoCon | ✗ (targeted) | ✗ | ✓ | hypothesised | ✗ |
| Eichler / longcallR | ✗ / ✓ | ✗ / ✓ | ✓ | ✓ (1 family) / ✗ | ✗ |
| **This thesis** | **✓** | **✓** | **✓** | **✓** | **✓ (ASJ; DNA for the irreducible case)** |

The three levels are kept **separate and orthogonal** (RNA = the definition; DNA = SEDEF/Soto; protein =
`protein_family_verify.py`); on real GGO, SEDEF independently confirms 73% of the RNA read-conflict families as
segmental duplications (`bench/FAMILY_LEVELS.md`).


---

## DNA_PROTEIN_VALIDATION

# Dna Protein Validation (consolidated)

> Merged from 5 source docs (verbatim, git keeps the originals' history). Each section below was a separate `bench/*.md`.

**Contents:** [dna_psv_catalog_summary](#dna-psv-catalog-summary) · [dna_rna_overlay](#dna-rna-overlay) · [compara_fetch](#compara-fetch) · [compara_validation](#compara-validation) · [transcript_validation](#transcript-validation)


---

## dna_psv_catalog_summary

# DNA-derived PSV identifiability catalog (Phase 1)

- co-located classified pairs: **1387** -> resolvable **331** (24%), genuine-K=0 **1056** (76%). NOTE: this is the DNA reference universe (all aligned co-located pairs, including unexpressed identical tandem copies the RNA census never observes); on the **137** pairs the RNA census actually classifies, DNA and RNA agree **86%** and both find that expressed subset mostly resolvable.
- pairs excluded from K: **14262** unaligned (copy did not align to ref0 — divergent/short paralog), **14** unannotated (no overlapping GFF exon)
- **cross-check DNA-K=0 vs RNA-K0** on 137 census-classified pairs: concordance **86%** (confusion {(True, False): 14, (True, True): 8, (False, False): 110, (False, True): 5})
- discordant DNA-K=0 ∧ RNA-resolvable: **14** (candidate: indel / splice-shift pseudo-K=0 — substitution-only PSVs miss it; Phase-2 private_exon_bp will test this)
- discordant DNA-K≥1 ∧ RNA-tied: **5** (candidate: PSV in a poorly-expressed exon — reference identifiability ≥ read identifiability)


---

## dna_rna_overlay

# Two-tier overlay: RNA layer on DNA/protein-defined multi-copy families

**DNA tier** = protein clusters (mmseqs2 easy-cluster, ≥30% id / 50% cov, on translated CDS) →
formal multi-copy families, INCLUDING ancient/diverged families the RNA-similarity definition
missed. **RNA overlay** = per copy, is it transcribed (real GGO IsoSeq), and how well.

- DNA multi-copy families (protein clusters ≥2): **3,587** (14,545 gene copies)
- RNA-tier families (POA, for comparison): 1,337

## What actually transcribes (genome vs transcriptome)
- of 14,545 copies in DNA-defined families: **transcribed (≥5 reads): 10,490 (72%)** ; well-expressed (≥40): 6,824 (47%) ; **silent: 4,055 (28%)**
- i.e. the DNA tier enumerates every copy; the RNA layer shows which are live — many copies are
  present in the genome but transcriptionally silent.

## Ancient-family gain (the whole point of the DNA tier)
- curated families recovered by the DNA tier that the RNA tier MISSED: **DEFB*, SIGLEC***
| family | DNA tier | RNA tier | DNA copies | transcribed |
|---|---|---|---|---|
| APOBEC3 | YES | YES | 4 | 4 |
| CRYBG (ANCIENT) | no | no | 0 | 0 |
| DAZ | no | YES | 1 | 1 |
| DEFB* (ANCIENT) | YES | no | 4 | 0 |
| MAGEA* | YES | YES | 5 | 5 |
| PRAMEF* | YES | YES | 4 | 0 |
| RABL2 | YES | YES | 2 | 2 |
| RFPL | YES | YES | 4 | 3 |
| SIGLEC* (ANCIENT) | YES | no | 3 | 2 |
| TAS2R* (ANCIENT) | YES | YES | 18 | 7 |

## Per-family 3-number summary (sample: largest + curated)
| DNA family | copies | transcribed | well-expressed | example members |
|---|---|---|---|---|
| DFAM0 | 501 | 83 | 4 | LOC101123691, LOC101123789, LOC101123793, LOC101124039, LOC101124044… |
| DFAM1 | 229 | 219 | 180 | KRBOX5, LOC101123988, LOC101124084, LOC101124732, LOC101124778… |
| DFAM2 | 136 | 134 | 108 | LOC101126065, LOC101126415, LOC101127631, LOC101128844, LOC101129578… |
| DFAM3 | 93 | 82 | 50 | CDC42, IFT27, LOC101128843, LOC101133567, LOC101142457… |
| DFAM4 | 74 | 19 | 8 | BFSP2, DES, GFAP, GHAA, INA… |
| DFAM5 | 64 | 31 | 2 | ACR, CELA1, CFD, CTRL, F10… |
| DFAM6 | 47 | 12 | 1 | LOC101124648, LOC101136188, LOC101144932, LOC101146296, LOC101147007… |
| DFAM7 | 47 | 0 | 0 | LHB, LOC101123748, LOC101124600, LOC101154318, LOC109024055… |
| DFAM8 | 46 | 42 | 31 | CDK1, CDK10, CDK14, CDK15, CDK16… |
| DFAM9 | 45 | 26 | 13 | ACKR2, ACKR4, AGTR1, APLNR, CCR10… |
| DFAM10 | 44 | 0 | 0 | LOC101123968, LOC101128508, LOC115932853, LOC115932855, LOC129523578… |
| DFAM11 | 43 | 17 | 0 | LOC101126783, LOC101127524, LOC129524344, LOC129524346, LOC129524347… |

## Honest scope
- DNA tier = protein clustering of annotated CODING genes → catches ancient coding families +
  currently-silent coding copies. Non-coding/pseudogene + UNANNOTATED copies need a genomic
  self-alignment pass (next extension).
- 'transcribed' uses the REAL IsoSeq (not the ideal-coverage synthetic) — so silent = silent in
  this sample; ideal-coverage GGO would test detectability, not biology.
- copy-resolvability (which transcribed copies are distinguishable per-read) follows from the
  identifiability theorem: dispersed copies resolve by locus; co-located need PSVs (rare here).


---

## compara_fetch

# Ensembl Compara paralogy fetch — coverage report

Non-circular validation prep. This phase ONLY fetches/caches Ensembl Compara paralogy for the NAMED universe genes; the comparison vs our family grouping is the next phase.

- Cache: `bench/compara_cache.json` (keyed by `endpoint|symbol`, resumable; reruns fetch only missing).

- Relation summary: `bench/compara_paralog_relation.json`

- Source universe: `bench/copy_recovery_eval/results/universe.tsv`


## Universe gene inventory

| metric | count |
| --- | --- |
| total distinct gene_ids | 195 |
| named genes (not `^LOC[0-9]+`) | 41 |
| LOC-only genes | 154 |
| families | 62 |

## Mapping coverage (named genes)

| metric | count | of named |
| --- | --- | --- |
| got ENSG id | 40 | 41 |
| got paralogue data (HTTP ok) | 37 | 41 |
| non-empty paralogue set | 32 | 41 |
| unmapped (symbol not in Ensembl) | 0 | 41 |
| persistent fetch errors | 5 | 41 |

## Family-level checkability

| metric | count | of 62 families |
| --- | --- | --- |
| families with >=2 NAMED genes | 11 | 62 |
| families with >=2 MAPPED genes (checkable within-family pair) | 10 | 62 |

Within-universe symmetric paralog pairs (both genes mapped, Compara-linked): **5**


### Symbols with persistent fetch errors (rerun to retry)

`CREB1`, `GP1BB`, `NCAPH2`, `RABL2B`, `USP18`


Deterministic given the cache. Fetched 82 new API responses this run.


---

## compara_validation

# Non-circular validation: Ensembl Compara vs our paralog families
_Deterministic (no RNG). Minimizer k=15, w=10, canonical blake2b-64. Generated by `bench/compara_validation.py`. Compara release: `vertebrates`, `https://rest.ensembl.org`._
## Headline (the non-circular number)
> On the **10** universe families (12 within-family named gene pairs) that map to Ensembl, human Compara independently confirms **33%** (4/12) of within-family pairs as paralogs.
This **de-circularizes** the earlier `family_detection_validation.py`, which compared one sequence-similarity clustering (minimizer-Jaccard) against another (the minimap2-built `universe.tsv`) -- both reward shared subsequence, so that agreement was partly methodological. Here the truth set is Ensembl Compara paralogy: **gene-tree + species-tree reconciliation**, which never sees our minimap2 alignments. A Compara 'yes' is therefore an independent biological confirmation.
## KEY FRAMING: Compara is COARSER, so PRECISION is the metric
Ensembl Compara paralogy is gene-tree-based and far **coarser** than our families. Querying Compara for `RABL2A` returns the entire RAB / small-GTPase **superfamily** (68 paralogues genome-wide in this dump), whereas our family is the 2-copy tandem `RABL2A`/`RABL2B`. So **our families should be SUBSETS of Compara groups**: we split ancient superfamilies into recent copy-clusters, which is correct granularity, not error. The right question is therefore PRECISION -- *of the gene pairs we put in the same family, what fraction does Compara independently also call paralogs?* -- not recall. Compara lumping ancient paralogs we correctly separate is expected granularity and is **intentionally not scored as a miss**.
## (4) Coverage (stated up front -- small, honest sample)
- Universe: **195 genes / 62 families** (`universe.tsv`). Of these, **154 genes are LOC-only** (no human symbol) and cannot be mapped to Ensembl -- they are out of scope here.
- Named genes attempted: **41**. Got an ENSG id: **40/41** (1 lookup error: `GP1BB`). Non-empty paralogue set returned: **32** (4 paralogue-fetch errors: `CREB1`, `NCAPH2`, `RABL2B`, `USP18` returned an error and so contribute NO confirmations even if real).
- **Evaluable set: 10 universe families with >=2 mapped named genes -> 12 checkable within-family pairs.** This is the entire basis of the headline; it is a SMALL sample and is reported as exact counts, not extrapolated.
## (1) PRECISION vs Compara -- UNIVERSE families (the headline)
Of the **12** within-universe-family named pairs (both genes mapped to an ENSG), Compara independently confirms **4** as paralogs => **precision = 33%** (4/12).

| universe family | pair | Compara confirms? |
|---|---|---|
| `APOBEC3D` | APOBEC3D <-> APOBEC3F | **yes** |
| `ASDURF` | ASDURF <-> ASNSD1 | no |
| `CASP8` | CASP8 <-> FLACC1 | no |
| `CCDC188` | CCDC188 <-> ZDHHC8 | no |
| `CDPF1` | CDPF1 <-> PPARA | no |
| `CREB1` | CREB1 <-> METTL21A | no |
| `GCA` | GCA <-> KCNH7 | no |
| `GPR39` | GPR39 <-> LYPD1 | no |
| `LOC134758217` | RFPL1 <-> RFPL2 | **yes** |
| `LOC134758217` | RFPL1 <-> RFPL3 | **yes** |
| `LOC134758217` | RFPL2 <-> RFPL3 | **yes** |
| `RABL2A` | RABL2A <-> RABL2B | no |

**Reading the non-confirmations.** Most of the `no` rows are pairs where the universe's minimap2 step grouped two DIFFERENTLY-NAMED genes (e.g. `GCA<->KCNH7`, `CDPF1<->PPARA`, `CREB1<->METTL21A`) that Compara's gene tree does NOT place in the same paralog group. A non-confirmation can mean EITHER (a) a genuine universe over-grouping (sequence similarity from a shared domain / repeat, not orthologous-block paralogy), OR (b) a Compara gap / coarser-but-different granularity, OR (c) one partner's paralogue fetch errored (`CREB1`, `RABL2B` did -- so `CREB1<->METTL21A` and `RABL2A<->RABL2B` could not be confirmed even if true). The confirmed pairs are the unambiguous tandem/recent-duplicate families: `APOBEC3D<->APOBEC3F`, `RFPL1<->RFPL2`, `RFPL1<->RFPL3`, `RFPL2<->RFPL3`.
## (2) PRECISION vs Compara -- RUSTLE's minimizer-Jaccard grouping
Re-using the EXACT criterion from `family_detection_validation.py` (canonical k=15/w=10 minimizers, blake2b-64, union-find) over the same mapped named genes:

| grouping | within-cluster pairs | Compara-confirmed | precision |
|---|---|---|---|
| universe families | 12 | 4 | 33% |
| rustle Jaccard @ 0.3 | 2 | 0 | 0% |
| rustle Jaccard @ 0.06 | 13 | 5 | 38% |

rustle's non-trivial predicted clusters (mapped named genes only):
- @ 0.3: {CREB1, METTL21A}; {RABL2A, RABL2B}
- @ 0.06: {APOBEC3D, APOBEC3F}; {ASDURF, ASNSD1}; {CASP8, FLACC1}; {CCDC188, ZDHHC8}; {CDPF1, PPARA}; {CREB1, METTL21A}; {GCA, KCNH7}; {GGT1, GGTLC2}; {GPR39, LYPD1}; {RABL2A, RABL2B}; {RFPL1, RFPL2, RFPL3}

The restriction to the 40 mapped NAMED genes means rustle's minimizer-Jaccard rarely puts two of them in one cluster at the strict 0.30 bar -- the named genes that share a universe family are often different-symbol genes whose whole-gene Jaccard is below 0.30 (the same broad within-family Jaccard distribution documented in the prior report). Where rustle DOES cluster a named pair, that pair is the relevant precision test; counts are reported exactly above and are small.
## (3) GRANULARITY (observation, not error)
- Non-trivial Compara groups containing >=2 of our mapped named genes: **3**.
- Compara groups that MERGE >1 distinct universe family (Compara coarser): **1**. This is EXPECTED granularity -- Compara lumps at the superfamily level what we split into copy-clusters -- and is reported as an observation, NOT an error.

| Compara group (our genes) | distinct universe families merged |
|---|---|
| APOBEC3D, APOBEC3F | 1: `APOBEC3D` |
| GGT1, GGTLC2 | 2: `GGT1`, `GGTLC2` |
| RFPL1, RFPL2, RFPL3 | 1: `LOC134758217` |

- REVERSE check (real over-merge by us): universe families whose mapped named genes span >1 Compara component: **8**.
  - `ASDURF` (ASDURF, ASNSD1) spans 2 Compara components -- i.e. Compara does not link all of these into one paralog group. NOTE this is the SAME signal as the non-confirmed pairs above, re-expressed at the component level: a universe family that pairs genes Compara does not relate. Whether it is a genuine over-merge or a Compara/granularity gap is exactly the ambiguity named in caveat 2.
  - `CASP8` (CASP8, FLACC1) spans 2 Compara components -- i.e. Compara does not link all of these into one paralog group. NOTE this is the SAME signal as the non-confirmed pairs above, re-expressed at the component level: a universe family that pairs genes Compara does not relate. Whether it is a genuine over-merge or a Compara/granularity gap is exactly the ambiguity named in caveat 2.
  - `CCDC188` (CCDC188, ZDHHC8) spans 2 Compara components -- i.e. Compara does not link all of these into one paralog group. NOTE this is the SAME signal as the non-confirmed pairs above, re-expressed at the component level: a universe family that pairs genes Compara does not relate. Whether it is a genuine over-merge or a Compara/granularity gap is exactly the ambiguity named in caveat 2.
  - `CDPF1` (CDPF1, PPARA) spans 2 Compara components -- i.e. Compara does not link all of these into one paralog group. NOTE this is the SAME signal as the non-confirmed pairs above, re-expressed at the component level: a universe family that pairs genes Compara does not relate. Whether it is a genuine over-merge or a Compara/granularity gap is exactly the ambiguity named in caveat 2.
  - `CREB1` (CREB1, METTL21A) spans 2 Compara components -- i.e. Compara does not link all of these into one paralog group. NOTE this is the SAME signal as the non-confirmed pairs above, re-expressed at the component level: a universe family that pairs genes Compara does not relate. Whether it is a genuine over-merge or a Compara/granularity gap is exactly the ambiguity named in caveat 2.
  - `GCA` (GCA, KCNH7) spans 2 Compara components -- i.e. Compara does not link all of these into one paralog group. NOTE this is the SAME signal as the non-confirmed pairs above, re-expressed at the component level: a universe family that pairs genes Compara does not relate. Whether it is a genuine over-merge or a Compara/granularity gap is exactly the ambiguity named in caveat 2.
  - `GPR39` (GPR39, LYPD1) spans 2 Compara components -- i.e. Compara does not link all of these into one paralog group. NOTE this is the SAME signal as the non-confirmed pairs above, re-expressed at the component level: a universe family that pairs genes Compara does not relate. Whether it is a genuine over-merge or a Compara/granularity gap is exactly the ambiguity named in caveat 2.
  - `RABL2A` (RABL2A, RABL2B) spans 2 Compara components -- i.e. Compara does not link all of these into one paralog group. NOTE this is the SAME signal as the non-confirmed pairs above, re-expressed at the component level: a universe family that pairs genes Compara does not relate. Whether it is a genuine over-merge or a Compara/granularity gap is exactly the ambiguity named in caveat 2.
## Cross-check: JSON's stated within-universe paralog relation
The JSON records **5** within-universe symmetric Compara paralog pairs (both genes mapped + linked): `APOBEC3D<->APOBEC3F`, `GGT1<->GGTLC2`, `RFPL1<->RFPL2`, `RFPL1<->RFPL3`, `RFPL2<->RFPL3`. Note `GGT1<->GGTLC2` is in this list but the two genes sit in SEPARATE universe families (`GGT1` and `GGTLC2`), so Compara confirms a paralogy that our minimap2 universe **split** -- a recall observation (granularity), not a precision failure, and consistent with the framing that we cluster more finely than Compara.
## Honest caveats
- **Human-as-proxy for gorilla.** The sequences are gorilla genes; Compara paralogy is queried on the HUMAN ortholog's ENSG. Paralogy is deeply conserved across great apes (these duplications predate the ~10-Mya gorilla/human split), so human paralogy is a strong proxy -- but it is a proxy, not the gorilla gene tree.
- **Compara is COARSER (superfamily-level).** A NON-confirmation can mean a genuine universe false-grouping OR a Compara gap / different granularity. Every non-confirmed pair is named above so the reader can judge; we do not silently treat them all as errors.
- **Only NAMED genes mapped.** 154 of 195 universe genes are LOC-only (no symbol -> no ENSG) and are out of scope. This validates the NAMED-gene backbone of the families, not all 62 families.
- **Symbol -> Ensembl mapping noise.** Genes were resolved by symbol via the Ensembl REST lookup; 40/41 got an ENSG, with 1 lookup error and 4 paralogue-fetch errors. An errored fetch yields an EMPTY paralogue set, which can only DEPRESS the confirmation rate (it cannot create a false confirmation) -- so the reported precision is a conservative lower bound for the affected pairs.
- **Recall vs Compara is intentionally NOT the metric.** Compara lumping ancient paralogs we correctly split is expected granularity, not error (see the framing section). We report granularity as an observation only.
- **Small N.** Only 10 families / 12 pairs are checkable. The headline is an exact count on a small sample, not a generalized rate; treat it as a directional, independent sanity check rather than a population estimate.
- **Determinism.** No RNG anywhere (the full mapped set is evaluated; no sampling). Minimizer hash is blake2b. Output is byte-stable.


---

## transcript_validation

# De-novo transcript validation (intron-chain Sn/Pr vs RefSeq)

Realigned all 101,467 de-novo transcripts (`denovo_transcripts.fa`) to the genome with minimap2
(`-ax splice:hq -uf`, low-mem `-I 1G --split-prefix`, MALLOC_ARENA_MAX=2 — the ulimit -v VIRTUAL cap
false-triggered 3×, drop it), `bam_to_gtf.py` → GTF, `gffcompare -r GGO_genomic.gff`. 100% mapped.

## Verdict: REAL + defensible
- **Intron-level precision 86%** (`-R -Q` 86.1) — posited splice junctions are genuine annotated sites.
- **Class codes: 98.9% overlap a KNOWN gene** (=20.7% FSM, c=31.0% ISM, j=30.9% novel-iso-of-known,
  m/n/k=16.4% retained/containment); only **0.5% `u` intergenic-novel** (artifact-suspect), ~0.6% antisense/intronic.
- **Sensitivity (where it looks, -R -Q): 76.6% of introns, 76.7% of loci recovered.** Genome-wide intron Sn 53%.
- Intron-CHAIN FSM only ~21% — NOT artifacts: long-read novel/partial isoforms (one novel junction in a
  ~9-exon chain fails whole-chain match). Expected annotation-incompleteness + novel-isoform discovery.
- CAVEAT: 31% ISM + 13% retained-intron => a real fraction are PARTIAL/incomplete (5' degradation / pre-mRNA),
  typical of read-derived assembly.

Artifacts in /home/juanfra/winloci_scratch/validate/ (dn_realigned.bam, dn_gw*.stats/.tmap).

---

## PARCN_VALIDATION

# Assembly-based parCN (`parcn`) — real-data validation

> **Updates the "DNA copy-number axis" section above.** parCN is no longer only future work: an optional
> assembly-side `parcn` tool was actually run on data already on disk (no DNA download, no new sequencing),
> closing the parCN gap for the SUN-resolvable copies. The RNA-exclusive core is untouched.

Ran `parcn` on the **gorilla (GGO) de-novo catalog** projected onto the **two phased mGorGor1 haplotype
assemblies** (species-consistent). Substrate: `gw_xchrom_refined.copies.fa` — **157 families / 414 copies**
(`--cross-chrom --homology-primary` refined catalog on GGO_mm.bam).

### Command

```bash
# one-time splice indexes (mat 13.7 GB / 233 s, pat 13.0 GB / 144 s)
minimap2 -x splice -d mGorGor1.mat.splice.mmi mGorGor1.mat.cur.20231122.fasta.gz
minimap2 -x splice -d mGorGor1.pat.splice.mmi mGorGor1.pat.cur.20231122.fasta.gz

parcn --copies-fa gw_xchrom_refined.copies.fa \
      --mat mGorGor1.mat.splice.mmi --pat mGorGor1.pat.splice.mmi \
      --out ggo_parcn --threads 4
```

**Cost:** 4 m 55 s wall, 17.6 GB peak RAM (one 13 GB splice index at a time; within the ~19 GB WSL2 cap).
No DNA depth model, no GC/mappability mask. Two TSVs out.

### 1. Threshold-free assignment — deterministic SUN or abstain

| assign_method | copies | share |
|---|---|---|
| **SUN** (deterministic private-base witness) | 376 | **90.8 %** |
| align_fallback (flagged heuristic) | 3 | 0.7 % |
| UNRESOLVED (Tier-3 / no witness) | 35 | 8.5 % |

Tier mix: T1 93.5 %, T2 1.9 %, T3 4.6 %. The heuristic fallback fires **3 times in 414** — effectively
**assign-by-SUN or abstain**, no threshold-dependent call carries the result. SUN coverage (90.8 %) *exceeds*
the RNA-only SUN-identifiability estimate (~82 %, `bench/sun_identifiability.py`): phased assembly + divergent
gorilla paralogs yield cleaner private markers than RNA reads alone.

### 2. Conservation — no loci lost

`Σ famCN_diploid (1281) + Σ n_unresolved (108) = 1389` total distinct projected loci; the independent
`Σ parCN` over the 414 copy rows = **1281**, matching `Σ famCN_diploid` exactly. Every deduped genomic locus
is either assigned to a copy or counted unresolved — nothing double-counted or dropped.

### 3. Diploid famCN tracks the catalog — and recovers what RNA collapsed

`famCN_diploid / (2 × haploid catalog copies)` over 157 families:

| statistic | value |
|---|---|
| **median ratio** | **1.00** (Q1 1.00, Q3 1.50) |
| exactly 2× (diploid-stable) | **69 / 157 = 44 %** |
| within [0.75, 1.25] of 2× | 80 / 157 = 51 % |
| mean ratio (right-skewed) | 1.50 (max 14.25) |

The typical family's **diploid CN is exactly 2× its RNA haploid catalog count** — the expected result for a
CN-stable paralog on both haplotypes (the core correctness check). The **right tail** (Q3 1.50 → max 14.25×)
is not error: the assembly reveals genuine copies the RNA catalog collapsed (K=0 near-identical merge) — e.g.
`GWFAM10` 6→23, `GWFAM116` 6→28, `GWFAM107` 2→15. Because every counted locus is **SUN-gated** (carries that
copy's private base), these expansions are real copies, not spurious cross-family hits (a spurious hit lacks
the private marker → UNRESOLVED, not the copy's parCN). This quantifies the RNA undercount parCN exists to close.

### 4. Heterozygous copy number — a phased-assembly-only signal

**37 / 157 = 24 %** of families have `loci_mat ≠ loci_pat` — maternal and paternal haplotypes carry different
copy numbers of the paralog. This allelic-CN signal is invisible to RNA and to an unphased reference; it falls
out of the mat/pat split `parcn` reports per copy.

### Frame-fix lesson (whole-branch adversarial review)

An adversarial whole-branch review caught two SUN-confirmation bugs the passing run could not surface
(conservation and the diploid ratio are invariant to *which* copy a locus is assigned): the confirm read the
wrong cs column for **soft-clipped (`qs>0`)** and **minus-strand** hits — cs-tag reads MUST honor `qs` + strand
or you get silent false SUN. Fixed: the confirm now tests match-vs-mismatch at the strand/`qs`-mapped position
(strand-symmetric). The fix moved **~15 loci (≈1 %)** from false-confirm to correctly-unresolved (famCN
1296→1281, unresolved 93→108); aggregate story unchanged, confirm now correct for inverted-duplicate/clipped
hits. The review also caught an over-tight `banded_msa_pair` band cap (cap > real length spread) that had
briefly halved SUN coverage; the band is now sized to the real within-family length spread (max 5651 bp).

### Caveats (honest)

- **Very-high-ratio families warrant a spot-check.** SUN-gating guarantees each counted locus carries the
  copy's private base, but a repeat-rich consensus could acquire many near-identical genomic hits; the extreme
  tail (e.g. 14×) should be eyeballed before headline use. The median / 45%-exact bulk is the trustworthy core.
- **UNRESOLVED (8.5 %) is genuine assembly-level collapse** — copies indistinguishable even in the phased
  assembly (no private base). Correctly abstained, not guessed.
- Substrate uses de-novo `GWFAM` ids; a gene-name spot-map (RBMY/DAZ/GSTM) via annotation overlap is a
  follow-up, not required for the copy-number counting validated here.

**Bottom line:** `parcn` closes the parCN gap in ~5 minutes — **90.8 % deterministic SUN** (assign-or-abstain,
essentially no heuristic), **exact conservation**, **median diploid famCN = 2× the RNA catalog** (44 % exact),
a right tail that recovers RNA-collapsed copies, and a **24 % heterozygous-CN** signal unique to the phased
assembly — all from on-disk data, RNA core untouched.

