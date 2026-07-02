# Single-exon exclusion: trade-off quantification

**Proposal (user):** define the RNA multi-copy family over **spliced / multi-exon transcripts only**
(exclude single-exon transcripts). Motivation: single-exon transcripts are the irreducible-from-structure
wall — exon-containment cannot separate a single-exon real paralog from a single-exon domain-bridge
(MAGE is intronless). This note quantifies the trade-off with numbers. **DNA truth = VALIDATION only.**

Reproduce (deterministic, `PYTHONHASHSEED=0`; identical md5 across runs):
```
PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/single_exon_exclusion.py
```
Outputs: `bench/single_exon_exclusion.{py,tsv,json}`. RNA-only: per-member exon count is
**read-structure-derived** = last field of `member_dn` (`DN_<chrom>_<start>_<n_exon>`, built in
`bench/denovo_assemble_gate.py`), i.e. the exact transcript fed to the POA family definition. The
annotation intron index (`gene_intron_index.json`) is used **only** as a validation overlay to
characterise the biological recall pool.

## Headline

The de-novo RNA pipeline **already applies single-exon exclusion implicitly** (intron-junction collapse +
chromosome-span drop). The transcript FASTA feeding the family definition has **1 / 102,455** single-exon
transcripts (the mitochondrion, `NC_011120.1`), and the 607-family default catalog
(`family_rna_refine.tsv`) has **0 single-exon members** (minimum de-novo `n_exon` = 2). So an **explicit
single-exon exclusion is a NO-OP on the shipped catalog**: it removes 0 families, 0 residual FPs, and
costs 0 recall. The definition is *already spliced-only*; single-exon families are a recall boundary
handled **upstream**, not a precision knob.

## [1] Population

| level | single-exon | multi-exon | note |
|---|---|---|---|
| de-novo skeletons (`denovo_skeletons.tsv`) | **26** | 164,769 | all 26 single-exon "loci" are **per-chromosome aggregation buckets** (spans 44–236 Mb + 16.4 kb mito; 0.7k–37k unspliced reads each) — the pipeline never resolves single-exon reads into individual loci; it dumps them per-chromosome and drops them as chromosome-spanning artifacts |
| transcript FASTA fed to family def | **1 / 102,455** | 102,454 | the single one = mitochondrion |
| catalog (607 fam / 1791 members) | single-only **0**, mixed **0** | multi-only **607** | catalog single-exon members = **0** |

**Exclusion removes: 0 families, 0 members.**

## [2] Residual FP removed by exclusion

Named residual-FP roster (`family_level_pr_current.json`) = **7 entries** (4 multifam + 3 oversize; 6
distinct DNA-confirmed over-merge blocks — the MAGE-A sub-cluster appears in both roles). **Removed by
exclusion = 0.**

| FP block | status | verdict |
|---|---|---|
| **MAGE-A sub-cluster LOC129529978 / LOC129529986** | annotation `n_intron` = **2 each**; de-novo n_exon 2–3 | **multi-exon → NOT removed (survives)** |
| **GSTM2 hub** | GSTM2 `n_intron` = **7** (co-members 15/27-intron) | **multi-exon → survives, still needs exon-containment** |
| **fam17** (block 17, 27-gene / 16-protein-family repeat-bridge hub) | only 1 of 27 genes is annotation-single-exon; `all_single_exon_denovo` = False | **survives** |

Net: exclusion **removes zero residual FPs → precision does not improve.** The MAGE double-role is
handled correctly: LOC129529978/986 sit in *both* the multifam and oversize roles, both are 2-intron
(multi-exon), so exclusion removes **neither** role.

## [3] Recall cost (the crux)

**Against the diploid CN oracle** (`diploid_cn_oracle.tsv`, `asm_hapCN>=2`): 57 named multi_copy genes;
**exactly 1** is single-exon in annotation — **LOC134758386** (hapCN 2, one 3184 bp exon, `n_intron` 0).
It is **recovered as family 80** via a **4-exon** de-novo transcript (the assembly recruited neighbors
over ~29 kb) → it **survives operational exclusion; it is lost only under a strict annotation-keyed rule**.

- **MAGE is NOT at risk against the oracle:** the only MAGE oracle gene, **MAGEA9, is 5-exon**
  (`n_intron` 4); MAGEA1/A4/A12 are 3-exon (spliced 5′UTR), MAGEA10 is 5-exon. On this
  reference/annotation MAGE-A is genuinely multi-exon and **dodges the wall**. The real MAGE-A family is
  catalog **family 550** (de-novo n_exon 2–3, all multi-exon).
- **Operational recall cost = [] (none). Strict-annotation counterfactual recall cost = [LOC134758386]
  (1 gene).**

### The concrete biological recall boundary (NOT in the oracle — honestly named)

The diploid oracle contains **no** intronless multi-copy family at CN ≥ 2 other than LOC134758386, so
the classic single-exon families are **not measurable against that oracle**. But they **are present in
this gorilla annotation** and are the concrete casualties a spliced-only definition loses **upstream**
(each carries no splice junction → bucketed per-chromosome → dropped before family definition, unless the
de-novo assembly recruits a neighboring exon). Enumerated from the intron index the script loads:

| family | intronless / total in annotation | in catalog | fully lost upstream |
|---|---|---|---|
| **IFNA** type-I interferon cluster | 4 / 6 (IFNA2, IFNA6, IFNA8, IFNA21) | **0** | **YES** |
| **TAS2R** bitter-taste GPCRs | 15 / 18 | **0** | **YES** |
| **H1** linker histones | 5 / 6 (H1-1, H1-2, H1-3, H1-5, H1-6) | **0** | **YES** |
| **MAGE-B** cancer-testis (MAGEB5, MAGEB17 intronless) | 2 / 9 | **0** | **YES** |
| **PCDHB** protocadherin-β cluster | 9 / 10 | 3 (PCDHB2/5/15, via multi-exon de-novo) | no (partial) |

This is the honest recall cost the user rightly worries about **on this testis IsoSeq substrate**: IFNA,
TAS2R, H1, and MAGE-B5/B17 are genuinely intronless multi-copy families that the RNA read-structure
pipeline **cannot** recover as multi-exon transcripts, and they are already absent from the 607-family
catalog — the loss is real and already paid upstream, **not** introduced by the proposed exclusion knob.
The PCDHB row shows the escape hatch: an intronless family can still be recovered *if* assembly recruits
a neighboring exon (3 of its genes made it in via multi-exon de-novo transcripts).

## [4] Retrocopy / repeat / artifact vs genuine split

- **De-novo single-exon resolved real loci = 0.** All 26 skeleton single-exon rows are chromosome-
  aggregation buckets; the "assembly-truncation" class is **never resolved** — it lives in the dropped
  bucket. This is exactly the irreducible-from-structure wall the user described.
- **Annotation single-exon genes = 9,821 / 40,703 (24.1 %)**; recovered in the catalog = **28**, all via
  multi-exon de-novo transcripts (de-novo n_exon 2–6, none = 1). Class split of the 28:
  **retrocopy_signature 14, genuine_single_exon_tandem 6, mixed 6, unknown 2**.
- **Genome-wide cDNA-homology edge split** (reused byte-for-byte from
  `family_def_retrocopy_filter.json`): **1931 one-side-intronless retrocopies dropped** (processed-
  pseudogene noise, already handled by the shipped filter) vs **3267 both-intronless genuine
  single-exon tandems kept** (IFNA/PCDHB/ELOA/HNRNPCL-like) vs 12,212 both-spliced. The genuine
  single-exon family pool exists genome-wide but is **not the recall-limiting factor for this catalog** —
  it is limited upstream by the read-structure drop, not by the retrocopy filter.

So the recall cost is only the **genuine single-exon families**, not the retrocopy/exonized-repeat/
assembly-truncation noise (which is correctly discarded).

## [5] Net family-level P/R (vs diploid oracle)

| scenario | P (dedup) | P (task-formula) | R | Δ |
|---|---|---|---|---|
| **current default** | 0.875 | 0.8542 | 0.8421 | — |
| **operational exclusion** (drop de-novo single-exon members) | **0.875** | **0.8542** | **0.8421** | **identical (no-op)** |
| strict-annotation counterfactual | 0.8723 | — | 0.8246 | P −0.003, R −0.018, **loses LOC134758386** |

Operational exclusion is byte-identical to the default (0 single-exon members to drop). The strict-
annotation counterfactual **lowers both precision and recall** (removes 0 FP, loses 1 real oracle
family) — a **net loss**.

## Verdict + recommendation

Excluding single-exon transcripts as proposed is **neither a win nor a loss on the shipped catalog — it
is a no-op**, because the de-novo assembly + intron-junction collapse **already enforces it**. It removes
**0 residual FPs** (the MAGE-A sub-cluster, GSTM2, and fam17 are all multi-exon and survive → precision
unchanged at **P 0.875 / R 0.8421**) and costs **0 recall** at the operational (RNA read-structure)
level. A stricter annotation-keyed rule would cost exactly **1 named oracle family (LOC134758386)** for
**0 precision gain** — a net loss.

**Recommendation: do NOT add an explicit single-exon exclusion as a precision lever — it buys nothing
here. State the scope honestly instead: the RNA family definition is *already spliced-only*, and
single-exon families are a recall boundary handled upstream, not a precision knob.**

**Honest recall cost:** the biological cost the user rightly worries about on this testis IsoSeq
substrate is **real but already paid and invisible in the catalog**. Genuinely intronless multi-copy
families in this gorilla annotation — **IFNA type-I interferon cluster (4 genes), TAS2R bitter-taste
GPCRs (15/18 intronless), H1 linker histones (5/6), MAGEB5/MAGEB17** — are all **absent from the
607-family catalog** because they carry no splice junction and were dropped upstream. On this specific
reference/annotation **MAGE-A dodges the wall** (annotated 3–5 exons; real family = catalog fam 550),
and the only intronless oracle multi_copy family (**LOC134758386**) is recovered via a multi-exon de-novo
transcript. If recovering intronless families like IFNA/TAS2R/H1 is a goal, that is a **separate,
sequence-homology (non-splice) recovery track** — it cannot be bought or lost by a single-exon
inclusion/exclusion knob at the family-definition stage.
