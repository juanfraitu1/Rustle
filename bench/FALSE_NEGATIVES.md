# False negatives: what the pipeline misses, and why

**Date:** 2026-07-11. Substrate: gorilla (GGO) testis Iso-Seq, `GGO_mm.bam` vs `GGO.fasta`, current binary
(`c394bfd`, refine-by-default). A companion to `bench/GW_CATALOG_FP_AUDIT.md` (false POSITIVES). This enumerates
the classes of real multi-copy families / copies the method does **not** recover, with the cause and a measured
size where possible. "FN" here means a real object we miss — not necessarily a fixable defect; several are the
correct behaviour of a sample or an information limit.

## Baseline — what we DO recover (refine-by-default)

The six flagship testis-expressed families all recover exactly under the new refine-by-default `copy_assign`, and
both single-copy controls stay silent:

| GSTM | MAGEA | DAZ | RBMY | TSPY | PCDHB | EEF1A1 | SRGAP2 |
|---|---|---|---|---|---|---|---|
| 3 | 2 | 2 | 6 | 5 | 5 | 0 | 0 |

Refine kept every real family (each `1 → 1`) and added no false negative to the flagships; it even cleaned an E_r
over-call at the SRGAP2 control (`3 → 0`). So the FN classes below are the *edges*, not the core.

## Detection FNs — real families / copies we never emit

### 1. Not transcribed in this tissue (the largest class — and correct)
Most annotated multi-copy families are silent in testis: HOX, TAS2R, S100A, MMP, KRT, WFDC, DEFB, SERPINA/B,
PRAMEF, XAGE, CEACAM, brain-specific expansions, etc. ~24 of 30 sampled named families carry 0 reads
(`FAMILY_ARTIFACT_AUDIT.md`). Returning nothing there is a **sample limitation, not a method failure** — a
multi-tissue panel would recover them.

### 2. Coverage floor — expressed, but below the assembly gate
A copy present and transcribed but too shallowly to clear `GATE_MIN_READS`(3)/`locus_support` is missed:
**RBMY proximal 2/6** (77 reads over 1 Mb — 4 copies not deep enough), **TSPY 6th copy** (c276 had 0 reads in
this sample; the ground-truth sim shows it *would* resolve if expressed), **CDY/HSFY** (0–10 reads)
(`YAG_CHECK.md`). This is the λ-floor; more depth recovers them.

### 3. Default-mode "globin problem" — unique-mappers need `--homology-primary`
Families whose copies each map **uniquely** (high MAPQ) form **zero** read-conflict edges, so the default E_c
oracle emits nothing: **GSTM, MAGEA, RFPL → 0 families under default**, recovered only with `--homology-primary`
(`FAMILY_SPOT_CHECK.md`). Refine does not change this (it is an oracle issue, not a gate issue). *Open question:
whether E_r (`--homology-primary`) should be the default now that refine cleans its domain-bridge over-calls.*

### 4. ⭐Refine recall cost (measured this audit): 13 real families dropped
Making refine the default removes false positives at a recall cost. Of the **42 families refine drops**
genome-wide (124 → 86), a 42-agent adversarial classification (align copies + check annotation) found **29
correctly dropped** (repeat-bridges + gene-splits) and **13 real paralog families wrongly lost** — ~10% of the
124-family raw catalog. Three causes:

| cause | n | examples | fixable? |
|---|---|---|---|
| **Partial transcript models → exon-sum coverage < 0.50** despite ~100% identity | 7 | EOLA1/EOLA2 (99.96% id, 42% cov), ZNF74/ZNF74-like (99.7%, 24%), RABGEF1, α2-macroglobulin, FRG1-like, GRAP | **yes** — the copies ARE near-identical segdups; the de-novo transcript models just don't overlap enough. Use genomic-span coverage (`--refine-introns`) or `max(exon-sum, genomic)` coverage for the floor. |
| **Genuine divergence below the identity floor** | 5 | ARMCX1/ARMCX6 (65% aa, 0 nt alignment), IFITM cluster, FRG1-like, KRAB-C2H2 ZNF cluster (ZNF677/761/665) | partly — the true precision/recall tradeoff; ancient paralogs below asm20 0.80 / sensitive 0.70. Protein-tier (`--protein-tail`) recovers some coding ones. |
| **Family-split logic edge case** — a real near-identical pair lost when a 3rd bridging copy is present | ~1–2 | ARHGAP23-like pair (99.2% id, 99.9% cov), PDPK1/PDPK2 (99.6%, 57% cov) | **yes** — a refine component/`distinct_locus_reps` bug: the good pair should survive on its own. |

**7 of the 13 (the coverage class) are a fixable metric artifact, not a real limit** — the copies are
homologous at 99%+ identity and only fail because the assembled transcript models are partial. The 5 divergence
FNs are the honest cost of a high-precision gate.

### 5. Under-merging / fragmentation — real family split, copy count understated
A real family emitted as two smaller families: **GBP (6 annotated → 4 + 2), TCEAL (6 → 3 + 2)**
(`FAMILY_ARTIFACT_AUDIT.md`). No false copy, but a copy-number readout of "4" understates the true 6.

### 6. Reference-absent copies (O4) — detected/flagged, not resolved
Copies absent from the reference assembly are detected and flagged only (collapsed-CNV signal or unmapped
reads), never emitted as resolved copies — copy-vs-allele needs DNA (`project_reference_absent_catalog`).

## Assignment FNs — copies found, but reads unassignable

### 7. K=0 frontier — exonically-identical copies
The copies are **detected** (χ_H correct) but their reads cannot be assigned to a specific copy — they tie.
**TSPY**: 4 of 6 copies are 100.000% identical, so a read carries zero copy-specific signal; the pipeline
correctly abstains rather than fabricating (`TSPY_SIMULATION.md`, confirmed by minimap2 *and* winnowmap). This is
an information limit, not a method defect; the escapes are DNA, aggregate quantification, or a copy divergent
enough to map uniquely.

## Closed FN classes (were misses, now fixed)
- **Inverted duplicates** — `colocated_families` used to split on strand, dropping every inverted pair; the
  chrom-only fix recovered MAGEA (0 → 2) (`o2_strand_fix`).
- **DAZ2** — a 5′-truncated collapsed copy the assembly gate lost; `locus_support` + chimeric-bridge exclusion
  recovered DAZ (0 → 2 copies) (`daz2_locus_support`).

## The bottom line — the precision/recall trade refine buys

Refine-by-default removes **11+ whole-family false positives** (gene-splits, repeat-bridges) at a cost of **13
real families** (recall), of which **7 are a fixable coverage-metric artifact** and ~2 are a fixable split-logic
edge case. So the *net* honest recall cost of the current gate is ~5 genuinely-divergent paralog families
genome-wide, and the coverage/split cases are worth fixing (use genomic-span coverage for the floor; fix the
component drop of a valid sub-pair). The dominant FN overall is **not** refine but **tissue silence** (families
not expressed in testis) and the **coverage floor** — both resolved by more data, not more method.

Reproduce: the 42-agent classification ran over the families in `gw_audit.copies.tsv` (raw) absent from
`gw_audit_refine.copies.tsv` (refined); flagship re-run is `copy_assign … --homology-primary` per region.
