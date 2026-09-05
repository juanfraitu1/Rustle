# PRE-REGISTRATION — RNA corroboration of the MCL catalog (written 2026-09-04, BEFORE the run completed)

**Run:** `mcl_families --paf allgenes.asm20.paf --gff GGO_genomic.gff --min-exonic-bp 1 --bam
fibroblasts/GCA_029281585.2_flnc_mm.bam` (FLNC IsoSeq, 23.2 GB, NC_ contigs verified present for all 3
pilot contigs). Corroborated = a member with >= 3 reads carrying an ALIGNED BLOCK inside it
(`corroboration_min_reads = 3`; ⚠ block, not span — `N` is an intron spliced out).

## The reasoning that makes this a real test

The project's framing is **"DNA PROPOSES, RNA DISPOSES"**. That requires RNA to be able to REJECT a
DNA-proposed cluster. ⚠ I do not think per-member corroboration can do that, and I am recording why
before seeing the numbers.

Corroboration asks **"is this gene expressed in fibroblast?"** — a property of the GENE. Family
membership is a property of the PAIR. The two are not the same question, and the adjudicated cases
predict the failure directly:

- **MCL2** is 33 unrelated but *well-characterised, mostly housekeeping* genes (ATP5F1A, CEP76, ELL,
  TRMT1, NDUFA13, NDUFAB1 …). Those ARE expressed in fibroblast ⟹ a **known-false** cluster should
  corroborate WELL.
- **MCL0** is a chorionic-gonadotropin / NTF4 **pseudogene** array. Those are NOT expressed in
  fibroblast ⟹ a real (if falsely-merged) tandem array should corroborate at ~**0**.

## Pre-registered predictions

1. ⭐ **PRIMARY: corroboration will NOT separate cohesive clusters from hairball slices.** Operationally:
   the mean corroborated fraction of clusters with `frac_in >= 0.90` will differ from that of clusters
   with `frac_in < 0.75` by **less than 0.20**, and/or the sign will be inconsistent with cohesion.
2. **MCL2's surviving members corroborate HIGH** (>= 0.5) despite being an adjudicated repeat clique.
3. **MCL0's members corroborate LOW** (<= 0.2) despite being a real tandem array.
4. **NPIP (MCL1) corroborates moderately-to-high** — it is expressed in fibroblast (prior: 36/43 = 0.837
   in the genome-wide run).

## What each outcome means

- **If 1 holds** ⟹ ⛔ **per-member RNA corroboration CANNOT dispose of DNA-proposed false merges.** The
  "RNA disposes" half of the framing is not implemented by the instrument currently called
  corroboration, and saying otherwise would be an overclaim. The honest role of this instrument is
  **necessary-condition evidence that a family is transcribed at all** — which is what makes this an RNA
  thesis — not adjudication of membership.
- **If 1 fails** (cohesion does predict corroboration) ⟹ a genuinely non-circular DNA/RNA agreement, and
  the strongest single result available for the MCL definition.

⚠ **Absence of expression is NOT evidence of a false family.** Tissue-specificity and pseudogene status
both produce zeros. Any use of a LOW corroboration rate as evidence AGAINST a cluster is forbidden by
this pre-registration.
⚠ Confounds not controlled by this run: cluster size, biotype (pseudogene vs protein-coding), gene
length, and the single-tissue substrate.

---

# OUTCOME (same day, after the run). 4/4 predictions HELD.

| # | prediction | outcome | result |
|---|---|---|---|
| 1 | cohesion will NOT predict corroboration (\|Δ\| < 0.20) | **HELD** | Δ = **−0.0792**, and the sign **INVERTED** |
| 2 | MCL2 (known-FALSE clique) corroborates HIGH (>= 0.5) | **HELD** ⚠ | **0.5758** (19/33) — ⚠ upper bound, see caveat |
| 3 | MCL0 (real pseudogene array) corroborates LOW (<= 0.2) | **HELD** | **0.0000** |
| 4 | NPIP corroborates moderate-high | **HELD** | **1.0000** (43/43) |

**Primary result, with the size control — the inversion GROWS with cluster size:**

| cluster size | cohesive (frac_in >= 0.90) | hairball slice (frac_in < 0.75) |
|---|---|---|
| 3-4 | 0.852 (n=50) | 0.864 (n=38) |
| 5-9 | 0.768 (n=13) | **0.903** (n=12) |
| 10-14 | 0.750 (n=4) | **0.971** (n=7) |
| 20-24 | **0.565** (n=6) | **0.970** (n=1) |

Overall mean corroboration **0.8478** over 142 clusters; only **12/142** have zero support.

⚠ **CAVEAT ON PREDICTION 2.** MCL2 was eliminated by the exonic clause (0/33 survive), so it has no
cluster to score; its members were scored as GENES with `samtools view -F 2308 -c`, which is **SPAN
overlap, not the aligned-block criterion the pipeline uses**. A read splicing straight over a gene counts
there and must not — the same error inflated a headline 3.4x in §6cm. **0.5758 is therefore an UPPER
BOUND** and the prediction is supported by a loose instrument, not the shipped one.
