# O3 — copies absent from the genome, present in the transcriptome: the evidence chain

**Status 2026-08-14.** O3 asks whether a gene copy can be detected that is missing from the assembly but
present in the individual's transcriptome. This document is the complete evidence chain, in the order an
examiner should read it: what the field accepts as evidence, what we found, what we retracted, the two
independent negative results, and — new on 2026-08-14 — the **whole-genome positive control** that
converts "the detector did not fire" into a measured sensitivity.

Read the limits in §7 before quoting any number. Every rate carries its 95% CI; where a number was
retracted the retraction travels with it.

---

## 1. What the field accepts as evidence — and what it does not

A literature sweep (53 searches, ~30 primary fetches) found **no paper that discovers a gene copy missing
from an assembly using transcriptome data as the evidence.** The "missing from the assembly" literature
is large, mature, and **entirely DNA-side**. Its three standards:

| | standard | representative work |
|---|---|---|
| **S1** | **produce the sequence** — re-assemble and show the copy (gold standard, 4 of 5 papers read) | T2T-CHM13 (Nurk 2022); Vollger 2022 SDs 5.4%→7.0%; Yoo/Rhie 2025 apes |
| **S2** | **read-depth excess + paralogous-sequence-variant structure** | WSSD/Bailey 2002; SDA (Vollger 2019); NucFreq/NucFlag; Flagger |
| **S3** | **molecules mapping uniquely once the copy exists** | WASH1/GPRIN2 — via **proteomics, not RNA** |

⚠⚠ **No collapse detector in the literature uses CLIPPING**, and transcript-alignment anomaly is not
accepted as evidence of a missing copy. Our original statistic (FARDIV ∧ FARCLIP) sits outside the
field's evidence standards, which is why §4 and §5 matter.

⚠ **Keep four things apart; the field blurs them.** (a) missing from the ASSEMBLY — our case ·
(b) absent from the REFERENCE but present in other individuals — pangenome PAV · (c) transcripts that
merely align poorly · (d) novel isoforms of sequence that is present.

⚠ **A mechanism correction we must carry.** "The orphan read leaves the locus" is **wrong for collapsed
loci**: FixItFelix (Behera 2023) measures depth in collapsed regions **~1.5× higher**. The correct
statement is that the orphan is **absorbed by its best-matching paralog wherever that is** — confirmed
independently in §5 at **1.75×**.

## 2. The transcript-side screen: 0/915, and the mechanism is the finding

On the matched arm the shipped detector fired on **0 of 915** loci. This is not a null result about the
phenomenon; it is a statement about the statistic:

* **FARCLIP is dead on a whole genome** — median **0.0006** (matched) / **0.0000** (cross) against a gate
  of **0.05**, so the FARDIV ∧ FARCLIP conjunction cannot fire.
* **Why**: on 3.6 Gb a read that does not fit locus *u* finds a better primary home elsewhere instead of
  being forced onto *u* with a clip. The orphan pile-up that every excision positive exercised is exactly
  what whole-genome alignment prevents. ⟹ **the mini-reference did not merely flatter the rate, it
  manufactured the signal.**

The published mini-reference performance (TPR **41.13%** [37.02, 45.37] per arm at a substrate-matched
FPR of **3.29%** [2.03, 5.27]) is therefore an **upper bound**, and was labelled as such.

## 3. The candidates: retracted, and fully explained

Three loci carried transcripts with no home above 96% anywhere in the individual's **diploid** genome
(`_pri` + `_pat` + `_mat`), against ordinary loci at **0.9987**. All are now withdrawn, with causes:

* **GWFAM488:0 and GWFAM382:2 are IMMUNOGLOBULIN.** Somatic hypermutation + V(D)J recombination is the
  one biological process that genuinely puts sequence in a transcriptome that is not in the genome.
  GWFAM488:0 is a single clonally expanded, productively rearranged λ light chain (IGLV3-1/IGLJ2/IGLC3);
  63/69 reads share an identical 37 bp non-templated V–J junction, and V is mutated **1.77–3.27×** more
  than the C exon. **Every germline segment is present in mGorGor1** (IGLJ2 identity **1.0000**), while
  the transcript sits 7–11% from both human and gorilla germline — no phylogeny puts a sequence 7% from
  two germlines that are 2% from each other.
* **GWFAM382:2** is IgG heavy chain whose entire support is **run-exclusive** (815/815 reads from 1 of 4
  pooled SRA runs, P = 10⁻⁶⁰·²) ⟹ index hopping / cross-library contamination. The substrate has no
  B-lineage cells (CD79A 6, MS4A1 13 vs COL1A1 183,670).
* **Two of three were killed by MIS-CHAINING alone.** ⚠ The shipped rule was blind by design:
  `is_giant_intron_mischain` fires only when the giant junction has **fewer than 3 reads**, and the
  827,011 bp mis-chain at GWFAM244:2 is carried by **97** primary reads. A *popular* mis-chain is
  evidence of *systematic* mis-chaining. Fixed in both languages with one source of truth.

⟹ **"Transcripts that align poorly" ≠ "copies missing from the assembly." O3 had demonstrated the first,
not the second.** Any future candidate must be screened against IG/TR loci **first** and run-exclusivity
**second**.

## 4. The DNA-side check: an independent zero, and it was predicted

Ran the field's own **S1** standard — assembly-vs-assembly, the method FixItFelix used on GRCh38 — over
three assemblies of the same gorilla (`_pri`/`_pat`/`_mat`).

* **Zero loci called collapsed at 816 of 817** with a flank-anchored orthologous interval. Zero tolerance:
  0 positives / 3 negatives, sign test p = 0.25. At the published ±100 kb tolerance: 6 up vs 5 down,
  p = 1.00. Six compare-phase candidates **withdrawn** as tolerance artifacts.
* **No absence call exists anywhere** — every interval present in both assemblies, **0 bp of N**, >100 kb
  from any sequence end ⟹ the fragmentation confound cannot manufacture a call.
* ⭐ **The zero was predicted.** Yoo/Rhie report *"1–2 Mbp of collapse per haplotype assembly"*; these
  windows span **1.1224%** of the assembly ⟹ the prior predicts only **0.47–0.94** collapses. Finding zero
  is what you expect. **The run is underpowered by construction — a statement about the compartment, not
  the instrument.**
* Instrument validated three ways: FP floor **0/817** [0.00, 0.47]% yet fires on **11.42%** of
  span-matched random intervals; recovers Yoo/Rhie's gorilla-specific MAPKBP1/SPTBN5/PLA2G4B expansion
  **exactly at 8/9/9** copies; synthetic control 0.923–0.985 detection given a visible probe.
* ⚠ Binding constraint is **probe visibility (16–21%)**, not detection: unconditional sensitivity ~0.17,
  and the product is flat whatever operating point is chosen. **The lever is a larger compartment, not a
  better instrument.**

⟹ **The transcript-side 0/915 is CONCORDANT with an independent DNA-side zero, not contradicted by it.**

## 5. The whole-genome positive control (2026-08-14) — the control every prior number lacked

**Design.** Delete one copy of **162 two-copy families** from the full 3.6 Gb genome and realign the
matched fibroblast IsoSeq. Two-copy families are 71% of the catalog and the biggest coverage hole: delete
one copy and the survivor has no family, so O1's definition cannot even be invoked on it.

*Substrate.* mGorGor1 v2.0_pri (KB3781 "Jim"); reads are **the same cell line as the assembly's DNA**
(hom-alt 13.50×, het control 1.096). Realigned with the BAM's **exact `@PG` command**, only `-t 28`→`-t 4`.

*Pre-registered mask rule.* Mask the copy with the **larger start coordinate** — deterministic and
uncorrelated with depth, span and divergence. ("Mask the copy with fewer reads" would tie the choice to
the outcome.)

*Gates, all passed before any readout.* **6,851,713 bp flipped to N == expected exactly**, with **0
pre-existing N** ⟹ every excised copy was fully assembled sequence, no gap confound; contig names and
lengths unchanged. ⭐ **32 of the original 194 families were EXCLUDED** because the copy to be masked
overlaps another catalog locus — **30 overlap their own sister**, and masking those would have destroyed
the survivor the comparison depends on.

### 5.1 A deleted copy has two fates, and both are measurable

| fate | families | signature |
|---|---|---|
| **ABSORBED** | **104/162 = 64.2%** | orphans concentrate on **one** paralog (median concentration **0.9667**) |
| **ORPHANED** | **54/162 = 33.3%** | median **92.73%** of that copy's reads become **UNMAPPED** |
| scattered | 4/162 = 2.5% | neither |

* ⭐ **The ghost is real and matches the literature.** Median depth ratio at the destination **1.75×**
  (q25 1.30, q75 2.83, n = 88) against FixItFelix's measured ~1.5×; **50/88 = 56.8%** fall in a
  ghost-shaped [1.2, 3.0]. ⚠ 16/104 destinations had zero baseline reads and 12/104 exceed 5× on a tiny
  baseline — **quote the median, never the mean**.
* ⭐⭐ **The missing-copy signature is UNMAPPED READS, not clipping.** Pooled **61,512 / 178,145 =
  34.53%** of orphans become unaligned, and they were **confidently mapped before** (baseline MAPQ-0
  fraction ≈ 0). Deleting a copy *destroys* multimapping rather than creating it: with one target left
  the read is no longer tied.

### 5.2 An unplanned, direct measurement of O1's incompleteness

**Of the 104 absorbed families, only 48 (46.2%) sent their orphans to the catalog's own sister; 56
(53.8%) sent them to a paralog the family does not contain.** Worked case — **GWFAM250**: the masked copy
(chr2:139 Mb, 16,882 reads) sent 16,377 reads to **chr6:138 Mb**, a real locus already carrying 36,569
baseline reads, while the catalog's designated sister sits on chr8.

⭐ **Why**: `E_r` admits an edge at identity ≥ 0.60 and coverage ≥ 0.50 — **far more permissive than what
minimap2 needs to place a read as primary.** Two loci can be one `E_r` family while their reads never
cross-map. ⟹ **the family graph and the read-competition graph are not the same graph**, and this run
measures the gap for the first time. Consistent with O1 reach **0.5500** [0.3983, 0.6929].

## 6. The S2 detector, built and held-out tested

**Statistic — the within-pile divergence mixture.** 1-D 2-means on the per-read `de` vector; fire iff the
higher-divergence cluster holds a fraction in [0.10, 0.50] **and** the cluster separation ≥ 0.01. The
separation term guards the "uniformly divergent pile" trap, where a 2-means split returns two clusters
regardless.

⭐ **Depth is deliberately NOT used.** The 1.75× ghost needs a "before" that a real case does not have,
and expression is not dosage ⟹ depth requires a DNA control (parCN/famCN). The divergence mixture is
self-contained. Depth is reported alongside, never used to decide.

⚠ **Confound cleared first.** Migrant reads came through minimap2 2.30, native reads through 2.31, so a
mixture detector could have learned the *version*. On 42,393 reads that stayed put in both runs,
**42,390 = 99.99% have byte-identical `de`**, median delta exactly 0.

*Positives* = destination pile (its own native reads + the reads that migrated there). *Paired negatives*
= the same loci, native reads only — within-locus pairing removes locus identity, depth, gene and
expression level as confounds. *Independent negatives* = panel loci that are neither masked nor a
destination. Thresholds fixed on TRAIN families; **all rates below are held-out TEST**.

| | rate | 95% CI |
|---|---|---|
| **TPR** | **0.2703** (10/37) | [0.1540, 0.4298] |
| FPR, paired negatives | 0.0541 (2/37) | [0.0150, 0.1770] |
| **FPR, independent negatives** | **0.0200** (1/50) | [0.0035, 0.1050] |

### 6.1 Sensitivity is set by the hidden copy's divergence, not its abundance

The **missed** positives carry **more** migrant reads than the detected ones (median m_true 0.4552 vs
0.3763) ⟹ abundance is not the limit.

| hidden copy's divergence from its host | TPR | 95% CI |
|---|---|---|
| **≥ 0.01** | **0.4500** (9/20) | [0.2582, 0.6579] |
| **< 0.01** | **0.0588** (1/17) | [0.0105, 0.2698] |

Bin detail: [0, 0.002) 0.111 · [0.002, 0.005) 0.000 · [0.005, 0.01) 0.000 · [0.01, 0.02) 0.500 ·
[0.02, 0.05) 0.500 · ≥ 0.05 0.417.

⚠⚠ **The ceiling, and it is the problem rather than the rule.** **27.71%** of real positives have a
hidden copy less than **0.002** diverged from its host; **45.78%** are below **0.01** (median true
separation 0.0144). A divergence mixture has nothing to separate there — and **collapses happen BECAUSE
copies are similar**, so the blind region is exactly where the phenomenon concentrates. This reproduces,
on a whole genome with a held-out split, the earlier mini-reference finding that clearing the gate
requires the collapsed copies to have been ~2% divergent.

## 7. Limits — state these before being asked

1. **"The reference is an average" does not apply to this substrate.** mGorGor1 is a haplotype-resolved
   assembly of **one animal**, and the IsoSeq is the **same cell line** as the assembly's DNA. That
   objection is true of **GRCh38** (~20 donors — which is why FixItFelix found collapses there) and
   remains live for the cross arm and for human work, but not here.
2. **The excision positives are synthetic.** A masked copy is a perfect deletion; a real collapse is a
   merge of two similar copies, which is harder. The 1.75× ghost and the 0.45 TPR are therefore upper
   bounds for real collapse of comparable divergence.
3. **Tissue restricts the panel.** Fibroblast (matched) and testis (cross) miss most classic
   copy-number-variable loci — CYP2D6 (liver), AMY1 (salivary), C4A/C4B (liver), DEFB4 (epithelium),
   FCGR3B (neutrophil), HBA (erythroid). SMN1/SMN2 is the one classic pathological CNV well expressed in
   fibroblasts.
4. **RNA cannot count copies.** Expression level is not dosage; a 2-copy gene can out-express a 6-copy
   one. RNA can *identify* which copies are transcribed, via PSVs — that is O2 — but copy **number**
   requires DNA (parCN/famCN, WSSD, ddPCR).
5. **The ORPHANED third is untouched by §6.** Where reads go unmapped there is no pile to run a
   within-pile statistic on; that needs a different detector (see §8).
6. **n is small.** 37 held-out positives; 162 families; one individual; one tissue.

## 8. What follows

* **A detector for the ORPHANED third.** 33.3% of deleted copies produce a coherent pile of unmapped
  reads (median 92.7% of that copy's reads). ⚠ Any such class must stay a **flagged candidate**, never a
  family: our object is nodes = genomic intervals, so a read-only cluster has no interval, no λ, and no
  seed-invariance — and **both retracted candidates in §3 were exactly such clusters.** IG/TR and
  run-exclusivity screens compulsory.
* **The 16 zero-baseline destinations** — orphans landing where nothing was expressed.
* **The Illumina depth check (S2 on DNA)** as a genuinely independent axis. ⚠ The §4 power argument
  applies to it too: this compartment predicts fewer than one collapse.

---

*Artifacts.* Pre-registration, panel, mask, scripts and results:
`/home/juanfra/winloci_scratch/o3_excise/` (`PREREG.md`, `panel.py`, `mask_genome.py`, `analyse2.py`,
`detector.py`, `RESULTS_DETECTOR.md`). Masked genome and BAMs:
`/mnt/linuxdisk/home/juanfraitu/o3_excise/`. DNA-side check: `winloci_scratch/o3_collapse/`.
