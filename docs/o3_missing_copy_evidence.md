# O3 — copies absent from the genome, present in the transcriptome: the evidence chain

> ⚠ **CATALOG PROVENANCE.** Figures quoted per **/494** or **/1415** describe the
> **superseded** 2026-07-17 catalog; the current default emits **627 families / 2,019 copies**.

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

⚠⚠ **BOTH of those figures have T1-clean replacements — do not quote them bare.** 41.13% is the
per-arm rate; for the mini-reference screen quote the **host-only arm yield 25.44% (189/743)**, and for
the false-positive rate quote **1.98%** (the fixed-denominator form of the 3.29%, which used a
called-only denominator conditioned on the method's own output). Five detection rates and four FPRs are
on file and are **not interchangeable** — see the *O3 RATE PROVENANCE* table in
`NEGATIVE_RESULTS_REGISTER.md` §6.4 before quoting any of them.

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

## 6. The divergence-mixture detector, built and held-out tested

⚠ **Naming, corrected 2026-08-14.** This was called "the S2 detector" in the first draft, and that is
wrong in a way a reviewer who knows the literature will catch: **S2 is depth excess *plus* PSV
structure**, and this detector deliberately uses **no depth at all** (§6, next paragraph). It
implements the PSV/variant half of S2 only. Call it a divergence-mixture detector.

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
   within-pile statistic on; that needs a different detector (see §9).
6. **n is small.** 37 held-out positives; 162 families; one individual; one tissue.

## 8. Possible O3 avenues — decision record (2026-08-16)

The original wording, "find copies not in the genome," is too strong for Iso-Seq alone. The observable
object is a set of **expressed read paths not adequately explained by the copies in the analysis
reference**. Whether that object is an extra genomic copy, a divergent allele, a novel isoform, or an
artifact is a second question. Therefore the proposed thesis wording is:

> **O3 — detect and flag expressed transcript paths not adequately explained by the copies represented
> in the analysis reference. Evaluate missing-copy recovery by leave-one-copy-out reference ablation.
> Do not call a natural candidate a confirmed genomic copy without independent DNA evidence.**

This wording does not assume that reference bias exists in the biological sample. It makes reference
bias a falsifiable stress condition in the controlled experiment and a hypothesis in natural data.
Also, a reference assembly is not an "average copy set": CHM13 and mGorGor1 are assemblies of particular
individuals, while a composite reference such as GRCh38 is still not guaranteed to contain every
population- or donor-specific copy.

### 8.1 Recommended main avenue: paired leave-one-copy-out ablation

This is the strongest O3 experiment available without importing another genome or using Liftoff.

1. Start from an O1 family whose members are known, expressed, non-overlapping enough to mask safely,
   and distinguishable by linked PSVs and/or splice structure. A Soto segmental-duplication family is
   eligible only after it passes the O1 gene-family rules; segmental duplication alone is not truth for
   a multi-copy gene family.
2. Make two otherwise identical references: the intact reference and one in which a known family copy
   is masked or removed. Select the hidden copy by a pre-declared rule, not from the outcome.
3. Align the same Iso-Seq molecules to both references with minimap2. Rustle receives only the degraded
   reference for discovery. The hidden sequence is evaluator-only truth.
4. Search both fates demonstrated in §5: **ABSORBED** reads that align to a sibling or outside-family
   paralog, including high-MAPQ primaries, and **ORPHANED** reads that become unmapped. A MAPQ-0-only
   screen is insufficient.
5. Cluster reads by linked residual variants and junctions and emit an auditable candidate signature.
   Reconstructing and adding a candidate transcript path to the local family graph is an optional
   stronger evaluation, not required by the detect-and-flag objective.
6. Split molecules into discovery and validation sets. On held-out molecules, require the candidate to
   improve alignment likelihood/read assignment over every represented sibling and over a
   length/composition-matched decoy.
7. Only after prediction, compare the candidate with the withheld copy. Run the same detector on the
   intact reference as the paired negative control.

Report at least: family-level detection recall; signature agreement with the withheld copy; held-out
read reassignment; false-positive rate on the intact reference; results by hidden-to-host divergence;
and the fraction on which Rustle correctly abstains. If optional path reconstruction is attempted, also
report sequence identity and covered fraction against the withheld copy. The intact/masked pair is the
causal test: removing one copy must create the signal and restoring it must remove the signal.

The current full-genome mask experiment in §5 already establishes the two read fates and gives a
held-out operating point for the absorbed/divergent stratum. A next experiment should evaluate the
complete candidate flag on held-out reads; candidate-path recovery can be reported separately as a
stretch result. Near-identical copies are an expected abstention stratum, not failures to be hidden by
the aggregate score.

### 8.2 RNA-only natural-candidate avenue: useful, but detection-and-flag only

On the unmodified reference, Rustle can report a **reference-divergent transcript candidate** when a
coherent set of molecules carries linked PSVs and/or splice structure that no represented family copy
explains. The output should include the best competing locus, number of independently supporting
molecules, phasing/linkage strength, run balance, remap identity, decoy result, and all rejection
reasons. It must screen IG/TR somatic rearrangement, run-exclusive contamination, giant-intron
mis-chaining, repeat-driven alignment, and ordinary novel isoforms before admission.

The interpretation has hard limits:

* an unexpressed extra copy is invisible to Iso-Seq;
* an expressed copy identical across observed bases is indistinguishable from increased expression;
* two coherent haplotypes at one represented locus are compatible with ordinary diploid alleles;
* three or more clean linked haplotypes can reject a simple one-locus diploid explanation under the
  stated assumptions, but do not locate or fully prove an extra genomic copy; and
* RNA abundance is not genomic copy number.

Consequently this avenue can produce ranked, auditable candidates and a transcript-haplotype lower
bound. It cannot by itself support the sentence "this individual has an additional genomic copy."

### 8.3 Independent-validation avenue: needed only for the stronger natural claim

If the thesis must confirm that a natural candidate is an additional genomic copy, an independent
DNA measurement is unavoidable. Suitable validation includes targeted genomic long reads, matched WGS
depth plus paralog-specific k-mers/PSVs, ddPCR, or a donor-matched phased assembly. These data may remain
**blinded evaluator-only evidence**: Rustle can make its prediction from the original reference and
Iso-Seq, after which the DNA assay tests it. This validates the method without turning the DNA result
into an input feature.

A second assembly is therefore optional validation, not the proposed solution. Feeding a better or
donor-specific genome to Rustle would largely remove the missing-reference challenge; it should not be
the main O3 experiment. Liftoff is neither required nor proposed. If projection is ever needed solely
to build evaluation truth, a minimap2-based, evaluator-side comparison is sufficient and must remain
separate from candidate discovery.

### 8.4 Claims supported by each avenue

| avenue | defensible conclusion | conclusion it does **not** support |
|---|---|---|
| paired intact/masked reference | Rustle detects the expressed consequence of a known omitted copy, at a measured sensitivity and FPR; optional reconstruction is scored separately | the unmodified biological sample naturally contains a reference-absent copy |
| unmodified Iso-Seq only | a transcript path/haplotype is not adequately explained by represented copies | the path is an additional genomic locus or a genomic copy-number increase |
| independent donor DNA validation | the flagged transcript is supported by an additional/different genomic copy | that RNA alone established the copy |
| alternate assembly used as method input | annotation/assignment against that assembly | robustness to a missing copy in the analysis reference |

**Recommended thesis order:** make §8.1 the primary O3 result; use §8.2 as exploratory candidate
generation; attempt §8.3 only if suitable donor DNA becomes available. A clean negative in natural data
and an explicit abstention region are valid outcomes. Do not claim natural reference bias unless §8.3
supplies evidence for it.

## 9. What follows

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

---

## ⛔ THE HAPLOTYPE ROUTE — NO-GO (2026-08-23)

*Can the non-chosen haplotype supply O3 a REAL positive panel of reference-absent copies?* 17 agents,
4 angles, 3 adversarial lenses each. **The premise half-holds; the panel does not survive one flag.**

### What is TRUE (verified at byte level)
`_pri` **is** a per-chromosome mosaic, verbatim: **16 PAT contigs (2,382,898,053 bp) + 9 MAT
(1,162,936,171 bp) + chrM (16,412 bp) = 3,545,850,636 bp, residual 0.** A 100 kb slice of `_pri`
NC_073229.2 is md5-identical to mat CM054587.2 at the same coordinates; the pat homolog is not.

### ⚠⚠ WHAT IS FALSE — AND IT IS THE SUBSTRATE, NOT THE IDEA
**`GGO_ds.bam` — the BAM the entire O1 catalog was built from — is a DIFFERENT ANIMAL from the
assembly.** Its `@PG` reads `raw/testis/long_reads/GGO/GGO_OR6737.IS...fastq`; 3,000/3,000 sampled
primaries carry movie `m64076_221110_210557`. Assembly = **KB3781 "Jim"** (SAMN04003007), RNA =
**OR6737**. Already recorded at `NEGATIVE_RESULTS_REGISTER.md:160` and `o1_read_evidence_repair.md:51`.
⭐ **Same-animal RNA DOES exist** — `/mnt/linuxdisk/home/juanfraitu/fibroblasts/GCA_029281585.2_flnc_mm.bam`
(SRR27438212/13, SRR27178662/63 → SAMN04003007). **Any same-individual argument must use it.**

### The route fails on its own measurement, independent of provenance

| kill test | measured | rule |
|---|---|---|
| **KT-1 parameter stability** — candidate set at `-p 0.1` vs `-p 0.8` | 56 vs 60, intersection 40, **Jaccard 0.5263** | ≥0.80 required; **FAILED** — 47% turnover on one flag |
| **KT-2 corrected instrument** — gene candidates vs composition-matched control | **37/3,853 = 0.96%** vs **82/3,862 = 2.12%**, CMH p = 0.5749 | genes are **DEPLETED**; **FAILED** |
| **KT-3 mirror arm** — excess in the *chosen* haplotype (truth value: not reference-absent) | 56 vs 51, binomial **p = 0.6992** | candidate supply is **symmetric with its own negative control**; **FAILED** |

Two instrument bugs found in the 08-19 code: `count.py` opens a new "copy" on a **1 bp** gap (ACACB's two
copies abut at CM054568.2:124,626,594/595), and counting genome-wide mat-vs-pat is wrong because
**`_pri` IS mat for 9 chromosomes**, so 26.91% of "surplus copies" are in the reference.

### ⚠⚠ RETRACT the "larger compartment" argument — it was mine and it is wrong
The non-chosen haplotype is 3,360,731,991 bp, but that is not the ABSENT set. Measured absent fraction:
**0.3295% of genic span** (297,879/90,406,212 bp) and 1.2760% of random span ⟹ **10.9–42.1 Mb
genome-wide against the 2026-08-14 compartment's ~37.0 Mb**. So **~1.1× at best, and 3.4× SMALLER in
the gene frame** — not the 84× I claimed. **The lever the 08-14 work named is NOT supplied.**

Expected detections after every repair: **≈ 0.9** — numerically identical to the 0.47–0.94 that made
the 08-14 collapse run underpowered by construction. Panel after the tandem/flank filter: **1 locus**
(TRANK1, 3,602 bp). Inserted sequence is **median 71.2% softmask repeat**; the two largest blocks are
one satellite array counted twice.

### The het-CNV objection, resolved — concede it, do not argue it
**A haplotype-specific copy is 100% in the ABSORBED stratum by construction**: its reads have a
near-perfect home in its own retained allele. It satisfies the letter of O3-restated but is a different
phenomenon, and it is present at 1× dosage in this individual and absent in ~half the population.

### What would have counted as confirmation (kept, for any future attempt)
Raw HiFi from `gorilla_hifi/SRR26152581` (verified HiFi, same animal, on disk) spanning **both**
breakpoints at 0.35–0.70× local diploid depth with zero such reads over the `_pri` interval — the only
thing separating "extra copy on haplotype A" from "collapse on haplotype B". Plus held-out molecules
carrying **≥2 linked PSVs**, scored against the third-allele null **51/19,365 = 0.0026** (the tightest
null in the design, 180× tighter than an allele-fraction null).


---

<!-- merged from o3_missing_copy_evidence.md on 2026-08-25 -->

### O3 — the unmapped-read route: is there ANY way to find a copy missing from the genome?

**Status 2026-08-15.** Companion to [`o3_missing_copy_evidence.md`](o3_missing_copy_evidence.md). That
document ends on a mechanism: on a whole genome the clipping statistic is structurally silent, so the
transcript-side detector could not fire and its 0/915 says nothing about the phenomenon. This document
asks the follow-up question directly — **is there any transcript-side route left, and can the reads that
fail to align at all be that route?** — runs it end to end, attacks it three ways, and states what it
bounds.

---

### READ THIS FIRST — the limits, before any number

1. **The single headline number is retracted.** "At most ~9 expressed reference-absent copies genome-wide
   (M ≤ 8.54)" is **not defensible as stated**. Two independent attacks broke the calibration that
   produced it, using the project's own data. Do not quote it.
2. **What replaces it is a two-stratum statement, and one of the two strata is a measured zero of power.**
   * **Unique / non-collapsible sequence** — the route works, and works better than the retracted
     headline claimed: **M ≤ 6.4 expressed reference-absent copies** (≤ 7.2 at the CI-conservative
     detection power).
   * **Collapsible sequence** — a copy that aligns ≥ 98% over ≥ 50% of its span to a surviving paralog,
     i.e. *every copy an assembler could actually collapse*: detection power measured at
     **1/35 = 0.0286 [0.0007, 0.1492]**, and **0/26 [0, 0.1323]** at coverage ≥ 0.8. **The route bounds
     this class not at all.** It is vacuous there, not merely loose.
3. **The class in (2b) is exactly O3's target.** O3 asks for a missing *copy of a family*. A copy goes
   missing from an assembly *because* it was collapsed with a near-identical paralog — that is the
   definition of collapse. So the route has high power precisely where absences are rare, and zero power
   precisely where the literature (Yoo/Rhie: 1–2 Mbp collapse per haplotype) says absences are expected.
4. **One real detection exists, and it is not a family copy.** STON1 (+ GTF2A1L), ~116.7 kb absent from
   mGorGor1 NC_073236.2, present in chimp and orangutan, 125 unmapped reads. **STON1 is single-copy.**
   It validates the O3 premise — an expressed gene absent from the assembly *is* recoverable from the
   unmapped pile — but it does not deliver a missing copy of a family, and it was scored in the
   π ≈ 0.74 stratum, which licenses no extrapolation to the π ≈ 0 stratum.
5. **Verdict (§7): VIABLE-BUT-CONDITIONAL for unique sequence, DEAD for collapsed paralogous copies.**

---

### 1. The question, and why the previous detectors could not answer it

The transcript-side detector (FARDIV ∧ FARCLIP) fired on **0 of 915** loci. §2 of
`o3_missing_copy_evidence.md` establishes that this is a property of the statistic, not of the biology:
**FARCLIP median 0.0006** (matched arm) against a gate of **0.05** — three orders of magnitude below
threshold, on a whole genome, so the conjunction is structurally unable to fire.

The mechanism is the one FixItFelix (Behera 2023) measures on the DNA side and we reproduced at 1.75×:
on 3.6 Gb, a read that does not fit locus *u* does not get forced onto *u* with a long soft clip — it
finds a **better primary home elsewhere**. Whole-genome alignment prevents the orphan pile-up that every
mini-reference excision positive exercised. That is why the mini-reference did not merely flatter the
rate; it manufactured the signal.

That leaves exactly one transcript-side stratum that whole-genome alignment cannot absorb: reads with
**no acceptable primary home anywhere** — SAM flag 4. If a copy is absent from the assembly and its
reads cannot be absorbed by a paralog, they must land there. This document tests whether that stratum
carries information, and how much.

**Scope note.** Only the flag-4 stratum was examined. The soft-clip stratum was never examined
(§4, exclusion iii) and remains open.

---

### 2. What the excision calibration established: a missing copy has two fates

The calibration masks one copy of an intact two-copy family in `GGO.fasta`, re-aligns with the BAM's own
settings, and records where the reads of the masked copy go. **n = 162 families.**

| quantity | value | 95% CI |
|---|---|---|
| P(orphan) — reads become unmapped | **54/162 = 0.3333** | [0.2613, 0.4116] |
| P(absorbed) — reads migrate to a surviving locus | 108/162 = 0.6667 | — |
| P(≥ 1 unmapped ≥ 500 bp read) | 117/162 = 0.7222 | [0.6465, 0.7896] |
| **P(≥ 8 unmapped ≥ 500 bp reads) = detection power ceiling** | **95/162 = 0.5864** | [0.5065, 0.6631] |
| families yielding **zero** unmapped reads | 45/162 = 0.2778 | — |

The dichotomy understates the route: **60 of the 104 nominally "ABSORBED" families still leak unmapped
reads.** Detection does not require orphaning, only leakage above the floor.

**Leakage is linear in expression.** On the 54 orphaned families, 61,398/61,512 = **99.81%** of orphan
reads are ≥ 500 bp (the length filter is free). log-log slope **0.9870** ⟹ linear:

* orphaned: `U = 0.7345 · X` (through origin); `U = 0.7731·X − 160.9`, R² = 0.8701
* unconditional over all 162 fates: `U = 0.2543 · X`

U quantiles at fixed class span two orders of magnitude: q0/q25/q50/q75/q100 = **5 / 120 / 357.5 / 822 /
7,103** (q10 37.7, q90 3,332.9).

**Detection function π(X) = P(U ≥ 8 | n_clean = X)** — saturates early and low:

| X (n_clean) | 10 | 25 | 100 | 431 (catalog median) | 1,000 | 5,000 | ∞ |
|---|---|---|---|---|---|---|---|
| π(X) | 0.0617 | 0.3086 | 0.4568 | **0.5494** | 0.6049 | 0.6852 | 0.5864 |

π ≥ 0.50 at X = 214. Buying 25× more expression (200 → 5,000) adds only **+0.19** to detection
probability. **Abundance is not the limiting variable; where the reads GO is.** This is the same
conclusion the divergence detector reached, by an independent route.

⚠ **Axis trap.** `n_clean` (catalog axis) and the primary-BAM read count over the same interval differ by
median factor **0.682** (r = 0.9796). Mixing them inflates every rate ~1.5×. **Every number in this
document is on the `n_clean` axis.** (An on-disk scratch script, `o3_unmapped/final.py`, uses the BAM
axis and yields π̄ = 0.5858 / M ≤ 8.10; that variant is wrong and was never shipped.)

⚠⚠ **This whole section is the part that the attacks broke.** The 162 families are *not* a random sample
of the sequence contexts in which assembly absence occurs. See §5.

---

### 3. What is actually in the unmapped pile

**Setup, re-verified in a single BAM pass:** 13,516,443 primary mapped + **959 flag-4**. Unmapped read
lengths: median 69 bp; **199 are ≥ 500 bp** (20.75% [18.23, 23.46] of 959), median 3,013 bp,
Σ = 634,416 bp. `GGO.fasta` = 3,545,850,636 bp vs NCBI's 3,545,834,224 bp over 25 scaffolds ⟹ **no
unplaced scaffolds**, so "absent from `GGO.fasta`" = absent from the whole assembly. Annotation is the
matching RefSeq `GCF_029281585.2-RS_2024_02` (isolate KB3781, male).

#### 3.1 It does not dissolve under a permissive preset

| preset | reads rescued of 199 | quality of the rescue |
|---|---|---|
| `map-hifi` | 0/199 | — |
| `asm20` | 0/199 | — |
| sensitive (`-k11 -w5 -p0.05`) | 75/199 | partial anchors |
| `-x sr` | 126/199 | **spurious**: aligned-fraction-of-read median **0.056** (q25 0.050, q75 0.080, max 0.331), **0/126 reach 50%**, median MAPQ 1 |

Hardest test: the 113-read core cluster cut into **2,232 × 400 bp tiles** → **0/2,232 tiles hit anywhere
in the gorilla genome.** The pile is not an alignment-preset artefact.

**Aligner verified** (the BAM's own command): `minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1
--secondary=yes`. Splice-aware, `-Y`, secondaries on ⟹ flag-4 is a genuine failure. The excision
calibration used byte-identical settings.

#### 3.2 It is FLNC-shaped, not junk

Against 800 mapped control reads ≥ 500 bp: mean-QV median **39.5 vs 39.6**; polyA absent in both
(refine-trimmed — the control proves absence is uninformative); **0/199 primer remnants**, 0/56 SMRTbell
adapter, 0 palindromic concatemers; 3-mer entropy median **5.79** of 6.0.

#### 3.3 It is not contamination and not run-exclusive

χ² GOF against the mapped run composition (6.56 / 12.54 / 77.36 / 3.54%):

* all 199: χ² = 15.1, df = 3, **p = 0.0017** — the pile *as a whole* is a run-quality stratum (per-run
  flag-4 rates 4.41 / 6.90 / 18.29 / 17.47 per 100k, max/min = 4.15×; whole 959-read pool
  χ² = 540.5, df = 3, p = 7.8e-117)
* **the 113-read core cluster: χ² = 4.74, p = 0.192, present in 4/4 runs** — it is the part that *tracks*
  the mapped composition. All 125 STON1 reads: χ² = 7.62, p = 0.055. Rate 125/13,516,443 =
  **0.92/100k [0.77, 1.10]**. ⟹ **index hopping refuted.**
* Contrast, this is what contamination looks like: the 33 EBV reads are **33/33 in SRR27438212,
  p = 2.1e-4**.

Species placement of the core cluster: identity **human 0.9806 / chimp 0.9812 / orangutan 0.9670** — the
gorilla outgroup position exactly. Contamination would give ~0.999 to human. NM decomposition:
substitution 1.53%, indel 0.33%, 0/125 with indel > 5%. **Human contamination refuted.**

**IG/TR = ZERO.** 0/133 human-aligned reads overlap any of 611 merged CHM13 IG/TR loci. Independently the
library is unambiguously **fibroblast** — COL1A2 310,210 / DCN 384,281 / COL1A1 183,670 / FN1 73,272
reads vs CD19 5 / MS4A1 13 / PAX5 2 / CD79A 6. No B cells ⟹ no V(D)J/SHM. (This is the trap that
retracted both earlier O3 candidates; it is closed here by cell type, not by argument.)

#### 3.4 Accounting of the 199

| class | n | share |
|---|---|---|
| **STON1** | **125** | **62.8% [55.7, 69.5]** |
| EBV (Human gammaherpesvirus 4 B95-8, identity 0.9991, all SRR27438212) | 33 | 16.6% |
| insertion-corrupted gorilla mRNAs | 23 | 11.6% |
| misc (Pongo mito, Pan, Symphalangus, human clone — all SRR27438212) | 18 | 9.0% |

The 23 BLAST to the **correct gorilla RefSeq mRNA** with 12–17% gaps and 0–2 mismatches (TUBA1A:
1660/1988 identities, 327 gaps ⟹ 1 mismatch) — genuine transcripts carrying ~15% spurious insertions
that break minimap2 seed chains (homopolymer runs ≥ 6: 1.26% vs 0.32% control; compression 0.608 vs
0.703).

**Residual after STON1 + EBV: 41 reads = one cluster of 4 + 37 singletons ⟹ ZERO additional candidate
loci with ≥ 3 reads.** The pile is essentially one locus.

#### 3.5 The locus

Gorilla NC_073236.2: LHCGR ends 91,170,740 → PPP1R21 starts 91,215,577 = **44,837 bp**. The same human
flanks span **161,543 bp** ⟹ a **~116.7 kb deficit containing STON1 + GTF2A1L**; both genes have
**0 lines** in the gorilla RefSeq GFF. The chromosome is **gapless** (0 N over 145,906,006 bp). STON1 is
**present in chimp mPanTro3** (NC_072410.2, 74,077,984–74,143,261) and **orangutan mPonPyg2**
(NC_072385.2, 85,089,179–85,154,643) — gorilla is the only one missing it.

125 reads align to CHM13 chr2 STON1 (48,536,051–48,601,303) near-full-length (aligned-fraction median
**0.999**, MAPQ 60), sharing one intron chain (junctions supported by **89, 71, 39** reads), with 5′/3′
ends matching the annotated TSS/end; longest read 6,131 bp = one full-length mRNA. Alignment geometry
kills the junk hypothesis outright: **longest contiguous CIGAR M-block median 1,230 bp** (q25 750, max
2,170) — nothing in 3.55 Gb hides from a 1.2 kb contiguous block under `-x sr -k11`. The STON1 reads are
also the *cleanest* thing in the pile (0.84 insertion events/kb, 0.26% inserted bases, vs 1.62/kb and
0.73% for the other human-hitting unmapped reads).

Expression context in the same BAM: PPP1R21 592, FOXN2 230, TTC7A 224 primary reads vs **STON1 125**;
11 reads fall inside the 44.8 kb gorilla gap. Against the excision calibration (median orphaned copy
≈ 549, q25 315, min 47), 125 sits between the minimum and q25.

**Interpretation.** A biological deletion is excluded by the reads themselves ⟹ this is an **assembly
absence, most likely a false join**, unconfirmed against raw mGorGor1 reads. And ⚠ **STON1 is
single-copy — genome-absence, not paralog-absence.**

---

### 4. THE BOUND — as it was stated, with its conditions

*This section records the bound as shipped. §5 breaks it. Read both.*

#### 4.1 The floor, and why it is not the binding constraint

**Floor: 8 unmapped reads ≥ 500 bp coalescing into one cluster.** Justified empirically, not chosen:
all-vs-all clustering of the 199-read pile gives **53 clusters — 120, 7, 6, 5, 5, 4, 3, 2, 2, 2, + 43
singletons**. Largest non-candidate cluster = **7** ⟹ a "≥ 8 in one cluster" rule fires on **0/52**
background clusters.

Converted to the catalog axis through the measured linear relation: **~12 `n_clean` if the copy orphans,
~31.5 (≈32) in expectation over all fates.** Adopt **32**. (The alternative whole-budget rule — a copy
must alone exceed the observed 199, or the Poisson-95 223.8 — gives 297–334 orphaned / 782–880 all-fates,
25× more conservative, and defensible only by refusing to cluster the background, which the 53-cluster
decomposition shows is the wrong model of it.)

**Coverage above the floor: 909/915 = 0.9934 [0.9858, 0.9976]** of catalog copies exceed 32 `n_clean`;
914/915 = 0.9989 [0.9939, 1.0000] exceed 12. Only **6 of 915** loci fall below. ⟹ **expression is not the
binding constraint; the 27.8% absorption sink is.**

#### 4.2 The bound as shipped

> M ≤ **8.5** expressed reference-absent copies. Poisson-95 upper U(1) = **4.744** on **D = 1** observed
> detection, divided by π̄ = **0.5553** (detection power integrated over the 915-locus `n_clean`
> distribution). Point estimate 1.80, 95% CI [0.09, 8.54].

Per-expression form: ≤ ~9 copies above 431 `n_clean`, ~8 above 1,000, ~7 above 5,000. **D = 1** counts
the STON1 cluster as one genuine absence; if it were a contaminant the bound would *tighten* to M ≤ 5.4,
so D = 1 is the conservative choice. If STON1 + GTF2A1L are two independent absences, D = 2,
U(2) = 6.296, **M ≤ 11.3**.

#### 4.3 What it explicitly did not cover

1. Copies **absorbed** onto a surviving paralog that leak fewer than 8 unmapped reads. π̄ < 1 corrects for
   these **only under the assumption that the excision panel's fate distribution transfers** to genuinely
   reference-absent copies. ← *this is the assumption that failed; §5.*
2. Copies expressed **below ~32 `n_clean`** (6/915 of catalog).
3. The **soft-clip stratum**, never examined — only flag-4 records were.

The 123,230-read SRA-minus-BAM gap was declared *not* an exclusion; §6 revises the arithmetic of that
claim but not its direction.

---

### 5. The attacks, and how each fared

Three independent attacks. **Two succeeded** (one fully, one partially), and they succeeded on the same
point by two different measurements.

#### 5.1 "The 199 reads are junk / zero information" — **REFUTED**

Three ways, all in §3: alignment geometry (1,230 bp median contiguous M-block, MAPQ 60, identity 0.9805,
two introns); the corruption confound is *absent* (STON1 reads are the cleanest in the pile); and the
junk stratum, which genuinely exists (χ² = 540.5 over the 959-read pool), **does not contain STON1**
(p = 0.192, 4/4 runs). Simulation confirms the run-composition screen is properly sized (pass rate
0.941–0.953 at n = 8…200 under the mapped-composition null), so passing it is informative. The axis trap
was checked and cleared: every shipped number reproduces exactly on the `n_clean` axis.

#### 5.2 ⚠⚠ "π̄ = 0.5553 does not transfer to reference-absent copies" — **SUCCEEDS. This breaks the headline.**

Two attackers reached this independently, by different measurements, using the project's own collected
data. Both stratified the 162-family calibration panel by **whether the excised copy has a
sequence-similar survivor left standing in the genome** — and detection power collapses in the stratum
where it does.

**Attack A — stratify by DNA homology to the retained partner**
(`minimap2 -c -x asm20 -k15 -w5 -N100 -p0.02`, plus a maximally sensitive `-x sr -k11` pass):

| calibration stratum | n | P(U ≥ 8) | π̄ | M ≤ (D = 1) |
|---|---|---|---|---|
| SHIPPED: all 162 | 162 | 0.5864 | 0.5553 | 8.5 |
| **no surviving homolog** | 105 | 0.7905 [0.700, 0.864] | 0.7624 | **6.2** |
| **homologous paralog survives** | 57 | 0.2105 [0.114, 0.339] | 0.1738 | **27.3** |
| … and de ≤ 5% | 29 | 0.1034 | 0.0788 | 60 |
| … and de ≤ 2% | 22 | 0.0909 | 0.0939 | 50 |
| **… and de ≤ 1% (collapse regime)** | 12 | **0.0000 [0.000, 0.265]** | 0.0099 | **478** |

Fisher on the homology split: **p = 8.0e-13**. Within the 57, detection is monotone in divergence in the
mechanistically expected direction (Spearman ρ = +0.254, p = 0.057). For the 105 with no surviving
homolog, a maximally sensitive pass finds merged homology of **median 416 bp inside ~30 kb intervals**
(q75 680, max 2,325) — excising them models **deleting unique sequence, not collapsing a copy**, and they
trivially orphan.

**The project's own flag reproduces it with no new alignment.** `per_family2.json:is_sister` —
True 14/52 = 0.2692 [0.156, 0.410] vs False 81/110 = 0.7364 [0.644, 0.816], **Fisher p = 2.3e-8**,
π̄ = 0.2292 ⟹ M ≤ 20.7. The stratification is not an artefact of anyone's alignment settings.

**The single decisive number.** The 12 families at de ≤ 1% carry **15,314 `n_clean` reads** (median 478 —
*above* the catalog median 431; GWFAM374 alone has 9,689 = 22× the catalog median). Deleting all 12
produced **exactly 1** unmapped ≥ 500 bp read. The shipped all-fates model `U = 0.2543·X` predicts
**3,894**. P(≤ 1 | 3,894) = 0. All 12 are classed ABSORBED with `unaln` 0.0000–0.0044. **A 3,900-fold
overestimate of leakage, in exactly the regime that causes assembly absence.**

**Attack B — stratify by pair coverage** (same panel, independent alignment):

| stratum | n | P(orphan) | π = P(U ≥ 8) | leak ΣU/ΣX |
|---|---|---|---|---|
| all 162 (what the bound uses) | 162 | 54/162 = 0.3333 | **95/162 = 0.5864** | 0.3446 |
| non-collapsible (cov < 0.5) | 127 | 54/127 = 0.4252 | **94/127 = 0.7402 [0.6549, 0.8139]** | 0.4024 |
| **collapsible (cov ≥ 0.5)** | 35 | **0/35 = 0.0000 [0, 0.100]** | **1/35 = 0.0286 [0.0007, 0.1492]** | 0.00272 |
| collapsible, ≥ 98% id | 18 | 0/18 | 1/18 = 0.0556 | 0.0041 |
| **collapsible (cov ≥ 0.8)** | 26 | 0/26 | **0/26 = 0.0000 [0, 0.1323]** | 0.000165 |

Fisher on U ≥ 8 × alignable: **p = 2.6e-13**. **Not an expression confound** — X distributions are
indistinguishable (Mann-Whitney p = 0.600), and **matched at X ≥ 431** the split is 1/16 vs 46/57,
**Fisher p = 7.0e-8**, per-read leak 0.00299 vs 0.398 (**133×**). Corroborated by a field the analysis
already had: 29/35 collapsible families' reads migrate to their own sister vs 23/127 non-collapsible
(p = 2.0e-12) — the collapsible stratum is **absorbed, by construction**. Extreme case GWFAM374: 60,360 bp
masked copy, **99.76% identity over 100% of its span** to a survivor 123 kb away, **7,120 reads,
0 unmapped**.

Attack B also shows **the panel is 78% non-duplications**: 66/162 pairs have *zero* alignment between
copies even at high sensitivity, only 35/162 align over ≥ 50% of the masked copy's span (median pair
coverage 0.009). These are RNA-defined (E_r) families whose members are not sequence-similar paralogs.

**Why this breaks the bound rather than loosening it.** The estimator M = U(D)/π̄ requires
π̄ = E[π | the copy is *genuinely reference-absent*]. The bound substitutes E[π | the copy is a *random
catalog copy*]. Solving the mixture backwards: **π̄ = 0.5553 is exactly the value obtained by assuming
the collapse-type share of absent copies equals the panel base rate (f = 0.260 vs base rate
35/162 = 0.216)** — i.e. by assuming **assembly absence is independent of sequence similarity**. That is
the load-bearing assumption; it is unstated, and it is false in the direction that inflates M:

`M ≤ 4.7439 / [(1−f)·0.7402 + f·0.0286]`

| f (collapse-type share of absent copies) | 0 | **0.26 (implied)** | 0.50 | 0.75 | 0.90 | 1.0 |
|---|---|---|---|---|---|---|
| M ≤ | 6.4 | **8.5 (published)** | 12.3 | 23.0 | 47.6 | **166.0** |
| M ≤ with π_collapsible at its CI upper 0.1492 | 6.4 | 8.1 | 10.7 | 16.0 | 22.8 | 31.8 |

And **masking is the best case for detection within the collapsible stratum**: masking leaves a target at
the full paralog divergence *d*, whereas a real collapse leaves a consensus/mosaic contig assembled from
the reads of *both* copies, at divergence ≤ *d* from the missing copy (0 if the assembler emitted the lost
haplotype). So π_collapsible ≤ 0.0286 is monotonically an **over**-estimate.

**⚠ Bonus: the remedy named in the bound's own `unsafe_if` would have falsely reassured.** That remedy was
"re-run the calibration deleting the copy *and* collapsing its flanks (a simulated false join) rather than
masking." It changes nothing measurable: reads of a masked copy lie inside its clean interval by
construction and never touch the flank junction; the available sequence (paralog only, at full divergence
*d*) is identical under N-mask and under deletion-with-join. That check would have returned ~the same
33.3/64.2 split and been read as confirmation. **The variable that matters is not the boundary — it is
what sequence remains standing in place of the missing copy**, which is only testable by making the
survivor a *consensus* of the two copies. No probe did that. **This is the experiment to run next.**

**Two smaller, independent hits from the same attack:**

* **Calibration sampling error was never propagated.** Bootstrapping the 162 families: π̄ 95% CI
  [0.4866, 0.6234] ⟹ **M ≤ 9.75** even on the shipped panel; double bootstrap (families + catalog)
  [0.4848, 0.6205] ⟹ M ≤ 9.79. "M ≤ 8.5" was a point estimate presented as a bound.
* **Expression-distribution sensitivity:** absent copies drawn from the catalog's lower quartile ⟹
  π̄ = 0.452, M ≤ 10.5; lower decile ⟹ π̄ = 0.414, M ≤ 11.5. The run-composition screen costs a further
  ~5% of power (its nominal α).

**Catalog representativeness is not the escape hatch.** All-vs-all of the 915 catalog intervals:
271/915 = 29.6% have a same-family DNA alignment (panel: 57/162 = 35.2%, comparable) and
64/915 = **7.0%** have a partner at de ≤ 1%. The panel *does* mirror the catalog's marginal composition —
which defends π̄ = 0.5553 as an average over *catalog-as-composed*, but **not** as detection power for the
*reference-absent subpopulation*, which the collapse mechanism selects into the de ≤ 1% tail.

**The expression-floor claim inverts with the stratification.** Converting the 8-read floor through each
stratum's own leak rate: **20 `n_clean` if non-collapsible, 2,941 if collapsible (cov ≥ 0.5), 48,421 if
cov ≥ 0.8.** Catalog coverage above those floors: **909/915 = 99.34% → 116/915 = 12.68% [0.1059, 0.1501]
→ 3/915 = 0.33% [0.0007, 0.0096].** "Expression is not the binding constraint" holds only for the stratum
where absences do not happen.

#### 5.3 "The 123,230-read gap hides the answer" — **premise REFUTED, arithmetic PARTIALLY BROKEN**

Attacked by download rather than inference: 300 MB byte-ranges of `SRR27178662/63_subreads.fastq.gz`,
split against the present-spot list, and mapped with the BAM's own command alongside a matched
present-read control. Sample **15,254 of 123,230 = 12.38%**, median length 2,933/2,610 bp, 98.6% ≥ 500 bp.

| | discarded (missing) | control (present, same untrimmed form) |
|---|---|---|
| n | 15,254 | 15,200 |
| map to `GGO.fasta` | **15,245 = 99.941% [99.888, 99.973]** | 15,198 |
| unmapped | 9 = 5.90e-4 [2.70e-4, 1.12e-3] | 2 = 1.32e-4 |
| aligned fraction, median | 0.9856 | 0.9762 |

The attack needs the missing reads to be ~100% unmapped; they are **0.059%** unmapped. The mapping is
conservative — reads were mapped **untrimmed** (primer + polyA on), strictly harder than what the BAM
saw. And **all 9 unmapped discarded reads are the artefact class `refine` exists to delete**: 9/9 carry
the IsoSeq primer `AAGCAGTGGTATCAACGCAGAGTAC` (0/199 of the real pile do), GC 0.016–0.28, 3-mer entropy
1.38–4.54 (real pile median 5.79), longest informative non-polyA/T block 17–174 bp; after the trimming
they would have received, **0/9 reach 500 bp**. **Candidate-grade additions: 0/15,254.**

*Strengthened, not weakened:* the GAP probe's claim that `isoseq refine` is genome-independent was an
assertion; it is now **measured** — discarded reads map at 99.94% with aligned-fraction median 0.9856 vs
0.9762 for retained reads. The discard is orthogonal to genomic origin.

*What the attack did break* — see §6.

---

### 6. ⚠ The gap probe's arithmetic, corrected — and what remains unexamined

#### 6.1 The gap itself reconciles exactly

| run | ENA read_count | BAM records | missing | miss % (Wilson 95%) | flag-4 |
|---|---|---|---|---|---|
| SRR27438212 | 10,457,208 | 10,457,208 | **0** | 0% (≤ 0.000029%, rule of 3) | 461 |
| SRR27438213 | 478,254 | 478,254 | **0** | 0% (≤ 0.000627%) | 33 |
| SRR27178663 | 1,770,087 | 1,694,745 | 75,342 | 4.2564% [4.2268, 4.2862] | 310 |
| SRR27178662 | 935,083 | 887,195 | 47,888 | 5.1213% [5.0768, 5.1661] | 155 |
| **total** | **13,640,632** | **13,517,402** | **123,230** | **0.9034% [0.8984, 0.9084]** | **959** |

75,342 + 47,888 = 123,230 exactly; 13,517,402 − 959 = 13,516,443 primary mapped.

**Cause:** SRA holds a **pre-`isoseq refine`** stage for exactly those two runs. `run_alias` correlates
perfectly: `Jim_GGO_MAS-1.isoseq.bam` / `Jim_GGO_IS-1.isoseq.bam` (already FLNC) → 0 missing;
`m64404e_*.hifi_reads.bc*.bam` (post-`lima`, pre-`refine`) → 4.26% / 5.12%, squarely inside `refine`'s
normal 3–8%. Alternatives are dead: both short runs still carry flag-4 records (no `-F 4` ran); a
BAM-creation filter cannot be run-selective (SRR27438212 lost 0/10,457,208 while carrying 461 unmapped);
max spot index = ENA `read_count` exactly in both short runs (no truncated tail); missing spots are
i.i.d.-scattered singletons (block-count observed/predicted 0.9975 and 0.9992, max block 5 and 4).

#### 6.2 What the attack broke, and the required edits

* **The projection was wrong and its CI excludes the truth.** The probe projected "22.1 new unmapped
  reads [19.4–25.2] = 0.040 copies-equivalent" from each library's *retained* unmapped rate. Expected on
  the sampled 15,254 under that model: 2.73; observed 9. **P(≥ 9 | 2.725) = 0.00204.** Corrected
  projection over the visible gap: **72.7 unmapped reads [33.3, 138.0] = 0.132 [0.061, 0.251]
  copies-equivalent** — 3.29× the claim.
* **A larger population no probe counted.** The gap exists because SRR27178662/63 were deposited
  *pre*-`refine`. By the same token SRR27438212/13 were deposited *post*-`refine`, so **their** discarded
  reads never reached SRA and are invisible. At the measured 4.26–5.12% loss that is **~527,000 reads**,
  making the never-examined population **~650,000 (≈4.5% of production), 5.3× the 123,230** the bound
  discusses. Projected: **383.5 unmapped [175.4, 727.8] = 0.70 [0.32, 1.33] copies-equivalent.** The
  bound's "adds at most ~22 reads = 0.040 copies-equivalent" understates this ~17×. *(Rate transferred
  from runs 62/63; MAS-Seq loss is unmeasurable from the deposit — an inference, flagged as such.)*
* **Soft-clip stratum (exclusion iii) is enriched in the discarded reads.** Query residue uncovered by
  primary + supplementary ≥ 500 bp: **11/7,866 (0.140%) missing vs 2/7,788 (0.026%) control, Fisher
  p = 0.0224, OR 5.45.** Over 650,000 that is ~910 such reads — but the examined BAM already carries
  ~3,500 by the control rate, so the gap is ~26% of an **already-unexamined** stratum, not a new dominant
  one.
* **Add the escape clause the bound lacks:** if any 8 of those additions coalesced into one cluster,
  **D = 2, U(2) = 6.296, M ≤ 11.3** (same magnitude as the STON1/GTF2A1L alternative). Measured
  occurrences: **zero**.
* **Add a low-complexity/primer filter to the floor rule.** The "≥ 8 in one cluster fires on 0/52
  background clusters" justification was derived on a primer-**trimmed** pile; the un-refined stratum's
  unmapped tail is 9/9 mutually similar polyA-primer artefacts that *would* coalesce into a cluster ≥ 8 if
  ever admitted. The filter is free (0/199 real reads carry a primer; 196/199 have a ≥ 500 bp
  non-polyA/T block) but it must be **stated**, not assumed.

#### 6.3 Net effect and closure cost

**Point estimate unchanged.** Candidate-grade additions measured **0/15,254**; Poisson-95 upper **24.2
reads** over the visible gap, **127.8** over the full ~650,000 — against a floor of ≥ 8 in one cluster.
The gap is **no longer the binding uncertainty**: §5.2 is.

**Cost to close what remains** (the full-depth re-run of the two pre-`refine` runs): 5.44 GB gz download,
~18 GB uncompressed, ~1–2 h wall clock on 5 cores; recovered population is by construction enriched for
concatemers/non-poly(A) — the exact chimeric-cDNA false-positive mode already burned here.
**Recommendation stands: not worth recovering.** The ~527,000 invisible reads on SRR27438212/13 cannot be
recovered at all from the deposit; only re-processing from the raw PacBio subreads would expose them.

---

### 7. Honest verdict

**The unmapped route is VIABLE-BUT-CONDITIONAL for unique sequence and DEAD for collapsed paralogous
copies.** It is not one route with one bound; it is two routes with opposite verdicts, and the wrong one
carries the O3 question.

#### 7.1 The replacement statement — quote this, not M ≤ 8.5

**Stratum 1 — unique / non-collapsible sequence** (78% of the calibration panel; STON1's own class:
single-copy, no paralog):
> π = **0.7402 [0.6549, 0.8139]**, D = 1 ⟹ **M ≤ 6.4 expressed reference-absent copies** (≤ 7.2 at the
> CI-conservative π; ≤ 8.7 if STON1 + GTF2A1L are counted as D = 2). Floor **~20 `n_clean`**;
> 909/915 = 99.34% of catalog copies clear it.

This is the part of the original bound that survives — and it is **tighter** than the retracted 8.5.

**Stratum 2 — collapsible sequence** (a copy aligning ≥ 98% over ≥ 50% of its span to a surviving
paralog, i.e. every copy an assembler could actually collapse):
> π = **1/35 = 0.0286 [0.0007, 0.1492]** at cov ≥ 0.5 and **0/26 [0, 0.1323]** at cov ≥ 0.8; D = 0 ⟹
> **M ≤ 105 at the point estimate, ≥ 20–23 even at the most defender-friendly CI upper limit, and
> formally unbounded at cov ≥ 0.8.** Implied floor **2,941–48,421 `n_clean`**, cleared by
> 12.68% [0.1059, 0.1501] → **0.33% [0.0007, 0.0096]** of the catalog.

Honest wording: **the unmapped-read route carries essentially no information about collapsed copies.**

**One-sentence quotable form:**
> *"≤ 6 reference-absent copies of unique sequence expressed above ~20 `n_clean`; this route bounds
> collapsed paralogous copies not at all (π = 0/26 measured, 95% CI upper 0.132 ⟹ M ≤ 23 at best, and
> only 0.33% of the catalog clears the implied floor)."*

#### 7.2 Why this is still a result worth writing

* **The premise is confirmed, once, on real data.** An expressed gene genuinely absent from mGorGor1
  (STON1 + GTF2A1L, ~116.7 kb, present in chimp and orangutan, 0 lines in the GFF) **was** recovered from
  the flag-4 stratum, survived every contamination, run-artefact, IG/TR and preset-dissolution check, and
  is not junk by any measurement applied. The transcriptome *can* see an assembly absence.
* **The negative is bounded where it can be bounded** (stratum 1: ≤ 6) and **explicitly unbounded where it
  cannot** (stratum 2: zero power). A bounded negative plus an admitted zero-power stratum is publishable;
  an unbounded negative is not, and a *falsely* bounded one is a retraction waiting to happen.
* **Three routes now converge on the same mechanism, independently.** Clipping is silent because orphans
  are absorbed (FARCLIP 0.0006 vs a 0.05 gate). Divergence saturates because abundance is not the limiting
  variable (π gains only +0.19 for 25× expression). The unmapped route dies in the collapsible stratum
  because the reads migrate to the survivor (29/35 vs 23/127, p = 2.0e-12; per-read leak 133× lower).
  **All three failures are the same fact: where the reads GO, not how many there are.** That is the
  finding, and it is a mechanism, not an excuse.

#### 7.3 ⚠ How this relates to the DNA-side zero — do not merge them

The DNA-side probe (assembly-vs-assembly, the field's S1 standard, three assemblies of one gorilla)
returned **0 collapses at 816/817 loci** with the instrument validated at an FP floor of 0/817, and that
zero was **predicted**: Yoo/Rhie's 1–2 Mbp collapse per haplotype over the 1.1224% of the genome those
windows span predicts only 0.47–0.94 collapses ⟹ underpowered by construction.

⚠ **The transcript-side zero and the DNA-side zero are concordant but not independent confirmations of
each other, and neither is evidence of absence-of-absences.** The DNA side says *"an instrument with
measured sensitivity found nothing in a compartment too small to expect anything."* The RNA side says
*"in the stratum where collapses live, this instrument has no measured sensitivity at all."* Reporting
them as "two independent zeros" would be exactly the denominator error this project has retracted before.

#### 7.4 The one experiment that would change the verdict

Re-run the excision calibration with the **survivor replaced by a consensus/mosaic of the two copies**
(what an assembler actually emits when it collapses), rather than by N-masking one copy and leaving the
other at full divergence *d*. That, and only that, measures π in the regime that matters. Predicted
direction: π_collapsible falls further, from 0.0286 toward 0, and stratum 2's "bounds nothing" hardens
from a CI statement into a mechanism.

Explicitly **not** worth doing: the flank-collapse variant named in the original `unsafe_if` (§5.2 shows
it cannot move the measurement), and the SRA gap recovery (§6.3).

---

### Artefacts

* `/home/juanfra/winloci_scratch/o3_unmapped/` — `CHARACTERISE.md`, `CALIBRATE.md`, `GAP.md`,
  `fates.json`, `perrun_counts.tsv`, `unm500.fa`, `ava.paf`, `cdhit90.clstr`, `hs_rows.json`,
  `hs_splice.sam`, `PTR.sam`, `PPY.sam`, `noany_tiles_ggo.paf`, `ston_tiles_ggo.paf`, `classification.tsv`
  ⚠ `final.py` on disk uses the wrong (BAM) axis — see §2.
* `/home/juanfra/winloci_scratch/o3_attack/` — `keep.fa`, `mask.fa`, `pairs.paf`, `pairs_sens.paf`,
  `pairdiv.json`, `pairdiv_sens.json`, `joined2.json` (per-family fam / cls / X / U / pair identity /
  pair coverage / is_sister)
* `/mnt/linuxdisk/home/juanfraitu/attack/` — `pairde.json`, `pairs2.paf`, `miss_sr.paf`, `cat_ava.paf`,
  `cat_pairde.json`, `maskR.fa`, `keepR.fa`
* `/mnt/linuxdisk/home/juanfraitu/o3_gapattack/` — `missing.fq`, `missing63.fq`, `tomap.bam`,
  `tomap63.bam`, `all6_unmapped.fa`, `unm63.fa`, `present_62.txt`, `present_63.txt`, `mm2.log`, `mm63.log`
* `/mnt/linuxdisk/home/juanfraitu/spots_short_runs.tsv`


---

<!-- merged from o3_missing_copy_evidence.md on 2026-08-25 -->

### Haplotype copy-number run: UNINFORMATIVE by its own pre-registered criterion

**Status 2026-08-19. Run complete. Verdict: the primary number must not be quoted.**
Pre-registration: [`o3_missing_copy_evidence.md`](o3_missing_copy_evidence.md).
Artifacts: `/mnt/linuxdisk/home/juanfraitu/o3_hapcnv/`.

### 1. What was run

10,019 probes (5,000 random autosomal protein-coding genes as centred windows capped at 30 kb,
plus named controls and a span-matched random-interval arm), 187.1 Mb, aligned to `mat` and `pat`
each restricted to its 24 chromosome-scale sequences (98.67% of mat, 100.00% of pat).
`d_hap = copies(pat) − copies(mat)` at the validated P3 floor (identity ≥ 0.973, ≥ 5,800 aligned bp).

**Amendment, declared:** the first pass used `-N 100` at minimap2's default `-p 0.8`. That failed
pre-declared control 2, and diagnosis showed why: `MAPKBP1_full` returned **one** record (58,172 bp
at identity 1.0000) because its paralogues, covering only part of a 58 kb probe, score under 0.8× the
self-hit and are discarded, while `PLA2G4B_full` (10.4 kb probe, paralogues spanning it at 0.984–0.986)
returned nine. **Copy counts were biased by probe length.** Re-run at `-N 200 -p 0.1`; records rose
28,623 → 97,930 (mat) and 21,983 → 78,278 (pat). This was fixing the instrument to meet a
**pre-declared** must-pass, not tuning to move the primary; all controls were re-evaluated.

### 2. Controls

| # | control | result | |
|---|---|---|---|
| 3 | self-recovery — a probe from a PAT-sourced chromosome must be found in pat | **2652/2652 and 1201/1201 = 1.0000** | PASS |
| 4 | single-copy panel | TBP, POLR2A, GTF2H1, SF1, TFRC, HMBS all **1 and 1** | PASS |
| 5 | probe-provenance stratification | PAT-sourced **0.0290 [0.0233, 0.0361]** vs MAT-sourced **0.0250 [0.0176, 0.0354]** | PASS |
| 2 | known answer (pat 8/9/9, mat 5/6/8) | SPTBN5 **9/8 exact**; MAPKBP1 **9/8** (expected 8/5); PLA2G4B **9/7** (expected 9/6) | **FAIL** |
| 1 | span-matched random intervals | **0.1512 [0.1403, 0.1629]** vs the gene rate **0.0278** | **FAIL — floor is 5.4× the signal** |

### 3. Verdict

The scope declared: *"If control 1 lands above the observed rate, the correct report is
'uninformative', not a number."* It landed **5.4× above**. ⟹ **The gene rate of 2.78% is not
quotable and no rate is reported.**

**Why control 1 failed is a design flaw, not noise.** The controls were **span-matched but not
composition-matched**: a random 30 kb interval is far more repeat- and segdup-rich than a gene body,
so it both varies more biologically and is harder to count. A signal *below* its own null is not a
weak signal — it is evidence the null is the wrong null. Any repeat attempt must match repeat
content and segdup overlap, not just span.

### 4. ⚠ A prior number is now in doubt — correction

`OBJECTIVES_AND_VERIFICATION.md` row 3.10 reports MAPKBP1/PLA2G4B/SPTBN5 at **8/9/9** on `_pri`/`_pat`
versus **5/6/8** on `_mat`, and that "several-copy difference between one animal's two haplotypes"
has been used as the standing refutation of *"copy number is stable; only SNPs and indels differ."*

**This run cannot reproduce it with either setting.** At `-p 0.8` MAPKBP1 gives 1/1; at `-p 0.1` it
gives 9/8. Neither is 8/5. SPTBN5 reproduces exactly (9/8) and PLA2G4B is off by one on mat (9/7 vs
9/6). With secondaries retained the mat deficits shrink from **3, 3, 1** to about **1, 2, 0**.

⟹ **Do not quote "1–3 whole gene copies differ between KB3781's haplotypes" until row 3.10's
instrument is re-derived and its `-p`/`-N` settings recorded.** The direction (mat ≤ pat at these
loci) survives; the magnitude does not. This does not rescue the "copy number is stable" proposition
— SPTBN5 and PLA2G4B still differ — but the difference is smaller than the number in circulation.

### 5. What does survive

**Control 5 is a clean positive result and is worth keeping.** `_pri` is a per-chromosome mosaic —
independently re-derived here from exact chromosome lengths as **16 PAT / 9 MAT**, reproducing the
documented split — so probes are paternal sequence on some chromosomes and maternal on others. The
`d_hap` rate agrees across the two strata (0.0290 vs 0.0250, overlapping CIs). **`d_hap` is not
driven by probe provenance.** That was the sharpest confound available and it is now excluded, which
any future attempt inherits.

### 6. Reproduction

```bash
cd /mnt/linuxdisk/home/juanfraitu/o3_hapcnv
python3 probes.py            # probe set + span-matched controls (span multiset asserted identical)
./align.sh mat && ./align.sh pat    # ~11 + ~8 min, peak 14.7 GB, ONE AT A TIME
python3 count.py             # controls first, primary last
```

---

# Appendix A — the PRE-REGISTRATION, written before the run

*Merged from `o3_missing_copy_evidence.md` on 2026-08-20. ⚠ Committed BEFORE the run (ed86742 then
1247d67). Its declared stop rule is what makes the verdict "uninformative" a result rather than a
disappointment.*

### Scope: how often do whole gene copies differ between one individual's two haplotypes?

**Status 2026-08-19: SCOPED, NOT RUN.** Pre-registration draft — thresholds and controls below are
declared before any result.

### The question and why it is the right one

The advisor's objection to O3 is that *"the genome is an average, so every individual differs from it
unless you compare to the individual that generated it."* For our substrate the premise is wrong on
two counts: **mGorGor1 is a haplotype-resolved assembly of ONE animal (KB3781 "Jim")**, not a mosaic
of donors, and **the fibroblast IsoSeq is that animal's own cell line**. The register already scopes
the objection as *"live only for the CROSS arm and human work."*

But the counter-proposition — *"copy number is stable; only SNPs and indels differ"* — is also false,
and our own data already refutes it. On this assembly, `_pri`/`_pat` vs `_mat`:

| probe | `_pri` / `_pat` | `_mat` |
|---|---:|---:|
| MAPKBP1 | 8 | 5 |
| PLA2G4B | 9 | 6 |
| SPTBN5 | 9 | 8 |

**One animal's two haplotypes, differing by 1–3 whole gene copies.** What is missing is the **rate**.
The existing measurement — `d_ortho` nonzero on **11/816 = 1.35% [0.75, 2.40]**, 6 up vs 5 down,
sign test p = 1.0000 — carries two stated limits: it is **tandem-only by construction** (`in_pri = 1`
for 762/817 = 93.3% of loci, so a dispersed event is structurally invisible), and the stratum spans
only **1.1224% of the genome**.

This run replaces the compartment and drops the tandem restriction.

### Design

**Probes.** 5,000 randomly sampled autosomal protein-coding genes (of 22,650 in `GGO_genomic.gff`),
**genomic spans, introns included** — the validated detection floor (P3: ≥97.3% identity over ≥5.8 kb)
was established with genomic queries. Extracted from `_pri`.

**Targets.** `mat.fa` and `pat.fa`, each restricted to its **24 chromosome-scale sequences (≥20 Mb)**.
This matters: mat has 225 sequences and pat has 24, but the 24 large ones are **98.67% of mat and
100.00% of pat**, so restricting makes completeness symmetric and discards only 48 Mb = 1.33%.

**Statistic.** `d_hap(g) = copies(g in pat) − copies(g in mat)`, counted genome-wide (not
window-restricted) at the P3 floors.

**Primary readout.** Fraction of genes with `d_hap ≠ 0`, with CI, plus the **sign test** — polymorphism
predicts symmetry, a systematic assembly deficit predicts skew.

### Controls, declared before any result

| # | control | must-pass |
|---|---|---|
| 1 | **Span-matched random intervals** — same SIZE distribution, not count-matched (count-matching alone has already killed a metric here) | the FP floor must sit well below the observed rate |
| 2 | **Known answer** — MAPKBP1 / PLA2G4B / SPTBN5 | pat 8/9/9 and mat 5/6/8, or this is not the instrument that produced the published result |
| 3 | **Self-null** — pat vs pat | `d = 0` on every probe; any nonzero is a bug |
| 4 | **Single-copy panel** — TBP, POLR2A, GTF2H1, SF1, TFRC, HMBS, PSMB6 | 1 and 1 |
| 5 | ⭐ **Probe-provenance stratification** | `_pri` is a per-chromosome MOSAIC — 16 chromosomes byte-identical to `_pat`, 9 to `_mat` — so the probe is paternal sequence on some chromosomes and maternal on others. **The rate must agree across the two strata.** If it does not, `d_hap` is measuring probe provenance, not biology |

Control 5 is the sharpest and has no analogue in the previous run, which was forced through `_pri`
as one of the two compared assemblies. Comparing mat to pat **directly** removes that entanglement
from the comparison, but not from the probe — hence the stratification.

### The main risk, stated in advance

Genome-wide counting is exactly where the previous module measured a **control artefact rate of
8.9–12.8%**, which is why it declined to treat `d_other` as its collapse statistic. If the true
polymorphism rate is ~1–2%, **the control floor could exceed the effect** and the run returns
nothing interpretable. Mitigation is the P3 floor (≥97.3% identity, ≥5.8 kb), which is what drove
`d_ortho`'s own FP floor to 0/817. **If control 1 lands above the observed rate, the correct
report is "uninformative", not a number.**

Second limit: this measures **polymorphism**, not collapse. The two are different questions and one
run cannot separate them. It answers the advisor's objection; it does not by itself advance O3's
detector.

### Cost

Two minimap2 index builds (~10 min each, ~12 GB RSS) plus ~130 Mb of query against each haplotype:
**≈2–2.5 job-hours, one at a time**, peak RSS ~14 GB of the 25 GB ceiling.

### The decision it produces

| outcome | reading |
|---|---|
| rate low (~1–2%), symmetric, control floor well below | the objection is quantitatively small; the matched design is sound and O3's negatives stand |
| rate high (>10%) | copy-level polymorphism is common ⟹ every O3 claim must be matched-individual — which we can do, and the cross/testis arm must be re-scoped |
| control floor ≥ observed rate | uninformative; report as such and do not quote a rate |


---

<!-- merged from o3_missing_copy_evidence.md on 2026-08-25 -->

### Over-collapse: is it happening, and could we see it?

**Status 2026-08-19.** Analysis of the original O3 premise, plus a scoped simulation that is
**specified, not run**. Companion to [`o3_missing_copy_evidence.md`](o3_missing_copy_evidence.md),
[`o3_missing_copy_evidence.md`](o3_missing_copy_evidence.md) and the O3 restatement in
[`THESIS_OBJECTIVES.md`](THESIS_OBJECTIVES.md).

### 0. The premise

*"Very similar copies may be collapsed into one in the assembly. Align RNA against that assembly and
the over-collapsed copies are missed."* That is O3's original motivation. It decomposes into two
questions with **different answers**, and conflating them is the main error to avoid.

### 1. Is it happening? — UNKNOWN, and our zero does not say no

| | |
|---|---|
| DNA-side S1, three assemblies of one gorilla | **0 collapses at 816/817** |
| transcript side | **0/915** — concordant |
| literature (Yoo/Rhie) | ~**1–2 Mbp of collapse per haplotype** |
| what that predicts over the **1.1224%** of the genome those windows span | **0.47–0.94 events** |

**The zero was predicted before it was observed.** An instrument with sub-one expected yield returns
zero whether or not collapse exists. This is underpowered *by construction*, not a negative result.

⚠⚠ **And the screen's validation now needs re-checking.** It rested on recovering a published
expansion exactly — MAPKBP1/PLA2G4B/SPTBN5 at 8/9/9. On 2026-08-19 that count proved
**`-p`-sensitive**: at minimap2's default `-p 0.8`, MAPKBP1 returns **1** copy, because its paralogues
cover only part of a 58 kb probe and score under 0.8× the perfect self-hit; at `-p 0.1` it returns
**9**. If the collapse screen ran at default `-p`, its copy counts were suppressed on long probes —
**biasing it toward zero**. See [`o3_missing_copy_evidence.md`](o3_missing_copy_evidence.md).
**Do not re-quote "0/817 with a validated instrument" until that screen's `-p`/`-N` are recorded.**

### 2. Could we see it? — the signature is NOT what O3 originally assumed

**A collapsed copy does not orphan its reads. It absorbs them.** Measured in the whole-genome excision
control: absorbed copies show **1.75× depth** on the surviving paralogue (FixItFelix reports 1.5×),
concordance 0.967. So the signature is **depth + PSVs** — the field-standard **S2** — and *not*
unmapped reads, and *not* clipping (no published collapse detector uses clipping).

#### 2a. In RNA, the depth half is dead

Depth conflates copy number with **expression level**, which varies over orders of magnitude, so 2×
coverage at a locus means nothing without a per-gene expectation. Already on the register as
human-review-flag-only. **That leaves only the PSV half.**

#### 2b. PSVs need divergence — and collapse is *caused by* the lack of it

| | TPR |
|---|---:|
| S2 detector, copies **≥ 0.01** diverged | **0.4500** |
| S2 detector, copies **< 0.01** diverged | **0.0588** |
| **fraction of true positives below 0.01** | **45.78%** |

held-out overall: **TPR 0.2703 / FPR 0.0200**.

**Copies collapse *because* they are near-identical.** The mechanism that creates the target is the
mechanism that destroys the evidence. This is **adverse selection, not bad luck**, and no threshold
choice escapes it.

At the limit it is a genuine impossibility rather than a power problem: **K = 0 — perfectly identical
copies — is unidentifiable, and that is entailed, not measured.** Between K=0 and 0.01 divergence it
is power; above 0.01 it is 0.45.

### 3. So it is not impossible — but not from RNA alone

| half | source | status |
|---|---|---|
| PSV | RNA | available, bounded, adverse-selected |
| depth | **RNA** | ⚠ **structurally dead** — expression confound |
| depth | **DNA** | ⭐ **available and unblocked** |

The gorilla WGS was confirmed usable — the pre-registered k-mer test rejected the "Y flow-sorted"
label at **17–1600×** (ENA's `sample_title` is simply wrong for every SAMN04003007 run). **The open
route is RNA PSVs + DNA depth on the matched individual, over a compartment large enough to expect
more than one event.** Both halves exist; neither has been run against the other.

### 4. Scoped simulation — a synthetic-collapse calibration ladder

**SPECIFIED, NOT RUN.** This is the calibration the collapse screen explicitly lacks: *"`d_ortho` has
never been calibrated by a synthetic collapse. The one such test on record is arithmetic, not
calibration — it proves the statistic is wired without an off-by-one, not that minimap2 plus the
clustering rule resolve two real adjacent tandem copies as two."*

#### 4a. Collapse is NOT deletion — simulate the right operation

The excision control **deleted** one copy: the survivor is copy A's sequence, and copy B's reads must
find a home. **Collapse merges**: the assembly carries **one** sequence where biology has **two**, and
reads from *both* copies land on it. The operations differ, and the existing panel only tested
deletion. **Both arms must be run on the same families** so the difference is measured, not assumed.

#### 4b. Design

Hybrid — **real reads, modified genome**. Simulating reads would forfeit the error and expression
structure that make the excision control credible.

1. Take the existing **162 two-copy family** panel and its matched IsoSeq.
2. Bin each family by its pair's **measured divergence** (0, ≤0.001, ≤0.002, ≤0.005, ≤0.01, ≤0.02, >0.02).
3. **Collapse arm:** in ONE modified genome, replace both copies of every family with a single
   sequence — as the excision arm masked all 162 in one genome, so it is **one alignment, not one per bin**.
4. **Deletion arm:** the existing excision genome, unchanged, for the contrast in §4a.
5. Realign the real matched IsoSeq; run the existing `detector.py`.
6. **Primary output: TPR as a function of divergence** — converting the two-bin summary
   (0.4500 / 0.0588) into a **located detection floor**.

Reuses `o3_excise/{mask_genome.py, align.sh, detector.py, panel.py}` and the frozen panel.

#### 4c. Controls, to declare before running

| # | control | must-pass |
|---|---|---|
| 1 | unmodified genome, same reads, same detector | the FP floor; must be far below the collapse-arm rate |
| 2 | single-copy loci | must not fire |
| 3 | deletion arm on the same families | establishes that collapse ≠ deletion rather than assuming it |
| 4 | `-p` / `-N` recorded and swept | ⚠ the 2026-08-19 lesson — a default `-p 0.8` silently discarded 8 of MAPKBP1's 9 copies |

⚠ **Known traps.** Identical copies **fake junctions** (the TSPY simulation). Simulations need `-N 50`
— `--secondary=no` yields **0 families**. And a planted-divergence experiment scored on planted
divergence is circular: the detector must be blind to the bin.

#### 4d. What it can and cannot establish

**Can:** where the detection floor actually sits, whether K=0's impossibility begins at 0.000 or bites
well above it, the FP floor under a known negative, and whether collapse and deletion really differ.

**Cannot:** ⚠ **whether over-collapse is happening in mGorGor1.** That is a *prevalence* question and
needs a larger compartment, not a calibration curve. A ladder tells you what the instrument can see;
it never tells you what is there.

#### 4e. Cost

One genome modification, one index, one IsoSeq realignment against the existing panel: **≈2–3
job-hours**, one at a time, comparable to the excision run.


---

<!-- merged from o3_missing_copy_evidence.md on 2026-08-25 -->

### Why family-similar reads can still be "unmapped"

A common worry when proposing to rescue unmapped reads against a gene-family
model: *if a read has any resemblance to known family copies, shouldn't the
aligner have found it already?* The answer is no — long-read aligners
(minimap2 in particular) use **seed-chain-extend with hard thresholds**, not
pairwise similarity. A read can carry substantial family-shared sequence and
still fail at any of several stages.

### Failure modes

| Failure mode | What happens in minimap2 | Why a k-mer / family-HMM scan still finds it |
|---|---|---|
| **1. Below-threshold chain score** | Read has matching k-mer islands separated by indels or divergence. Each island is too short to anchor a primary chain (minimap2's `-s` minimum chainable score, default ≈ 40 for HiFi). | Family scan sums *all* k-mer hits ignoring chain geometry. Fifty scattered 15-mer hits add up to a strong family signal but don't chain into a single alignment. |
| **2. Repetitive-seed masking** | Minimap2's `-f` flag drops the top X% most-frequent seeds. In a multi-copy family, those are exactly the *family-conserved* seeds — the most informative ones for family attribution. | The family k-mer index keeps every family k-mer, abundant or not. Profile-HMM emissions weight conserved columns *up*, not down. |
| **3. Divergence above the affine-gap budget** | A novel copy at ~85% identity with several copy-specific indels exceeds minimap2's score-per-base budget at the alignment-extension stage even after seeding succeeds. The candidate is dropped. | HMM forward scoring sums log-probabilities over all paths; there is no hard "give up" threshold, and indels are first-class transitions, not penalties against an affine model. |
| **4. Copy-specific structural difference** | A read spanning a copy-specific exon insertion gets fragmented into supplementaries with no primary alignment, or rejected by primary-alignment fraction filters. | A per-exon profile-HMM family graph allows the read to traverse a copy-specific branch at a bubble — the structural difference is modelled, not penalised. |
| **5. Reference is wrong (genome is an average)** | The true source copy is not in the linear reference at all (the assembly is a consensus over individuals; some paralogs exist in the donor but not the reference). Minimap2 has no positive anchor. | Cross-family conserved emissions still light up at conserved positions; the HMM lands the read on the closest paralog branch, and the divergence shows up as low per-column posterior — a *signature* of a missing copy. |

### The two questions, side by side

- **Aligner asks:** does a high-scoring chain exist between this read and the linear reference?
- **Family HMM asks:** what is the probability this sequence was emitted by *some* path through the family graph?

The second question is strictly more permissive in exactly the regimes that
matter for multi-copy families: high paralog-internal divergence,
copy-specific structural variation, and reference bias.

### Validation: per-read failure-mode classification

Rescue is only a defensible scientific claim if we can show *which* of the
five failure modes each rescued read fell into. Rustle's rescue pipeline
classifies every rescued read, distinguishing:

- **Buckets 1–2 (below-threshold / seed-masked):** the aligner could have
  found these with relaxed parameters. The HMM is acting as automatic
  re-parameterisation. Useful but not novel.
- **Buckets 3–5 (divergent / structural / reference-absent):** no minimap2
  parameter setting recovers these. The family HMM is doing irreducible
  work; these are the real novel-copy candidates.

For buckets 3–5 specifically, Rustle re-runs minimap2 as an external
subprocess with stepped-down parameters as verification — turning the
"this could not have been found by the aligner" claim into a reproducible
fact rather than a model assumption.
