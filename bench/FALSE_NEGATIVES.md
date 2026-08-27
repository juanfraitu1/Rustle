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

**⭐FIXED (commit pending): a genomic-span homology tier** now runs alongside the exon-sum core in
`refine_families_exon_sum` — a real segdup covers >=50% of its GENOMIC extent at high (gap-compressed) identity
even when its partial transcript models fail the exon-sum coverage floor, while a repeat-bridge covers <50% of the
genomic span regardless of the repeat's identity. Measured genome-wide: refined families **86 -> 100** (+14),
**8 of the 13 measured FNs recovered** (the coverage + split classes: EOLA1/2, ZNF74, ARHGAP23, PDPK1, alpha2M,
FRG1...), the structural audit stays CLEAN (0 giant-span, 0 cross-shared), and all spot-checked FP gene-splits/
bridges (PBX1, EBF1, CTNNA2, HS6ST2, NNT-GHR, GARRE1-ZNF540) stay ABSENT — no FP regression. The remaining **5 FNs
are the DIVERGENCE class** (ARMCX 65% aa, IFITM, KRAB-ZNF, RABGEF1, GRAP): below the identity floor, the genuine
precision/recall tradeoff, recoverable only with `--protein-tail` or a lower `--min-identity` (which risks FPs).

Original analysis (why the coverage class was a metric artifact): **7 of the 13 (the coverage class) are a fixable metric artifact, not a real limit** — the copies are
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

---

## ⭐ 2026-08-13 — O-4: THE TWO CALL SITES' OPPOSITE DEFAULTS, ADJUDICATED

This file was **deleted in `9b0814f`** ("docs consolidation", 104 bench files, *"git retains all"*). The
deletion repointed dangling cross-links across the bench docs but missed the one that lives in **Rust
source**, so `denovo_pipeline.rs`'s genomic-span tier — a **default-ON edge rule** — spent a month
citing a file that was not in the tree. It is restored above verbatim from `9b0814f^` (byte-identical to
the git blob). ⚠ **DELETED-WITH-A-DANGLING-CITATION, not never-existed**; an earlier audit line recording
"never committed in git history" is wrong (`git log --all -- '*FALSE_NEGATIVES*'` → `9b0814f`, `4586ba8`,
`0bb3f83`).

### The question
Two call sites ship OPPOSITE defaults for "does `E_r` see the genomic span as well as the exon-sum?"

| site | what it does | default |
|---|---|---|
| `homology_edges_all_reps_pooled` (O1 definition) | **SWAPS** substrate, globally, over the whole pooled rep graph; γ-quasi-clique on top | `homology_genomic_span` **OFF** |
| `refine_families_exon_sum` (O2 / legacy catalogs) | **UNIONS** a genomic-span tier in, **gated** on `!edges_connect_all`, per input family; connected components on top | additive tier **ON** |

### Decision: KEEP BOTH. They are not the same operation.
A **gated, family-local ADDITION** can only ever connect copies *inside* a copyset the upstream stage
already isolated — it cannot bridge two families. A **global SWAP** re-decides every pair. The gate is
the whole difference, and it is a *sufficient condition for harmlessness* on the connected case: refine
returns `Vec<Vec<DenovoTranscript>>` (copy sets, no edges), so the union's only downstream channel is
the PARTITION, and gate-true is a provable no-op.

### Refine site — additive tier stays ON (recall gain, no measured precision cost)
Denominator printed by `RUSTLE_UNION_AUDIT` before any rate: the **26** input families with ≥2 copies
that refine examined on the `--cross-chrom` human panel run.
* tier fired **and moved a partition: 11/26 = 0.4231** Wilson95 [0.2554, 0.6105]
* **block sets differing, SHIPPED(`E_x ∪ E_g`) vs DNA-only(`E_g`): 0/26 = 0.0000** [0.0000, 0.1287]
* emitted catalogs **IDENTICAL**: 26 families / 143 copies both ways, **ARI 1.0000**, forbidden pairs 0 and 0
* O2 `copy_assign`, 25 control regions, denominator **7** families examined: gate-false 3/7, 1 edge added,
  1 partition moved. **SHARP `CAFAM0`** (`NC_086017.1:207204670-207210179`, 2 copies, 339 reads, 114
  resolvable PSVs, 89 assigned reads) exists **solely** because of this leg: exon-sum 2148/2323 bp at
  identity **1.0000** but coverage **0.0480** → 0 edges; genomic span 5509/33521 bp at identity 0.9811,
  coverage **1.0000** → 1 edge. Its sibling `SHARP fam1` fired the tier and was correctly **REJECTED**
  (genomic coverage 0.1316) — the anti-bridge guard holds.
* ⚠ **"the tier moved a partition" is NOT informative alone**: a degree-sequence-matched null moves the
  partition with P = 0.66–1.00 per family. The informative statistic is WHERE it lands — and it lands on
  `E_g`'s partition **26/26**.
* This restores the original justification above: **+14 refined families genome-wide (86 → 100), 8 of 13
  measured FNs recovered, 0 giant-span / 0 cross-shared, no FP regression.**

### O1 site — `homology_genomic_span` stays OFF, on **unresolved sign**, NOT on measured harm
Denominator fixed before the run: **U = 18,528** annotation-labelled pairs among the **244** emitted reps
(193 labelled; 6,241 same-family / 12,287 different-family), from the 7-family panel roster
(`o1_closure/nodes/*.tsv`), which predates the run. **The node set is identical in both arms**, so only
the partition moves.

| | OFF (shipped) | ON | Δ |
|---|---|---|---|
| `E_r` edges | 696 | 1,575 | +879 (+126%) |
| pairwise **recall** on U | 0.1490 [0.1404, 0.1581] | 0.1780 [0.1687, 0.1877] | **+0.0290**, family-clustered P(Δ>0) = 0.9673 |
| pairwise **precision** on U | 0.8501 [0.8277, 0.8700] | 0.9876 [0.9792, 0.9926] | +0.1375 — but **P(Δ>0) 0.4070 / tie 0.1993 / P(Δ<0) 0.3937** |
| F1 | 0.2536 | 0.3017 | +0.0481, P(Δ>0) = 0.9673 |
| frozen-metric `gene_capture` (its headline recall) | 49/740 | 49/740 | **0** |
| frozen-metric `purity_strict` | 0.4140 | 0.3636 | −0.0504 |
| frozen-metric `contamination` | 0.0375 | 0.0685 | +0.0310 |
| **same formula, panel-roster truth** | 0.0645 | 0.0268 | **−0.0377 — OPPOSITE SIGN** |
| runtime | 1m52s | 10m36s | 5.7× |

⚠⚠ **THE PRECISION SIGN IS NOT STABLE.** Identical formula, identical two partitions, two different
truths → opposite signs. And the panel roster's own precision advantage does **not** survive
family-clustered resampling (40.7% positive / 19.9% tied / 39.4% negative). **A definition-level default
whose every downstream number is on record may not be moved on an axis whose sign flips with the choice
of truth.** The only statistic robust across every instrument is RECALL, which is neutral to positive,
never negative.

⚠ **Correcting the record**: the frozen-metric numbers that *looked* like they condemned genomic span do
not survive the frozen metric's own rules. `contamination`'s docstring reads *"⚠ its denominator is
prediction-conditioned; it is a DIAGNOSTIC, never the headline"* — the project's own banned metric shape
— and `purity_strict` keeps unlabelled nodes in the denominator **by design**, so both move whenever a
method merges more (unlabelled fraction 0.5699 → 0.6096). The frozen metric's headline recall statistic
was **identical**. It does not adjudicate this question in either direction.

### The two results that were said to be in tension — PARTLY resolved
1. **RESOLVED — they were never the same operation.** "recall 0.698 → 0.788 while precision 0.709 → 0.269"
   is a **locus-enlargement** result: it changes what a NODE IS (and 17.9% of RNA loci then span >1 DNA
   copy). The substrate question changes only *which sequence an edge is computed on*, with the node set
   **held fixed** (244 identical reps here; the 540-rep A/B likewise). The project's own rule — never
   judge a change to what a node is on node-level metrics — says these are different operations, and the
   data agree: at fixed node boundaries the added edges are **97.0% same-family (841 T / 26 F of 867)**.
2. **RESOLVED — the A/B was a RECALL result misreported as a PRECISION result.** Re-measured at the
   literal shipped tier, "precision 0.908 → 0.916 UP" is **+0.0053, CI [−0.0216, +0.0412], P(Δ>0) = 0.628**
   — a coin flip whose sign flips on the **coverage form**, not the tier. What survives is the recall
   half: **+176 true edges, CI [+67, +312], P(Δ>0) = 1.0000.** Quote it as *"a large recall gain at
   neutral precision"*, never as *"precision goes up"*.
3. ⚠ **NOT RESOLVED, and stated plainly**: whether genomic span improves or degrades **catalog-level**
   precision at the O1 site. Three independent measurements — the 540-rep A/B (P = 0.628), the 244-rep
   panel pairs (0.407 up / 0.394 down), the frozen metric under RefSeq-coarse truth (down) — give a coin
   flip, a coin flip, and a negative. **The sign is truth-definition-dependent and the available
   substrates cannot settle it.** That is a measurement limit, not a preference.

⚠ **Open, and NOT closed by this decision**: `docs/o1_investigations.md#the-joint-dna-rna-family-definition-retracted` §1 states *"`E_r` is
computed on `seq_g` only"* and §1.2 rejects exon-sum edges as *"strictly weaker"*. **The shipped O1
default is exon-sum.** The binary contradicts its own spec on this axis; O-4 keeps the binary's default
and records the contradiction rather than silently resolving it in either direction.
