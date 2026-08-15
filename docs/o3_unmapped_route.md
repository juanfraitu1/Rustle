# O3 — the unmapped-read route: is there ANY way to find a copy missing from the genome?

**Status 2026-08-15.** Companion to [`o3_missing_copy_evidence.md`](o3_missing_copy_evidence.md). That
document ends on a mechanism: on a whole genome the clipping statistic is structurally silent, so the
transcript-side detector could not fire and its 0/915 says nothing about the phenomenon. This document
asks the follow-up question directly — **is there any transcript-side route left, and can the reads that
fail to align at all be that route?** — runs it end to end, attacks it three ways, and states what it
bounds.

---

## READ THIS FIRST — the limits, before any number

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

## 1. The question, and why the previous detectors could not answer it

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

## 2. What the excision calibration established: a missing copy has two fates

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

## 3. What is actually in the unmapped pile

**Setup, re-verified in a single BAM pass:** 13,516,443 primary mapped + **959 flag-4**. Unmapped read
lengths: median 69 bp; **199 are ≥ 500 bp** (20.75% [18.23, 23.46] of 959), median 3,013 bp,
Σ = 634,416 bp. `GGO.fasta` = 3,545,850,636 bp vs NCBI's 3,545,834,224 bp over 25 scaffolds ⟹ **no
unplaced scaffolds**, so "absent from `GGO.fasta`" = absent from the whole assembly. Annotation is the
matching RefSeq `GCF_029281585.2-RS_2024_02` (isolate KB3781, male).

### 3.1 It does not dissolve under a permissive preset

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

### 3.2 It is FLNC-shaped, not junk

Against 800 mapped control reads ≥ 500 bp: mean-QV median **39.5 vs 39.6**; polyA absent in both
(refine-trimmed — the control proves absence is uninformative); **0/199 primer remnants**, 0/56 SMRTbell
adapter, 0 palindromic concatemers; 3-mer entropy median **5.79** of 6.0.

### 3.3 It is not contamination and not run-exclusive

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

### 3.4 Accounting of the 199

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

### 3.5 The locus

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

## 4. THE BOUND — as it was stated, with its conditions

*This section records the bound as shipped. §5 breaks it. Read both.*

### 4.1 The floor, and why it is not the binding constraint

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

### 4.2 The bound as shipped

> M ≤ **8.5** expressed reference-absent copies. Poisson-95 upper U(1) = **4.744** on **D = 1** observed
> detection, divided by π̄ = **0.5553** (detection power integrated over the 915-locus `n_clean`
> distribution). Point estimate 1.80, 95% CI [0.09, 8.54].

Per-expression form: ≤ ~9 copies above 431 `n_clean`, ~8 above 1,000, ~7 above 5,000. **D = 1** counts
the STON1 cluster as one genuine absence; if it were a contaminant the bound would *tighten* to M ≤ 5.4,
so D = 1 is the conservative choice. If STON1 + GTF2A1L are two independent absences, D = 2,
U(2) = 6.296, **M ≤ 11.3**.

### 4.3 What it explicitly did not cover

1. Copies **absorbed** onto a surviving paralog that leak fewer than 8 unmapped reads. π̄ < 1 corrects for
   these **only under the assumption that the excision panel's fate distribution transfers** to genuinely
   reference-absent copies. ← *this is the assumption that failed; §5.*
2. Copies expressed **below ~32 `n_clean`** (6/915 of catalog).
3. The **soft-clip stratum**, never examined — only flag-4 records were.

The 123,230-read SRA-minus-BAM gap was declared *not* an exclusion; §6 revises the arithmetic of that
claim but not its direction.

---

## 5. The attacks, and how each fared

Three independent attacks. **Two succeeded** (one fully, one partially), and they succeeded on the same
point by two different measurements.

### 5.1 "The 199 reads are junk / zero information" — **REFUTED**

Three ways, all in §3: alignment geometry (1,230 bp median contiguous M-block, MAPQ 60, identity 0.9805,
two introns); the corruption confound is *absent* (STON1 reads are the cleanest in the pile); and the
junk stratum, which genuinely exists (χ² = 540.5 over the 959-read pool), **does not contain STON1**
(p = 0.192, 4/4 runs). Simulation confirms the run-composition screen is properly sized (pass rate
0.941–0.953 at n = 8…200 under the mapped-composition null), so passing it is informative. The axis trap
was checked and cleared: every shipped number reproduces exactly on the `n_clean` axis.

### 5.2 ⚠⚠ "π̄ = 0.5553 does not transfer to reference-absent copies" — **SUCCEEDS. This breaks the headline.**

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

### 5.3 "The 123,230-read gap hides the answer" — **premise REFUTED, arithmetic PARTIALLY BROKEN**

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

## 6. ⚠ The gap probe's arithmetic, corrected — and what remains unexamined

### 6.1 The gap itself reconciles exactly

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

### 6.2 What the attack broke, and the required edits

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

### 6.3 Net effect and closure cost

**Point estimate unchanged.** Candidate-grade additions measured **0/15,254**; Poisson-95 upper **24.2
reads** over the visible gap, **127.8** over the full ~650,000 — against a floor of ≥ 8 in one cluster.
The gap is **no longer the binding uncertainty**: §5.2 is.

**Cost to close what remains** (the full-depth re-run of the two pre-`refine` runs): 5.44 GB gz download,
~18 GB uncompressed, ~1–2 h wall clock on 5 cores; recovered population is by construction enriched for
concatemers/non-poly(A) — the exact chimeric-cDNA false-positive mode already burned here.
**Recommendation stands: not worth recovering.** The ~527,000 invisible reads on SRR27438212/13 cannot be
recovered at all from the deposit; only re-processing from the raw PacBio subreads would expose them.

---

## 7. Honest verdict

**The unmapped route is VIABLE-BUT-CONDITIONAL for unique sequence and DEAD for collapsed paralogous
copies.** It is not one route with one bound; it is two routes with opposite verdicts, and the wrong one
carries the O3 question.

### 7.1 The replacement statement — quote this, not M ≤ 8.5

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

### 7.2 Why this is still a result worth writing

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

### 7.3 ⚠ How this relates to the DNA-side zero — do not merge them

The DNA-side probe (assembly-vs-assembly, the field's S1 standard, three assemblies of one gorilla)
returned **0 collapses at 816/817 loci** with the instrument validated at an FP floor of 0/817, and that
zero was **predicted**: Yoo/Rhie's 1–2 Mbp collapse per haplotype over the 1.1224% of the genome those
windows span predicts only 0.47–0.94 collapses ⟹ underpowered by construction.

⚠ **The transcript-side zero and the DNA-side zero are concordant but not independent confirmations of
each other, and neither is evidence of absence-of-absences.** The DNA side says *"an instrument with
measured sensitivity found nothing in a compartment too small to expect anything."* The RNA side says
*"in the stratum where collapses live, this instrument has no measured sensitivity at all."* Reporting
them as "two independent zeros" would be exactly the denominator error this project has retracted before.

### 7.4 The one experiment that would change the verdict

Re-run the excision calibration with the **survivor replaced by a consensus/mosaic of the two copies**
(what an assembler actually emits when it collapses), rather than by N-masking one copy and leaving the
other at full divergence *d*. That, and only that, measures π in the regime that matters. Predicted
direction: π_collapsible falls further, from 0.0286 toward 0, and stratum 2's "bounds nothing" hardens
from a CI statement into a mechanism.

Explicitly **not** worth doing: the flank-collapse variant named in the original `unsafe_if` (§5.2 shows
it cannot move the measurement), and the SRA gap recovery (§6.3).

---

## Artefacts

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
