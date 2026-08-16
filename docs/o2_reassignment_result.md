# O2 — REASSIGNMENT on real reads with structural-anchor truth: the result (2026-08-15)

**Status: RUN, ATTACKED, UNDERPOWERED — and the route itself is refuted for the population O2 cares
about.** Companion to [`o2_reassignment_ground_truth.md`](o2_reassignment_ground_truth.md) (the design)
and to `/home/juanfra/winloci_scratch/o2_truth_real/PREREG.md` (the pre-registration, sha256
`0154e72a077dfa002eb72c418e5d566b0766f88c5ea76551e8a99c748de79791`, frozen before any read was
selected).

The experiment was executed exactly as pre-registered, scored through the real `copy_assign` binary,
and then attacked from three independent angles. All three attacks landed. The verdict did not change;
three of the four quantitative statements in the first draft of the result did.

---

## READ THIS FIRST — the limits, before any number

1. **n = 53 reads, and the pre-registration abandons the head-to-head at n < 150.** Prereg §8, written
   before scoring, binds: *"If n < 150, the head-to-head is abandoned; only descriptive statistics are
   reported and the result is labelled UNDERPOWERED."* Expected n was ~680 (plausible range 250–1,500).
   **No significance claim is made in either direction.** The phrase "assigns better than minimap2" is
   not written anywhere in this document.
2. **n = 53 is really n_eff ≈ 3.** The 53 reads come from **2 families, 2 anchors, 13 distinct evidence
   vectors (11 informative)**. Family-clustered ICC on the argmax diagnostic is **0.694**, design effect
   **18.5**, **n_eff ≈ 2.9**. Every pooled rate below is quoted with its per-family decomposition
   because the pooled value routinely lies **outside** both family values.
3. **Scored through the pipeline, not a proxy** (T8 satisfied): the real `copy_assign` binary,
   `--release`, `--families`/`--copies-fa` from the shipped catalog (sha256 recorded before scoring),
   no detection / refine / rescue path. Nothing here carries the `PROXY` label.
4. **The truth is non-circular instance-wise, but the design does not guarantee it.** 29/53 = 54.7% of
   N rests on anchor **A00018, which has 354 genome-wide alignment records** (336 at ≥0.80 query
   coverage, ≥14 chromosomes, every non-self hit MAPQ 0) — a dispersed ~1 kb element. It is non-circular
   only because the nearest other instance is **6.6% divergent** (per-read exact-25-mer containment:
   self 0.710–1.000 median 0.974 vs nearest paralogous instance 0.210–0.255), and that was measured
   **after the fact**. No prereg filter (L1–L9) tests genome-wide anchor uniqueness. **L10 is now
   required** (§3.4).
5. **The truth set is the intron-retention minority.** Both surviving anchors are copy-specific
   **intronic** insertions. Of the primaries overlapping them at the carrier locus, **88–92% splice the
   anchor out** (`N`) and can never be labelled: A00018 387/421 = 91.9%, A00035 156/180 = 86.7%. Pooled,
   **N is 53/601 = 8.82% [0.0681, 0.1136]** of the reads present at the anchor. The result does not
   generalise to spliced transcripts, which are the population O2 actually operates on.
6. **⚠⚠ The selection effect is the opposite of the one pre-registered, and it is fatal to the route.**
   Prereg §10 hypothesised *"anchored families are plausibly more divergent"*. Measured: anchored
   families are **less** divergent (median `de` 0.0698 vs 0.1830, Mann-Whitney AUC 0.2414,
   p = 2.7e-07) and near-tie reads are **140× enriched** inside them (23.83% vs 0.17%). The
   anti-correlation with the contested population is created **one stage later, by filter L1**, and it
   is total: in the 6 contested L6-passing families, **29/29 candidate anchors were rejected by L1**.
   The two families that survived to be scored are **more divergent than every panel family with a
   near-tie rate ≥ 5%**. See §5.2 — this is the single most important finding in the document.
7. **The C2 sham control does not cleanly pass.** Pooled it does (4.9 pp, under the 10 pp threshold);
   stratified it **reverses**: GWFAM364 +16.1 pp, Fisher p = 0.00070, in the family supplying **100% of
   O2's commitments**. Reported as a caveat, not as a clean pass.

---

## 1. The gap this closes, and why the alternatives are not available

O2's abstention capability has a non-circular validation on real reads (whole-genome excision: held-out
TPR 0.5066 / FPR 0.0280, **AUC 0.7995**, against MAPQ's **0.4944**). **Reassignment has none.** Its only
real-data evidence is agreement with minimap2's primary flag — which *is* the "98.4% restatement"
complaint, so the validation and the criticism have always been the same measurement.

Every obvious alternative is closed:

* **PSV-derived truth is circular.** A PSV is a scoring decision, and O2 scores PSVs.
* **Simulation is planted.** `sim_ground_truth` (2-chromosome synthetic genome) and `k0_flank`
  (4-copy planted locus, 60 reads/copy) are both airtight and both synthetic.
* **Excision truth answers a different question.** There the true copy is *absent*; that is the
  abstention test, already done.

A **structural anchor** — sequence present in copy A and absent from copy B — is presence/absence rather
than scoring. Label the read by the anchor, **trim the anchor away**, hand the method only shared
sequence, score against the label. That is the construction. §5.1 shows the presence/absence premise is
true only *instance-wise*, and only because it happened to be checked.

---

## 2. Construction (as pre-registered, as executed)

**Substrate.** Exon-sum spliced representatives `q915_exon.fa` — not genomic intervals (rep-vs-rep
reproduces the shipped E_r rule 103/104; genomic-vs-genomic 1/3). Anchor discovery at the **shipped E_r
sensitive tier** (`minimap2 -c -X --no-long-join -k 11 -w 5 --cs`, identity ≥ 0.60, coverage ≥ 0.50 of
the shorter), never `asm20`. Read alignment at the source BAM's exact `@PG`
(`-ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`, only `-t` changed). No number crosses
between the two presets (T13).

**Reads.** `GCA_029281585.2_flnc_mm.bam`, the matched individual (KB3781 "Jim", same cell line as
mGorGor1). `-F 2308` before any per-read statistic. Anchor coverage ≥ 90% of **aligned reference blocks
excluding `N`** — a read whose intron merely spans the anchor does not carry it — and the overlap must
be **internal** (no clip within 20 bp of either boundary), read length ≥ 500 bp.

**Trim: TRUNCATE, never EXCISE.** Anchor ± 50 bp removed, **longer flank kept**, minimum residual
300 bp. Internal excision was rejected *a priori* in the prereg: a re-joined read must open a
~anchor-length deletion against the carrier and aligns contiguously to the sister — a systematic
manufactured push toward reassignment.

**Denominator (T1).** `N` = anchor-labelled ∧ selection-passing ∧ residual ≥ 300 bp, **fixed at freeze,
never conditioned on any scorer's output**. Abstentions count **against**. Out-of-scope reads: **0**.

### 2.1 Yield, filter by filter

| stage | count |
|---|---|
| within-family pairs with ≥ 1 PAF record | 144 / 162 |
| best-record identity ≥ 0.60 | **144 / 144** (identity never fails) |
| **L6** best-record coverage ≥ 0.50 of shorter ⟹ E_r edge | **18 / 144** (median best coverage 0.061) |
| raw anchors (indel, N-free run ≥ 100 bp) on E_r-edge pairs | 71 in 16 families |
| *diagnostic, L6 not applied* | *214 anchors in 52 / 162 families* |
| **L1** anchor absent from sister's **genomic** span ± 20 kb | **6 survive; 65 rejected** |
| certified anchors → reads | **2 anchors yield reads; 4 yield zero** |
| **N** | **53** (GWFAM215 29, GWFAM364 24) |

Two facts from this table matter more than the yield:

* **L6 is the binding filter, not the aligner preset.** The shipped tier finds anchors in 52 families
  vs `asm20`'s 34 (so 34 was indeed a lower bound), but requiring the pair to actually form an E_r edge
  cuts 52 → 16, because the 0.50 coverage-of-shorter floor almost never clears on concatenated-exon
  substrate.
* **L1 kills 65/71, and the rejected anchors are not marginal: median coverage in the sister genome is
  1.000.** ~92% of rep-vs-rep indels are **expression differences of a read-derived envelope, not
  presence/absence of sequence.** Without L1 this would have been an expression assay. Prereg design
  fact #2 was the right call.

The 4 zero-yield certified anchors have ~180 overlapping primaries each and **zero** at ≥ 90% aligned
coverage: those reads splice *over* the anchor. That is the `N`-exclusion rule working (and it is the
same fact as limit #5).

---

## 3. Leakage — what passed, what is not interpretable, and what was missing

| id | check | result | verdict |
|---|---|---|---|
| L2 | residual overlap with the anchor interval after realignment | 0 / 53 | pass |
| L3 | exact 25-mer shared with anchor or revcomp | 0 / 53 | pass |
| — | trimmed read aligns to its own anchor sequence | 0 / 53 (max cov 0.000) | pass |
| L4 | blinding: FASTQ names `T0000xx`, truth in a separate sha256'd file, absent from every scoring command line | — | pass |
| L5 | shipped family definitions, sha256'd **before** scoring, never rebuilt | — | pass |
| L7/L8/L9 | N-free anchors; exon-block arithmetic 0/324 fail; rep slice == genomic slice **71/71 = 100%** | — | pass |
| L1 positive control | 20 shared exon blocks pushed through L1 | **20 / 20 rejected** (min cov 1.000) | fires (T20) |
| (b) | junction structure by label class | 12 vs 11 junctions, 0 shared | **NOT INTERPRETABLE** |
| (c) | AUC of residual length alone | 0.3089 | **NOT INTERPRETABLE** |

**(b) and (c) must not be quoted as "no leak".** The two label classes are different loci in different
families, so any pooled class comparison measures the locus, not the copy. The intended within-family
version does not exist in this panel — the §6.4 one-sidedness problem, now empirical.

### 3.1 Four attacks on the label that FAILED (reported because they failed)

* **"`D` counts toward anchor coverage, so a sister read that opens a deletion gets labelled carrier."**
  A real code hole (`s5_reads.py`, `CONS_REF=set([0,2,7,8])`). Empirically closed: coverage by `M/=/X`
  alone is min **0.9890**, median 1.0000; `D` contributes max 0.0110, and no read has >11 bp of `D` over
  a 1001 bp anchor.
* **"The read never contains the anchor; only its alignment does."** Alignment-free exact-25-mer
  containment in the raw `SEQ`, both strands: min **0.6883**, median 0.9765, 0/53 below 0.5. The reads
  physically carry the sequence.
* **"Truth = the primary flag, because minimap2 put the read there."** Enumerating every primary read in
  both copies' spans ± 20 kb of both families (1,263 reads) and scoring anchor containment independent
  of placement: 57 anchor-carrying reads, **0/57 whose primary is not at the carrier**. *Scoped:* see
  §3.3 — this frame cannot see reads whose primary is outside both loci entirely.
* **"Trimming leaves residual signal."** L2 0/53, L3 0/53, 0/53 align to their own anchor; C2's pooled
  4.9 pp is consistent.

### 3.2 The reads are not pseudo-replicates at the molecule level — but the *evidence* is

27/29 and 23/24 distinct (start, end) pairs; 28/29 and 24/24 distinct lengths. Both copies are genuinely
expressed (primary `-F 2308`: GWFAM215 :0 = 197, :1 = 437; GWFAM364 :0 = 201, :1 = 147). But
`anchor.psv_reads.tsv` collapses the 53 reads to **13 distinct evidence vectors (11 informative)**, and
**27/29 GWFAM215 residuals start at the identical coordinate 84205633**. The replication is in the
molecules; it is absent from the evidence.

### 3.3 Leakage that SURVIVED: the sampling frame is conditioned on the primary flag

`s5_reads.py` selects reads via `bam.fetch(anchor)` + `-F 2308`. **A carrier read that minimap2
primary-placed outside both copies can never enter N.** Measured at A00035: **37 reads have an internal
≥ 90%-coverage alignment across the anchor as a SECONDARY record with their primary outside both
copies** (median 90 match / 11 mismatch over 101 bp) versus 24 admitted — up to **61% of anchor-spanning
candidates removed, every one of which scorer (b) would have scored INCORRECT**. (A00018: 0.) The §3.1
enumeration searched only *within* the two loci ± 20 kb, so it could not see these. Consequence: the
primary-flag number is not merely a C1 artefact — its denominator is partly the scorer's own correct
decisions.

### 3.4 Leakage that was NEVER CHECKED: L10, genome-wide anchor uniqueness

The prereg's non-circularity premise — *"a structural anchor is presence/absence, not scoring"* — is
**false as written** for A00018. Presence/absence of *this instance* is what the label needs, and
separating instance from element is exactly an identity/alignment decision. **L1 is pairwise-local**
(sister interval ± 20 kb only). Under the prereg as written, a read from any of the other 353 instances
placed at `GWFAM215:1` by minimap2 would have been labelled `GWFAM215:1` with no check firing.

**Required addition — L10:** every certified anchor must be aligned genome-wide, and the per-read
exact-k-mer containment margin against its nearest paralogous instance reported. Post hoc here:
A00018's best non-self identity is **0.9336**, self containment 0.710–1.000 vs nearest instance
0.210–0.255 (margin ≈ 2.8×); the other five certified anchors are genome-unique at ≥ 0.5 query coverage.
The label survives. The **design** does not.

---

## 4. The result

**Unit: the READ. Denominator: N = 53, fixed at freeze. Abstentions count against.**

| scorer | accuracy | 95% Wilson | note |
|---|---|---|---|
| **C1 constant baseline** (always predict the anchor-carrying copy) | **53/53 = 1.0000** | [0.9324, 1.0000] | **1.0000 by construction** — labels are 100% one-sided |
| **(b) minimap2 primary flag** | 53/53 = 1.0000 | [0.9324, 1.0000] | all 53 primaries MAPQ 60 — **= C1, and see below** |
| **(c) O2 divergence/PSV channel** | **4/53 = 0.0755** | [0.0297, 0.1786] | scored through `copy_assign` |
| S1 accuracy-when-committed | 4/4 = 1.0000 | [0.5101, 1.0000] | **descriptive only** (T1: denominator conditioned on the method's own output) |
| S2 abstention | **49/53 = 0.9245** | [0.8214, 0.9703] | 40 `ambiguous` + 9 `tied` + 4 `assigned` |

**The primary-flag 53/53 measures neither aligner skill nor accuracy.** Three reasons, in increasing
severity: (i) it exactly equals C1, which is 1.0000 by label imbalance; (ii) scorer (b) asks whether the
primary lands inside a copy's genomic interval, and those intervals are **read-derived envelopes built
from these reads** — 2/53 reads have aligned start exactly equal to their copy's `orig_start`, 4/53
within 5 bp of the start, 29/53 within 5 bp of the end, 46/53 fully inside; (iii) §3.3, the frame
excludes the primary flag's own errors. **Restate it as: the truncated read reproduces the untrimmed
primary placement that both selected it into N and defined the interval it is scored against — placement
stability under truncation, not accuracy.**

**S1 = 4/4 is 2 distinct evidence vectors at 1 locus.** Three of the four commits share the identical
allele string `GAAAATATTCCAGAAACTTAATATTTGTGATG` with **10 reads that did not commit**; the fourth
differs by one base. The gate is not even a function of the evidence vector. n_eff = 1.

### 4.1 The argmax diagnostic — pooled value RETRACTED

The first draft reported *"O2's raw likelihood argmax is correct on only 32/53 = 0.6038 [0.4694, 0.7241],
favouring the wrong copy on 21/29 GWFAM215 reads."* **That pooled number is withdrawn.**

| family | n | O2 correct | commits | **argmax correct** | abstain |
|---|---|---|---|---|---|
| GWFAM215 (A00018) | 29 | 0 | 0 | **8 (0.276)** | 29/29 |
| GWFAM364 (A00035) | 24 | 4 | 4 | **24 (1.000)** | 20/24 |

The pooled 0.6038 lies **outside both family values**. Cluster support is {0.2759, 0.6038, 1.0000};
ICC 0.694, deff 18.5, **n_eff 2.9**; the Wilson CI is anticonservative by ≈ 4.3×. Worse, **all 9 `tied`
reads have `n_decisive = 0`, `margin = 0.000`, and `assigned_copy = 0`** — every one. Truth is copy 0 in
GWFAM364 (7 free "correct") and copy 1 in GWFAM215 (2 free "incorrect"), so the tie-break is the **copy
index**. Removing them: **25/44 = 0.5682 [0.4222, 0.7032]** — a coin flip.

**The "21/29 wrong" is one haplotype class in a 23–27 bp window.** At the evidence-vector level the
GWFAM215 split is **3/6 = exactly 0.500**; the 21-vs-8 read split is how many FLNC molecules of each
haplotype were sequenced. 26/29 reads observe **all** their PSVs inside a **27 bp** window where the two
copies differ at **16 of 27 bases** (59% divergence, at an exon-block edge — column 25 → column 26 jumps
4.6 kb), while the **whole ~2 kb residual aligns 24% better to copy 1 for all 29 reads** (AS ≈ 1900 vs
1500). A read cannot be copy-0 at 16 consecutive discriminating columns and copy-1 over 2 kb; one of the
two is broken, and **the GWFAM215 label is not certified** — L1 checked the sister's ± 20 kb only, on a
**single-haplotype assembly of a diploid animal**, for an anchor that is a 351–354-instance dispersed
repeat with only **one** hit ≥ 95% identity genome-wide. A polymorphic insertion of that element at the
copy-0 locus predicts exactly the 19 observed copy-0-allele reads.

### 4.2 Resolvable fraction — and why it does not answer the K = 0 question

* Prereg §7 stratum **S_PSV+ = 53/53 = 1.0000** [0.9324, 1.0000]; **S_PSV0 is EMPTY**, so §6.2's S3
  stratification is degenerate (one stratum). The K = 0 stratum does not appear in this panel.
* The pipeline's own identifiability gate, `n_decisive ≥ 1`: **44/53 = 0.8302** [0.7077, 0.9080].
* **⚠⚠ None of these 53 reads is contested.** On the realigned trimmed BAM the minimum relative AS gap
  `(best − second)/best` is **0.103**, median 0.208, and **0/53 fall within 5% of the primary's AS**.
  Untrimmed the gaps are larger still (min 0.4116, median 0.5230, 0/53 near-ties).
* All 53 reads retain PSVs (median 70, min 38) yet O2 abstains on **92.45%** of them. That abstention is
  the **IsoCon significance gate**, not an absence of distinguishing sequence.

### 4.3 The C2 sham arm — pooled pass, stratified reversal

|  | anchor abstain | sham abstain | Δ | Fisher p |
|---|---|---|---|---|
| GWFAM215 | 29/29 = 1.0000 | 411/426 = 0.9648 | **−3.5 pp** | 0.613 |
| GWFAM364 | 20/24 = 0.8333 | 179/180 = **0.9944** | **+16.1 pp** | **0.00070** |
| pooled | 49/53 = 0.9245 | 590/606 = 0.9736 | +4.9 pp | 0.069 |

The pooled 4.9 pp is below the pre-declared 10 pp threshold, and O2/primary agreement-among-committed is
100% in both arms (0 pp). But the per-family differences point in **opposite directions** with different
arm composition (anchor 45% GWFAM364, sham 30%) — a **Simpson reversal**, the same failure mode as this
project's k-core p 0.0086 → 0.8125. In GWFAM364 the anchor arm commits **4/24 = 16.7%** vs sham
**1/180 = 0.6%**, a 30× difference, in the family supplying **100% of O2's commitments**. The arms are
also not matched: sham `n_decisive` median **251** vs anchor **33.5** (T20 — this control controlled for
nothing), and sham anchors are shorter than real (median 101 bp vs 265.5 bp), because the longest shared
aligned segments available are 1026 bp and 330 bp. (Bookkeeping: the sham table has 606 rows; the first
draft quoted 605. The rate is unchanged to 4 dp.)

**Do not write "NO TRIM ARTEFACT declared."** Write: *pooled pass, stratified +16.1 pp reversal in the
family that owns every commitment, on unmatched arms.* Prereg-legal, but it must be reported.

### 4.4 The head-to-head — ABANDONED BY PRE-REGISTRATION

Prereg §8 binds at n = 53 < 150. The following are **descriptive, non-inferential, prereg-void**, and
are recorded only so nobody recomputes them and claims more:

Δ = 4/53 − 53/53 = **−0.9245**; McNemar discordant b(primary-only) = 49, c(O2-only) = 0, n_disc = 49,
exact two-sided p = 3.55e-15; family-clustered bootstrap 95% CI on Δ = [−1.0000, −0.8333] — with 2
families the bootstrap has 3 possible resamples and **quotes nothing**.

**No significance claim is made in either direction.** Per prereg §6.5 this is an **absence of a
demonstrated improvement**, not evidence of parity, and not a success. Equally, the apparent primary-flag
win is not evidence of aligner skill (§4).

**What n would be needed.** The n ≥ 150 floor binds before power on discordance does (observed
discordance 92.5% would give n_disc ≈ 100 by n ≈ 108). At the observed 26.5 anchor-labelled reads per
anchored family, 150 reads needs ~6 L1-surviving families; at the measured L1 family survival of 2/16
(12.5%) and L6 pass rate of 18/162 (11.1%), that is **~430 two-copy families screened — ~2.7× the
current panel.** §5.2 shows why that arithmetic is nevertheless the wrong plan.

---

## 5. The attacks, and how each fared

Three independent attacks. **All three partially or wholly succeeded.** Nothing below is a defence.

### 5.1 Attack A — "the label is alignment-derived": PARTIALLY SUCCEEDS

Four routes to break the label failed empirically (§3.1). Three findings survived and are folded into
the document above: the **354-instance dispersed repeat** carrying 54.7% of N with no design-level check
(§3.4, new filter L10); scorer (b)'s **decision boundary built from the reads it scores** (§4); and
**N being the ~8.8% intron-retention fraction** of reads at the anchor (limit #5). Net effect: the label
holds, the *guarantee* does not, and the generalisation shrinks to unspliced/IR reads.

### 5.2 Attack B — "anchoring selects against the contested population": SUCCEEDS, and INVERTS the stated mechanism

The first draft explained the "0/53 near-ties" finding as *"structural anchors exist only where two
copies are structurally divergent, so the anchor-labelled population is anti-correlated with the near-tie
population."* **That mechanism is refuted.** Measured on the 144 panel families with a within-family
E_r-tier record:

| group | n fam | median `de` | IQR | median best_cov |
|---|---|---|---|---|
| carries a ≥ 100 bp N-free anchor | 52 | **0.0698** | 0.0231–0.1593 | 0.284 |
| no anchor | 92 | **0.1830** | 0.1264–0.2434 | 0.031 |

Mann-Whitney AUC **0.2414**, z = −5.15, **p = 2.7e-07**: anchor-carrying is a proxy for *alignability*
(best_cov AUC 0.8250, p = 1.0e-10), and alignability anti-correlates with divergence. And near-tie reads
live almost entirely inside anchored families — per-read relative AS gap from the source BAM
(`-F 2052`, 459,473 records over the 324 panel intervals, 319,606 reads with a panel primary):

* anchored families **15,229/63,901 = 23.83%** near-tie (≤ 5%) — essentially the 21.75% contested-population figure
* unanchored families **376/216,673 = 0.17%**; no-alignment families 0/39,032; panel pooled 4.88%

**Anchor discovery is enriched ~140× TOWARD the contested population.** The kill happens one stage later,
at **L1**, and L1 kept the wrong end of the distribution:

| family | reads | near-tie | `de` | fate |
|---|---|---|---|---|
| GWFAM409 | 1335 | 0.9985 | 0.0231 | L1 rejected |
| GWFAM262 | 182 | 0.9945 | 0.0022 | L1 rejected |
| GWFAM253 | 922 | 0.9837 | 0.0126 | L1 rejected |
| GWFAM404 | 109 | 0.9817 | 0.0034 | L1 rejected |
| GWFAM29 | 1049 | 0.2898 | 0.0130 | L1 rejected |
| GWFAM372 | 732 | 0.1544 | 0.0048 | L1 rejected |
| **GWFAM215** | 557 | **0.0018** | **0.0383** | **SCORED** |
| **GWFAM364** | 344 | **0.0000** | **0.1347** | **SCORED** |
| (8 others) | | ≤ 0.0145 | 0.0337–0.1437 | L1 rejected |

Mechanism, measured: **in the 6 contested L6 families, 29/29 candidate anchors were rejected by L1, with
sister-genome coverage median 1.000, min 0.881.** In the near-tie regime a rep-vs-rep indel is *always*
an expression difference, never presence/absence. **Zero structural anchors exist there.** Spearman
across 144 families: near-tie rate vs `de` **rho = −0.689, p = 1.7e-16**; the contested population is the
low-divergence population (read-weighted `de` of the 15,605 near-tie reads: median **0.0147**, p95 0.0231).

Selection is also **read-level, within a single locus** — immune to the "different loci" objection that
voided leakage checks (b) and (c): GWFAM215's 29 selected reads have median untrimmed AS gap **0.4380**
vs **0.2113** for the other 528 reads at the same locus (AUC 0.9988, p = 1.4e-19).

**Honesty on this attack's own limits.** At n = 2 survivors the family-level claim is not statistically
established (Fisher, 2 survivors from 16 families of which 6 contested: **p = 0.375**). The defensible
form is the **anchor-level 29/29**, not the family-level count. Near-tie rates use secondaries falling
inside the 324 panel intervals only, so they are lower bounds (applied identically to all groups); the
within-family paralog is always in-panel, so §5.2's within-locus comparison is unaffected. `de` is the
frozen E_r-tier per-family value and is never compared to a `splice:hq` number (T13).

**Consequence for §4.4's scaling paragraph: it is refuted, not merely expensive.** Rule-of-three upper
bound on contested-family L1 survival is **≤ 10.3%** (Wilson ≤ 11.7%), so 2.7× panel expansion buys an
expected **0** contested anchors. Loosening L6 does not help either — L6 rejects on *coverage*, and
low-coverage pairs are the *high*-divergence ones. **The structural-anchor route is barren for near-tie
reads by construction; no amount of panel expansion turns this design into the reassignment test O2
needs.**

### 5.3 Attack C — "the diagnostics are pooling artefacts": SUCCEEDS

Delivered §4.1 (n_eff ≈ 2.9, the 9 index tie-breaks, the 27 bp window vs the 2 kb residual, the
uncertified GWFAM215 label), §4.3 (the Simpson reversal in C2), §4 (S1 = 2 evidence vectors), and §3.3
(the frame conditioned on the primary flag). Its own attacks that failed are reported in §3.1–3.2 and
above: the sister copy **is** expressed, the reads **are** distinct molecules, and leakage L2/L3 is
clean. It did not overturn `underpowered`; it deepened it.

### 5.4 What survives all three attacks, unchanged

`verdict: underpowered`; the abandoned head-to-head; C1 = 1.0000; `divergence_channel_accuracy = 4/53`;
the abstention rate 49/53; **0/53 near-ties**; and the absence-of-improvement language of prereg §6.5.
The attacks removed positive claims. They created none, and the result makes none.

### 5.5 Transferability figure, to be attached to every number above

The 53 reads describe a stratum holding **1 of the 2,955** near-tie reads in its own 16-family L6 panel
(0.034%), and **0.10%** of the 15,605 near-tie reads in the 324-interval panel, at family divergence
**1.6×–9.2× the contested median** (0.0383 and 0.1347 vs 0.0147). This is stronger than prereg §10's
"does not generalise to the whole catalog": **it specifically excludes the population O2 was re-scoped
onto.**

---

## 6. Verdict

**On real reads with non-circular structural-anchor truth, O2's reassignment value could not be measured:
the pre-registered design yielded n = 53 reads from 2 families / 2 anchors / 13 evidence vectors
(n_eff ≈ 3), which its own §8 stop rule abandons at n < 150 — and the binding filter (L1, which is what
makes the truth structural) removes 29/29 candidate anchors in exactly the contested, low-divergence
families where O2 operates, so the route cannot be scaled into the test it was built to be.**

Examiner-checkable: `truth.final.tsv` sha256
`87d1119f06ef9ccac3ff7f750ceee4d558eb9e336eaa2c2d236723e35f556e80` (53 rows, frozen before scoring);
prereg §8 line *"If n < 150, the head-to-head is abandoned"*; `L1_coverages.json` (65/71 rejected, median
sister coverage 1.000); 29/29 anchor rejections in the 6 contested L6 families.

---

## 7. Does O2's objective statement need further amendment?

O2 was restated on 2026-08-15 as **abstention under alignment-score near-ties**, not reassignment. **That
restatement stands and is strengthened, but it needs two amendments and one prohibition.**

1. **The restatement is not contradicted — and this experiment is not evidence for it either.** The
   abstention scoping rests on the excision validation (AUC 0.7995 vs MAPQ 0.4944). This experiment
   contains **0/53 contested reads**, so it neither supports nor undermines that. It must not be cited
   as corroboration.
2. **AMENDMENT 1 — say plainly that the reassignment component has no real-read validation, and that the
   only designed route to one is closed.** Current wording implies reassignment was merely
   *de-emphasised*. It should say: reassignment is validated **only in simulation** (`sim_ground_truth`,
   `k0_flank`); the structural-anchor route was designed, pre-registered, executed, and is **barren for
   near-tie reads by construction** (§5.2); therefore no non-circular real-read reassignment number
   exists or is currently obtainable. A future route must generate truth that does **not** require
   presence/absence structure, because presence/absence structure and alignment-score near-ties are
   mutually exclusive in this data.
3. **AMENDMENT 2 — the near-tie population now has a second, independent measurement, and it should be
   quoted.** Independently of the 21.75%-at-multi-copy-loci figure, near-tie rate inside E_r-anchored
   two-copy families is **23.83% (15,229/63,901)** vs **0.17%** in unanchored ones. The contested
   population is the **low-divergence, high-alignability** stratum (near-tie vs `de` rho = −0.689,
   p = 1.7e-16; read-weighted median `de` **0.0147**). O2's target should be stated in those terms — it
   is a *divergence-defined* population, not a MAPQ-defined one, and this is the sharpest handle on it
   the project has.
4. **PROHIBITION (unchanged, now with a second independent basis).** "O2 assigns better than minimap2"
   remains unwritable. Prior basis: the re-scope analysis (a secondary fits better by `de` in 1.96% of
   reads-with-secondary; net headroom ~0.1%). New basis: the only non-circular real-read reassignment
   experiment ever run is underpowered by its own pre-registration and its route is refuted. **Defend O2
   on abstention.**

---

## Artifacts

* Prereg: `/home/juanfra/winloci_scratch/o2_truth_real/PREREG.md` (sha256 `0154e72a…`)
* Frozen truth: `/home/juanfra/winloci_scratch/o2_truth_real/truth.final.tsv` (sha256 `87d1119f…`);
  sham `sham.truth.final.tsv` (sha256 `ff189403…`)
* Anchors: `anchors.tsv` (71 raw) · `anchors.certified.tsv` (6) · `L1_coverages.json` ·
  `L1_poscontrol.json` · `build_counts.json` · `L5_shipped_families.sha256`
* Scoring: `/mnt/linuxdisk/home/juanfraitu/o2_truth/score/` — `scored_anchor.json`,
  `anchor.assignments.tsv`, `sham.assignments.tsv`, `anchor.psv_reads.tsv`, `anchor.psv_copies.tsv`,
  `anchor.psv_cols.tsv`
* Big: `/mnt/linuxdisk/home/juanfraitu/o2_truth/` — `panel324.fa`, `panel324.er.cs.paf`, `trimmed.fq`,
  `trimmed.bam(.bai)`, `sham.trimmed.*`, `L1/`, `mm2_*.log`
* Scripts: `s1_panel_fa.py` … `s8_freeze.py` (same dir as the prereg)
* Machine note: the realignment peaked at **21.6 GB RSS of 25 GB** — never run two concurrently.
