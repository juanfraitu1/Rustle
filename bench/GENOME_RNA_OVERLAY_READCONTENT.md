# A2 — Read-Content Cross-Modal Corroboration of O1 (SEDEF identity ↔ PSV resolvability)

Status: **PARTIAL corroboration with a major provenance caveat.** Honest report; every refuter
finding is in the Limitations section, not buried.

Artifacts:
- `bench/genome_rna_overlay_readcontent.py` (analysis)
- `bench/genome_rna_overlay_readcontent.json` (numbers)
- `bench/genome_rna_overlay_readcontent_per_family.tsv` (per-family table, n=154)
- Inputs: `bench/psv_graph_genomewide.json`, `/home/juanfra/winloci_scratch/validated_families.tsv`,
  `/mnt/c/Users/jfris/Desktop/final.bed` (SEDEF self-alignment)
- Prior (substrate-shared) overlay being differentiated against: `bench/genome_rna_overlay.json`,
  `bench/GENOME_RNA_OVERLAY.md`

---

## 1. What A2 set out to do, and why the prior overlay failed

The prior overlay (`GENOME_RNA_OVERLAY.md`) grouped de-novo "transcript" sequences and overlaid them
on SEDEF segdups. Adversarial review caught it as **circular at the substrate level**: the
"transcript" sequence it grouped on is a REFERENCE-genome substring (`denovo_assemble_gate.py:58`
fetches `GGO.fasta`), the exact bytes SEDEF self-aligns. Different code, same data — its 205×
"lift" measures the reference agreeing with itself.

A2's goal was a corroboration that uses a **read-content** signal: correlate, per multi-copy family,
two estimates of copy divergence computed from *different data streams* —
- (i) **SEDEF reference %identity** between copy intervals (`final.bed` col-21; reference FASTA
  self-alignment, never sees a read), against
- (ii) **read-derived PSV resolvability** (PSV density / K / class `no_psv|partial|fully_resolvable`)
  from `psv_graph_genomewide.py::column_alleles`, which pileups actual HiFi/IsoSeq read bases.

**Prediction:** near-identical copies → HIGH SEDEF identity → FEW PSVs → `no_psv`/K-frontier
(read-unresolvable); divergent copies → LOW identity → MANY PSVs → `fully_resolvable`.

## 2. The cross-modality argument — and its hard ceiling

Corroboration is two different instruments estimating one latent quantity (here, true paralog
divergence, a common *cause*) and agreeing. Circularity is one instrument's output being re-read by
the other (a common *substrate*). **Perfect orthogonality is impossible and is not the bar:** both
quantities reflect the same copy divergence, so if both are valid they *must* correlate. What can be
independent is the DATA + CODE provenance — and that is exactly where A2 turns out to be only
partially independent (Section 7).

## 3. PRIMARY result — SEDEF identity ↔ read-PSV density (continuous, load-bearing)

n = 154 families (all SEDEF-linked; dedup reproduces json `n_copies` 154/154; my independent
re-derivation of `sed_ident` from `final.bed` matches the TSV to 5 dp, so the join is genuine).

| axis (vs SEDEF %identity)            | Spearman ρ | p        |
|--------------------------------------|-----------:|----------|
| psvs (raw, length-inflated)          |     −0.329 | 3.0e-05  |
| **psv_per_kb (density — headline)**  | **−0.443** | **9.1e-09** |
| K_per_copy (saturates at n_copies)   |     −0.103 | 0.21 (NS)|

Direction is exactly the prediction: more reference identity → fewer PSVs per kb. Independently
reproduced byte-for-byte (ρ = −0.4426, p = 9.08e-09).

**Confound controls — all survive.** Raw-PSV confound structure: mean_copy_len ρ=+0.384 (p=9e-7),
reads/coverage +0.278 (p=5e-4), n_copies +0.239 (p=3e-3).
- (a) **Normalize** to density removes length: ρ strengthens −0.329 → **−0.443**.
- (b) **Partial** Spearman(identity, psv_per_kb | log10 reads, n_copies, mean_len) = **−0.412**, p=1.4e-7.
- (c) **Stratify** n_copies==2 only (n=111) = **−0.456**, p=5.0e-7 (removes K~n_copies saturation).
- (d) **Coverage tertiles** (SEDEF identity is depth-free, so this guards the depth→`no_psv` artifact):
  low −0.360 (p=9e-3), mid −0.507 (p=1.5e-4), high −0.423 (p=2e-3) — significant in all three.
- (e) **Join-permutation null** (shuffle family pairing, B=5000): mean|ρ|=0.065, max|ρ|=0.272,
  p(|perm|≥0.443) = **0.0002** — the signal is the family join, not rank autocorrelation.
- (f) **Robustness (refuter-added, reproduced):** dropping the 27 zero-PSV `no_psv` families
  *strengthens* ρ to **−0.585** (p=5.4e-13) — the signal is a real continuous interior gradient,
  not a `no_psv`-vs-`full` pile-up. Leave-one-out ρ ∈ [−0.466, −0.433] (no influential family).
- (g) **Coverage orthogonality (reproduced):** Spearman(SEDEF, reads)=0.014 (p=0.87),
  Spearman(psv_per_kb, reads)=0.083 (p=0.31) — "more reads → more PSVs" cannot manufacture this.

**Honest magnitude note:** Pearson r = **−0.134** (vs Spearman −0.443). The relationship is
rank-monotone but weakly *linear* — SEDEF is ceiling-saturated and heavy-tailed, so Spearman is the
correct statistic, but −0.44 is a rank-correlation, NOT a linear effect size.

## 4. K-frontier vs SEDEF (directional support only — NOT significant)

Thesis claim under test: the read-derived K-frontier (`no_psv`/K=1) is *genuinely near-identical
copies*, not a coverage/pileup artifact. Mann-Whitney one-sided (`no_psv` > `fully_resolvable`):

- `no_psv` (n=27) median SEDEF id **0.9920** vs `fully_resolvable` (n=124) median **0.9889**;
  U=1927, **p=0.110**, rank-biserial +0.151, AUC 0.576. K==1 vs K≥2: 0.9920 vs 0.9889, p=0.123.

Direction is correct (the K-frontier sits at the higher-identity edge) but **not significant**:
SEDEF identity saturates near a 0.99 ceiling (min 0.837, q25 0.983, median 0.989, q75 0.995), which
compresses the binary contrast. This is **supporting-directional only**; the continuous density
Spearman is load-bearing. Do not headline the Mann-Whitney.

## 5. SECONDARY — membership overlay (flagged reference-shared, NOT a headline)

Mirror of the prior overlay's machinery: within-family distinct-copy pairs vs SEDEF-link PRESENCE.
Precision **0.830** (409/493 pairs linked), distance/chrom-matched null 0.00706, **lift 117.5×**,
permutation p=0.001. Link-fraction by read-class: `no_psv` 0.967, `partial` 1.000,
`fully_resolvable` 1.000.

**Flag:** this is substrate-shared at the MEMBERSHIP level — family membership is reference-derived
(`family_def_vg_reinforce.py`) and SEDEF is reference, so both axes are reference homology. High
precision is **expected by construction**, and read CONTENT does not enter the grouping at all, so
this number carries **strictly less** read-independence than even the prior RNA-POA overlay. Reported
only for (a) sanity that the read-validated families are the reference-homologous ones, and (b) as
the differential baseline. Not the corroboration.

## 6. DIFFERENTIAL — how much does A2 escape the prior's L1 circularity?

Prior overlay: precision 0.270 (excl-DNFAM0), **lift 204.9×** — a *binary link-presence*
coincidence with NO graded axis. It assigns the same "linked" label to a 0.91-identity pair and a
0.999-identity pair; its sequence side was reference substrings.

A2's increment is operationalized by holding membership constant — **restricting to SEDEF-linked
families (all 154 link), so the prior's entire presence signal is a CONSTANT** and cannot drive any
within-set result. Within that set, Spearman(SEDEF identity, psv_per_kb) = **−0.443** (p=9.1e-09).

So A2 *does* deliver something the prior structurally could not: a graded correlation on the
%identity gradient the prior was blind to, computed by a **different aligner** (minimap2 asm20 vs
SEDEF/lastz), over a **different interval** (read-covered exonic columns vs whole segdup block).
**That increment is real and in CODE + INTERVAL.** What it is NOT, is an increment in DATA — see
Section 7.

## 7. LIMITATIONS — every refuter finding, stated up front

Two of three adversarial lenses returned **refuted = true (major)**; the third (statistics) returned
**refuted = false (minor)**. All findings are reproduced and accepted.

**L-1 (MAJOR — hostile-examiner + meta-provenance, the load-bearing flaw): the headline axis is
reference-SEEDED, not read-discovered. Circularity is relocated, not escaped.**
`psv_graph_genomewide.py:135-141` declares a column a PSV iff `len(set(hap.values()))>=2`, where
`hap` = `refseq[p]` (reference backbone) ∪ `copy_cols[p]`, and `copy_cols` comes from
`minimap2 -ax asm20 ref.fa copies.fa` (line 97) — the REFERENCE copy FASTA aligned onto the
REFERENCE backbone. **No read base contributes to whether a column is a PSV.** Reads enter the
headline axis ONLY as (a) a coverage-cardinality gate `len(read_cols)>=3` (line 140 — read base
VALUES are never read; a column where all reads match the backbone but the asm20 copy differs still
counts) and (b) an exonic/expression restriction. Reads cannot ADD a PSV the reference doesn't
already show. So `psv_per_kb ≈ (1 − SEDEF_identity)` measured by a *second reference self-aligner*
over a read-masked exonic subinterval. The −0.443 is **two aligners scoring divergence on the same
reference bytes**, i.e. the prior's exact trap relocated from SEDEF↔grouped-bytes to
SEDEF-self-align ↔ asm20-self-align. The build's own `honest_caveats(1)` concedes this verbatim
("NOT different DATA for the divergence content … a milder version of the same trap that sank the
prior"), which directly contradicts the build's own "genuinely cross-modal in DATA+CODE" wording.

**L-2 (MAJOR — the genuine read-base axis fails). ~0% of the significant signal is read-derived.**
The ONLY axis where a read BASE determines anything is the assignment rate `one_copy/total`
(`psv_graph:153-173`, uses the actual read base) — and it is **NULL**: ρ=−0.138 (p=0.087) all-linked,
ρ=−0.109 (p=0.22) non-no_psv (the script self-labels it "weak_nonsignificant_coverage_confounded").
The clean read-discovered corroboration A2 was tasked to produce **does not independently hold.** The
significant `psv_per_kb` number must not be allowed to stand in for it.

**L-3 (membership is reference-selected).** WHICH pairs are tested comes from
`validated_families.tsv`, whose membership is reference-derived. The primary mitigates by asking the
graded *within-set* question (the binary membership substrate cannot generate a gradient), but
membership selection can still bound the gradient's range.

**L-4 (SEDEF dynamic range is compressed).** SEDEF emits only segdups ≳90% identity; every tested
family sits in a narrow top end (median 0.989). This is why the K-frontier Mann-Whitney is
underpowered (p=0.11) and why the discrimination is carried almost entirely by the PSV side
(range 0–246), not by the SEDEF side. Honest framing: "read-PSV content tracks the *narrow top end*
of the SEDEF scale," not "two wide-range estimators agree."

**L-5 (shared blindspots at the extremes).** If the T2T assembly collapsed two near-identical copies,
SEDEF cannot represent them as two AND reads have one locus to map to — both pipelines are blinded
the same way, inflating agreement at the high-identity end. Symmetrically, too-divergent copies fall
out of SEDEF (correlated missingness; `no_psv` link-fraction 0.967 vs full 1.000). Both effects live
at the truncated ends; the interior gradient (−0.443, and −0.585 among PSVs>0) is over families both
methods see.

**L-6 (interval mismatch — a strength and a noise source).** SEDEF identity is over the whole segdup
block (introns + flanks + exons); read-PSVs sample only expressed exonic columns. Different bases
(*reduces* circularity) but adds noise (exonic divergence under purifying selection can decouple from
block divergence). This is part of why |ρ|=0.44 is *lower* than two aligners on identical bytes would
give — the attenuation is noise, not read-supplied divergence content.

**L-7 (statistics lens — refuted=false, minor).** Every fair confound control survives independent
recomputation (Section 3). Remaining honesty edit: report Pearson −0.13 alongside Spearman −0.44 so
the magnitude is not read as a linear effect; report the PSVs>0 subset (−0.585) as the strongest
robustness check. No statistical fix required — the residual issue is provenance (L-1/L-2), not
statistics.

**Decisive break-test — NOW RUN, read gate proven INERT (2026-06-30, `genome_rna_overlay_breaktest.json`,
92 min over n=154).** Recomputed `psv_per_kb` **without the read-coverage gate** (count ALL asm20-divergent
columns, a pure reference quantity) and re-correlated with SEDEF. The gate keeps only **36%** of asm20
columns (9,089/30,893), yet the correlation is essentially unchanged: **gated ρ=−0.443 → ungated ρ=−0.420**
(p=6e-08), gated-vs-ungated per-family PSV counts **ρ=+0.908**. The predicted outcome held: **the read gate
adds nothing; the −0.44 headline is a two-reference-aligner (SEDEF lastz vs minimap2 asm20) concordance,
carrying zero read information.** This settles the strongest axis against A2's read-content claim.

**Only true closure of L1:** PSV alleles called from a read-only consensus/POA of each copy's reads,
with NO asm20 reference seeding (or a test of read disagreement at columns where the reference copies
*agree* — a class the current code can never emit). That is not in production.

---

## VERDICT

A2 delivers a **bounded, partial** increment over the prior overlay — genuine in CODE + INTERVAL
(different aligner, graded gradient, exonic read-masked subset; ρ=−0.443, p=9.1e-09, confound- and
permutation-robust) — but it does **NOT** deliver the read-CONTENT cross-modal corroboration of O1 it
was tasked with: the load-bearing axis is read-GATED but reference-SEEDED (asm20 self-align = same
substrate family as SEDEF), and the only genuinely read-base axis (assignment rate) is non-significant
(−0.138, p=0.087). Two of three adversarial lenses correctly flag this as MAJOR; circularity is
relocated, not escaped.

**Must-fix status:** (1) re-headlined to "second-aligner reference concordance over a read-masked exonic
subset; cross-modal in CODE+INTERVAL only — divergence content is shared reference substrate" ✅; (2) the
decisive ungated-PSV break-test is **DONE and confirms the read gate is inert** (ρ −0.443 → −0.420 ungated)
✅; (3) the read-base axis (assignment rate) is NULL (−0.138), so no read-content corroboration holds ✅;
(4) membership is reference-selected ✅. **Final position:** both O1 overlay attempts (POA-substring,
`GENOME_RNA_OVERLAY.md`; and this PSV-resolvability one) reduce to shared reference substrate. The only
fully read-discovered closure — PSV alleles from a read-consensus POA with NO asm20 seeding — is not in
production; the genuinely read-independent O1 evidence we *do* have is (a) recovery of known families from
the genome alone (`genome_family_def`) and (b) structural witnesses like APOBEC3D/F, not either overlay.
