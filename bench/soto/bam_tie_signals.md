# BAM signals for multimapper ties: what AS can't do, what `de` fixes, and what's still unused

**Date:** 2026-07-28. Context: the advisor's objection that **`AS` is a coin toss for very similar copies**.
Rustle already switched the conflict-tie criterion to `de` (keeping `AS` reported alongside). This note
inventories what else the BAM carries, measures whether it helps, and states the ceiling.

## 1. The objection, restated precisely

`AS` is the alignment score of **one** alignment. It contains **no reference to the alternative placement**.
For two ~99%-identical copies the two alignments score within a few points of each other, so ranking by `AS`
is close to arbitrary — exactly the advisor's point. `de` (gap-compressed per-base divergence) is a genuine
improvement because it is *continuous* rather than a coarse integer, but it shares the structural flaw: it
is still a property of a single alignment, not a comparison.

## 2. Tags actually present (minimap2 `-ax splice:hq -N 50 -p 0.1`, verified on `soto_regions.bam`)

| tag | meaning | used by Rustle? |
|---|---|---|
| `AS` | alignment score of this alignment | reported |
| `de` | gap-compressed per-base divergence | **yes — the tie criterion** |
| `NM` | edit distance | — |
| `ms` | DP score of the max-scoring segment | — |
| **`s1`** | **chaining score of the best chain** | **no** |
| **`s2`** | **chaining score of the best SECONDARY chain** | **no** |
| **`rl`** | **length of repetitive seeds in the read** | **no** |
| `cm` | number of minimizers on the chain | — |
| `SA` | supplementary alignments (chimeric reads; ~3% of reads) | — |
| `ts` | transcript strand from splice sites | — |
| `nn` | ambiguous bases | — |

## 3. Measurement

Primary alignments only (`-F 2308`). `margin = (s1 − s2) / s1` = how much better the best chain is than the
runner-up. `repeat-frac = rl / read length`.

| locus | n | MAPQ0 | median margin | median `de` | median repeat-frac |
|---|---:|---:|---:|---:|---:|
| **NCF1** (K=0-hard) | 44 | 2 | **0.026** | 0.0045 | **0.112** |
| SRGAP2 (resolved) | 1328 | 6 | 0.036 | 0.0017 | 0.021 |
| PMS2P1 (resolved) | 1142 | 0 | **0.676** | 0.0016 | 0.022 |

**Two useful, unused signals:**

- **`s2` → the margin.** PMS2P1 (0.676: best chain 68% better than the alternative) vs NCF1 (0.026: 2.6%)
  is precisely the discrimination `AS` cannot express, because it is the *comparison* the advisor's critique
  says is missing.
- **`rl` → repeat-driven placement.** NCF1 reads are **11%** repetitive seeds vs **2%** at resolved loci
  (~5×). This is an axis orthogonal to both `AS` and `de`: it says "this placement rests on repeat anchors",
  which no score or divergence encodes.

## 4. A hypothesis that was REFUTED

Expectation: since MAPQ is derived from `s1`/`s2` and saturates at 0, the raw margin should retain
information where MAPQ has thrown it away — a "continuous MAPQ".

**It does not.** Among MAPQ = 0 reads the margin is ~0 and barely varies (min −0.006, median 0.001,
max 0.011). When minimap2 reports MAPQ 0, `s1` really does equal `s2`; the saturation is faithful, not lossy.
⚠Weak evidence: only 2 and 6 MAPQ-0 reads in these loci — suggestive, not settled.

**Second caveat:** SRGAP2 has a *low* margin (0.036) yet resolves cleanly. So the margin grades **chain
ambiguity**, which correlates with but does not determine **resolvability** — a chain-ambiguous read can
still carry a decisive PSV.

## 5. The hierarchy — and its ceiling

| signal | what it measures | tie-detection |
|---|---|---|
| `AS` | score of one alignment, no comparison | weakest — ties trivially |
| `de` | continuous divergence of one alignment | better (current criterion) |
| `s1`/`s2` | **explicit margin to the runner-up** | better still |
| `rl` | is the placement repeat-driven? | orthogonal axis |

**All four are scalar summaries of an alignment, and at true K=0 every scalar ties by construction** — that
is what K=0 *means*. So `s2` and `rl` improve **tie DETECTION** (and let reads be triaged by *type* of
ambiguity: score-tie vs repeat-driven vs chimeric), but they cannot **BREAK** a tie.

Tie-breaking requires **positional** evidence — which positions differ, which allele, at what base quality —
i.e. the PSV + junction likelihood in `copy_assign`. This is the correct division of labour and it is
already how the engine is built:

> The advisor's critique retires `AS` as a **tie-breaker**; `de` is a better **tie-detector**; `s2` + `rl`
> would be better still; and tie-**breaking** must live in the PSV layer regardless of which scalar is used.

## 5b. Reconciliation with the earlier tag-dig (IMPORTANT qualifier)

A previous analysis (`conflict_criterion_bakeoff`) already dug through these tags and concluded:
**"ms/cm/s1 are length-confounded — AS is too; `as_per_base` divides it out"**, and explicitly decided
**"Do NOT port nintron/ms/cm/s1/s2/ts"**. Two points reconcile that with §3-§4 here:

1. **The prior decision was about RAW scores, and it is correct.** `s1` alone scales with read length exactly
   like `AS`, so porting it as a raw score adds nothing. What §3 uses is the **ratio** `(s1 − s2)/s1`, which
   divides length out by construction — the same normalisation `as_per_base` applies to `AS`. The
   *comparison* is the new content, not the score.

2. ⚠**But the margin inherits the PARTIAL-RUNNER-UP confound, and this is serious.** The same prior work
   found that real Iso-Seq secondaries are often *partial*, which "manufactures a huge spurious margin" —
   raw `as_margin` on certified-unassignable reads had median **202** (vs Eichler's AS≥10 rule), so `AS≥10`
   would confidently assign 37/38 reads our gate certifies as UNASSIGNABLE. **`s2` is the best SECONDARY
   CHAIN's score, so a partial runner-up depresses `s2` and inflates `(s1 − s2)/s1` in exactly the same way.**
   The PMS2P1 0.676 vs NCF1 0.026 contrast in §3 may therefore partly reflect *partial vs full runner-ups*
   rather than genuine placement confidence.

   **Consequence:** `s2` must NOT be adopted as a tie-detector without first normalising for the runner-up's
   aligned length (a per-base margin), and validating on the simulation where K=0 reads have `as_margin`
   exactly 0. Until that is done, treat §3's margin numbers as *descriptive*, not as a validated criterion.

This does not affect `rl` (repetitive-seed fraction), which is a read-intrinsic property with no
runner-up dependence, or `SA`, which is a factual list of the read's other alignments.

## 6. Recommendations (not yet implemented)

1. **`s2` — only after length-normalising the runner-up (see §5b).** The margin is the comparison `AS`/`de`
   lack, but it inherits the partial-runner-up confound that already discredited raw `as_margin`. Required
   first: a per-base margin, validated on the sim (K=0 reads must give margin ~0). Do not ship the raw ratio.
2. **`rl`/read-length — the SAFEST of the three** (read-intrinsic, no runner-up dependence). Add as an
   ambiguity-type annotation. A read that is 11% repetitive seeds is
   ambiguous for a *different reason* than one with a genuine score tie; the abstention record should say which.
3. **Use `SA` for mis-chain detection.** ~3% of reads carry supplementary alignments; a read split across
   paralogous loci is exactly the NCF1 mis-chain signature, and `SA` names the partner locus directly.
4. **Do not expect any of these to raise copy-assignment recall** — see §5. They improve the *honesty and
   granularity of abstention*, not the identifiability limit.

**Reproduce:** `samtools view -F2308 <bam> <locus>` and parse `s1`, `s2`, `de`, `rl`.

---

## 7. TSS/TES as a copy discriminator — TESTED AND REJECTED (2026-07-28)

**The blind spot is real.** The exon-sum is TSS/TES-blind by construction: `pass1_skeletons_robust` groups
reads by **exact intron chain only**, and sets boundaries to the *k*-th smallest start / *k*-th largest end
(`k = min_terminal_support`, default 2) — a robustness quantile, not a TSS call. No TSS/TES/polyA-aware logic
exists anywhere in the assembly or family path (the only polyA fields live in `phasing.rs`).

This matters *a priori* because IsoSeq FLNC reads are selected for carrying both the 5' primer and the 3'
polyA, so their boundaries approximate **real** TSS/TES rather than random truncation. The hypothesis was
that copies might use different promoters — a **positional, RNA-only** discriminator that could break ties
where no PSV exists (DNA cannot tell you which promoter is used).

**Test (common coordinate frame).** Raw offsets within each copy are not comparable (copies differ in
length/strand), so for every family with ≥2 copies at ≥40 reads: align copy B's genomic span onto copy A's
(`minimap2 -a -x asm20`), build a B→A position map from the CIGAR (reverse-strand aware), project each read's
5' end into A's frame, and compare median shift against **within-copy MAD** (5' degradation is real even in
FLNC — the *k*=2 quantile exists to absorb exactly that noise). Script: `tss_common_frame.sh`.

**Result — 51 families:** **35 same TSS (69%)**, 9 ambiguous, 7 "distinct" (14%).

The 7 "distinct" do not survive inspection:

- **2 are trivial in magnitude** — GOLGA6L1/GOLGA6L22 shift **12 bp**, NPIPB6/NPIPB9 shift **9 bp**. They
  clear the ratio test only because within-copy MAD is ~1 bp; the exon-sums are effectively identical.
- **5 are non-equivalent copy pairs, not alternative promoters** — TRIM73/**NSUN5P1** (15 kb),
  GOLGA8A/**UBE2Q2P2** (24 kb), SHLD2/SHLD2P3 (72 kb), WASH8P/WASHC1, AC110079.1/AL669831.1. These are
  *different genes* sharing a duplicated block inside one Soto family region, so the projection lands far
  away. This is the §2 over-merge resurfacing as a spurious TSS shift — not a biological signal.

**Conclusion: sibling copies overwhelmingly share their TSS.** The exon-sum's boundary collapse is
*correct*, and TSS/TES is **not** a usable copy discriminator here. Do not build TSS-aware logic for
copy assignment or add a "boundary conflict" term to χ(H).

⚠**Limits of this test:** the median 5' end is a crude summary — a *bimodal* TSS, or the same peaks used in
different proportions, produces no median shift. The stronger distributional test was subsequently run: §8.

**Residual value:** the 5' end distributions are still informative for *biology* (genuine alternative
promoter usage) and possibly for making exon-sum boundaries more consistent — which feeds the `min_coverage`
floor implicated in family fragmentation (§7 of `merge_quality_analysis.md`). Just not for telling copies apart.


---

## 8. Distributional TSS test — the answer is still NO, but only because of a control (2026-07-28)

§7's median test is weak by construction. The stronger test (`tss_distribution_test.py`): compare the full
5'-end **distributions** in a common frame via **Wasserstein-1 (earth-mover)**, against an **empirical null**
= EMD between two random halves of the *same* copy (200 permutations), which captures sampling noise *and*
each copy's own 5'-end heterogeneity — exactly the noise the `k=2` boundary quantile is built to absorb.
It also **rejects non-equivalent pairs** (B→A alignment covering <50% of B), which caused 5 of §7's 7 false
"distinct" calls.

**Raw result: 28/40 homologous pairs "DISTINCT" (11 skipped as non-equivalent).** That REVERSES §7.

⚠**But the raw result is an artifact, caught by a control.** A real TSS is a *sharp peak*; if the 5' ends are
scattered across the locus they are not a promoter at all — they are wherever coverage happened to fall.
Scoring **sharpness** (fraction of 5' ends in the densest 400 bp window) in both copies:

| of the 28 "DISTINCT" calls | n | interpretation |
|---|---:|---|
| sharp peak in **both** copies (≥0.30) | **3** | plausible genuine TSS difference |
| **broad** scatter (0.06–0.29) | **25** | **differential partial COVERAGE, not a promoter** |

The three survivors are ID_8 (PMS2P1/PMS2P7, sharpness 0.94/0.72), ID_116 (GOLGA6L1/GOLGA6L22, 0.86/0.57)
and ID_443 (SHLD2/SHLD2P3, 0.52/0.68). The tell on the rest: EMDs of 20–50 kb (ANKRD36B 49 kb, SHLD2 50 kb)
are far too large to be promoter shifts, and within-copy nulls of 1–13 kb confirm the distributions are broad.

**Corrected verdict: 3/40 pairs (7.5%) show a genuine TSS difference** — consistent with §7's conclusion,
now established with a test that *can* detect shape differences and still doesn't find them.

### Answer to "is the exon-sum enough, or should it be extended?"

**The exon-sum is enough.** Extending it to encode TSS/TES would, in ~92% of families, encode **coverage
noise rather than biology** — and coverage noise is precisely what the `k`-th-quantile boundary rule already
exists to suppress. Do not extend the representation, and do not add a boundary-conflict term to χ(H).

**Two caveats kept honest:**
- For the 3 sharp families TSS *is* real and discriminating. A **conditional** extension (encode TSS only
  where the 5' distribution is sharp in all copies) would be defensible — but at 3/40 families the value is low.
- The 25 "broad" cases are not merely noise: they independently corroborate the **differential-coverage**
  finding in `merge_quality_analysis.md` §7 (copies of one family covering different parts of the gene).
  The same phenomenon that fragments families also fakes TSS differences.

## 9. TES is NOT TSS — the 3' end carries ~5x more copy-discriminating signal (2026-07-28)

§8 tested only **5' ends** and concluded the exon-sum needs no TSS encoding (3/40 pairs distinct). That
conclusion was then applied loosely to "transcript boundaries" in general. **It does not transfer to the 3'
end**, and polyadenylation-site choice is a different mechanism from promoter choice, so it never should have
been assumed to. Re-running the same machinery with `--tes` (`bench/soto/tss_distribution_test.py --tes`):

| | distinct pairs | same | BROAD (coverage artifact) |
|---|---:|---:|---:|
| **TES (3')** | **14/42 = 33%** | 12 | 16 |
| TSS (5') | 3/40 = 8% | 12 | 25 |

TES is ~4.7x more often distinct, and its BROAD (coverage-artifact) rate is lower — expected, since IsoSeq
reads are polyA-anchored, making the 3' terminus the most reliably observed boundary in the data, whereas 5'
ends suffer degradation.

**The magnitudes are large, not polyA wobble.** Median EMD **4,855 bp**, max 58.8 kb:

```
ANAPC1  / ANAPC1P2   58,792 bp     GTF2I    / GTF2IP1    10,000 bp
NOTCH2  / NOTCH2NLB  58,592 bp     RGPD2    / RGPD5       6,517 bp
SHLD2   / SHLD2P1    11,113 bp     SPDYE1   / SPDYE2      5,206 bp
```

These are **truncated paralogs terminating early** — NOTCH2NL is a partial duplication of NOTCH2, ANAPC1P2 of
ANAPC1. The 3' end is where that truncation is visible, and it is precisely the copy-distinguishing
information the exon-sum currently discards.

### Answers to the two questions

1. **Is adding boundary variability detrimental?** No — that claim was RETRACTED (`merge_quality_analysis.md`
   §9). The paired test found no detectable effect in either direction (sign test p = 0.69). The snap is off
   for absence of demonstrated benefit, which was measured on the 5' end where there is little signal to find.
2. **Can TES be used in the exon-sum?** **Yes, and it is the most promising untried lever in this document.**
   Encoding the terminal exon's observed 3' extent (rather than the k-th-read quantile) would add real
   differential sequence for 14/42 copy pairs, at a median of ~4.9 kb — orders of magnitude more than the
   sub-50 bp differences the 5' test dismissed.

⚠ **Why this is worth prioritising:** the families showing distinct TES are the ones currently FAILING.
ID_400 (NOTCH2/NOTCH2NL) and ID_395 (RGPD1-4) are both in the set of 11 families dropped by the isoform
requirement in `merge_quality_analysis.md` §17, and ID_14 (LRRC37A), ID_208 (GTF2I), ID_207 (SPDYE) are all
over-merge/split cases elsewhere in these analyses. A lever that adds kilobases of discriminating sequence
exactly where copies are currently indistinguishable is targeted at the known failure mode, not a generic
tweak.

**Not yet implemented.** Suggested next step: extend the terminal exon in `refine_copy_seq`'s exon-sum to the
observed TES where the 3' distribution is sharp (the `sharpness()` control already exists), then re-measure
the §17 family table — the prediction is that some of the 11 isoform-dropped families gain separability.

⚠ Also fixed while running this: `tss_distribution_test.py` crashed at its own summary line with
`NameError: b_` (the BROAD bucket was added to the verdicts but never tallied), so every run since that
verdict was added died before printing the summary.
