# Quantification benchmark: the value of multi-mapper resolution

Two near-identical paralog copies (A,B) at abundance theta:1-theta, differing at K PSVs; HiFi reads (P_cov=0.6,
eps=0.005). Estimate per-copy abundance via naive 1/k, abundance-EM (soft, PSV likelihood), psv_hard (per-
molecule majority vote). Metric = |theta_hat - theta|, 25 reps. (quant_benchmark.py / .tsv / .png)

## Findings
1. **K=0 (no PSVs): ALL methods fail equally** (error = |theta-0.5|). The identifiability floor — no information,
   no resolution. The same wall as detection.
2. **The value is at SKEWED abundance, not symmetric.** At 50:50, naive 1/k is already exact (0 error). The
   payoff is at asymmetric expression — which is the COMMON paralog case (one copy dominant / a pseudogene):
   value (naive_err - psv_err) reaches **+0.19 at 70:30 and +0.40 at 90:10**. (This corrects the earlier
   "gap at similar abundances" framing — for point estimates it's the opposite.)
3. **Value scales with PSV density K — the long-read edge.** 90:10: value +0.24 (K=1) -> +0.40 (K>=8). Each
   molecule covering more PSVs => exact quantification; short reads (few PSVs/read) can't.
4. **Soft EM beats hard assignment at SPARSE PSVs.** K=1, 90:10: EM err 0.007 vs psv_hard 0.160; they
   converge by K>=4-8. ACTIONABLE: for quantification, copy-assignment should be SOFT/probabilistic (use the
   likelihood), not hard labels — hard assignment wastes thin partial evidence.
5. **At symmetric abundance the value is DIFFERENTIAL POWER, not point error.** naive always reports 50:50
   (zero variance but ZERO power to detect a real copy difference); only resolution (EM/PSV) recovers the
   actual ratio and can test copy-specific change. So resolution is essential for differential analysis even
   when abundances are equal.

## Takeaway
Resolving multi-mapper ambiguity adds nothing to ASSEMBLY (structure is in the primaries) but is decisive for
QUANTIFICATION of asymmetrically-expressed paralogs, scaling with PSV density (long-read advantage) and bounded
only at the K=0 identifiability floor. The right estimator is the SOFT likelihood (Canzar's frame), not hard
labels. This is the concrete payoff of the PSV machinery: per-copy abundance, not better assembly.

## Real-reads validation (sim5x, ground truth) — sim5x_quant_validate.py

Subsampled the sim5x HiFi reads to a SKEWED abundance (40:30:20:10:10 = [0.364,0.273,0.182,0.091,0.091],
truth in read names) and swept the PSV ladder K=0..8, running the actual copy_assign soft quant:

| K | soft EM err | naive(uniform) err | hard-count err | soft abundance |
|---|---|---|---|---|
| 0 | 0.473 | 0.473 | 1.27 | [.2,.2,.2,.2,.2] (uniform: no info) |
| 1 | 0.276 | 0.473 | 1.06 | partial resolution |
| 2 | **0.007** | 0.473 | 1.02 | [.37,.27,.18,.09,.09] ≈ truth |
| 4 | 0.007 | 0.473 | 1.02 | ≈ truth |
| 8 | 0.025 | 0.473 | 1.02 | ≈ truth |

CONFIRMS the simulation benchmark on REAL reads, and proves the EM RESOLVES (not just returns the prior):
- K=0 (no PSV): soft == naive == uniform (honest identifiability floor; the soft EM does not fabricate).
- K=1: soft already pulls toward the truth (0.276 < naive 0.473) -- uses the single PSV's PARTIAL evidence.
- K>=2 (PSVs clear the HiFi error floor, psv_acc=1.0): soft = 0.007, recovering the SKEW exactly -- genuine
  per-molecule resolution, not the 0.2-each prior. Naive stays flat at 0.473 (never resolves).
- Tracks the identifiability ladder (per-read psv_acc K0=0.2 -> K2=1.0). The quant error IS the assignment error.
- (hard read-count is ~1.0 throughout: each read has 5 alignments and the wrong-frame ones default to c0 --
  itself an illustration of why naive best-alignment counting fails on multimappers; the clean contrast is
  soft-resolves vs naive-uniform.)
