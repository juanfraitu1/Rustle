# Row 1813: the chain score without the tie instability (pre-registered 2026-11-19)

**Problem (row 1812).** Chain score r1769 turns each feature into nine indicators "above the sample's own k-th decile". Where a
decile falls on a value many chains share, a small change of the chain set moves it across the tie and flips the indicator of
every tied chain. Seen: mouse ES, extension off + primary seeds, 16,480 chains dropped instead of 5,582 for a `PREFIX.gtf` that
lost 2,919 chains; OR6737 6,729 instead of 3,518. Reproduced on a rerun. User (2026-11-19): fix the filter.

**Change, fixed now.** One change to the design, everything else as r1769 (same 23 features, same extras, same 11 fit samples and
feature tables of workspace `r1769_multisample_rules`, same fit routine): a chain tied with others gets, for decile k, the share of
its tie's rank interval that lies above k/10 instead of 0 or 1 (untied chains: unchanged 0 or 1). The indicator then moves
continuously when the chain set changes. Model `tools/data/chainscore_r1813.json`, design name `decile9t`.
**Cut, rule fixed now:** the cut (grid 0.001) at which the 11 fit samples together keep the same number of correct default chains
as r1769 at its cut 0.02 (closest from above).

**Data.** Feature tables of `PREFIX.gtf` for arms `base` and `ext0_seed` of row 1811 on seven datasets (fly, mouse, zebrafish
simulated; mouse ES, mouse heart; OR6737, A119b), jobs 4482852-4482858, made by the unchanged feature code. Registered before any
of these tables or any refit was looked at.

**Bars.**
- **T1, stable.** On each of the seven datasets: chains dropped under `ext0_seed` differ from chains dropped under `base` by at
  most 25% of the latter (r1769: +195% mouse ES, +91% OR6737, by the same measure).
- **T2, not worse where it was fitted.** On each of the 11 fit samples and on HG002 (unseen), at the cut above: precision of the
  kept default chains at least r1769's minus 0.005, and correct chains kept at least 0.99 of r1769's.
- **T3, not worse on new data.** On each of the seven datasets, arm `ext0_seed` (today's default): correct chains kept at least
  0.99 of what r1769 keeps on arm `base` scaled by the correct chains present, i.e. kept share of correct chains not more than
  0.01 below r1769's share on `base`; precision shown.
**Decision.** T1, T2 and T3 pass: r1813 is put to the user as the default model. Otherwise the failed bars are reported; frozen
thresholds on the raw feature scale are the registered second attempt (they lose the adjustment to the sample's depth).

**Limits.** The fit data are the 11 primate samples at the defaults of 2026-11-16; the filter stays a primate-trained score (row
1793: its features reward truncated reads). T1 is one perturbation of the chain set, not a proof of stability.

## Result of the first attempt (2026-11-19): FAILS T1, T2 and T3

Model `r1813_filter_ties/chainscore_r1813.json` (design `decile9t`, cut 0.019 by the rule; the refit of the old design reproduces
the shipped r1769 weights to 2e-10). Tables `refit.txt`, `bars.txt`.
- **T2 FAIL.** Precision of the kept chains falls by 0.019 to 0.025 on A119b, PPA fibroblast, PPY, PAB and by 0.010 on HG002
  (it rises by 0.006 to 0.012 on the six others); correct chains kept 0.994 to 1.010 of r1769's.
- **T1 FAIL** on fly (-35%), mouse ES (+52%), A119b (+59%); pass on four. r1769 by the same measure fails on five of seven
  (fly -46%, zebrafish -25%, mouse ES +195%, OR6737 +91%, A119b -44%).
- **T3 FAIL** on mouse simulated, mouse ES, mouse heart (kept share of correct chains 0.915 / 0.936 / 0.934 against r1769's
  0.943 / 0.971 / 0.970 on `base`).
Splitting ties halves the jump on mouse ES but does not remove it: the deciles and the standardisation both move with the chain set.
**A flaw of T1 as registered, seen in these numbers:** the chains that extension added are chains the filter drops, so fewer dropped
chains under `ext0_seed` is expected and is not instability (A119b r1769: 43,157 -> 24,051).

## Second attempt, as registered (frozen thresholds), details fixed before it ran
Design `frozen9`: columns standardised by the pooled mean and sd of the 11 fit samples and cut at their pooled deciles, the
shortest intron likewise; nothing is taken from the sample being filtered. Same features, fit samples, fit routine and cut rule.
Bars T2 and T3 as above. T1 as above, and **added before the run because of the flaw above (T1b):** of the chains present in
both arms, the share whose keep / drop decision differs is at most 1% on each of the seven datasets (shown for r1769 too).

## Result of the second attempt (2026-11-19): STABLE (T1b PASS), T2 and T3 narrowly FAIL; not a default by the rule

Model `tools/data/chainscore_r1813.json` (design `frozen9`, cut 0.018 by the rule), opt-in: `RUSTLE_CHAINSCORE_MODEL=r1813`.
Tables `refit2.txt`, `bars2.txt`; `chainscore.py apply` reproduces the table's count on fly (115 chains dropped).

| dataset | decision differs between arms, r1769 | r1813 | r1813 on today's default: chains dropped, correct kept share, precision before -> after | r1769 there |
|---|---|---|---|---|
| fly simulated | 4.1% | 0.01% | 115, 0.998, 0.926 -> 0.938 | 748, 0.940, 0.964 |
| mouse simulated | 2.8% | 0.05% | 713, 0.990, 0.903 -> 0.915 | 2,609, 0.950, 0.937 |
| zebrafish simulated | 3.2% | 0.08% | 487, 0.993, 0.813 -> 0.834 | 2,939, 0.910, 0.913 |
| mouse ES | 20.5% | 0.13% | 2,056, 0.989, 0.331 -> 0.340 | 16,480, 0.895, 0.424 |
| mouse heart | 1.8% | 0.03% | 598, 0.994, 0.546 -> 0.557 | 2,253, 0.974, 0.586 |
| OR6737 | 4.8% | 0.24% | 6,380, 0.997, 0.321 -> 0.348 | 6,729, 0.995, 0.349 |
| A119b | 6.1% | 0.74% | 61,844, 0.981, 0.158 -> 0.207 | 24,051, 0.995, 0.174 |

- **T1b PASS** on all seven (r1769 fails on all seven: between 1.8% and 20.5% of the chains present in both arms change sides).
- **T2 FAIL** on 2 of 12: PPA fibroblast and PPY lose 0.010 precision (limit 0.005) while keeping 0.5% more correct chains; the
  other ten pass (HG002 unseen: +0.011 precision, 0.994 of the correct chains).
- **T3 FAIL** on A119b only: 0.981 of the correct chains kept against the limit 0.983 (precision 0.207 against r1769's 0.174).
- **T1 FAIL** (fly, zebrafish, mouse ES), the measure with the flaw named above.
**What the frozen model is.** A chain's score depends on its own features only. With thresholds fixed on the primate samples'
scale it filters hard where the library is deep (A119b) and hardly at all on the non-primate datasets, where it keeps 0.99 of
the correct chains and raises precision by about 0.01 to 0.02. r1769's higher precision there came with its instability.
By the registered decision r1813 is not put forward as the default on its own; the user rules.
