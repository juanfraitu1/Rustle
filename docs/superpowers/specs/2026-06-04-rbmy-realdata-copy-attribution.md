# Real-data copy attribution on RBMY — held-out PSV concordance (HONEST NEGATIVE)

**Date:** 2026-06-04. **Verdict:** the synthetic copy-attribution capability does **NOT** transfer
to the one real tandem array tested. Recorded as the honest negative the audit's "RBMY accuracy
UNVALIDATED" gate demanded. Analysis: `bench/paralog_secondary_scan/rbmy_heldout_concordance.py`.

## The question + why it's hard
Real RBMY1 reads (6-copy tandem array, chrY NC_073248.2 ~19.60–19.73 Mb) have **no per-read
source-copy label** and there is **no external per-copy annotation**. So attribution accuracy can't
be measured directly. Leakage-free substitute: a read genuinely from copy X carries X's paralogous
sequence variants (PSVs) **consistently** along its length. Test = **held-out PSV concordance**:
align each read to each copy separately (read-coordinate match vectors via pysam `aligned_pairs`),
take discriminating positions (copies disagree), split them into random disjoint halves, and check
whether the two halves pick the same best copy. Identifiability gated on the **normalized margin**
(net discriminating matches per position). (Method audited: an odd/even split manufactures a parity
artifact — spuriously perfect on one near-tie copy, below-chance on its sister — so averaged
**random-half** is used.)

## Result 1 — real RBMY reads are identifiable for only ONE of six copies
Overall held-out concordance 0.432 vs chance 0.167. Per copy:

| copy | gene | n reads | concordance | norm-margin | identifiable |
|---|---|---|---|---|---|
| c0 | LOC129530243 | 8 | 1.000 | **0.442** | **YES** |
| c1 | LOC129530239 | 10 | 0.966 | 0.040 | no (weak consistent, below gate) |
| c2 | LOC129530240 | 19 | 0.503 | 0.030 | no (near-tie) |
| c3 | LOC129530238 | 1 | — | — | noise |
| c4 | LOC129530241 | 45 | 0.207 | 0.037 | no (≈chance — the dominant pool) |
| c5 | LOC129530242 | 2 | 0.090 | 0.044 | no (divergent/starved, 2 reads) |

c0 is genuine, verified: its 8 reads align to all 6 copies (~97–99% to others, 100% to c0), gap
distributed across the whole read (not a private exon), reference c0 vs others NM 361–374. The
other copies are statistical near-ties (~1–2 net discriminating matches over ~48 positions). The
divergent c5 (the "starved copy") has only 2 reads and **is not identifiable** from them.

## Result 2 — the EM's confidence is MIS-CALIBRATED on real RBMY (anti-correlated)
Tandem-mode EM (`RUSTLE_VG_TANDEM=1`), copies mapped by transcript span:

| copy | PSV concordance | capacity_confidence |
|---|---|---|
| c0 (identifiable) | 1.000 | 0.439 |
| c4 (non-id, ≈chance) | 0.207 | **0.728** |
| c5 (non-id, 2 reads) | 0.090 | **0.956** (highest) |

- `capacity_confidence` is **anti-correlated** with real identifiability (Pearson −0.55): highest on
  the *least* identifiable copies (c4, c5) — it tracks coverage/expression, not attribution
  certainty. `copy_confidence` is uncorrelated (peaks at c5 where concordance = 0).
- Per-read (73/87 reads joined by name hash): the EM **collapses 96/110 multimappers onto a single
  sink copy index** and agrees with the PSV best-copy only **~1/72** times, while reporting
  prior-inflated weight gaps (median 0.54) — i.e. it reports confidence and does **not** abstain
  where the data is non-identifiable.

## Honest conclusion
"Copy attribution works on real RBMY" is **TRUE for c0 only** and **false / non-identifiable for the
other five copies**, including the dominant c4 read pool and the divergent c5. The synthetic 100%
decisive-accuracy result does **not** transfer to this real tandem array: the per-read EM does not
attribute reads to their copy of origin (it sinks them), and its scalar confidences are
anti-correlated with the real signal. Any claim of "recovering the starved c5" or resolving the c4
pool would be the DAZ3 failure mode (confident on non-identifiable data).

This validates the Tier-2 capability as **synthetic-only**; the real-data gap the audit named is now
characterized. Next directions: (1) **DONE — see below**; (2) error-aware per-site scoring may rescue
weakly-identifiable copies (c1); (3) `capacity_confidence` must not be presented as attribution
certainty (still open — it's a coverage proxy).

## Fix applied (2026-06-04): evidence-based per-read confidence
The per-read EM confidence was the **posterior weight gap**, which includes the prior — so a lopsided
RBMY prior drove ambiguous reads onto a sink with a confident-looking gap. Fixed (`vg.rs`
`run_fingerprint_em`): the reported decisiveness now gates on the **pre-prior evidence margin**
(`ev_gap` = fingerprint + junction + NM, EXCLUDING the prior) clearing the same `eff_gap` bar used for
anchoring; exposed as `evidence_gap` + `ev_decisive` columns in the FP-attr TSV. 0-site (structural /
DAZ) reads keep the prior path (DAZ byte-identical). Result — the EM's evidence-decisiveness now
**tracks real PSV identifiability** instead of anti-correlating:

| copy | PSV concordance | EM ev_decisive_frac (after fix) |
|---|---|---|
| c0 (identifiable) | 1.00 | **1.000** |
| c1 / c2 (near-tie) | 0.97 / 0.50 | 0.000 / 0.000 |
| c4 (non-id, dominant pool) | 0.21 | **0.026** (was the confident sink) |

362/458 RBMY reads now abstain. Validated no regression: DAZ3 unchanged (cov 4.04, low_confidence);
synthetic benchmark id 0.97 capability 1.0 + id 1.0 abstains; default headline 95.6/90.5; suite 222/0;
EM-correctness 5/5. (c5 n=2 is statistical noise either way.) **The per-read CONFIDENCE is now honest
on real data.**

## Fix #3 applied (2026-06-04): copy_confidence is now evidence-based attribution certainty
`compute_per_copy_confidence` was the **mean post-EM multimapper weight** — it includes the coverage
prior, so the sink copy of a lopsided family scored high. Replaced with the **fraction of the copy's
ambiguous reads that are evidence-decisive winners** (new per-read `em_ev_decisive` flag = winner AND
pre-prior margin clears eff_gap). The `copy_confidence` GTF attribute is now genuine attribution
certainty; on RBMY the family copies read **low** (c4 sink 0.011, was prior-inflated) — honest
abstention, anti-correlation gone. `capacity_confidence` is left as-is but its gtf comment now states
it is COVERAGE/flow adequacy, **not** attribution (it can be high on a non-identifiable dominant copy
— read `copy_confidence` for attribution). No regression: DAZ3 4.04/low_conf, obj5 1.0/abstain,
headline 95.6/90.5, suite 222/0, EM-correctness 5/5. Residual: **#2** (rescue weakly-identifiable copies like c1).

## Fix #2 investigated (2026-06-04): error-aware FAILS; coherence is the real (future) lever
The proposed #2 lever — error-aware per-site scoring (`RUSTLE_VG_FP_ERROR_AWARE`, uses the read's
`de` instead of a fixed 9:1) — does **NOT** rescue c1 and **HURTS** the genuinely-identifiable c0
on real RBMY:

| copy | PSV concordance | ev_decisive_frac default | ev_decisive_frac error-aware |
|---|---|---|---|
| c0 (identifiable) | 1.00 | 1.000 | **0.250** (worse) |
| c1 (real-weak) | 0.97 | 0.000 | 0.000 (no rescue) |
| c4 (non-id pool) | 0.21 | 0.026 | 0.051 (more false confidence) |

RBMY reads are clean (median `de` 0.0016; 90% clamp up to 0.01), so error-aware *strengthens*
per-site evidence — but it shifts the apportionment (winner), so c0's reads scatter. On synthetic it's
a saturated no-op (id 0.95/0.99 already 1.0). **Keep error-aware default-off.**

**The real distinguishing signal is within-read PSV COHERENCE, not margin.** Per-read held-out
concordance cleanly separates the real-weak c1 from the noise c4 *despite identical margins*:

| copy | mean per-read concordance | mean norm-margin |
|---|---|---|
| c1 (real-weak) | **0.966** | 0.040 |
| c4 (noise) | **0.207** | 0.040 |

So c1 carries a real, consistent signal (its few PSVs agree across the read) that the magnitude-based
`eff_gap` can't credit. A **coherence-gated** decisiveness (reusing `build_read_site_obs`' per-site
per-copy match bits to compute internal per-read PSV agreement, and crediting weak-but-coherent reads)
is the genuine lever. **NOT implemented** — c1 has only ~2-3 PSVs, so coherence over so few sites is
statistically thin and forcing c1 decisive flirts with the DAZ3 over-confidence failure mode. This is
scoped future research with real calibration risk, not a quick win; abstention on c1 remains the
conservative-correct default. The grounded recommendation: if pursued, gate on coherence over a
MINIMUM site count, and validate against the held-out concordance that it rescues c1 without crediting
c4/c2.
