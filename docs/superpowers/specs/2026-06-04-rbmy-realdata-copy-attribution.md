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
characterized, not closed. Concrete next directions (not yet done): (1) the EM should abstain where
the within-family per-read margin is ≈chance (it currently collapses-with-confidence); (2) error-
aware per-site scoring may rescue weakly-identifiable copies (c1); (3) `capacity_confidence` must not
be presented as attribution certainty — a margin-based per-read confidence is needed.
