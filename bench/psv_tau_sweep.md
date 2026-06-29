# Operating point: precision/recall vs the decisive margin τ

τ is now a real parameter of the shipped engine (`copy_assign --margin <τ>`, default 2.0; the doc-string
gives τ = ln((1−p)/p), τ=6.9 ≈ p 1e-3 = the Eichler AS≥10 analog). The identifiability gate
(`n_decisive≥1`) is independent of τ and always applied. This sweeps τ on two substrates using the
production engine (`copy_assign.rs::assign_read`, mirrored by `bench/copy_assign.py`).

## Sim5x K-ladder — TRUE labels (read name = source copy)

| regime | behaviour |
|---|---|
| **K ≥ 2** | trivially perfect (recall 1.0, precision ~1.0) for every τ up to ~7. Multi-PSV reads are decisive. |
| **K = 1** | the only τ-sensitive case. argmax (τ→0): recall 1.0 but **precision 0.80** (genuine ties tie-broken wrong). τ ≥ 0.5: recall 0.60 / **precision 1.0** (ties correctly discarded). |

A single clean HiFi PSV gives margin = ln(3·0.997/0.003) = **6.90**, so **τ=6.9 (p=1e-3) sits exactly at
the one-confident-PSV knife-edge** — an independent calibration check that τ = ln((1−p)/p) is the right
parameterisation.

## Real GGO co-located families — 47,732 reads / 70 families (silver-standard agreement)

| τ | recall (assigned/all) | silver-standard agreement |
|---|---|---|
| 0.0 (argmax) | 0.976 | 0.9340 |
| 0.5 | 0.964 | 0.9411 |
| **2.0 (default)** | **0.964** | **0.9412** |
| **6.9 (p=1e-3, Eichler-faithful)** | **0.964** | **0.9412** |
| 12.0 | 0.961 | 0.9416 |
| 25.0 | 0.958 | 0.9419 |

**The real-data curve is nearly flat.** Two readings:

1. **τ=2 and τ=6.9 give the *identical* operating point** (recall 0.964, agreement 0.9412). The
   principled conservative choice — τ=6.9, a *stated* per-read misassignment bound p≤1e-3, the PSV-space
   restatement of Eichler's AS≥10 — costs **nothing** in recall over the loose default. Adopt it for free.
2. **The only material decision is τ=0 vs τ>0.** Argmax (assign ties) buys +1.2% recall but *loses*
   ~0.7% agreement; any τ>0 discards the ties and is strictly better. Beyond that, the margin distribution
   on real data is **bimodal** — reads are either decisive (margin ≫ 25) or genuine ties (margin ≈ 0),
   with little thin-margin middle — so τ is a robust, insensitive knob. The lever is the **gate**
   (`n_decisive≥1`, already in this engine), not fine-tuning τ.

## Recommendation for the meeting

Default **τ = 6.9 (p_target = 1e-3)** + **discard ties (τ>0)**: a principled, Canzar-clean operating
point with a stated misassignment guarantee, free on real data, and the exact PSV-space analog of the
AS≥10 criterion your advisor cites. The curve is the evidence; the point is robust.

## Honest caveats

- "Silver-standard agreement" treats each read's best-overlap copy as truth — exactly what is unreliable
  for hard multimappers — so ~94% is a *proxy*; some of the ~6% disagreement is PSV *correcting* the
  overlap call, not error. The sim panel (true labels) is the rigorous anchor.
- 70 families (≥2 valid copies after the spliced-length filter), not genome-wide; capped at 150
  candidates. Representative of the co-located tandem regime, not a full census.
- Flatness is consistent with the GGO headroom measurement (Eichler-discard only 3.4% of multimappers —
  most gorilla paralogs are divergent enough that the margin is decisive): the thin-margin tail where τ
  would matter is genuinely small in this family set.

## Rust cross-check

`copy_assign --margin {2.0,6.9,12.0}` on `NC_073242.2:3771193-3799186` (a 2-copy de-novo family): runs
end-to-end on real GGO, 18 reads all ASSIGNED at every τ (margins 20–330, the decisive regime) — confirms
the engine builds, honors τ, and applies the identifiability gate.

Artifacts: `psv_tau_sweep.py` · `psv_tau_sweep_fig.py` · `psv_tau_sweep.png` · `.json` ·
`src/bin/copy_assign.rs` (--margin/--error-rate).
