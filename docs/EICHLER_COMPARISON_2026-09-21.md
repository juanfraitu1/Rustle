# Eichler's AS-margin rule, computed alongside ours

Tool: **`bench/eichler_compare.py`**. Substrate: the Y ampliconic genes (O2's own hard case),
`copy_assign --families` with **`--no-as-tied-only`** so the full pre-gate population is visible.

    EICHLER(T):  assign to the best alignment iff no other alignment scores within T units (T = 10)
                 otherwise discard as ambiguous
    OURS:        AS-tied gate -> PSV certificate -> assign-or-abstain, never 1/k

⭐ **No pipeline change was needed.** `copy_assign` already emits `as_best`, `as_second` and `as_margin`
beside our `status`, so both calls come out of one file. The "mode" is a reader, not a new computation.

## The population, which is the whole story

| AS margin | reads | share |
|---|---|---|
| **0** (aligner indifferent) | 10,894 | **69.9%** |
| 0 < m < 10 | 2,151 | 13.8% |
| ≥ 10 | 2,536 | 16.3% |

⭐ **Eichler's rule discards 83.7% of the multi-mapping population by construction** — and that discarded
set is exactly O2's subject. Running our comparison on our *gated* output shows this starkly: every one
of those 10,881 reads has margin 0, so Eichler discards 100% of them. The two rules barely share a
subject, and any "agreement" figure has to say which population it is on.

## Where both make a call, we agree with the aligner 97.4%

| | Eichler (T = 10) | ours |
|---|---|---|
| assigns | **2,536** (16.3%) | 653 (4.2%) |
| **both assign** | **508** | |
| **agreement on those 508** | | **495 / 508 = 97.4%** |
| we assign where Eichler discards | | **145** |
| Eichler assigns where we abstain | 2,028 | |

⭐ **This is the reassuring number for the advisor**: where the aligner is confident, our PSV certificate
reproduces its call 97.4% of the time. We are not doing something exotic — we agree with the aligner
wherever the aligner has an opinion, and we add calls only where it does not.

Per family, the 13 disagreements are not spread evenly:

| family | agree | disagree | |
|---|---|---|---|
| DAZ | 365 | 0 | 100% |
| RBMY | 59 | 0 | 100% |
| CDY | 53 | 0 | 100% |
| BPY | 1 | 0 | 100% |
| **TSPY** | 13 | **5** | 72.2% |
| **PRY** | 4 | **8** | **33.3%** |

⚠ **PRY is a flag, not a result.** On 12 reads we pick a different copy than the aligner two times out of
three. That is either real PSV signal the alignment score misses, or it is us being wrong on a hard
family — and with no ground truth on chrY this run cannot say which. It is the obvious next thing to
adjudicate.

## Threshold sensitivity

| T | Eichler assigns | |
|---|---|---|
| 1 | 4,687 | 30.1% |
| 5 | 4,089 | 26.2% |
| **10** | **2,536** | **16.3%** |
| 20 | 1,588 | 10.2% |

The rule is quite sensitive to T — halving from 10 to 5 adds 61% more assignments — which is worth
knowing when "AS ≥ 10" is quoted as though it were a constant of nature.

## What this does and does not establish

- ⭐ **Establishes**: the two methods are complementary, not competing. Eichler's rule answers the easy
  16.3%; ours answers 145 reads inside the 83.7% it discards, and reproduces its answer 97.4% of the
  time where they overlap.
- ⛔ **Does not establish that our 145 extra assignments are correct.** There is no ground truth here.
  Register 690 already measured the honest headroom over this gate at **28%**, not the ~50% a naive
  reading suggests, and that number stands.
- ⚠ The 4× gap in raw assignment count (2,536 vs 653) is not a deficit: 2,028 of Eichler's assignments
  are reads we deliberately abstain on, which is the assign-or-abstain discipline the advisor asked for
  (Q8), not a failure to decide.
