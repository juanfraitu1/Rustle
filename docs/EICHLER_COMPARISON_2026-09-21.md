# Eichler's AS-margin rule, computed alongside ours

⭐ **How to emit it:** `copy_assign --eichler-margin 10 --no-as-tied-only`. Two columns are appended
LAST — **`eichler_call`** (`assign`/`discard`) and **`eichler_same_copy`** (1/0/NA) — so every existing
column keeps its position and the schema is **byte-identical when the flag is unset** (verified: a
flag-off run with the new binary is `cmp`-identical to the pre-change output). `bench/eichler_compare.py`
summarises the two columns into the tables below.

⚠ **CORRECTION, 2026-09-21.** The first version of this document reported Eichler assigning **2,536**
reads. That was wrong: the Python tool silently skipped the **1,522** rows whose `as_margin` is `NA` — a
read with *no rival placement*, which has nothing within T and which Eichler therefore **assigns**. The
Rust implementation exposed the bug (4,058 − 2,536 = exactly 1,522). Every count below is the corrected
one; the agreement rate moves 97.4% → **97.1%** and is unaffected in substance.

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
| **0** (aligner indifferent) | 10,894 | **63.7%** |
| 0 < m < 10 | 2,151 | 12.6% |
| ≥ 10 | 2,536 | 14.8% |
| no rival placement (margin NA) | 1,522 | 8.9% |

⭐ **Eichler's rule discards 76.3% of the population by construction** (13,045 of 17,103) — and that discarded
set is exactly O2's subject. Running our comparison on our *gated* output shows this starkly: every one
of those 10,881 reads has margin 0, so Eichler discards 100% of them. The two rules barely share a
subject, and any "agreement" figure has to say which population it is on.

## Where both make a call, we agree with the aligner 97.4%

| | Eichler (T = 10) | ours |
|---|---|---|
| assigns | **4,058** (23.7%) | 657 (3.8%) |
| **both assign** | **512** | |
| **agreement on those 512** | | **497 / 512 = 97.1%** |
| we assign where Eichler discards | | **145** |
| Eichler assigns where we abstain | 3,546 | |

⭐ **This is the reassuring number for the advisor**: where the aligner is confident, our PSV certificate
reproduces its call 97.1% of the time. We are not doing something exotic — we agree with the aligner
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
| **10** | **4,058** | **23.7%** |
| 20 | 1,588 | 10.2% |

The rule is quite sensitive to T — halving from 10 to 5 adds 61% more assignments — which is worth
knowing when "AS ≥ 10" is quoted as though it were a constant of nature.

## What this does and does not establish

- ⭐ **Establishes**: the two methods are complementary, not competing. Eichler's rule answers the easy
  16.3%; ours answers 145 reads inside the 76.3% it discards, and reproduces its answer 97.1% of the
  time where they overlap.
- ⛔ **Does not establish that our 145 extra assignments are correct.** There is no ground truth here.
  Register 690 already measured the honest headroom over this gate at **28%**, not the ~50% a naive
  reading suggests, and that number stands.
- ⚠ The 6× gap in raw assignment count (4,058 vs 657) is not a deficit: 2,028 of Eichler's assignments
  (3,546 of them) are reads we deliberately abstain on, which is the assign-or-abstain discipline the advisor asked for
  (Q8), not a failure to decide.
