# Pre-registration — can the intron-chain certificate CREATE edges minimap2 never proposes?

**Written 2026-09-22, §6y2, before any chromosome is scored.** User: *"can we try one or more of these
approaches but in edge construction instead?"*

## Why this is the right place to ask, and what it is NOT

§6u4's FN decomposition: of 88 false negatives over 263 truth pairs — **MCL split 45.5% · NO ALIGNMENT
26.1% · identity/cov gate 22.7% · exon conjunct 5.7%** ⇒ **48.8% are edge-construction losses**, and §6o8
put the pairwise recall ceiling at 0.052 (human) with grouping SATURATED.

Edge construction has **two** failure modes and they need different tools:
- **gate rejection (22.7%)** — an alignment exists and is thrown away. Already addressed: `--min-cov-shorter`
  (r1006) admits the fully-contained pairs `cov_longer` evicts. **Not re-tested here.**
- **no alignment (26.1%)** — minimap2 proposes nothing at all. **No alignment-derived signal can reach
  this population by construction.** It is the only place an alignment-free structure can contribute, and
  it is the whole subject of this file.

⚠**r1024 already measured the certificate globally and it was 99.3% redundant** (139 of 140 true hits
already aligned). That is not the question here. The question is whether the residual — the pairs with no
alignment — is recovered at usable precision, measured on **four chromosomes** rather than the single
number r1024 produced.

## The rule under test

A pair is proposed as an edge iff its genes share a **3-intron length shingle at 5% tolerance**
(log-binned), the operating point r1024 measured at precision 0.979 / recall 13.5% overall. Chains come
from the longest annotated transcript, strand-oriented 5'→3'. **No alignment is consulted.**

## Substrates and truth

Development **chr16**; held out **chr2 / chr8 / chr10**, no re-tuning. Truth = the **protein-family
referee** (neutral; independent of both our gate and Soto's SD track). Restricted to genes with ≥3 exons —
⚠**a gene with <3 exons has no 3-shingle and is unreachable by this rule**, which is a ceiling on it, not
a property of the data, and is reported as such.

## Metrics — reported on the NO-ALIGNMENT population only

For each chromosome: truth same-family pairs · how many have **no alignment in the all-vs-all PAF** ·
how many of those the certificate recovers · and the **false-positive count among unaligned pairs**
(different-family pairs the certificate would wrongly join).

## The bar — committed now

| outcome | verdict |
|---|---|
| on held-out, recovers **≥10%** of the no-alignment truth pairs at **precision ≥0.70** among unaligned pairs | ⭐ **REAL EDGE SOURCE** — implement it |
| recovers ≥10% at precision 0.40–0.70 | ⚠ **FLAG ONLY** — usable as a certificate beside an edge, not as one |
| **<10% recovered, or precision <0.40** | ⛔ **NO** — the no-alignment population is out of reach of intron structure |

**Predicted, before looking — ⛔.** r1024's single-chromosome read of this exact population was **1 true
against 3 false (precision 0.25)**, and the mechanism it exposed is general: intron-length chains only
match when the sequences are similar enough that the aligner already found them, so the residual should be
dominated by chance collisions. If the prediction is wrong and precision clears 0.70 on held-out, that is a
genuine alignment-free edge source and the first thing this session has found that could move §6o8's
ceiling.

I will not change the rule, the operating point, the substrates, the truth or the bar after seeing any number.

---

# OUTCOME (2026-09-22) — ⛔ **NO. 0.1% of the target at precision 0.250. Held-out NOT spent.**

Development, human chr16, protein-family referee:

| | |
|---|---|
| referee genes | 306 (256 have a ≥3-intron chain) |
| truth same-family pairs | 1,474 |
| **of which NO ALIGNMENT — the target population** | **1,195 (81.1%)** |
| reachable by the rule (both genes have a chain) | 771 = 64.5% of target |
| ⭐**recovered by the certificate** | **1 = 0.1% of target** (0.13% of reachable) |
| false joins among unaligned pairs | 3 |
| ⭐**precision on unaligned pairs** | **0.250** |

⛔**Two orders of magnitude short of the bar on recall and well under it on precision.** Per the §6w6
precedent, **gorilla-style held-out substrates are not spent on a rule that fails development this
decisively** — chr2/chr8/chr10 were not run.

⭐**The ceiling is not chain availability.** 64.5% of the no-alignment truth pairs have a usable ≥3-intron
chain on both sides; the chains simply do not match. So the failure is not "most genes are too short for a
trie" — it is that **intron-length structure does not survive the divergence that already defeated the
aligner.** r1024's mechanism, confirmed on 771 pairs instead of one.

## ⚠ The number worth keeping is the denominator, not the verdict

**81.1% of referee same-family pairs have NO alignment at all** (1,195 of 1,474). That is a much starker
statement of §6o8's recall ceiling than the ceiling itself, and it independently corroborates §6t1's
divergence cliff from a different direction: §6t1 measured 69% of pairs below protein identity 0.60 where
RNA does not align; this measures 81% with no alignment, on a different population and a different
statistic. ⚠**It also bounds every structural idea in this session**: poset, trie, VG, junction strings and
DP chaining are all computed FROM alignments, so none of them can address four fifths of the truth pairs.

⭐**Scope, not defeat**: §6t1 already established there is no ceiling where the thesis operates (NPIP 90.6%,
FAM90A 100%) — the 81% is dominated by ancient protein families the method deliberately does not target.
The right response is to quote the scope, not to chase the 81%.

## Prediction scorecard

Predicted ⛔ with precision ~0.25 by the r1024 mechanism. **Correct, and the measured precision is 0.250 —
exactly r1024's single-chromosome value, now on a 771-pair reachable population.** The one thing I did not
predict is that reachability would be as high as 64.5%; I had assumed short genes would be the binding
limit, and they are not.

> **Generator (2026-09-22 consolidation):** `python3 bench/edge_probes.py intron-chain ...` — the original `bench/intron_chain_edges.py` was folded in verbatim (§6z2); the register rows above cite this file.
