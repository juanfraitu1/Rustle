# The length-asymmetry blind spot: real, characterised, and correctly rejected

Run 2026-09-20, following §6t5. Held-out chr2/chr8/chr10, Soto S1C truth, shipped edge set obtained
with `mcl_families --dump-graph` (identity ≥ 0.7, cov_longer ≥ 0.3, exon conjunct at 0.60 already
applied). Tool: ad-hoc analysis over `bench/pair_metric_sweep.py`'s loader.

## First, a correction to my own §6t5 framing

§6t5 reported that symmetric metrics recover **0%** of asymmetric true pairs. That was true of the
simplified arms (one global threshold, connected components). **It is not true of the shipped rule:**

| true pairs | present in the SHIPPED edge set |
|---|---|
| similar length (179) | 173 = **96.6%** |
| asymmetric, min/max ≤ 0.5 (59) | 38 = **64.4%** |

So the real gap is **21 pairs**, and asymmetric pairs are ~10× more likely to be dropped (35.6% vs
3.4%) — a genuine bias, but an order of magnitude smaller than "0%".

## What the 21 missed pairs are

Near-complete containments of a small gene inside a large one: **18 of 21 have containment ≥ 0.94**,
and 16 of 21 involve a pseudogene. Examples: `ANAPC1 ~ ANAPC1P3` (containment 0.995, length ratio
0.020), `PTPN20 ~ PTPN20CP` (0.995), `CTSLP3 ~ CTSLP1/4/6` (≈ 0.98), `SHLD2P1 ~ RHEBP1` (1.000). They
fail `cov_longer ≥ 0.30` because the fragment covers 2–30% of the *long* gene while covering ~100% of
*itself*.

## The obvious fix fails, and so does every pairwise rescue

Adding an OR-admission on containment, scored **only on pairs the shipped rule currently rejects**:

| containment ≥ | TRUE admitted | FALSE admitted | precision |
|---|---|---|---|
| 0.99 | 7 | 29 | 0.194 |
| **0.90** | **19** | **57** | **0.250** |
| 0.70 | 23 | 62 | 0.271 |
| 0.30 | 24 | 79 | 0.233 |

⛔ **Roughly 3 false edges per true one at every threshold** — the curve is flat, so there is no
operating point. This is r293 in a new guise: *"the domain-sharer range ENCLOSES the paralog range."*

And no second pairwise signal rescues it. Among the rejected pairs with containment ≥ 0.90:

| signal | TRUE median (n=19) | FALSE median (n=57) |
|---|---|---|
| **alignment identity** | **0.987** | **0.989** ← false pairs are marginally *more* identical |
| matched bp | 1,342 | 5,150 |
| shorter gene length | 1,345 | 5,238 |

An identity cut on top of containment moves precision from 0.250 to at best **0.268**. Length *does*
differ, but in the unusable direction (the false pairs are the long ones), and a length cut would be
fitted to this sample, not principled.

## Conclusion: the rejection is correct, not a defect

A short sequence sitting fully inside a long one at ~99% identity is **equally consistent** with "a real
duplicate" and "a shared repeat, domain, or unrelated processed pseudogene". The pairwise view does not
contain the information needed to tell them apart — which is precisely why the **exon conjunct** exists,
and these pairs fail it.

⭐ **The shipped rule is making the right trade: it gives up 21 true pairs to avoid 57+ false ones.**
The blind spot is the price of precision, and on this evidence the price is strongly favourable.
§6t3 measured that precision is where the definition's value sits (the conjunct buys +0.136); spending
it to recover 21 asymmetric pairs would be a bad exchange.

⭐ **This closes the asymmetry thread.** It is not a metric problem (§6t5: Ochiai, Dice and Jaccard all
fail), not a threshold problem (the precision curve is flat), and not an identity problem. Anything that
separates these pairs has to come from **structure** — matching intron chains, or exon-to-exon evidence,
i.e. the conjunct the rule already applies — not from a better pairwise scalar.

## Limits

- 76 rejected containment-≥-0.90 pairs with both endpoints Soto-labelled, over three chromosomes.
- "FALSE" means "different Soto family". Soto's universe is not exhaustive, so some of the 57 could be
  real relationships Soto does not record; that would raise the precision figure but not flatten the
  curve, which is what kills the rule.
