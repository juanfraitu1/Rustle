# Are Soto's families a REFINEMENT of ours? (each Soto family inside one of ours)

**Answer: yes at L1/L2, and it is exact for NPIP — all six Soto NPIP-side families sit entirely inside our
single NPIP family. It breaks at L3, where we become finer than they are.**

Test: for every Soto family with >= 2 genes present in our node set, count how many of OUR groups its
members span. Span = 1 means Soto's family is a subset of one of ours, i.e. Soto is a finer level.

## First, a fact that reframes the question

**Soto's families are themselves a COVER, not a partition.** Of 2,333 genes in their universe, **148 belong
to two or more Soto families** — 104 to two, 25 to three, 6 to four, 6 to five, **7 to six**. So "their
families as a subset of ours" cannot be a partition-refinement on their side either; the right statement is
containment of each Soto family in one of ours.

## Development region (NPIP / TBC1D3 lattice levels)

| level | Soto fams (>= 2 genes) | **contained in ONE of ours** | rate | **NPIP-side (6 fams)** | our NPIP parts | our TBC1D3 parts |
|---|---|---|---|---|---|---|
| L1a | 10 | 9 | **90.0%** | **6/6** | 2 | 2 |
| L1b | 10 | 9 | **90.0%** | **6/6** | 2 | 2 |
| **L2** | 10 | 8 | **80.0%** | **6/6** | 2 | 2 |
| L3 | 10 | **2** | **20.0%** | **0/6** | 2 | 4 |

**The user's case confirmed: NPIP is one family for us and six for Soto (ID_149, ID_151–ID_155), and at
L1a/L1b/L2 every one of those six lies entirely inside our NPIP group.** Their families are a subfamily
level of ours, exactly as the nested-lattice story wants.

**L3 is the boundary.** At L3 containment collapses to 2/10 overall and **0/6** on the NPIP side — we cut
finer than Soto does, so we split their families rather than refining them. L3 is past their granularity.

## Held-out region (chr5/7/21, MCL catalogs, 76 / 72 Soto families)

| catalog | Soto fams >= 2 | contained in one of ours | rate |
|---|---|---|---|
| E1 | 76 | 61 | **80.3%** |
| E1S | 72 | 55 | **76.4%** |

Independently reproduces the ~80% containment seen at L2, on a substrate that was never used to develop
any of this.

## Where containment fails — it is a size effect

| Soto family size | contained (E1) |
|---|---|
| <= 3 genes | **53/59 (89.8%)** |
| > 3 genes | **8/17 (47.1%)** |

Small Soto families are nested almost always; large ones are split about half the time. That is consistent
with §6kr's hand-picked result (chr1 all four nested; chr15/17 only CHRFAM7A nested, with GOLGA8 spanning 3
of ours and LRRC37A 4) — the failures there were the large families too.

## Reading

**Soto sits between our L2 and our L3.** Above L3 they are a refinement of us at ~80–90%; at L3 we overtake
them. The residual disagreement is concentrated in their large families, not spread evenly — which is a
much more tractable statement than "we disagree with Soto", and it is the form the nested-lattice
(§0★★★) framing needs.

⚠ The 10-family development table is descriptive (NPIP/TBC1D3 are development families). The chr5/7/21
numbers are the ones to quote: **80.3% / 76.4% on a held-out substrate.**
