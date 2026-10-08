# Intron placeholders — tested, and the reason they fail closes the containment question

Run 2026-09-20 against `docs/PREREG_intron_placeholder_2026-09-20.md` (committed `521c2ef9` before any
placeholder was aligned). Tool: `bench/intron_placeholder_edges.py`. Population: the identical
**19 TRUE / 57 FALSE** pairs as §6t7 and §6t8 — those the shipped rule rejects at containment ≥ 0.90 on
held-out chr2/chr8/chr10.

## The representation fails, and makes things worse

Each gene rendered four ways and aligned with the same `minimap2 -x asm20 -c -X -N 50 -p 0.1`:

| representation | TRUE pairs that align at all | FALSE | AUC (identity) |
|---|---|---|---|
| **S** — exons concatenated, no spacer | **3 / 19** | 3 / 57 | 0.556 |
| **P50** — 50 bp `N` spacer between exons | **3 / 19** | 3 / 57 | 0.555 |
| **Pprop** — spacer = min(true intron, 500) | **1 / 19** | 2 / 57 | 0.508 |

⛔ **Pre-registered verdict: NO** (bar was AUC ≥ 0.80; 0.70 to beat §6t8's multiplicity). Worse, the
secondary measure moves the wrong way: **94–95% of TRUE pairs stop aligning entirely.** The exons are
short and fragmented, so minimap2 loses the seeds and chains it had on the genomic span. Preserving
spacing does not rescue that — `Pprop`, which preserves geometry most faithfully, is the worst arm.

## But the underlying intuition was right, and measuring it directly says why

The hypothesis — that exon structure is what distinguishes a real duplicate — is testable without
re-rendering, by asking what fraction of the **genomic** alignment is exonic:

| | TRUE (n=19) | FALSE (n=57) | AUC |
|---|---|---|---|
| exonic fraction, the more exonic side | **0.367** | **0.107** | **0.655** |
| exonic fraction, **both** sides (min) | 0.000 | 0.000 | 0.582 |

⭐ **A 3.4× gap in medians**, where identity gave 0.987 vs 0.989. The signal is real — but at
**AUC 0.655** it lands *below* §6t8's multiplicity (0.681) and inside the pre-registered "no better"
band. It is another ~0.65 signal, the third in a row.

⭐⭐ **The second row is the decisive fact: the median pair has ZERO exonic overlap on at least one
side.** The containment in these pairs is in **intronic or otherwise non-exonic sequence** — a fragment
sitting inside an intron, not inside a coding structure. That is why re-rendering as exons deletes the
shared sequence outright, and it is why the placeholder idea could not have worked.

A two-sided exonic requirement is perfectly precise and nearly empty:

| min exonic fraction ≥ | TRUE | FALSE | precision | recall |
|---|---|---|---|---|
| 0.10 | 3 | **0** | **1.000** | 0.158 |
| 0.30 | 1 | 0 | 1.000 | 0.053 |

## What this settles

⭐ **The shipped exon conjunct is not merely a defensible trade — it is the correct test, and it is
already applying exactly this criterion.** `mcl_families` ships `--exonic-both-sides` on by default, and
these 76 pairs are the ones it rejects. This run shows what it is rejecting: pairs whose shared sequence
is not exonic on both sides. §6t7 concluded the rejection was a favourable trade; this explains the
mechanism.

⭐ **The containment question is now closed from four directions**: pairwise identity (§6t7, no
separation), containment thresholds (§6t7, flat curve), graph multiplicity (§6t8, AUC 0.681), and
exon structure (here, AUC 0.655 one-sided, empty two-sided). Every channel available at the pair level
lands at or below AUC 0.68.

⚠ **What would still be new** is not another pair statistic but **intron-chain identity** — do the two
genes share the same *junctions*, not merely exonic bases. That needs spliced RNA evidence at both loci,
which most of these pseudogene fragments do not have, so it is a question about node construction, not
about edge scoring.

## Limits

- 76 pairs, three chromosomes, Soto truth.
- Exonic fraction uses the transcript with most exonic bases per gene; a different transcript choice
  moves the fractions, not the zero-overlap finding.
- `Pprop` caps spacers at 500 bp; uncapped spacers would reproduce the genomic arm by construction.
