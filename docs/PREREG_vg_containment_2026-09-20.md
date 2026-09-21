# Pre-registration — does a variation-graph signal separate real containments from repeat/domain ones?

**Written 2026-09-20 before any multiplicity is computed.** User: *"determine if using variation graphs
helps with the containment problem."*

## The problem, as §6t7 left it

The shipped rule rejects 21 asymmetric true pairs. They are near-complete containments, and **no pairwise
signal separates them from the 57 false containments in the same regime**: containment-OR gives precision
0.250 with a FLAT curve from 0.30 to 0.99, and alignment identity is **0.987 TRUE vs 0.989 FALSE**.
§6t7's conclusion: *a short sequence fully inside a long one at ~99% identity is equally consistent with a
real duplicate and with a shared repeat/domain — the pairwise view lacks the information.*

## What a VG would add, stated precisely

A variation graph is a **multi-way** alignment, so it carries one thing a pairwise alignment cannot:
**how many other sequences traverse the same segment.** The hypothesis follows directly:

> **H:** a real duplicate's shared segment is traversed by its family and little else (low multiplicity),
> while a repeat/domain's shared segment is traversed by many unrelated genes (high multiplicity).

⚠ This does **not** require building a VG to test. Node multiplicity is exactly *"how many distinct other
genes align to this interval"*, which the existing all-vs-all PAF already determines. If multiplicity
does not separate the two classes, then the main extra information a VG carries does not separate them
either, and building one cannot help.

## What is already refuted, and why this is not that

- **r472** — *"use a variation graph to DEFINE families (family = one VG, copies = paths)"* — ⛔
  **circular: a VG presupposes its members.** Not proposed here. This uses a graph-derived SIGNAL to
  admit or reject an edge; it does not define the family.
- **r467** — minigraph for the DNA VG of family copies: SV-level, collapses near-identical copies. Not
  used; no graph is built at all.
- **r384** — the library-free VG repeat catalog (minimizer multiplicity) as a **general** separator:
  **AUC 0.686 vs aln_frac's 0.85**, cutting 54% of borderline-real edges, *"defensible only as a targeted
  mult ≥ rule"*. ⚠ **That is the unfavourable prior for this test, and it is a different question**: r384
  asked whether multiplicity separates edges in general; this asks whether it separates the specific
  76-pair containment population where every pairwise signal has already failed. r384's own verdict
  invites exactly a targeted rule.

## Method — frozen

Population: the **76 pairs the shipped rule rejects with containment ≥ 0.90**, both endpoints
Soto-labelled — **19 TRUE** (same family) and **57 FALSE** — on held-out chr2, chr8, chr10.

For each such pair, on the **long** gene's aligned interval, the multiplicity is

    mult = number of DISTINCT other genes on the chromosome whose PAF alignment overlaps
           that interval by >= 50% of the interval

computed from the same all-vs-all PAF, excluding the two genes of the pair.

## The bar — committed now

| outcome | verdict |
|---|---|
| a `mult ≤ k` cut reaches **precision ≥ 0.60** while keeping **≥ 10 of the 19** TRUE pairs | ⭐ **VG HELPS** — the multi-way signal is the missing information; building one is justified |
| precision 0.40–0.60 at ≥ 10 TRUE | ⚠ **PARTIAL** — real signal, not enough to admit edges on |
| best precision < 0.40, or < 10 TRUE retained | ⛔ **NO** — multiplicity does not carry it either, and a VG would not help |

Baseline to beat: **0.250**, the flat precision of containment alone.

⚠ I will report the full precision/recall curve over k, not just the best point, and I will report AUC so
it is comparable to r384's 0.686. I will not change the population, the truth, or the 0.90 containment
floor after seeing the numbers.
