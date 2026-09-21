# Pre-registration — intron placeholders: does preserving exon SPACING separate real containments?

**Written 2026-09-20 before any placeholder representation is aligned.** User: *"for some members the
exon sum still lacks introns, can we check if there is any way to add space before and after or put
placeholders for the introns."*

## Why this is the right question to ask now

Three independent results converge on **structure** as the only channel left:

- **§6t7**: among pairs the shipped rule rejects at containment ≥ 0.90 there are **19 TRUE and 57 FALSE**,
  and *no pairwise signal separates them* — containment precision is flat at ~0.25 from 0.30 to 0.99,
  and identity is **0.987 TRUE vs 0.989 FALSE**. Its conclusion: any fix must come from structure.
- **§6t8**: a variation graph's multiplicity channel is only **AUC 0.681**, so the multi-way signal does
  not carry it either. Its explicitly untested channel was **path topology**.
- **r505**: comparing **concatenated exons** with the single-record coverage rule *"penalises CONCATENATED
  EXONS specifically"* — aggregating over records recovered 53 RNA edges and 0 DNA edges. Splicing
  fragments the alignment because concatenation destroys the spacing.

The proposal sits exactly between the two failures: **genomic span** drowns the signal in intron
sequence (and its length asymmetry is what §6t7 is about), while **bare concatenation** destroys the
spacing (r505). A placeholder keeps exon order and approximate geometry without intron sequence.

## The hypothesis, and why it should DISCRIMINATE rather than merely help

> **H:** a true segmental duplicate retains its parent's **intron structure**, so with intron
> placeholders the two align as a colinear chain of exons at matching offsets. A false containment —
> a shared repeat or domain, or a processed (retrotransposed, intron-less) pseudogene of an unrelated
> parent — does not, because it either has no introns or has different ones.

So the placeholder representation is predicted to **separate** the classes, not just to raise everyone's
score. That is a falsifiable prediction and the reason this is worth running.

## The representations — every gene rendered four ways

| arm | representation |
|---|---|
| **G** | genomic gene body (the shipped substrate) — baseline |
| **S** | exons concatenated, introns deleted (r505's substrate) |
| **P50** | exons joined by a **fixed 50 bp** `N` spacer — order preserved, intron length neutralised |
| **Pprop** | exons joined by an `N` spacer of **min(true intron length, 500)** — order and approximate geometry preserved, long introns capped |

Each aligned all-vs-all with `minimap2 -x asm20 -c -X -N 50 -p 0.1`, unchanged from every other arm.

## Population and measure

The **76 pairs the shipped rule rejects at containment ≥ 0.90** on held-out chr2/chr8/chr10, both
endpoints Soto-labelled: **19 TRUE, 57 FALSE** — the identical population as §6t7 and §6t8, so the
numbers are directly comparable.

Primary measure: **AUC** for separating TRUE from FALSE, per representation, using the pair's best
alignment identity and its containment. Reported against the two standing baselines:
**identity on genomic ≈ 0.50 (0.987 vs 0.989, no separation)** and **multiplicity AUC 0.681 (§6t8)**.

## The bar — committed now

| outcome | verdict |
|---|---|
| a placeholder arm reaches **AUC ≥ 0.80** | ⭐⭐ **STRUCTURE CARRIES IT** — the channel §6t7/§6t8 pointed to is real and usable |
| AUC 0.70–0.80 | ⭐ **BETTER THAN MULTIPLICITY** — beats 0.681, worth a pre-registered follow-up |
| AUC 0.60–0.70 | ⚠ **NO BETTER THAN §6t8's multiplicity** |
| AUC < 0.60 | ⛔ **NO** — structure does not separate them either, and the containment problem is closed for good |

Secondary, committed now: **how many of the 19 TRUE pairs align at all** under each representation. If
placeholders make TRUE pairs alignable that were previously unalignable, that is a real gain even if the
AUC bar is missed, and it will be reported either way.

⚠ I will not tune the spacer length after seeing scores — 50 and min(intron, 500) are fixed now. I will
not change the population or the truth.
