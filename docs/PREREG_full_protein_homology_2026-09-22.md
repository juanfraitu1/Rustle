# Pre-registration — does requiring FULL protein homology rescue the sensitive mode?

**Written 2026-09-22, §6y5, before any floor is swept.** User: *"lets now try to require full protein
homology"* — the single option r1029 left standing.

## Why this is the only remaining form of the idea

§6y4/r1028 measured what r906 could not: of the protein-only edges the sensitive mode ADDS,
**90.0–97.5% have no nucleotide alignment at all**, held-out pessimistic false-merge is **0.900–0.975**,
and the additions are **domain sharing, not paralogy** — different symbol roots in 54–85% of cases
(`ADRA2A~HTR7` GPCR 7TM, `AAMP~CIAO1` WD40, `ACTA2~ACTR1A` actin fold). r1029 named the cause:
**§6ko's rule admits a pair at ≥0.30 coverage of the LONGER protein, and a shared domain easily reaches
0.30.** Raising that floor is the direct fix and the only untested one.

## The sweep

§6ko's rule otherwise unchanged (blastp `-evalue 1e-5`, greedy non-overlapping HSPs by bitscore projected
onto the longer protein, CDS translated in-house). **Only the coverage floor moves**:
**0.30 (r906's shipped value) · 0.50 · 0.70 · 0.80 · 0.90**.

## ⚠⚠ The cost this is expected to have, stated before measuring

Raising the floor removes domain sharers **and genuinely divergent paralogues in the same stroke** — a
real paralogue that has lost a terminal exon, or a partial retrocopy, also covers well under 0.90 of its
longer partner. **So a clean false-merge number at a high floor is not a win by itself**; it must still add
edges the nucleotide graph does not have. Both numbers are reported at every floor and neither is allowed
to stand alone.

## Metrics — identical to §6y4 so the floors are comparable

Per floor, per chromosome: protein-only edges (the additions) · TRUE / FALSE (both Soto-labelled) ·
UNKNOWN corroborated by a sub-threshold nucleotide record · UNKNOWN with no alignment at all ·
**pessimistic false-merge `(FALSE + bare)/all`** · and the **cross-symbol-root fraction** of the bare set
(⚠diagnostic only, r902) to confirm domain sharers are actually being removed rather than merely thinned.

Development **chr16**; held out **chr2 / chr8 / chr10**, no re-tuning, verdict taken there.

## The bar — committed now

| outcome | verdict |
|---|---|
| a floor gives held-out pessimistic false-merge **≤ 2.00%** (§6bt.1's rate) **AND still adds ≥ 50 protein-only edges per chromosome** | ⭐ **ADOPT AS AN OPT-IN SENSITIVE MODE** |
| ≤ 2.00% but adds < 50 edges per chromosome | ⚠ **VACUOUS** — the floor bought cleanliness by emptying the mode |
| no floor reaches ≤ 2.00%, or the cross-root fraction does not fall | ⛔ **NO** — the protein line is closed for good |

**Predicted, before looking — ⚠ VACUOUS.** The mechanism r1024 established should repeat: full-length
protein homology implies enough nucleotide similarity for minimap2 to align the pair, so the surviving
edges should be ones the nucleotide graph already has. I expect the false-merge rate to fall cleanly with
the floor — confirming r1029's diagnosis — while the protein-ONLY count collapses toward zero. If instead a
floor is both clean and non-empty, that is a genuine sensitive mode and the first thing in this line to
earn adoption.

I will not change the rule, the floors, the metrics or the bar after seeing any number.

---

# OUTCOME (2026-09-22) — ⛔ **NO. The floor does not clean the set, and it refutes my own r1029 diagnosis.**

Pessimistic false-merge `(FALSE + bare)/protein-only`, by floor:

| floor | chr16 (dev) | chr2 | chr8 | chr10 | protein-only edges (chr16/2/8/10) |
|---|---|---|---|---|---|
| 0.30 | 0.913 | 0.975 | 0.919 | 0.900 | 1009 / 2652 / 418 / 368 |
| 0.50 | 0.898 | 0.973 | 0.898 | 0.873 | 640 / 2023 / 236 / 205 |
| 0.70 | 0.878 | 0.975 | 0.887 | 0.815 | 353 / 1444 / 141 / 124 |
| 0.80 | 0.897 | 0.977 | 0.869 | 0.802 | 213 / 1257 / 99 / 101 |
| **0.90** | 0.912 | 0.976 | 0.968 | **0.794** | 113 / 991 / 62 / 68 |

⛔**No floor comes near 2.00% on any chromosome** — the range across the whole sweep is **0.79–0.98**. And
the cross-symbol-root fraction of the bare set does **not** fall as required: chr16 **54.3% → 64.1%** and
chr2 **85.3% → 88.9%** (it RISES), chr8 65.1% → 33.3% and chr10 71.9% → 53.7% (falls). Both ⛔ conditions
fire. **Edges are removed roughly uniformly across classes; the floor thins the set without cleaning it.**

## ⛔ This refutes r1029's diagnosis, which was mine

r1029 said §6ko's 0.30 floor *"lets a shared domain qualify"*, implying a higher floor would exclude domain
sharers. **It does not, and the reason is visible in the survivors at 0.90**: `ACTA2~ACTR1A` (377/376 aa),
`CYP17A1~CYP2C8` (508/490 aa), `ADRA2A~HTR7`. **For these families the shared fold IS the whole protein** —
a GPCR is 7TM end to end, a P450 is one fold, an actin is one fold. **A coverage floor cannot separate
"same family" from "same fold" when the fold spans the protein.**

## ⭐⭐ But the protein-only set contains GENUINE families the nucleotide graph misses

The chr10 floor-0.90 survivors are a **mixture**, not noise:

- **genuine paralogues**: `AKR1C1~AKR1C4`, `AKR1C2~AKR1C4`, `AKR1C3~AKR1C4` (the real AKR1C tandem cluster),
  `CALHM1~CALHM3`, `CALML3~CALML5`, `CYP26A1~CYP26C1`, `ARL3~ARL5B` — **with no nucleotide alignment at all
  in our graph.** This is exactly what a sensitive mode is supposed to find.
- **fold sharers**: `ACTA2~ACTR1A`, `CYP17A1~CYP2C8`.

⚠⚠**So the blocker is the TRUTH, not the method.** The "bare" class lumps both together, and the
pessimistic rate therefore counts `AKR1C1~AKR1C4` — a real family — as a candidate false merge. **No truth
on hand separates an ancient paralogue from a fold sharer**: Soto is SD-scoped (recent ≥98% duplications
only, so it labels neither) and the protein referee is circular for a protein arm (§6y3). ⭐**That is why
this line cannot be adopted, and it is a different reason from "the edges are wrong."**

⚠**A discriminator is visible but unusable**: symbol root separates the two cleanly here (`AKR1C1`/`AKR1C4`
share a root; `ACTA2`/`ACTR1A` do not). **It cannot be used** — r902 voided symbol-root truth, and it would
not transfer to the thesis's own substrate, where gorilla NPIP copies are annotated as *"titin-like"* /
*"NACA-like"* LOCs with no shared root at all.

## Verdict

⛔ **The protein line is closed as specified.** r906/r1028/r1029/§6y5 together: protein edges are real, they
find families nucleotide misses, and **there is no available truth that can license them**. Reopening it
requires a curated ancient-paralogue truth (gene trees, Ensembl Compara paralogues, or a domain database to
subtract folds) — a new external dependency and a separate question, not a floor.

## Prediction scorecard

Predicted ⚠ VACUOUS: *"the false-merge rate should fall cleanly with the floor while the protein-only count
collapses toward zero."* **Wrong on the first half and right on the second.** The count did collapse
(1009 → 113 on chr16), but **the rate did not fall at all** — which is the informative failure, because it
is what exposed that the fold spans the protein. I had the right suspicion about redundancy and the wrong
model of what the 0.30 floor was admitting.

> **Generator (2026-09-22 consolidation):** `python3 bench/edge_probes.py protein-false-merge --floors 0.30,0.50,0.70,0.80,0.90 ...` — the original `bench/protein_false_merge.py` was folded in verbatim (§6z2); the register rows above cite this file.
