# Pre-registration — assembly-only precision polish (2026-09-19)

Written BEFORE running gffcompare on the held-out chromosome.

## Rule (fixed on chr20; not re-tuned afterwards)

Two post-assembly filters on the `--assemble-only` GTF, using only the emitted
`reads "N"` attribute — no reference, no annotation.

1. **Support-aware ISM collapse.** Drop transcript *y* when its intron chain is a contiguous
   sub-chain of transcript *x*'s on the same contig and strand, UNLESS `reads(y) >= reads(x)`
   (`--support-ratio 1.0`). A mono-exonic transcript contained in a multi-exon transcript's span
   is treated the same way.
2. **Mono-exonic support floor.** A surviving single-exon transcript has no junction evidence at
   all, so it must carry at least as much read support as the **upper quartile of the multi-exon
   transcripts in the same run**: drop it when `reads < p75(reads over multi-exon transcripts)`.
   The threshold is computed per run, from the prediction itself.

## Development substrate (chr20, human testis Iso-Seq on CHM13)

| arm | mRNAs | chain Sn/Pr | transcript Sn/Pr | chains |
|---|---|---|---|---|
| raw | 976 | 8.0 / 44.6 | 7.6 / 35.6 | 345 |
| +ISM (ratio 1.0) | 842 | 7.9 / 50.8 | 7.4 / 40.3 | 337 |
| +mono floor (p75 = 8) | 682 | 7.9 / 50.8 | 7.4 / 49.6 | 337 |
| StringTie `-L -p 4` (⚠3.0.1 on chr20, 3.0.3 on chr11/7/14/5/9) | 712 | 7.7 / 47.4 | 7.3 / 47.1 | 331 |

## Held-out substrate: chr11, same BAM (`human_testis.t2t.bam`), same reference source

chr11 was chosen before looking at any chr11 result: it is a different, gene-denser chromosome
(10,534 RefSeq transcripts vs chr20's 4,574) from the same library, and nothing in this project
has been tuned on it.

## Hypotheses (decided in advance)

- **H1 — the filters cost almost no true transcripts.** Matching intron chains on chr11 fall by
  **at most 3%** from raw to fully filtered. FAIL ⟹ the ISM/mono signature does not generalise.
- **H2 — transcript-level precision rises by at least 8 points** from raw to fully filtered
  (chr20 moved 35.6 → 49.6, +14.0). FAIL ⟹ the precision gain was a chr20 artifact.
- **H3 — we still beat StringTie on matching intron chains** on chr11 after filtering.
  FAIL ⟹ the recall lead does not hold on a second chromosome.
- **H4 — transcript-level precision reaches StringTie's, within 3 points or better.**
  This is the goal's actual bar; H2 can pass while H4 fails.

Denominator guard: the scored reference set is fixed (`chr11_ref.gtf`, 10,534 transcripts) and is
never filtered by the prediction. Sensitivity denominators are therefore identical across arms.

---

## VERDICTS (recorded after the chr11 run; nothing above was edited)

| chr11 held-out result | bar | outcome |
|---|---|---|
| `full`: matching chains 715 → 690 = **−3.50%** | H1 ≤ 3% | ⛔ **FAIL** by 0.5 pt |
| `mono` alone: 715 → 715 = **0.00%** | H1 ≤ 3% | ✅ PASS |
| transcript Pr 34.8 → 47.2 = **+12.4** | H2 ≥ +8 | ✅ PASS (`mono` alone +6.7 fails) |
| 690 matching chains vs StringTie 648 | H3 > StringTie | ✅ PASS (+6.5%) |
| transcript Pr 47.2 vs StringTie 50.0 | H4 within 3 pts | ✅ PASS (−2.8) |

**Consequence of the H1 failure:** the two filters are reported and shipped separately. `mono` is the free
half (0 chains lost on both chromosomes); `full` is a deliberate recall/precision trade, not a free win.
No threshold was re-tuned on chr11 — the floor it self-selected (8, then 9 under `full`) came from chr11's
own read distribution, which is the rule working as written, not a fit.

Full numbers: `bench/ASSEMBLY_POLISH.md`. Ledger: §6p8.

---

# Addendum A — locus isoform fraction (§6p9), written before the sweep

The held-out chr11 failure is entirely class `j` ("novel junction combination"): we emit 589 against
StringTie's 481, which is the whole of the 119-transcript non-matching excess (773 vs 654). Every other
class code is at parity or better. A `j` transcript shares junctions with a reference transcript but its
full chain does not match — a minor alternative flow at a locus we already reconstruct.

**New filter: `--polish-isoform-fraction F`.** Drop a transcript whose `reads` is below `F ×` the
best-supported transcript at the same `gene_id`. The locus dominant is never dropped, so no locus is
emptied. This is StringTie's `-f` criterion, which our assembler has never applied (the existing
`--min-isoform-fraction` is a fraction of the locus TOTAL and only tags `low_confidence`).

**Selection rule, fixed now:** F is swept on **chr20 only** and chosen as the LARGEST value at which
matching intron chains fall by at most **1%** from the unswept arm. No other criterion. Whatever F that
rule returns is then applied unchanged to the held-out chromosomes.

**Held-out chromosomes: chr11 and chr7.** chr11 has so far only been used to *diagnose* the class
composition above — no threshold has been fitted to it — but because it has now been inspected, a
completely untouched third chromosome (chr7) is added and is the primary held-out test.

**Hypotheses:**
- **A1** — on BOTH held-out chromosomes, matching intron chains fall by ≤ 2% from the unswept arm.
- **A2** — on BOTH held-out chromosomes, transcript-level precision meets or exceeds StringTie's.
- **A3** — on BOTH held-out chromosomes, matching intron chains still exceed StringTie's.

A2 is the goal's bar. A1 and A3 are the guards that stop A2 being bought with recall.

---

# Addendum B — the (k, F) operating point (§6p9), written before any chr14 number exists

## What Addendum A returned, and why it was not enough

Its rule ("largest F with ≤1% chain loss on chr20") returned **F = 0.02**. Verdicts:

| | chr11 | chr7 |
|---|---|---|
| A1 chain loss ≤ 2% | −1.01% ✅ | −0.77% ✅ |
| A2 transcript precision ≥ StringTie | 49.2 vs 50.0 ⛔ **FAIL** | 43.6 vs 43.0 ✅ |
| A3 matching chains > StringTie | 683 vs 648 ✅ | 515 vs 515 — a tie, **not** "exceed" ⛔ |

So F = 0.02 is not an operating point that matches or outperforms everywhere. The sweep also showed no
single F does: F = 0.03 dominates StringTie on all four axes on chr20 AND chr11 but loses 4 chains on
chr7, while F = 0.02 matches chr7 and misses chr11's precision by 0.8.

## The corrected reading

F trades chains for precision along one axis, so it cannot fix both ends alone. The other knob —
`--read-isoform-k 3` with `RUSTLE_JUNCTION_MAJORITY=1` — moves the OTHER way: it adds candidate chains at
a precision cost. Pairing them (generate more, then filter by locus share) is the actual two-sided lever.

## Selection rule, fixed now, on chr20 ONLY

Over the (k ∈ {0, 3}) × (F ∈ {0, .02, .03, .05, .08, .12}) grid already measured on chr20: take the cells
where **all four gffcompare rates are ≥ StringTie's AND matching intron chains ≥ StringTie's**, and among
them pick the one with the **most matching intron chains**. Ties break toward the smaller F.

Applied to the chr20 grid this returns **k = 3, F = 0.05** (338 chains; 7.9/49.9 and 7.4/48.3 against
StringTie's 331, 7.7/47.4, 7.3/47.1). No held-out chromosome was consulted.

## Held-out test

chr11 and chr7 have both now been inspected, so the primary held-out test is a **fourth, completely
untouched chromosome: chr14**. chr11 and chr7 are reported as secondary replication.

**Hypotheses (the goal's bar, stated for each of the three test chromosomes):**
- **B1** — matching intron chains ≥ StringTie's.
- **B2** — all four gffcompare rates ≥ StringTie's.
- **B3 (primary)** — B1 and B2 both hold on **chr14**, the untouched chromosome.

Anything short of B1 ∧ B2 on a chromosome is reported as a miss on that chromosome, with the number.

---

# VERDICTS for Addenda A and B (recorded after chr7 and chr14; nothing above was edited)

## Addendum A — F = 0.02, tested on chr11, chr7 and chr14

| | chr11 | chr7 | chr14 |
|---|---|---|---|
| A1 chain loss ≤ 2% from F = 0 | −1.01% ✅ | −0.77% ✅ | −1.02% ✅ |
| A2 transcript precision ≥ StringTie | 49.2 vs 50.0 ⛔ | 43.6 vs 43.0 ✅ | 42.5 vs 41.9 ✅ |
| A3 matching chains > StringTie | 683 vs 648 ✅ | 515 vs 515 — a tie, not "exceed" ⛔ | 389 vs 388 ✅ |

Over the four chromosomes and the five gffcompare quantities, F = 0.02 matches or outperforms StringTie
in **19 of 20 cells**; the miss is chr11 transcript precision. On matching intron chains alone it matches
or beats on all four. A3's chr7 tie is a "match", which the goal's wording admits but the hypothesis as
written did not.

## Addendum B — k = 3, F = 0.05

| | chr11 | chr7 | chr14 (primary) |
|---|---|---|---|
| B1 matching chains ≥ StringTie | 679 vs 648 ✅ | 514 vs 515 ⛔ | 403 vs 388 ✅ |
| B2 all four rates ≥ StringTie | 4/4 → transcript Pr 48.6 vs 50.0 ⛔ | 1/4 ⛔ | ✅ |
| B3 B1 ∧ B2 on chr14 | — | — | ✅ **PASS** |

B3, the primary hypothesis, passes. B1/B2 fail on chr11 and chr7, and the cell is strictly worse held out
than Addendum A's. **Conclusion: the selection rule in Addendum B was the wrong rule** — maximising chains
on the development chromosome subject to beating StringTie there selected a cell that overfits chr20's
recall/precision balance. The shipped recommendation is Addendum A's `k = 0, F = 0.02`. Register row 855.

Full numbers: `bench/ASSEMBLY_POLISH.md` §6p9. Ledger §6p9.

---

# Addendum C — the shadow rule and the final setting (§6q0), written before chr5/chr9 exist

## Why a new rule

At `k=0, F=0.02` the only miss was chr11 transcript precision (49.2 vs 50.0), and it was measured to be
**single-exon transcripts, not chains**: chr11 emitted 36 to StringTie's 16. Removing all of them fixed
chr11 but cost chr14 two real matching transcripts, so a discriminator was needed, not a blanket policy.

## The discriminator, designed on chr20 alone

On chr20's 177 unfiltered single-exon predictions (2 of which match a reference transcript):

| feature | matching | non-matching |
|---|---|---|
| overlaps a same-strand multi-exon EXON | 0 / 2 | 24 / 175 |
| overlaps an ANY-strand multi-exon EXON | 1 / 2 | 41 / 175 |
| overlaps a same-strand multi-exon SPAN | 0 / 2 | 28 / 175 |
| overlaps an anti-strand multi-exon SPAN | **2 / 2** | 46 / 175 |

**`--polish-mono-shadow`**: drop a single-exon transcript that overlaps any multi-exon EXON on either
strand, or any same-strand multi-exon SPAN. A single-exon read pile has no splice motif, so its strand
label carries no evidence — which is why exon overlap is taken on either strand. Anti-strand SPAN overlap
is deliberately NOT a criterion: both chr20 matches have one. Of the chr20 single-exon predictions that
survive the read floor and match a reference, the rule removes none.

## The mono floor quantile

With the shadow rule on, `--polish-mono-quantile` was swept over {0.75, 0.80, 0.85, 0.90, 0.95} on the
four chromosomes already in use. The 20-cell scorecard: 0.75 → 19/20 (chr11 transcript Pr 49.8), **0.80
→ 20/20**, **0.85 → 20/20**, 0.90 → 19/20 (chr14 transcript Sn), 0.95 → 19/20 (same).

⚠**That window was chosen by looking at all four chromosomes, so it is a fitted value, not a validated
one.** The registered setting is the **midpoint of the passing window, q = 0.82**, which maximises the
distance to both observed failure edges, and it is tested on two chromosomes that do not yet exist in this
project.

## Registered setting

```
copy_assign --assemble-only --assembly-polish full \
            --polish-isoform-fraction 0.02 --polish-mono-shadow --polish-mono-quantile 0.82
```

## Held-out test: chr5 and chr9, both untouched

**Hypotheses:**
- **C1 (primary)** — on **both** chr5 and chr9, all four gffcompare rates and the matching-intron-chain
  count are ≥ StringTie's, i.e. 5/5 on each, 10/10 over the two.
- **C2** — matching intron chains ≥ StringTie's on both.
- **C3** — the setting is still 20/20 on chr20/11/7/14 at q = 0.82 (a consistency check on the midpoint,
  not a held-out result).

A miss on either new chromosome is reported as a miss, with the number, and the fitted nature of q is
reported whatever the outcome.

---

# VERDICTS for Addendum C (recorded after chr5 and chr9; nothing above was edited)

Registered setting: `--assemble-only --assembly-polish full --polish-isoform-fraction 0.02
--polish-mono-shadow --polish-mono-quantile 0.82`.

| | chr5 ★ | chr9 ★ |
|---|---|---|
| C1 all four rates and the chain count ≥ StringTie | ⛔ **3/5** | ✅ **5/5** |
| C2 matching intron chains ≥ StringTie | 473 vs 476 ⛔ | 405 vs 403 ✅ |

**C1 fails on chr5, passes on chr9.** C3 holds: the midpoint q = 0.82 is 5/5 on all four development
chromosomes. **Overall 28 of 30 cells** across six chromosomes; five of six chromosomes are 5/5.

**The chr5 miss is recall, not precision** (chr5 precision is 47.6/47.2 against StringTie's 45.5/45.3).
Two cells: matching intron chains 473 vs 476 (−0.6%) and transcript Sn 5.8 vs 5.9. Its causes were
measured separately:
- the 3 chains are lost in the **ISM collapse**, which is harsher at chr5's depth (519,887 records, the
  deepest of the six); the mono filters cost chr5 **zero** chains (485 → 485).
- the transcript-Sn cell is **single-exon recall**: StringTie matches 5 single-exon reference transcripts
  on chr5 (matching transcripts 481 vs chains 476) and we match 0.

Four attempts to close chr5 were measured and all made the six-chromosome scorecard worse or equal:

| attempt | result |
|---|---|
| `--polish-ism-escape` (self-tuned absolute escape for the ISM pass) | recovers chr5's chains 473 → 476 but costs chr11 its precision lead → **27/30** (row 858) |
| the same escape at quantile 0.86 / 0.90 / 0.94 | 27, 27, 26 / 30 |
| drop the ISM pass entirely and raise F to 0.05–0.16 | **9–12 / 30** — F is a blunt per-locus instrument and removes true minor isoforms (row 859) |
| mono floor taken from the single-exon distribution instead of the multi-exon one (q 0.82/0.90/0.95) | 27, 27, **28** / 30 — never better, and chr5's two cells are untouched because they are not mono cells |
| `--read-isoform-k 3` with the full polish | **18/30**; precision falls on every one of the six (row 855 confirmed a second time) |

No further search was run: with six chromosomes consulted, continuing to tune would be fitting the panel
rather than the problem.
