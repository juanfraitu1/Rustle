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
| StringTie 3.0.1 `-L -p 4` | 712 | 7.7 / 47.4 | 7.3 / 47.1 | 331 |

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
