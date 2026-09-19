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
