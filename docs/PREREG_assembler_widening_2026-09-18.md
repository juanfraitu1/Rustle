# Pre-registration — read-isoform widening in the assembler, single-exon strand, all-alignment pool

**Written 2026-09-18 BEFORE any number.** Follows §6m5, which found the reconstruction bottleneck is
EXACT-CHAIN COLLAPSE in `pass1_skeletons_robust_with` (groups keyed by exact intron chain, filtered
`n >= min_reads`), not ambiguity: dropped junctions have median largest-chain **2 reads** vs **32** kept.

## Baselines (measured, §6m5) — 9-copy co-located cluster, 106 annotated junctions

| | junctions | complete |
|---|---|---|
| assembler today | **64 (60.4%)** | 0/9 |
| primary-read union (the target) | 85 (80.2%) | 2/9 |
| ALL-alignment union (§6m4 ceiling) | 91 (85.8%) | 4/9 |

## Change 1 — read-isoform widening (PORT)

`shared_definition::widen_with_read_isoforms` admits a chain when EVERY junction has ≥ k read support,
instead of requiring one exact chain to clear a floor. Ported into `pass1_skeletons_robust_with`:
a group passes if `n >= min_reads` **OR** (`isoform_k > 0` AND every junction of its chain has
per-`(chrom, junction)` support ≥ `isoform_k` across all spliced reads in the region).
CLI `copy_assign --read-isoform-k`, **default 0 = OFF**.

- **W-1 (off ⇒ byte-identical).** `--read-isoform-k 0` reproduces the §6m5 GTF byte for byte. A failure
  here blocks everything else.
- **W-2 (the win).** Junction coverage on the cluster reaches **≥ 74/106** (i.e. ≥ +10) at some
  k ∈ {2, 3, 5, 8}. k is reported at every value; no k is chosen after seeing the result without saying so.
- **W-3 (no runaway).** Coverage must not exceed the ceiling **91/106**; exceeding it is a bug, not a win.
  Emitted-transcript count is reported at every k — a rule that triples the transcript count to buy
  junctions has bought fragmentation.

## Change 2 — single-exon isoforms (EXISTING FLAG, not new code)

Measured: **97 of 750 transcripts (12.9%) are single-exon and ALL 97 are `'+'`**, while the 653 spliced
split 327 `-` / 326 `+`. 93/97 sit inside an annotated NPIP locus and **91/97 overlap a spliced
transcript**. 15.3% of reads in the NPIP windows are unspliced. The `'+'` is the documented placeholder;
`RUSTLE_READ_STRAND=1` (with `RUSTLE_READ_STRAND_MARGIN`, default 0.90) replaces it with the FLAG-0x10
read-orientation majority, already validated at **386/400 = 0.965** against junction-determined strand
versus **0.4867** for the constant.

- **S-1.** With the flag on, the single-exon `'+'` fraction falls below **0.90**.
- **S-2.** Junction coverage must not DECREASE. (Single-exon models carry no junctions, so this is a
  no-harm check, not a recall bar. **Stated up front: this change cannot raise junction coverage** — it is
  a correctness fix, and register 843 already measured folding single-exon fragments as worth only +0.012.)
- **S-3.** Report how many single-exon models ABSTAIN under the margin rather than flipping.

## Change 3 — all-alignment pool (measured LAST)

- **A-1.** Coverage rises above whatever Change 1 reaches, toward the **91/106** ceiling, with the
  transcript count reported. Expected gain is small (**+6 junctions** is the entire headroom).

## Global

- The §6m4 ceiling is an integrity check throughout: **> 91/106 or > 4/9 complete = bug or truth leakage.**
- No annotation is an input to any arm; the annotation is used only for scoring.
- Every arm reports transcript count, junction coverage, complete chains and exact chains together.
