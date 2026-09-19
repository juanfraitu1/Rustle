# Assembly polish: matching StringTie in `--assemble-only` mode (§6p8, 2026-09-19)

Pre-registration: `docs/PREREG_assembly_polish_2026-09-19.md` (written before any chr11 number existed).
Shipped as `copy_assign --assembly-polish <none|mono|full>` (default `none` = byte-identical old output,
verified: the `none` GTF diffs clean against the published `ours_2026_09_19/base.gtf`).

## What the two filters are

Both use ONLY the `reads "N"` attribute the assembler already emits on every transcript line. No
reference, no annotation — legal in de novo mode.

1. **Mono-exonic support floor** (`mono`). A single-exon transcript carries no junction evidence at all,
   so it must reach the upper quartile of the multi-exon read support **in the same run**
   (`--polish-mono-quantile 0.75`). Self-tuning: chr20 and chr11 independently landed on a floor of 8.
2. **Support-aware ISM collapse** (`full` = 1 + 2). Drop a transcript whose intron chain is a contiguous
   sub-chain of another's on the same contig and strand, unless it carries at least as much read support
   as its container. Mono-exonic transcripts inside a multi-exon span are handled the same way.

Both passes are order-deterministic (containers sorted by chain length then transcript id); the Rust
implementation is byte-identical to the Python reference `bench/assembly_polish.py` on both chromosomes,
and run-to-run identical.

## Result

Human testis Iso-Seq on CHM13 (`human_testis.t2t.bam`), RefSeq reference, gffcompare v0.12.10.
`k0` = baseline flags, `k3` = `--read-isoform-k 3` + `RUSTLE_JUNCTION_MAJORITY=1`.

### chr20 — development substrate

| arm | mRNAs | intron-chain Sn / Pr | transcript Sn / Pr | **matching chains** | novel loci |
|---|---|---|---|---|---|
| k0 raw | 976 | 8.0 / 44.6 | 7.6 / 35.6 | 345 | 123/456 |
| k0 `mono` | 794 | 8.0 / 44.6 | 7.6 / 43.6 | 345 | 23/319 |
| ⭐**k0 `full`** | 682 | **7.9 / 50.8** | **7.4 / 49.6** | **337** | 23/319 |
| k3 raw | 1,275 | 8.4 / 33.4 | 7.9 / 28.2 | **358** | 123/458 |
| k3 `mono` | 1,130 | 8.4 / 33.4 | 7.9 / 31.8 | **358** | 48/350 |
| k3 `full` | 794 | 8.1 / 45.1 | 7.6 / 43.8 | 347 | 26/326 |
| StringTie 3.0.1 `-L -p 4` | 712 | 7.7 / 47.4 | 7.3 / 47.1 | 331 | 19/359 |
| FLAIR 3.0.0 | 820 | 6.2 / 35.1 | 5.8 / 32.3 | 264 | 61/335 |

⭐**On chr20 `k0 full` beats StringTie on all four gffcompare axes** — matching chains 337 vs 331,
intron-chain 7.9/50.8 vs 7.7/47.4, transcript 7.4/49.6 vs 7.3/47.1 — with 30 fewer emitted transcripts.

### chr11 — HELD OUT (chosen before any chr11 measurement; 10,534 reference transcripts)

| arm | mRNAs | intron-chain Sn / Pr | transcript Sn / Pr | **matching chains** | novel loci |
|---|---|---|---|---|---|
| k0 raw | 2,059 | 7.3 / 42.6 | 6.8 / 34.8 | 715 | 208/870 |
| k0 `mono` | 1,728 | 7.3 / 42.6 | 6.8 / 41.5 | 715 | 49/621 |
| k0 `full` | 1,465 | 7.0 / 48.3 | 6.6 / 47.2 | 690 | 46/616 |
| k3 raw | 2,720 | 7.8 / 32.5 | 7.2 / 28.1 | **761** | 208/876 |
| k3 `mono` | 2,437 | 7.8 / 32.5 | 7.2 / 31.3 | **761** | 77/671 |
| k3 `full` | 1,727 | 7.4 / 43.2 | 6.9 / 42.0 | 723 | 54/637 |
| StringTie 3.0.1 `-L -p 4` | 1,307 | 6.6 / **50.2** | 6.2 / **50.0** | 648 | 40/653 |

## Pre-registered hypotheses — verdicts

| | bar | chr11 result | verdict |
|---|---|---|---|
| H1 | filters cost ≤ 3% of matching chains | `full` 715 → 690 = **−3.50%** | ⛔ **FAIL** (`mono` alone: 0.00%, PASS) |
| H2 | transcript precision +≥ 8 points | 34.8 → 47.2 = **+12.4** | ✅ PASS (`mono` alone +6.7 would fail) |
| H3 | still beat StringTie on matching chains | 690 vs 648 (+6.5%) | ✅ PASS |
| H4 | transcript precision within 3 pts of StringTie | 47.2 vs 50.0 = **−2.8** | ✅ PASS |

⚠**H1 failed by half a point.** The ISM collapse genuinely costs real transcripts on chr11 (25 chains),
about three times its chr20 cost (8). The mono floor is the free half of the rule and the ISM collapse is
the paid half; they are reported separately for that reason and `mono` is the safer default of the two.

## Honest reading of the goal ("match or outperform the assembly tools")

- **Sensitivity: we outperform on both chromosomes, at every polish setting.** Matching intron chains
  337–358 vs 331 on chr20, 690–761 vs 648 on chr11. The recall lead is the robust result.
- **Precision: outperformed on chr20 (50.8/49.6 vs 47.4/47.1), NOT on chr11 (48.3/47.2 vs 50.2/50.0).**
  StringTie keeps a ~2-point precision lead on the held-out chromosome. "Match" is the right word there,
  "outperform" is not.
- The `mono` setting is the one free lunch: **zero matching chains and zero sensitivity lost on both
  chromosomes**, +8.0 (chr20) / +6.7 (chr11) transcript precision, and novel loci 123 → 23 and 208 → 49.
  Every one of our single-exon predictions on chr20 was junk against this reference — the filter removed
  182 of them and cost nothing.

⚠Both chromosomes are ordinary human autosomes measured against a single annotation; the
`docs/o1_ledger.md` §6kl/§6km ground-truth ceiling applies, and neither chromosome measures the
project's multi-copy contribution.


---

# §6p9 — the locus isoform fraction closes the gap (four chromosomes)

The chr11 precision deficit above was **entirely class `j`** ("novel junction combination"): 589 of them
against StringTie's 481, which is the whole 119-transcript non-matching excess (773 vs 654). Every other
gffcompare class code was at parity or better. A `j` transcript shares junctions with a reference
transcript but its chain does not match — a minor alternative flow at a locus we already reconstruct.

**`--polish-isoform-fraction F`** (new): drop a transcript whose `reads` is below `F ×` the best-supported
transcript at the same `gene_id`; the locus dominant is never dropped, so no locus is emptied. This is
StringTie's `-f` criterion, which the assembler had never applied. (`--min-isoform-fraction` is a
different thing: a fraction of the locus TOTAL, and it only tags `low_confidence`.)

F was chosen on chr20 alone by the Addendum-A rule — largest F with ≤1% chain loss — which returned
**F = 0.02**. Two further chromosomes were then built from scratch to test it: **chr7** and **chr14**,
neither previously touched by this project.

## Recommended setting: `--assemble-only --assembly-polish full --polish-isoform-fraction 0.02`

| | chr20 (dev) | chr11 | chr7 | chr14 |
|---|---|---|---|---|
| reference transcripts | 4,574 | 10,534 | 8,726 | 6,241 |
| **matching intron chains** | **335** / 331 | **683** / 648 | **515** / 515 | **389** / 388 |
| matching transcripts | **336** / 335 | **685** / 653 | **518** / 516 | 391 / **392** |
| intron-chain Sn | **7.8** / 7.7 | **7.0** / 6.6 | 6.4 / 6.4 | **7.1** / 7.1 |
| intron-chain Pr | **51.9** / 47.4 | **50.3** / 50.2 | **44.7** / 43.5 | **43.8** / 42.0 |
| transcript Sn | **7.4** / 7.3 | **6.5** / 6.2 | 5.9 / 5.9 | **6.3** / 6.3 |
| transcript Pr | **50.6** / 47.1 | 49.2 / **50.0** | **43.6** / 43.0 | **42.5** / 41.9 |
| emitted mRNAs | 664 / 712 | 1,393 / 1,307 | 1,188 / 1,199 | 921 / 935 |

(ours / StringTie 3.0.1 `-L -p 4`; bold = ours at least matches.)

⭐**Scorecard: 19 of 20 (chromosome × metric) cells match or outperform StringTie.** The single miss is
chr11 transcript precision, 49.2 vs 50.0.

⭐**On matching intron chains — "does it find real transcripts" — we match or beat StringTie on all four
chromosomes**, including both chromosomes built after the rule was fixed.

## The remaining chr11 miss, measured

It is mono-exonic transcripts, not chains. At this setting chr11 keeps 36 single-exon predictions to
StringTie's 16; they contribute 2 matches. Removing all of them gives chr11 transcript precision 50.3
(> 50.0, 5/5) — but costs chr14 two real matching transcripts and drops chr14 to 4/5. Single-exon
predictions are therefore mostly, but not always, junk, and no single-exon policy is 5/5 everywhere:

| policy | chr20 | chr11 | chr7 | chr14 |
|---|---|---|---|---|
| mono floor at p75 (shipped) | 5/5 | **4/5** | 5/5 | 5/5 |
| drop every single-exon transcript | 5/5 | 5/5 | 5/5 | **4/5** |

## Negative result: the (k, F) grid rule picked a worse cell

Addendum B registered a second selection — over `k ∈ {0,3} × F`, take the cell where all four rates and
the chain count beat StringTie on chr20 and maximise chains — which returned **k = 3, F = 0.05**. Held
out, that cell is *worse*: 5/5 on chr20 and chr14 but **3/5 on chr11 and 1/5 on chr7** (chr7: 514 chains
vs 515, chain Pr 42.7 vs 43.5, transcript Pr 40.9 vs 43.0). Its primary hypothesis B3 (chr14) passes and
B1/B2 fail on two chromosomes. Adding recall with `--read-isoform-k 3` and buying it back with a larger
F is a worse trade than not adding it: **k = 0 with F = 0.02 dominates the k = 3 arms on every held-out
chromosome.** Register row 855.

## Per-chromosome F sensitivity (post-hoc, NOT a validated selection)

| F (`full`, k=0) | chr20 chains/transPr | chr11 | chr7 | chr14 |
|---|---|---|---|---|
| 0.00 | 337 / 49.6 | 690 / 47.2 | 519 / 42.4 | 393 / 41.1 |
| **0.02** | **335 / 50.6** | **683 / 49.2** | **515 / 43.6** | **389 / 42.5** |
| 0.03 | 335 / 51.9 | 679 / 50.5 | 511 / 44.3 | 385 / 43.4 |
| 0.05 | 329 / 53.6 | 653 / 52.6 | 496 / 44.8 | 381 / 44.6 |
| 0.08 | 321 / 54.5 | 623 / 54.4 | 482 / 45.4 | 368 / 45.3 |

⚠ chr11 alone would prefer F = 0.03 (5/5) and chr7 alone F = 0.02; no F is 5/5 on all four. Do not quote
a per-chromosome best as if it were the rule — F = 0.02 is the one that was fixed in advance on chr20.

## Substrates built for this work

`bakeoff/human_chr{11,7,14}`, each: chromosome BAM sliced from `human_testis.t2t.bam` (**not** the deeper
`A119b.t2t.bam`, which the chr20 bakeoff does not use), chromosome FASTA from `chm13v2.0.fa`, and a
reference GTF from `chm13v2.0_RefSeq_full.gff.gz` via `/mnt/linuxdisk/tmp/gff2gtf.py` — validated to
reproduce gffread's `chr20_ref.gtf` transcript count exactly (4,574 = 4,574); `gffread` is not installed.
StringTie 3.0.1 `-L -p 4` on the same BAM in every case. FLAIR was run on chr20 only.
