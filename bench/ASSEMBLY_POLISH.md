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
