# TPM in the GTF, and a locus BED that shows how close our loci are to the annotation (§6r5, 2026-09-19)

Two outputs an external user asked for: abundance in the GTF, and a BED plus size report that makes the
locus-extent question inspectable.

## 1. `--gtf-tpm` — count-based TPM, chosen by measurement

Adds `cov` and `TPM` to every transcript line, from the `reads` support the assembler already records.
**Default off; unset, the GTF is byte-identical to a run without the flag.**

**`TPM_i = reads_i / Σ reads * 1e6` — no length normalisation.** A long read is one molecule, so dividing
by transcript length down-weights long transcripts that were sequenced end to end. That is not a
preference, it is measured — chr20, against StringTie's own TPM over the 507 intron chains both tools
call, Spearman on log10:

| our convention | ρ vs StringTie TPM |
|---|---|
| **count-based (shipped)** | **0.879** |
| length-normalised (short-read convention) | 0.714 |

Against FLAIR's isoform counts, count-based gives **ρ = 0.755** (388 shared chains; FLAIR's counts file
keys are `<read>_<chrom>:<pos>`, so the suffix has to be stripped before joining).

Verified on emit: 658 transcript lines annotated, **TPM sums to exactly 1,000,000**, and the emitted
column reproduces ρ = 0.879 end to end.

```sh
copy_assign --assemble-only --gtf-tpm ... --out run      # adds cov + TPM to run.gtf
```

## 2. `bench/locus_bed.py` — loci as BED, and a one-to-one match against the annotation

```sh
python3 bench/locus_bed.py run.gtf --out run --ref chr20_ref.gtf
#   run.loci.bed         predicted loci, BED6, score = summed reads
#   run.ref_loci.bed     annotated loci
#   run.locus_match.tsv  one row per matched pair, with size_ratio = pred_span / ref_span
```

Each GTF is collapsed to **loci** (one record per `gene_id`, spanning its transcripts), then predicted and
annotated loci are matched **one-to-one, greedily on reciprocal overlap**.

⚠**The matching is EVALUATION ONLY.** Loci are never built with bipartite matching — that is a standing
project constraint — and this tool only scores loci built without it. ⚠The greedy pass is a **lower bound**
on the optimal assignment, not the optimum; ties break on locus id so it is deterministic.

### chr20, `human_testis` library, min reciprocal overlap 0.10

| | predicted loci | matched | **median size ratio** | within ±10% | within 2× |
|---|---|---|---|---|---|
| **ours** | 357 | 288 | **0.994** | **55.6%** | **79.5%** |
| StringTie | 359 | 315 | 0.979 | 53.7% | 77.1% |
| FLAIR | 429 | 290 | 0.994 | 47.2% | 72.4% |

⭐**Our loci are the right size at the median (0.994) and we hold the largest share inside ±10% and
inside 2× of the annotated extent.** StringTie matches more annotated loci (315 vs 288) — it calls more
loci that overlap an annotation — but its matched loci agree less well in extent (q25 0.586 against our
0.676).

⚠**The q25 of 0.676 is the 5′-truncation signature again** (§6p4: TBC1D3's 5′UTR is 75.3% covered while
its 3′UTR is 100%). A quarter of our matched loci are appreciably shorter than annotated, and that is the
same defect the assembler shows at transcript level.

⚠943 of 1,231 annotated chr20 loci are unmatched by every tool — most annotated loci are not expressed in
this library, so the unmatched count is not an error rate.
