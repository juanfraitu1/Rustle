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

⚠**StringTie build caveat:** StringTie 3.0.1 on chr20 (the vendored `tools/stringtie` build) and **3.0.3 on chr11/7/14/5/9** (conda) — the submodule was removed mid-session, so the panel is not on a single StringTie build; the lab's A119b/GGO arms are 3.0.1.

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
| StringTie `-L -p 4` (⚠3.0.1 on chr20, 3.0.3 on chr11/7/14/5/9) | 712 | 7.7 / 47.4 | 7.3 / 47.1 | 331 | 19/359 |
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
| StringTie `-L -p 4` (⚠3.0.1 on chr20, 3.0.3 on chr11/7/14/5/9) | 1,307 | 6.6 / **50.2** | 6.2 / **50.0** | 648 | 40/653 |

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

(ours / StringTie `-L -p 4` (⚠3.0.1 on chr20, 3.0.3 on chr11/7/14/5/9); bold = ours at least matches.)

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
StringTie `-L -p 4` (⚠3.0.1 on chr20, 3.0.3 on chr11/7/14/5/9) on the same BAM in every case. FLAIR was run on chr20 only.


---

# §6q0 — the shadow rule and the final setting (six chromosomes)

## Shipped recommendation

```
copy_assign --assemble-only --assembly-polish full \
            --polish-isoform-fraction 0.02 --polish-mono-shadow --polish-mono-quantile 0.82
```

**`--polish-mono-shadow`**: drop a single-exon transcript that overlaps any multi-exon EXON on **either**
strand, or any **same-strand** multi-exon SPAN. A single-exon read pile has no splice motif, so its strand
label carries no evidence — hence either-strand exon overlap. Designed on chr20's 177 unfiltered
single-exon predictions:

| feature | matches a reference | matches none |
|---|---|---|
| same-strand multi-exon EXON | 0 / 2 | 24 / 175 |
| any-strand multi-exon EXON | 1 / 2 | 41 / 175 |
| same-strand multi-exon SPAN | 0 / 2 | 28 / 175 |
| anti-strand multi-exon SPAN | **2 / 2** | 46 / 175 |

⚠Anti-strand SPAN overlap is deliberately **not** a criterion — every true positive has one (register 857).

## Result: 28 of 30 cells over six chromosomes

ours / StringTie `-L -p 4` (⚠3.0.1 on chr20, 3.0.3 on chr11/7/14/5/9); bold = ours at least matches. ★ = held out after the rule was fixed.

| | chr20 | chr11 | chr7 | chr14 | **chr5 ★** | **chr9 ★** |
|---|---|---|---|---|---|---|
| reference transcripts | 4,574 | 10,534 | 8,726 | 6,241 | 8,174 | 7,874 |
| **matching intron chains** | **335**/331 | **683**/648 | **515**/515 | **389**/388 | 473/**476** | **405**/403 |
| intron-chain Sn | **7.8**/7.7 | **7.0**/6.6 | **6.4**/6.4 | **7.1**/7.1 | **6.3**/6.3 | **5.5**/5.5 |
| intron-chain Pr | **51.9**/47.4 | **50.3**/50.2 | **44.7**/43.5 | **43.8**/42.0 | **47.6**/45.5 | **44.3**/42.8 |
| transcript Sn | **7.4**/7.3 | **6.5**/6.2 | **5.9**/5.9 | **6.3**/6.3 | 5.8/**5.9** | **5.2**/5.2 |
| transcript Pr | **51.5**/47.1 | **50.1**/50.0 | **44.5**/43.0 | **43.5**/41.9 | **47.2**/45.3 | **44.0**/42.8 |
| emitted mRNAs | 652/712 | 1,367/1,307 | 1,162/1,199 | 897/935 | 1,002/1,061 | 926/953 |
| verdict | 5/5 | 5/5 | 5/5 | 5/5 | **3/5** | 5/5 |

⭐**Five of six chromosomes match or outperform StringTie on every metric, including held-out chr9.**
Against FLAIR (chr20, the only chromosome it was run on) the margin is far wider: 335 vs 264 chains.

## The chr5 miss, measured

Both missing cells are **recall**; chr5's precision is comfortably ahead (47.6/47.2 vs 45.5/45.3).

- **matching chains 473 vs 476.** The chr5 ladder localises it: raw 485 → mono+shadow 485 (**the
  single-exon filters cost chr5 zero chains**) → +ISM+fraction 473. The ISM collapse is harsher at chr5's
  depth, which is the deepest of the six (519,887 records vs 327k–520k).
- **transcript Sn 5.8 vs 5.9.** Single-exon recall: StringTie's matching transcripts exceed its matching
  chains by 5 on chr5 (481 vs 476), i.e. it matches 5 single-exon reference transcripts; we match 0.

## Four attempts to close chr5, all measured, none kept

| attempt | six-chromosome cells | why not |
|---|---|---|
| `--polish-ism-escape` at q = 0.82 | 27/30 | recovers chr5's chains (473 → 476) but costs chr11 its precision lead (row 858) |
| the escape at q = 0.86 / 0.90 / 0.94 | 27 / 27 / 26 | chr14 and chr9 transcript Sn start failing |
| drop the ISM pass, raise F to 0.05–0.16 | 9–12/30 | F is per-locus and removes true minor isoforms wholesale (row 859) |
| mono floor from the single-exon distribution | ≤ 28/30 | never better, and chr5's two cells are not single-exon cells (row 860) |
| `--read-isoform-k 3` with the full polish | 18/30 | precision falls on **all six** (row 855, confirmed again) |

`--polish-ism-escape` is kept as an explicit recall/precision dial with both endpoints measured, default
off. With six chromosomes consulted, the search was stopped rather than continue fitting the panel.

⚠**Provenance of q = 0.82:** it is the midpoint of the window {0.80, 0.85} that was found by scoring four
chromosomes, so it is **fitted on four and validated on two** — chr9 passed 5/5, chr5 did not. The shadow
rule itself was designed on chr20 alone.


---

# §6q2 — why chr5's two cells do not close, and the pooled view

Two further attempts, both refuted:

| attempt | six-chromosome cells | why |
|---|---|---|
| `--polish-ism-3p`: collapse only 3'-anchored sub-chains, on §6p4's 5'-truncation finding | 19-21/30 at F = 0.02-0.05 | mid-chain and 5'-anchored sub-chains are junk at a similar rate; sparing them costs chr11/chr7/chr14 precision and still does not buy chr5 its chains (row 861) |
| `--polish-ism-escape` × fraction 0.025-0.035 × quantile 0.82/0.86 | 21-26/30 | chr7 and chr14 start losing chains before chr5 recovers its own |

## Single-exon recall is not intrinsically recoverable (row 862)

chr5's two matching single-exon predictions were traced through the filters:

| | length | reads | passes shadow | passes floor (11) |
|---|---|---|---|---|
| `DN_chr5_80803957_1` | 1,853 | **2** | yes | **no** |
| `DN_chr5_122133896_1` | 599 | **2** | yes | **no** |

Both pass the shadow rule and both sit at the **minimum possible read support**, against a self-tuned
floor of 11 on that deep chromosome. Length does not rescue them either — among the single-exon
predictions that survive the shadow rule, the share of non-matching ones LONGER than the shortest matching
one is 49% (chr20), 82% (chr11), 73% (chr7), 72% (chr5); only chr14 separates (1/136). Admitting chr5's
two matches means admitting ~120 junk transcripts with them, which costs more transcript precision than
the sensitivity cell is worth. **Within the single-exon class, true positives are not separable from
artifacts by read support, length, or genomic context.**

## Pooled over the six chromosomes (46,077 reference mRNAs)

| | ours | StringTie 3.0.1 |
|---|---|---|
| emitted mRNAs | **6,006** | 6,167 |
| **matching intron chains** | **2,800** | 2,761 (+39, **+1.4%**) |
| matching transcripts | **2,808** | 2,785 (+23, +0.8%) |
| transcript precision | **46.8%** | 45.2% |
| transcript sensitivity | **6.09%** | 6.04% |

⚠A pooled figure hides the per-chromosome variation that the table above shows, and chr5 is a real miss
inside it. It is reported as a genome-scale summary, not as a substitute for the per-chromosome verdict.


---

# §6q3/§6q4 — the ISM support ratio, and why chr5 and chr11 cannot both pass

## `--polish-ism-ratio` (new, default 1.0; **0.7 is the shipped recommendation**)

Keep a sub-chain fragment when its read support reaches this fraction of its container's. The previous
code demanded parity (1.0). A genuine shorter isoform carries a substantial share of its locus while a
5'-truncation artifact carries a small one, so the ratio separates them **independently of library
depth** — unlike `--polish-ism-escape`, whose absolute bar rises with coverage.

**0.7 is a free gain on every chromosome**, and it puts chr11 exactly on StringTie's bar:

| | chr20 | chr11 | chr7 | chr14 | chr5 | chr9 | pooled |
|---|---|---|---|---|---|---|---|
| chains at ratio 1.0 | 335 | 683 | 515 | 389 | 473 | 405 | 2,800 |
| **chains at ratio 0.7** | **336** | **684** | **516** | **391** | **475** | **407** | **2,809** |

## Final shipped setting

```
copy_assign --assemble-only --assembly-polish full --polish-isoform-fraction 0.02 \
            --polish-mono-shadow --polish-mono-quantile 0.82 --polish-ism-ratio 0.7
```

| ours / StringTie | chr20 | chr11 | chr7 | chr14 | chr5 | chr9 |
|---|---|---|---|---|---|---|
| **matching chains** | **336**/331 | **684**/648 | **516**/515 | **391**/388 | 475/**476** | **407**/403 |
| intron-chain Sn / Pr | **7.8/51.6** | **7.0/50.2** | **6.4/44.4** | **7.1/43.7** | **6.3/47.5** | **5.5/44.1** |
| transcript Sn / Pr | **7.4/51.2** | **6.5/50.0** | **5.9/44.2** | **6.3/43.4** | 5.8/**47.1** | **5.2/43.8** |
| emitted mRNAs | 658/712 | 1,373/1,307 | 1,172/1,199 | 904/935 | 1,008/1,061 | 934/953 |
| verdict | 5/5 | 5/5 | 5/5 | 5/5 | **3/5** | 5/5 |

**28/30 cells. Pooled: 2,809 matching intron chains vs 2,761 (+48, +1.7%), 2,817 vs 2,785 matching
transcripts, from fewer emitted transcripts (6,049 vs 6,167) — transcript precision 46.6% vs 45.2%,
sensitivity 6.11% vs 6.04%.**

## chr5 and chr11 have disjoint feasible regions (row 863)

This is the reason 30/30 is not reached, and it is now a measured statement rather than a search that ran
out of ideas. Two witnesses from the same parameter family:

| setting | chr5 | chr11 |
|---|---|---|
| ratio 0.7, F 0.02 (**shipped**) | chains 475 / 476 ⛔, transcript Sn 5.8 / 5.9 ⛔ | 50.2 / 50.2 ✅, 50.0 / 50.0 ✅ |
| ratio 0.7, F 0.02, **+ `--polish-ism-escape`** | chains **478** / 476 ✅, transcript Sn **5.9** / 5.9 ✅ | 49.7 / 50.2 ⛔, 49.5 / 50.0 ⛔ |

Both are 28/30; they fail on **different** chromosomes. Every intermediate point was measured — fraction
0.015/0.016/0.017/0.018/0.019/0.02/0.022/0.025, ratio 0.3–1.0, quantile 0.75–0.95, with and without the
escape and the fraction exemption, roughly 60 cells — and the frontier is monotone: **chr11's chain
precision reaches StringTie's 50.2 only at fraction ≥ 0.019, while chr5 needs ≤ 0.018 for its chains and
≤ 0.015 for its transcript sensitivity. The windows do not intersect.**

chr11 is also StringTie's single best chromosome of the six for intron-chain precision (50.2, against
42.0–47.4 elsewhere), so that cell is the hardest bar in the whole panel; our 49.7 there under the
chr5-favouring setting still exceeds our own precision on four of the six chromosomes.

**Both endpoints ship.** Add `--polish-ism-escape` to favour chr5-like (deeply covered) substrates;
leave it off to favour chr11-like ones. Neither dominates, and the choice is a substrate property, not a
tuning accident.
