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

## 2026-09-23 addendum — two precision levers became the `--assemble-only` defaults (§6za, register 1074-1076)

`--assembly-junctions strict` (the transcript product no longer inherits the §6m8 family-recovery tolerance
for non-canonical junctions) and `--polish-retained-ratio 10` (drop a chain whose exon contains another
transcript's junction carrying ≥ 10× its reads — the aligner's short-exon read-through, r1070). Held-out
gorilla, 26 contigs: intron-chain precision 33.1 → 35.6 for −0.46% matching chains; human chr20-22 16.0 →
18.7 for −1.6%. Pre-registration and tables: `docs/PREREG_assembly_precision_levers_2026-09-23.md`. The
2026-09-22 output is `--assembly-junctions majority --polish-retained-ratio 0`, byte-for-byte.

## 2026-09-27 addendum — `--polish-subchain off|tag|drop`: sub-chain tag, opt-in drop (register 1128-1130)

**What it flags.** A multi-exon transcript y is flagged when it is **an end-compatible contiguous sub-chain of a
longer emitted transcript of the same locus with >= 1/2 its reads**. That means some x satisfies all of:
- the same `gene_id`, contig and strand;
- a strictly longer intron chain;
- y's chain is the exact contiguous block `x.chain[k..k+m]`;
- y's first exon starts no more than 10 bp before x's exon k, and y's last exon ends no more than 10 bp after x's
  exon k+m;
- `2*reads(x) >= reads(y)`.

Among the qualifying x, the one with the most reads is named, with ties going to the first `transcript_id`.
- Mono-exon transcripts and the longest chain of a locus are never flagged.
- Every flag is decided on the input set.
- It uses reads only, so it is legal de novo. It runs on the output of `--assembly-polish`.

**How it differs from the shipped ISM collapse.** `--polish-ism-ratio 0.7` already removes contiguous sub-chains that
have < 0.7 of their container's reads, whatever their ends. This pass adds the **end condition**, which is what
separates the two cases on dev:
- a 5′-truncated copy ends inside the container's exons (co-emitted `c`: median overhang 0 bp);
- a real shorter isoform starts in the container's intron (`=` transcripts with a container: median 5′ overhang
  33-80 bp).

**Modes.** The prototype is `compat_collapse.py` e14c5646, and this is its port with the ½-reads condition kept and no
further guard.

| mode | effect |
|---|---|
| `off` (default) | byte-identical to a run without the flag (394/394 dev products against f480847a) |
| `tag` | output only. The flagged transcript line ends with `subchain_of "<x>"; subchain_missing "5p"\|"3p"\|"both";`, where `subchain_missing` names the ends of x's chain that y lacks, in transcript orientation (a `-` strand swaps them). Every other byte is unchanged, and no consumer (families, gffcompare, the Python loaders) parses the tag. **It does not mean "incomplete":** on dev, 63% of tagged transcripts are `c`, 30% `j`, and 8 of 864 are `=`. |
| `drop` | removes every flagged transcript, before `--gtf-tpm`, so TPM stays normalised. No locus and no mono-exon transcript is lost. **A deliberate trade; see the table.** |

When not `off`, the pass writes one log line (`SUB-CHAIN (tag|drop): N of M multi-exon transcripts ...`) and 8
`params.tsv` rows (`polish_subchain*`).

**How to enable.**
```
RUSTLE_POLISH_SUBCHAIN=tag  tools/rustle_pipeline.sh assemble --bam B --fasta G --out PREFIX   # or =drop
copy_assign --assemble-only ... --polish-subchain tag|drop                                      # direct
```
When `RUSTLE_POLISH_SUBCHAIN` is unset, the driver command is unchanged. The driver accepts only `off`, `tag` or `drop`
and exits with status 2 on any other value (`tools/rustle_pipeline.sh`). Exporting the variable also costs a needless catalog-cache miss; it never gives a stale
hit.

**Drop trade-off: DEV, IN-SAMPLE.** The rule was selected on these three contigs, and no held-out test was run.
Species are never pooled. The data are human A119b and gorilla OR6737, dedup, gffcompare 0.12.10 against the
contig-restricted RefSeq.

| contig | `c` share of multi-exon queries | partial share c+k+m+n (rel.) | `=` precision | `=` queries lost | genes losing their only `=`/`c` query: drop vs matched random drop (binding seed; 5-seed mean) | `j` dropped | loci / mono-exon lost |
|---|---|---|---|---|---|---|---|
| human chr20 | 5.7 → 2.8% | −18.2% | +1.02 pt | 2 | 7 (1.46%) vs 3; 3.2 | 74 | 0 / 0 |
| human chr16 | 6.3 → 3.1% | −18.2% | +1.07 pt | 3 | 9 (1.15%) vs 7; 6.6 | 171 | 0 / 0 |
| gorilla NC_073244.2 | 5.4 → 3.1% | −16.9% | +1.23 pt | 3 | 7 (0.79%) vs 5; 3.4 | 21 | 0 / 0 |

- **Precision.** A matched random drop from the same sub-chain pool (NULL_S) gets 78-91% of the precision gain, so
  most of it is the smaller denominator.
- **Genes.** The genes that lose their only `=`/`c` query keep every junction. Their surviving container is `j`, `k`
  or `m` for the same reference.
- **What `drop` also changes (not measured).**
  - It changes the families-stage representative of 15-19 loci per contig.
  - It shrinks 8-11 gene spans per contig by ≤ 10 bp, and `flag` reads those spans.

**Status: opt-in only, not a default candidate.**
- The held-out pre-registration (`docs/PREREG_complete_transcripts_2026-09-27.md`) was **withdrawn on dev**. The drop
  loses more intron-correct gene labels than the matched random drop on 14 of 15 seed × contig cells (clause C5).
- 14 reads-only guards failed to fix that (r1129). An annotation oracle shows that the limit is the reads'
  selectivity, not the clauses.
- The SIRV E0 truth is unspent. A default flip needs its own held-out prereg first.

## 2026-09-27 addendum 2 — `--polish-tss` and `--polish-tes`: read-proven 5′ and 3′ ends (opt-in; register 1131-1135)

Two opt-in passes on the output of `--assembly-polish` (the `--gtf` / `--assemble-only` path). Both default to `off`,
which is byte-identical to a run without the flags. **Neither is a default candidate.** Every number below is **DEV,
IN-SAMPLE**: the rules and every review fix were designed on human A119b chr20 and chr16 and gorilla OR6737
NC_073244.2. No held-out contig and no SIRV read was used. Species are shown apart and never pooled.

### Shared evidence

- **Records.** Every primary, spliced record (not secondary, supplementary, unmapped or QC-fail), deduplicated on
  (strand, start, end, intron chain). `--keep-coordinate-duplicates` is honoured, and the first record of a key wins.
- **Per record.** The oriented 5′ end, the oriented 3′ end, the oriented intron chain, and the cap flag. A record is
  **capped** when its 5′-most CIGAR op in read orientation is a 1-3 bp soft clip of G only. That is the non-templated G
  that template-switching reverse transcription adds opposite the m7G cap.
- **When.** The streaming reader collects them in pass 1, before the pool's home, AS-tie, window and dedup rules. Under
  `--materialize-reads` there is one lazy pass per window. Nothing is collected unless one of the flags is on. The
  memory cost is about 180 B per record, a few GB genome-wide.

### `--polish-tss off|tag|rescue|split` (driver `RUSTLE_POLISH_TSS`)

**The proof uses reads only.** The constants are α = 0.05 and TOL = 10 bp (the `--polish-subchain` tolerance), so the
window is W = 21 bp.
- **Null.** A constant 5′-truncation hazard h. The window count is negative binomial with mean h·ΣN(t) and variance
  μ + aμ², where N(t) is the number of reads crossing t on their way to the junction. h and a are fitted per contig on
  the internal exons of the polished multi-exon transcripts.
- **Acceptor stratum.** Within TOL of an acceptor through which reads enter the exon, the hazard is h·r(d), where r(d)
  is the observed start excess at offset d from internal-exon acceptors. Truncated molecules whose upstream-exon bases
  were soft-clipped pile at d = −2: r(−2) = 97.6 (chr20), 83.1 (chr16), 114.5 (gorilla). Before this stratum was added,
  66 / 63 / 33 % (chr20 / chr16 / gorilla) of the first design's rescues were exactly such piles.
- **Case A** (a short form y that is the 3′-flush intron chain of a longer x of the same locus):
  - proven iff the best W-window of y's first-junction 5′ ends has p·⌈L/W⌉·|candidates| < α;
  - in scope iff that window lies downstream of x's upstream exon.
- **Case B** (one chain, several TSSs): significant windows over the chain's own first exon, merged, and kept apart
  only across a valley.

**The cap-signal dependence.** Whether the rule may act is decided per contig, without labels. A contig has a cap
signal iff capped reads are a **majority** at its proven first-exon clusters **and** a **minority** among the starts
inside internal-exon bodies. Each half is a one-sided exact binomial test against ½ at α.

| contig | capped share of all records | capped at the proven clusters | capped in internal-exon bodies | cap signal |
|---|---|---|---|---|
| human chr20 | .391 | 30,897 / 43,997 (70%) | 1,656 / 42,192 (4%) | yes |
| human chr16 | .364 | 46,971 / 67,636 (69%) | 4,383 / 80,364 (5%) | yes |
| gorilla NC_073244.2 | .034 | 1,726 / 41,932 (4%) | 121 / 11,762 (1%) | **no** |

- **With a signal,** a window or cluster counts only when most of its reads are capped.
- **Without one,** `rescue` and `split` do nothing on that contig. The output is exactly `tag`'s. The log line says
  `NO CAP SIGNAL … (output = tag; withheld: N rescues, M splits)`, and `polish_tss_applied` records it.
- **Why.** The negative-binomial count excess alone also proves truncation piles (tss_critique §2). On gorilla, before
  this gate, the rule rescued 38 forms for +1 chain against a random-rescue null maximum of 1, which is not selective.
- **So the TSS rule's usefulness is a property of the library, not of the method or the species.** Why the OR6737
  library lacks the G clip is not established. Check the signal on any new library before claiming anything for it.

**Modes.**

| mode | effect |
|---|---|
| `off` (default) | byte-identical; no evidence collected |
| `tag` | output only. Each multi-exon transcript with ≥ 1 proven cluster ends with `tss_clusters "<pos>:<capped>/<n>,…";` (5′→3′; pos = the cluster's mode, 1-based genomic). Families are byte-identical |
| `rescue` | `tag`, plus case A. A proven, in-scope short form that the polish dropped at the ISM, isoform-fraction or retained-intron step is kept, with `tss_rescued "ism\|frac\|ret";`. Every other polish decision is unchanged, so every `off` transcript is present unchanged (monotone). Its 5′ end is the most-started position among **its own** exact-chain reads inside the proven window (ties 5′-most). A form with no own start inside the proven window is neither rescued nor protected: the window proves another chain's TSS (G1 fix, 2026-09-27) |
| `split` | `rescue`, plus case B. A chain with ≥ 2 cap-proven clusters is emitted as `<tid>` plus `<tid>_tss<i>`, with `tss_split "<i>/<k>";`. Each piece's 5′ end is its cluster's mode. Reads are divided among the pieces, so their sum and the TPM normalisation are unchanged. **Not monotone:** the chain itself moves its 5′ end and changes `reads` |

`--polish-subchain` is decided on the **final** set, after the rescue, the split and any TES end move. Under `drop`, the
proven in-scope forms are spared, but only where rescue or split acted.

### `--polish-tes off|tag|pas-end` (driver `RUSTLE_POLISH_TES`)

**The proof uses reads plus the genome.** It has no p-value and no fitted constant.
- For every multi-exon transcript of the output (after `--polish-tss`), take its own exact-chain, same-strand,
  deduplicated primary 3′ ends over its last exon.
- Cluster them by single linkage: a gap > 21 bp starts a new cluster. A cluster needs ≥ 2 reads, the pass-1 floor. Its
  mode is the most-ended position, with ties going 3′-most.
- A cluster is **PAS-proven** when a canonical AATAAA or ATTAAA lies wholly inside [mode−35, mode−10] (transcript
  orientation) **and** the mode is not internally primed. Primed means ≥ 60% A in the 20 bp downstream, or an A6 run
  there. That is r1064's instrument.
- The reads carry no poly(A) tail (they are FLNC-trimmed; the 3′ clip is 0 in 78% of reads), so the priming test has to
  be genomic.

**Why the genome and not a scan.**
- The 3′ background is extremely clumped. Its NB dispersion is a = 114 (chr20), 254 (chr16) and 368 (gorilla), against
  1.2 / 1.4 / 6.6 on the 5′ side.
- Without labels, "canonical PAS and not primed" holds at:
  - .585 / .543 / .650 of emitted 3′ ends (chr20 / chr16 / gorilla);
  - .022 / .074 / .125 of internal-exon 3′ piles;
  - .007 / .007 / .004 of random internal-exon positions.
- The likelihood ratio of a real end over a pile is therefore **27 / 7.3 / 5.2**.
- **No library property is needed, so the TES proof runs on both species, but it separates about 5× less well on
  gorilla.**

**Modes.**

| mode | effect |
|---|---|
| `off` (default) | byte-identical; no evidence collected for it |
| `tag` | output only. Each multi-exon transcript ends with `tes_clusters "<n PAS-proven>"; tes_pas "yes\|no"; tes_primed "yes\|no";`. The last two describe the **emitted** 3′ end, which pass 1 takes from the single most-3′ read (r1135). Families are byte-identical |
| `pas-end` | `tag`, plus an end move. When the emitted end is primed **and** the transcript has ≥ 1 PAS-proven cluster, the 3′ end moves to the most-3′ proven mode, with `tes_end_moved_from "<old 1-based genomic>";`. The move always goes upstream and stays within the transcript's own reads. There is no rescue and no split: the transcript set and `reads` are unchanged, and `cov` follows the new length |

A TES rescue of 3′-shorter forms and a tandem-APA split were designed and measured, but **not built** (r1132).

### How to enable

```
RUSTLE_POLISH_TSS=tag|rescue|split  RUSTLE_POLISH_TES=tag|pas-end \
  tools/rustle_pipeline.sh assemble --bam B --fasta G --out PREFIX
copy_assign --assemble-only ... --polish-tss tag|rescue|split --polish-tes tag|pas-end      # direct
```

- **Driver.** An unset variable leaves the command unchanged, `off` passes an explicit `off`, and any other value exits
  2.
- **Refusals.** `copy_assign` refuses either flag without `--gtf` or with `--assembly-polish none`.
- **Output when on.** `params.tsv` gets `polish_tss*` / `polish_tes*` rows, and each contig gets one log line
  (`[copy_assign] TSS PROOF (<mode>) <contig>: …`, `… TES PROOF …`).
- **Cache.** Exporting either variable costs a needless catalog-cache miss, never a stale hit.

### Dev effect: DEV, IN-SAMPLE

The metric is multi-exon queries against the contig-restricted RefSeq (Gnomon for gorilla), scored with gffcompare
0.12.10, on the seeded driver cell. Annotation is used as a label only.

**TSS `rescue`.** The null is a depth-matched random rescue from the same dropped pool, 20 seeds.

| contig | rescued | Δ chains (null mean, max) | annotated TSS ≤ 250 bp (null mean) | Δ query precision | Δ `c` | Δ ISM-like |
|---|---|---|---|---|---|---|
| human chr20 | 59 → **47** after the G1 fix | **+3** (0.4, 2) | 13 of 59 (1.3) | −0.18 pt | +4 | +37 |
| human chr16 | 116 | **+7** (1.0, 3) | 23 of 116 (3.2) | −0.16 pt | +3 | +67 |
| gorilla NC_073244.2 | 0 (38 withheld: no cap signal) | 0 | – | 0 | 0 | 0 |

- **It is selective on human chains, but it is a trade, not an improvement.** The rescued forms match a reference chain
  at .05 / .06, below the output's own precision. That is the r1064 / r1072 regime.
- **Provenance of the numbers.** chr20's 13 is after review fix F4, which moved 26 of the 59 5′ ends. Before F4 it was
  16. The chr16 row and the chain, precision, `c` and ISM-like columns come from tss_impl's binary (fbe4fa86, before
  the review fixes). The rescued set and the intron chains do not depend on the 5′ end, but `c` was not re-scored after
  F4, and chr16's on-modes were not rerun.

**TSS `split`.**

| contig | chains split (records added) | pieces with an annotated TSS ≤ 250 bp | chains whose pieces sit at ≥ 2 distinct annotated TSSs | families (O1) |
|---|---|---|---|---|
| human chr20 | 49 (+54) | .68 | 4 | 26 → 27, copies 97 → 99 |
| human chr16 | 88 (+94), pre-fix binary | .64 | 2 | 2 representatives change; one becomes a `_tss2` twin |
| gorilla NC_073244.2 | 0 (64 withheld) | – | – | unchanged |

- **Split's query-precision "gain" (+0.35 / +0.39 pt on human) is a metric trap.** A twin of an `=` chain counts as a
  second `=` query; distinct-chain precision equals rescue's.
- **RefSeq rarely lists same-chain TSS variants,** so the "≥ 2 distinct TSSs" count is a lower bound, not a refutation.
- **Split twins are not exempt from the families-stage representative choice.**

**TES `pas-end`** (chr20 and gorilla; chr16 was not run). The label is an annotated TES with the same last intron,
within 50 bp.

| contig | multi-exon | emitted end primed | moved (median shift upstream) | annotated TES ≤ 50 bp: before → after | genes recovered / lost | `=`, chains, precision | families (O1) |
|---|---|---|---|---|---|---|---|
| human chr20 | 5,048 | 1,027 | **106** (868 bp) | 6 → **32** (of 73 labelled) | 28 / 2 (of 51) | unchanged (1,059 `=`, 0 membership changes) | 26 → 27, copies 97 → 99 |
| gorilla NC_073244.2 | 3,899 | 455 | **118** (728 bp) | 8 → **66** (of 98 labelled) | 43 / 5 (of 67) | unchanged (1,595 `=`, 0 membership changes) | 29 → 29, copies 93 → 95 (a 2-copy and a 3-copy family merge; a new 2-copy family) |

- **The rest of the primed ends stay.** 921 (chr20) and 455 − 118 = 337 (gorilla) primed ends have no PAS-proven
  cluster.
- **The PAS is what selects.** With the same trigger and a different target, the new-near rate of the labelled moved
  ends is:

  | target | chr20 | gorilla |
  |---|---|---|
  | PAS (shipped) | .438 | .673 |
  | any unprimed cluster | .317 | .574 |
  | any upstream cluster | .236 | .473 |

  The unprimed target recovers slightly more in absolute count (37 vs 32; 77 vs 65), but it also loses more (9 vs 6;
  13 vs 7).
- **Classes.** The only gffcompare class changes are 2 `n→j` (chr20) and 1 `n→c` (gorilla).
- **Sub-chain drop.** Under `--polish-subchain drop`, 2 shortened transcripts per contig become end-compatible
  sub-chains and are dropped.
- **Gorilla's label gain is against Gnomon models,** which may be circular. Its label-free PAS separation is the weaker
  one (5.2×).

### Why both stay opt-in

1. **In-sample.** Both rules, and every critique and review fix, were designed on these three contigs. A default flip
   needs a held-out prereg on a fresh development contig outside V1-V6 (tss_critique §7, r1129). SIRV E0 is unspent,
   and a PAS-gated rule is inert on SIRV by construction (r1132).
2. **The TSS rescue trades precision for recall.** Its forms are annotated below the output's precision.
3. **The TSS rule depends on a library property.** The only gorilla library at hand has no cap signal, so on gorilla
   `--polish-tss` is only `tag`.
4. **`split` and `pas-end` change O1** (families, copies and representatives), and no truth on the dev contigs can
   score the change.
5. **The labels are lower bounds.** RefSeq and Gnomon are blind to same-chain TSS variants and to tandem APA.
6. **Open review items** (tss2_review):
   - **G1, medium — FIXED 2026-09-27.** A proven window now counts for a form only when the form has an own
     (exact-chain) start inside it; otherwise the form is neither rescued nor protected (`params.tsv` row
     `polish_tss_proven_no_own_start`). chr20 rescue 59 → 47 (the 12 fallbacks removed; `drop`+`rescue` 5,450 → 5,437
     transcripts, `DN_chr20_61676274_23.4` no longer spared); gorilla 3 proven forms lack an own start. `off` stays
     byte-identical (chr20 gtf/quant/families/assignments/params); `cargo test --release` 933 passed, 0 failed.
     Binaries `/mnt/linuxdisk/tmp/rustle_figures/tss3_bin_frozen/` (`copy_assign` b1709a96). The other table
     columns of this section were not re-scored after the fix.
   - **G2, low.** `pas-end` ignores read support. 7 of 106 chr20 moves go from a larger cluster to a smaller one, the
     largest from 24 reads to 4.
   - **G3, low.** The 21-bp linkage boundary is not unit-tested; the Python parity covers it.
   - **G4, low.** A contig missing from `--fasta` is silent: every PAS and priming test is false, and nothing moves.

**Verification.**
- **Off.** Byte-identical to the frozen pre-change products:
  - 818/818 products over 118 runs with the flags omitted;
  - 150/150 over 24 runs with an explicit `off`.
  - In substance these are the GTF, `params.tsv` and readthrough products; in `--assemble-only` the quant and families
    files are header-only.
- **On.** The GTFs are byte-identical to an independent Python reference in 38/38 chr20 and gorilla cells.
  Independent pysam re-derivations reproduced:
  - the rescue sets (59/59, 38/38 before F1);
  - chr20's `pas-end`: 106/106 moves, and 0 attribute mismatches over 5,631 transcripts.
- **Paths.** Streaming equals `--materialize-reads`. `--genome-wide` equals the concatenated per-contig runs for TSS;
  this was not rerun for TES.
- **Tests.** `cargo test --release`: 933 passed, 0 failed, 13 ignored.
- **Frozen** at `/mnt/linuxdisk/tmp/rustle_figures/tss2_bin_frozen/` (`copy_assign` aebfcf96).

## 2026-09-29 addendum 3 — `--bridge-regroup off|f1|f1v2`: bridge-aware regrouping (opt-in; register 1145-1147)

One opt-in pass over the final GTF of `--assemble-only`: the Rust port (`vg_family::bridge_regroup`) of two frozen
post-processors, **F1** = `bench/f1_bridge.py --mode full` (37ee8e77, `docs/PREREG_f1_bridge_locus_2026-09-28.md`)
and **F1v2** = F1 plus `f1v2.py --rule min` (b4e788ad, `docs/PREREG_f1v2_readshare_2026-09-29.md`). `off`, the
default, is byte-identical to a run without the flag. The rule, verbatim in the module header:

- **Bridge.** A transcript whose intron J is the only link between the other transcripts of its `gene_id` upstream and
  downstream of J (no component straddles J). The upstream side must have a PAS-proven 3′ cluster of deduplicated
  spliced primaries ending inside J's intron: `--polish-tes`'s cluster rule and constants, called, not copied. The
  downstream side must have V1 ≥ 1: reads starting inside J's intron in a real start cluster, reaching beyond the
  acceptor, with their own first exon inside the intron. This is the readthrough filter's `rt_v1`, now shared.
- **F1v2** keeps a bridge only when its transcripts carry fewer reads than EACH side: share = reads(link) /
  (reads(link) + min(UP, DOWN)) < 1/2. A tie abstains.
- **Output.** Bridges become `<gene_id>.fus<k>` with `fusion_of` / `fusion_junction`. The other transcripts split into
  exon-overlap pieces named as `--gtf-regroup` names them, so without a bridge the output is RG3's. No line is added,
  removed or moved, and intron chains are unchanged. `PREFIX.families.gtf` is the GTF without the bridges: it is the
  families input, as it was in the held-out runs.
- **Evidence.** The script's reads, collected by the pass-1 reader: every primary record with an `N`, QC-fail
  included, strand = `ts` XOR reverse. This is not `--polish-tes`'s evidence, which takes the alignment's strand and
  skips QC-fail; 1,230 of testis's 1.10 M spliced primaries carry `ts:A:-`.

### How to enable

```
RUSTLE_BRIDGE_REGROUP=f1v2 tools/rustle_pipeline.sh all --bam B --fasta G --out PREFIX   # assemble + families
copy_assign --assemble-only --genome-wide ... --bridge-regroup f1|f1v2                     # direct
```

- **Products when on.** `PREFIX.families.gtf`; `PREFIX.bridge_junctions.tsv` (F1's `junctions.tsv`: every structural
  junction with U, its 3′ clusters, V1 and the decision); under `f1v2` also `PREFIX.bridges.tsv` (`f1v2.py`'s table: the
  reads of each F1 bridge, its share and `keep`). `params.tsv` gets `bridge_regroup*` rows (the scripts' `stats.json`
  counts), and the log gets one `BRIDGE REGROUP` line.
- **Order.** It runs last: after the polish, `--polish-tss` / `--polish-tes`, the sub-chain drop and the attribute
  passes, as the scripts post-processed the emitted GTF. It includes `--gtf-regroup`'s split, so the two are exclusive:
  `copy_assign` and the driver refuse both. The held-out runs used it without `--gtf-regroup`.
- **Refusals.** `copy_assign` refuses it without `--assemble-only`, with `--families`, or with two regions of one
  contig (its V1 counts each record once). The driver refuses `families` on a bridge-regrouped GTF when the variable is
  unset.
- **Row order.** The side tables follow a single call of the script: contigs by name. The held-out tables merged
  contig batches, so they hold the same rows with the contig blocks in batch order.

### Read before quoting

- F1 fails on human A119b (bridge splits 219 SEP / 451 FRAG), and cuts about 40 annotated genes per gorilla sample.
- F1v2 was EFFECTIVE on both human libraries, but still leaves 149 single-gene cuts on A119b.
- The fused-locus gain comes from moving fusions into explicit `fusion_of` relation records, not from removing them.
  Counting each bridge as its own locus, F1v2 is no better than `--gtf-regroup` (A119b 1,829 vs 1,786).
- Opt-in; a default flip is the user's call.

### Verification (2026-09-29)

Inputs: the stored driver-default BAM runs of 2026-09-25 (`rustle_figures/runs/<s>/`), re-assembled by the port with
the driver's flags. `off` reproduced every stored BASE GTF byte for byte. The frozen outputs are those of the held-out
tests (`f1_heldout/`, `f1v2/<s>/`, and the gorilla F1v2 dev arm `f1v2/dev/<s>/min.*`).

| sample | `off` vs HEAD (6 products) | F1: GTF, families GTF, junctions | F1v2: GTF, families GTF, bridges, junctions | peak RSS off → on |
|---|---|---|---|---|
| human_testis | identical | identical | identical | 0.79 → 0.88 GB |
| gorilla_OR6737 | identical | identical | identical (dev arm) | 1.18 → 1.82 GB |
| gorilla_KB3781 | — | identical | identical (dev arm) | 1.50 (2026-09-25) → 2.28 GB |
| human_A119b | — | identical | identical | 2.58 (2026-09-25) → 3.09 GB |

- **Junction tables.** They are identical after the batch-order permutation above (`reorder.py`), and identical as
  sorted sets.
- **Evidence.** The evidence counts equal `samtools view -F 2308` spliced records (testis 1,104,114).
- **Buffered path.** On the synthetic fixture, `--materialize-reads` (its own indexed evidence pass) equals streaming.
- **Tests.** `cargo test --release`: 998 passed, 0 failed, 13 ignored (979 before; 16 unit tests and 3 end-to-end
  tests in `tests/copy_assign_bridge_regroup.rs`). The fixture's `off` GTF run through the frozen scripts gives the
  asserted files byte for byte.
- **Cost.** The pass takes 5-16 s per sample (A119b: 16 s, 490 s for the whole run). The memory is about 40 B per
  spliced primary plus the contig sequences for the PAS test, bounded by the genome cache.
