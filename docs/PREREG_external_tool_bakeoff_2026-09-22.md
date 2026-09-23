# Pre-registration — is our `--assemble-only` step as good as isoseq collapse, StringTie and FLAIR, on the lab's own genome-wide runs?

**Written 2026-09-22 (§6z5), before any new number is computed.** User: *"I have external data from
stringtie, flair and isoseq collapse, can we check the results and use that info to ensure our assembly
step is as good as the ones for these tools?"* — data in `../benchmark_collapse/` (StringTie 3.0.1 and
FLAIR 3.0.1, both species) and `../isoseq_upload/` (isoseq collapse, both species).

## What is already known, and therefore NOT a prediction

- Human A119b, ours vs **StringTie 3.0.3 run here** (`bench/GENOME_WIDE_BAKEOFF_2026-09-22.md`): ours wins
  every chain-level metric (intron chain 22.2/16.8 vs 18.8/15.8; 38,859 vs 32,871 matching chains),
  StringTie wins base sensitivity (59.1 vs 48.6) and misses fewer loci (39.9% vs 44.4%).
- Human chr20/21/22 panel, FLAIR (our bamToBed route): 23.1/6.4 intron chain vs ours 23.3/16.0.
- The lab's own tool-vs-tool gffcompare (no annotation, no ours) already shows the shape of the three
  tools: isoseq emits 3.2 M transcripts (human) / 551 k (gorilla), FLAIR 1.1 M / 153 k, StringTie
  250 k / 68 k. So isoseq will have the highest sensitivity and the lowest precision of any arm, by
  construction.
- The six-chromosome testis panel (§6p8-§6q6) is on a **different library** (r867) and is not quoted here.

## Substrates — same BAMs the lab used, no exceptions

| species | BAM (identical `@PG` line to the lab's) | genome | annotation |
|---|---|---|---|
| human | `A119b.t2t.bam` (68 M records, 25 contigs) | `chm13v2.0.fa` | `chm13v2.0_RefSeq_full.gff.gz` → `ref_genome.gtf` (188,903 mRNAs) |
| gorilla | `GGO_mm.bam` (10.7 M records, 26 contigs) | `GGO.fasta` (GCF_029281585.2) | `GGO_genomic.gff` → GTF via `gff_to_gtf` (**106,508 transcripts**; the converter keeps every transcript type, as it did for human) |

⚠ Human and gorilla are reported in **separate tables and never pooled.** ⚠ The human ours-vs-StringTie
number above is on **StringTie 3.0.3**; the lab's file is **3.0.1**. Both are reported; they are not the
same arm.

**Ours** = the genome-wide `--assemble-only` sweep already on disk for human
(`/mnt/linuxdisk/tmp/gw22/ours_genome.gtf`, shipped polish, `RUSTLE_JUNCTION_MAJORITY` default ON) and
the same recipe run now for gorilla. **No flag is changed for this measurement.**

## Metrics (gffcompare v0.12.10, query = tool, reference = annotation)

Per species, per tool: intron-chain sensitivity, intron-chain precision, transcript sensitivity,
transcript precision, matching intron chains — **5 metrics × 3 tools = 15 cells**. A cell is "won" if ours
≥ the tool on that metric. Also reported (not scored): base/exon/intron level, missed/novel loci, and the
lab-style pairwise tool-vs-tool agreement with ours as the fourth arm (`summarize.sh`'s framework).

⚠ Metric-trap guard: sensitivity is denominated on the FULL annotation (188,903 / 80,997 mRNAs) in every
arm — no universe intersection, no per-tool restriction. ⚠ The gorilla annotation has 80,997 mRNAs but 106,508 transcripts after conversion (lncRNAs etc.); the human 188,903 is the same converter's output, so the two are like-for-like. Precision is per emitted transcript, so the
3.2 M-transcript isoseq arm is penalised for volume exactly as it should be. No pooling across species.

## The bar — committed now

| outcome | verdict |
|---|---|
| ≥ 12/15 cells in **both** species, and ours has the most matching intron chains of any arm in both | ⭐ **AS GOOD OR BETTER** |
| ≥ 12/15 in one species only, or most matching chains in one only | ⚠ **SPECIES-DEPENDENT** — say which and why |
| < 12/15 in both | ⛔ **NOT AS GOOD** |

**Predicted, before looking:** ⚠ at best. Ours should win every precision cell against all three tools and
the sensitivity cells against StringTie (already known on human), but **isoseq will beat us on both
sensitivity cells** because it emits 12× more transcripts, and FLAIR will tie or edge us on sensitivity
(panel: 23.1 vs 23.3). That is 11-13/15 on human. Gorilla is **unmeasured for every arm including ours**
and is the held-out species for any change this work motivates.

## The second question — "use that info to improve the assembly step"

Only asked if a sensitivity gap to isoseq/FLAIR exists. Diagnosis rule, fixed now: take the reference
transcripts that isoseq (or FLAIR) matches exactly (`=` in its tmap) and ours does not, and split them into
(a) **absent from our raw pre-polish output** — an assembly-stage loss — versus (b) **present raw, removed by
the polish** — a filter loss. Only (a) is "the assembly step". Any remedy is then pre-registered separately,
developed on human, and judged on gorilla. ⛔ Already-refuted dials are not re-proposed: fuzzy junctions
(r866), `--read-isoform-k 3` (r855), polish fractions/escapes (r858-r864), relaxing canonicity (§6m8).

I will not change the metrics, the cells, the substrates or the bar after seeing any number.

---

# OUTCOME (2026-09-22) — gorilla measured here; human FLAIR/isoseq arms go to the cluster (user's call)

⚠ The two human files (FLAIR 1.1 M transcripts / 1.0 GB GTF, isoseq 3.2 M / 1.8 GB) are scored on the cluster,
not on this box (user, 21:26). What was shipped for that: `../benchmark_collapse/assembler_{A119b,GGO}/*.assembler.gtf.gz`,
`ref/*.gtf.gz` (the same converter output as here, md5s in `ref/md5.uncompressed.txt`), `score_vs_annotation.sh`
(the 15-cell table) and `ASSEMBLER_ARM.md`. Until those land, the human row is ours-vs-StringTie only.

Gorilla `--assemble-only` sweep: 26 contigs, 196 s, 0 failures, 78,602 transcripts. The 18:17 binary reproduces this
morning's human chr21 GTF byte-for-byte (74 s), so the human arm on disk is the current code.

## Gorilla (`GGO_mm.bam`, RefSeq GCF_029281585.2, 106,508 transcripts) — query = tool, reference = annotation

| arm | transcripts | loci | intron SN/PR | **intron chain SN/PR** | **transcript SN/PR** | **matching chains** | matching loci | missed loci | novel loci |
|---|---|---|---|---|---|---|---|---|---|
| **ours** | 78,602 | 17,398 | 59.6 / 81.7 | **27.0 / 33.1** | **24.3 / 32.9** | **25,829** | 12,422 | 60.2% | 9.4% |
| StringTie 3.0.1 (lab) | 68,249 | 19,406 | 63.9 / **83.1** | 24.4 / **34.6** | 22.0 / **34.3** | 23,320 | 12,601 | 57.1% | 11.3% |
| FLAIR 3.0.1 (lab) | 153,402 | 23,077 | 62.1 / 61.5 | 25.5 / 17.6 | 23.1 / 16.0 | 24,351 | 11,634 | 58.5% | 29.6% |
| isoseq collapse (lab) | 551,342 | 86,860 | **75.6** / 40.5 | **30.8** / 6.7 | **28.0** / 5.4 | **29,426** | 13,236 | **40.2%** | 73.1% |

**Cells won by ours: 10/15** — vs StringTie 3/5 (both precision cells lost by 1.4-1.5 pts), vs FLAIR 5/5, vs isoseq
2/5 (both sensitivity cells and the chain count lost). **Most matching chains: isoseq (29,426 vs our 25,829).**
By the bar that is **< 12/15 on gorilla**, so the overall verdict cannot be ⭐; it is ⚠ or ⛔ depending on the human row.

## Where isoseq's extra chains come from — read support, not assembly

Exact-matched reference transcripts (`=` in each tmap): ours 25,869 · StringTie 23,395 · FLAIR 24,565 · isoseq 29,802;
union of all four arms 36,077 (33.9% of the annotation).

| isoseq's exact matches | n | median FL reads | FL = 1 | FL ≤ 2 | FL ≥ 5 |
|---|---|---|---|---|---|
| shared with ours | 21,332 | **6** | 10.2% | 24.6% | 58.5% |
| **isoseq-only** | 8,470 | **1** | **77.8%** | 89.3% | 5.1% |

⭐ **The 8,470 reference transcripts isoseq matches and we do not are single-read isoforms in 77.8% of cases**
(6,588 at FL = 1; isoseq's whole output is 68.8% FL = 1). Our assembler needs 2 reads at each terminal
(`min_terminal_support = 2`), so a one-read chain is unreachable by construction, and r975/§6w6 measured that
k = 1 is a no-op on the polished output. 84.9% of them sit at a gene where we do emit a transcript (996 are class
`j`, a different junction combination), so this is isoform-level, not locus-level.

The reverse set is the mirror image: of the 4,537 references **we** match and isoseq does not, isoseq has a
**contained fragment (`c`) at 4,119 (90.8%)** — it emitted a truncated chain where we emitted the full one.
StringTie-only 2,813 vs ours-only 5,287; FLAIR-only 2,861 vs ours-only 4,165.

## Human chr20/21/22 (development contigs; the whole-genome FLAIR/isoseq rows come from the cluster)

Same three contigs, same reference subset (10,851 transcripts), every arm cut to those contigs.

| arm | transcripts | intron SN/PR | **intron chain SN/PR** | **transcript SN/PR** | **matching chains** | novel loci |
|---|---|---|---|---|---|---|
| **ours** (shipped polish) | 15,876 | 61.0 / 56.3 | **23.3 / 16.0** | **21.8 / 14.9** | **2,295** | 55.2% |
| StringTie 3.0.1 (lab) | 14,169 | 62.7 / 55.9 | 18.9 / 14.8 | 17.4 / 13.3 | 1,863 | 54.5% |
| FLAIR 3.0.1 (lab) | 54,765 | 66.5 / 28.2 | 23.1 / 6.3 | 21.3 / 4.2 | 2,268 | 79.7% |
| isoseq collapse (lab) | 169,422 | 76.9 / 15.0 | 27.7 / 2.4 | 26.3 / 1.7 | 2,723 | 88.1% |

**12/15** (loses only isoseq's two sensitivity cells and its chain count). Human genome-wide, ours vs the
lab's StringTie 3.0.1 (250,500 transcripts): 22.2/16.8 vs **18.1/15.0** intron chain, 38,859 vs 31,792
matching chains — the lab's 3.0.1 file scores slightly below the 3.0.3 run here (18.8/15.8, 32,871).

## The second question, answered — what the tools match that we do not

Diagnosis rule as pre-registered (tool `=` matches absent from ours, split by our RAW pre-polish output):

| species / tool | tool-only exact matches | absent from RAW (assembly-stage) | present RAW (polish removed) | tool's read support, assembly-stage losses |
|---|---|---|---|---|
| gorilla / isoseq | 8,470 | **84.0%** | 16.0% | median FL **1**; FL = 1 in 88.6% |
| gorilla / FLAIR | 2,861 | 74.5% | 25.5% | median 5 (FLAIR's post-hoc count) |
| gorilla / StringTie | 2,813 | 90.4% | 9.6% | `longcov` ≤ 1 in 75.0% |
| human chr20-22 / isoseq | 867 | 72.7% | 27.3% | median FL 1; FL = 1 in 75.1% |
| human chr20-22 / FLAIR | 453 | 66.9% | 33.1% | median 7 |
| human chr20-22 / StringTie | 250 | 79.2% | 20.8% | `longcov` ≤ 1 in 56.1% |

Counting reads in the BAM that carry the EXACT reference chain (controls: 99-100% of shared and
polish-removed refs have ≥ 2): of the FLAIR-only / isoseq-only refs absent from raw, **86.8% / 88.8% have
exactly one such read** (gorilla FLAIR-only: 94.5%). Pass-1 needs two, so those are unreachable by
construction and FLAIR's "median 5-7 reads" is quantification, not chain evidence. The **13.2% / 11.0%
residue with ≥ 2 exact-chain reads (75 human cases) is entirely the coordinate de-duplication defect of
`docs/PREREG_primary_dedupe_2026-09-22.md`** — confirmed, fixed behind `--keep-coordinate-duplicates`, and
refuted as a default because the shipped polish was fitted on the de-duplicated counts (polished chains fall
on both substrates; gorilla precision rises 3.2 pts, cells 10 → 12/15).

**Answer:** the assembly step is not losing anything the tools assemble from ≥ 2 exact-chain reads except
through the de-duplication key; everything else the tools add is single-read isoforms (isoseq: 68.8% of its
gorilla output is FL = 1) and the price of that is their precision (isoseq 6.7%, FLAIR 17.6% vs our 33.1%
on gorilla).

## What is missing (2026-09-23) — taxonomy of every missed multi-exon reference transcript

For each multi-exon reference transcript our polished output does not match exactly: the relation of our
closest transcript on the strand, whether the exact chain exists in our RAW output, how many primary BAM reads
carry the exact chain, and which lab tools match it (`miss_taxonomy.py`; gorilla is a 3,000-transcript sample).

| substrate | class of our closest transcript | share of misses | exact-chain reads 0 / 1 / ≥2 | in RAW | matched by any tool |
|---|---|---|---|---|---|
| human chr20-22 (7,446 misses) | other (different junction set) | 40.2% | 86 / 9 / 5% | 3.4% | 13.0% |
| | sub-chain (ours shorter) | 29.5% | 85 / 11 / 4% | 2.9% | 14.9% |
| | none (no spliced transcript there) | 18.8% | 96 / 4 / 0% | 0.1% | 3.7% |
| | super-chain (ours longer; ref is a contained isoform) | 4.2% | 53 / 11 / **36%** | 34.0% | 34.0% |
| | retained intron | 3.3% | 83 / 12 / 5% | 3.2% | 16.7% |
| | junction shift ≤ 10 bp | 2.3% | 59 / 29 / 11% | 10.3% | 33.1% |
| gorilla (3,000 of 69,987) | other | 41.6% | 80 / 17 / 3% | 2.3% | 18.5% |
| | none | 31.4% | 94 / 6 / 0% | 0.0% | 5.0% |
| | sub-chain | 18.7% | 82 / 16 / 2% | 1.4% | 18.3% |
| | super-chain | 4.4% | 64 / 11 / **25%** | 24.8% | 26.3% |
| | junction shift | 2.9% | 61 / 32 / 7% | 4.6% | 23.0% |

**Misses with ≥ 2 exact-chain reads: 384 human (5.1%), 83 of 3,000 gorilla (2.8%).** Of those, 303 / 74 are in
RAW (the polish removed them: fraction 170, ISM 110, both 23 on human) and the human remainder (81) is the
r1058 de-duplication. Everything else — ≥ 95% of what we miss — has at most one read carrying the chain,
and the tools match it only at singleton level (isoseq). Among the ≥ 2-read super-chain cases the container's
extra introns lie at the 5′ side in 54 and the 3′ side in 56 (human), 19 / 12 (gorilla); the ≥ 2-read junction
shifts are the 3-bp NAGNAG case (16 of 20). Priced and refuted on this residue: an absolute-support exemption
from the fraction rule (r1063) and an internal-priming-guarded ISM exemption for 3′-shorter isoforms (r1064).
The one lever left is the polish re-fit on true counts (`docs/PREREG_polish_refit_true_counts_2026-09-23.md`).

## Two better universes than the whole annotation (2026-09-23, user's point)

The annotation holds transcripts that are not in this library. Two denominators that are:

**(a) Annotation restricted to what the reads carry** (≥ 2 / ≥ 1 primary reads with the exact chain; ours-matched
transcripts counted as ≥ 2, verified 99–100% in the r1057 controls; gorilla scaled from the 3,000-miss sample):

| universe | ours | isoseq | FLAIR | StringTie |
|---|---|---|---|---|
| human chr20-22, ≥ 2 reads (2,750) | **86.0%** | 83.0% | 74.3% | 62.4% |
| human chr20-22, ≥ 1 read (3,457) | 68.4% | **80.6%** | 65.9% | 52.8% |
| gorilla, ≥ 2 reads (27,805) | **93.0%** | 82.6% | 80.7% | 75.7% |
| gorilla, ≥ 1 read (37,160) | 69.6% | **80.7%** | 65.3% | 62.5% |

The tools also match reference transcripts with ZERO supporting reads (StringTie 43 human / 397 gorilla, FLAIR
1 / 140, isoseq 2 / 140) — annotation leakage through fuzzy ends, not recovery.

**(b) The tools' own consensus, annotation-free** (strand-aware intron chains, multi-exon):

| substrate | chains in all 3 tools | ours has | in our RAW | we miss: in RAW / absent from RAW | our chains in no tool | tool-exclusive chains (ST / FLAIR / isoseq) |
|---|---|---|---|---|---|---|
| human chr20-22 | 4,178 | 3,132 (75.0%) | 77.4% | 103 / 943 | 1,502 (10.5%) | 21.6% / 25.3% / 68.9% |
| gorilla genome-wide | 31,657 | 29,073 (**91.8%**) | 93.2% | 431 / 2,153 | 4,592 (**5.9%**) | 10.6% / 22.2% / 67.3% |

What we lack from the all-3 consensus, absent from raw:

| | human 943 | gorilla 2,153 |
|---|---|---|
| exactly one supporting read | 273 | 1,794 (83%) |
| non-canonical junction (CCAG, CTAG, GTCT, … shared by all tools because they share the alignments; §6m7's artifact class) | 388 | 239 |
| canonical with ≥ 2 exact reads ("buildable") | 389 | 237 |
| … of which ALL exact reads share one coordinate (the r1058 de-duplication signature) | 389 (100%) | 237 (100%) |
| … present once duplicates are kept (raw / shipped polish) | 389 / 276 | – / 218 |

So the tool-consensus gap is singletons, non-canonical junctions, and coordinate-identical reads that our key
counts as one molecule and the tools count as several. Whether the last class is PCR duplicates or distinct
molecules with identical ends is exactly what the polish re-fit (`PREREG_polish_refit_true_counts_2026-09-23.md`)
arbitrates: if no true-count setting recovers the de-duplicated operating point, the key is a legitimate
duplicate filter and the tools' agreement on those chains is duplicate-inflated.
