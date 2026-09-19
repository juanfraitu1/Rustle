# isoseq collapse and FLAIR: as comparison arms, and as mechanisms (§6q6, 2026-09-19)

Two separate questions, answered separately.

## 1. As comparison arms — FLAIR now runs on all six chromosomes

Previously FLAIR existed only on chr20 and the six-chromosome tables were StringTie-only.
FLAIR 3.0.0, the documented recipe: `samtools fastq -F 2308` → `flair align -g chrN.fa -r reads.fq` →
`flair collapse -g chrN.fa -q flair.bed -r reads.fq --generate_map`.

⚠**A second, broken FLAIR install on `/mnt/linuxdisk/.../_from_wsl/miniforge3/envs/flair` shadows the
working one** — `flair collapse` shells out to `filter_transcriptome_align.py` by bare name and picked up
the broken copy, dying with `ModuleNotFoundError: No module named 'flair'`. Fix: shims in
`/mnt/linuxdisk/tmp/flair_shims` that exec each `site-packages/flair/*.py` with the active env's python,
prepended to `PATH`. Without them only the chr20 run (from 09-15, before the shadowing) works.

### gffcompare, all six chromosomes

| chrom | tool | mRNAs | chains | chain Sn/Pr | transcript Sn/Pr |
|---|---|---|---|---|---|
| chr20 | **ours** | 658 | **336** | **7.8 / 51.6** | **7.4 / 51.2** |
| | StringTie | 712 | 331 | 7.7 / 47.4 | 7.3 / 47.1 |
| | FLAIR | 820 | 264 | 6.2 / 35.1 | 5.8 / 32.3 |
| chr11 | **ours** | 1,373 | **684** | **7.0 / 50.2** | **6.5 / 50.0** |
| | StringTie | 1,307 | 648 | 6.6 / 50.2 | 6.2 / 50.0 |
| | FLAIR | 1,671 | 589 | 6.0 / 38.6 | 5.6 / 35.5 |
| chr7 | **ours** | 1,172 | **516** | **6.4 / 44.4** | **5.9 / 44.2** |
| | StringTie | 1,199 | 515 | 6.4 / 43.5 | 5.9 / 43.0 |
| | FLAIR | 1,326 | 436 | 5.4 / 37.0 | 5.0 / 33.2 |
| chr14 | **ours** | 904 | **391** | **7.1 / 43.7** | **6.3 / 43.4** |
| | StringTie | 935 | 388 | 7.1 / 42.0 | 6.3 / 41.9 |
| | FLAIR | 1,120 | 306 | 5.6 / 30.1 | 5.0 / 27.8 |
| chr5 | ours | 1,008 | 475 | 6.3 / **47.5** | 5.8 / **47.1** |
| | StringTie | 1,061 | **476** | 6.3 / 45.5 | **5.9** / 45.3 |
| | FLAIR | 1,200 | 383 | 5.1 / 36.5 | 4.8 / 32.3 |
| chr9 | **ours** | 934 | **407** | **5.5 / 44.1** | **5.2 / 43.8** |
| | StringTie | 953 | 403 | 5.5 / 42.8 | 5.2 / 42.8 |
| | FLAIR | 1,142 | 331 | 4.5 / 32.0 | 4.2 / 29.2 |

**Pooled (46,077 reference mRNAs):**

| | emitted | matching chains | matching transcripts | transcript Pr | transcript Sn |
|---|---|---|---|---|---|
| **ours** | **6,049** | **2,809** | **2,817** | **46.6%** | **6.11%** |
| StringTie | 6,167 | 2,761 | 2,785 | 45.2% | 6.04% |
| FLAIR | 7,279 | 2,309 | 2,332 | 32.0% | 5.06% |

**Pooled SQANTI3:**

| | n | FSM | ISM | artifact cats | rules-filter PASS |
|---|---|---|---|---|---|
| **ours** | 6,049 | **2,819 (46.6%)** | 735 (12.2%) | **354 (5.9%)** | **5,705 (94.3%)** |
| StringTie | 6,153 | 2,794 (45.4%) | 770 (12.5%) | 384 (6.2%) | 5,759 (93.6%) |
| FLAIR | 7,279 | 2,393 (32.9%) | **546 (7.5%)** | 858 (11.8%) | 5,227 (71.8%) |

FLAIR is behind on every gffcompare quantity on every chromosome and on the SQANTI3 pass rate
(71.8% vs our 94.3%); its one lead is the lowest ISM fraction (7.5%), which is its `--no_redundant`
subset removal being more aggressive than ours — bought with 858 artifact-category transcripts.

## 2. `isoseq collapse` could NOT be run on this substrate

isoseq 26.2.0 (`mamba create -n isoseq -c bioconda isoseq`). `isoseq collapse` emits 0-1 transcripts from
our 25,341 chr20 primary alignments regardless of `--min-aln-coverage`/`--min-aln-identity`/
`--keep-non-real-cells`. It parses PacBio ZMW read names and PacBio read-group metadata; renaming reads to
`m64011_000000_000000/<n>/ccs` moved it from 0 to 1 transcript, and supplying an `@RG` with
`DS:READTYPE=SEGMENT` made it fail outright (`segment read group is missing SOURCE type`). **The substrate
is an SRA-derived minimap2 BAM, not a native PacBio BAM, so the tool is not runnable here** — this is a
substrate limitation, not a result about isoseq. Recorded so nobody re-attempts it (register row 865).

Its *mechanisms* are fully specified by its own CLI defaults, so they were implemented and measured
instead.

## 3. The mechanisms: what we already had, and the one we did not

| mechanism | tool | our equivalent |
|---|---|---|
| collapse 5'-shorter transcripts missing 5' exons (default; `--do-not-collapse-extra-5exons` disables) | isoseq | **`--assembly-polish full`** — the support-aware ISM collapse (§6p8) |
| 3'/5' redundancy removal (`--no_redundant`) | FLAIR | same ISM collapse |
| minimum isoform read support (`-s`) | FLAIR | `--polish-mono-quantile` floor, `--polish-ism-ratio` |
| minor-isoform fraction of the locus (`-f`) | StringTie | `--polish-isoform-fraction` (§6p9) |
| **fuzzy junction tolerance (`--max-fuzzy-junction 5`)** | isoseq | **nothing — implemented here as `--polish-fuzzy-junction`** |
| 5'/3' end windows (`--max-5p-diff 50`, `--max-3p-diff 100`) | isoseq | subsumed: our loci already collapse by exact chain, so same-chain transcripts are one record before any end comparison |
| splice-junction correction against annotation/short reads | FLAIR (`flair correct`) | not applicable — de novo by construction, and it was skipped in the FLAIR arm too |

### `--polish-fuzzy-junction` — implemented, measured, REFUTED at isoseq's default

Two chains count as the same chain when they have the same number of junctions and every corresponding
donor/acceptor is within N bp; near-duplicates merge into the best-supported member.
`--polish-fuzzy-ism` additionally extends the tolerance to the ISM sub-chain test (what isoseq does).

| tolerance | cells vs StringTie | emitted | matching chains | transcript Pr |
|---|---|---|---|---|
| **0 bp (shipped)** | **28/30** | 6,049 | **2,809** | 46.57% |
| 2 bp, merge only | 28/30 | 6,039 | 2,805 | 46.58% |
| **5 bp, merge only (isoseq default)** | **15/30** | 5,845 | 2,728 | 46.81% |
| 5 bp, merge + ISM (exactly isoseq) | 15/30 | 5,841 | 2,727 | 46.82% |
| 10 bp, merge + ISM | 15/30 | 5,772 | 2,707 | 47.04% |

⛔**5 bp costs 81 matching intron chains for +0.24 points of transcript precision, and halves the
scorecard.** The merge and the containment test cost the same, so **the damage is the merge itself**.

**Why:** our pass-1 already enforces canonical GT-AG motifs (§6m8), so a sub-5 bp junction difference that
survives into our GTF is a *real tandem splice site*, not alignment noise. On chr20, of the 25 chain pairs
that 5 bp merges but 0 bp keeps apart, the modal offset is **3 bp — the NAGNAG acceptor signature** (32.1%
of differing coordinates; the rest are 2, 4 and 5 bp). Merging them deletes whichever variant has less
read support, and that is not always the annotated one. `isoseq collapse` needs the tolerance because it
collapses raw alignments with no motif constraint; we do not. Register row 866.

Both flags ship, default off (0 bp = our historical behaviour, byte-identical).
