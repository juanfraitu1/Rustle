# The reference-absent-copy chain on the Y ampliconic genes — HELD-OUT: gorilla OR6737 testis on mGorGor1 `_pri`, 2026-10-02

Prereg: `docs/PREREG_yag_isocon_chain_2026-10-01.md` (158795c0). Nothing was tuned on this substrate; the rules, delta_Y's definition,
the masking rule and the flag floor were frozen on the human DEV run (`docs/YAG_CHAIN_HUMAN_2026-10-01.md`). Script
`bench/rna_allele/yag_test.py --sub gorilla`; work dir `/mnt/linuxdisk/tmp/rna_allele/yag_ggo/`; outputs in
`docs/YAG_CHAIN_GORILLA_score.out.txt`. Annotation `GGO_genomic.gff` (families by description; chrY NC_073248.2, 67.4 Mb).

## What the gorilla Y offers (the panel rule applied blind)

- Protein-coding copies on chrY: TSPY 6, RBMY 13, DAZ 2, HSFY 5, BPY2 2; CDY 1, PRY 0, VCY 0, XKRY 0 (dropped: < 2 protein-coding copies).
- OR6737 testis primaries per copy: DAZ1 200, DAZ3-like 20, RBMY best copy 45 (next 22, 14, ...), TSPY <= 8, HSFY <= 4, BPY2 <= 1. The
  masking rule (>= 20 primaries, most-expressed copy; a second if >= 4 eligible) yields **two masked copies: DAZ1 and RBMY LOC129530241**;
  TSPY, HSFY and BPY2 stay in the panel unmasked. n = 2: the held-out is as thin as this library's germ-cell expression.
- **delta_Y,ggo = 0.00383** (p99 of d over 99 IsoCon transcripts of 7 X-degenerate single-copy genes: DDX3Y 55, KDM5D 18, ZFY 9, RPS4Y1
  8, UTY 6, USP9Y 2, AMELY 1; median 0.00043, p95 0.00129). No chimeric tail this time; the registered statistic is 10x smaller than the
  human one.
- d_min: **DAZ1 0.1019** to the surviving DAZ3-like copy (prediction: FLAG); **RBMY LOC129530241 0.00108** to its nearest surviving RBMY
  copy (prediction: silent).

## Registered result (delta_Y,ggo = 0.00383)

- 34 IsoCon outputs (DAZ net 220 reads, RBMY 177, TSPY 41, HSFY 13; BPY2 1 read, no run) -> 19 flagged -> 7 linked back -> 12 new-copy
  contigs, all DAZ1-derived -> one flag (>= 2 transcripts): **DAZ1**.
- **Y1 HOLDS: 2/2** — DAZ1 flagged as predicted; the RBMY copy not flagged as predicted (its transcripts link to the surviving RBMY copies
  within delta_Y; its 45 reads tie among them in both arms).
- **Y2 HOLDS: false moves 10/207 S reads = 4.83%** (bar 5%; close, on 207 reads).
- DAZ1's reads (arm R -> arm M, right / wrong / unplaced): **0/200/0 -> 185/9/6** — without the chain all 200 sit on the DAZ3-like copy
  at 10% divergence; with it 185 are on DAZ1's candidate.
- Reported beside: at the autosomal delta (0.00958) the same (2/2); at p95 (0.00129) an RBMY transcript 0.13% away becomes a flag (1/2,
  fails) — the same ordering as in the human run.

## Reading

- The identifiability prediction transfers: on a second species, library and annotation, the chain flags the masked Y copy exactly when
  its divergence from the nearest survivor exceeds the individual-vs-reference divergence of the chromosome, and stays silent otherwise.
  Human 10/11, gorilla 2/2, under rules frozen before either substrate was looked at.
- What the Y adds to the autosomal story: the human palindrome arms put 7 of 11 copies at d_min = 0 — below any floor — and the one
  human exception (DAZ3) is structural; the gorilla DAZ pair, 10% apart, is the easy case and behaves like an autosomal family.
- Power: the gorilla held-out has two masked copies, because OR6737's testis library carries <= 8 reads per TSPY / HSFY copy. A deeper
  gorilla testis library (or KB3781 testis, which would also bring the matched Y assembly) is what a stronger held-out needs.
