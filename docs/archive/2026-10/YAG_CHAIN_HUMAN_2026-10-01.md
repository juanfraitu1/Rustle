# The reference-absent-copy chain on the human Y ampliconic genes (DEV: A119b on CHM13v2.0), 2026-10-01

Prereg: `docs/archive/2026-10/PREREG_yag_isocon_chain_2026-10-01.md` (commit 158795c0, before any run). Script `bench/rna_allele/yag_test.py`; work dir
`/mnt/linuxdisk/tmp/rna_allele/yag_hsa/`; outputs copied to `docs/YAG_CHAIN_HUMAN_score.out.txt`. Gorilla OR6737 testis is the
held-out substrate (separate document).

## delta_Y (computed before the panel ran)

IsoCon on the nets of 13 X-degenerate single-copy Y genes with >= 20 primaries (DDX3Y 11,358 reads ... RPS4Y2 22) gave 383 transcripts
(12 genes aligned); d = 1 - matches / length against the unmasked genome: **median 0.00037, p95 0.00464, p99 = 0.03893 = delta_Y as
registered.** The p99 is carried by partially aligned outputs (ZFY max 0.22, UTY 0.11, USP9Y 0.05, KDM5D 0.04: chimeras and unaligned
ends counted as divergence by the whole-length d), not by donor-vs-reference substitutions, which the median puts at ~0.04%. The
registered statistic is used for the verdict; p95 and the autosomal delta (0.00958) are reported beside.

## The deletion panel (11 copies, 8 families, one masked genome; 240,583 bp masked)

| masked copy | primaries | d_min to the nearest surviving copy (canonical transcript) |
|---|---|---|
| TSPY3 (9.61 Mb) | 30 | 0.00053 |
| RBMY1A1/RBMY1B (22.39 Mb) | 259 | 0.00000 (= RBMY1D) |
| RBMY1F (23.00 Mb) | 123 | 0.00000 |
| DAZ3 | 354 | 0.00000 (= DAZ1) |
| DAZ4 | 201 | 0.00307 |
| CDY1/CDY1B (26.43 Mb) | 108 | 0.00086 |
| CDY2B | 23 | 0.00000 (= CDY2A) |
| BPY2 | 30 | 0.00000 (= BPY2B) |
| HSFY1 | 1,107 | 0.00000 (= HSFY2) |
| PRY2 | 31 | 0.00000 (= PRY) |
| VCY1B | 26 | 0.00185 |

Seven of eleven masked copies have a surviving palindrome arm with an IDENTICAL canonical transcript; the other four are 0.05-0.31% away.
Every d_min is below delta_Y under every definition (registered 0.039, p95 0.0046, autosomal 0.0096), so the registered prediction is
**no flag anywhere**: the human YAG copies sit below the chain's identifiability floor by construction (register 572 / 647 said so for
TSPY and DAZ; this says it for all eight families, as a prediction tested).

## Registered result (delta_Y = 0.03893)

- 282 IsoCon outputs -> 119 flagged (not in the masked genome at 0.999) -> 111 linked back -> 8 new-copy contigs (7 from DAZ3, 1 from a
  TSPY pseudogene) -> flags (>= 2 transcripts): **1, DAZ3**.
- **Y1 HOLDS: the prediction is right for 10/11 masked copies (90.9%, bar 80%).** The miss is DAZ3: d_min = 0 (its canonical transcript
  equals DAZ1's), yet five full-length DAZ3 transcripts (2.1-3.9 kb, support 2-4 reads) align to DAZ1 at **81-87% identity in the block**
  — transcripts with a DAZ-repeat exon arrangement no surviving copy carries. Structural copy differences, not substitutions, make DAZ3
  visible (the same lesson as the gorilla NPIP loss, `project_ggo_npip_loss`: exonic deletions / SVs, not SNVs). The d_min covariate,
  one canonical transcript, cannot see it.
- **Y2 HOLDS: false moves 34/2,303 S reads = 1.48%** (bar 5%).
- Reads of the masked copies (arm R -> arm M, right / wrong / unplaced): DAZ3 0/55/299 -> **171/28/155**; every other masked copy
  unchanged, its reads on the identical surviving arm (HSFY1 500 on HSFY2, CDY1 108 on CDY1', VCY1B 26 on VCY) or tied (TSPY3 30, BPY2 30,
  RBMY1A1 239 of 259).

## Reported beside: the smaller deltas

| delta | new-copy contigs | flags (>= 2 transcripts, D-derived) | Y1 |
|---|---|---|---|
| 0.03893 (registered p99) | 8 | DAZ3 | 10/11 HOLDS |
| 0.00958 (autosomal) | 19 | DAZ3, DAZ4 | 9/11 HOLDS |
| 0.00464 (p95) | 32 | DAZ3, DAZ4, RBMY1A1/1B | 8/11 FAILS |

At the autosomal delta DAZ4 joins (transcripts at 97-99% in-block identity, 2-3% whole-length d); at p95 an RBMY1A1 transcript 1.3% away
(support 3) becomes a flag. The ordering is the point: the flags appear exactly where the deleted copy's transcripts are structurally
or substantially different from every survivor, and the palindrome arms (HSFY, PRY, BPY2, CDY2, RBMY1F) never produce one.

## Reading

- On the human Y the chain behaves as predicted: it stays silent where copies are identical to a surviving arm and speaks only where a
  deleted copy's transcripts differ in structure (DAZ3 at the registered delta; DAZ3 + DAZ4 at the autosomal one). Y1 and Y2 hold.
- IsoCon's own YAG claim (separating transcripts a few substitutions apart) is not what this chain uses: the link step treats any output
  within delta of a reference locus as the individual's variant of that locus. On a haploid chromosome with palindrome arms > 99.9%
  identical, that is the right abstention for copy assignment (register 877) and the wrong tool for counting copies a few SNVs apart;
  the chain's unit of detection on the Y is structural divergence.
- The registered delta_Y is a weak estimator (p99 of a whole-length d over few genes is set by chimeric outputs). For the gorilla
  held-out it is applied as registered; a better definition (d over the aligned block, or the median plus a consensus-error allowance)
  belongs in a new prereg, not here.
- Caveats: one donor, one library; d_min from one canonical transcript per copy; A119b's tissue is not documented (it expresses the
  germ-cell YAGs); CAT and RefSeq name the TSPY copies differently.
