# If NPIP and TBC1D3 were ideally expressed, would the current default pipeline find them? (2026-10-06)

Protocol: `docs/PREREG_ideal_expression_2026-10-06.md` (v1 56637089, Amendment 1 d745c6a9). Tools: `bench/ideal_expression/`. Products: `/mnt/linuxdisk/tmp/ideal_expression_2026-10-06/{NPIP,TBC1D3}/rep{1,2}/`.
HEAD release binaries (`copy_assign` 87824d91, `mcl_families` a6308244, `as_table` 48786e8f, `family_score` 542923fd), human CHM13 v2.0, CAT/Liftoff v2.0. **DEV, human only, circular by construction** (the reads come from the annotation that scores them):
the numbers say whether the default pipeline can find these copies when expression, truncation, readthrough and depth do not limit; they do not predict recovery from a real library.
Two independent read replicates per family (seeds 20261006 and 20261007): NPIP 20,660 reads from 2,066 simulated transcripts of 474 genes (547 genes and 2,147 transcripts lie in the windows), TBC1D3 16,680 reads from 1,668 transcripts of 464 genes (519 genes, 1,725 transcripts in the windows) (every CAT/Liftoff transcript within +-500 kb of a copy, 10 full-length jittered reads each, error .001).

## Answer

| family | copies | reachable (R) | IDEAL-FOUND on R, rep1 / rep2 | rule | all copies | E2 / E3 / E4 on R | cluster precision |
|---|---|---|---|---|---|---|---|
| **TBC1D3** | 16 | 13 | 12 / 12 | **YES** (both replicates) | 14 / 16 | 13 / 12 / 13 | .938 |
| **NPIP** | 25 | 14 | 12 / 12 | **NO** (both replicates; one copy short of the 13 needed) | 16 / 25 | 12 / 13 / 14 | .909 |

IDEAL-FOUND = the copy owns a locus (both purities >= .5, no other copy shares it), that locus's representative has exactly one of the copy's chains, and the locus is in the family's cluster. **R** = copies that are not entangled with another annotated gene (>= 100 shared exonic bp),
have a canonical multi-exon chain, are not readthrough images, and have an observable chain (>= 3 reads of the seeding pool): fixed from the truth tables and the BAM before the pipeline ran (R sha1 83ab9229 for NPIP, e43717f7 for TBC1D3; identical in both replicates).
YES needs >= ceil(0.9 n) on R and cluster precision >= .5. All four arms pass G1, G2, G3, G4 and G6 (VALID). Predictions: P1 (TBC1D3 YES) held, P2 (NPIP not YES, and <= 21 of 25 over all copies) held, P3 (replicates agree) held.

The excluded copies, printed beside (identical in both replicates): **NPIP** entangled 10 (IDEAL-FOUND 4: h05, h08, h10, h11; three of the ten are also readthrough images), chainless 1 (NPIPA3, found 0), R 14 (12). **TBC1D3** entangled 2 (found 2), chainless 1 (TBC1D3P7, found 0), R 13 (12).
A control that gives the same scorer PERFECT loci (the canonicalized annotation through `mcl_families`, G6) reaches 14 / 14 on NPIP's R (17 / 25 over all copies) and 13 / 13 on TBC1D3's R (15 / 16): relative to perfect loci the pipeline loses 2 reachable NPIP copies and 1 reachable TBC1D3 copy.

## Why the reachable copies are missed (NPIP two, TBC1D3 one per replicate)

- **NPIPA2 (h01), both replicates:** its two true chains are assembled (2 / 2 recovered) but the locus also holds three artifact chains built from reads of paralogous copies, and the most-reads representative is one of those (6 exons, representative purity .454): E2 and E3 fail. The representative problem of `docs/SPLICED_COPY_SUPPORT_2026-10-04.md` exists with ideal reads.
- **NPIPA7 (h06), both replicates:** the representative is exact (purity 1.0, 7-intron chain), but the locus (gene_id) also contains nine transcripts, most of them 21-30 exon chains of the neighbouring pseudogenes (AC138969.2, PKD1P6): locus purity .112, so E2 fails. What the families stage and O2 receive (the representative) is right.
- **TBC1D3P4 / TBC1D3P3 (h26, h27), one per replicate:** near-identical copies with one 11-intron chain each. The representative is a coin toss (10 reads against 10) between the copy's own chain and a chain with one intron displaced by about 25-31 bp. The only reads that carry the displaced chain are the paralog's (h27's reads at h26 in replicate 1, h26's reads at h27 in replicate 2), all 10 as secondary alignments at this copy: the aligner put one junction off in the paralog's sequence, and the good-secondary seeding admitted them.
  The assembled transcripts contain every chain of every TBC1D3 copy (31 / 31 chains recovered, both replicates), so this is representative choice, not assembly. Which of the two copies fails flips between replicates.

## Same copies, same registered instrument, real reads against ideal reads

The real-read column is the current default on A119b (CAT/Liftoff copies, instrument of `docs/PREREG_spliced_copy_support_2026-10-04.md` Amendments A and B, `bench/copy_support.py`). Only these annotation-anchored measures are like-for-like; the registered strict FOUND needs the cap signal and cannot run on simulated reads.

| instrument | NPIP real | NPIP ideal (rep1 / rep2) | TBC1D3 real | TBC1D3 ideal (rep1 / rep2) |
|---|---|---|---|---|
| copies with a locus overlapping them | 25 / 25 (own node 24) | 25 / 25 (own node 24) | 13 / 16 (no node file) | 16 / 16 (own node 15) |
| annotation-anchored FOUND (Amendment A) | 12 | 23 / 23 | 12 | 16 / 16 |
| exact-chain FOUND (Amendment B) | 6 | 21 / 21 | 9 | 15 / 15 |
| locus level, annotation-anchored | 23 | 24 / 24 | 13 | 16 / 16 |
| locus level, chain | 14 | 24 / 24 | 12 | 16 / 16 |

On the same 25 NPIP copies the annotation-anchored FOUND goes from 12 to 23 and the exact-chain FOUND from 6 to 21 when the reads are ideal; for TBC1D3 12 to 16 and 9 to 15. The real-read shortfall is therefore mostly the data, with the representative problem left over at NPIP.

## Observability and chains (E0, E1)

NPIP: 24 of 25 copies have an observable chain; 142 of 148 chains are observable (>= 3 reads of the seeding pool; 131 with primary alignments alone); the primary alignment of a read lands on its own copy for 84% of reads (12 copies below 90%: the near-identical paralogs tie), yet junction-level observability is intact.
Chains recovered exactly as an assembled transcript: 133 / 148 (rep1), 129 / 148 (rep2); the unrecovered chains are at NPIPA1 (6 / 8), NPIPB5 (17 / 21), NPIPB4 (17 / 21) and AC138894.1 (7 / 12), all entangled copies. TBC1D3: 15 of 16 copies have an observable chain, 31 / 31 chains observable and recovered, 98% of reads land on their own copy.

## Family level (beside; NOT a ceiling)

`family_score` is keyed by RefSeq genes while the reads are CAT/Liftoff (RefSeq NPIPB3 has no CAT gene, and several NPIP symbols differ), and its universe here is the simulated windows, so these numbers are not the real-read run's and not an upper bound.
NPIP family (Compara CF153, 19 genes in the universe): 11 in the NPIP cluster, sens .579, prec .917, F .710 (U2 ID_154: 11 / 19, F .611 and .595); pooled over the nine Compara families F .733 (.667 without size-2 clusters). TBC1D3 family (Compara CF185, 9 genes): 9 / 9 in the cluster, prec .900, F .947; pooled F .753 and .737 (.720 and .685 without size-2 clusters).

## Controls

G1 0 violations in 2,066 and 1,668 simulated transcripts (12,857 and 8,093 introns), re-read from the genome by a separate script. G2 100% (NPIP) and 99.9% (TBC1D3) of the reads of 275 and 304 single-copy genes land on their gene. G3 one record per FASTQ read, driver exit 0, `families.gtf` newer than `gtf`. G4 scorer tables identical under two hash seeds.
**G5** an independent scorer written from the pre-registration text alone (it never opened `score.py`) reproduces E2, E3, E4 and IDEAL-FOUND for all 82 copy-arms (25 + 25 + 16 + 16) with no difference; a line-by-line code review found no blocker (one chromosome-check bug, found and fixed during the first scoring, and several consistency fixes made after the results: unmapped reads in the G2/E0 denominators, unrounded purity thresholds, E3p over all locus transcripts,
cluster precision without a strand clause, the folded-holder count; every per-copy flag was identical before and after). **G6** perfect loci give E2 = E3 = E4 = 100% on R in all four arms. **G7 (single-copy genes, >= 95% with exactly one locus) is NOT attainable as registered:** the annotation-as-loci control reaches 89.5% (NPIP) and 94.1% (TBC1D3) because genes overlapping other genes cannot each get one locus,
and the pipeline arm reaches 63-65% and 62% (many single-copy genes are single-exon and receive no locus at ten reads); it is reported against the control, not as a gate.

## Robustness of the NPIP verdict

The independent scorer computed twelve alternative readings of definitions the text leaves open (copy exon union against copy span, folded-record handling, strand clauses, cluster-precision variants, tie order, ...). Eleven keep NPIP at NO; one (purity measured against the copy's span instead of its exon union) gives PARTLY; none gives YES.
Post hoc, and not a verdict: if E2 judged the representative only (the representative is what the families stage and O2 receive), NPIPA7 would pass and NPIP would be 13 / 14, the YES threshold. NPIP's verdict therefore sits one copy below the bar and TBC1D3's exactly on it (12 of 13 needed 12).

## Held-out contrast (no bars)

SMG1P (7 genes in the NPIP windows): 2 reachable, 1 found; KRTAP (33 genes in the TBC1D3 window): 30 are single-exon, 2 reachable, 0 found (not grouped in one family). Too few reachable copies to say more.

## Per-copy results (replicate 1; replicate 2 differs only where the last column shows it)

### NPIP

| copy | name (CAT/Liftoff) | strata | chains (observable) | chains recovered | reads on own copy | holder purity rep / locus | E2 | E3 | E4 | IDEAL-FOUND rep1 / rep2 | if not found |
|---|---|---|---|---|---|---|---|---|---|---|---|
| h00 | NPIPB2 | R | 8 (8) | 8/8 | 1.00 | 1.0 / 1.0 | 1 | 1 | 1 | 1 / 1 |  |
| h01 | NPIPA2 | R | 2 (2) | 2/2 | 1.00 | 0.454 / 0.535 | 0 | 0 | 1 | 0 / 0 | locus |
| h02 | NPIPA1 | E with CHM13_G0020712 (3881 bp, 53%) | 8 (8) | 6/8 | 1.00 | 1.0 / 0.674 | 1 | 0 | 1 | 0 / 0 | annotation:E |
| h03 | PKD1P6 | E+X with CHM13_G0020726 (1888 bp, 31%) | 6 (6) | 6/6 | 1.00 | 1.0 / 0.716 | 1 | 1 | 0 | 0 / 0 | annotation:E |
| h04 | NPIPA5 | R | 3 (3) | 3/3 | 1.00 | 1.0 / 0.835 | 1 | 1 | 1 | 1 / 1 |  |
| h05 | AC138969.1 | E+X with CHM13_G0020764 (5447 bp, 68%) | 12 (12) | 12/12 | 0.83 | 1.0 / 0.837 | 1 | 1 | 1 | 1 / 1 |  |
| h06 | NPIPA7 | R | 2 (2) | 2/2 | 0.65 | 1.0 / 0.112 | 0 | 1 | 1 | 0 / 0 | locus |
| h07 | NPIPA8 | E with CHM13_G0020801 (1227 bp, 76%) | 3 (3) | 3/3 | 1.00 | 0.906 / 0.129 | 0 | 1 | 1 | 0 / 0 | annotation:E |
| h08 | AC138969.1 | E+X with CHM13_G0020802 (1568 bp, 25%) | 9 (9) | 9/9 | 0.86 | 1.0 / 0.648 | 1 | 1 | 1 | 1 / 1 |  |
| h10 | NPIPB5 | E with CHM13_G0020899 (2061 bp, 22%) | 21 (19) | 17/21 | 0.56 | 1.0 / 0.606 | 1 | 1 | 1 | 1 / 1 |  |
| h11 | NPIPB4 | E with CHM13_G0020920 (386 bp, 5%) | 21 (20) | 17/21 | 0.70 | 1.0 / 0.537 | 1 | 1 | 1 | 1 / 1 |  |
| h12 | NPIPB3 | E with CHM13_G0020934 (386 bp, 10%) | 13 (12) | 13/13 | 0.70 | 1.0 / 0.234 | 0 | 1 | 1 | 0 / 0 | annotation:E |
| h13 | NPIPB6 | R | 2 (2) | 2/2 | 1.00 | 1.0 / 0.864 | 1 | 1 | 1 | 1 / 1 |  |
| h14 | NPIPB8 | R | 2 (2) | 2/2 | 0.65 | 1.0 / 1.0 | 1 | 1 | 1 | 1 / 1 |  |
| h15 | AC138894.1 | E with CHM13_G0021066 (2480 bp, 59%) | 12 (10) | 7/12 | 1.00 | 0.944 / 0.477 | 0 | 0 | 1 | 0 / 0 | annotation:E |
| h16 | NPIPB9 | R | 2 (2) | 2/2 | 1.00 | 1.0 / 0.928 | 1 | 1 | 1 | 1 / 1 |  |
| h17 | NPIPB10P | R | 1 (1) | 1/1 | 1.00 | 1.0 / 1.0 | 1 | 1 | 1 | 1 / 1 |  |
| h18 | NPIPB11 | R | 2 (2) | 2/2 | 1.00 | 1.0 / 1.0 | 1 | 1 | 1 | 1 / 1 |  |
| h19 | NPIPB12 | R | 8 (8) | 8/8 | 0.67 | 1.0 / 0.746 | 1 | 1 | 1 | 1 / 1 |  |
| h20 | NPIPB13 | R | 5 (5) | 5/5 | 0.46 | 1.0 / 0.656 | 1 | 1 | 1 | 1 / 1 |  |
| h21 | NPIPA3 | C | 0 | - | 1.00 | 1.0 / 0.702 | 1 | 0 | 1 | 0 / 0 | annotation:C |
| h22 | NPIPB14P | E with CHM13_G0022124 (1678 bp, 97%) | 3 (3) | 3/3 | 1.00 | 0.373 / 0.27 | 0 | 0 | 1 | 0 / 0 | annotation:E |
| h23 | NPIPB15 | R | 1 (1) | 1/1 | 0.70 | 1.0 / 1.0 | 1 | 1 | 1 | 1 / 1 |  |
| h24 | NPIPB15 | R | 1 (1) | 1/1 | 0.60 | 1.0 / 1.0 | 1 | 1 | 1 | 1 / 1 |  |
| h25 | NPIPB15 | R | 1 (1) | 1/1 | 0.70 | 1.0 / 1.0 | 1 | 1 | 1 | 1 / 1 |  |

### TBC1D3

| copy | name (CAT/Liftoff) | strata | chains (observable) | chains recovered | reads on own copy | holder purity rep / locus | E2 | E3 | E4 | IDEAL-FOUND rep1 / rep2 | if not found |
|---|---|---|---|---|---|---|---|---|---|---|---|
| h26 | TBC1D3P4 | R | 1 (1) | 1/1 | 1.00 | 0.971 / 0.972 | 1 | 0 | 1 | 0 / 1 | representative |
| h27 | TBC1D3P3 | R | 1 (1) | 1/1 | 1.00 | 1.0 / 0.977 | 1 | 1 | 1 | 1 / 0 |  |
| h28 | TBC1D3P5 | R | 3 (3) | 3/3 | 1.00 | 1.0 / 1.0 | 1 | 1 | 1 | 1 / 1 |  |
| h29 | TBC1D29P | R | 5 (5) | 5/5 | 1.00 | 1.0 / 1.0 | 1 | 1 | 1 | 1 / 1 |  |
| h30 | TBC1D3B | E with CHM13_G0024042 (1383 bp, 49%) | 4 (4) | 4/4 | 0.90 | 1.0 / 1.0 | 1 | 1 | 1 | 1 / 1 |  |
| h31 | TBC1D3I | R | 3 (3) | 3/3 | 1.00 | 0.993 / 0.92 | 1 | 1 | 1 | 1 / 1 |  |
| h32 | TBC1D3G | R | 1 (1) | 1/1 | 1.00 | 0.974 / 0.727 | 1 | 1 | 1 | 1 / 1 |  |
| h33 | TBC1D3H | R | 1 (1) | 1/1 | 1.00 | 0.979 / 0.818 | 1 | 1 | 1 | 1 / 1 |  |
| h34 | TBC1D3B | R | 4 (4) | 4/4 | 0.78 | 1.0 / 1.0 | 1 | 1 | 1 | 1 / 1 |  |
| h35 | TBC1D3K | R | 1 (1) | 1/1 | 1.00 | 0.998 / 0.781 | 1 | 1 | 1 | 1 / 1 |  |
| h36 | TBC1D3D | R | 1 (1) | 1/1 | 1.00 | 0.998 / 0.781 | 1 | 1 | 1 | 1 / 1 |  |
| h37 | TBC1D3E | R | 1 (1) | 1/1 | 1.00 | 0.98 / 0.767 | 1 | 1 | 1 | 1 / 1 |  |
| h38 | TBC1D3 | R | 1 (1) | 1/1 | 1.00 | 0.973 / 0.817 | 1 | 1 | 1 | 1 / 1 |  |
| h39 | TBC1D3P7 | C | 0 | - | 1.00 | 1.0 / 1.0 | 1 | 1 | 0 | 0 / 0 | annotation:C |
| h40 | TBC1D3P1 | E with CHM13_G0024955 (182 bp, 11%) | 1 (1) | 1/1 | 1.00 | 1.0 / 1.0 | 1 | 1 | 1 | 1 / 1 |  |
| h41 | TBC1D3P2 | R | 3 (3) | 3/3 | 1.00 | 1.0 / 1.0 | 1 | 1 | 1 | 1 / 1 |  |


## Departures from the plan, and limits

Stratum X (readthrough images) is not exclusive: NPIP h03, h05, h08 are X and also E. The scorer was changed once after the first scoring (a missing chromosome check, which made G7 return 0 and could have affected cluster precision; every E2, E3, E4 and IDEAL-FOUND flag was unchanged) and again after the independent review (consistency fixes above, flags unchanged).
The legacy 2026-09-21 reads were not used (register 1253). The held-out families, the family-level scores and the like-for-like instruments were computed after the verdicts and are reported beside them. Not shown: real reads, O2 (about 16% of NPIP reads have their primary alignment on a paralogous copy even when every copy is expressed, so contested reads exist in the ideal case and O2 could be tested there),
O3, gorilla, other families, and anything about the annotation's correctness. A YES or NO here is about the algorithm at these loci: 10 of the 25 NPIP copies overlap another annotated gene and cannot each get a locus from any assembler.
