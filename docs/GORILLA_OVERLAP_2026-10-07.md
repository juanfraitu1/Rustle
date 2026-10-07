# Overlapping genes on real reads, gorilla OR6737 testis: our assembler against StringTie, FLAIR and isoseq collapse (2026-10-07)

**Status: descriptive real-read check on a partly spent library (cross-individual: reads OR6737, assembly and annotation KB3781), RefSeq annotation as truth, one library. Protocol: `docs/PREREG_gorilla_overlap_2026-10-07.md` (committed 120fc627, before any truth table or arm existed). Tools: `bench/entangled/{real_truth.py, real_score.py, run_gorilla.sh}` (14 unit tests, `python3 -m unittest test_real_overlap`, commits 9d8f6b16 / 89beddd1). Products: `/mnt/linuxdisk/tmp/gorilla_overlap_2026-10-07/` (truth tables, the three arms, scores).**
User direction: 'lets do 1' (overlapping genes on real gorilla reads). Arms: ours D (default), P (`--no-seed-secondaries`), PC (P + `--polish-subchain drop`), HEAD release binaries, genome-wide `assemble` (D 4 min with the AS table, P and PC 1.7 min); the lab's S (StringTie 3.0.1 `-L`), F (FLAIR 3.0.1 collapse) and I (isoseq collapse) on the same BAM. D reproduces the 09-25 transcript count exactly (74,045), P the 09-25 primaries-only count (72,690).

## Answer

1. **On real reads our default recovers more of the expressed annotated chains than any of the three tools, at overlapping genes and elsewhere.** Primary denominator (chains with >= 3 primary reads, what the tools see): same-strand overlapping genes with both genes expressed (E_both, 108 genes, 202 chains) D 170 (84.2%), I 159 (78.7%), S 140 (69.3%), F 140 (69.3%); non-overlapping genes (N, 10,079 genes, 19,771 chains) D 94.6%, F 86.2%, I 81.2%, S 78.8%. `gffcompare` agrees on E_both (intron-chain sensitivity D 84.1, F 69.7, S 69.2).
2. **Overlap costs every method, FLAIR most**: E_both against N: D 84.2 vs 94.6, S 69.3 vs 78.8, I 78.7 vs 81.2, F 69.3 vs 86.2. Antisense overlap (A, 1,583 genes): D 93.8%, F 85.9%, I 84.8%, S 75.2%.
3. **Gene separation (E_both genes resolved)**: D 74 / 108, F 70, S 52, I 13; without shared junctions (E_both_x, 100 genes) D 71, F 67, S 47, I 13. StringTie fuses (146 gene_ids carrying exact chains of two genes against D's 61 and FLAIR's 4), isoseq collapse gives one id per cluster (951). Only 8 of the 108 E_both genes share a junction with their partner: the RefSeq gorilla annotation has none of the CAT / Liftoff readthrough and exon-reusing models that made 84% of the human ideal-window entangled genes inseparable.
4. **The two precision options cost almost nothing here and gain nothing visible**: P and PC change E_both by 0 and 1 chains (P1), N by 2 and 18; the annotation-match share moves from 42.4% to 42.5% / 43.8% (N). On the secondary denominator (>= 3 reads of the seeding pool) P loses 2 chains at E_both and 102 at N (0.5%): the seeding gain is visible genome-wide, small at overlaps.

## Results (chains = expressed annotated chains recovered exactly; P1 denominator unless stated)

| stratum (genes, chains) | D | P | PC | S | F | I |
|---|---|---|---|---|---|---|
| E_both (108, 202) | 170 (84.2%) | 170 | 169 | 140 (69.3%) | 140 (69.3%) | 159 (78.7%) |
| E_both_x (100, 188) | 162 (86.2%) | 162 | 161 | 131 (69.7%) | 136 (72.3%) | 148 (78.7%) |
| E_both_j (8, 14) | 8 | 8 | 8 | 9 | 4 | 11 |
| E_one (176, 329) | 299 (90.9%) | 299 | 298 | 251 (76.3%) | 251 (76.3%) | 266 (80.9%) |
| A (1,583, 3,491) | 3275 (93.8%) | 3275 | 3274 | 2625 (75.2%) | 3000 (85.9%) | 2960 (84.8%) |
| N (10,079, 19,771) | 18702 (94.6%) | 18700 | 18684 | 15586 (78.8%) | 17034 (86.2%) | 16062 (81.2%) |
| transcripts | 74,045 | 72,690 | 69,449 | 68,249 | 153,402 | 551,342 |
| annotation-match share, E_both / N | 43.2 / 42.4% | 43.3 / 42.5% | 43.2 / 43.8% | 40.4 / 43.4% | 19.6 / 19.5% | 16.1 / 14.4% |

Secondary denominator (E_both 110 genes, 207 chains; N 10,216 genes, 19,953 chains): E_both D 175 (84.5%), P 173, PC 172, S 142, F 141, I 162; N D 94.4%, P 93.8%, PC 93.8%, S 78.3%, F 85.4%, I 80.8%. Evaluated genes: 11,946 (P1) and 12,094 (P2) of 41,193; 93,564 distinct valid annotated chains of 29,312 genes; 23,793 (P1) and 23,993 (P2) expressed. Strata (P1): N 10,079, A 1,583, E_one 176, E_both 108 (E_both_j 8).

## Predictions

**G1 held** (D above S and F at E_both and at N, both denominators). **G2 failed narrowly at E_both** (P loses 2 chains on the secondary denominator, 3 were required) **and holds at N** (102). **G3 failed**: the annotation-match share of D is not below S or F (E_both 43.2% against 40.4% and 19.6%; N 42.4% against 43.4% and 19.5%), and PC (43.8% at N) is above both; the over-emission seen on the ideal reads (precision 85% against 92%) does not show in this proxy, whose denominator holds every unannotated isoform. **G4 held** (E_both_x D 71 against S 47). **G5 held** (I 159 against D 170 at E_both).

## Limits

One library, cross-individual (private variants drop exact-chain reads: the denominators are lower bounds for every arm alike); the RefSeq annotation is partly Gnomon and may cite long reads (if OR6737 is among them the truth is partly circular, unverified); the annotation-match share is a proxy, not a precision; the lab tools are in their default annotation-free recipes (versions 3.0.1, local ones differ) and were never tuned; E_both_j has 8 genes (14 chains), too few for a conclusion; no E_both result is a pass or fail of anything except G1 to G5; OR6737 is not held-out (the seeding rule and F1/F1v2 were decided with it).
