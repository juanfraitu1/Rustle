# PREREG r1634: opt-in profiles that move the initial assembly toward StringTie or FLAIR, and their price (2026-11-11)

Goal (user, 2026-11-11): match either tool "almost exactly" to know the precision / recall trade. Step 0 = r1633. Written
before any arm is run. Only existing switches; defaults untouched; no new code.

Measure (`bench/famcmp/tool_match.py`): multi-exon intron chains, Jaccard with the tool's set, share of the tool's chains
reached, share of ours inside the tool's set. Price = RefSeq-correct chains and precision, with chains new and lost against
the default counted apart. Development: gorilla OR6737. Held out: A119b (only if the stop rule lets it run).

Arms (assemble only), from the row-1267 difference model and r1633:
- S1 `RUSTLE_READ_ISOFORM_K=2` (StringTie's single-read choice)
- S2 `RUSTLE_READ_JUNCTION_SNAP=25` (StringTie's 25 bp junction correction)
- S3 S1 + S2; S3f = S3 scored on `PREFIX.final.gtf` (sub-chains dropped)
- F1 `RUSTLE_MIN_CHAIN_READS=1` (r1616's arm, the closest to FLAIR's one-read chains); F1 + `RUSTLE_CHECK_SPLICE=1` = F2
- P `--seed-pool primary` (does it remove the 1,032 ours-only chains with < 2 full-chain primary records?)

"Almost exactly", fixed here: Jaccard >= 0.80 with the tool AND >= 0.90 of the tool's chains reached AND >= 0.90 of ours
inside the tool's set. Default today: 0.442 / 0.639 / 0.588 (StringTie), 0.309 / 0.361 / 0.680 (FLAIR).
- PASS: an arm meets all three on gorilla and on A119b.
- Stop rule: if no arm reaches Jaccard 0.60 on gorilla, A119b is not run; the profile closest to each tool is reported with
  its price, the residual is split by class, and the user is asked once how to continue.

## Result (2026-11-11, `r1634/`): STOP at the registered rule (best Jaccard 0.519 < 0.60); A119b arms not run
Gorilla OR6737, against StringTie (Jaccard / share of its chains / share of ours inside it) and RefSeq:
| arm | chains | correct | precision | Jaccard | reached | inside | new (correct) | lost (correct) |
|---|---|---|---|---|---|---|---|---|
| default | 73,348 | 25,794 | 0.352 | 0.442 | 0.639 | 0.588 | | |
| S1 read-isoform-k 2 | 91,746 | 27,520 | 0.300 | 0.439 | 0.720 | 0.530 | 18,472 (1,730) | 74 (4) |
| S2 junction snap 25 | 66,272 | 23,538 | 0.355 | 0.494 | 0.656 | 0.668 | 3,922 (567) | 10,998 (2,823) |
| S3 = S1 + S2 | 81,144 | 24,840 | 0.306 | 0.494 | 0.729 | 0.606 | 18,834 (1,870) | 11,038 (2,824) |
| S3f = S3, final.gtf | 76,003 | 24,785 | 0.326 | 0.519 | 0.727 | 0.645 | 17,235 (1,858) | 14,580 (2,867) |
| F1 1-read floor | 127,014 | 29,006 | 0.228 | 0.339 | 0.729 | 0.387 | 53,773 (3,217) | 107 (5) |
| P seed-pool primary | 72,274 | 25,634 | 0.355 | 0.446 | 0.638 | 0.596 | | |
FLAIR: no arm beats the default's Jaccard by more than 0.002 (default 0.309, P 0.311); F1 reaches 0.412 of FLAIR's chains at
Jaccard 0.273; F2 (check-splice) equals F1.
Closest to StringTie = S3f: -1,009 correct chains and -0.026 precision for +0.077 Jaccard. Only 6,529 of its 17,235 new
chains are StringTie's.

Residual, gorilla (`tool_match_drop.py`, default vs StringTie): ours-only 30,198 = other combination of shared introns
11,213 (0.143 correct), sub-chain of a StringTie chain 7,745 (0.068), within 25 bp of one 7,485 (0.326), super-chain 2,115,
no shared intron 1,640. StringTie-only 24,341 = other combination 11,197, super-chain of ours 6,380 (0.165), no shared
intron 4,041, sub-chain 2,296, near 427.
A119b read split (`r1633/step0b.a.tsv`): StringTie-only 125,409 = 0 full-chain reads 28,099, 1 read 38,595, 2+ 58,715;
FLAIR-only 507,176 = 1 read 345,191 (0.012 correct). Ours-only vs StringTie 133,061, 121,667 of them with 2+ reads.

## A119b price of the closest StringTie-like profile (2026-11-11, run after the stop, at the advisor tool's suggestion; `r1634/hS3`)
| arm | chains | correct | precision | Jaccard | reached | inside | new (correct) | lost (correct) |
|---|---|---|---|---|---|---|---|---|
| default | 219,517 | 38,715 | 0.176 | 0.251 | 0.408 | 0.394 | | |
| S3 | 243,862 | 36,826 | 0.151 | 0.273 | 0.461 | 0.401 | 56,125 (2,291) | 31,780 (4,180) |
| S3f | 225,096 | 36,591 | 0.163 | 0.286 | 0.459 | 0.432 | 52,390 (2,236) | 46,811 (4,360) |
Price on A119b: -2,124 correct chains and -0.013 precision for +0.035 Jaccard. The stop stands on both samples.
