# PREREG row 1655: multiple representatives per locus, transfer to gorilla (OR6737 and KB3781)

## Why
Row 1654: under the exon-sharing default, multiple representatives pass all four bars on HG002 and fail the precision bar
on A119b (0.896 vs 0.918; recall +0.038). User 2026-11-12: test it on gorilla to see if it is useful. Neither gorilla
sample has been run through `locus_reps.py`; nothing in the procedure is changed from row 1654.

## Arms (binary `target_sddefault/release/mcl_families`, GGO.fa, damped with `sedef_dl/gorilla/final.bed`)
Flags: `--min-exonic-bp 1 --emit-units --min-shared-exon-frac 0.0 --min-exon-share 0.5 --inflation 2.0 --sd-damp-bed BED`.
- base: `or.families.gtf` / `kb.families.gtf` unchanged (rerun here: the stored tables predate the default).
- F2: `locus_reps.py --disjoint` -> families -> piece_merge -> piece_attach --keep-all onto base. Fa descriptive.
Scoring: pairscore --species gorilla (symbol-inherited truth, simiiformes, PAIRSCORE_SAME_SYMBOL=ignore), universe =
the sample's families.gtf. Descriptive: the same with DENOVO_GGO_ENS_MAP (row 1611 protein-hit ids for LOC genes).

## Bars: F2 vs base, each sample (as row 1654)
1. Precision >= base.  2. Recall >= base + 0.02.  3. Base hits lost <= 5.  4. Copies <= 1.05 x base.
All four on both samples -> proposed as an opt-in families post-step. Otherwise FAIL as registered.

## Prediction
OR6737 (testis, like A119b): bars 2-4 pass, bar 1 borderline FAIL. KB3781 (fibroblast, fewer expressed copies): bars
1, 3, 4 pass, bar 2 uncertain.

## Result (2026-11-12, job 4467154, 603 s, `/anvme/workspace/v138da13-rustle/r1655_multirep_gorilla/`)
**FAIL as registered on both samples; the registered truth is too small to decide much.** Fa = F2 on both samples.

| sample | truth | arm | pred pairs | true | false | precision | recall | hits lost | copies |
|---|---|---|---|---|---|---|---|---|---|
| OR6737 | symbol (88 pairs) | base / F2 | 35 / 35 | 18 / 18 | 17 / 17 | 0.514 / 0.514 | 0.205 / 0.205 (bar 2 FAIL) | 0 | 1,339 / 1,357 |
| KB3781 | symbol (68 pairs) | base / F2 | 19 / 21 | 12 / 13 | 7 / 8 | 0.632 / 0.619 (bar 1 FAIL) | 0.176 / 0.191 (bar 2 FAIL) | 0 | 694 / 708 |
| OR6737 | + protein-hit ids (1,012), descriptive | base / F2 | 2,095 / 2,138 | 753 / 768 | 1,342 / 1,370 | 0.359 / 0.359 | 0.744 / 0.759 | | |
| KB3781 | + protein-hit ids (238), descriptive | base / F2 | 464 / 480 | 97 / 107 | 367 / 373 | 0.209 / 0.223 | 0.408 / 0.450 | | |

Prediction wrong on OR6737 (no recall gain at all on the registered truth). The layer changes little on gorilla: 18 and
14 copies added (1.3%, 2.0%), against 287 on A119b (4.4%). pairdiag (symbol truth): the absorbed class is 1 pair on each
gorilla sample (A119b 71, HG002 42); the gorilla misses are unlinked (OR 64 of 70, KB 50 of 56), which multiple
representatives do not reach. With the protein-hit truth the added pairs are OR 15 true / 28 false, KB 10 true / 6 false.

## Amendment 1 (2026-11-12, post hoc, user "can we check that"): why the protein-hit precision is low
The descriptive protein-hit truth counts as false every pair of gorilla genes that best-hit the SAME human gene: human
Compara has no pair for one gene, and the same-symbol exemption (row 1530) keys on names, which LOC genes do not share.
OR6737 base: 1,238 of 1,342 false pairs (0.92) are such pairs; KB3781 315 of 367 (0.86). The rest: OR 47 older Compara
paralogs + 57 other; KB 9 + 43. With same-human-gene pairs set aside (undecidable by this truth, not counted either way):

| sample | arm | true | other false | precision |
|---|---|---|---|---|
| OR6737 | base / F2 | 753 / 768 | 104 / 116 | 0.879 / 0.869 |
| KB3781 | base / F2 | 97 / 107 | 52 / 56 | 0.651 / 0.656 |

Multi-rep on gorilla is then +15 true / +12 false (OR) and +10 / +4 (KB): the registered FAIL stands. pairscore is not
changed; an exemption for pairs sharing a mapped human gene needs its own row (the same-gene pairs are probably real
gorilla duplicates, but nothing here verifies that).
