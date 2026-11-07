# PREREG row 1576: held-out of the locus-piece attach rules on HG002 (user option 1)

## Why
Rows 1573 (positional pieces, attached) and 1575 (exon-disjoint multiple representatives, attached) were designed on
A119b after seeing its 25 absorbed genes. HG002 (Kinnex, paternal haplotype, `~/runs/hg002_kinnex/hg.*`) is a held-out
sample: different individual, technology and genome; its absorbed genes are not known yet. Registered before row 1575's
A119b result and before any HG002 run.

## Arms (HEAD mcl_families, minimap2-gate, `--fasta ~/data/hg002/hg002.pat.fa`, the families stage flags)
- A: `hg.families.gtf` unchanged.
- D (row 1573's procedure): B = locus_pieces.py on all loci -> families -> piece_merge (Bm); S = gene_ids with >= 1
  piece in Bm; C = A's input with only S loci replaced by their pieces -> families -> piece_merge (Cm);
  D = piece_attach.py(A copies, C, Cm).
- F (row 1575's procedure): locus_reps.py --disjoint -> families -> piece_merge (Fm); F = piece_attach.py(A, F, Fm).
Scoring: pairscore --species hg002 (simiiformes, same-symbol ignored), pairdiag --species hg002, universe `hg.families.gtf`.

## Bars, per arm vs A (as row 1571)
1. Of A's absorbed missed pairs, <= half still missed.  2. Precision >= A - 0.01.  3. Recall >= A.  4. Copies <= 1.05 x A.
An arm that passes all four is proposed as an opt-in families post-step (`RUSTLE_*` switch, default off).

## Prediction
HG002's truth is small (414 pairs at base) and its absorbed class unknown: D bars 2-4 as on A119b (bar 2 borderline),
bar 1 FAIL; F unknown (depends on row 1575).

## Result (2026-11-07, jobs 4452575 / 4452604, `/anvme/workspace/v138da13-rustle/r1576/`)
**Both arms FAIL bar 1 only.** Arm A at HEAD = the October HG002 run within one pair (0.466 / 0.345, 143 hits; 414 truth
pairs); 42 absorbed missed pairs.

| arm | precision | recall | hits | absorbed still missed (bar <= 21) | A hits lost | copies (bar <= 4,269) |
|---|---|---|---|---|---|---|
| A | 0.466 | 0.345 | 143 | 42 | - | 4,066 |
| D (r1573 procedure) | 0.507 | 0.457 | 189 | 25 FAIL | 9 | 4,139 |
| F (r1575 procedure) | 0.514 | 0.490 | 203 | 25 FAIL | 0 | 4,155 |

Bars 2-4 PASS for both, with gains far above A119b's (F: precision +0.048, recall +0.145). Neither arm meets the
registered opt-in rule (all four bars). Prediction right on bar 1 for D; wrong on D's bar 2 (no precision cost here).
Reading: the multiple-representative layer (F) is a large, collateral-free gain on the held-out sample; it does not
recover most absorbed pairs (17 of 42), whose remaining losses are the two mechanisms named in row 1575.
