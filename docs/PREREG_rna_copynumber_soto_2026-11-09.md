# PREREG row 1594: does RNA see Soto's copy number? copy_phasing extra-copy calls vs S1E population parCN classes

Written 2026-11-09, before any script output, copy_phasing run or count of hits. Advisor context: Soto 2025 defines duplicated
genes by SD98 + shared exons and measures their copy number from DNA read depth (QuicK-mer2 parCN). The exon-sharing half
transfers to RNA (rows 1591, 1593). RNA depth is expression, so depth cannot carry copy number; what RNA can count is distinct
expressed sequences. This row asks whether the existing extra-haplotype detector sees the copies that Soto's DNA counts.

**Not a new method.** The detector is the frozen `copy_phasing` certificate of rows 1582-1584 (OR6737 5/7 vs Jim 0/7; GRCh38
panel 10/13 hit but labelled 5/13; HG002 held-out 2/6 FAIL, depth-limited). New here is the truth set: Soto's S1E population
classes, about 50x more labelled genes than those panels.

## Truth (`bench/soto/soto_parCN_S1E.tsv`, Soto Table S1E, 2,504 1KGP genomes)

No RNA individual of ours is in S1E, so only genes whose copy number is (nearly) fixed in the population are labelled:
- **Positive** (CHM13 holds fewer copies than nearly everyone): `parCN class` Fixed or Nearly-Fixed AND `Median parCN` >= 2.5.
  (S1E: 92 genes, all Nearly-Fixed.) parCN counts the gene's paralog-specific k-mers in a diploid genome, so a CHM13 gene at
  median parCN 4 has, in most genomes, a second copy per haplotype carrying its specific sequence that CHM13 lacks.
- **Negative**: Fixed or Nearly-Fixed AND 1.5 <= `Median parCN` < 2.5 (S1E: 1,085 genes).
- Polymorphic genes and Median parCN < 1.5: descriptive only.

Genes are placed on CHM13 v2.0 by S1E `Gene ID` through the CAT table
(`~/work/jfim/rustle_data/annotation/chm13v2.0_CAT_Liftoff.genes.tsv`, all 1,833 IDs present; S1E coordinates are CHM13 v1.0
and are not used), with the gene's exon blocks.

## RNA side (A119b testis Iso-Seq on CHM13 v2.0; blind to the labels)

- Families: A119b's shipped de novo families `~/runs/human_A119b/a1.fam.copies.tsv`, restricted to the families with a copy
  whose exons overlap the exons of any S1E gene (any class). Restriction only saves time: copy_phasing runs per family.
- `bench/o3_graph/real_net.py A119b.t2t.bam SUBSET.copies.tsv cp` (the >= 2-copy families and their net, as rows 1582-1583).
- The frozen binary `~/o2ties_real/copy_phasing.port run --bam ... --copies cp.copies.tsv --genome chm13v2.0.fa --net
  cp.net.tsv --out cp.A --threads 48`, defaults (as row 1583).

## Definitions (fixed now; `bench/soto/rnacn.py`)

- **Testable gene**: a copy of `cp.copies.tsv` (>= 2-copy family) on the same chromosome has exon overlap >= 1 bp with the
  gene's CAT exons.
- **Hit (extra-copy call)**: a `flagged=1` row of `cp.A.candidates.tsv` whose home copy (`nearest_locus`) or a row of
  `cp.A.putative.tsv` whose home copy (`copy`) is a copy overlapping the gene's exons (same rule).
- Sensitivity = hits / testable positives; false-positive rate (FPR) = hits / testable negatives. A hit on a negative is an
  RNA count above the DNA count.

## Bars

- **STOP, no verdict**, if fewer than 15 positives are testable (underpowered).
- **PASS** iff sensitivity - FPR >= 0.15 AND one-sided Fisher exact p < 0.01 (hits vs testable, positives vs negatives).
- FPR is reported against row 1583's single-origin null (7/30 = 0.23); no bar on it alone.

There is no dev/held-out split: nothing is fitted (frozen binary, defaults, labels fixed above). HG002 is not used: its DNA
copy number is not in S1E and only its paternal assembly is on disk.

## Reported, no bar

- Hit rate in the Polymorphic class and at median parCN < 1.5.
- Confounds: median net reads and family size (copies) of testable positives vs negatives; the PASS test repeated within
  families of >= 3 copies.
- Per positive gene: symbol, Soto family, median parCN, testable, hit.

## What follows

- PASS: RNA extra-copy calls carry Soto's copy-number signal at population level; next would be an individual-matched test
  (an RNA sample with its own DNA copy number).
- FAIL: on RNA, Soto's copy-number half is not measurable with the current detector; families stay defined by exon sharing +
  SD, and copy number stays a DNA quantity.

## Outcome (2026-11-09, job 4453941): FAIL

Net 2,171,512 rows over 432 families (1,829 copies, every family >= 2 copies); copy_phasing 2,575 s.

| label | S1E genes | testable | hit | rate | median family copies | median copy reads |
|---|---|---|---|---|---|---|
| positive (CHM13 short of nearly everyone) | 92 | 43 | 18 | **0.419** | 4 | 14 |
| negative (fixed at 2) | 1,085 | 596 | 186 | **0.312** | 4 | 23 |
| Polymorphic | 515 | 233 | 39 | 0.167 | 7 | 12 |
| median parCN < 1.5 | 141 | 68 | 18 | 0.265 | 6 | 10 |

Sensitivity 0.419 - FPR 0.312 = **+0.107** (bar +0.15), one-sided Fisher p = 0.10 (bar 0.01): **FAIL** as registered.
Descriptive: within families of >= 3 copies, 13/36 vs 121/396 (+0.056, p 0.30). FPR 0.31 is above row 1583's single-origin
null (0.23). Positives carry fewer reads per copy than negatives (14 vs 23), so depth does not inflate their hit rate.

Reading: the extra-haplotype certificate fires on about a third of genes Soto calls fixed at 2 copies, and only slightly
more often on genes where CHM13 is short of the population. At population level, with this detector, RNA does not carry
Soto's copy-number signal. Families stay defined by exon sharing + SD (rows 1580, 1593); copy number stays a DNA quantity.
An individual-matched test (an RNA sample with its own DNA copy number) would remove the population-label noise, but the
detector's false-positive rate on fixed genes is the first problem. Per gene: `r1594/per_gene.tsv`.
