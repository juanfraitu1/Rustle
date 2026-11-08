# PREREG: combined collapse rule (certificate on the copy's own catalog copy OR S/H pile), held out on HG002 (register 1584, 2026-11-08)

Follows rows 1538 and 1583. The certificate catches collapses on a heterozygous survivor. The S/H pile catches a pure second
haplotype on a homozygous survivor. This row combines them, ties the certificate to the copy instead of the landing window, and
tests the result on a second individual.

## Rule (blind, `bench/nont2t/combined_rule.py`)
- **Units.**
  - These are r1538's loci (`psv_pile.py` on the sample): absorbing loci of expressed scorable missing copies, and single-origin
    null loci.
  - A null locus is single-origin when its own gene is the top labelled origin of the window's primaries and no other gene reaches
    10 % of them.
  - Identical absorbing windows count once.
  - Stratum autosomal; chrX/Y is reported only.
- **Certificate component.**
  - The unit's own copy is the de novo catalog copy (`cp.copies.tsv`) whose exons hold the most labelled primary reads of the
    gene inside the window.
  - It fires when a `flagged=1` candidate or a putative copy of `copy_phasing.port` (frozen binary, default flags) is homed on
    exactly that copy.
- **S/H component.** S or H of the locus exceeds the 95th percentile over the autosomal single-origin null loci of the same sample
  (own null).
- **Combined.** Either component fires.

## Development (A119b, all parts already seen; recorded before the held-out exists)
`/anvme/workspace/v138da13-rustle/nont2t/human/cp/combined_dev.txt`. The own-null thresholds are S 6.13, H 0.454.

| A119b autosomal | n | certificate | S/H | combined | certificate with purity >= 0.5 (labelled) |
|---|---|---|---|---|---|
| missing | 13 | 10 | 3 | 10 (0.77) | 5 |
| single-origin null | 32 | 6 | 2 | 8 (0.25) | 1 |

- With the frozen r1538 post hoc thresholds (S 4.64, H 0.431), null is 10/32 (0.31): FAIL.
- Tying the certificate to the copy did not remove the neighbours' haplotypes. They are homed on the same absorbing copy.
- The labelled purity column separates missing from null: 5/13 against 1/32. The blind certificate does not.

## Held-out (HG002 Kinnex on GRCh38)
- **Data prep.** Slurm 4452739, `/anvme/workspace/v138da13-rustle/nont2t/hg002h/`, launched before this file and blind to it.
  - It follows the r1504 recipe on `HG002_kinnex.t2t.bam` (minimap2 2.31 with the same flags; `de` and `=`/`X` present).
  - The panel and `copies.paf` are reused, so the truth status is CHM13 copies on GRCh38, as for A119b.
  - Expression and origin are recomputed from HG002.
- **Differences from A119b, disclosed.**
  1. Panel reads are the primaries in panel exons, fetched by index; A119b also took reads whose only panel hit was a secondary.
     The BAM is 147 GB and a full name scan would take about 9 h.
  2. Only the assemble and families stages run.
  3. There are 15,613 panel primaries against A119b's 79,226, so no read cap is applied. HG002 Kinnex is not testis, and fewer
     missing copies may be expressed.
- **Confound.** HG002 is not CHM13. A copy CHM13 has may be absent or extra in HG002, so some "missing" and null labels can be
  wrong for this individual. Null certificates can be checked post hoc against `hg002.pat.fa`; that check is not barred.

## Bars (HG002 autosomal)
1. Combined hits >= 50 % of the missing units.
2. Combined hits <= 25 % of the single-origin null units.
3. If either has fewer than 5 units, the verdict is UNDERPOWERED and nothing is closed.

**Reported, not barred:**
- each component alone;
- frozen A119b thresholds (S 6.13, H 0.454);
- labelled-correct certificates (purity >= 0.5), missing against null;
- X/Y.

**PASS:** the combined rule is a blind collapse detector across individuals. Next it would be wired as an opt-in report in
`missing_copy_flag`.

**FAIL on bar 2:** the certificate fires on present tandem clusters in a second individual too. The blind rule is not a detector,
and the remaining question is a blind proxy for the labelled purity.

**FAIL on bar 1:** the collapsed regime stays below half on RNA.

## RESULT r1584 (2026-11-08): FAIL as registered (bar 1); depth-limited
- **Run.** Slurm 4452739, 19 min. 15,613 panel reads; 12 expressed missing copies, of which 6 are scorable.
  - The de novo catalog is small: 12 families with >= 2 copies, against 60 on A119b.
  - Own-null thresholds: S 2.27, H 0.448, from n = 7. Output `/anvme/workspace/v138da13-rustle/nont2t/hg002h/cp/combined_heldout.txt`.

| HG002 autosomal | n | certificate | S/H | combined | labelled-correct certificate |
|---|---|---|---|---|---|
| missing | 6 | 2 | 2 | 2 (0.33) | 1 (NPIPB11, purity 1.00) |
| ... of which collapsed | 3 | 1 (NPIPB4, purity 0.43) | 1 | 1 | 0 |
| ... of which absent | 3 | 1 | 1 | 1 | 1 |
| single-origin null | 7 | 0 | 0 | 0 (0.00) | 0 |

- **Bar 1 FAILS** (2/6 < 0.5). **Bar 2 PASSES** (0/7). Frozen A119b thresholds give the same counts.
- **The misses are coverage misses.**
  - NPIPB9, NPIPB7 and NPIPB15 have no position at DP >= 20 (S = nan).
  - NPIPB7 and NPIPB15 also have no catalog copy at the locus ("-").
  - GOLGA6L9 is again the pure second haplotype on a homozygous survivor (H 0.375, under 0.448), the same miss as on A119b.
- **Reading.** No false call in a second individual: 0/7 nulls. Reach stays at about a third, as on A119b's labelled 5/13.
  - Neither component reaches half of the collapsed copies in either sample.
  - The collapsed regime remains below half from RNA.
  - HG002 Kinnex expresses too few of these genes to say more. The n is small (6 / 7).
