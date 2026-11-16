# PREREG row 1756: `sd_core` label on family pairs (SD-held and exon-matched through the SD alignment), opt-in (2026-11-15)

User (2026-11-15): "yes lets add that rule as opt in and test". Source: row 1755 post hoc, annotation level on CHM13 (pairs that
pass both the span rule and CORR >= 0.5 were young-true 0.777, any-age 0.903). A label on pairs of copies inside one family; it
changes no family. Rule frozen as in row 1755, nothing tuned here.

**Rule.** Two copies of one family are `sd_core` when one SEDEF pair (a) holds each of them by the `--sd-damp` span rule (>= half
the exon bases and >= half the span inside its copy) and (b) projects, through its CIGAR, exon bases of one onto exon bases of the
other for >= 0.5 of the smaller copy's exon bases. Tool `bench/denovo_eval/sdcore.py`, driver stage `sdcore` (opt-in, named only)
-> `PREFIX.fam.sdcore.tsv`.

**Test.** Scored as `pairscore.py` scores (copies mapped to genes, protein-coding pairs in the universe, same-symbol pairs set
aside, truth = Compara simiiformes among expressed genes). A gene pair is flagged when any copy pair behind it is `sd_core`.
Within-family scored gene pairs are split into flagged / not flagged.

- **Decisive: HG002** on its own paternal genome with its own SEDEF build (`r1741_rdna_hg/d/hg.fam.copies.tsv`). The label was
  derived on CHM13 annotation; HG002's genome, SD calls and annotation are separate.
  - bar 1: precision of flagged pairs >= 0.90;
  - bar 2: precision of flagged pairs >= precision of unflagged pairs + 0.05;
  - bar 3 (reach): flagged pairs hold >= 0.20 of the family pairs that are true.
  PASS = all three. FAIL: the stage stays opt-in as a descriptive column and is not recommended for the advisor report.
- Supporting, same three numbers, no verdict: A119b (`r1749_fam_darklist/base`; same genome as row 1755, not independent), gorilla
  OR6737 and KB3781 (`r1733_shiftdrop/{or,kb}/base`, gorilla SEDEF; truth is small, counts reported as they are).

Also reported, deciding nothing: precision against related pairs (any-age Compara / HGNC group), and the count of flagged pairs
per sample.
