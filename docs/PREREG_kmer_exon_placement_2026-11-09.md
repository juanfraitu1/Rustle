# PREREG row 1599: k-mer counter with exon-based unique placement (row 1598 rerun, one definition changed)

Written 2026-11-09, after row 1598 (FAIL: reach 0.42-0.43 in all nine arms, 75% of unassigned reads with no SPEC hit),
before any count under the new definition. User: "ok lets try once more".

Post hoc reading of row 1598, tested here: rows 1596 and 1598 called a read uniquely placed when its MAPQ 60 primary
overlapped a copy's genomic span (`locus_start..locus_end`). That includes reads in introns, flanks and nested genes,
which hold no copy exon sequence and so cannot carry a SPEC k-mer.

## The one change

- **Exon overlap** of a primary record with a copy = aligned bases (`get_blocks`) inside the copy's exons.
- **Unique-placement read** (truth): primary MAPQ 60 with exon overlap >= 100 bp on exactly one copy of the family, and
  0 bp on the family's other copies.
- **Tied read**: primary with exon overlap >= 100 bp on >= 1 copy of the family, and MAPQ < 60 or >= 2 such copies.
- **Span-only**: the primary overlaps a copy span but has < 100 bp of exon overlap on every copy. Reported, out of every
  bar.

Everything else is row 1598's, unchanged:
- reads in scope;
- SPEC k-mers at each k;
- the arms (k {21, 25, 30} x m {1, 2, 3});
- the chromosome split;
- the dev rule (highest reach with agreement >= 0.98; ties larger k, then larger m);
- the held-out bars: agreement >= 0.98, reach >= 0.70, Spearman >= 0.90. Spearman uses each copy's k-mer assigned reads
  (all in-scope reads) against its unique-placement reads under the new definition;
- the STOP rule and the diagnostic (none / weak / conflict).

## What follows

- PASS: the counter is valid where SPEC k-mers exist; the KB3781 tied-read test becomes the next row.
- FAIL: the counter is closed for good on reference-derived k-mers. Families keep exon sharing + SD, and copy number
  stays a DNA quantity.

Outputs: `/anvme/workspace/v138da13-rustle/r1599/`.

## Outcome (2026-11-09, job 4454075, 4,927 s; read pass 75 min, slower filesystem): FAIL

Dev: every arm has agreement >= 0.988 and reach 0.668-0.677. Dev choice k = 21, m = 1 (reach 0.6767).
Held-out: agreement 0.9903, reach **0.6733** (bar 0.70), Spearman **0.7706** (bar 0.90): **FAIL** as registered.
(The log line says "VERDICT 1598": a hardcoded label in `rnakmer.py arms`; this is row 1599's verdict.)

Unassigned unique-placement reads at the chosen arm:

| half | none | conflict | weak |
|---|---|---|---|
| dev | 31,474 | 73,048 | 0 |
| held | 13,383 | 51,264 | 0 |

Buckets, (read, family) pairs:

| half | unique | tied | span-only |
|---|---|---|---|
| dev | 362,780 | 636,403 | 850,969 |
| held | 220,680 | 157,589 | 238,479 |

Reading:
- Row 1598's post hoc reading holds. Moving unique placement from span to exons cut the "none" class about 7-10x and
  raised reach from 0.43 to 0.67.
- It still misses the bar by 0.027, and Spearman stays at 0.77.
- What is left is mostly **conflict**: about 75% of unassigned reads carry another copy's reference SPEC k-mer. That
  fits A119b alleles that differ from CHM13 at the paralogous sites (gene conversion, or PSVs that are polymorphic).
- As registered, the counter is closed on reference-derived k-mers. Families keep exon sharing + SD, and copy number
  stays a DNA quantity.
- Not tried: k-mers taken from the individual (A119b's own reads or assembly) instead of CHM13.
