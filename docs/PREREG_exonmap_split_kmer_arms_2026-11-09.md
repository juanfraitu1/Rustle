# PREREG rows 1597-1598: does the exon map split false families, and does a looser k-mer counter reach more reads?

Written 2026-11-09, after rows 1595 (exon map PASS) and 1596 (k-mer counter FAIL: agreement 0.985, reach 0.43, Spearman
0.67), before any split family, arm or diagnostic is computed. User: "lets try the exon map and counter and see if any of
them make a difference".

## Row 1597: split families by the exon map (pairscore)

Input: row 1595's columns (`bench/denovo_eval/exonmap.py`). In a family, two members are **linked** when they share
columns (both hold an exon in the same non-private column). Each family is split into the connected components of its
linked members, and every component becomes a family (`FAM.sN`; one-member components give no pairs).

- **S1**: linked iff >= 1 shared column.
- **S50**: linked iff shared columns >= 0.5 x the column count of the member with fewer distinct columns.

Scored with `pairscore.py` (defaults: simiiformes, `PAIRSCORE_SAME_SYMBOL=ignore`), the same command as rows 1592-1593.

- **Dev (A119b, `port1593/Ab`)**: baseline = Ab unsplit. Choose the arm with the larger precision gain among arms with
  recall >= baseline - 0.01. Go to held-out only if that gain is >= +0.01.
- **Held-out (HG002, `port1593/Hb`, exon map built the same way)**: **PASS** iff precision >= Hb + 0.01 AND recall >= Hb -
  0.01 (Hb rescored in the same job).
- What follows: a PASS makes the split an opt-in families-stage step; defaults stay the user's call. FAIL or STOP: the
  exon map stays a description, not a family filter.

## Row 1598: k-mer counter arms on a chromosome split (row 1596's definitions otherwise unchanged)

- **Arms**: k in {21, 25, 30} x minimum hits m in {1, 2, 3}. SPEC at each k is defined as in row 1596 (no other copy,
  genome count <= 1). Resolved = |SPEC| >= 100 at that k.
- **Split** by the chromosome of each family's `copy_idx` 0:
  - dev = odd-numbered autosomes + chrX;
  - held-out = even-numbered autosomes, chrY and any other contig.
- **Dev**: choose the arm with the highest reach among arms with agreement >= 0.98 (ties: larger k, then larger m).
- **Held-out bars** (row 1596's three, unchanged): agreement >= 0.98, reach >= 0.70, Spearman >= 0.90 (over all
  in-scope reads, as registered). STOP as row 1596 (resolved copies < 200 or assigned unique reads < 1,000 in held-out).
- **Diagnostic** (descriptive, both halves, every arm): unique-placement reads on resolved copies that stay unassigned,
  split into three classes:
  - none: no SPEC hit on the family;
  - weak: hits only on the primary's copy, fewer than m;
  - conflict: a hit on another copy of the family.
- What follows: a PASS makes the counter usable at that arm, and the KB3781 tied-read test (row 1596's follow-up)
  becomes the next row. A FAIL with conflict dominant points to individual alleles (A119b is not CHM13); a FAIL with
  none or weak dominant points to too few SPEC k-mers per read. Either closes the counter on reference-derived k-mers.

Outputs: `/anvme/workspace/v138da13-rustle/r1597/`, `r1598/`.

## Outcome row 1597 (2026-11-09, job 4453988): STOP, no held-out run

The exon map split 48 (S1) and 50 (S50) of A119b's 1,395 families. The HG002 exon map was also built: collision 0.0025,
closure 0.979.

| A119b dev | pairs | hit | false | precision | recall |
|---|---|---|---|---|---|
| Ab unsplit | 1,860 | 1,707 | 153 | .918 | .782 |
| S1 | 1,856 | 1,703 | 153 | .918 | .780 |
| S50 | 1,850 | 1,702 | 148 | .920 | .780 |

Best gain +0.002 (S50) < +0.01: **STOP** as registered. The split removes 4 or 5 true pairs and at most 5 false ones.
Reading: the shipped exon-sharing gate already keeps only exon-linked members, so the family-level map has almost nothing
left to cut. The exon map stays a description, not a filter.

## Outcome row 1598 (2026-11-09, job 4453989, 2,205 s): FAIL

SPEC k-mers / resolved copies:

| k | SPEC k-mers | resolved copies |
|---|---|---|
| 21 | 1,931,635 | 2,179 |
| 25 | 2,288,885 | 2,375 |
| 30 | 2,671,098 | 2,571 |

All nine arms on dev:
- reach 0.423-0.433;
- agreement 0.978-0.986 (k21 m1 and k25 m1 fall below 0.98);
- Spearman 0.64-0.67.

Dev choice k = 21, m = 2 (reach 0.4310). Held-out: agreement 0.9860, reach **0.4282** (bar 0.70), Spearman **0.7117**
(bar 0.90): **FAIL** as registered. Full table: `r1598/ra.arms.tsv`.

Diagnostic (unassigned unique-placement reads on resolved copies, chosen arm):

| half | none | weak | conflict |
|---|---|---|---|
| dev | 222,346 | 461 | 73,980 |
| held | 136,261 | 467 | 51,503 |

So about 75% none, 25% conflict, weak about 0.

Reading:
- Neither k nor m moves reach: every arm sits at 0.42-0.43. The counter is not limited by k-mer length or hit count.
- Most unassigned reads carry no SPEC k-mer of the family at all. Post hoc, not tested: unique-placement was defined on
  the copy's span (`locus_start..locus_end`), so it includes reads in introns, flanks and nested genes that have no copy
  exon sequence to match.
- The conflict quarter is consistent with A119b alleles carrying another copy's reference k-mer.
- As registered, this closes the counter on reference-derived k-mers. The KB3781 tied-read test is not run.
