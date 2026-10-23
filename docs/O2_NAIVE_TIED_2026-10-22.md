EXPLORATORY (post hoc, not pre-registered)
# Naive tied-read policies (NAIVE-BEST, NAIVE-ALL) against minimap2's primary: read level and copy level, 2026-10-22

Question (user): as a naive O2, impute that the right copy is "the one within 10 AS", and test whether that improves on the primary placement and on loci determination. Read-only analysis of the frozen audit `/mnt/linuxdisk/tmp/o3_graph/o2vc/audit_frozen_eccbc4e0/{hg002,hg00733}` (calls.tsv, `<family>.reads.tsv`, `<family>.case.tsv`, score.C.json, score.E.json); no decision rule was fixed beforehand, so every number is descriptive. Scripts (not committed): `/tmp/claude-1000/-mnt-c-Users-jfris-Desktop/e542c6ba-28cf-42e7-becf-fd99d54f830f/scratchpad/naive_tied.py` (analysis) and `.../build_md.py` (this layout).

## Definitions and checks
- Tied read: alignments within 10 AS of the best on >= 2 copies of one family. On every read the AS spread across its tied copies is <= 10 (checked: 0 violations).
- **NB-abst**: place on the tied copy with the strictly highest AS, abstain on an exact tie at the max (scored reads: 622 of 6333 in HG002, 5132 of 12015 in HG00733). **NB-prim**: same, but at an exact tie keep minimap2's primary (N1, even when it is not one of the max-tied copies: 156 / 25 tied reads, HG002 / HG00733).
- **NAIVE-ALL**: the read is credited to every tied copy (`candidates`). Not a single placement, so it has no wrong/placed; it is scored by set containment and, at copy level, by credits.
- Scoring as in `cmd_score`: scored = tied reads of class unique, unique_elsewhere, elsewhere, no_counterpart; right iff the placement is the truth copy (elsewhere / no_counterpart have none, so any placement is wrong); fixes = right where N1 is wrong, breaks = wrong where N1 is right; PAR families (par_families.json) excluded (15 HG002 families with tied reads, none in HG00733). Of the 20100 non-PAR tied HG002 reads 6333 are scored (not scored, no definite truth: own_tie 9642, poor_match 3800, no truth row 325); all 12015 HG00733 tied reads are.
- Checks: my N1, N2, N4, VC placed/right/wrong/abstained/fixes/breaks equal score.C.json and score.E.json exactly (asserted); NB-abst / NB-prim / N1 were re-derived a second way from the `case.tsv` truth columns and agree to the read.

## 1. Read level

**HG002**, left block all four scored classes (N = 6333), right block unique + unique_elsewhere only, i.e. without elsewhere / no_counterpart (N = 6083); abstained = N - placed.

| policy | placed | right | wrong | wrong/placed | fixes | breaks | placed | right | wrong | wrong/placed | fixes | breaks |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| N1 (minimap2 primary) | 6333 | 3875 | 2458 | 0.388 | - | - | 6083 | 3875 | 2208 | 0.363 | - | - |
| N2 (MAPQ >= 1) | 5725 | 3647 | 2078 | 0.363 | 0 | 0 | 5507 | 3647 | 1860 | 0.338 | 0 | 0 |
| VC | 4579 | 4308 | 271 | 0.059 | 1137 | 24 | 4393 | 4308 | 85 | 0.019 | 1137 | 24 |
| N4 (corrector) | 6333 | 5096 | 1237 | 0.195 | 1252 | 31 | 6083 | 5096 | 987 | 0.162 | 1252 | 31 |
| **NB-abst** | 5711 | 3645 | 2066 | 0.362 | 36 | 38 | 5505 | 3645 | 1860 | 0.338 | 36 | 38 |
| **NB-prim** | 6333 | 3873 | 2460 | 0.388 | 36 | 38 | 6083 | 3873 | 2210 | 0.363 | 36 | 38 |

**HG00733** (N = 12015, every read class unique):

| policy | placed | right | wrong | wrong/placed | fixes vs N1 | breaks vs N1 |
|---|---|---|---|---|---|---|
| N1 (minimap2 primary) | 12015 | 7420 | 4595 | 0.382 | - | - |
| N2 (MAPQ >= 1) | 6953 | 5223 | 1730 | 0.249 | 0 | 0 |
| VC | 2801 | 2625 | 176 | 0.063 | 279 | 30 |
| N4 (corrector) | 12015 | 7788 | 4227 | 0.352 | 408 | 40 |
| **NB-abst** | 6883 | 5410 | 1473 | 0.214 | 342 | 112 |
| **NB-prim** | 12015 | 7650 | 4365 | 0.363 | 342 | 112 |

wrong/placed is not comparable across arms that abstain (N2, VC, NB-abst place fewer reads): compare NB-abst with N2, which abstains on the same exact ties, and the others by `right` out of N.

**NAIVE-ALL** (the set of tied copies as the answer):

| dataset | scored reads | reads with a truth copy | set contains truth | fraction of scored | fraction of truth-copy reads | mean set size | credits on the truth copy / all credits | N1 right / scored |
|---|---|---|---|---|---|---|---|---|
| HG002, four classes | 6333 | 6083 | 5892 | 0.930 | 0.969 | 2.06 | 0.470 | 0.612 |
| HG002, unique + unique_elsewhere | 6083 | 6083 | 5892 | 0.969 | 0.969 | 2.06 | 0.470 | 0.637 |
| HG00733 | 12015 | 12015 | 11607 | 0.966 | 0.966 | 2.37 | 0.407 | 0.618 |

Truth outside the tied copies (never containable): 191 reads in HG002, 408 in HG00733. Set sizes in HG002: {2: 5973, 3: 337, 4: 15, 5: 5, 6: 1, 7: 2}; HG00733: {2: 9566, 3: 1270, 4: 850, 5: 110, 6: 49, 7: 71, 8: 42, 9: 56, 10: 1}.

## 2. Copy level (loci determination)
Credits per family over the scored tied reads. TRUTH = truth copy of unique / unique_elsewhere reads; N1 = primary; VC = placed reads only; VC+N1 = VC, else N1; NB-prim; NAIVE-ALL = +1 to every tied copy. Expressed = credited >= 2 reads; sens / prec pooled over (family, copy) against the copies with >= 2 TRUTH reads. A primary on no copy counts as a credit (false expression, TVD mass) but not as a called copy. TVD = 0.5 * sum |policy fraction - TRUTH fraction| per family, over families with a TRUTH read; weighted by the family's TRUTH reads; VC (and NB-abst) are undefined where they place nothing, so the family count is shown.

**HG002, tied reads only** (TRUTH: 6083 reads on 117 copies of 61 families, 87 expressed copies); left block four classes (6333 reads, 250 credited but with no truth copy), right block unique + unique_elsewhere (6083 reads):

| policy | called | sens | prec | zero-truth credits | TVD med | TVD wtd | called | sens | prec | zero-truth credits | TVD med | TVD wtd |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| N1 | 105 | 0.874 | 0.724 | 743/6333 = 0.117 | 0.187 | 0.347 | 97 | 0.874 | 0.784 | 700/6083 = 0.115 | 0.187 | 0.347 |
| VC | 48 | 0.471 | 0.854 | 17/4579 = 0.004 | 0.083 | 0.126 | 45 | 0.471 | 0.911 | 7/4393 = 0.002 | 0.049 | 0.129 |
| VC+N1 | 105 | 0.885 | 0.733 | 728/6333 = 0.115 | 0.100 | 0.166 | 97 | 0.885 | 0.794 | 683/6083 = 0.112 | 0.060 | 0.165 |
| NB-abst | 99 | 0.839 | 0.737 | 530/5711 = 0.093 | 0.153 | 0.338 | 91 | 0.839 | 0.802 | 514/5505 = 0.093 | 0.132 | 0.338 |
| NB-prim | 110 | 0.897 | 0.709 | 731/6333 = 0.115 | 0.180 | 0.345 | 102 | 0.897 | 0.765 | 700/6083 = 0.115 | 0.164 | 0.344 |
| NAIVE-ALL | 147 | 0.954 | 0.565 | 1658/13062 = 0.127 | 0.500 | 0.469 | 135 | 0.954 | 0.615 | 1598/12543 = 0.127 | 0.500 | 0.468 |

TVD on the 35 four-class families where every policy has credits (median / read-wtd mean): N1 0.182/0.327, VC 0.083/0.126, VC+N1 0.079/0.130, NB-abst 0.125/0.315, NB-prim 0.164/0.324, NAIVE-ALL 0.499/0.464. VC places a scored read in 35 of the 61 families with a truth read, NB-abst in 59.

**HG00733, tied reads only** (12015 reads; TRUTH: 12015 reads on 352 copies of 174 families, 294 expressed copies):

| policy | copies called (>= 2 reads) | sens | prec | credits on zero-truth copies | TVD median (families) | TVD read-wtd mean |
|---|---|---|---|---|---|---|
| N1 | 330 | 0.939 | 0.836 | 762/12015 = 0.063 | 0.136 (174) | 0.248 |
| VC | 158 | 0.531 | 0.987 | 71/2801 = 0.025 | 0.064 (116) | 0.304 |
| VC+N1 | 331 | 0.956 | 0.849 | 729/12015 = 0.061 | 0.093 (174) | 0.236 |
| NB-abst | 272 | 0.833 | 0.901 | 398/6883 = 0.058 | 0.086 (168) | 0.226 |
| NB-prim | 329 | 0.932 | 0.833 | 772/12015 = 0.064 | 0.143 (174) | 0.227 |
| NAIVE-ALL | 399 | 0.997 | 0.734 | 2741/28533 = 0.096 | 0.373 (174) | 0.358 |

TVD on the 114 HG00733 families where every policy has credits (median / read-wtd mean): N1 0.138/0.254, VC 0.069/0.304, VC+N1 0.090/0.234, NB-abst 0.085/0.173, NB-prim 0.143/0.213, NAIVE-ALL 0.327/0.327.

**Expressed copies that N1 misses versus copies NAIVE-ALL adds** ("beyond N1" = called by NAIVE-ALL, not by N1; "expressed" = >= 2 TRUTH reads):

| analysis | TRUTH-expressed copies | N1 misses | NAIVE-ALL finds | VC+N1 finds | NAIVE-ALL calls beyond N1 | of which expressed | with 1 truth read | with 0 truth reads |
|---|---|---|---|---|---|---|---|---|
| HG002 tied, four classes | 87 | 11 | 7 | 1 | 44 | 7 | 2 | 35 |
| HG002 tied, unique + unique_elsewhere | 87 | 11 | 7 | 1 | 39 | 7 | 4 | 28 |
| HG00733 tied | 294 | 18 | 17 | 5 | 69 | 17 | 17 | 35 |
| HG002 tied + untied, unique + unique_elsewhere | 407 | 22 | 4 | 0 | 15 | 4 | 1 | 10 |
| HG00733 tied + untied | 656 | 2 | 1 | 0 | 14 | 1 | 0 | 13 |

**Tied + untied reads** (untied = 205524 HG002 reads of unique / unique_elsewhere in truth_untied.tsv, placed on minimap2's primary, which is right on 98.9% of them; HG00733: 71685 reads, placed on the truth copy by assumption). TRUTH HG002: 211607 reads on 486 copies of 222 families, 407 expressed; HG00733: 83700 reads on 657 copies of 274 families, 656 expressed.

| policy | HG002 sens | prec | zero-truth credits | TVD wtd | HG00733 sens | prec | zero-truth credits | TVD wtd |
|---|---|---|---|---|---|---|---|---|
| N1 | 0.946 | 0.941 | 1927/211607 = 0.009 | 0.020 | 0.997 | 0.965 | 717/83700 = 0.009 | 0.036 |
| VC | 0.919 | 0.964 | 1267/209917 = 0.006 | 0.015 | 0.980 | 0.998 | 70/74486 = 0.001 | 0.039 |
| VC+N1 | 0.946 | 0.944 | 1915/211607 = 0.009 | 0.014 | 0.997 | 0.965 | 699/83700 = 0.008 | 0.034 |
| NB-abst | 0.941 | 0.950 | 1752/211029 = 0.008 | 0.020 | 0.994 | 0.979 | 381/78568 = 0.005 | 0.034 |
| NB-prim | 0.951 | 0.942 | 1921/211607 = 0.009 | 0.019 | 0.997 | 0.965 | 725/83700 = 0.009 | 0.033 |
| NAIVE-ALL | 0.956 | 0.920 | 2514/218067 = 0.012 | 0.021 | 0.998 | 0.947 | 2522/100218 = 0.025 | 0.072 |

With the 260772 untied no_counterpart / elsewhere reads kept (HG002 four classes), every policy puts ~55% of its credits on zero-truth copies (labels partly liftover artifacts) and prec differs by only 0.034 between N1 and NAIVE-ALL, so that variant is dominated by those reads.

## Caveats
- HG002 `elsewhere` / `no_counterpart` labels partly come from liftover-chain artifacts; hence both variants above. HG00733 is a simulation from the donor haplotypes, aligned to CHM13, not a real-read result.
- The >= 2-reads rule is lenient for NAIVE-ALL: two reads tied on the same pair make both copies "expressed" by construction. Tied-only copy sets are small (HG002 61 families with a truth copy, HG00733 174).
- Untied reads dilute every difference (the policies agree on them); HG00733 untied reads are credited to the truth copy because the frozen outputs hold no minimap2 placement for them, an upper bound for N1.

## Reading
- NAIVE-BEST gives no gain over the primary in HG002 (NB-prim wrong/placed 0.388 vs N1 0.388, 36 fixes / 38 breaks; NB-abst 0.362 vs N2 0.363) and a modest one in the simulation (0.363 vs 0.382, 342 / 112; NB-abst 0.214 vs N2 0.249), far from VC (0.059 / 0.063): the highest tied AS carries little beyond the primary's own pick.
- NAIVE-ALL does recover expressed copies the primaries miss (HG002 7 of 11, HG00733 17 of 18; VC+N1 only 1 and 5), its set holds the truth for 96.9% / 96.6% of reads, at 2.06 / 2.37 copies per read.
- It also invents expression: 39 / 69 extra copies called beyond N1, only 7 / 17 expressed (18% / 25%) and 28 / 35 with zero truth reads (~4 / ~2 per recovered copy); prec 0.615 vs N1 0.784 and 0.734 vs 0.836, TVD median 0.500 / 0.373 vs 0.187 / 0.136: a recall-only lever (HG002 / HG00733).
