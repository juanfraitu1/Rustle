# PREREG r1612: tie-group loci cut by aligned-block overlap, not span overlap (2026-11-10)

Audit finding (2026-11-10, after r1611): `tie_groups.py build` joined a record to a locus when the two SPANS overlapped by half
of the shorter. A record lying inside a long seed's intron always passed: on the A119b panel records 489 of 5,541 loci had a seed
more than 5x the median member span (largest: chr1:146.15-146.28 Mb, a 128 kb seed holding 10,592 records of median 9.5 kb). Such
a locus gets one model, so "one model per locus" is not true there. Written before any re-scoring.

## Change (one rule, nothing else)
A record joins the current locus when its aligned blocks overlap the seed's aligned blocks by >= 50% of the smaller of the two
block totals. Recipe otherwise frozen as r1607: `build PREFIX as98 0 1 RMSK_BED:0.5`, stranded, scorer `full`, SEDEF co-duplicated.
Self-check: `bench/famcmp/test_tie_groups.py` (fails on the span rule).

## Runs
A119b `r1603/recs.slim.tsv` and gorilla OR6737 `r1611/recs.slim.tsv`, same scoring inputs as r1607 and r1611.

## Bars (those of r1607 and r1611, unchanged, both samples)
B1g grouped reach (bar printed by the scorer), B2s purity with co-duplicated groups >= 0.90, B3 TAF11L models <= 5 kb >= 20 (A119b),
B5r largest group on one chromosome or >= 0.90 one family, B4 <= 1 h.
- PASS on both: the block rule is the recipe; r1607 and r1611 figures are replaced by these.
- FAIL on either: reported as is; the span rule's figures stay on record with the audit caveat, and the locus rule is the open item.
Descriptive: loci, groups, loci set aside as repeat, strict purity.

## Outcome (2026-11-10, `r1612/`, job 4456476; human 1,406 s, gorilla 222 s)
| | A119b span rule (r1607) | A119b block rule | gorilla span (r1611) | gorilla block |
|---|---|---|---|---|
| loci / groups | 2,509 / 1,141 | 20,008 / 9,766 | n/a | 1,579 / 952 |
| B1g grouped reach | 70/72 | 69/72 | 124/128 | 124/128 |
| B2s purity (co-duplicated) | 220/224 = 0.982 | 4,839/4,840 = 1.000 | 129/129 | 404/404 |
| B2f strict | 0.884 | 0.979 | 0.907 | 0.941 |
| B3 TAF11L <= 5 kb | 31/31 | 31/31 | - | - |
| B5r largest group | 86 loci | 1,582 loci, 1 chromosome, 1.000 one family | 33 | 58, 1 chromosome |

**Bars: PASS on both samples as registered. NOT adopted.** The bars do not see the new defect: the rule compares a record only
with the last locus's seed, so records of two interleaved shapes at one place open a new locus each time. 15,791 of 20,008 A119b
models (669 of 1,579 gorilla) overlap the previous model on the same strand by half: the same place is modelled many times, the
group counts are inflated and the purity figures are not comparable with r1607. The span rule under-splits (489 loci), the block
rule against one seed over-splits. Open: compare a record with every open locus (or cluster the blocks), with a registered bar on
duplicate models per place.
