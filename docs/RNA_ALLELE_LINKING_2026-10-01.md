# Linking IsoCon transcripts to their source locus, held-out on 53 multi-copy families (Amendment 7), 2026-10-01

Prereg: `docs/PREREG_rna_allele_haplotype_count_2026-10-01.md` Amendment 7 (commit e66d1c04, before any run). Script
`bench/rna_allele/link_test.py` (+ `iso_batch.sh`); work dir `/mnt/linuxdisk/tmp/rna_allele/linktest/` (`score.out`).

- **Linking rule:** an IsoCon output flagged "not in the reference" is an allele of its best masked-genome locus when its whole-length
  divergence d <= 0.00958 (the 99th percentile of exonic divergence between KB3781's two haplotypes over 28,541 single-copy genes); else a
  new copy.
- **Substrate (never used before):** 53 families with >= 3 copies (2 of 55 dropped: masked interval overlaps another locus), 201 copies;
  in each the copy last by position hard-masked (1.30 Mb); 17,286 reads of the deleted copies (D) and 41,727 of the surviving copies (S).
- IsoCon: 2,309 outputs; 1,372 flagged; the rule keeps **565 as new copies** (497 from the deleted copies, 68 from survivors) and links
  **807** back (750 from survivors, 57 from deleted copies that sit within 0.96% of a survivor).

## Registered result

| | R (masked) | R+I (all flagged, Amendment 6) | R+I+L (linked) |
|---|---|---|---|
| D right | 0 | 1,525 | 1,616 |
| D wrong | 10,218 | 3,439 | 3,584 |
| D unplaced | 7,068 | 12,322 | 12,086 |
| S stay | 41,051 | 15,521 | 38,595 |
| S unplaced | 670 | 26,191 | 3,115 |
| S false moves | 0 | 12 | 14 |

- **Linking WORKS:** S unplaced 26,191 -> 3,115 (-88%), D right 1,525 -> 1,616 (+6%).
- **Overall HELP (R vs R+I+L):** wrong D 10,218 -> 3,584 (-65%), false moves 14 / 41,727 = 0.03%.

## Post hoc (not the verdict): a copy's transcripts counted as one locus

IsoCon returns several transcripts per copy (isoforms, ends), so most remaining "unplaced" D reads tie between two transcripts of the same
deleted copy. Counting all contigs derived from one source copy as one locus (source from the unmasked genome, so this uses truth):

| | R | R+I+L |
|---|---|---|
| D right | 0 | **12,879 (74.5%)** |
| D wrong | 10,218 | 3,510 |
| D unplaced | 7,068 | 897 |
| S stay | 40,736 | 39,926 (95.7%) |
| S moved to another copy or its contigs | 15 | 273 (0.65%) |

## Reading

- The procedure (family-scoped IsoCon -> flag outputs not in the reference -> link outputs within allele divergence of a reference locus
  -> add the rest as new copies) passes both registered rules on held-out multi-copy families: it gives a deleted copy's reads a home and
  stops the ties the unlinked version created, at a 0.03% cost to the surviving copies.
- Missing step, untested: merging the several new-copy transcripts that come from one missing copy into one copy without truth (e.g.
  link new-copy outputs to each other at d <= delta). The post hoc grouping shows what it would buy: 9% -> 75% of D reads placed right.
- 57 transcripts of deleted copies are within 0.96% of a survivor and get linked to it: a missing copy closer than allele divergence to a
  paralog is indistinguishable from an allele, by construction.
- Caveats: one individual (the reference animal); fibroblast expression; IsoCon capped at 1,000 reads per family; S-derived new copies
  (68) are alleles beyond the 99th percentile or other variants.
