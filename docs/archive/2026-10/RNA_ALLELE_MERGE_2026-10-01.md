# Merging a missing copy's new-copy transcripts into one candidate copy, without truth (Amendment 8), 2026-10-01

Prereg: `docs/archive/2026-10/PREREG_rna_allele_haplotype_count_2026-10-01.md` Amendment 8 (commit b44934d6, before any run). Script
`bench/rna_allele/merge_test.py`; work dir `/mnt/linuxdisk/tmp/rna_allele/linktest/merge/` (`components.out`, `score.out`, copies in
`docs/RNA_ALLELE_MERGE_score.out.txt`). Same 53 families, same reads and same R+I+L alignments as Amendment 7
(`docs/archive/2026-10/RNA_ALLELE_LINKING_2026-10-01.md`): nothing was realigned, only which contigs count as one locus changed.

- **Merge rule:** within a family, the 565 new-copy contigs (46 families) aligned all-vs-all (`minimap2 -c -x asm20 --cs -N 200 -p 0.1`,
  9,805 alignments); two contigs joined when the best alignment of the pair covers >= 50% of the shorter contig and its gap-compressed
  divergence `de` <= delta = 0.00958 (Amendment 7's allele-divergence cut, unchanged); components = candidate copies.
- **Result of the merge:** 565 contigs -> **75 candidate copies**: 53 hold only deleted-copy transcripts, 21 hold only survivor
  transcripts (spurious new copies, in 14 families, at most 4 per family), 1 mixed (GWFAM28). 44 families have >= 1 deleted-copy
  candidate: **36 have exactly one, 6 have two, 2 have three** (over-split).

## Registered result

| | R (masked) | R+I+L (Amendment 7) | **R+I+L+M (merged)** | T (grouped WITH truth) |
|---|---|---|---|---|
| D right | 0 | 1,616 | **12,787 (74.0%)** | 12,797 |
| D wrong | 10,218 | 3,584 | 3,599 | 3,592 |
| D unplaced | 7,068 | 12,086 | 900 | 897 |
| S stay | 41,051 | 38,595 | 40,421 | 39,883 |
| S unplaced | 670 | 3,115 | 1,278 | 1,816 |
| S false moves | 0 | 14 | **25 (0.06%)** | 25 |

- **M1 (merging works): PASSES** — D right 1,616 -> 12,787 (bar >= 6,440), false moves 0.06% (bar <= 5%). The truth-free merge reaches
  the truth-grouped ceiling (12,787 vs 12,797; the ceiling recomputed here groups D- and S-derived contigs by source and is 82 reads
  below Amendment 7's post hoc figure of 12,879, which grouped every source).
- **Overall (Amendment 5's rule, R vs R+I+L+M): HELP** — wrong D 10,218 -> 3,599 (-64.8%), false moves 0.06%.
- **M2 (one missing copy, one candidate): HOLDS** — the deleted copy's transcripts fall in ONE component in 33/41 families with >= 2 such
  transcripts (80.5%, bar 2/3); mixed components 1/75 (1.3%, bar 10%).
- **Sensitivity (reported, not the verdict):** at delta/2 and 2 x delta the read-level numbers are the same to the read (D right 12,787
  in all three; false moves 25 / 25 / 29); components 76 / 75 / 71; one-component families 33 / 33 / 35 of 41. The merge is not a
  tuned threshold: within a factor of 4 around delta nothing moves.

## Why the 8 over-split families stay split

Pairs of deleted-copy contigs that ended in different components: 62 pairs have an alignment covering < 50% of the shorter contig and 41
have no alignment at all (non-overlapping fragments of one gene with no full-length transcript to bridge them), 60 pairs align but diverge
beyond delta (two transcripts of the same copy more than 0.96% apart: errors, alleles beyond the 99th percentile, or chimeric outputs).
The split costs almost nothing at the read level (12,787 vs 12,797 right) because the reads of each fragment tie only within their own
component; it costs at the count level: 10 extra candidate copies over 44 families.

## Reading

- The full chain is now truth-free end to end: family-scoped IsoCon -> flag outputs not in the reference -> link outputs within allele
  divergence to their reference locus -> merge the rest into candidate copies -> O2 with the candidates as loci. On held-out multi-copy
  families it gives a missing copy's reads a home 74% of the time (9% before the merge, 0% without the chain), moves 0.06% of the
  surviving copies' reads, and names one candidate per missing copy in 80% of the families where it names any.
- **O3 reading:** a candidate copy is the RNA-only flag; the count is a lower bound with a known inflation (over-split fragments: 54 deleted-copy candidates for 44 deleted copies)
  and a known false-flag source (21 survivor-derived candidates in 14 families: alleles beyond delta or variants). The false-flag rate
  WITHOUT any deletion is not measured here (Amendment 8 lists it as the separate control).
- **O2 reading:** unchanged from Amendment 7 — the gain is through the copy set the aligner sees, not through the assignment rule.
- Caveats as before: one individual (the reference animal), fibroblast expression, IsoCon capped at 1,000 reads per family, 9 of the 53
  deleted copies produced no new-copy transcript at all (their outputs linked to a survivor within delta, or IsoCon returned none).
