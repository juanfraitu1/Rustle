# The no-deletion control of the IsoCon chain (Amendment 9), 2026-10-01

Prereg: `docs/PREREG_rna_allele_haplotype_count_2026-10-01.md` Amendment 9 (commit 41dace23, before any run). Script
`bench/rna_allele/control_test.py`; work dir `/mnt/linuxdisk/tmp/rna_allele/control/` (`contigs.out`, `classify.out`, `score.out`,
`candidates.tsv`; copies in `docs/RNA_ALLELE_CONTROL_score.out.txt`). Same 53 families and 59,013 reads as Amendments 7-8, NOTHING
masked: the chain ran against the full `_pri`, so every candidate copy it names is a flag raised without a deletion.

- **Chain:** R0 alignment -> IsoCon inputs 44,746 reads (53 families, median 1,000) -> 2,441 outputs -> 1,152 flagged (not in `_pri` at
  0.999) -> 1,076 linked back (d <= delta: alleles / variants of their own locus) -> **76 new-copy contigs in 16 families -> 28 candidate
  copies** (delta/2: 30, 2 x delta: 26).
- **Diploid adjudication (KB3781's own mat / pat assemblies):** 24 candidates match NEITHER haplotype at 0.999 (best haplotype hit
  0.90-0.999, mostly their own locus), 3 are alleles beyond delta (GWFAM247, GWFAM28, GWFAM331: inside the lifted B interval of a copy of
  the family), and **1 is a haplotype-only locus: GWFAM175, 12 transcripts, paternal chr5 `CM054563.2:40,028,172-40,031,100` at identity
  1.0000**, outside the lift of all six maternal (= `_pri`) copies of this tandem array — a genuine reference-absent expressed copy, found
  with no deletion.
- 20 of the 28 candidates (18-21 across the three deltas) reproduce the deletion run's survivor-derived candidates (same family, identity
  >= 0.999): the false flags are locus-specific and reproducible, not sampling noise.

## Registered result

| | R0 (`_pri`) | C (`_pri` + candidates) |
|---|---|---|
| stay on own copy | 58,097 | 55,306 |
| on a candidate derived from own copy | – | 2,018 (3.4%) |
| **false move** (candidate of another copy) | – | **27 (0.05%)** |
| other copy of the family | 13 | 11 |
| unplaced | 903 | 1,651 |

- **C1 (specificity): FAILS** — families with >= 1 false candidate (classes b + c) **16/53 = 30.2%**, bar <= 28% (<= 14 families). The
  count against the haploid reference alone is the same 16/53 (GWFAM175's true flag shares its family with four false ones). Family-level
  likelihood ratio of a flag: 0.830 / 0.302 = **2.75** (bar 3).
- **C2 (cost without a deletion): HOLDS** — false moves 27/59,013 = 0.05% (bar 5%); 2,018 reads move onto candidates derived from their
  own copy (harmless for assignment; they are the reads behind the flags); unplaced +748 (+1.3%).

## Post hoc (not the verdict): candidate support

| candidate support | deletion run: families detected | control: families with a false candidate | LR |
|---|---|---|---|
| any (>= 1 transcript) | 44/53 = 0.830 | 16/53 = 0.302 | 2.75 |
| **>= 2 transcripts** | 41/53 = 0.774 | **5/53 = 0.094** | **8.2** |

- In the deletion run 45 of the 54 deleted-copy candidates hold >= 2 transcripts; in the control 22 of the 27 false candidates are single
  transcripts (IsoCon read support per transcript: min 2, median 4, max 25 — a single transcript is typically 3-4 reads). A floor of two
  transcripts per candidate is the analogue of the assembler's floor of 2 and is the rule to PRE-REGISTER for the next tests (real
  reference-absent copies, YAGs); it is not applied here after the fact.

## Reading

- As registered the flag misses the LR-3 bar by two families (16 vs 14): a 2.75-fold enrichment at a 0.05% read cost is real, but it is
  not the claim we wanted at any support.
- What the false flags are: 24/27 are transcripts that exist in neither haplotype at consensus accuracy while matching their own locus at
  0.93-0.999 — single IsoCon outputs 1-7% away from the genome (few-read consensus, chimeras, one 150 kb spliced hit); only 3 are alleles
  beyond the 99th percentile, so the delta link is doing its job.
- The control produced the first reference-absent copy found without a deletion (GWFAM175, paternal chr5, 12 transcripts): the test case
  for item 2 (real reference-absent copies, truth = the haplotypes).
- Caveats: one individual (the reference animal), fibroblast expression, 1,000-read cap per family; the allele / haplotype-only split
  rests on the asm5 lift (200/201 copies lifted; GWFAM175:2 did not lift, consistent with a copy-number difference between the haplotypes
  in that array).
