# Liftoff copy pairs against matched decoys: O1's independent-truth test (2026-10-14)

**Pre-registered before any decoy statistic is computed.**
- Design: agreed with the user in chat on 2026-10-14 ("close F1, then a new test": count every read-supported copy
  pair; use a null matched for genomic distance; human for development, gorilla held out). It implements the family
  bar of `docs/superpowers/specs/2026-10-13-o1-layered-hierarchy-design.md`, section 3.
- Code: `bench/liftoff_pairs/` at commit `97f831b4`, run with the defaults (1,000 null draws, seed 7):

      bash tools/rlock.sh light python3 bench/liftoff_pairs/run_pairs.py human_A119b --out /mnt/linuxdisk/tmp/o1_liftoff_pairs/human_A119b
      bash tools/rlock.sh light python3 bench/liftoff_pairs/run_pairs.py human_testis --out /mnt/linuxdisk/tmp/o1_liftoff_pairs/human_testis
      bash tools/rlock.sh light python3 bench/liftoff_pairs/run_pairs.py gorilla_KB3781 --out /mnt/linuxdisk/tmp/o1_liftoff_pairs/gorilla_KB3781

## Question

Do our default de novo families put a gene's Liftoff copies together, more often than they put the gene together with
comparable loci that are not its copies? Liftoff's copy relation is independent of our edge rule: it comes from
projecting the annotation onto the assembly with `-copies`. Human is for development; gorilla is held out.

## Inputs (frozen)

| What | Path |
|---|---|
| Liftoff loci (self-lift with `-copies -sc 0.95`, figure 8) | `/mnt/linuxdisk/tmp/rustle_figures/liftoff/{human,gorilla}/liftoff_loci.tsv` |
| Read support per sample (≥ 2 primary reads with an aligned block on the exon union, `-F 2308`) | `/mnt/linuxdisk/tmp/rustle_figures/liftoff/support/<sample>.tsv` |
| Our default families (shipped, exact-mcl; one copy per member locus) | `/mnt/linuxdisk/tmp/rustle_figures/runs/<sample>/<sample>.fam.copies.tsv` |
| Figure 8's F1 rows (the wiring check) | `figures/data/fig8_samples.tsv` |

## Known before running (recorded now)

- **F1** (`docs/archive/2026-09/PREREG_liftoff_loci_2026-09-25.md`; verdict recorded today, register 1283) asks whether
  copy pairs whose two loci are both covered by family members share a family.
  - Human A119b: 131/148 (0.885), failing the 0.90 bar by three pairs.
  - Human testis 34/34; gorilla KB3781 11/12; chimp 18/18.
  - So among pairs we already place in families, most share one. F1 has no null and conditions on our own
    catalog.
- **The same denominator was defined but never scored.** Supplementary figure 6s defined the unconditioned
  sensitivity (every read-supported pair; `docs/archive/2026-09/PREREG_genome_wide_families_2026-09-25.md`,
  Amendment 1, item 3). It was never scored: no table of it exists, and no number of it has been seen.
- **Why not a family-size null.** The layered spec named "a size-matched null". A null that only keeps family sizes
  is passed by design: plan 1's E1 (register 1277) and the Soto study's T1a (register 1281) both passed this way. This
  test therefore matches decoys on coverage and on genomic distance.
- **Development contigs.** The family rule was developed on human chr16 and gorilla chr20 (NC_073244.2). Results
  without those contigs are reported.
- **Nothing is tuned.** The distance bins, the widening rule and the decoy definition are fixed in code commit
  97f831b4. They were written before any decoy statistic existed.
- **Input counts** (`inputs.json`; inputs, not outcomes). Read-supported copy pairs at identity ≥ 0.95 / 0.98 / 0.99 / 1.00:

  | sample | genome | without development contigs | decoy pool (genome / without) |
  |---|---|---|---|
  | human A119b | 514 / 338 / 246 / 53 | 469 / 307 / 227 / 47 | 43,703 / 42,018 |
  | human testis | 38 / 31 / 24 / 3 | 33 / 26 / 19 / 1 | 19,454 / 18,606 |
  | gorilla KB3781 | 78 / 49 / 35 / 5 | 75 / 46 / 33 / 5 | 18,551 / 17,380 |

## Definitions

- **Copy pair (a, c).** a is a source gene's annotated placement (in place or moved) and c is one of its Liftoff extra
  copies.
  - The extra copy has `sequence_ID` ≥ the floor.
  - Both exon unions are ≥ 200 bp.
  - Both loci are read-supported in the sample (`_liftoff.copy_pairs` with support).
- **Covered.** A locus is covered by a family when a member of that family (its representative's exons) covers ≥ 50%
  of the locus's exon bases. This is figure 8's rule (`_liftoff.pair_families`).
- **Shared.** A pair is shared when one family covers both loci.
  - A locus covered by no family member cannot be shared.
  - There is no conditioning on our catalog: every read-supported pair is in the denominator.
- **Decoy.** A decoy for c is a read-supported Liftoff locus of another source gene that does not overlap a, matched to
  c on two things:
  1. whether a family member covers it;
  2. its distance bin from a: another contig, or the same contig at < 100 kb, 100 kb-1 Mb, 1-10 Mb, or ≥ 10 Mb
     (gap between spans).
- **Widening.** When no decoy matches:
  - level 1 keeps the contig relation at any distance;
  - level 2 allows any contig;
  - a pair without a decoy at level 2 leaves the test, on both sides.

  The number of pairs at each level is reported.

## The test

- **Statistic.** k = the number of shared pairs among the n pairs in the test.
- **Null.** Each draw replaces every pair's copy by a uniformly drawn decoy and counts the shared decoy pairs. There
  are 1,000 draws, seed 7, and p = (1 + draws ≥ k) / 1,001.
- **Bar:** p < 0.001 on the genome at identity ≥ 0.95, in each of human A119b, human testis and gorilla KB3781.
- **Reported, no bar:**
  - the identity floors 0.98, 0.99 and 1.00;
  - the genome without the development contigs;
  - k / n with a 95% interval (2,000 bootstrap resamples of source genes, figure 8's seed);
  - the null's mean rate and maximum;
  - the widening levels;
  - the pairs' coverage: both loci covered, one, or none.

## Readings (fixed now)

| human A119b | human testis | gorilla KB3781 | Reading |
|---|---|---|---|
| pass | pass | pass | O1's independent-truth claim: "a gene's Liftoff copies share one of our families far more often than distance- and coverage-matched loci of other genes (rate x against null y), in human and held out in gorilla". |
| pass | pass | fail | The same claim for human only; the gorilla result is reported. |
| fail on either human sample | — | — | No claim; the per-sample tables say where it fails. |

**What a pass does and does not show.**
- Decoys are matched on whether a family covers them, so the test asks whether our families group a gene with its own
  copies rather than with other covered loci at the same distance.
- It would fail if our families lumped neighbouring genes duplicated with the copy.
- A decoy can be a paralog annotated as a separate gene, which makes the null stricter, not looser.
- It is not a precision measure for loci without a Liftoff relation.

## Wiring check (the run stops, and nothing is read, if it fails)

F1 recomputed by the code must reproduce `figures/data/fig8_samples.tsv`'s n and k, for every scope and floor of the
sample. It did for all three samples in the input reports and in the tests.

## Results

Added after the runs. Nothing above this heading changes once a statistic has been seen.

### Results (2026-10-14)

Run: commit `2cd1446c` (code `97f831b4`). Wall time and peak memory: human A119b 4.56 s, 355 MB; human testis 2.17 s, 273 MB; gorilla KB3781 1.89 s, 233 MB.
Wiring check passed: F1 reproduced the figure 8 table for every scope and floor of each sample.

**The bar (genome, identity ≥ 0.95)**

| sample | pairs n | shared k | k / n [95% interval] | null mean rate (max k) | p | bar |
|---|---|---|---|---|---|---|
| human A119b | 514 | 75 | 0.1459 [0.1025, 0.1955] | 0.0214 (16) | 0.000999 (1/1001) | **PASS** |
| human testis | 38 | 12 | 0.3158 [0.1702, 0.5000] | 0.0321 (3) | 0.000999 (1/1001) | **PASS** |
| gorilla KB3781 | 78 | 4 | 0.0513 [0.0112, 0.1098] | 0.0006 (1) | 0.000999 (1/1001) | **PASS** |

**Coverage of the pairs and widening (genome, identity ≥ 0.95)**

| sample | both loci covered | one | none | shared / both covered | decoys at level 0 / 1 / 2 / none |
|---|---|---|---|---|---|
| human A119b | 91 | 133 | 290 | 75 / 91 | 507 / 7 / 0 / 0 |
| human testis | 12 | 2 | 24 | 12 / 12 | 35 / 3 / 0 / 0 |
| gorilla KB3781 | 4 | 8 | 66 | 4 / 4 | 78 / 0 / 0 / 0 |

**Every floor and scope (reported, no bar)**

| sample | scope | identity | n | k | k / n | null mean rate (max) | p |
|---|---|---|---|---|---|---|---|
| human A119b | genome | 0.95 | 514 | 75 | 0.1459 | 0.0214 (16) | 0.000999 (1/1001) |
| human A119b | genome | 0.98 | 338 | 62 | 0.1834 | 0.0284 (15) | 0.000999 (1/1001) |
| human A119b | genome | 0.99 | 246 | 51 | 0.2073 | 0.0383 (15) | 0.000999 (1/1001) |
| human A119b | genome | 1.00 | 53 | 18 | 0.3396 | 0.0855 (8) | 0.000999 (1/1001) |
| human A119b | minus_dev | 0.95 | 469 | 57 | 0.1215 | 0.0210 (13) | 0.000999 (1/1001) |
| human A119b | minus_dev | 0.98 | 307 | 48 | 0.1564 | 0.0291 (13) | 0.000999 (1/1001) |
| human A119b | minus_dev | 0.99 | 227 | 39 | 0.1718 | 0.0390 (14) | 0.000999 (1/1001) |
| human A119b | minus_dev | 1.00 | 47 | 15 | 0.3191 | 0.0942 (8) | 0.000999 (1/1001) |
| human testis | genome | 0.95 | 38 | 12 | 0.3158 | 0.0321 (3) | 0.000999 (1/1001) |
| human testis | genome | 0.98 | 31 | 10 | 0.3226 | 0.0392 (3) | 0.000999 (1/1001) |
| human testis | genome | 0.99 | 24 | 9 | 0.3750 | 0.0503 (3) | 0.000999 (1/1001) |
| human testis | genome | 1.00 | 3 | 1 | 0.3333 | 0.3333 (1) | 1.000000 (1001/1001) |
| human testis | minus_dev | 0.95 | 33 | 11 | 0.3333 | 0.0367 (4) | 0.000999 (1/1001) |
| human testis | minus_dev | 0.98 | 26 | 9 | 0.3462 | 0.0465 (3) | 0.000999 (1/1001) |
| human testis | minus_dev | 0.99 | 19 | 8 | 0.4211 | 0.0624 (3) | 0.000999 (1/1001) |
| human testis | minus_dev | 1.00 | 1 | 1 | 1.0000 | 1.0000 (1) | 1.000000 (1001/1001) |
| gorilla KB3781 | genome | 0.95 | 78 | 4 | 0.0513 | 0.0006 (1) | 0.000999 (1/1001) |
| gorilla KB3781 | genome | 0.98 | 49 | 3 | 0.0612 | 0.0000 (0) | 0.000999 (1/1001) |
| gorilla KB3781 | genome | 0.99 | 35 | 3 | 0.0857 | 0.0000 (0) | 0.000999 (1/1001) |
| gorilla KB3781 | genome | 1.00 | 5 | 0 | 0.0000 | 0.0000 (0) | 1.000000 (1001/1001) |
| gorilla KB3781 | minus_dev | 0.95 | 75 | 4 | 0.0533 | 0.0006 (1) | 0.000999 (1/1001) |
| gorilla KB3781 | minus_dev | 0.98 | 46 | 3 | 0.0652 | 0.0000 (0) | 0.000999 (1/1001) |
| gorilla KB3781 | minus_dev | 0.99 | 33 | 3 | 0.0909 | 0.0000 (0) | 0.000999 (1/1001) |
| gorilla KB3781 | minus_dev | 1.00 | 5 | 0 | 0.0000 | 0.0000 (0) | 1.000000 (1001/1001) |

### Verdict (2026-10-14)

**Reading 1.** The bar passes in human A119b, human testis and gorilla KB3781 (held out), each at p = 0.000999
(1/1001).
- **Claim:** a gene's Liftoff copies share one of our families far more often than distance- and coverage-matched loci
  of other genes.

  | sample | copy pairs shared | decoys (null mean) |
  |---|---|---|
  | human A119b | 75 of 514 (0.146) | 0.021 |
  | human testis | 12 of 38 (0.316) | 0.032 |
  | gorilla KB3781 | 4 of 78 (0.051) | 0.0006 |

  It holds in human, held out in gorilla, without the development contigs, and at the floors 0.98 and 0.99 in
  every sample. At 1.00, testis (3 pairs) and gorilla (5 pairs) do not pass.

**What this does and does not say.**
- **The test is effectively conditioned on both copies being covered** (post-hoc note below). Only those pairs can
  count, on either side: 91 (A119b), 12 (testis) and 4 (gorilla).
  - Among them, the copies share a family in 75 / 91 (0.824), 12 / 12 and 4 / 4.
  - The coverage-matched decoys give null rates of 0.121, 0.102 and 0.012.
  - A119b's 0.824 is below the 0.90 that F1 asks for, and its 16 split pairs are F1's.
- **What a pass shows.** A gene's covered copies are grouped with it much more often than covered loci of other genes at
  a similar distance.
  - It is a specificity check against unrelated loci, not proof against an independent family definition.
  - The copies are ≥ 0.95-identity sequence matches and our families are similarity clusters, so agreement is expected.
- **Recall is low, mostly because of coverage.** In most read-supported copy pairs neither copy is covered by any family
  member: 290 of 514 (A119b), 24 of 38 (testis), 66 of 78 (gorilla).
  - The other unshared pairs are grouping failures, F1's split pairs among them.
  - This is O1's known recall limit (node completeness and edge recall), and it is not measured as a bar here.
- **The gorilla pass rests on 4 pairs from 4 genes.** It passes because no random decoy draw reached 4 (null maximum 1).

**Post hoc (added after the final review): the denominator does not do what Definitions says.**
- **Uncovered pairs carry no information.** Decoys are matched on the copy's coverage, so a pair with either copy
  uncovered adds 0 to k and to every null draw.
  - k, the null and p are therefore exactly those of the same test run on the both-covered pairs.
  - "There is no conditioning on our catalog" (Definitions) holds for the denominator only.
  - The headline rates (0.146 against 0.021 in A119b) are diluted by 423, 26 and 74 uninformative pairs.
- **The bar is weak.** It passes from k ≥ 17 in A119b (null maximum 16), k ≥ 4 in testis and k ≥ 2 in gorilla.
  - A119b would pass with about 80% of its covered copy pairs split, so the test cannot detect F1's failure mode.
  - With F1 already known (131 / 148), a pass was nearly certain.
- **These numbers were already defined elsewhere, and are now seen.** k / n is figure 7's Liftoff-pair sensitivity of
  the de novo families (`docs/archive/2026-09/PREREG_genome_wide_families_2026-09-25.md`, Amendment 1, item 4). It is
  also the default arm of F6.6 (item 3).
- **For a reuse.** To put recall into the bar, do not match decoys on coverage, or test the unconditioned rate against
  an external expectation.
