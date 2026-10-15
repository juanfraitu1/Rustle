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
