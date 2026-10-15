# O3 as the advisor defines it: multi-mapping reads that no reference copy explains — design

Date 2026-10-15. Branch `laptop/o3-psv-multimappers`. Status: implemented as `copy_assign --psv-clusters` (`rustle::missing_copy::psv_clusters`, v2 with the linkage test); held-out PASS in Rust on base and hard worlds (`docs/PREREG_o3_psv_clusters_v2_2026-10-15.md`: absent copies 425/426, controls 0/152, separation 41/44, labels 135/135 and 262/262; v1 `docs/PREREG_o3_unexplained_2026-10-15.md`). Next: real families on the main machine.

## 1. The objective (advisor, relayed 2026-10-15)

O3 was never about recovering unaligned reads. It is: **cluster the multi-mapping reads of multi-copy loci by their PSVs; when
multi-mappers have no adequate original copy in the reference, that is evidence of copies absent from the genome — without WGS and
without copy-number estimates.** Experiments add and remove copies of those loci.

In the record's terms this is the **ABSORBED stratum** (whole-genome excision control: 64.2% of deleted copies' reads land on the
nearest paralogue; `docs/THESIS_OBJECTIVES.md`, O3 RESTATED), which unmapped-read methods cannot see. The unmapped-rescue chain and
`candidate_copies` serve the other stratum (ORPHANED) and are secondary here.

## 2. What already exists (O2's machinery) and what was already tried

- **"No adequate copy" per read** = O2's origin certificate (`copy_assign`, `src/family.rs` ~24937-25006): the read's edits against
  its best candidate tested one-sided against Binomial(n, 0.003) at alpha 1e-3; output `origin_rejected` (0/1) in
  `<out>.assignments.tsv`. Real reads already show it is common: 606 of 607 abstaining exact-tied gorilla molecules are
  origin-rejected (register 775); 285 of 433 low-MAPQ NPIP reads are explained by no unit (register 702).
- **Clustering the unexplained reads** = `copy_assign --flag-missing-copies` (`detect_missing_copy_pairs`, ~71782): origin-rejected
  reads grouped by best copy Y, re-aligned to Y's locus, consistent mismatch sites (≥ 3 reads and ≥ 50% of covering reads) per kb
  against Y's own reads.
- **The excision experiments of 2026-09 removed the copy from the CATALOG ONLY** (`docs/archive/2026-09/PREREG_excision_2026-09-05.md`:
  "removed from the family's candidate set; nothing else changes"). The genome and the BAM kept the copy, so its reads stayed at its
  own locus. NPIP / ZNF875 (siblings 0.966-0.991): the consistent-site signature held (§6fl/§6fn). Five tandem arrays (siblings
  0.909-0.969, §6kv): 5/5 "orphaned" — the reads touch no remaining candidate, because they are still placed on the removed copy's own
  locus and have no tied placement on a sibling. That outcome is a property of catalog-only removal; a copy truly absent from the
  reference sends its reads to the siblings (absorbed) or leaves them unmapped. The genome-level control (`o3_excise`, 162 TWO-copy
  families, genome masked, reads re-mapped) has one survivor per family, so no tie and no abstention test is possible there.
  ⟹ **The advisor's experiment in its proper form — remove the copy from the GENOME, re-map, with ≥ 3 copies — has never been run.**
- k+1 MEC clustering of the PSV matrix is closed (§6ai, four grounds). The design below is a goodness-of-fit test per read plus a
  consistency test per group, not model selection over k; it inherits §6ai's lesson that reference-PSV columns cannot show a
  missing copy's own variants, so the consistency test works on read-vs-best-copy mismatches, not on PSV columns.

## 3. Known limits, stated before any run

- **Absorption floor.** A removed copy within ~0.6% of a sibling is absorbed with a valid certificate (register 711: 12 of 32 NPIP
  reads; register 796). That is also the allelic range (delta = 0.00958, p99 allelic divergence): below it, a missing copy and an
  allele absent from the reference are the same signal from RNA alone.
- **Orphan ceiling.** A copy far from every sibling leaves its reads unmapped (the other stratum).
- **The read set.** Exact AS ties (the O2 default gate) are NM-identical twins; reads of an absent copy that is closer to one sibling
  than to the others map UNIQUELY to it, with extra mismatches. Two arms are therefore measured: **T** = the shipped contested set
  (exact ties), the advisor's literal "multi-mappers"; **P** = every molecule placed on a family copy (primary or tied), so the origin
  certificate sees absorbed reads that are not ties.

## 4. Experiments (famsim worlds first, then real families on the main machine)

Per family, the reference holds copies A1..Ak (k = 2-4, some pairs near-identical so that ties exist); the reads come from:
| scenario | truth | expected |
|---|---|---|
| NULL | every expressed copy in the reference | no read unexplained beyond the error null; no cluster |
| REMOVED@d | copy E expressed, absent from the GENOME, d = 0.3-8% from its nearest sibling | E's reads origin-rejected and clustered (above the floor) |
| ADD-BACK | REMOVED's world with E put back in the reference | E's reads assigned; no cluster |
| DECOY | an extra unexpressed copy added to the reference | nothing changes |
| ALLELE@d | a heterozygous allele of A1 absent from the reference, ~half of A1's reads | a cluster labelled allele, not copy |

Measured per family: the fate of the absent copy's reads (silently assigned to a sibling / origin-rejected / tied / unmapped); whether
the family is flagged and whether the flagged group is E's reads; flags in NULL / ADD-BACK / DECOY; how ALLELE is labelled.

## 5. Plan

1. World generator (`bench/o3_unexplained/world.py`, famsim) and harness: map reads (`sim.MM2`), `copy_assign --families` (+ the P-arm
   read set; to be implemented if no existing flag widens the gate), `--flag-missing-copies`, scorer against read-name truth.
2. Dev worlds (seeds 31-33): get the machinery right, freeze every rule.
3. Pre-registration with held-out bars; held-out seeds; outcome.
4. Main machine: the same with real families — remove real copies from the genome, re-map the real reads, add them back.
