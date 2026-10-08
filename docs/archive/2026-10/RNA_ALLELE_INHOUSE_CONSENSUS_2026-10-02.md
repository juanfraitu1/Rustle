# The in-house consensus in place of IsoCon (Amendment 11): FAILS, and why — 2026-10-02

Prereg: `docs/PREREG_rna_allele_haplotype_count_2026-10-01.md` Amendment 11 (commit 092e2c8b, before any run). Work dir
`/mnt/linuxdisk/tmp/rna_allele/linktest_ih/` (shares Amendment 7's panel, reads, masked genome and R-arm alignments; `ih/` holds the
`missing_copy_flag` scan and verdicts). Outputs copied to `docs/RNA_ALLELE_INHOUSE_CONSENSUS_score.out.txt`.

## What ran

`missing_copy_flag` (defaults: `--delta-min 0.01 --min-reads 10 --min-sub 3`) on the 148 surviving copies of the 53 families, over the
same reads IsoCon saw, aligned to the same masked genome: **148 loci scanned, 23 fired** (15 `fired`, 5 `fired_both`, 3
`fired_structural`), 125 `no_mixture`; **23 patched consensus sequences** (verdicts: 15 `reference_absent_candidate`, 8 `scattered`).
Chain unchanged from Amendments 7-8: 19 flagged -> 15 new-copy contigs (10 D-derived, 5 survivor-derived) -> 14 components (9 D-derived).
IsoCon at the same point: 497 D-derived new-copy contigs -> 54 D-derived components in 44 families.

## Registered result

| | R (masked) | **in-house, arm M** | IsoCon, arm M (Amendment 8) |
|---|---|---|---|
| D right | 0 | **2,534 (14.7%)** | 12,787 (74.0%) |
| D wrong | 10,218 | 8,165 | 3,599 |
| D unplaced | 7,068 | 6,587 | 900 |
| S false moves | 0 | 2 (0.00%) | 25 (0.06%) |

- **IH1 FAILS**: D right 2,534 against the bar of 10,230 (80% of IsoCon); false moves pass (0.00%). Overall rule: wrong D -20.1% -> MIXED.
- **IH2 FAILS**: one candidate per deleted copy in 1/2 families with >= 2 D contigs (IsoCon 33/41).
- Post hoc, counting any own-family candidate — hybrid or not — as a home for the deleted copy's reads: 3,122 / 17,286 (18.1%) in 12 of 53
  families; 8,907 sit on a reference locus and **5,257 are unmapped in the masked genome**.

## Why, by cause (deleted copies grouped by where their reads land in arm R)

| deleted copies | n | survivor locus fired | D-derived candidate |
|---|---|---|---|
| reads on a surviving copy at `de` >= 0.01 | 18 | 14 | 8 |
| reads on a surviving copy at `de` < 0.01 | 4 | 1 | 0 |
| **no read with its primary on any surviving copy** | **31** | 3 | 0 |

- **The locus-anchored scan cannot see most missing copies.** 31 of 53 deleted copies have no read whose primary alignment lands on a
  surviving copy of the family — their reads go to another locus, or are unmapped (5,257). `missing_copy_flag` scans piles at given loci;
  IsoCon's input was the family's read NET (any record on a family copy, secondaries included, plus unmapped reads), which is where these
  reads were. This is the first-order reason (31 of 45 missing families).
- **The patched consensus is a host-backbone hybrid, not the hidden copy.** Of the 14 fired loci in the >= 0.01 stratum, 6 produced a
  consensus whose best hit in the unmasked genome is still the SURVIVOR (reference patched at the consistent sites only), so it labels as
  survivor-derived and the deleted copy's reads placed on it count as wrong; the five survivor-only components (GWFAM104, GWFAM407,
  GWFAM331 x2, GWFAM490: sub-piles of 47-398 reads at 3-16% divergence) are exactly these. IsoCon's transcripts are the reads' own
  consensus and label as the deleted copy.
- **Below the floor, as designed:** 4 copies at `de` < 0.01 fall under `--delta-min` and under delta alike — undetectable by construction
  in both arms.
- The 3 `fired_structural` loci and the `n_bigins` / `n_rearr` columns show the scanner notices structure, but the consensus does not carry
  it.

## What this decides for the `candidates` stage

The in-house replacement for IsoCon cannot be a per-locus consensus. It has to reproduce the two things IsoCon's architecture gets right,
in our own code: (1) **a family-scoped read net** — reads with any record on a family copy, secondaries included, plus unmapped reads
attributed to the family by sequence (a minimizer index of the family's copies; the test above attributed them by truth); (2) **a de novo
consensus per read cluster** — reads clustered within the net by pairwise divergence (the same rule the merge step applies to contigs),
one POA consensus per cluster, clusters merged when the real-vs-error significance test (already O2's certificate) cannot separate them.
The chain downstream (flag, link, merge, floor, union representative) is unchanged. `missing_copy_flag` keeps its role as the per-locus
screen (verdict classes, editing / immunoglobulin / contamination screens, DNA-depth expectation); it is not the consensus source.

Register rows 1213-1215.
