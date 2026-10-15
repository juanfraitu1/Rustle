# O3: a read counts only for the family holding its best alignment (2026-10-15)

Pre-registration, written before any `--psv-family-best` run and before the fresh held-out inputs were built. Code under
test: branch `o3-family-best`, commit `7c9d2e8b` (`copy_assign --psv-family-best`, opt-in). User goal (2026-10-15): "keep
only reads whose best alignment is in the family, or cluster related families together or require the unexplained alleles
to be carried by most of the cluster's reads". This test is the first option.

## Question

O3 calls a PSV cluster a copy missing from the genome when nothing in the genome explains it better. The real-data test
(register 1294, `docs/PREREG_o3_realdata_2026-10-15.md`) left `copy` clusters in about 30% of families whose copies were
all in the genome. Its follow-up (register 1295) traced most of them to reads from other loci, placed on the family through
secondary alignments by `--no-as-tied-only`.

For each read, the aligner already gives part of O3's test: a read whose best alignment lies at another locus is better
explained there. `--psv-family-best` lets a read into family F's clustering only when F holds one of its best alignments:
- a primary or secondary record carrying the read's highest AS (ties count wherever the primary went; supplementary records
  are segments, not placements);
- with an aligned base inside a copy of F.

Does this:
1. remove most `copy` clusters from intact families, and
2. keep the detection of copies removed from the genome?

## Known before running (recorded now)

- **Register 1295** classified the `copy` clusters of the added-back arms of heldout2 (57) and dev (27).
- **`primary_local` tallies** were computed while designing the rule, from the existing assignments. `primary_local` is the
  primary-only form of the rule, without ties. It is the share of a cluster category's reads whose primary alignment
  overlaps a copy of the family:

  | substrate | category | reads with `primary_local` |
  |---|---|---|
  | heldout2 | other family's copy | 34 / 2,292 |
  | heldout2 | not in the genome | 9 / 2,276 |
  | heldout2 | catalog lacks the exons | 165 / 320 |
  | heldout2 | own-copy variant | 27 / 129 |
  | heldout2 | reads in all other clusters | 16,829 / 26,169 |
  | dev | GWFAM175, absent from the reference | 469 / 469 |
  | dev | minority segment | 219 / 930 |
  | dev | catalog lacks the exons | 11 / 16 |
  | dev | own-copy variant | 4 / 39 |
  | dev | not in the genome | 19 / 59 |
  | dev | reads in all other clusters | 13,112 / 56,545 |

  So **heldout2 and dev are development substrates for this rule**. Their old-rule results are also known: recovered 4 / 12
  and 15 / 25; families with a `copy` cluster 17 / 60 and 16 / 53.
- **Without the flag**, this branch's `copy_assign` reproduces dev's added-back arm byte for byte (9 output files). The
  binary is frozen at `/mnt/linuxdisk/tmp/o3_fb_build/copy_assign.fb`.
- **Not seen:** any `--psv-family-best` output, and which families the fresh held-out rule picks.

## Method (fixed)

**Fresh held-out families** (`heldout3`):
`prepare.py heldout --min-copies 2 --min-expressed 2 --exclude-picked .../heldout2/picked.tsv --seed 1 --cap 60`.
- Families of the shipped KB3781 catalog with >= 2 copies.
- At least 2 copies with >= 20 primary reads on their exons. Register 1294's rule required every copy; that is now relaxed.
- Every locus hull <= 200 kb, and no copy overlaps a copy of A13.
- Not one of heldout2's 60 families. The unrun 14-family build is not excluded: no output of it was ever made.
- The last expressed copy by (chrom, start) is masked. Silent copies stay in the catalog as survivors.
- A family whose masked hull overlaps any other catalog copy is dropped.
- More than 60 families: a seeded sample of 60.

**Reads, masking and re-mapping** are as in register 1294:
- up to 500 reads per copy;
- `bedtools maskfasta`;
- `minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`, the masked genome's own index, `GGO.splice.mmi`
  for the added-back arm.

**Runs.** Four on heldout3: {removed, added-back} × {old rule, new rule}. The old rule is register 1294's command:
`copy_assign --bam BAM --fasta GENOME --regions regions.txt --families COPIES.tsv --copies-fa COPIES.fa --no-as-tied-only
--psv-clusters --region-threads 1 --skip-poa-diagnostic --out PREFIX`. The new rule adds `--psv-family-best`. Both use the
frozen binary. On heldout2 and dev, only the new rule is run; their old-rule outputs exist.

## Scoring (fixed now; `bench/o3_realdata/score_fb.py`)

- A `copy` cluster is an `unexplained` cluster with est_div > delta = 0.00958.
- **Placed**, **eligible** (>= 6 masked-copy reads placed on the family) and **recovered** are register 1294's rules:
  - one `unexplained` cluster holds >= 0.8 of the masked copy's PLACED reads, with purity >= 0.9;
  - the denominator is the placed reads, which the filter does not change, so a read the filter drops counts against
    recovery;
  - the eligible set is the same under both rules.
- `cleared` is register 1294's rule.

| id | bar (fresh held-out, heldout3) |
|---|---|
| S1 | under the new rule, <= 0.10 of the families have a `copy` cluster in the added-back arm |
| S1b | the new rule's number of such families is <= 0.5 x the old rule's (paired, same families) |
| S2 | eligible masked copies recovered under the new rule >= those recovered under the old rule - 1 |
| G | denominator guard: families with an `explained` cluster in the added-back arm, new rule >= 0.9 x old rule |

- S2 needs >= 5 eligible masked copies. With fewer it is reported as not measurable, and the reading uses S1, S1b and G.
- S1b with 0 old-rule families is not measurable.
- G is there because a filter that deleted every read would clear every `copy` cluster. G checks that the families' own
  reads still form clusters the catalog explains.

## Readings (fixed now)

1. **S1, S1b, S2 and G pass.** Claimed, quoted with the counts: "keeping each read in the family that holds its best
   alignment removes most `copy` clusters of intact families on real gorilla reads, without losing detections". Recommended
   as the default inside `--psv-clusters`. That is the user's decision; the flag stays opt-in until then.
2. **S1, S1b and G pass, S2 fails.** Specificity gained at a cost in detection. Reported and not recommended; the lost
   detections are described.
3. **S1 or S1b fails, or G fails.** Not claimed. The remaining `copy` clusters are classified with register 1295's tools.

**Reported, not barred:**
- R1 and R2 of register 1294 under both rules;
- hosted (family, read) pairs kept by the filter;
- the masked copy's placed reads kept;
- `copy` clusters per family;
- the development substrates' numbers (heldout2, dev).

## Results

Added after the runs. Nothing above this heading changes once a statistic has been seen.
