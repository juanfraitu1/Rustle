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

## Results (2026-10-15)

Run with the frozen binary (`7c9d2e8b`) and `bench/o3_realdata/` on this branch. Outputs:
`/mnt/linuxdisk/tmp/o3_realdata/heldout3/`, plus `*_fb` and `score_fb.json` in `heldout2/` and `dev/`.

**Build.**
- `prepare.py` picked 40 families out of 326 candidates, after excluding heldout2's 60: 40 masked copies and 153 survivors.
  Reads: 7,510 masked-copy (D) and 15,692 survivor (S).
- Masked arm: 22,909 primary, 293 unmapped. Added-back arm: 23,202 primary.
- Cost: the masked index took 251 s at 22.6 GB; mapping, 8 parts of about 330 s; each `copy_assign` run, 226-293 s at
  <= 2.7 GB.

**The filter is the rule as specified.**
- `family_best.py`, an independent Python implementation, keeps exactly the pairs the Rust runs keep:
  - heldout3: 23,180 / 25,537 (added back) and 15,069 / 18,349 (removed);
  - heldout2: 26,793 / 31,798 (added back);
  - dev: 58,352 / 59,105 (added back).
- No clustered read is outside its kept set.
- The prototype re-clusters every remaining `copy` cluster at Jaccard 1.0.

| heldout3 (40 fresh families) | old rule | `--psv-family-best` | bar | |
|---|---|---|---|---|
| families with a `copy` cluster (added back) | 10 | 6 = 0.15 | S1 <= 0.10 | **FAIL** |
| the same, paired | 10 | 6 | S1b <= 5 | **FAIL** |
| `copy` clusters | 35 | 11 | | |
| families with an `explained` cluster | 36 | 36 | G >= 32.4 | pass |
| eligible masked copies recovered (24 eligible) | 7 | 5 | S2 >= 6 | **FAIL** |
| cleared once added back | 5 of 7 | 4 of 5 | | |

**Lost detections.** MCL117 (recall 0.879 -> 0.515) and MCL31 (0.841 -> 0.305). Once their copy is masked, their reads'
best alignment lies outside the family.

| development | heldout2 (60) old -> new | dev (53) old -> new |
|---|---|---|
| families with a `copy` cluster | 17 -> 10 | 16 -> 16 |
| `copy` clusters | 57 -> 18 | 27 -> 27 |
| recovered | 4 -> 3 of 12 (MCL126 lost) | 15 -> 8 of 25 (7 lost) |

### Verdict (2026-10-15)

**Reading 3: not claimed.**
- S1 and S1b fail: 6 of 40 families keep a `copy` cluster.
- S2 fails: 2 of 7 detections are lost.
- G passes.

The filter removes every cluster made of another family's reads (0 of the 11 left, against 35 of 57 in heldout2 under the old
rule). But it costs detections: when a copy is missing, its reads' best alignment often lies in another family.

**Known before running, corrected.** `primary_local` is always 0 for families assigned through `copy_assign`'s
cross-chromosome pass:
- dev, multi-chromosome families: 44,367 of 44,367 rows read 0, though the filter keeps them;
- one-chromosome families: it agrees with the filter on 14,076 of 14,091 rows.

So the tallies above understated what the rule keeps on dev. The column's doc comment says every per-family read rate must be
taken over it; that is a defect of its own, not fixed here.

### The remaining `copy` clusters (register 1295's tools, as reading 3 requires)

heldout3 under the new rule: 11 clusters in 6 families, 270 reads. None is another family's reads.

| what the cluster is | clusters (reads) |
|---|---|
| chimpanzee transcripts, MCL45 (identity to chimp 0.997-0.998) | 2 (51) |
| MCL45, closest to chimpanzee (0.97) or nothing (0.85); MCL125 (0.92 to every ape) | 4 (89) |
| gorilla sequence outside both assemblies, MCL79 (0.994) | 1 (4) |
| the catalog lacks the exons: MCL63, MCL79 | 2 (66) |
| a minority segment, MCL31 | 1 (54) |
| an own-copy variant, MCL32 | 1 (6) |

### Post hoc: would a variation graph with every isoform help? (user question; descriptive)

**More reference content, same exact rule** (`genomic_explain.py`). Register 1294's explained rule was asked against
three references for heldout3's 136 unexplained clusters (added back, old rule):

| reference | clusters explained |
|---|---|
| G0, each copy's spliced representative (as the method does) | 0 |
| G1, each copy's genomic locus +- 2 kb: every isoform, intron and PSV of the family | 1 |
| G2, G1 plus every locus where the cluster's reads align | 4 |

G1's one explanation is MCL0's own copy 3, an isoform its representative lacks, in an `allele_or_copy` cluster (est_div
0.0005). G2 adds three clusters explained by another family's copy: MCL51 twice by MCL107, and MCL95 by MCL23. The exact
rule needs one locus to carry all of a cluster's
alt sites with no conflict: 20-338 per `copy` cluster, over 100 in 27 of 35. Reads 0.1% away from their own source locus almost always miss somewhere. So the
representative is not the bottleneck; the rule's exactness is.

**A tolerant rule** (`tolerant.py`). An unexplained cluster counts only when its reads' median best-alignment divergence
(minimap2 `de` at the primary, i.e. to the nearest locus anywhere in the genome) exceeds delta:

| old-rule clusters, as run -> tolerant | heldout3 (40) | heldout2 (60) | dev (53) |
|---|---|---|---|
| families with a `copy` cluster | 10 -> **2** | 17 -> 6 | 16 -> 9 |
| `copy` clusters | 35 -> 7 | 57 -> 14 | 27 -> 17 |
| masked copies recovered | 7 -> 5 of 24 | 4 -> 3 of 12 | 15 -> 11 of 25 |

- Held out, the 2 remaining families are MCL45 (chimpanzee transcripts) and MCL79.
- Dev keeps GWFAM175, the genuine reference-absent copy.
- The lost recoveries are mostly `allele_or_copy` copies, already within delta of another locus (MCL22, MCL97, 3 of dev's
  4). MCL10's masked copy is labelled `copy`, but its reads still align at 0.9985 somewhere else in the masked genome:
  by O3's definition it is not absent.
- This was measured on substrates already seen, so it is a hypothesis for a new prereg on a fresh substrate. The gorilla
  pool is exhausted under these rules; human A119b is next.

**For a variation graph**, this says:
- The graph is the natural structure for the tolerant rule: a cluster is a missing copy when its consensus is more than delta
  from every path.
- Its content gain alone is small (G1 explains 1 of 136). The gain comes from the distance rule, which the aligner's best
  alignment already measures genome-wide.
- What a graph would add: both haplotypes, so allelic reads are explained; related families in one structure; and a clean
  path definition for the thesis.
