# Held-out family test re-run on CAT/Liftoff v2.0 (2026-10-01)

Step 3 of `docs/CAT_RERUN_PROTOCOL_2026-10-01.md`. The benchmark is Test 2 of `docs/PREREG_heldout_families_2026-09-20.md`
(result: `docs/HELDOUT_FAMILIES_RESULT_2026-09-20.md`, registered against NCBI RefSeq; it stays as registered). Test 1
(symbol-root truth) was declared void in the pre-registration and is not re-run. Nothing was re-tuned (R6). The canonical
repo, the original outputs in `/mnt/linuxdisk/tmp/heldout/` and `/mnt/linuxdisk/tmp/regress/`, `score.py` and `lib.py`
are unchanged.

## Result

**The pre-registered bar holds on CAT, as it did on RefSeq.** The bar is: held-out pooled F (chr2 + chr8 + chr10)
minus development pooled F (chr5 + chr7 + chr21) of at least -0.10 means HOLDS; -0.25 to -0.10 means PARTIAL; below
-0.25 means FAILS.

| | truth universe | held-out: families, touched, **F**, sens, prec, exact | development: families, touched, **F**, sens, prec, exact | delta F | verdict |
|---|---|---|---|---|---|
| **RefSeq (registered)** | Soto by `Gene Name` | 20, 20, **0.7016**, 0.808, 0.722, 3 | 11, 11, **0.7183**, 0.819, 0.730, 4 | **-0.017** | HOLDS |
| **CAT** (primary, R5) | Soto by `Gene ID` | 41, 38, **0.7523**, 0.780, 0.776, 17 | 23, 22, **0.7194**, 0.718, 0.763, 6 | **+0.033** | HOLDS |
| CAT, registered families only | the same 20 + 11 families, members by `Gene ID` | 20, 20, **0.7668**, 0.824, 0.787, 6 | 11, 11, **0.7924**, 0.805, 0.831, 3 | -0.026 | HOLDS |
| decomposition: RefSeq clusters, CAT truth | `Gene ID`; RefSeq members re-keyed by R1 | 41, 36, 0.6359, 0.631, 0.688, 10 | 23, 20, 0.5710, 0.534, 0.639, 2 | +0.065 | (not a test) |
| same, registered families only | | 20, 20, 0.6726, 0.660, 0.754, 2 | 11, 11, 0.7183, 0.685, 0.796, 2 | -0.046 | (not a test) |

Pooled F, sensitivity and precision are the means over all of the arm's truth families (the registered definition; an
unmatched family scores 0 and is kept). The secondary bar (at least one exact recovery on a zero-exposure chromosome,
chr2 or chr8) is met on CAT: chr2 has 4 and chr8 has 12 (2 and 3 on the registered families).

**The two headline numbers are not on the same families.** The `Gene ID` join roughly doubles the truth: held-out 20 -> 41
families (97 -> 255 members), development 11 -> 23 families (58 -> 150 members), and chr21 becomes scoreable (3 families).
The CAT truth contains every registered family and every registered member; nothing was lost. The "registered
families only" row compares like with like: CAT gives higher F on both arms (+0.065 held-out, +0.074 development) and the
gap between the arms is about the same as on RefSeq (-0.026 vs -0.017).

The decomposition rows score the ORIGINAL RefSeq clusters against the CAT truth (each RefSeq cluster member is replaced
by the CAT gene it maps to under R1, `chm13v2.0_CAT_Liftoff.refseq_map.tsv`; a member with no CAT gene can only cost
precision). On the registered families the larger truth alone moves held-out F 0.7016 -> 0.6726 and development
0.7183 -> 0.7183; the CAT clustering then gives 0.7668 / 0.7924. Averaged over the arm, re-keying the truth alone does
not raise F: it lifts a few families sharply (FAM90A, TAF11L, whose RefSeq truth missed most copies) and lowers others
whose truth gains members RefSeq never annotated. The net gain comes from clustering on the CAT loci and gene models
(which include loci RefSeq lacks). R1 maps are approximate for "partial" and "weak" matches, so these two rows are a
decomposition, not a benchmark.

## Sanity check: the scorer reproduces the registered numbers

`bench/annotation/heldout_cat.py score --join name` (the registered path: genes keyed on the GFF `Name=`, Soto joined by
`Gene Name`, built from `score.heldout_load_genes`, `score.heldout_predicted_clusters`, `lib.soto_gene_family`,
`lib.families_on` and `lib.bipartite_families`, all unchanged), run on the original RefSeq cluster files, writes
`chrN_soto.json` files that are **byte-identical** to the registered `/mnt/linuxdisk/tmp/heldout/chrN_soto.json` for all
seven scored chromosomes (chr2, chr8, chr10, chr5, chr7, chr16, chr17). chr21 is reported as not scoreable, as
registered. Pooled: held-out 0.7016 (20 families, 3 exact), development 0.71836 (11 families, 4 exact). The
registered 0.7183 is the same number computed from the rounded per-chromosome means ((9 x 0.7825 + 2 x 0.4296) / 11 =
0.71834); the per-family values are identical. The check was repeated after the script's last edit, with the same
result.

**Binary check.** The current `mcl_families` (md5 `baee52a4...`) differs from the one that ran on 2026-09-20 in one
default: `--min-cov-shorter` became 0.70 on 2026-09-29. With `--min-cov-shorter 0` it reproduces all 8 original cluster
files byte for byte (chr2, chr8, chr10, chr5, chr7, chr21 from the original PAFs; chr16, chr17 from
`/mnt/linuxdisk/tmp/regress`); the `.loci.tsv` files are also identical. The CAT runs therefore pass
`--min-cov-shorter 0`. Every other option is the binary default, as in the original (the `params.tsv` files differ from
the originals only in paths, counts and the new `min_cov_shorter` / `admitted_by_containment` rows).

## Per chromosome

Same rule (Sep-20 options), same scorer. "CAT reg." = CAT scored on the families the registered run scored.

| chromosome | arm | RefSeq: fam, touched, F, sens, prec, exact | CAT: fam, touched, F, sens, prec, exact | CAT reg.: F (exact) |
|---|---|---|---|---|
| chr2 | held out | 7, 7, 0.7295, 0.770, 0.792, 2 | 12, 12, 0.7762, 0.766, 0.846, 4 | 0.7320 (2) |
| chr8 | held out | 5, 5, 0.6377, 0.950, 0.583, 1 | 19, 17, 0.7892, 0.856, 0.771, 12 | 0.8158 (3) |
| chr10 | held out | 8, 8, 0.7172, 0.752, 0.747, 0 | 10, 9, 0.6534, 0.653, 0.701, 1 | 0.7668 (1) |
| chr5 | development | 2, 2, 0.4296, 1.000, 0.281, 0 | 7, 6, 0.7365, 0.752, 0.741, 2 | 0.8445 (0) |
| chr7 | development | 9, 9, 0.7825, 0.778, 0.830, 4 | 13, 13, 0.7200, 0.686, 0.804, 3 | 0.7808 (3) |
| chr21 | development | not scoreable (0 families) | 3, 3, 0.6768, 0.778, 0.639, 1 | not scoreable |
| chr16 | reference (NPIP) | 8, 7, 0.6698, 0.731, 0.691, 0 | 37, 30, 0.7067, 0.681, 0.769, 15 | 0.7779 (2) |
| chr17 | reference (TBC1D3) | 14, 13, 0.7328, 0.758, 0.733, 4 | 26, 21, 0.6047, 0.596, 0.654, 5 | 0.6796 (3) |

Reference arm pooled (not part of the pre-registered comparison; the registered document gave only the per-chromosome
rows): RefSeq 22 families, F 0.7099, 4 exact; CAT 63 families, F 0.6646, 20 exact; CAT on the registered 22 families,
F 0.7154, 5 exact.

On the registered families the worst development chromosome changes: on RefSeq chr5 (TAF11L, GUSBP) was the worst in the
panel at 0.430; on CAT it is 0.845. On the full CAT universe the worst held-out chromosome is chr10 (0.653) and the worst
overall is chr17 (0.605).

## Why the universes differ

`bench/annotation/heldout_cat.py universe` puts the two truths side by side. For every member of the CAT truth it records
why the name join did or did not count it (`universe.summary.tsv`, `universe.genes.tsv`).

| chromosome | families (name / id) | members (name / id) | id members also in the name truth | symbol absent from RefSeq on this chromosome | name shared by several Gene IDs (collapsed) | name carries > 1 Family ID (excluded) | family < 3 under the name join |
|---|---|---|---|---|---|---|---|
| chr2 | 7 / 12 | 32 / 61 | 32 | 23 | 2 | 1 | 3 |
| chr8 | 5 / 19 | 19 / 115 | 19 | 65 | 9 | 0 | 22 |
| chr10 | 8 / 10 | 46 / 79 | 46 | 13 | 0 | 19 | 1 |
| chr5 | 2 / 7 | 11 / 65 | 11 | 12 | 33 | 0 | 9 |
| chr7 | 9 / 13 | 47 / 76 | 47 | 25 | 0 | 0 | 4 |
| chr21 | 0 / 3 | 0 / 9 | 0 | 7 | 0 | 0 | 2 |
| chr16 | 8 / 37 | 43 / 186 | 43 | 105 | 14 | 1 | 23 |
| chr17 | 14 / 26 | 53 / 116 | 53 | 49 | 3 | 0 | 11 |
| chr6 (not run) | 0 / 1 | 0 / 4 | 0 | 2 | 0 | 0 | 2 |

- Every name-join member has a CAT counterpart (the "also in the name truth" column equals the name-join member count on
  every chromosome), and no name-join member is a RefSeq symbol whose Soto gene lies on another chromosome. The CAT truth
  is a strict superset of the registered truth.
- **Symbol absent from RefSeq** (299 of the panel's CAT members): Soto's names are GENCODE/Ensembl names. RefSeq lacks the
  clone-based lncRNA and pseudogene names (AC..., AL..., FO..., `MSTRG.*` models) and names many copies `LOC...`.
- **Name shared by several Gene IDs** (61): CAT/GENCODE gives several paralogs one name (e.g. about 30 copies called
  `TAF11L5` on chr5; `FAM90A13P`, `DUX4L24`; `DEFB103A`/`DEFB103B` for three chr8 copies). The name join counts each name
  once. This is why 8 small chr8 families (DEFB103, DEFB104, DEFB105, DEFB106, DEFB107, DEFB4, SPAG11, PRR23D) enter the
  CAT truth: each has 3 Gene IDs, mostly Liftoff records of the 8p23.1 defensin repeat, but only 2 names. They are
  distinct loci (none of their member pairs overlap on exons).
- **Name carrying more than one Family ID** (21): the name rule excludes the name. On chr10 `DUX4L24` names copies that
  Soto puts in two families (ID_330 and ID_347), so 19 members were dropped; by Gene ID each copy has one family.
- **Family below 3 under the name join** (75): the members' names are in RefSeq but the family had fewer than 3 distinct
  names on the chromosome.

## Per-family changes that matter

Registered families, RefSeq vs CAT (truth size / predicted cluster size / hit, F):

| family | arm | RefSeq | CAT | what changed |
|---|---|---|---|---|
| chr8 ID_356 FAM90A | held out | 4/46/4, 0.160 | 56/56/56, **1.000** | Truth. The RefSeq cluster held 46 FAM90A records (19 named FAM90A*, 27 named `LOC...` with a FAM90A description), so the name truth saw only 4 and called the cluster imprecise. The RefSeq clusters scored against the CAT truth already give 0.784. |
| chr5 ID_210 TAF11L | development | 8/43/8, 0.314 | 43/47/40, **0.889** | Truth, the same pattern (29 `LOC...` members; RefSeq clusters vs CAT truth 0.837). |
| chr10 ID_330 DUX4L | held out | 11/28/11, 0.564 | 28/34/28, **0.903** | Truth grows 11 -> 28; the CAT cluster contains all 28 (RefSeq clusters vs CAT truth 0.500). |
| chr16 ID_2 ABCC6 | reference | 3/0/0, 0.000 | 3/2/2, 0.800 | Clustering: CAT pairs ABCC6P1 with ABCC6P2; none of the three RefSeq records was in any cluster. |
| chr17 ID_136 CCDC144 | reference | 3/4/2, 0.571 | 5/4/4, 0.889 | Truth and clustering. |
| chr5 ID_163 GUSBP | development | 3/8/3, 0.545 | 6/9/6, 0.800 | Truth grows; the cluster carries the new members. |
| chr17 ID_114 RDM1P | reference | 3/2/2, 0.800 | 3/3/1, **0.333** | Clustering: under CAT, RDM1P1 and RDM1P4 are in no cluster and RDM1P2 sits with AC091132.4 and LRRC37A3. |
| chr17 ID_247 KRT17P | reference | 3/3/3, 1.000 | 6/4/3, 0.600 | Truth grows by three processed-pseudogene copies named `AL353997.4`, which are in no cluster; the cluster adds KRT42P. |
| chr2 ID_65 FAR2P / CYP4F | held out | 7/4/4, 0.727 | 9/5/3, 0.429 | Truth grows 7 -> 9 (two lncRNA models, unclustered); the CAT cluster keeps 3 FAR2P, loses the CYP4F pseudogenes and adds PLEKHB2 and AC073869.5. |
| chr16 ID_39 RRN3 | reference | 4/3/3, 0.857 | 7/3/3, 0.600 | Truth grows by two lncRNA models (clustered separately) and one pseudogene (unclustered); the same 3-member cluster. |
| chr10 ID_223 ANXA8 | held out | 3/2/2, 0.800 | 5/2/2, 0.571 | Truth grows by two lncRNA models (unclustered); the same pair. |

Of the 53 registered families, 22 improve, 16 are unchanged and 15 get worse. Restricting the truth to the registered
families changes no family's matched cluster or F (the matching on the full CAT universe gives the same values for them).

**New CAT families at F = 0** (held out 3, development 1, reference 11): chr8 ID_96 and ID_97 (clone-named pseudogenes
and lncRNAs, no cluster), chr5 ID_313 (CDH12P1-3, no cluster), and on chr16/chr17 families of processed or unprocessed
pseudogenes and clone-named lncRNAs (`AC140658.*`, `AC243829.*`, `AC133555.*`, IGHV1OR16-3/4, WASH6P with MIR6859-4,
TRIM16/TRIM16L, YWHAEP2/3). **chr10 ID_347** is the NPIPA/NPIPB shape again (advisor Q9): Soto splits the DUX4L24 /
MIR8078 array into a second family; our single DUX4L cluster is matched to ID_330, so ID_347 is unmatched.

**Exact recoveries** on CAT (17 held out) are dominated by chr8 (12): 7 of the 3-ID / 2-name families above (all but
PRR23D, F 0.667), the registered DEFB108, FAM90A and ZNF705, and the new ID_105 (clone-named lncRNAs) and ID_277
(ALG1L11P/12P/14P). They are real copies, but most are 3-member families of near-identical short genes; do not quote the
exact count without that.

**Folded records.** `mcl_families` folds annotation records that overlap on exon bases inside a cluster into one locus,
and the folded record is not a row of `clusters.tsv`, so a truth member folded away cannot be hit. This is a property of
the registered scorer, and it bites slightly harder on CAT: 10 of 255 held-out CAT truth members (RefSeq 5 of 97),
3 of 150 development (2 of 58), 16 of 302 reference (5 of 96). In the CAT truth, 28 families contain at least one pair of
members whose exons overlap (38 pairs, e.g. a gene and a lncRNA model over it).

## Inputs and graphs

CAT annotates 7-33% more gene bodies per chromosome than RefSeq gene + pseudogene records, so the graphs are larger:

| chromosome | bodies RefSeq / CAT | PAF records RefSeq / CAT | graph nodes / edges, RefSeq | graph nodes / edges, CAT | clusters (members), RefSeq | clusters (members), CAT |
|---|---|---|---|---|---|---|
| chr2 | 4,053 / 4,333 | 487,933 / 563,771 | 639 / 714 | 771 / 964 | 168 (441) | 181 (512) |
| chr8 | 2,347 / 2,622 | 306,343 / 360,555 | 393 / 2,375 | 472 / 4,619 | 70 (282) | 66 (313) |
| chr10 | 2,270 / 2,419 | 261,372 / 285,484 | 352 / 696 | 440 / 1,124 | 81 (235) | 91 (270) |
| chr5 | 2,742 / 3,074 | 472,482 / 573,677 | 433 / 1,400 | 507 / 1,625 | 87 (274) | 92 (292) |
| chr7 | 2,928 / 3,194 | 361,096 / 404,757 | 619 / 1,466 | 806 / 2,831 | 148 (456) | 181 (583) |
| chr21 | 1,003 / 1,123 | 32,338 / 59,562 | 270 / 4,242 | 314 / 7,944 | 27 (229) | 37 (263) |
| chr16 | 2,081 / 2,773 | 79,161 / 116,915 | 380 / 627 | 771 / 1,419 | 90 (271) | 144 (477) |
| chr17 | 2,538 / 3,155 | 77,703 / 110,554 | 430 / 478 | 679 / 705 | 90 (257) | 120 (330) |

## Supplementary: the current default `--min-cov-shorter 0.70`

This is not the registered rule. It was run because 0.70 has been the binary default since 2026-09-29. Same PAFs, same scorer:

| | held-out F (fam, exact) | development F (fam, exact) | delta |
|---|---|---|---|
| RefSeq | 0.6936 (20, 3) | 0.7227 (11, 4) | -0.029 |
| CAT | 0.7637 (41, 17) | 0.7047 (23, 6) | +0.059 |
| CAT, registered families | 0.7820 (20, 6) | 0.8114 (11, 3) | -0.029 |

Per chromosome on RefSeq, 0.70 raises chr2 (0.7295 -> 0.7415) and chr7 (0.7825 -> 0.7878), leaves chr5 unchanged, and
lowers chr8 (0.6377 -> 0.5927), chr10 (0.7172 -> 0.7147), chr16 (0.6698 -> 0.5908) and chr17 (0.7328 -> 0.7009). The
verdict is HOLDS in every row. I did not check how this relates to the measurements behind the 09-29 default (register
1006/1014).

## Deviations from the original run

1. `--min-cov-shorter 0` is passed explicitly, because the binary default changed after the original run (see the binary
   check). Without it the CAT run would not be the registered rule.
2. The all-vs-alls were split with `tools/mm2_shard.sh` (10 Mbp query shards against one shared index, resumed until
   exit 0); the original ran one minimap2 process. On chr21 the sharded PAF is byte-identical to a single
   `minimap2 -x asm20 -c --eqx -P -t 4` run (checked here). minimap2 is the same version (2.30-r1287).
3. Regions: the original awk filter (`gene` or `pseudogene`, `sort -u`) is applied to the CAT slim GFF, which has only
   `gene` features (all biotypes). Records with identical spans collapse to one body, as in the original (chr2: 4,336
   genes, 4,333 bodies).
4. `--gff` is the whole uncompressed CAT slim GFF (all contigs), as the original used the whole RefSeq GFF.
   `mcl_families` attaches exons by `gene=` and genes by `Name=`; under CAT both are the CAT gene id, and every node had
   an exon-union length (0 fell back to the span).
5. Soto is joined by `Gene ID` (R5) through a new script, `bench/annotation/heldout_cat.py`; `score.py` is unchanged. The
   member lookup keeps the registered probe order ((start+1, end), then (start, end)). `--exact-only` (the B4 fix) gives
   identical pooled numbers on both annotations.
6. chr16 and chr17 (the registered "reference" rows, clusters built in `/mnt/linuxdisk/tmp/regress` by the same recipe)
   were re-run as well. chr6 was not re-run: the original could not score it. Under the `Gene ID` join it would have
   one scoreable family (ID_306, 4 members).
7. The FASTA path is `winloci_data/Reference/chm13v2.0.fa`, a hard link of the `winloci_data/chm13v2.0.fa` the original
   used (same inode).

## Exact commands

```sh
W=/mnt/linuxdisk/tmp/heldout_cat
zcat winloci_data/gencode_chm13/chm13v2.0_CAT_Liftoff.slim.gff3.gz > $W/cat_slim.gff3     # md5 e0de3789...
cp rustle_target/release/mcl_families $W/mcl_families                                     # md5 baee52a4...
# per chromosome (chr21 chr16 chr17 chr10 chr8 chr7 chr5 chr2), $W/run_cat.sh:
bash $W/run_cat.sh prep chrN   # awk -F'\t' '$1=="chrN" && ($3=="gene"||$3=="pseudogene"){print $1":"$4"-"$5}' cat_slim.gff3 | sort -u > chrN.regions
                               # samtools faidx Reference/chm13v2.0.fa -r chrN.regions > chrN.bodies.fa
RLOCK_WAIT=150 BUDGET=400 bash $W/run_cat.sh mm2 chrN     # repeated until exit 0 (75 = call again):
  # MM2_SHARD_BUDGET_S=$BUDGET tools/rlock.sh heavy bash tools/mm2_shard.sh paf chrN.paf -x asm20 -c --eqx -P -t 4 chrN.bodies.fa chrN.bodies.fa
bash $W/run_cat.sh fam chrN
  # tools/rlock.sh heavy $W/mcl_families --paf chrN.paf --gff cat_slim.gff3 --min-exonic-bp 1 --min-shared-exon-frac 0.60 --min-cov-shorter 0 --out chrN_fam
bash $W/run_cat.sh fam070 chrN # supplementary: the same without --min-cov-shorter (default 0.70)

# scoring (worktree root; G = winloci_data/gencode_chm13, O = /mnt/linuxdisk/tmp/heldout, R = /mnt/linuxdisk/tmp/regress)
ARMS="--arm heldout=chr2,chr8,chr10 --arm development=chr5,chr7,chr21 --arm reference=chr16,chr17"
# sanity: the registered path on the original outputs
python3 bench/annotation/heldout_cat.py score --gff $R/chm13.gff --soto bench/soto/soto_famCN_S1C.tsv --join name \
  --clusters chr2=$O/chr2_fam.clusters.tsv ... --clusters chr16=$R/chr16_guided.clusters.tsv --clusters chr17=$R/chr17_guided.clusters.tsv \
  $ARMS --json-dir $W/sanity --families-tsv $W/sanity/families_refseq_name.tsv
# CAT
python3 bench/annotation/heldout_cat.py score --gff $G/chm13v2.0_CAT_Liftoff.slim.gff3.gz --soto bench/soto/soto_famCN_S1C.tsv --join id \
  --clusters chr2=$W/chr2_fam.clusters.tsv ... $ARMS --json-dir $W/score --families-tsv $W/score/families_cat_id.tsv \
  --names $G/chm13v2.0_CAT_Liftoff.genes.tsv
#   + --only-families $W/registered_families.txt (the 53 registered chrom/family pairs)   -> CAT, registered families
#   + --exact-only                                                                         -> B4 check
# decomposition: RefSeq clusters re-keyed to CAT ids, CAT truth
python3 bench/annotation/heldout_cat.py score --gff $R/chm13.gff --truth-gff $G/chm13v2.0_CAT_Liftoff.slim.gff3.gz \
  --translate $G/chm13v2.0_CAT_Liftoff.refseq_map.tsv --soto bench/soto/soto_famCN_S1C.tsv --join id --clusters (RefSeq) $ARMS
# universes
python3 bench/annotation/heldout_cat.py universe --refseq-gff $R/chm13.gff --cat-gff $G/chm13v2.0_CAT_Liftoff.slim.gff3.gz \
  --soto bench/soto/soto_famCN_S1C.tsv --chroms chr2,chr8,chr10,chr5,chr7,chr21,chr16,chr17,chr6 --out $W/universe
```

Inputs: CAT slim GFF md5 `a3051ee1a0aa698b52afd49c3841eb50`; refseq_map md5 `7c8e13f97aacfe6bfa526fed936d4a20`;
S1C md5 `afd30b7381cc53b1a67afd560212017e`. Wall time: all-vs-alls 48 s (chr21) to about 15 min (chr2, 3 resumed calls);
each `mcl_families` run under 15 s.

## Files (`/mnt/linuxdisk/tmp/heldout_cat/`)

- `run_cat.sh`, `cat_slim.gff3`, `mcl_families` (the copy used); `chrN.{regions,bodies.fa,paf}`; `chrN_fam.*` (registered
  rule), `chrN_fam070.*` (supplementary).
- `sanity/` (name-join JSONs, byte-identical to the registered ones; `families_refseq_name.tsv`); `binchk/` (current binary
  on the original PAFs: `*_mcs0.*` byte-identical to the originals, `*_mcs070.*` supplementary).
- `score/` (`chrN_soto.json` and `families_cat_id.tsv` for CAT; `families_cat_id_registered.tsv`;
  `families_refseqclusters_cattruth*.tsv`; `families_*_mcs070.tsv`); `universe.{summary,families,genes}.tsv`;
  `registered_families.txt`.
- The shard cache for these runs (about 5.2 GB, duplicating the PAFs) is under `/mnt/linuxdisk/tmp/mm2_shard_cache/`
  in the key directories named by the md5 of each `chrN.bodies.fa`; it can be deleted.
