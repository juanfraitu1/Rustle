# Pre-registration: the O3 -> O2 sample-copy step, a second held-out test (HG002 pseudoautosomal families, 2026-10-19)

**User (2026-10-19):** "Second held-out test (Recommended)", on the plan to run the frozen original rule with the same bars
on the six unused pseudoautosomal families.

## Why

**The first held-out test failed** (register 1333). On 19 families, 82 placements were right and 4 wrong (0.047 against the
0.02 bar), with clear wins in six families.

**Development found no rule that separates the 4 errors from the clear cases:**
- cluster purity (the register 1333 correction);
- teaching labels only with reads resolved by >= 10 AS (register 1334).

So the original rule is tested again, unchanged, for an independent estimate of its precision.

**The framing** (user, the same day): an O3 step feeding O2. Reads that CHM13 cannot explain are clustered by their own
variants and handed back to a CHM13 copy, or left as candidate copies CHM13 lacks. Phased data is never an input; HG002's
haplotypes are the truth only.

## The method

`bench/o3_graph/samplecopy.py` as in the first held-out test (2f23fe47), run without `--margins`. Its rule:
- the anchor is the copy most untied reads resolve to;
- the allele floor is 0.5 / kmax;
- the error rate is the reads' background mismatch rate;
- tied reads need two decisive sites and a posterior of at least 0.99 on their own evidence;
- labels need a confident majority of all untied reads.

## The held-out families (by the same fixed rule, never run by this method)

The pseudoautosomal HG002 families with at least 5 AS-tied reads that have a donor copy.

**6 families, 622 tied reads with a donor copy, 6,681 donor-tied:** AKAP17A 10, CD99 8, CD99P1 24, DHRSX 544, PLCXD1 20 and
SLC25A6 16.

**Here the copies are chrX and chrY,** HG002's maternal and paternal alleles of one recombining locus. The clusters are
those two haplotypes. CHM13 v2.0's chrY is HG002's own Y.

## The run

```
H=/mnt/linuxdisk/tmp/o3_graph/hg002proof; O=/mnt/linuxdisk/tmp/o3_graph/samplecopy_par
for F in S19AKAP17A S19CD99 S19CD99P1 S19DHRSX S19PLCXD1 S19SLC25A6; do
  python3 bench/o3_graph/samplecopy.py run /mnt/linuxdisk/tmp/o3_graph/r19seq/families19cal_full2_all.json seq19 $F $H/chm13.bam \
    $H/tied.tsv $H/untied.tsv /mnt/linuxdisk/home/juanfraitu/npip_ladder/idx/target.fa $O
done
python3 bench/o3_graph/samplecopy.py combine $H/truth.tsv $O /mnt/linuxdisk/tmp/o3_graph/fast/hg002/par/INTACT <families> $O/combine.json --tied $H/tied.tsv
```

## Bars (pooled; as the first held-out test)

- **S1, precise:** wrong / placed <= 0.02. Floor: 50 placed.
- **S2, calibrated:** `own_tie` placed <= 0.05.
- **S3, adds without harm:** v14 + sample copies has more right placements than v14 and no more wrong ones. v14 is the PAR
  world's `o2.integv.m2.tsv`.

**Verdict:**
- **PASS:** S1 is powered and holds, and S2 and S3 hold.
- **FAIL:** a powered bar fails.
- **UNDERPOWERED:** otherwise.

**Rules:** the code is frozen at the commit that adds this file. A bug is recorded as an amendment before any rerun.

## Results (2026-10-19, run once; `/mnt/linuxdisk/tmp/o3_graph/samplecopy_par/`)

**Verdict: PASS.** S1, S2 and S3 hold, and S1 is powered.

| pooled, 6 families (622 tied reads with a donor copy) | right | wrong | abstain | own_tie placed |
|---|---|---|---|---|
| sample copies | **508** | **0** | 114 | 0 / 6,681 |
| v14 | 1 | 0 | 621 | 0 |
| v14 + sample copies | 509 | 0 | 113 | 0 |

- **S1:** 0 of 508 placements are wrong.
- **S2:** 0 of the 6,681 donor-tied reads are placed.
- **S3:** v14 + sample copies makes 509 right placements against v14's 1, with no wrong one.

**One family carries it.**
- **DHRSX:** 852 reads, 41 informative columns and 2 clusters, both labelled. 508 of its 544 tied reads with a donor copy
  are placed right; v14 places 1.
- **The other five abstain throughout:**
  - AKAP17A, CD99 and SLC25A6 have a single cluster (SLC25A6 has 4 informative columns over 7,535 reads, CD99 has 35);
  - CD99P1's clusters are unlabelled;
  - PLCXD1's 4 labelled clusters place nothing (none of its tied reads passed the two-site and posterior rules).

**What it claims.** On a held-out stratum of real reads, the O3 -> O2 step resolves AS-tied reads without error when the
sample's own variants separate the copies. Here the copies are HG002's chrX and chrY alleles of DHRSX. Where nothing
separates them, it abstains.

**Limits.**
- The pass rests on one family.
- Separating chrX from chrY in the PAR is read-backed phasing of one recombining locus, not paralog assignment. It tests the
  method's machinery: clustering by the sample's own alleles and labelling by confidently placed reads.
- Together with register 1333 (82 / 4 on paralog families), the step's real-data precision is 590 right / 4 wrong over the
  two held-out tests.
