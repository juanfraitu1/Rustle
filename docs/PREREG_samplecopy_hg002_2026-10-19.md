# Pre-registration: the sample's own copies place O2's ambiguous reads (HG002, held-out families, 2026-10-19)

**User (2026-10-19):**
- "can we show clear cases of O2 working?"
- "Can we try to find families with 95-99% similarity members, those should have tied AS reads"
- "lets stop this and focus on building on those families that are more likely to work"
- On the design: "yes".

## Why

**Where ties live.** On HG002's reads, AS-tied reads are 0.0-0.1% of reads in families whose closest copies are 95-99%
identical, 3.5% at 99.5-99.9%, and 28% at 99.9% or more.

**Why O2 placed so few.** The donor-truth test (register 1329) found that most tied reads carry the donor's own alleles,
which both CHM13 copies lack. O2's PSVs, taken from CHM13, cannot place them.

**What `samplecopy.py` does instead.** It clusters each family's reads by the alleles they share, read-backed phasing over
the sample's own copies, after SDA / Chaisson 2017 / IsoCon in the thesis vault. A cluster takes the copy its untied reads
resolve to, by a confident majority. A tied read takes its cluster's copy when it carries two decisive sites; otherwise it
abstains.

## Development (spent): TBC1D3B and PDPR

- **PDPR:** 14 right, 0 wrong of 60 tied reads with a donor copy. v14 places none.
- **TBC1D3B:** abstains on all 97. Its reads hold no within-sample variation that separates the copies: only 4 informative
  columns, and the donor separates its tied reads at 2 sites that every one of its reads shares.
- **Rules fixed in development** (each was a development finding):
  - The error rate is the reads' background mismatch rate to the anchor.
  - The allele floor is half of one copy-haplotype's expected share.
  - A tied read needs two decisive sites.
  - Labels need a confident majority: the 95% Wilson lower bound of the share above 1/2.

## The held-out families (by a fixed rule, before any held-out run)

Every HG002 family with at least 5 AS-tied reads that have a donor copy (`hg002proof/truth.tsv`: `unique` or
`unique_elsewhere`), not pseudoautosomal, with no chrM copy, and other than the two development families.

**19 families, 399 tied reads with a donor copy, 7,228 donor-tied (`own_tie`):** AC0064533 14, AC1275021 6, ANKRD36C 29,
DMRTC1 5, DUSP22 7, EIF3C 24, GOLGA8A 9, GTF2H2C 5, GTF2IRD2 35, HSPA6 34, IKBKG 47, MSTRG3150 60, NBPF1 27, SERF1B 13,
SLC35E2B 18, SMG1 5, SMN2 22, SPDYE1 34, SULT1A1 5.

## The run

```
H=/mnt/linuxdisk/tmp/o3_graph/hg002proof; O=/mnt/linuxdisk/tmp/o3_graph/samplecopy_heldout
for F in <the 19 families>; do
  python3 bench/o3_graph/samplecopy.py run /mnt/linuxdisk/tmp/o3_graph/r19seq/families19cal_full2_all.json seq19 $F $H/chm13.bam \
    $H/tied.tsv $H/untied.tsv /mnt/linuxdisk/home/juanfraitu/npip_ladder/idx/target.fa $O
done
python3 bench/o3_graph/samplecopy.py combine $H/truth.tsv $O /mnt/linuxdisk/tmp/o3_graph/fast/hg002/nonpar/INTACT <families> $O/combine.json
```

- **Truth:** the donor's own haplotypes (register 1329, `hg002truth.py`).
- **v14:** the HG002 world's `o2.integv.m2.tsv`.
- **The combined system:** v14's placement stands, and the sample copies fill v14's abstentions.

## Bars (pooled over the 19 families)

- **S1, precise:** the sample copies' wrong / placed <= 0.02 on tied reads with a donor copy. Floor: 50 placed; below it,
  UNDERPOWERED.
- **S2, calibrated:** the sample copies place <= 0.05 of the `own_tie` reads.
- **S3, adds without harm:** v14 + sample copies has more right placements than v14 and no more wrong ones.

**Verdict:**
- **PASS:** S1 is powered and holds, and S2 and S3 hold.
- **FAIL:** a powered bar fails.
- **UNDERPOWERED:** otherwise.

**Reported:** per family, recall (right / tied reads with a donor copy), and the clusters (columns, sizes, labels).

**Rules:** the code is frozen at the commit that adds this file. A bug is recorded as an amendment before any rerun.

## What a PASS would and would not claim

**A PASS would claim:** on real reads of 19 held-out families, clustering the sample's reads by their own alleles places
AS-tied reads that CHM13's PSVs cannot place, at >= 98% precision. It abstains on donor-tied reads, and it adds to v14
without errors.

**It would not claim:** testis-only families (the YAGs, absent from this library), power where the sample's reads carry no
separating variation (TBC1D3B), or other donors.

## Amendment 1 (2026-10-19, 13:05, during the run and before any score was computed)

**The bug.** The scorer counted only the tied reads present in `<family>.calls.tsv`. A tied read that never aligned to the
anchor copy has no row there and so dropped out of every count. MSTRG3150, for example, has 94 tied reads and 30 rows.

**The fix.** `samplecopy.py combine ... --tied hg002proof/tied.tsv` scores every tied read of the family that has a truth
row. A read with no call counts as an abstention, for the sample copies and in v14 + sample copies alike.

**What does not change.** No call. S1 (wrong / placed) is unaffected. S2's and S3's denominators now cover all the
family's tied reads. The pre-registered score is the `--tied` one.

## Results (2026-10-19, run once; `/mnt/linuxdisk/tmp/o3_graph/samplecopy_heldout/`, scored with Amendment 1's `--tied`)

**Verdict: FAIL.** S1 and S3 are powered and fail; S2 holds.
- **S1:** the sample copies are wrong on 4 of 86 placements (0.047) against the 0.02 bar.
- **S2:** they place 7 of the 7,228 `own_tie` reads (0.001).
- **S3:** v14 + sample copies makes 126 right placements against v14's 48, but 11 wrong against 7.

| pooled, 19 families (399 tied reads with a donor copy) | right | wrong | abstain | own_tie placed |
|---|---|---|---|---|
| sample copies | 82 | 4 | 313 | 7 |
| v14 | 48 | 7 | 344 | 52 |
| v14 + sample copies | 126 | 11 | 262 | 59 |

**Clear cases** (sample copies against v14, on the same tied reads):

| family | sample copies right / wrong | v14 right / wrong |
|---|---|---|
| HSPA6 | **34 / 0** | 0 / 0 |
| GTF2IRD2 | **31 / 0** | 0 / 0 |
| MSTRG3150 | 7 / 0 | 0 / 0 |
| GTF2H2C | 4 / 0 | 3 / 0 |
| SMG1 | 2 / 0 | 0 / 0 |
| AC1275021 | 2 / 0 | 0 / 2 |

The development family PDPR adds 14 / 0 against v14's 0 / 0.

**Other families:**
- **The 4 errors:** GOLGA8A 3 and AC0064533 1.
- **Abstained on every tied read:** IKBKG (no informative column; v14 is right on 42 by structure), SMN2, SPDYE1, SERF1B,
  SULT1A1, NBPF1, DUSP22, SLC35E2B, ANKRD36C, DMRTC1 and EIF3C.

**The errors, read post hoc.** Each sits in a cluster labelled by a majority that is not pure:
- **GOLGA8A:** 829 of 907 untied reads (91%) and 152 of 162 (94%).
- **AC0064533:** 35 of 37.
- **Their donor margins:** 5-15, i.e. one to three sites.

A tied read assigned by its cluster inherits that cluster's impurity. The majority test (Wilson lower bound above 1/2) admits
clusters that are 6-9% another copy. A purity test is the post hoc candidate. It would need its own pre-registration on
fresh families.
