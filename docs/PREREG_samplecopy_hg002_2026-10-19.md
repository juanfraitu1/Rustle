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
