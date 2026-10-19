# Pre-registration: the haploid rule of the O3 -> O2 sample-copy step, held out on HG002's chrX families (2026-10-19)

**User (2026-10-19):** "yes lets implement the rule", after register 1336. On A119b testis (male), two RBMY1 clusters were
labelled with one RBMY1F copy, and CDY1's two copies formed one cluster. On a haploid chromosome that cannot be one copy
per cluster, yet every agreed YAG call followed such a label.

## The rule (code 2adc91d7, `samplecopy.py run ... --male`)

A male sample carries its chrX and chrY copies once outside the pseudoautosomal regions (`hg002truth.PAR`), so each of
those copies gives one group of reads. A read tied only on such copies is placed only when each copy it ties on labels
**exactly one** cluster (`one_copy_one_group`):
- **two or more clusters** on one copy: at least one of them is not that copy;
- **no cluster**: the copy's reads may sit unseen inside another copy's cluster.

Reads tied on any diploid copy are untouched, so autosomal and pseudoautosomal families are unchanged. That was checked
byte-identical on PDPR, and on CDY2B without `--male`. The flags name, for O3, the haploid copies that label several
clusters or none.

**Development (A119b YAGs, spent):** the rule blocks all 208 placements (RBMY1 134, CDY 59, DAZ 15); split-half then calls
nothing. On HG002's two haploid families that have donor-copy reads (IKBKG 47, DMRTC1 5; spent, register 1333), the step
placed none, so the rule has nothing to remove there.

## The held-out families (by a fixed rule, never run by the step)

Every HG002 family whose copies are all haploid in a male, with at least 5 tied reads in `hg002proof/tied.tsv`, other
than the spent IKBKG and DMRTC1. Tied = equal best AS, the tables of registers 1329-1336.

| family | copies | tied | own_tie | elsewhere | poor_match | untied |
|---|---|---|---|---|---|---|
| F8A2 | F8A2, F8A3 (chrX) | 2,335 | 1,043 | 918 | 374 | 5 |
| FAM156B | FAM156B, FAM156A (chrX) | 634 | 532 | 0 | 102 | 36 |
| MAGED4B | MAGED4B, MAGED4 (chrX) | 50 | 39 | 0 | 11 | 2 |
| H2AB2 | H2AB2, H2AB3 (chrX) | 6 | 4 | 2 | 0 | 3 |

**No tied read here has a donor copy in its family** (no `unique`):
- own_tie reads tie in the donor too, so placing one is unsupported;
- elsewhere reads come from a donor copy off the family, so placing one is wrong;
- poor_match reads match the donor poorly and are reported but not scored.

The scored denominator is **2,538** (1,618 own_tie + 920 elsewhere).

## The run

```
H=/mnt/linuxdisk/tmp/o3_graph/hg002proof; O=/mnt/linuxdisk/tmp/o3_graph/samplecopy_haploid
for F in S19F8A2 S19FAM156B S19MAGED4B S19H2AB2; do
  python3 bench/o3_graph/samplecopy.py run /mnt/linuxdisk/tmp/o3_graph/r19seq/families19cal_full2_all.json seq19 $F $H/chm13.bam \
    $H/tied.tsv $H/untied.tsv /mnt/linuxdisk/home/juanfraitu/npip_ladder/idx/target.fa $O/plain
  python3 bench/o3_graph/samplecopy.py run <same arguments, output $O/male> --male
done
python3 bench/o3_graph/samplecopy.py score $H/truth.tsv $O/plain <families>; ... score $H/truth.tsv $O/male <families>
```

## Bars

- **HS, calibrated on haploid families (with the rule):** own_tie placed + elsewhere placed <= 0.05 x 2,538, so at most
  126 reads.
- **HB, the rule's benefit:** the step without the rule places at least 50 own_tie + elsewhere reads (the floor), and
  the rule places fewer.

**Verdict:**
- **PASS:** HS and HB hold.
- **UNDERPOWERED:** HS holds, but the step without the rule places fewer than 50 such reads, so there was little to remove.
- **FAIL:** HS fails.

**Reported:** per family, placements by class with and without the rule, the clusters (labels, untied counts) and the
O3 flags.

**Rules:** code frozen at 2adc91d7. A bug is recorded as an amendment before any rerun.

## What a PASS would and would not claim

**A PASS would claim:** on real reads where the donor itself cannot place them, the haploid rule stops the sample-copy
step from placing reads that the step without it places.

**It would not claim:**
- that the rule keeps right placements: these families have none to keep. The no-harm side rests on construction, since
  diploid families are untouched;
- anything about the YAGs, which have no truth;
- anything about female samples.
