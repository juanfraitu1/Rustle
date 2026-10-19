# Pre-registration: the O3 -> O2 sample-copy step on the Y-ampliconic families in testis, by split-half agreement (2026-10-19)

**User (2026-10-19):** "ok lets do the next step". That step is DAZ and the other YAGs in testis, checked by split-half,
since no donor assembly exists.

## Why split-half, and what HG002 says it means

There is no truth for a testis donor's Y. The check is reproducibility, using `samplecopy.py splithalf`:
- the family's reads are split into two disjoint halves by a seeded hash (seed 1);
- each half learns its own clusters and labels;
- every AS-tied read is called by both models.

**Calibration on HG002** (spent families, the donor's haplotypes as truth; `fast/../samplecopy_split_cal/`), tied reads with
a donor copy:
- **Both halves call the same copy:** 519 right, 0 wrong.
- **Only one half calls:** 38 right, 3 wrong.
- **The halves disagree:** never.

So agreement is the calibrated signal. The calibration comes from HG002's families, not from YAGs, and that limit is stated
below.

## Data

- **Reads:** A119b testis (`winloci_data/A119b.t2t.bam`, CHM13 v2.0, minimap2 splice:hq `-N 50 -p 0.1`). They are the
  4,953 reads whose primary lies within 1 kb of a copy of the 10 Y-ampliconic families of the adopted pool, with all their
  chrX / chrY alignments (`yag_testis/yag.bam`).
- **Limit:** secondaries off chrX / chrY are not seen. DAZL and the TSPY-like genes are autosomal but divergent.
- **Tables:** `hg002truth.py tied` / `untied`. 3,270 reads are AS-tied within one family and 1,021 are resolved by minimap2.

| family | copies | tied | untied |
|---|---|---|---|
| HSFY1 | 2 | 1,797 | 30 |
| RBMY1A1 | 6 | 590 | 293 |
| DAZ1 | 4 | 579 | 346 |
| TSPY2 | 48 | 136 | 123 |
| CDY2B | 4 | 105 | 84 |
| BPY2 | 3 | 47 | 5 |
| PRY | 4 | 10 | 11 |
| XKRY / VCX3A / RBMY2BP | 2 / 6 / 3 | 3 / 2 / 1 | 0 / 123 / 4 |

## The run (the frozen original rule; no change since be1edda6)

```
Y=/mnt/linuxdisk/tmp/o3_graph/yag_testis
for F in S19HSFY1 S19RBMY1A1 S19DAZ1 S19TSPY2 S19CDY2B S19BPY2 S19PRY; do
  python3 bench/o3_graph/samplecopy.py run /mnt/linuxdisk/tmp/o3_graph/r19seq/families19cal_full2_all.json seq19 $F $Y/yag.bam \
    $Y/tied.tsv $Y/untied.tsv /mnt/linuxdisk/home/juanfraitu/npip_ladder/idx/target.fa $Y/out
  python3 bench/o3_graph/samplecopy.py splithalf <same arguments>
done
```

The families with at least 10 tied reads are run. XKRY, VCX3A and RBMY2BP have 1-3 each and are left out.

## Bar and reports

- **Y1, reproducible:** pooled over the families, the tied reads both half-models call agree on >= 0.95 of them. Floor: 50
  such reads; below it, UNDERPOWERED.
- **Reported:**
  - per family, the full model's placements, both-halves calls and agreements, and the agreed calls per copy (DAZ1-4, ...);
  - the clusters (labels, unlabelled clusters as candidate copies CHM13 lacks: O3's flag).
- **The claim (if Y1 holds):** testis tied reads of the YAGs are placed on CHM13 copies reproducibly from the sample's own
  variants.
- **Not claimed:** correctness. The HG002 calibration suggests agreed calls are rarely wrong, but YAG copies, testis
  expression and a donor whose Y is not CHM13's are untested by it.
