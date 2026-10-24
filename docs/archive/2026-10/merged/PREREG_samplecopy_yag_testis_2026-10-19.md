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

## Results (2026-10-19, run as frozen; `yag_testis/run.log`, `out/*.calls.tsv`, `*.split.tsv`, `*.clusters.json`)

**Y1 holds: 187 of 187 tied reads called by both half-models agree (1.000), above the floor of 50.** The claim the bar allows
stands: the YAG tied reads that are placed are placed reproducibly. Correctness is not claimed, and the post-hoc reading
below gives a reason to doubt that these placements are provenance.

| family | reads / informative columns | clusters (labelled) | tied reads | full model placed | A / B / both / agree |
|---|---|---|---|---|---|
| HSFY1 | 2,191 / 6 | 1 (1) | 1,797 | 0 | 0 / 0 / 0 / 0 |
| RBMY1A1 | 942 / 101 | 7 (5) | 589 * | 134 | 163 / 218 / 134 / 134 |
| DAZ1 | 1,053 / 326 | 3 (1) | 579 | 15 | 0 / 19 / 0 / 0 |
| TSPY2 | 263 / 275 | 12 (6) | 136 | 0 | 0 / 0 / 0 / 0 |
| CDY2B | 207 / 38 | 2 (1) | 105 | 59 | 59 / 53 / 53 / 53 |
| BPY2 | 57 / 0 | 1 (0) | 47 | 0 | 0 / 0 / 0 / 0 |
| PRY | 61 / 214 | 2 (0) | 10 | 0 | 0 / 0 / 0 / 0 |

\* 590 in the table above; one read has only 73 bp inside the copy window (its other 1,442 bp lie in a downstream exon)
and does not reach the model: an abstention.

**The agreed calls, by copy:** RBMY1F@chrY:23,239k 132, RBMY1HP 2, CDY1@chrY:26,434k 53. DAZ1-4, HSFY1/2, TSPY, BPY2 and PRY: none.

**What each family's reads show** (untied reads = the copy minimap2's unique best AS picks):
- **HSFY1/2:** every read ties on the pair (1,797 of 1,827 reads), and the reads have 6 variable columns and form one group.
  The two copies are the same sequence over what is expressed. Nothing can be assigned, and nothing is.
- **DAZ1-4** (DAZ1 / DAZ4 and DAZ2 / DAZ3 are inverted pairs): 393 of the 579 tied reads tie on all four copies. The reads
  form three groups:
  - a DAZ1 group, labelled (untied DAZ1 79, DAZ4 6, DAZ2 3, DAZ3 2);
  - a DAZ2 + DAZ3 group, unlabelled (untied DAZ3 78, DAZ2 76): the two copies are not separable by this sample's reads;
  - a mixed group, unlabelled (DAZ1 32, DAZ3 27, DAZ4 18, DAZ2 11).
  The full model places 15 tied reads on DAZ1, and no read is placed by both halves. **DAZ is not resolved.**
- **RBMY1:** 7 groups. Five are labelled: RBMY1HP, RBMY1E, RBMY1A1, and **two groups (115 and 78 reads) both labelled
  RBMY1F@23,239k**. Two groups are unlabelled (73 and 44 reads). All 132 RBMY1F calls are reads tied between the two RBMY1F
  copies (23,004k on the minus strand, 23,239k on the plus strand). Of the family's untied reads, 62 resolve to 23,239k and
  8 to 23,004k.
- **CDY:** two groups.
  - A CDY1 group, labelled 26,434k by 76 of its 84 untied reads; the other 8 resolve to the other CDY1 copy (24,814k).
    All 53 agreed calls are CDY1-pair tied reads placed on 26,434k.
  - A CDY2A + CDY2B group of 47 reads, none untied (the pair is one sequence over what is expressed), unlabelled. Its 46
    tied reads abstain.
- **TSPY** (48 copies, a tandem array): 12 groups, 6 labelled, and no tied read meets the placement bar. **BPY2** (no
  variable column) and **PRY** (no labelled group) abstain.

**Post-hoc reading (not a bar; it does not change Y1).** The Y is haploid, so one copy gives one group of reads. In RBMY1,
two groups carry one label: the donor has two RBMY1F-like sequences that both look like 23,239k at the sites where CHM13's
two RBMY1F copies differ. Each read's call then follows its group's label, and nothing in the read tells the two CHM13
copies apart. Two explanations fit:
- an arm-to-arm conversion that left the donor's 23,004k with 23,239k's alleles, in which case about one group's reads are
  placed on the wrong copy;
- a donor duplication of RBMY1F.

CDY1 shows the converse: one group for two copies, labelled by a majority of 76 to 8. Its calls are the majority copy
for every CDY1 tied read, which is right only if the donor's other CDY1 is silent or absent. (Losing one CDY1 copy is
common: the gr/gr deletion removes one.)

Split-half agreement cannot see either problem, because both halves learn the same groups and the same labels. **The 187
agreed calls are reproducible, not provenance.** No YAG call is offered as a clear case of O2 working. The clear cases
remain those with donor truth:
- HSPA6, GTF2IRD2 and PDPR (register 1333);
- DHRSX (register 1335).

**O3 flags (groups no CHM13 copy labels):** RBMY1 has two unlabelled groups plus the duplicate RBMY1F label; DAZ has the
DAZ2 + DAZ3 group and the mixed group; CDY has the CDY2A + CDY2B group; PRY has two groups. These are candidates for donor
copies or alleles that CHM13 lacks. On these RNA reads alone, a group that is an extra copy cannot be told from one that
is a copy pair identical over what is expressed.
