# Round 19: integrated v14 (unchanged) on fresh families whose non-coding members are copies like any other, expressed as real testis expresses them (2026-10-18)

Pre-registration, written before any world of round 19's families was built: no read of these families has been simulated.
Code: branch `o3-fast`, frozen at the commit that adds this file.

**Why.**
- **Round 18 failed the two seed-803 cells on B1** (register 1322; 4 of 6 cells passed). Nearly every lost read was OR4F21's: O1's
  genome search found five OR4F unprocessed pseudogenes, the simulator never expressed them, and v14 abstained on the reads
  tied with them. The current method compares catalogued copies only, so in that simulation it was always right.
- **Real A119b testis expresses those pseudogenes** (41 primary reads against 59 on the two coding copies).
- **The user's decision (2026-10-18):** non-coding paralogs are family members like any other. They go into the catalogue for
  both methods, and the simulator expresses them at the level the real libraries show. Retrocopies the family rule admits (those
  of single-exon genes, where structure cannot decide) are kept.
- **Development on round 18's families with their non-coding members** (register 1323, spent): 6 of 6 cells pass. The thinnest
  B1 margin is +93 reads (round 18: -107) and the thinnest B2 margin +7.

## The two systems

**Current:** copy_assign's O2 and `--psv-clusters` O3, unchanged since round 1.

**Integrated v14:** unchanged since round 18 (register 1321), with the same flags:
```
--o2-distance gapcomp --member-genome-check --member-attribution --union-exact --o2-genome veto --per-read-detectability
--twins --error-floor --o1-hosts --branch-genome-novel --mosaic-shared --haplotype-check --haplotype-check-genome
--branch-linear --aln-cache --veto-abstain --locus-pairwise --o2-significance --structure-check --struct-rescue
--short-path-align --ref-copies --sim-het-draw --heavy-current --step-timeout 1800 --het 0 --jobs 2 --threads 2
```

**Only the simulation changed:**
- **Families:** `seqfamilies.py build --paralogs` (d8ea648d). A gene of any biotype joins a family when its locus carries a coding
  gene's canonical transcript by the family rule (identity >= 0.90 over >= 0.8, exon chain carried; that direction only).
  - A retrocopy of a multi-exon gene fails the structure test.
  - None may overlap a coding gene, and a pseudogene wins over a lncRNA that overlaps it.
- **Expression:** profiles v4 (19a1b482) = v2 plus the measured expression of non-coding members (`calibrate.py paralogs`,
  uniquely placed reads, families whose coding copies average >= 10 reads):
  - A119b: 3.7% silent, the rest expressed at a median 0.34 of the family's coding copies;
  - public testis: 37.5% silent, the rest at 0.14.
  - A silent member stays catalogued (kind `silent`).
- **Scenario copies** (REMOVED, OUTSIDE, BRANCH) and SILENT's silenced third are drawn among the coding copies, as before
  (d8ea648d; BRANCH's extra copy fixed in b2502ca4).

## Families (held out; drawn by the fixed rule)

**The pool:** `seqfamilies.py build --paralogs` with every earlier simulated set spent, round 18's included (23 files,
`ls families*.json | grep -v _all.json`).
- 341 families by sequence and structure; 130 are unspent (301 copies, 106 of them non-coding).
- **File:** `r19seq/families19seq_all.json`.

**The draw:** the rule of rounds 17-18 (51f9b3a6) with the round's number as seed: the unspent families sorted by name,
`random.Random(19).shuffle`, taken in that order until the copy total reaches 90, a family skipped if it would pass 150.
- **35 families, 90 copies** (54 coding, 36 non-coding); copies per family {2: 23, 3: 7, 4: 3, 5: 1, 6: 1}.
- **Non-coding members:** 19 unprocessed pseudogenes, 4 transcribed unprocessed, 8 processed, 3 transcribed processed, 1
  polymorphic, 1 unitary.
  - 9 of the 11 retrocopies are in families of single-exon genes.
- **File:** `/mnt/linuxdisk/tmp/o3_graph/families19.json` (set `r19`, sha1 73ee6bdb294c).

| family | copies | coding genes | non-coding members | chromosomes |
|---|---|---|---|---|
| S19ADAM21 | 2 | ADAM21 | ADAM21P1 (tr. proc.) | chr14 |
| S19ANXA8L1 | 2 | ANXA8L1, ANXA8 | - | chr10 |
| S19AZGP1 | 2 | AZGP1 | AZGP1P1 (tr. unproc.) | chr7 |
| S19DDX11 | 3 | DDX11 | AC009533.1 (unproc.), DDX12P (unproc.) | chr12 |
| S19DEFB106A | 3 | DEFB106A, DEFB106A, DEFB106B | - | chr8 |
| S19DEFB130A | 3 | DEFB130A | DEFB130C (tr. unproc.), DEFB130D (unproc.) | chr4, chr8, chr11 |
| S19GPAT2 | 2 | GPAT2 | GPAT2P1 (tr. unproc.) | chr2 |
| S19GPR32 | 2 | GPR32 | GPR32P1 (unproc.) | chr19 |
| S19H2BC9 | 4 | H2BC9, H2BC13, H2AC13, H2BC17 | - | chr6 |
| S19H4C11 | 2 | H4C11, H4C12 | - | chr6 |
| S19HSD3B2 | 4 | HSD3B2, HSD3B1 | HSD3BP2 (unproc.), HSD3BP4 (unproc.) | chr1 |
| S19IKBKG | 2 | IKBKG | IKBKGP1 (unproc.) | chrX |
| S19IL3RA | 2 | IL3RA, IL3RA | - | chrX, chrY |
| S19NACA2 | 3 | NACA2 | NACA3P (proc.), NACA4P (tr. proc.) | chr4, chr8, chr17 |
| S19NAT8 | 2 | NAT8 | NAT8B (polymorphic) | chr2 |
| S19NUDT4B | 2 | NUDT4B | NUDT4P2 (proc.) | chr1 |
| S19OR4A5 | 2 | OR4A5 | OR4A42P (unproc.) | chr11 |
| S19OR4K15 | 3 | OR4K15 | OR4K8P (unproc.), OR4K11P (unproc.) | chr14, chr18, chr21 |
| S19OR5H14 | 3 | OR5H14, OR5H6 | OR5H7P (unproc.) | chr3 |
| S19OR5M8 | 2 | OR5M8 | OR5M6P (unproc.) | chr11 |
| S19OR8I2 | 2 | OR8I2 | OR8I1P (unproc.) | chr11 |
| S19PABPC1L2B | 2 | PABPC1L2B, PABPC1L2A | - | chrX |
| S19PDZK1 | 2 | PDZK1 | PDZK1P1 (tr. unproc.) | chr1 |
| S19PGAM4 | 4 | PGAM4 | PGAM1P8 (tr. proc.), AC012085.1 (proc.), AC078819.1 (proc.) | chrX, chr11, chr12 |
| S19PRAMEF10 | 5 | PRAMEF10, PRAMEF33, PRAMEF17 | PRAMEF31P (unproc.), PRAMEF31P (unproc.) | chr1 |
| S19PRR23D1 | 6 | PRR23D1, PRR23D2, PRR23D1 | AC134684.4 (unproc.), PRR23D3P (unproc.), AC134684.4 (unproc.) | chr8 |
| S19PTTG2 | 2 | PTTG2 | PTTG3P (proc.) | chr4, chr8 |
| S19RHOXF2B | 2 | RHOXF2B, RHOXF2 | - | chrX |
| S19SAA2 | 2 | SAA2, SAA1 | - | chr11 |
| S19SERF1B | 2 | SERF1B, SERF1A | - | chr5 |
| S19SLC25A52 | 3 | SLC25A52 | SLC25A51P1 (proc.), SLC25A51P3 (proc.) | chr6, chr8, chr18 |
| S19TGIF2LX | 2 | TGIF2LX, TGIF2LY | - | chrX, chrY |
| S19VN1R4 | 2 | VN1R4 | VN1R6P (unitary) | chr19 |
| S19ZC3H11B | 2 | ZC3H11B | AL365217.1 (proc.) | chr1, chr6 |
| S19ZNF181 | 2 | ZNF181 | ZNF807P (unproc.) | chr19 |

## Worlds and runs

- **Worlds:** six world types (INTACT, CROSS, OUTSIDE, REMOVED, BRANCH, SILENT) under `profiles/A119b_testis_v4.json` and
  `profiles/public_testis_v4.json`.
- **Seeds:** 901, 902 and 903, fixed with no seed selection: 6 cells, 36 worlds.

```
F=/mnt/linuxdisk/tmp/o3_graph/families19.json; O=/mnt/linuxdisk/tmp/o3_graph/fast/heldout19
for S in 901 902 903; do for P in A119b public; do
  python3 bench/o3_graph/fastset.py run $O/${P}_s$S $F r19 $S --local $O/local $V14 --sim-profile bench/o3_graph/profiles/${P}_testis_v4.json
done; done
```

## Bars

**Correctness (as rounds 15-18).** Per-read detectability less each world's error floor, on the veto table
(`--o2 o2.integv.m2.tsv`). Bars per world:
- **B1:** integrated correct >= current correct.
- **B2:** integrated errors <= current errors.

**A cell PASSES when all six of its worlds pass both bars. Round 19 PASSES only when all six cells pass** (72 bar checks).

**Cost (C1), in CPU seconds, as rounds 16-18,** on the mainstay set (seed 306, both profiles), nothing else running:
```
M=/mnt/linuxdisk/tmp/o3_graph/fast/r14dev/mainstay_npip306/local
for P in A119b public; do
  python3 bench/o3_graph/fastset.py run /mnt/linuxdisk/tmp/o3_graph/fast/r19cost/v14/$P \
    /mnt/linuxdisk/tmp/o3_graph/families_mainstay_npip.json mainstay 306 --local $M $V14 \
    --sim-profile bench/o3_graph/profiles/${P}_testis_v4.json
done
```
- The mainstay families carry no biotypes, so v4 simulates them as v2 did.
- **Bar:** v14's integrated CPU seconds (O1 + O2 + O2_genome + O3, summed over the six worlds) at most the current method's, in
  each profile.
- **Reported, no bar:** wall times, peak RSS, the cells' CPU seconds, and the genome copies O1 found per world.

**Rules:**
- Nothing is changed after the first run.
- A bug is recorded as an amendment, with its reason, before any rerun.
- A run that crashes is rerun once, and the crash is recorded.

## What a PASS would and would not claim

**A correctness PASS** would claim that v14 does at least as well as the current method in every world type, on 35 held-out
name-free families with their non-coding members expressed at the measured level, under both profiles and three seeds.

**A cost PASS** would claim no more CPU than the current method on NPIP-sized families.

**Neither would claim:**
- anything on real reads;
- anything about arrays of FAM90A's or DUX4's size beyond the spent runs;
- that the family definition is settled: the Compara study (in design, 2026-10-18) may change it;
- anything about paralogs no annotation holds: loci found only by O1's genome search still express nothing.

## Results (2026-10-18 16:13-16:52; nothing changed after the first run, no crash)

**Launch note:** the first launch (16:12) failed before any run started. The batch's log directories (`fast/r19cost/`,
`fast/heldout19/`) did not exist, so every run's output redirect failed and nothing was simulated. The directories were
created and the same batch relaunched unchanged.

**Verdict: PASS -- 6 of 6 cells, all 72 bar checks. C1 PASSES** in both profiles.

| cell | INTACT | CROSS | OUTSIDE | REMOVED | BRANCH | SILENT |
|---|---|---|---|---|---|---|
| A119b 901 | 2,785 / 0 vs 2,682 / 0 | 2,828 / 0 vs 2,725 / 0 | 2,340 / 1 vs 722 / 780 | 2,375 / 0 vs 1,355 / 24 | 3,664 / 516 vs 2,743 / 1,035 | 1,461 / 0 vs 1,402 / 0 |
| public 901 | 2,359 / 0 vs 2,213 / 21 | 2,382 / 0 vs 2,237 / 21 | 1,856 / 0 vs 433 / 853 | 1,827 / 0 vs 1,116 / 7 | 3,145 / 631 vs 2,149 / 1,212 | 1,063 / 0 vs 974 / 21 |
| A119b 902 | 2,725 / 0 vs 2,652 / 0 | 2,778 / 0 vs 2,705 / 0 | 2,342 / 0 vs 726 / 534 | 2,391 / 0 vs 1,144 / 25 | 3,295 / 256 vs 2,742 / 491 | 1,336 / 0 vs 1,285 / 0 |
| public 902 | 2,275 / 0 vs 2,086 / 0 | 2,294 / 0 vs 2,105 / 0 | 1,636 / 1 vs 461 / 587 | 1,856 / 1 vs 890 / 11 | 2,746 / 267 vs 2,174 / 424 | 874 / 0 vs 836 / 0 |
| A119b 903 | 3,149 / 0 vs 3,025 / 0 | 3,213 / 0 vs 3,089 / 0 | 2,846 / 0 vs 724 / 336 | 2,758 / 1 vs 1,043 / 8 | 3,918 / 324 vs 3,025 / 694 | 1,159 / 0 vs 1,121 / 0 |
| public 903 | 2,703 / 0 vs 2,530 / 0 | 2,739 / 0 vs 2,566 / 0 | 1,927 / 6 vs 434 / 295 | 2,261 / 1 vs 715 / 11 | 3,417 / 350 vs 2,553 / 638 | 717 / 0 vs 688 / 0 |

(integrated correct / errors vs current correct / errors)

- **Over the 36 worlds:** v14 makes 85,440 correct calls against the current method's 62,070, and 2,355 errors against 8,028.
- **B1's thinnest margin is +29 reads** (public 903 SILENT), then +38 (A119b 903 SILENT, public 902 SILENT).
- **B2's thinnest margins are ties at zero:** in eight complete-catalogue worlds both methods make no error.
- **Genome copies O1 found per world:** 4-5 in INTACT, CROSS, REMOVED, BRANCH and SILENT, and 31-34 in OUTSIDE (seed 901).

**C1, the cost on the mainstay set:** v14 used 995 CPU seconds (A119b) and 825 (public), against the current method's 1,032 and
943 (-4% and -13%).
- Reported: round 19's cells used 849 CPU seconds in all against 750. On these small families O1's genome search (about 20
  CPU seconds a world) costs more than it saves.
- The set's one-time local-reference build peaks at 14.4 GB.

**What it claims, as registered:** v14 does at least as well as the current method in every world type, on 35 held-out name-free
families with their non-coding members expressed at the measured level, under both profiles and three seeds, and uses no more
CPU on NPIP-sized families. It claims nothing on real reads, on arrays of FAM90A's or DUX4's size, or on paralogs no annotation
holds.
