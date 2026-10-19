# Round 20: integrated v15 (O1's genome search under the carried intron chain, every chain kept) on fresh families built by the same rule (2026-10-18)

Pre-registration, written before any world of round 20's families was built: no read of these families has been simulated.
Code: branch `o3-fast`, frozen at the commit that adds this file.

**Why.**
- **Round 19 passed** (register 1324): v14 held every bar on fresh families whose non-coding members are copies.
- **The family rule missed known families** (register 1325):
  - GOLGA6L failed the exact structure test, because copies split a long exon with short introns of their own and indels come
    out of the spliced alignment as introns.
  - TSPY2 was lost to minimap2's default chain selection among about 45 near-identical targets.
- **The user's decisions (2026-10-18):**
  - build the next family pool with the carried intron chain (`--intron-chain carried`) and every chain kept
    (`--retain-chains`);
  - apply both to O1's genome search.
- **Development of v15 on spent runs** (register 1326):
  - Round 19's six cells and the mainstay set, both profiles: 8 of 8 pass.
  - Cold-cache cost with the search cached per set (O1_ALN_CACHE): 1,019 / 848 CPU s against the current method's 1,031 / 955.
    A119b's margin is only 12 s.

## The two systems

**Current:** copy_assign's O2 and `--psv-clusters` O3, unchanged since round 1.

**Integrated v15** = v14 + O1 `--ref-carried` (05cb9c6d, cached through O1_ALN_CACHE under `--aln-cache`, 5b3d50ac). O1's genome
search takes a locus as a copy when the catalogued transcript's intron chain is carried there (each junction has its intron; the
locus may hold more), with every minimap2 chain kept (-P).

**Integrated flags (`$V15`):**
```
--o2-distance gapcomp --member-genome-check --member-attribution --union-exact --o2-genome veto --per-read-detectability
--twins --error-floor --o1-hosts --branch-genome-novel --mosaic-shared --haplotype-check --haplotype-check-genome
--branch-linear --aln-cache --veto-abstain --locus-pairwise --o2-significance --structure-check --struct-rescue
--short-path-align --ref-copies --ref-carried --sim-het-draw --heavy-current --step-timeout 1800 --het 0 --jobs 2 --threads 2
```

**The simulation, as round 19:**
- profiles v4: non-coding members silent at the measured fraction, else expressed at the measured ratio;
- scenario copies drawn among the coding copies.
- **Fixed since round 19:** a world's regions end at their contig's end (49323d8b). Round 20's draw holds two mitochondrial
  families (chrM is 16,569 bp). No region of rounds 15-19 ran past an end.

## Families (held out; drawn by the fixed rule)

**The pool:** `seqfamilies.py build --paralogs --intron-chain carried --retain-chains`.
- Every earlier simulated set is spent, round 19's included (24 files, `ls families*.json | grep -v _all.json`).
- 341 families by sequence and structure; 96 are unspent (218 copies, 73 of them non-coding).
- **File:** `r20seq/families20seq_all.json`.

**The draw:** the rule of rounds 17-19 (51f9b3a6) with the round's number as seed: the unspent families sorted by name,
`random.Random(20).shuffle`, taken in that order until the copy total reaches 90, a family skipped if it would pass 150.
- **38 families, 90 copies** (60 coding, 30 non-coding); copies per family {2: 30, 3: 5, 4: 1, 5: 1, 6: 1}.
- **Non-coding members:** 11 unprocessed pseudogenes, 6 transcribed unprocessed, 10 processed, 1 polymorphic, 2 StringTie
  models.
- **File:** `/mnt/linuxdisk/tmp/o3_graph/families20.json` (set `r20`, sha1 fcf105530369).

| family | copies | coding genes | non-coding members | chromosomes |
|---|---|---|---|---|
| S19BX6646152 | 2 | BX664615.2 | FGF7P6 (tr. unproc.) | chr9 |
| S19CFC1B | 2 | CFC1B, CFC1 | - | chr2 |
| S19FAHD2A | 3 | FAHD2A, FAHD2B | FAHD2CP (tr. unproc.) | chr2 |
| S19FKBP9 | 2 | FKBP9 | FKBP9P1 (tr. unproc.) | chr7 |
| S19GK3P | 2 | GK3P | AC110373.1 (unproc.) | chr4 |
| S19GLYATL1B | 2 | GLYATL1B | GLYATL1P4 (unproc.) | chr11 |
| S19HCAR2 | 2 | HCAR2, HCAR3 | - | chr12 |
| S19HNRNPA1L2 | 4 | HNRNPA1L2 | HNRNPA1P7 (proc.), HNRNPA1P10 (proc.), HNRNPA1P35 (proc.) | chr2, chr13, chr18, chr19 |
| S19HP | 2 | HP, HPR | - | chr16 |
| S19HPF1 | 2 | HPF1 | AL592293.1 (unproc.) | chr4, chr9 |
| S19HSD17B7 | 2 | HSD17B7 | HSD17B7P2 (tr. unproc.) | chr1, chr10 |
| S19HSFX2 | 2 | HSFX2, HSFX1 | - | chrX |
| S19IFNA13 | 2 | IFNA13, IFNA1 | - | chr9 |
| S19IFNL3 | 2 | IFNL3, IFNL2 | - | chr19 |
| S19KRT81 | 2 | KRT81, KRT83 | - | chr12 |
| S19LAIR1 | 2 | LAIR1, LAIR2 | - | chr19 |
| S19MTCO2 | 2 | MT-CO2 | MTCO2P22 (proc.) | chr5, chrM |
| S19MTND5 | 5 | MT-ND5, MT-ND6 | MTCYBP18 (proc.), MTND6P4 (proc.), MTND5P11 (proc.) | chr5, chrM |
| S19MYH1 | 2 | MYH1, MYH2 | - | chr17 |
| S19OR10A5 | 2 | OR10A5, OR10A2 | - | chr11 |
| S19OR10H3 | 2 | OR10H3, OR10H4 | - | chr19 |
| S19OR2T27 | 2 | OR2T27 | OR2T7 (polymorphic) | chr1 |
| S19OR4C12 | 3 | OR4C12 | OR4C2P (unproc.), OR4C50P (unproc.) | chr11 |
| S19OR4C3 | 2 | OR4C3 | OR4C49P (unproc.) | chr11 |
| S19OR5AK2 | 2 | OR5AK2 | OR5AK1P (unproc.) | chr11 |
| S19OR5L1 | 2 | OR5L1, OR5L2 | - | chr11 |
| S19PDPR | 3 | PDPR | MSTRG.2242 (StringTie), MSTRG.2243 (StringTie) | chr16 |
| S19PNMA6A | 2 | PNMA6A | PNMA6B (unproc.) | chrX |
| S19RAB40AL | 2 | RAB40AL, RAB40A | - | chrX |
| S19REG3G | 2 | REG3G, REG3A | - | chr2 |
| S19RNFT1 | 2 | RNFT1 | RNFT1P3 (unproc.) | chr17 |
| S19RPL10L | 2 | RPL10L | RPL10P11 (proc.) | chr7, chr14 |
| S19SDHA | 2 | SDHA | SDHAP2 (tr. unproc.) | chr3, chr5 |
| S19TRIM74 | 3 | TRIM74, TRIM50, TRIM73 | - | chr7 |
| S19TUBA1B | 2 | TUBA1B, TUBA1A | - | chr12 |
| S19UBTFL1 | 3 | UBTFL1 | UBTFL7 (proc.), UBTFL2 (proc.) | chr11 |
| S19ZNF479 | 6 | ZNF479, ZNF722P, ZNF735, ZNF679 | ZNF733P (tr. unproc.), ZNF734P (unproc.) | chr7 |
| S19ZNF736 | 2 | ZNF736 | AC095038.1 (unproc.) | chr7 |

## Worlds and runs

- **Worlds:** six world types (INTACT, CROSS, OUTSIDE, REMOVED, BRANCH, SILENT) under `profiles/A119b_testis_v4.json` and
  `profiles/public_testis_v4.json`.
- **Seeds:** 1001, 1002 and 1003, fixed with no seed selection: 6 cells, 36 worlds.

```
F=/mnt/linuxdisk/tmp/o3_graph/families20.json; O=/mnt/linuxdisk/tmp/o3_graph/fast/heldout20
for S in 1001 1002 1003; do for P in A119b public; do
  python3 bench/o3_graph/fastset.py run $O/${P}_s$S $F r20 $S --local $O/local $V15 --sim-profile bench/o3_graph/profiles/${P}_testis_v4.json
done; done
```

## Bars

**Correctness (as rounds 15-19).** Per-read detectability less each world's error floor, on the veto table
(`--o2 o2.integv.m2.tsv`). Bars per world:
- **B1:** integrated correct >= current correct.
- **B2:** integrated errors <= current errors.

**A cell PASSES when all six of its worlds pass both bars. Round 20 PASSES only when all six cells pass** (72 bar checks).

**Cost (C1), in CPU seconds, as rounds 16-19,** on the mainstay set (seed 306, both profiles), nothing else running:
```
M=/mnt/linuxdisk/tmp/o3_graph/fast/r14dev/mainstay_npip306/local
for P in A119b public; do
  python3 bench/o3_graph/fastset.py run /mnt/linuxdisk/tmp/o3_graph/fast/r20cost/v15/$P \
    /mnt/linuxdisk/tmp/o3_graph/families_mainstay_npip.json mainstay 306 --local $M $V15 \
    --sim-profile bench/o3_graph/profiles/${P}_testis_v4.json
done
```
- **Bar:** v15's integrated CPU seconds (O1 + O2 + O2_genome + O3, summed over the six worlds) at most the current method's, in
  each profile.
- **Reported, no bar:** wall times, peak RSS, the cells' CPU seconds, and the genome copies O1 found per world.

**Rules:**
- Nothing is changed after the first run.
- A bug is recorded as an amendment, with its reason, before any rerun.
- A run that crashes is rerun once, and the crash is recorded.

## What a PASS would and would not claim

**A correctness PASS** would claim that v15 does at least as well as the current method in every world type, on 38 held-out
families built by the carried-intron-chain rule with their non-coding members expressed at the measured level, under both
profiles and three seeds.

**A cost PASS** would claim no more CPU than the current method on NPIP-sized families.

**Neither would claim:**
- anything on real reads;
- anything about arrays of FAM90A's or DUX4's size beyond the spent runs;
- anything about paralogs no annotation holds.
