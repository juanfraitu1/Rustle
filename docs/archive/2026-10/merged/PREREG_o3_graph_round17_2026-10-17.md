# Round 17: integrated v13 on name-free families (linked by sequence identity and gene structure), three seeds per profile, and the cost of full families in CPU seconds (2026-10-17)

Pre-registration, written before any world of round 17's families was built: no read of these families has been simulated.
Code: branch `o3-fast`, frozen at the commit that adds this file.

**Why.**
- **Round 16 failed only on B1** (register 1318): in 3 of 36 worlds, 4-7 correct reads short of the current method, with no
  error where the current method has none.
- **v13 fixes the mechanisms on 52 spent runs** (register 1319), round 16's included:
  - thin and tied path calls are decided again on structure-aware distances;
  - reads with no path alignment are realigned.
- **The families are name-free** (user, 2026-10-17: "we should go by families linked by sequence identity, also if possible
  by intron/exon structure"). The advisor does not want the annotation leaned on, and the name-stem pools of rounds 11-16
  missed unnamed paralogs.

## The two systems

**Current:** copy_assign's O2 and `--psv-clusters` O3, unchanged since round 1.

**Integrated v13** = v12 (register 1317) + two O2 genome-step rules (68722272, 67ba1cd2). On 52 spent runs: 44 -> 46 runs pass,
no world newly fails, correct +4,289, errors +95 (OUTSIDE 84, BRANCH 11).
- **`--struct-rescue`:** a thin or tied path call is decided again on gapstruct path distances (a long indel at its length).
- **`--short-path-align`:** reads the path step left without any alignment are aligned again with minimap2 -s 40.

**Integrated flags (`$V13` below):**
```
--o2-distance gapcomp --member-genome-check --member-attribution --union-exact --o2-genome veto --per-read-detectability
--twins --error-floor --o1-hosts --branch-genome-novel --mosaic-shared --haplotype-check --haplotype-check-genome
--branch-linear --aln-cache --veto-abstain --locus-pairwise --o2-significance --structure-check --struct-rescue
--short-path-align --sim-het-draw --heavy-current --step-timeout 1800 --het 0 --jobs 2 --threads 2
```

## Families (held out; name-free; drawn by a rule fixed before the build)

**The pool** (`seqfamilies.py`, 9b985ec5; structure test amended before the draw):
- **Genes:** every CAT protein-coding gene, whatever its name, with span <= 60 kb and overlaps resolved: 12,583 genes.
- **Sequence link:** gene A's canonical transcript, spliced-aligned to gene B's locus (+- 1 kb), covers >= 0.8 of A at identity
  >= 0.90, in either direction (families.py's locus rule, genome-wide). That gives 4,738 links among 1,119 genes.
- **Structure link:** one gene's exon chain is carried by the other's genome. In that spliced alignment every intron lies within
  10 bp of one of A's junctions, and none of A's junctions is left without one, in either direction. 4,384 of the 4,738 links
  pass; known retrocopies fail (SETSIP and SET, NUDT4B and NUDT4).
- **The amendment, before the draw:** the first version compared A's alignment with B's annotated introns. That failed identical
  genes whose canonical models differ (OR4F21 and OR4F16; ATAD3A's extra exon): a choice of the annotation, not of the genomes.
- **Families:** connected components of the links that pass both tests, 244 of them; 128 are unspent (no gene in any earlier
  simulated set), 282 copies.

**The draw** (fixed in 9b985ec5, before the build): the unspent families sorted by name, `random.Random(17).shuffle`, taken in
that order until the copy total reaches 90, a family skipped if it would take the total past 150.
- **40 families, 91 copies:** 33 pairs, 6 trios and one family of 7.
- **Four are X/Y pseudoautosomal pairs** (ASMTL, PPP2R3B, SHOX, SLC25A6), identical by construction.
- **File:** `/mnt/linuxdisk/tmp/o3_graph/families17.json` (set `r17`, sha1 b869562e5628), and `r17seq/families17seq_all.json`
  (+ `.links.tsv`).

| family | copies | genes | chromosomes |
|---|---|---|---|
| S17ADH1A | 2 | ADH1A, ADH1B | chr4 |
| S17AL6272301 | 7 | AL391987.2, AL627230.1 | chr9, chr13, chr21 |
| S17ASMTL | 2 | ASMTL | chrX, chrY |
| S17CBWD2 | 2 | CBWD1, CBWD2 | chr2, chr9 |
| S17CXorf49 | 2 | CXorf49, CXorf49B | chrX |
| S17DEFB104A | 3 | DEFB104A, DEFB104B | chr8 |
| S17DEFB107A | 3 | DEFB107A, DEFB107B | chr8 |
| S17DEFB4B | 3 | DEFB4A, DEFB4B | chr8 |
| S17DHRS4 | 2 | DHRS4, DHRS4L2 | chr14 |
| S17FCGR3A | 3 | FCGR3A, FCGR3B | chr1 |
| S17H2BC21 | 2 | H2AC20, H2BC21 | chr1 |
| S17H4C14 | 2 | H4C14, H4C15 | chr1 |
| S17KRT6B | 3 | KRT6A, KRT6B, KRT6C | chr12 |
| S17LCE3E | 2 | LCE3D, LCE3E | chr1 |
| S17LGALS7 | 2 | LGALS7, LGALS7B | chr19 |
| S17LIMS3 | 2 | LIMS3, LIMS4 | chr2 |
| S17LYZL1 | 2 | LYZL1, LYZL2 | chr10 |
| S17METTL2A | 2 | METTL2A, METTL2B | chr7, chr17 |
| S17MT1E | 2 | MT1E, MT1M | chr16 |
| S17NBPF6 | 2 | NBPF4, NBPF6 | chr1 |
| S17OR1L4 | 2 | OR1L4, OR1L6 | chr9 |
| S17OR1S2 | 2 | OR1S1, OR1S2 | chr11 |
| S17OR2T3 | 2 | OR2T3, OR2T34 | chr1 |
| S17OR52E6 | 2 | OR52E6, OR52E8 | chr11 |
| S17OR52K2 | 2 | OR52K1, OR52K2 | chr11 |
| S17OR5H1 | 2 | OR5H1, OR5H15 | chr3 |
| S17OR6B2 | 2 | OR6B2, OR6B3 | chr2 |
| S17OR9G1 | 2 | OR9G1 | chr11 |
| S17PPP2R3B | 2 | PPP2R3B | chrX, chrY |
| S17PRAMEF7 | 2 | PRAMEF7, PRAMEF8 | chr1 |
| S17RFPL4A | 2 | RFPL4A, RFPL4AL1 | chr19 |
| S17SFTPA2 | 2 | SFTPA1, SFTPA2 | chr10 |
| S17SHOX | 2 | SHOX | chrX, chrY |
| S17SLC25A6 | 2 | SLC25A6 | chrX, chrY |
| S17SMIM10L2B | 2 | SMIM10L2A, SMIM10L2B | chrX |
| S17SPACA5 | 2 | SPACA5, SPACA5B | chrX |
| S17SPATA31D4 | 3 | SPATA31D1, SPATA31D3, SPATA31D4 | chr9 |
| S17SPRR1A | 2 | SPRR1A, SPRR1B | chr1 |
| S17ULBP2 | 2 | RAET1L, ULBP2 | chr6 |
| S17XCL2 | 2 | XCL1, XCL2 | chr1 |

## Worlds and runs

- **Worlds:** six world types (INTACT, CROSS, OUTSIDE, REMOVED, BRANCH, SILENT) under `profiles/A119b_testis_v2.json` and
  `profiles/public_testis_v2.json`.
- **Seeds:** 701, 702 and 703, fixed with no seed selection: 6 cells, 36 worlds.

```
F=/mnt/linuxdisk/tmp/o3_graph/families17.json; O=/mnt/linuxdisk/tmp/o3_graph/fast/heldout17
for S in 701 702 703; do for P in A119b public; do
  python3 bench/o3_graph/fastset.py run $O/${P}_s$S $F r17 $S --local $O/local $V13 --sim-profile bench/o3_graph/profiles/${P}_testis_v2.json
done; done
```

## Bars

**Correctness (as rounds 15-16).** Read outcomes are per-read detectability less each world's error floor, scored on the veto
table (`--o2 o2.integv.m2.tsv`). Bars per world:
- **B1:** integrated correct >= current correct.
- **B2:** integrated errors <= current errors.

**A cell PASSES when all six of its worlds pass both bars. Round 17 PASSES only when all six cells pass** (72 bar checks).

**Cost (C1), in CPU seconds, as round 16.** Measured on the mainstay set (NPIP as one family, TBC1D3, GOLGA8, GOLGA6L, CT47A,
XAGE, SPATA31; seed 306; both profiles), nothing else running:
```
M=/mnt/linuxdisk/tmp/o3_graph/fast/r14dev/mainstay_npip306/local
for P in A119b public; do
  python3 bench/o3_graph/fastset.py run /mnt/linuxdisk/tmp/o3_graph/fast/r17cost/v13/$P \
    /mnt/linuxdisk/tmp/o3_graph/families_mainstay_npip.json mainstay 306 --local $M $V13 \
    --sim-profile bench/o3_graph/profiles/${P}_testis_v2.json
done
```
- **Bar:** v13's integrated CPU seconds (O1 + O2 + O2_genome + O3, summed over the six worlds) at most the current method's, in
  each profile.
- **Reported, no bar:** wall times, peak RSS, and the CPU seconds of round 17's cells.

**Rules:**
- Nothing is changed after the first run.
- A bug is recorded as an amendment, with its reason, before any rerun.
- A run that crashes is rerun once, and the crash is recorded.

## What a PASS would and would not claim

**A correctness PASS** would claim that v13 does at least as well as the current method in every world type, on 40 held-out
families chosen by sequence and structure rather than by name, under both profiles and three seeds.

**A cost PASS** would claim that, on NPIP-sized families, the integrated pipeline uses no more CPU than the current method.

**Neither would claim:**
- anything on real reads;
- that near-identical uncatalogued twins are handled (register 1314-1315);
- anything about arrays of FAM90A's or DUX4's size beyond the spent runs.

The pipeline under test still starts from the annotation's copies and isoforms (guided O1); round 17 removes the gene names from
the choice of test families, not the annotation from the method.

## Results (2026-10-17 23:34 - 2026-10-18 00:14; nothing changed after the first run, no crash)

**Verdict: FAIL, 0 of 6 cells. C1 PASSES** in both profiles.

| cell | failing worlds (integrated vs current: correct / errors) |
|---|---|
| A119b 701 | OUTSIDE 3,727 / 668 vs 893 / 191 |
| public 701 | OUTSIDE 3,671 / 781 vs 931 / 173 |
| A119b 702 | OUTSIDE 3,690 / 367 vs 725 / 163 |
| public 702 | OUTSIDE 3,487 / 376 vs 769 / 120; REMOVED 2,988 / 2 vs 803 / 1 |
| A119b 703 | OUTSIDE 3,792 / 540 vs 531 / 109 |
| public 703 | OUTSIDE 3,641 / 549 vs 530 / 112 |

- **B1 holds in all 36 worlds.** Round 16's shortfall is gone: v13 makes 124,264 correct calls against the current method's
  82,008.
- **Every failing bar but one is OUTSIDE's B2.** Over the 36 worlds v13 makes fewer errors, 6,986 against 8,489, but in OUTSIDE
  it makes 3-5 times the current method's.
- **The mechanism, post hoc (A119b 701 OUTSIDE):** 639 of the 668 errors are plain path calls on a family whose uncatalogued
  twin O1 did not seed as a new copy, so the catalogue keeps only the other copy:
  - LGALS7: 223 errors; CXorf49: 197; SPACA5: 27; SLC25A6: 18. The current method makes none in these families.
  - The twin's reads go to the only catalogued copy as a sole candidate. The current method's genomic candidates tie the read
    with the twin's locus, and it abstains.
- **This is the O2 proof's open failure (register 1314-1315),** made dominant by a set of pairs: take one copy of a pair out of
  the catalogue, and its near-identical twin has no rival in O2.
- **The guard measured for it in register 1315** (abstain when a locus off every seeded copy explains the read as well) cut these
  errors to 0-2. It was rejected because it also abstained the reads of an expressed copy tied with a silent catalogued copy
  that O1 had dropped.

**C1, the cost on the mainstay set:** v13 used 830 CPU seconds (A119b) and 675 (public) against the current method's 1,042 and
952 (-20% and -29%). Wall time, reported: 469 s and 378 s against 313 s and 281 s, at the flags' threads. Round 17's cells
used 591 CPU seconds in all against 1,033.
