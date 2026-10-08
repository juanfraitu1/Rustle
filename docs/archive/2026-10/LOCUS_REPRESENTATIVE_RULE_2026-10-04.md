# The locus representative, most-junctions (R_J) vs most-reads (R_M), on pre-f1v2 de novo GTFs: R_J carries the copy's structure, held-out families are mixed — chr6 fails clause (a), R_J stays opt-in (2026-10-04)

> ⚠ **A verdict on the configuration as run, not on the shipped one.** Both arms ran on the same pre-f1v2 GTFs (assembled 2026-09-25)
> with `RUSTLE_BRIDGE_REGROUP=off`. That setting was forced: the driver refuses f1v2 on a GTF that has no bridge relations. The chain
> behind the chr6 failure is a fused chain of the kind f1v2 may turn into a relation. Whether it would is unmeasured.

Coordinates are 1-based, closed (GTF / GFF).

Prereg: `docs/archive/2026-10/PREREG_locus_representative_rule_2026-10-04.md` (3bf6aca1, written before any run). Flag: `mcl_families --representative
most-junctions` / driver `RUSTLE_REPRESENTATIVE=most-junctions` (559649a9). Runner `bench/rep_rule/run.sh`, scorer `bench/rep_rule/score.py`;
work dir `/mnt/linuxdisk/tmp/rep_rule/<species>_<contig>/` (`R_M.*`, `R_J.*`, `h3.*`, `score.json`; `decision.json`). Per-gene tables:
`docs/LOCUS_REPRESENTATIVE_RULE_2026-10-04_h3_changed.tsv` (every gene whose H3 call differs between the arms, 7 contigs),
`docs/LOCUS_REPRESENTATIVE_RULE_2026-10-04_npip.tsv` (the 25 NPIP copies on chr16).

## Outcome

**R_M stays the default; R_J stays opt-in.** The decision is held-out only and needs all three clauses on all five held-out contigs. One
violation: **human chr6, clause (a): Compara bipartite F 0.372 under R_J vs 0.444 under R_M (bar 0.439)**, precision 1.000 in both. Every
other held-out clause passes: (a) chr2 (F 0.609 vs 0.588), chr8 (0.388 = 0.388), chr10 (0.542 = 0.542), and gorilla NC_073234.2 has no
primary reference to test (0 Liftoff pairs on the contig); (b) H3 found genes go up on all five (+16 to +70); (c) the H2 median junctions per
copy ties on all five (0, 1 or 2 under both rules).

What the run shows:
- **R_J does what it was built for, with losses.** On every contig the representative carries more of the gene's read-supported
  structure.
  - Held-out found genes: chr2 960 → 1,030, chr6 757 → 815, chr8 500 → 531, chr10 562 → 598, NC_073234.2 734 → 750. That is close to the
    locus-level ceiling (chr2 1,037).
  - On the 25 NPIP copies, strict FOUND goes from 10 to 20, equal to the locus-level reading.
  - The gain does not depend on fused representatives: dropping every gained gene that a two-gene R_J representative overlaps still leaves
    R_J ahead on every contig.
  - The 3-14 losses per contig have two causes (see H3): the representative moves to a neighbour gene in a multi-gene locus, or it stays on
    the gene but carries fewer well-supported junctions.
- **Held-out families are mixed, and chr6 loses one.**
  - chr2 improves: Compara F 0.588 → 0.609, Soto 0.600 → 0.634, and the LIMS pair joins one family.
  - chr8 and chr10 do not move.
  - chr6 loses one 2-gene Compara family, BTN3A (BTN3A2 + BTN3A3). R_J represents the BTN3A2 locus by `DN_chr6_26200364_15`, a 2-read,
    15-exon chain from the lncRNA LOC124901288 through BTN3A2 into BTN2A2 (26,200,365-26,251,178). R_M uses BTN3A2's own 5-read, 11-exon
    body.

  At chr6, the pair's one PAF record is the same in both arms (`loci.paf` is byte-identical). What changes is the shared-exon fraction the
  edge rule computes on that record: the record's exon-to-exon bases divided by the smaller locus's exonic length, with threshold 0.60.
  - R_M: 1,524 / 1,707 = 0.893.
  - R_J: 1,300 / 2,174 = 0.598, which is 0.002 under the bar.

  This is a re-derivation outside the binary, from `loci.gff3` + `loci.paf`, and it matches the reviewer's independent figure. Two things
  move the fraction:
  - The chain adds exons outside BTN3A2 (the lncRNA's and BTN2A2's), so the exonic length grows 1,707 → 2,174.
  - Its BTN3A2 exon at 26,244,169 ends at 26,244,363 instead of 26,244,587, 224 bp shorter, so the exon-to-exon bases fall 1,524 → 1,300.

  Contig-wide, `rejected_low_shared_exon` goes 9,585 → 9,824 and the edges 2,082 → 2,046. **This is the fused-locus caveat** (see the ⚠ at
  the top).
- **On the development contig, R_J's full-length representatives split NPIP by subfamily.** This is a post hoc reading of a registered dev
  readout, not a pre-registered test. Under R_M, 18 of the 25 annotated chr16 NPIP copies map to one 27-copy family whose representatives
  are mostly 2-exon fragments. Under R_J, 19 annotated copies map to four families: the same 18 minus NPIPA7, plus NPIPA6 and NPIPB9.
  - 11 annotated copies (NPIPB7-B11, B13, B14P, B15, LOC124907834/807/808) in a family of 18 de novo copies;
  - 4 (NPIPA1/A6/A8/A9) in a family of 5;
  - 2 (NPIPA2/A5) in a family of 2;
  - 2 (NPIPB4/B5, which share one de novo locus) in a family of 3.

  The family-level references score that split as a loss: Compara 0.615 → 0.541, the NPIP set 0.576 → 0.481, Liftoff pairs 15/36 → 8/36.
  Development results are reported only and decide nothing.

## Setup (as registered, with one forced input setting)

- Inputs: the Figure 7 per-contig de novo GTFs (`/mnt/linuxdisk/tmp/rustle_figures/fig7/current/<species>_<contig>.denovo.gtf`, read-only,
  copied). Each arm ran the driver's `families` stage once per contig, with one binary: `mcl_families` sha1
  `96728e963222db9d90ad49249828fc8a8decb265` and `family_score` `7723029bb3af134da0d8bc5606afde743d7474b0`, both identical before and after
  every call. R_M had `RUSTLE_REPRESENTATIVE` unset (`env -u`). R_J had `RUSTLE_REPRESENTATIVE=most-junctions` on its one command, and its
  `params.tsv` ends with `representative most-junctions` (that row is absent in R_M's, as checked). Shipped stage defaults:
  `--min-exonic-bp 1 --min-shared-exon-frac 0.60 --emit-units`, and `--min-cov-shorter` is the binary's 0.70.
- ⚠ **`RUSTLE_BRIDGE_REGROUP=off` on both arms (forced, not chosen).** The prereg lists `--bridge-regroup f1v2` among the shipped defaults.
  But the Figure 7 GTFs were assembled on 2026-09-25, before f1v2 became the default (2026-09-29). They hold 0 `fusion_of` relations and no
  `PREFIX.families.gtf`, and the driver's guard refuses f1v2 on them. Bridge transcripts are therefore still loci members here. Under the
  current default some of them would become relations; whether that includes the chr6 chain is not measured. f1v2 can regroup these base
  GTFs from the BAM without a re-assembly (the frozen `f1_bridge.py` + `f1v2.py`; `docs/PENDING_2026-10-04.md` item 5), but that is a new
  test with its own pre-registration, not a re-scoring of this one.
- The R_M numbers are not Figure 7's printed dev numbers. The family binary's defaults changed after 2026-09-25 (`--min-cov-shorter 0.70`
  since 09-29, the copy table), so chr8 Compara is 0.388 here vs 0.364 in Fig. 7. The comparison is between the two arms of this run.
- References and scorers are Figure 7's own, imported: `_o1_recovery.family_score` (the binary), `score_arm` and
  `check_against_family_score` (every row agrees). The GFF slice is Figure 7's. References:
  - Compara at Primates, restricted to the contig. It was built by `_o1.filter_rows` from the cached `compara.Primates.families.tsv`, because
    `compara_contig_truth` refuses today: the raw Compara download it checks first, `human_e116.tsv`, is no longer on disk. All five
    restrictions are byte-identical to Figure 7's cached `human_<c>_compara.tsv`.
  - Soto 2025.
  - The NPIP reference set (chr16).
  - The contig's protein-homology families (Figure 7 cache).
  - Liftoff self-lift pairs via `_liftoff.copy_pairs(loci, 0.95, support, both loci on the contig)` + `pair_families`, on each arm's copy
    table. Read support comes from the fig8 products of human_A119b / gorilla_OR6737.
- H2 counts junctions as gaps >= 50 bp between a copy's exons (the rule's own junction). `n_exon - 1` (the prereg's wording) is shown in
  brackets. The two disagree on 0-5 copies per contig and never change a median or a fraction.
- H3: `bench/copy_support.py`, unchanged. Copies are the contig's annotated protein-coding genes: human from the CAT/Liftoff v2.0 slim GFF
  (`gene_id` = the CAT ID), gorilla from the RefSeq `GGO_genomic.gff`. `score.py h3-inputs` writes the copies table and a truth GTF of every
  transcript of those genes. The strict FOUND rule:
  - Spliced-expressed: >= 2 same-strand primary reads, each with >= k read-supported junctions, where k = min(2, annotated introns) and a
    supported junction is >= 3 reads. For k = 0, a read counts when its blocks cover >= 50% of the exon union.
  - Found: a same-strand locus whose representative carries >= k of those junctions. For k = 0, the representative's exons must cover
    >= 50% of the exon union.
- Check of the binary on real data: a Python re-derivation of both rules from each GTF reproduces every locus's `loci.gff3` exons under both
  arms (0 mismatches on 7 contigs, 25,027 loci). `loci.fa` and `loci.paf` are byte-identical between the arms on every contig, as the prereg
  said.

## H1 — families (bipartite F (sensitivity / precision); Liftoff = pairs in one family / pairs with both loci read-supported)

| contig | status | arm | Compara F (s / p) | Soto F (s / p) | NPIP set F | protein homology F (s / p) | Liftoff pairs |
|---|---|---|---|---|---|---|---|
| chr16 | development | R_M | 0.615 (0.444 / 1.000) | 0.617 (0.465 / 0.917) | 0.576 | 0.225 (0.127 / 0.975) | 15/36 |
| chr16 | development | R_J | 0.541 (0.370 / 1.000) | 0.615 (0.451 / 0.970) | 0.481 | 0.199 (0.111 / 0.971) | 8/36 |
| NC_073244.2 | development | R_M | — | — | — | 0.086 (0.045 / 1.000) | 0/1 |
| NC_073244.2 | development | R_J | — | — | — | 0.089 (0.046 / 1.000) | 0/1 |
| chr2 | held out (reused verdict set) | R_M | 0.588 (0.417 / 1.000) | 0.600 (0.436 / 0.960) | — | 0.157 (0.089 / 0.661) | 0/11 |
| chr2 | held out (reused verdict set) | R_J | 0.609 (0.438 / 1.000) | 0.634 (0.473 / 0.963) | — | 0.161 (0.091 / 0.678) | 0/11 |
| chr6 | held out (untouched) | R_M | 0.444 (0.286 / 1.000) | 0.800 (0.667 / 1.000) | — | 0.119 (0.064 / 0.812) | 1/2 |
| chr6 | held out (untouched) | R_J | **0.372** (0.229 / 1.000) | 0.800 (0.667 / 1.000) | — | 0.116 (0.062 / 0.893) | 1/2 |
| chr8 | held out (reused verdict set) | R_M | 0.388 (0.245 / 0.929) | 0.280 (0.163 / 1.000) | — | 0.179 (0.102 / 0.762) | 0/7 |
| chr8 | held out (reused verdict set) | R_J | 0.388 (0.245 / 0.929) | 0.280 (0.163 / 1.000) | — | 0.184 (0.105 / 0.767) | 0/7 |
| chr10 | held out (reused verdict set) | R_M | 0.542 (0.394 / 0.867) | 0.557 (0.415 / 0.844) | — | 0.157 (0.087 / 0.759) | 0/6 |
| chr10 | held out (reused verdict set) | R_J | 0.542 (0.394 / 0.867) | 0.547 (0.400 / 0.867) | — | 0.150 (0.083 / 0.750) | 0/6 |
| NC_073234.2 | held out (untouched) | R_M | — | — | — | 0.020 (0.010 / 1.000) | 0/0 |
| NC_073234.2 | held out (untouched) | R_J | — | — | — | 0.020 (0.010 / 1.000) | 0/0 |

Reference families whose match changes:

| contig | reference family | R_M | R_J |
|---|---|---|---|
| chr2 | Compara CF250 LIMS (2) | 1 member | 2 members |
| chr6 | Compara CF314 BTN3A (2) | 2 members | 0 members |
| chr16 | Compara CF153 NPIP (19) | 10 members | 7 members |
| chr16 | Compara CF154 SLX (2) | 1 member | 0 members |
| chr16 | NPIP set ID_154 (19) | 12 hits / 17 predicted | 8 / 11 |

The 2 Liftoff pairs on chr6 and the 11 / 7 / 6 on chr2 / chr8 / chr10 are too few to read. Gorilla NC_073234.2 has none with both loci on
the contig.

## H2 — copy table structure (junction = gap >= 50 bp; `n_exon - 1` in brackets)

| contig | arm | copies | families | median junctions / copy | mean junctions / copy (gap >= 50 bp) | copies with >= 2 junctions | total exon bp | loci whose rep changed |
|---|---|---|---|---|---|---|---|---|
| chr16 | R_M | 366 | 110 | 1 [1] | 2.08 | 0.249 [0.249] | 1,525,952 | |
| chr16 | R_J | 366 | 105 | 1 [1] | 3.11 | 0.287 [0.287] | 1,468,532 | 547 / 2,802 |
| NC_073244.2 | R_M | 94 | 29 | 3 [3] | 3.11 | 0.968 [0.968] | 270,596 | |
| NC_073244.2 | R_J | 93 | 29 | 4 [4] | 3.71 | 0.946 [0.946] | 258,350 | 323 / 1,142 |
| chr2 | R_M | 456 | 131 | 1 [1] | 1.71 | 0.250 [0.250] | 2,841,771 | |
| chr2 | R_J | 452 | 131 | 1 [1] | 2.62 | 0.265 [0.265] | 2,771,563 | 973 / 7,163 |
| chr6 | R_M | 250 | 58 | 0 [0] | 1.56 | 0.240 [0.240] | 1,874,626 | |
| chr6 | R_J | 241 | 55 | 0 [0] | 2.01 | 0.253 [0.253] | 1,807,827 | 705 / 5,012 |
| chr8 | R_M | 275 | 53 | 0 [0] | 0.93 | 0.164 [0.164] | 1,660,964 | |
| chr8 | R_J | 272 | 52 | 0 [0] | 1.38 | 0.180 [0.180] | 1,633,380 | 523 / 4,076 |
| chr10 | R_M | 215 | 56 | 1 [1] | 1.63 | 0.251 [0.251] | 1,365,969 | |
| chr10 | R_J | 219 | 58 | 1 [1] | 2.84 | 0.292 [0.292] | 1,349,825 | 601 / 3,876 |
| NC_073234.2 | R_M | 15 | 7 | 2 [2] | 2.60 | 0.600 [0.600] | 56,928 | |
| NC_073234.2 | R_J | 15 | 7 | 2 [2] | 3.27 | 0.600 [0.600] | 61,059 | 353 / 956 |

Clause (c) ties on all five held-out contigs: the median does not move. Most human family members have a 0-1 junction representative under
both rules (single-transcript loci have no alternative). The mean, though, rises on every contig (chr2 1.71 → 2.62, chr10 1.63 → 2.84), and
so does the >= 2 fraction everywhere except NC_073244.2 (0.968 → 0.946) and NC_073234.2 (0.600 = 0.600). Total exon bp falls on 6 of 7
contigs: the junction-richest transcript is often shorter in exon bp than the most-read one. 13-37% of loci change representative.

## H3 — annotated protein-coding genes FOUND (strict rule)

| contig | genes | spliced-expressed | found R_M | found R_J | R_J gains (a two-gene R_J rep overlaps) | R_J losses (a two-gene R_M rep overlaps) | k = 2 genes: R_M / R_J | locus level R_M / R_J | two-gene reps R_M / R_J, all loci (in the copy table) |
|---|---|---|---|---|---|---|---|---|---|
| chr16 | 857 | 761 | 636 | 673 | 51 (11) | 14 (0) | 616 / 652 | 675 / 685 | 12 / 27 (3 / 5) |
| NC_073244.2 | 1,520 | 1,005 | 884 | 909 | 29 (14) | 4 (0) | 843 / 866 | 898 / 913 | 12 / 30 (1 / 1) |
| chr2 | 1,243 | 1,111 | 960 | **1,030** | 78 (14) | 8 (0) | 926 / 996 | 1,018 / 1,037 | 6 / 25 (0 / 0) |
| chr6 | 1,047 | 900 | 757 | **815** | 62 (14) | 4 (0) | 707 / 763 | 805 / 818 | 14 / 28 (0 / 0) |
| chr8 | 698 | 574 | 500 | **531** | 39 (3) | 8 (0) | 481 / 511 | 529 / 534 | 4 / 9 (0 / 0) |
| chr10 | 729 | 650 | 562 | **598** | 41 (13) | 5 (1) | 547 / 583 | 590 / 605 | 4 / 17 (0 / 3) |
| NC_073234.2 | 1,119 | 818 | 734 | **750** | 19 (7) | 3 (0) | 704 / 719 | 744 / 752 | 1 / 7 (0 / 0) |

- A **two-gene representative** has exons overlapping the exon unions of >= 2 same-strand protein-coding genes of the GFF slice that are
  exon-disjoint from each other (RefSeq readthrough genes excluded). R_J makes 2-7 times more of them at the locus level (chr2 6 → 25,
  chr10 4 → 17).
- In the copy table they stay rare (chr10 0 → 3, chr16 3 → 5, otherwise unchanged), but that count is survivorship-biased. A fused
  representative can drop its locus out of its family, and so out of the copy table: the chr6 BTN3A2 locus (BTN2A2 + BTN3A2 under R_J) is
  one.
- The bracketed gain counts are an upper bound on what fused representatives explain: they count overlap, not the junction-carrying locus.
  R_J stays ahead without them: chr2 1,016, chr6 801, chr8 528, chr10 585, NC_073234.2 743, each vs R_M's count above.
- **The 46 losses (3-14 per contig; 28 held-out) have two causes.** The table below comes from a recount over the lost genes with
  `bench/copy_support.py`'s own functions (read filter, junctions, support). "Moves off" means R_J's representative of the locus that found
  the gene under R_M no longer overlaps the gene's exons; the split is the same with territory overlap.

  | cause | all | held-out | of which |
  |---|---|---|---|
  | moves off the gene | 30 | 17 | 28 (held-out 17) in loci whose span holds another same-strand RefSeq protein-coding gene's exons |
  | stays on the gene, fewer supported junctions | 16 | 11 | 14 (held-out 9) in single-gene loci |

  The reviewer's independent recount, with its own overlap definition, gave 28 moved / 18 stayed (held-out 16 / 12). Both counts find the
  same two causes and differ by 2 genes (1 held-out).
  - **Moves off the gene:** R_J picks a transcript with more junctions that does not overlap the gene, in 28 of 30 cases a neighbour's in
    a multi-gene locus.
    - HSP90AB1 (chr6, 9,966 support reads): its locus `DN_chr6_44080652_12` (44,053,598-44,088,620) also holds an upstream gene. R_M picks
      HSP90AB1's 606-read, 11-junction transcript; R_J picks the neighbour's 61-read, 12-junction transcript (44,053,598-44,070,604).
    - ATF2 (chr2): locus 175,288,066-175,657,658. R_J takes a 16-junction transcript at 175,288,135-175,494,189.
    - PPM1B (chr2): locus 44,172,945-44,777,912. R_J takes a 10-junction transcript at 44,255,365-44,327,759.
  - **Stays on the gene but carries fewer of its >= 3-read-supported junctions:**
    - ZNF174 (chr16): R_M's representative, 3 exons / 2 junctions, carries 2; R_J's, 5 exons / 4 junctions, carries 1.
    - CUTC (chr10): 8 → 1.

    The junction count is blind to per-junction support: the assembler admits a junction at 2 reads, H3 calls it supported at >= 3.
- **The 25 NPIP copies (chr16; the page's question).** R_M reproduces the chr16-wide reading of the GOOD arm exactly: strict 10, locus
  level 20, overlap 22. Those figures are in `/mnt/linuxdisk/tmp/readpool_npip/support_hsa.json`.
  `docs/archive/2026-10/SPLICED_COPY_SUPPORT_2026-10-04.md` quotes 9 / 19 instead, because it counts only the page's NPIP-cluster nodes.

  | arm | strict found | locus-level found | any same-strand overlapping locus |
  |---|---|---|---|
  | R_M | 10 | 20 | 22 |
  | R_J | **20** | 20 | 22 |

  - Gained (11): NPIPA6, NPIPA8, NPIPA9, NPIPB5, NPIPB7, NPIPB9, NPIPB12, NPIPB14P, NPIPB15, LOC124907808, LOC124907807. At NPIPA9 the
    representative goes from 1 supported junction to 22.
  - Lost (1): NPIPB4. Its locus (22,350,725-22,422,849) also holds a 13-junction transcript downstream of the copy
    (22,396,489-22,419,976), and R_J picks that one.
  - Not found under either rule:
    - NPIPB2 and NPIPB6: no same-strand locus representative overlaps them.
    - NPIPB13: not spliced-expressed.
    - NPIPA7: R_M's representative carries 1 junction. Its locus (16,329,977-16,406,194) also spans NPIPA6, and R_J's representative is
      NPIPA6's 22-junction transcript (16,329,977-16,359,006).

## Run times (`/usr/bin/time -v`)

| contig | families R_M, s (max RSS GB) | families R_J, s (max RSS GB) | H3 copy_support, s (max RSS GB) |
|---|---|---|---|
| chr16 | 43 (2.36) | 43 (2.37) | 67 (0.24) |
| NC_073244.2 | 40 (2.14) | 40 (2.18) | 30 (0.07) |
| chr2 | 390 (3.18) | 392 (3.06) | 175 (0.17) |
| chr6 | 286 (3.04) | 282 (2.85) | 87 (0.16) |
| chr8 | 163 (2.70) | 160 (2.70) | 63 (0.23) |
| chr10 | 123 (2.39) | 121 (2.39) | 65 (0.14) |
| NC_073234.2 | 55 (2.48) | 43 (2.62) | 18 (0.08) |

The representative rule costs nothing measurable: the all-vs-all dominates and is the same alignment under both rules.

## Decision (held-out contigs only, as registered)

| contig | (a) primary reference | (b) H3 found | (c) H2 median junctions / copy |
|---|---|---|---|
| chr2 | PASS: Compara F 0.609 vs 0.588 (bar 0.583), prec 1.000 vs 1.000 (bar 0.990) | PASS: 1,030 vs 960 | PASS: 1 vs 1 [n_exon-1: 1 vs 1] |
| chr6 | **FAIL: Compara F 0.372 vs 0.444 (bar 0.439)**, prec 1.000 vs 1.000 (bar 0.990) | PASS: 815 vs 757 | PASS: 0 vs 0 [0 vs 0] |
| chr8 | PASS: Compara F 0.388 vs 0.388 (bar 0.383), prec 0.929 vs 0.929 (bar 0.919) | PASS: 531 vs 500 | PASS: 0 vs 0 [0 vs 0] |
| chr10 | PASS: Compara F 0.542 vs 0.542 (bar 0.537), prec 0.867 vs 0.867 (bar 0.857) | PASS: 598 vs 562 | PASS: 1 vs 1 [1 vs 1] |
| NC_073234.2 | not testable: 0 Liftoff pairs with both loci on the contig | PASS: 750 vs 734 | PASS: 2 vs 2 [2 vs 2] |

**R_M stays the default; R_J stays opt-in** (`RUSTLE_REPRESENTATIVE=most-junctions` / `--representative most-junctions`, unchanged since
559649a9).
- ⚠ This is the verdict on the configuration as run. Both arms ran on the same pre-f1v2 GTFs with `RUSTLE_BRIDGE_REGROUP=off` (forced).
  The chain that fails chr6 is a fused chain of the kind f1v2 may turn into a relation; whether it would is unmeasured.
- The prereg is silent on an empty primary reference. NC_073234.2's clause (a) is reported as not testable, and the verdict does not
  depend on it: chr6's violation alone keeps R_M, so counting NC_073234.2 as a pass or as a fail changes nothing.
- The development results (chr16 families down, NPIP strict found 10 → 20) do not enter the decision.

## What follows (not run here)

- **The lead for any follow-up rule is per-junction support.** R_J counts junctions, not how well each is supported. 16 of its 46 H3
  losses stay on the gene but carry fewer >= 3-read-supported junctions: the assembler admits a junction at 2 reads, H3's "supported" is
  >= 3 reads. A support-aware representative rule is the next arm. It must be pre-registered, held-out first, and never tuned on chr16.
- The other 30 losses move off the gene, 28 of them to a transcript in a locus that also holds another protein-coding gene. In such a
  locus, any one-transcript representative stands for one gene only.
- **Two consumers.** The representative serves two consumers:
  - The copy table (O2, per-copy structure) is better under R_J on every contig (H3).
  - For the family edge rule, the held-out result is mixed (chr2 up, chr6 down, chr8 and chr10 equal) and the development contig is worse.

  Hypothesis, untested: the alignments are the same under both rules (`loci.paf` is identical), and what changes is which exon bases the
  shared-exon rule scores over them. 2-exon fragments would share most of their exon bases with their paralogs' fragments. Full-length
  representatives add subfamily-specific or fused exons to the denominator; BTN3A2 is one measured case (0.893 → 0.598). Separating the two
  consumers would be a new rule, to be pre-registered: the current representative for the family graph, R_J's (or a support-aware) transcript
  for the copy table O2 reads.
- **The chr6 failure is a fused chain; how f1v2 would classify it is unmeasured.** A re-test on f1v2 GTFs is a new pre-registration, with
  the same clauses. f1v2 can regroup the existing base GTFs from the BAM without a re-assembly, using the frozen
  `/mnt/linuxdisk/tmp/rustle_figures/f1_frozen/f1_bridge.py` and `/mnt/linuxdisk/tmp/rustle_figures/f1v2_frozen/f1v2.py`.
- O2 was not re-run, as registered. R_J changes the copy table O2 reads; its effect is a separate measurement.

Register rows 1234-1236.

## Amendment A of the spliced-support prereg (2026-10-04 11:45): H3 re-scored with support anchored on the ANNOTATED introns

The user's correction to `docs/archive/2026-10/PREREG_spliced_copy_support_2026-10-04.md` (Amendment A, 8735ceeb, before this re-score): a read or a
representative supports a gene only through the gene's own annotated introns (exact splice sites), not through any read-supported junction.
H3 re-scored from the saved arms (`h3A.support.*` in each contig's work dir; no families re-run; `bench/copy_support.py` `ann_*` columns):

| contig | genes | spliced-expressed (annotated introns) | (read-defined) | FOUND R_M | FOUND R_J | read-defined found R_M | R_J | locus level R_M | R_J |
|---|---|---|---|---|---|---|---|---|---|
| human_chr2 | 1,243 | 1,127 | 1,111 | **943** | **1,026** | 960 | 1,030 | 1,017 | 1,035 |
| human_chr6 | 1,047 | 909 | 900 | **741** | **800** | 757 | 815 | 798 | 810 |
| human_chr8 | 698 | 586 | 574 | **498** | **524** | 500 | 531 | 529 | 534 |
| human_chr10 | 729 | 660 | 650 | **555** | **589** | 562 | 598 | 589 | 604 |
| gorilla_NC_073234.2 | 1,119 | 855 | 818 | **738** | **750** | 734 | 750 | 750 | 759 |
| human_chr16 (dev) | 857 | 766 | 761 | **619** | **653** | 636 | 673 | 662 | 673 |
| gorilla_NC_073244.2 (dev) | 1,520 | 1,075 | 1,005 | **893** | **918** | 884 | 909 | 911 | 927 |

- **Clause (b) holds under Amendment A as well: found genes rise under R_J on all five held-out contigs** (chr2 943 -> 1,026, chr6 741 ->
  800, chr8 498 -> 524, chr10 555 -> 589, NC_073234.2 738 -> 750). The verdict is unchanged: R_M stays the default (chr6 clause (a)).
- Annotation anchoring lowers the found counts by 0-3% on the human contigs for both arms (the read-defined rule credited reads spliced at
  junctions the annotation lacks) and raises the gorilla "spliced-expressed" denominators (RefSeq gorilla models carry introns the
  >= 3-read rule had not yet confirmed at low depth).

## Amendment B of the spliced-support prereg (2026-10-04 12:00): H3 under the CHAIN rule — the sign of clause (b) flips

Amendment B (prereg 1e2801b0, before this re-score): support and FOUND require the read's / representative's in-span junction chain to be a
contiguous sub-chain of an annotated transcript's intron chain (>= 2 junctions; single annotation on these contigs: CAT human, RefSeq
gorilla). Re-scored from the saved arms (`h3B.support.*`), no families re-run.

| contig | genes | chain-expressed (>= 2 ISM/FSM reads) | CHAIN FOUND R_M | R_J | locus level R_M | R_J | (Amendment A coordinate rule R_M / R_J) |
|---|---|---|---|---|---|---|---|
| human_chr2 | 1,243 | 1,091 | **827** | **564** | 971 | 987 | 943 / 1,026 |
| human_chr6 | 1,047 | 876 | **638** | **469** | 754 | 766 | 741 / 800 |
| human_chr8 | 698 | 561 | **436** | **295** | 502 | 506 | 498 / 524 |
| human_chr10 | 729 | 637 | **473** | **337** | 561 | 575 | 555 / 589 |
| gorilla_NC_073234.2 | 1,119 | 808 | **680** | **567** | 724 | 733 | 738 / 750 |
| human_chr16 (dev) | 857 | 735 | **538** | **363** | 622 | 633 | 619 / 653 |
| gorilla_NC_073244.2 (dev) | 1,520 | 1,004 | **823** | **738** | 880 | 896 | 893 / 918 |

- **Under the chain rule R_J finds FEWER genes than R_M on all seven contigs** (held-out chr2 827 vs 564, chr6 638 vs 469, chr8 436 vs 295,
  chr10 473 vs 337, NC_073234.2 680 vs 567): the junction-maximal transcript carries, more often than the most-read one, a junction the
  annotation lacks (an alternative donor/acceptor, an extra exon, a readthrough junction), and one such junction fails the whole chain. The
  locus level barely moves (R_J still slightly ahead: the loci hold an annotated sub-chain either way); it is the single representative that
  the chain test judges.
- **Decision unchanged, now on two clauses:** R_M stays the default; under Amendment B clause (b) is violated on every held-out contig as well.
- What this measures: agreement of ONE representative with the annotated isoform chains (ISM/FSM). A representative that is a real, novel
  isoform fails it; the metric is the user's standard for "does the node support the annotated gene", not a count of real transcripts. Any
  future representative rule is judged on this chain rule and the family metrics together.

## Amendment C (2026-10-04 12:25): H3 against the reads' own EXPRESSED chains (annotation-free)

Expressed chain = an identical >= 2-junction chain of >= 3 uniquely placed reads; a representative is FOUND when its in-span chain equals or is
a contiguous sub-chain of one (re-scored from the saved arms, `h3C.support.*`).

| contig | genes | expressed (>= 1 expressed chain) | FOUND C R_M | R_J | locus level R_M | R_J | (Amendment B, annotated chains, R_M / R_J) |
|---|---|---|---|---|---|---|---|
| human_chr2 | 1,243 | 1,047 | **920** | **891** | 980 | 998 | 827 / 564 |
| human_chr6 | 1,047 | 843 | **718** | **699** | 771 | 783 | 638 / 469 |
| human_chr8 | 698 | 539 | **475** | **439** | 506 | 511 | 436 / 295 |
| human_chr10 | 729 | 609 | **536** | **513** | 564 | 580 | 473 / 337 |
| gorilla_NC_073234.2 | 1,119 | 717 | **674** | **613** | 684 | 692 | 680 / 567 |
| human_chr16 (dev) | 857 | 684 | **587** | **538** | 627 | 637 | 538 / 363 |
| gorilla_NC_073244.2 (dev) | 1,520 | 866 | **807** | **753** | 816 | 830 | 823 / 738 |

- **R_M finds slightly more than R_J on every contig under C as well** (held-out chr2 920 vs 891, chr6 718 vs 699, chr8 475 vs 439, chr10 536 vs
  513, NC_073234.2 674 vs 613; the gap is a third of Amendment B's): the junction-maximal transcript often carries a junction that fewer than 3
  unique reads share, so it is not a sub-chain of any expressed chain; the locus level stays slightly in R_J's favour (chr2 980 vs 998).
- Verdict unchanged (R_M default); clause (b) is against R_J under B and under C, for R_J under A and the first rule — the four readings
  disagree on what a representative should carry, which is itself the finding: a representative should be an expressed chain, neither the
  most-read fragment nor the junction-maximal transcript. That rule ("the most-read EXPRESSED CHAIN of the locus") is the one to pre-register
  next; its H3 under C is bounded above by the locus-level column.

## Amendment D (2026-10-04 12:40): H3 with TSS-anchored support — the sign of clause (b) flips back to R_J

D' = the read / representative starts within ±150 bp of an expressed chain's TSS (modal 5' end of its unique carriers) and carries its first
min(3, n) introns; D = the same against the annotated models' TSS (re-scored from the saved arms, `h3D.support.*`).

| contig | genes | expressed (D': >= 2 reads from the expressed TSS with the first 3 introns) | FOUND D' R_M | R_J | locus level R_M | R_J | D (annotated TSS): expressed, found R_M / R_J |
|---|---|---|---|---|---|---|---|
| human_chr2 | 1,243 | 1,019 | **541** | **696** | 911 | 926 | 1,003, 442 / 537 |
| human_chr6 | 1,047 | 778 | **410** | **497** | 675 | 686 | 776, 347 / 372 |
| human_chr8 | 698 | 528 | **294** | **355** | 473 | 477 | 518, 239 / 283 |
| human_chr10 | 729 | 591 | **331** | **410** | 522 | 535 | 582, 282 / 301 |
| gorilla_NC_073234.2 | 1,119 | 690 | **446** | **465** | 604 | 609 | 598, 437 / 400 |
| human_chr16 (dev) | 857 | 672 | **350** | **412** | 566 | 575 | 669, 320 / 349 |
| gorilla_NC_073244.2 (dev) | 1,520 | 825 | **599** | **585** | 729 | 739 | 734, 546 / 513 |

- **Under a TSS-anchored test the junction-maximal representative wins on the human contigs** (D': chr2 541 vs 696, chr6 410 vs 497, chr8 294
  vs 355, chr10 331 vs 410, chr16 350 vs 412; gorilla NC_073234.2 446 vs 465, NC_073244.2 599 vs 585): reaching the 5' end is what the test
  rewards, and the most-read representative is the 3' fragment. Under D (annotated TSS) the same, with lower counts (the annotated starts).
- The four readings now split 2-2 (A, D for R_J; B, C for R_M) and the family clause (a) still fails on chr6. The decision stays as registered
  (R_M default, R_J opt-in); the next representative rule is pre-registered against the user's standard (D') together with the family
  metrics: the most-read expressed chain that starts at the locus's expressed TSS.

## Amendment E (2026-10-04 14:00): H3 anchored on CAPPED starts

| contig | genes | expressed (>= 2 reads from a capped start with the first 3 introns) | FOUND E R_M | R_J | locus level R_M | R_J |
|---|---|---|---|---|---|---|
| human_chr2 | 1,243 | 963 | **572** | **693** | 867 | 882 |
| human_chr6 | 1,047 | 713 | **431** | **490** | 628 | 636 |
| human_chr8 | 698 | 495 | **316** | **365** | 451 | 452 |
| human_chr10 | 729 | 551 | **354** | **404** | 494 | 506 |
| human_chr16 (dev) | 857 | 617 | **365** | **405** | 539 | 547 |
| gorilla_NC_073234.2 — NO cap signal in this library: E undefined | 1,119 | 275 | **170** | **190** | 247 | 248 |
| gorilla_NC_073244.2 (dev) | 1,520 | 280 | **183** | **178** | 255 | 259 |

- Human: the junction-maximal representative is ahead on every contig (held-out chr2 572 vs 693, chr6 431 vs 490, chr8 316 vs 365, chr10 354 vs
  404), as under D/D' — a TSS-anchored test rewards reaching the 5' end. **Gorilla OR6737 carries no cap signal** (`project_read_proven_ends`,
  09-27), so E is undefined there (the 280 / 275 "expressed" genes are chance G clips); D' stands for gorilla.
- Decision unchanged; the next representative rule is pre-registered against E on human and D' on gorilla, with the family metrics.
