# The locus representative, most-junctions (R_J) vs most-reads (R_M): R_J carries the copy's structure but costs families — chr6 fails clause (a), R_J stays opt-in (2026-10-04)

Prereg: `docs/PREREG_locus_representative_rule_2026-10-04.md` (3bf6aca1, written before any run). Flag: `mcl_families --representative
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
- **R_J does what it was built for.** On every contig the representative carries more of the gene's read-supported structure. Held-out found
  genes: chr2 960 → 1,030, chr6 757 → 815, chr8 500 → 531, chr10 562 → 598, NC_073234.2 734 → 750. That is close to the locus-level ceiling
  (chr2 1,037). On the 25 NPIP copies, strict FOUND goes from 10 to 20, equal to the locus-level reading. The gain does not depend on fused
  representatives: dropping every gained gene that a two-gene R_J representative overlaps still leaves R_J ahead on every contig.
- **It does not improve families, and on chr6 it costs one.** The whole chr6 loss is one 2-gene Compara family, BTN3A (BTN3A2 + BTN3A3). R_J
  represents the BTN3A2 locus by `DN_chr6_26200364_15`, a 2-read, 15-exon chain from the lncRNA LOC124901288 through BTN3A2 into BTN2A2
  (26,200,365-26,251,178). R_M uses BTN3A2's own 5-read, 11-exon body. The alignment does not change (`loci.paf` is byte-identical between the
  arms), so the pair is lost in the exon-dependent edge rule. **This is the fused-locus caveat.** These input GTFs predate the f1v2 bridge
  regroup (below).
- **On the development contig, R_J's full-length representatives split NPIP by subfamily.** Under R_M, 18 of the 25 annotated chr16 NPIP
  copies map to one 27-copy family whose representatives are mostly 2-exon fragments. Under R_J those copies map to four families:
  - the NPIPB copies (B7-B15, B14P, LOC124907834/807/808): one family of 18;
  - NPIPA1/A6/A8/A9: one of 5;
  - NPIPA2/A5: one of 2;
  - NPIPB4/B5: one of 3.

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
  current assembly default some of them would become relations; whether that includes the chr6 chain is not measured, because it needs a
  re-assembly. That would be a new test, not a re-scoring of this one.
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

| contig | arm | copies | families | median junctions / copy | mean | copies with >= 2 junctions | total exon bp | loci whose rep changed |
|---|---|---|---|---|---|---|---|---|
| chr16 | R_M | 366 | 110 | 1 [1] | 2.08 | 0.249 [0.249] | 1,525,952 | |
| chr16 | R_J | 366 | 105 | 1 [1] | 3.11 | 0.287 [0.287] | 1,468,532 | 547 / 2,802 |
| NC_073244.2 | R_M | 94 | 29 | 3 [3] | 3.11 | 0.968 [0.968] | 270,596 | |
| NC_073244.2 | R_J | 93 | 29 | 4 [4] | 3.71 | 0.946 [0.946] | 258,350 | 323 / 1,142 |
| chr2 | R_M | 456 | 131 | 1 [1] | 1.72 | 0.250 [0.250] | 2,841,771 | |
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
both rules (single-transcript loci have no alternative). The mean, though, rises on every contig (chr2 1.72 → 2.62, chr10 1.63 → 2.84), and
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
  chr10 4 → 17), but few reach the copy table (chr10 0 → 3, chr16 3 → 5, otherwise unchanged). The bracketed gain counts are an upper bound
  on what fused representatives explain: they count overlap, not the junction-carrying locus. R_J stays ahead without them (chr2 1,016,
  chr6 801, chr8 528, chr10 585, NC_073234.2 743, each vs R_M's count above).
- **The losses (3-14 per contig) come from over-merged multi-gene loci, where R_J moves the representative to the neighbour gene.** Example:
  HSP90AB1 on chr6, 9,966 support reads. Its locus `DN_chr6_44080652_12` (44,053,598-44,088,620) also holds an upstream gene. R_M picks
  HSP90AB1's 606-read, 11-junction transcript; R_J picks the neighbour's 61-read, 12-junction transcript. Two chr2 losses are the same
  case, in larger loci:
  - ATF2: locus 175,288,065-175,657,658. R_J takes a 16-junction transcript at 175,288,135-175,494,189.
  - PPM1B: locus 44,172,944-44,777,912. R_J takes a 10-junction transcript at 44,255,365-44,327,759.

  Under either rule a multi-gene locus is represented by one gene. The locus over-merge (`project_node_overmerge_reversal`) is the defect,
  not the representative rule.
- **The 25 NPIP copies (chr16; the page's question).** R_M reproduces the chr16-wide reading of the GOOD arm in
  `docs/SPLICED_COPY_SUPPORT_2026-10-04.md` exactly: strict 10, locus level 20, overlap 22. That doc's 9 / 19 count only the page's
  NPIP-cluster nodes.

  | arm | strict found | locus-level found | any same-strand overlapping locus |
  |---|---|---|---|
  | R_M | 10 | 20 | 22 |
  | R_J | **20** | 20 | 22 |

  - Gained (11): NPIPA6, NPIPA8, NPIPA9, NPIPB5, NPIPB7, NPIPB9, NPIPB12, NPIPB14P, NPIPB15, LOC124907808, LOC124907807. At NPIPA9 the
    representative goes from 1 supported junction to 22.
  - Lost (1): NPIPB4. Its locus (22,350,724-22,422,849) also holds a 13-junction transcript downstream of the copy
    (22,396,489-22,419,976), and R_J picks that one.
  - Not found under either rule:
    - NPIPB2 and NPIPB6: no same-strand locus representative overlaps them.
    - NPIPB13: not spliced-expressed.
    - NPIPA7: R_M's representative carries 1 junction. Its locus (16,329,976-16,406,194) also spans NPIPA6, and R_J's representative is
      NPIPA6's 22-junction transcript.

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
559649a9). The development results (chr16 families down, NPIP strict found 10 → 20) do not enter the decision.

## What follows (not run here)

- The representative serves two consumers, and they want different things:
  - The copy table (O2, per-copy structure) is better under R_J on every contig.
  - The family edge rule is better under R_M. Fragment representatives align fragment-to-fragment, so they keep coarse families together.
    Full-length representatives resolve subfamilies (NPIPA vs NPIPB) and are exposed to readthrough chains.
  
  One representative per locus cannot serve both. Separating them would be a new rule, to be pre-registered: R_M for the family graph,
  R_J's transcript for the copy table O2 reads, with the families unchanged.
- The chr6 failure is a bridge chain in a GTF assembled without f1v2. A re-test on f1v2 GTFs would be a new pre-registration: the contigs
  re-assembled with the current default, and the same clauses.
- O2 was not re-run, as registered. R_J changes the copy table O2 reads; its effect is a separate measurement.

Register rows 1234-1236.
