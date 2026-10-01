# Pre-registration: are NPIP fusion transcripts duplicon-boundary crossings inside co-duplicated blocks? (2026-09-30, KEY=npipfusion)

Written before any statistic below was computed. Human CHM13 only, NPIP only. Nothing in `src/` or `tools/` changes.
Script: `bench/soto_m2/npip_fusion_duplicons.py`, written after this file is committed; the commit hash of this file is recorded
in the result section. After the first run the script may change only to fix a crash, and every such change is listed as a
deviation.

## 1. Question

A descriptive look at the four annotated NPIP read-throughs (section 6) suggests that NPIP fusions are not two families joined by
chance: each fused transcript switches from a partner's duplicons (PKD1's, PDXDC1's) to the NPIP core at a duplicon boundary,
and the partner's duplicons travel with the NPIP core from one segmental duplication to the next. **Does this hold for the NPIP
fusion transcripts seen in long reads, annotated or not?** If yes, an NPIP fusion is a property of the duplicated block, and the
right container is a path across one duplicon boundary. If the partners are single-copy sequence, NPIP fusions are copy-specific
events instead.

## 2. Inputs (frozen)

- **Reads.** Two human Iso-Seq libraries aligned to CHM13 v2.0 with minimap2 `-ax splice:hq`:
  - development: `/mnt/linuxdisk/home/juanfraitu/winloci_data/A119b.t2t.bam`;
  - held-out: `/home/juanfra/human_val/human_testis.t2t.bam`. **The verdict is the held-out library's.**
- **NPIP genes.** Genes on chr16 and chr18 whose name starts with `NPIP` and contains no `-`, from RefSeq
  (`winloci_data/Reference/chm13v2.0_RefSeq_full.gff.gz`: NPIPA1, A2, A5-A9, NPIPB2-B15 including B10P and B14P, NPIPB1P) or from
  the CAT v4 gene list used by the meeting page (`families_cn.json` `genes`; it adds copies RefSeq lacks, such as NPIPP1).
  `T_N` = the union of their exons; `B_N` = the union of their gene spans.
- **Duplicons.** DupMasker colours BED, Vollger et al. 2022 (`winloci_data/duplicons/chm13.draft_v1.0_plus38Y_dupmasker_colors.bed`),
  every overlapping record counts. **SD98 regions:** `winloci_data/soto_replication/sd98_v1.bed`.
- **Coordinates.** The duplicon and SD98 files are CHM13 v1.0; reads and RefSeq are v2.0. The script asserts that every SD98 region
  on chr16 and chr18 has the same sequence in `sd98_regions.fa` (cut from v1.0) and in `Reference/chm13v2.0.fa`, and stops otherwise.

## 3. Units and measures (fixed now)

- **Reads:** primary alignments (`-F 2308`) overlapping `B_N`, each read counted once. **Blocks:** the read's aligned reference
  segments split at `N` (`M`, `=`, `X`, `D` extend a block; `I`, `S`, `H` do not).
- **Block class:** `N` if it overlaps `T_N` by at least 1 bp; `O` if it overlaps no `B_N`; otherwise `I` (NPIP intron or a new NPIP
  exon). `I` blocks never form a switch.
- **Junction:** (chromosome, intron start, intron end) between two consecutive blocks. **Units are distinct junctions supported by at
  least 2 reads** (the project's floor). For each unit, each flank is its most frequent adjacent block across the supporting reads
  (ties: the longer block).
- **Switch unit:** one flank `N`, the other `O`. **Internal unit:** both flanks `N`.
- **Duplicons of a block:** the DupMasker IDs overlapping it; the **dominant** one covers the most bases (ties: ID order).
- `b(j)` **boundary:** the two flanks' duplicon sets share no ID (an empty set on either side counts as no shared ID).
- **Core duplicons `K`:** duplicons overlapping the exons of at least half of the NPIP gene records (RefSeq and CAT records counted
  separately).
- **Co-duplicated duplicon:** one that occurs in at least 2 SD98 regions each holding at least one segment of a core duplicon.
- `c(j)` **co-duplicated partner** (switch units only): the `O` flank's dominant duplicon exists and is co-duplicated.

## 4. Tests and decision rule (fixed now)

- **H1 boundary.** Fraction of units with `b = 1`, switch vs internal; one-sided Fisher exact test (switch greater). Passes if
  p < 0.01 and the switch fraction is the larger.
- **H2 co-duplicated partner.** Fraction of switch units with `c = 1`; one-sided exact binomial test against 1/2 (most partners
  co-duplicated). Passes if p < 0.01.
- **Verdict** (held-out decides; development reported beside it; if they disagree the verdict is SPLIT with both):
  - **EXPLAINED:** H1 and H2 pass.
  - **BOUNDARY ONLY:** H1 passes, H2 fails.
  - **NOT EXPLAINED:** H1 fails.
  - **UNDERPOWERED:** fewer than 10 switch units in the held-out library (reported, no verdict).

## 5. Secondary (descriptive, no verdict)

- MAPQ >= 1 reads only: H1 and H2 recomputed.
- **Location-matched null for H2:** for each switch unit with intron length g, a block of the `O` flank's length placed on the `O` side
  with its near edge at a distance drawn uniformly from [g/2, 2g] from the `N`-side splice site, redrawn (up to 100 times) until it
  overlaps no `B_N`; mean `c` over 10,000 replicates, seed 20260930. Says whether partners are more co-duplicated than sequence at a
  similar distance (enrichment beyond locality), which H2 does not ask.
- **Partner table:** per switch unit, the `O` flank's dominant duplicon, the RefSeq or CAT gene(s) it overlaps, the NPIP gene(s) of
  the `N` flank, read support.
- **Recurrence:** per partner duplicon, the number of distinct NPIP genes it is joined to (placement among identical copies is
  ambiguous, so this is not a test).
- Fraction of switch units whose `O` flank lies in no SD98 region (single-copy partners).
- Whether the junctions of the four RefSeq read-throughs are among the switch units.

## 6. Seen before writing this file (disclosed)

- RefSeq longest transcripts, 5' to 3' dominant duplicon per exon: PKD1P3-NPIPA1, PKD1P4-NPIPA8 and PKD1P6-NPIPP1 share one string,
  SD9474 (SD9605 or SD9609) SD9607 SD9611 SD9613 then SD9449 or SD9450, SD9443, SD9622; PDXDC2P-NPIPB14P is SD9585, SD9479 ... SD9456,
  SD9450, SD9449, SD9443, SD9445. PKD1 (chr16:2.1 Mb) exonic duplicons SD9613, SD9606, SD9605, SD9607, SD9474, SD9609; PDXDC1's main
  duplicon SD9479. The 4 SD98 regions other than PKD1's that carry PKD1 duplicons all carry NPIP duplicons. RefSeq NPIPA1 exonic
  duplicons SD9443, SD9622, SD9449; NPIPA8 SD9443, SD9613, SD9622, SD9450.
- A119b: the PKD1P6-NPIPP1 fusion junction chr16:15,120,015-15,126,650 has 110 MAPQ-60 reads (2026-09-18). A119b is therefore the
  development library.
- human_testis: no NPIP read has been examined for this question.
- Expectation written now: H1 will likely pass (gene ends tend to sit at duplicon ends); H2 is the informative test.

## 7. Result

(Filled in after the run, below this line, without editing anything above.)

**Amendment 1 (2026-09-30, before any statistic was computed).** The first run stopped at the section 2 coordinate check, as the
rule says it should: on chr16, CHM13 v1.0 and v2.0 differ by a 5 bp indel in the first telomeric repeat, so every v1.0 coordinate
on chr16 sits 5 bp to the right of v2.0 (47 of 64 chr16/chr18 SD98 regions differed unshifted). After shifting chr16 by -5, all
chr16 SD98 regions are identical except the two that touch the chromosome ends (chr16:0-12,460 and chr16:96,183,274-end, 539 and 151
mismatches); chr18 needs no shift (one mismatch, in chr18:0-217,988, which touches the chromosome start). NPIP lies at
chr16:11.9-30.7 Mb and 75-80 Mb and on chr18 at 11.8 Mb, away from all three. Changes: (1) every v1.0 coordinate on chr16 (CAT
exons and spans, DupMasker segments, SD98 regions) is shifted by -5 before it meets v2.0 data (reads, RefSeq); chr18 is not shifted;
the co-duplicated set is computed within v1.0, where duplicons and SD98 regions already agree; (2) the coordinate check becomes:
after the shift, every SD98 region on chr16 and chr18 that does not touch a chromosome end has identical sequence, else stop.
Nothing else changes. No read had been read and no unit counted when this was written.

---

Run on 2026-09-30 after this file was committed (`329ce4bc`, sha1 of the file at that commit
`de9224a059054f4d97b886a471f71d06f8a0be4f`) and after Amendment 1 (`f2b123fd`, with the script); 19 s, heavy lock. Partner table:
`docs/NPIP_FUSION_PARTNERS_2026-09-30.tsv`.

**VERDICT (section 4): none. The held-out library is UNDERPOWERED** (human_testis: 135 primary reads at NPIP genes, 3 switch units,
fewer than 10). The script printed "SPLIT"; by section 4 an underpowered held-out library gives no verdict, so SPLIT does not apply.
**The development library alone passes both tests (EXPLAINED).**

| library | reads at NPIP | switch units | internal units | boundary, switch vs internal | H1 p | co-duplicated partners | H2 p |
|---|---|---|---|---|---|---|---|
| development (A119b) | 17,948 | 103 | 725 | 0.689 vs 0.279 | 1.2e-15 | 89 / 103 = 0.864 | 8e-15 |
| development, MAPQ >= 1 | 13,567 | 94 | 686 | 0.691 vs 0.278 | 1.3e-14 | 80 / 94 = 0.851 | 1.1e-12 |
| held-out (human_testis) | 135 | 3 | 54 | 0.667 vs 0.259 | 0.19 | 1 / 3 | 0.88 |
| held-out, MAPQ >= 1 | 94 | 2 | 53 | 1.000 vs 0.245 | 0.071 | 0 / 2 | 1 |

Core duplicons found (section 3): SD5887, SD9443, SD9445, SD9449, SD9450, SD9456, SD9621, SD9622; 103 co-duplicated duplicons.

Secondary (section 5), development: location-matched null 0.785 vs observed 0.864, p = 0.0016 (partners are more co-duplicated than
sequence at a similar distance, by a modest margin: most sequence near NPIP is co-duplicated anyway). Single-copy partners (outside
every SD98 region): 18 of 103. Each of the four RefSeq read-throughs has exactly one intron among the switch units (PKD1P6-NPIPP1
chr16:15,120,014-15,126,650, the junction already seen on 09-18). Partners by reads: SD9526 (EIF3C) 2,744; SD9613 (PKD1P5, PKD1P6,
PKD1P1, PKD1P4 fusion parts) 884 over four loci; SD9527 (SMG1P1, SMG1P4) 350; SD9534 253; SD9559 (OTOAP1) 15.

Descriptive breakdown written after the run (not pre-registered): 35 of the 103 switch partners sit on core duplicons, i.e. NPIP
sequence outside every annotated NPIP gene (co-duplicated by definition); 54 on co-duplicated non-core duplicons; 6 on duplicons not
co-duplicated; 8 on no duplicon. Without the core-duplicon partners, 54 of 68 (0.79) are co-duplicated. Without the four RefSeq
read-through junctions seen before the run, 85 of 99 (0.86) are co-duplicated and 70 of 99 (0.71) cross a boundary.

**Reading.** In the one library with power, NPIP fusion transcripts, annotated or not, leave NPIP at a duplicon boundary and mostly
land on duplicons that travel with the NPIP core: EIF3C, the PKD1 pieces and SMG1P (LCR16u), the known neighbours in the 16p
mosaic. That supports the block reading (a fusion as a path across one duplicon boundary inside a co-duplicated block) for NPIP in
A119b. It is not held-out validated: the second human library has almost no NPIP expression, and no other human long-read library is
on disk. A third of the "partners" are NPIP-core sequence outside annotated NPIP genes, which is an unannotated-copy question (O3),
not a fusion with another family.
