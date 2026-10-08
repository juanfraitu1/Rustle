# Hierarchy study: results (2026-10-07). Y block, concordance with the Yoo et al. 2025 gorilla expansions

**Status: DESCRIPTIVE; no verdict class.** Protocol: `docs/PREREG_hierarchy_duplicon_expansion_2026-10-07.md`, Amendment 3 (git 8fdff163, committed before the Y2 and Y3 numbers existed; Y1 was looked at before and is labelled spent). Results of T1 (Phases R0 to R2) and T2 (Phase E) will be appended to this file under their own headings. Truth files: `/mnt/linuxdisk/home/juanfraitu/winloci_data/yoo_2025/` (workbook sha256 865d4e868568b32d...; `table_38.tsv`, `yoo_GGO_table35_novel_paralogs.tsv`, `yoo_GGO_table40_lineage_specific.tsv`, `GGO_chr_map.tsv`, `yoo_GGO_unit_loci.tsv`, `y1_lrpap1_concordance.tsv`). Match rule: same chromosome, overlap of at least 50% of the shorter locus. Gorilla only; one animal (Jim, mGorGor1).

## T1-d, the gorilla SEDEF arm (Phase R2; descriptive, no verdict class)

**Protocol.** `docs/PREREG_hierarchy_duplicon_expansion_2026-10-07.md` section 5 and Amendment 1. Code `bench/hierarchy/t1_gorilla_pairs.py` with `bench/dna_sd_atoms.py` (restored unchanged, sha256 9d96d976...); 132 tests (`test_t1_gorilla_pairs.py`, `test_t1_gorilla_release.py`, `test_t1_gorilla_review_extras.py`), also under `PYTHONHASHSEED` 0, 1 and 2. Products in `/mnt/linuxdisk/tmp/hier_run/r2/` (`report.json`, `report_pairs_0.90|0.95|0.98.tsv`, `report.stdout`, `gate0_r2.tsv`, `gate0_r2_env.tsv`, `counts.json`, `counts.stdout`); the second hash-seed copy is in `hier_run/r2_report_hs1/`.

**Gates (all pass).** Gate 0: the five inputs and the atoms script match their registered sha256 prefixes (python 3.14.4, numpy 2.4.2, scipy 1.17.1, pysam 0.23.3). Gate 3: `report.json`, the three pair tables and the printed report are byte-identical under `PYTHONHASHSEED` 0 and 1. Gate 6: the verdict-set depth-matched pairs and families are 38/33, 15/15 and 5/5 as registered; under the registered rules the atoms, classes and edges differ from the recon values by -6, +22 and -186 (identical tuples collapsed, coverage >= 0.5). Positive control: the 8 LRPAP1 copies each touch at least one class at tau 0.90 and form one connected component (identity 0.957 to 0.990).

**Review and one disclosure.** Two independent reviewers of the first version found no defect that changes a registered number and asked for a validity coupling in `report` (INVALID and no rate if the control or a Gate 0 fingerprint fails), a fuller printout, hermetic guard tests and the interpreter in Gate 0; a delta reviewer then confirmed that no computed number changed (counts and report byte-identical on five made-up worlds, `counts` on the real inputs equal to the earlier record) and found no fail-open path (verdict: safe to release). After that review I made the reviewer's optional fail-closed hardening (all six inputs required in `validity`, a wrong self-lift table is INVALID with a reason, the registry is validated before anything is built, the registry source and the scorer's sha256 are recorded, stale pair tables are removed on INVALID, the B_viol caveat no longer calls the pool small); these were covered by new tests and were not reviewed again. **Disclosure:** during the first review a mutation run of one reviewer made the old guard tests run `report` on the real inputs once (18:12) and write `report.json` and three pair tables into the real output directory; the reviewer never opened them, deleted them by exact path at 18:16, and no number was seen by anyone; the report is deterministic, so nothing was lost. The leaking tests were replaced by hermetic ones.

**Result (released 2026-10-07T19:25:47 PDT, scorer sha256 e48c4c29...; one run, no tuning; tau 0.90 is the primary level).**

| tau | depth-matched pairs / families | S-rate, pair-weighted [Wilson] | family-weighted | A_viol | bootstrap, pair-weighted | path pairs | aligned < 0.5 | B_viol | UNDERPOWERED |
|---|---|---|---|---|---|---|---|---|---|
| 0.90 | 38 / 33 | 0.789 [0.637, 0.889] | 0.773 | 0.211 | 0.641 to 0.907 | 23 of 38 | 9 | 9 of 11,422 | no |
| 0.95 | 15 / 15 | 0.733 [0.480, 0.891] | 0.733 | 0.267 | 0.467 to 0.933 | 6 of 15 | 4 | 9 of 4,250 | yes |
| 0.98 | 5 / 5 | 0.800 [0.376, 0.964] | 0.800 | 0.200 | 0.400 to 1.000 | 2 of 5 | 2 | 4 of 1,419 | yes |

B_viol at tau 0.90 by the registered distance classes (the same count, split): cross-contig 0 of 10,506, 100 kb to 1 Mb 0 of 27, under 100 kb 5 of 11, over 1 Mb 4 of 878. B_viol is a ceiling-type quantity (the pool is not identity-matched and 92% cross-contig).

**Reading (descriptive; no verdict class is issued).** At tau 0.90, 30 of 38 depth-matched same-family gene pairs lie in a shared SEDEF atom class, and 8 do not. Different-family pairs almost never share a class (9 of 11,422). The 8 pairs without a shared class (post hoc listing from `report_pairs_0.90.tsv`, no statistic): five are cross-contig pairs of a multi-exon parent and a compact copy with one or two exon blocks (GLUD1 13 blocks and GLUD2 1; RBMX 12 and RBMXL1 2; UTP14A 15 and UTP14C 1; RPE 9 and RPEL1 2; FKBP1A 4 and FKBP1C 1), the signature of retrocopies, which are not segmental duplications; three are tandem ancient paralogs on one contig (GON4L and YY1AP1, APOBEC3D and APOBEC3F, FCGR2A and FCGR2B), two of them aligned over less than 0.56 of the shorter transcript. So in the gorilla the share of same-family pairs that sit in a common SD atom is about four fifths in this roster, and most of the rest are family relations that are not SD-derived, which fits the layering SD -> duplicon -> family with retrocopies and old tandem paralogs entering the family level by another route.

**Limits (registered, repeated).** The S-rate measures SEDEF and atom sensitivity at the genes as much as unit boundaries, because depth-matching on transcript identity >= tau selects the pairs for which a SEDEF alignment at fracMatch >= tau is expected. The layer is a human Compara family projection by exact symbol onto the gorilla RefSeq annotation (depleted of LOC-named recent copies, not gorilla truth). The within-contig control pool is tiny (11 and 27 pairs below 1 Mb). The gorilla roster is SD-rich by construction (the development slice and the most-junctions contig are excluded); nothing here is compared with the human T1 cells, and no species comparison is claimed. Path genes (a gene touching two or more classes) are 23 of the 38 pairs, so the pair-level S depends on how classes are cut.

## Y1, LRPAP1 (spent)

Our 11 loci (8 full-length copies and 3 5' fragments; `docs/LRPAP1_FAMILY_2026-10-04.md`) against Yoo's Table VIII.38 (the ancestral copy and 9 copies, 10 loci) and the 5 AMRP rows of Table VIII.35 (AMRP is the UniProt mnemonic of LRPAP1):

| our locus | class | chromosome | Yoo locus | table | overlap of ours / of Yoo's |
|---|---|---|---|---|---|
| LRPAP1 | full-length | chr3:12,090,719-12,110,719 | LRPAP1-anc | VIII.38 | 1.00 / 1.00 |
| LOC134756753 | full-length | chr12:22,553,726-22,573,741 | LRPAP1-1 (also AMRP 358 aa) | VIII.38, VIII.35 | 0.97 / 1.00 |
| copyB (unannotated) | full-length | chr12:24,785,460-24,805,469 | LRPAP1-3 | VIII.38 | 0.99 / 1.00 |
| copyC (unannotated) | full-length | chr12:30,205,371-30,225,550 | LRPAP1-4 | VIII.38 | 0.99 / 0.71 |
| LOC129526389 | full-length | chr14:23,392,677-23,412,705 | LRPAP1-5 (also AMRP 389 aa) | VIII.38, VIII.35 | 0.99 / 1.00 |
| LOC134757218 | full-length | chr16:15,822,038-15,842,048 | LRPAP1-7 | VIII.38 | 0.99 / 1.00 |
| LOC129523574 | full-length | chr22:11,292,022-11,312,055 | LRPAP1-8 | VIII.38 | 0.99 / 1.00 |
| LOC129530227 | full-length | chrY:45,277,394-45,297,429 | LRPAP1-9 | VIII.38 | 0.99 / 1.00 |
| LOC134756368 | 5' fragment | chr12:23,070,423-23,086,604 | LRPAP1-2 (solitary; AMRP 258 aa) | VIII.38, VIII.35 | 1.00 / 1.00 |
| LOC115932954 | 5' fragment | chr14:25,250,978-25,267,164 | LRPAP1-6 (solitary; AMRP 258 aa) | VIII.38, VIII.35 | 1.00 / 1.00 |
| LOC115932756 | 5' fragment | chr16:17,176,045-17,192,323 | AMRP 258 aa | VIII.35 only | 1.00 / 1.00 |

All 10 loci of Table VIII.38 match one of our loci, and all 11 of our loci are in Yoo's supplementary catalogue (10 in VIII.38; the chr16 fragment only in VIII.35, as a novel paralogous gene). Sensitivity 10/10 on Table VIII.38, precision 11/11 against the union of VIII.38 and VIII.35. The text count of 10 copies (1 ancestral and 9) and our 11 differ only by the chr16 fragment that Yoo lists in another table. Yoo's solitary copies 2 and 6 are our chr12 and chr14 fragments with identical coordinates. Both Yoo tables for every locus, recomputed by `python3 -B bench/hierarchy/yoo_concordance.py y1 ...`: `winloci_data/yoo_2025/y1_lrpap1_all_matches.tsv` (sha256 b5e729097c973abd...). Our O1 and O2 results for these loci (de novo family of the expressed copies, no contested read, no haplotype-only copy) are in `docs/LRPAP1_FAMILY_2026-10-04.md`.

## Y2, the chr1 MAPKBP1 / JMJD7-PLA2G4B / SPTBN5 unit

Truth: the 27 loci of Table VIII.38 (8 chr1 triplets, the ancestral chr16 triplet; 9 MAPKBP1, 9 JMJD7-PLA2G4B (the 8 chr1 loci plus the ancestral chr16 record, which Table VIII.38 labels PLA2G4B), 9 SPTBN5 loci). Ours: copy rows of the BASE family tables (genome-wide, pre-f1v2 snapshot, one library at a time). The KB3781 table is the registered one; the OR6737 table is an addition made after the first result and is labelled so.

| library | gene | Yoo loci | matched by a BASE copy | families touched | primary reads at the loci with at least 2 reads |
|---|---|---|---|---|---|
| KB3781 (fibroblast, the animal's own) | JMJD7-PLA2G4B | 9 | 8 | MCL11 | 9 of 9 (26 to 138 reads) |
| | MAPKBP1 | 9 | 1 | MCL11 | 2 of 9 |
| | SPTBN5 | 9 | 1 | MCL11 | 3 of 9 |
| OR6737 (testis) | JMJD7-PLA2G4B | 9 | 6 | MCL31 | 9 of 9 (12 to 132 reads) |
| | MAPKBP1 | 9 | 0 | none | 4 of 9 |
| | SPTBN5 | 9 | 0 | none | 8 of 9 (2 to 69 reads) |

KB3781 family MCL11 has 8 members: 7 match a Yoo JMJD7-PLA2G4B locus (6 on chr1, the ancestral chr16 locus), one member (chr1:19,235,748-19,243,916) lies 2.5 kb from Yoo's chr1:19,225,471-19,233,170 without overlapping it, and the two single matches of MAPKBP1 and SPTBN5 come from one 24-exon member (chr1:15,783,756-15,912,704) that runs through three adjacent gene copies (a readthrough locus). So the unit's expressed gene is recovered as one family with 8 of Yoo's 9 loci (KB3781) and 6 of 9 (OR6737), and the MAPKBP1 and SPTBN5 copies are not recovered because the libraries carry few or no reads at most of their loci (the loci with at least 2 reads are in the last column); reads are present at every JMJD7-PLA2G4B locus, so the 1 and 3 missed loci are losses of the BASE snapshot, not of data. The prediction of Amendment 3 (at least 5 of 9 MAPKBP1 loci covered) failed (1 of 9 and 0 of 9). The member-level count differs from the locus-level count by design: in KB3781 the 24-exon member overlaps four Yoo loci (two JMJD7-PLA2G4B, one SPTBN5, one MAPKBP1), so 8 Yoo loci are covered by 7 members (the 8th member is the one 2.5 kb away). The object is the copy rows (transcript spans) of the BASE `copies.tsv`, as registered; a post hoc sensitivity with the locus extents (`--use-locus`) raises the JMJD7-PLA2G4B count to 9 of 9 (KB3781) and 8 of 9 (OR6737) and leaves MAPKBP1 and SPTBN5 at 1 of 9 in both libraries, so the conclusion does not move. Reproduce: `python3 -B bench/hierarchy/yoo_concordance.py y2 --unit-loci winloci_data/yoo_2025/yoo_GGO_unit_loci.tsv --copies /mnt/linuxdisk/tmp/rustle_figures/rt_arms/<library>/<library>.BASE.fam.copies.tsv --out PREFIX` (copies.tsv sha1 64abad727435 for KB3781, a2de88e1097a for OR6737; unit-loci sha256 19a30600139d667c...; tests `bench/hierarchy/test_yoo_concordance.py`, 15). Not done: the current-default (f1v2) families on the two contigs, and a testis library of Yoo's own animal.

## Y3, PSMA5 (count only)

The gorilla RefSeq annotation (`GGO_genomic.gff`) has the ancestral PSMA5 (NC_073224.2 = chr1:130,657,639-130,682,285) and three 'proteasome subunit alpha type-5-like' gene records upstream of it on the same chromosome (LOC129525664 at 116.67 Mb, LOC129525329 at 125.73 Mb, LOC129532613 at 129.28 Mb), which is Yoo's ancestral copy plus 3 duplicated copies upstream. A fourth record, LOC109028113 (NC_073231.2, 1 kb, 'proteasome subunit alpha type-5'), is a short retro-like record that Yoo does not count. The prediction (the ancestral gene and at least one duplicate as separate records) held (3 of 3). The coordinates of Yoo's PSMA5 copies are in no table, so the comparison is by count and chromosome only.

## Y4, genome-wide E-table against Tables VIII.40 and VIII.35

Not yet: it needs the E-table (Phase E).
