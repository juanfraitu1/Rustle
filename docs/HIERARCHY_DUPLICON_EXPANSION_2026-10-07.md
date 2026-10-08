# Hierarchy study: results (2026-10-07). Y block, concordance with the Yoo et al. 2025 gorilla expansions

**Status: DESCRIPTIVE; no verdict class.** Protocol: `docs/PREREG_hierarchy_duplicon_expansion_2026-10-07.md`, Amendment 3 (git 8fdff163, committed before the Y2 and Y3 numbers existed; Y1 was looked at before and is labelled spent). Results of T1 (Phases R0 to R2) and T2 (Phase E) will be appended to this file under their own headings. Truth files: `/mnt/linuxdisk/home/juanfraitu/winloci_data/yoo_2025/` (workbook sha256 865d4e868568b32d...; `table_38.tsv`, `yoo_GGO_table35_novel_paralogs.tsv`, `yoo_GGO_table40_lineage_specific.tsv`, `GGO_chr_map.tsv`, `yoo_GGO_unit_loci.tsv`, `y1_lrpap1_concordance.tsv`). Match rule: same chromosome, overlap of at least 50% of the shorter locus. Gorilla only; one animal (Jim, mGorGor1).

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
