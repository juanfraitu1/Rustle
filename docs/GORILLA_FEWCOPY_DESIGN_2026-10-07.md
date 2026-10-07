# Gorilla OR6737 testis few-copy class test: design sweep result (2026-10-07)

**Status: DESIGN ONLY, nothing of the pipeline was run on gorilla. The class that the inputs define is too small for the drafted bar: |K| = 3, so rule 6 of `docs/PREREG_fewcopy_class_DRAFT_2026-10-06.md` (n_reachable < 8 = UNDERPOWERED, no verdict) applies before any arm exists.**
User direction (2026-10-06): "lets try with the gorilla testis". Workflow wf_fbf8afe6-9cb (exposure ledger and symbol-projected truth in parallel, then criteria and class K blind to every pipeline product, then an adversarial review that recomputed 5,050 values). Outputs `/mnt/linuxdisk/tmp/fewcopy_gorilla_2026-10-07/` (tables, scripts, the four sweep results). Frozen gorilla criteria: `bench/fewcopy/criteria_gorilla.md`, sha256 `2b17d529dc62710ade954e67363118cd9685e31780b2c2cb612f52196c19086e` (frozen 00:07 PDT, first count 00:12 PDT, hash re-checked by me and by the reviewer; human criteria 89363117 unchanged).

## Label to use

**Not 'held-out'.** "Second-species, inputs-frozen (outcome-blind) class check on a partly spent library (gorilla OR6737 testis)"; KB3781 is reported separately as dev (same assembly, annotation and truth source). OR6737 is DEV for the seeding rule (NC_073244.2), F1 / F1v2 (both gorilla libraries) and the DAZ2 gate; held-out and spent for strict junctions + retained-intron ratio 10 (26 contigs), F1 (minus NC_073244.2) and R_J (NC_073234.2); the RefSeq annotation (the universe and the truth source) fixed the MCL constants; no contig is unexposed. Missed by the ledger and found by the review: NC_073230.2 + NC_073228.2 were the held-out substrate of Addendum AH (09-14, CF343 sits there), the MHC (C4A / C4B) is in the 06-25 reference-absent catalog, ZNF600 is in the 09-05 O2 / O3 gap-locus rows.

## Truth: human Compara families projected by exact RefSeq symbol

Human Ensembl Compara 2-4-gene families (377) projected onto the gorilla RefSeq annotation (GCF_029281585.2, individual KB3781) by exact symbol of a protein-coding gene record: **342 of 864 members project 1:1 (39.6%); 90 families whole, 135 partial, 152 none**; the 522 absent members are 305 with LOC-named protein-coding candidates by name evidence and 217 without. All rows recomputed against the GFF by the reviewer, 0 mismatches. Selection effect: old dispersed pairs project (whole in 54.0% of human cross-contig families against 17.8% of same-contig ones); the recent duplicates (LOC-named) do not. This is a third truth method, neither Compara gorilla nor the 07-29 projection of register row 861; it is human truth projected, not gorilla truth.

## The class K under the frozen criteria

Funnel (families): 377 -> **90** whole projection -> 38 (C1 annotation) -> 30 (C2 not entangled) -> 29 (C3 separate loci) -> 7 (C4 identity [0.90, 0.999), coverage >= 0.5) -> **3** (C5 >= 3 exact-chain reads per member). Independent passes among the 90: C1 38, C2 73, C3 86, C4 33, C5 17. Thresholds are the human ones; the forced changes are the single RefSeq annotation (no CAT mapping tier) and the longest mRNA as reference transcript.

| family | members | contig | identity | min reads | verdict on the case |
|---|---|---|---|---|---|
| CF343 | STEAP1, STEAP1B | NC_073230.2 (68.7 Mb apart) | .987 | 28 | the only clean O1 case: dispersed pair, no read ties |
| CF315 | C4A, C4B | NC_073229.2 (MHC) | .998 | 134 (same reads at both copies, P1 65 / 70) | an O2 tie control: cross-link .91 / .90, read divergence ~ the inter-copy .0020; hyperpolymorphic locus |
| CF236 | ZNF600, ZNF611 | NC_073244.2 | .947 | 23 | unreachable by listed size: ZNF808 / ZNF888 / ZNF578 are unlisted relatives in a 43-gene KRAB-ZNF cluster; its human outcome ('merged') is dev-known |

Near misses (12, fail one criterion): C4 x8 (FCGR2A / FCGR2B coverage .30, ALG10 / ALG10B .47, RABL2A / RABL2B .38, ZNF155 / ZNF230 identity .78, DGCR6 / DGCR6L .87, BEX1 / BEX2 .87, EOLA1 / EOLA2 identity .9997, ZNF195 / ZNF429 cross-contig), C5 x4 (KRT37 / KRT38, USP32 / USP6, PF4 / PF4V1, UGT2B15 / UGT2B17 one read short).

## What follows from it

- The drafted class-level test cannot return a verdict on gorilla testis with these inputs: |K| = 3 (one clean). Admitting every LOC-named family at the observed 3 / 90 rate would give about 7 before the reachability control.
- Read support is the binding constraint, not annotation: **17 of the 90 whole families have every member read-supported (15 pairs, 2 triples)**, but 7 of the 17 are KRAB-ZNF families on NC_073244.2 (ZNF181 / ZNF302, ZNF585A / B, ZNF224 / ZNF225 / ZNF284, ZNF600 / ZNF611, ZNF195 / ZNF429, ZNF324 / ZNF324B, ZNF155 / ZNF230), whose listed size hides relatives (not independent, exact recovery unreachable by listed size). The other ten: C4A / C4B, TCEAL2 / TCEAL4, BEX1 / BEX2, STEAP1 / STEAP1B, EOLA1 / EOLA2 (identical), ALG10 / ALG10B, FCGR2A / FCGR2B, DGCR6 / DGCR6L, RABL2A / RABL2B, APOBEC3C / D / F.
- Risks specific to gorilla: cross-individual reads (OR6737 reads on the KB3781 assembly; private variants drop exact-chain reads, so P2 is a lower bound), RefSeq-only annotation (C2 'not entangled' may mean under-annotated: the human CAT contrast fails C4A / C4B at 520 bp shared with antisense lncRNAs and STEAP1B at 986 bp), one movie and one library, 57 of 342 projected members on chrX, Gnomon models may cite long SRA reads (if OR6737 is among them, C5 is partly circular; unverified).
- Overlapping genes on gorilla (input-only, exploratory, not frozen): 2,202 gene pairs share >= 100 exonic bp; 1,362 have a valid chain in both genes; 423 have both genes with >= 3 exact-chain reads (35 on the same strand); at >= 1 kb shared 59 (10 same strand). This is a real-read test bed for the overlapping-gene question (the lab's StringTie, FLAIR and isoseq outputs for OR6737 exist genome-wide), independent of the few-copy class.

## Unsupported claims withdrawn by the review

'All three K families have 0 unlisted relatives' is name-based (true by sequence for CF315 and CF343, false for CF236); 'family-level exactness of 2-4-member families was never scored on any gorilla contig' is unverified (r1232 and r1237-r1245 are family-level default runs on OR6737 testis, not 2-4-member Compara families); the quoted '46.6% of genes are LOC ids' does not reproduce (55.7% of gene records, 23.4% of protein-coding genes); 'no better truth source on disk' misses `GGO_sedef_final.bed` (DNA segmental-duplication calls, a sequence-based label for LOC-named recent duplicates, with a circularity caveat); CF315 ranks first only because min P2 counts the same 135 reads at both copies.
