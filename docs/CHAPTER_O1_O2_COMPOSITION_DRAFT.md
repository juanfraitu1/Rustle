# O1 → O2 by composition: families of transcribed loci, and assignment among exactly those loci (draft, 2026-09-05)

Compact version; prose draft verbatim at `archive/2026-11/CHAPTER_O1_O2_COMPOSITION_DRAFT.md`.

## Current state (2026-11-01)
- This is the 2026-09-05/06 draft on the 3-contig substrate; its numbers are of that date. The composition (O1 units -> O2 assigns among exactly them) still holds; the thesis model is the variation graph (sketch families = member discovery only).
- O2's domain since 2026-10-19: AS-tied multimapping reads only (rival within 10 AS). At sample scale (LEDGER_O2 `O2 AMBIGUOUS DECOMPOSITION`): HG002 17,805 of 141,265 tied reads assigned; discards are outside-catalog origins and junk.
- Newer rows: 1477 (`--tie-ms` worth porting), 1483 (`--tie-ms` port FAIL, opt-in), 1481 (O2 LEVERS, amended to a 57-family panel), 1484 (O2 = pileup + caller; bcftools agrees 54.6% as run / 93.9% with MAPQ reset). Family-stage rules 1469-1476 FAIL or opt-in; defaults unchanged.

## 1. The two objects, and the one that is not the thesis's
A duplication block (SEDEF) is not a gene family (a set of loci sharing a transcribed core); on 16p LCR16a (NPIP) and LCR16u (SMG1 pseudogenes) differ. The catalog defines the second and must say why it did not build the first. (full: archive/2026-11/CHAPTER_O1_O2_COMPOSITION_DRAFT.md:3)

## 2. O1: the family as a set of units
Node = locus (overlapping exon unions, core linked by SD pairs to >= half the family, read-supported chain: >= 3 reads). Edge = asm20 homology (identity >= 0.70, coverage >= 0.30, >= 300 bp, exon-on-exon), inside a measured plateau (walls at 0.85 / 0.60 / 1 kb). Family = MCL (inflation 2.8, prune 1e-9, min 2) proposals certified by the core rule, SD, repeat library, CHM13 landing. Old γ-quasi-clique definition: precision 0.974 vs 0.954, recall 0.874 -> 0.940 (`O1_DEFINITION_SWITCH.md`). Certificates per unit: SD depth, repeat fraction, nearest-paralogue identity, block (`blocks.tsv`), CHM13 landing. (full: archive/2026-11/CHAPTER_O1_O2_COMPOSITION_DRAFT.md:11)

## 3. What the certificates decided on the 3-contig substrate (46 clusters with n ≥ 5)
39 SD-corroborated (3 anchored), 3 repeat artefacts (L1, ERV1, L1), 1 element-welded (MCL4); genome-wide at prune 1e-9 the tandem array is recovered. (full: archive/2026-11/CHAPTER_O1_O2_COMPOSITION_DRAFT.md:41)

## 4. O2: assignment among exactly O1's units
Domain (user, 2026-10-19): AS-tied reads only. Each molecule aligned splice-aware to every candidate locus; competitors rejected at `min_p < α/(n−1)` on shared distinguishing columns, else tie/abstain; origin certificate abstains when no candidate explains the read. 35 paired families: 56% assigned, 9,236 of 9,236 MAPQ-60 agreeing; 62 valid NPIP anchors (55 of 117 were annotation gaps), no wrong confident call on the spliced path. (full: archive/2026-11/CHAPTER_O1_O2_COMPOSITION_DRAFT.md:49)

## 5. The result on NPIP, and catalog-wide
NPIP (32 loci): 0 wrong on 62 anchors, MAPQ-60 201/201; of 433 low-MAPQ reads 115 assigned, 64 K = 0, 254 unexplained; on NPIP, O2 is an abstention certificate. NPIPA/B is a human 16p label, not a cut in this catalog (§6ew). Catalog-wide (35 families): reads-following units 28% assigned / 8% tied / 63% ambiguous, agreement 99.4% (4,674/4,704); identity does not predict abstention (ρ −0.05; row 688). (full: archive/2026-11/CHAPTER_O1_O2_COMPOSITION_DRAFT.md:68)

## 6. What is a limitation and what is a defect
Limitations: copy-inconsistent annotation boundaries, library-silent elements (MCL4), slice-conditioned cohesion, SD corroborates the duplication not the gene (row 669). Defects fixed (rows 689-691): PSV read-filter (ZSCAN5, 11 wrong; off is better), SD-hull clipping (ZNF569-like, 190 wrong; chain follows reads). Open then: secondary records scored as observations (row 691). (full: archive/2026-11/CHAPTER_O1_O2_COMPOSITION_DRAFT.md:86)

## 7. The sentence
DNA defines the edge, RNA defines the node where the annotation is absent, SEDEF and the repeat library audit both, and O2 assigns — or, on NPIP, declines to assign — among exactly the units O1 defined. (full: archive/2026-11/CHAPTER_O1_O2_COMPOSITION_DRAFT.md:101)

## 8. Figures
Figs 1-4 `docs/figures/fig_sweep_status.png`, `fig_sweep_agreement.png`, `fig_znf569_locus.png`, `fig_npip_anchors.png`, rendered by `bench/chapter_figures.py`. Update 2026-09-06 (§6fh–§6fj): candidate set = every cluster locus; sole candidate = assigned; 79.8% assigned, agreement 18,994/18,995, NPIP anchors 16/16; `n_reads_soft`, `--dump-star`. (full: archive/2026-11/CHAPTER_O1_O2_COMPOSITION_DRAFT.md:105)
