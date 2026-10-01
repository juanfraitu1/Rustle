# Pre-registration: does our genome-only self-alignment + quasi-clique grouping reconstruct Soto's families? (2026-09-30, KEY=quasiclique)

Written before the self-alignment was run or any family was built. DNA only. Human, Soto's 2,334-gene table (S1C). Nothing in `src/` or
`tools/` changes. Scoring script: `bench/soto_m2/soto_m2_quasiclique.py`, written after this file is committed.

## 1. Question

The July detection page reported that our genome-only mode (`gw_family_catalog --from-genome`: self-alignment, then the γ-quasi-clique
homology core) places 361 of 362 members of 83 Soto families in some multi-copy family. That is member detection on 83 of 491 families.
**Does the same method, run on all of Soto's genes, reconstruct Soto's families as groups?**

## 2. Method (frozen)

- **Windows = Soto's 2,334 gene loci** (CAT v4 gene spans), the mode's own design ("each window is a locus to group … this is Soto's
  setup"). The page's coordinates are CHM13 v1.0; the run uses CHM13 v2.0 (`winloci_data/Reference/chm13v2.0.fa`, whose prebuilt splice
  index `/home/juanfra/human_val/chm13v2.0.splice.mmi` is passed as `RUSTLE_PROJECT_MMI`), so each gene span is lifted v1.0 -> v2.0 by
  sequence: the gene's v1.0 sequence must occur at the shifted position in v2.0 with identical bases (offset searched within ±1 Mb);
  genes that fail are listed and left out of both truth and prediction.
- **Run:** `gw_family_catalog --from-genome windows.bed --fasta chm13v2.0.fa` with every other option at its default (binary built
  2026-09-30 from main, copied before use). The self-alignment (`minimap2 -c -x splice -N 50 -p 0.01` of the windows against the genome)
  is run in query chunks with that exact command and supplied to the binary by a wrapper that first checks, with `cmp`, that the binary's
  query FASTA equals the chunks concatenated; the all-vs-all of the representatives goes through `tools/mm2_shard.sh` (byte-identical
  to one run by its own check).
- **Families -> genes:** each window's representative is its gene; a gene in no family of >= 2 members is a singleton. Discovered extra
  loci (paralogs outside every window) take part in the grouping but are not scored.

## 3. Arms and scoring (`soto_m2_families.classify`, the scorer behind the page's numbers)

- **Q:** edges = every pair of genes in the same quasi-clique family (the families are γ-quasi-cliques, so near-complete). Scored with no
  copy-number cut, with S1C famCN under Soto's pair rule, and with our famCN.
- References from KEY=ourhomology: H0 (Soto's map-back) held-out ARI 0.9681 with S1C famCN; H1 (our MCL homology) 0.8844.
- Measures: ARI and exact families (all, dev, held-out by the frozen split), nesting (Soto families with >= 2 clean members inside one
  family, no copy-number cut), bipartite sensitivity / precision, genes placed in a family of >= 2 (member detection).

## 4. Decision rule (primary: Q with S1C famCN, held-out half; same margins as KEY=ourhomology)

- **RECONSTRUCTS:** held-out ARI(Q) >= 0.9681 - 0.05.
- **PARTIAL:** held-out ARI(Q) >= 0.9681 - 0.15.
- **DOES NOT:** otherwise.

## 5. Seen before (disclosed)

KEY=ourhomology (MCL on Soto's loci: 249 / 491 exact, held-out ARI 0.8844, nesting 257 / 444; 158 non-nested families have members
with no edge, many short genes); July DNA-mode detection (361 / 362 members of 83 families, raw precision 47%); the 106 bp gene
AC243829.6 needs a minimap2 score floor below asm20's 200 to align at all. No `--from-genome` run on Soto's 2,334 loci has been made.

## 6. Result

(Filled in after the run, below this line, without editing anything above.)

**Amendment 1 (2026-09-30 23:15, before any family was built).** The genome-wide self-alignment (paralog discovery) needs ~40 min
(1 of 14 chunks done when stopped for a meeting). A **windows-only arm (W)** is run first: the same command, but the self-alignment
returns no hits (wrapper `PROJ_EMPTY=1`), so the representatives are exactly Soto's 2,334 loci and the γ-quasi-clique grouping runs on
them alone. Arm W is reported as a deviation and decides nothing; the full arm of section 2 still decides when it is run. Everything
else (scoring, decision margins) unchanged.
