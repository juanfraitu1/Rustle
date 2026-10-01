# Pre-registration: does OUR homology method find Soto's families when given Soto's copy numbers? (2026-09-30, KEY=ourhomology)

Written before the all-vs-all alignment was run or any family was built. DNA only. Human CHM13 v1.0, Soto's 2,334-gene table (S1C).
Nothing in `src/` or `tools/` changes. Script: `bench/soto_m2/soto_m2_our_homology.py`, written after this file is committed.

## 1. Question

Soto's families = a homology graph (their exon map-back) cut by copy number (per-pair MAD < 1 on famCN). Every reproduction so far used
**their** homology graph (map-back + S1C famCN: 479 / 491 exact, ARI 0.9698). This test swaps in **our** homology method and keeps
everything else Soto's: their genes, their famCN values (S1C "Median famCN", used as given — "imputed" into our pipeline), their pair
rule and grouping. If our homology finds the same families, the difference between the two methods is not in finding homologs.

## 2. Our homology method (frozen)

- Loci: each S1C gene's span (CAT v4, `winloci_data/soto_replication/cat_v4_soto2334.gff3`), sequence from CHM13 v1.0
  (`t2t-chm13-v1.0.fa.gz`); exons = the union of the gene's CAT exons.
- All-vs-all: `minimap2 -x asm20 -c -X -N 50 -p 0.1 --secondary=yes` (the command `mcl_families --from-gtf` runs).
- Families: `mcl_families` (binary built 2026-09-30 from main, copied before use) with the shipped DNA-level settings
  `--min-exonic-bp 1 --min-shared-exon-frac 0.60` and every other option at its default (min-cov-shorter 0.70, inflation 2.8,
  prune 1e-9, min-size 2), `--dump-pairs` to obtain the homology edges inside each family.

## 3. Arms (all scored with `soto_m2_families.classify`, the scorer behind the meeting page's numbers)

- **H0 (reference):** Soto's exon map-back edges (`bench/soto/shared_exons_5154_exon_mapback.tsv`).
- **H1 (ours):** the within-family edges of our MCL families (`<out>.pairs.tsv`).
- For each: sequence only (no copy-number cut), S1C famCN with Soto's pair rule (MAD < 1), and our recomputed famCN (secondary).
- Measures: ARI and exact families against S1C (all; dev / held-out by the frozen split, a family's half = the half of most of its
  clean members, as in `soto_m2_loosen.py`), nesting (Soto families with >= 2 clean members inside one sequence-only component),
  bipartite sensitivity / precision.

## 4. Decision rule (primary: H1 with S1C famCN, held-out half)

- **FINDS THEM:** held-out ARI(H1) >= held-out ARI(H0) - 0.05.
- **PARTIAL:** held-out ARI(H1) >= held-out ARI(H0) - 0.15.
- **DOES NOT:** otherwise.
The 0.05 / 0.15 margins are fixed here without data; they are not fitted. Exact families, nesting and bipartite are reported beside it.

## 5. Seen before (disclosed)

- H0 numbers: sequence only 0.7307 / 345 exact, S1C famCN 0.9698 / 479 (held-out 0.9681 / 263), our famCN 0.9277 / 411; nesting
  440 / 444 on map-back edges, 346 / 444 on SEDEF-projected edges.
- Our pipeline's catalogs (different gene sets, RNA and DNA levels) nest 80-90% of Soto's families (`bench/SOTO_AS_A_REFINEMENT.md`).
- Nothing has been computed with our homology on Soto's 2,334 loci.

## 6. Result

(Filled in after the run, below this line, without editing anything above.)
