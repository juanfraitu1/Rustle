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

Run on 2026-09-30 after this file was committed (`ffc7fb13`, sha1 `94dd0a83a411ea099eede3c846c39198ee08a41b`). All-vs-all 164,616 PAF
records (4 query chunks, 2.5 min, heavy lock); `mcl_families` (binary sha1 e4cc13b9) with the shipped DNA settings: 1,819 nodes, 11,031
edges, 26,778 pairs dropped for no exonic evidence, 3,053 for a shared-exon fraction < 0.60, 93 annotation records folded into 86 loci;
380 families covering 1,709 of 2,334 genes; 10,045 within-family edges -> 10,731 gene pairs. Scoring: `soto_m2_our_homology.py`, 6 s.
H0 reproduces every earlier number exactly (0.7307 / 345; 0.9698 / 479, held-out 0.9681 / 263; ours 0.9277 / 411).

| homology | copy numbers | ARI all / dev / held-out | exact all (dev / held-out) | nesting | bipartite sens / prec |
|---|---|---|---|---|---|
| H0 Soto map-back | none | 0.7307 / 0.6418 / 0.8693 | 345 | 440/444 | - |
| H0 Soto map-back | S1C | 0.9698 / 0.9708 / 0.9681 | 479 (216 / 263) | 440/444 | 0.980 / 1.000 |
| H0 Soto map-back | ours | 0.9277 / 0.9227 / 0.9343 | 411 (177 / 234) | 440/444 | 0.935 / 0.975 |
| H1 ours | none | 0.7648 / 0.7146 / 0.8242 | 186 | 257/444 | - |
| **H1 ours** | **S1C** | **0.8235 / 0.7675 / 0.8844** | **249 (113 / 136)** | 257/444 | **0.690 / 0.966** |
| H1 ours | ours | 0.7856 / 0.7323 / 0.8446 | 221 (91 / 130) | 257/444 | 0.656 / 0.947 |

**VERDICT (section 4): PARTIAL.** Held-out ARI 0.8844 against the map-back reference 0.9681 (within 0.15, not within 0.05).

Descriptive, after the verdict: our families rarely over-merge Soto's (bipartite precision 0.966) but miss part of them (sensitivity
0.690). Of the 185 Soto families (>= 2 clean members) not inside one of our components, 158 have members with no edge at all in our graph
and 27 are split across our components. The unlinked members are unprocessed pseudogenes (126), protein-coding (76), processed
pseudogenes (59), miRNAs (51), transcribed unprocessed pseudogenes (39) and lncRNAs (34). 8,657 of the 12,231 map-back pairs are also our
edges; we add 2,074 pairs the map-back does not have.

**Reading.** The two homology definitions differ in what counts as related: Soto links two genes that share one >= 98%-identical exon;
our shipped DNA edge asks the alignment to cover a substantial part of a gene (30% of the longer gene's exons, or 70% of the shorter's)
over >= 300 bp, with >= 60% of exons shared. Short genes (miRNAs) and genes that share only a piece fall below ours. Given Soto's copy
numbers, our homology recovers most of Soto's structure without merging their families, but not the single-exon links that make their
families larger.
