# How faithfully does each Soto family's genes get reproduced by our de novo loci?

**§6w5, 2026-09-22.** User: *"we have an analysis of how close our predicted loci are vs ground truth for
NPIP... can we run this for all the families in Soto and find the worst approximations."*
Tool `bench/soto_family_locus_fidelity.py`, run on the **genome-wide** de novo assembly built in §6w4
(`ours_genome.gtf`, 264,996 transcripts, A119b) — the first time this question could be asked outside chr16.

⚠**HUMAN A119b only.** Never pooled with gorilla.

## Population and guards

Soto genes carrying a Family ID: **1,078**. Matched one-to-one to a de novo locus by exonic overlap
(greedy, a locus already claimed by a better gene cannot be reused and the loser is a MISS):
**765 = 71.0%**. Of those, **258 (33.7%) sit on a locus covering ≥2 annotated genes** and are reported
SEPARATELY — a readthrough-fused locus's "boundary error" is over-merge, not boundary error (§6v8), and
pooling them reproduces the exact confound §6w3's population filter existed to avoid.

## Headline: the limiting factor is EXPRESSION, not assembly

| biotype | n clean | median exonic Jaccard | med \|d5\| | med \|d3\| | **match rate** |
|---|---|---|---|---|---|
| protein_coding | 287 | **0.621** | 755 | **18** | 79.3% |
| transcribed_unprocessed_pseudogene | 106 | 0.493 | 480 | 360 | **85.9%** |
| unprocessed_pseudogene | 74 | 0.315 | 2,807 | 962 | **46.6%** |
| processed_pseudogene | 15 | 0.225 | 5,998 | 159 | **42.1%** |
| lncRNA | 17 | 0.222 | 1,116 | 887 | 85.4% |

⭐⭐**A clean natural experiment sits inside this table.** `unprocessed_pseudogene` and
`transcribed_unprocessed_pseudogene` are the SAME sequence class; they differ only in whether the gene is
transcribed — and the match rate splits **46.6% vs 85.9%** with median Jaccard **0.315 vs 0.493**. The
genes we fail to represent are overwhelmingly the ones with no RNA to represent them. **This is not an
assembler defect and no node-construction change can reach it.**

⭐**§6w3's "the 3′ end is exact" is a PROTEIN-CODING statement, not a general one**: median |d3| is
**18 bp** for protein_coding and **159–962 bp** for every other class. The 5′/3′ asymmetry the NPIP work
described generalises only where the gene is well expressed.

## The worst families

| family | n | medJac | med\|d5\| | med\|d3\| | width ratio | members |
|---|---|---|---|---|---|---|
| ID_481 | 3 | **0.027** | 1,447 | 529 | **0.24** | UBE2Q2P11, UBE2Q2P12, UBE2Q2P6 |
| ID_339 | 2 | 0.068 | 9,955 | 3,180 | **12.65** | DEFB109A, DEFB109B |
| ID_98 | 2 | 0.086 | 18,204 | 13,668 | 6.85 | ENPP7P1, ENPP7P12 |
| ID_451 | 2 | 0.100 | 11,420 | 29,494 | 0.11 | SLC25A24P1, SLC25A24P2 |
| ID_352 | 2 | 0.120 | 43,794 | 33,632 | 1.91 | FAM21EP, FAM21FP |

Best, for contrast: **RSPH10B/B2 0.990 · NOMO1/2/3 0.989 · TP53TG3B/D 0.974 · ZNF322/P1 0.930 ·
SPATA31A×6 0.930** — all protein-coding, with med |d5| 6–230 bp and med |d3| 2–5 bp.

## Two opposite failure modes, one mechanism

- **ID_481 UBE2Q2P\*** — loci capture only **11–28%** of the gene, on **2, 3 and 2 reads**.
- **ID_339 DEFB109A/B** — loci are **10.0× and 15.3× TOO BIG**, on 3 and 18 reads (DEFB109A's locus
  starts 12,838 bp past the gene's own 5′ end).

Across all 507 clean genes the two modes separate cleanly **by depth, in opposite directions**:

| reads | n | ratio < 0.5 (under-capture) | ratio > 2 (over-extension) |
|---|---|---|---|
| 2 | 142 | **28.2%** | 12.7% |
| 3-4 | 100 | 17.0% | 12.0% |
| 5-9 | 70 | 15.7% | 12.9% |
| 10-29 | 116 | 11.2% | 19.0% |
| 30-99 | 47 | 12.8% | **25.5%** |
| ≥100 | 32 | 12.5% | 15.6% |

⭐⭐**This is §6w3's "k is a fixed RANK, so it is a moving QUANTILE" reproduced on an independent metric
(width, not the 5′ boundary) and an independent population (Soto families genome-wide).** The boundary is
the k-th most extreme read end with **k=2 fixed**: at n=2 that is the *inner* read, so the locus shrinks
(under-capture); at n=30-99 it sits deep in the tail, so one long molecule drags the locus out
(over-extension). Median exonic Jaccard is nearly flat with depth (0.504 at ≤4 reads vs 0.548 at ≥10)
**precisely because the two errors trade off.**

⚠**47.7% of matched Soto genes rest on ≤4 reads**, so the low-depth regime is not a tail here, it is half
the population.

## Reproduce

```sh
python3 bench/locus_probes.py soto-fidelity --gtf ours_genome.gtf \
  --gff chm13v2.0_RefSeq_full.gff.gz --soto bench/soto/soto_famCN_S1C.tsv --out soto_fid
```

> **Generator (2026-09-22 consolidation):** `python3 bench/locus_probes.py soto-fidelity ...` — the original `bench/soto_family_locus_fidelity.py` was folded in verbatim and verified identical on its documented inputs (§6z3); the register rows above cite this file.
