# Family-merge threshold: minimizer-Jaccard sensitivity sweep

**Question (advisor):** "What is the threshold to merge graphs? Is there an exon
similarity metric in place?"

**Answer:** Yes. The exon-similarity metric is **minimizer-Jaccard** (k=15, w=10)
computed per exon between paralog copies, with merge bar
`RUSTLE_VG_FAMILY_MERGE_JACCARD` (**default 0.30**).
`vg_family/family_graph.rs`: `minimizers()` (line ~226),
`refine_by_minimizer_jaccard()` (line ~243), `merge_singletons_by_sequence()`
(line ~293). Two exons go in the same `ExonClass` iff
`|min(A)∩min(B)| / |min(A)∪min(B)| ≥ bar`.

This sweep characterises that exact metric (calls the production `minimizers()`)
to locate where the bar should sit. Harness: `#[ignore]`d, env-guarded test
`merge_jaccard_sensitivity_sweep` in `family_graph.rs`. Reproduce:

```
MERGE_SWEEP_DATA=/tmp/merge_sweep_exons.tsv \
  cargo test -p rustle merge_jaccard_sensitivity_sweep -- --ignored --nocapture
```

Inputs: 6 real LRPAP1 exons (NC_073227.2, lengths 119–471 bp, mean ~240) as the
homolog source; 6 real exons from two unrelated genes (NC_073236.2) as the
non-homolog control. Divergence axis is controlled (substitution-only SNPs,
seeded, 40 trials/rate); the **sequences are real**.

## Result 1 — minimizer-Jaccard vs sequence divergence (homologs)

| divergence | mean SNPs | mean Jaccard | % still merging @0.30 | @0.10 |
|-----------:|----------:|-------------:|----------------------:|------:|
| 0.0%       | 0.0       | **1.000**    | 100% | 100% |
| 0.5%       | 1.3       | 0.867        | 100% | 100% |
| 1.0%       | 2.4       | 0.747        | 100% | 100% |
| 2.0%       | 4.7       | 0.578        |  96% | 100% |
| 3.0%       | 7.1       | 0.444        |  82% | 100% |
| **5.0%**   | 12.3      | **0.283**    |  37% |  96% |
| 8.0%       | 19.3      | 0.156        |   7% |  74% |
| 10.0%      | 24.1      | 0.100        |   2% |  40% |

The **0.30 bar is crossed at ~5% divergence (~10–12 SNPs over a ~240 bp exon)**.
This validates the code comment at `family_graph.rs:304` ("~9 SNPs in a 300 bp
exon push Jaccard to ~0.3") — directionally exact; measured crossing is ~10–12
SNPs. k=15 minimizers are SNP-sensitive: a single SNP destroys every k-mer (and
candidate minimizer) spanning it, so Jaccard decays steeply with divergence.

## Result 2 — non-homologous baseline (false-merge ceiling)

51 unrelated exon pairs (cross-gene + family×unrelated): **mean Jaccard 0.0000,
max 0.0000**. Unrelated exons share **zero** k=15 minimizers. The false-merge
ceiling is a hard floor at 0 — there is no value above which unrelated exons
start merging within the tested set.

## Result 3 — threshold sweep (recall vs false-merge)

| bar  | recall@0.5% | @1% | @2% | @3% | false-merge |
|-----:|------------:|----:|----:|----:|:------------|
| 0.30 | 100% | 100% | 96% | 82% | clean |
| 0.20 | 100% | 100% | 99% | 95% | clean |
| 0.15 | 100% | 100% | 100% | 96% | clean |
| **0.10** | 100% | 100% | 100% | 100% | **clean** |
| 0.05 | 100% | 100% | 100% | 100% | clean |

## UPDATE (2026-06-06) — the operative bar was 0.05; raised it to 0.30 (+2 RBMY)

The separation analysis above (and a first indel-inclusive read) chased the wrong
number. Instrumenting the resolver (`RUSTLE_MJ_DEBUG`) revealed that the family-
merge constant I was tuning (`vg.rs:480/637`) is a **dead path** for assembly: the
real VG assembly path is **pipeline.rs:10955**, which **hardcoded the cross-copy
merge bar at 0.05** (refine 0.0). So the operative default was 0.05 — far too low.

**Shipped change:** pipeline.rs:10955 `0.05 → family_merge_jaccard()` (default
**0.30**, env-overridable). Benchmark, old(0.05) vs new(0.30) over the 15-locus
panel: **RBMY (LOC129530243) 5→7 TP**, all 14 other loci unchanged, total
**64→66 TP / +1 pred** (ΔTP +2 at essentially no precision cost, 0 regressions).

The substitution-only separation suggested *lowering* to ~0.10; the indel model +
the assembly benchmark show the opposite — the bar was already **too low** and
needed RAISING.

**Indel-inclusive curve (copy-vs-copy, the quantity union-find thresholds):**

| divergence (incl. indels) | Jaccard | Containment |
|---:|---:|---:|
| 1% | 0.50 | 0.67 |
| 2% | 0.28 | 0.44 |
| 3% | 0.16 | 0.28 |

**Assembly benchmark — threshold ladder, recovered RefSeq transcripts (TP/pred):**

| locus | 0.40 | 0.30 | 0.25 | 0.20 | 0.10 |
|---|---|---|---|---|---|
| **RBMY** (LOC129530243) | **8/13** | 7/12 | 7/12 | 5/11 | 5/11 |
| LOC129523503 / LOC101124683 / LOC115932683 | flat | flat | flat | flat | flat |

Recovery on the flagship dispersed paralog is **monotone in the bar — higher is
better.** Lowering to 0.10 *regresses* RBMY (7→5) and improves nothing; raising
toward 0.40 helps (8). 14/15 panel loci are unchanged (the cross-copy merge only
fires for *dispersed* paralogs).

**Why:** the merge objective is **not pure recall.** A lower bar over-merges
near-identical-but-DISTINCT copies — union-find folds a copy on a *single*
above-bar edge and chains transitively, and the completion gate (vg.rs:2803)
then fabricates phantom exons on coverage asymmetry with no similarity-magnitude
guard — collapsing copy-distinct isoforms. An adversarial sweep also found
false-merge modes that appear only **below ~0.25** and that 0.30 refuses:
embedded shared TE/domain (J≈0.22), poly-A/low-complexity cores (J≈0.17). Pure-
satellite exons collapse at J≈1.0 at *any* bar (a metric limitation, mitigated
upstream by repeat-masking). All encoded as permanent guards in `family_graph.rs`.

**Conclusion for the meeting:** the similarity metric (minimizer-Jaccard k15/w10)
is sound and containment is correctly rejected (length-disparity false merges).
The *separation-optimal* bar (~0.10) ≠ the *assembly-optimal* bar (~0.30) because
over-merging is costly. The operative default was a hardcoded **0.05** (over-
merging); it is now **0.30** (+2 RBMY, 0 regressions). This **corrects the earlier
`project_o5_already_built` framing**: forcing more cross-copy merges does NOT
recover starved copies — it collapses distinct ones; the real lever is coverage/
flow plus a similarity-magnitude guard on the completion gate (vg.rs:2803).
0.40 (+3 RBMY on this panel) is a candidate pending genome-wide validation.
