# Pre-registration — one edge rule instead of three tiers, and the 60–90% band

**Written 2026-09-24 (§6zi), before the rule is run on any catalog.** Advisor's two objections, via the user: the
family definition has too many modes to follow, and the 60–90% identity band — old enough to have diverged, young
enough to still cross-map partly — may be the most interesting biology, which the spectrum (§6zh, row 1096) shows
we recover at 0.24–0.36.

## What the spectrum established (fixed inputs to this design)

- The sensitive pass (`-k11 -w5`) subsumes asm20: over 116 Compara pairs at ≥ 50% protein identity on chr16, asm20
  finds exactly one pair the sensitive pass does not (NPIPB12–NPIPB3). Two nucleotide tiers are one.
- The 80–90% band is lost to the coverage clause (7 of 8 misses have a record at 0.86–1.00 identity covering 8–44%
  of the shorter representative), the 60–70% band half to coverage (9/22) and half to seeds (13/22).
- Post hoc, union-of-records coverage at a 0.30 floor gave 80–90% recall 0.75 and 70–80% 0.50 at pair precision
  0.87 (shipped clause: 0.97). That number was looked at; this document is the test it now has to pass.

## The rule (one sentence)

**Two expressed loci are joined when their representative transcripts share ≥ 30% of the shorter one in aligned
sequence at ≥ 60% nucleotide identity**, where "aligned sequence" is the union of every local alignment between
them (one `minimap2 -x asm20 -k11 -w5 -c -X --no-long-join` pass), and the identity is the best record's. The
protein extension (`--protein-tail`) is not part of the rule; it is the separately named step for < 60% nucleotide
identity, where membership is a gene-tree question.

Implementation: the union coverage already exists in the builder as `RUSTLE_ER_SUM_COVERAGE=1` (blocks ≥ 200 bp,
`RUSTLE_ER_SUM_MIN_BLOCK`); the floor is `RUSTLE_ER_MIN_COVERAGE=0.30` (new env, default 0.50). Unset, every output
is byte-identical. The asm20 pass is left running (it adds one pair); the one-sentence rule is what the definition
SAYS, the second pass is an implementation detail it does not depend on.

⚠ **Prior result on this switch (register, "Union-of-records coverage is the biggest win of the session", class H):**
on the internal within/cross-family audit, single-record coverage gave 319 within-family and 0 cross-family pairs
(precision 1.00) and the union 423 within and 45 cross (0.90) — "it buys recall by admitting exactly what the
per-record rule exists to prevent". That verdict was against the catalog's own families; this test re-asks the
question against two EXTERNAL truths and for the specific band the advisor cares about, and its bar (precision
≥ 0.85 vs Compara, held-out drop < 0.05) is set knowing that 0.90 is the expected level.

## Test — committed now

Catalogs built twice, shipped rule vs single rule, on the family-rich development chromosome (human chr16, A119b)
and on the held-out gorilla contig NC_073244.2 (OR6737 testis), scored at the PAIR level against the two
independent truths the project has for them:
- chr16: Ensembl Compara same-chromosome paralogue pairs with both genes expressed (§6zh truth) — recall by band,
  precision of catalog pairs Compara can judge;
- NC_073244.2: the protein-referee families (`sec/ref/NC_073244.2.tsv`, 113 families / 775 genes) — pair recall
  and precision, plus the family count and the largest family (the over-merge sentinel).
Loci map to genes by RefSeq exon overlap; a catalog pair is two copies of one family whose loci map to two genes.

| outcome | verdict |
|---|---|
| chr16 recall in the 60–90% bands rises by ≥ 0.15 (absolute, pooled over the three bands) AND precision vs Compara stays ≥ 0.85 AND NC_073244.2 pair precision vs the referee drops by < 0.05 with the largest family growing by < 50% | ⭐ **adopt the single rule as the default** |
| recall rises but precision falls below 0.85 on chr16 or by ≥ 0.05 on the held-out | ⚠ **keep it as the named option; the band costs precision** |
| recall does not rise ≥ 0.15 | ⛔ **the rep-level gain does not survive the catalog; report and leave the tiers** |

**Predicted:** ⭐ on chr16 (the rep-level numbers are the catalog's own edges); the held-out is the risk — the gorilla
catalog's families are near-identical SD copies where a 0.30 floor could let repeat-mediated partial alignments
glue neighbours (the reason the repeat gate exists), so the largest-family sentinel is the number to watch.

I will not change the floors, the truths, or the bar after seeing any number.

---

# OUTCOME (2026-09-24) — catalogs in `/mnt/linuxdisk/tmp/gw22/rule/`, scored with `bench/identity_spectrum.py --catalog`

Two corrections to the design text, found while running it, neither touching the bar:
- The catalog builder already runs the sensitive pass ALONE (`RUSTLE_ER_SENSITIVE_ONLY` defaults on; both logs say
  "asm20 run SKIPPED"), so the shipped catalog is already one nucleotide tier at floor 0.60. The single rule changes only
  the coverage clause (best record ≥ 0.50 → union of records ≥ 0.30).
- The catalog scorer's RefSeq branch mapped every copy to one chromosome-wide span built from exon lines with no gene
  attribute (0 families scored); fixed to `gene_name` → `gene_id`, skip empty. Recall is reported over the §6zh universe
  (275 Compara pairs with both genes expressed, `chr16.truth_pairs.tsv`), the denominator the catalog cannot move.

## Results

| substrate | rule | copies | families | largest (genes) | pair recall | pair precision |
|---|---|---|---|---|---|---|
| gorilla NC_073244.2, referee truth | shipped | 357 | 52 | 37 | 0.117 | 0.983 |
| | single rule | 527 | 70 | **81 (+119%)** | 0.205 | 0.970 |
| human chr16, Compara truth | shipped | 1,418 | 258 | 25 | see bands | 0.268 all copies · **0.864 multi-exon** (70/81) |
| | single rule | 2,576 | 430 | 32 (+28%) | see bands | 0.119 all · 0.791 multi-exon (72/91) |

chr16 recall of the 275-pair universe, by Compara protein identity (shipped → single rule):
≥ 90 **0.903 → 0.806** (28 → 25/31) · 80–90 0.917 → 0.917 (11/12) · 70–80 0.786 → 0.857 (11 → 12/14) · 60–70 0.552 → 0.517
(16 → 15/29) · 50–60 0.767 → 0.867 · 30–50 0.197 → 0.246 · < 30 0.020 → 0.020.
**Pooled 60–90: 38/55 → 38/55, Δ = 0.000.** The rule gains CHST5–CHST6, MT2A–MT3, MT1X–MT3 and loses NPIPB12/B13–NPIPB3,
NPIPB12–NPIPB5, NPIPA5–NPIPB12/B13: with more edges the γ-quasi-clique blocking re-partitions the NPIP-B family rather
than extending it.

## Verdict — ⛔ by the pre-registered table (recall does not rise ≥ 0.15)

The rep-level gain (§6zh post hoc: 80–90% 0.25 → 0.75) does not survive the catalog, for a reason the design missed:
**the family definition is transitive and the shipped catalog already recovers the band.** The rep-level 0.333 at 80–90%
counted DIRECT edges; in the shipped families 11 of the 12 pairs are already together through a third copy (the NPIP-A/B
cores), and the one miss (CHST5–CHST6) has neither gene in any family. 70–80% is at 0.786, 60–70% at 0.552 — the losses
that remain are the seeding losses (MT/NPIPA×NPIPB cross pairs at 65–70% nt) and node absences, not the coverage clause.
On the held-out gorilla contig the rule buys +0.088 pair recall by doubling the largest family (37 → 81 genes), the
over-merge sentinel the design named as the number to watch. **The tiers stay; the rule is not adopted.** The one-sentence
description of the shipped definition is still true, with the shipped clause: *two expressed loci are joined when their
representatives share ≥ 50% of the shorter one in one colinear alignment at ≥ 60% nucleotide identity, and families are the
transitive closure (γ-quasi-clique blocks) of those joins.*

## What the chr16 catalog's precision actually is (post hoc, labelled)

Precision over ALL judgeable within-family pairs is 0.268 under the shipped rule, far below the 0.96 of the rep-level
edges. The gap is not the edges: **85% of the chr16 catalog's copies are single-exon loci** (gorilla: 35%), and the families
Compara rejects are 95% single-exon, 76% overlapping no annotated exon, at median identity 0.83 — 73 repeat-glued groups of
unspliced pileups inside big genes (CDH8, FTO, RBFOX1, HYDIN, WWOX, GRIN2A…), r1043's class. Restricted to multi-exon copies
the catalog is 0.864 precise (0.791 under the rule). This is a node-admission fact about the human chr16 substrate, already
registered as a dead end for rule-making (r356, r1049: dropping single-exon loci costs precision and recall on gorilla);
it is reported here so the 0.268 is not read as an edge-rule failure.
