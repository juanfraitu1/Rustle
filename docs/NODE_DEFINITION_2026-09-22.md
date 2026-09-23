# Node definition — what limits it, what was fixed, and what is closed (§6x3–§6x7)

**2026-09-22.** Session goal: *improve node definition to avoid false positives and false negatives in all
modes.* This is the consolidated result. Every number is reproducible from the commands at the end.

## TL;DR

| | before | after | state |
|---|---|---|---|
| node **false negatives** (loci with homology that never become graph nodes) | 1,026 of 1,892 = 54.2% | `--min-cov-shorter 0.70` recovers 37 (0 lost) | ⭐ **improved, shipped behind a flag** |
| node **false positives** (over-merged loci) | 228 fused loci = 28/306 scored genes (9.2%) | — | ⛔ **CLOSED: perfect repair buys +0.021 F, of which the flag already has +0.016** |

**The two halves are independent (r1013).** They were treated as one problem for four operators.

## 1. The false-negative half — diagnosed, and the cause was not what §6w2 said

99.8% of the 1,026 dropped loci fail `cov_longer < 0.30` **alone**, and **679 align along a median 1.00 of
their own length** at passing identity and `alen` (r1002). Nothing is wrong with the victim, the alignment
or the boundary — `cov_longer` divides by the **longer** locus's span.

⛔**§6w2/r971's causal story is wrong (r1013).** It attributed 406 of those evictions (39.6%) to partners
holding ≥2 whole genes, i.e. to over-merge. But the **median fused locus is 125,450 bp and the constituent
gene its alignment lands in is 97,026 bp — a 1.02× denominator shrink available.** There is no giant to
shrink. The eviction is *a small locus against a genuinely large gene*, not *a small locus against a
readthrough fusion*.

That single mis-attribution produced **four independent refutations of a cause that was not there**:

| remedy | outcome |
|---|---|
| boundary pull-in (depth-adaptive k) | §6w6 ⛔ — shrinking costs the locus its own alignments |
| read-level splitting (§6w0 triggers) | r1001 ⛔ — needs a median **4.83×** shrink; a binary split gives 2×; 49.3% need >5× |
| node cut at the parent boundary | r846 ⛔ — doubles rather than separates; short pieces become hubs |
| per-pair local denominator | r1012 ⛔ — ceiling 13.0%, because only 1.02× is available |

⭐**The correctly-shaped fix normalises by the SHORTER locus**, which is the one lever that moved anything.

## 2. What shipped — `mcl_families --min-cov-shorter <C>`

A pair whose `cov_longer` fails also passes if the alignment covers ≥ C of the **shorter** gene's exonic
length; the edge weight then uses `cov_shorter`. Default `0.0` = OFF, **verified byte-identical when
unset**, 843 lib tests pass, unit test included.

- Held out chr2/chr8/chr10 (protein referee): **F up on 2 of 3, none down, precision up-or-equal in every
  cell**, largest cluster 0.96–1.13× baseline.
- Development: de novo F .214 → **.230**, guided .312 → **.322**, precision up in both.
- NPIP de novo: .727 → **.750** at **precision .923 → 1.000**.
- ⭐**C is not a fitted threshold**: C = 0.40/0.50/0.60/0.70/0.80 are identical to three decimals (r1014).
  **Use C = 0.70.**

⚠⚠**Never enable it without the exon conjunct** (`--min-exonic-bp 1 --min-shared-exon-frac 0.60`).
Unguarded this is register 913's refuted `min(la,lb)`: the largest component runs **5.7–16.3× baseline**
with the guard off versus **2.1×** with it on.

### Known regressions — both real, neither patchable

- ⛔**NPIP in GUIDED mode**: Soto F .833 → .800 at every C tested (0.70 and 0.90 alike).
- ⛔**semi-guided mode**: referee F .210 → .175 (C=0.90) / .148 (C=0.70), precision .973 → .833. The SD
  region node set was **built from self-alignment**, so containment there is generic — the escape fires
  35,728 times versus de novo's 53. A positional-overlap guard was built to fix this and **refuted**
  (r1010): it did nothing for semi-guided and cost de novo .230 → .199, below baseline.

⟹ **the flag is valid only where a node is a gene-like unit** (assembled locus or annotated gene body).

## 3. The false-positive half — CLOSED by an oracle ceiling

The over-merge population, characterised (r1015/r1016): **228 fused de novo loci**, costing **28 of 306
referee-scored genes (9.2%)**, concentrated on the thesis family — `NPIPB4`→`LOC112268174`,
`NPIPB5`→`SMG1P1`, `NPIPB12`→`SMG1P2`. Bimodal: 49.6% have a passenger <20% of the dominant gene (a trim),
34.2% hold two comparable genes (a split, refuted). The dominant gene is only 42.2% of the locus span.

**The ceiling (r1017).** All 233 fused loci replaced by one node per constituent gene — a perfect oracle
repair, same aligner and config as the baseline:

| arm | referee F | collapsed | NPIP Soto F | NPIP sens |
|---|---|---|---|---|
| shipped | 0.214 | 12 | 0.727 | 0.600 |
| `--min-cov-shorter 0.70` | **0.230** | 12 | **0.750** | 0.600 |
| ORACLE perfect split | 0.235 | **5** | **0.687** | **0.550** |

⛔**Perfect repair buys +0.021 and makes NPIP worse** (sensitivity 0.600 → 0.550).
⭐⭐**The flag captures 76% of that (+0.016) with no oracle, and beats it on the thesis family.**
**≤ +0.005 F remains in the FP half and reaching it costs NPIP. Do not build another trim/split/cut rule.**

⭐**Why**: the oracle *does* recover NPIPB4/B5/B12 (collapsed 12 → 5, cleanest of any arm), but the pieces
then fail the coverage gate the fused locus passed. **Over-merge is simultaneously the cause of NPIP's
collapse and the reason its members are admitted at all** (r1018).

## 4. Two open decisions — both the user's

1. **Flip `--min-cov-shorter` default to 0.70?** Clears its pre-registered bar, threshold-free,
   held-out positive — but regresses NPIP-guided, which is the thesis's own family. Same shape as
   `RUSTLE_JUNCTION_MAJORITY` (§6m8), which is also recorded as the user's call.
2. **Commit.** Nothing in §6x3–§6x7 is committed, per the standing rule that commits happen only on request.

## Reproduce

```sh
cd /mnt/linuxdisk/tmp/regress
# FN diagnosis and the split-arm re-scoring
python3 bench/node_graph_admission.py --paf dn16.paf --graph dn16.graph.tsv --gff chr16.genes.gff
# the flag
mcl_families --paf dn16.paf --gff dn16.gff3 --min-exonic-bp 1 --min-shared-exon-frac 0.60 \
  --min-cov-shorter 0.70 --out ARM
# scoring (protein referee truth built by bench/soto_vs_us_referee.py's protein_referee)
python3 bench/mode_family_score.py --clusters ARM.clusters.tsv --gff chr16.genes.gff \
  --soto REFEREE.tsv --chrom chr16
```

Pre-registrations: `PREREG_split_graph_admission_2026-09-22.md` ·
`PREREG_guarded_containment_admission_2026-09-22.md` · `PREREG_cov_shorter_adoption_2026-09-22.md` ·
`PREREG_overmerge_ceiling_2026-09-22.md`. Register rows **1000–1018**.

## 5. ⚠ Correction to §§1–4 above — the flag is an EDGE-WEIGHT change, not an admission fix (r1019)

Sections 1–2 present `--min-cov-shorter` as repairing the **admission** gate for the 679 fully-aligned
loci that `cov_longer` evicts. **That framing is wrong.** Setting the escaped edge's weight to the pair's
real mutual coverage (`identity × max(cov_longer, 0.30)`) instead of `identity × cov_shorter`:

| arm | OFF | escape, w = cov_shorter | escape, w = mutual |
|---|---|---|---|
| de novo chr16 | 0.214 | **0.230** | 0.214 |
| guided chr16 | 0.312 | **0.322** | 0.312 |
| chr2 / chr8 / chr10 | .236 / .436 / .197 | **.243 / .440** / .197 | .236 / .436 / .197 |
| **NPIP guided (Soto)** | **0.833** | 0.800 | **0.833** |

⭐⭐⭐**Every arm returns to baseline to three decimals.** The nodes still enter the graph; the partition
does not move. **Putting a locus in the graph is worth 0.000 — ranking its edge highly is the entire
result.** The NPIP-guided regression and the chromosome-wide gain are therefore *the same effect measured
on two populations*, and cannot be separated by tuning the weight: the one setting that removes the
regression removes the gain with it.

⚠ Semi-guided is harmed under both weights (0.159 vs 0.175), so the weight is not that regression's cause.
