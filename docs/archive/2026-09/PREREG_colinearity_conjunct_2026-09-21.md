# Pre-registration — does block colinearity separate the 76 rejected containment≥0.90 pairs?

**Written 2026-09-21, §6w1, before any colinearity score is computed.** User: *"should we lower the
coverage metric so more NPIP members are found? I know this increases false positives, but we can
enforce colinearity checks."*

## Why not just lower the coverage/identity floor

Already closed. r919 swept containment/identity thresholds from 0.30 to 0.99 on exactly this population
and found the precision curve **flat at ~0.25 the whole way** (true 0.987 vs false 0.989 identity — no
separating value). Lowering a scalar floor is not proposed here; this tests a **structural** signal
instead, per §6t7's own conclusion ("any fix must come from STRUCTURE, not a better scalar") and to give
Canzar a non-arbitrary-threshold story if it works.

## What "colinearity" means here, stated precisely

`vg_multiplicity_containment.py` and every earlier pairwise-metric script (`pair_metric_sweep.py`,
`chained_pair_score.py`) keep only the **single best** PAF record per gene pair. Checked today on
chr2.paf: **95.6% of all 90,765 gene pairs on chr2 have more than one PAF record** (minimap2 was run
`-P`, retaining every chain). Discarding all but the best record throws away exactly the information a
colinearity check needs — this is a genuinely untested signal, not a rerun of §6t4 (DP chaining was
refuted as *redundant with the best single record*; this uses the *other* records, which §6t4 never
touched).

For a gene pair (A, B):

1. Collect every PAF record with `{qname,tname} == {A,B}`, **excluding self-pairs and any block with
   `nmatch < 30bp`** (a low bar to drop single-minimizer noise, fixed now, not tuned on the labels below).
2. Normalize each record into a block `(a_start, a_end, b_start, b_end, strand)` in A's and B's own
   body-relative FASTA coordinates (`qs/qe` or `ts/ts` depending on which side A is), `strand` = the
   record's own relative strand.
3. If the surviving blocks disagree in `strand`, the pair is **NOT colinear** (score 0) — a real
   duplicate should sit at one consistent relative orientation across every block.
4. Otherwise sort blocks by `a_start` and take the sequence of `b_start` values. **Colinearity score**
   = fraction of concordant pairs among all `C(n,2)` block-pairs (Kendall-tau concordance): for
   `strand +`, a comparison is concordant if `a_start` and `b_start` order the same way; for `strand -`,
   concordant if they order oppositely. Score is undefined (reported separately, not imputed) when a
   pair has **fewer than 2 surviving blocks** — colinearity cannot be evaluated and this is itself a
   result to report, not a value to guess.

This is a genuine synteny/concordance measure, distinct from everything already refuted: not identity/
containment (r919), not the exon conjunct (already applied upstream), not exact intron-junction identity
(§6u1, which failed specifically because ~50% of members are intronless — colinearity needs ≥2 blocks
of *any* kind, not spliced junctions), not node multiplicity (§6t8).

## Predicted failure mode, stated before looking

If most of the 76 pairs have 0-1 surviving blocks, colinearity is **undefined or trivially satisfied**
for most of the population and this signal is underpowered — I will report the block-count distribution
for TRUE vs FALSE regardless of whether it helps, exactly as §6t8 reported the multiplicity distribution.

## Population — frozen, identical to `docs/PREREG_vg_containment_2026-09-20.md`

The same **76 pairs** the shipped rule rejects with containment ≥ 0.90, both endpoints Soto-labelled,
on held-out chr2/chr8/chr10 (19 TRUE / 57 FALSE) — same truth (`bench/soto/soto_famCN_S1C.tsv`), same
shipped-edge exclusion (`famgraph/chrN.graph.tsv`), same PAFs (`/mnt/linuxdisk/tmp/heldout/chrN.paf`).
Reusing the frozen population makes this directly comparable to r919 (baseline 0.250) and r906
(multiplicity, AUC 0.681, best precision 0.324).

NPIP itself is on chr16, a **dev** chromosome (`chr16 NPIP cluster` in memory) — this test runs
held-out first per the "hold a substrate back" rule. Only if colinearity clears the bar below does it
get checked against the actual NPIP case on chr16, and chr16 numbers do not set the adoption threshold.

## The bar — committed now, same shape as the VG prereg for comparability

| outcome | verdict |
|---|---|
| a `colinearity ≥ t` cut reaches **precision ≥ 0.60** while keeping **≥ 10 of the 19** TRUE pairs | ⭐ **HELPS** — worth building as a real conjunct |
| precision 0.40–0.60 at ≥ 10 TRUE | ⚠ **PARTIAL** — real signal, not enough alone |
| best precision < 0.40, or < 10 TRUE retained, or >50% of the population has <2 blocks (underpowered) | ⛔ **NO** |

Baseline to beat: **0.250** (containment alone, r919). Comparable prior art: multiplicity best precision
**0.324** at AUC 0.681 (r906/§6t8).

I will report the full precision/recall curve over `t`, the AUC, and the block-count distribution
(median blocks, % with <2) for TRUE vs FALSE. I will not change the population, the truth, the 0.90
containment floor, or the nmatch≥30bp block filter after seeing the numbers.

## OUTCOME (2026-09-21, `bench/colinearity_conjunct.py`)

Population reproduced exactly: **76 pairs, 19 TRUE, 57 FALSE** — matches r919/the VG prereg bit for bit,
confirming the reimplementation is sound.

⛔ **NO.** Best precision at ≥10 TRUE retained = **0.247** (below baseline 0.250 — worse than doing
nothing), AUC = **0.490** (chance). The threshold sweep is flat at precision 0.247 across every `t` from
0.50 to 1.00.

The reason is not the predicted failure mode (too few blocks). **Both TRUE and FALSE populations have a
block-count median of exactly 2.0, with 0% under 2** — the population is fully scoreable, not
underpowered by block count. The actual reason: with exactly 2 blocks, there is only **one** pairwise
comparison, and minimap2's own secondary-chain placement almost always reports it in spatial order
regardless of whether the relationship is a real duplicate or a repeat/domain match — **both classes hit
colinearity score 1.000 at the median.** The extra PAF records are not independent structural evidence;
they are mostly small satellite fragments of the same local alignment, and fragments of one contiguous
region are trivially "in order" whether or not the two genes are actually related. This is a variant of
§6t7's own diagnosis — *"the pairwise view lacks the information"* — extended to the multi-block view:
having a second data point does not add information if that second point is not independent.

**Closes this specific structural avenue.** Combined with r919 (scalar thresholds: NO), §6u1/§6u2
(exact intron-junction identity: NO, biological ceiling), and §6t8 (node multiplicity: best precision
0.324, PARTIAL at best) — **four independent structural/scalar signals have now been tried against this
exact 76-pair population and none clears even a 0.40 precision bar.** The shipped rejection of these 21
asymmetric true pairs remains the correct call: no known signal recovers them without an unacceptable
false-positive cost. Lowering the coverage/identity floor is not recommended, with or without a
colinearity guard — the guard has no separating power to contribute.

> **Generator (2026-09-22 consolidation):** `python3 bench/locus_probes.py colinearity ...` — the original `bench/colinearity_conjunct.py` was folded in verbatim and verified identical on its documented inputs (§6z3); the register rows above cite this file.
