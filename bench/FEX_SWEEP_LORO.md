# Raising `--min-shared-exon-frac` from 0.30 to 0.60 — leave-one-region-out validated

**Toward the goal "an RNA-level multi-copy family definition with high precision, recall and bipartite
matching": the single highest-value change found so far is raising the shipped shared-exon floor from
0.30 to 0.60.** It improves bipartite F in every region and is validated leave-one-region-out.

Scoring: FIXED, node-independent truth set (Soto families with >= 2 genes in the node set), unclustered
genes get singleton predictions so recall is penalised. This is the §6o2 scorer, not §6o1's (which was
retracted for conditioning the truth on the prediction). All arms rebuilt with the CURRENT binary and
`--min-exonic-bp 1`, so the only thing that varies is `--min-shared-exon-frac`.

## The sweep

| region | f_ex | R | P | **F** | pairwise sens | pairwise prec |
|---|---|---|---|---|---|---|
| chr5/7/21 (82 fams) | 0.00 | 0.772 | 0.845 | 0.807 | 0.848 | 0.822 |
| | **0.30 (shipped)** | 0.758 | 0.862 | **0.807** | 0.834 | 0.851 |
| | 0.40 | 0.768 | 0.881 | 0.821 | 0.837 | 0.854 |
| | **0.60** | 0.783 | 0.918 | **0.845** | 0.830 | **0.898** |
| | 0.70 | 0.780 | 0.925 | 0.846 | 0.829 | 0.894 |
| chr15/17 (131 fams) | 0.00 | 0.645 | 0.749 | 0.693 | 0.599 | 0.367 |
| | **0.30 (shipped)** | 0.651 | 0.814 | **0.723** | 0.555 | 0.508 |
| | 0.40 | 0.655 | 0.841 | **0.737** | 0.538 | 0.544 |
| | 0.60 | 0.643 | 0.851 | 0.733 | 0.524 | 0.568 |
| chr1 (55 fams) | 0.00 | 0.649 | 0.774 | 0.706 | 0.572 | 0.475 |
| | **0.30 (shipped)** | 0.694 | 0.843 | **0.761** | 0.520 | 0.637 |
| | 0.60 | 0.702 | 0.874 | **0.779** | 0.505 | 0.716 |
| | 0.70 | 0.691 | 0.877 | 0.773 | 0.477 | 0.728 |

| f_ex | mean F | min F | mean pairwise precision |
|---|---|---|---|
| 0.00 | 0.735 | 0.693 | 0.555 |
| **0.30 (shipped)** | 0.764 | 0.723 | 0.665 |
| 0.40 | 0.772 | 0.737 | 0.685 |
| 0.50 | 0.774 | 0.734 | 0.707 |
| **0.60** | **0.786** | 0.733 | 0.727 |
| 0.70 | 0.784 | 0.733 | **0.760** |

## Leave-one-region-out (the value is chosen on two regions, reported on the third)

| held-out region | picked on the other two | F @ 0.30 (shipped) | F @ picked | delta |
|---|---|---|---|---|
| chr5/7/21 | 0.60 | 0.807 | **0.845** | **+0.039** |
| chr15/17 | 0.60 | 0.723 | **0.733** | **+0.009** |
| chr1 | 0.70 | 0.761 | **0.773** | **+0.012** |

**Never worse in any fold**, and 0.60 is the pick in 2 of 3. This is not a dev-set optimum: the value was
never chosen on the region it is reported on.

## What it costs

Recall and pairwise sensitivity fall slightly (chr15/17 sens 0.555 → 0.524; chr1 0.520 → 0.505). The gain
is concentrated in precision, which is where §6o2 located the remaining gap. chr5/7/21 gains on BOTH axes
(R 0.758 → 0.783, P 0.862 → 0.918).

## Caveats that must travel with this

- **Three regions is a small LORO.** Real, but weak; a fourth region would strengthen it materially.
- **§6km's ground-truth ceiling.** GENCODE-vs-RefSeq on the same construction scores F 0.79-0.81, so
  chr5/7/21's **0.845 is at or above the ceiling for agreement with a single-annotation truth**. Treat that
  region's number as saturated rather than as headroom.
- Soto is one annotation's view and is itself a cover (§6n7), so "precision against Soto" is not the same
  as correctness.

**Recommendation: propose `--min-shared-exon-frac 0.60` as the new default**, reporting the recall cost
alongside. The decision is the user's; nothing has been changed.
