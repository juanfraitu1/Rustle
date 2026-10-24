# Post-MCL strongest-edge reassignment of family members (2026-10-15)

Pre-registration, written before the rule was run on any catalog. Code: `bench/o1_reassign/` at `cb237430`
(16 tests). Branch `post-mcl-reassign`.

## Question

Human A119b fails claim F1 of figure 8 (131 / 148 Liftoff copy pairs share a family; the bar is 0.90). Exploring the
17 split pairs (2026-10-15) showed that 7 of them come from one MCL decision in a tandem array. The user asked to try a
rule that undoes such decisions: after MCL, a member moves to the family that holds its strongest edge when that edge
is stronger than every edge it has inside its own family. Does the rule fix F1 without harming the held-out catalogs,
and do its moves agree with independent truth?

## Known before running (recorded now)

- **F1 as shipped** (genome, sequence_ID ≥ 0.95): A119b 131 / 148, testis 34 / 34, gorilla KB3781 11 / 12.
- **The decoy test as shipped** (register 1284, genome, ≥ 0.95): A119b 75 / 91, testis 12 / 12, KB3781 4 / 4 against
  null means 0.121, 0.102 and 0.012 (p 0.000999 each).
- **A119b's split pairs, traced through the family graph** (`bench/` not involved; session scratch):
  - 7 CYKILR pairs: the chr19 38.8-40.1 Mb array is one dense block (edge weights 0.95-0.999) that MCL cut into
    MCL42, MCL451 and MCL566. The CYKILR source node's strongest edge (0.986) goes to MCL42, which holds all 7 copies;
    its best edge inside MCL451 is 0.963. Two other MCL451 members also have their strongest edge into MCL42 (0.994 vs
    0.976; 0.985 vs 0.976), and the fourth ties (0.916 and 0.916).
  - 8 pairs (5 genes) have no edge: the 0.60 shared-exon floor rejects them. The rule cannot join them.
  - RPL23AP88's copy node has edges of 1.000 inside its family and 0.999 to the source's family, so the rule leaves
    it. LOC124906865's cross edge is 0.391, below both nodes' best inside edges.
  - Expected from this: the CYKILR source node moves to MCL42 and A119b's F1 rises to about 138 / 148. That is why A119b
    is development: its F1 is not evidence for the rule.
- 14 of the 17 split pairs already share a superfamily of the emitted hierarchy, but 11 of them through SF166, which
  holds 2,991 members in 435 families.
- **Not seen**: how many members move in any catalog, any F1 or decoy number after the rule, any Soto comparison after
  the rule.

## Method (fixed)

**The rule** (`reassign_lib.reassign`; no parameter). For a member v of family F:
- in(v) is its strongest edge to another member of F;
- out(v, G) is its strongest edge into family G ≠ F.

If max over G of out(v, G) > in(v), v moves to that G; when two families tie, the smaller family name wins. Members are
visited in sorted order and see the moves already made (asynchronous passes), and passes repeat until one moves
nothing, up to 50.

Only members move:
- an edge to a node in no family is ignored, and no node joins or leaves the catalog;
- a member left alone in its family is reported (lost).

The edges are the shipped graph: `mcl_families --dump-graph` from each catalog's cached all-vs-all
(`bench/o1_hier/recluster.sh`). Weights are written to 6 decimals, so edges closer than 1e-6 count as ties.

**Gates** (a sample whose gate fails has no verdict):
- the re-clustered `clusters.tsv` equals the shipped one byte for byte;
- every row of the shipped copy table maps (gene → span node → folded representative) to a node of its own shipped
  family.

**Scoring** (`run_reassign.py score`):
- F1 before and after, every scope and identity floor, with figure 8's code (`run_pairs.f1`);
- the Liftoff pairs (genome, ≥ 0.95) the moves join and split;
- the decoy test's primary row (genome, ≥ 0.95, 1,000 decoy sets, seed 7) before and after;
- human only: Soto gains and losses. Per moved member reached by a gene of Soto's 491 families (soto_subset's de
  novo universe), count its Soto family-mates in its old family and in its new one; more in the new one is a gain.
  Also reported: pairwise precision and recall against Soto over the universe's members, before and after.

## Bars

| id | sample | bar |
|---|---|---|
| D | A119b (development) | F1 ≥ 0.90 after the rule |
| H1 | testis, KB3781 (held out) | F1 stays ≥ 0.90, and the decoy test still passes (p < 0.001) |
| H2 | testis, KB3781 | Liftoff pairs joined − split ≥ 0 in each; in testis, Soto gains > losses when gains + losses ≥ 5 (fewer: not powered) |
| guards | every sample | no member lost; F1's covered pairs change by ≤ 2%; the largest family grows by ≤ 10%; the passes converge |

**Evidence of benefit** (held out): Liftoff joined > split in testis or KB3781, or a powered Soto H2 in testis.

## Readings (fixed now)

1. **Adopted as an opt-in flag** (then ported to `mcl_families`): D, H1, H2 and the guards hold, and there is evidence
   of benefit.
2. **Harmless but unproven**: D, H1, H2 and the guards hold, but nothing held out changes enough to count as evidence.
   It is recorded as a rule that fixes A119b's F1 without held-out support, and it is not adopted.
3. **Not adopted**: D, H1, H2 or a guard fails.

A119b's moves outside the chr19 array (Soto gains and losses, Liftoff joined and split) are reported as a description,
not as a bar, since the array motivated the rule.

## Results

Added after the runs. Nothing above this heading changes once a statistic has been seen.

## Results (2026-10-15)

Run: code `cb237430`, outputs `/mnt/linuxdisk/tmp/o1_reassign/<sample>/` (`apply.json`, `score.json`, `moves.tsv`).
Both gates held in all three samples: the re-clustered `clusters.tsv` equals the shipped one byte for byte, and every
copy row maps to a node of its shipped family. Re-clustering took 6 s (testis), 17 s (KB3781) and 155 s at 4.7 GB
(A119b).

| | A119b (development) | testis (held out) | KB3781 (held out) |
|---|---|---|---|
| members moved (moves, passes) | 1,609 of 9,907 (2,006, 6) | 5 of 1,115 (5, 3) | 11 of 918 (11, 3) |
| families | 2,150 → 1,944 | 342 → 340 | 326 → 322 |
| largest family | 260 → 349 | 47 → 47 | 34 → 34 |
| members lost (left alone) | 9 | 1 | 0 |
| **F1** (genome, ≥ 0.95) | 131 / 148 → **137 / 148 = 0.926** | 34 / 34 → 34 / 34 | 11 / 12 → 11 / 12 |
| F1 without the development contigs | 113 / 128 → 119 / 128 | 33 / 33 → 33 / 33 | 11 / 12 → 11 / 12 |
| Liftoff pairs joined / split | 7 / 1 | 0 / 0 | 0 / 0 |
| decoy test (k / n, null mean; p) | 75 / 514, 11.0 → 80 / 514, 8.4; 0.000999 | 12 / 38, 1.2 → unchanged | 4 / 78, 0.05 → unchanged |
| Soto gains / losses (same, unlabelled) | 50 / 14 (21, 1,524) | 1 / 0 (0, 4) | — |
| Soto pairwise precision | 0.7344 → **0.6881** | 0.8717 → 0.8721 | — |
| Soto pairwise recall | 0.3008 → 0.3150 | 0.8386 → 0.8418 | — |

Bars and guards:
- D holds (A119b F1 0.926 ≥ 0.90).
- H1 holds: held-out F1 and the decoy test are unchanged.
- H2 holds trivially but without evidence: no Liftoff pair is joined or split in either held-out sample, and testis's
  Soto comparison is not powered (1 gain, 0 losses; the bar needs 5).
- **Guards fail.** A119b's largest family grows by 34% (260 → 349; bar 10%), and members are lost in A119b (9) and
  testis (1). KB3781's guards hold.

### Verdict (2026-10-15)

**Reading 3: not adopted.** The guards fail in two of three samples. Even with the guards set aside, nothing held out
shows a benefit, which is reading 2 at best. The rule does fix A119b's F1 as predicted, but in the large catalog it
merges much more than the tandem array: Soto pairwise precision falls by 0.046 there.

### Post hoc (descriptive, A119b)

- **The array is a small part of it.** Only 9 of the 2,006 moves are inside the chr19 38.8-40.1 Mb array. The other
  moves spread over the genome; chr13 has the most, 283 moved members.
- **Members MCL placed without a direct edge move first.** 718 moves (36%) start from in = 0: MCL put the member in its
  family through flow, with no edge to any member of it, so any edge to another family pulls it out.
- **Large families attract.** The largest family (MCL3) keeps 234 of its 260 members and gains 115. The top five
  destination families take 517 of the moves (MCL3 124, MCL121 121, MCL60 106, MCL617 94, MCL12 72).
- **The deciding gap is often small.** The median out − in is 0.171, but 381 moves (19%) are decided by less than 0.01
  and 685 (34%) by less than 0.05.
- **Outside the array** (1,601 of the 1,609 moved members; pre-registered as a description):
  - every Soto-labelled move is outside it: 50 gains, 14 losses (21 same, 1,516 unlabelled);
  - Liftoff pairs: 0 joined and 1 split, so all 7 pairs the rule joins are the array's CYKILR pairs.
  - Pairwise precision against Soto still falls, because the absorbed members bring their neighbourhoods.
