# PREREG — two candidate improvements to the family definition

**Written 2026-09-21 before any scoring of either rule.** Two independent tests, one document, one
held-out arm. Either may be adopted without the other.

## Shared setup

- **Graph**: the SHIPPED DNA gene-body graph via `mcl_families --dump-graph` (edge set and exon conjunct
  are the shipped ones). Nodes are `chrom:start-end`.
- **Operator baseline**: MCL at I=2.8 through `bench/mcl_port.py`.
  ⚠**`mcl_port` is NOT bit-identical to the Rust MCL** (245 vs 168 clusters on chr2, register 917), so
  the comparator is mcl_port-MCL **on the same graph**, never the shipped Rust F of 0.7016/0.7123.
- **Development**: chr2, chr8, chr10, chr16 — all already exposed this session.
- **HELD-OUT: chr4 + chr9**, run LAST. chr4 and chr9 have 3 and 16 doc mentions respectively (vs chr1
  194, chr15 89) and no ledger result. Together: 9 cover families, 52 members, 8 multi-family genes.
  ⚠chr19 was considered and REJECTED as held-out: it has **1** Soto family with >= 3 members.

## Truth — and a correction to how it has been read

Soto's S1C carries an explicit **`No. Assigned Families`** column, and 149 of 2,334 CHM13 gene IDs
(6.4%) are assigned to more than one family. The derived family-set size equals that column for
**2,334 / 2,334** records, so **the published truth is a COVER, authored as one.**

Every scorer in this project applies `if len(fids) != 1: continue` — *"a partition needs one label per
gene"*. That drops the multi-label genes, and because dropping a member can push a family below the
>= 3-member floor it also removes whole families: on chr16 the truth goes **15 families / 70 members
(cover) -> 8 / 43 (partition)**, on chr7 12 -> 9. So the restriction suppresses ~47% of chr16's families.

- **TRUTH USED HERE = the cover**: a gene's truth label set is the union of its Family IDs.
- Matching stays **one-to-one bipartite** on overlap (as shipped), so a truth family may only claim one
  cluster; a cover PREDICTION gains only when two truth families that share a gene can each match a
  cluster containing it.
- ⚠Keying is by **Gene Name** because no CHM13 CAT annotation (`CHM13_G*` IDs) is available to reach
  loci; 249 names map to >1 gene ID and 16 names are multi-family PURELY by name collision. Every
  headline is reported a second time with those 249 ambiguous names dropped; if the two disagree in
  sign, neither rule is adopted.

## Test 1 — the prediction may emit a COVER

O1 currently emits a strict partition: 0 of 2,670 loci sit in >1 cluster, so a fusion gene is
inexpressible. Register 845 refuted dual membership on the **truth** side and closed with an explicit
re-open condition — *"do not re-propose without a family whose chimera fraction is large"* — now met
(chr10 AGAP, 4 of 10 members curated readthroughs = 40%). **The prediction side has never been tested.**

    RULE (purely combinatorial, no weight threshold):
      run MCL as shipped, then for every node v and every cluster C with v not in C:
          add v to C  iff  v has >= k edges into C
      k swept over {2, 3, 4}.  k = 2 is triangle support, the project's existing §6kd notion.

Register 846 (cutting a chimera at its parent boundary DOUBLES rather than separates, and short pieces
become hubs) is about SPLITTING a node and does not bind: nothing here changes any node.

## Test 2 — neighbourhood Jaccard as an MCL edge weight, size-gated

`J_N(u,v) = |N(u) ∩ N(v)| / |N(u) ∪ N(v)|` is the strongest separator measured in this project
(AUC 0.924 at component size >= 10, 0.872 at 3-4) and the only one that SURVIVED size-residualisation
(0.739 -> 0.740, corr +0.49, where r523's betweenness collapsed 0.683 -> 0.531).

It has only ever been judged as a **post-clustering merge** (register 933, refuted). Register 916/917
established that this is the wrong seat for a scalar: guarded containment was +0.013 in connected
components and +0.025 fed into MCL — *"a scalar can't be judged at the operator's weakest setting."*

    RULE:  w'(u,v) = w(u,v) * (1 + J_N(u,v))   for edges inside a component of size >= 5
           w'(u,v) = w(u,v)                    otherwise
    then MCL at I = 2.8 as shipped.

The size gate is part of the rule, not a knob: J_N is AUC **0.500 — exact chance** at component size 2,
and 60% of truth families are pairs (register 931). Sweep: gate in {3, 5, 8}.

## Bars — stated before looking

Pooled bipartite F over the development chromosomes, against the cover truth, versus mcl_port-MCL on the
same graph:

1. **PRIMARY:** pooled F must rise by **> +0.02** — strictly outside the project's stated tied band
   (register 922, where label-prop's +0.0172 was declared not a switch).
2. **GUARD A (precision):** pooled precision must not fall by more than **0.03**. Test 1 can only add
   members to clusters, so precision is its exposed flank.
3. **GUARD B (small families):** the number of matched 2-member truth families must not fall at all.
   Register 921: every triangle operator dissolved ALL 362 two-member groups, and 2 is the modal size.
4. **GUARD C (coverage):** node coverage must not fall by more than 2 percentage points.
5. **HELD-OUT (chr4+chr9), run LAST at the single (k or gate) chosen on development:** pooled F must
   rise by any positive amount. **A held-out failure means NOT ADOPTED**, whatever development showed
   ([[feedback_hold_a_substrate_back]]).

Parameters are frozen on development before chr4/chr9 is touched. Reporting is sensitivity, precision
and bipartite F together, never F alone.

## What a negative means

Test 1 negative ⟹ dual membership is not worth its precision cost even where the truth is a cover, and
the partition is the right shape for the definition — which is itself an answer to the advisor.
Test 2 negative ⟹ J_N is a scorer, not a definitional ingredient, closing the last live lead from §6u3.
