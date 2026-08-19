# O1 fresh-emission validation

The current Rustle binary rebuilt nodes from regional BAMs extracted from the original
full-genome alignments. Frozen node ids, family ids, and dispositions were not inputs to Rustle;
they were used only after emission for coordinate matching.

Intervals were fixed from the 19-family audit with ±10,000 bp padding. Fresh loci match
**124/133** frozen nodes at ≥0.50 overlap of the shorter interval, including
**72/75** independently named target nodes; **69/75** land in
the modal fresh family for their test case. Of **16** previously
flagged conflicting-gene nodes, **14** are independently re-emitted and
**1** re-enter the same fresh component as the named target family.
Separately, **1/1** documented broad-family outgroups are re-emitted,
and **1/1** remain in the broad RNA family.

| species | family | frozen→fresh matched | named emitted/in target | related emitted/in target | conflicts emitted/in target | fresh named-family ids |
|---|---|---:|---:|---:|---:|---|
| GGO | C4 | 2→2 (2 unique) | 2/2 | 0/0 | 0/0 | GWFAM2 |
| GGO | TSPYL | 2→2 (2 unique) | 2/2 | 0/0 | 0/0 | GWFAM3 |
| GGO | GSTM | 4→4 (4 unique) | 4/4 | 0/0 | 0/0 | GWFAM0 |
| GGO | RGPD | 6→6 (6 unique) | 1/1 | 0/0 | 0/0 | GWFAM4 |
| GGO | ANKRD18 | 4→4 (3 unique) | 1/1 | 0/0 | 0/0 | GWFAM5 |
| GGO | HERC2 | 2→2 (2 unique) | 1/1 | 0/0 | 0/0 | GWFAM8 |
| GGO | MAGEA | 11→11 (11 unique) | 5/5 | 0/0 | 1/1 | GWFAM9 |
| GGO | DAZ | 2→2 (2 unique) | 2/2 | 0/0 | 0/0 | GWFAM12 |
| GGO | APOBEC3 | 2→2 (2 unique) | 2/2 | 0/0 | 0/0 | GWFAM13 |
| HSA | GOLGA6_8 | 22→19 (19 unique) | 10/9 | 1/1 | 0/0 | GWFAM16,GWFAM17 |
| HSA | TBC1D3 | 6→6 (6 unique) | 5/4 | 0/0 | 0/0 | GWFAM21,GWFAM22 |
| HSA | RANBP2_RGPD | 9→8 (8 unique) | 8/8 | 0/0 | 0/0 | GWFAM26 |
| HSA | RABL2 | 2→2 (2 unique) | 2/2 | 0/0 | 0/0 | GWFAM28 |
| HSA | AMY1 | 6→6 (6 unique) | 2/2 | 0/0 | 0/0 | GWFAM4 |
| HSA | PCDHB | 2→2 (2 unique) | 2/2 | 0/0 | 0/0 | GWFAM31 |
| HSA | MAGEA | 10→10 (10 unique) | 7/7 | 0/0 | 3/0 | GWFAM34 |
| HSA | RBMY1 | 3→2 (2 unique) | 2/2 | 0/0 | 0/0 | GWFAM39 |
| HSA | NBPF_core | 21→19 (19 unique) | 13/12 | 0/0 | 1/0 | GWFAM0,GWFAM1 |
| HSA | NBPF_repeat_bridge | 17→15 (15 unique) | 1/1 | 0/0 | 9/0 | GWFAM1 |

## Emitted evidence

The `*.fresh.copies.tsv`, `*.fresh.families.tsv`, and `*.fresh.er.*.tsv` files are direct
snapshots from the new Rustle runs. The two `graphs/*.fresh_emitted.gfa` files contain only
freshly emitted family copies and the actual fresh `E_r` edges between representative nodes;
no audit edge is synthesized. Colours are post-hoc labels: green named target, blue RNA
candidate, purple documented broad-family outgroup, red conflicting gene still in the target
family, orange conflicting gene emitted in another family, grey review/unscored.

| species | fresh family copies | mapped to E_r nodes | fresh E_r edges in GFA |
|---|---:|---:|---:|
| GGO | 44 | 44 | 73 |
| HSA | 165 | 165 | 349 |

Run-certificate disclosure: HSA: 1 co-located same-strand pair(s) were decided by reads_distinguish (O2's chi(H) predicate), not by sequence. O1's node set is NOT a function of sequence alone for this run — disclose it with any number derived from it.
Consequently, this fresh HSA node set is not a function of sequence evidence alone.

## Interpretation boundary

A re-emitted conflicting node is not an audit artifact: current upstream locus construction
produced it again. If it also enters the target's fresh component, the false merge is an O1
family-edge/typing problem. Failure to re-emit is evidence of representative-selection or
regional-substrate instability, not proof that the old node was biologically false.
A related outgroup is different: its broad-family edge is biologically expected, while its
exclusion applies only to the recent-copy subfamily view.

This remains a regional reconstruction. `samtools -M -L` preserves records from the original
whole-genome alignment but removes records outside the panel intervals. Therefore the experiment
validates node reproducibility and local regrouping, not byte-identical whole-genome partitioning.
