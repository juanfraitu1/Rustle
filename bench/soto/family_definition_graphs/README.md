# Rustle family-definition graphs

These graphs are the correct companion to the compact similarity plots: every
node is an annotated Soto copy and every link is an **E_r family-definition
edge**. They implement the nucleotide rules used by Rustle's
`homology_edges_all_reps` engine:

- asm20 tier: identity >= 0.80 and coverage of the shorter member >= 0.50; or
- sensitive tier (`-k11 -w5`): identity >= 0.60 and coverage of the shorter
  member >= 0.50.

The final family condition is a connected E_r block whose internal density is
at least gamma = 0.20. Both exported families exceed that condition maximally:

| family | copies | E_r edges / possible edges | density | components | certificate |
|---|---:|---:|---:|---:|---|
| NPIP | 14 | 91 / 91 | 1.000 | 1 | PASS |
| FAM72 | 4 | 6 / 6 | 1.000 | 1 | PASS |

Thus these are **cliques**, which is stronger than the gamma-quasi-clique
criterion. The matching `.certificate.tsv` is the machine-readable statement;
`.edges.tsv` records the qualifying tier, identity, coverage, and score for
every edge. Load the `.gfa` and matching `.colours.csv` in Bandage.

Regenerate from the repo root:

```bash
python3 bench/soto/export_family_definition_gfa.py
```
