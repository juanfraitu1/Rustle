# Compact family-member similarity graphs

These are the small graphs requested for Bandage: **one node is one annotated
copy**, and one link represents a strong whole-locus sequence relationship.
They are intentionally distinct from `../dna_graphs/`, where segments are
base-level variation-graph pieces and copies are paths.

| file stem | nodes | retained links | colours |
|---|---:|---:|---|
| `ID_154_NPIP` | 14 | 25 | NPIPA blue, NPIPB red, NPIPP yellow |
| `ID_354_FAM72` | 4 | 6 | FAM72 green |

For every pair, minimap2 (`-x asm20`) reports alignment identity and reciprocal
coverage. A displayed link must have identity >= 0.80, reciprocal coverage >=
0.30, and `score = identity × reciprocal coverage >= 0.80`. All pairwise
values, including filtered pairs, are in `<stem>.pairwise_similarity.tsv`.
Every GFA link carries `ID:f`, `CV:f`, and `SC:f` tags.

In Bandage, open a `.gfa`, then load its matching `.colours.csv` and choose the
`Colour` column for custom node colours. Node names are copy names. The
`.legend.tsv` gives locus, class, length, and colour.

Regenerate from the repo root with:

```bash
python3 bench/soto/export_member_similarity_gfa.py
```
