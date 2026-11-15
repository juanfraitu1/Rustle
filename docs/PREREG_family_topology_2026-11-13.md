# Row 1681 — family topology as the witness for a refused pair (2026-11-13)

Goal: "find empirical rules to improve edge recall in the multi-copy gene families definition". User's choice after
row 1680 stopped ("Test family topology"). Families layer (O1), opt-in only if found, not for `advisor_share`. No read
is assigned to a copy; families and `fam.copies.*` are not changed by this row. No score: rows 1675-1680 closed it.

**Mechanism tested.** A gene that belongs in a family should align to several of its members, not to one. Registered
risk: large old families (KRAB-ZNF) align to each other widely and satisfy this too.

**Rule (frozen, `bench/denovo_eval/edgetopo.py`, on the row-1674 `edgefeat.py` tables).** family(b) = b plus the genes
linked to b in the table (same default family). For a refused pair (a, b): support = members of family(b) with which
a has an alignment row; frac = support / |family(b)|. The pair is accepted when, in either direction, |family| >= 2,
support >= K and frac >= F. Arms: K in {2, 3} x F in {0.5, 0.67, 1.0}.

**Go/no-go (A119b + HG002 pooled).** Arm = most true pairs among arms with precision >= 0.75 and >= 20 true. None ->
STOP, gorilla not scored.
**Test (frozen arm; rebuilt OR6737 + KB3781 pooled):** >= 20 new true pairs at precision >= 0.60.
**If it passes:** row 1674's step 2 (family level through `--extra-edges`, same bars).
**If it stops or fails:** the row closes.
**Limit stated now:** family(b) is read from table rows, so members without a row to b are not counted.

## Result (2026-11-13, workspace `r1681_topology/topo.out`) — STOP

No arm comes near the bar. A119b + HG002 pooled, 370,077 refused pairs, 277 truth:

| K | F | true | false | precision |
|---|---|---|---|---|
| 2 | 0.5 | 160 | 5,983 | 0.026 |
| 2 | 0.67 | 157 | 4,828 | 0.031 |
| 2 | 1.0 | 146 | 4,580 | 0.031 |
| 3 | 0.5 | 98 | 2,644 | 0.036 |
| 3 | 0.67 | 95 | 2,112 | 0.043 |
| 3 | 1.0 | 68 | 1,829 | 0.036 |

Reading: topology has reach (over half of the refused truth pairs touch several members of the partner's family) but
no precision: thousands of unrelated or older-relative genes align to several members as well. Gorilla not scored.
Row closed.
