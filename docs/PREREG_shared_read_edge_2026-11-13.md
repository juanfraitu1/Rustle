# PREREG row 1664: shared reads as a second witness for a weak alignment edge (2026-11-13)

**Why.** Member funnel (`bench/denovo_eval/member_funnel.py`): on A119b 237 expressed spliced truth members have a model
but no family. For their best truth partner the edge rule stops at: no record 121, SHARE < 0.5 83, passes 23, other 10;
223 of the 237 share reads with that partner (>= 5 molecules and >= 10% of the smaller gene's molecules).

**Already closed:** read sharing alone (row 1512: precision <= 0.032, unrelated genes dominate) and SHARE over the union
of records (row 1662 step 0c). **Seen before registration:** the truth side above. Not seen: any count on false pairs.

**Step 0 (A119b, `port1593` tables, no pipeline change).** Rule: a gene pair whose loci have an alignment record that
passes identity and length but fails SHARE (`edgemiss.py` stage `low_share`) gets an edge when the two genes also share
reads: molecules with a record on the exons of both, >= 5 and >= F of the smaller gene's molecules, F in {0.1, 0.3, 0.5}.
Pairs counted as pairscore counts them (expressed coding genes, different symbols, unlinked today).
**Stop rule (frozen):** go on only if some F gives >= 30 new true pairs at new-pair precision >= 0.75. Otherwise STOP.
A go = prototype behind a flag, gorilla KB3781 as the transfer test, families layer only.

## Result: STOP

2,639 unlinked low-SHARE pairs, 166 true. Shared reads do not discriminate: they are present for almost every pair.

| F | new true | new false | precision | same, gene spans that overlap set aside |
|---|---|---|---|---|
| 0.1 | 166 | 1,742 | 0.087 | 164 / 1,022, 0.138 |
| 0.3 | 163 | 1,411 | 0.104 | 161 / 767, 0.173 |
| 0.5 | 158 | 1,212 | 0.115 | 156 / 623, 0.200 |

Two loci similar enough to give an alignment record also attract the same reads, whether or not the genes are young
paralogs (false examples: ATAD3A/B/C, HNRNPCL1 x HNRNPC, ANKRD65 x ZNF888). Read sharing adds nothing to the record.
