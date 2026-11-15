# Row 1682 — step 0: the SD pair's own identity inside the SD-held, low-SHARE cell (2026-11-13)

Goal (user, 2026-11-13): "explore what other ways can we include the segmental duplications so they help the
pipeline". Families layer (O1), a measurement on existing tables; nothing wired; not for `advisor_share`. No read is
assigned to a copy.

**Why this cell.** Counted before this registration on the row-1674 tables: nearly all refused truth pairs that have
an alignment record are SD-held and refused for low exon sharing (A119b 133 of 936, HG002 79 of 489, OR6737 rebuilt
105 of 661, KB3781 81 of 492). The false pairs there are mixed: about a quarter unrelated, the rest older Compara
paralogs and HGNC-group relatives. Rows 1674-1681 already searched SHARE, identity and `sd_held` together, so only an
SD attribute the alignment features do not carry can add anything.

**This reopens dating, on purpose.** Row 1679 closed dating for alignment-derived features and row 1680 said I would
propose no further score variant. This is not a score: it asks whether SEDEF's own identity for the holding SD pair
(fracMatch, over the whole duplicated block, not over our locus alignment) dates the pair. Row 996 found SEDEF to be a
lossy copy of our PAF, so a control arm is registered.

**Measure (frozen, `bench/denovo_eval/edgesd.py`).** In the cell only. `sd_id` = highest fracMatch among SD pairs that
hold the two loci (same holding rule as `sd_held`: 300 bp over each locus). Thresholds 0.90, 0.95, 0.98, 0.99.
Control = the same sweep by rank on the pair's PAF `identity`, read at the same number of picked pairs.

**Go (all three, else STOP):** a threshold chosen on A119b with >= 20 true pairs at precision >= 0.75; the same
threshold on HG002 with >= 20 true at >= 0.75; and on both samples precision at least 0.10 above the control.
**If GO:** I report and ask; a gorilla test and the family-level step are the user's decision.
**If STOP:** SD identity is closed as a pair-level gate.

## Result (2026-11-13, workspace `r1682_sdid/sdid.out`) — STOP

| sample | sd_id >= | picked | true | precision | control true (PAF identity, same picks) | false: unrelated / Compara / HGNC-only |
|---|---|---|---|---|---|---|
| A119b (cell 936, truth 133) | 0.90 | 184 | 39 | 0.212 | 18 | 102 / 28 / 15 |
| | 0.95 | 99 | 16 | 0.162 | 11 | 61 / 13 / 9 |
| | 0.98 | 42 | 5 | 0.119 | 5 | 24 / 4 / 9 |
| | 0.99 | 23 | 2 | 0.087 | 1 | 15 / 1 / 5 |
| HG002 (cell 489, truth 79) | 0.90 | 113 | 35 | 0.310 | 29 | 60 / 15 / 3 |
| | 0.95 | 58 | 25 | 0.431 | 26 | 29 / 3 / 1 |
| | 0.98 | 20 | 7 | 0.350 | 8 | 10 / 2 / 1 |
| | 0.99 | 8 | 1 | 0.125 | 1 | 5 / 1 / 1 |

No threshold reaches 0.75; precision falls as the SD gets younger. Reading: inside this cell a young SD mostly holds
unrelated genes that were duplicated together in one block (the "unrelated" column dominates the false picks), and 94
of the 133 A119b truth pairs sit in SD pairs below 0.90 identity. SD identity dates the block, not the gene pair, and
it is no better than our own PAF identity on HG002 (row 996 again). Closed as a pair-level gate. One implementation
note: fracMatch is column 21 in the CHM13 and gorilla files and column 22 in the HG002 build; `edgesd.py` fixes the
column per file from its first row.

## Follow-up count (2026-11-13, not registered, no row): a lower SHARE threshold for SD-held pairs

Asked after rows 1683-1685: would the family curve gain from admitting SD-held pairs at SHARE below 0.5? Pair-level
count on the same tables before building anything (refused low-SHARE pairs, true / all):

| sample | SD-held, SHARE >= 0.2 | >= 0.3 | >= 0.4 | not SD-held, SHARE >= 0.2 |
|---|---|---|---|---|
| A119b | 54 / 260 (0.208) | 20 / 104 (0.192) | 7 / 46 | 5 / 131 (0.038) |
| HG002 | 33 / 139 (0.237) | 14 / 41 (0.341) | 1 / 9 | 0 / 94 |
| OR6737 rebuilt | 42 / 214 (0.196) | 15 / 76 (0.197) | 2 / 15 | 1 / 110 |
| KB3781 | 27 / 145 (0.186) | 3 / 36 (0.083) | 1 / 7 | 1 / 81 |

SD-held pairs are about five times richer in true pairs than the rest, but at about 0.2 they are below what lowering
MCL inflation already buys on the default curve (A119b, I 2.0 -> 1.6: about 22 more true pairs for 31 more false,
0.42). An SD-conditional SHARE gate would sit under the existing curve, so no binary flag was built and no row taken.
