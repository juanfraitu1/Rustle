# PREREG row 1703 (2026-11-14): the informative-secondary pool for candidate copies, re-registered

User go 2026-11-13 ("Start r1513 re-registration"). Row 1513 failed on its wall bar only (447 vs 416 s; the bar asked for
half the wall) and passed both quality bars on the r1509 40-family panel (absent_from_ref 6 vs 3, false flags 19 vs 21).
Its outcome named the follow-up: re-register without the speed bar, on more families, wall reported, not barred.

## Run (everything else as row 1513, `docs/archive/2026-11/PREREG_secondary_pool_2026-11-03.md`)
- Sample KB3781, catalog `~/runs/gorilla_KB3781/kb`, R = 0.5 (the r1511 bound, unchanged).
- Families: HELD OUT from row 1513: `panel.py --top 80` on the same molecules table, minus the 40 families of
  `~/o2contested/panel/families.txt` (ranks 41-80 of the same ranking).
- Arms, same binary (today's `candidate_copies`, copied to the workspace), wall limit 2 h each: A default pool,
  B `--molecules ... --min-ms-ratio 0.5`.
- Scoring `bench/o2_contested/pool_score.py`, unchanged classes. The catalog copies' haplotype loci for these families
  are made as `haplotype_truth.py` makes them (`splice:hq -uf --eqx -N 10`, best NM per haplotype;
  `bench/o2_contested/copy_hits.py`), since the old hits table covers the old panel only.

## Bars (written before the run)
1. B keeps at least 80% of A's absent_from_ref loci (same haplotype contig, within 1 kb).
2. B's false flags (allele_of_known_copy + not_in_genome) are no more than A's.
Fewer than 3 absent_from_ref contigs in arm A: UNDERPOWERED, report only.
Both hold = PASS: with row 1513's quality result this is two panels, and the informative pool is proposed as the pool
of the opt-in `candidates` stage. The stage's default is the user's decision and is not changed by this row. One fails
= FAIL: the flag stays opt-in and the question is closed.
Reported without a bar: wall per arm; contigs per class; contigs whose haplotype locus recurs in two or more families
(row 1513's "repeat-driven clusters" caveat), per arm.

## Outcome
**FAIL as registered, on both bars** (job 4470984, workspace `r1703_secpool`; 40 held-out families, ranks 41-80).

| arm | flagged | absent_from_ref contigs (loci) | allele_of_known_copy | not_in_genome | between | wall |
|---|---|---|---|---|---|---|
| A default pool | 46 | 12 (7) | 6 | 21 | 7 | 499 s |
| B `--min-ms-ratio 0.5` | 39 | 7 (3) | 8 | 20 | 4 | 423 s |

- Bar 1: B keeps 2 of A's 7 absent_from_ref loci (0.29; bar 0.80). B has one locus A lacks.
- Bar 2: false flags B 28 vs A 27 (bar: no more than A).
- Not barred: B drops 307,685 of 407,516 secondary records and is 15% faster. Contigs on a haplotype locus shared with
  another family's contig: A 22 of 46, B 13 of 39. One maternal region (CM054600.2: 6.73-6.80 Mb) holds 6 of A's 12 and
  6 of B's 7 absent_from_ref contigs, from 5 families: the "repeat-driven clusters" caveat of row 1513 again. Both
  loci B keeps are in that region. Outside it B keeps 0 of A's 5 loci: at 4 B has no contig of any class, at 1 a
  not_in_genome contig. So B loses them; it does not reclassify them.
- Copy loci check: 123 copies in the 40 families, all placed on a haplotype, 121 on both.
- Scoring note: the copy loci cover these 40 families only (as in row 1513, its own panel only), so a contig on another
  family's copy counts as absent_from_ref in both arms alike.

Reading: row 1513's quality result (6 vs 3) does not repeat on held-out families; here the default pool finds more (outside one repeat region, 5 loci vs 1).
`--min-ms-ratio` stays an opt-in flag, the `candidates` stage is unchanged, and the question is closed.
