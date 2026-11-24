# PREREG: the seeding pool on `ms` instead of AS (row 1942, 2026-11-24)

User: "the near-tied multimappers are those at 0.98 x AS, right? why not 0.8 x AS or 0.98 x ms". Two different rules are in
the pipeline: the assembler's seeding pool `good` (a secondary alignment helps build transcripts when its AS >= 0.98 x
the read's genome-wide best AS) and the near tie of the labels, the EM and the sequence test (`ms` within 10 of the
best). A wider AS band was tested and rejected (row 1767: additions 1-16% correct, loci on secondaries alone 70-75%
shadows; row 1773). The pool on `ms`, the score minimap2 ranks by, was never tested, although the user's rule of
2026-11-18 asks for the `ms` counterpart of every AS band.

## How (no change to the assembler)
`bench/o2_families/pool_bam.py` writes a BAM with only the secondaries a rule admits (best scores from the BAM's score
table); it is assembled with `--seed-pool all`. Arms:
- A0 = today's `good` pool, the normal command.
- C = `as:0.98` through the emulation: must reproduce A0 (the check that the emulation is faithful).
- M1 = `ms:0.98`: the same relative width on `ms`.
- M2 = `msgap:10`: the near-tie band of the labels.
Each arm: `assemble`, then stage `multicopy` (today's default, one floor-1 assembly with the bypass), then the labels
(`em_filter.py --label --weighted 1`) on the ORIGINAL BAM and score table, so the definition of a tie is the same in
every arm.

## Data and measures
- Planted human named families, fresh seeds 1951, 1952, 1953, scenarios E and S: copies and isoforms (`score.py`) of the
  seeded assembly and of the mode, units (`units_score.py`), labels (`label_score.py`).
- Real: OR6737 and KB3781 against the runs of row 1937: transcripts, rescued loci, the named families' annotated copies
  and junction precision (`toolcompare.py --species gorilla`), and how many secondaries each rule admits.

## Claims
- E0: C reproduces A0 (same transcripts by chain, >= 0.999). If not, the arms are not comparable and the row stops.
- An `ms` arm is BETTER than A0 when, on planted reads, expressed units found are not lower in any set, silent units
  labelled `expressed` or `shared` are not higher by more than 1 per set, and the mode's isoform precision or isoform
  recall is higher by >= 0.01 in the sum without the other falling by more than 0.005; AND on both real samples annotated
  copies are not lower and junction precision is not lower by more than 0.01.
- EQUIVALENT when every one of those is within the same tolerances in both directions: then `ms` is the cleaner
  definition at no cost, and the choice is the user's.

## Predictions (mine)
E0 holds. M2 (within 10) admits far fewer secondaries (the census of row 1829: near ties are a minority of the 0.98 pool)
and should raise isoform precision; whether it loses copies at low coverage is the open point. M1 should be close to
A0, differing where a spliced and an unspliced alignment of one read compete (an intron costs 24 on AS, nothing on
`ms`). I have been wrong on most predictions this week.

## Result (jobs 4503952 planted, 4503953 OR6737, 4503954 KB3781)
**E0 HOLDS.** The emulation reproduces today's pool exactly: assemblies identical on all six planted sets and on
OR6737 (70,181 = 70,181 transcripts), the mode identical on planted and 0.99937 shared on OR6737.

Secondaries admitted (share of all secondary alignments): planted AS 0.98 0.43-0.48, `ms` 0.98 0.55-0.63, within 10 `ms`
0.37-0.43; OR6737 0.018 / 0.010 / 0.005 of 6.27 M; KB3781 (M1, M2) 0.013 / 0.008 of 21.2 M.

Planted, the mode (seeds 1951-1953; E, S):

| | expressed units found | silent units labelled `expressed` or `shared` (S) | isoform recall | isoform precision, six sets |
|---|---|---|---|---|
| A0, AS 0.98 | 104 / 57, 104 / 55, 104 / 56 | 2, 4, 3 | 0.828 / 0.815, 0.828 / 0.846, 0.824 / 0.756 | 0.639, 0.432, 0.645, 0.433, 0.644, 0.421 |
| M1, `ms` 0.98 | the same | 2, 4, 3 | the same | 0.619, 0.423, 0.635, 0.424, 0.631, 0.421 |
| M2, within 10 `ms` | the same | 2, 4, 3 | the same | 0.664, 0.445, 0.676, 0.409, 0.679, 0.413 |

Real (`toolcompare.py`: the named families' annotated copies with a spliced transcript / junction precision / spliced transcripts at those copies; corrected 2026-11-25, the third number was first labelled "junctions"):

| | OR6737 | KB3781 | rescued loci OR / KB | transcripts with the mode OR / KB |
|---|---|---|---|---|
| A0 | 37 / 0.67 / 176 | 28 / 0.63 / 108 | 1,688 / 514 | 73,173 / 61,772 |
| M1 | 38 / 0.67 / 196 | 28 / 0.63 / 109 | 1,710 / 492 | 73,091 / 61,558 |
| M2 | 37 / 0.70 / 135 | 29 / 0.68 / 98 | 1,303 / 375 | 71,746 / 61,072 |

- **M2 (within 10 `ms`, the near-tie band) is BETTER by the registered rule.** Units and labels unchanged, isoform recall
  unchanged, isoform precision +0.072 over the six planted sets (mean +0.012); on both real samples annotated copies
  not lower (37 = 37, 29 against 28) and junction precision higher (+0.03, +0.05). What it costs, said plainly: precision
  is LOWER in two of the three half-silent planted sets (0.433 -> 0.409, 0.421 -> 0.413); it rescues fewer loci (1,303
  against 1,688; 375 against 514) and reports about 2% fewer transcripts.
- **M1 (`ms` 0.98) is not better and not equivalent:** isoform precision lower in five of six planted sets (-0.061 in the
  sum), real samples the same as A0. On planted reads the `ms` ratio admits MORE secondaries than the AS ratio, on real
  reads fewer.
- Predictions: E0 right; M2 raises precision, right, and loses no copies; M1 "close to A0", right on real data, slightly
  worse on planted.
- Not decided and not built. The pool rule lives in the assembler; this row used a filtered BAM. Adopting `within 10 ms`
  as the pool would (a) make the pool and the labels use ONE definition of a near tie, (b) need the rule in the Rust code
  (the score table already holds the best `ms`) and a rebuild, (c) by the working rules, one more unseen sample first
  (A119b or HG002; OR6737 and KB3781 are now used for this lever). A default: the user's call.

## Confirmation on unused samples (registered before the run, 2026-11-24)
A119b (human testis) and HG002 (human cell line), neither used for this choice. Arms: A0 = today's pool with today's
multi-copy stage (the stage rerun with the bypass on the runs of rows 1920 / 1931, whose assemblies are today's), M2 =
`ms within 10` through `pool_bam.py`. Measures as for the gorillas (`toolcompare.py --species human`).
CONFIRMED when on both samples M2's annotated copies with a spliced transcript are not lower than A0's and its junction
precision is not lower by more than 0.01. Prediction: confirmed, with junction precision up and fewer rescued loci.
