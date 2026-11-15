# PREREG r1660: 3' end = the commonest full-chain read end (registered 2026-11-12, before the A119b run)

User 2026-11-12: "yes, lets register the 3' end rule". Source: `docs/STEP0_tool_match_2026-11-11.md` §10 (gorilla OR6737,
discovery): the default emits the most extended full-chain read end on both sides; against RefSeq ends that is the best 5'
rule (0.351 within 10 bp) and the worst 3' rule (0.532; IsoQuant 0.661, the mode of our own reads 0.668).

## Rule (frozen)
A multi-exon transcript's 3' end = the commonest 3' end among the primary alignments carrying its exact intron chain; ties
go to the more extended end. The 5' end is not touched. (FLAIR's windowed peak is reported beside it, not eligible.)

## Test: A119b on CHM13 (`stage/A119b.t2t.bam`, default `audit/g2/r76e/a1.gtf`, RefSeq), run once, `bench/famcmp/end_rule.py`
Population: every RefSeq-correct chain of the default with a full-chain primary read. Truth: nearest RefSeq transcript
with the same chain. A119b ends were never scored before.
- T1: 3' ends within 10 bp, rule minus emitted >= +0.05 (gorilla +0.136).
- T2: 3' ends within 50 bp and within 100 bp, rule >= emitted (no loss at the coarse scale).
PASS = both. PASS -> ask the user whether to port it (default or opt-in; it changes every multi-exon model's 3' coordinate,
no chain). FAIL -> closed, no second rule.
Limits: RefSeq 3' ends are one annotated site per transcript; an expressed alternative polyA site counts as a miss for
every rule alike. The measure is on correct chains only.

## Result, A119b (2026-11-12, job 4468217, `r1660/a.txt`; first submission crashed on a NameError before any score): FAIL T1 and T2
38,533 correct chains of the default with a full-chain read.
| ends | 5' within 10 / 50 / 100 | 3' within 10 / 50 / 100 |
|---|---|---|
| ours emitted (= most extended) | 0.153 / 0.458 / 0.620 | 0.495 / 0.556 / 0.580 |
| rule: mode, ties to longer | 0.332 / 0.690 / 0.794 | 0.459 / 0.503 / 0.527 |
| median | 0.324 / 0.708 / 0.803 | 0.424 / 0.473 / 0.501 |
| FLAIR window peak | 0.336 / 0.717 / 0.807 | 0.457 / 0.500 / 0.526 |
| 10th percentile inward | 0.244 / 0.624 / 0.757 | 0.519 / 0.564 / 0.586 |
| IsoQuant emitted (where it has the chain) | 0.126 / 0.441 / 0.624 | 0.580 / 0.635 / 0.657 |
T1: 3' within 10 bp 0.459 vs 0.495 = -0.036 (bar +0.05). T2: lower at 50 and 100 bp. Closed as registered, nothing ported.
Post hoc reading: the gorilla pattern does not transfer and partly reverses. On A119b the commonest end is the better 5'
rule (+0.18 within 10 bp) and the worse 3' rule; on gorilla the opposite. IsoQuant's 3' ends are ahead of ours on both
samples (0.661 vs 0.532, 0.580 vs 0.495), but no simple statistic of our full-chain read ends reproduces that on A119b.
