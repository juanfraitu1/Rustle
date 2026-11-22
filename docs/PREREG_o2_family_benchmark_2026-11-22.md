# PREREG: O2 as a benchmark of known gene families across tools, on planted reads of known origin (row 1903, 2026-11-22)

User (2026-11-22, goal): "replace the O2 machinery by this newly discovered modality that looks at what do known gene
families actually perform under different tools, use simulation for planting copies and having some provenance truth".
Written before any read of this row is planted. Seen: rows 1855-1894 (real-sample family tables per tool; planted reads
scored for read sharing, joins and aligner settings, with Rustle and StringTie only; never the other assemblers on
planted reads, never isoforms, never two isoforms per copy).

## What O2 is from this row on
O2 was: assign an AS-tied read to a copy on PSVs, or abstain (`assign` stage, `copy_assign`, truths on haplotype
assemblies; `LEDGER_O2.md`). Its record: no rule beats the primary alignment, nothing reads the assignments.
O2 is now: **how completely and how correctly each tool reports the copies and isoforms of known multi-copy families**,
measured (a) on planted reads, where the copy and isoform of origin of every read is known, and (b) on real samples
against the annotation (`bench/denovo_eval/toolcompare.py`, rows 1859-1869). Code: `bench/o2_families/`
(`plant.py`, `run_tools.sh`, `score.py`, `tied_share.py`). The old machinery stays in the repository, runnable, as legacy.
Amendment (2026-11-22, after the planted run, on the goal's wording "replace"): in the driver `all` runs the new stage `o2`
(family report of the run) in place of `assign`; `RUSTLE_O2_LEGACY=1`, `--candidates` or `--legacy-catalog` keep `assign`;
stage `o2bench` runs this planted benchmark. Dispatch tests changed with it.

## This run (workspace `r1903_o2_family_benchmark`)
- Planted: the 133 protein-coding copies of `MAINSTAYS` + `YAGS` on CHM13, up to 2 annotated isoforms per copy with
  distinct intron chains, reads 2 : 1; seed 1903; errors 0.1% / 0.2% / 0.2%; 70% full length. Scenarios E (60 reads per
  copy), D (one copy 600, the others 20), S (half the copies of each family silent).
- One BAM per scenario (`minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1`), the same for every tool.
- Tools: Rustle with primaries only, Rustle with near-tie secondaries, Rustle's multi-copy mode, StringTie -L, FLAIR,
  IsoQuant, bambu (NDR 1), each de novo with the recipe used on the real samples. A tool that fails on such a small BAM
  is reported as failed, not replaced.
- Measures, in pairs: copy recall (expressed copies with a transcript) and copy precision (called copies that are
  expressed; meaningful in S only); isoform recall (planted isoforms reproduced exactly at their copy) and isoform
  precision (transcripts on planted copies that are such an isoform); plus joined transcripts and transcripts elsewhere.
  Read level: reads placed on the true copy / another copy / tied.

## Claims under test
- C1 (the multi-copy mode's claim, rows 1867-1869, now with truth): the mode's copy recall is >= every other tool's in
  E, D and S. FAILS if any tool is higher in any scenario.
- C2 (its cost): in S the mode's copy precision is not lower than StringTie's by more than 0.05. FAILS otherwise.
- C3: the mode's isoform recall is >= StringTie's in all three scenarios, with isoform precision not lower by more
  than 0.05.
All three are reported whatever they say; none tunes anything.

## Predictions (mine)
C1 holds, but narrowly: with every read's primary tossed among identical copies, every tool reaches most expressed
copies in E (recall > 0.9 for StringTie too). C2 FAILS: the mode calls more silent copies than StringTie (near-tie
secondaries land on silent identical copies). C3 holds. FLAIR and IsoQuant: lower isoform precision than Rustle.
Simplifications as rows 1874-1875: the sample equals the reference, simple error and truncation models, three
scenarios; a scenario, not a mimic of real data.

## Result (2026-11-22, job 4495069): C1 FAILS, C2 FAILS, C3 holds on recall and FAILS on precision
All seven tools ran on all three scenarios. Tables: `O2_FAMILY_BENCHMARK.md`, `r1903_o2_family_benchmark/{E,D,S}/score.tsv`.
- C1: copy recall, mode against the best other tool: E 0.970 = StringTie and bambu 0.970; D 0.970 against StringTie
  0.977 (129 against 130 copies); S 0.971 against StringTie 0.986 (68 against 69). FAILS in D and S by one copy.
- C2: S copy precision mode 0.618 (42 of 63 silent copies called), StringTie 0.734 (25): lower by 0.116. FAILS.
- C3: isoform recall mode 0.791 / 0.741 / 0.750, StringTie 0.649 / 0.636 / 0.648: holds. Isoform precision 0.677 /
  0.592 / 0.409 against 0.791 / 0.776 / 0.606: lower by 0.11-0.20. FAILS.
- Every tool calls 18-25 silent copies with primaries alone; near-tie secondaries add 17.
- Reads: placed alone on the true copy 4,585 of 4,665 (E), 5,868 of 5,906 (D), 2,490 of 2,503 (S); tied 3,315, 5,454, 1,697.
Predictions: "C1 holds narrowly" wrong (it fails narrowly); every tool above 0.9 copy recall in E wrong for FLAIR
(0.887); C2 fails right; C3 holds wrong on precision; FLAIR and IsoQuant lower isoform precision than Rustle wrong (they
are higher than Rustle with secondaries).
