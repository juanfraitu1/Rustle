# PREREG r1821: resolution-aware output of multi-copy families, judged once on the last HG002 quarter (2026-11-19)

Bars accepted by the user on 2026-11-19 ("lets accept the bars, spend the residues and commit changes"), before residues
75-99 were selected. Workspace `~/ws/r1821_resolution_call` (PSU; claimed after the register tail on origin/main was 1820).
Goal: `docs/GOAL_multicopy_families_2026-11-19.md`. Development record: workspace `goal_multicopy_dev/LEDGER.md` (seen HG002
quarter, residues 0-24; simulated arrays). Judging set: HG002 residues 75-99, the last unseen quarter, opened once.

## Why
The per-copy output's main error on the seen quarter is NARROW: one copy claimed where the reads fit several (whole-sample
GTF 187 of 1,613 judged one-copy claims). The `famwise` layer fixes that but joins separable copies (copy-set claims 0.42
EQUAL, mostly WIDE). A call per family transcript from the reads alone (`bench/famwise_eval/resolution_call.py`: keep a
copy set only when >= half its supporting molecules tie across its places) raised both precision and resolved units on
development, and makes identical tandem arrays one copy set over all units.

## Measure (fixed)
`bench/famwise_eval/copyset_score.py --match contained` (a read chain supports a transcript holding it as a contiguous run
of introns): claim verdicts EQUAL / NARROW / WIDE / MIXED / OFF by the supporters' strict plurality (a tie = MIXED);
"origin" = the HG002 copies a read fits equally (the truth's origin set), not a literal source; precision = EQUAL / judged
claims over one-copy and copy-set claims together, with the molecule-weighted share and the claim counts beside (grouping
copies shrinks the claim count); resolved units = true units (>= 2 molecules) held by one claim with exactly their copies. The
exact-match version is printed beside, no bar.

## Arms (one selected BAM; the frozen code of the development ledger)
- P: `assemble --seed-pool primary`. D: `--seed-pool good`. S: `--seed-pool good` with `RUSTLE_SEED_SPLICED=support`
  (spliced secondaries as support only; development: same as the post-filter emulation, 0.850 vs 0.851). StringTie `-L` on P's BAM: the trusted baseline.
- For each of P and S: `eq_classes --famwise` -> `famwise` -> `resolution_call --tied-share 0.5`, window 0 (W0)
  and 10 (W10). Four candidate outputs: P_W0, P_W10, S_W0, S_W10.

## Readings (accepted 2026-11-19)
Primary candidate, named on development before the quarter is opened: **P_W10** (P -> famwise -> call, window 10 ms,
share 0.5). The other three candidates are printed beside, no bar.
- **R1 (the call adds to the layer).** Against P_fw (P -> famwise, no call): P_W10's precision interval lower bound >=
  P_fw's point estimate minus 0.01, AND resolved units >= P_fw's x 1.01, AND one-copy resolved units >= P's x 0.99.
  Development (corrected scorer, contained): P_fw 0.913 / 1,056; P_W10 0.914 / 1,075 (+1.8%); one-copy 1,037 vs P 1,039.
- **R2 (against the per-copy output).** P_W10's precision >= P's and resolved units >= P's x 1.02. Development: 0.914 vs
  0.896, 1,075 vs 1,039 (+3.5%).
- **R3 (copy sets are honest).** P_W10's copy-set claims: molecule-weighted EQUAL share >= 0.8 (development 0.96), claim
  share printed with its interval (development 17/30). Fewer than 30 judged copy-set claims: R3 undecidable on claims.
- **Beside, no bar:** StringTie P; the exact-match table; P_W0, S_W10, S_W0, S_fw; the six array simulations rerun with
  the frozen code, scored on the fits-equally and on the literal-source truth.

## Prediction (mine)
R1 narrow PASS or BETWEEN (+1-2% over the layer), R2 PASS, R3 weighted PASS with too few claims for the claim share.
Development after the scorer fix: S_W10 0.901 / 1,085, P_W0 0.908 / 1,068, S_W0 0.870 / 1,072.

## Data budget
Residues 75-99 are the last unseen HG002 quarter; the user chose to spend them on this row (questions 1 and 3). After it
no unseen HG002 quarter remains; the parked cross-copy test's ruling (residues 50-99 at its own call loci) now meets only
spent reads, said in its record when it runs.

## Question 3 beside (the goal's rule: new and lost chains apart)
`copyset_score.py --new-lost P`: S's and D's claims absent from P (new) and P's claims absent from them (lost), by verdict.
Printed, no bar (development: S new EQUAL 28 / NARROW 105 / OFF 16, D new 50 / 335 / 103, both lose 25).

## Freeze
Snapshot `~/ws/r1821_resolution_call/snap/` (tools, bench, the S binaries of `goal_multicopy_dev/target/release`, copied);
hashes in `snap/FROZEN.sha256`, written before the judging job is submitted.

## Not in this row
Question 4 (borrowed junction columns): development ceiling ~8% of units, mostly 2-read isoforms; not built.

## Result (2026-11-19, job 57497220 on PSU; workspace `judge/`)
Snapshot hashes checked by the job before selection. 136,129 reads selected (residues 75-99); true units one 1,127 / set 103.
Contained matching (bars):
| output | precision [Wilson 95%] | resolved one + set | set claims EQUAL / judged (weighted) |
| P | 1894/2132 = 0.888 [0.874, 0.901] | 1019 + 0 = 1019 | - |
| P_fw | 1879/2068 = 0.909 [0.895, 0.920] | 995 + 39 = 1034 | 19/62 (0.850) |
| **P_W10** | 1903/2093 = 0.909 [0.896, 0.921] | 1011 + 33 = 1044 | 17/37 (0.939) |
| P_W0 | 0.901 [0.888, 0.913] | 1019 + 29 = 1048 | 13/17 (0.993) |
| S_fw | 0.894 | 995 + 46 = 1041 | 33/107 (0.830) |
| S_W10 | 0.893 [0.879, 0.905] | 1011 + 43 = 1054 | 29/69 (0.918) |
| S_W0 | 0.869 | 1018 + 30 = 1048 | 17/29 (0.981) |
| StringTie P | 0.861 | 919 + 0 | - |
- **R1 FAILS** (narrowly, on both counts): precision lower bound 0.8962 < 0.8986 (P_fw 0.9086 - 0.01); resolved 1,044 <
  1,044.3 (P_fw x 1.01); one-copy 1,011 >= 1,008.8 holds. Over the uncalled layer the call adds 10 units (+1.0%).
- **R2 PASSES**: precision 0.909 >= 0.888; resolved 1,044 >= 1,039.4 (+2.5% over P).
- **R3 PASSES** on the weighted share (0.939 >= 0.8); 37 judged copy-set claims, claim share 17/37 = 0.46.
- My prediction (R1 narrow pass or between, R2 pass, R3 weighted pass) was wrong on R1.
Beside, question 3, new / lost against P: D new EQUAL 51 / NARROW 298 / OFF 60 (lost 36, EQUAL 17); S new EQUAL 23 /
NARROW 88 / OFF 15 (lost 32, EQUAL 15): support-only keeps ~half the right new chains and cuts wrong-copy ones 4x, as on
development. Arrays (frozen code): as on development (identical units: every arm puts chains on every silent unit, the call
makes one copy set; S removes silent-unit chains at 0.1-0.2% divergence where D does not).
Reading (mine): grouping copies is what raises precision and resolution over the per-copy output (R2); the call's own
contribution over the layer is +1% and not distinguishable from the layer here. The exact-tie call (W0) gives the most
honest copy sets (weighted 0.99) and the most resolved units of the P arms (1,048) at precision 0.901.
