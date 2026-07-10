# Design: Multi-copy gene family benchmark oracle

Status: DESIGN (approved 2026-05-31). Next: implementation plan (writing-plans).
Grounding: VG objectives/state (`VG_OBJECTIVES_AND_ROADMAP.md`, `project_vg_wiring` memory); existing
benchmark assets (`bench/multi_copy_eval/OBJECTIVES_ASSESSMENT.md` + scorers; `test_data/synthetic_family/`);
the DAZ verification + canonical-k-mer fix (`docs/experiments/DAZ_vg_verification.md`).

## 1. Problem & insight

The thesis's core objective is multi-copy gene family assembly — what StringTie cannot do. Unlike the
StringTie-parity work (measurable via gffcompare), this has no off-the-shelf ground truth, so progress
risks being anecdotal. GROUNDING REVEALED the benchmark is already ~80% built — `bench/multi_copy_eval/`
has `OBJECTIVES_ASSESSMENT.md` (all 4 objectives, Rustle-vs-StringTie numbers), `no_absent_copy_eval.py`
(Obj 2, re-runnable: GOLGA8 13/14, YAG 21/21 vs StringTie 71%/81%), per-family scorers, and a
deterministic synthetic fixture (`test_data/synthetic_family/`, `truth.gtf` + read-name copy-truth,
`copy_confidence` GTF attribute emitted). BUT it is scattered + partly static markdown, predates today's
DAZ canonical-k-mer fix, and has measurement GAPS: no formal synthetic per-copy isoform Sn/Pr (Obj 3
stays "descriptive") and no consolidated assignment-accuracy score (Obj 4).

**Insight:** the work is NOT to build a benchmark — it is to consolidate the existing pieces into a
single re-runnable, all-objective REGRESSION ORACLE, fill the two measurement gaps, and refresh with
DAZ. This makes every objective provable with numbers (for the advisor) and turns regressions into a
non-zero exit code (like the parity bench tools).

## 2. Goal & non-goals

- **Goal:** `bench/multi_copy_eval/run_oracle.py` measures all 4 objectives, scores actuals vs an
  expectations file, prints one advisor-ready table + per-objective PASS/FAIL, exits non-zero on
  regression, and regenerates `OBJECTIVES_ASSESSMENT.md` from live data.
- **Decision rule:** analysis-only; no production/VG-algorithm changes; default 95.6/90.5 reported
  separately as an isolation guard.
- **Non-goals:** the Obj-1 low-depth simulation (future), the intra-bundle splitter, any VG algorithm
  change. This is the MEASUREMENT layer only.

## 3. Architecture — two tiers + expectations

- **FAST tier (synthetic, deterministic, ~seconds; the everyday regression gate):** `rustle --vg` on
  `test_data/synthetic_family/` scoring:
  - Obj 3 — per-copy isoform Sn/Pr (NEW) vs `truth.gtf`.
  - Obj 4 — read-to-copy assignment accuracy (NEW) vs read-name truth.
- **SLOW tier (real families, ~minutes, on-demand):** wraps `no_absent_copy_eval.py --check` (Obj 2),
  adds the DAZ canonical-fix run (DAZ3 ≥5 isoforms + EM-ran), and GOLGA6L7 as an explicit expected-FAIL
  (the same-bundle Obj-1 blocker, tracked so it flips to PASS when the splitter lands).
- **Expectations file** `bench/multi_copy_eval/expectations.json`: per objective/dataset, expected value
  + tolerance (exact for synthetic counts; `≥` floors for recovery %; `±` band for the de-novo
  headline; DAZ3 `≥5`; GOLGA6L7 `expected-fail`). `--check` diffs actual vs expected → PASS/FAIL/
  REGRESSION, non-zero exit on any FAIL (SKIPPED rows don't fail).

## 4. Components

### `bench/multi_copy_eval/run_oracle.py` (new — the driver)
- Modes: `--fast` (synthetic only), `--full` (+ real panel), `--check` (compare to expectations.json,
  exit non-zero on regression). Orchestrates scorers, collects results, diffs, prints the table, and
  (with `--write-report`) regenerates `OBJECTIVES_ASSESSMENT.md`.

### `bench/multi_copy_eval/score_synth_isoforms.py` (new — Obj 3)
- Run `./target/release/rustle --vg test_data/synthetic_family/reads_sorted.bam -o /tmp/synth_vg.gtf`
  (+ a non-VG baseline). Partition `truth.gtf` and predictions by copy via gene span (copy A: 1000–4500,
  copy B: 10000–13500). Per copy compute intron-chain Sn/Pr (reuse `bench/gtf_chain_diff.py` chain
  logic). Output `{copyA:{Sn,Pr}, copyB:{Sn,Pr}, overall, vg_vs_baseline_delta}`.

### `bench/multi_copy_eval/score_synth_assignment.py` (new — Obj 4)
- Run with `RUSTLE_VG_FP_ATTR_TSV=/tmp/attr.tsv --vg-snp`. Parse rows
  `(family_id, read_name_hash, placement_copy, n_sites_covered, final_weight, weight_gap, weight_sum)`.
  Build `hash→true_copy` by FNV-1a-hashing the BAM read names (`uniq_A*`/`uniq_A2*`→A, `*B*`→B,
  multimap-origin per the generator) the SAME way the dump hashes them. Score: of multi-mappers,
  argmax-weight copy == true copy → accuracy; bucket decisive (gap>0.8)/uncertain (≤0.5); report
  overall accuracy + accuracy-among-decisive.

### Reused (existing)
- `no_absent_copy_eval.py --check` (Obj 2); `classify_per_family.py` / `per_paralog_breakdown.py`
  (real-family per-copy tables); a small DAZ wrapper (the committed subset command, count DAZ3 isoforms,
  assert the `HMM-EM ... reweighted` log line).

## 5. Data flow

```
expectations.json ─┐
synthetic fixture ─▶ rustle --vg ─▶ score_synth_isoforms (Obj3) ─┐
 (truth.gtf, name-truth)        └─▶ score_synth_assignment (Obj4)─┤
real panel ─▶ no_absent_copy_eval (Obj2); DAZ wrapper ───────────┼─▶ run_oracle.py ─▶ table + PASS/FAIL + exit
GOLGA6L7 ─▶ expected-FAIL (Obj1 blocker) ────────────────────────┘                  └─▶ regenerate OBJECTIVES_ASSESSMENT.md
default de-novo GGO_19 ─▶ headline (isolation guard, reported separately)
```

## 6. Validation (of the oracle itself)

- **Scorers validated against hand-checked synthetic values:** the deterministic fixture has a knowable
  correct Obj-3 Sn/Pr and Obj-4 accuracy on first run; the plan captures and asserts them.
- **Determinism check:** FAST tier run twice → identical output (catches flakiness).
- **Expectations seeded** from the first clean run; `--check` then flags drops below floor as REGRESSION.
- **Default isolation:** the de-novo headline (95.6/90.5) is reported each run so VG scoring never masks
  a default regression.

## 7. Safety & exit

- Analysis-only — no production/VG code changes; zero risk to the default operating point.
- Graceful degradation: FAST tier needs only the in-repo synthetic fixture (always runs / CI-able);
  SLOW-tier rows print `SKIPPED (data missing)` when the large real BAMs are absent, never hard-error;
  per-invocation exit codes checked (a crash → `ERROR` row, surfaced not swallowed).
- `run_oracle.py --check` returns 0 iff every measured objective meets its expectation. Output is the
  advisor-ready table + a regenerated `OBJECTIVES_ASSESSMENT.md`.
- This is the measurement layer that makes the next VG improvements (intra-bundle splitter, absent-copy
  validation, Obj-1 low-depth) provable.
