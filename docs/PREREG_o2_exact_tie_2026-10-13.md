# PREREG — O2 exact ties decided by rounding noise: the fix (2026-10-13)

**Found by** the anchor-trim truth test (`docs/PREREG_o2_anchor_2026-10-13.md`, human held-out). Three MCL1623
molecules match both candidates at exactly 5 of 10 PSV columns. They were `assigned`, with margin 0.000 and p 2.5e-13,
and all three were wrong.

**Root cause.** Every O2 certificate decides on log-likelihood differences built as floating-point sums in column order:
- the read-star pair LLR, which drives both the maximin best-copy pick and the "LLR > 0" assign rule;
- the population `fixed_only_certificate`;
- the legacy per-read `read_copy_evidence` logl, which drives the argmax and the `margin > 0` rule.

When two copies' evidence is exactly balanced, the true difference is 0. The computed difference can be off by one
rounding error: the order `+++--+--+-` of ±lr sums to 1.78e-15. So "> 0" and the best-copy pick are decided by rounding
noise. The design itself says exact twins must fall through to the tie key, in the read-star maximin comment: "twins
(0 both ways) fall through to the old key".

In the shipped runs, the assigned contested calls with margin 0.000 are: gorilla_KB3781 7 of 140, chimp_PTR 2 of 46.

**Fix.**
- **Read-star certificate (and the maximin pick) and `fixed_only_certificate`:** the pair LLR is computed from integer
  counts, LLR(a, b) = (k - m) x lr + (kj - mj) x lrj. Here k / m are the decisive columns matching a / b, and kj / mj
  the junction terms for and against a. It is exactly 0 iff the evidence is balanced.
- **Legacy `read_copy_evidence`:** the decision uses each copy's terms summed in a canonical order (sorted), so copies
  with the same multiset of terms get bit-equal totals. The printed numbers keep the old sums.

No statistical rule changes: alpha, e/3, lr, lrj and the tie keys are untouched.

**Bars.**
- **F1 (unit, exactness).** An exactly balanced read is not Assigned and has margin exactly 0.0, in every column order.
  This holds for the read-star pair LLR helper, `fixed_only_certificate` and the legacy `assign_read`. Each test must
  fail on the unfixed code (red first).
- **F2 (real data, gorilla_KB3781).** The full assign stage is run with the before binary (`copy_assign.before`, md5
  15ba26eb, HEAD b84ddf84) and the fixed binary, same command (the pipeline's), and `assignments.tsv` is compared
  row-for-row.
  - Every row whose status changes must have had an old margin with |margin| < 5e-4 (printed `0.000`), and must move
    from `assigned` to `ambiguous` or `tied`.
  - Expected: those 7 rows.
  - Rows whose `assigned_copy` changes without a status change are reported. These are tie-key picks for exact twins.
  - Any other difference is a FAIL: stop and explain.
- **F3 (the anchor test).** The human arm-1 batches are re-run with the fixed binary and re-scored. The 3 MCL1623 wrong
  calls must no longer be `assigned`. A1 / A2 are reported, descriptive only (A1 stays UNDERPOWERED).
- **F4.** `cargo test --release --no-fail-fast` passes, except the 3 known `pipeline_driver_dispatch` failures. Golden
  tests that pin bit-level sums are updated only if their decisions are unchanged, and each such update is listed.

Outputs `/mnt/linuxdisk/tmp/o2_tie/`.

## Results (2026-10-13; fix 7ac557e7; outputs `/mnt/linuxdisk/tmp/o2_tie/`)

**F1 PASS.** The three tests were red first: legacy `assign_read` assigned the 5:5 tie, and `fixed_only_certificate`
had margin 1.78e-15. Now green; lib 1,047 / 0.

**F2 (gorilla_KB3781 full assign, before vs after): FAIL as registered, explained.** 60,214 rows; 60,189 identical.
- **4 status changes,** all assigned -> ambiguous, with old margin `0.000`, 8-10 decisive columns, p 7e-11 to 2.5e-13,
  two candidates. These are the fix.
- **The registered expectation of 7 was wrong.** The other 3 zero-margin assigned rows are sole-candidate calls
  (n_candidates 1, 0 decisive columns), whose margin is 0 by design. The fix correctly leaves them.
- **3 rows changed only the printed margin,** `-0.000` -> `0.000`. These are exact ties whose old rounding residue
  was negative. This is the same mechanism, with no status or copy change, but the bar did not foresee it, hence the
  FAIL.
- **18 rows changed only the copy:** exact twins now fall through to the tie key, as the prereg allows.
- `families.tsv` and `quant.tsv` change only through those rows.
- Nothing outside the exact-tie mechanism.

**F3 PASS (human anchor test, held-out).** The 3 MCL1623 wrong calls are now `ambiguous`.
- Anchor-trimmed reads: A1 (amended) 1 assigned, 1 right, 0 wrong (UNDERPOWERED). A2 0 / 42.
- Whole human run (129 families, 44,491 rows): assigned 804 -> 675; contested assigned **691 -> 562 (-19%)**.
  - 127 tie-noise calls removed;
  - 185 exact-twin copy changes;
  - 82 margin-sign rows;
  - 5 other status changes:
    - 3 tied -> ambiguous: exact twins; the tie key picks another best copy and both statuses abstain;
    - 2 assigned -> tied (MCL23, real margins 6.9 / 13.8): **a separate, pre-existing race**, explained below.

**F4 PASS.** `cargo test --release --no-fail-fast`: 1,286 passed. The only failures are the 3 known
`pipeline_driver_dispatch` ones.

**Second bug found (not fixed here): the tie-outside gate depends on region order.**
- **The mechanism.** A molecule with an AS-tied placement outside every catalog unit must never be Assigned. The
  registry of such molecules is filled while each region is swept, but read when a family's own region is
  certified. If the region holding the outside placement is swept later, the gate misses and the call stays
  `assigned`. The write-time column `tie_outside_catalog` still says 1.
- **The proof.** The before binary on 1 region thread gives `tied` for both MCL23 rows; on 4 threads it gave
  `assigned`. So "`--region-threads` is byte-identical by construction" is false here.
- **In the shipped runs,** assigned rows with `tie_outside_catalog` = 1: gorilla_KB3781 1 / 624, chimp_PTR 1 / 54, and
  the human anchor run 1 / 675.
- **Fix candidate:** apply the demotion after all sweeps, at write time, where the column is computed.

**Reading.** O2's "LLR > 0" rule was decided by rounding noise on exactly balanced molecules. In human testis these
are common: mosaic molecules matching both copies at several columns each. They made up about a fifth of the
contested assigned calls in these families. Every human O2 number from before 7ac557e7 overstates assignments by that
much.
