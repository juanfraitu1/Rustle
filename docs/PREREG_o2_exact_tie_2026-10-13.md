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

## Results

(pending)
