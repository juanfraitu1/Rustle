# PREREG — the read-star "best" must be scored the way the verdict is (2026-09-09, before the code change)

## Finding (audit of every `tied`/`ambiguous` molecule, `bench/audit_ambiguity.py`, human MCL0 + gorilla MCL1)
Of 343 human `tied` molecules with ≥ 2 catalog candidates: **188 are absolutely ambiguous** (≥ 2 copies match the
read at every column it covers), **140 are equally imperfect** (best two share the same mismatch count), and
**15 have a unique perfect candidate** with the runner-up 6–9 substitutions worse — yet were `tied`.
Cause for 12 of the 15: `bk` (the best) is `argmax(aligned − X)` over the WHOLE read; a copy that aligns 11 more
bases with 11 substitutions ties a copy that gaps those bases with 0 substitutions, and candidate order breaks
the tie. With the wrong best, a copy identical to it sets `k0` ⟹ `Tied`, while the LLR toward the true best is
−48 … −62 (the `margin` column). The other 3 are identical at every co-observed substitution column and differ
only by an indel — a policy case (the certificate excludes indels by design), not a bug.

## Rule
`bk` = argmax over candidates of the **PSV score** Σ_columns (+1 if the candidate's allele equals the read's
base, −1 if both are known and differ) — the same evidence the pairwise LLR uses — with `aligned − X` and then
index as tie-breaks. Behind `p.best_by_psv` (default ON); escape `--best-by-alignment` restores the old choice.
Full escape `--no-as-tied-only --best-by-alignment` must be byte-identical to the pre-09-09 outputs.

| # | prediction | refuted by |
|---|---|---|
| P1 | the 12 bug cases become `assigned` to their perfect candidate (copy 2 or 11) | < 9 of 12 |
| P2 | the 3 policy cases stay `tied` | any of them assigned |
| P3 | human contested set: assigned rises from 14 by ≥ 9 and by ≤ 40; tied falls by the same; ambiguous unchanged within ± 5 | outside those bands |
| P4 | gorilla: 0 molecules change verdict (its 21 tied-inside audited as 18 identical + 3 equally imperfect) | any change |
| P5 | full escape md5s unchanged: human `91081887`, gorilla `ff0b8f16` | any differs |

---
## Outcome (2026-09-09, after the runs — scored against `ours_fix4`, the fully-§6gz-fixed pre-PSV state, NOT `ours_gate`)
⚠ The first scoring pass diffed against `ours_gate.assignments.tsv`, which predates all three §6gz leak
fixes (contested 69/374/316). Corrected to `ours_fix4` (contested 14/429/316 — the true pre-PSV baseline).

| # | prediction | verdict |
|---|---|---|
| P1 | the 12 bug cases become `assigned` | ⛔ **REFUTED as written.** `assigned` stayed flat at 14 on human, in every intermediate run. The mechanism DID fire correctly — for all 12, `bk` now picks the PSV-perfect candidate (was copy 6, now copy 2), the LLR margin flips from −41…−62 to the SAME MAGNITUDE **positive**, and `n_decisive` goes from the k0-forced 0 to 6–9 real decisive columns. But the corrected certificate finds the family has 19–20 near-identical candidates and not all are rejected at α/(n−1) — the honest verdict is `tied → ambiguous`, a real correctness improvement (the molecule is no longer falsely reported as spanning zero decisive features), not a new `assigned`. |
| P2 | the 3 policy cases stay `tied` | ✓ **held** — margin 0.000, `n_decisive 0` unchanged, indel-only difference, as designed |
| P3 | contested set: assigned +9..+40, tied −same, ambiguous ±5 | ⛔ **REFUTED.** assigned +0, tied −20, ambiguous −14, contested total −34 (759→725). The band was built on the wrong (P1) expectation |
| P4 | gorilla: 0 molecules change verdict | ⛔ **REFUTED as written**, but the SUBSTANTIVE claim holds: **`status` is byte-identical for all 57,646 rows**; 38 rows have a diagnostic-only change (`origin_rejected` flip against the corrected `bk`, same mechanism as human) that never crosses a status boundary |
| P5 | full escape (`--no-as-tied-only --best-by-alignment`) byte-identical | ✓✓ **held on both species** — human `91081887`, gorilla `ff0b8f16` |

### ⭐⭐ A finding the PREREG did not anticipate: the pairwise and origin certificates can disagree under a corrected `bk`
37 human molecules moved OUT of the contested set (`origin_rejected 0→1`) and 3 moved IN (`1→0`), net −34,
fully accounting for the P3/P4 gap. Mechanism: the origin certificate is an ABSOLUTE test (total substitutions
+ indels + unaligned bases over the WHOLE read, against `a.best_copy`), independent of the PAIRWISE test that
`bk`-selection now optimises (PSV-column agreement specifically). A candidate can win the pairwise comparison
(best PSV-column match) while losing the absolute one (more total edits elsewhere in the read) — so correcting
`bk` to be internally consistent with the pairwise evidence can make it INCONSISTENT with the absolute one for
the same molecule. ⚠ Not chased to a root cause beyond this (would need per-column edit-distance instrumentation
beyond what `--dump-star` records); flagged as an open question, not resolved.

### Conclusion
The fix is CORRECT and SHIPS: it repairs a real defect (wrong `bk`, sign-flipped LLR margins, false `k0` on
19–20-candidate families) verified against real data, changes no molecule's status on gorilla, and both
escapes are exactly byte-identical. Its measurable effect on the O2 headline is **not** what was predicted —
it reclassifies some `tied` to a more honest `ambiguous` and, via the certificate interaction above, shrinks
the contested population rather than growing `assigned`. Every number in §6gz/§6gy that quoted human contested
as **14/429/316** should be read alongside the corrected **14/409/302** (total 725, not 759).
