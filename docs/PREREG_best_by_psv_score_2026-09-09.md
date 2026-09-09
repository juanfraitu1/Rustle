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
