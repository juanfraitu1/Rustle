# PREREG — candidates are ALL tied placements (user's rule, 2026-09-09; written before the code change)

**Rule.** For an AS-tied molecule, the candidate set is the union of (a) this family's catalog copies overlapping any of
its tied placements and (b) an ad-hoc locus for every tied placement that overlaps NO catalog copy of this family.
Read-star scores every candidate. Outcomes: `assigned` (one catalog copy wins), `tied`, `ambiguous`, and the new
**`assigned_outside_catalog`** (an ad-hoc locus wins — never counted as a copy assignment; O3/O1 material).
The L3 sole-candidate label is retired: a tied molecule is never `assigned` on one candidate.

| # | prediction | refuted by |
|---|---|---|
| P1 | Human MCL0: the 4,706 EIF3C/EIF3CL reads stop being `assigned`; they become `assigned_outside_catalog` or `tied` | > 100 of them remain `assigned` |
| P2 | Human contested set is **unchanged**: 759 / 69 assigned / 374 tied / 316 ambiguous (these already had ≥ 2 catalog candidates) | any of the four numbers moves by > 2 % |
| P3 | Gorilla MCL1 output is **unchanged** (0 sole candidates there) | any row changes |
| P4 | Escape `--no-as-tied-only` stays byte-identical (`ff0b8f16` gorilla, `91081887` human) | any md5 differs |
| P5 ⚠ | Between the 4,706, EIF3C and EIF3CL are separable by PSVs for a MAJORITY (they are ~99.9 % identical; ≥ 1 decisive column per read is plausible) | < 50 % separable ⟹ report as tied; that is also a legitimate result |

---
## Outcome (2026-09-09, after the runs)
| # | verdict |
|---|---|
| P1 | ✓ **4,705 of 4,706** EIF3C leaks are `tied` with `tie_outside_catalog = 1`; human `assigned` 4,958 → 15 |
| P2 | ⛔ refuted as written — contested human 69 → **14** assigned, 374 → 429 tied: 55 contested reads had a THIRD tied placement outside the family. The rule applied consistently; the prediction was wrong |
| P3 | ⛔ refuted as written — gorilla 6 → **2** assigned (597 tie-outside molecules); the md5 form of the test was ill-posed since the gated output gained a column |
| P4 | ✓ both escapes byte-identical: human `91081887`, gorilla `ff0b8f16` |
| P5 | ⛔ refuted — O2 run on the EIF3C/EIF3CL family itself: **4,685 of 4,706 `tied`**, 19 ambiguous, 2 assigned (6 decisive columns). Not separable over these reads; the correct outcome is abstention |
`assigned_outside_catalog` (scoring the outside locus as a candidate) is **deferred**: the EIF3C experiment
says its payoff for that class is 2/4,706; the 55 + 4 contested-with-outside reads are untested.
