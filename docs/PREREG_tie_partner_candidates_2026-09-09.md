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
