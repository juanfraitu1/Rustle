# PREREG — excision control on the origin-drop-indels assignments (2026-09-09, before the runs)

**Question.** Are the 216 human MCL0 molecules newly `assigned` under `--origin-drop-indels` (§6hc) real?
**Instrument.** §6ff's excision control: remove the assigned copy from the catalog, re-run, follow its reads.
A real assignment has no substitute — its reads must abstain. A leak jumps to the next-best copy confidently.
**Arms.** Copy 2 excised (153 reads) and copy 22 excised (42 reads), `--origin-drop-indels` on, everything
else at today's defaults, same BAM/regions.

| # | prediction | refuted by |
|---|---|---|
| P1 | ≥ 90 % of copy 2's 153 reads are NOT `assigned` to any remaining copy (tied / ambiguous / origin-rejected / absent) | < 75 % |
| P2 | ≥ 90 % of copy 22's 42 reads likewise | < 75 % |
| P3 | Of any that DO assign elsewhere, the margin is small (median < 20) — a near-identical sibling absorbing them, not a confident wrong call | median margin ≥ 40 |
| P4 | The rest of the family is unchanged: reads never assigned to the excised copy keep their status (≥ 99 %) | < 97 % |
⚠ Human only; never pooled with gorilla. Escapes untouched (no code change).

---
## Outcome (2026-09-09)
| # | verdict |
|---|---|
| P1 copy 2 | ✓ **149/153 = 97.4 %** not assigned elsewhere (107 ambiguous+origin-rejected, 41 tied, 1 ambiguous) |
| P2 copy 22 | ✓✓ **42/42 = 100 %** — every one origin-rejected: with copy 22 gone, NO remaining copy explains these reads |
| P3 | ⚠ marginally missed: the 4 copy-2 reads that do assign go to **copy 8** (a near-identical sibling, locus identity 0.98) at margin median **20.7** vs the < 20 predicted — the boundary case, not a confident wrong call |
| P4 | ✓ rest of family unchanged 99.60 % / 99.59 % |
⟹ The origin-drop-indels assignments are real: their copy has no substitute. Precision gate passed.
