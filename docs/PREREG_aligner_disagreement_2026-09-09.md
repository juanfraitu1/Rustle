# PREREG — `--admit-aligner-disagreement`: widen the gate to aligner self-disagreement (2026-09-09, before any run)

**Finding it follows from** (§6hc addendum 2): 206 of the 216 origin-drop-indels assignments had the TRUE copy
(the primary) scoring BELOW the AS tie because a real insertion costs gap penalties; when the wrong copy's AS is
unique there is no tie and the gate skips the read. Human MCL0: **2,278** molecules whose primary unit ≠ best-AS
unit; 64 % carry ≥ 50 bp of insertion in the primary placement; `--as-tie-ratio 0.98` reaches only 42 %.

**Rule.** With the flag, a molecule is admitted to the certificate if its PRIMARY record's unit differs from
every best-AS record's unit (both by unit-span overlap), in addition to the AS-tie rule. Marked in a gate-only
column `aligner_disagreement`. Default OFF. It is a SCOPE widening of the user's "same AS" definition — the
"coin toss" in a stronger form (the aligner's two stages disagree) — and is measured, not shipped.

| # | prediction | refuted by |
|---|---|---|
| P1 | human: ≈ 2,278 admitted by disagreement (± 10 %; the count here uses the binary's own unit spans) | outside 2,000–2,600 |
| P2 | of those, **≥ 60 % `assigned`**, and **≥ 80 % of the assigned go to their PRIMARY copy** (chaining was right; scoring was fooled by the SV) | < 40 % assigned, or < 60 % to the primary |
| P3 | the assigned ones carry strong evidence: `n_decisive` median ≥ 8, `margin` median ≥ 40 | either below half |
| P4 | the pre-existing 1,143 contested molecules are unchanged (status byte-identical) | any change |
| P5 | gorilla MCL1 / MCL7: the same rule admits far fewer (both families have small SV-driven populations) and creates no new false calls; no directional count predicted | — |
| P6 | escapes untouched: `--no-as-tied-only --best-by-alignment` stays `91081887` / `ff0b8f16` | any md5 differs |
⚠ Runs with `--origin-drop-indels` ON (the two interact by design). Human/gorilla never pooled.

## Outcome (2026-09-09) — runs `bakeoff/human/ours_dis`, `bakeoff/mcl1_dis`, `bakeoff/mcl7_dis`; scorer `bench/o2_disagreement_score.py`
| # | verdict |
|---|---|
| P1 | ✓ **2,278** admitted by disagreement (the binary's unit spans reproduce §6hc's count exactly) |
| P2 | ⛔ **REFUTED: 759/2,278 = 33.3 % assigned** (< 40); of the assigned **477/759 = 62.8 % go to the PRIMARY** — inside the 60–80 band, neither pass nor refuted |
| P3 | ⛔ **REFUTED: margin median 13.8** (< 20); n_decisive median 5 (4–8 band). Against the same run's AS-tied contested assignments: n_decisive 9 / margin 55 — the disagreement population carries **a quarter of the evidence** |
| P4 | ✓ status of all 7,482 base `contested` rows identical (0 changed); over the whole 7,657-row base file, 1 already-`ambiguous`, origin-rejected row reports a different `bk` (16→18) because PSV profiles are built from the admitted read set — status unchanged |
| P5 | gorilla MCL1: 64 admitted → **12 assigned (7 primary / 3 best-AS / 2 other)**, n_decisive 12.5, margin 76 — the contested set grows 33→46 and `assigned` 4→16; MCL7: 10 admitted → 0 assigned (all ambiguous, all origin-rejected). ⚠ gorilla `assigned_copy` is the binary's internal unit id; the catalog index is `catalog_copy_idx` (identical in human, not in gorilla) — the scorer uses `catalog_copy_idx` |
| P6 | ✓ human escape `91081887`; gorilla escape `ff0b8f16` |

### What the 759 human assignments are (class by the assigned copy's relation to the read's own placements)
| class | n | n_decisive med | margin med | AS(assigned) − AS(primary) | AS(best) − AS(assigned) | ≥ 50 bp ins |
|---|---|---|---|---|---|---|
| **primary** (chaining was right, scoring fooled — the predicted mechanism) | 477 | 6 | 13.8 | 0 | 52 | 298 |
| best-AS (scoring was right, chaining wrong) | 43 | 6 | 20.7 | +45 | 0 | 19 |
| **other** (neither — a placement 211 AS BELOW the primary and 306 below the best) | **239** | 5 | **6.9** | −211 | 306 | 207 |
The `other` class is the warning: 31 % of the flag's assignments go to a copy both aligner stages ranked far
down, on the weakest margins in the file (92 of them are copy 20 → copy 21, two 0-read catalog copies 80 kb
apart). The 1,490 `ambiguous`: 660 origin-rejected even with `--origin-drop-indels` (526 carry ≥ 50 bp
insertion), 830 origin-pass with margin median **0** (a PSV dead heat between two copies; `min_p` ≤ 1e-3 in
all 830 — the read is certainly IN the family, not certainly in one copy).

### Reading
The mechanism predicted from §6hc's 206 (true copy under-scored by an SV gap penalty) is real for at most
477/2,278 = 21 % of the blind spot, and even there the evidence is a quarter of the AS-tied contested
assignments'. Widening the gate to aligner self-disagreement buys 759 human assignments at margin 13.8 with a
31 % `other` class, versus 230 at margin 55 without it. **Not recommended as a default**; the flag stays
(default OFF) as the measured answer to "why not admit the primary/best disagreement too?".
