# PREREG — the held-back family for the O2 refactor: gorilla MCL58 (2026-09-09, before any run)

**Why.** The AS-tied gate, `best_by_psv`, the tie-outside rule and `--origin-drop-indels` (§6gv–§6hf) were
developed and their defaults decided on human MCL0 and gorilla MCL1 / MCL7. Trap 15: a rule validated only
where it was developed is not validated. The nominal hold-out MCL2 is uninformative for O2 (1 AS-tied molecule).

**Selection rule (outcome-blind, `scratchpad/count_astied.py`, BAM only):** among the 32 sweep_v19 families
never used for O2 development, the one with the MOST molecules whose best AS is tied across ≥ 2 of its own
non-dropped units. Result: **`fam_MCL58_073242`** — 9 copies (3 kept / 6 partner), region
NC_073242.2:29933824-35695970, 7,847 molecules in its units, **166 AS-tied**. Runner-up MCL54 (83).

**Runs.** Shipped defaults, no flags (`copy_assign --bam npip3.bam --fasta npip3_contigs.fa --families
copies.tsv --copies-fa copies.fa --regions regions`); then the excision of the copy receiving the most
assignments (catalog without it, `bench/o2_excision_sweep.py` logic, single copy); then
`--no-as-tied-only --best-by-alignment --no-origin-drop-indels` for the record.

| # | prediction (transferred from the development families: human 230/1,143 = 20 %, MCL1 4/33 = 12 %, MCL7 0/11) | refuted by |
|---|---|---|
| P1 | contested set (AS-tied ∧ origin-pass ∧ ≥ 2 catalog candidates) between **20 and 166**; assigned fraction **between 3 % and 35 %** of it | assigned > 50 % (the gate leaks) or contested < 10 |
| P2 | assigned molecules: n_decisive median **≥ 4**, margin median **≥ 14** (≥ 2 decisive columns) | either below half of that |
| P3 | excision of the top copy: **≥ 90 %** of its assigned molecules abstain, 0 move with margin ≥ 40 | < 80 %, or ≥ 2 silent-confident |
| P4 | the origin certificate rejects **< 60 %** of the AS-tied in-catalog molecules (human 429/4,115 ≈ 10 % of gate rows; MCL1 higher) | ≥ 80 % |
| P5 | nothing in the shipped pipeline is changed by this family; if a defect is found it is recorded, not fixed here | — |

## Outcome (2026-09-09) — `bakeoff/mcl58/{base,no2,esc}`
Gate: 276 of 22,105 molecules AS-tied, 110 with a tie outside the catalog; 169 rows, **129 origin-rejected**,
114 with < 2 candidates. **Contested 40 = 2 assigned (5.0 %) / 38 tied / 0 ambiguous.**
| # | verdict |
|---|---|
| P1 | ✓ contested 40 (20–166), assigned 5 % (3–35 %) — no leak |
| P2 | ⚠ missed, not refuted: both assignments rest on **2 decisive columns, margin 13.8** (copy 2, p 1e-6, 2 candidates) |
| P3 | ✓ excise copy 2: **2/2 abstain** (the contested set empties: 0/0/0), 0 silent |
| P4 | ⚠ missed, not refuted: origin-rejected **129/169 = 76 %** of the rows (predicted < 60; refuted at ≥ 80). 6 of the 9 copies are `partner` rows; the certificate rejects most reads against them |
| P5 | ✓ nothing changed |
The 38 tied molecules have **n_decisive median 0** — K = 0 twins: the family's copies are identical over the
read footprints; abstention is the only honest answer and the pipeline gives it. Escape
(`--no-as-tied-only --best-by-alignment --no-origin-drop-indels`) `5b6402ee`: 2 / 35 / 0 at width 1.0.
**Transfer verdict:** the refactored O2 behaves on the held-back family as on the development families —
abstention dominates, the few assignments survive excision, the gate leaks nothing.
