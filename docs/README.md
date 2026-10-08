# docs/ — what is where

Reorganised 2026-10-07: the undated living documents, this week's open studies and the current handoff stay here;
closed studies and superseded notes moved, names unchanged, to `archive/2026-09/` and `archive/2026-10/`. Every file,
one line each, is in [`INDEX.md`](INDEX.md) (regenerate with `python3 tools/docs_index.py > docs/INDEX.md`).
Grepping a filename still finds it; only the directory changed.

## Read these first

| document | what it is |
|---|---|
| [`HANDOFF_2026-10-04.md`](HANDOFF_2026-10-04.md) | **Start here.** Everything done 10-03/04 and what is left, with paths |
| [`PENDING_2026-10-04.md`](PENDING_2026-10-04.md) · [`PENDING_2026-09-23.md`](PENDING_2026-09-23.md) | Parked work: O3 Amendment 15 run, row-1220 run, HG002, chapter drafts; item 3 of 09-23 (O3 DNA step) still open |
| [`THESIS_OBJECTIVES.md`](THESIS_OBJECTIVES.md) | The three objectives (O1 define, O2 assign-or-abstain, O3 detect+flag), scope, what is WON / DEAD / OPEN |
| [`NEGATIVE_RESULTS_REGISTER.md`](NEGATIVE_RESULTS_REGISTER.md) | **Every dead end, one row each. Consult before proposing anything.** Rows cite the study file |
| [`seeded_family_definition.md`](seeded_family_definition.md) | O1: the shipped family definition (§0★, 2026-09-25) and its full development record; 248 KB, use its INDEX |
| [`copy_assignment_definition.md`](copy_assignment_definition.md) | O2: assignment and abstention; read its "State as of 2026-10-07" preamble first |
| [`O3_STATUS.md`](O3_STATUS.md) | O3: where the reference-absent-copy work stands |
| [`o1_ledger.md`](o1_ledger.md) | The O1 ledger, section by section (§3–§6z); 1.8 MB, search a § number |
| [`ADVISOR_QUESTIONS.md`](ADVISOR_QUESTIONS.md) | The advisor's standing questions and what we concede; status table at the top |
| [`TERMINOLOGY_FAMILY_SD_DUPLICON_EXPANSION_2026-10-07.md`](TERMINOLOGY_FAMILY_SD_DUPLICON_EXPANSION_2026-10-07.md) | Family vs SD / duplicon / expansion: one object, four levels |
| [`DATA.md`](DATA.md) | Datasets: what they are, where they live, how to rebuild them |
| [`REFERENCE.md`](REFERENCE.md) | Glossary, disk mount procedure, the DAZ worked example |
| [`MODULE_STATUS.md`](MODULE_STATUS.md) | Which Rust modules exist and their state (checked by a test in `src/lib.rs`) |
| [`ACTIVE_WORKING_SET.md`](ACTIVE_WORKING_SET.md) · [`MEMORY_DIGEST.md`](MEMORY_DIGEST.md) · [`REGISTER_DRAFTS_machine2.md`](REGISTER_DRAFTS_machine2.md) | Session bookkeeping: the active file set, compacted memory entries verbatim, register rows drafted on machine 2 |

## Open studies (dated 2026-10-05 to 10-07)

Each study is a `PREREG_<name>_<date>.md` (decision rules, committed before looking) and, once run, a
`<NAME>_<date>.md` results file; some have a `_REVIEW_DISPOSITION` file. This week: seed pool on real reads,
entangled baseline, locus units and levels, ideal expression through the default, few-copy ideal cases, default
re-score of NPIP, O2 default roster, gorilla overlap, hierarchy (duplicon boundary vs family; expansions in
families; Yoo 2025 concordance), inversion aligner blindness, O3 candidates Amendment 15. See `INDEX.md`.

## Conventions

- Pre-register before looking (`PREREG_*.md`); hold a substrate back; report sensitivity, precision and bipartite
  matching; never pool human and gorilla numbers; check the register before proposing.
- When a study closes, move its files to `archive/YYYY-MM/` with `git mv`, regenerate `INDEX.md`, and make sure
  the register row that records the verdict cites the archived path.
- `figures/` holds the publication figures' sources; `superpowers/` holds agent plans and specs.
