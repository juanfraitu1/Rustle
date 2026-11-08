# docs/ — what is where

Compacted 2026-11-01 and merged 2026-11-03: fourteen living documents stay here, each opening with a `## Current state` block (now through register row 1514). Closed studies and superseded notes live in `archive/`, names unchanged; every file, one line each, is in [`INDEX.md`](INDEX.md) (regenerate with `python3 tools/docs_index.py > docs/INDEX.md`). Grepping a filename still finds it; only the directory changed.

## Read these first

| document | what it is |
|---|---|
| [`NEGATIVE_RESULTS_REGISTER.md`](NEGATIVE_RESULTS_REGISTER.md) | **Every dead end and every verdict, one row each (rows to 1514 on 2026-11-03). Consult before proposing anything.** Rows cite the study file |
| [`THESIS_OBJECTIVES.md`](THESIS_OBJECTIVES.md) | The three objectives (O1 define, O2 assign-or-abstain on reads minimap2 could not place, O3 detect+flag), scope, what is WON / DEAD / OPEN |
| [`LEDGER_O1.md`](LEDGER_O1.md) · [`LEDGER_O2.md`](LEDGER_O2.md) · [`LEDGER_O3.md`](LEDGER_O3.md) · [`LEDGER_GENERAL.md`](LEDGER_GENERAL.md) | The per-objective ledgers since 2026-10 (pre-registrations, results), each opening with its current state. The O3 do-not-quote list is in `LEDGER_O3.md`; the O1 ledger to 2026-10 (§3-§6z) is `archive/2026-11/o1_ledger.md` |
| [`seeded_family_definition.md`](seeded_family_definition.md) | O1: the shipped family definition (§0★) and its development record |
| [`copy_assignment_definition.md`](copy_assignment_definition.md) | O2: assignment, abstention, the gate (`--tie-ms`, `--primary-on-family`) and the truth |
| [`ADVISOR_QUESTIONS.md`](ADVISOR_QUESTIONS.md) | The advisor's standing questions and the current best answer with register rows |
| [`DATA.md`](DATA.md) | Datasets and products: what they are, where they live on Fritz and in the workspace, how to rebuild them |
| [`PERFORMANCE.md`](PERFORMANCE.md) | Where pipeline time goes and the speed switches, with measured effect and whether products stay identical |
| [`REFERENCE.md`](REFERENCE.md) | Glossary, the working set (binaries, driver, bench, cleanup audit), the DAZ worked example |
| [`MODULE_STATUS.md`](MODULE_STATUS.md) | Which Rust modules exist and their state (checked by a test in `src/lib.rs`) |

Merged on 2026-11-03 into the documents above (verbatim copies in `archive/2026-11/`): `STATUS_O2_O3_2026-10-24.md` and `O3_STATUS.md` -> `LEDGER_O2/O3.md`; `MEMORY_DIGEST.md` -> `LEDGER_GENERAL.md`; `ACTIVE_WORKING_SET.md` -> `REFERENCE.md`; `CHAPTER_O1_O2_COMPOSITION_DRAFT.md` and the compact `o1_ledger.md` -> archive only.

## Archive

- `archive/2026-09/`, `archive/2026-10/` : closed studies and superseded notes (`archive/2026-10/merged/` removed 2026-11-08, see `archive/2026-10/MERGED_REMOVED.md`).
- `archive/2026-11/`: the verbatim pre-compaction copies (2026-11-01) of the ledgers, the register and the living docs (the compact files above end entries with `(full: archive/2026-11/X.md:LINE)`), plus every pre-registration of 2026-11-01..03 (`PREREG_*.md`, rows 1485-1514, each with its RESULT section).

## Conventions

- Pre-register before looking (`PREREG_<name>_<date>.md` in `docs/`, a register row with the bars, committed before any result); hold a substrate back; report pair precision and recall (`bench/denovo_eval/pairscore.py`: no F, no bipartite matching, pseudogene pairs ignored); never pool human and gorilla numbers; check the register before proposing.
- When a study closes, write its RESULT into the prereg, fill the register row, `git mv` the file to `archive/YYYY-MM/`, repoint the row's citation, regenerate `INDEX.md`, and update the `## Current state` block of the ledger it touches.
- `figures/` holds the publication figures' sources; `superpowers/` holds agent plans and specs.
