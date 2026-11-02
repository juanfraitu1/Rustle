# docs/ — what is where

Reorganised 2026-10-07: the undated living documents, this week's open studies and the current handoff stay here;
closed studies and superseded notes moved, names unchanged, to `archive/2026-09/` and `archive/2026-10/`. Every file,
one line each, is in [`INDEX.md`](INDEX.md) (regenerate with `python3 tools/docs_index.py > docs/INDEX.md`).
Grepping a filename still finds it; only the directory changed.

## Read these first

Start with the `## Current state` block at the top of each ledger and living document (added 2026-11-01).

| document | what it is |
|---|---|
| [`LEDGER_O1.md`](LEDGER_O1.md) · [`LEDGER_O2.md`](LEDGER_O2.md) · [`LEDGER_O3.md`](LEDGER_O3.md) · [`LEDGER_GENERAL.md`](LEDGER_GENERAL.md) | **Start here.** The per-objective ledgers since 2026-10 (pre-registrations, results), each opening with its current state |
| [`NEGATIVE_RESULTS_REGISTER.md`](NEGATIVE_RESULTS_REGISTER.md) | **Every dead end, one row each (rows to 1484 on 2026-11-01). Consult before proposing anything.** Rows cite the study file |
| [`THESIS_OBJECTIVES.md`](THESIS_OBJECTIVES.md) | The three objectives (O1 define, O2 assign-or-abstain on AS-tied multimappers only, O3 detect+flag), scope, what is WON / DEAD / OPEN |
| [`seeded_family_definition.md`](seeded_family_definition.md) | O1: the shipped family definition (§0★, 2026-09-25) and its development record; use its INDEX |
| [`copy_assignment_definition.md`](copy_assignment_definition.md) | O2: assignment and abstention |
| [`O3_STATUS.md`](O3_STATUS.md) · [`STATUS_O2_O3_2026-10-24.md`](STATUS_O2_O3_2026-10-24.md) | O3 evidence summary; the dated 10-24 O2/O3 snapshot |
| [`o1_ledger.md`](o1_ledger.md) | The O1 ledger to 2026-10, section by section (§3–§6z); search a § number |
| [`ADVISOR_QUESTIONS.md`](ADVISOR_QUESTIONS.md) | The advisor's standing questions and the current best answer with register rows |
| [`CHAPTER_O1_O2_COMPOSITION_DRAFT.md`](CHAPTER_O1_O2_COMPOSITION_DRAFT.md) | Chapter draft: O1 units, O2 assignment among exactly them (2026-09-05) |
| [`PREREG_sqanti_single_copy_2026-11-01.md`](PREREG_sqanti_single_copy_2026-11-01.md) | Open study: SQANTI single-copy comparison (register 1482) |
| [`DATA.md`](DATA.md) | Datasets: what they are, where they live (Fritz and laptop), how to rebuild them |
| [`PERFORMANCE.md`](PERFORMANCE.md) | Where pipeline time goes, the speed switches (`RUSTLE_STAR_THREADS`, `RUSTLE_FAMILY_THREADS`, `RUSTLE_TIMING`) with measured effect and whether products stay identical |
| [`REFERENCE.md`](REFERENCE.md) | Glossary, laptop disk mount procedure, the DAZ worked example |
| [`MODULE_STATUS.md`](MODULE_STATUS.md) | Which Rust modules exist and their state (checked by a test in `src/lib.rs`) |
| [`ACTIVE_WORKING_SET.md`](ACTIVE_WORKING_SET.md) · [`MEMORY_DIGEST.md`](MEMORY_DIGEST.md) | Session bookkeeping: the active file set; compacted memory entries |

Moved to `archive/2026-10/merged/`: `HANDOFF_2026-10-04.md`, `PENDING_2026-10-04.md`, `PENDING_2026-09-23.md`, `TERMINOLOGY_FAMILY_SD_DUPLICON_EXPANSION_2026-10-07.md`, `REGISTER_DRAFTS_machine2.md` (and later handoffs).

## Archive

- `archive/2026-09/`, `archive/2026-10/` (and `archive/2026-10/merged/`): closed studies and superseded notes, names unchanged.
- `archive/2026-11/`: **verbatim pre-compaction copies (2026-11-01)** of the ledgers (`LEDGER_GENERAL/O1/O2/O3.md`, `o1_ledger.md`), the register (`NEGATIVE_RESULTS_REGISTER.md`) and the definition and living docs (`seeded_family_definition.md`, `copy_assignment_definition.md`, `THESIS_OBJECTIVES.md`, `ADVISOR_QUESTIONS.md`, `MEMORY_DIGEST.md`, `ACTIVE_WORKING_SET.md`, `O3_STATUS.md`, `STATUS_O2_O3_2026-10-24.md`, `CHAPTER_O1_O2_COMPOSITION_DRAFT.md`, `REFERENCE.md`). The living files in `docs/` are compact; each entry ends with `(full: archive/2026-11/X.md:LINE)` pointing to the full text.

## Open studies (dated 2026-10-05 to 10-07)

(Historical list; those studies are closed and archived. Current open work: the ledgers' `## Current state` blocks.)

Each study is a `PREREG_<name>_<date>.md` (decision rules, committed before looking) and, once run, a
`<NAME>_<date>.md` results file; some have a `_REVIEW_DISPOSITION` file. This week: seed pool on real reads,
entangled baseline, locus units and levels, ideal expression through the default, few-copy ideal cases, default
re-score of NPIP, O2 default roster, gorilla overlap, hierarchy (duplicon boundary vs family; expansions in
families; Yoo 2025 concordance), inversion aligner blindness, O3 candidates Amendment 15. See `INDEX.md`.

## Conventions

- Pre-register before looking (`PREREG_*.md`); hold a substrate back; report pair precision and recall
  (`bench/denovo_eval/pairscore.py`, since 2026-10-29: no F, no bipartite family matching, pseudogene pairs ignored);
  never pool human and gorilla numbers; check the register before proposing.
- When a study closes, move its files to `archive/YYYY-MM/` with `git mv`, regenerate `INDEX.md`, and make sure
  the register row that records the verdict cites the archived path.
- `figures/` holds the publication figures' sources; `superpowers/` holds agent plans and specs.
