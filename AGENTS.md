# Notes for AI agents (and human collaborators)

Read this first. It says what is safe to change, what is not recoverable, and where the project context lives.

## 1. Git is the safety net
- Remote `git@github.com:juanfraitu1/Rustle.git`; the only branch is `main` (the former working branch
  `dna-from-genome` was fast-forwarded into it on 2026-09-24; retired branch tips are kept as `archive/*` tags).
- `git status` before any non-trivial edit. Never `git reset --hard`, `git clean -fd` or `git checkout .` without
  first seeing what would be discarded. Never force-push.
- Commit or push only when the user explicitly asks.
- Removed files stay recoverable through the `notebook-YYYY-MM-DD*` tags (`git show notebook-2026-09-24:<path>`) and
  the attic at `~/Desktop/Rustle_attic/<date>/` (one `MANIFEST.tsv` per wave).

## 2. What is not recoverable
Input BAMs/FASTAs and bench outputs under `/mnt/linuxdisk` (and `~/_from_wsl/winloci_scratch`) have no backup. Do not
delete them without the user's confirmation.

## 3. Start here
1. `README.md` — the three stages (family definition, copy assignment, missing copies) and the binaries.
2. `REPRODUCE.md` — every reported number with its exact command.
3. `tools/rustle_pipeline.sh` — the whole pipeline, one command per stage (`assemble|families|candidates|catalog|assign|flag|all`;
   `assign` reads the families' copy table, the legacy catalog only with `--legacy-catalog`; `candidates` is opt-in,
   `--candidates`, since its first acceptance failed (`docs/O3_CANDIDATES_ACCEPTANCE_2026-10-02.md`; the re-run passed: A13, `docs/O3_CANDIDATES_ACCEPTANCE_A13_2026-10-03.md`) and Amendment 14's no-deletion
   control failed: `docs/O3_CANDIDATES_CONTROL_A14_2026-10-03.md`), with a default-on
   intermediate cache in `PREFIX.cache/` and `--inspect` for analyst dumps; `merged` is the genome-wide entry point (assemble + families + assign on the families copy table, byte-identical, resumable; `bench/MERGED_PIPELINE.md`).
4. `docs/NEGATIVE_RESULTS_REGISTER.md` — every refuted idea; check it before proposing anything.
5. `docs/MODULE_STATUS.md` — what each Rust module is (shipped, opt-in, other binary); a test keeps it in sync.
6. `bench/README.md` — the analysis scripts and the old-name → new-command table.

## 4. Build / verify
```bash
CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo build --release --bins > build.log 2>&1
CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo test --release > test.log 2>&1
```
Always `--release`; send cargo output to a file (a pipe hides the exit code). Behaviour-preserving changes are proven
by `cmp` of the products against a build of the previous commit on the same inputs; `tools/identity_check.sh golden|check`
(scripted version: full suite + an end-to-end MCL cmp; `FULL=1` adds a real-data slice) is the harness used for
the 2026-10-05 consolidation and stays the one command to run after any refactor.

**Profiles for iterating (debug mode).** The tree is debug-insensitive (no `debug_assertions` blocks, no timing
tests). Plain `dev` (opt 0) is 10-30x too slow for the alignment-heavy runs; use the pre-defined profiles:
- `cargo test --profile dev-opt` — opt 2, no LTO: the full suite in minutes. THE iterate-on-a-decision loop:
  edit → `cargo test --profile dev-opt` → real numbers still come from `--release` (unchanged rule).
- `cargo build --profile dev-opt --bins` — a run/debug build tolerable on real slices (roughly 1.5-2x release).
- `cargo check --profile quick` — fastest compile for type-check iterations (opt 0, 256 CGUs).

## 5. Machine rules (WSL2, 5 cores, crashes under load)
One heavy process at a time, in the foreground; big outputs and `TMPDIR` under `/mnt/linuxdisk`; never `pkill -f`
(kill by PID). Pre-register a decision rule (`docs/PREREG_*.md`) before looking at any result it decides.

## 6. Style
- Treat StringTie and other tools as trusted baselines in user-facing text, never as competitors.
- The user prefers one cohesive arc, one example, one metric over parallel concepts; terse answers.
