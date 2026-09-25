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
3. `tools/rustle_pipeline.sh` — the whole pipeline, one command per stage (`assemble|families|catalog|assign|flag|all`),
   with a default-on intermediate cache in `PREFIX.cache/` and `--inspect` for analyst dumps.
4. `docs/NEGATIVE_RESULTS_REGISTER.md` — every refuted idea; check it before proposing anything.
5. `docs/MODULE_STATUS.md` — what each Rust module is (shipped, opt-in, other binary); a test keeps it in sync.
6. `bench/README.md` — the analysis scripts and the old-name → new-command table.

## 4. Build / verify
```bash
CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo build --release --bins > build.log 2>&1
CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo test --release > test.log 2>&1
```
Always `--release`; send cargo output to a file (a pipe hides the exit code). Behaviour-preserving changes are proven
by `cmp` of the products against a build of the previous commit on the same inputs.

## 5. Machine rules (WSL2, 5 cores, crashes under load)
One heavy process at a time, in the foreground; big outputs and `TMPDIR` under `/mnt/linuxdisk`; never `pkill -f`
(kill by PID). Pre-register a decision rule (`docs/PREREG_*.md`) before looking at any result it decides.

## 6. Style
- Treat StringTie and other tools as trusted baselines in user-facing text, never as competitors.
- The user prefers one cohesive arc, one example, one metric over parallel concepts; terse answers.
