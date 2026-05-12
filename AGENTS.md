# Notes for AI agents (and human collaborators)

If you are an AI agent newly opening this repo, **read this file first**. It tells you what is safe to change, what is not recoverable if you break it, and where to find the rest of the project context.

## 1. This repo is git-tracked. Use git as your safety net.

- Remote: `git@github.com:juanfraitu1/Rustle.git`
- Active branches: `main` and `master` — both point at the same commit by convention.
- The current safe checkpoint is `8877923` (or `HEAD` if you have not made changes yet — check with `git log -1 --oneline`).
- Run `git status` before any non-trivial edit. Commit or stash before risky changes so you can undo with `git reset --hard <SHA>`.
- **Do not force-push.** The remote at `8877923` (or whichever SHA was on `origin/main` when you arrived) is the single recovery anchor. Overwriting it is the one thing that can lose work that has actually been saved.
- **Never run `git reset --hard`, `git clean -fd`, `git checkout .`, or similar destructive commands without first checking `git status` to see what would be discarded.** If there is uncommitted work, commit or stash it first.

## 2. What is NOT recoverable

These are outside the git tree and have no automatic backup. **Do not delete them without confirming with the user.**

| Location                                                        | What it is                                          |
|-----------------------------------------------------------------|-----------------------------------------------------|
| `bench/` (most subdirs — see `.gitignore`)                      | Multi-GB benchmark artifacts; some are slow to regenerate (24 min+ per HMM-EM run) |
| `/tmp/jct_bench/`                                               | Regional BAMs + truth GTFs for divergent-paralog testing |
| `/tmp/exonlen_bench/`                                           | Partial bench runs from the previous session         |
| `~/.claude/projects/-scratch-jxi21-Assembler/memory/`           | Persistent memory notes (if you are Claude Code; other agents won't see this) |
| `/storage/home/jxi21/scratch/Assembler/stringtie/`              | StringTie source we are instrumenting in parallel    |
| `/scratch/jxi21/Assembler/GGO*.bam`, `GGO.fasta`                | Input data — large, irreplaceable from this machine  |

## 3. Start by reading these

In order:

1. **`docs/CONTINUATION_2026_05_11.md`** — full handoff note: workspace layout, build/run commands, parity instrumentation state, what is still missing for single-locus parity and for the VG approach, interrupted work to resume.
2. **`README.md`** — user-facing project overview, with theory references into `docs/ALGORITHMS.md`.
3. **`docs/ALGORITHMS.md`** — derivations (network flow, EM, HMM scoring, variation graphs, SNP rule).

If you are Claude Code specifically, **also** read `MEMORY.md` at `~/.claude/projects/-scratch-jxi21-Assembler/memory/MEMORY.md` for the persistent project-specific notes accumulated over many sessions.

## 4. Build / run / verify

```bash
# Build (~2 min cold, far less incremental)
cd /scratch/jxi21/Assembler/Rustle && cargo build --release

# Primary regression: chr19 benchmark (current ceiling = 1710 matching, 100% Pr)
cd /scratch/jxi21/Assembler && \
  ./Rustle/target/release/rustle -L -p 8 \
    --vg --vg-solver em-hmm --genome-fasta GGO.fasta \
    -o /tmp/rs.gtf GGO_19.bam && \
  gffcompare -r Rustle/bench/parity_gap_analysis/GGO_19_stringtie.gtf -o /tmp/cmp /tmp/rs.gtf && \
  grep -E "Matching transcripts|Transcript level" /tmp/cmp.stats
```

If those numbers drop below ~1710 matching / 100% precision unexpectedly after a code change, you have introduced a regression. Investigate before committing.

## 5. Style notes for working with this user

- The project advisor is skeptical of the pipeline; the user has asked us to **avoid combative or dismissive language toward StringTie or other tools** in user-facing artifacts (README, slides, papers). StringTie is the reference implementation we port; treat it as the trusted baseline, not a competitor.
- The user has explicitly flagged that they lose focus on "word salad concepts" — for user-facing material, prefer **one cohesive arc, one example, one metric** over parallel concepts.
- Prefer terse responses, no trailing summaries the user can read from the diff.
- Only commit and push when the user explicitly asks.

## 6. If something goes wrong

1. Stop and tell the user what happened. Do not try to "fix" by piling on more changes.
2. `git status` and `git log -5 --oneline` so they can see the state.
3. Offer `git reset --hard <SHA>` (with the SHA from §1) as the simplest recovery path.

The remote on GitHub is the single source of truth. As long as it is not force-pushed over, local work can always be restored.
