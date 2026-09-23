# Cleanup candidates — likely dead / likely superseded files

Generated 2026-09-22 at `dna-from-genome@f0dc0d33` by `tools/audit_cleanup_candidates.py` (re-run it; this file is overwritten). **Read-only: nothing was moved, edited or deleted.** Full per-file table: `docs/cleanup_candidates.tsv` (filter on `class` and `confidence`).

⚠ A mark is a *candidate*, not a verdict. Before deleting anything: (1) grep the path once more, (2) check the `ledger_sections` / `anchor_citers` columns, (3) prefer `git mv` into an archive directory over `rm` for anything tracked, (4) remember `bench/` data can be slow to regenerate (AGENTS.md §2).

## How a file is kept

A file is **anchored** if a path-like token naming it (full path, unique path suffix, unique basename, output prefix, a glob matching exactly 1 file, an enclosing directory of ≤40 files, or a Python `import` resolved against the importing script's directory and its parents) appears in an anchor source: `docs/o1_ledger.md`, `docs/NEGATIVE_RESULTS_REGISTER.md`, `docs/PREREG_*.md`, top-level `docs/*.md`, `docs/experiments/*.md`, `README.md`, the auto-memory directory, Rust `src/`/`tests/`, or `Cargo.toml`. Anchoring then propagates: a file named by an anchored script or anchored markdown write-up is anchored, index sidecars (`.fai`, `.bai`) follow their parent, and a script that names anchored data (its generator — the naming line looks like a write) is anchored; a script that merely reads anchored data is not. An ambiguous bare basename (e.g. `reads.fa`, found in many directories) never anchors on its own. Citations from `docs/archive/`, `docs/superpowers/`, `AGENTS.md` and un-anchored scripts are *weak*: recorded, not counted.

## Rules (applied in this order; stale cutoff = 2026-08-23)

| class | confidence | rule | files | bytes |
|---|---|---|---:|---:|
| **PROTECTED** | - | build/config files, Rust sources and tests (Rust reachability is docs/MODULE_STATUS.md's job, enforced by module_status_tests), and the anchor docs themselves. Never a candidate. | 185 | 7M |
| **TEMP** | high | Python/pytest caches; untracked or git-ignored files at the repo root; untracked *.log / *err* / *out* / *.patch.txt / checkpoint files that nothing cites. | 1 | 0B |
| **REFUTED-MODULE** | medium | Rust module whose `//! **STATUS:**` header is REFUTED and that no other file names (a REFUTED module that is still imported, e.g. collapse_gate.rs, stays PROTECTED). | 0 | 0B |
| **SUPERSEDED-PORTED** | medium | Python script that a Rust source line declares it ports ('Port of', 'Faithful Rust port of', 'Mirrors', 'migration'); 'low' when only a function or part is ported (`x.py::f`, `x.py loaders`). The Python may still serve as a parity oracle or golden-fixture generator -- check tests before deleting. | 0 | 0B |
| **SUPERSEDED-VERSION** | medium | Older member of a version series in the same directory (_v1.._vN, foo/foo2/foo3, dated _YYYY-MM-DD copies, foo vs foo_fix/_final/_new) that nothing anchors. 'high' when the newest member IS anchored. | 0 | 0B |
| **SUPERSEDED-CITED** | low | Older member of a version series that IS anchored: provenance for a recorded result -- archive, don't delete. | 0 | 0B |
| **LEGACY-ASSEMBLER** | medium | Not anchored, and its path or first 200 lines name StringTie-era assembler machinery (bundle/transfrag/parity/gffcompare/...) -- the assembler layer was retired (docs/RETIREMENT_AND_MIGRATION.md). | 0 | 0B |
| **AMBIGUOUS-CITE** | low | Not anchored; an anchor source names it only by a bare basename shared by several files, a directory too large (>40 files) or a glob matching several files -- may or may not mean this copy. | 13 | 189K |
| **ORPHAN** | medium | Not anchored; named only by files that are themselves not anchored (e.g. a figure named only by its un-cited plotting script). | 0 | 0B |
| **UNCITED-STALE** | medium | Not named by anything (wide globs / big directories / scripts' ambiguous basenames ignored), last touched more than 30 days ago. | 4 | 50K |
| **PROBABLE-PROVENANCE** | - | Not cited by name, but it sits in an experiment directory (below bench/, docs/, ...) holding anchored files, or under a folder whose anchored README/write-up covers it -- usually an output of that experiment written under a computed name and cited as a folder. Verification judged 5/6 such files provenance: not a candidate. | 3 | 55K |
| **UNCITED-RECENT** | low | Not cited, touched within 30 days: may be work in progress. | 16 | 1M |
| **KEEP-CITED** | - | Anchored directly or transitively. Not a candidate. | 468 | 22M |

**34 candidates** of 690 files (1 high, 4 medium, 29 low confidence).

## Candidates by directory

| directory | TEMP | REFUTED-MODULE | SUPERSEDED-PORTED | SUPERSEDED-VERSION | SUPERSEDED-CITED | LEGACY-ASSEMBLER | AMBIGUOUS-CITE | ORPHAN | UNCITED-STALE | UNCITED-RECENT | kept |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `docs/` |  |  |  |  |  |  | 3 |  | 1 | 16 | 163 |
| `bench/` |  |  |  |  |  |  | 9 |  |  |  | 239 |
| `(root)` |  |  |  |  |  |  | 1 |  | 1 |  | 6 |
| `test_data/vg_hmm/` |  |  |  |  |  |  |  |  | 2 |  | 0 |
| `bench/__pycache__/` | 1 |  |  |  |  |  |  |  |  |  | 0 |

## TEMP (1)

- `bench/__pycache__/` [ignored, high] (cache-dir:31-files)

## Overlay: files tied to the dropped ASJ objective (4)

Not classified as dead by this flag alone — ASJ was dropped as an objective (memory, 2026-08-07), but the binaries still build. Decide as a scope question, not a cleanup one.

- `bench/ASJ.md` — AMBIGUOUS-CITE
- `src/bin/asj.rs` — PROTECTED
- `src/bin/asj_verify.rs` — PROTECTED
- `src/rustle/vg_family/allele_specific_junctions.rs` — PROTECTED

## Known blind spots

- Paths built at run time (`f"bench/{name}.tsv"`, shell loops) are invisible to a token scan: such outputs land in UNCITED-IN-ANCHORED-DIR or UNCITED-*, never KEEP.
- A citation proves a file was *named*, not that the naming text is still true; an anchor in a retracted ledger section still anchors (see `ledger_sections`).
- Data outside the repo (`/mnt/linuxdisk`, BAMs) is not audited; `.worktrees/`, `tools/stringtie` (submodule), `.remember/`, `.superpowers/` are excluded.
- Rust module reachability is not re-derived here — see `docs/MODULE_STATUS.md`.
