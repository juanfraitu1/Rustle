# Agent handoff — Rustle gorilla multi-copy gene-family thesis

You are picking up a bioinformatics research codebase mid-thesis. Read this first; it tells you what
the project is, where the data lives, the conventions, and the ground rules. **Work on the branch
`collaborator` (see "Branch discipline" below) so your changes never undo the committed progress.**

> Also read **`AGENTS.md`** (general git-safety notes for this repo). Note: its "safe checkpoint"
> SHA/branch info is **stale** — the live thesis work is on `vg/flow-capacity-apportionment` (and now
> your `collaborator` branch), far ahead of the SHA it names. Its *safety conventions* still hold; its
> *specific pointers* are superseded by this document.

---

## 0. Branch discipline (do this first)

All committed progress is on **`vg/flow-capacity-apportionment`**. A branch **`collaborator`** has been
created from that HEAD for you.

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
git checkout collaborator      # <-- work here, NOT on vg/flow-capacity-apportionment
```

- Do **all** your work, commits, and experiments on `collaborator`.
- Do **not** commit to `vg/flow-capacity-apportionment` or `main` — that's how progress gets undone.
- Do **not** `git push`, force-push, rebase, or delete branches unless the human explicitly asks.
- Commit often on `collaborator`. Small, described commits.

---

## 1. What this project actually is

Rustle **began** as a Rust port of StringTie (a transcript assembler), but the **thesis has moved on**.
The thesis is **NOT an assembler**. It is about **defining multi-copy gene families at the RNA level and
assigning reads to copies under mapping ambiguity**, in gorilla (GGO) HiFi/IsoSeq long-read transcriptomes.
Advisor: **Stefan Canzar** — he values clean combinatorial structure, provable theorems, and **no
arbitrary thresholds / no `1/k` heuristics**.

Four objectives:
- **O1 — family definition.** Which loci are paralogous copies of one family (homology graph).
- **O2 — copy assignment.** Assign multi-mapping (MAPQ-0) reads to specific copies via PSVs/SUNs, with a
  significance gate that **abstains** rather than guessing when copies are unresolvable.
- **O3 — allele-specific junctions (ASJ).**
- **O4 — reference-absent copies** (copies missing from the primary assembly).

⚠ **We are moving AWAY from the StringTie-like engine, but NOT removing it yet.** The Rust transcript-
assembly engine (`src/rustle/`, ~65 `.rs` files) is legacy/parity-era. Leave it in place — do not delete
or gut it — but new thesis work happens in the family/copy-assignment layer, mostly under `bench/`
(Python) and `src/rustle/vg_family/` (Rust). If a task seems to call for extending the assembler, check
with the human first; it's probably the wrong layer.

---

## 2. Where everything lives

**Repo:** `/mnt/c/Users/jfris/Desktop/Rustle` (git; remote `origin` = github.com:juanfraitu1/Rustle).

**Python / tools** (a mamba/conda env under miniforge3 — always use these absolute paths):
- Python: `/home/juanfra/miniforge3/bin/python`
- minimap2: `/home/juanfra/miniforge3/bin/minimap2`
- vg: `/home/juanfra/miniforge3/bin/vg`  ·  GraphAligner: `/home/juanfra/miniforge3/bin/GraphAligner`
- Install new deps with **mamba** (`/home/juanfra/miniforge3/bin/mamba install ...`), not pip-into-system.

**Data (big files) live in `/home/juanfra/winloci_scratch/`, NOT `/mnt/c`** (see disk note below):
| file | what it is |
|---|---|
| `GGO.fasta` (3.4 G) | the reference genome = **GCF_029281585.2** (the RNA reference the reads are aligned to) |
| `GGO_mm.bam` (11 G) | the aligned IsoSeq reads. Alignment cmd: `minimap2 -ax splice:hq -uf --eqx -Y -N50 -p0.1 --secondary=yes` (the `-N50 -p0.1 --secondary=yes` forces secondaries for copy-assignment; note it depresses MAPQ — use the AS-margin for a clean read-confusability signal) |
| `mGorGor1.mat.cur.20231122.fasta.gz` (926 M) + `mGorGor1.pat.cur...pat...gz` (880 M) | the **phased diploid genome** (maternal + paternal haplotypes) of the same individual. DNA is for **validation/calibration only** — the thesis is RNA-primary; DNA never gates a family/copy decision. |
| `denovo_skeletons.tsv` (37 M) | per-de-novo-locus exon/intron structure |
| `validated_families.tsv` (72 K) | per-copy loci of the validated families |
| `/mnt/c/Users/jfris/Desktop/final.bed` (103 M) | **SEDEF segmental-duplication map** (WGAC ≥90% id); DNA duplication truth |

**Disk situation (important):** `/mnt/c` is a 477 GB Windows drive, **~80% full (~96 GB free)**. Put all
large intermediates in `/home/juanfra/winloci_scratch/` (a larger Linux disk), **never** in `/mnt/c`.
Small analysis outputs (TSV/JSON/PNG/MD, a few MB) go in `bench/` in the repo.

---

## 3. The shipped state (what "the definition" is right now)

**Family definition = `bench/family_rna_refine.py`.** It builds the multi-copy family catalog from the
RNA read-conflict/homology graph with **three composable, default-on VG-native over-merge gates**, each
with a `--no-*` opt-out, each byte-identical when disabled:

1. **repeat-hub gate** (`--no-repeat-gate`): cuts single extreme-multiplicity repeat/Alu-bridge edges.
2. **recombinant-split gate** (`--no-split-recombinants`): splits mosaic/exon-shuffled bridge families.
3. **multi-repeat-bridge gate** (`--no-repeat-bridge-gate`): family-level "repeat-bridged AND no full
   shared exon" conjunction (the dominant FP class).

Regenerate + verify:
```bash
PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/family_rna_refine.py
# default catalog md5 must be dca64cbd... (573 families). Golden is in test_family_rna_refine.py.
PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/test_family_rna_refine.py   # full self-check
PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/family_level_pr_current.py  # P/R metrics
```
Current metrics vs the diploid-CN gold oracle: **R_oracle 50/57 = 0.877, P_oracle(dedup) 0.917,
E_p protein-purity 0.892, distinct-FP 4.**

**Copy assignment (O2):** the IsoCon-style significance gate (`min_p` certificate; per-read Poisson-
binomial p-value; abstains on unresolvable reads) + the **recombinant-read abstain gate**
(`bench/recombinant_abstain.py`, default-on) that refuses reads carrying two copies' SUNs. Shipped
Rust path: `src/rustle/vg_family/copy_assign.rs`. VG visualization: `bench/o2_vg_visualization.py`.

**Two-regime finding (read this — it frames the whole scope):** `bench/IDENTITY_REGIME_MAP.md`. The
family DEFINITION (O1) is a whole-copy-identity homology problem, clean; the copy ASSIGNMENT (O2) is a
per-segment read-confusability problem (~0.99 identity). They're orthogonal. The apparent near-identical
"FP rate" is mostly the protein truth (E_p) over-splitting real duplications, not RNA over-merges.

**Per-family validation showcase:** `bench/KNOWN_FAMILY_SHOWCASE.md` (RABL2, APOBEC3, MAGEA, ANKRD18,
RGPD8, ZNF92, GSTM, HERC2 — P/R + graphs).

There are **~54 `bench/*.md`** design/finding docs and **~335 `bench/*.py`** analyses. When starting a
task, grep the `bench/*.md` for the topic first — most questions have already been investigated and the
finding (often a *negative* result) is documented. Don't re-run a settled experiment.

---

## 4. Conventions & gotchas (these bite)

- **Determinism:** run everything with `PYTHONHASHSEED=0`. Catalog builds and tests assert byte-identity.
- **Big data → `winloci_scratch`, never `/mnt/c`** (disk is ~80% full).
- **WSL filesystem lag:** `/mnt/c` (drvfs/9p) can lag file existence/content after a child process writes.
  Tests poll (`_settle`) to handle this. If a freshly-written file "isn't there," retry/poll, don't assume failure.
- **NEVER `pkill -f rustle`** — the cwd is `.../Rustle`, so it matches and kills the session itself. Also
  beware `pkill -f <pattern>` where the pattern appears in your own command line (it self-matches). Kill by PID.
- **Concurrent writers:** don't run two things that both write `bench/family_rna_refine.tsv` at once — it
  can leave a mixed-md5 transient on disk (this happened; it caused a stale-artifact bug).
- **Commit footer:** end commit messages with
  `Co-Authored-By: <your model> <noreply@anthropic.com>` (adjust to your agent).
- **Slides:** `three_slides.pptx` may be OPEN/LOCKED in PowerPoint — write to `family_slides.pptx` /
  `family_targeted_slides.pptx` instead. Don't overwrite a file that's open.
- **DNA is validation-only.** The thesis is RNA-primary; a family/copy decision must not be gated by DNA
  (SEDEF, mGorGor1 CN). DNA calibrates and validates; it never decides. Keep this firewall.
- **Honesty discipline:** this project keeps negative results and states caveats plainly. When something
  doesn't work, say so and commit the finding; don't spin a wash as a win.
- **⚠ `.gitignore` is a WHITELIST** (`*` then `!/bench/`, `!/src/`, `!/docs/`, … + individual root
  files). New files **inside** whitelisted dirs (`bench/`, `src/`, `docs/`, `scripts/`, `tests/`,
  `tools/`) are tracked normally. But a new **root-level file, or a file in a non-whitelisted dir, is
  silently ignored** — `git add newfile` does nothing and no warning is obvious. Fix by adding a
  `!/path` line to `.gitignore`, or `git add -f newfile`. (This is why `AGENT_HANDOFF.md` needed a
  whitelist entry.) **Prefer creating new work under `bench/`.**

---

## 5. Good first moves

1. `git checkout collaborator`
2. `PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/test_family_rna_refine.py` — confirm green.
3. Read `bench/IDENTITY_REGIME_MAP.md` and `bench/KNOWN_FAMILY_SHOWCASE.md` for the current framing.
4. For any task, `grep -l <topic> bench/*.md` before coding — the finding may already exist.
