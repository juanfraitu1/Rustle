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

⚠ **The StringTie-like assembler + network-flow + HMM machinery is DEAD CODE — already marked, do NOT hunt
for it or extend it.** As of 2026-07 the whole legacy transcript-assembly engine (`src/rustle/*.rs` at the
top level: `pipeline.rs` 952 KB, `max_flow.rs`/`global_flow.rs` = the network flow, `transcript_filter.rs`,
`path_extract.rs`, the `*_st.rs` StringTie modules, `bundle*`, `transfrag*`, the `rustle` binary
`src/rustle/main.rs`, `treepat.rs`/any HMM-ish path scoring, `tracing/`, `parity/`, `graph.rs`,
`junction*`, etc.) carries a `//! DEAD CODE — StringTie-era …` banner at the top of each file, and
`src/rustle/lib.rs` tags every `pub mod` line `[KEEP]` / `[DEAD]` / `[SHARED-tendril]`. **Do not read it to
"understand the pipeline", do not fix its warnings, do not extend it — it is slated for removal.** The full
map + removal order is in **`docs/RETIREMENT_AND_MIGRATION.md`** (+ `bench/RETIREMENT_MAP.tsv`, 89 modules
classified). The only reason it still compiles is that the LIVE code borrows a few symbols from 4
"SHARED-tendril" files (`vg.rs`/`graph.rs`/`path_extract.rs`/`bundle.rs`) — those get extracted before the
assembler is deleted; don't build new work on them.

✅ **The LIVE thesis code is `src/rustle/vg_family/` (Rust) + `bench/` (Python analyses).** The whole O1
family-definition pipeline is now **migrated to Rust** (see §3). If a task points you at the assembler,
it's the wrong layer — the answer is in `vg_family/` or `bench/`.

---

## 2. Where everything lives

**Repo:** `/mnt/c/Users/jfris/Desktop/Rustle` (git; remote `origin` = github.com:juanfraitu1/Rustle).

**Python / tools** (a mamba/conda env under miniforge3 — always use these absolute paths):
- Python: `/home/juanfra/miniforge3/bin/python`
- minimap2: `/home/juanfra/miniforge3/bin/minimap2`
- vg: `/home/juanfra/miniforge3/bin/vg`  ·  GraphAligner: `/home/juanfra/miniforge3/bin/GraphAligner`
- Install new deps with **mamba** (`/home/juanfra/miniforge3/bin/mamba install ...`), not pip-into-system.

**Big files live in `/home/juanfra/winloci_scratch/`, NOT `/mnt/c`** (see disk note below). All absolute
paths below are exact.

**Genomes** (`/home/juanfra/winloci_scratch/`):
| file | what it is |
|---|---|
| `GGO.fasta` (3.6 G) + `GGO.fasta.fai` | the reference genome = **GCF_029281585.2** (the RNA reference the reads map to). **Soft-masked** (lowercase = repeats/TEs — the family gates read this). |
| `GGO_genomic.gff` (693 M) | NCBI RefSeq **gene annotation** (coords, strand, biotype). Source of `bench/gene_meta_strand.tsv` (the antisense gate) + the annotation loader. |
| `mGorGor1.mat.cur.20231122.fasta.gz` (971 M) + `mGorGor1.pat.cur.20231122.fasta.gz` (923 M) | the **phased diploid genome** (maternal + paternal) of the same individual. DNA is **validation/calibration ONLY** — the thesis is RNA-primary; DNA never gates a family/copy decision. |

**Read alignments — BAMs** (`/home/juanfra/winloci_scratch/`):
| file | what it is |
|---|---|
| `GGO_mm.bam` (11.7 G) + `GGO_mm.bam.bai` | the aligned IsoSeq reads. Cmd: `minimap2 -ax splice:hq -uf --eqx -Y -N50 -p0.1 --secondary=yes` (`-N50 -p0.1 --secondary=yes` forces secondaries for copy-assignment; it depresses MAPQ, and expresses multimapping as **secondary** alignments not MAPQ-0 — use the AS-margin / the 0x100 flag). `-N50` is the shipped choice (bumping `-N` only enriches KRAB-ZNF repeat-fragmentation — `bench/MINIMAP2_N_SWEEP.md`). |
| `GGO.bam` + `.bai` | **symlink → `GGO_mm.bam`** (same 11.7 G file). Some scripts (`psv_graph_genomewide.BAM`, the O2 read-fetch) reference `GGO.bam`; identical data. |

**De-novo + derived data** (`/home/juanfra/winloci_scratch/`):
| file | what it is |
|---|---|
| `denovo_skeletons.tsv` (38 M) | per-de-novo-locus intron chain (exons = complement) |
| `denovo_transcripts.meta.tsv` (6.6 M) | per-DN-locus metadata (chrom/start/end) — the coordinate source for `dn_exons` / `load_meta` |
| `validated_families.tsv` (71 K) | per-copy loci of the validated families (the `FAM_TSV` input) |
| `o2vg/` (154 × `fam<N>.json`) | cached **materialized O2 variation graphs** (`o2_vg_visualization` / the Rust `materialize_family` reproduce these) |

**Precomputed family-definition inputs — committed in `bench/`** (small; the def READS these, does not recompute them here):
`denovo_family_edges.tsv` (core_recip edges) · `ri_sharedlen_universal.tsv` (aln_frac) · `vg_repeat_catalog.tsv`
(canonical-minimizer multiplicity — repeat/multi-repeat gates) · `a1_read_consensus_o1.tsv` (allele signal, demote) ·
`gene_meta_strand.tsv` (gene strand — antisense gate) · `denovo_families.tsv` (raw single-linkage families).

**DNA duplication truth:** `/mnt/c/Users/jfris/Desktop/final.bed` (107 M) — **SEDEF segmental-duplication map**
(WGAC ≥90% id). Validation-only.

**Disk situation (important):** `/mnt/c` is a 477 GB Windows drive, **~80% full (~96 GB free)**. Put all
large intermediates in `/home/juanfra/winloci_scratch/` (a larger Linux disk), **never** in `/mnt/c`.
Small analysis outputs (TSV/JSON/PNG/MD, a few MB) go in `bench/` in the repo.

---

## 3. The shipped state (what "the definition" is right now)

**Family definition (spec) = `bench/family_rna_refine.py`.** It builds the multi-copy family catalog from
the RNA homology graph (edge oracle `core_recip≥0.19 AND aln_frac≥0.24`) → γ-quasi-clique refine (γ=0.20)
→ **FOUR composable, default-on over-merge gates** (each `--no-*` opt-out, byte-identical when disabled)
→ allele-demote:

1. **repeat-hub gate** (`--no-repeat-gate`): cuts single extreme-multiplicity repeat/Alu-bridge edges.
2. **recombinant-split gate** (`--no-split-recombinants`): splits mosaic/exon-shuffled bridge families.
3. **multi-repeat-bridge gate** (`--no-repeat-bridge-gate`): family-level "repeat-bridged AND no full
   shared exon" conjunction.
4. **antisense/reciprocal-overlap gate** (`--no-antisense-gate`): NEW (2026-07). A *genome-architecture*
   axis (coords+strand) — cuts edges between same-contig, opposite-strand, ≥50%-reciprocally-overlapping
   genes (can't be two copies of one gene), mega-span-guarded. Strand from `bench/gene_meta_strand.tsv`.

Regenerate + verify:
```bash
PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/family_rna_refine.py
# default catalog md5 = 548029ad... (566 families). Golden is in test_family_rna_refine.py.
# --no-antisense-gate reproduces the PRIOR default dca64cbd... (573 fam) byte-identical.
PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/test_family_rna_refine.py   # full self-check
PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/family_level_pr_current.py  # P/R metrics
```
Current metrics vs the diploid-CN gold oracle: **R_oracle 51/57 = 0.895, P_oracle(dedup) 0.917,
E_p protein-purity 0.894, distinct-FP 4.** (Recall 51 not 50 = a genomic-linkage relabel-credit fix for a
captured-but-mislabelled gene; measurement-only, precision byte-identical.)

⭐ **The O1 pipeline is now MIGRATED to Rust** (byte-parity per module, `PYTHONHASHSEED=0` golden fixtures):
`src/rustle/vg_family/{family_definition,family_loaders,edge_oracle,minimizers,repeat_catalog,
bridge_detector,recombinant_split,multi_repeat_bridge}.rs` + the **`family_define` binary**
(`src/bin/family_define.rs`) which reproduces `dca64cbd` **570/573 families byte-identical** (the 3 differ
only in the one NP-hard γ-quasi-clique-splitter blob — a documented, accepted non-parity: networkx-Louvain
witness vs the deterministic Rust splitter). networkx is fully retired from the O1 core. See
`docs/RETIREMENT_AND_MIGRATION.md`. The Python `family_rna_refine.py` remains the canonical spec + is where
the antisense gate currently lives (not yet ported).

**FP boundary is PROVEN, not asserted** (`bench/FP_EXCLUSION_DISCRIMINATORS.md`): the residual domain-share
over-merges (~21 blocks, ~3–4%) are separable ONLY by DNA copy-number — nucleotide / protein-E_p / TE /
VG-topology / (architecture, the antisense corner) all tested adversarially; the dispersed bulk is
DNA-CN-bound because an FP domain-share ≡ a real conserved-domain paralog at the RNA level. Do NOT re-run
those exclusion axes on domain-shares — it's settled.

**Copy assignment (O2) — also migrated to Rust.** The `min_p` significance gate (`copy_assign.rs`) + the
recombinant-read **abstain gate** (`recombinant_abstain.rs`, refuses reads carrying two copies' SUNs) + the
full **VG-materialization** (`o2_materialize.rs`: copy-extraction, minimap2/samtools orchestration,
`o2_columns.rs` = a pysam-pileup byte-parity port, `o2_margin_gate.rs` = the margin gate) — the Rust
`materialize_family` reproduces a fresh Python `o2_vg_visualization.materialize_family` end-to-end (zero
field divergence, no minimap2 drift). The Python `bench/o2_vg_visualization.py` / `recombinant_abstain.py`
remain the specs. **O3 (ASJ):** core is Rust (`allele_specific_junctions.rs` + `asj` bin); the analysis
layer `asj_strand_bias.rs` / `asj_verify.rs` / `asj_genetic_core.rs` (the shipped 54-call SOR-clean core)
are ported byte-parity; `asj_motif_check`/`aggregate`/`evidence` remain Python. **O1 copy-NUMBER:**
`bench/family_copy_number.py` now reports a 3-probe cardinality stack `n_loci → χ(H) → depth_cn` — the
depth leg (`rna_copy_number_depth.py`, expressed-copy estimate) counts the K=0 identical copies χ(H) is
blind to (validated Spearman 0.907 vs DNA CN on collapsed families; a noisy lower bound, silent copies
invisible). See `bench/NONEXPRESSED_MULTIMAPPING.md` + `bench/DNA_VS_RNA_SPECIFICITY.md` for two recent
advisor-question findings.

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
