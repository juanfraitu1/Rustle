# Continuation note — 2026-05-11

Pick up where the previous session left off. This document is for whoever (human or agent) opens this repo next and needs full context without scrolling through chat history.

Memory files at `/storage/home/jxi21/.claude/projects/-scratch-jxi21-Assembler/memory/MEMORY.md` index project-specific facts in greater detail; this note is the navigation map.

## 1. Workspace layout

```
/scratch/jxi21/Assembler/                       (mirror of /storage/.../scratch/Assembler/)
├── GGO.bam, GGO.bam.bai                         Full Gorilla gorilla PacBio IsoSeq alignment
├── GGO_19.bam, GGO_19.bam.bai                   chr19 subset (primary benchmark)
├── GGO.fasta                                    Reference genome
├── GGO_genomic.gff                              NCBI gorilla annotation (NOT used as `-r` truth)
├── Rustle/                                      THIS REPO (cargo project)
│   ├── src/rustle/                              Rust source
│   ├── bench/                                   Benchmark dirs (many gitignored)
│   ├── docs/                                    Algorithms + parity guides + this file
│   └── target/release/rustle                    Compiled binary
└── stringtie/                                   StringTie C++ source — we instrument this in parallel
    ├── stringtie                                  Compiled binary
    ├── rlink.cpp, rlink.h                         Core algorithm — JFINAL/JKILL traces emitted from here
    └── tmerge.cpp, …
```

**StringTie is the reference implementation we port and parity-test against**. It lives at `/storage/home/jxi21/scratch/Assembler/stringtie/`, builds via `make`, and has been patched to emit parity traces on the same JSONL contract as Rustle (see §3).

## 2. Build / run / benchmark

```bash
# Rustle (this repo)
cd /scratch/jxi21/Assembler/Rustle && cargo build --release          # ~2 min cold

# Single chr19 reference benchmark (current ceiling: 1710 matching, 100% Pr)
cd /scratch/jxi21/Assembler && \
  ./Rustle/target/release/rustle -L -p 8 \
    --vg --vg-solver em-hmm --genome-fasta GGO.fasta \
    -o /tmp/rs.gtf GGO_19.bam && \
  gffcompare -r Rustle/bench/parity_gap_analysis/GGO_19_stringtie.gtf -o /tmp/cmp /tmp/rs.gtf

# StringTie equivalent (reference output is precomputed in bench/parity_gap_analysis/)
cd /storage/home/jxi21/scratch/Assembler/stringtie && \
  ./stringtie -L -p 8 -o /tmp/st.gtf /scratch/jxi21/Assembler/GGO_19.bam
```

Most session work runs through `cargo build --release` and uses `target/release/rustle` directly. Tests are `cargo test` but the real signal is the gffcompare metric on chr19 (and increasingly the regional BAMs under `/tmp/jct_bench/`).

## 3. Parity instrumentation — both sides

We are instrumenting **both Rustle and StringTie** with matched JSONL traces so each stage's decisions can be diff'd. The contract:

| Side          | Trace name      | What it records                                            |
|---------------|-----------------|------------------------------------------------------------|
| Rustle        | `JFINAL_TRACE`  | Per-junction KEEP/KILL + reason at the canonical decision point |
| Rustle        | `parity_decisions.jsonl` | Per-stage decisions (longtrim, predcluster, checktrf, etc.) |
| Rustle        | `GNODE_TRACE`, `RUSTLE_SEED_STATS` | Splice-graph nodes, path-extraction outcomes |
| StringTie     | `JFINAL_TRACE`  | Same junction events (we patched `rlink.cpp`)              |
| StringTie     | `parity_decisions.jsonl` | Same per-stage emission (patched)                  |

Diff scripts under `scripts/` (e.g. `compare_junction_parity.py`, `parity_fidelity_report.py`, `compare_live_parity_to_backward.py`) consume both sides and report first-divergence.

**Current parity status on chr19 (commit `6c590a7`, `main` + `master`):**
- Matching: **1710 / 1839** ref transcripts
- Tr / Intron / Locus precision: **100%** (with `RUSTLE_CHIMERA_FILTER_GTF` + `RUSTLE_NOVEL_JX_MAX_USAGE=3` filter chain)
- Sn: 92.9 % (de novo, no -G)

## 4. What is still missing for **parity** with StringTie

Tracked in detail across these memory files (under `~/.claude/projects/-scratch-jxi21-Assembler/memory/`):

- **`project_parity_ceiling_2026_05_08.md`** — analysed remaining 131 ref misses: isoform-diversity dominates; no readthr threshold helps; BACK guard recovers 1; **1708 is the practical ceiling without changing algorithm**, current 1710 is at that ceiling.
- **`project_parity_dive_2026_05_04.md`** — RSTL.135-class j-class transcripts: real bottleneck is `collapse_high_overlap_variants` intentionally disabled in StringTie; parameter sweep on that function is the next move if you want to push past 1710.
- **`project_m_class_dive_2026_05_04.md`** — 25 m-class artifacts come from rustle's flow extracting both spliced + merged paths where ST only extracts spliced. Fix is upstream path-enum or k-class truncation; deferred.
- **`project_st_port_b_2026_05_05.md`** — stage-by-stage decision-diff methodology validated; one WIN landed (orphan_junc_rescue OFF, +0.28pp F1, all 13 paths were FPs). Next stage targets: `flow` and `checktrf_rescue`.

**Concrete TODOs the next session can pick up:**
1. Parameter sweep on `collapse_high_overlap_variants` — the only candidate path past 1710.
2. Run the parity_decisions diff over the `flow` stage on a single high-divergence locus (RSTL.135 still the canonical test case).
3. Direction A intron-assignment "identity bonus" at alpha != 1.0 (currently a no-op at default α=1.0) — see `project_intron_assign_2026_05_05.md`.

## 5. What is still missing for the **VG approach**

VG = variation-graph mode for multi-copy gene families. Five objectives, status:

| # | Objective                                      | Status                       | Where to look next                       |
|---|------------------------------------------------|------------------------------|------------------------------------------|
| 1 | StringTie parity (single-locus)                | ✅ 1710/1710 matching        | §4 above                                  |
| 2 | Multi-mapping resolution (EM, HMM-EM)          | ✅ wired, validated on AMY  | `project_hmm_em_validation_2026_05_05.md` |
| 3 | Novel-copy rescue from unmapped reads          | ✅ scaffolded, validated chr19 | `project_novel_copy_rescue.md`           |
| 4 | EM-based copy assignment for similar paralogs  | ✅ + SNP rule + k-mer Jaccard | `project_snp_assignment.md`, `project_kmer_jaccard_signal.md` |
| 5 | Discrimination of divergent paralogs           | ⚠️ this session's work       | §5.1 below                                |

### 5.1 Divergent-paralog discrimination — current state

The advisor flagged that HMM-EM works for high-identity paralogs (GOLGA6L7-style, >95%) but the per-base emission scoring may collapse for divergent paralogs (GOLGA8, TBC1D3, NBPF — <85% identity, differ mainly in indels/cassette exons).

Three additive signals added in recent sessions, all gated behind flags (defaults preserve byte-identical chr19 behaviour):
- `--vg-junction-bonus λ` — per missed junction (`project_junction_bonus_2026_05_11.md`)
- `--vg-em-uniform-prior` — diagnostic: disables EM iteration prior to test if it's load-bearing (`project_em_uniform_prior_2026_05_11.md`)
- `--vg-exon-len-penalty μ` — per bp of spliced-length mismatch (`project_exon_len_penalty_2026_05_11.md`) — **just added this session, build verified, default 0**

**Validation status:**
- chr19 metric is invariant at default 0 for all three (correct, target case isn't there)
- Regional BAMs `/tmp/jct_bench/{GOLGA6L7,GOLGA8,TBC1D3}.bam` exist for divergent-case testing
- Bench runs at `/tmp/exonlen_bench/` got partially completed before the user paused:
  - GOLGA6L7: both μ=0 and μ=0.01 → 1/3 matching (sample too small, no HMM-EM engaged)
  - TBC1D3 default `--vg-em-max-work 2000`: family skipped (work=10×373=3730 > 2000), no signal
  - TBC1D3 with `--vg-em-max-work 8000`: was killed mid-run (>1h, CPU contention with GOLGA8)
  - GOLGA8 μ=0: 1/6 matching (single iteration, no μ=0.01 sibling to compare)
- **No conclusive divergent-paralog signal yet** — restart the bench runs sequentially (one at a time, no contention) when you have ~2-3h compute. See §6.

### 5.2 VG-side TODOs (not yet started)

1. **Per-position Bayesian SNP prior**: currently `--vg-snp` is a hard rule. Promote to a soft prior fed into HMM emission, so SNP evidence is naturally additive with the HMM signal.
2. **Indel-aware emission**: profile HMM currently models substitutions cleanly but not indels. For divergent paralogs that mostly differ by indels, this is the dominant signal — affine gap penalties or a pair-HMM extension.
3. **Family-graph alignment as scoring**: replace per-path forward scoring with an alignment against the family-VG (sharing nodes where copies agree). Both expensive and architecturally invasive — deferred.
4. **Phase 2 of module reorg**: `pipeline.rs` is 646 KB. Splitting into `pipeline/bam_ingest.rs`, `pipeline/bundle_build.rs`, `pipeline/vg_resolve.rs`, `pipeline/assemble.rs`, `pipeline/emit.rs` is the advisor-visible refactor. Plan in `project_module_reorg_phase1_2026_05_11.md`.
5. **Cross-validation extension to GOLGA8 + TBC1D3**: complete the Compara-style table for these families. `bench/compara_e/` has the methodology and the script can be extended (add ENSG IDs to `HUMAN_ENSG` dict, rerun).

## 6. Active / interrupted work

- **Bench runs at `/tmp/exonlen_bench/`** were killed mid-run at user request. Partial results retained. To resume:
  ```bash
  cd /tmp/exonlen_bench
  bash run_one_bigwork.sh GOLGA8 0.0  /tmp/jct_bench/GOLGA8.bam  /tmp/jct_bench/GOLGA8.truth.gtf
  bash run_one_bigwork.sh GOLGA8 0.01 /tmp/jct_bench/GOLGA8.bam  /tmp/jct_bench/GOLGA8.truth.gtf
  bash run_one_bigwork.sh TBC1D3 0.0  /tmp/jct_bench/TBC1D3.bam  /tmp/jct_bench/TBC1D3.truth.gtf
  bash run_one_bigwork.sh TBC1D3 0.01 /tmp/jct_bench/TBC1D3.bam  /tmp/jct_bench/TBC1D3.truth.gtf
  ```
  **Run them sequentially** — running two at once doubles wall time per run (the machine has only 2 effective cores).
- **Compara validation E** is complete and pushed (`bench/compara_e/`). Two raw FTP downloads (`*.tsv.gz`, ~180 MB total) are gitignored; the README has the curl commands to refetch.
- **Latest commit:** `6c590a7` on both `main` and `master`. Tree is clean of tracked changes. ~110 untracked files in `bench/` are deliberately ignored.

## 7. Advisor context (for tone/communication)

- Advisor is skeptical the pipeline works — wants reproducible external validation. Compara validation E delivers that.
- Advisor disengages on "word salad" — prefer **one cohesive arc per advisor-facing artifact** (single example, single metric). See `feedback_advisor_attention.md`.
- The README on GitHub was previously perceived as "snarky toward StringTie" — softened this session in §"Relation to StringTie". Avoid combative framing.
- Slide decks and meeting prep live in `/scratch/jxi21/Assembler/Rustle/presentation/` (not exhaustively reviewed; check there before starting any presentation work).

## 8. Recent memory files (most relevant for picking up)

- `project_compara_validation_e_2026_05_11.md` — completed Compara validation E
- `project_exon_len_penalty_2026_05_11.md` — this session's main code addition
- `project_junction_bonus_2026_05_11.md`, `project_em_uniform_prior_2026_05_11.md`, `project_em_perf_gates_2026_05_11.md` — earlier-in-session HMM-EM signal additions
- `project_module_reorg_phase1_2026_05_11.md` — file layout change, Phase 2 not yet done
- `project_perf_findings_2026_05_11.md` — chr19 bottlenecks (BAM I/O, duplicate FASTA loads, NOT HMM-EM itself)
- `project_compact_dispatch.md` — `--vg-solver auto` collapsed routing
- `project_meeting_2026_05_04_followups.md` — punch list of meeting-prep deferred items (a few still open)

If any of those mention a file path, function name, or flag that no longer exists, treat the memory as stale and trust the current code over the note.
