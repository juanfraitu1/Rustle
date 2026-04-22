# Session continuation — 2026-04-21

Capture of the state at end of the KRAB-ZNF S2-falsifier session so the work
can resume cleanly.

## What landed this session

### 1. Zero-width graph node compaction — DEFAULT-ON

- `Graph::compact_zero_width_nodes()` in `src/rustle/graph.rs:581`
- `Graph::compact_transparent_nodes(include_pass_through: bool)` helper
- Wired in `src/rustle/graph_build.rs` `create_graph_with_longtrim` after
  `create_graph_inner`
- Transitive-closure edge rebuild: for every path `u → [zw…] → v` in the old
  graph, add `u → v` in the compacted graph. Safe because zero-width nodes
  have `start == end`, so no read ever maps to them.
- **Impact on GGO_19:** 15,333 → 14,031 nodes (−1,302), 17,324 → 12,813 edges
  (−26%). Final tx 1,926 → 1,911. Matching tx loss 1637 → 1637 (zero valid
  chain loss — FP rate on removed tx is 20%, well below the 40–60% of generic
  filters).
- **Opt-out:** `RUSTLE_COMPACT_ZW_NODES_OFF=1`
- **Trace:** `RUSTLE_COMPACT_ZW_NODES_TRACE=1`

### 2. Chain-pass-through compaction — experimental, DO NOT ENABLE

- `Graph::compact_chain_pass_through()` in `src/rustle/graph.rs:604`
- Wired opt-in only: `RUSTLE_COMPACT_CHAIN=1`
- Merges contiguous pass-through nodes (1 parent + 1 child, non-source/sink,
  not hardstart/hardend/longtrim_cov) by extending head's `end` to absorb
  descendants.
- **Regresses catastrophically:** matching tx 1637 → 729 (−908), transcript-
  level F1 87 → 36 on GGO_19.
- Root cause: pathpat bit length changes for transfrags spanning the merged
  region, which breaks downstream compatibility heuristics. A working version
  would need simultaneous recalibration of the pathpat-based filters — multi-
  session refactor.
- Kept in the tree as a clearly-marked experiment. See
  `memory/project_chain_compaction_negative.md`.

### 3. VG family discovery — strand-mirror filter

- `is_strand_mirror()` in `src/rustle/vg.rs:230`: same chrom, opposite strands,
  ≥90 % coord overlap → not linked as family.
- Removes the 67/76 "families" on NC_073244.2 that were just same-coordinate
  opposite-strand bundle pairs (per-strand bundling artifact, not paralogs).
- Result: 76 real-looking families → **10 real families**, 1126 mirror pairs
  skipped.
- Trace: `RUSTLE_VG_TRACE=1`.

### 4. VG family report — per-gene sub-loci refinement

- `refine_bundle_into_gene_loci()` in `src/rustle/vg.rs`; `write_family_report_
  with_em()` now takes `transcripts: Option<&[Transcript]>` and adds an
  `n_gene_loci` column.
- For each bundle in a family, clusters assembled transcripts by coord
  proximity (gap ≥ `RUSTLE_VG_REPORT_GAP`, default 5 kb). Each cluster becomes
  one "copy" in the report.
- On GGO KRAB-ZNF chromosome: the Family containing ZNF383 paralogs collapses
  from 5 mega-bundle spans (each 120–430 kb wide) to **12 per-gene sub-loci**
  that align 1-to-1 with reference ZNF gene annotations.

### 5. S2 falsifier — LOC115930085 ZNF383 paralog recovery

This is the clean demonstration for the advisor that multi-mapping-aware
assembly is needed for gene family completeness.

- GGO KRAB-ZNF-dense chromosome is **NC_073244.2**, not NC_073243.2 (which
  GGO_19.bam targets — only 18 ZNFs there). NC_073244.2 has 299 annotated
  zinc finger genes.
- LOC115930085 at 48,772,798–48,792,019 ("ZNF383", − strand). 25 reads, **21
  secondary, 4 primary.** Primary-elsewhere reads have their primary at
  LOC101153727 (48,384,759–48,415,507, + strand, also "ZNF383") — the paralog.

| pipeline | tx at LOC115930085 |
|---|---:|
| StringTie `-L` | 0 |
| Rustle default | 0 |
| Rustle `--vg --vg-solver none` | 6 |
| Rustle `--vg` (EM) | 2 (matches 4/5 ref exons, 3 exactly) |

See `docs/experiments/05_loc115930085_s2_falsifier.md` for the full write-up
with figures.

### 6. Figures + docs

- `docs/figures/make_loc115930085_figs.py` — generates all three S2 figures
  from the pipeline outputs in `/tmp/ks_*.gtf` and the family TSV.
- `docs/figures/fig_loc115930085_reads.png` — primary vs secondary alignment
  pattern at the locus.
- `docs/figures/fig_loc115930085_recovery.png` — per-pipeline tx count +
  exon-match with NCBI ref.
- `docs/figures/fig_family1_refinement.png` — before/after family region
  decomposition.
- `docs/experiments/05_loc115930085_s2_falsifier.md` — narrative doc for
  the advisor story.

### 7. Earlier session work consolidated in this commit

The modified-but-uncommitted files also include the accumulated work from the
prior StringTie-parity arc: full-chain witness filter, TF_COALESCE default-on,
end-RI strict filter, micro-node tail trim, alt_tts_end decoupling, and the
StringTie-parity scaffold (`src/rustle/stringtie_parity.rs`). Each has a
corresponding memory file — see the project's `memory/` directory for the
full audit trail.

## Current state

- **Default baseline on GGO_19**: 1911 tx (zero-width compact on; StringTie
  1839). Parity arc plateau; remaining +72 gap is architectural.
- **Default baseline on GGO KRAB-ZNF chromosome (NC_073244.2)**: StringTie
  3733, Rustle default 3826, Rustle `--vg` 3829. The story is per-locus
  paralog recovery, not whole-chromosome sensitivity.
- **Binary**: `./target/release/rustle` (rebuild with `cargo build --release`).
- **Active input**: `GGO_KRABchr.bam` (NC_073244.2 subset) for paralog work;
  `GGO_19.bam` for parity work.

## The 5-step advisor story

1. **S1 — parity.** Rustle ≈ StringTie on standard metrics. In progress:
   1911 vs 1839 on GGO_19; remaining +72 gap is architectural.
2. **S2 — multimap advantage.** KRAB-ZNF paralog recovery (this session).
   LOC115930085 is the falsifier — replicate across more loci before pub.
3. **S3 — family definition.** Graph-structure clustering **and** HMM-
   alignment clustering should agree with each other on known Pfam families.
   HMM side not yet built.
4. **Novel copies.** Align unmapped/supplementary-only reads to family
   HMM/graph. Scaffold exists at `--vg-discover-novel`; synthetic bundle
   creation not done.
5. **SNP/indel disambiguation.** `--vg-snp` flag exists (MD-tag parsing
   planned); no implementation yet.

## Resume points — pick one

### A. Broaden the S2 falsifier sweep
Sweep NC_073244.2 for more LOC115930085-style paralog pairs where
StringTie = 0 but Rustle-VG ≥ 1. The VG family report (10 families,
≥3 copies) is the candidate pool. Expected outcome: a short list of 3–10
loci to make the S2 claim robust.

### B. Build the HMM-on-VG scoring module (`src/rustle/vg_hmm.rs`)
Pfam HMMs are amino-acid. Recipe:
1. Write assembled-transcript FASTA from GTF + reference genome.
2. Six-frame translation per transcript.
3. Run `hmmscan` against Pfam (KRAB PF01352, zf-C2H2 PF00096 for this work).
4. Parse per-transcript best-frame bit-score.
5. Family-score matrix rows × Pfam hits — cluster agreement vs
   graph-structure clustering = the S3 falsifier.

Tradeoff: depends on `hmmer` + Pfam-A. Novel-copy detection requires the
same pipeline on reads with no mapping or only supplementary hits.

### C. Tighten the pre-assembly bundle split
Family-report refinement (step 4 above) is cosmetic — under the hood the
bundles are still 400 kb wide. A real fix splits bundles at coverage gaps
≥ N bp (N kb region with zero coverage between genes). Currently the
runoff distance is 200 bp, not triggering here.
This would also lift downstream quality (smaller graphs, tighter flow
decompositions).

### D. Novel-copy discovery (Module 1 of original VG plan)
Complete `discover_novel_copies()` — currently scaffolded. Group junction-
matched unmapped candidates by region; create synthetic bundles; feed into
the standard assembly loop. Prerequisite for the "find copies not in the
genome" part of the advisor story.

## Key files

| File | Purpose |
|---|---|
| `src/rustle/graph.rs` | `compact_zero_width_nodes`, `compact_chain_pass_through` |
| `src/rustle/graph_build.rs` | wiring (default-on zw, opt-in chain) |
| `src/rustle/vg.rs` | strand-mirror filter, gene-loci refinement, family report |
| `src/rustle/pipeline.rs` | passes `&all_transcripts` to family-report writer |
| `src/rustle/stringtie_parity.rs` | `RUSTLE_STRINGTIE_EXACT` meta-flag scaffold |
| `docs/experiments/05_loc115930085_s2_falsifier.md` | S2 narrative |
| `docs/figures/make_loc115930085_figs.py` | figure generator |
| `docs/figures/fig_loc115930085_*.png` | S2 figures |
| `docs/figures/fig_family1_refinement.png` | family refinement before/after |
| `memory/project_chain_compaction_negative.md` | chain-compact DO-NOT-ENABLE log |
| `memory/project_zero_width_compaction_shipped.md` | zw compact shipped log |
| `memory/project_S2_loc115930085_falsifier.md` | S2 audit trail |
| `memory/BENCHMARK_multicopy_gene_families.md` | 10-family benchmark plan |

## Reproducing the S2 run

```bash
# one-time: extract the right chromosome
samtools view -b GGO.bam NC_073244.2 -o GGO_KRABchr.bam
samtools index GGO_KRABchr.bam

# pipelines
stringtie -L GGO_KRABchr.bam -o /tmp/ks_stringtie.gtf
./target/release/rustle -L GGO_KRABchr.bam -o /tmp/ks_rustle.gtf
./target/release/rustle -L --vg --vg-solver none GGO_KRABchr.bam -o /tmp/ks_vg_none.gtf
./target/release/rustle -L --vg --vg-report /tmp/ks_vg_refined.tsv \
    GGO_KRABchr.bam -o /tmp/ks_vg_refined.gtf

# figures
python3 docs/figures/make_loc115930085_figs.py

# view the falsifier story
less docs/experiments/05_loc115930085_s2_falsifier.md
```
