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
- **StringTie parity experiment (cross-strand RI)**: StringTie can suppress
  opposite-strand low-coverage models via `retainedintron()` inside
  `print_predcluster` (because it clusters both strands together). Rustle runs
  predcluster per bundle+strand, so this interaction is normally missing.
  An opt-in post-pass (`RUSTLE_CROSS_STRAND_RI=1`) now applies the same
  retained-intron predicate across opposite-strand transcripts, gated by a
  StringTie-style significant-overlap matrix (the `overlaps.get(n1,n2)` gate).
  On `GGO_19.bam` vs `GGO_19_stringtie.gtf`:
  - baseline: 1911 query tx, matching tx 1634, transcript-level Sn/Pr 88.9/85.5
  - with `RUSTLE_CROSS_STRAND_RI=1`: 1856 query tx, matching tx 1627, Sn/Pr 88.5/87.7
  This fixes the STRG.309-like antisense explosion, but still costs 7 matches,
  consistent with lingering `pred->cov` (flow-derived coverage) drift at some loci.
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

## 2026-04-22 trace follow-up: next StringTie deviation after checktrf audit

### What was tested

I added a narrow `checktrf` change in [src/rustle/path_extract.rs](../src/rustle/path_extract.rs)
so that a failed long-read seed with the same outer boundaries as an already-kept
path, but a different internal splice pattern, prefers independent rescue over
immediate redistribution.

Result on `GGO_19.bam` vs `GGO_19_stringtie.gtf`:

- query transcripts: `1911 -> 1915`
- matching transcripts: `1634 -> 1638`
- transcript sensitivity / precision: `88.9 / 85.5 -> 89.1 / 85.5`

So the patch appears net-positive and is currently left in place.

### What the STRG.110-class trace showed

The targeted downstream-gene locus did **not** change under that patch.
The decisive difference is earlier than `best_trf_match`.

StringTie trace (`bundle=22118027-22157884`):

- dominant downstream-gene seeds `idx=4,8,10,11,15,16,18` all hit
  `back=fail fwd=fail outcome=checktrf`
- then only a few are rescued:
  - `t=4` rescued
  - `t=8` rescued
  - `t=18` rescued
  - `t=11`, `t=16` matched afterward

Rustle trace at the analogous locus (`SEED_STATS NC_073243.2:22149137-22155355`):

- `stored=8`, `checktrf_redist=3`, `fwd_fail=1`
- the dominant downstream-gene seed is stored directly instead of being deferred:
  - `[TRACE_LOCUS] idx=4 ... -> STORED`
- this happens because endpoint support is accepted:
  - `[TRACE_BACK_SOURCE] tf=23 first=0 last=9 ... exact=true accepted=true`
  - `[TRACE_FWD_SINK] tf=27 first=20 last=21 ... accepted=true`

That is the next major deviation. Rustle is certifying complete downstream-gene
paths from exact source-edge / sink-edge support transfrags, while StringTie is
still failing `back_to_source` / `fwd_to_sink` for the same seed family and
deferring them to `checktrf`.

### Next warranted target

Focus next on `back_to_source_fast_long()` and `fwd_to_sink_fast_long()` in
[src/rustle/path_extract.rs](../src/rustle/path_extract.rs).

The evidence-backed question is:

- when should exact `source -> first_seed_node` or `last_seed_node -> sink`
  support be considered sufficient to store a full path immediately?

Current trace evidence says this gate is too permissive for the
`22152353-22155355` family. The over-permissive endpoint acceptance inflates the
pre-`checktrf` keep-set, and that in turn makes later `checktrf` redistribution
look too eager.

### 2026-04-22 checkpoint: disable-checktrf audit

I also ran the traced loci with `RUSTLE_DISABLE_CHECKTRF=1` to isolate the base
long-recursion behavior from rescue logic.

Key result:

- the `22149137-22155355` / `22152353-22155355` class is already diverging
  **before** `checktrf`
- with the broad endpoint-stub experiment enabled, Rustle started sending the
  dominant downstream-gene seed into `LONGREC_FAIL -> checktrf`, which matches
  StringTie qualitatively
- but that experiment was too broad: it also forced many legitimate seeds in
  the same merged bundle into `checktrf`, increasing this locus from
  `seeds=17` to `seeds=28`, `fwd_fail=1` to `fwd_fail=14`, and
  `checktrf_rescued=0` to `checktrf_rescued=8`

Conclusion:

- the endpoint-certification hypothesis is real
- but the correct fix is **not** to globally reject all 2-node
  `source -> node` / `node -> sink` stubs
- the next patch must be narrower and topology-aware, likely only blocking
  endpoint-only certification when the seed is being pulled across a merged
  adjacent-gene structure with competing internal children/parents

### 2026-04-22 checkpoint: endpoint-stub gate was real but not yet keepable

I wired trace-only provenance for long-rec endpoint certification:

- `back_used_pure_source_stub`
- `fwd_used_pure_sink_stub`

and tested a narrow force-`checktrf` experiment for long-read direct paths that
were certified by pure source/sink stubs.

What held up:

- At the hotspot merged bundle, the target downstream-gene seed did move in the
  desired direction:
  - `idx=4` changed from direct store to
    `LONGREC_FAIL -> checktrf -> RESCUED_INDEP`
- With `RUSTLE_DISABLE_CHECKTRF=1`, the base long-rec layer moved closer to the
  intended StringTie shape without increasing hotspot seed count:
  - previous no-`checktrf` hotspot:
    `seeds=14`, `stored=8`, `fwd_fail=1`, `unwitnessed=2`
  - force-`checktrf` hotspot:
    `seeds=14`, `stored=6`, `fwd_fail=3`, `unwitnessed=1`

What failed:

- The same endpoint-stub gate regressed full `GGO_19` before rescue/depletion
  parity was fixed.
- Direction-agnostic asymmetric gating was too broad:
  - full run: `1946` query tx, `1605` matching tx,
    transcript-level `87.3 / 82.5`
- Restricting the gate to the observed hotspot direction
  (`hardstart=true && hardend=false`) was better but still regressive:
  - full run: `1923` query tx, `1621` matching tx,
    transcript-level `88.1 / 84.3`
- Both variants were worse than the current checked-in baseline and the earlier
  same-boundary `checktrf` patch.

Decision:

- keep the provenance/trace scaffolding
- **revert** the force-`checktrf` behavioral gate for now
- treat this as evidence that endpoint certification is upstream of the hotspot,
  but that the next actionable parity target is now **seed-failure abundance
  depletion / rescue-time depletion**, not another broad endpoint gate

Practical implication for the next pass:

- use the endpoint-stub diagnostics only to classify which seeds are suspicious
- then audit how much abundance those deferred seeds keep, when they are
  reintroduced, and whether StringTie is depleting or suppressing competing seed
  families earlier in the same bundle

### 2026-04-22 checkpoint: depletion trace says the hotspot is about residuals, not split count

I added hotspot-only abundance lifecycle tracing in:

- [src/rustle/path_extract.rs](../src/rustle/path_extract.rs)
- [src/rustle/max_flow.rs](../src/rustle/max_flow.rs)

This traces:

- max-flow depletion on overlapping transfrags
- `checktrf` matched redistribution
- independent `checktrf` rescue
- second-pass redistribution

Hotspot command:

```bash
RUSTLE_TRACE_LOCUS=22118027-22157884 \
RUSTLE_SEED_STATS=1 \
RUSTLE_DEPLETION_DIAG=1 \
target/release/rustle -L GGO_19.bam -o /tmp/rustle_trace_depletion_audit.gtf \
  2> /tmp/rustle_trace_depletion_audit.err
```

Key readout:

- the hotspot still runs in the same baseline branch:
  `SEED_STATS NC_073243.2:22149137-22155355 strand=+ seeds=17 checktrf_redist=3 fwd_fail=1 skipped_hard_boundary_low_abund=2 skipped_single_exon_from_multinode=1 stored=8 unwitnessed=2`
- strongly supported direct seeds are being depleted normally:
  - `idx=7`: `ab_before=8 -> flux=8 -> ab_after=0`
  - `idx=11`: `ab_before=2 -> flux=2 -> ab_after=0`
- the suspicious class is different:
  small deferred seeds keep enough residual abundance to survive until
  `CHECKTRF_MATCH`, then contribute fractional mass into already-kept paths

Concrete hotspot examples from `/tmp/rustle_trace_depletion_audit.err`:

- `t=10`
  - deferred by `hard_boundary_low_abund`
  - reaches `checktrf` with `abund=4.0`
  - matched only to kept path `3`
  - all `4.0` is redistributed there
- `t=12`
  - deferred by `hard_boundary_low_abund`
  - reaches `checktrf` with `abund=2.0`
  - split across kept paths `0` and `2`
  - redistribution is `1.8711 + 0.1289`
- `t=19`
  - deferred by `longrec_fail`
  - no later max-flow depletion on this transfrag was observed before `checktrf`
  - reaches `checktrf` still at `abund=1.0`
  - gets split across four kept paths:
    `0.7421 + 0.1485 + 0.0387 + 0.0707`

Interpretation:

- for this hotspot, the next deviation does **not** look like a primary
  bundle-over-segmentation / seed-count problem
- earlier endpoint-force experiments changed the downstream behavior without
  increasing hotspot seed count, which already pointed away from "wrong split
  count" as the main local blocker
- the stronger current explanation is:
  **some deferred seeds retain too much residual abundance before
  `CHECKTRF_GATE`, then get redistributed into multiple kept paths**

What this rules out:

- this is not simply a failure to zero seeds after `checktrf`
- the new traces show matched and rescued seeds being zeroed correctly after
  redistribution / rescue

What this points to next:

1. compare which deferred StringTie seeds reach `CHECKTRF_MATCH` with non-zero
   abundance versus which ones have already been drained to zero by the time
   the `CHECKTRF_GATE` runs
2. focus especially on `longrec_fail` residuals like hotspot `t=19`
3. keep bundle color-break / fragmentation on the list as a separate
   architectural target, but not as the main explanation for this hotspot
