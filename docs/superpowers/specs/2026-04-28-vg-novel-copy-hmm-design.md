# Novel-copy discovery via family RNA-HMM — design

Status: design, not yet implemented.
Date: 2026-04-28.
Supersedes: the flat k-mer rescue currently in `src/rustle/vg.rs::discover_novel_copies`.

## 1. Goal

Replace the flat k-mer rescue with a **hybrid per-exon profile-HMM family
model** that:

1. Scores unmapped BAM reads against a family-aware generative model.
2. Places each scoring read on the family graph (which exon classes, which
   copy-specific branch).
3. Clusters surviving reads into **synthetic bundles** that flow through
   Rustle's existing assembly pipeline as first-class citizens.
4. Classifies every rescued read into one of five **failure-mode buckets**
   explaining why minimap2 missed it, with external minimap2 verification on
   the load-bearing buckets.

Non-goal: replace the EM/Flow solvers, replace the SNP scaffold, or change
phased-assembly behaviour. Those are separate specs.

## 2. Why this approach

The companion document `docs/UNMAPPED_FAMILY_RESCUE.md` lays out the five
reasons a read with family resemblance can still be unmapped by minimap2.
The HMM is the smallest model that can absorb all five reasons in a single
inference pass: it does not have a hard chain-score threshold, it does not
mask repetitive seeds, indels are first-class transitions, copy-specific
exons are explicit branches, and reference-absent copies show up as flat
per-column posteriors over the closest paralogs.

Architectural choice (vs full splice-graph HMM or single profile HMM over
the whole gene) is **C: hybrid per-exon profile HMMs glued by a junction
transition layer**. Rationale, as agreed in brainstorming: it matches
Rustle's existing per-exon / per-junction data model, sidesteps the "MSA of
a whole gene with copy-specific exons" problem (each exon is its own MSA,
which is tractable), and produces a per-exon emission profile that is
itself the "summary of divergence" deliverable for downstream questions.

## 3. Module layout

All new code lives under `src/rustle/vg_hmm/`, exposed through a single
`pub mod vg_hmm;` declaration in `lib.rs`. The existing `vg.rs` is left
intact except for `discover_novel_copies`, which becomes a thin dispatcher
(see §8). This quarantines the new code until proven and lets the old
flat-k-mer path remain selectable as a fallback.

```
src/rustle/
├── vg.rs                         (existing; minimal edits)
└── vg_hmm/
    ├── mod.rs                    (re-exports)
    ├── family_graph.rs           (FamilyGraph builder; §4)
    ├── profile.rs                (per-exon profile HMM; §5)
    ├── scorer.rs                 (family-graph forward / Viterbi; §6)
    ├── rescue.rs                 (read selection + synthetic bundles; §7)
    └── diagnostic.rs             (failure-mode classifier; §8)
```

## 4. `family_graph.rs` — family graph

### 4.1 Inputs

- An existing `&FamilyGroup` (already produced by `vg.rs::discover_family_groups`).
- The per-copy `&[Bundle]` slice referenced by `family.bundle_indices`.
- A `&genome::GenomeIndex` for sequence retrieval.

### 4.2 Output

```rust
pub struct FamilyGraph {
    pub family_id: usize,
    pub nodes: Vec<ExonClass>,            // exon-equivalence classes
    pub edges: Vec<JunctionEdge>,         // family-wide observed junctions
    pub copy_index: HashMap<CopyId, Vec<NodeIdx>>, // node sequence per copy
}

pub struct ExonClass {
    pub idx: NodeIdx,
    pub representative_chrom: String,
    pub representative_span: (u64, u64),
    pub per_copy_sequences: Vec<(CopyId, Vec<u8>)>, // raw genomic seq per copy
    pub copy_specific: bool,              // true if only one copy contributes
}

pub struct JunctionEdge {
    pub from: NodeIdx,
    pub to: NodeIdx,
    pub family_support: u32,              // number of copies whose splice graphs use it
    pub strand: Strand,
}
```

### 4.3 Algorithm

1. For each copy, walk its existing splice graph and emit (chrom, start, end,
   strand) per exon node.
2. Cluster exons across copies into `ExonClass` instances using two-stage
   matching:
   - Stage 1: genomic-position overlap (≥30% reciprocal overlap on the
     same chromosome and strand) groups exons.
   - Stage 2: within each Stage-1 group, refine by minimizer Jaccard
     similarity (k=15, w=10 minimizers, Jaccard ≥ 0.3) — this splits a
     position-overlap cluster into multiple classes when copies have
     diverged enough that they are not the same exon any more.
3. Exons unique to one copy after Stage 2 become `copy_specific = true`
   singleton nodes — these are the bubble branches.
4. Edges are emitted from each copy's junctions, then deduplicated; the
   `family_support` count is the number of distinct copies emitting that
   edge.

### 4.4 Edge cases

- **Tandem repeated exons within a copy:** treat each occurrence as a
  separate copy-of-copy node; do not collapse. Documented in module
  doc-comment.
- **Strand mismatch in a family:** already handled by `discover_family_groups`'s
  strand-mirror filter; `FamilyGraph` builder asserts single strand and
  errors otherwise.

## 5. `profile.rs` — per-exon profile HMM

### 5.1 Per-class model

For each `ExonClass` with ≥2 contributing copies, fit a Krogh-style profile
HMM with Match / Insert / Delete states over an MSA of the per-copy
sequences.

```rust
pub struct ProfileHmm {
    pub n_columns: usize,
    pub match_emit: Vec<[f64; 4]>,        // log-prob per A/C/G/T per column
    pub insert_emit: Vec<[f64; 4]>,       // background-like per column
    pub trans: Vec<TransRow>,             // M->M, M->I, M->D, I->M, I->I, D->M, D->D
    pub n_seed_copies: usize,             // copies that informed the MSA
}
```

### 5.2 MSA strategy

- For up to 10 copies (covers all four families benchmarked in
  `MULTI_COPY_FAMILY_PROOF.md`): banded pairwise progressive alignment,
  guided by pairwise edit distance. Implementation in pure Rust to avoid
  FFI; banded Needleman-Wunsch with band width = 2 × max pairwise edit
  distance.
- Future: replace with POA (spoa) via FFI if a family ever exceeds 10
  copies and progressive alignment quality drops.

### 5.3 Emission distributions

- Per-column match emissions are smoothed with **family pseudo-counts**:
  the prior is the family-wide base composition for that exon class
  (computed across all copies before per-column refinement). Smoothing
  weight `α = 1.0 / (1 + n_seed_copies)` so 2-copy classes lean on the
  prior and 10-copy classes are nearly empirical.
- Insert emissions use the family-wide background.
- This is the "cross-family information" channel: a column conserved in
  8/10 copies has a peaked distribution; a column varying across copies
  has a flat one.

### 5.4 Singleton classes

`copy_specific` classes (only one copy) get a degenerate ProfileHmm where
match emissions are 1-hot per column with a small smoothing floor; this
still scores reads but does not pretend to family generality.

## 6. `scorer.rs` — family-graph forward / Viterbi

### 6.1 State space

A scorer state is `(node_idx, hmm_column, M|I|D)`. Transitions are:

- Within a node: standard profile-HMM transitions (M↔I, M→D, etc.).
- Between nodes: at the last `M` state of node `u`, transition to the
  first `M` state of node `v` for each junction edge `u→v`. Edge cost is
  `−log(family_support / total_family_junctions)`.

Source/sink are virtual nodes connected to all class nodes whose copies
include them as first/last exons.

### 6.2 Two scores per read

- **Forward log-likelihood** `log P(read | family)` — used for the
  rescue threshold.
- **Viterbi best path + per-column posteriors** — used to (a) place the
  read on a specific exon-class sequence, (b) detect copy-specific
  branches it took, (c) flag low-posterior columns as candidate
  divergence sites. Posteriors come from forward+backward, run only on
  reads that pass the forward threshold (cost control).

### 6.3 Threshold

A read is **rescue-eligible** iff:

- `forward_loglik − null_model_loglik ≥ τ` where the null model is the
  family-wide background composition; default `τ = 30 nats` (subject to
  tuning on GOLGA6L7).
- AND: prefilter passed (≥ 3 family k-mer hits, as in the current
  scaffold). Prefilter is kept to bound the per-read cost; the HMM only
  evaluates reads with a coarse signal.

### 6.4 Cost notes

Forward over a banded profile HMM is O(L × W) per node where L is read
length and W is the band width; total per read is O(Σ L × W) over visited
nodes. For typical IsoSeq reads (~3 kb) and an 8-exon family this is
~10⁵ ops — cheap. The expensive piece is family-graph construction
(§4), done once per run, not per read.

## 7. `rescue.rs` — synthetic bundles

### 7.1 Pipeline

1. Iterate unmapped BAM records (re-using existing iterator in
   `discover_novel_copies`).
2. Apply k-mer prefilter against family k-mer index (re-use existing
   k-mer code; do not duplicate).
3. Run §6 forward scorer; keep reads passing threshold.
4. Run Viterbi + posteriors on survivors → record per-read best path.
5. Cluster survivors by `(family_id, implied_chrom_region)` where region
   is the union of representative spans of the path's exon classes,
   widened by ±1 kb. Cluster minimum: 3 reads (matches existing
   `min_kmer_hits`).
6. For each cluster, emit a `Bundle`:
   - `source: "vg_rescue"` (new variant in `BundleSource`).
   - `start, end` = cluster span.
   - `chrom` = path's representative chromosome.
   - `junction_stats` = junctions implied by the consensus path.
   - `reads` = the rescued reads, with their (synthesized) exon coordinates
     mapped from path Viterbi traceback.
   - A new `synthetic: bool` flag on `Bundle`, default false; true here.
7. Return `Vec<Bundle>` to the pipeline; pipeline appends them to the
   bundle list before path-extract runs.

### 7.2 Filter exemptions

Synthetic bundles bypass:

- `transcript_isofrac` filter relative to neighbouring loci (their
  coverage is not comparable to mapped-read bundles).
- Pairwise-contained filter (§7.1.3 of the assembly pipeline) when the
  containing transcript is on a different bundle — a common false-positive
  for novel paralogs that overlap an existing copy.

They do NOT bypass:

- Splice-consensus checks (GT-AG / GC-AG via genome).
- Minimum-coverage threshold (`-c`).

Documented in `bundle.rs` next to the `synthetic` flag.

### 7.3 GTF attributes for emitted novel transcripts

```
copy_status "novel"; family_id "FAM_<n>"; rescue_class "<bucket>";
rescue_kmer_hits "<n>"; hmm_loglik "<f>"; n_diagnostic_columns "<n>";
```

## 8. `diagnostic.rs` — failure-mode classifier

### 8.1 Five buckets

Per `docs/UNMAPPED_FAMILY_RESCUE.md`:

1. `below_threshold_chain` — k-mer islands exist but cumulative chain
   score below minimap2's `-s` default.
2. `seed_masked` — k-mer hits exist but they are dominated by seeds that
   would be filtered by minimap2's `-f` repetitive-seed cap.
3. `divergent` — relaxed-parameter alignment yields < 85% identity to any
   reference copy.
4. `structural` — relaxed-parameter alignment yields a primary with
   indels >= 50 bp, indicative of copy-specific structural variation.
5. `reference_absent` — even with maximally relaxed parameters
   (`-s 10 -f 0 -N 100`) no primary alignment is produced.

### 8.2 Hybrid implementation

- **Internal stage (always run, all rescued reads):** Rust reproduction of
  - minimap2's chain-score formula (sum of seed weights minus gap
    penalties, with the documented default scoring),
  - minimap2's `-f` seed-frequency cap on the family's k-mer index.
  Output: an initial classification into buckets 1, 2, or "needs external
  verification."
- **External stage (subprocess; only on reads not classified 1/2):** call
  `minimap2` as a subprocess, piping the read sequence to a small FASTA
  containing only the family's genomic region(s) ± 50 kb. Stepped
  parameters: default → `-s 20` → `-f 0.0001` → `-N 50 --secondary=yes`
  → `-s 10 -f 0 -N 100`. Record at which step (if any) a primary
  alignment appears. Buckets 3, 4, 5 are derived from the external
  stage's identity / indel-length / no-alignment outcome.

External stage is gated by `--vg-rescue-diagnostic` (default off) so that
production runs do not pay the subprocess cost. Off-mode still emits
buckets 1 and 2 from the internal stage.

### 8.3 Outputs

- `rescue_class` GTF attribute on each novel transcript (per §7.3).
- A per-family TSV alongside `--vg-report`:

```
family_id  n_rescued  n_below_thresh  n_seed_masked  n_divergent  n_structural  n_ref_absent
```

This is the table that closes the advisor's "how can it be unmapped?"
question with data.

## 9. CLI surface

Existing flags reused:
- `--vg-discover-novel`: enables novel-copy discovery (now HMM-mode by
  default).
- `--genome-fasta`: required, as today.
- `--vg-report <TSV>`: family report; gains the rescue-class columns.

New flags:
- `--vg-discover-novel-mode {kmer|hmm}` (default `hmm` once stable;
  initially `kmer` until stage-by-stage validation in §11 passes, at
  which point the default flips). Lets us A/B compare on the same input.
- `--vg-rescue-diagnostic` (default off): turn on the external minimap2
  verification stage.
- `--vg-rescue-min-loglik <FLOAT>` (default 30.0): override τ in §6.3
  for tuning experiments.

## 10. Test fixtures

In order of integration testing:

1. **Unit: family-graph construction** — fixture: 2 hand-crafted bundles
   with 3 exons each, two shared and one copy-specific. Assert 4
   `ExonClass` nodes, correct edges, correct `copy_specific` flags.
2. **Unit: profile HMM emissions** — fixture: 3 toy sequences over an
   8-column MSA with 1 fully conserved column. Assert that column's
   match emission is peaked (`P(observed base) > 0.95`) and a divergent
   column is near-uniform.
3. **Integration: GOLGA6L7 chr19 triple smoke** — full Rustle run on
   `GGO.bam`, `--vg-discover-novel-mode hmm`. Assert: family graph for
   GOLGA6L7 has 9 exon classes, no copy-specific singletons, ≥ 1
   synthetic bundle emitted, `rescue_class` populated. Cross-check
   diagnostic counts against `tools/golga6l7_snp_probe.py`.
4. **Synthetic ground truth** — take one annotated GOLGA copy, perturb 5%
   of bases, insert a copy-specific exon, simulate 50 reads via a small
   read simulator (or even fixed perturbations of real reads). Assert:
   ≥ 80% of simulated reads rescued, classified as bucket 3 or 4,
   external minimap2 confirms (does not align under default params,
   does align under relaxed).
5. **Negative control** — random low-complexity unmapped reads,
   contam-like sequences. Assert: rescue rate ≤ 1%.

Tests live under `tests/regression/vg_hmm/` to match existing test
layout. The synthetic fixture's small inputs check into git; the GOLGA
fixture references the existing `GGO.bam` test data.

## 11. Rollout

1. Land §4–§6 with unit tests; run alongside but do not invoke from
   `discover_novel_copies`.
2. Land §7 behind `--vg-discover-novel-mode hmm`; default remains `kmer`.
   Validate on GOLGA6L7 fixture.
3. Land §8 internal stage; gate external via `--vg-rescue-diagnostic`.
4. Run full benchmark on `GGO.bam`; if novel-copy count and
   `MULTI_COPY_FAMILY_PROOF.md` exact-match counts hold or improve, flip
   default mode to `hmm`.
5. Update `MULTI_COPY_FAMILY_PROOF.md` and `ALGORITHMS.md §10` with the
   new architecture and a rescue-class table.

## 12. Risks and mitigations

| Risk | Mitigation |
|---|---|
| MSA-per-exon-class fails when copies are too divergent in a class — produces nonsense profiles | Detect at build time (mean pairwise identity < 60%): split into per-copy singletons (i.e. an explicit bubble) instead of one class. |
| Synthetic bundles disrupt downstream filters and inflate noise | `synthetic: bool` flag with documented filter exemptions (§7.2); emit `source "vg_rescue"` so all downstream code can identify them. |
| External minimap2 subprocess is slow on large unmapped pools | Already gated behind `--vg-rescue-diagnostic`; only buckets 3/4/5 candidates reach it; production runs leave it off. |
| HMM rescue produces phantom novel copies in regions covered by reference | `family_support` weighting in §6.1 deprioritises paths through reference-supported exons unless emissions strongly prefer them; synthetic bundles still go through splice-consensus and min-coverage gates (§7.2). |
| τ threshold (§6.3) is hard to set globally | `--vg-rescue-min-loglik` flag for tuning; default chosen on GOLGA6L7 ground-truth ROC. |

## 13. Out of scope (deferred to later specs)

- `--vg-snp` completion (indel-aware diagnostic alleles). Will reuse the
  per-column emission distributions from §5 as its likelihood term.
- Replacing the EM solver's `P(r | k)` with HMM forward score. Reuses
  §6's machinery but is a much larger refactor of `run_pre_assembly_em`.
- Phased / HP-tag interaction with the family graph. The current
  scaffold splits bundles by HP tag before VG processing; that ordering
  may need to flip if we want per-haplotype HMMs.
