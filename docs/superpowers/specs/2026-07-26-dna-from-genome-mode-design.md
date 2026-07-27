# `--from-genome` mode: reproduce Soto's families from the genome alone — Design

**Date:** 2026-07-26
**Status:** design (approved 2026-07-26)

## Goal

Add a genome-only family-detection front-end to `gw_family_catalog` that discovers segmental-duplication
gene families **from the CHM13v2.0 genome sequence alone — no annotation, no reads** — and reproduces
(almost all of) the Soto 2025 benchmark families. The DNA front-end feeds the **exact same** family-grouping
core (`homology_edges_all_reps` → `gamma_quasi_clique_partition`) that the RNA `--homology-primary` path uses.
Running both modes on the same benchmark regions then produces one apples-to-apples table — same engine, two
substrates — showing that RNA's lower recall is the **substrate/question**, not the method.

## Why (the argument this serves)

The advisor's question is "does the method work?" — proven, in his terms, only by reproducing Soto. But
Soto's own inputs were the **finished T2T-CHM13 assembly** (copies already separated by the assembler) plus
**WGS read-depth** (famCN/parCN). Soto classified already-separated loci; our RNA pipeline must *deconvolve
ambiguous reads into copies*. So:

- **Give the method the genome** (pre-separated copies) → it reproduces Soto.
- **Give it RNA** (must deconvolve) → it hits the identifiability ceiling (~85%).
- **Same grouping engine both times.** The only thing that changes is the front-end that discovers loci.

That is the concrete, undeniable form of "he is asking a different question."

### What this proves — and what it does NOT (written here so it is not over-claimed)

- **Proves:** (1) the family-grouping engine is correct — it reproduces Soto from the genome; (2) the RNA
  shortfall is the substrate (splicing discards the divergent intron/flank sequence that separates copies),
  not a broken method.
- **Does NOT prove:** that RNA copy-detection equals Soto. That is impossible by construction — Soto consumed
  an assembled genome + DNA depth; RNA has neither. No experiment can deliver it, so the spec does not promise it.
- **Mild, disclosed circularity:** DNA mode uses Soto's family *regions* as search windows (scope decision
  below). It discovers the copies and families *within* those windows de novo (from sequence, not from the
  member list), and the Soto BED is used only to score. Genome-wide discovery (removing even the window prior)
  is a named later extension, not this spec.

## Scope

**The 83-family / 362-member Soto benchmark** (`bench/soto/80_fams.chr.bed` on CHM13v2.0), scored the same way
the RNA side is scored, so the DNA-vs-RNA comparison is apples-to-apples. NOT the full 213-family catalog
(separate, much larger project).

## Architecture — one engine, two front-ends

```
                    ┌─ RNA front-end (EXISTS): reads → skeletons → assemble_gate → collapse → reps
 input ─────────────┤
                    └─ DNA front-end (NEW): genome FASTA + windows → self-align → duplicated-locus reps
                                                                                        │
                                          reps: Vec<DenovoTranscript>  ◄────────────────┘
                                                                                        ▼
   SHARED CORE (UNCHANGED): homology_edges_all_reps(&reps) → gamma_quasi_clique_partition(reps.len(), edges, 0.20)
                            → family catalog (families.tsv / copies.tsv) → Soto scorer
```

The RNA and DNA reps differ in exactly one meaningful way, and it *is* the scientific claim:

- **RNA rep:** `seq` = spliced exon-sum (introns removed), `introns` = the splice chain.
- **DNA rep:** `seq` = the **full genomic locus sequence (introns included)**, `introns` = **empty**.

## Components

### New: `src/rustle/vg_family/from_genome.rs`

A single new module with one public entry:

```rust
/// Discover duplicated genomic loci within `windows` by self-alignment, return them as reps
/// (genomic sequence, empty intron chain) for the shared homology-grouping core. Read-free, annotation-free.
pub fn genome_reps(
    fasta_path: &str,
    windows: &[(String, u64, u64)],   // search neighborhoods (benchmark family regions)
    p: &GenomeRepParams,
) -> Result<Vec<DenovoTranscript>>;
```

Three steps:

1. **Locus discovery (SD detector).** For each window, extract its sequence and align it against the genome
   FASTA with minimap2 (`-c --eqx -N100 -p0.02`, the recipe validated in the 2026-07-25/26 member map-back
   that recovered 98–100% of members). Keep alignments with identity ≥ `min_identity` (default **0.90**) over
   an aligned block ≥ `min_block` (default **1000 bp**). Each retained target segment that is **non-self**
   (does not overlap the query window) marks a duplicated locus; the query window itself is also a locus when
   it has ≥1 such duplicate. This is the annotation-free, read-free analog of Soto's SD98 map-back.
2. **Rep construction.** Per window, merge overlapping duplicated target segments into candidate loci
   (single-linkage by genomic overlap — the same rule `cluster_unspliced` uses). For each candidate locus
   `(chrom, start, end)`, fetch its genomic sequence and emit
   `DenovoTranscript { tid: format!("DN_{chrom}_{start}_1"), chrom, start, end, seq: <genomic>, introns: vec![],
   strand: '+', n_reads: <#self-align hits supporting the locus>, distinguishing_uniq: 0 }`.
   Strand `'+'` is a placeholder — homology grouping is strand-agnostic (minimap2 tries both), matching how
   the RNA reps are compared.
3. **Return** the rep vector. No grouping here — that is the shared core's job.

`GenomeRepParams { min_identity: f64 (0.90), min_block: u64 (1000), max_locus_span: u64 (MAX_SPAN=3_000_000) }`,
`from_env` reads `RUSTLE_GENOME_MIN_IDENTITY` / `RUSTLE_GENOME_MIN_BLOCK` so the SD-identity floor can be
swept (e.g. 0.98 for a strict "SD98" run) without recompiling.

**Two distinct identity thresholds — do not conflate them:**
- `GenomeRepParams.min_identity` (default 0.90) is the **locus-discovery** floor: "is this segment duplicated
  enough to be admitted as a *candidate locus*." It governs step 1 only.
- The homology-grouping `--min-identity` (E_r floor, default from `RefineParams`; 0.98 = SD98 mode) is the
  **family-membership** floor inside the shared core: "are two admitted loci homologous enough to be *one
  family*." It governs the reused `homology_edges_all_reps`.
  They are separate steps and may hold different values (a locus can be admitted at 0.90 yet only join a
  family when two loci clear the grouping floor).

### Modified: `src/rustle/bin/gw_family_catalog.rs`

Add a `--from-genome <windows.bed>` flag (mutually exclusive with the BAM input). When set:
- `--bam` is not required; `--fasta` and `--from-genome <bed>` are.
- The binary calls `genome_reps(fasta, windows, &params)` for the reps, then runs the **same** homology
  grouping + catalog emission the `--homology-primary` path already runs on RNA reps. `--min-identity 0.98`
  continues to mean the SD98 identity floor (already wired for the homology edges).
- Emits `<out>.families.tsv` / `<out>.copies.tsv` in the identical format to the RNA path.

The RNA code path is untouched and remains the default; `--from-genome` is a separate, opt-in entry.

**Where the DNA path joins and where it stops.** It **skips** the read-based front-end entirely
(`primary_reads_from_bam`, `pass1_skeletons_robust`, `assemble_gate`, the readthrough/mis-chain filters, and
`collapse_loci_span_aware` — all of which consume reads). It **joins** at the reps → homology-grouping step
(`homology_edges_all_reps` → `gamma_quasi_clique_partition`) and **stops** after emitting the family catalog
(`families.tsv` / `copies.tsv`, one copy row per grouped rep). It does **not** run the within-family
read-conflict / PSV / copy-assignment (O2) steps — those require reads and are irrelevant to reproducing the
family *definition*, which is the whole claim. A DNA family is emitted when its grouped component has ≥2
distinct loci (the same `>=2 copies` requirement the RNA catalog applies).

### Reused unchanged (the "same engine")

- `homology_edges_all_reps` (`denovo_pipeline.rs:2810`) — minimap2 all-vs-all over `rep.seq`.
- `gamma_quasi_clique_partition` (`family_split.rs`) with γ = 0.20 — the exact call at
  `denovo_pipeline.rs:1583`.
- The families/copies TSV writers.
- The Soto scorer (`bench/soto/soto_cache_score.py` / the per-member overlap scorer) — a DNA family "recovers"
  a Soto member when a ≥2-copy DNA family contains a locus overlapping that member, identical to the RNA rule.

## Data flow

```
80_fams.chr.bed (windows) ─┐
CHM13v2.0.fa ──────────────┴─► genome_reps ─► reps ─► homology_edges_all_reps ─► gamma_quasi_clique_partition
                                                                                          │
                                                                    families.tsv/copies.tsv
                                                                                          │
                                                        soto scorer (vs 80_fams.chr.bed)  ▼
                                                                    member recall / precision
```

## The deliverable

1. `gw_family_catalog --from-genome 80_fams.chr.bed --fasta chm13v2.0.fa --min-identity 0.90 --out dna_mode`
   → member recall vs Soto (expected ~90–98%).
2. The existing RNA result (85%, `bench/soto/member_attribution.final.tsv`).
3. A comparison table + short honest-framing note (`bench/soto/dna_vs_rna_mode.md`): same engine, DNA vs RNA,
   with the "genome hands you separated copies; RNA must deconvolve" framing and the caveats above.

## Testing

1. **Unit — `genome_reps` on a synthetic 2-copy FASTA:** a tiny in-repo reference with one sequence present
   twice at ≥95% identity in two windows → `genome_reps` returns 2 reps with non-empty genomic `seq` and empty
   `introns`; downstream `gamma_quasi_clique_partition` groups them into 1 family / 2 copies.
2. **Unit — rep fields:** emitted reps have `introns.is_empty()` and `seq.len() == end - start` (genomic, not
   spliced) — guards the core scientific invariant.
3. **Byte-identity guard:** a default (no `--from-genome`) run on an existing RNA fixture is byte-identical to
   the pre-change output — the new front-end must not perturb the RNA path.
4. **Integration — known family:** `--from-genome` on the SRGAP2 (`ID_462`) and PMS2P (`ID_8`) windows recovers
   each as one family with its copies (these are clean, well-separated families).
5. **End-to-end:** the full 83-benchmark run produces `families.tsv`/`copies.tsv` and a scored recall; recorded
   in `dna_vs_rna_mode.md`.

## Non-goals

- The full 213-family / 1002-paralog genome-wide catalog (separate project).
- Genome-wide de-novo discovery without window priors (named extension, not this spec).
- Copy-number (famCN/parCN) — Soto's read-depth leg; out of scope, we reproduce the *family grouping*, not CN.
- Any change to RNA-path behavior.
