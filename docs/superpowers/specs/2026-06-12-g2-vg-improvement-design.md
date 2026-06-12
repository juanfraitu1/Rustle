# `-G2` Annotation Flag for VG-Mode Improvement — Design

**Status:** DESIGN (approved 2026-06-12) · Branch `vg/flow-capacity-apportionment`
**Goal:** With `-G` = StringTie's GTF anchoring the base assembly (rustle base ≈ StringTie), add `-G2` = the real reference annotation, which (1) **defines paralog families** that guide VG copy recovery, and (2) emits a **per-run head-to-head scorecard** showing which annotation copies rustle-VG recovers that StringTie misses. This operationalizes the project's north star: demonstrate rustle-VG > StringTie on paralog/copy recovery.

**Relationship to the parity track:** The StringTie-parity work is banked at HEAD 41f1d71 (see `2026-06-12-parity-milestone-and-continuation.md`); we assume parity is eventually attainable and pivot here. `-G` feeding StringTie's GTF as the base guide is the operational stand-in for "StringTie is the baseline" until full de-novo parity lands.

## 1. Context (existing machinery this builds on)

- **`-G/--guide`** (`main.rs:44`, `guide: Option<String>`): one guide GTF → guided assembly via `reference_gtf::parse_reference_gtf` + `process_refguides`; sets `RunConfig.guide_mode` (`types.rs:821`).
- **`--trace-reference`** (`main.rs:346`): already reports, per reference transcript, why it is missing (NotExtracted vs filter) — the matching logic the scorecard reuses.
- **VG family machinery** (`vg_family/*`): family discovery (multimapper or `--vg-family-loci` coordinate-overlap) → copy-partition → phasing DP → fingerprint-EM. Similarity thresholds already exposed as CLI flags: `--vg-family-min-kmer-jaccard`, `--vg-family-min-poa-identity`, `--vg-family-min-primitive-jaccard`, etc. (`main.rs:455-490`).
- **`RUSTLE_VG_RESCUE_REPORT`** → `<out>.rescue.tsv`: per-VG-family copy certificate (primary/own/strict/tied/frac + `rescued` flag). The scorecard is a sibling artifact (annotation-anchored rather than evidence-anchored).
- **Validated showcase = RABL2** (multi-locus paralog family rustle-VG recovers, StringTie misses). ⚠ DAZ is retired as a false positive (HMM manufactured isoforms from DAZ1 reads; see `project_vg_wiring` memory) — do NOT use DAZ as a test oracle.

## 2. Architecture & data flow

```
reads.bam ─► bundles ─► base assembly guided by -G (stringtie.gtf)   ≈ StringTie baseline
                         │
                         └─► VG mode (--vg), when -G2 present:
                               families DEFINED from -G2 annotation
                               ├─ cluster -G2 transcripts into paralog families
                               │   (reuse vg_family k-mer/POA/jaccard similarity)
                               ├─ each annotated transcript = a copy structure
                               ├─ copy-partition → phasing DP → fingerprint-EM
                               │   (unchanged; runs against annotated copies)
                               └─► recovered copies
                         ▼
              out.gtf  +  out.vg_eval.tsv   (scorecard)
```

When `-G2` is present, family discovery sources families from the annotation instead of multimapper/coordinate-overlap discovery. The downstream copy-partition / phasing / EM machinery is unchanged — only the family *definition* input changes.

## 3. Components

### 3.1 CLI & config
- New flag `-G2 / --guide2 <PathBuf>` in `main.rs` (clap). Field `RunConfig.guide2_path: Option<String>`.
- **Constraint:** `-G2` requires `--vg` (error out otherwise), consistent with the other `--vg-family-*` flags. `-G` and `-G2` are independent paths (StringTie GTF and annotation respectively).
- `-G2` absent ⇒ zero behavior change (hard invariant; byte-identical output).

### 3.2 `annotation_families` module (new: `src/rustle/annotation_families.rs`)
- **Input:** parsed `-G2` transcripts (`reference_gtf::RefTranscript`).
- **Logic:** group transcripts into paralog families: transcripts merge into one family **only if** their shared-exon similarity ≥ `--family-exon-similarity` (the §3.5 gate) AND they occupy distinct genomic loci; each transcript is a copy structure (intron chain + locus). The similarity gate is the methodological control your advisor asked for — no family graph is formed without it.
- **Output:** `Vec<AnnotationFamily>` where `AnnotationFamily { family_id, copies: Vec<CopyStructure> }`, `CopyStructure { copy_id, chrom, strand, exons, intron_chain }`.
- **What it depends on:** `reference_gtf` (loading), `vg_family` similarity primitives (clustering). No dependency on reads/bundles.

### 3.3 VG family-source switch (`vg_family` discovery entry)
- When `RunConfig.guide2_path.is_some()`, the family-discovery step consumes `annotation_families` output instead of running multimapper/coordinate-overlap discovery. Reads are bundled and routed to the family whose copies overlap their locus; phasing DP + fingerprint-EM assign multimapped reads to specific annotated copies.
- Copies with no read support emit nothing (no fabricated isoforms — the DAZ failure mode).

### 3.4 Eval scorecard (new: emitted post-output)
- For each `-G2` copy, compute by intron-chain match (±2 bp tolerance, reusing `--trace-reference` matching):
  - `in_stringtie`: a matching chain exists in `-G` (StringTie's GTF).
  - `in_vg`: a matching chain exists in rustle-VG's `out.gtf`.
- **Verdict per copy:** `WIN` (in_vg & !in_stringtie), `tie` (both), `miss` (neither), `regression` (in_stringtie & !in_vg).
- **Output:** `<out>.vg_eval.tsv` — columns `copy_id family in_annotation in_stringtie in_vg verdict`, plus a `SUMMARY` line with counts and a `mode=guided` tag.

### 3.5 Family-graph merge similarity control (methodological rigor)

**Motivation:** a family graph must not merge two copies on position/coverage overlap alone — it needs an explicit, tunable, *reported* sequence-similarity gate, or the method is not defensible (copies that co-locate but are not homologous could be fused). rustle already gates merging by shared-exon similarity (`vg_family/family_graph.rs`: `family_merge_jaccard()` = "similarity bar for unifying homologous exons across paralog copies"; `refine_by_minimizer_jaccard()` splits position-overlap clusters by minimizer-Jaccard of the shared exons; `merge_singletons_by_sequence()`). This control exists but is not surfaced or reported. This section makes it first-class.

- **Metric — shared-exon similarity:** for two candidate copies, take the exons at homologous positions (overlapping coordinates / aligned positions in the family graph) and compute the minimizer-Jaccard of those exon sequences (the existing `refine_by_minimizer_jaccard` primitive). This is "similarity of shared exons" — it ignores position agreement and measures sequence agreement.
- **Gate:** two copies merge into the same family graph **only if** their shared-exon similarity ≥ `T`. Below `T`, they stay separate families even if they overlap on the genome.
- **CLI:** new first-class flag `--family-exon-similarity <0..1>` (default = current `family_merge_jaccard()` value; env override `RUSTLE_VG_FAMILY_MERGE_JACCARD` retained for back-compat). It sets the merge gate for **both** the `-G2` annotation-family clustering (§3.2) and the general VG family merging — one principled control for the whole method.
- **Reporting (auditability):** for each emitted family, write the **achieved** shared-exon similarity (min and mean across merged copy pairs) into a `<out>.vg_families.tsv` report (`family_id  n_copies  merge_threshold_T  achieved_min_sim  achieved_mean_sim`). This lets the user state "every family was merged at ≥ T identity, with achieved identity X" — directly answering the no-threshold suspicion.
- **`-G2` interaction:** in §3.2, `annotation_families` clustering uses this same gate: annotation transcripts merge into one family only if shared-exon similarity ≥ T. Distinct-locus check still applies (a family spans ≥2 loci).

### 3.6 Honesty caveat (claim integrity)
Because `-G2` *guides* the VG output, `in_vg` is **guided recall**, not a de-novo win — the scorecard tags `mode=guided`. A clean "rustle-VG beats StringTie de-novo" claim requires the same run with `-G2` used eval-only (no guiding). That is **iteration 2**: a `--g2-eval-only` flag that loads `-G2` for the scorecard but does NOT feed family definition, so VG recovers de-novo and the scorecard (`mode=denovo`) is a fair benchmark. Iteration 1 ships guided + the labeled scorecard; iteration 2 adds the fair-benchmark mode.

## 4. Testing

- **Unit — clustering:** synthetic 2-copy paralog GTF (two near-identical transcripts at distinct loci) → `annotation_families` yields 1 family / 2 copies; two dissimilar transcripts → 2 families / 1 copy each.
- **Unit — similarity gate (§3.5):** two copies with shared-exon similarity just above `T` merge (1 family); the same pair with `T` raised above their similarity stay separate (2 families). The `<out>.vg_families.tsv` report records `merge_threshold_T` and the achieved min/mean similarity. This is the test that demonstrates the threshold is real and reported.
- **Unit — scorecard verdicts:** synthetic `{annotation, stringtie, vg}` chain sets exercising each verdict (WIN/tie/miss/regression).
- **Integration — RABL2:** real RABL2 inputs; `-G2` with RABL2 copies → VG assigns multimapped reads to the copies → `out.vg_eval.tsv` shows a `WIN` on the copy StringTie misses, no fabricated copies on empties.
- **Invariant:** `-G2` absent → output byte-identical to current `--vg` run. `-G2` without `--vg` → clean error.

## 5. First-iteration scope (YAGNI)

Iteration 1 = §3.1–3.5 + §4:
1. `-G2` flag + config + `--vg` constraint.
2. `annotation_families` module (load + cluster, gated by the §3.5 similarity threshold).
3. VG family-source switch.
4. **Family-graph merge similarity control (§3.5):** surface `--family-exon-similarity` as a first-class flag wiring the existing `family_merge_jaccard`/minimizer-Jaccard gate, applied to both `-G2` and general VG merging, with the `<out>.vg_families.tsv` achieved-similarity report.
5. Eval scorecard (`out.vg_eval.tsv`, `mode=guided`).
6. The tests above.

**Deferred:** `--g2-eval-only` fair-benchmark mode (iteration 2); clustering-threshold auto-tuning; genome-wide scorecard aggregation; integration with `rescue.tsv`.

## 6. Open questions (resolve during planning, not blocking)
- Exact clustering thresholds for `annotation_families` (start from the existing `--vg-family-*` defaults).
- Whether the base (non-VG) output should also be scored against `-G2` (default: scorecard covers the VG output only).
