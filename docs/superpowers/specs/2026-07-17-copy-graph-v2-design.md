# Copy-Graph v2 — Exon Presence/Absence Arms + MI Tag — Design Spec

**Status:** approved (brainstorming), ready for implementation plan
**Date:** 2026-07-17
**Builds on:** v1 (`copy_graph.rs` + `copy_assign --phase`, shipped `bf47732`).
**Objective served:** O4 (copies-not-in-annotation made visible) — a copy that differs by EXON CONTENT becomes an arm nothing else walks; plus the `MI` corroboration signal v1 left honesty-omitted.

## Goal

Two additions to `copy_assign --phase`, sharing one upstream plumbing change:
1. **`MI` tag** — fill the genome-map-identity corroboration tag (currently always `None`) for reference-absent copies.
2. **Exon presence/absence graph** — a NEW sibling GFA (`<out>.exon.gfa`) where each copy is a walk over its genomic EXON classes; a copy-specific exon (esp. one carried by a reference-absent / unannotated copy) is an arm only that copy enters, and the shared **REFERENCE walk** skips it via a skip edge.

## Motivation / grounding (what a 4-agent code audit established, 2026-07-17)

- **"Serialize `FamilyGraph`" is dead.** `FamilyGraph`/`ExonClass`/`copy_specific`/`per_copy_cov` (`family_graph.rs`, `consensus.rs`) are **test-only**, never populated in the shipped pipeline; `per_copy_cov`'s doc references `annotate_per_copy_exon_coverage`, a function that does not exist. v2 must BUILD the exon structure, not serialize that layer.
- **The data we need is computed upstream then discarded.** Each copy's genomic exon chain lives on `DenovoTranscript.introns` (genomic donor/acceptor pairs, `family_detect.rs:63-76`) and is retained on `ColocatedFamily.copies` during assignment (`denovo_pipeline.rs:229`, used at `:1474/:1544/:1592`), but at the `FamilyAssignment` build sites only the whole-copy envelope `copy_spans=(chrom,start,end)` (`:1641`) and the SPLICED-length `copy_junctions` (`:1650`) survive. `copy_junctions` (cumulative exon lengths, `copy_assign_pipeline.rs:90-103`) can COUNT an extra exon — which is why `chi_h_with_junctions` (`readonly_copy_number.rs:54-82`) already separates such copies — but cannot PLACE or align exons across copies.
- **A copy-specific exon has ZERO PSV columns** (a PSV needs a multi-copy aligned column, `copy_assign_pipeline.rs:385`), so it is invisible to v1's PSV-column graph; forcing it in is a module rewrite (v1 `walk_tokens` visits every column, no skip edges). The clean shape is a SEPARATE exon-level graph at genomic (donor/acceptor) coordinates.
- **The absent-copy remap identity IS computed then thrown away.** `absent_copy::remap_identity_minimap2` (`absent_copy.rs:221-297`, from the minimap2 `de:f:` tag) is bound at the admission gate (`absent_copy.rs:153`), compared to `remap_max_identity`, and discarded — `Admission::Copy(DenovoTranscript)` (`absent_copy.rs:72-77`) carries no identity. In-genome copies are never remapped (identity ≈1.0 by construction; honest choice is to omit, not fake).

## Design

### §1 — Shared plumbing: thread genomic per-copy data into `FamilyAssignment`

- **Preserve the remap identity.** Change `Admission::Copy(DenovoTranscript)` → `Admission::Copy(DenovoTranscript, Option<f64>)` (or a named struct). At `absent_copy.rs:153`, pass the bound `id`. Update the match arms at `denovo_pipeline.rs:1233` and `:1514` and the test literal at `:3417`.
- **Two new per-copy fields on `FamilyAssignment`** (`denovo_pipeline.rs:377`), each parallel to and length-aligned with `copy_tids`:
  - `copy_map_identity: Vec<Option<f64>>` — absent copies → `Some(remap_id)`, in-genome copies → `None`.
  - `copy_introns: Vec<Vec<(u64,u64)>>` — genomic donor/acceptor chains (`all_copies[k].introns.clone()`).
- **Populate at EVERY `copy_tids.push` / construction site** to keep alignment (the audit flagged a pre-existing gap where `admit_novel_pools_with_admitter`, `denovo_pipeline.rs:1240-1244`, pushes `copy_tids` without `copy_spans` — the new fields must be pushed here too):
  - main build `denovo_pipeline.rs:1623-1653`: `copy_introns` from `all_copies[k].introns`; `copy_map_identity` = `None` for ref copies, `Some(id)` for the absent tail (capture the admitted ids from the stage-2 admission at `:1508-1537`).
  - `FamilyAssignment::empty()` (`:449`): both fields `Vec::new()`.
- **Constraint:** these fields are consumed ONLY by the v2 code below. No existing writer reads them, so default `copy_assign` output stays byte-identical.

### §2 — MI tag (small)

In `build_copy_graph` (`copy_assign.rs:360`, `Corrob` at `:420`), replace `map_identity: None` with `map_identity: fa.copy_map_identity.get(ci).copied().flatten()`. Renders `MI:f:{:.3}` on absent-copy P-lines in BOTH `<out>.phase.gfa` and the new `<out>.exon.gfa`. In-genome copies keep `None` → tag omitted (honesty rule, `copy_graph.rs:49,192`). No fake `1.0`.

### §3 — Exon presence/absence graph (`<out>.exon.gfa`)

New sibling types in `copy_graph.rs` (the v1 `CopyGraph` is UNTOUCHED):
```rust
pub struct ExonClass { pub chrom: String, pub start: u64, pub end: u64 }   // a node = one genomic exon interval
pub struct CopyExonPath {                                                    // a copy = an ordered walk over exon classes
    pub id: String, pub exon_nodes: Vec<usize>, pub status: CopyStatus, pub corrob: Corrob,
}
pub struct ExonGraph { pub family: String, pub nodes: Vec<ExonClass>, pub copies: Vec<CopyExonPath> }
impl ExonGraph {
    pub fn to_gfa(&self, exon_seq: impl Fn(&ExonClass) -> Vec<u8>) -> String;  // exon_seq fetches the reference bases
    pub fn colours_csv(&self) -> String;
    pub fn legend_tsv(&self) -> String;
}
```

**Build (`build_exon_graph` in the wiring, from `fa` + the genome):**
1. Per copy, reconstruct genomic exons: `exons_of` (`copy_assign_pipeline.rs:43-53`) from `copy_introns[ci]` + the copy's outer bounds `copy_spans[ci]` (first exon start / last exon end).
2. **Cluster exons into shared classes:** reciprocal-overlap union-find (same chrom+strand, overlap fraction ≥ 0.30) over all copies' exons — implement in the v2 code (the proven logic is `family_graph::cluster_by_position`, `:170`, but that module is dead/test-only; lift the algorithm, do not call it). Sort clusters by genomic position → backbone `E0..En`.
3. Each `CopyExonPath.exon_nodes` = the sorted class indices that copy contributes an exon to. A class with exactly ONE contributing copy = a **copy-specific exon** (owned by that copy).

**Emit `<out>.exon.gfa` (GFA 1.1):**
- `S` per exon class: `{fam}_E{k}`, sequence = reference bases over `[start,end)` via `GenomeIndex::fetch_sequence` (`genome.rs:113`, reached through the `genome_for` closure), tag `PO:i:{start}`. A per-class read-support tag `RC:i:` computed inside `to_gfa` as the SUM of `corrob.reads` over the copies whose `exon_nodes` include that class (no field on `ExonClass`; derived from `self.copies`).
- `L`: the UNION of consecutive-class adjacencies across the REFERENCE walk and every copy walk. This covers both directions of structural difference uniformly — a copy that INSERTS an exon contributes the detour `Ea→Eb`,`Eb→Ec` while the reference contributes the skip `Ea→Ec`; a copy that SKIPS a reference exon contributes `Ea→Ec` while the reference contributes `Ea→Eb`,`Eb→Ec`. Because every walk step is thus a real L-line, nothing dangles.
- `P` **REFERENCE walk** `{fam}_REFERENCE` = the exon classes present in ≥1 in-genome (non-absent) copy, in genomic order (tag `ST:Z:reference`); if there are NO in-genome copies, fall back to the classes present in ALL copies. The anchor a copy-specific arm departs from.
- `P` per copy `{fam}_copy{ci}[_ABSENT]` = its `exon_nodes` walk, tags `RC:i:`(reads) `MI:f:`(map id) `ST:Z:`(status). (`SU` is a PSV-column concept — omitted here.)
- `--gff` (existing, §Task-8 v1): overlay `ST` on exon classes / copies — an exon class covered by no annotated feature that a copy walks = `in-genome-unannotated`; this is how "not in the annotation" is labelled. Without `--gff`, classes carry `annotation-unknown`.

**`<out>.exon.gfa.colours.csv`** (segment-keyed): a copy-specific arm node → its owner copy's `CopyStatus::colour()` (red for an absent copy, per v1); shared/reference exon classes → grey (`#9aa0a6`). **`<out>.exon.gfa.legend.tsv`** as v1.

Triggered by the SAME `--phase` flag (one run emits both `.phase.gfa` and `.exon.gfa`); the exon graph builds only if `copy_introns` is non-empty. Individual read W-lines are OUT OF SCOPE for the exon graph (they clutter at exon scale; per-exon read support is the `RC:i:` node tag; the v1 `.phase.gfa` keeps read W-lines at PSV scale).

### §4 — Error handling / degradation

- Family with no exon differences → exon graph is a linear chain (all copies same classes); still valid, no arms.
- `copy_introns[ci]` empty (single-exon copy) → that copy walks one class; fine.
- No absent copies → no `MI` tags, no red arms (all shared); valid.
- A class whose `fetch_sequence` returns `None` (out of contig) → emit the node with an `N`-run of the interval length + log; do not panic.
- Reference/other copies genuinely disjoint (dispersed paralogs, no genomic overlap) → each copy's exons are its own classes; the graph shows parallel chains rather than a shared backbone (honest — there is no shared exon to anchor). Log this case.

### §5 — Testing

- **ExonGraph builder unit tests** (synthetic, injected `exon_seq`): 3 copies where copy1 has an extra internal exon → assert (a) a copy-specific class node exists with only copy1 in it, (b) copy1's walk includes it while REFERENCE + copy0/copy2 take the skip edge, (c) both skip and detour L-lines exist (no dangling — reuse a walk-parsing check like v1's `assert_no_dangling`), (d) the arm node is coloured by copy1's status, (e) `RC:i:` per-class support is correct.
- **MI test:** a `FamilyAssignment` with `copy_map_identity=[None, Some(0.952)]` → `build_copy_graph` puts `MI:f:0.952` on copy1's P-line and omits it on copy0.
- **Plumbing round-trip:** a `detect_and_assign` (or `admit_novel_pools`) test asserting `copy_introns`/`copy_map_identity` are populated length-aligned with `copy_tids` (incl. the admit-novel-pools site).
- **Real-data render (data-gated):** GWFAM61 — human CNTNAP3 (9 exons) / CNTNAP3B (8 exons), Soto ID_245, chr9; BAM `/mnt/linuxdisk/home/juanfraitu/winloci_data/soto_reads.bam`, ref `.../chm13v2.0.fa`. Assert `<out>.exon.gfa` has a REFERENCE walk and one copy-specific arm. ⚠Verify the differing exon is INTERNAL (cassette) not a partial-assembly terminal exon (compare the two copies' genomic intron chains). On-thesis alternative: a gorilla reference-absent copy from `bench/denovo_rescued_copies.tsv` (retains genomic `introns`) carrying the arm.

## Out of scope (v2.1+)

- Individual read W-lines in the exon graph (per-exon `RC:i:` depth stands in).
- Computing a real self-identity for in-genome copies (would add a minimap2 pass to `--phase`); they stay `MI`-omitted.
- Reviving the dead `FamilyGraph`/`consensus.rs` layer (v2 builds its own focused types).

## Global constraints

- Default `copy_assign` (no `--phase`) MUST be byte-identical; the new `FamilyAssignment` fields are consumed only by v2 code.
- New `FamilyAssignment` fields MUST be length-aligned with `copy_tids` at every construction/push site.
- Valid GFA 1.1; every P walk step backed by an L-line (skip + detour edges both emitted) — no dangling.
- Honesty rule: `MI` omitted when unknown, never faked (in-genome copies stay omitted).
- No new k-mer computation. v1 `CopyGraph` and `.phase.gfa` output unchanged.
