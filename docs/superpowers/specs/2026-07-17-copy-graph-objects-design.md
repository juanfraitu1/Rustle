# Copy-Graph Objects (v1) — Design Spec

**Status:** approved (brainstorming), ready for implementation plan
**Date:** 2026-07-17
**Objective served:** O4 (reference-absent + unannotated copies) unified with O1/O2 into one family variation graph.

## Goal

Make every copy of a multi-copy family a **first-class, tagged, corroborable path** in one Bandage-loadable
family variation graph emitted by `copy_assign --phase`. Add a **REFERENCE walk** so a copy *not in the genome*
is visibly a path the reference does not take, and carry each copy's **corroboration evidence as GFA tags** so the
same file is simultaneously the picture and the queryable record.

## Motivation / context

An advisor who "only trusts what he sees in Bandage/IGV" needs the reference-absent and unannotated copies to be
*visible, checkable graph objects*, not table rows. A code audit (2026-07-17) established the current state:

- The **only** GFA emitter is `copy_assign --phase` (inline at `src/bin/copy_assign.rs:916–996`, writer `:1178–1204`,
  flag `:194`). It already emits GFA 1.1 with PSV columns as allele bubbles, copies as P-lines, reads as W-lines
  (the Canzar shared-evidence flip), and a `.colours.csv`.
- Three gaps block the target visualization: (1) **no REFERENCE walk** is ever emitted, so "the reference does not
  take this arm" is unshowable; (2) **reference-absent copies carry no distinguishing marker** — their P-line is
  identical to a reference-present copy's; (3) nodes are **single-base** (1 bp specks, no backbone), and read W-lines
  reference **adjacencies with no backing L-line** (dangling edges).
- A richer in-memory `FamilyGraph` (`src/rustle/vg_family/family_graph.rs:16–52`, with `copy_specific` presence/absence
  bubble arms + `per_copy_cov`) exists but is **never serialized**; it is the home for v2 (see Out of Scope).

This spec is **v1**: allele bubbles + backbone + reference walk + tags. **v2** (exon presence/absence arms, ORF/protein
tags) is a separate spec.

## The unifying model

One family graph; every copy is a **path**; copies differ only in the **tags** their path carries, across two
existence axes and one corroboration axis:

|              | annotated              | unannotated                 |
|--------------|------------------------|-----------------------------|
| **in genome**    | `in-genome-annotated`  | `in-genome-unannotated`     |
| **not in genome**| (rare)                 | `absent-collapsed` / `absent-divergent` |

The **reference walk** makes the vertical axis visible. A "copy not in the annotation" is the *same* object as a
"copy not in the genome" — a copy-path — differing only in the genome/annotation tags.

## Architecture

New isolated module **`src/rustle/vg_family/copy_graph.rs`**, a *pure* graph-builder and GFA serializer.
`copy_assign --phase` assembles the inputs it already has and calls into it, replacing the inline block at
`copy_assign.rs:916–996`. Isolation buys unit-testability and is where v2's exon-arms slot in.

Rejected alternatives: inlining into `copy_assign.rs` (already large); a bench post-processor (drifts from the engine —
the user chose the proper feature).

### Module interface (target shape; the plan may refine names)

```rust
pub enum CopyStatus {
    Reference,
    InGenomeAnnotated,
    InGenomeUnannotated,
    AnnotationUnknown,      // in-genome, no --gff supplied
    AbsentCollapsed,
    AbsentDivergent,
}

pub struct Corrob {
    pub reads: Option<u32>,        // RC — assigned reads (computed from ReadWalk.assigned_copy)
    pub suns: Option<u32>,         // SU — private-column count (computed here from the allele matrix)
    pub map_identity: Option<f64>, // MI — genome-map identity, passed in from the O4 projection
}

pub struct PsvColumn {
    pub col: usize,
    pub genome_pos: Option<u64>,   // from FamilyDetail.psv_col_pos
    pub ref_allele: Option<u8>,    // reference base at genome_pos (read from FASTA by copy_assign); None => unresolvable
}

pub struct CopyPath {
    pub id: String,                // e.g. "DSFAM26_copy8"
    pub alleles: Vec<u8>,          // per-column base; b'.' = gap (routes through the reference allele node)
    pub status: CopyStatus,        // set by caller (absent subtype from O4; annotation axis from --gff)
    pub corrob: Corrob,            // map_identity set by caller (O4 projection); reads/suns computed by builder
}

pub struct ReadWalk {
    pub name: String,
    pub obs: Vec<Option<u8>>,      // observed allele per column; None = unobserved (walk backbone across it)
    pub assigned_copy: Option<usize>, // for colour; None = tied / K=0 (grey)
}

pub struct CopyGraph {
    pub family: String,
    pub columns: Vec<PsvColumn>,
    pub backbone: Vec<Vec<u8>>,    // shared reference sequence between columns; len == columns.len()+1
    pub copies: Vec<CopyPath>,
    pub reads: Vec<ReadWalk>,
}

impl CopyGraph {
    pub fn to_gfa(&self) -> String;     // GFA 1.1
    pub fn colours_csv(&self) -> String; // keyed on SEGMENT names
    pub fn legend_tsv(&self) -> String;  // status -> colour, one row each
}
```

The **REFERENCE walk is derived inside the builder** from `columns[*].ref_allele` + `backbone` (it is not passed as a
copy). `reads`/`suns` are computed by the builder; `map_identity` and `status` are set by the caller (which holds the
O4 projection and any `--gff` context).

## The graph object (GFA 1.1)

- **Backbone S-nodes** — `{fam}_bb{i}`, sequence = `backbone[i]` (shared reference sequence between column i−1 and
  column i; `bb0` = left flank, `bb[N]` = right flank). Sized nodes give visible arm length.
- **PSV allele S-nodes** — `{fam}_c{col}_{base}` (sequence = the single base), one per distinct allele present at that
  column among {reference ∪ all copies}. Tagged `PO:i:{genome_pos}`.
- **L-lines** — for **every** allele present at a column: `bb{i} -> c{col}_{allele}` and `c{col}_{allele} -> bb{i+1}`.
  Because reads and copies only ever walk backbone + allele nodes, **every** P/W step is backed by an L-line → no
  dangling edges (fixes audit gap 3).
- **REFERENCE P-line** — `P  {fam}_REFERENCE  bb0+,c0_{ref}+,bb1+,...  *  ST:Z:reference`. Built from the genome's own
  bases. This is the anchor for "not in the genome."
- **Copy P-lines** — `P  {fam}_copy{ci}[_ABSENT]  <walk>  *  RC:i:_ SU:i:_ MI:f:_ ST:Z:_`. Walk uses each copy's allele
  per column (gap `.` routes through the reference allele node). Absent copies get an `_ABSENT` suffix in the name in
  addition to the `ST` tag, so they are legible in Bandage's path list without opening tags.
- **Read W-lines** — `W  {read}  {hap}  {fam}  0  {len}  >bb0>c0_{obs}>...` over the read's observed span only.

### Colouring (keyed on segment names, as Bandage expects)

`colours.csv` rows `{segment_name},{#RRGGBB}`, painting the divergence directly:

- a PSV allele node **walked by an absent copy but NOT by the REFERENCE walk** → **red** (`absent` colour) — these are
  the arms that light up as "the copy the genome doesn't have";
- a node on the REFERENCE walk → **grey**;
- other copy-divergent allele nodes → that copy's palette colour (reuse `copy_color`, `copy_assign.rs:343`).

`legend_tsv` lists each `ST` status and its colour so the picture is self-documenting.

## Corroboration tags (honesty rule)

Every tag is derived from data `copy_assign` / the O4 pipeline already computes. **If a signal is unavailable for a
copy, the tag is omitted — never faked.**

- `RC:i:` assigned read count (from `ReadWalk.assigned_copy`).
- `SU:i:` SUN / private-column count = number of columns where this copy's allele is unique among the family's copies
  (reference excluded; computed in the builder from the allele matrix).
- `MI:f:` genome-map identity (1.0 for reference-derived copies; <1 for `absent-divergent`, from the O4 projection,
  e.g. DSFAM26 `genome_id = 0.952`). Omitted if no projection ran.
- `ST:Z:` status (the taxonomy above).

## Annotation axis (`--gff`, optional)

`annotated` vs `unannotated` requires a gene model. `copy_assign` gains an **optional `--gff <path>`**:

- with `--gff`: an in-genome copy whose locus overlaps an annotated feature → `in-genome-annotated`, else
  `in-genome-unannotated`;
- without `--gff`: in-genome copies are tagged **`annotation-unknown`** (we never claim "unannotated" unchecked).

The DSFAM26 demo lands without a GFF because its headline is *absent vs reference*, not the annotation axis.

## Inputs / outputs

**Inputs** (all already available to `copy_assign --phase`): the BAM, the reference FASTA (for backbone + reference
alleles), the region, the assembled family (`copy_psv_alleles` — `copy_assign_pipeline.rs:880–885`; `psv_col_pos`;
`read_psv_obs`), the per-copy absent/collapsed status (indices ≥ `n_ref` and the `discovery_coupled` flag,
`copy_assign.rs:64–66`; collapsed/divergent subtype from the O4 hidden-copy classification), the O4 projection identity
for `MI`, and optionally `--gff`.

**Outputs:** `<out>.phase.gfa` (now with the REFERENCE walk + tagged copy paths + backbone + L-lines),
`<out>.phase.gfa.colours.csv` (segment-keyed), `<out>.phase.gfa.legend.tsv`.

## Data flow

`copy_assign --phase` → build `CopyGraph { columns (with ref_allele read from FASTA), backbone (reference sequence),
copies (alleles + status + MI), reads }` → `to_gfa()` / `colours_csv()` / `legend_tsv()` → write three files.

## Error handling / graceful degradation

- No `--gff` → in-genome copies tagged `annotation-unknown` (not an error).
- No O4 projection / `MI` unavailable for a copy → omit its `MI` tag.
- A column with no resolvable `genome_pos` / `ref_allele` → **skipped** from the graph with a logged count (do not
  guess a reference base).
- A copy with a `.` gap allele at a column → routes through the reference allele node (no evidence of divergence there).
- A family with **no** absent copies → still a valid graph (REFERENCE + in-genome copies); the feature is not
  absent-copy-only.
- Empty family / zero PSV columns → emit a header-only valid GFA and log; do not panic.

## Testing

**Unit (no external data):** construct a `CopyGraph` in memory — 3 copies (two in-genome, one `absent-divergent`
diverging at columns {2,5,7}), reference alleles, a handful of reads. Assert:
1. `to_gfa()` contains `P\t{fam}_REFERENCE`.
2. the absent copy's P-line walks allele nodes at cols 2,5,7 that differ from the reference walk's nodes there.
3. parse every P/W walk step: each adjacency `(a,b)` has a matching `L` line (no dangling edges).
4. copy P-lines carry `RC`/`SU`/`MI`/`ST`; `SU` equals the count of columns where that copy's allele is unique.
5. `colours.csv` marks the absent copy's divergent nodes red and reference-walk nodes grey.
6. backbone S-nodes exist and have length > 1.

**Integration smoke (data-gated, may be `#[ignore]`):** render DSFAM26 on the narrow region and assert a valid GFA
containing `{fam}_REFERENCE` and one `_ABSENT` copy P-line.

Demo command:
```
copy_assign --bam /home/juanfra/winloci_scratch/GGO_mm.bam \
  --fasta /home/juanfra/winloci_scratch/GGO.fasta \
  --region NC_073229.2:49040000-49075000 \
  --out /home/juanfra/winloci_scratch/dsfam26_o4 \
  --phase --dump-psv
# -> dsfam26_o4.phase.gfa  (+ .colours.csv, .legend.tsv); load in Bandage, apply colours.csv.
```
(Run FOREGROUND, serial, output under winloci_scratch — WSL2 crash rule.)

## Out of scope (future specs)

- **v2 — exon presence/absence arms:** represent a copy-specific exon as a real branch only that copy walks, by
  serializing the existing `FamilyGraph` (`copy_specific` + `per_copy_cov`). Makes structurally-divergent unannotated
  copies visible as an arm nothing else walks.
- **ORF / protein tag (`OR`)** — needs a translation/`mmseqs` step; belongs with v2.
- Cross-species graph-to-graph alignment (out of the sequence-to-graph scope entirely).

## Global constraints

- Default `copy_assign` behaviour (no `--phase`) must be **byte-identical** to today.
- No new k-mer computation (consistent with the "corner & name" k-mer position).
- GFA must be valid GFA 1.1 and open in Bandage; every walk step backed by an L-line.
- Honesty rule: no fabricated tag values; omit when unknown.
