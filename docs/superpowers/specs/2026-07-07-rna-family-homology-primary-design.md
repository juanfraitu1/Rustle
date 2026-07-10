# RNA multi-copy family definition: homology-primary (E_r)

**Date:** 2026-07-07
**Status:** design approved (brainstorming); implementation plan pending
**Working dir:** `/mnt/c/Users/jfris/Desktop/Rustle` (branch `collaborator`)
**Companion:** `2026-07-07-rna-family-crossmodal-validation-design.md` (Liftoff `-copies`, SEDEF SD98, Soto famCN/parCN) — separate spec.

## 1. Problem

The shipped genome-wide RNA family catalog (`gw_family_catalog --cross-chrom --refine`,
`detect_conflict_catalog_genome_wide_xchrom` in `src/rustle/vg_family/denovo_pipeline.rs`) defines a
family as a connected component of the **read-conflict graph** `E_c` (loci a read cannot de-tie),
optionally validated by exon-sum homology in `refine_families_exon_sum`. A recall audit against the
phased-assembly oracle (`bench/diploid_cn_oracle.tsv`, `asm_hapCN` = true copy number) found the
conflict-primary definition misses real, expressed, homologous families by **two** mechanisms:

- **refine-strict** (e.g. ZNF727/736, KRAB-ZNF ZNF14/682/429, MAGEA10, GSTM3): the conflict family
  *does* form, but the copies sit at exon identity **0.62–0.78**, below the `--refine` `id ≥ 0.80` gate,
  so refine drops them. `refine` can only prune within a component; it never creates a link.
- **E_c-no-edge / "globin"** (e.g. GSTM1/2/4/5, the ZNF92 cov-602 hub): paralogs ~80–88% identical but
  divergent enough that minimap2 confidently picks one copy (`mapq0_prim = 0`), so no de-tie edge forms
  and no family is built at all — the homology shows up only in secondary alignments.

Using the read-*conflict* (indistinguishability) signal as the *definition* oracle is inverted for the
goal: it biases the catalog toward copies we cannot resolve and leaks the copies we can. The thesis
already reframed this — RNA family = **transcript-homology component `E_r`** (uniform with DNA `E_b` and
protein `E_p`); read-conflict belongs *inside* the family as copy-number `χ(H)` and O2 assignment. This
spec brings the shipped catalog in line with that reframe.

## 2. Goal / non-goals

**Goal.** Define an RNA multi-copy family as a homology component so the catalog captures **both**
regimes the pipeline is meant to model: families whose reads *distinguish* the copies (O2-resolvable)
and families whose reads are *indistinguishable* (K=0 / χ(H)). Recover the audit's ancient-paralog and
no-edge misses. Keep the family *definition* on a single nucleotide (RNA) axis, with protein and DNA as
**orthogonal** confirmation only.

**Non-goals.** (a) Not changing O2 copy-assignment (`copy_assign` / `colocated_families`) — it stays as
the within-family layer and keeps its tight same-strand/disjoint-loci input. (b) Not *separating reads*
onto K=0 near-identical copies (indistinguishable by definition) — but their genomic loci ARE enumerated
by genome projection (§7 = famCN), and per-copy read assignment (parCN) is attempted where PSVs exist.
(c) ~~Not using protein as a definition edge~~ SUPERSEDED 2026-07-08: protein IS an opt-in
definition edge (`--protein-tail`) for the divergent coding tail below the ~0.65 nt limit — see §2 tier 3. (d) Cross-modal validation (SEDEF/Soto tables, real Liftoff
cross-check) lives in the companion spec — §7 *implements* the projection mechanism in-engine; the
companion spec *validates* it against the published tool.

## 3. Design

### §1 Architecture

```
reps  ─►  exon-sum homology graph E_r (over ALL reps)  ─►  γ-quasi-clique families (≥2 distinct loci)
             (nucleotide only)                                │
                                                              ├─►  WITHIN each family:
                                                              │      conflict graph → χ(H) copy-number
                                                              │      → O2 assignment (PSV/SUN)
                                                              └─►  genome-projection (§7, Liftoff-style):
                                                                     consensus → genome @≥0.98
                                                                     → famCN (incl. K=0 collapsed) → parCN
```

Homology is the primary/definition edge (regime-agnostic: links copies whether or not reads distinguish
them). Read-conflict moves inside the family. `E_r` is a superset of today's conflict families (any
conflicting pair is homologous), so nothing currently found is lost. Rep construction
(`pass1_skeletons_robust → assemble_gate → collapse_loci_span_aware`) is unchanged.

### §2 Edge oracle (RNA nucleotide homology only)

An `E_r` edge between two reps' exon-sum sequences, **all tiers gated by coverage-of-shorter ≥ 0.50**
(the repeat-bridge guard — a shared Alu is < 50% of a real transcript):

| tier | method | identity floor | recovers |
|---|---|---|---|
| 1 | minimap2 `-x asm20` (nt) | `id ≥ 0.80` | recent paralogs (current behaviour) |
| 2 | minimap2 `-k11 -w5` (nt) | `id ≥ ~0.60` (from 0.70) | ancient paralogs (ZNF14↔429=0.62, ZNF727↔736=0.74, GSTM3=0.77) |

Reuses `nucleotide_edges` in `denovo_pipeline.rs`. The tier-2 floor is the recall lever; its exact value
is **fixed by the existing family precision/recall sweep** against the annotation oracle, not hand-picked
— so the threshold is pinned to a measured operating point. Coverage gate + γ-quasi-clique operator (§4)
absorb the added repeat-bridge risk.

**Tier 3 — protein (UPDATED 2026-07-08, opt-in `--protein-tail`).** Originally protein was QC-only; a
proof (`bench/sim_allprimary.py`) showed nt alignment gives out at ~0.65 pairwise identity — beyond that,
copies map all-primary (no read-conflict possible) *and* nt homology fails, so a coding family in that tail
is invisible to every nt/read signal. Protein stays conserved there, so protein is now an OPTIONAL
**definition edge**: `homology_edges_all_reps` unions `batch_protein_edges` (longest-ORF → mmseqs, fident ≥
0.50, cov ≥ `min_coverage`) over all reps when `--protein-tail` is set. Default OFF (needs mmseqs; nt-only
otherwise). This reverses the earlier "protein = QC only" scope for the divergent coding tail; the separate
orthogonal per-family `protein_coheres` QC flag (§6) still exists and is unaffected.

### §3 Tractable all-vs-all (two-source candidate prefilter)

All-vs-all over genome-wide reps is O(n²); shortlist candidate homologous pairs first, from the union of:

- **(a) secondary-alignment co-mapping** — free from the BAM: for each read, (primary locus, each
  secondary locus) is a candidate pair. This is exactly the Mechanism-2 homology the conflict graph
  discards, surfaced directly.
- **(b) minimizer / MinHash index** over exon-sum sequences (existing `minimizers` / `vg_repeat_catalog`)
  — catches homologous pairs that share no reads.

Each candidate pair is then **confirmed** by the tier alignments (§2). Only confirmed pairs become edges.

### §4 Operator: γ-quasi-clique (GAMMA = 0.20), uniform with E_b / E_p

Family = a **γ-quasi-clique** (internal edge density ≥ 0.20) with **≥ 2 spatially-distinct loci**
(`family_definition::distinct_loci`), via the density-gated recursive splitter from
`bench/genome_family_def.py` (`_refine_component` / `_split_once`, seed-fixed NP-hard witness). Port to
Rust in `family_split.rs`, reusing the existing `louvain_communities`. This is the same operator the DNA
(`E_b`) and protein (`E_p`) families use — the thesis's "one criterion, four axes." It fixes the two
structural miss modes the audit found:

- a tandem **array is a dense clique → stays whole** (kills the ANKRD18/ZNF92 fragmentation that Louvain
  on the *conflict* graph caused, where across-subgroup conflict = 0 chopped one array into 3 families);
- **repeat-bridge chains are sparse → split out**, which is the guard that makes dropping the nt floor to
  ~0.60 safe.

### §5 Conflict / χ(H) / O2 → within-family

Inside each `E_r` family, build the conflict graph (`conflict_edges` over the family's reps) → compute
**χ(H)** = number of mutually-unresolvable copies (the copy-number, including K=0 collapsed copies) →
run O2 assignment (PSV/SUN: resolvable → assign, else abstain). Reuses `family_copy_number` and
`copy_assign`, now scoped inside the homology family. The distinguishable/indistinguishable split becomes
a within-family property — the stated modelling goal.

### §6 Orthogonal validation hooks (definition stays single-axis)

- **Protein QC (orthogonal, never an edge).** After families are built on nt homology, one batched
  `mmseqs` run (`batch_protein_edges`, longest-ORF) emits a per-family **`protein_coheres` flag**:
  does the nt-defined family also form a protein-homology cluster? Confirms coding families and red-flags
  families with no protein coherence (candidate spurious nt merges / non-coding). Optional (`mmseqs`
  absent → flag = `NA`, no effect on membership).
- **Catalog-diff interface.** Emit the catalog in a form the companion validation spec can diff against
  DNA (SEDEF SD98) and Liftoff `-copies` catalogs (stable `family_id`, per-copy `chrom/start/end`,
  exon-sum FASTA — already produced by `gw_family_catalog`).

### §7 Genome-projection copy enumeration (famCN / parCN) — Liftoff-inspired, reference-free

**Motivation.** Near-identical copies (≥ ~98%) are indistinguishable at the RNA read level (K=0): the
assembler collapses them into fewer reps (audit: GSTM2 `asm_hapCN` 19 → ~6 reps). E_r groups whatever reps
exist, so the family is found but *undercounted*. Liftoff (Shumate & Salzberg, *Bioinformatics* 2021)
proves that projecting a **gene model onto the genome** enumerates near-identical copies reliably — exactly
the regime where genome alignment beats RNA. We adopt the **mechanism**, in-engine via minimap2 (no Liftoff
dependency), seeded by our **de-novo family consensus instead of reference annotation** — Liftoff's proven
enumeration *without* its reference-annotation circularity (the reason Liftoff cannot find novel families).

**Mechanism (per E_r family):**
1. Consensus = exon-sum sequence of the family's best-supported copy (POA consensus later if needed).
2. Project onto the genome with splice-aware minimap2 (`-x splice -N <k>`), all secondary hits.
3. Keep a hit as a **genomic copy locus** iff: identity ≥ `projection_identity` (default **0.98** — the
   Soto/very-similar regime, a separate operating point from the §2 definition floor); colinear and covers
   ≥ 0.90 of the consensus (structure-preserving — the Liftoff discipline that rejects partial/repeat
   hits); and disjoint from every other kept locus and from the family's known RNA copies.
4. **famCN** = |RNA-known copies ∪ projection-found loci|. Projection-only loci are flagged
   `collapsed/projection_only` (present in the genome, merged or silent in RNA).
5. **parCN** = the within-family read layer (§5 conflict/PSV/SUN) assigns RNA reads across the loci; loci
   with no distinguishing PSV stay K=0 (counted in famCN, abstained in parCN).

**Output (per family):** `famCN`, and per copy locus `(chrom, start, end, source ∈ {rna, projection},
rna_reads, parCN_assignable)` — the Soto famCN/parCN pair, produced reference-free.

**Wiring:** new `project_family_copies(consensus, genome, ProjectionParams) -> Vec<CopyLocus>` (new
`src/rustle/vg_family/genome_projection.rs`), called after E_r families are built. Gated by
`--enumerate-copies` (default ON in Soto `--min-identity 0.98` mode). `projection_identity` default 0.98.

### Modes

- `--min-identity <f>` on the E_r definition. Default = the P/R-tuned operating point (~0.60 tier-2).
  **`--min-identity 0.98` = "Soto SD98 mode"**: tier-1 only, a high-precision recent-duplication catalog
  directly comparable to Soto et al. SD98 (family ≈ famCN, χ(H) ≈ parCN).

## 4. Wiring

- New `detect_homology_catalog_genome_wide` in `denovo_pipeline.rs`, parallel to
  `detect_conflict_catalog_genome_wide_xchrom`: build reps → `E_r` edges (§2–§3) → γ-quasi-clique
  families (§4) → within-family χ(H)/O2 (§5). Returns `Vec<Vec<DenovoTranscript>>` (+ per-family flags).
- `gw_family_catalog` gains `--homology-primary` (and `--min-identity`). The conflict path
  (`--cross-chrom`) stays available for A/B comparison and is not removed.
- `colocated_families` and the O2 `copy_assign` path are untouched.

## 5. Test plan (TDD — write failing test first)

Unit (in `denovo_pipeline.rs` / `family_split.rs`):
1. `e_r_edge_keeps_ancient_paralog`: two reps at nt id 0.62, cov 0.87 form an edge; at cov 0.20 (Alu-like)
   do **not**.
2. `gamma_quasi_clique_keeps_array_whole_splits_repeat_chain`: a 5-node dense clique stays one family; a
   sparse A–B–C bridge chain splits.
3. `within_family_chi_h_counts_collapsed_copy`: a family with a K=0 collapsed locus reports χ(H) > n_reps.
4. `homology_primary_superset_of_conflict`: every conflict-linked pair is also E_r-linked.

Integration / regression:
5. On a fixture (or the chr-scoped real subset), the audit misses **ZNF727/736, KRAB-ZNF trio, MAGEA10,
   GSTM1/2/4/5** are now grouped (each in one family with its siblings).
6. `project_family_copies` on the GSTM consensus against NC_073224.2 returns **> the RNA-rep count** of
   disjoint ≥0.98 loci (famCN recovers collapsed copies); on a single-copy gene it returns 1.
7. Family P/R sweep vs annotation oracle: precision does not regress below the agreed floor at the chosen
   tier-2 identity (KRAB-ZNF is the stress case).
8. `bench/sim_run.sh` ground-truth stays 100% (planted families still recovered, 0 over-merge).

## 6. Risks & mitigations

| risk | mitigation |
|---|---|
| precision loss at nt floor ~0.60 (KRAB-ZNF over-merge) | coverage ≥ 0.50 gate + γ-quasi-clique split + P/R sweep pins the floor |
| O(n²) all-vs-all | two-source prefilter (§3) bounds to candidate pairs |
| `mmseqs` dependency | protein QC is optional/orthogonal; absent → `NA`, no membership effect |
| K=0 collapse undercount (GSTM2) | out of scope by design: family found, count = χ(H)/DNA depth |
| γ-quasi-clique port correctness | port `genome_family_def._refine_component` with the seed-fixed witness; unit-test against the Python operator on shared fixtures |

## 7. Accepted limitations

- Truly synonymous-divergent coding paralogs (nt < ~0.60 but high protein id, extreme RABL2B-class) are
  **not formed** by the nt-only definition; the protein QC flag *surfaces* them and DNA can confirm.
- K=0 near-identical copies cannot have their *reads separated* (parCN abstains), but their genomic loci
  are enumerated by §7 projection (famCN). Projection identity floor (0.98) means copies below 98% but
  above the definition floor are found by E_r, not projection — the two layers meet in the middle.

## 8. Mode coexistence (decided) + open questions

**Decided:** `--homology-primary` ships **alongside** the existing `--cross-chrom` conflict path — both
modes are kept. The P/R sweep + cross-modal validation (companion spec) decide which becomes the default
*later*, on evidence; neither is removed in this work. This gives a direct A/B on the same substrate.

Open questions (resolved during implementation):
- Exact tier-2 identity floor (P/R sweep) and whether tier-1 `asm20` stays at 0.80 or merges into one
  tuned tier.
- Whether the within-family conflict/χ(H) layer (§5) runs by default in homology-primary mode or behind a
  flag for the first cut.
