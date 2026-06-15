# Design: family variation super-graph, amended by multimappers (`RUSTLE_VG_RECOVER_COPIES`)

**Date:** 2026-06-14. **Status:** design for review. **Supersedes** the earlier per-copy
read-reassignment draft (that was structurally baseline + assignment, not a variation graph).

**One-liner:** For each multi-copy gene family, build **one shared variation graph** (common
backbone + variant nodes at the PSVs that distinguish copies; each copy is a *path*). Amend that
graph with the **multimapper reads' allele evidence** — a read's PSV genotype, pooled across all
its alignments, phases it onto its true copy-path, overriding minimap2's arbitrary primary flag.
This recovers paralog copies the `-G` baseline starves. Gated, default-off, byte-identical when
off. This is the advisor's exact interest — *multimappers improving graphs and assigning reads to
very similar copies* — realized as a genuine variation graph, not read-clustering.

## 1. Why the super-graph (not per-copy reassignment)

Per-copy read-reassignment recovers the same copies but produces independent StringTie-style
graphs — structurally baseline + assignment (Canzar's "resolve the ambiguity" frame). The
**family super-graph** is the genuine contribution: copies are *paths through one shared graph*,
and multimappers are *shared evidence amending that graph* (the FLIP of the resolve-ambiguity
frame). It also pools backbone-only (allele-ambiguous) reads on the shared structure instead of
discarding them, so thin copies that strict per-read assignment drops can still phase.

## 2. The super-graph substantially EXISTS — this is reposition + complete, not build

| piece | exists | file |
|---|---|---|
| family super-graph (shared exon-class nodes, junction edges, **copies = paths**) | ✅ `FamilyGraph`, `recover_paralog_path(copy)` | `vg_family/family_graph.rs` |
| family construction (homology clustering, jaccard, minimizer, POA/MSA backbone) | ✅ `build_family_graph`, `poa_msa` | `family_graph.rs` |
| PSV variant nodes on the graph | ✅ `psv_columns_for_family`, `PsvColumn`/`PsvCopyAllele` | `psv_linkage.rs` |
| read genotyping over **all alignments incl. multimapper secondaries** (molecule-level) | ✅ `genotype_family_reads` (keeps primary+secondary+supplementary; votes global per read) | `psv_linkage.rs:318` |
| PSV-based copy assignment (phase read → copy by alleles + margin) | ✅ `assign_read_to_copy(g, min_psv, margin)` | `psv_linkage.rs` |
| per-copy isoform assembly on the shared graph | ✅ `assemble_psv_isoforms`, `link_junctions` | `psv_linkage.rs` |
| identifiability gate | ✅ `family_identifiability(fg, error_rate, min_psv)` | `psv_linkage.rs` |
| MEC haplotype phasing (bounded DP) | ✅ `mec_dp(matrix, n_sites, cap)`, `phase_reads` | `vg_family/phasing.rs` |

**Why it is currently inert (the thing to fix):** PSV-linkage is wired as an *additive channel
alongside Part A* — it re-emits chains Part A already produced (union-by-chain dedups them → net
0), and the `(E)` identifiability gate confines it to the regime where Part A is already correct.
It never targets the copies the `-G` baseline **misses**.

## 3. The redesign: target the MISSED copies, amend the graph, decompose into paths

Reposition the existing machinery from "additive side-channel" to "**core recovery of copies
absent from the `-G` output**", driven by multimapper allele evidence.

1. **Build the family super-graph** (`build_family_graph`) for each multi-copy family (homology
   clustering of Layer-1 graphs). Shared backbone + junction edges; copies are paths.
2. **Place the PSV variant nodes** (`psv_columns_for_family`): the positions where copies differ
   (derived from the copies' reference sequences + read evidence).
3. **Amend with multimappers** (`genotype_family_reads`): ingest **all** alignments of every
   family read — crucially the **secondary** alignments of reads whose primary minimap2 placed at
   the sister. Each read's PSV genotype is the union over its alignments (molecule property).
4. **Phase reads onto copy-paths** (`assign_read_to_copy`): a read pins to the copy its alleles
   match decisively (≥`min_psv` supporting PSVs, `margin` over the next copy); allele-ambiguous
   (backbone-only) reads stay pooled on the shared structure (coverage for all copies).
5. **Decompose into copies** (`assemble_psv_isoforms`/`link_junctions`): each family copy-path,
   carrying its phased reads + pooled backbone coverage, yields its isoform(s).
6. **Recover only the MISSED** : emit a copy-path iff it is **absent from the `-G` output** AND
   passes the phantom guard (§5). This is the change from the inert wiring — recover what `-G`
   lacks, not re-emit what it has. Strictly additive (VG ⊇ baseline preserved).

**Advisor framing (no clustering):** the family graph is built by homology (his bundles/graphs,
collapsed across copies by sequence); multimappers *amend* that graph via their allele evidence;
copies are *paths*; reads are *assigned to copies* by their variants. Nothing here clusters reads
to invent loci — the copies are the homologous loci, and phasing is standard allele assignment.

## 4. Flag, gating, invariants

- **Flag:** `--vg-layer2-recover-copies` → `RUSTLE_VG_RECOVER_COPIES=1`. Requires `--vg`.
- **Default OFF**, output **byte-identical** to current `--vg --vg-layer2` when off.
- **`RUSTLE_PRECISE`-exempt**; **strictly additive** (only appends copy-paths absent from output).

## 5. Phantom guard (the DAZ3 lesson)

A copy-path is emitted only if it is **PSV-decisive**: ≥`K` reads phased to it by `assign_read_to_copy`
with the required margin (variant evidence, not coverage alone). Allele-tied loci (homology
shadows) accumulate no decisively-phased reads → never emitted. Gated additionally by
`family_identifiability` (the `(E)` gate) so we only attempt families where PSVs can separate the
copies above the error floor. This is stronger than the AS-tie guard — it uses the actual variants.

## 6. Components / files

- **Reuse (the bulk):** `family_graph.rs`, `psv_linkage.rs`, `phasing.rs` — as above.
- `src/rustle/vg_family/layer2.rs` (rewire): add a recovery path that, per multi-copy family,
  runs build-graph → PSV columns → genotype (incl. secondaries) → phase → assemble, then keeps
  **only copy-paths absent from `all_transcripts`**; merge additively (before identical-dedup).
- `src/rustle/vg_family/secondary_index.rs` (reuse): enumerate the family's multimapper reads/loci.
- `src/rustle/main.rs`, `types.rs`, `pipeline.rs`: flag plumbing (mirror `--vg-layer2-*`).
- Likely small extensions: ensure `genotype_family_reads` is invoked with the multimapper loci in
  guided mode (today's inert wiring may scope it differently); a "is this copy already emitted?"
  overlap check before keeping a recovered path.

## 7. Validation methodology

- **Target:** the 23 multimapper-recoverable copies (`/tmp/cre_guided/headroom.py`), plus any
  thin copies the pooled backbone lets phase.
- **Realness:** SQANTI3 (FSM) + authenticity guard (decisive/tied/phantom) via the copy-recovery
  eval harness. **Success:** material fraction of the 23 recovered as SQANTI3-FSM + `authentic`,
  **0 phantom** emissions. Flag-off vs flag-on diff = the recovered copies (additive).

## 8. Testing (TDD)

- Synthetic 2-copy family fixture (`bench/sim/` IsoSeq 2-copy paralog, commit 2ec1745): copy B
  starved (reads' primary at A, secondaries at B carry B's PSV alleles). Flag-on: B's copy-path
  phases from the secondary allele evidence and emits; flag-off: B absent. Byte-identical off.
- Unit: `genotype_family_reads` aggregates a read's PSV votes across its primary+secondary alignments.
- Unit: `assign_read_to_copy` pins a PSV-decisive read; leaves an allele-tied read unassigned.
- Unit: recovery keeps a copy-path absent from output; never re-emits an already-present copy.
- Phantom guard: an allele-tied (homology-shadow) copy is never emitted; `RUSTLE_PRECISE` unaffected.

## 9. The hard core, and what we ship first

The thesis-elegant endpoint is a **joint constrained flow-decomposition** of the super-graph —
infer copy-count + per-copy isoforms together under allele-linkage constraints (one path per copy,
PSV-consistent), à la min path-cover with guarantees (the advisor's taste). The **tractable
version we ship first** uses the existing pieces: copies are the *known* homologous loci, their
PSV signatures are reference-derived, so phasing is per-read assignment (`assign_read_to_copy`) +
per-copy assembly on the shared graph — no combinatorial blow-up (this sidesteps the previously
hanging full DP; `mec_dp`'s cap is available if joint phasing is needed). The joint-decomposition
upgrade is then a localized replacement of steps 4–5, on the same super-graph.

## 10. Risks

- **Identifiability bound (honest ceiling ~23–50):** copies with no PSV-spanning reads cannot be
  phased on any graph; the ~364 allele-tied copies stay capped. The super-graph pools ambiguous
  reads to rescue *thin* copies but does not beat the data limit — and that limit is itself a
  clean, defensible thesis statement.
- **Re-inerting:** if scoped like the old side-channel it nets 0 again; the load-bearing change is
  targeting **missed** copies via multimapper (secondary) allele evidence. Tests pin this.
- **Phantoms (DAZ3):** PSV-decisive emission + `(E)` gate + SQANTI3/authenticity validation.
- **Family-graph construction cost** on large families: bounded by the existing jaccard/POA paths;
  run only for multi-copy families with a candidate missed copy.
