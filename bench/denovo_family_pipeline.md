# De-novo Family + Copy-Assignment Pipeline (the "read-coherence way")

Genome-wide, **annotation-free, minimizer-free** pipeline on gorilla (GGO) HiFi/IsoSeq data, built to
serve the advisor's two objectives — **(1) identify multi-copy gene families at the RNA level** and
**(2) assign reads to copies under ambiguity (hard multimappers)** — plus a third lever that fell out of
the data: **(3) copy-specific / allele-specific junctions** as an assignment signal.

Branch: `vg/flow-capacity-apportionment`.
Commits: Python arc `b8abbbf`; Pass-1 `2698457`; Rust copy-assign core `f5def11`.
Big data lives in `/home/juanfra/winloci_scratch/` (NOT `/mnt/c`): `GGO.bam`, `GGO.fasta`,
`denovo_skeletons.tsv`, `denovo_transcripts.{fa,meta.tsv}`, `sim5x/`.

---

## Architecture: a general assembler with a copy-aware specialization

It is a **general-purpose assembler** (assembles *all* transcripts, family or not) whose "areas of focus"
are multi-copy families and PSVs. Pipeline:

```
BAM ─► Pass-1 read-coherence skeletons ─► gate ─► GENERAL ASSEMBLY (101k→93k transcripts)
                                                        │
                                          intron-junction collapse ─► 15,676 gene loci
                                                        │
              canonical-kmer pre-filter ─► contiguous-span ─► POA contiguous-core ─► family edges
                                                        │
                                weighted-modularity decomposition ─► 1,193 families (1,190 + 3 webs)
                                                        │
                         family-aware RESCUE (borrow strength) ─► +118 under-assembled copies
                                                        │
                  COPY ASSIGNMENT: per-read PSV + copy-specific-junction likelihood + ID gate
```

---

## Stages (script · input · output · key params)

### 1. Pass-1 read-coherence skeletons — `twopass_denovo_gw_pass1.py` (committed `2698457`)
Group **primary** reads by **exact intron chain** → de-novo transcript skeletons. `MIN_READS=2`.
Out: `/home/.../denovo_skeletons.tsv` (165,994 skeletons; 4.4M reads scanned, ~66s).
NB: exact-chain grouping is correct for HiFi — splice sites are exact (see Rejected: fuzzy collapse).

### 2. Gate → general assembly — `denovo_assemble_gate.py`
Keep skeletons with `>=3` reads + **all-canonical junctions** (GT-AG/GC-AG/AT-AC, consistent strand) +
spliced length in `[100, 300_000]` and **genomic span `<= 3 Mb`**.
⚠ The span/spliced caps were a *bug fix*: without them the gate emitted chromosome-spanning artifacts
(max "transcript" 236 Mb → 3.65 GB FASTA). With caps: 101,467 transcripts, max 16 kb, 346 MB.
Out: `/home/.../denovo_transcripts.{fa,meta.tsv}`.

### 3. Family detection — `denovo_families.py`  → `denovo_families.tsv`, `denovo_family_edges.tsv`
1. **Collapse isoforms → gene loci** by **intron-junction union-find** (two transcripts share a locus iff
   they share an identical `(donor,acceptor)`). NOT raw span overlap — span overlap single-linkage-chains a
   whole dense chromosome into ~1 "locus" (the original 93k→27 bug). → 15,676 loci.
2. **Pre-filter ("bloom filter" idea, done exactly):** `KMER=18` **canonical** (RC-symmetric) k-mers;
   `np.unique(return_counts)` gives exact distinct-gene ownership; family-informative = owned by `2..40`
   reps. A rep needs `>=6` informative k-mers to be a candidate (single-copy genes rejected → never POA'd).
   This cut C(15676,2)=123M all-pairs to 237,833 pairs and fixed the original OOM (old code indexed *all*
   k-mers incl. millions of unique singletons). KMER=18 chosen so coincidental sharing is negligible
   (~0.07% load); KMER=13/16 were saturated.
3. **Contiguous-span filter:** keep a pair only if its shared k-mers SPAN `>= T_CORE * min(len)` (cheap
   proxy for POA's core_recip, same threshold). 237,833 → 25,540 POA pairs (drops domain-sharers).
4. **POA contiguous-core = THE decision** (`poa_family_definition.poa_pair_stats`, `core_recip >= T_CORE=0.13`).
   **Strand-aware:** canonical k-mers (step 2) + POA tries **both orientations** (`_poa` RC fallback) →
   recovers copies assembled on inconsistent strands. ~11,400 confirmed edges → `denovo_family_edges.tsv`.
   POA is the slow part: ~30 min for ~25k pairs, 5-proc pool, O(L²) (190ms@2kb, 2.2s@8kb).

### 4. Decomposition — `denovo_family_split.py` → `denovo_families_split.tsv`
Connected components transitively close a **superfamily** (KRAB-ZNF) into one blob. A real recent-duplicate
family is a **dense, mutually-homologous** subgraph; an over-merge is a **sparse chain bridged by
lower-core_recip domain edges**. Decompose components `>= 6` by **weighted modularity** (networkx
`louvain_communities`, weight=`core_recip`, resolution=1.0). Robust: resolution 1→8 holds concordance
96.0–96.7% while max family 174→65. Per-family structural columns (density, avg_core_recip, n_articulation);
**flag webs** (size>=10 & density<0.15) — non-discretizable homology continua, not claimed as families.
Result: **1,193 families = 1,190 discrete + 3 webs**, name-concordance 97.3% (names never used in the decision).

### 5. Family-aware RESCUE — `family_rescue.py` → `denovo_rescued_copies.{tsv,fa}`
**Borrow strength / partial pooling:** a single confident canonical-junction read forming a multi-exon
chain that **POA-confirms against a confirmed family** is a real copy, even below the `>=3` gate. Scan
merged family neighbourhoods (`WIN=1Mb`), **canonical-kmer pre-filter** (`>=20` shared k-mers with the
best member), POA vs that **single best member** (both orientations), rescue if `core_recip >= 0.13`.
**118 rescued** (88 single-read, 25 two-read), genome-wide, without lowering the global gate.
⚠ Perf gotcha: POA against *all* members of large ZNF families is the cost — POA only the best k-mer match.

### 6. COPY ASSIGNMENT — `copy_assign.py` (modes: `sim5x`, `real`)
Per-read **PSV + copy-specific-junction likelihood** behind an **identifiability gate**:
- Discover PSVs from copies' spliced sequences (`discover_psvs`, all-pairs vs copy[0]; a copy that MATCHES
  the ref inherits that allele — NOT `None`, else the likelihood favours copies it shares no columns with;
  this was a real bug). Each column = `{copy: (gpos, seq_base)}`.
- `logL(read|copy) = Σ_PSV log P(obs base | allele, e) + Σ_junctions (±jw if copy has/lacks the boundary)`.
- **Resolvable iff** the read spans `>= 1` feature where copies differ (the theorem). Margin over runner-up
  → Assigned / Ambiguous; 0 decisive → Tied ("N equally good places").
- Junction boundaries are computed strand-agnostically via `gen2off` (`max(off(d-1), off(a))`).
Rescued copies merged into rosters via `load_rescued`.

---

## Key results

- **General assembler:** 93,373 transcripts → 15,676 gene loci.
- **Families:** 1,190 discrete + 3 flagged webs. De-novo recovers textbook families with names never used:
  MAGEA(15), MAGEB(13), GSTM(10), PCDHB, KRT, TRIM, HSPA, IFITM, and dense ZNF subfamilies.
- **Known-family validation:** RABL2 ✓, APOBEC3 ✓, **RFPL ✓ (unified after strand fix + RFPL2 rescued)**.
- **Copy assignment — synthetic (sim5x, ground truth):** K=0 → 100% tied; K>=2 → 100% accurate;
  error-floor sweep 1.0→0.805 as e→0.15. The realistic BAM-driven assigner reproduces the oracle curve.
- **Copy assignment — real GGO (25 co-located families, 25,940 reads):** **99.9% silver-standard agreement**
  (vs minimap2's confident mappings), 95% resolvable, 82% confidently assigned.
- **Junctions are decisive for near-identical families:** DSFAM43 (5 copies, 95% MAPQ-0) goes 10% → 99%
  resolvable when copy-specific junctions are added (+167 reads PSVs alone could not resolve).
- **Rescue:** 118 under-assembled copies recovered; RFPL2 (a single read) rescued into the RFPL family.

---

## Rejected approaches (don't re-try without new evidence)

- **Minimizers as the criterion** — never used; only a *contrast* panel in `poa_family_definition.py`. The
  decision is exact canonical k-mers + POA core. Advisor dislikes minimizers.
- **Fuzzy (splice-site-tolerant) intron-chain collapse** — MEASURED net-NEGATIVE (−3,988 transcripts at
  `>=3`). HiFi splice sites are exact; `<=8bp` differences are real isoforms, not noise. Does not recover
  under-assembled genes.
- **Lowering the `>=3` gate** — RFPL2 has 4 reads / 4 distinct support-1 chains; no gate helps; only
  family-aware rescue (with the family prior) recovers it.
- **Forward-only family detection** — under-merges any mixed-strand family (RFPL). Fixed by strand-awareness.

## Determinism / scale / gotchas

- Determinism: fixed seeds; no `Date.now`/`random` in derivations; networkx Louvain seeded.
- WSL2 box: **5 cores, 19 GB**. POA pools use 5 procs. Watch memory; prior runs OOM'd the candidate dict
  before the bloom pre-filter. `kmer_hashes` uses a **vectorized Horner** rolling hash, NOT int64 matmul
  (numpy integer matmul is a non-BLAS slow loop).
- The de-novo FASTA has some **duplicate transcript IDs** (same `DN_chrom_start_nexon` for distinct chains);
  handled by a length-guard in `copy_assign`; the 99.9% agreement shows no material harm.
- Honest residual: RFPL1's reads are a messy readthrough — *represented* (LOC134758217) but not cleanly
  delineated. A data-quality limit, not an algorithm gap.

## Re-run (order)

```bash
PY=/home/juanfra/miniforge3/bin/python
cd bench
$PY twopass_denovo_gw_pass1.py     # -> denovo_skeletons.tsv (~66s)
$PY denovo_assemble_gate.py        # -> denovo_transcripts.{fa,meta.tsv}
$PY denovo_families.py             # -> denovo_families.tsv + denovo_family_edges.tsv (~30 min POA)
$PY denovo_family_split.py         # -> denovo_families_split.tsv (fast, reads edges)
$PY family_rescue.py               # -> denovo_rescued_copies.{tsv,fa} (~2 min)
$PY copy_assign.py sim5x           # synthetic identifiability validation
$PY copy_assign.py real            # real-family copy assignment (merges rescued copies)
```

---

## Rust port (status & plan)

Goal: port the advisor-aligned cores into `src/rustle/vg_family/`. Rust has `poasta` (POA, no pairwise
aligner), `noodles` (BAM), FNV-1a minimizers + exact k-mer Jaccard, clap CLI, rayon, BTreeMap determinism.
Existing `copy_split.rs` does PSV *discovery* (`ReadObs`, `AlignedRead`, `allele_at` CIGAR bridge, 19 tests).

- ✅ **Copy assignment** — `src/rustle/vg_family/copy_assign.rs` (commit `f5def11`). Pure `assign_read(read,
  copies, params) -> {best_copy, log_lr_margin, n_decisive, resolvable, status∈{Assigned,Ambiguous,Tied}}`.
  PSV likelihood + junction term + identifiability gate; reuses `allele_at`; 13 tests incl. the sim5x
  K=2/K=1 logic and the ref-allele-inheritance guard. Full crate compiles (454 existing tests intact).
- ⬜ **Family-aware rescue** — POA-homology borrow-strength. Use `poasta::poa_msa([a,b])` for a 2-seq MSA →
  derive contiguous-core/`core_recip`; try both orientations. Medium.
- ⬜ **Strand-aware family detection** — add a **canonical** exact-k-mer counter (existing code uses
  minimizers); contiguous-span; POA core via poasta; **weighted-modularity decomposition** needs a
  community-detection approach (petgraph + Louvain, or a small implementation). Largest.
- ⬜ **Integration** — junction-boundary mapping (per-copy exon context, `gen2off`-style) + BAM/FASTA
  orchestration to run copy assignment inside the pipeline (bundle/family/transcript structs in
  `types.rs`/`vg.rs`/`pipeline.rs`).

Build note: test builds are memory-heavy on the 19 GB box — use `CARGO_PROFILE_TEST_DEBUG=0 -j4`.

## File index

Scripts: `twopass_denovo_gw_pass1.py`, `denovo_assemble_gate.py`, `denovo_families.py`,
`denovo_family_split.py`, `family_rescue.py`, `copy_assign.py`, `poa_family_definition.py` (shared POA core),
`build_sim5x.py` (5-copy synthetic + oracle), `sim_reads.py`.
Tables: `denovo_families{,_split,_annotated}.tsv`, `denovo_family_edges.tsv`, `denovo_rescued_copies.tsv`.
Rust: `src/rustle/vg_family/copy_assign.rs` (+ `copy_split.rs` it builds on).
