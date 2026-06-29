# Denovo Pipeline (consolidated)

> Merged from 7 source docs (verbatim, git keeps the originals' history). Each section below was a separate `bench/*.md`.

**Contents:** [denovo_family_pipeline](#denovo-family-pipeline) · [denovo_families.SUPERSEDED](#denovo-familiessuperseded) · [twopass_denovo](#twopass-denovo) · [integrate_end_to_end](#integrate-end-to-end) · [intronchain_discovery](#intronchain-discovery) · [readcoherence_finding](#readcoherence-finding) · [readcoherence_psv_headroom](#readcoherence-psv-headroom)


---

## denovo_family_pipeline

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
- ✅ **Family-aware rescue** — `src/rustle/vg_family/family_rescue.rs`. Pure decision core:
  `canonical_kmer_set` (exact base-4 KMER=18, `min(fwd,rc)`, N-drop — a faithful port of
  `denovo_families.py::kmer_hashes`, the missing canonical-k-mer counter), best-member pre-filter by
  k-mer overlap (`K_RESCUE=20`, strict-`>` earliest-wins), then `rescue_thin_locus` POA-confirms against the
  single best member via the already-ported `family_graph::contiguous_core_coverage` (`core_recip`), trying
  BOTH orientations (RC fallback), rescue iff `core_recip >= T_CORE=0.13`. 14 tests (TDD). The BAM-neighbourhood
  scan that builds thin loci is deferred to Integration (below). Adversarially verified faithful vs the python
  (encoding bit-for-bit over a 2000-seq sweep; decision logic confirmed); fixed one latent defect the review
  caught — lowercase/soft-masked input passed the case-insensitive k-mer pre-filter but `reverse_complement`
  mapped lowercase→N and voided the RC fallback, so POA operands are now uppercased (mirrors python
  `poa_pair_stats.upper()`).
- ✅ **Strand-aware family detection** — `family_detect.rs` (`denovo_families.py` core) + `family_split.rs`
  (`denovo_family_split.py` core). `family_detect`: intron-junction union-find loci collapse, canonical
  exact-k-mer counter (`canonical_kmer_first_pos`, position-compacted over N-dropped windows), ownership
  pre-filter + contiguous-span, both-orientation POA edges via `contiguous_core_coverage`. `family_split`:
  connected components + **hand-rolled deterministic weighted-modularity Louvain** (local-moving with the
  python-louvain gain formula + multi-level aggregation), iterative-Tarjan articulation points, structural
  metrics, web flagging. Chose hand-rolled over petgraph+a community crate (no Rust Louvain dep exists;
  networkx-byte-identical is infeasible; the target is same modularity DEFINITION + determinism). 43 tests.
  Adversarially verified: Louvain gain proven equal to networkx 3.6.1 (sympy) + identical modularity on
  300+ planted graphs, articulation 0/2000 mismatches, fully deterministic. Review caught + fixed one real
  bug — `canonical_kmer_first_pos` used absolute window indices vs python's N-compacted positions (shifts
  the contiguous-span filter on N-containing transcripts).
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


---

## denovo_families.SUPERSEDED

# ⚠ `denovo_families.tsv` is SUPERSEDED — do NOT cite it as the O1 family catalog

`bench/denovo_families.tsv` (1,130 "families") is the **OLD catalog built by the arbitrary similarity threshold
`core_recip ≥ 0.13`** — exactly the "with no arbitrary thresholds" claim the defense audit retired. It
**over-merges**: e.g. `DNFAM0` is a single 728-member family spanning chr1→chrY (unrelated genes bridged by
shared Alu/domain sequence).

**The principled O1 catalog is the threshold-free de-tie READ-CONFLICT-GRAPH catalog**
(`gw_family_catalog → detect_conflict_catalog_genome_wide`, persisted as
`/home/juanfra/winloci_scratch/gw_conflict_catalog.{families,copies}.tsv`): **82 same-chrom families / 207
copies**, no similarity threshold. This is the catalog O2 was run on (P1: 63.9% assigned / 35.7% certified-tied
genome-wide — `bench/P1_P4_RESULTS.md`).

`denovo_families.tsv` is kept only because a few legacy scripts (`copy_resolution_census.py`,
`denovo_families.py`, `dna_psv_catalog.py`, `gw_family_catalog.rs`) still read it for backward comparison. It is
**not** the family-definition result; cite `gw_conflict_catalog` (and, for external precision, the deferred
SEDEF segdup map — `bench/SEDEF_BUILD.md`).


---

## twopass_denovo

# Annotation-free two-pass: read-coherence preserves per-read PSVs → recovers copies flow loses

Prototype of the architecture: **don't use StringTie flow** (which collapses reads and discards the
per-read PSV linkage); use **read-coherence / direct transfrags** as Pass 1 so per-read PSVs survive to
Pass 2, where families are identified and reads assigned to copies. Demonstrated on the *collapsed*
5-copy regime (all copies' reads piled onto a single-copy reference), where the contrast is starkest.
No annotation, no StringTie. `bench/twopass_denovo.py`.

- **PASS 1 (read-coherence):** map all 5 copies' reads to the single-copy reference, group reads by
  exact intron chain → transcript skeletons, keeping each read's sequence. → **ONE skeleton** (the
  copies share the chain), with all reads + their PSVs.
- **PASS 2 (no annotation):** call PSVs *de novo* from the skeleton's collapsed pileup → split reads by
  PSV allele-vector into copies (copy_split) → declare the multi-copy family → assign each read to a copy.

## Result (PSV ladder K = identifiability axis)
| K | Pass-1 skeletons | de-novo PSVs | copies recovered (of 5) | read→copy assignment acc |
|---|---|---|---|---|
| 0 (identical) | 1 | 0 | 0 | 0.0 — UNASSIGNABLE (no info) |
| 1 | 1 | 1 | 4 | 0.80 (1 column < 5 copies) |
| 2 | 1 | 2 | **5** | **1.00** |
| 4 / 8 | 1 | 2 | 5 | 0.99 |

(de-novo PSVs saturate at 2 = ⌈log₄5⌉, the true # of distinguishing columns for 5 copies.)

## What this shows
- **The decisive contrast:** Pass 1 yields ONE skeleton at every K — exactly where **StringTie flow
  would stop and emit ~1 transcript, losing the 5 copies entirely** (their differences collapsed away).
  Read-coherence keeps the skeleton *with its reads*, so Pass 2 can call PSVs and **recover the 5 copies**.
- **The pipeline is annotation-free end-to-end** — skeletons and PSVs come from the reads, not a GTF.
- **It hits the same identifiability boundary** as the annotation-anchored version: 5/5 copies at K≥2,
  partial at K=1, impossible at K=0 (identical) — confirming the substrate change doesn't move the ceiling,
  it just *preserves the signal* needed to reach it.

## Honest scope
- This is the COLLAPSED regime (the copy_split target) — the case where read-coherence's advantage over
  flow is clearest. For dispersed copies (separate loci) the coordinate already resolves them and flow's
  collapse is less harmful.
- Identifiability ceiling unchanged: a read spanning no distinguishing PSV (tied) stays unassignable.
- Prototype on one synthetic family with truth; a production read-coherence Pass 1 would also handle
  noise (chimeras, minor splice variants) that this clean simulation doesn't stress.

## Verdict
Confirms the architectural point: **read-coherence (not flow) is the correct Pass 1 for a copy-aware,
annotation-free transcript caller** — it preserves the per-read PSVs that let Pass 2 recover and assign
copies that flow would have collapsed. The two-pass `(read-coherence skeleton) → (family + PSV split)`
is `copy_split`'s `(intron chain ⊕ PSV)` realized as a full pipeline.

## Genome-wide: not restricted to a handful (bench/twopass_genomewide.py)
The synthetic prototype above is one locus only because it needs a constructed collapsed reference +
ground truth. The pipeline itself is genome-wide — every stage already ran at genome scale this session
(read-coherence: 25 chroms; family graph definition: 1,337 families; de-novo PSV / co-segregation:
10,178 loci in the hidden-collapse scan). Composed over ALL graph-defined families:

- **1,337 families processed genome-wide** (603,267 reads at family loci) in ~5 s.
- **848 (63%) dispersed** — copies at separate loci → coordinate-resolved (Pass-1 separate skeletons).
- **489 (37%) co-located** — copies share/near one frame → de-novo PSV split applies.
- **hard multimappers (MAPQ-0): 3,803 (0.6%)**; **46 co-located families have ≥5** — those are where the
  PSV-split is genuinely decisive.

So: genome-wide-capable and fast; most families resolve by coordinate, the PSV-split is decisive at the
~46 co-located+hard families (the sparse hard regime — consistent with every prior finding; abundant
only in deep co-located data like testis HiFi).

## Reproduce
- `MINIFORGE python bench/twopass_denovo.py` ; `python3 bench/twopass_fig.py`  (synthetic, with truth)
- `MINIFORGE python bench/twopass_genomewide.py`  (genome-wide tally over all families)


---

## integrate_end_to_end

# Integrated end-to-end: identify families → assign reads to copies (hard multimappers)

Runs the full pipeline and answers, honestly: can the method identify multi-copy families accurately
AND assign reads to copies — especially hard multimappers (MAPQ-0)? Validated WITH ground truth on the
synthetic 5-copy benchmark; censused on real GGO (where there is no per-read truth).

## End-to-end on the 5-copy benchmark (WITH truth) — bench/integrate_end_to_end.py
For each divergence level K: (1) POA variation graph of the 5 copies → graph core-score (family
detected?); (2) PSVs recovered from the graph's variant columns; (3) the MAPQ-0 *hard multimappers*
assigned by their PSV allele-vector → accuracy vs truth.

| K | family score | detected | PSVs from graph | hard MM (MAPQ-0) | PSV acc on hard MM |
|---|---|---|---|---|---|
| 0 (identical) | 1.00 | ✓ | 0 | 200 | 0 — UNASSIGNABLE (no info) |
| 1 | 1.00 | ✓ | 1 | 80 | 0.50 (1 column < 5 copies) |
| 2 | 1.00 | ✓ | 2 | 0 | — (minimap2 resolved) |
| 4 / 8 | 1.00 | ✓ | 2* | 0 | — (minimap2 resolved) |

- **Family detection works at every K** — the 5 near-identical copies are correctly called one family
  (score 1.0 ≫ T), and the PSVs are recovered from the graph (the true # of *varying* columns; the
  base-4 design saturates at ⌈log₄5⌉=2 distinguishing columns, so "K=4/8" add no new variation — the
  graph correctly reports 2).
- **Hard multimappers**: at K=0 (identical) they are information-theoretically unassignable (both the
  aligner and PSV); at K=1 PSV resolves 50% of the MAPQ-0 reads the aligner left tied; at K≥2 there are
  **no** MAPQ-0 reads because **minimap2 itself resolves copies once ≥2 PSVs are spanned by full-length
  reads**.

**The honest catch:** for full-length HiFi reads the aligner already resolves non-marginal cases, so
PSV-assignment's *distinctive* edge over minimap2 on hard multimappers is the **marginal (single-PSV)
regime** — plus it provides a principled, identifiability-bounded assignment and an explicit
"unassignable" verdict where the aligner just emits MAPQ-0.

## Real GGO census — how much hard-multimapper signal exists
Over the graph-defined multi-copy families (sampled 400; 174,459 primary reads at family loci):
- **hard multimappers (MAPQ-0) = 0.7%** of family reads (1,156); **6% of families** have ≥5; a few
  carry most (FAM156: 358).

So the hard-multimapper regime the method most distinctively addresses is **sparse in GGO** — consistent
with every prior finding (most paralogs are coordinate-separated; collapsed/co-located real copies are
rare). The 26 hard families are where it applies; there is no per-read truth on real data, so real-GGO
is a *demonstration* + census, and the accuracy proof comes from the synthetic benchmark.

## Bottom line (answering the question)
- **Identify multi-copy families accurately: yes** — the per-family POA-graph definition (validated,
  fixes over-merge) detects them, including on real data.
- **Assign reads to copies, incl. hard multimappers: yes, up to the identifiability boundary** —
  validated with truth on the synthetic hard case; it resolves the marginal hard multimappers PSVs can
  resolve and provably declares the rest unassignable.
- **Caveat:** real GGO under-exercises the hard regime (0.7%), and for full-length reads the aligner
  already handles ≥2-PSV cases — so the method's biggest distinctive payoff is in a **deep co-located
  dataset** (testis HiFi for DAZ/RBMY) where near-identical copies and hard multimappers are abundant.

## Reproduce
- `MINIFORGE python bench/integrate_end_to_end.py` ; `python3 bench/integrate_fig.py`


---

## intronchain_discovery

# Minimizer-free copy discovery via intron-chain alignment (structural axis)

**Question (user):** the family identification uses minimizers, but we're not restricted to them — a
full alignment of intron chains should also work. Try it without minimizers.

**Answer:** intron-chain alignment is a viable, discriminative, **minimizer-free** family criterion.
It recovers the flagship cross-chromosome dup, rejects domain-sharers by construction, and at high
precision is **94% concordant with the sequence pipeline** — i.e. a strong *independent confirmation*
axis. It is not a standalone recall expander in this data (it is blind to retrocopies, which were a
large fraction of the sequence finds).

## Method (bench/extract_intron_chains.py + family_intronchain_discovery.py)
1. Per gene, the **intron chain** = ordered units `(exon_len, intron_after)`, multi-exon only
   (20,502 genes; 2,565 single-exon are the explicit retrocopy blind spot).
2. **Candidate generation (no sequence minimizers)** — two modes (a documented tradeoff, below).
3. **Full Needleman-Wunsch alignment of the two intron chains.** A unit matches iff exon lengths
   agree (±max(6bp,10%)) AND the flanking intron lengths agree (±max(40bp,20%)). Gaps model intron
   gain/loss / exon skip. Structural score = matched units / shorter chain; **gate: score≥0.6 AND
   ≥4 matched units** (RABL2 matches 6; few-exon and chance matches fall below 4).

**Coupling exon+intron is essential.** Exon lengths alone are NOT discriminative (internal exons
cluster ~120–150 bp) → matching them gave 173,283 chance "copies". Adding intron lengths (which span
4–5 orders of magnitude and are preserved between copies) collapsed that to a precise set.

## Recall — the criterion works
- **RABL2A (9 exons, NC_073235.2) ↔ RABL2B (10 exons, NC_086018.1): matched 6, score 0.67** —
  recovered despite the 9-vs-10-exon intron gain/loss (the NW gap absorbs it).
- Discriminative against domain-sharers: the universe "family" LOC101144552 (19 vs 4 vs 16 exons) is
  **correctly rejected** (score 0.06–0.25) — structure is naturally hostile to the FP modes (domain
  share, shared transposon) that needed extra filtering in the sequence pipeline.

## The headline: independent confirmation (structure ⟂ sequence)
At the high-precision gate (matched≥4, length-window): **360 cross-chromosome structural copies,
340 (94%) of which the POA-sequence pipeline also found.** Two fully independent axes — intron-chain
structure vs sequence alignment — agree on the same copies. A pair confirmed by BOTH is robust;
this is a provable, minimizer-free second line of evidence (advisor-aligned).

| axis overlap (cross-chrom) | count |
|---|---|
| confirmed by BOTH structure + sequence | 340 |
| structure-only (sequence missed) | 20 |
| sequence-only (structure blind) | 7,964 |

The 7,964 sequence-only are dominated by **retrocopies / single-exon** copies (no intron chain) plus
copies whose structure diverged past tolerance — the structural axis cannot see these by construction.

## Candidate-generation tradeoff (a finding; both modes kept, `--cand`)
| mode | candidate pairs | cross-chrom confirmed | recall (univ 5) | structure-only | runtime |
|---|---|---|---|---|---|
| length-window, matched≥4 (default) | 5.98M | 360 | 1/5 | 20 | ~5 min |
| length-window, matched≥3 | 5.98M | 993 | 1/5 | 527 | ~5 min |
| 2-intron-shingle, matched≥3 | 14.5M | 38,588 | 3/5 | 37,138 | ~50 min |

- A **single intron length is not discriminative** (common sizes shared by thousands of genes → an
  intron-length index over-collapses). **2-intron shingles** (consecutive intron-size pairs) are
  specific enough to index — they recover more (partial / large-indel copies the length-window
  excludes, e.g. LOC129529434's 5-exon↔10-exon contained match) but are broad → need a tighter NW
  gate and ~10× the compute. The length-window's precision partly comes from its length constraint,
  which is also why it misses partial copies — an inherent tradeoff.

## Honest limitations
- **Retrocopies / single-exon invisible** (no intron chain) — the structural axis's hard blind spot;
  the sequence axis covers them. The two are complementary, neither subsumes the other.
- **Few-exon genes** (≤3 exons) are structurally under-determined (can't reach matched≥4); universe
  2-exon families (LOC101127159, LOC129529611) are missed for this reason.
- **Universe recall 1/5** at high precision, but the misses decompose as: 2 few-exon (under-determined),
  1 structurally-different (correctly rejected = precision win), 1 partial-copy (candidate-gen excluded);
  RABL2 (the one structurally-tractable non-partial case) is recovered.
- **Input = RefSeq gene reps**, one representative isoform per gene.

## Verdict
Minimizer-free intron-chain alignment is a **viable, discriminative, independent structural axis**,
best used as a **confirmation layer** alongside sequence (94% concordant; a both-axes copy is robust),
not a standalone recall expander here. Strengths the sequence axis lacks: domain-sharer resistance and
sequence-divergence robustness. Blind spot the sequence axis covers: retrocopies. **Best system = both
axes.** Next: the user's planned use is exactly this — structure as a second, minimizer-free signal.

## Reproduce
- `python3 bench/extract_intron_chains.py` (gene chains → /tmp/gene_chains.tsv)
- `python3 bench/family_intronchain_discovery.py --cand length` (default; `--cand shingle` for recall mode)


---

## readcoherence_finding

# Read-coherence + degradation-aware extraction — finding (2026-06-14)

**Status: DOCUMENTED, PARKED.** Strong empirical result, but parked in favor of the
multimapping/VG thesis direction (the advisor's interest). Resume from here.

## Result (genome-wide GGO, guided `-G stringtie.gtf`, SQANTI3-validated)

Read-coherence = `rustle -G st.gtf --read-chain` (per-molecule path extraction), compared
additively against the `-G` flow baseline (`gd`). The realness verdict (SQANTI3, paralog
universe + structural categories + canonical junctions + RT-switch):

| metric | value |
|--------|-------|
| additive multi-exon extras over flow baseline | **10,144** (drops only 85 flow chains) |
| canonical junctions | 9,650 / 10,144 (**95%**) |
| **STRICT real** (FSM/NIC + canonical + non-RT-switch) | **1,857 (18%)** |
| REAL (+ NNC canonical) | 4,685 (46%) |
| of real extras, cov≥5 | 3,607 / 4,685 (**77%** — well-supported) |

Of the 1,857 strict-real: **503 FSM** (exact RefSeq matches *both* StringTie and rustle-flow
missed — indisputable recall) + ~1,354 NIC (real novel combinations of known junctions).

**Contrast with VG-multimapper (same harness, guided): 0** FSM paralog-copy recoveries over
baseline. Read-coherence is the bigger recall lever by ~3 orders of magnitude — but it is
**not** the advisor's interest (multimapping is), hence parked.

## VALIDATED vs StringTie (2026-06-15, rigorous cut — survives the PSV-grade scrutiny)

gffcompare rc_$C vs RefSeq minus st_$C vs RefSeq, genome-wide (25 chroms):
- StringTie FSM(=)=23,371; read-coherence FSM(=)=23,951.
- **FSM (exact) read-coherence finds that StringTie MISSES: 580** (≈ the parked 503 estimate, confirmed).
- broad (=/c/j) StringTie misses: **2,784**.
- FSM = exact annotated-transcript match ⇒ canonical+real by construction (NOT the PSV non-canonical-
  artifact failure mode). This is the real beat-StringTie lever (~100–200× the PSV margin of 2–3 exact).
- COST: the 580 sit inside ~10,144 raw extras (~half noise) ⇒ need the annotation-free realness gate.

## SHIPPED — gated `--read-coherence` layer VALIDATED (2026-06-15, genome-wide 25 chroms)

The gated layer (degrade-collapse 3′+internal + annotation-free realness gate + holdout/union-back
additivity) is built, wired (default-off byte-identical, 422 lib tests), and validated. Rigorous
FSM(=) cut, `rcg`/`rc`/`gd`/`st` vs RefSeq (`bench/rcg_gen.sh` + `bench/rcg_validate.sh`):

| arm | FSM(=) total | beat-ST cut (\\ st_FSM) | total tx |
|---|---|---|---|
| StringTie (`st`) | 23,372 | — | — |
| flow baseline (`gd`) | 23,549 | 177 | — |
| ungated read-chain (`rc`) | 23,952 | **580** (reproduces the documented baseline ✓) | 78,559 |
| **gated read-coherence (`rcg`)** | **24,373** | **1001** | **76,680** |

- `st ⊆ rcg` exactly (guides preserved) ⇒ the 1001 is **pure additive recall over StringTie**.
- **Degrade-collapse is the hero:** gated finds +589 *new* FSM over ungated (folding degraded
  fragments back into full-length parents promotes ISM→FSM), at the cost of 168 ungated FSM the
  realness gate over-kills (71% retention). Net +421 over ungated, with ~1,900 FEWER total tx.
- **CRITICAL-1 cost (flow replacement):** `--read-coherence` REPLACES flow extraction per-bundle
  (pipeline.rs:7283 `RUSTLE_READCHAIN` else-if), so it gives up 96 flow-only FSM (all ST-misses).
  Truly-additive ceiling (`rcg ∪ flow`) = **1097** (+96). Option B = produce flow AND read-chain
  per bundle + merge (stronger "⊇ flow baseline" claim for the adversarial advisor).
- **Gate over-kill (168) partly fixable:** `is_canonical_junction` only accepts GT-AG/CT-AC, missing
  the minor-canonical GC-AG / AT-AC (U12) classes → some of the 168 are real GC-AG isoforms.

## OPTION B (truly-additive) + GC-AG fix — VALIDATED (2026-06-15, user chose "max")

Re-architected the layer from read-chain-REPLACES-flow to truly-additive `flow ∪ read-chain`:
per-bundle dual extraction (read-chain on a transfrag clone so flow sees pristine input) +
a per-bundle hold-out (`rc_layer_buf`) + a re-introduced global hold-out/union-back
(`read_coherence_holdout`, mirrors `union_baseline_holdout`) so read-chain NEVER displaces a
flow/guide transcript at any filter stage. Also fixed `is_canonical_junction` to accept the
minor canonical classes GC-AG / AT-AC (+ revcomps), recovering real isoforms the gate over-killed.

Genome-wide rigorous FSM(=) cut (25 chroms):

| arm | FSM(=) | beat-ST (\\ st) | total tx |
|---|---|---|---|
| StringTie | 23,372 | — | 68,157 |
| flow (gd) | 23,549 | 177 | 70,373 |
| ungated read-chain (rc) | 23,952 | 580 | 78,559 |
| replacement gated (prev) | — | 1,001 | 76,680 |
| **OPTION B (rcg)** | **25,620** | **2,248** | **94,145** |

- **+2,248 exact annotated isoforms StringTie misses** (~4× ungated, ~2.2× replacement). The jump:
  read-chain now bypasses the pairwise/isofrac/predcluster filters (gated only by canonical +
  RT-switch + read-depth), which were culling real FSM in the replacement path.
- **⊇ flow guarantee holds genome-wide: `gd \ rcg = 0` AND `st \ rcg = 0`** (never loses a flow
  find or a StringTie guide) — the provable additivity claim for the adversarial advisor.
- **Precision cost:** 94,145 tx vs StringTie 68,157 (+38%); = 70,373 flow floor + 23,791 gated
  read-chain extras. The gate `min_cov` (default 2, env RUSTLE_READ_COHERENCE_MIN_COV) is the
  tunable precision lever; raising it trades read-chain recall for fewer transcripts.
- Default-off BYTE-IDENTICAL (chr19 off==gd exactly); RUSTLE_PRECISE airtight (BOTH read-chain
  extraction branches precise-gated; precise+rc-env == precise-plain, diff 0); 424 lib tests.

### Precision of the 23,791 read-chain extras (SQANTI3, genome-wide — answering "is the +38% real?")
| structural_category | count | % |
|---|---|---|
| novel_not_in_catalog (NNC) | 10,059 | 42% |
| novel_in_catalog (NIC) | 6,275 | 26% |
| incomplete-splice_match (ISM) | 2,427 | 10% |
| full-splice_match (FSM) | 2,066 | 8% |
| intergenic | 1,134 | 4% |
| fusion | 949 | 4% |
| antisense | 709 | 3% |
| genic_intron / genic | 172 | 1% |

- **100% canonical** (all 23,791 — the gate guarantees it). RT-switch (SQANTI3 RTS): 3,971 (17%).
- **REAL (FSM/NIC/NNC + canonical + non-RT) = 14,999 (63%)**; STRICT (FSM/NIC + …) = 6,731 (28%).
- FSM+NIC+NNC = 18,400 (**77%** are matched/novel isoforms); likely-noise (intergenic+antisense+
  fusion+genic_intron) = 2,909 (**12%**); ISM/genic (degradation/partial) = 2,482 (10%).
- The gate IMPROVED composition vs the ungated 10,144: real 46%→63%, strict 18%→28%, canonical 95%→100%.
- Residual precision lever: rustle's narrow 8bp RT-switch heuristic misses the 3,971 SQANTI3 calls;
  a wider/two-orientation RTS detector (or higher gate min_cov) would shed those.

## Noise anatomy (the ~5,459 non-real extras) — filter-vs-generate split

| artifact | count | fix path |
|---|---|---|
| ISM truncations (3′ 794 / 5′ 698 / internal 636 / IR 95) | 2,223 | **better generation** — degradation-aware collapse (fold fragments into full-length parent). rustle's existing 5′-degrade collapse (`RUSTLE_READCHAIN_DEGRADE_TOL`, default-on) is incomplete: leaves 3′ + internal + residual 5′. |
| RT-switching | 1,538 | **filter** (data artifact; the false junction is genuinely in the read — annotation-free detectable via genome sequence around the junction) |
| non-canonical junctions | 494 | **filter** (canonical GT-AG check, annotation-free) |
| intergenic / antisense (~55% single-read) | 1,331 | read-depth + locus-context gate |
| fusion (read-through) | 510 | read-through logic |

**Verdict:** hybrid. ~2,128 (ISM degradation) avoidable by better *generation* (the
"degraded-transcript modelling" idea); ~2,032 (RT-switch + non-canonical) are *data* artifacts
needing a *filter* (no extractor avoids them). Neither alone suffices.

## Proposed design (when resumed)

Shippable #1 = **degradation-aware read-coherent extraction** (extend the degrade-collapse to
3′ + internal fragments; the biggest, cleanest noise win, doubles as the 3′-degradation
feature) **+ an annotation-free realness gate** (canonical junctions + RT-switch detection +
read-depth). Additive over the byte-exact `-G` flow floor, opt-in/gated like the other
"better decisions", `RUSTLE_PRECISE`-exempt. SQANTI3 stays the validator each iteration.
Expected yield: 10,144 → ~4,000–4,700 real isoforms genome-wide over the StringTie-exact floor.

## Reproduce
- generate: `bench/rc_gen.sh` (per-chrom `rustle -G st --read-chain`, OOM-safe) → `/tmp/gw/rc_*.gtf`
- validate: `/tmp/cre_guided/rc_validate.sh` (merge + SQANTI3 + realness verdict)
- harness reuses the copy-recovery SQANTI3 stages + cached paralog universe.


---

## readcoherence_psv_headroom

# Read-coherence × PSV — HEADROOM probe (go/no-go)

Cheap geometric+structural probe: do read-coherence's recall wins and PSV copy-resolution
act on the SAME loci (→ unified molecule-threaded graph worth building), or disjoint loci
(→ not worth the rebuild; keep them as two separate levers)?

## Inputs
- read-coherence recall set: transcripts tagged `source "read_coherence"` in rcg_*.gtf (25 chroms)
- multi-copy families (universe.tsv, n_copies>=2): **58** families, 556/614 member transcripts located in RefSeq
- family regimes: TANDEM_NEAR=27, DISPERSED=17, COLLAPSED=13, SINGLE_LOCUS=1

## Recall set
- read_coherence transcripts total: **23791**
- multi-exon (the recall lever): **23791**
- of those, FSM (intron chain == a RefSeq transcript = real recall): **1921** (8.1%)

## Tier 1 — are the recall wins even at multi-copy family loci?
- recall isoforms overlapping a multi-copy family locus: **88** / 23791 (**0.37%**)
- of those, FSM: **7** (of 1921 FSM = 0.36% of all read-coherence real recall)

## Tier 2 — of the family-locus hits, which regime (does PSV actually SPLIT?)

| regime | recall isoforms | FSM | PSV leverage |
|---|---|---|---|
| COLLAPSED | 35 | 4 | **YES — copies share one frame, PSVs split** |
| TANDEM_NEAR | 32 | 2 | no — distinct coords separate them |
| DISPERSED | 21 | 1 | no — distinct coords/contigs (RABL2 boundary) |

## Tier 3 — of the COLLAPSED hits, are the 'copies' REAL paralogs or domain-sharers?
(PSVs can only copy-resolve genuine paralogous copies. Two different genes that merely share a
protein domain/repeat and happen to overlap by annotation are NOT copies — PSVs are meaningless there.)

- COLLAPSED recall isoforms at CONFIRMED real paralog families (APOBEC3/RFPL/RABL2): **0**
- COLLAPSED recall isoforms at DOMAIN-SHARER false families (CDPF1/PPARA, CREB1/METTL21A, GCA/KCNH7, TTN/NCAPH2/MIEF1-LOC): **35**
- COLLAPSED recall isoforms at unclassified families: **0**

- read-coherence recall isoforms landing on ANY confirmed-real paralog family (any regime): **3**

## Verdict
Geometric PSV-resolvable headroom (recall isoforms at COLLAPSED loci) = **35** (0.147% of the recall set, 4 of them FSM).
TRUE headroom after removing domain-sharer false families = **0** (confirmed-real collapsed paralog recall isoforms = 0).

**GO/NO-GO: NO-GO** for a molecule-threaded graph justified by *combining* recall + copy-resolution.
Read-coherence recall (99.6% at single-copy loci) and PSV copy-resolution (needs collapsed real
paralogs) are disjoint in this data. The recall lever already ships additively (`--read-coherence`);
threading PSVs through the same graph adds copy-resolution to ~0 real recall isoforms.

Top families receiving recall isoforms:
  - LOC134756662 (COLLAPSED, 2 copies): 14 isoforms [DOMAIN-SHARER]
  - CREB1 (COLLAPSED, 2 copies): 7 isoforms [DOMAIN-SHARER]
  - LOC109024534 (TANDEM_NEAR, 2 copies): 6 isoforms
  - CDPF1 (COLLAPSED, 2 copies): 5 isoforms [DOMAIN-SHARER]
  - IGLL1 (DISPERSED, 2 copies): 5 isoforms
  - LOC101126655 (TANDEM_NEAR, 11 copies): 5 isoforms
  - GCA (COLLAPSED, 2 copies): 4 isoforms [DOMAIN-SHARER]
  - GGT1 (DISPERSED, 2 copies): 4 isoforms
  - LOC101136027 (TANDEM_NEAR, 2 copies): 4 isoforms
  - LOC129529456 (COLLAPSED, 2 copies): 4 isoforms [DOMAIN-SHARER]
  - CCDC188 (DISPERSED, 4 copies): 3 isoforms
  - LOC101129569 (TANDEM_NEAR, 10 copies): 3 isoforms
  - LOC101134642 (DISPERSED, 2 copies): 3 isoforms
  - LOC115931965 (TANDEM_NEAR, 2 copies): 3 isoforms
  - LOC129529434 (DISPERSED, 4 copies): 3 isoforms

![funnel](readcoherence_psv_headroom.png)

## What this does and does NOT settle
- **Settles:** the *PSV-unification* motivation. Read-coherence's recall win does not co-occur
  with PSV-resolvable copies, so building one graph to deliver BOTH on the same loci pays ~nothing.
- **Does NOT settle:** the molecule-threaded graph as a *pure recall* architecture (gate-not-kill,
  keep per-molecule evidence through paths). That stands or falls on recall alone — and the recall
  is *already* delivered additively by the shipped `--read-coherence` layer, so a rebuild would buy
  cleaner architecture / less noise, not new recall. It would NOT add copy-resolution.

## Honest caveats (which way they push)
1. **Geometric undercount (pushes headroom UP a little):** COLLAPSED is detected only when *annotated*
   copies overlap. Cross-mapping can pile a divergent/unannotated copy's reads onto a *single-copy*
   annotated locus, creating PSV signal the annotation hides. But prior copy_split real-data work
   (RABL2 / DAZ / co-located tests) showed that regime is rare and coverage-limited in GGO; even a
   3–5× undercount keeps the headroom sub-1% and the *real-paralog* count near zero.
2. **Small family universe (neutral):** 58 multi-copy families is the genome's reality, not a sampling
   gap. Multi-copy paralogs are a small slice of loci; read-coherence's recall is dominated by
   single-copy alt-splicing/novel-junction isoforms regardless of family-set size.
3. **Domain-sharer contamination (pushes headroom DOWN, decisively):** all 35 COLLAPSED hits land on
   families the Compara + contiguous-core validation already flagged as domain-sharing FALSE families
   (two different genes sharing a repeat, overlapping by annotation accident) — not paralog copies.
   PSVs cannot copy-resolve non-paralogs, so the TRUE headroom is 0, not 35.

## Reproduce
- `python3 bench/readcoherence_psv_headroom.py`  (reads /tmp/gw/rcg_*.gtf + ref_*.gff3 + universe.tsv)
- `python3 bench/readcoherence_psv_headroom_fig.py`  (renders the funnel figure)
- COLLAPSED-regime hit loci: `bench/headroom_loci.tsv`

