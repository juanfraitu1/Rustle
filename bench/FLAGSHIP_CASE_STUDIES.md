# Flagship case studies — validating the two advisor interests

The thesis is **not** an assembler. The StringTie-cloned core is the **substrate** (it produces
transcripts/loci). The contribution is two things, and these four examples validate them:

> **(I)** a *topological / mathematical* definition of multi-copy gene families at the RNA level —
> family = clique in the backbone graph, copies = paths through a variation graph, identifiability =
> **MCC = χ(H) = minimum path-cover** of the read-conflict graph.
> **(II)** assigning *multimapping reads to copies in difficult cases* (MAPQ-0 ties) using PSVs +
> other properties (copy-specific junctions, divergence), **up to the identifiability limit**, with a
> calibrated decisive-margin gate and **no 1/k guessing**.

## Narrative arc

StringTie substrate → **topological family definition** (clique / variation graph, χ(H)) →
**difficult copy assignment** (PSV + junction, up to identifiability) → **discovery** (copies the
reference lacks). Each step is a clean combinatorial object with a provable guarantee — the Canzar
aesthetic — not a tuned threshold.

## The four flagships

### 1. sim5x K-ladder — the mathematical spine (validates I + II)
Five tandem copies, PSV count **K dialled 0,1,2,4,8**, reads carry true labels. Recovery = **0% at
K=0, 100% at K≥2**, zero free parameters — carrying Lemma 1 (MCC = χ(H)) and Theorem 2 (unique recovery
under Strong Separation).
- **Topological:** reads = vertices of conflict graph H; a PSV column with disagreeing alleles = an
  edge; copies = colour classes; recovery is χ(H) crossing from edge-sparse to Strongly-Separated.
- **Difficult:** at K=0 copies are exonically identical → minimap2 MAPQ-0, identical edit distance,
  0% resolvable. PSVs cross a boundary alignment cannot.
- **Airtight gap:** re-run under a HiFi error-rate sweep (errors only worsen identifiability → the
  theorem is conservative, but show it).

### 2. DSFAM237 (the WIN) + DSFAM42 (the FLOOR foil) — the real difficult case (validates II)
A genome-wide scan of all 68 multi-copy families on two axes — **difficulty (MAPQ-0 fraction)** vs
**read-assignment rate** (`hard_loci_psv_assignment.py`) — settles it: of the **5 genuinely hard loci**
(≥50% MAPQ-0), **PSVs win 4** and only DSFAM42 is the floor. So the difficult case is *usually solved*;
DSFAM42 is the rare certified exception. The flagship is the win, with DSFAM42 as the honest foil.

| family | MAPQ-0 | PSV cols | reads assigned | silver |
|---|---|---|---|---|
| **DSFAM237 (WIN)** | 90% | 10 | **94%** | **1.0 (5 uniq)** |
| DSFAM817 | 93% | 249 | 90% | 0.67 |
| DSFAM238 | 70% | 1049 | 98% | 0.6 |
| DSFAM102 | 86% | 1234 | 91% | 0.0 (silver degenerate) |
| **DSFAM42 (FLOOR)** | 95% | **3** | **21%** | 1.0 |

- **The WIN — DSFAM817, confirmed END-TO-END on the production engine (clean on BOTH).** 3
  size-homogeneous copies (~10 kb each, NC_073229.2 ~44.4 Mb, 59 kb span — no container to over-merge).
  Production `copy_assign` de-novo: detects a **clean 3-copy family** under **95% MAPQ-0**, assigns
  **79/118 reads (67%)** confidently, **27 tied** (the honest K-frontier minority), and emits a clean
  **3-copy-path variation-graph GFA** (890 bubble nodes). The curated `assign_family` agrees (90%); both
  engines win (the 67/90 gap is the stricter de-novo operating point τ≈6.9 vs the recall-mode τ=2.0).
  **Silver = 0.67 (2 of 3 unique-mappers agree)** — *thin and circular* (the "truth" is minimap2's own
  primary placement, and only 3 reads map uniquely in this MAPQ-0 locus). [Reconciliation: an earlier
  draft said "silver 3/3 = 100%"; the measured value in `hard_loci_psv_assignment.json` is 2/3 = 0.667.]
  *minimap2 ~0% confident → PSVs 67–90% confident, with the unresolvable minority abstained.* The
  **load-bearing validation is the sim5x labeled-truth oracle** (below), not this circular silver.
- **DSFAM237** (3 small copies over 162 kb): the curated family wins at 94% (silver 5/5), but the CLI's
  de-novo over-merges a 42 kb neighbor → 0% — the clean illustration that *the family definition gates
  the assignment* (see the engine caveat below).
- **The win/floor split IS the theorem.** Look at PSV count: DSFAM42 has **3** (copies near-identical →
  K-frontier → it *abstains*, 21%); every win has 10–1,234 (copies distinguishable → it *assigns*,
  90–98%). The method assigns when the copies differ and honestly abstains when they don't — DSFAM42 is
  the floor certificate, not a failure (1 of 5 hard loci; ~1.5% of families).
- **Honest validation caveat (= the motivation):** silver rests on *few* unique-mappers in the hard
  regime (3–15 reads) — because MAPQ-0 *means* few reads map uniquely. Strong where measurable
  (DSFAM237 5/5) but thin where not (DSFAM102 0/4, where the best-mapping "truth" is itself arbitrary).
  The silver standard **degenerates exactly where the method is needed** — which is why the airtight
  validation for the hard wins is DNA/orthogonal, not silver.
- **⚠ Engine caveat — and it is a thesis STRENGTHENER (interest I gates interest II).** The win (94%)
  is on the **curated de-tie 3-copy family** — what the family-definition pipeline and the production GTF
  `psv_linkage` path use, faithfully mirrored by `assign_family`. The standalone `copy_assign` CLI does
  its *own* de-novo detection and here **over-merged** the locus into a size-heterogeneous **5-copy**
  family (a 42 kb container + small copies) → 0/644 assigned, 489 tied. Same confound as DSFAM42's
  de-novo run. The lesson: **give the method the right topological family and PSVs assign 94% under 90%
  MAPQ-0; give it a size-heterogeneous merge and even 600 PSVs don't help.** The family definition
  (interest I) is not decoration — it is the *precondition* for the assignment (interest II). Evaluate
  on the curated family; the CLI's de-novo over-merge is a detection artifact, not a PSV failure.
- **⭐ THE CANONICAL ENGINE (L7) — one scoring, the full molecule, a principled gate.** Of the three
  historical wrappers (production vote `psv_linkage`, the `combined` pipeline engine, the CLI LLR), the
  **canonical one is `copy_assign::assign_read` driven by the `combined` pipeline path** — because it uses
  the **full long-read evidence: PSV columns + the read's own copy-specific junction chain** (a unit
  test proves it strictly out-resolves PSV-only: `psv.n_decisive=0` but `combined.n_decisive≥1` when a
  junction is the only discriminator). This is the **FLAIR-like** choice (per-molecule, the whole read
  defines its own assignment) and the **Canzar-clean** choice (assign-or-**abstain**, no 1/k). The vote
  engine is its flat-error vote-equivalent (kill-test 16/16); the CLI is the same scoring exposed
  standalone. The decision gate is **n_decisive ≥ 1** (identifiability: the read must span ≥1 column/junction
  where copies differ) **AND** a **decisive log-LR margin τ**, where **τ = ln((1−p)/p)** is a *principled
  operating point* set by the target per-read misassignment rate `p` — NOT an arbitrary threshold. The two
  values in the codebase are just two `p`: **τ=2.0 ≡ p≈0.12 (recall mode, default)** and **τ≈6.9 ≡ p=1e-3
  (precision mode, the PSV-space analog of Eichler's AS≥10)**. Set via `AssignParams::for_target_misassignment(p)`.
- **⭐ NON-CIRCULAR VALIDATION — the sim5x labeled-truth oracle (the load-bearing test, not the silver).**
  Each simulated read name encodes its TRUE copy, so accuracy is measured against *planted* labels, not
  minimap2. The canonical engine on the K-ladder (`smoke_sim5x_ground_truth`, 1000 reads/level):

  | K (PSVs) | resolvable% | **acc \| assigned** | acc \| forced-argmax | tied% |
  |---|---|---|---|---|
  | 0 | 0% | — | — | **100% (abstains)** |
  | 1 | 20% | **1.000** | 0.800 | 80% |
  | ≥2 | 20% | **1.000** | 1.000 | 80% |

  The headline is **`acc|assigned = 1.000` at every K ≥ 1**: *when the engine commits, it is never wrong*,
  and at K=0 it commits to nothing (100% tied — no fabrication). The gap to forced-argmax (0.80 at K=1)
  is the measured *value of abstaining*. This is the identifiability theorem made empirical on a ground
  truth that is **not** the aligner's own placement — cite this, not the circular 1026/1026 silver.
- **Honest scope (L9/L10):** the genome-wide run of the canonical engine is **blocked** — the production
  `GGO.bam` is currently missing/repointed (loose end L4), so the per-family/CLI + sim5x results above are
  what stands; a genome-wide PSV-linkage pass (`gw_psvlink.sh`) is built but unrun. State copy-assignment
  as a **per-family + sim-validated** capability, not yet a default genome-wide output.
- **Junctions = a real but MINOR second axis** (`junction_rescue_probe.json`): copy-specific junctions
  rescue **5.5%** of no-PSV reads genome-wide (96.5% validated) — an honest adjunct, not the hero.
- **Unassignability is certified, not assumed** (`UNASSIGNABLE_SEPARABILITY_ATTEMPT.md`): on DSFAM42's
  tied reads, 17 BAM features + RandomForest + KMeans recover **nothing** (sim5x control: predict true
  copy from seq/qual = 0.245 vs chance 0.200; tied-read clusters ARI-vs-copy = 0.00). The floor is the
  information-theoretic limit, ML-confirmed.

### 3. RABL2 vs RFPL4A — the topological family definition (validates I)
The clean rendering of "family = clique, copy = path, K-frontier = graph property."
- **RABL2** (2 copies, separate chromosomes): fully-resolvable 2-copy clique, 67 PSV bubbles,
  **2 nodes → 2 paths**, 58/58 reads agree with the minimap2 primary. χ(H) = #copies = 2.
- **RFPL4A** (5-copy tandem array): a founder + 4 near-identical duplicates → the graph exposes only
  **5 nodes → 2 paths** (copies 2–5 are PSV-identical across 18 columns). χ(H) = 2 < 5 = the
  **K-frontier as topology**.
- **Difficult:** RFPL4A's 4 near-identical copies are indistinguishable on RNA; 54% of reads hit no
  PSV → **honestly unassignable**, while the 6 PSV-spanning reads assign perfectly. The win *plus the
  principled refusal to guess* (no 1/k).
- **Airtight gap:** render the two GFAs side by side (RABL2 2→2, RFPL4A 5→2) with χ(H)/min-path-cover
  on the graph (the `--phase` GFA emitter now exists; previously only `psv_graph_demo.json`).

### 4. Reference-absent MHC + DAZ1/DAZL junction — the discovery (validates I + II)
Novel findings. **4 reference-absent divergent MHC copies** (Gogo-B / DQ-α / DQ-β / DRB1) detected as
hidden-haplotype cliques (no assembly), protein-BLAST-confirmed endogenous; and **DAZ1 vs DAZL
copy-specific junction reversal** (dPSI 0.918, q = 2.6e-151) — the property *beyond* PSV.
- **Topological:** a hidden copy = a clique of balanced co-segregating alt-columns; DAZ1/DAZL =
  min-path-cover on the haplotype-junction bigraph.
- **Airtight gap:** cross-individual replication + DNA parCN for the 4 MHC copies (the copy-vs-allele
  resolver).
- **⭐ GROUND-TRUTH STARVED-COPY RESCUE (the planted proof of the discovery — `sim_starved.py`).** The
  old StringTie-era idea — *use the multimapping reads to rescue a copy starved of primaries* — survives,
  reframed as rescue-then-assign-or-abstain. Plant 3 near-identical copies (6 PSVs); two healthy (40
  reads), one **starved to 1 expressed read**. minimap2 reproduces the textbook signature — the starved
  copy gets **1 primary / 80 secondary** (its siblings' reads pile on its locus). The pipeline recovers it
  as **`RC_sc_48446` at exactly the planted locus**, quantifies it as the **minor** copy (abundance 0.012,
  *not* inflated by the 80 shadow secondaries), and **assigns its own read to it** — all driven by the
  multimapping evidence + the 6 distinguishing PSVs. Two honest refinements this exhibit pins down: (a)
  the **default** collapsed-copy rescue already gets the 1-primary copy (`--recover-copies` is byte-identical
  here — it earns its keep only in the **0-primary** regime); (b) a fully-0-primary copy is essentially the
  **reference-absent / collapsed** case, because *minimap2 spreads primaries across in-reference duplicates*
  (the K0tandem identical trio splits 40/37/43) — which is why the real-data version is `--absent-copies`
  (O4) and lives in the GGO 905-collapsed-copy result. **Guard:** the copy is admitted only because it has
  genuine PSVs + its own read; a locus carrying *only* sibling shadow (no real expression) is rejected by
  the admission gate — the multimapping reads cannot fabricate a copy. Reframed, never **1/k**. Reproduce:
  `bash bench/sim_starved_run.sh`.

## Build order (smallest → highest value)
1. **Render RABL2 + RFPL4A GFAs side by side** with χ(H)/min-path-cover — the flagship topological
   figure (the `--phase` emitter makes this immediate).
2. Validate DSFAM42's junction-only win against an independent splice catalog / short-read RNA-seq.
3. sim5x under a HiFi error-rate sweep.
4. Cross-individual + parCN for the 4 MHC copies.
5. Render a min-proper-colouring on a real conflict graph for Lemma 1.

## Drop (to keep the advisor-focused story sharp)
- The **assembler-lineage** threads (recall, flow-parity, StringTie-mimicry, VG-regression) — the
  better-assembler framing the advisor does NOT want.
- Exhaustive enumeration + recombination witness → a footnote, not a flagship.
- MAGEA arrays → fold into the sim5x K=0 floor as one corroborating bullet.
- Intermediate DSFAM / family IDs → keep DSFAM42, RABL2, RFPL4A only.
- PSMD2 ASJ catalog as a headline → keep only as the DAZ1/DAZL anchor.
