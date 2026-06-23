# Bridge: from the family definition to copy assignment

**PSV-on-backbone *are* the columns.** One object — the family's variation graph — is both the detection output and the
resolution input; identifiability is the single thread through both.

This note connects two formal pieces:
- **Detection** — `family_definition_note.md`: a multi-copy gene family = a $\sim_B$-connected component within a
  read-coupled ($R^\*$) group of de-novo loci. Output: the family's member loci **and** their shared **backbone** (the
  $\sim_B$ alignment of the all-isoform exon-union copy models).
- **Resolution** — `copy_assignment_theory.md`: reads as partial allele-vectors over **columns** $[m]$ (PSVs); the
  **conflict graph** $H$; $\mathrm{MCC}(R) = \chi(H)$ (Lemma 1); **Strong Separation** $\Rightarrow$ unique recovery
  (Thm 2); the **K-frontier** (recombination / collapse).

Detection produces exactly what resolution consumes. The map is below.

---

## 1. The handoff

Detection hands resolution a family $C \subseteq V$ and the backbone that $\sim_B$ already computed. Resolution needs
three things from its model (§2 of the theory): the **columns** $[m]$, the **copies** as allele-vectors, and the
**reads** as partial allele-vectors. All three are read off the family's variation graph.

Build the family's graph once: align the copy models $\{S(v) : v \in C\}$ (the same exon-union sequences $\sim_B$
aligned pairwise, now taken multi-way / POA). The graph is a shared **backbone** (consensus) with **bubbles** where the
copies diverge.

---

## 2. The bridge map (three identifications)

> **Columns = bubbles.** A **PSV column** $j \in [m]$ is a backbone position at which the copy models do not all carry
> the same base — a **bubble** of the variation graph — kept only if it is **read-supported** (recurs across reads, so
> it is a paralog difference, not a one-read sequencing error; the error-vs-recurrence test). The allele alphabet
> $A_j$ is the set of bases the copies present at $j$.

> **Copies = paths.** Each copy $v \in C$ realises an allele-vector $c_v = \big((S(v) \text{ aligned})_j\big)_{j\in[m]}$
> — its **path** through the bubbles. This is the theory's "gene copy" $c_v \in \prod_j A_j$.

> **Reads = partial paths.** A read $r$ threaded through the graph observes the bubbles it spans:
> $\mathrm{obs}(r) = \{\,j : r \text{ covers column } j\,\}$, and $r(j)$ is the allele it carries there. This is the
> theory's partial allele-vector.

With these three identifications, every object of `copy_assignment_theory.md` is defined directly from the family
graph, and the theory applies **verbatim**: the conflict graph $H$ on the reads, $\mathrm{MCC} = \chi(H)$ (the copy
number), Strong Separation (Thm 2) for unique recovery, and the K-frontier for the limit.

---

## 3. The two relations, in resolution's terms

The detection relations are not discarded — they are the *preconditions* that make the columns meaningful:

- **$\sim_B$ (shared backbone) $\;\Longrightarrow\;$ the columns exist.** Two loci are $\sim_B$ iff their copy models
  align reciprocally — i.e. they share a backbone on which the bubbles $D_{ij} \subseteq [m]$ (the distinguishing
  columns of the pair, §5 of the theory) are *defined*. Without $\sim_B$ there is no common coordinate to call a
  "column," so resolution is ill-posed. $\sim_B$ is precisely the well-posedness guarantee resolution assumes.
- **$\sim_R$ (read-confusability) $\;\Longrightarrow\;$ the problem is non-trivial.** Two loci are $\sim_R$ iff reads
  cross-map at tied divergence — i.e. the conflict graph $H$ is non-empty and the reads genuinely must be assigned.
  If the loci were not $\sim_R$, the reads would already place uniquely and $H$ would be edge-free.

So **a family ($\sim_R \cap \sim_B$) is exactly the regime where copy assignment is both well-posed (columns exist) and
non-trivial (reads conflict)** — which is the entire point of detecting it first.

---

## 4. One object, two readings — and the identifiability thread

The family's variation graph is read **twice**:
- **Detection** reads its *macro* structure: do the copies form long **parallel paths** (a shared backbone)? That is
  $\sim_B$; it certifies "these loci are copies."
- **Resolution** reads its *micro* structure: do the **bubbles** distinguish the paths, and do reads **link/cover**
  them? That is Strong Separation; it certifies "these reads can be assigned to copies."

Both are governed by the same bubbles $D_{ij}$, so a single quantity — the distinguishing-column structure — controls
the whole arc. This yields the precise containment:

> **Detection $\supseteq$ resolution.** The definition can *detect* a family that resolution *cannot* resolve. The
> family needs only $D_{ij} \neq \varnothing$ for the pair to be a paralog *and* reads to cross-map; resolution needs,
> in addition, that the bubbles are **covered and linked** so no recombinant or collapsed cover competes (Strong
> Separation). The **K-frontier** names exactly the gap:
> - $D_{ij} = \varnothing$ for some true pair $\Rightarrow$ those copies are one path (**collapse**, $K=0$ between them)
>   — detectable as a family, *provably unresolvable* into individual copies;
> - $D_{ij} \neq \varnothing$ but reads do not link/cover them $\Rightarrow$ recombinant covers compete (the $K\!\ge\!3$
>   phasing failure) — resolvable only with longer reads.

So the identifiability theorem is not a separate result bolted on after detection: it is the *same* graph, read at the
bubble level. Detection establishes the backbone on which the K-frontier is even stateable.

---

## 5. Empirical instantiation (`bench/psv_graph_demo.py`)

The two regimes, on real great-ape families, are the two sides of §4:
- **RABL2 (2 copies)** — the bubbles separate the two paths and reads cover them: Strong-Separation-like. **58 reads
  thread to a single copy, 100 % agreeing** with the best-mapping copy; the rest stay on the shared backbone (cover no
  bubble) — the *coverage* face of the K-frontier.
- **RFPL4A array (5 copies)** — 18 bubbles, but copies 2–5 take the **same allele at every bubble** ($D_{ij}=\varnothing$
  within the cluster): a founder (RFPL4A) plus **4 collapsed paths**. The graph *shows* the collapse — 6 reads resolve to
  RFPL4A, 5 only to the cluster (which of the 4 is unidentifiable) — the **$K=0$ collapse** face of the K-frontier.

The variation graph makes $D_{ij}$ visible (the bubbles), so "resolvable vs not" is read directly off the picture — the
identifiability dichotomy as graph structure.

---

## 6. Exact vs. pipeline

The **column abstraction is exact**: given the family graph, the bubbles, copy-paths, and read-paths are defined with no
free parameter, and the theory's results (Lemma 1, Thm 2, the K-frontier) hold as stated. The **pipeline** supplies one
calibration — the read-support threshold that separates a paralog bubble from a one-read error (the recurrence test,
§2-bubbles) — and inherits the de-novo locus vertices from detection. Both are the same disclosed dependencies the
family definition already carries; the bridge adds none.

---

## 7. In the code — the bridge is the production pipeline

Every object above already exists in the Rust path (`src/rustle/vg_family/`), and the copy-assignment driver `run_layer2`
(`layer2.rs:537`) walks exactly this chain:

| formal object (this note) | Rust function | location |
|---|---|---|
| **$\sim_B$ family** (copies share a backbone) | `contiguous_core_coverage` $\ge$ `T_CORE` (= the family criterion) | `family_detect.rs:15,31`; `family_graph.rs:933,1046` |
| **columns = bubbles** (PSVs from the copy backbone) | `psv_columns_from_reference` (POA-aligns the copies' exon sequences, emits a column per divergent position) | `psv_linkage.rs:159` |
| (node-shared PSVs, the co-located case) | `psv_columns_for_family` | `psv_linkage.rs:53` |
| **copies = allele-vectors / paths** | `PsvColumn.per_copy` (`PsvCopyAllele`) | `psv_linkage.rs:21,31` |
| **reads = partial allele-vectors** | `genotype_family_reads` → `ReadGenotype.psv_votes` (2nd BAM pass) | `psv_linkage.rs:517,420` |
| **assignment** (thread reads → copies) | `assemble_psv_isoforms` | `psv_linkage.rs:798` |
| **K-frontier gate** (skip if $<$ `min_psv_columns`) | `family_identifiability_gated` / `psv_columns_and_identifiable` | `psv_linkage.rs:340,350,372` |
| **entry point** (end to end) | `run_layer2` | `layer2.rs:537,985` |

So the abstraction and the implementation are the *same object*: `psv_columns_from_reference` POA-aligning the copies' exon
sequences **is** the $\sim_B$ backbone yielding the bubbles; `genotype_family_reads` + `assemble_psv_isoforms` **is**
threading reads to paths; `family_identifiability` **is** the K-frontier cut. The PSV-aware variation graph is not a future
plan — it is what the `--vg` layer-2 path runs, and it is unit-tested (`family_graph.rs` tests:
`family_merge_default_merges_near_identical_paralogs`, `fragment_in_larger_exon_does_not_merge_under_jaccard_default`,
`shared_low_complexity_tract_does_not_merge_at_default`).

**The one open lever (a deliberate flip, not a gap).** The copy-assignment *family-merge* step still gates the $\sim_B$
backbone-coverage check **opt-in** — `family_min_core_coverage()` defaults to `0.0` (off), enabled by
`RUSTLE_VG_FAMILY_MIN_CORE_COVERAGE` (`family_graph.rs:416,447,1066`) — for byte-identical-default discipline, even though
`family_detect.rs` already uses the same $\sim_B$ test as *the* criterion. Making it the production default is a single
threshold flip ($\tau$), but it changes output, so it is owned by a genome-wide validation pass (the byte-identical-off
discipline), not a blind change. That validation is the last step to make the formalized $\sim_B$ definition the *default*
substrate of assignment, rather than the available one.

---

## Summary

The family definition does not merely *precede* copy assignment — it **constructs its input**. The $\sim_B$ backbone is
the coordinate system; its bubbles are the columns $[m]$; the copies are paths; the reads are partial paths; and the one
identifiability structure (the distinguishing columns $D_{ij}$) governs both whether the loci are a family and whether
its reads can be assigned. Detection and resolution are two readings of one variation graph — which is why making the
graph PSV-aware is the whole pan-transcriptomic copy-resolution method in a single object.
