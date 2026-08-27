# The seeded multi-copy family: definition, properties, certificate

Status 2026-08-10. Numbers below are measured on CHM13 (human) and NHGRI_mGorGor1-v2.0 (gorilla);
scripts in `bench/crossspecies/`. Anything not yet measured is marked **OPEN**.

---

## ⚠⚠ TIER NOTICE — READ BEFORE QUOTING ANY NUMBER IN THIS DOCUMENT (2026-08-10)

**Every `E_r` number is a statement about ONE alignment tier, and until 2026-08-10 this document never
said which.** Two tiers were in play and they are not interchangeable:

| name used below | argv | who ran it |
|---|---|---|
| **SHIPPED** | `minimap2 -c -X --no-long-join -t N -k 11 -w 5` | the binary (`ER_TIER_FLAGS`, `denovo_pipeline.rs`) |
| **PANEL** | `minimap2 {-x asm20 \| -k 11 -w 5} -c --eqx -N 200 -p 0.02 -t N` (**no `-X`**) | the eight `bench/crossspecies/` scripts, incl. `multifamily_p1p4.py` — the script behind P1–P4 |

`-N` and `-p` are **inert** at this tier. `-X` is the operative difference: it implies `--dual=no`, so
exactly ONE orientation per pair is emitted. The PANEL scripts additionally kept an `-x asm20` leg that
the shipped default **skips** (`RUSTLE_ER_SENSITIVE_ONLY` has defaulted true since 2026-08-07).

Measured over 14 panels on **byte-identical FASTA**: the **partition** differs on **4/14** and the
**edge count** on **10/14**. This is not a rounding difference.

**FIXED 2026-08-10.** The tier now exists once, as `ER_TIER_FLAGS` / `er_tier_command()`
(`denovo_pipeline.rs`), mirrored once in Python as `bench/soto/rustlib.py::ER_TIER_FLAGS` / `ava()`;
a unit test fails if the flag string is ever re-typed anywhere in shipping Rust. All eight panel
scripts now call the shipped command. **A `RUSTLE_ER_EDGE_DUMP` run writes `<prefix>.params.tsv`
carrying `mm_args_sensitive`, `coverage_form`, `coverage_denominator` and `identity_metric` — settle
"same rule?" by diffing two `params.tsv` files, never by reading code.**

⚠⚠ **THAT SENTENCE WAS AN O1-ONLY CAPABILITY UNTIL 2026-08-11, AND THE TIER HAD ALREADY DRIFTED
BEHIND IT.** `params.tsv` was written only by `homology_edges_all_reps_pooled`, which `copy_assign`'s
`refine_families_exon_sum` never calls (measured 25/25 written on the O1 side, **0/25** on the O2
side), and refine called `primary_seed_args()` unconditionally, so the 2026-08-07 sensitive-only
default never reached it: O2 ran `-x asm20 @ 0.80` while O1 ran `-k 11 -w 5 @ 0.60`.

**FIXED (X.4).** One predicate `er_sensitive_only()` and one selector `er_primary_tier()` now serve
both call sites, and both emit a **`<prefix>.rule.tsv`** — the edge-deciding knobs *only*, with no
counts, lengths or paths — from the shared `er_rule_rows()`. The parity question is now literally

```
diff <prefix>.rule.tsv <prefix>.refine.rule.tsv     # O1 catalog edge  vs  O2 refine edge
```

⚠⚠ **SUPERSEDED 2026-08-13 (O-3) — THE EMPTY DIFF WAS AN ARTEFACT OF A BLIND FILE.** The measurement
below stands as made, but `rule.tsv` had **no row for the substrate the run actually used**: it wired
`substrate` to `params.include_introns` and so printed `exon-sum` for refine runs whose edge set is
`E_x ∪ E_g`. With `core_substrate` + `additive_genomic_tier` in the file, the shipped-default diff is
**one line**, not zero — `absent (single-substrate site)` (O1 catalog, which *swaps* substrate) vs
`armed (genomic-span …)` (O2 refine, which *unions* one in). Quote the diff as **one line on the
substrate axis, empty on every other edge-deciding knob**, never as "the same rule".

**Measured 2026-08-11 (gorilla MAGEA control region, 5 loci): EMPTY DIFF, on all 5 refine calls.**
Under `RUSTLE_ER_SENSITIVE_ONLY=0` — which reproduces refine's pre-fix tier — the same diff prints
five lines (`-k 11 -w 5`/`0.600000`/`sensitive` vs `-x asm20`/`0.800000`/`asm20`), so the instrument
demonstrably *can* fail. `params.tsv` remains a superset of `rule.tsv` on both sides; **only
`rule.tsv` is the parity object** — `params.tsv` carries input-dependent rows that will always differ.

⚠ **Scope.** This concerns the ALL-VS-ALL (`E_r`) pass only. The seed→**genome** pass that locates
copies legitimately keeps `-N 200 -p 0.02` and no `-X`, and was **not** changed. **P1, P3 and P4's
locus counts therefore do not move**; P2's densities and P4b's edge sets do.

> ⭐ **RESOLVED 2026-08-11 — every §2 number below is now a SHIPPED-tier number and names its species;
> the prediction in the paragraph above is confirmed on both counts.** P1's 19 seeds, P3's two species
> and P4's three probes reproduce to the digit (they never touch `E_r`). **P2 moved by exactly one edge
> in total** — human NPIP is byte-identical between the tiers (351 edges, `d` = 1.000, Jaccard 1.0000);
> gorilla NPIP loses one of 261 to a **0.0011 coverage near-miss** (`c` = 0.4989 against the 0.50 floor)
> and keeps its single component. **P4b moved a lot** — RNA 338 → 260 edges, Jaccard 0.963 → 0.741 —
> **but `E_RNA ⊂ E_DNA` survives (0 RNA-only edges) and all 91 DNA-only edges are coverage misses at
> identity ≥ 0.9709.** ⚠⚠ **THE CONTAINMENT IS PANEL-SPECIFIC — see §1★.3.** It holds on this 27-node
> `o1_joint` NPIP panel and **fails on the 80-node `o1_closure` panel: 33 RNA-only edges of 481 (6.9%)**.
> What survives everywhere is the PARTITION-level claim: exactly **1 of those 33** joins two DNA
> components. Never quote "0 RNA-only edges" without naming the panel.
> ⚠ **§3b's N0 table is at a THIRD tier and stays OPEN.**

### Coverage form (defect M1), fixed the same day

`c` was computed as `(qend-qstart) / min(qlen,tlen)` — a **query-axis numerator over a possibly-target
denominator**, which is not an aligned fraction of anything. With `-X` the query is the **longer**
sequence in **54.5%** of records (25-region gorilla control panel), so the statistic exceeded 1.0 on
**21/255 records (8.2%), maximum 1.142** — a "fraction" above 1 is the defect made visible.

Now: **the numerator's axis follows the denominator** —
`qlen ≤ tlen ? (qend-qstart)/qlen : (tend-tstart)/tlen`. This is the literal reading of the rule this
document already stated ("aligned-fraction-of-shorter"); it is bounded by 1.0 and **symmetric under
exchanging query and target**, which the old form was not.

End-to-end cost on the 25-region gorilla control panel (75 files): **1 region changes, HERC2**, which
gains a third copy at `NC_073240.2:23,459,856-23,492,691`. That copy lies inside annotated
`LOC109023386` **"E3 ubiquitin-protein ligase HERC2"** and joins its siblings at identity **0.9940**
(coverage 0.4509 → 0.6220, floor 0.50) and 0.9518. Strict single-copy controls stay **0/7**; no new
family, no cross-family fusion, **0 edges lost and 2 gained** panel-wide. The change is a **recall
gain at a true paralog**, not merely a difference.

---

## 0. What the definition must survive

The objection this is written against is not "does it work" but "is it *arbitrary*". Three specific
forms of arbitrariness have to be closed:

1. **Annotation-dependence** — if the answer changes with which annotated gene you happened to start
   from, the family is an artifact of the annotation, not of the genome.
2. **Threshold-dependence** — if the answer changes with γ, the family is an artifact of a tuning knob.
3. **Circularity** — if the evaluation consumes the same labels that produced the answer, nothing was
   demonstrated.

Each is addressed by a property below, and each property is a measurement, not an argument.

## 1★. THE OBJECT, RESTATED — ONE OBJECT AT TWO LEVELS (2026-08-14)

Written to answer one question directly: *is there a single consistent object that describes a gene
family at both the DNA and the RNA level, or are there two objects that must be reconciled?*
**There is one object.** This section states it, and every clause carries the measurement that earns
it. Nothing here changes the shipped rule — §1 below is unchanged and still normative. This is a
**restatement**, and where it corrects an earlier claim it says so.

> **A multi-copy gene family is a connected block of the graph whose nodes are genomic intervals and
> whose edges are contiguous, high-coverage homology computed on assembly sequence. RNA determines
> which intervals are nodes and how far they extend. RNA does not determine an edge.**

### 1★.1 Why this is one object and not two

The DNA/RNA "disconnection" dissolves once the node is stated correctly. An RNA node's **bases are
fetched from the assembly** (`genome.fetch_sequence`); RNA supplies *coordinates*, the assembly
supplies every base compared. So `E_r` is computed on genomic sequence at both levels, and DNA-vs-RNA
is not two relations but **one relation over two node extents** — the full span `S_g`, or the exon
subset `S_x`. The minigraph result states the same thing structurally: *the spliced path is a strict
sub-path of the genomic path; introns are deletion-only bubbles.* A path and a sub-path, not two
graphs.

The decisive measurement is that **node length contributes nothing**: a contiguous genomic window cut
to *exactly* the RNA node's length reproduces the stored DNA edge set — **351 edges, symmetric
difference 0, 91/91 of the DNA-only pairs recovered.** What separates the levels is **contiguity**,
not length, and DNA has contiguity for free. This is why *"some single alignment record"* in §1 is
the operative clause: **the single-record requirement IS the contiguity requirement**, and the spec
has never said so in those words.

### 1★.2 What graph class a family actually is — measured, and it is NOT one class

Recomputed 2026-08-14 on the `o1_closure` panel (7 families, **80 nodes**, human/CHM13, shipped tier,
frozen `rustlib.er_edges`, `c = 0.50`, denominator `min`). Largest component of each family:

| family | \|V\| | dens | δ | degeneracy | **λ** (edge conn.) | bridges | artic. | shape |
|---|---|---|---|---|---|---|---|---|
| NPIP | 30 | 0.871 | 5 | 25 | **5** | 0 | 0 | dense |
| TBC1D3 | 12 | 0.864 | 2 | 10 | **2** | 0 | 0 | dense |
| MAGEA | 13 | 0.744 | 6 | 7 | **6** | 0 | 0 | dense |
| GSTM | 3 | 1.000 | 2 | 2 | **2** | 0 | 0 | complete |
| RABL2 | 2 | 1.000 | 1 | 1 | **1** | 1 | 0 | single edge |
| APOBEC3 | 3 | 0.667 | 1 | 1 | **1** | 2 | 1 | fragile |
| HERC2 | 12 | **0.348** | 1 | 4 | **1** | 1 | 1 | fragile |

⚠⚠ **This RETRACTS the stated reason for γ's inertness.** §1's parameter table says γ is inert
*"on real families — they are complete graphs."* **They are not.** Only 2 of 7 are complete, and both
are trivially small (n = 2, n = 3). Real densities run **0.348 – 1.000, median 0.864**. The
*conclusion* survives — γ = 0.20 is still inert, since the observed **minimum** density 0.348 clears
it with margin 0.148 — but the mechanism is *"0.20 is far below the sparsest real family"*, not
*"families are cliques."* ⭐ Note also that **HERC2 at 0.348 would SPLIT at γ = 0.40**, which is
exactly the documented "γ = 0.40 binds on exactly one component" observation, now localised.

⟹ **"Quasi-clique" oversells the object.** What ships is the **connected component** (γ inert on
seeded nodes); γ = 0.20 is a hairball floor that only activates when the node set is noisy — which is
why it is inert seeded and load-bearing blind (blind chr1 recovery 0.30 → 0.45; the transitive-closure
variant `ER_CC` scores **below its own null**, 0.3000 vs 0.3408, and its failure mode is lumping).

### 1★.3 The DNA↔RNA relation, measured on the same nodes

| quantity | value |
|---|---|
| `\|E_dna\|` / `\|E_rna\|` / shared | 526 / 481 / 448 |
| DNA-only edges | 78 |
| **RNA-only edges** | **33 (6.9% of `E_rna`)** |
| RNA-only edges that **join two DNA components** | **1 of 33** |
| Partition equality `P(dna) == P(union)` | **6/7** |
| `V_rna(big) ⊆ V_dna(big)` | **6/7** |
| internal DNA edges retained by RNA on `V_rna` | **448/491 = 0.9124** |

⚠ **This RESTATES the stored "RNA-only edges = 0."** That figure is **panel-specific** — it holds on
the 27-node `o1_joint` NPIP panel and **fails here**: 33 RNA-only edges across 7 families. The claim
that survives is at the **partition** level, and it is sharper for being true rather than absolute:
**exactly ONE RNA-only edge in 7 families changes the structure** (1 of 33; 1 of 481 = 0.21%), and it
is in APOBEC3 — the family already known to be pathological (articulation point, β-sensitive, DNA
fragmented into 3 components). ⟹ **RNA is not epiphenomenal, and it is not a second opinion either:
it contributes one structural edge per seven families, and the union lands on DNA 6/7.**

⚠ **Provenance**: the panel's stored `summary.json` (2026-08-09) reports NPIP RNA = 346 edges;
`rustlib.py` changed **2026-08-11** and the current frozen rule gives **339** (346 corresponds to a
coverage bar of 0.45). Numbers in this section are on the **current** rule. Do not mix them with
pre-08-11 stored artifacts.

### 1★.4 ⚠ TWO CANDIDATE FRAMINGS TESTED AND REJECTED — record these, they look right

**(a) "DNA is a quasi-clique, RNA is merely connected components." REFUTED — it is inverted.** RNA's
largest component is **denser** than DNA's in 4 of 7 families (NPIP 0.897 > 0.871, TBC1D3 1.000 >
0.864, MAGEA 0.848 > 0.744, HERC2 0.455 > 0.348), sparser in 2, equal in 1. The reason is not that
RNA is tighter: **RNA sheds nodes** (NPIP 1 component → 3), so the survivor is denser by
construction. Density comparisons across differing node sets are not comparisons.

**(b) "The RNA family is a k-core of the DNA family." NOT ESTABLISHED — and the first test of it was
a Simpson's-paradox artifact.** The description *fits* the big dense families exactly — NPIP's RNA
component **is precisely the 11-core** of the DNA graph, TBC1D3's the 3-core — but the mechanism it
implies (RNA trims the low-degree periphery) does not hold. Pooling raw degree across families gave
AUC 0.7668, permutation **p = 0.0086**; that pooling is invalid because degrees range from ~3 (HERC2)
to ~27 (NPIP), so the test measured family size, not trimming. ⭐ **On the correct within-family unit
(normalised degree), the dropped node is lower-degree than the kept median in only 2 of 5 families,
exact sign test p = 0.8125.** k-core is a *description* of two families, not a mechanism. **Do not
put it in the definition.**

### 1★.5 ⭐ WHAT THE GRAPH THEORY DOES BUY: λ AS A CERTIFICATE, NOT A DEFINITION

The one invariant that separates the families that have never caused trouble from the ones that
always have is **edge connectivity λ** — and it is already implicit in the project's own conclusion
*"use CUT EDGES, not margin."* Measured above: **λ ≥ 2 for NPIP, TBC1D3, MAGEA, GSTM; λ = 1 with a
bridge for RABL2, APOBEC3, HERC2** — and the λ = 1 set is exactly the set of families that has
produced every β-sensitivity, every articulation-point split, and the only γ = 0.40 casualty.

> **SHIPPED 2026-08-14 — the certificate, and it is not a change to the definition.** Every emitted
> family carries λ, or minimally the predicate `λ ≥ 2`. It is exact, cheap (Stoer–Wagner global min
> cut, `family_split::edge_connectivity`), needs **no threshold**, and states a guarantee an examiner
> can check: **at λ ≥ 2 no single alignment record's loss can split this family.** At λ = 1 the family
> hangs on one edge and is reported as such.

This is the honest version of "a cleaner combinatorial object": not a new family class, but a
**theorem-shaped confidence statement attached to each family it emits.** It costs nothing, changes
no partition, and is falsifiable per family.

⭐ **PURPOSE — the question it answers, stated before the mechanism.** The catalog previously asserted
that a family is a connected block of `E_r` and gave the reader no way to ask *how strongly* it is
connected. Two rows reading `n_copies = 7` could be a 7-clique and a 7-node chain; nothing
distinguished them. λ is the smallest statistic that separates those cases, and it is the one that
predicts which families have historically been fragile — the λ = 1 set {RABL2, APOBEC3, HERC2} is
exactly the set that produced every β-sensitivity, every articulation-point split, and the only
γ = 0.40 casualty (§1★.2).

⚠⚠ **WHY IT IS NOT A MEMBERSHIP CRITERION — the decisive argument, not a preference.** A 2-copy
family has `λ = 1` necessarily. Enforcing `λ ≥ 2` would therefore delete **every 2-copy family**, and
2 is the modal family size. `cut_certified = false` on a 2-copy row is a statement about what a
2-node graph can be, **not a defect flag**, and must never be read or filtered as one. This is pinned
by `certificate_two_copy_family_is_always_lambda_one`.

⚠ **What λ does NOT tell you.** It is a property of `E_r` on the emitted copies, so it inherits every
limit of `E_r`: it says nothing about whether the family is biologically real, whether a copy is
missing (that is O3), or whether reads can be assigned to copies (that is O2). A false-merged family
can have high λ — λ measures how firmly the block holds together, not whether it should exist.

### 1★.6 What the restatement concedes — state these first

1. **RNA cannot be shown to add edges in general** — 1 of 481 changes structure. The relation is
   DNA-computed at both levels; the joint object is a **property**, not a definition (§P4c).
2. **The node set is an input, not derivable** — settled by the blind arm (delineation, not
   detection, is what fails: 245/245 lost edges fail COVERAGE).
3. **Reach is bounded by DIVERGENCE, not by dispersal** — 22/40 = 0.5500 [0.3983, 0.6929] on chr1,
   and chr1 is representative (Fisher p = 0.6090). See §4a.
4. **Families are not one graph class.** Density spans 0.348–1.000 and λ spans 1–6 across seven
   families. Any single-class claim ("clique", "quasi-clique", "k-core") is false as stated.

---

## 1. Definition

### ⭐⭐ 1.0 THE OPERATING ASSUMPTION (stated 2026-08-12, `/home/juanfra/winloci_scratch/o1_minann/O1_MINIMAL_ANNOTATION.md`)

> **O1 assumes the genome, a minimal annotation, and reads aligned to that genome.** "Minimal" is
> bounded on three independently measured channels — INTERVAL (which loci have coordinates), SYMBOL
> (which carry an informative name), FAMILY (which carry a copy assertion). **O1's task is to supply
> the copy relation, for members the annotation never named.** It is not to discover loci ab initio.

**Granted.** A reference good enough to align reads against ships with a gene model — reads exist only
relative to the reference they were aligned to — and the annotation supplies exactly the thing
measured to bind: because M1 coverage divides by the **shorter** sequence, an `E_r` edge needs **ONE**
well-delineated endpoint, not two (122/122 of the blind boundary arm's surviving pairs; 4 exact nodes
of 27 recover 79 of the 80 pairs touching them).

**⚠ CONCEDED, PERMANENTLY: the fully annotation-free claim (row O1.13).** The concession is the honest
one because the blind failure was never a homology failure — detection was fine (**0/27 loci absent**;
chr1 gene capture **120/162 = 0.7407**) and **delineation** bound (blind nodes median **5.97×** too
long; of 245 lost edges **0 lack a record, 0 fail identity, 245/245 fail COVERAGE**). Independently
reproduced under this reframe on a different substrate and node construction: **152/152 = 100%** of
within-family pairs that have a record and still fail, fail the coverage clause — at *every*
completeness level tested. Delineation is an ab-initio gene-finding sub-problem O1 never claimed.

**⭐ WHAT THE ASSUMPTION REQUIRES OF AN ANNOTATION — measured, human/CHM13, shipped tier, chr1 (K=1),
40-family GFF denominator fixed before the run.** Two axes ablated independently; both knee criteria
pre-declared before any surface existed.

| axis | tolerated | breaks at | evidence |
|---|---|---|---|
| **interval completeness `c`** | gain present at **c = 1.00** | **gone by c = 0.50** (Δ +0.0250, P = 0.3625); **exactly +0.0000, CI [0,0], at c ≤ 0.25** | oracle ceiling 40 → 23 → 5 → 1 |
| **boundary error `d`** | **±2 kb costs literally nothing** (recovery, Δ and bootstrap bit-identical at d = 0/1000/2000) | **cliff at 5 kb** — false-merge rate **0.0806 → 0.5000**, Wilson disjoint, Fisher p = 2.7e-11 | γ = 0.20 goes inert → load-bearing |

⚠ The **c = 0.75** row and the **d ∈ {250, 500}** interior were **NOT RUN** — the completeness knee is
bracketed to `(0.50, 1.00]` and **not located**; cells not run are never interpolated. **OPEN.**

⭐ **The rule does not break as the annotation degrades — it is STARVED.** Within-family pairs collapse
277 → 58 → 8 → 2 while the **fire rate holds at 0.12–0.18 down to c = 0.25**. And the ceiling
decomposition separates two failures that must not be conflated: at c = 0.50 the annotation still
permits 23 families and the method takes 1 of 10 headroom (**a method limit**); at c ≤ 0.25 headroom is
1 then 0 and the method takes it (**an annotation limit — no edge rule could do better**).

⭐ **~5 kb of boundary error is what converts a seeded arm into a blind one** — annotation quality and
the blind failure mode are the same axis seen from two ends. ⚠ **The panel's "±1 kb buys ~90%" does
NOT transfer to this scope: here ±1 kb buys 100%.** chr1 truth genes have median length **18,278 bp**,
so ±1 kb is a ~2.7% perturbation — the panel's pricing is about **relative**, not absolute, boundary
error. Confirmed by the one quantity with power: `fire SHORT` slides monotonically to 34% of baseline
while `fire LONG` is flat and ends *above* it.

**⚠ WHAT THE ASSUMPTION DOES NOT BUY — state these first.** Gain over the unimproved annotation at the
pre-declared τ = 0.50 headline is **+0.1000 [+0.0000, +0.2250], P(Δ≤0) = 0.0635 — NOT significant**;
the surviving claim is the m ≥ 3 stratum (**5/22 vs 0/22 for both degenerate controls**, size-matched
null 0.003, p_perm = 0.0000, and it survives a position-respecting null whose mean is **124× higher**).
**Dispersed families: 0/8 at every rung, p = 1.0000 — and 0 of the 8 have even ONE within-family `E_r`
edge, so this is a RULE-NEVER-FIRED failure, not a grouping failure.** Below ~0.75 identity the method
is **exactly the all-singletons floor** against a declared floor of 0.60. The gain is **not monotone**
(one m=2 family is lost that all-singletons recovers). chr1 is 80.0% clustered vs 36.3% genome-wide
⟹ **no number here is genome-wide**, and scopes may never be pooled or differenced (the `-k 11` tier is
not scope-invariant). ⚠ **Robustness to a WRONG annotation is untested** — supplied labels and truth
come from the same T4 parse, so the input can be incomplete but never wrong. **OPEN.**

---

Rewritten 2026-08-09 to state only what measurement supports. The previous version listed **five**
constants (`τ, c, β, γ, d`) as equal parts of the definition. Three of them have now been measured
against real data and cannot stand where they were: `τ` cannot bind, `γ` cannot bind, and `β`/`d`
belong to NODE CONSTRUCTION, not to the family relation. **The edge rule has one operative constant,
`c`, and it has a measured plateau: the emitted 25-region catalog is byte-identical over
`c ∈ [0.50, 0.60]`, with the shipped 0.50 at that interval's lower edge.** Every demotion below is a
measurement, and each carries the number that produced it — including the one that went the wrong way
for the story (`β` is operative, and the claim that it does not change the partition is refuted).

Fix a genome `G` and an alignment tier `T`.

> **Node set.** `V` = the loci proposed **from the data alone**: read-supported transcribed loci
> assembled from the BAM (`pass1_skeletons` → `assemble_gate`), or, under `--from-genome`, duplicated
> genomic loci found by genome self-alignment. **`V` does not depend on any seed** and does not depend
> on annotation.
>
> **Substrate `S` — a term of the definition, not an implementation detail (added 2026-08-13, O-2).**
> `S(v)` = the bases an edge is computed **on**: either `S_x(v)` = the **exon-sum** (the copy's spliced
> model, `refine_copy_seq(v, None)`) or `S_g(v)` = the **genomic span** `G[start(v), end(v))`
> (`refine_copy_seq(v, Some(G))`). ⚠⚠ **`E_r(S_x)` and `E_r(S_g)` are different relations on the same
> nodes, and neither contains the other**: measured on the binary over 244 human reps, `|E_x| = 696`,
> `|E_g| = 1575`, shared 630, **exon-only 66**, genomic-only 945, block sets differing on **6 of 7**
> panel families, **ARI 0.4702**. A number computed under one substrate may never sit beside one
> computed under the other.
>
> **Support relation.** `E_r` = the unordered pairs `{u,v} ⊆ V` for which **some single alignment
> record** between `S(u)` and `S(v)`, at tier `T`, covers at least `c` of the shorter of the two.
>
> *(The shipped code also requires that record's identity to clear `τ`. That clause is retained in the
> code and omitted from the definition here because it is measurably incapable of rejecting anything —
> see `τ` below. It is written as `id ≥ τ ∧ cov ≥ c`; it BEHAVES as `cov ≥ c`.)*
>
> **Family.** The families are the blocks of the γ-quasi-clique refinement of the connected components
> of `(V, E_r)`. **The blocks PARTITION `V`** — the refinement starts from `all_components`, singletons
> included (`family_split.rs:477`), and only ever splits.
>
> **Seed — a QUERY, not a term.** `F(s)` = the block containing `s`. Because the blocks partition `V`
> and `V` is seed-free, this is defined for every node, and `F(s') = F(s)` for every `s' ∈ F(s)`.
> **P1 is a theorem about a partition, not a measurement.**
>
> **Certificate — part of what is EMITTED, never of what is DECIDED (added 2026-08-14).** Every
> emitted family is written as a pair `(B, λ(B))`, where `λ(B)` is the **edge connectivity** of `E_r`
> induced on `B`'s emitted copies: the minimum number of alignment records whose removal would
> disconnect the family. `λ ≥ 2` certifies that **no single record's loss can split this family.**
> ⚠⚠ **λ is not a term of the definition and must never become one**: a 2-copy family has `λ = 1`
> NECESSARILY — one edge is all a 2-node graph can hold — so gating membership on `λ ≥ 2` would delete
> every 2-copy family, the most common size in the catalog. It reports confidence; it never decides
> membership. **Nothing in the pipeline branches on it, and it carries no threshold.**

**Emission ≠ definition.** The partition covers all of `V`; the catalog *writes out* only blocks with
at least `--min-copies` (default 2) spatially distinct loci. A block below that bar exists in the
partition and is simply not printed. Say "emitted family" when that is what is meant.

⭐ **What the emitted family now carries, and why each field is there.** `<out>.families.tsv` gains
four columns, appended last so existing header-keyed readers are unaffected:

| column | meaning | why it is emitted |
|---|---|---|
| `n_edges` | `E_r` edges induced on the emitted copies, de-duplicated | the evidence count the family rests on — a reader can divide it by the copy count |
| `density` | `2\|E\|/(n(n−1))` | the γ statistic, so **γ's inertness is checkable per family** rather than asserted |
| `lambda` | edge connectivity λ | the certificate: how much of that evidence is load-bearing |
| `cut_certified` | `λ ≥ 2` | the one-bit form of the guarantee |

⚠ Three honesty constraints are enforced in code, not by convention. (i) λ is computed on the
**EMITTED** node set — after the coverage split and after `distinct_locus_reps` merges co-located
copies — because that merge is exactly where the node set changes, and a λ measured before it would
describe a different object than the row printed (`certificate_is_measured_after_the_locus_merge_not_before`).
(ii) Several alignment records for one copy pair are **one** edge, since only one edge's loss can
actually split the family. (iii) On paths that build no `E_r` graph (the conflict catalogs) and under
`--refine` (which re-clusters over a different edge set), the columns print **`NA`, never `0`** — a
missing certificate and a genuinely disconnected family must not read the same.

Shipped values (2026-08-09), with the measured role of each:

| symbol | shipped | where it lives | measured role |
|---|---|---|---|
| `c` | **0.50 of the SHORTER locus** | the edge rule | **OPERATIVE — the only one.** Plateau [0.50, 0.60], shipped value at its lower edge |
| `S` | ⚠ **TWO DIFFERENT VALUES SHIP.** O1 catalog: `S_x` (exon-sum) — `homology_genomic_span` default **OFF**, an ungated global SWAP. `refine` (`copy_assign`, conflict catalogs): `S_x` **∪** `S_g`, a **gated family-local ADDITION** | the edge rule | **OPERATIVE, and the axis on which the two call sites genuinely differ** — see §1.1 |
| `T` | `minimap2 -c -X --no-long-join -k11 -w5` | the edge rule | fixed input; ⚠ see *T is not one invocation* below |
| `τ` | 0.60 (`sensitive_identity`) | the edge rule | **INERT** on `[0, 0.65]` — cannot reject anything |
| `γ` | **0.20** (`RUSTLE_GENOME_GAMMA`) | the refinement | **INERT** on seeded real families. ⚠ **NOT because they are complete graphs — that reason is RETRACTED, see §1★.2.** Only 2/7 are complete (both n ≤ 3); densities run 0.348–1.000, median 0.864. γ = 0.20 is inert because it sits **0.148 below the sparsest observed family** (HERC2, 0.348 — which *would* split at γ = 0.40). ⭐ γ is **load-bearing when the node set is noisy** (blind chr1 0.30 → 0.45) |
| `β`, `d` | 5 kb, 10 kb | **NOT in the shipped pipeline** — the seeded probe's node leg | **OPERATIVE where they exist**; report them with any probe result |

### Why each demoted constant is demoted (and what remains free)

**`τ` = 0.60 — INERT. It sits below the aligner's own emission floor.**
Three independent lines, all negative. (i) A **634,502-record census** of non-self sensitive-tier
records over six corpora (including 106,568 deliberately CROSS-family genomic pairs and 161,818 NPIP
self-alignment records): minimum reported identity anywhere = **0.6430**, zero records below 0.60;
restricted to the 19,596 records that clear the coverage floor — the only records `τ` could ever
block — the minimum is **0.6582**. (ii) The direct counterfactual: recomputing `E_r` with the floor
**removed entirely (τ = 0)** gives edge-for-edge identical DNA, RNA and union edge sets and identical
component counts on two node constructions and on an 18,418-transcript panel (9,468 pairs at τ = 0 and
at τ = 0.60). `τ` is inert across `[0, 0.65]`; the first bite is at **0.70**. (iii) 99 adversarial
pairs built at known true identity 0.30–1.00 (a substitution ladder plus mosaics of conserved islands
in unrelated filler) never produced a reported identity below **0.6469**.
*Mechanism, which is the reportable part:* identity = `1 − de` is measured **on the aligned block**, so
rising divergence does not lower the reported identity — it SHORTENS the block. True 0.66 → coverage
0.975 at identity 0.686; true 0.62 → coverage 0.063 at identity 0.701; true 0.60 → **no record at
all**. **Divergence exits through coverage, not identity**, which is also why 100% of cross-family
rejection on the panel is coverage and 0% is identity, and why 9 of 9 misses on the external paralog
panel failed on coverage.
⚠ This is empirical + mechanistic, **not a theorem**. The natural derivation (default `A=2/B=4` makes
an ungapped alignment break even at identity 2/3) is contradicted by the data: minimap2 emits
**negative-scoring secondaries** (3.35% of the corpus, down to `AS = −7,184` at identity 0.6833). The
conditional statement is the one that holds everywhere looked: every `AS > 0` record in 634,502 has
identity ≥ 0.6822. A future substrate or minimap2 version could in principle emit a sub-0.60 record.
⚠ And disclose the upstream help: the panel's node constructor applies its own per-record filter
`nmatch/blocklen ≥ 0.80` *before* anything reaches `E_r` — a stricter identity floor one stage
earlier. The all-vs-all census corpora are not subject to it, so the aligner-floor claim survives, but
panel-internal evidence for `τ` is not independent of it.

**`γ` = 0.20 — INERT, by structure and then by measurement.**
Structurally, `gamma_quasi_clique_partition` keeps a block iff `block.len() <= 2 ||
induced_density >= γ`, so γ can only bind on a component with `|C| ≥ 3` whose density is below γ.
Measured **on the β/d panel of this section**: every family but GSTM is effectively a complete graph, so
γ = 0.20 and γ = 0.40 produce **byte-identical** output, and all 21 components clear 0.40 with minimum
density **0.4722**.
⚠⚠ **DO NOT GENERALISE THAT SENTENCE — corrected 2026-08-14, see §1★.2.** On the 80-node `o1_closure`
panel only **2 of 7** families are complete (both n ≤ 3) and densities run **0.348–1.000, median 0.864**.
The two minima are different panels and must not be cross-quoted (0.4722 here, **0.348** there). The
conclusion that γ = 0.20 is inert survives — it sits 0.148 below the sparsest observed family — but the
REASON is "0.20 is far below the sparsest real family", **not** "families are cliques", and at γ = 0.40
**HERC2 (0.348) would split**. Human NPIP is `27·26/2 = 351` edges on 27 nodes — density exactly
**1.000**, so that family is a γ-quasi-clique for *every* γ ∈ (0,1] and no threshold was chosen.
⚠ **Inert is not harmless.** γ *is* the only stated defence against over-merge, and it catches
over-merges roughly half the time: chr15's 61-locus GOLGA8 blob sits at density 0.210 and is correctly
rejected, but chr7 merges SPDYE with PMS2P at density **0.660**, comfortably above γ, and is admitted
(§3b). "γ never binds on real families" and "γ cannot be relied on to split false ones" are the same
fact seen from two sides.
⚠ **Four values for γ across four documents** — this file said 0.40, the binary defaults to **0.20**
(`RUSTLE_GENOME_GAMMA`), `bench/DEFINITIONS_FORMAL.md` says 0.20 plus a theoretical γ_core = 0.13, and
`bench/VALIDATION_AND_STATUS.md` says 0.30. **0.20 is what ships**; 0.40 described the seeded probe.
The P2 margins in §2 are computed against 0.40 and are LARGER at 0.20 (5.0× and 4.35×), so this was a
wrong-number problem, not a wrong-conclusion one. Never quote a margin without naming its γ.

**`β` = 5 kb and `d` = 10 kb — NOT DEMOTED. They are node-construction choices, and they were tested
and found operative.** The earlier reading "β changes the node count but not the partition" is
**REFUTED on both halves.** Node set: at the documented `d`, β = 0/2k/5k/7.5k/10k/15k/20k/30k gives
**73/73/69/64/62/56/52/45** nodes — **−38%** over the swept range, strictly nested, so β only ever
deletes. What it deletes are annotated members *with read support*: at β = 10 kb MAGEA loses MAGEA1
(79 spliced reads), MAGEA8 (34), MAGEA10 (51), MAGEA12 (69), GSTM loses GSTM3 (1,867 reads), APOBEC3
loses APOBEC3A. Partition: **APOBEC3A is a cut vertex** — at β ≤ 5 kb the union partition is one
component `{A, B, CDF, GH}`; deleting A at β = 10 kb splits it into `{CDF,GH}` and `{B}`. That is 1 of
35 (β, family) cells on each tier, but one is enough to kill a universal claim. The old verdict was
measured over a **two-point window {2 kb, 5 kb}** that sits inside the region where β is dominated by
an **undocumented 1 kb per-record block floor** in the same constructor (β = 0 with the span floor also
removed gives 74 nodes vs 73). And β is an absolute-bp threshold compared against gene size: MAGEA
genes span 3.3–10.7 kb, so β ≥ 15 kb annihilates MAGEA **by construction**, while NPIP (10–50 kb)
survives to 30 kb — the same non-transferability that killed eight earlier absolute-threshold rules
here. `d` is decisive too and equally unjustified: union == DNA is **6/7 at d = 10 kb but 7/7 at
d = 9.5 kb**, γ-legal at both.
⭐ **The resolution is not to tune β and d — it is that they are not in the shipped definition at
all.** They parameterise `V(s)`, the *seeded* node leg of `bench/crossspecies`. Once `V` is seed-free
(above), the shipped pipeline has no β and no d. **This does not make the problem disappear; it names
where it lives.** The shipped node leg has its own constants and they must be reported the same way:
RNA — `pass1_min_reads = 2`, `min_terminal_support = 2`, `GATE_MIN_READS = 3`; DNA `--from-genome` —
`min_identity = 0.90`, `min_block = 1,000 bp`, `max_locus_span = 3 Mb`. **Node construction, not
homology, remains the binding constraint on the whole definition**, and any result quoting a node set
must quote its constants. ⚠ **WITHDRAWN 2026-08-11: the numbers that used to appear here — "blind
self-alignment tops out at purity 0.136, node purity at 0.237" — are no longer quotable.** No artifact
survives (the one candidate file is 0 bytes); their substrate is recoverable only in prose (HUMAN, CHM13
v2.0 self-alignment of chr16:11.0–31.5 Mb + 79.7–80.7 Mb, 20,007 records, 677 nodes → 510 components),
and **three** different purities circulate for the one experiment (0.136 = 19/140 best component,
0.156 shipped rule at id ≥ 0.90, 0.237 = 19/80 gene-level ceiling), the last of which is
**oracle-selected** — it presumes you already know which of 377–510 components to look at. The
qualitative claim stands on its own and needs no number. See `OBJECTIVES_AND_VERIFICATION.md` row 1.13.

**`c` = 0.50 — OPERATIVE. It does 100% of the discrimination, and it is the only constant that has to
be justified on its own.** The honest headline is:
> **the definition is indexed by `c`. The shipped operating point is `c = 0.50` of `min(|u|,|v|)`,
> and it sits at the LOWER EDGE of the one interval over which the emitted catalog is invariant.**

Swept end-to-end through the shipped binary (`RUSTLE_GENOME_MIN_COVERAGE`, 27 values from 0.10 to
0.90) over the 25-region gorilla control panel, and over an external 15-positive / 15-negative human
paralog-pair panel:

| `c` | strict single-copy FPs (7 controls, external label) | 25-region panel | paralog panel: recall / false merges |
|---|---|---|---|
| 0.10 | **7 / 7** | 7 fam / 33 copies | 13/15 / **9** |
| 0.30 | 2 / 7 | 7 fam / 30 copies | 8/15 / 6 |
| 0.40 | 1 / 7 | 6 fam / 26 copies | 7/15 / 3 |
| 0.46 | 1 / 7 | GTF2H1 (a single-copy control) emits a family | — |
| 0.48 | 0 / 7 | 6 fam / **25** copies (HERC2 admits an extra copy) | — |
| **0.50 – 0.60** | **0 / 7** | **6 fam / 24 copies — 75/75 files BYTE-IDENTICAL across the whole interval** *(PRE-M1; post-M1 this reads 6 fam / **25** copies — see the correction below)* | 6/15 / 1 → 4/15 / **0** |
| 0.65 | 0 / 7 | 5 fam / 20 copies (GSTM loses members) | 4/15 / 0 |
| 0.90 | 0 / 7 | 4 fam / 18 copies | — |

⭐ **There IS a plateau, and it is one-sided.** The emitted 25-region catalog is byte-identical at
every tested value in **[0.50, 0.60]** (0.50, 0.52, 0.55, 0.56, 0.60 — 75/75 files each, verified by
`diff -rq` against `o1_famctl/out_fixed`), and changes at both **0.48** and **0.65**. So the invariance
window contains [0.50, 0.60] and is contained in (0.48, 0.65): at least **0.10 wide**, with the shipped
value at its lower end. The two walls are of *opposite kind*, which is what makes the interval
meaningful rather than arbitrary: **below it the method starts inventing families** (0.48 gives HERC2
a spurious extra copy; 0.46 makes the single-copy control GTF2H1 emit a family; 0.10 makes all 7
strict controls emit), **above it the method starts losing real members** (0.65 costs GSTM). The
single-copy-control column is the non-circular half of this — those 7 genes are labelled single-copy
by external annotation, and their FP count falls monotonically 7 → 2 → 1 → 0 and then stays 0.
⚠ **The "6 families / 24 copies" column is partly circular** and must not be quoted as independent
evidence: that target was *defined* from the c = 0.50 run, so 0.50 passes it by construction. The
non-circular content is the *width* of the byte-identical interval above 0.50 and the identity of the
two things that break it.

> ⚠⚠ **CORRECTION, 2026-08-11 — this whole table is PRE-M1, and the M1 fix removed its lower wall.**
> Two things must be said before any figure above is quoted.
>
> **(1) What "6 fam / 24 copies" counts.** It is the 25-region panel **excluding SHARP**, which this
> document reports separately as "SHARP 1f/2c" (see § refine-flip). Including SHARP the same run is
> **7 families / 26 copies**. Both are recounts of
> `/home/juanfra/winloci_scratch/o1_famctl/out_fixed/`: APOBEC3 2, GSTM 4, HERC2 **2**, MAGEA 9 + 2,
> RABL2 5, SHARP 2.
>
> **(2) Post-M1 the numbers and the walls both move.** The M1 tier fix adds HERC2's third copy
> (`DN_NC_073240.2_23459856_11`), so the post-M1 reference
> `/home/juanfra/winloci_scratch/o1fix_o2audit/fix/out_m1fix/` is **25 copies / 6 families** excluding
> SHARP, **27 / 7** including it. And the **c = 0.48 wall no longer exists**: `diff -rq` of the post-M1
> sweep arms against `out_m1fix` gives **0 of 75 files differing at 0.48, 0.55 and 0.60**, with 6
> differing at **0.46** and 6 at **0.65**. The byte-identical interval on this panel is therefore
> **[0.48, 0.60]**, not [0.50, 0.60], and the "extra copy" that defined the old wall is the copy the
> pre-registered truth list scores as **correct** (HERC2 recall 2/5 → 3/5, precision 3/3).
>
> **(3) None of that rehabilitates `c = 0.50`.** A *wider* byte-identical interval on this panel is
> evidence the panel is **insensitive**, not that the parameter is flat: on 35 gene-tight windows the
> output changes at 0.50 → 0.51 and again at 0.51 → 0.55, and the widest interval containing 0.50 over
> which nothing changes is **(0.4721, 0.5000]** — 0.50 on the **upper boundary**. See
> `docs/OBJECTIVES_AND_VERIFICATION.md` row 1.7, which records this sweep as **REFUTED**. Treat the
> table above as a record of what the control panel can and cannot see, not as a plateau argument.
⚠ **The plateau does not make 0.50 optimal — it makes it defensible.** Inside the window the two
paralog-panel axes move in opposite directions monotonically: 0.50 gives recall 6/15 with 1 false
merge (ATP1A1/ATP4A), 0.60 gives 4/15 with 0. No value is best on both. What can be said without
circularity is narrower: **0.50 is the smallest value at which no single-copy control emits a family
AND the panel catalog stops changing with `c`.** Calling the 25th HERC2 copy at 0.48 "spurious" is
already the circular reading; the non-circular statement about 0.48 is only that the catalog is still
moving there.
On the edge-level substrate measurement the operating point also carries **0 of 1,383 cross-family
edges on any substrate**, but the safety margin below it is substrate-dependent — on genomic span
cross-family edges appear at c = 0.30 (1 edge), 0.20 (22), 0.10 (141); on the rep transcript only at
c = 0.10 (29), zero at 0.20 — **0.30 of headroom vs 0.10**. Attempts to replace the fraction with an
absolute aligned-base floor were **refuted** (smallest safe `N` is 5,000 bp / 2,000 bp / 500 bp on the
three substrates — a 10× spread with no shared plateau; §1a). Two further properties of `c` that must
travel with it:
*(i)* it is **scale-dependent** — the same 0.50 demands a median 6,338 bp on genomic span and 876 bp on
a rep transcript (§1a, recomputed 2026-08-10), and 5 of 7 families change partition when the substrate
changes (§1a in full);
*(ii)* the statistic **was not a fraction, and now is — defect M1, FIXED 2026-08-10.** It was
`(qe−qs)/min(|u|,|v|)`: a **query-axis numerator over a possibly-target denominator**, which could
exceed 1 (on the gorilla SDHA probe **52 of 68 accepted edges scored "coverage" > 1.0, up to 2.019** —
a 12.9 kb node clearing the floor against a 6.1 kb node by denominator artifact, not by reciprocal
coverage; on the 25-region control panel, 21/255 records, max 1.142). The form is now
**`qlen ≤ tlen ? (qe−qs)/qlen : (te−ts)/tlen`** — the numerator's axis follows the denominator, so the
value is bounded by 1.0 and symmetric under exchanging query and target. **This was the "fix that first"
this paragraph asked for**; the level `c = 0.50` was not touched.
⚠ Every pre-2026-08-10 number in this file that involved `c` was computed under the old form. Where it
mattered it is now shown in both forms (§1a); the measured effect is **+0.22 to +2.29 bite points**,
and **one region of 25** changes end to end.
⚠ **Do not change `c` as a side effect of this rewrite.** 0.50 is what every number in this file was
measured at.

**`T` was not one invocation. RESOLVED 2026-08-10 — it is now, in both languages.**
The shipped binary runs `minimap2 -c -X --no-long-join -t N -k 11 -w 5`. The closure panel and the
documented panel were built with `minimap2 -k 11 -w 5 -c --eqx -N 200 -p 0.02 -t 2`. These are **not
the same tier**: `-X` implies `--dual=no`, so one orientation per pair is emitted; `-N`/`-p` are inert
here. On **byte-identical FASTA** the two give HERC2 `E_dna` **29 vs 17** (density **0.8056 vs
0.4722**), HERC2 `E_rna` 9 vs 20, APOBEC3 `E_rna` 0 vs 1. The **0.4722** quoted above as the panel's
minimum density and the single APOBEC3 RNA edge carrying the "union ≠ DNA, bridging is 6/7 not 7/7"
headline **exist only under the panel invocation**. Direction is favourable — under the shipped tier
the γ-vacuity claim gets *stronger* and the bridging break disappears.

**`T` is now defined exactly once per language and cannot drift:**
`ER_TIER_FLAGS` + `er_tier_command()` in `denovo_pipeline.rs` (it was hardcoded at **four** sites — two
`Command` builders, the `.args` dump, and the `params.tsv` row — which is *how* it drifted), mirrored by
`ER_TIER_FLAGS` + `ava()` in `bench/soto/rustlib.py`. A unit test
(`er_tier_flags_appear_as_a_literal_exactly_once_in_the_source`) fails if a fifth copy is ever typed.
All eight `bench/crossspecies/` panels now issue the shipped command; **`-N`/`-p` dropped, `-X` added,
the `-x asm20` leg removed** (the shipped default skips it). ⚠ Numbers computed before 2026-08-10 still
carry whichever tier produced them — **every table in this document now names its tier**, and any that
could not be re-run is marked OPEN rather than silently relabelled.

**Constants that are in the code but were not in this section.** Stating "one operative constant" is
only honest if the rest are named. Beyond the node-leg constants above: the panel node constructor
applies `blocklen ≥ 1,000 bp` and `nmatch/blocklen ≥ 0.80` per record plus a 2 kb span floor
(`mknodes_bd.py:12`), and the shipped engine carries a second, stricter coverage application,
`RUSTLE_COVERAGE_SPLIT` — **default 0.0, i.e. off**, verified a no-op (`coverage_edges_all_reps`
returns an empty set at ≤ 0.0, so `coverage_split_block` returns the block whole).

### The seed: what changed, and what `--seed` now is

The old §1 made the node set `V(s)` a function of the seed. Measured on the human NPIP panel, seeding
independently from each of the 19 annotated NPIP genes:

* **strict seed-invariance FAILS.** Only **4 distinct `F(s)` as sets of loci** (sizes 26/27/27/27),
  agreeing on **64 of 171** seed pairs (37.4%). The three loci that move are the three SHORTEST, all
  within 1–5 kb of β, and all **unannotated**;
* the mechanism is that **`V(s)` inherits the seed's own length** — seed NPIPB8 (10.6 kb) returns 27
  loci of ~9.9–10.6 kb, seed NPIPA1 (14.6 kb) returns 26 loci of ~14.4–14.8 kb;
* **the seed never acts through the edge rule.** For 19/19 seeds `F(s)` = all of `V(s)`, induced
  density 1.000 for 16 seeds (min 0.9544), and two different edge rules give byte-identical family
  position sets. The rule's contribution to the seeded/unseeded gap is **0**;
* over a fixed 29-position node universe the unseeded partition is **ONE component of 29** under all
  six rule × tier combinations tested, `F(s) ⊆ Comp(s)` for 19/19, `⋃ₛ F(s) = Comp(s)` exactly, and
  `F(s) = Comp(s)` for **0/19**;
* and the seeded construction has a false-positive mode the shipped binary does not: seeding from the
  gorilla single-copy control **SDHA** yields `V(s)` = 14 loci and one component at density **0.747**,
  legal at γ = 0.20 and 0.40 — a **14-copy "family" where the shipped binary correctly emits nothing.**

⟹ The seeded node set is not a defensible part of a definition, and it is not what ships. With `V`
seed-free, `F(s)` = "the block containing `s`" and **P1 stops being a claim in need of a panel.**

`--seed <chrom:start-end>` (repeatable) is implemented on `gw_family_catalog` accordingly, as a
**projection over the emitted catalog**: the run is byte-identical with and without it except for the
added `<out>.seed.tsv`, and the flag then reports the block the seed lands in by maximum overlapping
bp, printing **every member of that block**. Coordinates are 1-based inclusive. It **abstains**
(`ABSTAIN_NO_OVERLAP`) rather than snapping to a nearest locus. Its one honest limit: at
`--min-copies 2` a singleton block is not written, so an abstention does not distinguish *"the seed's
locus is a singleton"* from *"there is no transcribed locus there"* — re-run with `--min-copies 1` to
separate them. Implementation and rationale: `src/rustle/vg_family/seed_projection.rs`; behaviour
fixed by `tests/gw_family_catalog_seed_projection.rs`.

Demonstrated on the gorilla RABL2 panel region — two seeds taken at two *different* members of the
same family, one seed off any locus:

```
seed                              status              family  n  member  chrom         is_seed  overlap
NC_073231.2:151170200-151185900   HIT                 GWFAM0  5  0..4    (all 5 copies)  copy 1   15,701
NC_073235.2:15131700-15147500     HIT                 GWFAM0  5  0..4    (all 5 copies)  copy 2   15,801
NC_073227.2:1-1000                ABSTAIN_NO_OVERLAP  .       0  —       —               —        0
```

The two seeds return the **identical five-member component** — P1 exhibited rather than argued — and
`families.tsv`/`copies.tsv`/`copies.fa` are byte-identical to the same region run without `--seed`.
Note also what the run's own log says: *"2 γ-quasi-clique blocks → 1 families (≥ 2 distinct loci)"* —
the partition has two blocks and one is emitted. That is the emission-≠-definition gap, visible in a
single line.

⚠ **What survives of "annotation is a seed, not a definition".** Annotation supplies `s` and nothing
else: it does not supply `V`, `E_r`, membership, or the family's name. That claim is now *stronger*
than before, because `s` no longer touches the node set either. What is **withdrawn** is any suggestion
that the seeded probe in `bench/crossspecies` computes the same object as the shipped binary — on
gorilla it does not (MAGEA4 seed → 1 locus vs 9 shipped copies; GSTM3 seed, 4,057 bp < β → **0** loci
vs 4 shipped copies; HERC2 → 164 vs 2; RABL2 is the one clean match, 5 of 6 seeded loci coinciding with
shipped copies to within 3–46 bp on every boundary).


⚠ **The asm20 tier was retired from this path on 2026-08-07** (`RUSTLE_ER_SENSITIVE_ONLY` default
flipped to true; `=0` restores it). asm20 is a subset of the sensitive run — structurally, because
whenever `sensitive_identity ≤ min_identity` the sensitive run has both the lower bar and the denser
seeding; and measured: **0 unique edges** on genomic sequence (NPIP 2162/2210, TBC1D3 55/55, Soto
2894 edges 0 unique). ⚠ **1 unique edge on the SPLICED substrate** (NPIP RNA: asm20 150, sensitive 269,
union 270) — say "0 on genomic, 1 on spliced", not "0". Verified `cargo test --release --all-targets`:
711 passed, 0 failed.

⭐ **The shipped catalog is now ONE definition (fixed 2026-08-09).** Until this date the RNA path ran a
second, undocumented clustering stage — `refine_families_exon_sum` — after γ-quasi-clique(E_r), while
the DNA `--from-genome` path returned before it. RNA therefore emitted `refine(γ-QC(E_r))` and DNA
emitted `γ-QC(E_r)`. Refine is not a pass-through: it recomputes edges with its own **asm20 @ 0.80** and
**genomic-span @ 0.80** tiers (running asm20 *after* the engine deliberately skipped it), **without** the
core-coverage denominator or the stub guard, and clusters by **connected components**
(`denovo_pipeline.rs:3345`) rather than by the γ-quasi-clique rule above. None of those tiers or rules
appear in §1. **Refine is now OPT-IN (`--refine`) on the homology (default) catalog** and stays on by
default only for the legacy conflict catalogs (`--window-catalog` / `--cross-chrom`), which it was
written for and where it was measured to bite (APOBEC3 1f/2c → 0f/0c, SHARP 2f/4c → 1f/2c).

Cost of the flip at shipped defaults: **zero**. 0 of 50 output files changed over the 25-region gorilla
panel; strict single-copy controls stayed **0/7**, the family panel stayed **6 families / 24 copies**
(⚠ that count **excludes SHARP**, reported separately in the next clause, and is **pre-M1**; post-M1 the
same panel reads **6 families / 25 copies** excluding SHARP and **7 families / 27 copies** including it —
the M1 tier fix adds HERC2's third copy. See the 2026-08-11 correction under the `c`-sweep table),
SHARP 1f/2c. This is provable, not just observed: with `RUSTLE_ER_SENSITIVE_ONLY` defaulting true,
refine's edge set within any raw family is a *superset* of `E_r`, and connected components over a
superset of a connected graph returns the same block. ⚠ **The no-op is conditional on shipped defaults.**
It breaks whenever `E_r`'s substrate or denominator moves — `--homology-genomic-span` (measured: the
family panel reads **7f/28c** with refine vs **6f/29c** without; refine splits two recovered MAGEA copies
into a spurious 2-copy family and drops GSTM4), `RUSTLE_ER_CORE_COVERAGE` (refine passes `cores=None` and
cannot reproduce it), `RUSTLE_SHARED_EXON`, `RUSTLE_SHARED_EXON_ISOFORMS`. Anyone quoting "0 of 50 files
differ" must say **"at shipped defaults"**.

**Annotation supplies `s` only.** It does not supply `E_r`, it does not supply membership, and it does
not supply the family's name. This is the whole content of the claim "annotation is a seed, not a
definition."

### 1a. ⚠⚠ `c` IS SCALE-DEPENDENT — the edge SUBSTRATE is a free parameter of the definition

Status 2026-08-09. **Deliberately left in the code and documented instead of "fixed".** Read this
before quoting any number that came out of `E_r`.

⭐⭐⭐ **§1.1 — WHAT THE BINARY ACTUALLY EMITS (O-2, 2026-08-13, verdict
`/home/juanfra/winloci_scratch/o2_shipped_union/O2_SHIPPED_UNION.md`).** The substrate is not "a free
parameter" in the sense of being unset: **two different values ship, at two call sites, and they are
two different operations.**

| call site | who runs it | substrate | operation |
|---|---|---|---|
| `homology_edges_all_reps_pooled` | the shipped O1 catalog | **`S_x` (exon-sum)**, `homology_genomic_span` default **OFF** | ungated **global SWAP** — never both |
| `refine_families_exon_sum` | `copy_assign` (default-ON), the legacy conflict catalogs | **`S_x ∪ S_g`** | **gated, family-local ADDITION** — the same primary tier is re-run on the genomic span and unioned in, but **only** for an input copyset the exon-sum core left disconnected (`!edges_connect_all`) |

⚠⚠ **THE SHIPPED O1 `E_r` IS THE EXON-SUM RELATION** — which is the object the joint corpus calls
**"RNA"**, not the genomic span it calls "DNA". `refine` is **opt-in** on the O1 homology catalog
(`refine_enabled = refine_flag || !o1_homology`, `gw_family_catalog.rs:423`), so the additive genomic
leg is **unreachable** there: 0 `.refine.*` artifacts over 25 gorilla control regions and over the
244-rep human panel, and the emitted certificate reads `core_substrate = exon-sum`,
`additive_genomic_tier = absent (single-substrate site)`.

**Measured on the binary, human, 244 reps, one index, node sets identical across arms:**
`|E_x| = 696` vs `|E_g| = 1575`; shared 630; **exon-only 66**; genomic-only 945 ⟹ **neither contains the
other**. Block sets differ on **6 of 7** panel families [0.4869, 0.9743]; **ARI 0.4702** against a
size-matched null whose mean is −0.0001 (97.5th pct 0.0213, n = 2000).
⟹ **`E_r(S_x)` and `E_r(S_g)` are different relations, and "the substrate" is a term of the definition,
not a flag.** The `homology_genomic_span` default stays **OFF** on *unresolved sign* (O-4): recall
0.1490 → 0.1780 (P(Δ>0) = 0.9673) but precision 0.8501 → 0.9876 does **not** survive family-clustered
resampling (40.7% up / 19.9% tied / **39.4% down**, U = 18,528 labelled pairs).

**Where the union fires, and what it does.** Firing is entirely determined by what feeds `refine`.
Fed γ-quasi-clique blocks built with the *same* tier, the gate cannot open — a **selection effect in
the code**, and the reason every earlier look found a no-op: 0/31 families on the human `--refine`
panel, 0/7 on the gorilla control panel, 1/14 gate-opens on `GGO_19_mm.bam` adding **0** edges. Fed
co-location windows with no homology prefilter, it fires hard: `--cross-chrom`, human, **D1 = 26**
input copysets, gate-false **14**, **288 edges added on 322 core edges**, components 52 → 32,
partitions moved **11/26 = 0.4231** [0.2554, 0.6105]. ⭐ **But it lands on the DNA-only partition:**
block sets shipped(`E_x ∪ E_g`) vs DNA-only(`E_g`) differ **0/26** [0.0000, 0.1287] and the emitted
catalogs are **identical** (26 families / 143 copies, ARI 1.0000, 0 forbidden pairs both ways).
⚠ A degree-sequence-matched null moves the partition with **P = 0.66–1.00** per family, so *"the tier
moved a partition"* is not informative on its own — **where it lands is.**

⭐ **THE CASE THAT PROVES IT IS LOAD-BEARING.** Gorilla, `copy_assign` at shipped defaults, SHARP
`NC_086017.1:207204670-207210179`: exon-sum 2148/2323 bp → identity **1.0000**, coverage **0.0480** →
**0 edges**; genomic span 5509/33521 bp → identity 0.9811, coverage **1.0000** → **1 edge**. The family
is emitted as `CAFAM0` with **339 read assignments**, and exists solely because of the additive leg
(`[provenance] tiers=genomic-span sole-support[genomic-span=1]`). Sibling `fam1` fires the tier and is
correctly **rejected** at coverage 0.1316 — the anti-bridge guard holds. **Identity never fails;
coverage is the only failure mode. Sixth substrate, same result.**

⚠ **Scope, not decoration.** APOBEC3 carries **1** genomic edge at 3-node scope and **0** at 61-node
scope (best 61-node record identity 0.9069 at coverage **0.4971** against the 0.50 floor). ⚠ `E_r` is
also **not invariant to node NAMING**: `-X` implies `--dual=no`, so renaming reps to integer indices
moves the 27-node RNA edge set by **8 of 260 = 3.1%** on byte-identical sequence.

`c = 0.50 of min(|u|,|v|)` is a **completeness filter, not a homology criterion**, and it is not
scale-invariant. Its absolute demand is set by whatever sequence the node contributes, and the two
substrates the two shipped paths use differ by roughly an order of magnitude in length. **The DNA
catalog and the RNA catalog are therefore not the same criterion applied twice.**

Measured on ONE fixed 61-node set (7 human families — NPIP 26, MAGEA 12, TBC1D3 9, HERC2 5, GSTM 4,
APOBEC3 3, RABL2 2 — identical intervals across substrates, identical rule `id ≥ 0.60 & cov ≥ 0.50 of
min`). **Recomputed 2026-08-10 at both tiers and both coverage forms** (`o1fix_o2audit/fix/s1a/`):

| substrate | tier | cov form | id-passing within-fam pairs | survive `c` | **pair-bite** | cross-fam edges |
|---|---|---|---|---|---|---|
| genomic span | SHIPPED | old | 446 | 409 | **8.30%** | 0 |
| genomic span | SHIPPED | **M1-fixed** | 446 | 408 | **8.52%** | 0 |
| genomic span | PANEL | old | 425 | 407 | 4.24% | 0 |
| genomic span | PANEL | M1-fixed | 425 | 406 | 4.47% | 0 |
| pooled exons | SHIPPED | old | 444 | 320 | **27.93%** | 0 |
| pooled exons | SHIPPED | **M1-fixed** | 444 | 316 | **28.83%** | 0 |
| pooled exons | PANEL | old | 436 | 413 | **5.28%** | 0 |
| pooled exons | PANEL | M1-fixed | 436 | 403 | 7.57% | 0 |
| **single rep transcript — the SHIPPED RNA substrate** | SHIPPED | old | 383 | 264 | **31.07%** | 0 |
| single rep transcript | SHIPPED | **M1-fixed** | 383 | 263 | **31.33%** | 0 |
| single rep transcript | PANEL | old | 326 | 214 | 34.36% | 0 |
| single rep transcript | PANEL | M1-fixed | 326 | 214 | 34.36% | 0 |

Median `0.50·min(len)` demand, same node set: genomic span **6,338 bp**, pooled exons **3,548 bp**,
single rep transcript **876 bp** (median node length 17,413 / 10,755 / 2,071 bp).

> ⚠⚠ **THE TWO POOLED-EXON FIGURES — 28.31% AND 5.71% — WERE NEVER A CONTRADICTION. THEY ARE THE TWO
> TIERS.** 28.31% was the SHIPPED tier; 5.71% was the PANEL tier. The recomputation above reproduces
> both to within a pair or two (**27.93%** shipped vs **5.28%** panel — a **5.3×** gap, against the
> 5.0× the two archived figures implied; ⚠ **both of those are `old`-form rows — on this same panel's
> `M1-fixed` rows the gap is 28.83% vs 7.57% = 3.81×, and 3.81× is the quotable one**), and the shipped genomic-span row reproduces **exactly**
> (409/446 = 8.30%). Neither number was wrong; **both were unlabelled**, and the document presented
> them as if one measurement had given two answers.
>
> **The tier is the dominant term; M1 is not.** Across the six substrate×tier cells the M1 fix moves
> the bite by **+0.22 to +2.29 points**, while switching tier moves the pooled-exon bite by **22.7
> points**. M1 is a correctness fix (the statistic is now a fraction); B1 is what moved the numbers.

Independently on the shipped gorilla catalog (`o1_famctl/out_fixed`, 26 emitted copies) the exon-sum /
genomic-span ratio has median **0.2625** ⟹ 4,918 bp vs 898 bp = **5.5×**. A separate fixed-node run
put the ratio at ~7.5× (median spliced `seqlen/span` **0.1292**).

**That second panel, at BOTH tiers (recomputed 2026-08-11, `closure_speed/scale3_both.py`).** HUMAN,
CHM13 v2.0; 59 nodes common to all three substrates, 434 within-panel pairs, same rule
(`1−de ≥ 0.60 & cov ≥ 0.50 of min`, pre-M1 coverage form, single record). Its PANEL series is the
`4.34% → 5.71% → 29.78%` quoted above; here is what it becomes at the shipped tier:

| substrate | tier | pairs | id-passing | survive `c` | **pair-bite** | median `0.5·min(len)` |
|---|---|---|---|---|---|---|
| genomic span | PANEL | 434 | 415 | 397 | **4.34%** | 8,728 bp |
| genomic span | **SHIPPED** | 434 | 433 | 400 | **7.62%** | 8,728 bp |
| pooled exons | PANEL | 434 | 420 | 396 | **5.71%** | 4,522 bp |
| pooled exons | **SHIPPED** | 434 | 432 | 317 | **26.62%** | 4,522 bp |
| single rep transcript | PANEL | 434 | 366 | 257 | **29.78%** | 940 bp |
| single rep transcript | **SHIPPED** | 434 | 427 | 305 | **28.57%** | 940 bp |

⚠⚠ **THE TABLE ABOVE IS ON THE PRE-M1 COVERAGE FORM, AND ITS MULTIPLIER IS WITHDRAWN (2026-08-11,
adversarial D1).** `scale3_both.py` reimplements the coverage statistic **inline** as
`(qe−qs)/min(ql,tl)` instead of importing `bench/soto/rustlib.py` — the form defect **M1** condemned,
applied to a `-X` PAF, which is precisely the combination M1 was raised for. It shows: under that form
the "coverage" **exceeds 1.0 on 106/415, 173/420 and 47/366** id-passing PANEL pairs — 41% of the
pooled-exon pairs — an aligned fraction above one, in the column the delta is measured on. Recomputed
on the **same PAFs and the same 59 nodes** with the M1 axis-following form (`close_o1o2/scale3_m1.py`
→ `scale3_m1.out`; `cov > 1.0` count is **0/415, 0/420, 0/366, 0/433, 0/432, 0/427**):

| substrate | PANEL (M1) | **SHIPPED (M1)** | × |
|---|---|---|---|
| genomic span | 4.58% | **7.85%** | 1.71× |
| pooled exons | 7.62% | **26.62%** | **3.49×** |
| single rep transcript | 30.33% | **28.34%** | 0.93× |

⭐ **What survives is the finding, not the figure.** The middle substrate is still where the tier
bites, on two independent node sets — but **the pooled-exon mover is 3.49×, not 4.66×**, and the
61-node panel's own M1 rows (7.57% → 28.83%, table above) give **3.81×, not 5.3×**. The *levels* barely
move under M1 (≤ +2.3 points; the shipped pooled-exon cell is **identical**, 26.62% either way) — it is
the **ratio** that moves, because a ratio of two small rates amplifies a small absolute correction.
**Quote the level table, not the ×.** ⚠ And the companion sentence *"the panel series was not monotone"*
is **FALSE and withdrawn**: the panel series is monotone in length under **both** coverage forms
(old 4.34 → 5.71 → 29.78; M1 4.58 → 7.62 → 30.33), so monotonicity separates neither tier.

⚠ Note also that the tier raises the **id-passing denominator** on every substrate (415→433,
420→432, 366→427): the shipped tier finds *more* pairs at `id ≥ 0.60` and then loses more of them on
`c`. Both halves must be reported — a bite rate alone hides that the denominators differ.

**It is the LENGTH, not the splicing.** Control: genomic sequence truncated to rep length bites
*harder* (61.70%, and 44.41% in an earlier length-matched control) than the spliced substrate.
⚠ Centred truncation also breaks positional correspondence, so 61.70% over-states; the direction is
what is established, not the value.

Four things that must be said with it:

1. **A coverage-rule fix can reach at most 56% of the substrate effect.** Of the 144 within-family
   edges the rep substrate loses relative to genomic span, **81 (56%) die on the coverage clause** and
   **63 (44%) never produce a record at `id ≥ 0.60` at all**. No coverage rule reaches the latter.
2. **The floor buys different PRECISION on each substrate.** At the shipped 0.50 there were
   **0 of 1383 cross-family edges on any substrate** — no measured false merge either way. But the
   safety margin differs: on genomic span cross-family edges appear at `c = 0.30` (1), 0.20 (22),
   0.10 (141); on the rep transcript only at `c = 0.10` (29), zero at 0.20. **0.10 of headroom vs 0.30.**
3. **Only PAIR-level bite may be quoted.** Record-level bite on genomic span is 99.66% (132,140
   records, dominated by short repeat hits) and means nothing.
4. **"RNA and DNA return the same families given the same loci" survives, but must be qualified by
   substrate.** On this panel: 0/1383 cross-family edges either way, TBC1D3 and RABL2 partitions
   identical — but **5 of 7 families change partition when the substrate changes, in both directions**
   (NPIP rep 3 components (24,1,1) d = 0.576 → genomic 1 component (26) d = 0.923; HERC2 rep 5
   singletons, i.e. no family at all → genomic 1 component of 5 d = 0.700; APOBEC3 3 singletons → 1
   component of 3; GSTM rep 2 comps (3,1) → genomic 3 comps (2,1,1), **worse**; MAGEA rep d = 1.000 →
   genomic d = 0.939, **worse**). Every surviving component passes both γ = 0.20 and γ = 0.40, so this
   axis is γ-insensitive on this panel.

**The substrate is RECORDABLE, and quoting a catalog without it is an error.** Set
`RUSTLE_ER_EDGE_DUMP=<prefix>` (opt-in — it is **not** written by default) and the run emits
`<prefix>.params.tsv` carrying `homology_genomic_span`, `genomic_span_active`, `min_coverage`,
`coverage_denominator`, and — added 2026-08-09 for this section — `substrate_median_len_bp` and
`coverage_floor_median_bp_demand`, the run's own absolute number. Added 2026-08-10: **`coverage_form`**
(the exact numerator/denominator expression) and a **`mm_args_sensitive` row now DERIVED from
`ER_TIER_FLAGS`** rather than re-typed, so the dump can no longer describe a command the binary did not
run. State them alongside any `E_r` result.

Worked example, same BAM, same 0.50 floor, same 16 gorilla MAGEA reps, only the substrate flag changed:

| run | `substrate_median_len_bp` | `coverage_floor_median_bp_demand` |
|---|---|---|
| default (exon-sum) | 1,749 | **875 bp** |
| `--homology-genomic-span` | 4,033 | **2,017 bp** |

One coverage *fraction*, **2.31× different in bases**, from the flag alone.

#### Why this is documented rather than changed

Three alternatives were swept and none is shippable today:

- **Make `--homology-genomic-span` the default.** It does remove the length term — on identical
  intervals the DNA and RNA bite rates become the same 8.30% by construction. But (i) the residual
  **moves into locus boundaries**, which the flag does not touch: on the same 7 regions the RNA leg
  yields 56 loci to the DNA leg's 80, and **10 of those 56 (17.9%) span more than one DNA-defined
  copy** (TBC1D3: one 146,022 bp RNA locus over 4 copies; APOBEC3: 56,826 bp over 3; HERC2: three
  77–88 kb loci over 2 each), so the edge would be computed on a **multi-copy chimera** — converting a
  recall problem into exactly the false-merge mode that is still untested (D4). (ii) End-to-end on
  gorilla it is not uniformly better: controls stay 0/7, RABL2 1f/5c and APOBEC3 1f/2c unchanged,
  HERC2 1f/2c → **1f/6c** (gain), but GSTM 1f/4c → **1f/2c** (loses 2 copies) and MAGEA 2f/11c →
  **3f/13c** (family fragments). (iii) Cost is ~**450×** on the 61-node all-vs-all (1.04 s → 467.6 s,
  43 MB PAF). (iv) Of the 56 edges genomic span makes that no RNA substrate makes, the exemplar
  block's median read-supported-exon fraction is 0.747 and **5/56 (8.9%) are carried by <50% exonic
  sequence** (3 of the 5 are HERC2) — a **lower** bound, measured in the orientation most favourable
  to genomic span.
- **Add an absolute aligned-base floor (`aln ≥ N bp`).** ⚠ **REFUTED — do not re-propose.** The
  smallest `N` reaching zero cross-family edges is 5,000 bp on genomic span, 2,000 bp on pooled exons,
  500 bp on the rep transcript — a **10× spread with no shared plateau**. At genomic span's safe
  N = 5,000 within-family edges fall 446 → 375; on the rep substrate N = 500 already costs 53 edges and
  N = 2,000 leaves 67 of 383. The `cov ≥ 0.50 OR aln ≥ N` hybrid degenerates at both ends. This is the
  eighth absolute-threshold rule swept and killed on this project, and it failed the same way:
  *the hidden oracle is not the number, it is the assumption that one number exists.*
- **Union-of-records coverage** (`RUSTLE_ER_SUM_COVERAGE=1`, already implemented, default off).
  ⚠ Note this is the rule P4b's table below calls "the shipped rule" — it is not; it is default OFF.
  ⚠ **M1 changed this path's output and only this path's.** On the `sc` parity fixture, with a
  **byte-identical PAF**, the Rust edge set moved by **4 edges in each direction (95 shared of 99)**
  when the coverage form was corrected, because the unioned intervals were being measured on the query
  axis while the denominator came from the target. The other four parity fixtures (`f78`, `mix`, `rna`,
  `lg`) were unchanged. The pre-M1 fixture is archived at `erdump/sumcov/PRE_M1_2026-08-07/` as the
  evidence. **Any earlier `RUSTLE_ER_SUM_COVERAGE` measurement predates this and should be re-run
  before it is quoted.**
  Recomputed at the shipped tier on the 61-node panel it adds DNA and RNA edges at comparable rates
  (NPIP 301→325 DNA, 157→183 RNA), not the RNA-only rescue the panel-tier run reported. Keeps
  the 0.50 fraction but unions aligned query intervals across records. Of the pairs `c` kills it would
  rescue 92.1% on pooled exons and 83.8% on genomic span but only **40.7% on the shipped rep
  substrate** — on the substrate that ships, most bitten pairs are genuinely partial. Left opt-in.

⚠ **And the honest reason for the delay:** with the two-definition defect (D1) open, a code change to
this clause was not attributable — the emitted catalog was `refine(γ-QC(E_r))`, so a change to `E_r`
could be silently re-split downstream. That defect is closed as of 2026-08-09 (refine is now opt-in on
the homology path), which is what makes a future change to `c` measurable. Documenting first removes an
**undeclared** free parameter from every published number without changing any number.

## 2. Properties

### P1 — Seed-invariance (the anti-arbitrariness property)

> **Claim.** For all `s' ∈ F(s)`, `F(s') = F(s)`.

If this holds, the family is a **fixed point of the seeding operator** and does not depend on which
annotated member you happened to possess. This is the property that answers objection 1, and it is the
one that makes "seeded" different from "hand-picked".

> ⭐ **Since the §1 rewrite (2026-08-09) P1 is a THEOREM, not a panel result.** With `V` seed-free the
> blocks partition `V`, so every node lands in exactly one block and any two nodes in the same block
> return the same block. Everything below measures the **retired seeded probe**
> (`bench/crossspecies`), whose `V(s)` is a function of `s` — and on that probe strict P1 **fails**
> (4 distinct `F(s)` as sets of loci, 64/171 seed pairs agreeing). Read this section as the evidence
> that motivated removing `V(s)`, not as the current status of the property. The membership-level
> result below is what survives of it and is still the right thing to quote about annotation.

> ⚠ **Two different claims live in this section and must not be merged.** The **theorem** is about the
> **partition** of a seed-free `V` (above). The **panel** below measures **membership agreement between
> 19 seeded probes**, each of which builds its own `V(s)`. The panel is the weaker statement, it is
> about the retired probe, and it is not evidence for the theorem — the theorem needs no evidence.

**TIER: not applicable, and that is a structural fact rather than a result.** Everything in this
subsection is computed from the seed→**genome** pass (`allseeds_hits.paf`,
`minimap2 -x asm20 -c --eqx -N 200 -p 0.02`) by `seed_invariance_report.py`, which never opens an
all-vs-all PAF and never evaluates `E_r`. The O1.11 tier question (`-X` or not, asm20 leg or not)
therefore **cannot** reach any number below. Re-run **2026-08-11 on HUMAN (CHM13 v2.0)**: every figure
reproduces to the digit.

Status **of the seeded probe**: **HOLDS ON MEMBERSHIP; boundaries are not invariant.** Measured by seeding independently from
each of the 19 annotated NPIP genes (all 19 queries in one genome pass):

- **19/19 seeds recover all 19 annotated NPIP genes.** Not one seed misses a member.
- **19/19 seeds produce the identical gene set**, equal to `F(NPIPB11)`'s.
- `|F(s)|` = 26 or 27 across every seed (median 27).
- Interval Jaccard against `F(NPIPB11)` under ≥50% reciprocal overlap: median **0.963**, max **1.000**
  (7 seeds reproduce the reference intervals exactly), min **0.256**.

⭐ So `F` is a fixed point of the seeding operator **at the level of membership** — the question the
definition is actually about — and the family does not depend on which annotated member you happen to
possess.

⚠ The single low-Jaccard seed is **NPIPB8, and it is a boundary effect, not a membership failure**:
NPIPB8's annotated span is 10,633 bp, the **shortest of all 19**, and its 27 loci have median span
9,894 bp against NPIPB11's 20,314 bp. Same 27 places, same 19 genes, intervals about half as long — so
50%-reciprocal matching fails while the gene set stays identical.

> **The refined statement: membership is seed-invariant; locus extent inherits the seed's extent.**

That residual is the node-extent problem already documented elsewhere in this project, not a new defect
in the definition — but it does mean `F(s)` should be reported as a set of members, and any downstream
quantity that depends on locus *length* must not be treated as seed-independent.

### P2 — γ-independence

> **Claim.** If `F(s)` has induced density `d`, then `F(s)` is a γ-quasi-clique for every `γ ≤ d`.

Trivially true; the content is the measured margin.

**RECOMPUTED AT THE SHIPPED TIER 2026-08-11** (`bench/crossspecies/seed_family.sh`, one invocation per
species, each on its own genome; artifacts in `close_o1o2/`). Both rows were OPEN; both are now closed.

| panel | species / genome | tier | nodes | edges | density `d` | γ=0.40 | **γ=0.20 (shipped)** |
|---|---|---|---|---|---|---|---|
| human NPIP, seed NPIPB11 | **HUMAN**, CHM13 v2.0 | **SHIPPED** | 27 | **351** | **1.000** | 2.50× | **5.00×** |
| human NPIP, seed NPIPB11 | HUMAN, CHM13 v2.0 | PANEL | 27 | 351 | 1.000 | 2.50× | 5.00× |
| gorilla NPIP, seed NPIPB11 | **GORILLA**, mGorGor1 | **SHIPPED** | 25 | **260** | **0.867** | 2.17× | **4.33×** |
| gorilla NPIP, seed NPIPB11 | GORILLA, mGorGor1 | PANEL | 25 | 261 | 0.870 | 2.17× | 4.35× |

⭐ **The tier moves P2 by one edge in total, and moves no partition.** Human is **byte-identical**
between the tiers — same 351 edges, edge-set Jaccard **1.0000**, still the complete graph. Gorilla
loses exactly **one** of 261 edges, Jaccard **0.9962**, one component of 25 either way.

⚠ **And that one edge is a `c` near-miss, not an alignment failure.** The pair is
`NC_073242.2:16,322,714-16,530,307` — `NC_073242.2:28,370,269-28,506,025`: identity **0.9849 (panel)
vs 0.9833 (shipped)**, coverage **0.8709 (panel) vs 0.4989 (shipped)** against the floor `c = 0.50`.
It fails by **0.0011**. The tier changed the *extent* minimap2 chained on that pair, not whether the
two loci align. Do not read "the shipped tier is stricter" into a single 1.1‰ miss.

⚠ The **node set** (27 / 25 loci) is identical because it comes from the seed→**genome** pass, which is
`-x asm20 -c --eqx -N 200 -p 0.02` in BOTH columns and was deliberately never changed. Only `E_r`
moves. The `|V|` figures in P1, P3 and P4 are tier-invariant **by construction**, not by luck.

⚠ The two rows were also recomputed under the **legacy report rule** (`nmatch/blocklen ≥ 0.80`, pre-M1
coverage) that `seed_family_report.py` applied until 2026-08-11: it gives 351 / 1.000 for human and
257 / 0.857 for gorilla on the shipped PAF. The published PANEL rows reproduce **exactly** under both
rules (351 / 1.000 and 261 / 0.870), so the historical numbers were never rule-sensitive — **the tier
was the only open variable, and it is now closed.** That report re-implemented a PANEL-era rule and has
been switched to `bench/soto/rustlib.py`; it now prints SHIPPED and LEGACY side by side with the tier.

Human's 351 edges on 27 nodes is exactly `27·26/2`: the family is a **complete graph**, so the human
result holds for *every* γ ∈ (0,1] and no threshold was chosen. That is the strongest available form of
objection 2's answer. ⚠ It is also the reason γ-invariance is **vacuous on a complete graph** — a
density of 1.000 clears every γ trivially, so the property is only informative on the families whose
density is *below* 1.

#### Recomputed 2026-08-10 at BOTH tiers, on the §1a 61-node panel (`o1fix_o2audit/fix/s1a/`)

This panel is not a complete graph, so the margins below are the informative ones.

**Genomic-span substrate:**

| family | \|V\| | SHIPPED `d` | γ=0.40 | γ=0.20 | PANEL `d` | γ=0.40 | partition same? |
|---|---|---|---|---|---|---|---|
| NPIP | 26 | 0.926 | 2.32× | 4.63× | 0.960 | 2.40× | yes |
| TBC1D3 | 9 | 1.000 | 2.50× | 5.00× | 1.000 | 2.50× | yes |
| RABL2 | 2 | 1.000 | 2.50× | 5.00× | 1.000 | 2.50× | yes |
| APOBEC3 | 3 | — (3 singletons, 0 edges) | — | — | 1.000 (comp of 2) | 2.50× | **NO** |
| MAGEA | 12 | **0.939** | 2.35× | 4.70× | **0.697** | 1.74× | yes |
| GSTM | 4 | 1.000 (comps 2,1,1) | 2.50× | 5.00× | 1.000 (comps 3,1) | 2.50× | **NO** |
| HERC2 | 5 | 0.700 | 1.75× | 3.50× | 0.700 | 1.75× | yes |

**Single-rep-transcript substrate (the shipped RNA substrate):**

| family | \|V\| | SHIPPED `d` | γ=0.40 | γ=0.20 | PANEL `d` | γ=0.40 |
|---|---|---|---|---|---|---|
| NPIP | 26 | 0.569 (comps 24,1,1) | **1.42×** | 2.84× | 0.435 | **1.09×** |
| TBC1D3 | 9 | 1.000 | 2.50× | 5.00× | 1.000 | 2.50× |
| RABL2 | 2 | 1.000 | 2.50× | 5.00× | 1.000 | 2.50× |
| APOBEC3 | 3 | 1.000 (comp of 2) | 2.50× | 5.00× | 1.000 | 2.50× |
| MAGEA | 12 | **1.000** | 2.50× | 5.00× | **0.818** | 2.05× |
| GSTM | 4 | 0.667 (comps 3,1) | 1.67× | 3.33× | 0.667 | 1.67× |
| HERC2 | 5 | — (5 singletons, no family at all) | — | — | — | — |

#### The component-density DISTRIBUTION, both tiers, three substrates (2026-08-11, `close_o1o2/dens_both.py`)

A second, larger fixed node set — the `o1_closure` 7-family panel, **HUMAN**, CHM13 v2.0, 20 family ×
substrate cells with ≥ 2 nodes — read from cached PAFs so the two tiers see byte-identical FASTA. Both
coverage forms were run; the table shows the shipped (M1) form, with the pre-M1 value in brackets where
it differs. `d` is the induced density of the **largest** component.

| family | substrate | \|V\| | PANEL `d` | PANEL comps | SHIPPED `d` | SHIPPED comps | partition same? |
|---|---|---|---|---|---|---|---|
| NPIP | genomic span | 30 | 0.8713 | 1 | 0.8506 | 1 | yes |
| NPIP | pooled exons | 28 | 0.8968 [0.9153] | 1 | 0.6878 [0.6905] | 1 | yes |
| NPIP | rep transcript | 25 | 0.5000 | 1 | 0.6433 [0.6367] | 1 | yes |
| TBC1D3 | genomic / pooled / rep | 12 / 11 / 11 | 0.8636 / 1.000 / 0.8545 | 1 / 1 / 1 | identical | 1 / 1 / 1 | yes |
| RABL2 | all three | 2 | 1.000 | 1 | identical | 1 | yes |
| APOBEC3 | genomic span | 7 | 0.6667 [0.4000] | 3 [2] | 0.6667 [0.4000] | 3 [2] | yes |
| APOBEC3 | pooled exons | 4 | 0.5000 | 1 | — (0 edges) | 4 | **NO** |
| APOBEC3 | rep transcript | 7 | 0.6667 [1.000] | 4 | 0.6667 | 4 | yes |
| MAGEA | genomic span | 13 | 0.7436 | 1 | 0.9487 | 1 | yes |
| MAGEA | pooled exons | 12 | 0.8485 | 1 | 0.9091 | 1 | yes |
| MAGEA | rep transcript | 12 | 0.8182 [0.8333] | 1 | 0.9394 [0.9545] | 1 | yes |
| GSTM | genomic span | 4 | 1.000 | 2 (3,1) | 1.000 | 3 (2,1,1) | **NO** |
| GSTM | pooled exons | 4 | 0.6667 [1.000] | 2 (3,1) | 1.000 | 3 (2,1,1) | **NO** |
| GSTM | rep transcript | 4 | 1.000 | 3 | 1.000 | 3 | yes |
| HERC2 | genomic span | 12 | **0.3485** | 1 | 0.5303 | 1 | yes |
| HERC2 | pooled exons | 11 | **0.4545 [0.4727]** | 1 (11) | **0.3333 [0.3778]** | 2 (10,1) | **NO** |

⭐ **γ = 0.20 — the SHIPPED value — is INERT under BOTH tiers, both coverage forms, in all 20 cells:
zero components of size ≥ 3 fall below it.** The refinement cannot split anything here whichever
aligner ran. That is the γ claim that ships, and the tier does not touch it.

⚠ **γ = 0.40 — the DOCUMENTED value — binds on exactly one component per tier, and NOT THE SAME ONE.**
PANEL: HERC2 genomic span, 12 nodes, `d` = 0.3485. SHIPPED: HERC2 pooled exons, 10 nodes, `d` = 0.3778.
**`HERC2` on the spliced substrate crosses γ = 0.40 with the tier — 0.4727 → 0.3778 under the pre-M1
coverage form, 0.4545 → 0.3333 under the shipped M1 form** — i.e. it *clears* γ = 0.40 at the panel tier
and *fails* it at the shipped tier, in both coverage forms. **A γ = 0.40 verdict is not tier-portable.
γ = 0.20 is.** This is the single strongest reason to keep the shipped γ where it is.

⚠ The tier changes the **partition** on 4 of these 20 cells (APOBEC3 pooled, GSTM genomic, GSTM pooled,
HERC2 pooled) and the density on 12 — but never on the two panels P2 actually reports (§ above), where
the partition is identical on both species.

#### P2 on families OTHER than NPIP, at the shipped tier (2026-08-11, `multifamily_p1p4.py`)

**HUMAN**, CHM13 v2.0; nodes from the same (unchanged) seed→genome pass, so `|V|` is tier-invariant.

| family | \|V\| | PANEL edges | PANEL `d` | SHIPPED edges | SHIPPED `d` | γ=0.40 (P→S) | γ=0.20 (P→S) | partition same? |
|---|---|---|---|---|---|---|---|---|
| PCDHB | 22 | 119 | 0.6263 | **127** | **0.6684** | 1.57× → 1.67× | 3.13× → 3.34× | yes |
| TUBA | 11 | 10 | 1.0000 | **10** | **1.0000** | 2.50× → 2.50× | 5.00× → 5.00× | yes |
| TBC1D3 | 205 | 8,033 | 0.6390 | **8,125** | **0.6511** | 1.60× → 1.63× | 3.20× → 3.26× | yes |
| GOLGA | 226 | 7,347 | **0.3078** | **14,442** | **0.6106** | **0.77× → 1.53×** | 1.54× → 3.05× | **NO** (8 → 9 comps) |
| NBPF | 168 | 7,926 | 0.5996 (2 comps, 163+5) | **OPEN** | **OPEN** | 1.50× → OPEN | 3.00× → OPEN | **OPEN** |

⚠ **NBPF is OPEN.** Its shipped-tier all-vs-all did not finish inside the 2026-08-11 session (>28 min
of CPU on a 7.6 Mb node set and still running; the panel-tier PAFs took ~40 min in 2026-08). **Do not
fill this row from the panel-tier archive** — every other row on this table moved.

⭐ **GOLGA is the second γ = 0.40 tier-flip, and it runs the OTHER WAY from HERC2's.** At the panel tier
GOLGA's 219-node component sits at `d` = 0.3078, **below γ = 0.40, so the refinement would split it**;
at the shipped tier the same 226 nodes give `d` = 0.6106 and γ = 0.40 does not bind. Together with
HERC2 (clears → fails) this is the case for the shipped γ: **γ = 0.40 verdicts flip with the aligner in
both directions; γ = 0.20 was never crossed by any component under either tier.**

⚠ Direction, stated per panel rather than summarised: PCDHB **+8 gained / 0 lost**, TUBA **0 / 0**,
TBC1D3 **+92 / 0**, GOLGA **+7,345 gained / 250 lost** (net +7,095, edge-set Jaccard **0.4831**) —
against gorilla NPIP's **0 gained / 1 lost**. **There is no "the shipped tier is stricter/looser"
summary to be had; it is per-panel, and on GOLGA it is not even one-directional.**

⚠⚠ **The `TBC1D3` and `GOLGA` node counts are inflated by a SEED-SET defect, not by the tier.**
`multifamily_p1p4.py` selects members by `symbol.startswith(family)`, which puts **TBC1D30, TBC1D31 and
TBC1D32** — unrelated TBC-domain genes — into a 9-copy family's seed set, and antisense/cluster records
(`TUBA1B-AS1`, `PCDHB1-AS1`, `PCDHB@`, `GOLGA4-AS1`, `GOLGA7B-DT`) into the others. That is why `|V|`
reaches 205 and 226 and why P1 reads NO on all four. **It is an annotation-string artifact, it is
independent of the alignment tier, and a P1 failure on these panels must not be read as a failure of the
definition.** The script now prints every member symbol so no row can be quoted blind.

Three things this establishes:

1. **γ-independence survives at the shipped γ=0.20 with the smallest margin 2.84×**, on every family
   that forms a component at all, under BOTH tiers and BOTH substrates. Nothing here fails γ=0.20.
2. **At the DOCUMENTED γ=0.40 the margin is thin and tier-dependent**: NPIP on the RNA substrate reads
   **1.42× shipped but 1.09× panel** — a family that clears γ=0.40 by 9% under one aligner and 42%
   under another. **Do not quote a γ=0.40 margin without its tier.**
3. **The tier changes the PARTITION, not just the density**, on 2 of 7 families (APOBEC3 and GSTM, both
   genomic-span). APOBEC3 is the sharper case: the SHIPPED tier makes **no family at all** where the
   PANEL tier makes a 2-node one.

### P3 — Non-circular recovery

> **Claim.** Members recovered were not supplied.

**Re-run at the SHIPPED tier 2026-08-11; every number below is unchanged.** P3's evidence comes from two
passes that are **not** the `E_r` all-vs-all tier and were deliberately not changed — the seed→genome
pass (`-x asm20 -c --eqx -N 200 -p 0.02`) and the cross-species gorilla-loci→human-NPIP pass
(`-x asm20 -c --eqx -N 50 -p 0.05`) — so P3 is tier-invariant by construction. The one place the tier
*could* enter is "in the largest component", and the largest component is 25/25 gorilla loci at both
tiers (§ P2), so it does not.

- Human (**CHM13 v2.0**): 1 seed in, **19/19 annotated NPIP genes** out, in one component of 27. The
  other 18 were never given. 8 of the 27 loci overlap no annotated NPIP gene.
- Gorilla (**mGorGor1**): the NCBI mGorGor1 annotation names **exactly one** NPIP gene. From it, 25 loci; all 25 overlap
  annotated genes spanning **37 distinct symbols of which exactly 1 is named** (the seed, NPIPB11) — the other 36
  are `LOC…` placeholders. Validation is against **human** annotation, a source independent of the
  gorilla annotation the seed came from: 23/25 loci carry a qualifying edge to a human NPIP gene
  (median identity 0.956, range 0.901–0.973), and all 19 human NPIP genes have ≥1 qualifying gorilla edge.

### P4 — Substrate convergence

> **Claim.** `F(s)` computed from the genomic substrate and from the RNA substrate agree.

Three probes built at the same single seed locus. **HUMAN, CHM13 v2.0.**

**TIER: not applicable.** All three columns are outputs of the probe→**genome** pass
(`-x asm20 -c --eqx -N 200 -p 0.02`), which is unchanged in both tiers; no all-vs-all PAF and no `E_r`
evaluation enters this table. Confirmed by re-running `seed_family.sh` at the shipped tier on
2026-08-11: the genomic probe still returns **27** candidate loci (gorilla still **25**), unchanged.

| probe | size | genes recovered | loci |
|---|---|---|---|
| RefSeq NPIPB11 transcript | 6,117 bp (8 exons) | 9/19 | 21 |
| **IsoSeq read-supported exons** (≥3 reads) | **40,632 bp (41 blocks)** | **19/19** | **27** |
| full genomic span | 25,154 bp | 19/19 | 27 |

The RNA probe and the genomic probe return the **same 27 loci**, including the same 8 non-NPIP-annotated
extras at near-identical coordinates. The failing probe is the *annotated transcript*, which represents
6.6× less of the locus than the reads do. **RNA and the T2T genome converge; annotation is the outlier.**

### P4b — DNA graph vs spliced-RNA graph, same node set

Both graphs over the identical 27 nodes of `F(NPIPB11)`; only the sequence each node contributes
changes (genomic span vs read-supported exons at ≥3 IsoSeq reads, from the full A119b BAM).

⚠⚠ **THE RULE NAMED IN THIS TABLE IS NOT THE SHIPPED RULE — corrected 2026-08-10.** The line below
used to read "Under the shipped two-tier rule (asm20 ≥ 0.80 ∪ sensitive ≥ 0.60, coverage ≥ 0.50
aggregated over records)". Three of those four clauses are wrong about what ships:

- **"two-tier"** — the shipped default is **sensitive-only** (`RUSTLE_ER_SENSITIVE_ONLY`, true since
  2026-08-07). The asm20 column below is a tier the binary skips.
- **"aggregated over records"** — that is `RUSTLE_ER_SUM_COVERAGE`, **default OFF**, and it carries a
  measured cross-family precision cost (0.9038; 45 cross-family edges vs **zero** for single-record).
  The shipped rule is **ANY ONE record** clearing both floors.
- the run was **PANEL-tier** (no `-X`).

The table is retained as the record of what was measured, **relabelled**, not deleted:

**PANEL tier, two-tier floors, coverage aggregated over records — NOT the shipped rule:**

| | asm20 | sensitive | union | density |
|---|---|---|---|---|
| DNA | 351 | 351 | **351/351** | 1.000 |
| RNA | 231 | 321 | **338/351** | **0.963** |

**Jaccard 0.963; 13 DNA-only edges; ZERO RNA-only edges** — *under the rules tested there*, which did
not include the shipped one. The residual 13 are coverage near-misses (0.039–0.481, most ≥0.28), not
identity failures.

#### The 27-node set AT THE SHIPPED TIER + SHIPPED single-record rule (2026-08-11) — this row is now CLOSED

`bench/crossspecies/graph_vs_graph.sh`, **HUMAN** (CHM13 v2.0, RNA nodes from the **full** A119b BAM,
`-F 2308`, exons at ≥3 reads). Same 27 nodes; **all 27 have read support**, so nothing is dropped.

| | tier | rule | DNA E | RNA E | shared | DNA-only | **RNA-only** | Jaccard | RNA density |
|---|---|---|---|---|---|---|---|---|---|
| 27-node human NPIP | **SHIPPED** | single-record, 1−de ≥ 0.60, M1 coverage | **351/351** | **260/351** | 260 | **91** | **0** | **0.741** | **0.741** |
| 27-node human NPIP | PANEL | two-tier, coverage aggregated | 351/351 | 338/351 | 338 | 13 | 0 | 0.963 | 0.963 |

⭐ **`E_RNA ⊂ E_DNA` SURVIVES the shipped tier on this node set: RNA-only edges = 0.** So the
qualification recorded below for the 61-node panel (13 RNA-only) is a **node-set** effect, not a tier
effect — on the original 27 nodes the containment is exact under the rule that ships. Both graphs are
still ONE component of 27 with no singletons; the **partition is unchanged**, only the density is.

⭐ **Identity never fails; coverage is the only failure mode — at BOTH tiers.** Of the 91 DNA-only
edges, **91/91 are coverage misses at identity 0.9709–1.0000** (best id-passing coverage 0.2928–0.4989,
median 0.4384; 70 of 91 are ≥0.40, and the maximum is again **0.4989**, 0.0011 under `c`). Zero are
identity failures and zero lack a record. This is the same shape as the panel's 13, 7× more of it.

⚠ **What moved is the count, and it moved a lot: RNA 338 → 260 edges, Jaccard 0.963 → 0.741.** The
headline "Jaccard 0.963" is a PANEL number computed with coverage **aggregated over records** (default
OFF). Quote **0.741** for the shipped rule at the shipped tier.
⚠⚠ **THE CONTAINMENT IS TIER-ROBUST BUT NOT PANEL-ROBUST — this instruction is CORRECTED (2026-08-14).**
It previously read "quote the containment (0 RNA-only) as the tier-robust part", which invited quoting a
figure that does not generalise: on the 80-node `o1_closure` panel there are **33 RNA-only edges of 481
(6.9%)**. Quote instead the claim that survives both tier and panel — **only 1 of those 33 joins two DNA
components**, i.e. containment holds at the PARTITION level (6/7 families), not at the EDGE level. §1★.3.

⚠⚠ **A defect in the diagnostic was found and fixed en route (2026-08-11).** `graph_vs_graph_report.py`
fed `why()` the output of `paf_pairs(..., min_identity, min_coverage)`, whose floors are applied *per
record inside* the parser — so every failing pair was simply absent from the dict and the breakdown
printed "**no alignment record at all: 91**" while in fact **all 351 pairs have records**. Its other two
branches were dead code and would have raised (they indexed a dict as a tuple). The breakdown now comes
from an unfiltered `diagnose()` scan. **Any earlier "no alignment record at all" line from this script
is void.**

#### The same comparison at the SHIPPED tier + SHIPPED single-record rule (2026-08-10)

Run on the §1a 61-node panel (a different, overlapping node set — 26 NPIP nodes, not 27), so this
**replaces nothing above**; it is the first measurement of this shape under the rule that ships:

| family | \|V\| | DNA E | RNA E | shared | DNA-only | RNA-only | Jaccard | partition same? |
|---|---|---|---|---|---|---|---|---|
| NPIP | 26 | 301 | 157 | 150 | 151 | 7 | 0.487 | no |
| TBC1D3 | 9 | 36 | 36 | 36 | 0 | 0 | 1.000 | **yes** |
| RABL2 | 2 | 1 | 1 | 1 | 0 | 0 | 1.000 | **yes** |
| APOBEC3 | 3 | 0 | 1 | 0 | 0 | 1 | 0.000 | no |
| MAGEA | 12 | 62 | 66 | 62 | 0 | 4 | 0.939 | **yes** |
| GSTM | 4 | 1 | 2 | 1 | 0 | 1 | 0.500 | no |
| HERC2 | 5 | 7 | 0 | 0 | 7 | 0 | 0.000 | no |
| **TOTAL** | **61** | **408** | **263** | **250** | **158** | **13** | **0.594** | 3/7 |

> ⚠ **"ZERO RNA-only edges under every rule tested" does not survive contact with the shipped rule on
> this panel: there are 13.** They are concentrated in NPIP (7) and MAGEA (4). This is a *qualification*
> of P4b, not a refutation of it — different node set, and the 27-node claim stands unretested — but
> the sentence "under every rule tested" must now be read as "under the rules tested in 2026-08, none
> of which was the shipped single-record rule at the shipped tier."

The aggregated rule's effect at the shipped tier, for reference (single → aggregated): NPIP 301→325
DNA and 157→183 RNA; APOBEC3 0→3 DNA; GSTM 1→3 DNA and 2→3 RNA; TBC1D3, RABL2, HERC2 unchanged.
**Aggregation adds DNA edges roughly as fast as RNA edges here**, which is the opposite of the
"aggregating recovered 53 RNA edges and 0 DNA edges" result recorded below — another number that
belongs to the PANEL tier.

### P4c — IS THE DNA/RNA JOINT OBJECT DEFINITIONAL? **NO — IT IS A PROPERTY** (closed 2026-08-13)

Full verdict `/home/juanfra/winloci_scratch/o1_joint/O1_JOINT_VERDICT.md`; the definition that was
tested, its pre-declared falsifiers and the phase-3 record are in `docs/o1_investigations.md#the-joint-dna-rna-family-definition-retracted`
(§9). The question asked was whether the DNA/RNA joint object should become **a clause of O1's
definition**. It should not.

⭐⭐⭐ **THE HONEST SENTENCE.** At the shipped tier the joint object's partition is identical to the
DNA-only partition on **0 of 5** disagreement families (Wilson 95% [0.000, 0.434]; 0 of 7 overall),
its only clause with teeth is a **label** (`κ` fires DISCORDANT on 4 of 6 opportunities but scores
**Fisher p = 0.40** against an outcome it did not produce, and at n = 6 the test **cannot reach**
0.05), and the node-level complementarity that motivated the proposal collapses on re-derivation —
**RNA-only nodes 24 → 8 of 362** Soto members (ceiling **13 = 3.6%**), with **all 13 carrying a
genome-wide DNA paralog at `1 − de` of 0.9703–1.0000 against a 0.60 floor and 0 below it** ⟹ **report
the DNA/RNA relation as a PROPERTY of this definition, never as a clause of it.**

⚠⚠ **RE-QUOTE THE FIRST CLAUSE (O-2, 2026-08-13).** *"0 of 5 / 0 of 7"* is an **entailment**, not a
measurement: the joint object's `E_r` is a DNA object by construction, so its partition equals DNA's by
construction (§9.1 says so). Apply the **shipped gated-union rule** to the joint run's own edge sets and
the forbidden set is **not** empty — **2/5 and 2/7 (APOBEC3, GSTM)** at the rep-transcript substrate,
**0/5 and 0/7** at pooled read exons. On the binary's own catalogs the union is a no-op **relative to
DNA-only** (0/26 [0.0000, 0.1287], catalogs identical, ARI 1.0000) and **never relative to the shipped
core** (11/26 partitions moved). The verdict stands; this sentence does not, as written. Full
re-measurement: `/home/juanfra/winloci_scratch/o2_shipped_union/O2_SHIPPED_UNION.md`; substrate
consequences in §1.1 above.

⚠⚠ **RETRACTED, DO NOT QUOTE AGAIN.** *">= 95%: DNA 90.3% / RNA 84.5% / union 97%, both 282, DNA-only
45, RNA-only 24, neither 11 — different questions, NOT a ranking; **RNA is NOT a subset of DNA**"*
(2026-07-25/26). It is the product of the **banned `nmatch/blocklen` identity metric** thresholded at a
post-hoc 0.95, an unsymmetrised query-axis coverage test, a pre-centralisation tier carrying the
condemned `-p 0.5`, an RNA leg that is a four-way union of detection legs (only 237/306 = 77% the O1
catalog) computed on a **`-M -L` SUBSET BAM**, and a denominator containing 7 members for which the DNA
predicate is vacuously false. On all 1,644 old-tier non-self records the metric bias is **strictly
one-sided** — 361 (22.0%) pass `1 − de` at 0.95 and fail `nmatch/blocklen`, **reverse count 0**; at the
shipped 0.60 floor **1,644/1,644** pass. There is **no divergence-rescue stratum**.

⭐ **THE SHIPPED-TIER RE-DERIVATION RAN** (`o1_joint/lens_vacuity/`): `minimap2 -c -X --no-long-join
-t 4 -k 11 -w 5 members.fa members.fa`, 1,659 s wall / 5,761 s CPU, 4.00 GB peak RSS, 1,306,867 usable
non-self records, ONE index over 11,264,772 bp / 588,643 distinct minimizers / `mid_occ` 752 — **single
scope; may not be set beside any genome-scale number.**

| leg (n = 362 Soto members, HUMAN / CHM13 v2.0) | DNA | RNA | union | both | **RNA-only** | DNA-only | neither |
|---|---|---|---|---|---|---|---|
| PUBLISHED 2026-07-25 (banned metric, `-p0.5`) | 90.3% | 84.5% | 97.0% | 282 | **24** | 45 | 11 |
| old tier `mm.paf`, shipped rule, symmetrised | 97.0% | 84.5%\* | 98.6% | 300 | **6** | 51 | 5 |
| genome map-back `gmap.paf`, shipped rule | 95.0% | 84.5%\* | 98.6% | 293 | **13** | 51 | 5 |
| ⭐ **SHIPPED tier + shipped rule** | **96.4%** | 84.5%\* | **98.6%** | **298** | **8** | **51** | **5** |
| same, any-member predicate | 97.5% | 84.5%\* | 98.6% | 302 | **4** | 51 | 5 |

\* ⚠ **the RNA column was NEVER re-derived** — inherited 2026-07-25, `-M -L` subset BAM, four-leg
detection union. It is identical in every row because it never moved. Label it INHERITED on sight.

⚠⚠ **AND THE PREDICTED DIRECTION OF THE TIER CORRECTION WAS REFUTED.** Two independent phase-1 reports
asserted the shipped tier "can only ADD records, so DNA can only rise and RNA-only can only fall."
Measured: **DNA FALLS 351 → 349 and RNA-only RISES 6 → 8.** `-X`/`--no-long-join` is not a superset of
`-c --eqx -N50 -p0.1` — suppressing long joins shortens the aligned span, so pairs drop below the M1
0.50 coverage floor. GOLGA8A is rescued; AC137800.1, POM121 and POM121C are newly lost, **all three on
COVERAGE** (best id-passing coverage 0.4460 / 0.2875 / 0.2875) at identity 0.9314–0.9867. **Quote 8.**
Nobody may reuse a one-sided-correction argument for this tier change.

⭐ **THE CEILING — the RNA leg cannot rescue the claim.** Re-deriving it can only move members between
cells, so RNA-only is bounded by `RNA-only + neither` = **13 of 362 (3.6%)**. All 13, recomputed as
`1 − de` (the stored `member_genome_paralogs.tsv` used the banned `matches/alnlen`): AC119751.3 0.9970,
AC137800.1 0.9970, AC243829.6 1.0000, AL590399.2 1.0000, BOLA2B 1.0000, CASTOR2 0.9949, CR381670.1
0.9975, CTSLP3 0.9927, DEFB104B 0.9927, DUX4L50 0.9790, HERC2 0.9703, POM121 0.9832, POM121C 0.9881 —
**min 0.9703; below the 0.60 floor: 0; with no non-self record: 0.** ⟹ **the set "loci RNA sees and DNA
cannot" is EMPTY** by a margin of 0.37 in identity, against every reachable RNA leg.

⭐⭐ **THE 27-NODE REDUNDANCY NEGATIVE HAS MEASURED POWER ZERO, AND MAY NEVER AGAIN BE QUOTED AS
EVIDENCE THAT RNA'S EDGE LOSSES DO NOT MATTER.** Pre-declared null (`o1_joint/test/arm2/PREREG.md`,
hashed before scoring): delete `|E_g \ E_x|` edges uniformly at random from `E_g`, 10,000 draws, seed
20260812. On the 27-node panel — which is the **COMPLETE** graph, 351 = 27·26/2, edge connectivity 26 —
a random deletion of 91 edges changes the partition **0/10,000 times [0.00000, 0.00038]**, so "0 of the
91 are cut edges, the partition is identical" is true and **uninformative**. ⭐ On **61-node NPIP**,
where the test has power, the *RNA-determined* deletion of 151 of 301 edges **does** disconnect the
family while a random deletion of the same size does so with probability **0.0002** (permutation
p ≈ 2 × 10⁻⁴): **RNA's losses are not a random thinning** — the two nodes intersection cuts off are
`L~chr16_22777929_22814315` (NPIPB5, 293 spliced reads, DNA degree 23/25, **RNA degree 0/25**) and
`L~chr16_75785926_75819336` (1,243 reads, DNA 23/25, **RNA 0/25**). ⚠ **Exactly ONE of seven family
cells is an informative null** (HERC2's is degenerate at P = 1.0000 by identity; five have zero lost
edges). ⟹ **the P4b line "the partition is unchanged, only the density is" is TRUE and CARRIES NO
WEIGHT on the 27-node set**, and at 61 nodes the intersection fold **destroys HERC2** (5 nodes → five
singletons) and costs NPIP two members.

⭐⭐ **THE 91 LOST EDGES ARE NOT ABOUT RNA, AND NOT ABOUT NODE SIZE — THEY ARE ABOUT CONTIGUITY.**
Build `T|locus` = a **contiguous genomic window of the DNA node whose length is EXACTLY the RNA node's
length** (asserted 27/27): same locus, same length, **DNA content, no RNA anywhere**. At the **literal
shipped tier**, `T × T` gives **351 edges, symdiff 0 vs stored `E_dna`, 91/91 of the lost edges
recovered**, coverage on the 91 of **0.8032 / 0.9913 / 1.0000**. A joint `R × D` arm on a *more
permissive* translation gets only 344 / 88 of 91 **and destroys 4 of the 260 edges `E_rna` already
had**. Per-node Spearman of median T-coverage against the rna/dna length ratio is **+0.118** (n = 27);
the most-shortened node (ratio 0.598, 16,428 → 9,828 bp) has median T coverage **1.0000**; an 18% loss
of total sequence costs **zero** edges. Controls: `D.fa` ava at the literal shipped tier reproduces
`E_dna` (351, symdiff 0), `R.fa` ava reproduces `E_rna` (260, symdiff 0). ⟹ **the corrected unifying
statement for the whole coverage-clause literature: an `E_r` edge needs at least ONE endpoint that is
CONTIGUOUS over the shared region — over-long/fused nodes are free (they hand the denominator to the
partner), short-but-content-complete nodes still pass iff the partner is contiguous, and RNA×RNA is the
one configuration where NEITHER endpoint supplies contiguity.** DNA has contiguity for free.

⚠⚠ **A DEFECT IN THIS SPEC, FOUND WHILE ENUMERATING THE SURFACE — §1 UNDERSPECIFIES THE SUBSTRATE, AND
THE BINARY USES THREE.** `denovo_pipeline.rs:3388-3425` (`family_refine`, the O1 refinement path) runs
its core tier on the **exon-sum (spliced)** sequence (`include_introns` defaults false, `:3231`) and
then, gated only on `!edges_connect_all(...)`, **unions in an additive GENOMIC-SPAN tier, DEFAULT ON**
(`TIER_GENOMIC`, `:3418-3421`) — so the shipped `E_r` at that site **is `E_x ∪ E_g`**. Its only cited
justification (`:3409`) is `bench/FALSE_NEGATIVES.md`, **which was DELETED from the tree in `9b0814f`**
(⚠ corrected 2026-08-13 by O-4 — the earlier wording "does not exist anywhere in the tree" wrongly read
as never-written; it was committed with the tier in `4586ba8`, is now RESTORED verbatim from the git
blob, and a test asserts the citation resolves).
The other call site, `homology_edges_all_reps_pooled` (`:4520-4535`), makes the same decision the other
way: `homology_genomic_span` is **default OFF** (`:3243`).
⭐⭐ **O-4 (2026-08-13) DECIDED BOTH DEFAULTS: KEEP THEM.** They are **not the same operation** — refine's
is a **gated, family-local ADDITION** (it fires only where the exon-sum core leaves a copyset
disconnected and is confined to that copyset, so it can never bridge two families), the catalog's is an
**ungated GLOBAL SWAP** that re-decides every pair (+879 edges, +126% on the human panel). Refine's leg
moved **11/26** examined families' partitions yet the emitted catalogs are **identical** to DNA-only
(block sets differ **0/26 [0.0000, 0.1287]**, ARI 1.0000, 0 forbidden pairs). Flipping the catalog's
default was **rejected on unresolved sign, not measured harm**: recall rises (P(Δ>0) = 0.9673) but
precision does not survive family-clustered resampling (40.7% up / 39.4% down), and ⚠⚠ **the same
contamination formula on the same two partitions flips sign with the choice of truth** (+0.0310 worse
under RefSeq-coarse, −0.0377 better under the panel roster). Denominators and CIs:
`bench/FALSE_NEGATIVES.md`. And the O1↔O2 rule certificate
(`er_rule_rows`, `:3701-3740`) wires `substrate` to `params.include_introns`, so it **prints
`substrate = exon-sum` on a run that unioned genomic-span edges in**, and has **no
`additive_genomic_tier` row at all**. ⟹ **§1 must name the substrate as part of the definition**, and
every "union is a no-op" measurement in this project was made on externally built `gvg_dna.paf` /
`gvg_rna.paf` graphs, **not on the shipped union**. Recorded, **not fixed**.

⭐ **THE CERTIFICATE HALF OF THAT DEFECT IS NOW FIXED (O-3, 2026-08-13) — AND FIXING IT FALSIFIED A
"DONE" ROW.** `er_rule_rows` no longer wires the substrate to `params.include_introns` alone:

* the key `substrate` is **renamed `core_substrate`** (a run whose edge set is `E_x ∪ E_g` has no single
  substrate, and the old name invited exactly the misreading it produced);
* a new **RULE** row `additive_genomic_tier` states the leg's armed state and, when off, **why** —
  `armed (genomic-span, …)` / `off (core substrate is already genomic-span)` / `off (no genome reachable
  at this call site)` / `absent (single-substrate site)`;
* new **DATA** rows in `<prefix>.refine.params.tsv` say whether it FIRED and what it contributed:
  `n_families_core_unconnected`, `n_families_genomic_tier_ran`, `n_families_genomic_tier_added_edges`,
  **`n_edges_genomic_tier_added`** (the quantity every "the union is a no-op" claim is about, on the
  object the binary actually emits), `n_edges_core_tier`, `additive_genomic_tier_fired`;
* the arming decision is taken **ONCE**, by `additive_genomic_tier(params, genome.is_some())`, and the
  gate and the certificate both read that one value (the `ER_TIER_FLAGS` precedent; a source-level test
  fails if a second derivation appears);
* the per-alignment dump is disambiguated too — the core and additive runs use the same command on
  different sequence, so their `.paf` / `.args` were **byte-identical but for a counter**. They are now
  `…refine.core.exon-sum.*` / `…refine.additive.genomic-span.*` and `.args` carries a `tier` line.

⚠⚠ **CONSEQUENCE — RE-QUOTE ROW X.2.** `diff <prefix>.rule.tsv <prefix>.refine.rule.tsv` was measured
EMPTY on all 5 refine calls; it is empty **no longer**, and was empty only because the file could not see
the substrate. At shipped defaults with a genome reachable the diff is now **exactly one line**,
`additive_genomic_tier`: `absent (single-substrate site)` on the O1 catalog vs `armed (genomic-span …)`
on O2's refine. The two sites agree on every other edge-deciding knob — **and disagree on the substrate**,
which is the F6 finding, now visible in the instrument built to find it. ⚠ Default **unchanged**: this is
an observability fix only; deciding `homology_genomic_span` / the additive default is still open (O-4).

⚠⚠ **AND `O1 ⊥ O2` IS ALREADY FALSE AT THE NODE ON THE DEFAULT CATALOG PATH.**
`detect_homology_catalog_genome_wide` (`:2964`) calls `locus_unique_mapper_counts` — defined in
**`read_conflict.rs:267`, the `E_c` module** — and feeds `distinguishing_uniq` to
`distinct_locus_reps(..., cfg.locus_min_reads())` (`:3116`, **not** behind `--refine`), where
`reads_distinguish` is documented as *"the χ(H) edge predicate restricted to a co-located pair"*;
`cfg.locus_min_reads()` (`:172`) literally `return self.conflict.min_reads` — **O2's
`RUSTLE_CONFLICT_MIN_READS`, default 3, renamed.** The guard test
`homology_catalog_never_touches_the_conflict_graph` (`:5756`) bans four strings and `cfg.conflict.`
and **passes while all of that flows** — it certifies spelling, not semantics, exactly as
`--cross-chrom` did on 2026-08-11. ⚠ Near-inert on chr1 (the same-strand MAPQ rule fired **0 times over
451 loci**), but that is a property of the data. **`E_r` itself is clean** — sequence-only, reference
bases at read-derived coordinates, no read base and no assignment in any comparison. **The violation is
entirely in `V`.** Prerequisites for any joint definition are listed in `o1_investigations.md`
§9.7 (O-5).

### ⭐ RESOLVED 2026-08-13 — SCOPED, MEASURED AND MONITORED, NOT REMOVED

**The claim is now stated with its exception rather than overstated.** O1 decides membership by
**SEQUENCE alone, EXCEPT** for a co-located, **same-strand** pair of reps that the junction-asymmetry
rule did not settle, where the unique-mapper count decides. Everything else — the identity floor, the
coverage clause, `E_r`, the γ-refinement — is sequence-only.

**Why it was not removed.** The branch's real domain is a pair of overlapping same-strand
**junction-less** reps: both advisor-flagged regression tests (`..._keeps_distinguishable_colocated_
copies_separate`, `..._two_junctionless_copies_still_decided_by_reads`) build their inputs with the
`rep()` helper, which sets `introns: vec![]`. Such a pair has **no splice structure to compare**, and
its spans overlap so the genomic sequence is **shared by construction** — nothing in the sequence
separates the two. This is the true K = 0 case, and the unique-mapper count is the only signal that
exists. Deleting the branch re-introduces the advisor-flagged over-merges.

⚠ **And it must not be "fixed" by comparing intron chains on the both-spliced side.** Two overlapping
same-strand reps with different chains are usually **isoforms of one gene**; splitting them would
over-split loci. Because the branch never fires on any substrate on record, such a change **cannot be
validated** — it would be judged on node-level metrics, which has failed end-to-end three times here.

**MEASURED — the exception is reachable but inert, on two species.** `RUSTLE_LOCUS_AUDIT=1`:

| substrate | co-located pairs | settled by junction asymmetry | settled by antisense ratio | **decided by the read leg** |
|---|---|---|---|---|
| 25-region control panel (gorilla, all 6 families) | 15 | 15 | 0 | **0** |
| 19-region family panel (gorilla, 14 fam / 39 copies) | 94 | 89 | 5 | **0** |
| chr1 (human, 451 loci) | — | — | — | **0 firings** |

The audit is behaviour-neutral: 75/75 control-panel outputs byte-identical to `out_m1fix`.

**MONITORED — silence is now an active claim.** `distinct_locus_reps` counts every verdict that came
from `reads_distinguish` and, if the count is non-zero, prints an **unconditional** stderr line
(`[o1-perp-o2] WARNING: …`) stating that the node set is not a function of sequence alone for that run
and must be disclosed with any number derived from it. Zero on every substrate to date.

**THE GUARD NOW CHECKS SEMANTICS, NOT SPELLING.** `homology_catalog_never_touches_the_conflict_graph`
was rebuilt: its checked set is **derived from the `use super::read_conflict::{…}` import** instead of
a hand-picked four-string list, so a new `E_c` symbol cannot appear unseen; it scans **transitively**
into `distinct_locus_reps`, because the body-only scan was structurally blind to `reads_distinguish`
(which is not in the catalog body at all); the two known uses are pinned in a `DISCLOSED_EC_USES`
ledger whose comment records that **adding to it is a thesis edit, not a code edit**; and
`locus_min_reads()` is pinned as an **alias of `conflict.min_reads`**, which is precisely why the
pre-existing `cfg.conflict.` assert passed while the dependency flowed. Red-before verified: injecting
a `family_mapq0_support` call into the catalog fails the guard. **771 passed / 0 failed / 11 ignored**,
real exit code captured to a file; control panel 75/75 byte-identical.

⚠ **"THE RNA GRAPH" IS NOT ONE OBJECT — the biggest caveat on every RNA edge number in this file.**
Same 61 nodes, same rule, same tier: **rep transcript 263 edges (the shipped substrate) vs pooled read
exons 316**, and the RNA partition differs between the two on **4 of 7 families**. Under pooled exons,
RNA's partition equals DNA's on **6 of 7** and union is a **clean no-op (0 of 7)**; intersection differs
on **1 of 7 (HERC2 only)**. **The only finding surviving both substrates is that intersection damages
HERC2.** Every RNA edge count in P4b and P4c is a **rep-transcript** number.

⚠ **What the joint object is still good for.** For **O3** it is *necessary for the objective to be
stateable at all* — a missing copy is a property of the (reads × genome) **product**, and neither a
genome-free nor a read-free statistic can express it — while buying **no measured power**. For **O2** it
is byte-identically inert: `K` enters only as `thr = alpha/max(K−1,1)`, and replaying five shipped
gorilla panels with evidence held fixed, **K ± 1 moves 0 of 7,691 reads (0.000%)** because the band
`1e-4 < p ≤ 1e-3` that `K = 2..11` can move `thr` across contains **0 reads** (evidence is bimodal:
977 at `p = 1.0`, 6,506 at `p ≤ 1e-20`). **What should be adopted is not a fold but a SUBSTRATE
decision** — assembly supplies every base and therefore the whole edge relation; RNA supplies which
loci exist as expressed objects and their transcribed extent; minimal annotation supplies candidate
boundaries — which is already half-implemented at two call sites with opposite defaults, and is
measured at **+176 true / +13 false edges at id ≥ 0.90, precision 0.908 → 0.916 (UP)** on 540 reps.

⚠⚠ **CORRECTION.** An earlier version of this comparison applied a uniform 0.80 floor to BOTH tiers and
reported RNA density 0.547 (single-record) / 0.698 (aggregated). That was wrong:
`homology_refine_params` forces both tiers equal **only when `--min-identity` is passed explicitly**;
with `--refine` alone the floors are asm20 0.80 / sensitive 0.60 (`denovo_pipeline.rs:2673-2683`).
**Do not quote 0.547 or 0.698.** The sensitive tier does the work — 321 of 338 RNA edges.

⚠ Two separate artifacts were found and fixed en route, both of which inflated the apparent gap:
the **single-record rule** penalises concatenated exons specifically (aggregating over records
recovered 53 RNA edges and **0** DNA edges), and the **uniform identity floor** above.

### Variation graph: RETRACTED (2026-08-08)

⚠⚠ **This section previously reported that the 27 loci "collapse to a single ~16.1 kb minigraph path
with no structural variation", that the spliced path is a "strict sub-path" of the genomic path, and
that this is *why* `E_RNA ⊂ E_DNA`. All of it is RETRACTED.** minigraph reported `inserted 0 events`:
nothing from the 26 sample loci entered the graph at all, so the "16,118 bp path" was simply the
REFERENCE LOCUS ECHOED BACK, and the "no structural variation between copies" finding was a graph
containing exactly one copy. The downstream use of it as an independent confirmation of the
size-invariant ~16 kb cassette is withdrawn with it.

The cassette result itself does not depend on this — it rests on the pairwise block analysis, which
stands on its own evidence (see `project_absolute_block_floor`). What is gone is the claim that a
variation graph confirmed it independently, and the path/sub-path account of `E_RNA ⊂ E_DNA`.

⟹ Establishing anything about intron bubbles needs a **base-level `vg construct`**, not minigraph:
minigraph only records variants above a length threshold, so it is structurally incapable of speaking
to the SNVs and small indels where the identity signal actually lives.

### On lowering the identity floor below 0.80

The shipped pipeline **already does**: the sensitive tier is 0.60 by default. Within this family, RNA
density by floor (coverage ≥ 0.50, aggregated): 0.80 → 0.698, 0.75 → 0.840, 0.70 → 0.920,
0.60 → 0.983, 0.50 → 0.997.

⚠⚠ **This sweep CANNOT justify lowering the global floor.** All 27 nodes are one family, so there are
no out-of-family negatives and rising density is guaranteed by construction — the tautology trap.
Lowering the global floor was tested genome-wide and produced **one 776-copy hairball at id 0.70**
([[project_rna_grouping_identity_floor]]).

### None of P1–P4 uses read conflict — on purpose

P1, P2, P3 and P4 are properties of `E_r` alone. No multimapping information enters any of them. P4 uses
IsoSeq reads only to *construct a probe sequence*; that probe is then aligned and thresholded by the same
`E_r` rule, and no read-conflict relation is consulted.

The reason is stated positively rather than as a prohibition:

> **`E_c` is a held-out test of `E_r`, not a term in it.**

`E_r` is built from sequence. `E_c` — which loci a read cannot be placed between — is then measured
independently and asked whether it agrees. Dense `E_c` inside an `E_r` family is confirmation. Dense
`E_c` *between* two `E_r`-separate families is a warning that `E_r` may have split something (and is
exactly what surfaced NBPF–NOTCH2NL and NPIP–SULT1A below). Feeding `E_c` back into `E_r` would be
training on the test set: the agreement becomes guaranteed and stops being evidence of anything.

This is the same separation as O1 ⟂ O2, stated as a measurement discipline rather than a rule.

## 3. The certificate (and why it is not the definition)

`E_c`, the read-conflict relation — pairs of loci a read cannot be placed between — is reported
**alongside** the family, never inside it.

Human NPIP, 21,775 reads: 64.6% place alignments on more than one NPIP gene (median 7, max 19); 97.4%
carry a MAPQ=0 alignment; **all 171 of 171 possible gene pairs are linked by ≥1 shared read**, with the
top pairs at 8.2–9.5k reads.

> ⚠ **`E_c` must never enter `E_r`.** O1 (identify, support) and O2 (assign, ambiguity) are independent by
> construction; using read conflict to build families collapses that independence and destroys the K=0
> argument, which reasons about ambiguity *within* a family already defined without it.

The certificate's job is to answer "does the method know these belong together?" with evidence rather
than assertion — and its honest form is **abstention**: these reads could not be placed, and we say so.

**Specificity is the load-bearing control.** A certificate that lights up between unrelated families as
readily as within one is vacuous. Measured on `human_testis.t2t.bam` (9.1M reads, full-genome alignment,
no region restriction — an *independent sample* from A119b), 10 families, 135 genes, 1,312 within-family
pairs:

| shared-read threshold | within-family linked | between-family linked | enrichment |
|---|---|---|---|
| ≥1 | 72.2% | 1.1% | **67×** |
| ≥3 | 69.6% | 0.4% | **179×** |
| ≥10 | 63.2% | 0.2% | **326×** |

Enrichment *rises* with the threshold, which is the correct direction: spurious repeat-driven links are
low-count.

Per family at ≥3 shared reads: **NPIP 171/171 (100%)**, NBPF 91/91 (100%), NOTCH2NL 6/6, SULT1A 6/6,
PCDHB 123/153 (80%), TUBA 67%, GOLGA 59%, SPANX 55%, GSTM 50%, TBC1D3 42%.

⭐ **NPIP's 171/171 replicates across two independent samples** (A119b chr16 and testis) — the
certificate is complete for the flagship family in both.

⭐⭐ **Every between-family link is a known co-duplication.** The 1.1% resolves into exactly two family
pairs and nothing else: **NBPF–NOTCH2NL** (573 reads over 23 gene pairs; NOTCH2NL and NBPF are
co-duplicated on chr1, and NOTCH2NLB–NBPF14 fusion transcripts were observed here directly) and
**NPIP–SULT1A** (289 reads over 60 gene pairs; SLX1B-SULT1A4 is one of the two loci that carried all 31
RNA/DNA partition violations in the 08-06 audit). There are **zero** between-family links outside these
two pairs. The certificate's off-diagonal signal is biology, not noise.

⚠ Within-family coverage is a **lower bound set by expression**: a pair cannot be linked if one member is
silent in the sample. TBC1D3's 42% and GSTM's 50% are therefore not necessarily method failures. Testis
was chosen for provenance, not for expression breadth.

### Certificate claims C1–C3 (properties of the EVIDENCE, not of the definition)

Numbered separately from P1–P4 so they can never be read as part of `F(s)`.

- **C1 — completeness.** Within a family, what fraction of member pairs is linked by a shared read.
  NPIP: **171/171 (100%)**, replicated in two independent samples. Bounded above by expression: a silent
  member cannot be linked, so C1 is a lower bound, never a failure claim on its own.
- **C2 — specificity.** Within-family linkage must exceed between-family linkage by a wide margin, and
  the margin must *grow* with the read threshold (spurious repeat links are low-count). Measured 67× /
  179× / 326× at ≥1 / ≥3 / ≥10 shared reads.
- **C3 — abstention semantics.** A link records that the aligner **could not separate** two loci. It is
  an abstention, not a merge instruction. This is what makes the certificate safe to publish alongside a
  family without contaminating it.

## 3a. THE THREE USES OF MULTIMAPPING — keep them apart

Multimapping enters the method in three places. They are easy to conflate and must not be, because two
of them sit on opposite sides of the O1 ⟂ O2 line.

| # | use | object | side | role |
|---|---|---|---|---|
| **N0** | **presence**, `sec_frac(locus)` | per-locus **scalar** | **O1** | proposes candidate NODES |
| **C** | **conflict**, shared reads between loci | pairwise **relation** on a built family | held out | certifies the family |
| **K** | **copy number**, χ(H) of the conflict graph | pairwise relation, then a colouring | **O2** | lower-bounds the number of copies |

The distinction that makes N0 legal is that it is a **scalar, not a relation**. `sec_frac` says *reads
placed elsewhere also fit here*; it never says **where**, so it carries no information about which loci
belong together. K and C both do, which is why they must come after the family exists, never before.

> **The one-line rule: a scalar over loci may build the graph; a relation between loci may only test it.**

⚠ K is the established copy-number leg — `max(n_loci, χ(H)) ≤ true copy number`, with MCC = χ(H). It is
NOT what N0 measures and the two must never be reported as one quantity. N0 says *this locus is
duplicated*; K says *this family needs at least k copies to explain its reads*.

## 3b. N0 — Multimapper PRESENCE as a node signal (OFFICIAL STAGE)

**Status: promoted from observation to a named stage of the pipeline, 2026-08-07.** It is the only
annotation-free node signal measured to work, and node construction is the binding constraint on the
whole definition. ⚠ **The "purity 0.136 / node purity 0.237" figures that used to close this sentence
were WITHDRAWN on 2026-08-11** (no surviving artifact — see § 1 and `OBJECTIVES_AND_VERIFICATION.md`
row 1.13). The binding-constraint claim does not depend on them.

---

### ⭐⭐⭐ 3b-0. THE BLIND / SEED-FREE ROUTE, MEASURED END TO END (2026-08-12) — read this before the N0 table below

**HUMAN, CHM13 v2.0, `A119b.t2t.bam` (full BAM, no `-M -L`), at the SHIPPED tier**
(`minimap2 -c -X --no-long-join -k11 -w5`, identity `1 − de ≥ 0.60`, M1 coverage on the shorter side
`≥ 0.50`, then γ-quasi-clique at γ = 0.20). Verdict, full numbers table and adversarial ledger:
`/home/juanfra/winloci_scratch/o1_blind/O1_BLIND_VERDICT.md`. **Status: OPEN, with a measured ceiling
and a named mechanism.**

> **THE HONEST SENTENCE. Without a seed the method still FINDS the loci but cannot DELINEATE them, and
> the failure goes through exactly one clause of `E_r`.**

**Finding is not the problem.** Loci **ABSENT from the blind node set = 0/27 [0.000, 0.125]** at panel
scope (19 windows, 3.04 Mb, chr16 + chr18, **27 loci of ONE gene family — ⚠ there is no TBC1D3 among
them; the "NPIP + TBC1D3" scope line used elsewhere is withdrawn**). The blind node build then
completed **genome-wide**: 20,706,985 primary reads → **38,130 nodes → 1,016,370,811 exonic bp spanning
2.163 Gbp**, and on chr1 **gene_capture = 120/162 = 0.7407**.

**Delineation is.** Blind exonic bp ÷ seeded exonic bp per matched locus: **median 5.97×, range
2.45–16.62×**; 14/27 loci FUSED; 12 blind nodes swallow more than one locus.

**The mechanism, and it is the same one this document already asserts.** Of the 260 seeded `E_r` edges
(251 testable) the blind node set survives **6/251 = 0.0239**, and the 245 losses decompose as
**no alignment record 0 · identity < 0.60 on every record 0 · identity passes and coverage < 0.50 =
245/245 = 1.0000 [0.9846, 1.0000]** (best identity median 0.9889; best identity-passing coverage median
0.1945). ⟹ **identity never fails, coverage always does, and coverage is set by node length.** This is
now the **third** substrate carrying that signature, alongside the **91/91** DNA-only edges on the
seeded 27-node panel (identity 0.9709–1.0000, § P4b) and the **129/129** tier-lost edges on the
byte-identical 6-chromosome `sec_frac` node sets recorded later in this section.
⚠ **`TRIM` = 250/251 is ORACLE-INFORMED — it clips the blind node to the seeded interval, so it is a
CEILING, not a method; `FILL` = 251/251 is a TAUTOLOGY CHECK (FILL ≡ SEEDED on 27/27).** Their value is
the decomposition: **TRIM restores 244 of the 245**, so the blind node's exon content *inside* the right
interval is fine and the whole damage is the material it drags in from outside.

**The one end-to-end score, with its degenerate control.** ATTEMPTED SCOPE = genome-wide
(pre-registered); **ACHIEVED SCORING SCOPE = chr1 only.** ⚠⚠ **The reason is a property of the rule in
this document, not of the machine: at 488 Mbp of node sequence `-k 11` yields 1,520,563 distinct
minimizers, average occurrence 109.7, `mid_occ` = 8456 — the shipped `E_r` seed is too short to index a
genome-scale node set, and because `mid_occ` is computed from the whole index the shipped tier is NOT
SCOPE-INVARIANT: the same two nodes can pass at one scope and be filtered at another, so a larger
scope's PAF may never be sub-set to score a smaller one.** On chr1, against a **40-family / 162-gene
coarse denominator built from the GFF before the run** and a frozen sha-verified metric:

| arm / partition | rec@τ=0.50 | matched NULL | **m ≥ 3** | τ=1.00 | purity | null | ARI | gene_capture |
|---|---|---|---|---|---|---|---|---|
| **SEEDED ER_QC ≡ ER_CC** (ceiling) | **22/40 = 0.5500** | 0.4464 | **5/22** | 0.1000 | 0.1645 | 0.0655 | **0.3231** | 0.8272 |
| SEEDED SINGLETONS (degenerate) | 18/40 = 0.4500 | 0.4500 | 0/22 | 0.0000 | n/a | n/a | 0.0000 | 0.8272 |
| **BLIND ER_QC (shipped)** | **18/40 = 0.4500** | 0.4178 | **1/22** | 0.0250 | 0.0733 | 0.0400 | **0.0467** | 0.7407 |
| BLIND ER_CC (no γ) | 12/40 = 0.3000 | **0.3408** ⚠ below its own null | 1/22 | 0.0000 | 0.0140 | 0.0057 | 0.0002 | 0.7407 |
| BLIND SINGLETONS (degenerate) | 17/40 = 0.4250 | 0.4250 | 0/22 | 0.0000 | n/a | n/a | 0.0000 | 0.7407 |

⭐ **The shipped blind partition beats "every node alone" by exactly ONE family** (fine truth 19/38 vs
18/38), the recovered sets are **nested**, and the one extra family is **FAM72**. ⭐ **The honest cell is
the stratum where a single node cannot reach τ = 0.50 by arithmetic: BLIND 1/22 · SEEDED 5/22 · both
degenerates 0/22.** Also **0/8 dispersed** and **0/8 size ≥ 5** for *every* partition. Purity excess
+0.0333 decomposes as contamination **0.0976** plus unlabelled **0.9188**.

⭐⭐ **γ = 0.20 BINDS in blind mode and is INERT in seeded mode — the sharpest statement of γ's role this
project has.** On the blind arm γ moves recovery 0.3000 → 0.4500 and contamination **0.7473 → 0.0976**
(1,660 → 2,887 components); on the seeded arm `ER_QC ≡ ER_CC` to four decimals, the **21st** seeded cell
in which γ is measured inert. At this scale the refinement is **not** a tidy-up — it is the only thing
between the blind edge set and a partition that scores below its own null.

⚠ **AND THE CEILING IS NOT HIGH.** The seeded arm still recovers **0/8 dispersed** and **0/8 size ≥ 5**,
leaves 81.6% of its component nodes unlabelled, and its 1,936 nodes carry **80.6 Mbp = 92%** of the
blind arm's 87.6 Mbp (median node span 24.7 kb) because unspliced IsoSeq reads fill the gene body ⟹
**the extent problem is a property of the read-to-block rule, not only of blind delineation.**

**What was tried and failed — volunteer all of it.**
* ⚠⚠ **READ-THROUGH IS REFUTED, NOT UNSUPPORTED.** Ported verbatim from the in-tree Rust, the filter
  flags 1,492 of 39,752 reads (3.8%) and removes **0 of 6,824 block-links**; the node set is
  **bit-identical** at both link thresholds. `recombinant_split` is **0/251 analytically** and is not
  even seed-free (it gates on `both_annotated`). ⭐ **Max-flow says why, and it bounds a class of rules:
  the read min-cut between two fused loci is min 1, MEDIAN 12, MAX 273.**
* The mis-chain gate is the only mover and **does not clear its null** (13/256, null median 11 [8,14],
  p = 0.1534). Dropping unspliced reads outright is **strictly harmful** — it lowers even the oracle
  TRIM ceiling, 250 → 237.
* ⚠⚠ **A depth split-and-trim rule reaching 178/251 IS NOT QUOTABLE.** Its constant `c = 0.20` is the
  **argmax of a sweep against the truth denominator**; the selectors that read only its own histogram
  give **150–151/251 with 38 components, largest 352 of 407**; and a **size-matched random-cut control**
  (same node count, same bp, same piece count, cut points uniform at random) scores **152/251**, so the
  depth statistic buys **26 of its 172 edges**. ⚠⚠ **Two of three zero-information random cutters clear
  the run's own matched null at p = 0.0000 ⟹ "clears its matched null" is not evidence of delineation;
  the null must be size-matched.** This is the eighth absolute threshold to fail in this project, and it
  failed in the specific way this document warns about: **it is not pickable without the answer.**

**What does work, and what to build next.** The one seed-free rule that survives is a **junction-chain
split** introducing **no new constant** (junction floor and link floor are both the shipped
`MIN_READS = 3`, block floor 20 bp): **119/251 = 0.4741** [0.4132, 0.5358], locus-clustered bootstrap
[0.3156, 0.6429], **null median 106 [98, 114], p = 0.0011**, while both of its own controls sit exactly
on theirs; it beats a 20-replicate size-matched **extent** null; it takes FUSED **14/27 → 0/27** and the
extent ratio 5.97× → 1.89×. ⚠ but its losses stay **132/132 coverage**, extent precision is only
**0.381**, and under the coverage-of-the-**LONGER** certificate it is **25/251** against a seeded 159.

⭐⭐ **An oracle ablation ladder prices the target and names the object to build** (⚠ every cell is
oracle-informed; none is a blind result): handed the **perfect split** — the count *and* the grouping of
every fused node — recovery is **29/251, gain 0.094, and it FAILS its null (p = 0.0795)**; handed only
the **boundary**, with no un-fusing at all, it is **170/251, gain 0.672**; and the boundary tolerance is
narrow — **±1 kb already costs 7–8% of the ceiling, ±2 kb costs 8–18%**, on loci whose exonic content is
8–24 kb. ⟹ **SPLITTING IS NOT THE LEVER; THE BOUNDARY IS, to about ±1 kb.** ⭐ And because the M1
coverage form measures the **shorter** side, **an `E_r` edge needs ONE well-delineated endpoint, not
two**: 4 exactly-delineated nodes out of 27 recover **79 of the 80** pairs that touch them (0.9875)
against wholly blind partners.

**Exposures, stated here rather than waited for.**
1. ⚠⚠ **The panel-scope truth IS the seeded pipeline's own output** (`gvg_graphs.json["rna"]`; the
   seeded arm scores **251/251 = 1.0000** against it — identity, not measurement; the same artifact's
   `["dna"]` is the **complete graph** on those 27 nodes, so the 91 "non-edges" are RNA-method misses,
   not negatives) ⟹ **every 19-window survival number is concordance, never accuracy.** Only the chr1
   arm is scored against annotation the method never saw.
2. ⚠⚠ **The panel's scope windows truncate the read set they are blind inside.** 298,847 / 2,035,590 =
   **14.68%** of blind exonic bp lies outside its own fetch window; refetching the same territory takes
   exonic bp **2.05×**, with the gradient tracking the 50 kb pad exactly (≤5 kb 1.11× · 5–50 kb 1.22× ·
   **> 50 kb 11.20×**), independently confirmed windowless at 1.71×. ⟹ panel "blind" means blind
   **delineation** inside a seed-derived scope. The chr1 arm is windowless and does not inherit this.
3. ⚠ **The blind arm's published null was in the wrong units** (observed = locus-pair survivals, null =
   node-pair hits). Recomputed in locus-pair units: **null median 10 [5, 16], p(≥ 6) = 0.9644** ⟹ blind
   edge survival is at the **3.6th percentile — significantly BELOW chance, not "at chance"**. The
   correction strengthens the negative, and the corrected form is what to say.
4. ⚠ **The 251 denominator is conditioned on the blind arm's own output** (260 − 9 fused-internal, the
   9 pairs it failed worst); the same convention can be exploited for a **3.0× inflation by deliberately
   destroying delineation** (6/85 = 0.0706). Arms reported on 256 or 257 are cross-incomparable.
5. ⚠ **Panel survival has no precision term** (Pearson r = 0.993 with the raw edge count among matched
   nodes), and the scoring-time `match` argmax is an oracle — now **quantified at ~12 edges (~5 pp)**,
   so it is recorded but is not doing the work.
6. ⚠ **chr1's denominator is 80.0% clustered against a pre-registered genome-wide 36.3%** ⟹ 0.4500 is
   biased **upward** and is not an estimate of the genome-wide number. And **18 of the 40 families have
   exactly two members**, which is why a degenerate partition reaches 0.4250 ⟹ **never quote τ = 0.50
   without the all-singletons number beside it; quote the m ≥ 3 cell.**
7. ⭐ **The standing worry that `min()`-coverage forgives fragmentation is REFUTED** — a shatter arm
   (every exon block its own node, near-perfect extent 0.9618) scores **102/251 and does not clear its
   null (p = 0.059)** — **but survival is therefore non-monotone in extent** (0.167→6 · 0.411→100 ·
   0.522→152 · 0.734→178 · **0.962→102**), so *"extent is the whole story"* holds on the way up and fails
   past the optimum.
8. ⚠⚠ **"All 27 panel loci land in ONE blind component" is true and useless** — that component holds
   **39** nodes (61 · 69 · 119–135 · 202 under the other arms). **Quoting 27/39 as a purity is exactly
   the oracle-selection shape that killed the withdrawn 0.237. Do not.**
9. ⚠ **HUMAN ONLY.** Gorilla cannot adjudicate this — `GGO_genomic.gff` names 1 of 19 NPIP genes and 0
   TBC1D3, and 46.6% of its named genes are `LOC` ids.
10. ⚠ **The withdrawn 13.6% / 0.237 are NOT restored by this run.** Nothing above re-derives them; they
    stay withdrawn permanently. **No Rust was changed and nothing is committed.**

**What would close it, costed.** ⭐ Cheapest and already built: **run the declared K-ladder** —
chr1..chr2 (177.0 Mbp, 6,209 nodes) and chr1..chr3 (250.6 Mbp, 8,556 nodes) node sets, FASTAs and
indexes exist and only the all-vs-all is missing (**~1.2 h + ~2.5 h**, one stream at a time); it moves
the clustered fraction **0.800 → 0.532 → 0.407** toward the pre-registered 0.363 and attacks the largest
caveat on the headline. Then **re-score the panel arms against the GFF** instead of against the seeded
pipeline's own output (**~2 h, no new alignment**). Then, high-risk: **a seed-free BOUNDARY rule** — not
a splitter — with a constant pickable from its own histogram or already shipped, scored against a
size-matched null from the start. **What will NOT close it, measured:** node *filters* (with
ABSENT = 0/27 the route already finds every locus, so a filter picks *which* loci become nodes and
cannot touch extent — this bounds N0 below); read-dropping filters of any kind (median-12 / max-273
min-cut); threshold sweeps on `c`; more panel windows of the same family.

---

> **`sec_frac(locus)` = secondary records / (primary + secondary) records overlapping the locus.**

**Stage contract.**
* INPUT: an aligned read set. No annotation, no seed, no catalog.
* OUTPUT: candidate intervals, ranked. **Candidates, not loci** — `E_r` is the filter that follows.
* GUARANTEE: nothing pairwise is emitted, so `E_r` and the certificate stay independent of it.

### N0 measured, seed-free and annotation-free (5 chromosomes, fixed parameters)

Candidates from `sec_frac` alone → `E_r` → components. The curated families are read ONLY to score.
The top-1000 cut, 5 kb bins and ≥10 kb interval floor were fixed on chr16 and applied **unchanged**
elsewhere — retuning per chromosome would prove nothing.

⚠⚠ **TIER + SPECIES, added 2026-08-11 — this whole table is HUMAN (CHM13 v2.0) and it is at NEITHER
of the two tiers named in the tier notice. It is OPEN.** `seedfree_multichrom.py` (and its two
siblings `seedfree_adaptive.py`, `seedfree_secfrac.py`) run
**`minimap2 -x asm20 -k11 -w5 -c --eqx -N 200 -p 0.02`** and threshold on
**`nmatch/blocklen ≥ 0.60`** — an asm20 preset carrying the sensitive tier's `k`/`w`, no `-X`, and the
identity metric the binary uses only as a *fallback*. The shipped tier is
`minimap2 -c -X --no-long-join -k11 -w5` with `1−de ≥ 0.60`. **Every `recovered`, `purity` and
`density` cell below therefore names a tier the binary does not run**, and `density` is the cell most
exposed (it is a pure `E_r` quantity). This was **not** recomputed in the 2026-08-11 shipped-tier pass:
it needs the 5-chromosome seed-free candidate build re-run end to end, which is on the expensive list.
**Quote these rows only with "third tier, human, OPEN" attached.** The *qualitative* claims —
`sec_frac` discovers nodes, γ catches roughly half of the over-merges, families come out as components
rather than a hairball — do not turn on the tier and are what to say.

⭐ **THE SCRIPTS ARE FIXED (2026-08-11, later the same day); THE TABLE IS STILL OPEN.** All three now
import `bench/soto/rustlib.py` for the aligner (`ava()`), the identity metric (`1 − de`) and the M1
coverage form, so there is one definition of the rule and no fourth tier. Their PAFs are written as
`*.er.paf` — a **different filename**, because `rustlib`'s disk cache can only see the tier through
`(path, mtime_ns, size)` and `st_mtime_ns` granularity on the scratch filesystem measured ~3.6 ms
(two writes 2 ms apart returned the *identical* `mtime_ns`); reusing `*.paf` would have let
`if not os.path.exists(paf)` serve panel alignments to the shipped rule. **The table below was NOT
regenerated** — its cells are per-family `max(comps, key=hits)`, i.e. **oracle-selected**, the shape
that killed the 0.237 blind-DNA purity, so re-running it at the right tier would only produce a
right-tier version of a metric that cannot be quoted. It is superseded, not refreshed. The scripts now
also write a whole-partition TSV (`*.partition.tsv`, every component including singletons) as the
substrate for a pre-registered whole-partition scorer.

What the tier change does, measured on the **byte-identical** `sec_frac` node sets (6 chromosomes,
876 nodes), decomposed so the aligner and the metric are not confounded:

| | edges | Δ from the metric alone | Δ from the aligner alone |
|---|---|---|---|
| old rule, panel PAF | 2,736 | — | — |
| shipped rule, panel PAF | 2,750 | **+35 / −21** | — |
| shipped rule, shipped PAF | 3,067 | — | **+446 / −129** |

⟹ **the aligner moves ~13× more edges than the metric**, and the partition changes on **6/6**
chromosomes (ARI 0.659–0.966 vs the old partition). γ = 0.20 is **inert on both sides**: 0 emitted
components of size ≥ 3 fall below density 0.20 at either tier, so "components" *is* the shipped
partition here — but the "γ catches roughly half of the over-merges" line above is a **panel-tier**
observation and is not reproduced at the shipped tier.

⭐ **All 129 tier-lost edges are lost on COVERAGE. Zero on identity, zero to the pair vanishing from
the PAF** (median best-record coverage of a lost pair 0.298–0.492, all below the 0.50 floor) — the
mechanism prediction "identity never fails, coverage always does", reproduced on a third substrate.
⚠ **COST:** the shipped tier is far more permissive per record (chr11: 19,656 records vs 3,856;
minimum block length 55 bp vs 290; median 304 bp vs 2,908) and on segdup-dense chr16 the all-vs-all
took **515 s for a 62.6 MB PAF against the panel's 2.8 MB** — budget for it.

| chrom | family | members | recovered | purity | density |
|---|---|---|---|---|---|
| chr7 | ID_207 SPDYE | 17 | 15/17 | **0.833** | 0.660 |
| chr15 | ID_116 GOLGA6L | 6 | 6/6 | **0.857** | 1.000 |
| chr17 | ID_14 LRRC37 | 7 | 7/7 | 0.714 | 0.297 |
| chr17 | ID_468 TBC1D3 | 9 | 7/9 | 0.636 | 0.600 |
| chr16 | ID_154 NPIP | 19 (RefSeq) | 18/19 | 0.556 | 0.530 |
| chr7 | ID_8 PMS2P | 9 | 7/9 | 0.389 | 0.660 |
| chr15 | ID_113 GOLGA8 | 16 | 16/16 | 0.279 | 0.210 |

⭐ **Purity 0.28–0.86 against blind self-alignment's 0.136** — a 2–6× improvement on the ceiling that had
not moved. Families come out as components, not as a hairball (chr16: 157 candidates → 78 components).

### N0's two failure modes — one fixable, one not

**(1) FIXED CANDIDATE BUDGET — engineering, fixable.** The top-N cut is a constant per chromosome
regardless of how much duplicated sequence that chromosome carries. chr2 got 97 candidates for 31
curated members and returned ID_63 1/8, ID_65 **0/7**. ⚠ This is NOT a signal failure: those families
have `sec_frac` 0.84–0.93 and 64–1047 reads per locus. The signal is present and the budget excludes
them. **This is the seventh absolute threshold to fail in this project — scale N to the chromosome's
duplicated content.**

**(2) CO-DUPLICATED FAMILIES — the real limit, NOT fixable by a threshold.** chr7 merges SPDYE and
PMS2P into ONE component at density **0.660**, comfortably above γ = 0.40, so the quasi-clique
requirement does not catch it. ⭐ Note γ *does* catch the other case: chr15's 61-locus blob holding
GOLGA8 + UBE2Q2P + fragments sits at density **0.210**, below γ, and is correctly rejected. So γ
discriminates over-merges roughly half the time, and cannot be relied on.

> ⟹ **N0 solves node DISCOVERY, not family SEPARATION.** The seed's remaining job is not finding
> copies — it is saying WHICH family they belong to when two families share a duplicated block.

### What N0 does not do

* it does not replace the seed (purity 0.56 seed-free vs 1.000 seeded on NPIP)
* it does not measure copy number — that is K, `χ(H)`, and it is O2
* it does not define edges, and must never be allowed to

## 3c. Multimapper PRESENCE — supporting measurements

A third use of multimapping, and the only one that touches O1 without breaking it:

> **`sec_frac(locus)` = secondary records / (primary + secondary) records overlapping the locus.**

A secondary alignment here means "a read whose best placement is elsewhere also fits HERE". It is a
**per-locus scalar**, never a pairwise relation — it does not name the other locus — so no `E_c`
information enters `E_r`. It answers *is this locus duplicated at all*, which is a NODE question.

Measured on chr16 (A119b, `chr16_sub.bam`), 81 loci:

| | n | median `sec_frac` | range |
|---|---|---|---|
| duplicated | 24 | **0.94** | 0.016–0.989 |
| single-copy | 57 | **0.05** | 0.000–0.960 |

AUC sits on a **plateau of 0.94–0.98 across 11 of 15 DNA-labelling settings**, peaking at 0.9786 when
"duplicated" = a paralogous block ≥3 kb at ≥95% identity.

⚠ **The plateau is the claim; the 0.9786 peak was selected post hoc.** AUC falls to 0.887 at the loosest
label (≥1 kb / ≥90%) because that label calls ordinary genes duplicated on an incidental fragment, and
to 0.900 at ≥20 kb because that label excludes real family members. `sec_frac` tracks *multi-copy gene*,
not *any paralogous fragment exists* — the two predicates are not the same, and the AUC dip measures the
label, not the signal.

⭐ Consistency: **SULT1A1** has the lowest positive `sec_frac` (0.417) *and* no DNA paralog at ≥1 kb/90%.
Signal and DNA independently agree it is the odd member; prefix-based labelling had wrongly called it a
positive.

**Scope, stated narrowly:** this is *evidence for node selection*, not a node builder. It scores a
candidate locus; it does not draw a boundary. No cut is proposed — the sweep is precisely the reason
not to propose one.

⚠ **Multimapper WIDTH does not work, and the reason is structural.** Using the footprint of the
multimap track as locus extent gives CV 0.482 against annotation's 0.412, and adjacent copies return
*identical* footprints (NPIPA6 ≡ NPIPA7, NPIPB6 ≡ NPIPB7, NPIPB8 ≡ NPIPB9). The track marks the
duplicated REGION, which contains several copies; it cannot resolve a per-copy boundary. Do not retry.

⚠ **First attempt used MAPQ=0 as the multimapping proxy and was wrong.** At NPIPB11, 234 of 271
*primaries* carry MAPQ 60 while the window holds 25,287 *secondaries* — minimap2 breaks ties
confidently, so the evidence is in the secondary records, which `-F 2308` removes. Also: minimap2
emits **no NH tag**; `sec_frac` is the NH-free equivalent.

## 4. Known exposures (state these first, do not wait to be asked)

- **PDXDC1 is a genuine false positive** in the human panel — a real co-duplicated passenger, not a
  random error. Human precision is **19/27 = 0.70 to 25/27 = 0.93** depending on how many unnamed `LOC`
  copies are counted as real; the spread is itself the annotation-is-a-hypothesis point.
- **Gorilla returns 25 loci against a published morpheus copy number of ~17** (human 15, chimp 25–30,
  orangutan 9; 1–2 copies is Old World monkeys). Unresolved.
- **Node extent, not homology, is the recurring failure.** Identity never fails (171/171 pairs clear it
  on both substrates); coverage is the only failure mode; gorilla's largest "locus" is a 207 kb merge of
  a tandem cluster.
- **The two-definition defect survives in ONE place: `copy_assign --homology-primary`.** The catalog
  binary is fixed (§1), but `detect_and_assign` still runs refine unconditionally. That is correct for
  copy_assign's default (`homology_primary` is off ⇒ the families are E_c conflict-derived, which is the
  catalog refine was written for); only the `--homology-primary` combination reproduces the D1 shape. It
  was **not measured**, and was therefore left alone rather than changed blind.
- **The coverage clause `c` is scale-dependent and the edge substrate is a free parameter** — §1a. Its
  absolute demand differs 5.5×–9.8× between the two shipped paths (6,338 bp on genomic span vs 876 bp on
  a rep transcript), and 5 of 7 families change partition when the substrate changes. Never quote an
  `E_r` number without its substrate. ⚠ **This exposure is UNCHANGED by the M1 fix** — M1 made `c` a
  fraction; it did nothing about the fact that a fraction of *what* is a free parameter.
  ~~And its denominator is unbounded above~~ — **CLOSED 2026-08-10 (M1).** `c` was
  `(qe−qs)/min(|u|,|v|)`, a query-axis numerator over a possibly-target denominator, which exceeded 1.0
  on **52 of 68** accepted edges (max 2.019) in the gorilla SDHA probe and on 21/255 (max 1.142) on the
  25-region control panel. It is now `qlen ≤ tlen ? (qe−qs)/qlen : (te−ts)/tlen` — bounded by 1.0 and
  symmetric under q/t exchange. Cost: **1 of 25 control regions changes** (HERC2 gains a true annotated
  HERC2 paralog at identity 0.9940), controls stay 0/7, 0 edges lost panel-wide.
- **`τ` and `γ` cannot bind, and `β`/`d` are node constants, not edge constants** — §1. The edge rule
  has **one** operative constant. This is a stronger claim than the old five-constant §1 and is stated
  there with the census, counterfactual and sweep that support each demotion — including the one that
  went the other way (`β` is operative and the "β does not change the partition" claim is refuted).
- **`T` WAS two different minimap2 invocations — CLOSED 2026-08-10 (B1), but the numbers it produced
  are not all re-run.** The shipped `-c -X --no-long-join -k11 -w5` and the panels'
  `-k11 -w5 -c --eqx -N200 -p0.02` differ on byte-identical FASTA (HERC2 DNA density 0.8056 vs 0.4722;
  partition differs on 4/14 panels, edge count on 10/14). The tier is now a single constant in the Rust
  and a single constant in `bench/soto/rustlib.py`, guarded by a test, and all eight panels issue it.
  ⭐ **BACKLOG CLEARED 2026-08-11.** The P2 NPIP density rows and the P4b 27-node table have been
  recomputed at the shipped tier on both genomes (one `seed_family.sh` invocation per species, plus
  `graph_vs_graph.sh`). **P2: human 351 / 1.000 unchanged (Jaccard 1.0000); gorilla 261 → 260 edges,
  0.870 → 0.867, one component either way.** **P4b: RNA 338 → 260 edges, Jaccard 0.963 → 0.741, still
  0 RNA-only edges** *(on this 27-node panel — ⚠ 33/481 on the 80-node panel, §1★.3)*. P1, P3 and P4 are tier-invariant by construction and reproduce to the digit.
  ⚠ **What remains open is §3b's N0 table**, which is at a THIRD tier again
  (`-x asm20 -k11 -w5 -c --eqx -N200 -p0.02`, `nmatch/blocklen ≥ 0.60`) and needs the 5-chromosome
  seed-free candidate build re-run end to end. ⚠ The HERC2 figures quoted in this bullet
  (0.8056 vs 0.4722) are from the 14-panel B1 sweep and are **not** the same node set as the HERC2 rows
  in §2; do not cross-quote them.
- **The over-merge mode is still UNTESTED (D4).** The single-copy control panel could only ask whether
  the method invents homology between *unrelated* genes; the identity floor was never exercised, because
  no control pair had coverage ≥ 0.5 at identity in (0.60, 0.999). "No over-merge anywhere" is void as
  stated — it was suppressed upstream by how nodes were pre-selected. The next panel must be
  **paralog pairs**.

## 4a. ⭐⭐⭐ REACH — THE BOUND, AND WHAT SETS IT (2026-08-14)

Verdict `/home/juanfra/winloci_scratch/o1_dispersed/O1_DISPERSED_VERDICT.md`; prereg
`PREREG_DISPERSED.md` sha256 `e8691ca1…0847b0`; frozen metric `blind_metric.py` sha256
`7da48616…4c18b710`, verified and unforked. Scope of every number below unless stated:
**HUMAN CHM13 v2.0, chr1, RefSeq full GFF, shipped tier.** Denominators fixed before any method ran
(`DENOMINATORS.txt` 22:25 < prereg 22:29 < first alignment 22:37): **W = 73** within-family gene pairs
over the 8 dispersed coarse families · **WC = 460** clustered within-family pairs · **X = 12,508**
cross-family pairs · 73 + 460 + 12,508 = 13,041 = C(162,2), exhaustive.

**THE BOUND.** O1 as shipped recovers **22 of 40 = 0.5500 [0.3983, 0.6929]** chr1 multi-copy families
(frozen `bm.recovery` at τ = 0.50, coarse, recall AND precision), and that panel is **statistically
indistinguishable from the genome by family composition** — the same 40 chr1 families scored on their
FULL genome-wide membership are **16/40 = 40.0% [0.2635, 0.5540]** clustered against genome-wide
**157/432 = 36.3%**, Fisher **p = 0.6090**. So ≈ 0.55 is a defensible genome-wide estimate, **not an
inflated one**.

**WHAT SETS IT IS DIVERGENCE, NOT DISPERSAL.** Families with median within-family **protein** identity
≥ 0.50 are edge-bearing **9/14 = 0.6429 [0.3876, 0.8366]**; below it, **1/26 = 0.0385 [0.0068, 0.1889]**.
Binning within-family pairs by protein identity (DNA edge set): < 0.40 → dispersed **0/66** vs clustered
**1/37** (rate 0.027), **Fisher p = 0.3714 — the coordinate stratum has NO measurable effect at matched
divergence**; 0.40–0.60 → 0/7 vs 5/176; 0.60–0.80 → 0/0 vs 24/152 (0.158); ≥ 0.80 → 0/0 vs 38/91 (0.418).
Genomic distance does not explain it either: clustered pairs > 5 Mb apart fire at **16/83 = 0.193**,
*higher* than clustered pairs < 100 kb apart (**15/120 = 0.125**).

**THE FAILURE, MEASURED.** On the 8 dispersed families the rule fires **0/73 [0.0000, 0.0500]** on the
RNA edge set (cross-family 10/12,508) and **0/73** on the DNA edge set (cross-family 17/12,508); on the
ideal exon-union substrate at the shipped tier the aligner emits a record for only **2/73** and an edge
for **0/73** (survivors SLC25A33/SLC25A34 id 0.8368 cov 0.0755, SLC35E2B/SLC35A3 id 0.8525 cov 0.0434 —
both die on COVERAGE). Genome-wide over all 235 members of these families: **16 edges of 4,287**
within-family pairs, with identity ≥ 0.60 on **194/194** records that exist.

**IT IS NOT AN ALIGNER-TUNING FAILURE — THE LADDER SETTLED THIS BEFORE ANY ROUTE RAN.** L1 seed census
(chr1 exon-union, the shipped `-k 11 -w 5` seed shape, canonical minimizers, low-complexity masked,
occurrence-filtered at `mid_occ` = 21): **3 of 73 = 0.0411 [0.0141, 0.1140]** pairs exceed the 99th
percentile of BOTH nulls (N1 = 100 dinucleotide-preserving shuffles per pair; N2 = the 12,508
cross-family pairs), against a clustered positive control of **207/460 = 0.4500 [0.4051, 0.4957]**.
Positives' usable-shared median is **7**; the *non-homologous* cross-family median is **3**.
**minimap2 cannot align what it cannot seed.** Confirmed downstream: over 20 sensitivity settings,
**57/73 = 78.1%** of pairs never received a single record at ANY setting, and the union ceiling is
14/73 edges / 16/73 records.

**AND FOR THE PAIRS THAT DO ALIGN, THE LIMIT IS THE SHAPE OF THE RULE, NOT THE ENGINE.** The 28
shipped-tier within-family records (8 distinct pairs) are real contiguous alignments — identity
0.746–0.891, matched/span 0.749–0.904, **zero** ≥ 100 bp gap-leaping — spanning **126–312 bp** of genes
**3,249–11,271 bp** long. Max M1 coverage across all 28 is **0.0856**: **6–15× below `c` = 0.50**. There
is only ~250 bp of alignable material per pair, so no chaining or sensitivity setting could ever
assemble a ≥ 0.50-coverage alignment. A **local**-homology rule would fire where a ≥ 0.50-coverage rule
provably cannot — that is a definitional question for O1, not an engine question.

**THE HOMOLOGY IS THERE — the sequences are ancient paralogues, not spurious.** Protein Smith-Waterman
(parasail 1.3.4, BLOSUM62, open 11 / extend 1) over the same 73 pairs: **70/73 = 0.9589
[0.8860, 0.9859]** at E < 1e-3 against a **7.01%** cross-family background; 62/73 at E < 1e-6; 52/73 at
E < 1e-10. **0 of 73** fall below the median cross-family pair by SW score. Median within-pair protein
identity is **0.2809**. ⟹ **A rule whose first clause is a ≥ 0.60 NUCLEOTIDE identity floor cannot, even
in principle, reach objects defined at ~0.28 protein identity.**

**THE SENTENCE FOR THE THESIS.** *O1 reaches recently-duplicated multi-copy families — roughly 55% of
annotated multi-copy families on the only panel where it is measured, and that panel is representative
of the genome by family composition (p = 0.61). It does not reach anciently-duplicated families: below a
median within-family protein identity of 0.50 it recovers 1 of 26, and on the 8 chr1 families that are
dispersed by coordinates it recovers 0 of 73 pairs [0.0000, 0.0500]. That boundary is set by sequence
divergence, not by copy dispersal, and it is not moveable by any nucleotide aligner.*

### 4a.1 What is WITHDRAWN or RESTATED by this pass (retractions travel together)

1. ⚠⚠ **WITHDRAWN — "chr1 is 80.0% clustered against a genome-wide 36.3%, so every headline O1 number
   lives on the stratum that flatters the method."** The two figures come from **two different rules**.
   `strata()` = "≥ 80% of members on one chromosome AND ≥ 50% inside one 5 Mb window", but
   `build_truth(scope=('chr1',))` puts **100% of members on chr1 by construction**, so the chromosome
   clause is **unfalsifiable in chr1 scope**: chr1's 80.0% is a **window-only** number and 36.3% is a
   **two-clause** number. Genome-wide, **271 of 432** dispersed families are dispersed because of the
   chromosome clause and only **4** on the window clause alone. Matched, chr1 is 40.0% vs 36.3%,
   **p = 0.6090**. ⟹ also delete the "0.4500 is biased upward" caveat attached to the blind chr1 score.
2. ⚠ **RESTATE "0 of 8 dispersed families" as "1 of 24 genome-wide-dispersed chr1 families"**
   = 0.0417 [0.0074, 0.2024]. The chr1-cut contrast (clustered 10/32 vs dispersed 0/8) is **p = 0.0761,
   NOT significant**; 0/8 has Wilson [0.0000, 0.3244] and is compatible with a true 32% rate. The
   genome-wide-membership stratification is significant: clustered **9/16** vs dispersed **1/24**,
   **p = 0.000333**. **This correction makes O1 look WORSE, and it must be published with (1) — making
   only the flattering-panel correction would be a thumb on the scale.**
3. ⚠ **RESTATE "dispersed" as "anciently diverged".** The coordinate stratum does not survive control
   for divergence (p = 0.3714). "8" is not a robust object either: window sweep 0.1 Mb → 14 dispersed ·
   1 Mb → 10 · **5 Mb (shipped) → 8** · 10 Mb → 3 · 50 Mb → 0, while the genome-wide figure moves 6 pp
   across the same 500× sweep.
4. ⚠ **"The declared identity floor is 0.60, the EFFECTIVE floor is ~0.75–0.80" does not describe this
   stratum.** An effective floor implies alignments proposed at 0.62 and rejected; there are none.
   **71 of 73 pairs produce NO RECORD AT ALL.** There is no floor; there is silence.
5. ⚠ **The expectation "dispersed families are the FIRST place identity binds" is REFUTED — identity has
   still never been the failure mode, on a SEVENTH substrate.** Zero records below identity 0.60 across
   all 20 settings (min observed **0.6437**); 194/194 genome-wide. On protein it is worse than useless:
   shuffle-null identity p99 **0.6364**, max **0.8571**, and median cross-family protein identity
   **0.3111 EXCEEDS** median within-dispersed **0.2809** — protein identity is **anti-informative** here,
   and the shipped rule transposed verbatim to protein (id ≥ 0.60, cov ≥ 0.50) gives **0/73**.
6. ⚠ **The truth's COARSE key over-merges THIS stratum specifically**: coarse families that split at the
   frozen FINE level are **dispersed 4/8 vs clustered 1/32** (SLC35 = four RefSeq subfamilies —
   SLC35**E**2B / **D**1 / **A**3 / **F**3 — in one key; kinesin 6 members / 5 fine keys; "tnf receptor
   **superfamily**" says it in the name). Any claim about these 8 must replicate at the fine level.
   The fine dispersed denominator is **27**, not the prereg's 25 (arithmetic slip; families unchanged).
7. ⚠ **`div_strata`'s per-family medians are GENOME-WIDE.** "kinesin family 0.7418" and "solute carrier
   family 35 0.8525" invite the reading that the chr1 members are 74–85% identical. Restricted to chr1
   members, **every** kinesin and **every** kelch pair has **no record at all**. Never pair a `div_*`
   label with a chr1 result without naming the label's scope.

### 4a.2 The three routes that were tried, and why the definition does not change

Recall is never quoted without the specificity from the **same setting**. Full tables in the verdict.

* **A more sensitive NUCLEOTIDE engine — DEAD, and it is the important negative.** 20 settings.
  Best gene-arm `S09_k7w1_max`: **14/73** [0.1178, 0.2966] **with 1,172/12,508** cross-family
  (ΔX = 1,142), lift 2.101, p = 6.045e-3, γ-QC largest block **127 of 162 genes**, frozen recovery
  **0.6500 → 0.2500**. Node-arm `N_splice_k9w3`: **48/73** [0.5433, 0.7561] **with 5,092/12,508 = 40.7%
  of ALL cross-family pairs**, ONE block of **121 genes spanning 39 of 40 truth families**, recovery
  **0.5500 → 0.1000**. The set of settings with W > 0 **and** recovery ≥ the 0.6500 baseline is **EMPTY**.
  ⭐ **Mechanism: every gain is a gap-spanning artefact.** Matched-residues/span median **0.0496** (S09)
  and **0.0363** (N_splice), with **0.9360** / **0.9684** of span consumed by ≥ 100 bp D/N operations
  (224/224 and 1000/1000 records below 0.50), against shipped-tier clustered **0.8556** and cross
  **0.9115** at long-gap fraction **exactly 0.0000**. A repaired coverage M1′ (matched residues /
  len(shorter) ≥ 0.50) takes W **14 → 0** and **48 → 0**. Both pre-declared controls **FAIL**: P4a
  length-matched control 171/1,114 = 0.1535 vs 14/73 = 0.1918, p = 0.235; P4b size-distribution-matched
  null observed 14 vs p97.5 = 14–15 (null mean ≈ 7.4). And S09's p = 0.006 does **not** survive the
  matched instrument (73 within vs 593 between-dispersed, same 37 genes): 0.1918 vs 0.1535, ratio 1.250,
  **p = 0.2439**.
* **PROTEIN — recovers the most and FAILS the hairball guard; adopt it as a DIAGNOSTIC, not a
  definition.** E < 1e-10 & cov ≥ 0.50: **34/73 = 0.4658** [0.3559, 0.5790] over 6/8 families **with
  634/12,508** cross-family (624 new), lift 9.19, Fisher p = 2.08e-24, largest component **S = 39 > 20**
  ⟹ **P3 FAILS** (P1, P2, P4a, P4b all pass). γ = 0.20 does **zero** work — it splits **0 of 27** raw
  components and the 39-gene block has density **1.0000**. Attribution: **606/634 = 95.6% OR↔OR** plus
  **28** nbpf↔"neuroblastoma breakpoint family" (one family under two RefSeq wordings), and **0** of the
  5,218 cross pairs touching a dispersed gene. ⚠⚠ **LEASH: hand-deleting OR↔OR gives S = 20 and may NOT
  be used to re-score P3** — that is the recorded "selecting WHICH component to score" death shape.
  Costs that must be written, not appended: **coding-only** (the truth omits ≥ 4/166 = 2.4%, a FLOOR —
  chr1 carries 1,462 `pseudogene` FEATURE records `build_truth` never reads; external anchor: Soto's
  families include pseudogenes and only 455/1002 = 45% are brain-expressed) and **domain, not
  duplication** (**29/73 = 0.3973** [0.2929, 0.5119] of the recovered relations rest on sub-half-length
  blocks). Median block purity 1.0000 is a **median trap** — size-weighted purity is **0.8219**.
  ⚠ Exposed to the tautology objection: the truth families are RefSeq **description strings** and those
  strings are largely domain-derived, so the honest sentence would be *"protein recovers the
  annotation's own naming criterion"*, not *"protein recovers the family"*.
* **SHARED CODING-INTRON POSITION + PHASE — the only route clearing all four bars, but DEFERRED.**
  At (a) W = 100, s_max ≥ 0.90, SIC ≥ 3: **5/73 = 0.0685** [0.0296, 0.1505] over 2 families **with
  ΔX_new = 0/12,508**, largest component **S = 14** at density 0.2637 ≥ γ = 0.20 (kept whole),
  Fisher p = 5.72e-12, P4a and P4b both **PASS**; frozen recovery 0.5500 → **0.6250**, dispersed-only
  ablation **+0.0250**; taste 1 receptor recovered whole (3/3) ⟹ **1 of 8 families**, not 0.
  Zero-hairball ceiling: of the **113** settings with ΔX_new = 0 the maximum ΔW is **5**. Strand trap did
  not fire (1,312/1,312 phase agreements, 651 minus-strand, 0 disagreements). ⚠ **Deferred for three
  reasons:** (i) transcript-representative policy is a **chosen, unswept** knob — longest-CDS instead of
  RefSeq Select gives ΔW **2/73** at (a), which **fails P1**; (ii) the route reaches only paralogs that
  **retained every intron** (it fires on 5 of 10 *structurally eligible* positives); (iii) it is
  **annotation-dependent**, a step away from "genome + minimal annotation + reads", and cannot be carried
  to gorilla without the recorded circularity. Its (b) recovery gain of +0.1000 is **WITHDRAWN** as a
  gene→node projection artefact (all-overlapping-nodes projection gives +0.0500 **and loses** a family);
  quote the dispersed-only ablation, never the headline.

### 4a.3 Limits of §4a, before being asked

**n = 8 families / 73 pairs, and the family-level instrument is nearly powerless** (0/8 = [0.0000,
0.3244]); only the pair count is powered, and pairs are not what O1 delivers. **There is no orthogonal
truth** — the families are RefSeq description strings, largely domain-derived, and the pair route to a
certified negative set is structurally dead (n = 0), which is why specificity is a cross-family **edge
count**. These are **superfamilies, not duplication classes**: gene lengths span **138×** (4.0 kb to
558 kb) and member spacing 8 kb to 235 Mb. **Human chr1 + RefSeq only**; the gorilla leg is untouched and
a gorilla truth built from human annotation is circular. The protein route's worst cost is **invisible on
this denominator by construction** (all 162 truth genes are protein_coding). The divergence explanation
is **partly circular** — the non-circular form is the screening-off statement (p = 0.3714). The chr1 →
genome extrapolation of 0.5500 is licensed only by a family-**composition** test; whether chr1's
**divergence** composition is representative was **never measured — OPEN**. And the "no nucleotide
engine can do it" claim is bounded to **this seed shape, this chaining and 20 settings around it**: a
different engine class (spaced seeds, profile HMMs, translated search) was not built — the L2 k-mer
**ruler** rung was never run, so a calibrated nucleotide local-identity band for the 73 pairs does not
exist.

## 5. Three false-positive filters that were tried and failed

Recorded so they are not re-proposed.

1. **Expected-copy-number prior per species** — premise false (gorilla has ~17 NPIP copies, not 1) and it
   is an oracle: it requires the answer per family per species before the method can run.
2. **Seeding with the mRNA instead of the genomic span** — human recall 19/19 → **9/19**, component
   splits. It looks good on gorilla, the panel with no ground truth, and fails on the one that has it.
3. **Coding fraction as a per-locus certificate** — runs *backwards*: annotated NPIP loci have median
   mRNA coverage **0.080**, non-annotated loci **0.443**; NPIPB2/A1/B8/B15 score **0.000**. Consequently
   the inference "gorilla loci with no coding sequence are flank, not gene" is **withdrawn**.
