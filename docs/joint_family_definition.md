# The joint (DNA + RNA) multi-copy family: definition, what it forbids, and its falsifier

## ⚠⚠⚠ VERDICT — 2026-08-13, AFTER PHASE 3. READ THIS BEFORE ANYTHING ELSE.

> **PROPERTY, NOT DEFINITION.** The DNA/RNA relation is a **measured property of the shipped O1
> definition** and must be reported as one. It may **not** be made definitional.
>
> **The honest sentence.** At the shipped tier the joint object's partition is identical to the
> DNA-only partition on **0 of 5** disagreement families (Wilson 95% [0.000, 0.434]; 0 of 7 overall),
> its only clause with teeth is a **label** (`κ` fires DISCORDANT on 4 of 6 opportunities but scores
> **Fisher p = 0.40** against an outcome it did not produce, and at n = 6 the test **cannot reach**
> 0.05), and the node-level complementarity that motivated the proposal collapses on re-derivation —
> **RNA-only nodes 24 → 8 of 362** Soto members (ceiling **13 = 3.6%**), with **all 13 carrying a
> genome-wide DNA paralog at `1 − de` of 0.9703–1.0000 against a 0.60 floor, and 0 below it**.
>
> **What it forbids, at the level of the partition: NOTHING — the set is empty.** That is why the
> verdict cannot be "definitional". The full argument, every adversarial finding, the limits and the
> costed open items are in **§9** below and in
> `/home/juanfra/winloci_scratch/o1_joint/O1_JOINT_VERDICT.md`.

>
> ⚠ **Sections 0–8 were written BEFORE phase 3 and are retained as the pre-declaration**
> (sha256 `f8996c55…1692a7` freezes F-J1..F-J6). **Four of their claims are corrected in §9.4.**
> Where §0–§8 and §9 disagree, **§9 is the record.**

---

### Implementation note — 2026-08-16

`gw_family_catalog --joint-dna-rna` now materializes the property/certificate verdict
without changing the partition. On the emitted RNA-locus universe it runs the same
nucleotide identity/coverage tiers on two typed substrates:

- RNA: spliced exon-sums in transcription orientation (and, when requested,
  `--rna-forward-only`);
- DNA: complete genomic spans of those loci, transcript-normalized but accepting both
  PAF orientations so inverted structural duplications remain visible.

It writes `<out>.joint_edges.tsv` (`RNA_DNA`, `RNA_ONLY`, `DNA_ONLY`),
`<out>.joint_families.tsv` (connectivity, edge Jaccard, and κ), and
`<out>.joint_rule.tsv` (the typed semantics). The files are reporting-only:
`membership_effect=none`. DNA-only cross-family edges are exposed as hypotheses, not
automatic merges, and RNA-only edges expose possible repeat/extent disagreement.

The node universe is deliberately stated in the rule file: this comparison covers
RNA-detected emitted loci. It cannot discover wholly unexpressed DNA-only nodes;
`--from-genome` and `--dna-family-fallback` remain the complementary mechanisms for
that part of O1.

There is a second type boundary: **a Soto family is an SD cluster, not automatically
a gene family**. DNA-only homology may corroborate a gene-family assignment, and it
may rescue an RNA-null locus that has independent same-family annotation. It must not
recruit an anonymous locus into a gene family merely because both loci share an SD
block. [`bench/o1_gene_family_audit.py`](../bench/o1_gene_family_audit.py) enforces
this distinction in the validation graphs; the production joint files remain
reporting-only and do not pretend to solve gene-family typing.

Example:

```bash
gw_family_catalog \
  --bam reads.bam --fasta genome.fa --out sample \
  --rna-forward-only --joint-dna-rna
```

Status **2026-08-12** (§0–§8), **2026-08-13** (§9). Companion to `docs/seeded_family_definition.md`, whose notation, tier notice and
constants this document inherits without restating. Written **before phase 3 measures anything**; every
number below is either (i) re-derived this session from artifacts on disk, or (ii) quoted from
`seeded_family_definition.md` §P4b with its provenance named. Anything not yet measured is marked
**OPEN**, and the falsifiers of §4 are **pre-declared**, not post-hoc.

Work dir for the phase-1 artifacts cited here: `/home/juanfra/winloci_scratch/o1_joint/strata/`
(`anatomy91/`, `mixed/`, `redundancy/`, `rederive/`, `rna_only/`).

---

## ⚠⚠ RETRACTION BANNER — READ BEFORE QUOTING THE MOTIVATION FOR THIS DOCUMENT

The proposal that DNA and RNA be made **definitionally joint** rested on one published number:

> "at ≥95%: DNA 90.3% / RNA 84.5% / **UNION 97%**, both 282, DNA-only 45, **RNA-only 24**, neither 11 —
> different questions, NOT a ranking; RNA is NOT a subset of DNA" (2026-07-25/26)

**It does not survive re-derivation at the shipped metric forms and is retracted.** Two independent
re-derivations this session, from two different surviving PAFs, agree on sign and magnitude:

| source | predicate | DNA | RNA | union | both | **RNA-only** | DNA-only | neither |
|---|---|---|---|---|---|---|---|---|
| PUBLISHED 2026-07-25 | `nmatch/blocklen ≥ 0.95`, `qcov ≥ 0.9`, `-N50 -p0.5` | 90.3% | 84.5% | 97.0% | 282 | **24** | 45 | 11 |
| `strata/rederive/` (`mm.paf`, sibling-roster) | shipped forms, symmetrised | 97.0% | 84.5%\* | 98.6% | 300 | **6** | 51 | 5 |
| same, DNA link to ANY member | shipped forms | 97.8% | 84.5%\* | 98.6% | 303 | **3** | 51 | 5 |
| `strata/rna_only/` (`gmap.paf`, genome-wide) | shipped forms | 95.0% | 84.5%\* | 98.6% | 293 | **13** | 51 | 5 |

\* **the RNA column was NOT re-derived** in either run — it is the inherited 2026-07-25
`member_attribution.final.tsv` value, computed on a `-M -L` **subset BAM** and itself a union of four
detection legs (only 237/306 = 77% of it the O1 catalog). Label it inherited wherever it appears.

The two re-derivations differ (6 vs 13) **because they use different DNA predicates** — "aligns to a
listed Soto sibling" on `mm.paf` versus "aligns anywhere in the genome" on `gmap.paf` — not because one
is wrong. What matters is that both are far below 24, and that:

⭐ **ZERO of the residual RNA-only members lack a DNA paralog at the shipped identity floor.** Minimum
genome-wide best-paralog identity across the 9 examined per-locus is **0.921**; all 13 in the other run
sit at DNA identity **≥ 0.9661**. As a set of *"loci that exist in RNA and cannot be seen in DNA"*,
RNA-only is **EMPTY**. The residual decomposes into defects of the DNA leg and of the benchmark roster:
singleton Soto families where the sibling predicate is vacuously false (3/9), the M1 coverage clause
against a mis-extended or partial node (4/9 and 6/13), the one-record clause (HERC2: union coverage
0.8919, best single record 0.4887), and roster restriction (1/9).

⚠ The published "divergent copies" story — *"RNA rescues via exon/protein homology BELOW the
sequence-identity floor (median 88%)"* — is an artifact of the **banned `nmatch/blocklen` metric**.
Under `1 − de` the same members read POM121 0.877→0.9686, UBE2Q2P1 0.854→0.9581, PMS2P4 0.892→0.9369,
GUSBP1 0.904→0.9765, ANAPC1 0.942→0.9875, HERC2 0.922→0.9703. On all 1,644 non-self records the bias is
**strictly one-sided**: 361 records (22.0%) pass `1 − de` at 0.95 and fail `nmatch/blocklen`, reverse
count **0**; at the shipped 0.60 floor **1,644/1,644** pass. There is no divergence-rescue stratum.

**Consequence for this document.** The node-level complementarity claim was the last standing empirical
support for making the DNA/RNA joint object definitional **at the level of the family relation**. It is
gone. What survives, and is re-derived, is the *reverse* asymmetry — DNA-only **51**, of which 21 are
`MISS:not-expressed` — which supports a **division of labour**, and forbids a symmetric fold.

---

## 0. THE ANSWER IN ONE PARAGRAPH

**The joint object is not a new partition, and any definition that says it is would be near-vacuous
here.** At the shipped tier `E_RNA ⊆ E_DNA` on the panel where it has been checked most carefully
(27-node human NPIP: 351 / 260 / 260 / **91** / **0**, Jaccard 0.7407), so a definitional **union of
edges is a measured no-op** and a definitional **intersection is exactly `E_RNA`** — and the 91 edges it
would discard are **all redundant** (0 of them are cut edges; 0 of all 351 are; the partition is
identical either way). The mechanism is near-structural: an RNA node's **bases are fetched from the
assembly**, so the edge relation is DNA all the way down. Therefore the definition below keeps the
**edge relation on DNA**, and makes the joint content live where it is real and falsifiable: at the
**node** (RNA decides which loci exist as *expressed* objects and supplies a second extent) and in a
**concordance certificate** on each family (RNA is a held-out test of the DNA partition, never a term in
it). That certificate is the only clause with teeth, its teeth are currently **7 families' worth**, and
§4 states in advance what would take them out.

---

## 1. Definition

Fix a genome `G`, an alignment tier `T` (the shipped `ER_TIER_FLAGS`), and a read set `R` aligned to `G`.

> **Node set `V`.** As in `seeded_family_definition.md` §1: the loci proposed **from the data alone**.
> A node `v ∈ V` is a genomic interval on `G` and carries **two extents and one label**:
>
> * `g(v)` — the **genomic extent**, the interval itself. `seq_g(v) = G[g(v)]`.
> * `x(v)` — the **transcribed extent**, the union of exon blocks supported by `≥ 3` reads of `R` at
>   `v`. `seq_x(v) = G[x(v)]` — the **concatenation of assembly bases at those coordinates**.
>   `x(v) ⊆ g(v)` by construction, and `x(v)` may be `∅`.
> * `e(v) ∈ {expressed, unobserved}` — `expressed` iff `x(v) ≠ ∅`.
>
> ⭐ **`seq_x(v)` contains not one base of read sequence.** RNA supplies **coordinates**; the assembly
> supplies **every base** (`genome.fetch_sequence`). This is why the containment below is near-theorem
> and not luck, and it is the single most load-bearing sentence in this document.
>
> **Support relation `E_r`.** The unordered pairs `{u,v} ⊆ V` for which **some single alignment record**
> between **`seq_g(u)` and `seq_g(v)`**, at tier `T`, has identity `1 − de ≥ τ` and M1 coverage
> `≥ c` of the **shorter** of the two. Shipped `τ = 0.60` (inert), `c = 0.50` (operative).
>
> ⚠ **`E_r` is computed on `seq_g` only.** Not on `seq_x`, and **never on a mixed pair
> `(seq_x(u), seq_g(v))`** — see §3.4 for why mixed pairs are excluded on tier grounds as well as
> definitional ones.
>
> **Family.** The blocks of the γ-quasi-clique refinement (γ = 0.20) of the connected components of
> `(V, E_r)`. Unchanged operator, unchanged constants. **The blocks partition `V`.**
>
> **Joint family.** A family `F` together with:
> * its **expression profile** `(|{v ∈ F : e(v) = expressed}|, |F|)`;
> * its **concordance certificate** `κ(F) ∈ {CONCORDANT, DISCORDANT, UNTESTABLE}` (§2), computed from
>   `E_x`, the same edge rule applied to `seq_x` restricted to `F`.
>
> **The certificate is a LABEL ON a family, never a TERM IN it.** Removing `κ` changes no block.

### 1.1 The division of labour, stated as a table

| supplies | what exactly | enters the definition as |
|---|---|---|
| **assembly / DNA** | every base of `seq_g` **and** of `seq_x` | **the edge relation `E_r`** — the whole of it |
| **RNA reads** | which loci exist as *expressed* objects; the transcribed extent `x(v)` | **node attributes** `e(v)`, `x(v)`, and the **certificate** — never an edge |
| **minimal annotation** | candidate boundaries for `g(v)` | **node construction**, upstream of everything (§5) |

### 1.2 What was considered and rejected, with the measurement that rejected it

| candidate joint rule | measured outcome | verdict |
|---|---|---|
| `E = E_DNA ∪ E_RNA` | `E_RNA ⊆ E_DNA` ⟹ `E = E_DNA` identically (27-node NPIP; 0 RNA-only in 30/30 degraded cells; joint partition = DNA partition in 44/44 cells) — ⚠ **27-NODE-LOCAL, see §9.4-(1)**: at 61 nodes with the rep-transcript substrate union differs from DNA on **2 of 7** | **NO-OP on 27 nodes. Rejected.** |
| `E = E_DNA ∩ E_RNA` | `= E_RNA` identically. Discards 91 edges at identity 0.9709–1.0000. **0 of the 91 are cut edges**; partition **identical**; γ still inert (min RNA degree 11 vs floor 5.2) — ⚠⚠ **and that negative has MEASURED POWER 0** (§9.3): a random deletion of 91 of 351 edges changes the partition **0/10,000** times. At 61 nodes intersection **destroys HERC2** (5 → five singletons) | **Rejected — and the outcome argument IS available at 61 nodes, §9.4-(1)** |
| node rule "member must be expressed" | would delete the **21** re-derived `MISS:not-expressed` DNA-only members | **HARMFUL. Rejected.** |
| edges on `seq_x` instead of `seq_g` | RNA edge count 351→260 on the same nodes; 91/91 losses are coverage, 0 identity | **Rejected — strictly weaker** |
| mixed edges `(seq_x(u), seq_g(v))` | **not expressible at the shipped tier** (`-X` voids two-file runs; `mixed_X.paf` is empty) | **Rejected on tier grounds, §3.4** |

⚠⚠ **CORRECTION to the framing this document was commissioned under.** The brief states that
intersection *"DISCARDS the 91 bridging edges … strictly worse."* That is right about the edges and
**wrong about the outcome**: intersection returns the **same partition** — density 1.000 → 0.741 with
families unmoved. The argument against intersection is a **margin** argument, not an outcome argument,
and the outcome argument is not available on this panel. Do not restate it as one.

---

## 2. What the definition FORBIDS

A definition that rules nothing out is a description. This section is deliberately unflattering.

### 2.1 At the level of the family relation, on the 27-node panel: it forbids NOTHING. Say this first.

| question | answer, measured |
|---|---|
| configurations the joint object declares **not** a family that DNA-only accepts | **none — zero blocks** |
| configurations the joint object accepts that DNA-only rejects | **none** (union is a no-op) |
| edges the strictest joint fold would remove | 91 of 351 — **0 cut edges**, partition unchanged |
| families whose membership changes under any fold tried | **0 of 7**, and 0 of 44 (family × degradation) cells |

⭐ **The vacuity trap was NOT escaped on this panel, and the brief predicted exactly that.** Given the
same loci, DNA and RNA return the identical partition. Every stratum nominated as a *disagreement*
stratum decided nothing: the 91 DNA-only edges are redundant, RNA-only edges are 0 by construction, and
the RNA-only **node** stratum collapsed 24 → 6/13 with none of the residual invisible to DNA.

### 2.2 Where it does forbid something: the certificate, and the panels with heterogeneous node sets

`E_RNA ⊆ E_DNA` is **a node-set effect, not a law.** On the 61-node panel at the **shipped** tier and
the **shipped** single-record rule (`seeded_family_definition.md` §P4b, 2026-08-10):

| family | \|V\| | DNA E | RNA E | shared | DNA-only | **RNA-only** | partition same? | `κ` |
|---|---|---|---|---|---|---|---|---|
| NPIP | 26 | 301 | 157 | 150 | 151 | **7** | no | **DISCORDANT** |
| TBC1D3 | 9 | 36 | 36 | 36 | 0 | 0 | yes | CONCORDANT |
| RABL2 | 2 | 1 | 1 | 1 | 0 | 0 | yes | CONCORDANT |
| APOBEC3 | 3 | **0** | **1** | 0 | 0 | **1** | no | **DISCORDANT** |
| MAGEA | 12 | 62 | 66 | 62 | 0 | **4** | yes | **DISCORDANT** |
| GSTM | 4 | 1 | 2 | 1 | 0 | **1** | no | **DISCORDANT** |
| HERC2 | 5 | 7 | **0** | 0 | 7 | 0 | no | **UNTESTABLE / DISCORDANT** |
| **total** | **61** | **408** | **263** | **250** | **158** | **13** | **3/7** | **4 fail** |

So the certificate is **not vacuous where node sets are heterogeneous**: it fires on **4 of 7** families,
and APOBEC3 is the sharp case — **DNA emits 0 edges and RNA emits 1**, i.e. RNA asserts a copy relation
DNA's own sequence does not support at `c = 0.50`. Under the definition above that is **not** a merge:
the family stays as DNA drew it and is **stamped DISCORDANT**. That is the whole of what the joint
object forbids — it forbids *emitting such a family as certified*.

> **The forbidding clause, stated once:**
> A family may not be emitted as **CONCORDANT** if `E_x|F ⊄ E_g|F`, or if the partition induced on `F`
> by `E_x` differs from the trivial one. Discordance is **flagged, never acted on** — it may not add,
> delete, split or merge a block.

### 2.3 What an RNA-primary definition would forbid that this one does not

Stated so the asymmetry is on the record and cannot be read as a preference:

* it would **delete 21** re-derived `MISS:not-expressed` members (untranscribed true paralogs);
* it would lose the **91** DNA-only edges at identity 0.9709–1.0000, and the **151** on NPIP@61;
* it would lose **HERC2 entirely** (DNA 7 edges, RNA 0).

An RNA method cannot detect an untranscribed feature. That is a statement about the substrate, not a
counterexample to containment — it is *consistent* with it, and it is why the fold is asymmetric.

### 2.4 The mechanism, corrected

The brief attributes the 91 lost edges to *"the RNA node is a smaller object"*. **Re-derived this
session, that is refuted.** Decomposing `cov_dna/cov_rna = (aln_dna/aln_rna) × (len_rna/len_dna)`:

| factor | min | median | max | direction |
|---|---|---|---|---|
| **length** `len_rna/len_dna` | 0.598 | **0.808** | 0.914 | **< 1 ⟹ HELPS the RNA node** (lower absolute bar) |
| **span** `aln_dna/aln_rna` | 2.221 | **2.788** | 3.929 | **> 1 ⟹ kills it** (91/91 pairs) |

The counterfactual the brief asks for — RNA's aligned span scored against the **DNA** node's length —
gives coverage 0.2183/0.3441/0.4313 and passes **0/91**: it is *strictly worse*, because the DNA
denominator is larger. There is no reading of "give the RNA node the DNA node's length" that recovers
anything.

⭐ The real mechanism is the **single-record clause meeting two spliced endpoints**. Lifting only that
clause (union of identity-passing RNA records, RNA denominator unchanged) gives coverage
0.7448/**0.9205**/1.0000 and passes **91/91**, from a median of 79 identity-passing records per pair.

⚠⚠ **THAT IS A DIAGNOSIS, NOT A PROPOSAL, AND IT IS THE MOST MISQUOTABLE NUMBER IN THIS FILE.** Summing
coverage across records is `RUSTLE_ER_SUM_COVERAGE`, **default OFF**, condemned in `bench/soto/rustlib.py`
("two loci sharing 60% in four scattered blocks are sharing fragments, not a gene") and carrying a
measured cross-family precision cost. Nothing in this definition relaxes it.

The corrected unifying statement — which makes the blind result and the DNA/RNA result the same theorem:

> **An `E_r` edge needs at least one endpoint that is CONTIGUOUS over the shared region.** Over-long or
> fused nodes are *free* (they hand the denominator to the partner); short-but-content-complete nodes
> (RNA nodes) are always the denominator yet still pass **iff the partner is contiguous**; RNA×RNA is
> the one configuration where **neither** endpoint supplies contiguity.

---

## 3. Where the O1 ⊥ O2 line falls, and the proof that this object stays on the O1 side

**The rule.** Reads as **SUPPORT** for a node is O1. Reads **ASSIGNED** to a copy is O2. `E_c` (read
conflict) and any per-read copy decision may not enter the definition of membership.

**Where the line falls here.** Reads enter this definition in exactly two places, both *within a single
locus*:

1. **node existence / `e(v)`** — "≥ 3 reads support exon blocks at `v`";
2. **`x(v)`** — the coordinates of those blocks.

Neither is a statement about *which of several copies* a read came from. `E_r` is then computed on
assembly bases; `E_c` appears nowhere; no matching, no facility location, no bipartite assignment.

**The invariance that makes it a proof rather than a claim.** Define the support predicate over
**alignment RECORDS, not assignments**: a read supports `v` if it has *an* alignment record at `v`, and
the same read may support many loci simultaneously. Under that predicate `x(v)` and `e(v)` are functions
of the alignment *set* alone, hence **invariant under every assignment consistent with the alignments** —
including the identity of the primary flag. An object invariant to all assignments cannot be encoding
one. That is the O1-safety argument, and it is a proof only of the record-counting predicate.

⚠⚠ **THE ONE PLACE THIS COULD BREAK, DISCLOSED RATHER THAN BURIED.** The shipped RNA node construction
filters `-F 2308` — **primary records only** — and the primary flag is itself a choice among tied
placements (MAPQ 0). So the *shipped implementation* of `x(v)` is assignment-dependent in principle even
though the *definition* above is not. This is not hand-waved away by the standing 95.9%-provably-invariant
and 99.10% max-weight-matching results: those are about *loci*, not about exon blocks. **F-J4 (§4) is
the pre-declared test, and until it is run the O1-safety argument holds for the definition and is OPEN
for the implementation.**

### 3.4 Mixed edges are excluded — two independent reasons

A mixed edge `{seq_x(u), seq_g(v)}` was measured this session and **recovers 85–88 of the 91**. It is
nonetheless excluded:

1. **It is not expressible at the shipped tier.** `-c -X --no-long-join -k 11 -w 5` on a two-file run
   emits **zero** records (`-X` ⟹ `--dual=no` skips self *and dual* mappings); `mixed_X.paf` is empty.
   The 88/91 comes from a translation `-P -N 50 -p 0` which gives **269** edges vs the stored **260** on
   the validation control (symdiff 11) — *more permissive than shipped*. The default-secondary
   translation gives **19/91**. The result is tier-sensitive across a 19–88 range and **88/91 must never
   be quoted as a shipped-tier number.**
2. **It is not the relation `E_r` is defined as.** `E_r` is symmetric on sequences of the same kind. A
   mixed edge asserts a relation between one locus's *transcribed* extent and another's *genomic* extent;
   its meaning under P1–P4 / C1–C3 is undefined.
3. And operationally it buys nothing: every one of the 91 it recovers **was already a DNA edge**, and all
   91 are redundant.

---

## 4. THE FALSIFIERS — pre-declared, before phase 3

Each names the observation that refutes it, the stratum it must be sought on, and where available the
prediction this document is committing to.

**F-J1 — the vacuity falsifier (the definition's own content).**
The certificate `κ` is the only clause with teeth. *If, on a panel of ≥ 5 families with ≥ 2
read-supported nodes each, `κ` returns **CONCORDANT** for every family, then the joint object forbids
nothing and must be **withdrawn as a definition** and retained as a diagnostic label only.* Opportunities
must be counted and reported (`families with ≥ 2 expressed nodes`), because a certificate that never
fires for lack of opportunity is not the same as one that never fires.
**Prediction, from §2.2: 4 of 7 fail on the 61-node panel.** If a re-run at the shipped tier returns
0 of 7, this document is refuted.

**F-J2 — the containment falsifier.**
*An RNA-only edge joining two nodes that `E_DNA` places in **different blocks**, reproduced at the
shipped tier, refutes "the edge relation is DNA all the way down".* Containment would then be false in
the strong sense and the intersection fold would no longer equal `E_RNA`.
**Prediction: this exists — APOBEC3 on the 61-node panel (DNA 0 edges, RNA 1).** It must be re-measured
at the shipped tier on that exact node set before it is quoted. If it reproduces, §1.2's row
"intersection = `E_RNA` identically" is *false off the 27-node panel* and must be narrowed there.

**F-J3 — the harm falsifier for the rejected RNA-required node rule.**
*If requiring `e(v) = expressed` for family membership deletes zero true members on a scored panel, the
rejection in §1.2 is unjustified.* **Prediction: it deletes ≥ 21** (`MISS:not-expressed`).

**F-J4 — the O1/O2 falsifier (the one that would be fatal).**
*Replace the `-F 2308` primary-only support predicate with the record-counting predicate of §3. If any
family's **membership** changes, then `V` depends on read assignment, the O1-safety proof fails as
stated, and the joint node must be reported as straddling the O1/O2 line.* A change in `x(v)` extents
alone is expected and is not a failure; a change in **blocks** is. **This document commits to reporting
a positive result here as fatal to itself, not as a nuance.**

**F-J5 — the redundancy falsifier (does any of this matter?).**
The Q4 result ("the 91 decide zero families") was obtained on a panel whose DNA graph is **complete**
(351/351, density 1.000) — the worst possible substrate, where redundancy is guaranteed a priori rather
than discovered. *If, on a family whose DNA graph is **sparse** (the standing genome-wide figure is 18%
of co-family pairs aligning against the ~10% a spanning forest needs — a 1.8× margin), deleting the
RNA-lost edges **disconnects** a family, then the intersection fold changes an outcome and §2.1's "it
forbids nothing" is local to NPIP.* This is the single most informative measurement phase 3 can make.

**F-J6 — the certificate's own precision.**
*If `κ = DISCORDANT` fires on families that are independently correct as often as on families that are
independently wrong, the certificate is noise.* It must be scored against a truth that did not produce
it. Note that on the 61-node panel `κ` fires on MAGEA, whose partition is nonetheless the **same** under
both graphs — so discordance at the **edge** level and at the **partition** level must be reported as two
different outcomes, never merged.

---

## 5. What is INHERITED, not new

### 5.1 P1 (seed-invariance) holds for the joint object, and the proof carries over verbatim

The P1 theorem in `seeded_family_definition.md` §1 has exactly three premises: (a) `V` is **seed-free**;
(b) the refinement starts from `all_components` **including singletons** (`family_split.rs:480`) and only
ever **splits**; (c) therefore the blocks **partition** `V` and `F(s)` is a lookup.

Each premise survives the joint construction, and the reason is that **nothing added is a term in the
relation**:

* `x(v)`, `e(v)` and `κ(F)` are functions of the data at `v` (or of `F` after the fact). None of them
  mentions a seed. `V` remains seed-free ⟹ (a) holds.
* The refinement operator is **unchanged** — same γ = 0.20, same `all_components` entry point, same
  split-only behaviour ⟹ (b) holds.
* `κ` is computed **after** the blocks exist and cannot move a node ⟹ (c) holds.

⟹ **P1 is inherited as a theorem, not re-measured.** For all `s' ∈ F(s)`, `F(s') = F(s)`. This is worth
more than any measurement in this run, and it is worth exactly as much as the guarantee that `κ` never
feeds back. **The moment discordance is allowed to split or merge a block, the refinement is no longer
split-only from `all_components`, and this proof is void.** That is the structural reason the certificate
is a label.

### 5.2 What is inherited as a *limitation*

* **Extent is not seed-invariant** (P1's NPIPB8 residual), and the joint node makes this **worse, not
  better**: there are now *two* extents per node, and `x(v) ⊆ g(v)` inherits every boundary error in
  `g(v)` plus its own.
* **`τ` is inert; `c` is the only operative constant** — unchanged, and the joint object adds no new
  constant. The `≥ 3 reads` support threshold is inherited from node construction, is not new here, and
  is **not** an edge constant.
* **The tier is not scope-invariant.** Every number in this document is panel-scope. Do not set any of
  them beside a genome-scale number. (`mid_occ` nearly halves between an 80.6 Mbp and a 38.6 Mbp scope.)
* **`E_c` remains a held-out test of `E_r`, never a term in it** — the joint object changes nothing here,
  and the certificate of §2.2 is the *same discipline applied to a second sequence channel*: RNA is to
  the DNA partition what `E_c` is to `E_r`.

---

## 6. Reconciliation with the minimal-annotation reframe

Both this document and the A0–A3 ladder (`/home/juanfra/winloci_scratch/o1_minann/`) act on the **same
node object**, from opposite sides: annotation proposes `g(v)`, reads carve `x(v)` out of it. They must
agree, and they do, on one mechanism and one ordering.

**They agree on the mechanism.** Of within-family pairs that have an alignment record and still fail:

| corpus | coverage-only failures | identity failures |
|---|---|---|
| blind run (`o1_blind`) | 245/245 | **0** |
| A3, chr1, `c=1.00, d=0` | 152/152 | **0** |
| A3, every completeness cell `c = 1.00 … 0.10` | 38/38, 5/5, 1/1 | **0** |
| this run, 27-node NPIP DNA-only edges | 91/91 | **0** |
| Soto member re-derivation | invariant from id 0.90 → 0.60 | **0** |

⭐ **Identity never fails. Extent binds, through the coverage clause, on every substrate looked at.** The
DNA/RNA gap and the blind delineation gap are the same failure mode; §2.4 gives the corrected form of
that unification (*contiguity*, not *length*).

**They agree on the ordering: completeness binds, accuracy is nearly free.** A3, chr1 coarse, τ = 0.50,
denominator fixed pre-run:

| axis | cells | effect on recovery |
|---|---|---|
| **completeness** `c` (fraction of gene intervals retained) | 1.00 / 0.50 / 0.25 / 0.10 | **0.5500 → 0.3500 → 0.1000 → 0.0250**; gain over the annotation's own partition +0.1000 → +0.0250 → **exactly 0.0000** (CI `[+0.0000,+0.0000]`) |
| **accuracy** `d` (endpoint displacement, bp) at `c = 0.50` | 0 / 1000 / 5000 | **0.3500 / 0.3500 / 0.3500** (QC); raw CC 0.3500 → 0.3250 only at `d = 5000` |

⚠ Two qualifications that must travel with that ordering, both from the A3 run itself:

* **"accuracy is free" has a ceiling.** At `d = 1000` the graph is *starved* (50 → 35 edges, within-family
  share holds 14.0% → 17.1%); at `d = 5000` it **floods** — edges rise 3.8× (50 → 189) while
  within-family edges *fall* (7 → 5) and the within-family share collapses to **2.6%**, with `CC ≢ QC`
  for the first time. Displacement past the kilobase scale does not merely lose edges, it **manufactures**
  them out of flanking repeat.
* **The knee is a hand-off, not a point.** Between `c = 1.00` and `c = 0.50` the *method's* contribution
  decays (gain capture 0.182 → 0.100); below `c ≈ 0.25` the *annotation* stops permitting an answer at
  all (oracle ceiling 40 → 23 → 5 → 1). Only the first is a method limit.

**Where the two runs meet, concretely.** The minimal-annotation reframe says: an `E_r` edge needs **one**
well-delineated endpoint, and an annotation hands over precisely that. This run says: what "well
delineated" means is **contiguous over the shared region**, and `x(v)` — the transcribed extent — is the
one channel that is *content-complete but never contiguous*. ⟹ **RNA is the right instrument for deciding
that a locus is a real expressed object, and the wrong instrument for supplying the interval that gets
aligned.** That is the division of labour of §1.1, arrived at from the other direction.

---

## 7. Known exposures (state these first, do not wait to be asked)

1. **The RNA column of every node-level table is INHERITED, not re-derived** (2026-07-25,
   `-M -L` subset BAM, four-leg detection union). Only the DNA column moved.
2. **`E_RNA ⊆ E_DNA` is a 27-node result.** On the 61-node panel there are **13** RNA-only edges. The
   containment is a *node-set* effect; the definition above does not depend on it (it never takes a
   union or an intersection), but §1.2's arithmetic does.
3. **"The 91 decide zero families" was measured on a complete graph** (density 1.000) — redundancy is
   near-entailed by that panel choice. F-J5 is the live version of the question.
4. **Control B is a standing warning about the tier, unresolved.** An RNA node aligned to *its own* DNA
   node is definitionally a spliced subsequence and should score 1.0000 at 27/27. It scores median
   **0.9344** and passes only **25/27**, minimum 0.4498. Any RNA-vs-genomic claim anywhere in the thesis
   sits on top of that floor.
5. **`κ` has never been scored against an independent truth** (F-J6). It is defined here, not validated.
6. **Single panel, one gene family, one species** for everything in §2.1 and §2.4.
7. **The `-F 2308` exposure of §3** — the definition is assignment-invariant, the implementation is not
   yet shown to be.
8. **Nothing in this document was measured at genome scale**, and `-k 11` exhausts its seed alphabet
   there.

---

## 8. Summary judgement, for the reader who reads one section

The DNA/RNA joint object **can** be defined cleanly, it **inherits P1 as a theorem**, and it **stays on
the O1 side of the line** — but on the panel where it has been measured most carefully it **forbids
nothing about families**, because `E_RNA ⊆ E_DNA`, the 91 differing edges are redundant, and the
node-level complementarity that motivated the whole proposal **did not survive re-derivation** (RNA-only
24 → 6/13, none of the residual invisible to DNA). The honest form of the object is therefore **not** a
joint *relation* but a **two-channel node plus a falsifiable concordance certificate**, whose entire
content is the 4-of-7 families it currently flags on a heterogeneous node set. If F-J1 returns 0 of 7 at
the shipped tier, the correct action is to withdraw this as a definition and keep it as a diagnostic —
and that outcome would be a better result than dressing up a no-op.

---

# 9. PHASE 3 — THE MEASUREMENT, AND THE VERDICT (2026-08-13)

Everything above this line was written **before** phase 3 and is the pre-declaration. This section is
the record. Full verdict, adversarial ledger and costed opens:
`/home/juanfra/winloci_scratch/o1_joint/O1_JOINT_VERDICT.md`.

## 9.0 VERDICT: **PROPERTY, NOT DEFINITION** — ⭐ **RE-TESTED ON THE SHIPPED OBJECT 2026-08-13 (O-2): STANDS, RESTATED**

`E_RNA ⊆ E_DNA` and the DNA/RNA node relation are **properties of the shipped O1 definition** — real,
measured, worth stating in the thesis. They are **not clauses of it**, and the joint object as defined
in §1 must be **retained as a diagnostic**, exactly as §8 anticipated.

⭐⭐⭐ **O-2 DONE — the folds have now been measured on the object the binary emits**
(`/home/juanfra/winloci_scratch/o2_shipped_union/O2_SHIPPED_UNION.md`). The verdict's **content**
survives; **three of its sentences must be re-quoted and one number is retracted.**

1. ⭐ **THE UNION IS REAL, IT FIRES, AND IT LANDS ON DNA.** The shipped `E_r` at the `refine` site is
   `E_x ∪ E_g`, confirmed by counter: on a 61-node single-family call the additive leg adds
   **101 edges = exactly `|E_dna \ E_rna|`**. On the human `--cross-chrom` panel (**D1 = 26** input
   copysets, printed before any rate) it added **288 edges on 322 core edges**, took components
   **52 → 32** and moved **11/26 = 0.4231** [0.2554, 0.6105] partitions — yet block sets
   shipped(`E_x ∪ E_g`) vs DNA-only(`E_g`) differ **0/26** [0.0000, 0.1287] and the **emitted catalogs
   are identical** (26 families / 143 copies, ARI 1.0000, 0 forbidden pairs in both directions).
   ⟹ **a no-op relative to DNA-only, never relative to the shipped core.**
2. ⚠⚠ **THE DIRECTION WAS WRONG.** This document measured *adding RNA to DNA* (9 RNA-only edges of 316
   at 61 nodes on the binary, 0/7 partitions). The shipped union runs the **other way** — it adds DNA
   to an RNA core — and that direction moves **101 edges and 2/7 families** at the same node set.
3. ⚠⚠ **"0 / 5, 0 / 7" IS AN ENTAILMENT, NOT A MEASUREMENT** — §9.1 says so itself ("JOINT's partition
   equals DNA's **by construction**"), but the figure has been read as a measured union result. Apply
   the **shipped gated-union rule** to this run's own `gate.json` edge sets and the forbidden set is
   **not empty**: **2/5 disagreement, 2/7 overall (APOBEC3, GSTM)** at the rep-transcript substrate,
   **0/5, 0/7** at pooled read exons — independently reproduced by a third parser on the stored sets.
   §9.4-(1) already carried the 2/7; the two numbers are different objects and must stop travelling as one.
4. ⚠⚠ **`E_RNA ⊆ E_DNA` IS FALSE ON THE BINARY past 27 nodes**: exon-sum-only edges **0 of 351** at 27,
   **9 of 316** at 61 (pooled) / 13 of 263 (rep), **66 of 696** at 244 human reps. It degrades
   monotonically with the node set and may never be quoted without one.
5. ⚠ **NEW DEFECT — `E_r` IS NOT INVARIANT TO NODE NAMING.** `-X` implies `--dual=no`; renaming reps to
   integer indices moves the 27-node RNA edge set by **8 of 260 = 3.1%** on byte-identical sequence
   (`-t 2` vs `-t 4` symdiff 0, so it is not threading). Every edge count in this document carries ~3%
   naming noise; the corpus's 351 / 260 / 91 read **351 / 258 / 93** off the binary.

## 9.1 What it forbids — the set is EMPTY at the level of the partition

Pre-registered test, `o1_joint/test/arm2/` (`PREREG.md` hashed **before** any score; acceptance gate
reproduced `seeded_family_definition.md`'s 61-node table **bit-for-bit** — DNA 408 / RNA 263 /
shared 250 / DNA-only 158 / RNA-only 13 and all seven per-family rows — and `s91.py` asserted
351/260/260/91/0 on the 27-node panel, so artifact drift would have failed loudly).

| comparison | families whose BLOCK SET differs | Wilson 95% |
|---|---|---|
| **JOINT vs DNA, disagreement stratum** | **0 / 5** (0 / 7 overall) | **[0.000, 0.434]** |
| JOINT vs RNA | 4 / 5 (NPIP, APOBEC3, GSTM, HERC2) | [0.376, 0.964] |
| UNION vs DNA | 2 / 5 (APOBEC3, GSTM) | [0.118, 0.769] |
| INTER vs DNA | 2 / 5 (NPIP, HERC2) | [0.118, 0.769] |

⚠⚠ **PREREG §6 forbids reading "JOINT = DNA" as a win**, and it is not one: JOINT's partition equals
DNA's **by construction**. The measurement establishes only that the construction is not accidentally
violated on the stratum where it could have been. **Verdict: `no-difference`, measured on the
disagreement stratum** — the strongest form of that statement available.

**The one non-empty clause is `κ`, and its precision is untestable.** `κ` has 6 opportunities (HERC2
UNTESTABLE, one RNA-bearing node) and fires **DISCORDANT on 4 of 6** [0.300, 0.903] — F-J1's
prediction of "4 of 7" **held**. Scored against an outcome it did not produce (DNA partition ≠ RNA
partition): fired&differ 3 · fired&same 1 (MAGEA, edge-level only) · notfired&differ 0 · notfired&same
2 ⟹ **Fisher exact two-sided p = 0.40**, and at n = 6 the **smallest attainable p is 0.10**. F-J6 is
**not answerable at this panel size**. `κ` is defined, not validated.

## 9.2 Where they disagree — and the re-derivation that removes the proposal's basis

⭐⭐⭐ **THE SHIPPED-TIER JOB THAT PHASE 1 DEFERRED WAS RUN** (`o1_joint/lens_vacuity/`):
`minimap2 -c -X --no-long-join -t 4 -k 11 -w 5 members.fa members.fa`, 1,659 s wall / 5,761 s CPU,
peak RSS 4.00 GB, 1,348,077 records / 1,306,867 usable non-self, ONE index over 11,264,772 bp,
588,643 distinct minimizers, `mid_occ` 752 — single scope; it may sit beside the `mm.paf` number and
beside **no genome-scale number**.

| leg | tier | DNA | RNA | union | both | **RNA-only** | DNA-only | neither |
|---|---|---|---|---|---|---|---|---|
| PUBLISHED 2026-07-25 | `-N50 -p0.1/-p0.5`, banned metric | 90.3% | 84.5% | 97.0% | 282 | **24** | 45 | 11 |
| old tier `mm.paf`, shipped rule | `-c --eqx -N50 -p0.1` | 97.0% | 84.5%\* | 98.6% | 300 | **6** | 51 | 5 |
| genome map-back `gmap.paf`, shipped rule | `-N100 -p0.02` | 95.0% | 84.5%\* | 98.6% | 293 | **13** | 51 | 5 |
| ⭐ **SHIPPED tier + shipped rule** | `-c -X --no-long-join -k11 -w5` | **96.4%** | 84.5%\* | **98.6%** | **298** | **8** | **51** | **5** |
| same, any-member predicate | SHIPPED | 97.5% | 84.5%\* | 98.6% | 302 | **4** | 51 | 5 |

\* INHERITED, never re-derived (2026-07-25, `-M -L` SUBSET BAM, four-leg detection union of which only
237/306 = 77% is the O1 catalog). It is identical in every row because it never moved.

⚠⚠ **THE PHASE-1 PREDICTION WAS REFUTED.** Both phase-1 reports stated the shipped tier "can only ADD
records, so DNA can only rise and RNA-only can only fall." Measured: **DNA FALLS 351 → 349, RNA-only
RISES 6 → 8.** `-X`/`--no-long-join` is not a superset of `-N50 -p0.1` — suppressing long joins
shortens the aligned span, so pairs drop below the M1 0.50 floor. GOLGA8A is rescued; **AC137800.1,
POM121 and POM121C are newly lost, all three on COVERAGE** (0.4460 / 0.2875 / 0.2875) at identity
0.9314–0.9867. **Quote 8. Never 6, 13 or 24. Never reuse the one-sided-correction argument.**

⭐ **THE CEILING.** Re-deriving the RNA leg can only move members between cells, so RNA-only is bounded
by `RNA-only + neither` = **13 of 362 (3.6%)**. All 13 recomputed as `1 − de` from `gmap.paf`:
AC119751.3 0.9970 · AC137800.1 0.9970 · AC243829.6 1.0000 · AL590399.2 1.0000 · BOLA2B 1.0000 ·
CASTOR2 0.9949 · CR381670.1 0.9975 · CTSLP3 0.9927 · DEFB104B 0.9927 · DUX4L50 0.9790 · HERC2 0.9703 ·
POM121 0.9832 · POM121C 0.9881 — **min 0.9703, below the 0.60 floor: 0, with no non-self record: 0.**
⟹ **"loci RNA sees and DNA cannot" is EMPTY against every reachable version of the RNA leg**, by a
margin of 0.37 in identity. The retraction banner at the top of this document is upheld and
strengthened at the shipped tier.

**The disagreement strata and their sizes** (all HUMAN / CHM13 v2.0, panel scope):
RNA-only **edges** 0/351 at 27 nodes and 13/408 at 61 · DNA-only **edges** 91/351 and 158/408 ·
RNA-only **nodes** 0/80 on the family panel and 8/362 (ceiling 13) on Soto · DNA-only **nodes** 13/80
and 51/362 · families where the two sides disagree at all **5 of 7** (50 of 61 nodes).

## 9.3 ⚠⚠ THE 27-NODE NEGATIVE HAS MEASURED POWER ZERO

Pre-declared null: delete `|E_g \ E_x|` edges uniformly at random from `E_g`, 10,000 draws, seed
20260812.

| panel | \|E_g\| | lost to RNA | RNA-lost deletion changes the partition? | **P(random deletion of same size does)** |
|---|---|---|---|---|
| **27-node NPIP** (complete, edge connectivity 26) | 351 | 91 | **no** | **0/10,000 = 0.00000 [0.00000, 0.00038]** |
| **61-node NPIP** | 301 | 151 | **yes** | **0.0002 [0.0000, 0.0011]** |
| 61-node HERC2 | 7 | 7 | yes | 1.0000 — **degenerate**, no information |
| TBC1D3 / RABL2 / APOBEC3 / MAGEA / GSTM | — | 0 | no | n/a |

⟹ **"0 of the 91 are cut edges" is true and UNINFORMATIVE, and may never again be quoted as evidence
that RNA's edge losses do not matter.** §2.1 and Exposure 3 suspected this; it is now measured.

⭐⭐ **AND WHERE THE TEST HAS POWER, IT FIRES AGAINST RNA.** On 61-node NPIP the *RNA-determined*
deletion of 151 of 301 edges disconnects the family; a random deletion of the same size does so with
probability 0.0002 (permutation p ≈ 2 × 10⁻⁴). RNA's losses are **not a random thinning**. The two
nodes intersection cuts off are `L~chr16_22777929_22814315` (NPIPB5; 7 blocks, 26,506 spliced bp,
**293 reads**; DNA degree 23/25, **RNA degree 0/25**) and `L~chr16_75785926_75819336` (9 blocks,
20,921 bp, **1,243 reads**; DNA 23/25, **RNA 0/25**) — heavily expressed, aligning to nothing on the
RNA side. ⚠ **Exactly ONE of seven family cells is an informative null.** That is the honest
denominator, and one cell is not a result (see F-J5 below).

## 9.4 FOUR CORRECTIONS TO §0–§8

1. ⚠⚠ **§1.2's "union is a NO-OP" and its correction "the outcome argument against intersection is not
   available" are BOTH 27-NODE-LOCAL.** At 61 nodes with the shipped rep-transcript RNA substrate,
   union differs from DNA on **2 of 7** (APOBEC3 1,1,1 → 2,1; GSTM 2,1,1 → 3,1) and intersection
   differs on **2 of 7** (NPIP 26 → 24,1,1; **HERC2 5 → five singletons**). **The outcome argument IS
   available — intersection destroys a family.** Stop saying it is not. Both statements must always
   travel with their node set.
2. ⚠⚠ **"The RNA graph" is not one object, and this is the single biggest caveat in the run.** Same 61
   nodes, same rule, same tier: **rep transcript 263 edges (the shipped substrate) vs pooled read
   exons 316 (+53)**, and the RNA partition differs between the two on **4 of 7 families**. Under
   pooled exons, RNA's partition equals DNA's on **6 of 7** and union is a **clean no-op (0 of 7)**,
   while intersection differs on **1 of 7 (HERC2 only)**. ⟹ **the only finding that survives BOTH RNA
   substrates is that intersection damages HERC2.** Every RNA edge number in this document is a
   **rep-transcript** number.
3. ⚠ **§2.3's framing "an RNA method cannot detect an untranscribed feature" does not describe this
   panel's DNA-only loci.** All **13 of 13** carry spliced reads (0/13 have none); **8/13** have zero
   read-supported exon blocks, but **5/13** have ≥ 1 exon block and were still missed by the RNA
   catalog; only **8 of 80** panel nodes are `e(v) = unobserved`. On this panel a DNA-only locus is
   usually a **delineation failure of the RNA side**, not silence. F-J3 is confirmed in **sign, not
   magnitude**: on the fixed 80-node denominator the expressed-member rule deletes **7 true members**
   [0.043, 0.170]; the pre-declared "≥ 21" belongs to the **Soto 362-member** denominator and must
   never be set beside 80. ⚠ the 80-node graph is **PANEL-tier** (`-k11 -w5 -c --eqx -N 200 -p 0.02`,
   no `-X`), a permissive superset, and may not be quoted beside the 61-node shipped-rule numbers.
4. ⚠⚠ **§1's "`E_r` is computed on `seq_g` only" and §1.2's "edges on `seq_x` — Rejected, strictly
   weaker" are BOTH FALSE OF THE SHIPPED BINARY.** `denovo_pipeline.rs:3388-3425` (`family_refine`,
   the O1 refinement path) runs its core tier on the **exon-sum (spliced)** sequence
   (`include_introns` defaults false, `:3231`) and then, gated only on `!edges_connect_all(...)`,
   **unions in an additive GENOMIC-SPAN tier, DEFAULT ON** (`TIER_GENOMIC`, `:3418-3421`). The shipped
   `E_r` at that site **is `E_x ∪ E_g`**. ⚠⚠ Its only cited justification (`:3409`) is
   `bench/FALSE_NEGATIVES.md`, **which was DELETED from the tree** (⚠ corrected 2026-08-13 by O-4: this
   line read "does not exist anywhere in the tree", which wrongly implied it was never written — it was
   committed in `4586ba8` alongside the tier and deleted in `9b0814f`'s docs consolidation, whose message
   claims "git retains all"; it is now RESTORED verbatim and a test asserts the citation resolves) — the string occurs once,
   in that comment. ⟹ **the union fold this document rejects as a measured no-op is what the code
   already computes**, and every phase-1/2 no-op measurement was made on externally built
   `gvg_dna.paf` / `gvg_rna.paf` graphs, **not on the shipped union**. Re-measuring the folds against
   the shipped union is the single most consequential unrun item (§9.7 O-2).
   ⭐⭐ **CORRECTED 2026-08-13 BY O-2 — THIS ITEM IS HALF WRONG, IN THE DIRECTION THAT MATTERS MOST.**
   `family_refine` is **opt-in on the shipped O1 catalog** (`refine_enabled = refine_flag || !o1_homology`,
   `gw_family_catalog.rs:423`), so on the default O1 path **`refine` is never called and the additive
   genomic tier is unreachable** — 0 `.refine.*` artifacts over 25 gorilla control regions and over the
   244-rep human panel. It is default-ON for `copy_assign` and the legacy conflict catalogs, and that is
   where it fires. ⚠⚠ **AND THE REAL DEFECT IS THE OPPOSITE ONE:** the shipped O1 `E_r` is computed on
   the **EXON-SUM**, i.e. the object this document calls **"RNA"** — `homology_genomic_span` is default
   **OFF**, so §1's *"`E_r` is computed on `seq_g` only"* is false because the binary uses `seq_g`
   **never** at that site, not because it also unions it in. Measured on the binary over 244 human reps,
   node sets identical across arms: `|E_x| = 696` vs `|E_g| = 1575`, shared 630, **exon-only 66**,
   genomic-only 945 — **neither contains the other** — block sets differing on **6 of 7** panel families
   [0.4869, 0.9743], **ARI 0.4702** (size-matched null mean −0.0001, 97.5th pct 0.0213, n = 2000).
   ⟹ **"DNA gives the relation, RNA gives expression" is a statement about a graph O1 does not ship by
   default.** ⚠ A degree-sequence-matched null moves a family's partition with **P = 0.66–1.00**, so
   *"the additive tier moved a partition"* is uninformative on its own; **where it lands is** — on
   `E_g`'s partition, 26/26.

## 9.5 F1–F2 — the two findings that decide the verdict

⚠⚠ **F1 (FATAL to the O1-safety claim of §3): `O1 ⊥ O2` is ALREADY FALSE on the default O1 catalog
path, at the NODE.** `detect_homology_catalog_genome_wide` (`denovo_pipeline.rs:2964`) calls
`locus_unique_mapper_counts` — defined in **`read_conflict.rs:267`, the `E_c` module** — sets
`distinguishing_uniq` (`:3053`) and passes it to `distinct_locus_reps(..., cfg.locus_min_reads())`
(`:3116`, **not** behind `--refine`), where `reads_distinguish` is documented as *"the χ(H) edge
predicate restricted to a co-located pair"*; and `cfg.locus_min_reads()` (`:172`) literally
`return self.conflict.min_reads` — **O2's `RUSTLE_CONFLICT_MIN_READS`, default 3, renamed.** So whether
two overlapping candidates are ONE node or TWO, and hence whether a family clears `--min-copies`, is
decided by **how many reads the aligner placed with MAPQ > 0**. Two further leaks: `V` is built from
**primary alignments only** (`RUSTLE_FLAGFREE_SITES` default OFF; the primary flag *is* the bipartite
assignment at 99.10%), and this proposal's own extent leg is prototyped as `locus_confident_extent`
(`:609`) — *"extent defined by only the reads that CANNOT have come from a paralog"* — gating on
**`de`**, the literal O2 criterion (opt-in). The guard test
`homology_catalog_never_touches_the_conflict_graph` (`:5756`) bans four strings and `cfg.conflict.`;
re-running its exact 311-line slice shows `locus_unique_mapper_counts(`, `reads_distinguish` via
`distinct_locus_reps(`, `build_read_placements(`, `locus_min_reads(` and `aligned_reads_from_bam(` all
**PRESENT**. **The guard certifies spelling, not semantics** — the 2026-08-11 `--cross-chrom`
precedent recurring. ⚠ Materiality against the finding: on chr1 the same-strand MAPQ rule fired **0
times over 451 loci**. Inert-on-chr1 is a property of the data, not of the definition. **What survives:
the EDGE relation `E_r` is clean** — sequence-only, reference bases at read-derived coordinates, no
read base and no read assignment in any comparison. **The violation is entirely in `V`** — which loci
exist, how many, and how long. ⟹ **A definitional joint object would PROMOTE an implementation-level
violation into the definition**, where §3's proof can no longer be stated. Prerequisites are listed in
§9.7 O-5.

⚠⚠ **F2 (FATAL to "jointness buys the edges"): the gain is CONTIGUITY, not jointness, and not node
size.** Build `T|locus` = a **contiguous genomic window of `D|locus` whose length is EXACTLY
`len(R|locus)`** (asserted 27/27) — same locus, same length as the RNA node, **DNA content, no RNA
anywhere** — and run the **LITERAL shipped tier**:

| arm | tier | edges | of the 91 lost | of the 260 shared | symdiff vs `E_dna` |
|---|---|---|---|---|---|
| `E_dna` (stored) | shipped, literal | 351 | 91/91 | 260/260 | 0 |
| **ARM 1 `T × T`** (DNA-only, RNA-SIZED) | **shipped, literal** | **351** | **91/91** | **260/260** | **0** |
| JOINT `R × D` | translation `-P -N 50 -p 0`, *generous* | 344 | 88/91 | **256/260** | 7 |
| `E_rna` (stored) | shipped, literal | 260 | 0/91 | 260/260 | 91 |

Coverage on the 91: ARM 1 **0.8032 / 0.9913 / 1.0000** vs joint 0.4498 / 0.9817 / 1.0000 vs RNA×RNA
0.2928 / 0.4384 / 0.4989. **Node length contributes 0 of the 91** — per-node Spearman of median
T-coverage against the rna/dna length ratio is **+0.118** (n = 27), the most-shortened node (ratio
0.598, 16,428 → 9,828 bp) has median T coverage **1.0000**, and an 18% loss of total sequence
(542,101 → 445,487 bp) costs **zero** edges. ⭐ **CONVERSE: the joint object is not a superset of either
single-sided graph** — `R × D` **destroys 4 of the 260 edges `E_rna` already had**, all at identity
0.978–0.979 with coverage just under the floor and all touching the low-ratio nodes. Controls (both
exact): `D.fa` ava at the literal shipped tier → 351, symdiff 0; `R.fa` ava → 260, symdiff 0.
⟹ **§2.4's corrected unifying statement is confirmed and sharpened: an `E_r` edge needs at least one
endpoint CONTIGUOUS over the shared region — and DNA supplies contiguity for free.** ⚠ Self-criticism:
on a **complete** DNA graph ARM 1's 351/351 is close to entailed — but on that same
entailed-to-succeed panel the joint arm still **lost 7 edges**.

## 9.6 The falsifiers, answered

| falsifier | pre-declared prediction | measured | verdict |
|---|---|---|---|
| **F-J1** vacuity | 4 of 7 fail | **κ fires DISCORDANT 4 of 6 opportunities** [0.300, 0.903] | **HELD** — κ is not vacuous. It is also not validated (F-J6) |
| **F-J2** containment | APOBEC3 exists (DNA 0, RNA 1) | reproduced at the shipped tier on the 61-node set; **13 RNA-only edges** total | **HELD** — §1.2's "intersection = `E_RNA` identically" is **false off the 27-node panel** |
| **F-J3** harm of the expressed-member rule | deletes ≥ 21 | **7 of 80** panel nodes [0.043, 0.170]; 21 belongs to the **362-member Soto** denominator | **HELD IN SIGN**, magnitude is denominator-specific |
| **F-J4** O1/O2 | a block change would be fatal | **not run as declared** — but F1 finds the violation upstream, in `V`, on the **default** path, independent of `-F 2308` | **WORSE THAN PREDICTED. OPEN + FATAL.** |
| **F-J5** redundancy | the live question | **ANSWERED, YES**: on sparse families the RNA-lost edges change an outcome (HERC2 destroyed; NPIP@61 loses two members, permutation p ≈ 2e-4) ⟹ **§2.1's "it forbids nothing" is local to NPIP@27** | **FIRES** |
| **F-J6** κ's precision | must be scored against a truth it did not produce | **Fisher p = 0.40, n = 6, minimum attainable p = 0.10** | **UNTESTABLE at this size** |

## 9.7 What to adopt instead, and what is still open

Two things have been called "the joint object" and they have opposite answers.

| reading | what it is | verdict |
|---|---|---|
| **A. the SET fold** — `E_DNA ∪ E_RNA` or `∩` on the same nodes | a graph operation | **PROPERTY, NOT DEFINITION.** 0 blocks moved in O1; **0 of 7,691 reads** in O2; 0 in O3 |
| **B. the SUBSTRATE fold** — which extent supplies the BASES one `E_r` edge is computed on | a node-object decision | **ADOPT — as a spec + code fix, not as "jointness".** Already half-shipped at two call sites with **opposite defaults** (a gated family-local ADDITION vs an ungated global SWAP — not two spellings of one knob) |

⚠⚠ **THE SUBSTRATE-FOLD NUMBER IS RE-QUOTED (O-2, 2026-08-13).** The recorded *"+176 true / +13 false at
id ≥ 0.90, **precision 0.908 → 0.916 UP**"* was a **PANEL-tier** measurement. At the **shipped** tier with
the M1 coverage form, same 540 reps (chr1 325 + chr15 215) in ONE pooled index, **40,755 labelled pairs
pre-declared**: **+176 true / +14 false**, precision **0.9111 → 0.9164**. ⭐ **The true-edge gain SURVIVES
EXACTLY** — family-clustered bootstrap **[+67, +312], P(Δ>0) = 1.0000**. ⚠⚠ **The precision claim does
NOT**: Δ **+0.0053**, CI **[−0.0216, +0.0412]**, **P(Δ>0) = 0.628** — a coin flip that rests on 2–3 false
edges out of ~360, and whose sign flips on the **coverage form**, not the tier. At id ≥ 0.95 and 0.98
precision goes DOWN in all four cells. ⚠ **And the mechanism attribution is REFUTED**: of 195 gained true
edges, **144 (74%) are one-single-exon** — the "stub" class the 2026-08-03 retraction declared untouched —
51 both-spliced, 0 both-single-exon. The median table that produced the original attribution hid its own
tail. ⚠ 2 families supply 52% of the gain (the net survives family-clustered resampling anyway). ⚠ The
tier is nearly irrelevant here: panel → shipped moves **one** exon-sum edge and **zero** span edges.

**Downstream, measured.** O2 is **byte-identically unchanged**: `K` enters at exactly one arithmetic
place (`thr = alpha/max(K−1,1)`, `copy_assign.rs:451-460`); replaying five shipped gorilla panels with
evidence fixed, **7,690/7,691 reads reproduce the shipped TIED label**, **K ± 1 moves 0 of 7,691
(0.000%)**, K ± 3 moves 18 (0.234%) — because `min_p` is **bimodal** (977 reads at `p = 1.0`, 6,506 at
`p ≤ 1e-20`) and the band `1e-4 < p ≤ 1e-3` that the whole range `K = 2..11` can move `thr` across
contains **0 reads**. And the "member must be expressed" fold has nothing to prune: **no shipped
catalog contains a zero-read copy** (min `n_reads` = 2 across 1,415 / 1,220 / 245 / 212 copies, DNA
path included) and `copy_assign --families` **aborts** on one. For **O3** the joint object is
**necessary for the objective to be stateable at all** — *a missing copy is a property of the
(reads × genome) product* — and buys **no measured power**; ⚠ and the supporting "reference-free
clustering is bit-identical intact vs degraded" is a **design tautology** (the read-overlap statistic
never opens the genome), while "26/26" must be quoted as **24/26 strictly negative, 26/26 within
|Δ| < 0.002**, whose real content is that **removing a copy removes the competitor that produced the
secondary alignments** — the obvious DNA proxy is **anti-correlated** with absence.

**OPEN, costed.** O-1 re-derive the RNA leg from the full BAM at the shipped tier (~4–6 h; bounded —
RNA-only cannot exceed 13, so it **cannot reverse the verdict**, but "84.5%" stays INHERITED until it
runs) · ⭐⭐⭐ **O-2 DONE 2026-08-13 — the folds are now measured on the object the binary emits**
(`/home/juanfra/winloci_scratch/o2_shipped_union/O2_SHIPPED_UNION.md`; summary in §9.0 above).
**VERDICT: STANDS, RESTATED.** The shipped `E_r` *is* `E_x ∪ E_g` at the `refine` site (counter-verified:
101 added edges = exactly `|E_dna \ E_rna|` at 61 nodes); it fires on **11/26** families under
`--cross-chrom` and is **sole support for one emitted O2 family carrying 339 read assignments**
(SHARP `CAFAM0`: exon-sum identity 1.0000 at coverage **0.0480** → 0 edges; genomic identity 0.9811 at
coverage **1.0000** → 1 edge) — yet it lands on the DNA-only partition **0/26** [0.0000, 0.1287] with
emitted catalogs identical (ARI 1.0000). ⚠ **The union is a no-op relative to DNA-only and never
relative to the shipped core**, the published **0/5, 0/7** is an entailment (**2/5, 2/7** under the
shipped gated-union rule at the rep substrate), **`E_RNA ⊆ E_DNA` is false past 27 nodes** (66 of 696 at
244 reps), and the shipped O1 substrate is the **exon-sum**, not `seq_g`. ⚠ **Still panel/region scope
only** — chr16 abandoned at 35 min, chr1's 96 GB BAM does not subset; **n = 7** on the gorilla
`copy_assign` path, so `1/7` is not a rate · ⭐ **O-3 DONE 2026-08-13** —
`er_rule_rows` now emits `core_substrate` (renamed from the misleading `substrate`) and a RULE row
`additive_genomic_tier` whose value states the armed state **and the reason when off**; the refine
`params.tsv` gained the DATA rows `n_families_core_unconnected` / `n_families_genomic_tier_ran` /
**`n_edges_genomic_tier_added`** / `additive_genomic_tier_fired`; the arming decision is taken once by
`additive_genomic_tier(params, genome.is_some())` and the gate and the certificate read that same value
(source-level test forbids a second derivation); and the core vs additive `.paf`/`.args` — previously
byte-identical but for a counter — are now tagged `refine.core.exon-sum` / `refine.additive.genomic-span`
with a `tier` line in `.args`. ⚠⚠ **Making it truthful FALSIFIED the X.2 "EMPTY DIFF" certificate:** at
shipped defaults with a genome reachable, `diff O1.rule.tsv O2.refine.rule.tsv` is now **one line** —
`additive_genomic_tier`: `absent (single-substrate site)` vs `armed (genomic-span …)`. The two sites agree
on every other edge-deciding knob and **disagree on the substrate**, exactly as F6 said; the old empty
diff certified nothing about that axis. ⚠ **the default was NOT changed** — observability only ·
⭐⭐⭐ **O-4 DONE 2026-08-13 — BOTH DEFAULTS KEPT; full decision, denominators and CIs in
`bench/FALSE_NEGATIVES.md` (restored).** The two sites are **not the same operation**: refine's leg is a
**gated, family-local ADDITION** (fires only where the exon-sum core leaves a copyset disconnected, so it
can never bridge two families), the catalog's is an **ungated GLOBAL SWAP** (+879 edges, +126%). Refine
keeps it **ON** — on the 26 examined families it moved **11/26** partitions yet block sets vs DNA-only
differ **0/26 [0.0000, 0.1287]** with emitted catalogs **identical (ARI 1.0000, 0 forbidden pairs)**, and
it is sole support for SHARP `CAFAM0` (exon-sum coverage **0.0480** → 0 edges; genomic **1.0000** → 1)
while sibling `fam1` was correctly rejected at 0.1316. The catalog keeps `homology_genomic_span` **OFF**
on **unresolved sign, not measured harm**: on U = **18,528** annotation-labelled pairs over the same 244
reps, recall **0.1490 → 0.1780** (P(Δ>0) = 0.9673) but precision **0.8501 → 0.9876** does **not** survive
family-clustered resampling (40.7% up / 19.9% tied / **39.4% down**), and ⚠⚠ **the same contamination
formula on the same two partitions flips sign with the truth** (+0.0310 worse under RefSeq-coarse,
**−0.0377 better** under the panel roster). ⟹ **§1's "`E_r` is computed on `seq_g` only" and §1.2's
"edges on `seq_x` — Rejected, strictly weaker" remain CONTRADICTED BY THE SHIPPED O1 DEFAULT**; O-4
records the contradiction rather than resolving it in either direction, because the axis that would
decide it has no stable sign. ⚠ **Item 4 of §9.4 above is now half-wrong and is corrected here:**
`bench/FALSE_NEGATIVES.md` was **deleted in `9b0814f`**, not never-created — it is restored verbatim from
the git blob, and a test now asserts the citation resolves · O-5 settle O1 ⊥ O2 at the node — default
`RUSTLE_FLAGFREE_SITES=1`, or remove `reads_distinguish`/`distinguishing_uniq` from
`distinct_locus_reps`, or retire the claim; and extend the guard to ban `locus_unique_mapper_counts(`,
`reads_distinguish(`, `distinct_locus_reps(` and assert `locus_min_reads ≠ conflict.min_reads`
(~1–2 days) · O-6 run `strata/attrib/run2.sh` (arm 2 `T × D`, placement, length ladder; ~1 h, nothing
depends on it) · **O-7 repeat the redundancy and attribution arms on a SPARSE DNA graph** (~1 day; one
informative null cell is not a result) · O-8 score `κ` against an independent truth on ≥ 20 families
(~2–3 days) · O-9 look directly at **AC124944.2** (family size 2, genome-wide best paralog 0.921, 0
paralogs ≥ 95%), the one residual not cleanly explained as a leg artifact (~1 h) · **O-10 retract the
published complementarity line in memory** (`project_dna_rna_complementarity.md` and the `MEMORY.md`
index line still carry *"RNA is NOT a subset of DNA"*, now directly contradicted; ~15 min) ·
**O-11 measure the additive tier's firing rate at GENOME SCALE** (~1 day on a quiet machine; every O-2
number is panel/region scope — chr16 was abandoned at 35 min in the all-vs-all and chr1's 96 GB BAM does
not subset) · **O-12 does the union change an assignment INSIDE an already-existing family?** (~2 h:
`copy_assign` default vs `--no-refine`, read-level diff on the gorilla panel; O-2 answered only the
existential form — 339 assignments exist *because* of the tier — and the project's own rule forbids
judging a change to what a NODE IS on node-level metrics) · **O-13 a real node-set ladder for
`E_RNA ⊆ E_DNA`** (~1 day; 0/351 → 9/316 → 66/696 is three points, not a curve).

## 9.8 Additional exposures found in phase 3 (append to §7)

9. **CONTROL B, unresolved and independent of this task.** An RNA node aligned to **its own** DNA node
   is definitionally a spliced subsequence and must score coverage 1.0000 at 27/27. It scores median
   **0.9344**, passes only **25/27**, minimum 0.4498. Every RNA-vs-genomic claim in the thesis sits on
   that floor.
10. **The CF-B number is a loaded gun.** *"91/91 pass if the single-record clause is lifted, median
    0.9205, from a median of 79 identity-passing records per pair"* is a **diagnosis, not a proposal**.
    `RUSTLE_ER_SUM_COVERAGE` is default OFF and condemned in `bench/soto/rustlib.py` with a measured
    cross-family precision cost. §2.4 already warns; the warning is repeated here because this is the
    most likely thing to be misquoted out of the run.
11. **Mixed edges are tier-sensitive across a 19–88 range.** `mixed_X.paf` is **empty** at the shipped
    tier; the 88/91 uses `-P -N 50 -p 0` (269 vs 260 on the control, symdiff 11 — *more permissive*),
    and the default-secondary translation gives **19/91**. **88/91 is not a shipped-tier number.**
12. **Two of this run's own intermediate claims were refuted by its own measurement** and are on disk
    with banners (`strata/mixed/PREDICTION_CORRECTION.md`, `PREDICTION_UPDATE.md`,
    `strata/redundancy/ONE_ENDPOINT.md`). The hashed pre-registration (predicted 70/91, interval
    55–88) was **correct at 88/91**; the mid-course revision against it was **wrong**.
