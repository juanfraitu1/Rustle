# The advisor's standing questions — and what the codebase can answer

> **Audience.** Stefan Canzar, who assumes a good number is luck or overfitting until shown
> otherwise. This document is written for that reading. Every claim carries the number that
> earns it and the section that derives it; every claim we **cannot** defend is stated as a
> concession **before** he finds it.
>
> **Posture.** Do not open with results. Open with the **discipline**, because the discipline is
> the only thing that makes the results admissible to someone who starts from disbelief. The
> order below is deliberate: Part 1 is the cross-examination he will actually run; Part 2 is the
> per-question evidence; Part 3 is what we concede unprompted; Part 4 is what to put on screen.

**Provenance.** Derivations in [`o1_ledger.md`](o1_ledger.md) (120 sections), negatives in
[`NEGATIVE_RESULTS_REGISTER.md`](NEGATIVE_RESULTS_REGISTER.md) (**836 rows**), status in
[`OBJECTIVES_AND_VERIFICATION.md`](OBJECTIVES_AND_VERIFICATION.md). Test baseline **824 passed /
0 failed / 11 ignored**.

---

## Part 0 — The object (read first)

**`docs/O1_O2_COMPOSITION.md`** states the method as a composition: **O1's vertex set IS O2's
path set**. Two graphs at two granularities — never "one decision rule", never "a family is one
variation graph" (`NEGATIVE_RESULTS_REGISTER.md:472`, `:1086`, `:1090` kill both, and `:1090` is
annotated as exactly what would irritate him). It also carries the corrections this file needs:
⚠ **§1.2's "four free numbers" is wrong** — ~25 constants are default-reachable, +6 added
2026-09-03 in `mcl_families.rs`. ⚠ **The excision abstention result (§Q-abstention) is not a run
of the shipped gate** — it is the robust-z of `de`, not the α-certificate.

---

## Part 1 — The cross-examination: "this is luck or overfitting"

These are not hypothetical. They are the five moves that follow from his stated priors, and each
has an answer that is already in the tree.

### 1.1 "You tried hundreds of things. The winners are the tail of a null."

**This is the strongest attack and it must be answered first, with the accounting, not with a
denial.**

| | count |
|---|---:|
| Routes attempted and killed, each with the number that killed it | **836** |
| Ledger sections (each an attempt, an audit, or a retraction) | **120** |
| `RUSTLE_*` behaviour flags in `src/` | **135** |
| …of which the shipped default path turns **ON** | **7** |
| Defaults flipped in the last month of work | **1** (`NODE_MIN_READS` 3 → 2, §6ac) |

⭐ **The accounting is the answer.** ~950 hypotheses were tested; **one** default changed. If the
survivors were a chance tail we would have shipped dozens of them — the selection pressure that
produces overfitting is *adopting* winners, and the adoption rate here is under 0.2%.

⭐⭐ **And the audit that looked specifically for winners to adopt found none.** §6aa classified
every boolean behaviour flag on the O1 path by measured verdict and concluded ***"every measured
flag is negative or break-even; the rest are unmeasured. The shipped defaults are already correct
and NOTHING should be flipped."*** A pipeline tuned to its benchmark does not produce that
sentence.

⚠ **What he can still say.** The register records *outcomes*, not a pre-registered analysis plan.
Many hypotheses were formed after seeing data. **Concede this** — and point at the mitigation
below, which is real: **35 ledger passages record a criterion fixed before the result**, including
the read-strand run (§4o, *"every pre-registered criterion passes"* — and it was then **still not
shipped**, §4p) and the engulfment test (§6ar, pre-registered, returned **partially supported**).
The instrument that most often failed pre-registration is the one we most often refused to ship.

⭐⭐⭐ **THE CLEANEST INSTANCE, AND THE ONE TO SHOW HIM (09-02, §6bq → §6bt.2).** The two-sided
coverage clause was developed **entirely on gorilla NPIP**. The 150-window negative panel is
**human** CHM13/A119b, frozen since 2026-08-10, never used to tune anything. The prediction was
written, **committed to git as `5cbced4`, and pushed BEFORE the arm ran** — and it was
*mechanistic*, not a count: which window would die, why, and that the other **must survive**.

| | |
|---|---|
| pre-registered criteria | **5/5 pass**; none of the three falsifiers fired |
| edge outcomes called by `cov_longer < 0.30`, across two node floors | **9/9** |
| edges or windows **added** in any arm | **0** (monotonicity held) |
| the case it predicted it **could not** fix | W063 — and it did not fix it |

⚠ **Say the caveat in the same breath, because it was also fixed in advance:** the rate movements
are **2–3 events** and the intervals overlap almost entirely. **Do not say the rate halved.** The
claim is that a rule fixed before the run called **9/9 edge outcomes across a species boundary**.

### 1.1b "You pre-registered something you already knew would pass."

The honest reply is that **the thing which passed was then refused a default**, on criteria also
fixed in advance (§6bu → §6bv).

`RUSTLE_ER_COVERAGE_LONGER_FLOOR` passed the cross-species panel **5/5, 9/9**. It is still **OFF**,
because the adjudication found:

| criterion | outcome |
|---|---|
| D1 recall (NPIP 31) | met — ⚠ but **the entire gain is one 274 bp, 2-read, single-exon copy** |
| D2 specificity | met, strongly |
| D3 losses in defective strata | consistent, **declared weak in advance** (correlated with the selector) |
| **D4** no penalty on corroborated copies | ⛔ **NULL** — the copy loss is **annotation-neutral** |
| **D5** trims rather than deletes families | ⛔ **fails** — **18 families deleted outright** |

⭐⭐ **The sentence that does the work:** *at the **edge** level the clause is sharply discriminating
(9/9 across a species boundary); at the **copy** level it is **indiscriminate** — within exon strata,
deleted and retained copies carry a reciprocal RefSeq match at the same rate (0.013 vs 0.009
single-exon; 0.791 vs 0.784 multi-exon). It removes 27.7% of copies without enriching for the bad
ones.* A project optimising its own numbers ships that flag; this one did not.

⚠ **And the pre-registration itself failed in a way worth volunteering.** §6bu named **length** as
D4's confound. Within length quartiles the clause looked like it preferentially deleted
uncorroborated copies — in all four strata. The operative confound was **exon structure**, and
controlling for that erased the effect. **Naming a confound in advance does not protect against
naming the wrong one**, and that is now a register row (640), not a footnote.

### 1.2 "You tuned the thresholds until the answer appeared."

⭐ **2026-09-05 (§6ez, `docs/g3_threshold_grid_2026-09-05.tsv`):** on the thesis definition the three edge
thresholds (identity 0.70, coverage 0.30, 300 bp) were swept one at a time on two substrates (gorilla anchors,
Soto slice): every anchored family and every Soto score is unchanged over identity 0.60–0.80, coverage 0.10–0.50
and 100–500 bp; the walls are at identity 0.85 (Soto recall 0.949 → 0.925), coverage 0.60 and 1 kb. Inflation
(2.0–4.0, §6ec), prune (size-safe, §6ec), the exonic floor (the wall at 1 bp, §6dt), "half" (one constant) and
the 3-read support (the corroboration floor) were justified before. Say: *"every default is a point inside a
measured plateau whose walls I can show you."*

The shipped rule has **four** free numbers, and they are visible in the source, not in a config
that drifted:

| parameter | value | `src/` |
|---|---|---|
| edge identity floor (sensitive tier) | **0.60** | `denovo_pipeline.rs` |
| edge identity floor (asm20 tier) | **0.80** | `denovo_pipeline.rs:3765` |
| edge coverage floor, **of the shorter** | **0.50** | `denovo_pipeline.rs:3766` |
| quasi-clique density γ | **0.20** | `family_definition.rs:173` |

⭐ **γ is not sitting on a tuned optimum — it is sitting on a *measured* one, and the sweep was run
by an arm that was trying to beat it.** §6bg swept γ upward against the seeded catalog looking for
a precision lever; **F1 peaks at the shipped γ = 0.20**, and raising it is recorded as ⛔ *"wrong
lever."* We went looking for a better value on a different catalog and the shipped value won.

⭐ **The coverage floor was attacked twice and survived both times, and the second attack produced a
retraction of our own recommendation** (§6g tested it, §6h **retracts** §6g's recommendation
because the register had already refuted it). §5c lowered it on edge evidence and it **failed end
to end**.

⚠ **Where he draws blood — and he is right.** The 0.50 floor is **one-sided**: it charges coverage
only on the shorter sequence, so *"a 10% fragment that aligns fully into a complete sibling scores
1.00"* (a comment the source has carried since long before this was measured). The named failure is
concrete: a **2,037 bp NPIPB6 fragment reaches coverage 0.948** against a **38,653 bp** chimeric
read-through node while touching **5%** of it, dragging EIF3CL into NPIP.

✅ **This is now fixed as an opt-in, validated end to end, and then formally adjudicated** —
`RUSTLE_ER_COVERAGE_LONGER_FLOOR` (§6bp, §6bt.2, §6bv): the OFF arm is **byte-identical** to the
prior catalog, the params certificate distinguishes the arms, the human panel passes **5/5 / 9/9**,
and NPIP recall goes **14/31 → 15/31** ⚠ *on one 274 bp, 2-read, single-exon copy*.
**It remains default OFF** — see §1.1b. The 27.7% copy loss is **annotation-neutral** and cannot be
priced without a positive stratum; NPIP labels **21 of 678** copies.

### 1.3 "Everything you have is one family."

**Concede immediately — the number is worse than he will guess.**

**66 of 120 ledger sections mention NPIP or its 31-locus panel.** Every O1 decision between
2026-08-25 and 2026-09-01 was scored on that one panel, which is itself a **minimap2 projection of
human NPIP onto the gorilla assembly**. The clean control has **n = 3**.

**What survives the objection anyway** — results measured where NPIP cannot reach:

| evidence | number | § |
|---|---|---|
| **Cross-substrate replication** (different animal **and** different tissue) | **87.06%** of edges reproduce; clean corner **130/136 = 95.6%**, marginal **390/479 = 81.4%** | §4l |
| **One-seed closure on HUMAN families** | **65/65 converge** | §5p |
| **Haplotype CNV proven inside one animal** | direct proof of the phenomenon, not a projection | §6u |
| **False-merge rate** on gene-tight windows with demonstrated power | **2/150 = 1.33%** [0.37, 4.73] — **reproduced exactly** by the current binary at `RUSTLE_GATE_MIN_READS=3`, same two windows, same 3 edges. At the shipped floor 2: **3/150 = 2.00%**, **disjoint** set. One parameter, sets swap completely | O1.10 · §6bt.1 |
| ⭐**Self-overlap defect, cross-species** — a standing blind spot, measured 09-02 | present at **7.09%** (GGO) and **7.10%** (PTR) in the `refine`-built catalogs; **0/4,176 pairs [0, 0.09%]** on the current path | §6bs |
| Two-sided coverage gain reproduced on a **second substrate** | holds — ⚠ but **the mechanism did not** (§6bl retracts the shared-domain explanation) | §6bl |
| ⭐⭐⭐**Pre-registered CROSS-SPECIES test** of a gorilla-derived clause on a human panel | **5/5 criteria, 9/9 edge outcomes**, prediction committed before the run | §6bq · §6bt.2 |

⭐ **The strongest single item is the cross-substrate replication**, because the relation was never
tuned on that animal or that tissue. ⚠ **But state its weakness in the same breath**: **not one base
of read sequence enters `E_r`** — every base comes from `genome.fetch_sequence` — so segmental-
duplication corroboration **shares the substrate**. It is corroboration, **not independence**.

### 1.4 "Your validation is circular."

This is where the project is genuinely strong, because **we killed our own metrics repeatedly and
recorded each kill**. Present it as a list of self-inflicted retractions — it is far more
persuasive than any surviving number:

- ratio-to-truth ⟹ must read the **in-band fraction**, never the median
- **"bases explained"** — banned: it rewards unspliced models
- **prediction ⊆ its own truth** — tautological; **3 metrics killed**
- **a denominator conditioned on the prediction** — **7 metrics killed**
- **selecting which component to score** — killed a purity of 0.237
- an **edge-count-matched null proves nothing**; the size distribution must be matched
- never judge a change to *what a node is* on node-level metrics — **3 failed end-to-end**
- **§6bb's own counts retracted** because the script was not preserved and never recorded its
  locus-assignment rule
- ⭐**09-02, §6br** — an annotation-based corroboration of the weak identity band reached **39.9×**
  a size-matched null and **survived a proximity control**, then died when we checked *what* was
  agreeing: **776 of 788 agreements (98.5%) were the single string `"zinc finger protein"`**, which
  the annotation carries on **552 genes**. Two instrument defects were found on the way — a
  large-gene attractor putting `"titin"` on 26 endpoints, and a size-matched null that a *spatial*
  predictor walks straight through (the first headline was pure co-location: cross-chromosome edges
  in the target band scored **0.36×, below chance**)

⭐⭐ **The one non-circular O2 result, and it is the one to lead with.** In the excision run copy A
is deleted and its reads migrate to locus B, so their true origin is A **by construction** — labels
no aligner produced.

| stratum | TPR (foreign) | FPR (native) | AUC |
|---|---|---|---|
| all loci | 0.2404 | 0.0239 | 0.6918 |
| **< 50% migrants** | **0.5066** [0.4957, 0.5176] | **0.0280** | **0.7995** |
| **MAPQ, same task (control)** | — | — | **0.4944 — chance** |

⭐ **minimap2 is not merely wrong here, it is *confidently* wrong**: median MAPQ **60 vs 60**, and
MAPQ = 60 covers **96.07% of migrants vs 94.98% of natives**. Its confidence carries **zero**
information about whether the read belongs.

### 1.5 "You would not show me the failures."

Hand him [`NEGATIVE_RESULTS_REGISTER.md`](NEGATIVE_RESULTS_REGISTER.md) — **836 rows, each with
the number that killed it**, and the two admitted exception classes (**NO-POWER**, and killed-by-
argument) marked as such rather than hidden. Then hand him the ledger's own index note: an earlier
auto-derived verdict tag scored **11/22 = 50% — a coin flip — against sections whose outcome was
known first-hand, so the tags were removed rather than shipped.**

⭐ **That is the single most disarming artifact in the repository.** It is a record of a convenience
feature deleted because it was not reliable enough to be honest.

---

## Part 2 — The standing questions, with the evidence

### Q1. "Do you have a method that identifies multi-copy gene families?"

**Yes, with a stated scope.** A family is a connected component (≥2 loci) of a homology graph
`E_r`, refined by a γ-quasi-clique partition. The definition is **seed-free**, **P1 (domain-sharer
exclusion) is a theorem**, and it carries a falsifiable false-merge rate of **1.33%** [0.37, 4.73] —
**reproducible today** at the node floor that produced it, and **2.00%** [0.68, 5.71] at the shipped
floor, on a disjoint window set (§6bt.1). Quote the floor with the rate.

⚠ **The binding constraint is named and not solved: node construction, not the definition.** The
oracle ablation (§5e) prices it — **the definition costs ~3%, node construction costs ~58%**. Say
this before he finds it; it converts a hidden weakness into evidence of self-knowledge.

### Q2. "Are these real families, or artifacts and overfitting?"

See Part 1. The summary line: **defensible narrowly, not broadly.**

⚠⚠ **The concession that matters most — the evidence covers the wrong end.** Median edge identity
is **0.8287** and **86.31% of edges sit below 0.90**, but *all* external support is at **≥ 0.90**
(NPIP annotated median 0.9779, GOLGA6L7 0.9673). **Nothing external covers the ~0.83 band that is
most of the catalog.** ⛔ **09-02: we tried to close this with the gorilla-native annotation and
failed** (§6br). The failure is structural, not a missing effort: **61.4% of catalog nodes have no
reciprocal gene match at all**, and among those that do the agreement is **domain-level** — strip
zinc-finger stems and the target band scores **0/8**. ⚠ Say what this does and does not mean: it
shows **the annotation cannot adjudicate these edges**, not that the edges are wrong. Also: **88 of 121 families (72.7%) have no segmental-duplication containment
at any floor**, and only **5** have per-family external adjudication. **Zero experimental
validation** — no ddPCR, qPCR, or FISH.

### Q3. "Does the method borrow information across the family?"

**No — and the honest answer is better than a hedge.** Cross-copy borrowing was implemented and
measured **inert or dead** (§6bd). The information that *is* shared is structural: multimapping
reads are treated as **shared evidence** rather than a conflict to be resolved — the deliberate
inversion of his 2016 framing, and the thesis's actual position.

⭐ **One orthogonal lever does exploit shared reads and it works**: `E_c` (shared-multimapper
edges) **splits** large families that `E_r` fuses — **0.2406 → 0.5431** on top of the direct-edge
rule (§6bi). ⚠ Scoped: it reaches only near-identical arrays, and it is **not** a definition tier
(§6ae: no depth threshold exists).

### Q4. "What about isoforms? Could two extremely similar copies produce the same isoforms?"

Isoforms are **detected** (174 supported chains at NPIP) and then **collapsed to one representative
per locus** before family definition — deliberate, because for O1 a family is a set of **loci**, not
of transcripts. ⚠ **A representative is one intron chain, not the isoform union**, and **26.8% of
shipped copies are single-exon stubs**. Both are disclosed, both are live defects.

### Q5. "Two tandem near-identical copies share a read. What happens?"

**Measured.** **152 primaries** touch ≥2 copies of the same family = **0.306%** of the cases where
it is geometrically possible (**89/678 copies, 14/121 families**). The cause points **away** from
biology: the copy-joining intron is annotated **0.182** against a length-matched **0.818** — a
**4.5× depletion**.

⛔ **But half is undetermined and must be conceded: 73/153 = 47.7% are canonical-but-unannotated =
"we do not know."** Defensible statement: **≥31% minimap2 artifact, ~18% real, nothing about the
plurality.**
⚠ Trap: §6aq's discriminator is **degenerate below ~3 bridging junctions** (one intron ⟹ modal
share 1.00 *by construction*); a naive tally gave a garbage 39/59 = 66%. **Check the junction count
first.**

### Q6. "Does it port to other families, tissues, apes?"

**Family/tissue: yes** (§4l, above). **Apes: not today.** The ape BAMs are drop-in — each aligned
to its own reference, identical minimap2 line, `-N 50`, indexed — ⛔ **but all four catalogs were
built with `refine`, a default since removed, so the current binary reproduces none of them.**
⚠ **Do not quote "149 ancestral + 84 expansions"** — that is a superseded 3-species run overwritten
84 minutes later; the repo file is the 4-way split.

### Q7. "Aren't the TSS/TES/UTR boundaries too convenient?"

**The premise is false and that is measurable** (§6ay) — ⚠ but the coverage statistic used to
answer it **was one-sided** until `cov_longer` was emitted (§6ba). Answer the question, then
volunteer the flaw in the instrument that answered it.

### Q8. "1/k is a bad assignment for tied multimappers."

**Agreed, and never used.** O2 is **assign-or-abstain**; the objective **provably decomposes**, so
the shipped per-read gate **is** the optimum and no joint estimator can beat it. **K = 0 abstention
is entailed, not chosen.** Empirically the EM changes **zero** of the gate's decisions on reads
that carry evidence.

⚠⚠ **Restate O2's target before he does.** "MAPQ-0 ambiguity" names **0.04%** of reads inside the
multi-copy loci — **rarer there than genome-wide**. The real contested population is **alignment-
score near-ties: 21.75%**, ~**500×** larger. **Do not claim O2 assigns reads better than minimap2 —
on ~99.9% of reads it measurably does not.** Claim **abstention**, which is where the excision
result (Part 1.4) lives.

### Q9. "NPIPA and NPIPB should be distinct subfamilies."

⭐ **2026-09-05, on the thesis definition (SD-core family, `rna_units_v3` MCL3, §6ew):** the 29 gorilla NPIP loci
land on CHM13 as **17 NPIPB-only, 3 NPIPA2+NPIPB13, 9 ABCC1/SORL1 chimeras** (5 trimmed to the LCR16a core, 4
dropped) — none on NPIPA alone. Identity (0.9015 vs 0.9007) and exonic coverage (0.588 vs 0.533) do not separate
the A-landing from the B-landing loci, and MCL keeps all 20 in one cluster from inflation 1.4 to 4.0; its first
cut is exactly the four core-0 chimeras. Say: *"in gorilla the A/B split is not a cut; the catalog reports the
CHM13 landing per locus, and the family's first partition coincides with the core rule's drop arm."* ⚠ The
coverage-splits-A-from-B answer below is from the EARLIER E_r graph and does not transfer (row 695).

**He was right, and the method already agreed.** ⭐ **Identity cannot separate them; coverage can** —
A↔A median coverage **0.46**, A↔B **0.12**, B↔B **0.06**, while identity is **0.99 / 0.96 / 0.99**.
On identity NPIP is one clique; the **coverage floor** is what splits it. ⚠ Superseded detail: the
dominant failure is now **fragmentation** (14 members → 13 families), not contamination — and
**§5j** finds NPIP is **one family fragmented into 5–6**, not three real subfamilies.

⭐ **Our measured precision on NPIP is understated**, because Soto's set is CAT-bounded: a real copy
CAT missed scores as a false positive. The defensible exhibit is **chr16:28,659,994 — 21 exons,
1,327 reads, identity 1.000 over 24 kb, 21/21-exon match to NPIPB9**, unannotated.

### Q10 / Q11. Non-canonical junctions; PSV credibility.

**He is right on non-canonical junctions** (§6au), the sites **recur across substrates** so they are
real (§6av), and `RUSTLE_JUNCTION_MAJORITY` is measured, works, and is **not yet a default** (§6aw).
Both PSV objections are answered (§6aj) — **and do not build the VCF.**

---

## Part 3 — What we concede before he asks

State these unprompted. With an examiner who assumes overfitting, volunteering the ceiling is the
only move that buys credibility for what is below it.

1. **No experimental validation.** No ddPCR, qPCR, or FISH. Every number is computational.
2. **The external corroboration covers the wrong end** of the identity distribution (Q2).
3. **72.7% of families have no external adjudication of any kind**; 5 of 121 have per-family review.
4. **The project is NPIP-bound** — 66/120 sections, clean control n = 3.
5. **Half the tandem-read cases are undetermined** (Q5).
6. **No ape catalog is reproducible by the current binary** (Q6).
7. **O1's binding constraint is node construction (~58%), not the definition (~3%)** — the part we
   formalised is not the part that limits the result.
8. **O2 does not beat minimap2 at read assignment** and we do not claim it does (Q8).
9. **O3 is a bounded negative, not a null result**: 0 collapse-shaped deficits at 816 evaluable loci
   with a 0/817 [0.00, 0.47]% false-positive floor, against a literature prior of **0.47–0.94**
   expected — ⟹ **a perfect screen returns zero 39–63% of the time.** The instrument cannot
   currently distinguish "nothing there" from "we cannot see it," and that is the honest claim.

---

## Part 4 — What to put on screen, in this order

1. **The register** (635 killed routes) and the deleted verdict-tag note — establishes the
   discipline before any result is shown.
2. **The forking-paths accounting** (Part 1.1): ~950 tested, **1** default changed, and §6aa's
   audit finding nothing to flip.
3. **The pre-registered human panel** (§6bt.2) — a gorilla-derived clause, a human substrate it
   never saw, the prediction pushed to git before the run, **5/5 and 9/9**, and the one case it
   said it could not fix left unfixed. Show the commit timestamp.
4. **The excision O2 result** with the **MAPQ AUC 0.4944 control** — the only fully non-circular
   accuracy number in the project, and its control is what makes it one.
5. **Cross-substrate replication** 87.06% / 95.6% clean corner — with the shared-substrate caveat
   said aloud.
6. **The NPIP coverage table** (Q9) — it shows the method *disagreeing usefully* with a naive
   identity reading, on the family he raised himself.

7. **[`REPRODUCE.md`](REPRODUCE.md)** — one catalog pinned end to end: source SHA, binary md5,
   both input md5s, the command, and the md5s of the outputs, with the OFF-arm byte-identity that
   proves the current binary still emits it. Have it open when he asks whether the numbers move.

⚠ **Do not lead with catalog sizes.** No pre-08-30 catalog count is reproducible by the current
binary (`NODE_MIN_READS` 3 → 2), and a number he cannot reproduce is a number he will assume was
chosen. `REPRODUCE.md` pins the **one** that is — three contigs, ~40 minutes, one command.

---

## Part 5 — Where each question is answered

| question | verdict | evidence |
|---|---|---|
| Q1 method exists | ⭐ defensible, scope stated | §5e, O1.10 |
| Q2 real vs overfit | ⚠ **narrowly** defensible | §4l, §6ax, Part 1 |
| Q3 borrowing | ⛔ no (inert) / ⭐ `E_c` splits | §6bd, §6bi |
| Q4 isoforms | ⚠ by design, two live defects | §6be, stubs |
| Q5 tandem reads | ⚖️ half measured, half unknown | §6aq, §6ar |
| Q6 portability | ⭐ tissue/animal · ⛔ apes | §4l, §5p |
| Q7 boundaries | ⭐ premise false | §6ay, §6ba |
| Q8 1/k | ⭐ never used; target restated | P-2, blind-spot audit |
| Q9 NPIP subfamilies | ⭐ answered on his own example | §5j, fam72 review |
| Q10 non-canonical | ⭐ he is right; recurs | §6au, §6av, §6aw |
| Q11 PSV | ⭐ answered | §6aj |
