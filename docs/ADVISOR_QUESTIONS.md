# The advisor's standing questions — and what the codebase can answer

> **2026-09-06 pass.** Q1, Q2, Q4, Q5, Q8 and Part 3 items 7–9 rewritten to the SD-core definition and the read-star O2
> (§6ev–§6fo). Paragraphs marked as the OLD node's measurements are kept as history. Q6 (apes) is unchanged and
> still unanswered.

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

## Part 0b — What Soto 2025 is, and what it is not (READ BEFORE QUOTING ANY SOTO NUMBER)

⚠⚠ **Soto's catalogue is a 98 %-identity SEGMENTAL-DUPLICATION catalogue, not a gene-family catalogue**
(§6fz, register 741, measured 2026-09-07 against CHM13 RefSeq biotypes):

| what its 362 members overlap | n | share |
|---|---|---|
| protein_coding | 140 | **38.7 %** |
| transcribed_pseudogene | 110 | 30.4 % |
| pseudogene | 58 | 16.0 % |
| lncRNA | 39 | 10.8 % |
| no annotated gene at all | 11 | 3.0 % |

**Fewer than two in five members are protein-coding genes; 46.4 % are pseudogenes. Of the 76 multi-member
families, 60 MIX biotypes.** Within-family length spread: median **3.1×**, 24 families above 5×, 12 above 20×,
extreme `ID_14` at **583.7×** (0.1 kb to 65.4 kb) — one "family" holding a 100-bp fragment and a 65-kb gene.

⟹ **What this means for the thesis.** A recall number against Soto is partly a count of how many pseudogene
fragments a method admits. A method that defines a family as a set of *genes* will lose members there **by
construction**, and that loss is not an error. So Soto is the wrong instrument for the question "did you find
the family", and the right one for two narrower questions.

**Soto remains the best available instrument for:**
- **Paired comparisons** — both arms pay the identical price, which is why `O1_DEFINITION_SWITCH.md` §2 is a
  valid head-to-head even though its absolute levels are not interpretable as gene-family recall.
- **Precision in the adjudicable [0.90, 1) band** — a pair it asserts is a real duplication.
- **Independence** — it is not our instrument, and it is CHM13-based, so it corroborates across a substrate
  boundary.

**Soto is NOT a ceiling, a gold standard, or the final word on whether a family was found.** ⚠ Its detection
levels are also slice-conditioned: only its own neighbourhoods exist in the comparison, so genome-wide false
merges cannot occur in either arm. ⭐ Our own gorilla truth carries the same defect in miniature and it is
declared: of the 26 LCR16a loci, 25 are protein-coding and `ID_9` is the single lncRNA (7.1 kb against a
25.5 kb median, 0 reads).
⛔ One inference that does NOT follow, and was tested: this conflation does not explain which members we miss.
Detection on the Soto slice is class-flat — pseudogene 0.935, protein_coding 0.929, lncRNA 0.872 (register 742).

---

## Part 0c — Read-throughs: what is proven, what is excluded, and where the fix lives

**The claim to make:** *30 of 42 gorilla read-through junctions are used by chimpanzee reads too.* That is the
only POSITIVE evidence and it is the one to lead with — an independent individual, library and species cannot
be produced by a systematic of our preparation.

**Everything else EXCLUDES an artefact route without demonstrating biology.** Say it that way:
| route excluded | how |
|---|---|
| concatemer chimera | the reads are **FLNC** — 5′ primer, 3′ primer and polyA required upstream |
| reverse-transcription template switch | junctions are **canonical**, and microhomology shows **no long tail** (≥6 bp: 2 % vs 4 % in real introns) |
| cross-copy mis-chain between paralogs | only **4 of 46** junctions have duplicate-linked flanks; the §6fw guard removes those and 2 opposite-strand joins |
| low-evidence chaining noise | ≥ 3 independent molecules, one strand; for the top junctions **full-length molecules span the whole structure** (5,043 of 5,166) |

**Where the fix lives: downstream, never in the aligner.** All three aligner levers are measured and all three
fail. `-G 50k` touches **1 of 46** junctions while damaging **8.2 %** of gorilla transcripts. `--junc-bed` from
the annotation is an uncalibrated per-copy prior that favours better-annotated copies and suppresses what O3
looks for. Disabling the long join (`-r 500,500`) removes **20 of 31 conserved junctions and 95 % of their
reads** — and does so incidentally, by deleting **98 % of ALL introns over 500 bp**. ⟹ *An aligner option
applies globally and invisibly and cannot be inspected afterwards; a downstream certificate applies to a named
set and leaves a record.*

⚠ **Two things not to over-claim.** The 12 junctions that do NOT replicate in chimpanzee are **not** shown to
be artefacts — NPIP is fast-evolving and a junction may be lineage-specific or unexpressed there. And the
cross-species positive control was **not expression-matched** (test 0.71 vs control 0.39), so that gap is not
evidence that read-throughs are more conserved than ordinary introns.

---

## Part 0d — What he is actually asking for: a FAMILY-AWARE assembler, and its IGV file

⚠⚠ **CORRECTED 2026-09-08 (user).** An earlier reading of comment 1 had it as "justify keeping a substrate
component two tools already provide". That is **not** what he means. Across other emails he asks **which copy
an isoform belongs to** and **where the GTF/GFF of produced isoforms is so he can check it in IGV**. Taken
together with "this is not a transcriptome assembler", the request is coherent and specific:

> **He wants a transcriptome assembler — one that is family-aware.** It should produce isoforms, assign each
> to a COPY, borrow information across members of a family, and resolve statistically which reads support
> which copy when the alignment scores are equal.

"Not a transcriptome assembler" is a statement about the **novelty claim** (a plain assembler is not the
contribution), not about the **deliverable** (an isoform GTF is exactly the deliverable). Comment 1's
comparison to flair/StringTie is then the natural benchmark: *those tools produce isoforms and cannot tell you
which copy; show what yours adds.*

**Where we stand against that specification, item by item:**
| what he asks for | status |
|---|---|
| produce isoforms | ⭐ **exists**: `copy_assign --gtf` writes a FLAIR-style transcript+exon GTF of every de-novo isoform in the swept regions, IGV-loadable, annotation-free (intron-chain collapse + the canonical gate) |
| say which copy each isoform belongs to | ⚖️ **partly**: each transcript carries `family_id`, `copy_index` and `multicopy "true"` — but that tag comes from **where the isoform assembled**, i.e. position, not from the read evidence |
| resolve reads statistically when scores tie | ⭐ **exists and is the contribution**: the origin certificate + posterior, assign-or-abstain, never 1/k |
| borrow information across family members | ⛔ **assembly-side borrowing is dead/inert** (register: `consensus.rs` has no call sites; `rescue_thin_loci_iterative` yields 0 copies in both shipped catalogs). Borrowing is live **only** in assignment |
| view it in IGV | ⭐ `bench/igv_tracks.py` writes `<out>.tagged.bam` with `cp:Z:<family>_c<idx>` per read — IGV "Group by tag" / "Color by tag" |

⟹ ⭐ **The one real gap is the join.** Isoforms are tagged by POSITION; reads are assigned to copies by
CERTIFICATE. Nothing yet aggregates the second onto the first, i.e. *"this isoform is supported by N reads
assigned to copy 3 with a certificate, M that abstain, and K that were assigned elsewhere"*. That is a small
piece of work over two files we already emit (`<out>.gtf` and `<out>.assignments.tsv`), and it is precisely
the sentence he keeps asking for. **Do this before Wednesday.**

⚠ The substitution test against flair/StringTie is still worth running, but it is now a **secondary**
question — it asks whether our isoforms are as good as theirs, when the point is that ours carry a copy
assignment and theirs cannot.

---

## Part 0e — "Not a transcriptome assembler" AND "compare to flair/StringTie": the earlier reading

The two instructions look contradictory. They are not, and the resolution is in his own earlier words: on
**2026-06-25** the reframe he asked for was that *the StringTie-clone assembler is the SUBSTRATE that produces
transcripts and loci, NOT the contribution* — the contribution is the family definition and the read-to-copy
assignment. Nothing has changed since.

⟹ **Comment 1 is therefore not "is your assembler better than flair?"** — on that question the answer is "we
do not claim it is, and isoform-level metrics are the wrong test". It is the sharper question: **"you built a
substrate component that two published tools already provide. Justify keeping yours."**

**The comparison that answers it is a SUBSTITUTION test, scored on families, not on transcripts.** Steps 6–7
consume one FASTA of one sequence per locus and do not care where it came from. So: run the family definition
three times, once on our stages 1–5, once on StringTie's loci, once on flair's, everything downstream held
fixed, and compare the FAMILY output — the locus set, the components, the hierarchy, and the O1 metrics
(sensitivity, specificity, bipartite coverage) against a truth family.

⭐ **Both outcomes are good, which is why this is worth running.**
- **Families are the same** ⟹ drop stages 1–5, cite StringTie for the substrate, and the thesis gets smaller
  while the novel part (the definition and the assignment) stands alone and unencumbered. That is a *better*
  thesis, not a worse one, and it removes the reviewer's largest "you reimplemented a solved problem" target.
- **Families degrade** ⟹ there is now a measured reason the custom stages exist, which is exactly what he
  asked for, and it is stated in the currency he cares about.

⚠ **Two things already in the record point at the first outcome.** A direct comparison found the foundations
**near-identical**: introns 99.3 / 98.0, bundles 3,351 / 3,430 exact, and **nodes 96.3 % byte-identical**, with
zero junctions unique to our path — the measured gap was over-enumeration of real introns, not a different
object. And the assembler stack was **already retired as dead code**: its 41 modules are unreachable from the
five thesis binaries. The substrate we ship is a thin locus-builder, not a rival assembler, which is why the
substitution is plausible and cheap to test.

⚠ **What is NOT yet done:** the substitution itself. Neither StringTie nor flair has been run into steps 6–7,
so no family-level comparison exists. That is the single open item behind comment 1 and it should be run
before Wednesday if anything is.

---

## Part 0f — How we guard against tuning: hold a substrate back

⭐ **The strongest single answer to "you tuned until it worked" is a worked example where we did the opposite,
and lost.** 2026-09-08, cross-family exon-overlap rule (§6gn, register 760):

The rule was **designed** on the gorilla read-throughs and checked on three substrates — gorilla NPIP, the
gorilla catalog, the human Soto slice. **All three agreed**: membership metrics unchanged, Soto specificity up
0.649 → 0.663, manufactured read-throughs 32 → 1. On that evidence the default was flipped ON.

The pre-registration had named a **fourth** substrate, human chr16+18, with the sentence *"P4 is the real test:
the other three substrates were used to design the rule, chr16+18 was not."* It **failed**: two real members,
`NPIPA1` and `NPIPA6`, were stripped of the sequence that qualified them, and NPIP sensitivity fell
**26/26 → 24/26**. The default was **reverted in the same run** and verified byte-identical to the shipped
catalog; the rule remains an opt-in flag.

⟹ **The agreement of the three development substrates carried no information about the rule's correctness.**
Only the held-out one did. That is the procedure to state when he asks about overfitting: *the substrate that
decides a change is named in the pre-registration before the change, and it is one the change was not built
on.* We have the failed instance to show, which is worth more than a list of successes.

---

## Part 0g — Yes, this is assembler-shaped work. Say so, then say what full-length changes

⭐ **Concede the original point first, because it was right.** A StringTie-clone assembler WAS built here, and
it is now **retired as dead code** — 41 modules unreachable from the five thesis binaries. Re-implementing a
solved problem was a real mistake and pretending otherwise costs credibility for nothing.

⭐⭐ **Then make the distinction that actually holds: with full-length reads, assembly is not reconstruction.**
Short-read assembly is a genuine inference problem — transcripts must be rebuilt from fragments, which is why
StringTie carries a flow model. **FLNC IsoSeq reads are already transcript observations**: 5′ primer, 3′
primer and polyA are required for a read to exist at all, and the median read here is **2,884 bp**. What
remains is not reconstruction but three much shallower operations: **group identical structures, discard
artefacts, and decide which copy each came from.**

**Measured, on the NPIP family: we group, we do not infer.**
| | |
|---|---|
| molecules | 15,960 |
| **distinct spliced structures OBSERVED in reads** | **1,818** |
| spliced isoforms we emit | **703** |
| isoforms we emit that no read shows | **0** |
⟹ Every isoform is an observed read structure; we emit a **subset** of what the data contains and invent
nothing. StringTie's flow model, by contrast, **can and does emit junction combinations no single read shows**
— that is inference, appropriate to fragmentary coverage and unnecessary here. It is also why it reports **3**
transcripts at NPIP copy 27 where **32 distinct full-length structures** were observed.

⟹ **The sentence to use.** *"With full-length reads the assembly problem reduces to grouping and attribution.
The grouping is FLAIR's collapse and we claim nothing new for it — if an assembler's output is preferable,
steps 6–7 consume it unchanged. The attribution — which copy each isoform came from — is undefined for every
existing assembler, because none of them has a notion of copies. That is the contribution, and it is the part
you have been asking to see."*

⚠ Do not overstate the grouping: StringTie recovers **more read-supported junctions** than we do (58.6 % vs
42.1 %, §6gj) precisely because it will assert a combination no read carries. Conservative and exhaustive are
different virtues; ours is the one a per-copy assignment needs, not the one that maximises junction recall.

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

⭐ **2026-09-05 — the excision experiment on the shipped certificate (§6ff, `docs/PREREG_excision_2026-09-05.md`):**
remove an NPIP copy 98.5 % / 98.7 % identical to its neighbour: 100 % of its reads abstain ("no candidate explains
this read"), none is absorbed; the abstaining reads carry consistent mismatch sites 7–100× above the controls that
point at the missing copy; the ZNF569-like copy (99.1 %) is reconstructed at 0.994 from its orphaned reads. The
wall: NPIP 13's expressed segment is 99.85 % identical to copy 12 (units 0.966 overall), 12 of 32 reads are
absorbed with a VALID certificate — a missing copy inside sequencing error is an allele. Say the wall as local
identity, not as a percentage of the unit.

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

**Yes — the SD-core definition (accepted 2026-09-05, §6ev; `docs/THESIS_OBJECTIVES.md`).** A family is a set of
annotated loci, pre-clustered by MCL over their pairwise genomic homology, in which each member shares a
duplicated core segment with at least half of the others; a member's unit is its read-supported exon chain, its
locus the read-supported extent clipped at every other catalog unit (§6fh, §6fm). The core comes from SEDEF's
segmental-duplication calls where they exist and **can be derived from the catalog's own alignments where they
do not** (`--core-from-paf`, §6fo: on NPIP the two agree on 30 of 32 members and give the same LCR16a core).
The old E_r/γ-quasi-clique definition is kept opt-in (`gw_family_catalog`) so the change can be audited.

**On his own example it is the better representation.** Old: NPIP fragmented into 5–6 families, 14–30 of the
31 loci recovered depending on seeds (§5j, §6be). New: 31/31 loci, one family of 29 units on three contigs,
LCR16u (SMG1P/SLC7A5P/PDXDC) separate, the 9 ABCC1/SORL1 chimeric models trimmed to the 23-kb LCR16a core or
dropped and kept as candidates with `member_status` (§6eh–§6ei, §6fh). Block versus family is a certificate,
29 versus 48 loci (§6eg).

**On false merges it is not worse where the benchmark can see, and it sees more.** Paired on the Soto slice
(`docs/O1_DEFINITION_SWITCH.md`): pair precision in the adjudicable [0.90, 1) band **0.974 → 0.954**, CIs
overlapping; recall on pairs both methods detect **0.874 → 0.940**; recall on all Soto pairs **0.173 → 0.580**;
families exact 21/33 → 40/56. SEDEF corroborates 39 of 46 three-contig clusters, the repeat library names 3
artefacts (§6dy–§6dz). ⚠ The old 1.33 % false-merge rate was a different instrument (window sets) and is not
comparable; do not quote the two side by side. ⚠ The 21 genome-wide MCL cuts certified by shared reads at
≥ 0.90 identity (§6fn, row 721) are the known false SPLITS; say so before he asks.

⚠ **The binding constraint moved with the node.** Under the old definition node construction cost ~58 % of the
loss (§5e); under the new one the node is the read-supported chain and the remaining losses are the annotation
(loci the GFF has no model for, §6ev P3) and the special cases recorded on 2026-09-06 (partial paralogues below
the coverage threshold, nested units).

### Q2. "Are these real families, or artifacts and overfitting?"

See Part 1. The summary line: **defensible narrowly, not broadly** — and since 2026-09-05 with three external
instruments instead of one: SEDEF (the core rule is a refinement over an independent SD call, 39/46 corroborated),
the RepeatMasker library (`rep_frac` names the repeat-clique clusters), and the CHM13 landing of every anchored
family (real families land on one stem, artefacts scatter, §6eh). ⚠ **The paragraph below on the ~0.83 band is
about the OLD edge set**; under MCL the pre-clustering edges are the same PAF, but membership is decided by the
core rule, whose evidence is at ≥ 0.90 by construction (SEDEF) — the band objection now applies to the
pre-clustering, not to the family.

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

**Answered by the genomic read-star (§6fd).** Each molecule is aligned splice-aware to each candidate's LOCUS,
so isoform structure aligns (introns are `N`, not edits) and only sequence the locus lacks counts against the
origin certificate — the isoform question and the origin question are separated by construction. Two copies
producing the same isoforms therefore tie at K = 0 columns (26 NPIP reads) or are told apart by the columns the
read covers; the structure of the read never decides. The unit (O1's node) is the read-supported chain, i.e. the
expressed isoform union at that locus (§6el), not one intron chain. ⚠ The 2026-09-05 residue: 456 NPIP-proper
reads abstain because no candidate locus explains their whole sequence.

### Q5. "Two tandem near-identical copies share a read. What happens?"

**Three cases, each with a rule (2026-09-05/06).** (i) The read aligns identically to both: a K = 0 tie,
reported as such — 26 of NPIP's unit reads, 90 of 24,462 on the paired 35 families (§6fj). (ii) A read-through
molecule spans two copies: the locus of each copy is clipped at the neighbour's chain, so the read leaves
unaligned bases on both targets and abstains with `origin_rejected` (§6fh; the unclipped form turned 226 MCL4
assignments into false ties, row 714). (iii) Two families' units over one place (44 genome-wide pairs, row 721)
are recorded as a special case, not resolved. ⚠ The paragraph below is the OLD node's measurement (152
primaries, 0.306 %) and its intron-annotation argument; keep it as history, do not quote it as current.

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

⭐ **Sizes, his own instrument (bipartite 1:1 matching, §6fs, 2026-09-06).** On the Soto slice the new node matches
272 of 362 members with 78 % of pairs within 2× and a median ratio of 1.00 (old node: median 0.54, 104
truncated); at NPIP the unit is the size of the LCR16a core (median 0.89) and the locus extent the size of the
annotated gene (1.02). Size stays a measurement, never a filter: 73 % of Soto's own families are more than 2×
size-heterogeneous.

**The premise is false and that is measurable** (§6ay) — ⚠ but the coverage statistic used to
answer it **was one-sided** until `cov_longer` was emitted (§6ba). Answer the question, then
volunteer the flaw in the instrument that answered it.

### Q8. "1/k is a bad assignment for tied multimappers."

**Agreed, and never used.** O2 is **assign-or-abstain**; the objective **provably decomposes**, so
the shipped per-read gate **is** the optimum and no joint estimator can beat it. **K = 0 abstention
is entailed, not chosen.** Empirically the EM changes **zero** of the gate's decisions on reads
that carry evidence.

⭐⭐ **Two metrics (2026-09-06, §6fr).** O1: of the 26 LCR16a copies in the three contigs (the family's cores aligned
back to the assembly at ≥ 90 %, expressed or not), 25 are family members, 24 have units, 0 are unannotated or in
another family; the one member without a unit has a single read. O2, on the reads the aligner could not place
(known-origin simulation): the true copy is in the posterior's tie set for 100 % of 1,380 contested reads and is
its unique maximum for 98.3 %; for the 197 reads with evidence below the certificate's α the most likely copy is
the right one 197 times; K = 0 ties hold the truth in a set of two at P = 0.5 each. The posterior answers "which
copy", the certificate answers "may I claim it".

⭐⭐ **One sensitivity, one specificity for NPIP (2026-09-06, §6fp–§6fq; PREREG npip_known_origin).** The machinery
decides the contested reads (MAPQ < 60); an uncontested read keeps its certified call or takes its placement.
Reads of known origin (5,000 simulated from the 25 NPIP units at the substrate's own error rate): **95.5 %
right, 0.12 % wrong, precision 0.9987**; the contested subcategory (1,382 reads) 84.0 % right, 0 wrong; the
equal-best-score reads (39 at MAPQ 0) 16 right, 0 wrong, 23 K = 0 ties against minimap2's 46 % / 54 % coin flip;
minimap2 alone 92.3 % right / 7.6 % wrong. Real NPIP, 1,000 unit reads of the 25 kept members: **701 assigned
(70.1 %), 29 K = 0 ties, 270 abstain; the 62 audited anchors 33 right / 0 wrong / 29 abstain** (contested
anchors 5 / 0 / 28). What abstains on real data is decomposed, not hidden: read-throughs the catalog has no
object for (71), divergence above 1 % (51, O3's material), the certificate's 0.3 % constant (43). Row 724: MAPQ
60 is not proof of origin (4 % of simulated MAPQ-60 reads sit at the wrong copy; the certificate corrects 139
of 145), which is why the certified call comes first.

⭐ **O2's population is the reads minimap2 cannot place; quote it on those first (2026-09-06).** NPIP-proper
(`sweep_v14`, the 25 kept members): 1,000 unit reads, **472 (47 %) below MAPQ 60**, 34 at MAPQ 0. On the 472:
assigned 182 (38.6 %; 181 certified against 2–7 competitors), of which **16 to a different copy than minimap2's
primary**; K = 0 ties 25; abstain 265 (190 because no candidate explains the read — O3's material). On the 34
MAPQ-0 reads: 1 assigned, 19 K = 0 ties, 14 abstain. Truth on the contested reads = the 33 audited junction
anchors below MAPQ 60: **5 assigned, 5 right, 0 wrong, 28 abstain** — sensitivity 15 %, specificity 100 %. The
unique-mapper agreement (5,208/5,208) is only the sanity check that O2 never contradicts a unique placement.

⭐ **Since the read-star (§6fa–§6fj) the claim is stronger than abstention, and it is measured.** O2 assigns a
molecule to a copy only with a certificate: its own edits against the best candidate's locus are within
sequencing error (the origin certificate) and every competitor is rejected on the columns the read covers (the
pairwise certificate). NPIP, `sweep_v14`: 62 audited junction anchors → 14 assigned, 14 right, 0 wrong, 48
abstain; MAPQ-60 placement agreement 5,208/5,208; **587 reads assigned against 2–7+ candidates by the pairwise
certificate**, 196 of 504 MAPQ<60 reads assigned; 456 NPIP-proper reads abstain because no locus explains them.
Paired 35 families: 78.8 % assigned at 18,772/18,772 agreement. The abundance per copy is the sum of the
per-molecule posteriors (`n_reads_soft`), with a half-width under 0.05 for all 26 NPIP copies. ⚠ The old line
"O2 does not beat minimap2" described the PSV-column era; retire it. What O2 still cannot do is stated with a
number: a copy whose expressed segment is within sequencing error of another's is absorbed with a valid
certificate (NPIP 13 → 12 at 99.85 %, row 711).

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
7. **O1's remaining losses are the annotation and three recorded special cases**, not the definition: loci
   the GFF has no model for (§6ev P3), partial paralogues below the coverage threshold with hundreds of
   cross-mapping reads, nested units of two families, and members whose reads fall outside the chain (§6fn).
   The 2026-09-05 thresholds are a measured plateau (§6ez), not a proof.
8. **O2's wall is local identity within sequencing error**: a copy whose expressed segment is 99.85 %
   identical to another's is absorbed with a valid certificate (row 711). Above that wall the certificate
   is measured (0 wrong anchors, 1.0000 placement agreement); below it, the reads are an allele.
9. **O3 flags a sequence, not a copy.** A flag means ≥ 3 reads consistently at least 0.7 % from every locus
   the assembly holds (§6fl–§6fn); RNA alone cannot say whether that is a copy absent from the reference or
   a diverged haplotype. On the intact catalogs 34.7 % (3 contigs) and 13.9 % (genome-wide) of the candidate
   pairs carry it — above the pre-registered 25 % on the three contigs (row 718) — and the reference-absent
   class proper is small (4 and 23 unannotated loci); most "missing" origins are annotated loci without a
   unit or MCL cuts. The old collapse-deficit screen (0/816) is a different instrument and is retired.

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
| Q1 method exists | ⭐ SD-core definition; NPIP 31/31 one family; Soto precision within CI, recall ×3.4 | §6ev, §6fh, §6fo, `O1_DEFINITION_SWITCH.md` |
| Q2 real vs overfit | ⚠ **narrowly** defensible; three external instruments | §6dy–§6dz, §6eh, Part 1 |
| Q3 borrowing | ⛔ no (inert) / ⭐ `E_c` splits | §6bd, §6bi |
| Q4 isoforms | ⭐ structure and origin separated by the genomic read-star | §6fd |
| Q5 tandem reads | ⭐ K = 0 tie / clipped locus / nested special case | §6fh, §6fj, row 721 |
| Q6 portability | ⭐ tissue/animal · ⛔ apes | §4l, §5p |
| Q7 boundaries | ⭐ premise false | §6ay, §6ba |
| Q8 1/k | ⭐ never used; certified assignment measured (0 wrong anchors, 1.0000 agreement) | §6fa–§6fj, `sweep_v14` |
| Q9 NPIP subfamilies | ⭐ answered on his own example | §5j, fam72 review |
| Q10 non-canonical | ⭐ he is right; recurs | §6au, §6av, §6aw |
| Q11 PSV | ⭐ answered | §6aj |

## Part 0h — THE ADVISOR'S CLARIFICATION OF WHAT THE PROJECT IS (relayed by the user, 2026-09-09)

Verbatim intent, and it supersedes Parts 0d–0g where they differ:

> The project **is an assembler** — but one more similar to **flair or isoseq collapse than to StringTie**.
> The entire focus is **multicopy gene families and better use of multi-mapping reads**, because that is a
> common gap in all other assemblers. **In the end we produce a GTF.** It does not matter whether it contains
> the transcripts other assemblers get easily; **what sets it apart is that it should contain the transcripts
> other assemblers struggle to produce.** O2 is an extension: take **all the ambiguous reads (same AS — a coin
> toss which is the "best" alignment)** and infer which copy they actually came from. **Uniquely mapped
> transcripts are irrelevant to O2.**

### What this settles
1. **The deliverable is a GTF, judged at the hard loci only.** Junction recovery against StringTie (§6gj,
   42.1 % vs 58.6 %) is **not a loss** — it is a comparison on the easy transcripts, which the advisor has just
   declared out of scope. ⛔ Do not present it as a cost.
2. **The right comparison is the one §6gt already made by accident**: isoseq collapse drops 23.1 % of the
   family-region molecules and **30.2 % of the ones we abstain on** because it discards secondaries. The
   transcripts carried by those molecules are *exactly* "the transcripts other assemblers struggle to produce".
   The benchmark becomes: **at multi-copy loci, which transcripts do we emit that isoseq collapse / flair emit
   nothing for, and are they attributed to the right copy?**
3. **O2's scope is confirmed as §6gv shipped it**: AS-tied multimappers only (`--as-tied-only`), unique mappers
   irrelevant. The honest population is the CONTESTED one — under the 2026-09-09 defaults (§6hk): **33 molecules
   on gorilla MCL1 (4 assigned / 28 tied / 1 ambiguous), 1,118 on human MCL0 (262 / 512 / 344)**, held-back
   MCL58 40 (2 / 38 / 0); on the certificate-independent in-catalog AS-tied pool, 262 / 1,643 = 15.9 % (human),
   4 / 122 (MCL1) — and the "38.7 % assigned" headline is retired. Precision: excision of every human copy in
   turn, 247/262 = 94.3 % of the assignments abstain (§6hg/§6hk); the hard-locus bakeoff (§6hh) carries 85 % of the
   AS-tied molecules vs flair 47 % / StringTie 44 % / isoseq 73 %. Two widenings were measured and refuted
   (aligner disagreement §6hd, indel PSV columns §6he); the pairwise best (§6hj) shipped.
   **The deliverable (§6hn–§6hp, default since 2026-09-09):** `copy_assign --gtf` emits the GTF O2 believes —
   family isoforms grouped across copies by lifting their intron chains, placed only where a unique mapper or
   a certified read backs them (phantoms dropped), certified isoforms lifted to O2's copy, and coin-toss
   isoforms emitted ONCE with `copies "A,B[,outside]"`. Measured against the same rule on every tool's GTF
   (human MCL0): phantom transcripts at twin copies ours 2 / flair 4 / StringTie 15 / isoseq 45; coin-toss
   isoforms given one arbitrary address ours 0 (40 sets) / flair 76 / StringTie 20 / isoseq 200; certified
   isoforms at O2's copy ours 100 % vs 78–86 %. Isoform-level pooling was tested and refuted as a placer
   (71 % excision, row 801); "similar copies produce the same isoform" is measurable: 43 shared groups.
4. **flair-like, not StringTie-like**: group observed full-length structures (intron-chain collapse), do not
   run a flow model that asserts chains no read carries. That is already the design (§6gr).

### What it changes in the standing framing
- `THESIS_OBJECTIVES.md` "NOT an assembler" ⟹ **"an assembler for the loci other assemblers cannot resolve"**.
- The bakeoff metric (`bench/tool_bakeoff.py`) must gain a **hard-locus restriction**: score only molecules that
  are AS-tied, and count transcripts emitted per copy where the competitor emitted none.
- ⚠ The user decides the wording change to `THESIS_OBJECTIVES.md`; this Part records the instruction.
