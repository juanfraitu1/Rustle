# O1 false-positive hardening: rules that survived falsification

> ⚠ **CATALOG PROVENANCE.** Figures quoted per **/494** (families) or **/1415** (copies) describe the
> **superseded** 2026-07-17 catalog, which **no invocation of the current binary reproduces** — it was
> built with `refine` on, a default removed on 2026-08-20. The current default emits **627 families /
> 2,019 copies**. Re-measure before quoting: see [`NUMBERS.md`](NUMBERS.md) and
> [`o1_catalog_provenance.md`](o1_catalog_provenance.md).

**Status: 2026-08-16.** This note extends
[`o1_coverage_repair.md`](o1_coverage_repair.md) and
[`o1_read_evidence_repair.md`](o1_read_evidence_repair.md). It does not replace their
negative results.

## Decision

The best new **hard-rule candidate** is a transcript-orientation guard:

> When both RNA O1 representatives are normalized to transcript 5′→3′ orientation,
> an `E_r` edge exists only if at least one identity/coverage-passing PAF record has
> strand `+`.

This is pair-local, annotation-free, seed-order invariant, and adds no fitted numeric
threshold. A `-`-only match between two sense-oriented transcripts is
reverse-complement homology, which is evidence for an inverted repeat rather than for
two homologous transcripts. An inverted *genomic* duplication still aligns `+` after
each expressed copy is normalized to its own transcript direction.

**⭐ SHIPPED AS THE DEFAULT 2026-08-19** (`793 passed / 0 failed / 11 ignored`, baseline 792 + the new
lock-in test). Opt out with `gw_family_catalog --no-rna-forward-only`; `--rna-forward-only` is retained
so an explicit request still conflict-checks against `--from-genome`.

⚠⚠ **The flip is at the RNA ENTRY POINTS, never on `RefineParams::default()`.** That struct is
**substrate-agnostic** — it configures both the RNA exon-sum path and the reference-oriented DNA path,
where a `-` record is a **real inverted segmental duplication**. Flipping the type default silently
applied the RNA guard to DNA and dropped an inverted duplication; the existing
`genome_mode_grouping_keeps_an_inverted_duplication` test caught it. A new test,
`refine_params_default_is_orientation_agnostic`, now locks that in and says so in its message. It filters PAF records before both the ordinary
single-record edge and the optional summed-coverage edge are formed, and writes its
effective value as `alignment_orientation` in the E_r rule certificate. The current
evidence consists of the frozen-PAF test, the corrected seven-family panel, and a fresh
paired 12-region end-to-end CLI safety run. A fresh whole-genome GGO/HSA comparison is
still required before changing the default.

## Formal rule

Let `P(a,b)` be the PAF records between transcript-oriented representatives `a` and
`b`, and let `R0(r)` be the current per-record predicate:

```text
identity(r) >= 0.60
and aligned_span_on_shorter_axis(r) / min(len(a), len(b)) >= 0.50
```

The candidate edge is

```text
E_r_plus(a,b) iff exists r in P(a,b): R0(r) and strand(r) = '+'
```

Pair semantics remain **ANY passing record**. A reverse record does not veto a pair
that also has an independently passing forward record.

### Typed-substrate restriction

Apply this only where Rustle owns the orientation invariant:

- RNA exon-sum representatives;
- RNA-locus genomic spans after `refine_copy_seq` normalizes a minus-strand locus.

Do **not** apply it to `gw_family_catalog --from-genome`. That mode deliberately keeps
read-free DNA intervals in genomic-plus orientation, where a `-` PAF record can be a
real inverted segmental duplication. Do not apply it to an arbitrary imported FASTA
unless its sequences are explicitly declared transcript-oriented.

## ⭐ The whole-genome GGO/HSA comparison (2026-08-19) — the blocking experiment, done offline

**T8: offline re-derivation over the shipped reps, not the binary.** The γ step is deliberately NOT
reimplemented (that engine is Louvain-based and reimplementing it offline has already produced one
Simpson's-paradox error here). Every statement below is **γ-independent**: edge removal, disconnection
and isolation hold whatever γ does, because γ can only split a component further, never rejoin one.
Species are never pooled.

| | GGO | HSA |
|---|---:|---:|
| families / copies | 494 / 1,415 | 394 / 1,220 |
| within-family `E_r` edges | 2,474 | 2,402 |
| minus-only, removed by the guard | **55 = 0.0222** | **153 = 0.0637** |
| families losing ≥1 edge | 42 = 0.0850 | 27 = 0.0685 |
| families that DISCONNECT | 42 = 0.0850 | 27 = 0.0685 |
| copies left ISOLATED | 73 = 0.0516 | 42 = 0.0344 |
| families that DISSOLVE | **31 = 0.0628** | **17 = 0.0431** |

### The dissolutions are not collateral — they are ~90% provably-wrong families

**All 31 dissolved GGO families are 2-copy**, which is arithmetic: a 2-copy family rests on one edge.
The question is whether that edge should have existed. Splitting the removals:

| variant | edges | families touched | dissolved |
|---|---:|---:|---:|
| BROAD — all minus-only (GGO) | 55 | 42 | 31 = 0.0628 |
| **NARROW — minus-only AND overlapping antisense** (GGO) | 33 | 33 | **28 = 0.0567** |
| BROAD (HSA) | 153 | 27 | 17 = 0.0431 |
| NARROW (HSA) | 15 | 15 | 14 = 0.0355 |

**28 of the 31 GGO dissolutions are overlapping ANTISENSE pairs** — the same DNA read in opposite
directions. Those are **provably not two homologous transcripts**, so deleting them is correct.
⟹ **~5.67% of GGO "families" are a gene paired with its own overlapping antisense partner, not
paralogues at all.** That is a precision defect of the same order as the definitional hole's 8.30%
exposure — and unlike that hole, it is fixable by a rule that is pair-local, threshold-free, and
already implemented.

### It also splits every known false-merge family

`GWFAM210` (the MRPS17 **AluY** hub), `GWFAM264`, `GWFAM82`, `GWFAM85` all go 1 → 2 or 1 → 3 — i.e.
the guard independently splits the anti-FP characterisation's Groups 1–3.

### ⚠ What is NOT corroborated

* **Only 1 of the 31 GGO dissolutions** is independently flagged as repeat-bridged by the
  genome-anchored veto (`gmult ≥ 50`). The two guards target different mechanisms, so this is expected
  — but it means **3 non-antisense dissolutions rest on the strand observation alone**.
* ⚠ **Unexplained species asymmetry:** HSA removes **6.37%** of within-family edges against GGO's
  **2.22%**, nearly 3×, yet only **15/153 = 9.80%** of HSA's removals are antisense overlaps versus
  **33/55 = 60.00%** of GGO's. HSA's non-antisense removals mostly split families rather than dissolve
  them. This should be understood before the default moves on the human arm.
* This is **T8**. It bounds the change; it does not replace running the binary both ways.

## Frozen-arm result

Reproducer: [`bench/o1_orientation_gate.py`](../bench/o1_orientation_gate.py).

The historical seven-family builder concatenated plus-genome exon blocks but did not
reverse-complement minus-strand loci. The reproducer first repairs that panel using the
dominant strand of contained spliced reads. This changes orientation only; sequence
content and the R0 identity/coverage predicate are unchanged.

### False edges rejected

| stratum | rejected / scored FP |
|---|---:|
| GGO catalog nodes | 6 / 14 |
| GGO node × HSA curated transcript | 2 / 11 |
| HSA curated transcripts | 12 / 28 |
| HSA census curated transcripts | 0 / 2 |
| HSA catalog nodes | 9 / 19 |
| **total** | **29 / 74 (39.2%)** |

Those 29 rows cover **12 of 24 independent FP mechanism components**. Fifteen are
edges the shipped GGO/HSA catalogs actually use. This is materially stronger on the
thesis substrate than coverage-of-longer at `0.20`, which rejected only 1/14 GGO
catalog FPs and damaged many true families.

### Positive-edge audit

On the non-panel scored TP arms, the guard removes only four of 9,032 edges:

- `OR2G6–OR2L2` and `OR2G6–OR2M4`, whose only frozen truth is the broad HGNC
  “olfactory receptors family 2” group;
- `EIF4H–OSBP2`, whose HGNC truth is the generic “microRNA protein-coding host genes”
  group and whose proteins are in different clusters;
- `LOC101151425–LOC101135168`, which shares a protein cluster and is the one credible
  adverse case that still needs manual adjudication.

The corrected seven-family graph gives the more informative family-level result:

| panel | R0 edges | forward edges | component sizes before → after |
|---|---:|---:|---|
| NPIP | 259 | 235 | `28 → 27+1` |
| TBC1D3 | 55 | 55 | `11 → 11` |
| RABL2 | 1 | 1 | `2 → 2` |
| MAGEA | 60 | 60 | `12 → 12` |
| GSTM | 1 | 1 | `2+1+1 → 2+1+1` |
| HERC2 | 16 | 16 | `10+1 → 10+1` |

The isolated “NPIP” node is `chr16:15105563-15115512`; it overlaps
`LOC100505915/PDXDC1`, has no NPIP seed annotation, and is connected to NPIP only by
reverse-complement records. Its removal is therefore a panel-label correction, not a
named NPIP loss. APOBEC3 has no R0 edge in this panel and supplies no positive test.

### Fresh regional end-to-end safety check (2026-08-16)

The current release binary was run twice, with identical BAM/FASTA inputs and the
orientation flag off/on, over seven strict single-copy controls (TBP, POLR2A, GTF2H1,
SF1, TFRC, HMBS, PSMB6) and five positive regions (RABL2, APOBEC3, GSTM, MAGEA, HERC2).
All seven controls emitted zero families in both arms. Every positive retained the
same family and copy counts: RABL2 1/5, APOBEC3 1/2, GSTM 1/4, MAGEA 2/11, and HERC2
1/3 (families/copies). The MAGEA copy representative differed between runs at one
same-locus alternative model, but the partition counts did not; this is existing
representative-selection variability, not an orientation-edge loss.

This is a safety result, not an FP-kill result: the regional slices contain none of
the reverse-only contaminants found in the whole-catalog frozen PAFs. Do not use it to
replace the 29/74 scored-edge result or the pending whole-genome comparison.

### Expanded known-family graph audit (2026-08-16)

The original seven-family audit is now complemented by 19 post-hoc family-bearing
graphs from the frozen whole-catalog GGO/HSA PAFs. The expanded set includes C4,
TSPYL, GSTM, RGPD/RANBP2, ANKRD18, HERC2, MAGEA, DAZ, APOBEC3, GOLGA6/8, TBC1D3,
RABL2, AMY1, PCDHB, RBMY1, and two NBPF strata. It deliberately includes the known
NBPF/TTC6-DNAH14 repeat-bridge mixture rather than selecting only clean pairs.

Results:

* primary graphs retain **108/133** emitted nodes and all **75** independently named
  target-family members;
* **16** nodes independently annotated as unrelated genes, **1** documented broad-family
  outgroup, and **8** untyped/unreachable nodes are withheld from recent-copy primary
  membership but remain in audit graphs;
* all **15/15 testable** cases (at least two independently named target members) keep
  those named members connected; four further cases are annotation-limited, not scored
  as passes; and
* the orientation rule removes **15** within-family edges, all in the adversarial NBPF
  repeat-bridge graph. That graph contracts from 17 audit nodes to the single named
  NBPF4 node after typing, while its unrelated TTC6, DNAH14, CEP152, TTBK2, and other
  genes remain visible in the audit view.

This is a **purity/typing audit conditional on emission**, not discovery recall: the
family-bearing catalog components were selected after Rustle emitted them. Its value is
showing that the typed output removes known mixtures without deleting named positives;
it cannot measure how many families Rustle never found. Reproducer, tables, and paired
Bandage graphs: [`bench/o1_expanded_family_audit.py`](../bench/o1_expanded_family_audit.py)
and [`bench/o1_expanded_family_audit/`](../bench/o1_expanded_family_audit/).

### Fresh-emission falsification of the audit (2026-08-16)

The expanded audit's frozen nodes were challenged by going one step upstream. The current
Rustle release binary, with `--rna-forward-only`, was rerun on regional BAMs extracted from
the original whole-genome-aligned GGO and HSA BAMs. The old node ids, family ids, gene labels,
and dispositions were not inputs to Rustle. They were applied only after emission by genomic
coordinate overlap.

Fresh Rustle output recovers **124/133** frozen loci, including **72/75** independently named
target nodes. Crucially, **14/16** nodes labelled as unrelated conflicting genes are also
rebuilt: those loci are not inventions of the audit. But only **1/16** re-enters the modal
fresh family of the named target. One additional node, GOLGA2, is reported separately as a
documented broad-family outgroup and remains in the broad RNA family. The distinction is
therefore:

- most suspect *nodes* are reproducible RNA locus constructions;
- most suspect old *memberships* do not reproduce under the current orientation guard; and
- one residual mixture remains an actionable subfamily-precision failure: GGO `MAGEB16` in
  the MAGEA component.

The former `GOLGA2/SWI5` false-positive call was a truth-granularity error. The fresh 26-exon
minus-strand prediction is a near-full GOLGA2 model; SWI5 is only a boundary-overlap label.
GOLGA2 has five redundant forward RNA edges to the chr15 GOLGA6/8 block, and published
core-duplicon work describes GOLGA6/8 as closest to the ancestral GOLGA2 locus. It is therefore
retained as `REVIEW_RELATED_OUTGROUP` in the broad family and excluded from the recent-copy
subfamily view. A direct RNA identity-0.80 run separates it, but also loses two named MAGEA
members and further splits NBPF, so that threshold is not promoted to a global rule. Detailed
evidence: [`bench/o1_golga2_subfamily_audit.md`](../bench/o1_golga2_subfamily_audit.md).
Production emission of the two hierarchy levels is deferred under the explicit implementation
contract in
[`o1_duplication_provenance_model.md`](o1_duplication_provenance_model.md); current Rustle still
emits only the broad family partition.

The NBPF adversarial case is especially informative. All nine conflicting TTC6/DNAH14/CEP152/
TTBK2/etc. loci that are recovered by the fresh run remain outside the fresh NBPF target family.
Thus the old 17-node NBPF repeat-bridge component is not reproduced even though its contaminant
loci themselves are real emissions.

The safety cost is visible too: only **69/75** named target nodes land in the modal fresh family.
Three named nodes are not re-emitted, while one each from GOLGA6/8, TBC1D3, and NBPF is emitted
into a sibling component. This is the precision/recall trade-off that must be reported rather
than hidden behind a pure graph.

The archived `copies`, `families`, `E_r` nodes/edges/rule certificates, logs, match table, and
fresh GFA graphs are in
[`bench/o1_fresh_emission_validation/`](../bench/o1_fresh_emission_validation/); the reproducer
is [`bench/o1_fresh_emission_validation.py`](../bench/o1_fresh_emission_validation.py). The GFAs
contain only actual newly emitted representatives and actual dumped `E_r` edges. Their colours
are post-hoc truth labels and do not influence graph construction.

This is still a regional falsification, not a whole-genome recall estimate. `samtools -M -L`
preserves records from the original genome-wide alignments but removes records outside the
predeclared ±10 kb intervals. It tests whether nodes and local groupings recur without feeding
the audit partition back into Rustle; it does not establish byte-identical whole-genome
partitioning or discover families outside the panel.

The HSA run certificate also reports **one co-located, same-strand pair decided by O2's
`reads_distinguish` predicate rather than by sequence**. Therefore the HSA node set in this
experiment is not a sequence-only O1 construction. This dependency is retained in the archived
log and must accompany these counts; it is evidence that the intended DNA/RNA synergy still
needs an explicit contract rather than being treated as an invisible implementation detail.

### Whole-catalog structural attack

Restricting the frozen all-vs-all PAFs to within-current-family edges gives:

| catalog | R0 edges | forward edges | edges removed | current families disconnected |
|---|---:|---:|---:|---:|
| GGO | 2,445 | 2,396 | 49 | 37 / 494 |
| HSA | 2,392 | 2,239 | 153 | 27 / 394 |

This is a blast-radius measure, not a false-negative count: the current family label is
the object under test, not truth. Among removed edges with independent labels, the
breakdown is:

| catalog | scored FP | scored TP | grey | unlabelled |
|---|---:|---:|---:|---:|
| GGO | 6 | 2 | 33 | 8 |
| HSA | 9 | 0 | 43 | 101 |

The large HSA unlabelled block is dominated by `GWFAM148`: 27 unannotated, mostly
two-exon nodes on acrocentric p-arms. Other large splits are the already named
repeat-bridge mixtures `GWFAM33` (TTC6/DNAH14), `GWFAM48` (TDRD5/AGAP), and
`GWFAM97` (ATRNL1/KIAA1328). The structural attack is therefore consistent with the
guard removing false families, but it cannot by itself prove that every split is
correct.

## Other candidate rules and their verdicts

### Reject alignment blocks that miss the predicted ORF — **not a hard rule**

The strongest earlier composite rejected a passing record when it was reverse, or
when identity was below 0.90 and less than 20% of the aligned block intersected the
longest forward ORF. On the frozen non-panel arms it looked excellent: 13/14 GGO
catalog FPs removed for 3/1,795 GGO catalog TPs.

That result did not test the seven named families; it **abstained on the entire panel**
because the old panel FASTAs were not orientation-normalized. On the corrected panel,
the ORF clause deletes the sole GSTM4–GSTM5 edge and annihilates the tested two-copy
GSTM family. It also deletes one additional NPIP edge. The apparent panel safety was
vacuous. Keep `orf_void` as a risk flag only.

### Junction-spanning homology — **confirmation flag, not membership**

On the curated-HSA arm, an R0 alignment crosses at least one transcript exon junction
on both sides for only 2/28 false edges, but also fails on 595/5,821 true edges
(10.2%). It is strong evidence when present but excludes real single-exon,
retrocopy, partial-assembly, and structurally divergent families. Emit
`junction_witness=true/false`; do not require it on every edge.

### Block promiscuity — **risk flag, not membership**

As documented in the coverage repair, asking how many catalog nodes match the same
block strongly identifies dispersed-repeat bridges. It is not pair-local: adding an
unrelated sequence to the catalog can change the decision for `(a,b)`, violating the
seed/library invariance expected of the O1 definition. Emit the count and a
`repeat_promiscuous` flag; do not use it to decide `E_r`.

### Coverage/length/soft-mask/read-depth retuning — **rejected**

Coverage-of-longer, absolute aligned-base floors, soft-mask thresholds, and
read-tiling/read-support thresholds all failed the same ordering test: a setting that
removes a useful number of repeat bridges breaks named families first or introduces
library dependence. Those negative results remain valid.

## Recommended O1 output contract

Strengthen the thesis claim by separating **membership** from **confidence**:

1. `E_r_plus`: R0 plus the transcript-orientation guard, only on typed RNA substrates.
2. `family`: the existing gamma-quasi-clique over `E_r_plus`, with at least two
   spatially distinct loci.
3. Per-edge evidence: orientation, identity, coverage, junction witness, ORF
   intersection, repeat-block promiscuity, and tier provenance.
4. Per-family confidence:
   - `supported`: at least one non-repeat structural witness and no sole repeat bridge;
   - `repeat_suspect`: connectivity depends on an ORF-void, junctionless, or
     promiscuous block;
   - `unresolved`: the family survives but lacks enough independent evidence.

Only `supported` families should feed the primary O1 precision headline. The other
families remain visible rather than being silently deleted, preserving recall and
making the limitation auditable.

The validation exporter now makes this separation physical as well as textual:
`<family>.gene_family.gfa` contains only admitted members and membership edges, while
`<family>.audit.gfa` retains the entire SD-discovery universe with dispositions and
`MB` membership tags. Consequently a Bandage rendering of the primary graph no longer
looks impure merely because review and rejected nodes were deliberately retained for
audit.

## SD-to-gene-family typing gate

The Soto benchmark is a catalog of segmental duplications. An `ID_*` match therefore
establishes homologous duplicated sequence, **not** membership in one multi-copy gene
family. O1 validation now uses the following asymmetric gate:

1. forward transcript homology may recruit an unnamed, expressed locus into the
   component of an independently named gene-family member;
2. DNA homology corroborates those edges and may rescue an RNA-null locus only when
   that locus has independent same-family annotation;
3. DNA-only SD homology cannot recruit an anonymous locus into a gene family;
4. a Soto-missing candidate is counted as a gene copy only when SD homology, at least
   two genomic loci, and independent same-family gene annotation agree.
5. forward RNA homology alone cannot relabel a locus independently annotated as a
   different named gene. Such a locus is `REVIEW_CONFLICTING_GENE`: visible in the
   audit graph but absent from the primary graph. `LOC*` placeholders do not trigger
   this rule, and an explicit same-family symbol overrides it.

Reproducer and Bandage-readable graphs:
[`bench/o1_gene_family_audit.py`](../bench/o1_gene_family_audit.py) and
[`bench/o1_gene_family_audit/`](../bench/o1_gene_family_audit/).

On the corrected seven-family panel, NPIP, TBC1D3, RABL2, and MAGEA pass as connected
typed primary graphs. APOBEC3H and GSTM3 remain named but disconnected and are
reported as false negatives, not removed. HERC2 has a connected ten-node candidate
graph but only one independently named member in this annotation, so it is reported
as annotation-limited rather than counted as a validated pass. Three adjudicated
contaminants are excluded from primary membership: a PKD1 interval, the reverse-only
PDXDC1 interval, and an RNA-null NPEPPS interval. A second, forward-matching PDXDC1
locus is withheld as `REVIEW_CONFLICTING_GENE`; RNA homology is sufficient to keep it
reviewable but not to relabel an independently named gene as NPIP. NPIP's primary
graph remains connected with all 19 independently named NPIP members. DNA-only
anonymous intervals remain `REVIEW`, not automatic false positives.

Applying the gate to the 18 previously reviewed Soto-missing loci supports four
unique annotated gene copies (NBPF8, NBPF19, GOLGA6L1, and putative ANKRD20A2/
LOC128966611). Three unannotated loci remain novel-copy candidates. CHRFAM7A is
explicitly rejected as an ULK4P-family assignment despite matching the ID_179 SD
block, and the single-copy chr7 locus is rejected from the multi-copy claim.

## Promotion gate

Before making `E_r_plus` the default:

1. ~~implement it as an RNA-only switch and include the effective value in the O1/O2
   rule certificates~~ — complete as `--rna-forward-only`;
2. run `gw_family_catalog` end to end on GGO and HSA with identical inputs, flag off
   and on;
3. re-score external TP/FP pairs and emitted family partitions, with special review of
   GGO `GWFAM103`;
4. require the named TBC1D3, RABL2, MAGEA, GSTM, and HERC2 panels to remain intact and
   treat the PDXDC1-labelled NPIP singleton separately;
5. ~~keep `--from-genome` unchanged and add a regression test proving an inverted DNA
   duplication is still accepted there~~ — complete in the shared genome-mode grouping
   core; the CLI also rejects the RNA-only flag in `--from-genome` mode.

---

# Appendix A — the guard's confirmation on the human negative panel

*Merged from `o1_false_positive_rules.md` on 2026-08-20. It is the same guard, measured; keeping
it as a separate file only made the evidence harder to find.*

## The false-merge rate, re-measured under the 2026-08-20 defaults

**Status 2026-08-20. Re-run of the frozen 150-window negative panel with the current binary.**
Harness: `/home/juanfra/winloci_scratch/o1neg/` (`run2.sh`, `score.py`, seed 101, panel unchanged).
July outputs preserved as `out.jul17` / `dump.jul17` and re-scored for the comparison.

## 0. ⚠ A provenance correction

`ONE_METHOD.md` was updated on 2026-08-19 to say the rates "were measured on the shipped 494-family
catalog". **That is wrong for the false-merge rate.** It is measured on **HUMAN CHM13 v2.0**, over 150
gene-tight single-locus windows drawn from 1,630 eligible, with **A119b** IsoSeq — not on the GGO
catalog at all. It is also a **specificity and a LOWER bound**, not a precision: the panel has no
positive stratum, so there is no prevalence and no precision to compute.

## 1. The headline is unchanged

| | July (pre-guard) | 2026-08-20 defaults |
|---|---:|---:|
| **windows emitting ≥1 family** | **2/150 = 0.0133** [0.0037, 0.0473] | **2/150 = 0.0133** [0.0037, 0.0473] |
| the two windows | W063 `ZNF492`, W106 `ANKHD1` | **the same two** |
| co-membership assertions | 4 | 4 |
| self-overlap (false by coordinates) | 1 | 1 |

**Everything shipped on 2026-08-19/20 — the transcript-orientation guard as default, one path with
refine rejected, substrate typing, λ recomputation, the streamed PAF — leaves the false-merge
specificity exactly where it was.** The definition's measured error rate is stable under all of it.

## 2. But the spurious-edge burden collapsed

Same panel, same scorer, `dump/*.edges.tsv`:

| | July | now | change |
|---|---:|---:|---|
| E_r edges emitted on the negative panel | **28** | **3** | **−89%** |
| edges between OVERLAPPING spans | **26** | **1** | **−96%** |
| self-identity CERTIFIED false | **7** | **1** | −86% |

This is the orientation guard doing exactly what the coordinate analysis predicted: overlapping
antisense pairs are ~46× enriched among minus-only edges, and on this human negative panel they were
**25 of the 28** edges. The guard removes them at the edge layer.

## 3. Why the family rate did not move — and why that is the right outcome

Those 26 spurious edges were already being stopped downstream: only 2 of 28 ever reached `copies.tsv`.
The guard cleaned a layer that was not the binding constraint on the emitted catalog.

And the two survivors were **never orientation artifacts**, so no version of this guard could have
fixed them:

* **W063 `ZNF492`** — self-identity: the aligned block IS the 1,204 bp genomic intersection of the two
  node spans, identity **exactly 1.000000**. That is census **pathology (a)**, one locus emitted as two.
* **W106 `ANKHD1`** — the 206 bp linking node is **100% soft-masked** interspersed repeat.

Both were predicted not to move: pathology (a) is a coordinate signature, not a strand one.

## 4. ⚠ A rate that rose while the burden fell

The **certified-false edge rate** went **7/28 = 0.2500 → 1/3 = 0.3333** while the absolute count fell
**7 → 1**. The denominator shrank faster than the numerator. Quoting the rate alone would say the edge
layer got *worse*; quoting the count alone would miss that what remains is proportionally more
concentrated. **Report both, or neither.**

## 5. What this does and does not establish

* ✅ The false-merge specificity is **stable at 1.33%** across a substantial change of defaults — the
  strongest evidence yet that it is a property of the definition and not of an invocation.
* ✅ The orientation guard's precision benefit is now **measured on an independent human negative
  panel**, not only on the GGO FP arm it was derived from: **28 → 3 edges**.
* ❌ It is **still a specificity and a lower bound**. No positive stratum, no prevalence, no precision.
* ❌ It does **not** transfer to gorilla. This panel is human CHM13/A119b. The GGO 627-family catalog
  has no comparable negative panel.

---

# Appendix B — the other two 2026-08-19 candidates

*Merged from `o1_genome_anchored_repeat_gate.md` and `o1_junction_crossing_guard.md` on
2026-08-20. All three candidates answer one question — can `E_r`'s precision be improved without
breaking P1? — and they are only comparable side by side. **Verdicts: the orientation guard SHIPPED
as the default; the genome-anchored veto is a FLAG (a veto, never an admission criterion); the
junction predicate is REFUTED as a gate** (12.80% of edges genome-wide, 100× exon-count bias).*

## A P1-safe repeat gate: fix the universe, not the statistic

**Status 2026-08-19.** Candidate derived and measured offline. **Not pipeline-confirmed (T8).**
Extends [`o1_coverage_repair.md`](o1_coverage_repair.md) §5 and
[`o1_false_positive_rules.md`](o1_false_positive_rules.md); it does not replace their results.

## 1. The defect being fixed

Candidate **R5** — block promiscuity — was the strongest discriminator ever found for O1's false
merges, and it was rejected because it **breaks P1**: MRPS17's block scores **50** partners over the
whole catalog and **1** from a 4-node seed, so the same pair is rejected run-whole and accepted
run-from-seed. `E_r` would stop being a relation between two sequences and become a function of the
node set.

The shipped repeat-hub gate has the same disease. `bench/vg_repeat_catalog.py` states its universe
outright: *"Universe = the 3462 gene-assigned loci that participate in E_r families."* It counts
**catalog members**. (It is confined to `family_define`/`driver.rs`; `gw_family_catalog` never calls
it, so the shipped O1 definition is not currently exposed.)

**The diagnosis was wrong in one word.** R5's disease is not *"it is a repeat count"* — it is
*"it is counted over the NODE SET."* Change the universe and the statistic becomes pair-local.

## 2. The statistic

Count occurrences in the **fixed reference assembly** instead of in the catalog. `GGO.fasta` is not
a function of the seed, so the value depends only on `(a, b, genome)`.

```text
block_a, block_b       the aligned interval on each rep, from the one passing E_r record
S                      canonical k-mers (k=21) present in BOTH blocks -- their shared anchors
g(x)                   occurrences of anchor x in GGO.fasta, both strands
min_shared_gmult(a,b)  = min over S of g(x)
```

Reading: *does the sequence these two loci share contain **any** anchor that is rare in the genome?*
A real paralogue pair shares private sequence, so some anchor is near-unique. A pair bridged only by
a mobile element shares nothing but high-multiplicity anchors.

Library-free: no RepeatMasker, Dfam or RepBase enters the definition. Strand-agnostic: the query set
carries each anchor and its reverse complement, so scanning the assembly's forward strand yields the
two-strand count. `k = 21` is odd, so no anchor is its own reverse complement and none is double
counted.

**Substrate:** the exon-sum spliced reps in `GGO_gwcat.copies.fa` (verified seqlen/genomic-span
0.076–0.447), never genomic spans.

## 3. It passes the test R5 failed

Both statistics computed on the **same pairs** under two universes — the whole 1,415-rep catalog, and
a seed containing only the pair's own family. R5-analogue = number of distinct reps in the universe
carrying any of the pair's shared anchors.

| statistic | pairs whose value CHANGES between whole-catalog and seed |
|---|---:|
| R5-analogue (counted over the node set) | **94 / 147 = 0.6395** |
| **genome-anchored `min_shared_gmult`** | **0 / 147 = 0.0000** |

Largest R5 swings, with the genome-anchored value beside them:

| pair | arm | R5 whole → seed | g-mult (both universes) |
|---|---|---:|---:|
| TP18481 | TP | 228 → 8 | 5 |
| FP00050 | FP | 189 → 7 | 12,663 |
| FP00051 | FP | 178 → 7 | 12,973 |
| FP00052 | FP | 136 → 7 | 8,952 |
| FP00049 | FP | 114 → 6 | 1,347 |

The FP swings are the GWFAM210 MRPS17 **AluY** hub — the mechanism the anti-FP characterisation
identified as Group 1. Under the genome universe those pairs score 1,347–12,973 regardless of what
else is in the catalog.

## 4. Discrimination on the frozen arms

Unit = **pair**. GGO only. FP arm = the 14 gorilla catalog false merges; TP arm = the 150 true pairs.

| arm | n scored | median `min_shared_gmult` | range |
|---|---:|---:|---|
| FP | 12 | **182** | 1 – 12,973 |
| TP | 135 | **2** | 0 – 44 |

**AUC (FP scores higher than TP) = 0.9429.**

| cut `M` | FP rejected | TP lost | TP cost |
|---:|---:|---:|---:|
| 10 | 11/12 | 7/135 | 0.0519 |
| 20 | 10/12 | 2/135 | 0.0148 |
| **50** | **10/12** | **0/135** | **0.0000** |
| 100 | 9/12 | 0/135 | 0.0000 |
| 500 | 4/12 | 0/135 | 0.0000 |

At `M = 50` the rule rejects **10 of 12 scored false merges at zero cost to 135 true pairs**. For
comparison on the same FP arm: coverage-of-longer at 0.20 rejected 1/14, R2@0.05 rejected 2/74
overall, and the transcript-orientation guard rejects 6/14 GGO catalog FPs.

## 5. Where it abstains, and why that is honest

The statistic requires at least one **exact** shared 21-mer. It has none on **11 of 158 pairs = 7.0%**
— 2 FP and 9 TP — and every one of them sits at identity **0.6927–0.8031**, where exact 21-mers are
expected to be rare (0.78²¹ ≈ 0.005 per position). Abstention is not biased toward one arm.

Those pairs must be reported `GMULT_UNMEASURED` and fall through to the incumbent rule. Absence of a
shared anchor is missing data, never evidence of a repeat bridge.

## 6. Is it just softmask restated?

Largely the same signal, but not identical, and library-free.

| | count |
|---|---:|
| FPs with softmask ≥ 0.70 on BOTH sides | 10/12 |
| FPs rejected by `min_shared_gmult ≥ 50` | 10/12 |
| rejected by g-mult, MISSED by softmask | **1** — FP00048, softmask 0.689 / 0.643, g-mult 76 |
| caught by softmask, missed by g-mult | **1** — FP00055, softmask 0.929 / 0.735, g-mult 19 |

So it is neither a superset nor a restatement of the softmask gate. Its value over softmask is that
it takes no repeat library into the definition.

**FP00058 scores 1 and is correctly NOT rejected.** That case is the LAGE3 processed pseudogene
against its own parent — a **truth-label failure**, not a false merge. A statistic that rejected it
would be wrong.

## 7. What is NOT established

1. **T8.** Offline re-derivation. Nothing here has been through the shipped binary. The E_r rule is
   unchanged and no default moves on this evidence.
2. **The AUC is not held out.** The mechanism of these same 14 FPs was characterised as
   repeat-driven (Group 1 AluY hub, Groups 2–3 low-copy elements), so scoring them on repeat content
   is partly circular. **The TP half — 0/135 at `M = 50` — is the load-bearing number**, because the
   TP arm was not selected for repeat content. A held-out FP set is required before quoting the FP
   rejection rate as a rate.
3. **Coverage is 12/14 FP and 135/150 TP.** Six TP pairs have no single passing record under this
   rule and eleven pairs abstain; both are stated above rather than dropped silently.
4. `M = 50` was read off this table. It must be fixed on a held-out set, or the rule quoted as an
   ordering rather than a threshold.
5. **Gate or flag is undecided.** As a gate it changes `E_r` and every downstream certificate; as a
   flag it costs nothing and P1 is untouched either way.

## 8. Can it REPLACE coverage? No — refuted

The tempting move is to drop the coverage clause entirely. Coverage is scale-free (a ~1 kb Alu is
≥0.50 of any node under 2 kb), which is the whole named hole; a rare-anchor test is structurally
immune to that, because a repeat has no rare anchor at **any** node size. So:

```text
E_r_free(a,b)  iff  exists a record with identity >= tier floor
                    AND min_shared_gmult(record) < M
```

**This cannot be judged on the FP/TP arms** — both were built from pairs that already pass
coverage ≥ 0.50. A rule that drops coverage must be judged on what coverage is currently holding
back. All-vs-all over the 1,415 shipped GGO reps, both tiers, every identity-passing pair split by
whether it clears coverage:

| | pairs | within shipped family | cross-family |
|---|---:|---:|---:|
| `COV_PASS` (the shipped edge set) | 2,727 | 2,474 | 253 |
| `COV_FAIL` (what coverage rejects) | 14,111 | 830 | 13,281 |

(The 253 cross-family `COV_PASS` pairs are not an anomaly: families are γ-quasi-cliques *of* E_r
components, so a cross-family E_r edge is by definition an edge γ cut.)

| `M` | shipped edges kept | NEW edges from `COV_FAIL` | of which cross-family |
|---:|---:|---:|---:|
| 2 | 835/2,727 = 0.306 | 131 | 23 |
| 3 | 1,415/2,727 = 0.519 | 437 | 214 |
| 5 | 1,822/2,727 = 0.668 | 878 | 528 |
| 10 | 2,173/2,727 = 0.797 | 1,463 | 967 |
| 20 | 2,363/2,727 = 0.867 | 1,915 | 1,310 |
| 50 | 2,459/2,727 = 0.902 | 2,239 | 1,567 |

**There is no operating point.** Holding new cross-family edges at parity with γ's existing 253
puts `M` near 3, which discards **48% of the shipped edge set**. Recovering 90% of the shipped
edges costs **1,567 new cross-family edges, 6.2× the 253 that exist**. Recall loss starts before
merge suppression does — the same shape as the coverage-repair impossibility argument.

**Why, mechanistically.** The TP distribution has median `min_shared_gmult` = 2, so an admission
criterion strict enough to exclude repeats (`< 2`, i.e. a genuinely unique anchor) also excludes
most real paralogue pairs. The statistic's discriminative power lives at the **top** of its range —
it separates "definitely a mobile element" (100–13,000) from everything else, and does not separate
within the low end at all.

⟹ **The genome-anchored statistic is a VETO, never an admission criterion.** It belongs on top of
the coverage clause, not in place of it. Coverage stays.

## 9. Reproduction

```bash
cd /mnt/linuxdisk/home/juanfraitu/o1_gmult
python3 blocks.py     # recover each pair's aligned block from the shipped reps
python3 gmult.py      # one streaming pass over GGO.fasta, ~9 min, ~3 GB RSS
python3 eval.py       # seed-invariance, discrimination, softmask additivity
python3 covfree.py    # §8: can it replace coverage? (all-vs-all + one genome pass)
```

Outputs: `pair_blocks.tsv`/`.fa`, `gmult.tsv`, `seed_invariance.tsv`, `covfree.tsv`.

## A threshold-free edge predicate that works at n = 2

**Status 2026-08-19. ⚠⚠ VERDICT: DO NOT ADD TO THE DEFINITION — see §3a.** The frozen arms
materially understated the cost; measured genome-wide the rule rejects **12.80%** of shipped `E_r`
edges with a **monotone bias against low-exon models**. It remains a legitimate **flag**.
Offline (T8), nothing through the shipped binary, no default moved.
Companion to Appendix B below (the genome-anchored veto) and
[`o1_error_case_census.md`](o1_error_case_census.md).

## 1. The problem this addresses

At **n = 2** — **348 of 494 = 70.45%** of the gorilla catalog — the γ-quasi-clique machinery is inert
and the entire definition reduces to **one coverage number**. That number is provably non-separating:
the accepted true pair GFPT1×GFPT2 scores **0.5353** while the rejected false pair ATP1A1×ATP4A scores
**0.5689**. No threshold on it orders the classes correctly.

The census prescribed the shape of any fix: **change the denominator or the substrate, not the
threshold.** This is the first candidate that does.

## 2. The predicate

> **The homology must cross a splice junction.** Reject an edge iff the passing alignment lies
> entirely within a **single exon on both sides**.

`max_exon_frac(side) = max over exons of (alignment ∩ exon) / alignment length`; reject iff it is
**1.0 on both sides**.

**Why it escapes the named hole: it has no length denominator at all.** The scale-free defect exists
because a ~1 kb repeat is ≥ 0.50 of any node under 2 kb. A repeat confined to one exon spans one exon
**at every node size**. The statistic is structurally immune to the defect rather than tuned against it.

**It is threshold-free** — "crosses a junction" is discrete. No fitted number enters.
**It is pair-local**, depending only on the two nodes and their exon structures ⟹ **P1 is untouched**.
**It is an edge predicate**, so it works identically at n = 2, where no graph structure exists.

## 3. Result on the frozen arms (unit = pair, GGO)

| | FP rejected / 14 | TP lost / 150 |
|---|---:|---:|
| coverage-of-longer @ 0.20 | 1 | damages many |
| transcript-orientation guard | 6 | 4 / 9,032 edges |
| genome-anchored repeat veto @ M=50 | 10 (of 12 scored) | **0** |
| **junction-crossing** | **12** | **9 = 0.0600** |
| **union of the last two** | **13** | 9 = 0.0600 |

The single FP the union misses is **FP00058**, the LAGE3 processed pseudogene against its own parent —
a **truth-label failure, not a false merge**, so it *should* survive. Effectively **13/13 real FPs.**

The two guards are complementary, not redundant: junction-only 3, gmult-only 1, both 9.
The genome-anchored veto abstains on 2 FP and 15 TP pairs (no shared 21-mer at identity 0.69–0.80).

## 3a. ⚠⚠ GENOME-WIDE: the arms understated the cost, and the rule does not generalise

The 164-pair arms are 6% of the catalog's 2,727 shipped edges. Measured on **all** of them — with
exon boundaries recovered by splice-aligning every rep back to the genome (**control: recovered exon
count equals `n_exon` on 1,405/1,415 = 0.9929**):

| rule | edges rejected / 2,727 | rate | within a shipped family |
|---|---:|---:|---:|
| **junction-crossing** | **349** | **0.1280** | **268** |
| genome-anchored veto @ M=50 | 100 | 0.0367 | 61 |
| genome-anchored veto @ M=100 | 77 | 0.0282 | 50 |

On the arms these cost 6.00% and **0%** respectively. Genome-wide they cost **12.80%** and **3.67%**.
⚠ **Both were measured only where the FPs were, and both look far worse on the full distribution.**
Since the measured false-merge rate is ~1.33%, **a rule rejecting 12.80% of edges is rejecting mostly
collateral, not false merges.**

### The pseudogene / retrocopy exposure is real and monotone

Stratified by the **smaller** model's exon count:

| `n_exon` (min) | edges | rejected | rate |
|---:|---:|---:|---:|
| 2 | 211 | 75 | **0.3555** |
| 3 | 231 | 45 | 0.1948 |
| 4 | 1,373 | 186 | 0.1355 |
| 5–6 | 356 | 41 | 0.1152 |
| > 6 | 556 | 2 | **0.0036** |

**A 100-fold gradient.** The mechanism is structural: a 2-exon model has exactly **one** junction, so
the alignment must hit that specific junction or be rejected. Low-exon models are disadvantaged by
construction — and that is where retrocopies, processed pseudogenes and compact genes live.

⚠ **The intended pseudogene fix — "abstain when a model is intronless" — is INERT.** The shipped
catalog has **0/1415 single-exon nodes** (the node builder requires ≥2 exons), so nothing abstains.
The damage sits at 2–4 exons, which *do* have junctions, just few. **No abstention rule rescues it.**

### Why this disqualifies it as a membership condition

It would make `E_r` a function of **how many exons a model happens to have** — a property of the
node builder and the assembly, not of the homology relation between two sequences. That is the same
genus of defect as R5's non-locality: membership must not depend on a property of the model rather
than of the pair. **Flag, never gate.**

## 4. The cost is real — 9 true pairs, and they are not an artifact

⚠ **First reported here as "all 9 are single-exon genes, so the cost is zero". That was WRONG**, and the
error is worth recording: the TP arm's `a_nex` column is **exons TOUCHED**, not total exons (that is
`a_tot_ex`), and touched-count is 1 **by construction** whenever `max_exon_frac = 1.0`. The check was
circular. With the correct column there are **0 single-exon nodes in either arm** — consistent with the
shipped catalog having **0/1415** single-exon nodes.

The 9 lost pairs are genuine: multi-exon nodes (2–7 exons) at identity 0.7117–0.7617 and coverage
0.5132–0.6586 whose homology genuinely sits inside one exon on both sides. That is a real biological
class — paralogues sharing a single domain-bearing exon — and rejecting them is a **recall cost, not
noise**. No scoping removes it.

**T13 verified:** both arms compute `max_exon_frac` by the identical formula
(`max(exon overlap) / alignment length`), from two different generators.

## 5. What is NOT established

1. **T8.** Offline. Nothing through the shipped binary; `E_r` is unchanged.
2. ⚠ **The FP side is not held out.** This discriminator was *found* on these same 14 FPs (the
   "≥99% of the alignment inside one exon on both sides, FP 13/14 vs TP 9/150" characterisation), so
   **12/14 is a description of a known set, not a rate**. **The load-bearing number is the TP side,
   9/150**, because the TP arm was never selected for exon structure.
3. ⚠ **Partly entailed by the truth predicate.** True families are annotated multi-exon protein-coding
   genes, so "shares more than one exon" is somewhat implied by being a real family. This does not
   dissolve the result — the rule rejects 12 pairs that the coverage clause *accepted* — but it means
   the FP rate needs an independent FP set.
4. **Exon structure is read-derived**, so "exon" is a model feature. Not more circular than the
   incumbent (the reps are exon-sums), but not independent evidence either.
5. **Gate vs flag undecided**, and the choice is a real one:
   - **gmult as a gate, junction as a flag** — 10/14 rejected at **zero** TP cost, 3 more flagged.
     The conservative ship.
   - **union as a gate** — 13/14 at 6.00% TP cost. A deliberate recall-for-precision trade.
   - **both as flags** — 13/14 surfaced, flag precision 13/22 ≈ 0.59.

## 6. What this changes about O1's defensibility

At n = 2 the definition was one non-separating number. It can become **one coverage number plus two
pair-local structural predicates**, one of them threshold-free and immune to the scale-free defect by
construction. That does not repair the coverage clause — the named hole and its 8.30% exposure ceiling
stand — but it means the *edge* carries structural evidence rather than a single scalar, and it does so
at exactly the family size where the graph machinery offers nothing.

## 7. Reproduction

Columns are already in the frozen arms: `fp14_detail.tsv` (`a_maxexonfrac`/`b_maxexonfrac`,
`a_nexon` = total exons) and `tp150_detail.tsv` (`a_mx`/`b_mx`, `a_tot_ex` = total exons,
⚠ `a_nex` = exons **touched**). Generators: `o1_antifp/analyze.py` and `o1_antifp/tpstats.py`.
