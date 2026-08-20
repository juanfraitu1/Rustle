# O1 false-positive hardening: rules that survived falsification

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
[`o1_hierarchical_family_followup.md`](o1_hierarchical_family_followup.md); current Rustle still
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
