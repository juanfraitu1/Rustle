# REGISTER OF NEGATIVE RESULTS

*Last updated 2026-08-14. 641 adversarially verified negative results, consolidated into 631 entries
(near-duplicate results recorded in two source files are merged into one row). Unsupported claims
were dropped before this document was written; every entry below carries the number that killed it.*

## How to use this

1. **Read this before proposing an approach, not after.** Search the claim you are about to make.
   Entries are written as the ORIGINAL claim, in the words it was first believed in, so you recognise
   your own idea in it.
2. **Check the verdict kind before arguing.** REFUTED/RETRACTED means it was believed and is wrong.
   DEAD-END means it was never wrong, it just cannot be made to work. TRAP/TAUTOLOGICAL means the
   MEASUREMENT was invalid, so the result meant nothing either way. SUPERSEDED means it was right and
   a better number replaced it — quote the replacement.
3. **A re-run needs a new mechanism, not a new mood.** If you want to revisit a DEAD-END, state which
   of its killing numbers you expect to move and why.
4. **`redo_risk` (column R: H/M/L) ranks how appealing the dead idea still is.** High-risk entries are
   sorted first inside every section; they are the ones that get proposed again.
5. **Read §2 METRIC TRAPS first if you are designing an evaluation.** Those patterns generalise past
   the specific claim that died and have each killed several results.

## §1 Contents

| § | Section | Entries |
|---|---|---|
| 2 | Metric traps (rules) | 20 rules |
| 3 | The false-omission measurement (2026-08-14) | new |
| 4 | O1 — family definition | 232 |
| 5 | O2 — copy assignment | 93 |
| 6 | O3 — reference-absent / collapsed copies | 71 |
| 7 | Benchmarking & evaluation | 131 |
| 8 | Infrastructure & tooling | 67 |
| 9 | Framing & scope | 37 |

Counts by verdict: REFUTED/RETRACTED **324** · TRAP/TAUTOLOGICAL **145** · DEAD-END **112** ·
SUPERSEDED **50**. By redo risk: high **309** · medium **260** · low **62**.

---

# §2 METRIC TRAPS

Rules, not anecdotes. Each has killed more than one result.

| # | Rule | Origin cost |
|---|---|---|
| T1 | **Never condition the denominator on the prediction.** Scoring only units the method admitted means abstention can never lower the rate. | 7 kills; O3 rate 37.14% → 31.09%, O2 "unique-mapper agreement 99.9%" |
| T2 | **A prediction that is a subset of its own truth is tautological.** If truth is built from the method's own output, recall cannot fall. | 4 kills; showcase recall 0.786 → 0.524, boundary in-band 0.796 → LOO 0.5227 |
| T3 | **Never select WHICH component/unit to score.** Oracle selection ("the best NPIP-containing component of 377–510") is unobtainable blind. | purity 0.237 and 13.6% permanently withdrawn |
| T4 | **Ratio-to-truth ⟹ report the IN-BAND FRACTION, never the median.** A median improvement can be pure distribution shift. | any-locus bounding: median 0.55→1.03 but in-band 47%→37% |
| T5 | **Never score "bases explained".** It rewards UNSPLICED models: one giant exon contains every aligned base. | baseline 84.8% vs co-threaded 68.3% was an artifact of 76% stubs |
| T6 | **An edge-count-matched null proves nothing.** Match the SIZE distribution, and respect genomic position. | zero-information RANDCUT scored 152/251 vs Rule 3's 178; position-blind null understated difficulty 124× |
| T7 | **Never judge a change to what a NODE IS on node-level or pair-level metrics.** | 3 predictions passed offline and failed e2e (0.553→0.433; +0.023→−0.072; 25%→51% in-band → F1 0.704→0.401) |
| T8 | **Offline/proxy re-derivation is a hypothesis generator, never a test.** It omits the partition, the ≥2-loci gate, and the emission rules. | 4 offline-proxy errors; replay_er predicted F1 0.745, pipeline 0.595 |
| T9 | **Recall-only evaluation of an edge-ADDING change is not an evaluation.** Recall is partition-blind and hides over-merge. | PCDHB recall stayed 14/14 while pure families fell 12→9 |
| T10 | **Median-ratio is null-informative.** A constant predictor scores 1.0000 on it. | famCN "validated at 0.997"; true integer match 34.8% |
| T11 | **Check the panel CAN fail before quoting its pass.** γ=0.20 can only bind at n≥11 on a panel whose largest component is n=10; a complete graph has zero redundancy power. | γ "invariance" vacuous; 27-node panel null power 0 |
| T12 | **State the unit.** Locus vs family vs pair vs node vs molecule changes p-values and rates. | p=0.0312 (locus) → 0.125 (family); blind "6 vs 6" was locus-pairs vs node-pairs |
| T13 | **Numbers from different tiers/presets are not comparable.** Preset alone moves identity by up to 0.17. | GWFAM227 0.9629 → 0.8766 on byte-identical sequence |
| T14 | **Filter `-F 2308` before ANY per-read statistic.** Secondaries are clipped and MAPQ-0 by construction. | 47.7% clip rate → 1.1% (~40×); "70–100% MAPQ0 mangled" retracted; junction counts 154→56 |
| T15 | **Circular truth: never consume the comparator's own files and call it replication.** 3 of 4 Soto legs do. Say "concordance". | 100% precision → 76.2%/91.4%; ARI 0.720 → 0.514 |
| T16 | **Agreement is not accuracy.** The silver standard is minimap2's own placement, and it is empty in the hard regime. | silver 456/456 = 100% covers only 24.3% of reads |
| T17 | **Aggregate histograms mask per-copy redistribution.** Always check per-copy. | GSTM support 1196/1255 → 242/2209 with copy calls unchanged |
| T18 | **An oracle prize is a CEILING, not a realizable gain.** Verify the killers survive in your own population. | averted a −114 catastrophe; corrected mechanism attribution 4× |
| T19 | **Argmax/best-hit counting is not counting.** Count qualifying edges. | best-hit said 6/19 orthologs, edges say 19/19 (elsewhere it INFLATED 1.7×) |
| T20 | **A control that removes nothing has excluded nothing.** | the 475-ASJ masquerade control removed 0/475 and cannot be cited |

---

# §3 THE FALSE-OMISSION MEASUREMENT (2026-08-14)

Answer to "does the 53.8% of absorbed orphans landing on a paralog the family lacks mean O1 missed
53.8% of its members?" — **No. (b), with a bounded minority of (a). The 53.8% must be RESTATED, not
quoted as recall failure.**

Artifacts: `/home/juanfra/winloci_scratch/o1_falseomission/` — `RESULTS.md`,
`per_family_er_ALL.tsv` (728 tests), `SUMMARY.txt`, `nonsister_calls.json`, `testis_depth.json`.

### 3.1 The sanity control saved the result

The brief's assumed substrate was wrong. Shipped E_r nodes are **exon-sum spliced representatives**
(`GGO_gwcat.log`: "exon-sum +sensitive refine"; `GGO_gwcat.copies.fa` median 2,982 bp), **not**
genomic intervals.

| control | reproduces shipped E_r |
|---|---|
| genomic-vs-genomic | 1/3 |
| fibroblast-rebuilt | 3/15, 4/15 |
| rep-vs-genomic | 24/162 = 0.1481 |
| **rep-vs-rep (gold geometry)** | **103/104 = 0.9904 [0.9475, 0.9983]** |

Run as literally specified (genomic intervals + read-span fallback) the experiment returns **43%
false omission** — pure artefact: destination spans are read-*start* spans (one is 135 bp) and
`denom="min"` is structurally blind to truncation. **That number is retracted by its own control.**

### 3.2 Population

162 panel families (design, frozen before masking) · absorbed 104/162 = 0.6420 [0.5656, 0.7117] ·
non-sister 56/104 = 0.5385 [0.4430, 0.6312].

### 3.3 ARM 1 — is the destination even a node?

Non-sister destinations that are nodes in O1's catalog at all: **2/56 = 0.0357 [0.0098, 0.1212]**.
**No node at all: 54/56 = 0.9643 [0.8788, 0.9902].** Of the 2 scoreable, edge forms **0/2** — both
`no_record`, not near-misses.

### 3.4 The tissue excuse is REFUTED

Catalog substrate is `GGO_ds.bam` (OR6737 testis, written 26 min before the catalog); excision reads
are fibroblast. But at the 54 no-node destinations, ≥2 testis reads (the catalog's own minimum
`n_reads`) are present in **44/54 = 0.8148 [0.6916, 0.8962]**, median testis depth 36.5.
**O1's own reads saw these loci.**

### 3.5 The fork, bracketed (side A = the exact shipped node in both arms)

| arm | control | non-sister edge forms |
|---|---|---|
| ARM 2 migrated-chain — exon structure forced to match ⟹ **pro-edge biased, UPPER BOUND** | 46/48 = 0.9583 | **39/56 = 0.6964 [0.5666, 0.8010]** |
| ARM 3 testis-rebuilt — independent reads, **unbiased** | 89/104 = 0.8558 | **9/56 = 0.1607 [0.0869, 0.2781]**, ceiling-corrected 0.1878 |
| ARM 3B window variant | 83/104 = 0.7981 | 10/56 = 0.1786 [0.1000, 0.2984] |

ARM3/ARM3B agree 51/56; union 12, intersection 7. Only **GWFAM79** fires in both ARM 2 and ARM 3 —
the arms measure different objects, which is why one is a bound. **The gap is the finding:** the
destination's DNA is homologous to the masked copy *over the masked copy's own exons* (ARM 2,
identity 1.0000), but the destination's *own transcript* is not (ARM 3). **Reads migrate into a
paralogous SEGMENT, not a paralogous GENE.**

### 3.6 Clause census, all 728 tests

EDGE 488 · coverage 27 · no_record 109 · no-node 54 · rep-unbuildable 50 · **identity-alone 0 ·
both-fail 0 · split 0**. **Identity never fails — 0/728**, stronger than the project's prior (sixth
independent substrate). Coverage failures are not marginal: identity-passing records reach coverage
median **0.219** (p25 0.169, max 0.463) against a floor of 0.50.

### 3.7 Rate, both denominators

PRIMARY denominator is the **162 panel families** (frozen before masking, part of the design); 104
and 56 are experiment OUTCOMES and therefore secondary (see T1).

| | /162 PRIMARY | /104 secondary | /56 secondary |
|---|---|---|---|
| **(a) ARM 3, unbiased** | **9/162 = 0.0556 [0.0295, 0.1022]** | 0.0865 [0.0462, 0.1563] | 0.1607 [0.0869, 0.2781] |
| (a) ARM 2, upper bound | 0.2407 [0.1814, 0.3121] | 0.3750 [0.2880, 0.4709] | 0.6964 [0.5666, 0.8010] |
| **(b) no edge** | | | **47/56 = 0.8393 [0.7219, 0.9131]** |

**Quote 9/162 = 5.6% [3.0, 10.2] as O1's false-omission rate**, with 24.1% [18.1, 31.2] as the
pro-edge upper bound.

Blunt version: *"53.8% of absorbed families sent orphans to a locus the family lacks"* is true about
where reads go. It is **not** evidence that O1 missed 53.8% of its members. It is a statement about
**E_r's coverage clause and O1's node set** — that is the defensible finding.

### 3.8 Caveats that could move it

* ARM 3's control ceiling is 0.8558, and 29/56 destinations built no testis rep (counted as no-edge)
  ⟹ 9/56 is a **floor within ARM 3**; the conditioned 9/27 is deliberately NOT quoted (T1).
* Destination windows are read-defined — ARM 3B's wider window flipped 5/56.
* ARM 1's control is 103/104 (GWFAM406, identity 0.99 / coverage 0.4346): a pairwise re-run of the
  shipped rule does not perfectly reproduce O1's shipped partition, so ~1% slack is inherited.
* `run.py` is the disqualified v1, kept only as the record of the substrate error.

---
# §4 O1 — FAMILY DEFINITION

## 4.1 REFUTED / RETRACTED — was believed, is wrong

### High redo risk

| R | Claim | Killing number | Why |
|---|---|---|---|
| H | γ is inert on real families because real families are complete graphs | only 2/7 complete, both n≤3; densities 0.348–1.000, median 0.864 | γ=0.20 is inert for a different reason (0.148 below HERC2); HERC2 SPLITS at γ=0.40 |
| H | A family is a quasi-clique on DNA and connected components on RNA (DNA tight, RNA loose) | RNA's largest component is DENSER in 4/7: NPIP 0.897>0.871, TBC1D3 1.000>0.864 | INVERTED, and only because RNA sheds nodes so the survivor is denser by construction |
| H | The RNA family is a k-core of the DNA family | pooled AUC 0.7668 p=0.0086 → within-family normalised degree: 2/5, sign test p=0.8125 | The first test was my own Simpson's paradox — pooled degree measured FAMILY SIZE |
| H | chr1 is 80% clustered vs 36.3% genome-wide, so chr1 flatters the method | matched: 16/40 = 40.0% vs 157/432 = 36.3%, Fisher p=0.6090 | Two different rules; the chromosome clause is unfalsifiable in chr1 scope. **chr1 IS representative** |
| H | Dispersal is the variable that determines whether O1 recovers a family | matched protein identity <0.40: clustered 1/37 vs dispersed 0/66, Fisher p=0.3714 | DIVERGENCE is the variable. Real split: median within-family protein id ≥0.50 → 9/14 vs 1/26 |
| H | Dispersed families are where the identity floor finally binds | 194/194 records with identity ≥0.60 over 4,287 within-family pairs | Identity never fails — sixth independent substrate. Only 16/4,287 = 0.37% would be an edge |
| H | DNA 90.3% / RNA 84.5% / union 97.0% with 24 RNA-only members — RNA rescues divergent copies | re-run at shipped tier: 298 / 8 / 51 / 5 ⟹ **quote 8, not 24**; O1 alone is 65.5% not 84.5% | The RNA column was never re-derived — it is a 07-25 file on a `-M -L` SUBSET BAM |
| H | The RNA-only members are the DIVERGENT copies (median identity 88%) | 0.8770 under nmatch/blocklen vs **0.9661 under 1−de**; this alone moves 17 of 24 to "both" | Artifact of the banned identity metric, which charges indels against identity |
| H | The shipped tier can only ADD records, so DNA coverage can only rise | DNA FALLS 351→349, RNA-only RISES 6→8 | `-X --no-long-join` is not a superset: shorter spans drop pairs below the 0.50 coverage floor |
| H | The extra 91 edges come from jointness (RNA supplying what DNA cannot see) | DNA-only window cut to len(RNA node): 351 edges, symmetric difference **0**, recovers 91/91 | The gain is CONTIGUITY, and DNA has it free. Node LENGTH contributes 0 of the 91 |
| H | Blind (seed-free) mode can define families without an annotation seed | blind 18/40 vs all-singletons 17/40; m≥3 blind 1/22 vs SEEDED 5/22 | Blind FINDS every locus and DELINEATES none; beats a degenerate by one family |
| H | Read-through transcripts fuse blind nodes — filter them and blind will delineate | removes **0 of 6,824** block-links; node set BIT-IDENTICAL | Its 6/251 score is an identity, not a measurement |
| H | Node EXTENT is the whole story — get extent right and edges follow | survival by extent 0.167→6, 0.411→100, 0.522→152, 0.734→178, **0.962→102** | Non-monotone with an interior optimum; near-perfect extent scores worse |
| H | The pipeline's output is a non-trivial function of the minimal annotation given | on 26 of 39 chr1 families every seed returns only itself; 30/40 best component matches ≤1 gene | Advisor criterion (c) FAILS — the method is INACTIVE on three quarters of the substrate |
| H | The DNA-vs-RNA coverage gap is caused by LENGTH, not splicing | union lifts coverage above the best single record for 150/392 pairs; pool+union closes 28.3 pp → 0.68 pp | Two measured causes replace it: FRAGMENTATION and SUBSTRATE CHOICE (one arbitrary isoform) |
| H | Use coverage of the LONGER sequence to equalise the DNA/RNA scale difference | gap narrows to 5.3 pp but **DNA recall 0.930→0.492** | It equalises by DESTROYING |
| H | Use the genomic projection of the RNA node as the definitional substrate | a true DNA edge's covered bases are only 13.2% exonic; 134/440 true edges <10% exonic | BLIND to transcription (13.2% = the genome-wide exon fraction) and invents 11 cross-family edges |
| H | We need BOTH DNA and RNA because RNA carries copies the genome lacks | 44/44 cells joint == DNA-only; 0 RNA-only edges in 30/30 degraded cells | RNA node bases are fetched from the assembly, and degraded copies' reads relocate 81–99% |
| H | β (seed extension) is INERT — it changes nodes but not the partition | nodes 73→69→62→52→45 at β=0/5k/10k/20k/30k; β=10 kb SPLITS APOBEC3 (cut vertex) | Measured over a two-point window where β is dominated by two undocumented constants |
| H | γ is VACUOUS because real families are complete graphs | genome-wide gorilla: 125 γ-eligible components, only 47 (38%) complete; one 221-node component at density 0.0686 fusing 24 GWFAM ids | γ's real job is HAIRBALL SUPPRESSION at genome scale — keep it |
| H | c = 0.50 sits on a plateau, so it is derived rather than tuned | external panel recall 13,10,8,7,6,6,4,4,4 strictly monotone; bootstrap P(max gap ≥0.0707)=0.9994 | No plateau exists; c is a tunable OPERATING POINT (and c=0.48 changes the catalog) |
| H | Change the coverage denominator min()→max() to stop over-merges | GSTM 5 copies → **0 families**; catalog-wide 84 fam/381 copies → 72/281 (26% of copies lost) | Removes only 2 of 4 artifacts and destroys a real family; suite stayed green |
| H | minigraph shows the 27 NPIP copies collapse to one 16,118 bp path | minigraph reported "inserted 0 events" | The path was the REFERENCE echoed back; "no variation between copies" described a 1-copy graph |
| H | The DNA∪RNA union buys extra family connectivity | 32 RNA-only edges, **0** joining two DNA components; 0 families decided | Insurance that only pays out under tightening (strictly better in 4 swept cells) |
| H | Use the margin ratio (density over γ) as the family robustness statistic | NPIP 13.07×, TBC1D3 5.18×, GSTM 1.00× all have ZERO cut edges; HERC2 at 2.82× has one | Margin does not predict fragility — 2-edge-connectivity does. Use CUT EDGES |
| H | Cluster stability across k certifies a subfamily | the most k-stable cluster (S3) is a 2,048 bp node with 0 exon blocks and coverage exactly 1.0000 | Stability alone cannot certify anything; a 2 kb piece aligns fully into a 12.7 kb copy |
| H | NPIP splits into several families at the RNA level | full BAM: 61 reps ⟹ 2 families, GWFAM0 = 53 loci with all 19 NPIP genes | The premise was the SUBSET-BAM trap (third hit): NPIPB3 0 reads vs 1,445 in the full BAM |
| H | RNA components REFINE DNA families, so RNA can only under-merge | 31/210 pairs same-RNA-family but different-DNA-family; **0** the reverse | Measured backwards: DNA refines RNA, RNA OVER-merges through its nodes |
| H | Seed with the annotated mRNA transcript instead of the genomic span | human 19/19 → **9/19**, component splits; gorilla 25 → 14 | Kills recall; looks good only on the panel that has no ground truth |
| H | Use coding fraction (mRNA coverage) as a per-locus certificate | annotated NPIP loci median 0.080 vs NON-annotated 0.443 | Runs BACKWARDS. Also retracts "only 10/25 gorilla loci carry coding" |
| H | Seven of the 19 NPIP genes are lost by the pipeline (chain shattering) | NPIPB3 0 reads in the subset BAM vs 1,445 in the full BAM; families 13→4 | Read counts from the FULL BAM, copies from a `-L` SUBSET catalog. All 17 genes get reps |
| H | Use minimizer-Jaccard similarity to group loci into families | domain-sharer CREB1~METTL21A J=0.406 > every true paralog (max 0.313) | The domain-sharer range ENCLOSES the paralog range |
| H | Use annotated gene coordinates as the vertex set of the read-conflict graph | annotated MAGEA members have de = 0 between them; real conflicts are 24–207 kb away | Vertices must be de-novo expressed loci |
| H | Use naive all-column POA reciprocal coverage as the family criterion | true copies [0.471,0.842] OVERLAPS domain-sharers [0.268,0.843] | Global alignment harvests chance matches; use the longest CONTIGUOUS co-aligned block |
| H | Rustle's similarity grouping constitutes automatic family detection | precision vs independent Compara paralogy 4/12 = 33%; 154/195 universe genes unverifiable | Vindicates the advisor; the universe.tsv validation was itself built by sequence similarity |
| H | Clean the locus models S(v) (repeat-mask or strict-exonic) to remove bridge edges | repeat-mask −3 bridges/−12 true (net −9); strict-exonic −3/−16 (net −13) | Real LOC/ZNF paralogs share the same features as the bridges |
| H | Add a coverage condition ~C (fraction covered by cross-mapping reads) | pruned 5 of 196, only 1 a genuine bridge; 2 REAL paralogs pruned (GBA1/GBAP1 0.95) | A diverged paralog cross-maps only over its conserved region (cov_frac 0.10) |
| H | Combining orthogonal levers in a higher-dimensional classifier beats any single feature | in-sample 0.915 but grouped leave-one-family-out CV 0.80±0.05 vs best single feature 0.840 | Clear overfit: it learns family-specific mixtures that do not transfer |
| H | Use 2-edge-connectivity (bridge cutting) to refine the over-merged DNA blob | 1,547 genes / 24,431 edges / density 0.020 with only **0.3% bridges** | The over-merge is a dense transitive blob, not a bridge chain |
| H | Graph-structural signatures (density, Fiedler, min-cut) will flag over-merged families | n_genes raw AUC 0.75 → residual 0.52; density 0.36 / fiedler 0.35 / mincut 0.35 | Size dominates, and connectivity INVERTS once size is controlled: over-merges are DENSER |
| H | Union-of-records coverage (RUSTLE_ER_SUM_COVERAGE) is the biggest win of the session | single-record: 319 within-family, **0** cross-family (P 1.0000); union: 423 within, **45 cross** (P 0.9038) | It buys recall by admitting exactly what the per-record rule exists to prevent |
| H | Rebuild each locus rep from the union of its group's exons (RUSTLE_LOCUS_EXON_UNION) | Soto recall 65.5% → **44.8%**, pure families 135 → 72 | Widened reps inflate their own coverage denominator and fall below the 0.50 floor |
| H | Strip passengers with BLAST-style RECIPROCAL coverage (aln ≥ 0.50·max) | shipped (min) 464 edges, NPIP 19/19, 1 component; reciprocal 138 edges, NPIP 5/19 over 9 components | Ceiling is min(len)/max(len): only 134/171 NPIP pairs could EVER reach 0.50 |
| H | Buy purity by tightening the E_r identity floor 0.70 → 0.93 | precision 0.796→0.932 but detection 296→282, copies 1710→1425, **F1 0.438→0.433** | Tightening removes whole LOCI through distinct_locus_reps and the ≥2-loci gate |
| H | Identity 0.93 + no γ partition beats the shipped pipeline on all four axes | offline predicted F1 0.553 / ARI 0.549; two pipeline runs gave **0.433 / 0.434** vs shipped 0.438 | Pure precision/recall repositioning at constant F1; the offline prediction did not transfer |
| H | Relax the E_r coverage floor from 0.50 to 0.30 to recover lost true pairs | cov≥0.30: 402 T / 132 F (P 0.753) vs shipped 340/62 (P 0.846) — 62 true bought for 70 false | And 127/180 of the target pairs have coverage EXACTLY 0.00 — no floor can admit a missing record |
| H | Any-locus homology bounding fixes locus sizes (median 0.55× → 1.03×) | in-band (0.5–2×): shipped 47% vs any-locus **37%**; >2× goes 5% → 37% | Distribution-shift artifact: 21 pp of "too small" became 32 pp of "too big" (T4) |
| H | Replace the E_r coverage denominator with the read-supported core (depth ≥10) | edge-level predicted +0.023 F1; pipeline gave **−0.072**. F1 0.704 → 0.632 (floor 0.80) → 0.601 (0.60) | core/span spans 0.10–0.68, so one floor is a different demand per locus. Span is worse but STABLE |
| H | Allowing NON-CANONICAL splicing recovers the junctions the locus model misses | 12,526 of missed junctions (95%) are CANONICAL GT-AG | They are dropped by CHAIN SELECTION; relaxing the motif rule is a 3–4% correction |
| H | HP_sharedexon beats HP_sumcov, so the shared-exon rule is the better edge rule | identity 0.70→0.98 halves contamination under BOTH rules (225→120 and 173→124) | CONFOUNDED — the arms differed on rule AND identity floor; identity carried essentially all of it |
| H | The shipped E_r genomic-span union at `refine` is a measured no-op (0 blocks changed) | **11/26 = 0.4231 partitions moved** relative to the CORE (+288 edges on 322, components 52→32) | Inert only on the O1 homology path, by a code selection effect; it is sole support for one O2 family |
| H | `--tied-seed` answers the primary-flag coin-toss objection (recovers NPIPB12) | 33,415 tied secondaries → **2 tied reps, both abundance 0.0000**, inside IL27 and PLA2G10CP | Pays only costs: CAFAM0 byte-identical, no ID_154 member covered by either arm |
| H | Sequence coverage recovers ~80% of famCN's benefit, so no WGS depth is needed | direct re-measurement on the corrected universe: **31%** | The ~80% was never re-derived after the audit |
| H | famCN adds nothing to the partition (ARI 0.440 without vs 0.437 with) | repaired (tiled) graph: **0.757 with famCN vs 0.596 without** (+0.161) vs permuted null 0.173±0.008 | Measured on the crippled whole-region graph missing a third of its edges |
| H | Our famCN separates fewer over-merges only because we used 5 WGS samples | full 271-sample panel: 2/6 solid either way; GWFAM207 7.2/6.9 → 6.7/6.2 | The 35 GB download did not move the answer at all |
| H | Attribute a discovered candidate copy to the family whose members are nearest | **23/37 candidates landed on the wrong family**; the NCF1 candidate was a 52 kb GTF2I copy | Proximity is not homology; use id ≥0.80 & cov ≥0.50 |
| H | Admitting pseudogenes and lncRNAs into the definition is a recall lever | pseudogenes/lncRNAs are already **54% of RESOLVED members**; pseudogenes resolve at 66% vs coding 64% | E_r is biotype-blind and already admits them |
| H | Splice-junction / intron architecture can DEFINE families with precision | strict concordance breaks the 389-member bridge 389→32 but leaves APOBEC3/RFPL with **0 members** | Paralog intron lengths and counts DRIFT after duplication; no operating point does both |
| H | Drop BOTH-intronless edges (3,267, ~63% of the graph) as retrocopy noise | both-intronless would destroy real single-exon families at ~10:1 false:true; valid one-side form flags 2/2 retros, 0/5 reals | The panel hid it because all 3 panel reals were multi-exon |
| H | Switch the community weight from co-threading frac_ref to cov_min (max comp 238→44) | orthogroup-recovery F1 frac_ref **0.70–0.72** vs cov_min **0.62–0.64** | cov_min's smaller max component was OVER-FRAGMENTATION. Judging by max-component size is the trap |
| H | Take the connected component of the pairwise homology graph as the family | the 238-member and 196-member blobs at density 0.08; only 81% of structural families are cliques | Components over-merge via transitive bridges; the pure-clique fix over-splits tubulins/LILRs |
| H | minimap2's `ts:A:+` tag identifies a read's transcript strand | DSFAM817 copyA introns GT-AG (+) vs copyB CT-AC (−) while ts:A:+ on both | It is an artifact of `-uf` forcing forward; the GENOME intron motif is authoritative |
| H | The raw cross-chrom conflict graph yields valid cross-chromosome families | cross-chrom median exon-sum purity **0.50** (15.5% pure) vs same-chrom 1.00 (51% pure) | Over-merges via Alu SINE bridging: 150–300 bp covering 4–7% of the transcript, high entropy |
| H | Lower minimap2 asm20's identity gate to detect more divergent families | id 0.70 = id 0.80 = **152 families**, no change whatsoever | The 19-mer minimizer SEED is the wall; you cannot accept an alignment never proposed |
| H | Swap the asm20 nucleotide detector for the protein (6-frame ORF + mmseqs) detector | swapping drops 152 → **132 families**; as a UNION tier it adds +13 (12 real / 1 false) | 6-frame ORF breaks frame on introns, losing HBA1/2, GBA1, SHLD2, SORL1 |
| H | Sweeping minimap2 parameters recovers the mis-chained NCF1-class copies | 13 settings (−G 200k→3k, −r, −z, −p/−N, −n/−m): reads-local 68→13 but **copies_resolved = 0 in all 13** | Tuning provably cannot recover it |
| H | Split mis-chaining reads into local + paralog fragments to salvage local evidence | net **+0, 6/6** | The split creates a phantom node that over-connects the conflict graph genome-wide |
| H | Mis-chaining is the cause of the missing family members | −G 50k gives 0/76 giant introns and 96% aligned, yet NCF1 yields **0 families** under both modes | Mis-chain is a SYMPTOM of K=0 copy collapse; near-identical copies collapse to <2 reps regardless |
| H | c = 0.50 sits inside a robustness plateau, so the threshold is not arbitrary | widest no-change interval containing it is **(0.4721, 0.5000]** — c=0.50 is the ENDPOINT | It sits ON the boundary; the completeness knee is only bracketed to (0.50, 1.00] |
| H | Dropping non-spliced (single-exon) candidate loci raises O1 precision | drop all single-exon: P 0.709→0.664, R 0.698→0.560, F1 0.704→**0.608** (617→166 copies) | Deleting them costs PRECISION too; the ≥2-loci gate dissolves the small pure families |
| H | Filtering stub EDGES rather than stub NODES avoids the ≥2-loci gate | edges 2388→763 (−68%), copies emitted 617→**273**, F1 0.684 | "MY REASONING WAS WRONG" — dissolving edges dissolves families exactly as dissolving nodes does |
| H | Identity thresholds separate NPIP subfamilies and reject NPIP false partners | NPIP-NPIP vs NPIP-other: identity 0.979 vs **0.973**, coverage-of-shorter 0.912 vs **0.979** | On identity NPIP is one clique; COVERAGE splits it (A↔A 0.46 / A↔B 0.12 / B↔B 0.06) |
| H | Minimizer-Jaccard is a good precision axis for calling paralog pairs | EEF1A1 retrocopies score jaccard ~0.07 but core_ident **0.92** | It penalises short-copy-vs-long-parent length mismatch; use POA core-block IDENTITY |
| H | The genome-wide families are built WITHOUT arbitrary thresholds | the headline 152–157 catalog is Louvain + density 0.30 + asm20 refine | True only of the 82-family raw conflict object, which is not the headline |

### Medium / low redo risk

| R | Claim | Killing number | Why |
|---|---|---|---|
| M | The declared identity floor is 0.60 but the effective floor is ~0.75–0.80 | 71 of 73 dispersed pairs produce **NO alignment record at all**; the 2 that do sit at 0.84/0.85 | There is no effective floor on that stratum — there is SILENCE, and both die on coverage |
| M | The JOINT arm gives better family rosters than the RNA arm | Fisher p = 0.3881 (45/50 vs 41/50) and p = 0.5238 (3/5 vs 1/5), intervals overlap | Not significant; κ fires discordant 4/6 but at n=6 the test cannot reach 0.05 |
| M | The constant reduction succeeded — τ was eliminated | **12 of 13 (k,w) settings change the partition** on byte-identical node FASTAs | τ was renamed to minimap2's `k`, and the effect is non-monotone |
| M | minimap2's secondary cutoffs `-p` and `-N` materially affect E_r | `-p` ∈ {0…1.0} and `-N` 1…1000 give **byte-identical** results | `-X` retains all chains. Only `-N 0` bites (11/12 panels) |
| M | The binary's "≥2 distinct loci" certificate proves a family spans two loci | **4 of 7 strict single-copy controls emit families** (5 families / 10 copies); MAGEA 13 records at 11 loci | UNSOUND — it prints the certificate for two reps at ONE locus; survival is an expression ratio |
| M | Fix the split-locus false positive with a SPAN CONTAINMENT test | fixed **1 of 4**; HMBS overhangs by 1 bp, TFRC by 2 bp, TBP by 5 kb | Wrong shape: a junction-less rep's ends are read-derived and routinely overhang |
| M | The shipped catalog computes one definition on both the DNA and RNA paths | 1 family / 28 copies → **2 families / 27 copies** | `refine_families_exon_sum` is RNA-only, runs asm20 after the skip, and clusters by components not γ |
| M | TBC1D3's two extra loci are the genes TBC1D3F and TBC1D3G | the real extras are the distal-17q pair LOC105371848/LOC105371853 (227 and 229 spliced reads) | F and G are inside the NAMED core 9 |
| M | Use a taxon prior (gorilla should have ~1 NPIP copy) to filter false copies | gorilla has ≈17 morpheus copies (human 15, chimp 25–30) | Premise false and it is an ORACLE; a count-keyed rule would delete true copies |
| M | NPIP RNA graph density is 0.547 / 0.698 | correct two-tier: DNA 351/351 d=1.000, RNA 338/351 d=0.963, **ZERO RNA-only** | A uniform 0.80 identity was applied to both tiers; with `--refine` alone it is asm20 0.80 / sensitive 0.60 |
| M | POA core completion (`--complete-core`) extends the definition's divergence REACH | +55 copies, **all 93–100% identity**; reach unchanged at ~81–87% | A completeness gain, not a reach extension |
| M | θ=0.9 (all-copies support) is the right POA column-support threshold | APOBEC3 scores **0.019 at θ=0.9** (domain-sharer level) vs 0.06 at majority θ | One divergent copy or indel breaks the all-copies run |
| M | Use `contiguous_core_coverage` as the copy-level backbone relation ~B | 21 clean 2-copy families at 99% id / 80% cov score core **0.06–0.09** | It is the EXON-level gate; at copy level it fragments real copies (~20% of families lost) |
| M | Spatial spread of cross-mapping reads is the residual precision lever | spatial-spread AUC **0.52**, contiguous-core 0.54, vs rl_frac 0.899 | The unbiased feature-discovery method rejected both guesses |
| M | Coarse 50 kb genomic bins are an acceptable vertex definition | a spurious **336-node mega-component** | FP-robustness is CONDITIONAL on tight per-locus vertices |
| M | Apply the Bailey segdup filter (≥90% id / ≥1 kb / non-TE) to the definition | only **26% of pairs pass**; transitive single-linkage max still 317 | Drops CEACAM5/6/7, KRAB-ZNF, PRSS, IFITM, ULBP and only partly fixes over-merge |
| M | Exon containment is a better family edge oracle than aln_frac | exon containment AUC **0.900 < aln_frac 0.915**; single-exon wall inert at 1 of 102,455 loci | Mechanism also wrong: the GSTM2 hub shares EXTENSIVE exons, so it survives |
| M | The library-free VG repeat catalog (minimizer multiplicity) is a general separator | AUC **0.686** vs aln_frac 0.85; cuts 54% of borderline-real edges | Defensible only as a targeted mult≥20 repeat-hub gate (0/375 TP cost) |
| M | The tier change moves the pooled-exon coverage bite by 4.66× | under M1 the multiplier is **3.49×**; the 61-node cross-check is 3.81× not 5.3× | The script reimplemented coverage inline as the condemned pre-M1 form on a `-X` PAF |
| M | The 08-07 sensitive-only default retired asm20 everywhere | log prints "asm20 run SKIPPED" then emits **two `-x asm20` PAFs** from refine | The guard lives only in `homology_edges_all_reps_pooled`; refine calls `primary_seed_args()` |
| M | The direction of the tier discrepancy is favourable under the shipped tier | HERC2.rna density **0.4727 → 0.3778** — it stops clearing γ=0.40 | Refuted in direction |
| M | Both halves of the pooled-exon + union fix already exist in-tree | `RUSTLE_ER_SUM_COVERAGE` exists; there is **no pooled-exon E_r substrate** in the binary | "pooled" reaches a different any-one-shared-exon rule that REPLACES coverage |
| M | The cross-chrom catalog's low impurity shows it is the more precise O1 mode | member recall **65.5% vs 81.8%**; pure families 135 vs 227 | Its purity is bought by DELETION — the Web filter drops any ≥10-locus community at density <0.30 |
| M | `--cross-chrom` is the flag that makes family detection cross-chromosomal | "chrom" appears 3× in the function (contig set, a comment, a print) | It always selected the CONFLICT-GRAPH catalog; homology mode was always cross-chromosomal |
| M | NPIP separates cleanly under a 12 kb absolute block floor | min true block **9,859 bp < max passenger 11,400 bp** | The separation window is empty; 12 kb buys precision 1.000 only at recall 0.895 |
| M | Class-(b) missing edges are caused by FRAGMENTATION (one copy is a fragment) | length-ratio Spearman ρ = **0.191**; identity rejects 0 of 732; 127/180 have coverage exactly 0.00 | The cause is spliced-vs-unspliced representative SUBSTRATE MISMATCH |
| M | Giving co-threading more defined pairs will make it gain separation | co-thread AUC **0.866 → 0.743** while coverage AUC 0.198 → 0.623 | The earlier dominance was co-threading compensating for the stub mismatch, not signal |
| M | The genomic-span substrate's gain comes from fixing the stub-vs-spliced mismatch | both-spliced 0.259→0.914; one single-exon 0.000→0.022; **both single-exon 0.000→0.000** | 2,165 of 2,386 true pairs involve a stub and are UNCHANGED |
| M | Bound a predicted locus by its SIBLINGS (other members of its own family) | median 0.91 but **fallback 311/575**; SRGAP2's 4 members sit in 3 predicted families | Structural feedback loop: wrong families bound loci by the wrong siblings and cement the error |
| M | Use the read-supported core as a BOUNDARY (trim each locus to its core extent) | in-band 25% → 24% → 22% → **20%** as depth ≥3/10/25; truncation rises 64%→75% | Low depth is not a marker of "not part of this locus" (3′ bias, low-expressed exons) |
| M | Extending loci to the TES fixes the truncation problem | **0 verdict changes**; the ±10 kb retest is near-symmetric, 73% 5′ vs 76% 3′ | Implicates SD flank homology, not TES truncation |
| M | `RUSTLE_SPLICED_REP` is net-marginal and worth shipping | chr7 F1 **0.570→0.411**; chr16 **0.910→0.761**; loses NPIPB12 | Two real e2e runs regress; the earlier verdict was measured on a different mode |
| M | γ can break the GWFAM0 over-merge — it is a topological failure | GWFAM0 (37 copies / 7 truth families) **survives γ=0.40** and dies only at identity 0.98 | The blob's edges were genuinely dense, just diverged |
| M | Block sets differing between the JOINT and DNA objects: 0/5 and 0/7 | applying the SHIPPED gated-union rule gives **2/5 and 2/7** (APOBEC3, GSTM) at rep-transcript | The 0/5 figure was ENTAILED by construction, never measured |
| M | E_r edge counts are a property of the sequences — 351/260/91 is exact | renaming nodes alone moves 260 → 258 edges (symdiff 8 = **3.1%** on identical sequence) | `-X` implies `--dual=no`, and WHICH orientation survives depends on the sequence NAMES |
| M | Seed starved copies from AS-tied secondary reads (extend the skeleton set) | chr1 amylase 21 → **6 copies**, assignments 94k → 27k; 9-family sweep net **−2** | Tied reps poison conflict → refine → assignment with spurious edges |
| M | Tied-secondary seeding via shared-intron-chain agreement recovers Soto members | spliced-only = **0 Soto recovery** | Genuine-K=0 targets are DISPERSED; co-located misses fail the 0.98 AS-tie gate |
| M | Switching the tied-seed gate from AS-ratio to `de` fixes the phantom copies | skeletons 47 → 13 but **FINAL OUTPUT IDENTICAL** (same 17 copies, same 4,756 assignments) | `de` is more principled but is not the phantom fix — os1's ties legitimately pass |
| M | The locus-boundary correction from escaping junctions improves delineation | DNA in-band 23/44 → 24/44, McNemar exact **p = 1.000**; E_r precision 0.6775 → **0.4391** | Non-significant, concentrated in one family, and 8 new edges build a 48-node hairball |
| M | The E_r family relation is the same graph as the read-competition relation | **53.8% (56/104)** of absorbed families sent orphans to a paralog the family does not contain | E_r admits id ≥0.60 / cov ≥0.50, far more permissive than what minimap2 needs to place a read |
| M | SD regions are mosaics, so SD-level family calling is bounded at precision ≤0.456 | mean 3.43 → **1.33 families/unit**; ceiling loosens 0.455 → 0.985 on tiles | The mosaicism was an artifact of OUR merge of ~11k Soto units into 817 blocks |
| M | Soto's family splitter is parCN (QuicK-mer2 / 1KGP population copy number) | parCN separates **0/8** real over-merges (77% of families in [1.5,2.5]); famCN separates 6/8 | parCN measures population CN constraint and has no dynamic range |
| M | No sequence tuning substitutes for copy-number depth | coverage+famCN ARI **0.764 vs famCN alone 0.720** | False on three counts; the splitter is famCN, and coverage recovers part of it with no WGS |
| M | famCN is the fix for RNA-mode over-merges | famCN separates **2/6 = 33% in RNA mode vs 11/17 = 65% in DNA mode** | It must be measured over the DUPLICATION UNIT; RNA copies are transcript spans |
| M | The isoform guard's harshness shows NOTCH2NL/SRGAP2/RGPD are genuinely unspliced | NOTCH2's 490 reads give **240 distinct intron chains**, the biggest holding 5% | Rep-SELECTION artifact: spliced reads shatter in SDs so max(n_reads) picks the unspliced pool |
| M | Use the strand of an unspliced cluster to block merging with a neighbouring gene | **all 740 single-exon copies are '+'** (spliced copies are 149+ / 218−) | '+' is a PLACEHOLDER; an unspliced cluster can then never merge with a minus-strand gene |
| M | The χ(H) guard protects the `--cross-chrom` catalog from over-merging | `distinguishing_uniq` was never set, so `reads_distinguish(0,0)` was **always false** | Inert; hidden by a golden test that WROTE to the goldens it then compared against |
| M | RFPL is recovered with 4 copies in gorilla | CAFAM0 has 707 reads and overlaps **NO RFPL gene**; RFPL2 was missed | Two of four are 46–48 kb intergenic readthroughs in a gene desert where filter R4 cannot fire |
| M | An intron-count-concordance gate (count_conc ≥0.6) separates real families from bridges | held-out GGT1~GGTLC2: cDNA id 0.98 but count_conc **0.267** — cut by the gate | Panel-overfit |
| M | CAFAM0 is a 3-copy family with 213 assigned reads at 99.1% silver agreement | after the (chrom, STRAND) fix: **3-copy family → 0 families** | Over-merged: a 12-exon gene, its own unspliced version, and an antisense gene |
| M | Use `--min-copies 3` to keep only convincing multi-copy families | control: 2 same-strand disjoint copies → silver 100/100, dropped at min-3 | Real paralog PAIRS are common; use `--min-copies 2` |
| M | 82 > 58 confirms 58 was a lower bound, and 82/82 single-chrom is a finding | `colocated_families` is same-chrom **by design** | Softened to "consistent with"; 82/82 is a structural limit blind to DAZ/DAZL, RABL2A/B |
| M | Refined families are "pure by construction" | GWFAM18 (the GSTM array) carries **1 NBPF10 locus** | Construction certifies homology + distinct locus, not functional family membership |
| M | The 17 genomic-silent families are retrocopies | **0/17 intronless**; mean n_exon 9.26 vs confirmed 8.29 | They are partial/structurally-divergent high-identity paralogs failing on COVERAGE |
| M | The Rust `family_define` binary reproduces the frozen golden Python catalog | **574 families, md5 65ca0462** vs golden 573 / dca64cbd | The whole diff is one 213-node below-γ blob where Option B differs from seeded Louvain |
| M | No intron-size cap can separate NCF1's real introns from spurious paralog hops | NCF1 real introns are **all ≤3,450 bp** while the copies sit 368 kb and 1.54 Mb apart | A size gap DOES exist; the cap simply does not help (−G 3kb → 0 copies, coverage collapses 68→13) |
| M | Raising the coverage parameter to c = 0.51 fixes the observed false merges | one failing window has **1,204 bp of shared DNA at identity 1.000000, coverage 0.5571** | c=0.51 does not remove it; false merges stand at 2 of 150 windows |
| M | O1.13's blind/seed-free definition achieves 13.6% (purity 0.237) | both numbers **WITHDRAWN** — no artifact survives and the substrate was never named | Not merely unpinned |
| M | The catalog's emitted copies are spatially distinct loci | **30 of 194** two-copy families emit a copy overlapping its own sister | A node-definition problem that was never audited |
| M | Stub edges are predominantly false, so removing them is a targeted FP fix | TP 227→213 vs FP 93→85; precision moved **+0.006** | They are removed roughly proportionally from true and false pairs |
| M | Exon-length structure alone is discriminative enough to match intron chains | exon length alone produces **173k chance matches** | Only exon coupled with flanking-intron lengths (2-intron shingles) works |
| M | The intron-chain structural axis is a recall expander | at matched ≥4, **340/360** structural copies are also found by the sequence pipeline; 20 structure-only vs 7,964 sequence-only | Confirmation axis, not an expander; blind to retrocopies and single-exon copies |
| M | The DNA/protein two-tier definition dominates the RNA-level POA definition | DNA/protein misses DAZ1/DAZL (protein coverage **<50%**) which the RNA POA tier caught | Complementary both ways; both over-merge superfamilies |
| L | DAZ2's "31 vs 16 introns" shows it is a divergent second copy | DAZ2's rep covers **70.1%** of the annotated span; 5′ gap depth 0.17× with 1 intron-bearing read | It is 5′ TRUNCATION. DAZ2 is still a genuine second copy, for different reasons |
| L | GOLGA6L shows clean separation between members and passengers under a block floor | counting the two unnamed LOC accessions as passengers gives max passenger block **10,608 bp @99.8%** | The "clean separation" was a NAMING DECISION; the window closes |
| L | The genomic span is a stricter substrate than the exon-sum (introns diverge faster) | only **19 true pairs** are exon-sum-only; at id 0.90 span 503 T/46 F (P 0.916) vs exon-sum 327/33 | The span wins on both axes |
| L | Every badly-undersized locus reduces to the stub-representative problem | NBPF26: GFF span 114,464 bp, predicted 22,740 bp, ratio 0.20, **nexon 20** | A genuine exception — a spliced rep covering part of a very large gene |
| L | The genomic-span leg inside `refine` is a hidden FOURTH tier | same seed, same floors, same `cores: None` — only the SUBSTRATE differs | Historically correct (it WAS asm20 when introduced), currently obsolete |
| L | The original tied-seed design recovered TRIM64 and TRIM64B — 2 genuine recoveries | there was **no tied-seed skeleton at 89,893,563** | TRIM64's "recovery" was a rep-shift artifact; restated +1 genuine, −0 destroyed |
| L | DUX4L50's lack of exon recurrence is a DNA-block-only link our exon method misses | **99% softmask**; recurs only at 96–97% identity, below the 98% cutoff | It is a D4Z4 macrosatellite REPEAT, correctly excluded |

## 4.2 DEAD-END — never wrong, cannot be made to work

| R | Claim / proposal | Killing number | Why |
|---|---|---|---|
| H | Build a second, more sensitive nucleotide aligner to reach the dispersed families | 20 settings swept; the set with any gain AND recovery ≥ baseline is **EMPTY**. Best setting buys 14/73 at 1,142/12,508 cross-family pairs, recovery 0.65→0.25 | Every gain is a gap-spanning artifact producing a hairball; 57/73 pairs never get a record at any setting |
| H | Adopt protein alignment as O1's edge substrate to reach divergent families | 70/73 dispersed pairs are protein-homologous at E<1e-3 but median identity **0.2809**; gene lengths span 138× | It changes what O1 is: drops non-coding members and links anything sharing a domain |
| H | Make the joint DNA/RNA object DEFINITIONAL (union or intersection of E_dna and E_rna) | set fold: **0 blocks in O1, 0 of 7,691 reads in O2, 0 in O3**; matched DNA arm 351 edges vs joint 344 | Measured no-op everywhere and strictly dominated by a matched DNA-only arm. Property, not definition |
| H | A read-level filter (dropping a principled subset of reads) can cut two fused loci | read min-cut between fused loci: min 1, **median 12, max 273** | Max-flow bounds every read-level rule; the fusion is REDUNDANT, not chimeric |
| H | Tune γ, choose CC vs QC, or move a threshold to recover the dispersed families | **0 of 8** dispersed families carry a single within-family E_r edge | RULE-NEVER-FIRED: an edge-rule failure, not a grouping failure |
| H | Define families at the DNA level (blind self-alignment) to escape RNA's node problems | blind best NPIP component = 127 nodes = 19 NPIP + **121 other genes**, purity 0.136; 192-setting sweep found nothing | It MOVES the problem: RNA under-merges, blind DNA over-merges |
| H | Use the WIDTH of the multimapper (sec_frac) footprint to delineate copies | footprint CV 0.482 vs annotation 0.412; **A6≡A7, B6≡B7, B8≡B9 identical** | The track marks the duplicated REGION, not a copy; adjacent copies are indistinguishable |
| H | Recover NPIPA vs NPIPB subfamilies by thresholding identity or coverage | A/B merge at 0.9752 while within-B reaches 0.9782 — a **0.3%-wide window** | The cut would have to be placed inside it without knowing where it is; hierarchy finds it, a threshold cannot |
| H | Lower the POA/RNA similarity threshold to catch ancient diverged families | ancient recall TAS2R 2/18, DEFB 0, SIGLEC 0, CRYBG 0; over-merge gives 14 mega-components, largest 250 | All NUCLEOTIDE definitions top out at ~81–87% identity reach; lowering re-admits domain-sharers |
| H | Use graph topology (neighbourhood Jaccard: family = clique, bridge = star) | standalone AUC 0.826 but as a hard add it cuts **3 bridges at 9 paralogs** = base rate | A 2-copy family is an isolated edge with no common neighbours — and 57% of families are pairs |
| H | Start from DNA families (a DNA-first scaffold) to escape the RNA bridge problem | cDNA all-vs-all components: max component **389 members**; 19 families >20 members = 23% of 5,517 genes | The DNA scaffold is itself over-merged; it RELOCATES the bridge problem |
| H | Fix the O1⊥O2 leak by comparing INTRON CHAINS on the both-spliced side | the branch fires **0 times on every substrate on record** | Two overlapping same-strand reps with different chains are usually ISOFORMS; and it cannot be validated |
| H | Use an ABSOLUTE shared-block-length floor (e.g. ≥12 kb) as the E_r edge rule | lower bound 11,401 bp vs upper bound 726 bp = **EMPTY by 16×**; the invariant block spans 58.7× across families | No single number exists; macro-F1 0.952 @700bp → 0.300 @11,401bp |
| H | Blind-DNA NPIP purity 13.6% can be raised by a better E_r edge predicate | purity **ceiling 0.237** for ANY edge rule (shipped 0.156, oracle floor 0.164) | A NODE-construction problem: the best-overlapping nodes already swallow 61 non-NPIP genes |
| H | Pool EVERY isoform's exons into the shared-exon edge rule | rep-only baseline F1 **0.562**; best pooled variant 0.552, ≥2 exons 0.473, ≥8 exons 0.392 | Premise correct (12.8× more exons) but no pooling setting beats rep-only; both refinements fail too |
| H | Link every gene in the query REGION to the mapped exon's gene | literal region-level linking ARI **0.267** with a 1,102-gene component vs 0.514 for CIGAR projection | Merged SD98 regions average ~120 kb and hold many genes. NB: projection is OUR deviation, not their method |
| H | Measure per-copy famCN over the predicted copy's transcript SPAN | GWFAM42 within-family spread collapses **~60 → ~4.4** when measured over exons | Spans include introns, so read depth averages over intronic repeat content |
| H | Use `minigraph` to build the DNA variation graph of family copies | (structural) minigraph is SV-level and ignores SNPs | It would collapse near-identical copies into one path — exactly the copies we need separated |
| H | Use `--cross-chrom` as the mode for the known-gene-family benchmark | `--cross-chrom` **14%** vs `--homology-primary` **57%** family recovery | It builds families from read CONFLICT, which dispersed paralogs never generate — and it fails silently |
| H | Relax `GATE_MIN_READS` from 3 to 2 to recover more Soto members | ~**0** real members recovered, high mis-chain risk | |
| H | Define a family by its shared k-mer COMMON CORE (k=15), admitting by containment | core body is only **7–83% of gene length**, always UNDER (RFPL 24%, LILRA 7%) | The core boundary is set by sequence conservation, not gene structure; reproduces the user's earlier failure |
| H | Seed one copy, follow its multimapping reads to the other copies, each with its own core | RABL2A cores are **116–354% of gene length**; MAGEA9's core is **2244%** (the whole cluster) | Symmetric failure: it OVER-merges. Reciprocal whole-gene alignment is default because its boundary IS the gene |
| H | Use a variation graph to DEFINE families (family = one VG, copies = paths) | (structural) a VG presupposes its members | Circular. Keep the homology/quasi-clique graph for the definition; VG is representation |
| H | Re-align mis-chaining reads with a `-G 50k` intron cap to recover missing members | recovery **1/12** (and that copy is a 52%-softmasked single-exon artifact) vs regressions CDH12 4→2 and ID_481 6→4 | Net 0 genuine, −4 real copies: the cap truncates genes with real long introns (CDH12 510 kb) |
| H | Transitive closure / connected components over the pairwise homology graph | largest "families" **145 and 114 genes** = superfamilies; 428 size-2 families are clean | Chains subfamilies through domain hubs. It is OrthoFinder without MCL, and skipping MCL is why |
| H | Raise the POA contiguous-core coverage gate above 0.13 to remove transposon-sharers | IGFL3↔USP12 slips the gate at 16% core @98% id, but raising it kills RABL2A/B (**core_cov 0.17**) | Binding tension: the threshold that removes element-sharers removes genuine recent duplicates |
| M | Fixing the unbounded coverage denominator will change the DNA/RNA scale problem | DNA bite **6.79 → 6.79** (edge set byte-identical), RNA 34.78 → 35.05 | Mechanically a null: the bug can only INFLATE coverage, so it can only ADMIT edges. Ship it because it is wrong |
| M | Add the divergence floor as an O1 precision gate | E_p 0.89→0.98 but **R 0.877→0.51**, and it cuts 54–60% of the read-confusable O2 target | Opt-in only, not wired |
| M | Keep an identity floor ι in the backbone relation ~B | **0 of 393** τ-passing pairs fall below it (minimum 0.963) | Measured inert; ~B is reciprocal coverage ALONE |
| M | Avoid single-exon loci to remove false family members | **0** single-exon catalog members | No-op: the pipeline is already spliced-only and genuine single-exon families are dropped upstream |
| M | Derive the block-length floor blind (KDE antimode) so it is not annotation-fitted | blind floor 10,247 bp vs fitted 12,000; purity **0.162 vs 0.164 oracle vs 0.156 shipped** | Derivation was never the obstacle |
| M | Apply coverage-of-the-longer-locus in the SEEDED setting (filter seed hits) | **26 loci / 19-of-19 NPIP at every threshold ≤0.50**; at 0.90 it drops to 10 loci / 8-of-19 | Inert: contaminants cover the reference as fully as real members |
| M | Replace the shipped coverage rule with coverage-of-LONGER ≥0.50 as the FILTER | 129 T / 2 F, precision 0.985 but **recall 0.248, F1 0.396** vs shipped F1 0.833 | Usable only as an assert-or-abstain CERTIFICATE |
| M | Make homology-based locus bounding work by tuning the merge gap or best-partner | merge gap changes the result by **0.00** across 0→5000; best-partner in-band 36% still below shipped 47% | Intrinsic: read extent spills into duplicated NEIGHBOURS, and "is this duplicated?" cannot separate them |
| M | A single representative transcript (max-weight path) adequately represents a locus | junction recall of ≥3-read junctions **7.0% baseline / 11.9% co-threaded**; 98–99% of loci miss ≥1 | Wrong OBJECT: ~24 supported junctions per locus, ~2.8 carried by the representative |
| M | Drop the dominant anchor of a reciprocal tie (tie-orientation hub-suppression) | window-dependent (os1 flips to starved genome-wide); secondaries carry no SA tag ⟹ **68M-record index** | Adversarially marked DO-NOT-SHIP: it deletes the BETTER-supported copy |
| M | Blind DNA-level reconstruction of Soto's families from the assembly alone | **13/83** families mix genes | Irreducible: co-location of overlapping genes plus Soto's block-vs-gene curation |
| M | Locus-stitching (merging orphan de-novo fragments into neighbouring loci) | **ZERO** orphan de-novo fragments to merge | |
| M | Soto member recall is bounded by our algorithm, so better rules keep raising it | baseline 59.4%; 5 silent + ~75–90 genuine K=0 misses ⟹ RNA-side ceiling **~76%** | The ceiling is identifiability, not the algorithm |
| M | Stranded homology graphs (dropping antisense edges) will fix the domain-bridge blob | 482/15,479 = **3.1%** antisense edges dropped, but comp-238 has exactly **1** and is unchanged | Real precision axis that does not touch the resistant blob |
| M | Tune minimap2 (asm30, `-k`) to detect more diverged family edges | asm20 floor ~82–85% is a sharp cliff; only hand-tuned `-k9 -w3 -B2` reaches ~70% | There is no asm30; and the read-level de_max=.05 caps families at ~90% id anyway |
| M | Close the remaining `family_define` parity gap by matching Python exactly | the diff is confined to one 213-node blob | Closing it means re-porting networkx's seeded Louvain, undoing the RNG-free splitter the user chose |
| M | Fuzzy TSS/TES clustering is part of the thesis path | the thesis path uses only **±4 bp** internal splice-junction jitter | `tss_tts.rs` is marked DEAD |
| M | Feed both a default-G and a `-G 50k` BAM to the detector as parallel inputs | (design) double-counts reads and the mis-chains persist | Over-connects the conflict graph. The diff→per-read rescue variant caps at the +1 repeat artifact |
| M | Measure O1 precision using gene-group (symbol-pair) negatives | 133,677 aligned symbol pairs → 32 candidates, **0 of 32** certifiable non-homologous | The pair route is exhausted; resolution needed the UNIT changed from PAIR to WINDOW |
| L | `len_max` is the most bimodal edge feature, so it should separate real from bridge | bimodality coefficient 0.87 but **AUC 0.32** | Bimodal but non-orienting: it measures short-vs-long cDNA |

## 4.3 TRAP / TAUTOLOGICAL — the measurement was invalid

| R | Claim | Killing number | Why |
|---|---|---|---|
| H | Rule 3 (depth valley) closes 70.5% of the blind-vs-oracle gap, 178/251 | honest seed-free selectors (max-gap / Otsu) give 150/251 with a **SHATTERED** partition (38 components, largest 352 of 407) | MAY NOT BE QUOTED — c=0.20 is the ARGMAX against the truth denominator |
| H | The method proposes N new loci beyond the annotation | **60% of chr1 and 56% of chr16 assertions are SPLIT_RECORD**; identity median 1.0000; 78 of 102 unanchored components have min gap 0 | Off by ~2.5× until overlapping-interval nodes are merged: two nodes hold the SAME BLOCK SEQUENCE |
| H | γ-invariance of the shipped catalog is demonstrated (γ ∈ {0.20…0.75} identical) | **every emitted family except GSTM has density 1.000** | Vacuous, not passed — γ cannot bite on a complete graph (T11) |
| H | Numbers from the shipped `gw_family_catalog` describe the documented definition P1–P4 | the two differ on **five axes** (seed, γ 0.20 vs 0.40, τ 0.60 vs 0.80, identity metric, node set) plus a 6th | γ alone has four values in the tree (0.40, 0.20, 0.13, 0.30) |
| H | P4 shows the same 38 loci give an IDENTICAL partition on DNA and RNA | bite rate on one fixed 59-node set: DNA 4.34% / pooled exons 5.71% / **transcript 29.78%** — 6.9× from input length alone | The coverage clause is not scale-invariant; both legs ran where the floor barely bites |
| H | Compare spliced and genomic substrates with the single-record coverage rule | aggregating over records recovers **53 RNA edges and 0 DNA edges** | The single-record rule penalises CONCATENATED EXONS specifically |
| H | Protein E_p evidence discriminates genuine families from domain-share FPs | genuine ≡ `in_ep = 0` **by construction** | Circular. On the residual, nucleotide / E_p / TE / VG-topology all fail; only DNA copy number separates |
| H | γ=0.40 verdicts on real families are stable | HERC2 pooled exons 0.4727 → **0.3778** (loses γ); GOLGA genomic 0.3078 → **0.6106** (gains it) | TIER-DEPENDENT and flips in BOTH directions; P2's human row (d=1.000) is vacuous |
| H | O1 discovers 328 loci absent from the annotation | 370 → 328 in SD98 → 25 spliced → **3 with ≥10 reads** | DO NOT QUOTE 328 — 303 are single-exon, the class our own isoform guard treats as suspect |
| H | E_r coverage = aligned / min(qlen, tlen) is a sound membership criterion | a prediction at **10% of its true size passes at 1.00** | Structurally blind to truncation — which is why no coverage check ever caught the 0.55× problem |
| H | Raise the identity floor to 0.98 as the default — it halves contamination | costs **5.8 recall points** (81.8%→76.0%, 296→275 members) and 380 copies | Purity bought by EXCLUSION; report the other cell beside it |
| H | Make locus construction flag-free by solving max-weight bipartite assignment | the per-read argmax reproduces minimap2's primary site for **99.10%** of primary-anchored reads | The LP has only read-side constraints ⟹ separable per read ⟹ it IS the flag, solved exactly |
| H | Group SD98 regions into families by sequence identity | **999/999** region-region edges are ≥0.98 by construction; identity-only F1 0.089 < grouping by chromosome NAME 0.106 | A literal tautology on SD98 |
| H | Coverage splitting beats famCN on the repaired graph (ARI 0.861 vs 0.757) | 0.8611 placing **1,385/1,793** vs famCN 0.757 placing **1,637/1,793** | A different operating point, not an improvement — and this claim has flip-flopped FIVE times |
| H | Report a discovery's best homology hit ranked by IDENTITY | "100% identity to CNTNAP3" is a **63 bp fragment at coverage 0.008** | Below ~0.01 coverage the fragment bridges two cliques into a false merge; coverage is the guard |
| H | cDNA sequence identity separates real family members from homologous counters | panel AUC real-vs-homologous-counter = **0.029** | ANTI-discriminative: domain-sharers and name coincidences are MORE identical than diverged real paralogs |
| H | O1's 1.33% false-merge rate is the method's precision / error rate | 2/150 = 1.33% [0.37, 4.73] is a SPECIFICITY on single-locus windows with no positive stratum | No prevalence ⟹ no precision; and false OMISSION is separately unquantified (see §3) |
| M | Quote a family's div_strata median identity beside a chr1 result | kinesin genome-wide median 0.7418 vs chr1-restricted: **no alignment record at all** | div_strata identities are GENOME-WIDE; the label invites the wrong reading |
| M | The union substrate gives its own TBC1D3 subfamily structure | `TBC1D3.union.tsv` is **BYTE-IDENTICAL** to `TBC1D3.dna.tsv` | Every "union" subfamily number is a DNA number relabelled (quote only the byte-identity) |
| M | HERC2 bridges as one family | HERC2 DNA density **0.348 fails γ=0.40**, passes only at the code's 0.20 | Not a legal family under the DOCUMENTED γ; APOBEC3 union 0.400 sits exactly ON the boundary |
| M | The isoform guard (requiring spliced copies) is a valid family filter | PCDHB 15/16 loci passed **only because PCDHB2 picked up a 2nd exon** | It misjudges intronless families; they need the exon-overlap check |
| M | Exon-count concordance is a good real-vs-bridge discriminator (panel AUC 0.929) | corr **0.71 with gene size**; AUC 0.929 on panel → **0.755** genome-wide | A gene-SIZE proxy, not a homology signal |
| M | Edge betweenness is an independent validator of real family edges (AUC 0.683) | corr **−0.71** with log component size; AUC 0.683 → **0.531** after size-residualisation | Degree in disguise; and 0/12 Compara positives sit in components ≥50 |
| M | Family-definition precision is perfect (panel FP = 0, P = R = 1.00) | residual audit: ~10 low-density bridge suspects (**1.2%**) + 4 length-disparity among 826 coding families | That is the 17-candidate hand-checked panel, not a genome-wide guarantee |
| M | Coverage topup recovered NCF1, so topup is a working recovery mechanism | (no number) topup resampled to ideal depth | It brute-forces locally-chaining reads into existence; it corroborates the coverage diagnosis |
| M | The false-merge panel is an annotation-independent falsification of the definition | every eligibility filter is HGNC/RefSeq; **no coordinate or sequence test enters the funnel** | A genuine unannotated duplicate scores as a false positive — the exact thing O1 exists to find |
| L | TBC1D3 validates the absolute block-length floor as an independent criterion | block/min(span) = **1.000 exactly**; r(block, min len) = **+0.9999** | An affine relabelling of the already-shipped aln ≥ f·min(len) rule |

## 4.4 SUPERSEDED — right then, better number now

| R | Claim | Replacement | Why |
|---|---|---|---|
| H | There are zero RNA-only edges — E_rna ⊆ E_dna always | **33 of 481 = 6.9%** RNA-only on the 7-family/80-node panel | True only on the 27-node NPIP panel. Only the PARTITION claim survives: 1/33 joins two DNA components |
| H | E_rna ⊆ E_dna, 0/351 violations | 0/351 at 27 nodes, 9/316 at 61, **66/696 at 244 human reps** | Containment degrades monotonically with the node set; never quote it bare |
| H | Primary-only locus construction misses 66 truth members (66 of headroom) | true ceiling **27 members** (+7.5 pp), of which only 9 rest on a near-tied alignment | 61/66 have a gate-passing secondary but primary-only already builds a transcript at 35 — family-layer losses |
| M | O1 recovers 0 of 8 dispersed families — quote that as the failure rate | **1 of 24 = 4.2%** [0.0074, 0.2024] on genome-wide-dispersed chr1 families | The 8 were mis-stratified domain superfamilies under a chr1-scoped rule |
| M | Annotation accuracy of ±5 kb costs one family | recovery 22/22/22/21, but false-merge **0.0806 → 0.5000**, Fisher p = 2.7e-11 | Under-reports it 7×; the cliff is in false merges, not recovery |
| M | The effective E_r identity floor is 0.80 | **0.60** — `nucleotide_sensitive` is on by default and E_r is the UNION | asm20 contributed 0 sole edges |
| M | Raise the shared-exon FRACTION toward Soto's 0.99 to improve precision | id≥0.95 + 50%: **354 T / 17 F, P 0.954** beats id≥0.90 + 99%: 288/19, P 0.938 | Raising IDENTITY beats raising the fraction; at id ≥0.98 the fraction is irrelevant |
| M | Run the O1 catalog at γ=0.20 | γ=0.40 gives **identical recall 81.8%** with fused 55→47, trapped 249→225, pure 195→227 | γ=0.40 dominates. Never run 0.20 on that path |
| M | The shared-exon edge rule was tested and is not worth adopting | at matched identity 0.70: STRICT F1 **0.438 → 0.584** (+0.146); SOFT 0.579 → 0.661 | The original comparison was confounded on rule AND identity floor |
| M | Our own famCN separates 4 of 8 residual over-merges | only **2/8 SOLID** once the between-group gap is compared to within-group spread | GWFAM42's split is an artifact: within-group spread 47.6 vs between-group gap 38.9 |
| M | The exon-sum needs no TSS encoding (3/40 sibling pairs show a distinct TSS) | **TES 14/42 = 33%** distinct, median shift 4.9 kb, max 58.8 kb | True but 5′-only; polyadenylation choice is a different mechanism carrying ~5× more signal |
| M | The refined family count is 157 (54 cross-chrom) | 152/50 corrected (DNA-confirmed 90.1%) → 155/49 with the sensitive tier → 158 with `--protein-tail` | Three robustness fixes plus tier promotion |
| L | Confident-read (de ≤ 5e-4) extents put 58% of loci in the 0.5–2× band | **51%** over all 340 loci (the 58% excluded the 56 that fall back to the rep span) | |
| L | The 281-region core_recip ≥0.13 catalog / 58 families is the genome-wide set | 82 clean families / 207 copies from the de-tie conflict catalog | It under-counted (N-capped arrays) and contained a 728-copy over-merge; relabelled a LOWER BOUND |
| L | NPIP's 14 members resolve into FOUR predicted families | genome-wide COV_off at γ=0.40: the same 14 land in **13** predicted families | The dominant failure is FRAGMENTATION, not contamination |

---
# §5 O2 — COPY ASSIGNMENT

## 5.1 REFUTED / RETRACTED

| R | Claim | Killing number | Why |
|---|---|---|---|
| H | MAPQ = 0 is the multimapping proxy at duplicated loci | at NPIPB11 **234/271 primaries have MAPQ 60** while the window holds 25,287 SECONDARIES | minimap2 breaks ties confidently; the evidence lives in secondaries, which `-F 2308` removes |
| H | O2's real-data assigned rate is 0.3466, geometry fixed before scoring | on the external region BED: N 5,378 → 5,862, assigned **0.3308** [0.3188, 0.3429]; GSTM 0.8594 → 0.6939 | The denominator was INSIDE O1's control — swept regions equal the catalog's own extents ±1 kb |
| H | O2's commit rate demonstrates assignment beyond what the aligner provides | **1,806 of 1,821 MAPQ>0 commitments (99.18%)** name the primary-flag copy; novel = 30/5,378 = 0.56% | 98.4% a restatement of the primary flag; and 496 of 497 MAPQ-0 reads are one family's |
| H | O2's headline sim5x accuracy measures copy-assignment correctness | rewriting the BAM so every primary sits on the WRONG copy leaves planted-label accuracy unchanged while the headline goes **100.0% → 0.0%** | The number is a report on minimap2, not on O2 |
| H | O2 is max-weight FACILITY LOCATION with assign-or-abstain and no 1/k | `grep -rn facility --include=*.rs src/` returns **0 hits** | The shipped rule is a per-read independent argmax plus a Bonferroni gate. Different problems |
| H | Separation + linkage is sufficient for identifiability of all K copies | fails at **K ≥ 3** via cross-copy RECOMBINATION | Sufficient only at K=2 (phasing); multi-copy is strictly harder than pairwise |
| H | The per-read significance gate (Thm 4) never misassigns when δ ≥ 1 | **2,616 misassignment witnesses** when the completeness precondition is dropped | Soundness holds only if origin(r) ∈ C — i.e. it fails exactly in the reference-absent regime |
| H | Use intronic divergence (intron PSVs) to separate exonically identical copies | **99.6% of reads have an N in the CIGAR** | N means the intron was spliced OUT; intron PSVs need intron-RETAINED reads |
| H | Copy-specific splice sites resolve K=0 copies (18 junctions, ~3 of 11 families) | per-read **resolved_fraction = 0.00 for ALL pairs**, pair0 included | minimap2 snaps every junction to the nearest canonical GT-AG at EACH copy; reference-level only |
| H | Longer reads will escape the K=0 wall by spanning more of the transcript | FLNC reads are already the full-length mature transcript | The whole spliced transcript is identical between the copies |
| H | `tie_invariant` (≥3 unique-mapper/copy-specific-junction reads) certifies a copy is real | on the NPIP A/B panel **both intronic phantoms are tie_invariant = true** (2/2 and 4/4) | Invariance tests robustness to the arbitrary primary label, not biological reality |
| H | Wire the QV-weighted LLR engine in place of PSV vote counting | scoring change = **ZERO gain, identical in 16/16** configs | Under flat error the LLR margin is monotone-equivalent to votes. The lever is the GATE (min_psv 3 → n_decisive ≥1) |
| H | Winnowmap (or a minimap2 ∪ winnowmap merge) improves assignability | same 2,783 reads: mm2 mapped 2,783 vs winnowmap 2,782; **winnowmap-only recovery = 0** | mm2 `-N50 -p0.1` already saturates placement; the merge double-counts into false confidence |
| H | NCF1 is the flagship K=0 identifiability case (76 reads tie at MAPQ 0) | NCF1 has **22 uniquely-mapped MAPQ-60 reads**; 34/44 mis-chain over ~2 Mb with introns to 234 kb | A third bucket: paralog MIS-CHAINING, distinct from K=0 and from coverage |
| H | Abstention is the K=0 identifiability floor | CAFAM1 has **252 PSV columns** yet only 10% of reads span one; tied fraction swings 7%→90% | Two causes: K=0 (hard) and coverage-limited (under-sampling, improvable) |
| H | Exonic SNVs can rescue assignment among exonically-identical tandem copies | injecting 5 copy-specific exonic SNVs left **psv_cols = 1 in both arms** | Reads spread arbitrarily across identical positions, so per-position assembly mixes copies |
| H | A better aligner or more sensitive parameters will separate the TSPY copies | minimap2 default, minimap2 sensitive and winnowmap give **identical 40/240 unique-mappable, 1/6 invariant** | 4/6 TSPY copies are 100.000% identical over 2,782 bp; there is no signal to manufacture |
| H | The ~0.13% MAPQ-0 multimapper prior means re-alignment touches few reads | within-family candidate density **53–117%** | That prior is genome-wide; inside a multi-copy family the candidates are the bulk |
| H | Percent identity predicts identifiability — there is a ~90–95% "win band" | RABL2A/B at **99.9% identity give 0/200 tied**; DAZ1/DAZ3, less identical, give 41/220 = 19% | The axis is discriminative-minimizer density / max discriminative-free gap vs read length |
| H | Winnowmap will resolve the tied multimapping reads minimap2 cannot | minimap2 31 tied (14%), Winnowmap 32 tied (15%) — **resolved 0 of minimap2's 31** | The tie zone is byte-identical sequence, not a seeding artifact |
| H | Force a copy call on every read (naive-primary) rather than abstaining | naive-primary is **9/12 WRONG** on the information-theoretically unresolvable reads | Both graph and linear-AS honestly abstain. The lever is AS-tie ABSTENTION — a LINEAR technique |
| H | Use minimap2 as the PSV-discovery engine (270× faster) | **109 columns vs poasta 3,238**; all 8,461 assignments changed; GSTM 0 vs 745 | It clips divergent flanks so nearly all PSV columns vanish. It would move the O2 headline |
| H | NCF1 is the flagship K=0 example (second instance) | NCF1B 11/15 and NCF1 37/44 reads at MAPQ ≥20; 3 skeletons → 0 reps → **0 copies** | Coverage/assembly-limited; even MAPQ-60 reads are 100–150 kb chimeric chains |
| H | An iterative joint estimator (EM, longcallR-style) improves on the shipped gate | **0 disagreements among 3,081** co-committed reads; converted abstentions 0.2074 vs chance 0.2000 | Its only effect is turning abstentions into chance-level guesses; EM_CERT accuracy 0.0000 |
| H | O2's target population is reads with MAPQ 0 | MAPQ-0 fraction inside the 915 multi-copy loci is **0.0004** vs 0.0020 genome-wide | The stated premise describes 0.04% of the reads O2 is built for; near-ties are ~500× larger |
| H | O2 assigns multimapped reads to copies better than minimap2's primary flag | a secondary fits better by `de` in **463/23,675 = 1.96%**; net headroom ~0.1% of reads | Divergence and alignment score agree ~98% of the time. It measurably does not, on ~99.9% of reads |
| H | Re-scoping O2 from MAPQ-0 to alignment-score near-ties restores its novelty | tight ties (AS margin ≤1%) disagree 12.16% but are only **0.88% of primaries** ⟹ ~0.1% of reads | 6× enrichment as predicted, and still too rare to matter |
| H | O2's real-data accuracy is 99.28% | the headline reads **100% on a hand-made false family**, on a denominator conditioned on the assignment | O2 has NO accuracy number on real data |
| M | Theorem 2's identifiability condition is per-read C2 | C2 fails at **K=2** (phasing ambiguity) under exhaustive enumeration | Only STRONG SEPARATION survived (2M-instance attack) |
| M | A SUN makes a copy recombination-immune | machine-checked NOT-iff: S1 0/1.25M, S2 0/6.68M, S3 a K≥3 witness, S4 a no-SUN Strong-Sep K=4 instance | SUN is per-read sufficient but not cover-level immune; the K≥3 obstruction bites at the cover |
| M | O2's family definition (the E_c read-conflict graph) is adequate for assignment | planted copy B (MAPQ 60, no de-tie edge) is dropped ⟹ its **98 reads come back 100% TIED** | E_c drops uniquely-mappable members; `copy_assign` references homology-primary 0 times |
| M | The EM consistency theorem explains the observed 91–93% per-read accuracy | coverage 2→100: abundance-L1 0.10 → 0.0077 (13×) while certified accuracy stays **~100% throughout** | Per-read accuracy is governed by δ/e, not coverage; the coverage effect lives in abundance |
| M | Partition family members by (chrom, STRAND) so antisense transcripts do not merge | after the fix MAGEA goes **0 families → 1 family / 2 copies**, 895/931 reads assigned | An inverted duplicate is a real paralog pair; the antisense justification was a DEAD GUARD |
| M | A `tie_invariant = FALSE` copy is spurious and should be dropped | DAZ2 is tie=FALSE (1 unique mapper) yet junction-invariant, proven by an e2e relabel flip | The unique-mapper bound is conservative: FALSE means "not guaranteed by unique support alone" |
| M | Add reference-absent copies to copy assignment in a single stage | (mechanism at `copy_assign.rs:338`) it tightens Bonferroni α/(n−1) and adds a competitor | Single-stage addition REGRESSES O2 on both recall and silver; a two-stage freeze is required |
| M | The 70–100% MAPQ0 / 45–59% clipped rate at tandem loci shows linear alignment failing | primary-only: RCF_611 **15/218** MAPQ0, RCF_518 **3/53** | Secondary-alignment artifact — secondaries are MAPQ0 and clipped by definition (T14) |
| M | Tandem-repeat unit count ("structural PSV") moves the identifiability frontier | **10/10** K-frontier candidates refuted, zero frontier gain | Fatal dichotomy: where unit counts differ the unique regions already diverge 8–52%; where sequence collapses the counts are equal |
| M | Gate the VG correction leg to skip clean-linear and confident-PSV reads | sim accuracy **DROPPED 100% → 99.8%**; GSTM/MAGEA candidate counts unchanged, still 5.3×/51× slower | The gain comes precisely from confidently-PSV-assigned-but-WRONG reads |
| M | Lowering minimap2 `-p` surfaces a middle band of suppressed genuine copy reads | `-p` 0.8 → 0.6 → 0.3: **TIED stuck at 31, breaks 0** | `-p` governs VISIBILITY, not discrimination; a read's true copy is never a low-ratio secondary |
| M | Retuning minimizer k and w (e.g. k11 w1) will break the tie zone | DAZ tie zone **2,876–2,891 bp across k11w5 / k15w10 / k19w19 / k21w11 / k11w1** | It is byte-identical sequence; there is nothing to seed |
| M | DAZ's 2.9 kb identical core is below the HiFi read floor (intrinsic) | TIED reads median **2,410 bp**, reach 5′ end 3/41; DECISIVE reads median **4,224 bp**, 133/174 | The ties are 5′-TRUNCATED reads; 176 full-length reads resolve the copies. The lever is 5′-cap selection |
| M | Compute discriminative-minimizer density over the whole genomic locus | corrected exon-only: RABL2A/B 24%/20%; DAZ1→DAZ3 **14%**, DAZ3→DAZ1 **1%** | Minimizers on INTRONS are meaningless — spliced reads sample exons only |
| M | A T2T-grade reference improves paralog discriminability | the 2,875 bp identical core was measured **ON** the T2T assembly, 0-N gapless | T2T helps REPRESENTATION, not discrimination; a worse assembly's errors act as false PSVs |
| M | poasta produces the optimal PSV alignment | the banded DP finds strictly cheaper alignments under poasta's OWN cost model: GSTM **1181 < 1331**, PCDHB **3474 < 4152** | poasta garden-paths on repetitive transcript sequence, over-fragmenting gaps |
| M | The conflict-edge δ = 0.005 threshold is a meaningful tuned parameter | flipping to the IsoCon significance test at the assign gate's α is **BYTE-IDENTICAL** on all known families | δ=0.005 was doing nothing α does not |
| M | minimap2's MAPQ tells you whether a read belongs to the locus it was placed at | MAPQ **AUC 0.4944**; MAPQ=60 covers 96.07% of migrants vs 94.98% of natives | At chance. minimap2 is not merely wrong here, it is CONFIDENTLY wrong |
| M | Emit and count one copy-assignment result per BAM RECORD | contradicting molecules **4,616 → 0** of 7,121 after moving to one result per (molecule, family) | Feeding secondaries into the POOL is load-bearing; emitting per record was the bug |
| M | A de-novo two-pass (read-coherence + PSV) pipeline raises the copy-resolution ceiling | by K: K=0 → 0 PSV / unassignable; K=1 → 0.80; K≥2 → 1.00 — **the same identifiability boundary** | Its real win is a unified de-novo representation, not more resolution |
| L | χ(H) certifies the assignment it accompanies | on DAZ: n_copies = 2 with 2,213/2,353 reads assigned, yet famcn_readonly printed **1** | It was blind to `copy_junctions` — the COUNT used weaker evidence than the ASSIGNMENT |
| L | `align_traceback` can compute path_obs without handling strand | (fix 37c5af9) minus-strand copies — roughly **half of all families** — got garbage path_obs | Corrupted the EM for those families |
| L | At region r4, pure E_r loses a real family that E_c had found | the family that vanishes is CAFAM0: **816 reads, both copies collapsed, 90% tied** | It was never a family; the strand fix removed a readthrough/collapse artifact |
| L | The recovered copy catalog is invariant under primary/secondary relabeling | TSPY: **all copies 0 anchored, 30 flips, copies 5 → 2** under adversarial reversal | True for 6/7 named families; TSPY is the honest failure |
| L | Amylase AMY2A/AMY2B show ~21% tied reads, so the array is non-identifiable | true **0% tied (0/123)** after deduplicating alignments by reference position | The annotations overlap ~4 kb, so a naive per-window count double-counts one alignment |
| M | Use primary AND secondary alignments (not primary-only) to define the confident-read extent | primary+secondary **54% in-band vs primary-only 58%**; absolute `de` 58% vs per-locus relative trimming 46% | A good-`de` SECONDARY fits a paralog equally well, so it carries no paralog-discriminating information |
| L | The PSV-only vote engine is the most canonical copy-assignment engine | test L744: psv.n_decisive = 0 but **combined.n_decisive ≥ 1** | The COMBINED pipeline (PSV columns + copy-specific junction chain) strictly out-resolves it |

## 5.2 DEAD-END

| R | Claim / proposal | Killing number | Why |
|---|---|---|---|
| H | Run an EM (soft SDA relaxation) over the PSV-aware VG to improve assignment | the EM changes **ZERO evidenced decisions (0/3,081)** | A reframing, not a lever |
| H | Read-through into divergent 3′ flank rescues K=0 copies on real Iso-Seq | of 140 tied primary reads only **7 (5.0%)** are flank-bearing; 133 (95.0%) are exon-confined | Mechanism real in simulation (60/60), magnitude negligible: mature mRNA ends at the polyA site |
| H | Some RNA signal can resolve co-located copies with identical exon bodies | the two locus references are **byte-identical** (NM:i:0); junctions 0/489, ends 0/490, flanks 0/60 bp | The distinguishing-column set is information-theoretically EMPTY. Not a threshold limit |
| H | Calibrated per-base QV information meaningfully reduces misassignment | at realistic HiFi error: **~0** wrong assignments removed; only at 3× error does it halve them | The only thing votes structurally cannot do, which is why it keeps being proposed |
| H | Tuning τ (the LLR margin) is a meaningful operating-point lever | τ=2 and τ=6.9 give an **IDENTICAL** operating point (recall .964, agreement .9412) | The margin distribution is bimodal; the only material choice is argmax vs τ>0 |
| H | Make `--vg-realign-correct` the default (it improves ground-truth accuracy) | **4–50× slower** on real families (MAGEA 2 s → 101 s) for **+0.2%** accuracy | Inherent cost: within a family most reads ARE the MAPQ-0 candidates |
| H | Stop minimap2/winnowmap from producing multi-mapping ties | (structural) a genuine tie means the read is identical across copies over its span | Forcing a single primary only HIDES the tie behind an arbitrary flag |
| H | Formulate copy assignment as a global joint optimisation (facility location / LP) | the objective DECOMPOSES: L(j) = Σ_r log P(obs(r)\|copy_j(r)) has no read-read or copy-copy coupling | The per-read argmax IS the global optimum; empirically the EM changes 0/3,081 decisions |
| M | Formulate copy assignment as MIN-COVER = χ(H) and solve/approximate it | MIN-COVER is graph colouring — **inapproximable** | Keep only as the hardness boundary; MAX-ASSIGNMENT is submodular, (1−1/e) by greedy |
| M | Intron-retained reads will supply the intronic PSVs needed for K=0 pairs | intron retention is **~1–4%** of reads; pair3's introns are 0-divergent | Where it occurs the divergence sits in the largest introns, which no ~2 kb read retains |
| M | Coordinates/MAPQ can assign reads among 5 equally good near-identical copies | K=0: **100% MAPQ0, 0/5 identifiable**, PSV accuracy 0.20 (random) | Information-theoretically impossible; PSVs work only from K ≥ 2 = ⌈log₄ N⌉ |
| M | Use divergent intergenic FLANK sequence to break the TSPY tie | the whole ~8 kb tandem UNIT including intergenic sequence is **99.005% identical** | There is no flank lever |
| M | Replace the assignment gate with `o2_materialize`'s margin ≥ 2.0 | it drops real junctions (e.g. DAZ2) | Keep the threshold-free Poisson-binomial certificate; an arbitrary margin is what Canzar dislikes |
| M | Raising minimap2 `-N` improves copy resolution | n=24 array: copies/read 4 → 13 at `-N500` but **best-AS tie fan-out UNCHANGED**, breaks 0 ties | Visibility only. Useful just as tie-set COMPLETENESS when the tie set exceeds ~6 copies |
| M | Use depth-based copy-number estimation as a copy-count signal in the VG path | (qualitative) depth conflates copy number with expression level | Human-review flag only |
| M | The per-read robust-z of `de` against its own pile separates foreign from native reads | TPR **0.5066 → 0.2404** once loci where migrants exceed 50% are included (29/83 loci) | Structural: the robust centre then estimates the WRONG population |
| M | A direct pairwise DP can be made byte-identical to poasta for `discover_psvs` | per-read assignments differ **<1% (6,641 vs 6,707)** purely from co-optimal gap placement | Matching it means re-implementing poasta's exact tie-breaking — impractical and unnecessary |
| M | Use `depth_cn` as the copy-count leg for known gorilla families | depth_cn gives **DAZ 40.6 and GSTM 46** against true counts of 2 and 3; χ_H matched exactly (3,2,2,6,5,5) | depth_cn stays expression-inflated; χ_H is the trustworthy count leg |
| L | `RUSTLE_PSV_MINIMAP2=1` changes the assignment outcome | an **exact null** on this substrate — bit-identical output | |

## 5.3 TRAP / TAUTOLOGICAL

| R | Claim | Killing number | Why |
|---|---|---|---|
| H | `--phase` implements the N-copy generalisation of read-backed phasing (min-path-cover) | the implementation sets `hap = a.best_copy`: **2,298/2,298 rows match `.assignments.tsv`, 0 mismatches** | Any phasing accuracy scored against the assignment is trivially 100%; the real MEC module has 0 call sites |
| H | O2's unique-mapper agreement is 99.9% | `if assigned_j && region_mapq > 0 { fa.uniq += 1 }` | Denominator conditioned on the prediction (abstention cannot lower it), restricted to MAPQ>0 (empty in the target regime), truth = the pipeline's own copy list |
| H | Colour the RAW read-conflict graph H to get the copy number (Lemma 1) | colouring median **3× K**; a K=2 family needs **~58 colours**; only 23% of H are chordal | The raw graph counts SEQUENCING-ERROR edges. Lemma 1 holds only for the error-free graph |
| H | Build a consensus from the minority cluster's reads, then certify those same reads | (structural) the significance model assumes copy profiles are independent of the reads | Anti-circularity failure; needs hold-out or at minimum a discovery-coupled FLAG |
| H | Silver-standard agreement is 100%, so copy assignment is validated | the silver-checkable subset is **456 reads = 24.3%** of 1,874, and is uniquely-mapping by construction | Resolvable BY CONSTRUCTION. Do not headline the 100% |
| H | The method abstains on 44.9% of reads (842/1,874 certified TIED) | CAFAM1 alone is **44% of the sample at 90.1% tied**; the definitive broad run gives 24.8% | Read-weighted and skewed by one outlier family. DO NOT QUOTE 44.9% |
| H | MAGEA tied reads show RNA resolves near-identical copies at 64.9% | on the genuine K=0 subset (ref-vs-ref NM:i:0) the rate is **0/386 = 0.0%** | The "resolvable" subset is not K=0 at all — the reference loci actually differ |
| H | Recovery mechanisms (FAMILY_BOOST, topo_borrow, FAMILY_RESCUE) recover non-identifiable copies | DAZ at **99.97% identity**; RBMY held-out PSV concordance identifies **1/6** copies; TOPO_BORROW has 0 beneficiaries | The ceiling is DATA identifiability; all such recovery FABRICATES. Keep them OFF |
| H | The K=0 frontier result is a rate: 0 of 386 unresolvable cases resolved | the denominator is **1,692 under a looser rule and 430 under another** | Definition-sensitive, and K=0 is ENTAILED by byte-identity (3599/3599, 4030/4030), not measured |
| H | A column where reads differ from the reference is a copy-distinguishing PSV | Clair3-RNA found **5,677 PacBio / 2,549 ONT** A-to-I editing sites in HG004 alone | Editing mimics variants (A→G on +, T→C on −); must be flagged before counting as PSV |
| M | Any de-novo transcript spanning the family locus can be admitted as a copy | a **128 kb single-exon** de-novo transcript with 12 reads was admitted; RFPL assignment hangs >400 s | `prune_same_locus` clause (b) only drops a structureless span when CONTAINED |
| M | Per-record assignment (each BAM record scored independently) is an adequate unit | **73.8%** of multi-record molecules are confidently assigned to CONTRADICTORY copies | Defect O2.14 |
| M | The `psv_graph` "resolvable by K" count is how many copies a single read can be tagged to | {SUN} ⊊ {hap-vector-unique} | Vector-unique copies need multi-PSV CO-OBSERVATION and are exposed to the K≥3 recombination obstruction |

## 5.4 SUPERSEDED

| R | Claim | Replacement | Why |
|---|---|---|---|
| M | There is ~50% headroom over Eichler's AS ≥10 gate in discarded ambiguous reads | honest recoverable **28%** — linked-PSV recovers 410/1,487 | 333 discarded reads have a PSV-vs-AS conflict and 274 cover no PSV at all |
| M | 12.4% genuine K=0 from the PSV scan independently matches the 12% census exactly | exon-union re-alignment: **9.7% K=0**, 0 align_fail, 131/145 resolvable | The 12.4% was the contaminated genomic-span figure; the "exact match" was a coincidence |
| M | Exact-A* (poasta `AffineMinGapCost`) is a faster drop-in for the PSV aligner | 07-12 re-measurement: **byte-identical AND not faster** (PCDHB marginally slower) | The 06-27 "faster but different traceback" finding was itself overturned. Default OFF either way |
| M | minimap2-local PSV discovery is fine because CAFAM0 gives 4 clean PSVs vs poasta's 3,179 | CAFAM0 was an **OVER-MERGED FALSE family** (1-exon + 8,570 bp 12-exon + 3,221 bp 8-exon) | Not a speed/accuracy trade: poasta force-aligned length-mismatched copies into spurious columns |
| M | A direct banded pairwise DP is a free, output-identical replacement for poasta | the first version degraded PCDHB **assigned 6,707 → 5,816, ambiguous 72 → 442** | Stateless-traceback bug; salvaged only by a 3-matrix Gotoh with stateful traceback |
| L | DSFAM817 is the flagship copy-assignment case at 3/3 silver accuracy | **2/3 = 0.67** on re-check; cite the labelled-truth sim5x result instead | The silver standard is circular and thin |

---
# §6 O3 — REFERENCE-ABSENT / COLLAPSED COPIES

## 6.1 REFUTED / RETRACTED

| R | Claim | Killing number | Why |
|---|---|---|---|
| H | The missing-copy detector has 0 false positives among 23 single-copy controls | on real CHM13, 220 single-copy autosomal loci: it fires on **20.0% (44/220)** | The control is padded with hemizygous nodes that can never form a PSV column (7 have diff_sites exactly 0) |
| H | Route the unmapped reads — they carry the evidence for reference-absent copies | unmapped primaries = **0.12% (5,506)**, 92% fail to map even in DNA mode, 8% already IN the genome | A homology-sharing absent copy CANNOT be fully unmapped; it maps, mis-placed |
| H | The shipped two-cluster minority-divergence statistic detects missing copies | fires on **29.2% of intact FAMILY loci vs 18.4% of single-copy**; AUC vs NULL B only 0.619 | It is a PARALOGY detector, not an absence detector |
| H | A missing copy is detectable because the minority consensus fails to align back | positives **26.1% vs NULL B 40.0%/48.5%** | The nulls fire MORE; and a spliced consensus vs a genomic target breaks into per-exon records |
| H | A missing copy is detectable because the minority cluster forms no E_r edge to the host | positives **34.8% vs intact-arm control 55.6%** | The control fires more often — no discrimination in the correct direction |
| H | A compound structural/E_r-native rule will detect missing copies | positives **21.7% vs control 50.0%** (2.3× the positive rate) | The entire structural/E_r-native route is dead as stated |
| H | Escaping/private splice junctions are an identity-independent collapse signal | **58.6%** of matched INTACT-arm nodes carry ≥1 polymorphic junction site; 70–87% of intact single-copy loci already carry one | Premise true (59.0% private vs a 3.9% floor) but alternative splicing and a collapsed copy are indistinguishable |
| H | Adding a junction statistic to FARDIV (FARDIV OR jdiff) improves detection | FARDIV alone **32.43%** TPR at 3.23% FPR; **FARDIV OR jdiff 22.90%** — LOWER | The 2D grid gain (+2 to +4 arms) was killed by a 200-fold permutation |
| H | Gate 5 (reject at remap_identity ≥0.98) tests distinguishability from the reference | gate 5's reconstruction scores **0.9961** to the host while the TRUE copy0-vs-copy1 identity is **0.9754** | Its QUERY is the HOST-derived reconstruction ⟹ ~18% of real divergence; near-unpassable by construction |
| H | O3's shipped detector (FARDIV AND FARCLIP) will fire on real data | **0/915** in both arms; median FARCLIP **0.0006** against a gate of 0.05 (~80× below) | Whole-genome alignment gives a bad read a better home in 3.0 Gb; the mini-reference MANUFACTURED the signal |
| H | The missing-copy signature is CLIPPING (soft-clipped reads at the surviving locus) | whole-genome excision: **34.53% (61,512/178,145) of orphans become UNALIGNED**; FARCLIP median 0.0006 | The signature is UNMAPPED reads. No published collapse detector uses clipping |
| H | When a copy is missing the orphan read LEAVES the locus | FixItFelix measures collapsed regions at mean depth **51.5, ~1.5× HIGHER** | Wrong for COLLAPSED loci: the orphan is absorbed by its best-matching paralog wherever that is — possibly the collapsed locus itself. Excision: 64.2% ABSORBED vs 33.3% ORPHANED |
| H | GWFAM244:2, GWFAM382:3, GWFAM488:0 are missing-copy candidates | matched arm **6/915 → 1/915** after mis-chain removal, and that survivor was already a known paralog | All IMMUNOGLOBULIN: 63/69 reads share an identical 37 bp non-templated V-J junction; V mutated 1.77–3.27× more than C |
| H | Detector sensitivity is set by how abundant the hidden copy is | missed positives carry MORE migrant reads: m_true **0.4552 (missed) vs 0.3763 (detected)** | Sensitivity is set by DIVERGENCE: TPR 0.4500 at ≥0.01 diverged vs 0.0588 below |
| H | Use SDA-style depth-based collapse detection | (retired `collapse_gate`, default-OFF) a fixed background cannot tell absent-copy from generic paralogy | A matched decoy can. Only SDA's PSV discovery + read-correlation clustering carry over |
| H | Yoo/Rhie 2025 discovered unknown gorilla copies from Iso-Seq (so it is established) | 770–1,482 novel gene models discovered by **RE-ASSEMBLY**; Iso-Seq only confirmatory | ZERO papers discover a copy missing from an assembly using transcriptome data |
| H | The raw hidden-copy flag identifies reference-absent gene copies | fires on **7.39% (828/11,206)** of single-copy protein-coding genes vs background 7.93% | A non-specific, heterozygosity-dominated SCREEN; "CONFIRMED-COPY" relabelled "second-haplotype candidate" |
| M | The minority-cluster rebuild recovers the deleted copy | intact **22/26** vs deleted **19/26** — the degraded arm does WORSE | Fails its own control; zero discriminative power |
| M | Discriminator D4: the rebuilt node joining the family is evidence of a missing copy | positives **84.2% vs intact-arm 77.8%** | The host's own consensus joins its own family trivially; both arms saturate |
| M | Five V4b families are too divergent to align (blank medians) | the shipped `sensitive` tier aligns all five: 6/3/3/10/2 pairs at 0.762/0.741/0.707/0.728/0.727 | An asm20 seeding-floor artifact, not biology |
| M | GWFAM382:2 is a missing-copy candidate (FARDIV 0.0726) | **815/815** IGH-constant and 748/748 IGH-V reads from 1 of 4 pooled SRA runs, P = 10^−60.2 | Run-exclusive ⟹ index hopping; and the substrate has no B-lineage cells (MS4A1 13 vs COL1A1 183,670) |
| M | The O3 candidates are heterozygous alleles of the host locus | candidates are **30–40×** gorilla nucleotide diversity (~0.0015–0.0020) | Killed by magnitude. Does NOT kill a heterozygous DUPLICATION — that needed the haplotype test |
| M | The O3 candidates are driven by repeats / transposable elements | softmask **0.000 / 0.011 / 0.166 / 0.287** vs random controls 0.189–0.494 | Backwards: they are LESS repetitive than ordinary loci, all below the 0.30 gate |
| M | DSFAM213 (×3) is a reference-absent copy — the strongest raw signal in the scan | still **59–68% genome identity** after per-isoform ava-cluster + POA retry | Raw flag strength is not copy evidence; the assembly + homology + protein gate is the filter |
| M | Confirming a gene-conversion event on RECURRENCE ALONE is sufficient | a planted recurrent 8 bp direct-repeat RT-switch is mis-called by recurrence alone | Microhomology hotspots make RT template-switch artifacts recur; the microhomology leg is required |
| M | `--collapse-enumerate` recovers the known category-B collapsed families | **0 of 9** pre-labelled category-B families fired | Conservative by design: no balanced witness ⟹ the gate abstains rather than guess (2 real families re-admitted at 100% precision) |
| M | Derive absent-copy status from the `collapsed_copies` / `rescued_copies` counts | GSTM: **collapsed_copies = 9 with n_copies = 3** ⟹ every copy wrongly suffixed _ABSENT | It is a diagnostic count that can exceed n_copies. Absent iff an assigned read is discovery_coupled |
| M | The published 54-call ASJ number is reproducible from the shipped tool | it comes only from UNWIRED modules (asj_genetic_core + asj_strand_bias + asj_verify with SOR) | The shipped `asj`/`asj_verify` binaries emit a different TSV without SOR |
| M | Against the same-species T2T gorilla there is still orphan signal to mine | Yoo 2025 on gorilla: **0.7%** additional reads, **0.7%** soft-clipped >200 bp (vs 2% on the old assembly) | Measured on our exact species and data type; the authors treat it as a mappability statistic |
| M | Closing the envelope oracle leak COST detection performance ("quote 37.1%, not 41.1%") | the fix is **+1.62 pp [+0.27, +3.20]** on the fixed denominator vs the published −3.87 pts | The apparent drop was an artifact of the moving (called-only) denominator — closing the leak IMPROVED it |
| M | The O3/gate-5 path can admit a real multi-exon gorilla copy as configured | **zero** real GGO copies admitted; the demo is synthetic | The asm20 preset at gate 5 structurally cannot chain real multi-exon candidates |
| L | GWFAM227 failed on read support (MIN_TARGET_READS / read_count = 3) | its own artifact records the decline reason as **"≥98% remap identity"** | DIVERGENCE, not read support; read_count=3 is a cluster at the HOST locus (coincidence) |
| L | The corrected screen would have re-targeted GWFAM227 to copy 1 | copies 1 and 2 score **0.9629 bit-identically**; deleting copy 1 produced **ZERO** candidates | The screen cannot discriminate them; the re-target is strictly worse |
| L | GWFAM382 produced 5 admissions (a recovery) | PSV identity to the target **0.53–0.65** against a 0.90 threshold; the paired control admitted 3 anyway | VOID — none sits at the deleted copy's locus |
| L | The V4b shortlist was selected on GENOMIC identity while gate 5 uses SPLICED | `hunt_v4b.py` reads the exon-concatenated `copies.fa`: median seqlen/span **0.1292**, 0/1415 exceed 0.99 | It was always spliced; two agents reproduced the shortlist exactly, 24/24 |
| L | GWFAM34/GSTM's admitted candidate is a recovery of the deleted copy | candidate scores **0.7210** to the deleted copy while the untouched host scores **0.7265** | The admission is worse than doing nothing (0 hard reads, 0 anchored) |
| L | GSTM hits an architectural stop (it cannot reach the required cluster count) | it reaches **n_clusters = 3 at both floors** | A false prediction never executed; the clusters are within-locus PSV haplotypes, not a copy count |
| L | The L1 read-level identity lens is a CEILING on the spliced-consensus measurement | consensus scores **exceed** reads by +0.0005 to +0.0036 | Backwards — L1 is a FLOOR (magnitude tiny; verdicts hold) |
| L | Same-individual data does not exist for O3 | 4 SRA runs → **SAMN04003007**, the assembly's own BioSample; hom-alt 13.50× vs a het control at 1.096 | The matched fibroblast IsoSeq is the SAME CELL LINE as the assembly's DNA |
| L | The O3 candidates are explained by a heterozygous duplication on one haplotype | identity essentially unchanged on pat and mat; best gain any read achieves anywhere is **+0.011 to +0.014** | 668/691 = 96.7% of candidate reads sit <0.96 to the ENTIRE DIPLOID assembly (ordinary loci 0.9986) |
| L | GWFAM175:1/2 is a missing-copy candidate (FARDIV 0.067) | FARDIV **0.0667** vs the in-family pair distance **0.0629**; cleaning moved it 0.0001 | A known paralog explains it |
| L | Six collapse candidates were found by assembly-vs-assembly comparison | zero tolerance: **0 positives, 3 negatives**, sign test p = 0.25; at ±100 kb 6 up vs 5 down, p = 1.00 | Tolerance artifacts, some non-independent. ZERO loci called collapsed at 816 of 817 |

## 6.2 DEAD-END

| R | Claim / proposal | Killing number | Why |
|---|---|---|---|
| H | Divergent reference-absent copies show up in the UNMAPPED read pile | **5,519 unmapped reads = 0.13%**; 79% already present at 99.7% identity; 1 genuine hit | Dry on a complete T2T assembly — divergent copies FORCE-MAP onto a paralog. Real only against GRCh38 |
| H | A clipping-based transcript-alignment anomaly is a recognised way to detect an absent copy | 53 searches / ~30 primary papers: **not one** collapse detector uses clipping | The field's standards are S1 re-assemble, S2 depth+PSV, S3 unique peptides |
| M | The family definition can hold a copy with no genomic coordinates | (structural) V's elements are genomic intervals | The definition as written cannot hold a coordinate-free copy |
| M | Genomic reciprocal-overlap exon-class clustering renders copy-specific "arms" | GSTM `.exon.gfa`: end-to-end with 0 dangling but **3 PARALLEL CHAINS** | Only works for CO-LOCATED copies; dispersed paralogs have no genomic overlap — and dispersed is most families |
| M | The detector will transfer to genuine collapses, including near-identical ones | pairs below 0.98 fire **15/20**; pairs ≥0.98 fire **0/7**. FARDIV = 0.29 × (pair divergence) + 0.015 | Clearing the 2% gate requires ~2% divergence, but collapses happen BECAUSE copies are similar |
| M | S2 orphan splice endpoints discriminate a collapsed copy from an intact locus | null p95 reaches **1.0000** (human) / 0.8803 (gorilla) ⟹ headroom exactly 0.00 | The null already saturates, and it is largely a DEPTH statistic (ρ −0.37 to −0.56) |
| M | S3 the non-threadable read fraction discriminates a collapsed copy | on a synthetic exon-skip collapse the non-threadable fraction is **exactly 0** | Dead by proof: the leave-one-out reference admits any ≥3-read junction, so a collapsed copy validates its own junctions |
| M | Re-screening families with a corrected ranking statistic will produce O3 recoveries | moved **1 of 24** verdicts and produced **0** recoveries | The per-copy gate is honest but changes nothing; the lever is candidate CONSTRUCTION |
| M | Lowering `min_clusters` from 3 to 2 will recover reference-absent copies | **0 of 2** recovered; false-positive cost ZERO on 73 intact loci; 40 files byte-identical | The knob is neutral; the blocker is elsewhere (gate 5 / candidate construction) |
| M | Tune the collapse detector's operating point to improve unconditional sensitivity | unconditional sensitivity 0.153–0.190 (union 0.216); the product is **FLAT at 0.17–0.18** | Binding constraint is PROBE VISIBILITY (only 16–21% of masked copies are visible). The lever is a larger compartment |
| M | The intronless signature distinguishes an absent retrocopy from a second haplotype | only **2/58** candidates intronless, and both are single-exon genes | Uninformative by construction: an absent retrocopy's reads map to the spliced PARENT |
| M | The `--absent-copies` pipeline will admit real reference-absent copies on gorilla RNA | 1 region / 2,203 reads → **0 AC_ copies**; all 25 candidates → DnaNeeds (15× ≥98% remap identity) | Real headroom is ~0 on GGO with a complete T2T reference and the conservative gate |
| M | A transcriptome-first search of unmapped reads yields missing protein-coding COPIES | Fan 2015: 2,409 individuals, ~300 M uncharacterised reads → 2,550 transcripts, **>95% lncRNA** | The one large-scale attempt returns almost no protein-coding copies |
| M | Copy-vs-allele can be settled on RNA alone | the raw flag fires on **7.39% of single-copy genes** vs 7.93% background; verdict relabelled "second-haplotype candidate" | An information boundary: distinguishing an absent copy from a hyperdivergent allele needs DNA |
| L | DSFAM26 can be reconstructed by running `copy_assign --region` on its MHC window | it does not form via plain `copy_assign --region` | It is a `gw_family_catalog` family, not reconstructable from one window; GSTM was substituted |

## 6.3 TRAP / TAUTOLOGICAL

| R | Claim | Killing number | Why |
|---|---|---|---|
| H | Reference-free read clustering is BIT-IDENTICAL in the intact and degraded arms | (structurally invariant — no number exists) | True BY CONSTRUCTION: `-x ava-pb -X` never reads the genome and the excision edits only the ASSEMBLY. Do not quote sides.json |
| H | Structural statistics (junction diff, junction PSVs, orphan endpoints, threading) carry information | a pure read-depth proxy **outscores every one**: ndistinct AUC 0.585 vs jdiff 0.549, njpsv 0.548, orph_a 0.487 | Whatever they carry is smaller than "how many reads landed here" |
| H | FARDIV (1 − min of two cluster medians) is valid on whole-genome alignments | at GWFAM382:2 all 37 local reads sit at ~0.95, max 0.9641 — the reported **0.0726 WAS the mis-chained reads** | Undefined by construction when the read pile is UNIFORMLY divergent: there is no minority cluster |
| H | Use read-depth excess at the surviving locus as the missing-copy statistic | **16/104** destinations had ZERO baseline reads (ratio undefined); 12/104 exceed 5× on a tiny baseline | The 1.75× ghost needs a "before" a real case does not have; expression is not dosage. Quote the MEDIAN, never the mean |
| H | Use genome identity < 0.97 as the absence criterion in a genome-wide scan | below 0.80, chimeras and repeats pass: **42** such candidates had to be removed | Divergence needs a BAND (0.80–0.97), not a ceiling |
| H | The 4 promoted MHC copies are confirmed reference-absent paralog COPIES | all four cluster in **NC_073229.2:49–50 Mb, the MHC region** (divergence 3.9–20.4%) | MHC is hyperpolymorphic: "ref-absent copy" vs "hyperdivergent allele" is settled only by DNA |
| H | Our IG-locus candidates showing SHM-like divergence are a novel finding | Rodriguez 2021: V(D)J somatic deletions average **464 kb** (74.3–937.2); 29% of alleles differ by >0.05 | A named, published, guarded-against trap (IGLoo removes somatic-haplotype reads). Cite it, never claim it |
| H | O3 detection rates should be computed on units the method admitted | called-only **37.14%** arm vs FIXED **31.09%** [27.87, 34.51] (16.29% abstention); matched FPR 3.29% → 1.98% | 7th instance of denominator-conditioned-on-the-prediction (T1) |
| M | The 3 defensible unannotated human loci are novel gene families | the chr14 locus sits **5.4 kb from POTEG on the same strand**, and POTE is itself multi-copy | Most likely an unannotated paralog of a known family; needs ORF + cross-species protein evidence |
| M | Use DNA-catalog "absence" of a mosaic as negative evidence to veto a conversion | ref0 intervals cover only **2.88%** of the genome | The DNA leg must be POSITIVE-ONLY (Some(true)/None) or it downgrades real conversions |
| M | O3's measured TPR estimates real-substrate sensitivity | positives scored on a family-only mini-reference; hard-null FPR **0.18% vs matched 1.98%** | An UPPER BOUND — an orphaned read has nowhere else to go. Real sensitivity is UNMEASURED |

## 6.4 SUPERSEDED

| R | Claim | Replacement | Why |
|---|---|---|---|
| H | O3's detection rate is 31.09% (arm TPR on the fixed denominator) | **HOST-ONLY ARM YIELD = 189/743 = 25.44%** — quote 25.4%, not 31.1% | 42/231 fired arms (18.2%) fire only on a non-host node, which the same file prints as an FP stratum |
| M | O3's V4b removal-recovery shows the detector fails to recover deleted copies | superseded twice: substrate scarcity → **gate 5 is structurally near-unpassable** | Only 2 of 19 computable families fall in the recoverable [0.96, 0.979] window; GWFAM227's 0.9962 was an ablation artifact |
| M | O3's zero recoveries are explained by substrate scarcity (bimodal gorilla duplicates) | gate 5's reconstruction reproduces only **~18%** of real divergence; 0 of 4 ablations validated the screen | Every "declined at ≥98% remap identity" is an artifact of the reconstruction, not a biological verdict |

---
# §7 BENCHMARKING & EVALUATION

## 7.1 TRAP / TAUTOLOGICAL — the measurement was invalid

| R | Claim | Killing number | Why |
|---|---|---|---|
| H | The depth-valley rule clears its null at p = 0.0000, so it demonstrates delineation | a zero-information RANDCUT with matched node count, bp and piece count scores **152 / 100 / 146 of 251** vs Rule 3's 178 | The null was matched on EDGE COUNT only; the depth statistic is worth 26–78 edges (T6) |
| H | The blind panel result (251/251 for the seeded arm) validates the method | the panel truth **IS** the seeded pipeline's own output ⟹ an identity; FILL 245/245 with FILL≡SEEDED on 27/27 | Concordance, not measurement (T2) |
| H | Blind-mode panel survival (n/251) measures delineation quality | Pearson **r = 0.993 with raw edge count**; exploitable for 3.0× inflation by DESTROYING delineation | No precision term, and the denominator is conditioned on the arm's own match map |
| H | Blind mode's chr1 recovery of 0.4500 is comparable to seeded genome-wide numbers | the tier is not scope-invariant (`-k 11` exhausts its alphabet: 1.52 M minimizers, mid_occ 8456 at 488 Mbp); **18 of 40 families have exactly 2 members** | τ=0.50 is gameable at m=2; always quote all-singletons beside it |
| H | The unimproved annotation's own labels are the right baseline for the minimal-annotation gain | truth-free proximity baseline: +0.1000 (P 0.063) → **+0.0500 (P 0.266)** at 10 kb; m≥3 +0.2273 → +0.0909 | The input also supplies 2,056 gene INTERVALS, so a proximity rule is a strict subset of the method's input |
| H | τ = 0.50 recovery of 22/40 is a substantive result | **17 of the 22 = 77.3%** of the numerator is the m=2 stratum; τ=0.55 → 6 | The NUMBER is 77% pre-paid (the DELTA is not: −1 on m=2, +5 on m≥3) |
| H | Gorilla TBC1D3 has 9 members by symbol match | gorilla TBC1D3 is **0/9** by symbol | Symbol-prefix false positive: the three "TBC1D3*" names are TBC1D30/31/32 |
| H | The control panel is byte-identical for γ ∈ [0, 0.50], so γ is inert | γ=0.20 can only bind at **n ≥ 11**; the panel's largest component is **n = 10** | Structurally powerless: no data could make that panel reject anything (T11) |
| H | The shipped c = 0.50 passes a "must not regress" criterion | the panel's 78 E_r edges have **ZERO coverage in [0.50, 0.55)** (nearest 0.5984) | Self-referential — the baseline IS the c=0.50 output, so any shipped value passes |
| H | `known_family_showcase` validates O1's copy recall (published 0.786) | honest re-measurement: **0.524 (22/42 distinct loci)**; 5 of 8 rows read exactly 1.000 | Truth = the shipped catalog's own grouped members, so an over-merged member joins the denominator (T2) |
| H | The single-copy control panel establishes O1's precision (0/7 false families) | achieved N 7, exact one-sided **95% upper bound 34.8%**; not one control pair has cov ≥0.5 at id in (0.60, 0.999) | The identity floor was never exercised; the false-MERGE mode is untested by construction |
| H | Per-locus boundary accuracy is 0.796 in-band with ZERO truncations | LEAVE-ONE-OUT gives DNA **0.5227 (44/80)** and RNA 0.3611 | Every truth gene was also a seed, so truncation was structurally impossible (T2) |
| H | O1 shows no over-merge anywhere on the closure panel | (no number — the mode was suppressed upstream by node pre-selection) | VOID. The paralog-pair panel is still the missing control |
| H | The closure panel's 80 nodes exercise the documented node parameters | at the documented β=5 kb, d=10 kb the panel is **69 nodes, not 80** | 13.8% of the nodes the bridge result rests on do not exist under the definition being defended |
| H | An identity-floor sweep on the NPIP panel justifies lowering the global floor | all 27 nodes are ONE family with **zero negatives**; genome-wide the same floor gives one **776-copy hairball** | Rising density is guaranteed (T11) |
| H | Count cross-species orthologs by BEST HIT | best-hit says **6/19** gorilla orthologs; counting by qualifying edge says **19/19** | Argmax trap; it had INFLATED 1.7× elsewhere (T19) |
| H | The Compara-labelled panel validates the discovered read-level feature | **AUC = 1.0** mirage; replaced by external DNA-paralogy labels giving 0.829 on 467 edges | The decoys are NESTED genes, so one record double-counts and every bridge reads aln_ratio 1.0 |
| H | sim5x shows assignment accuracy is FLAT from K=2 to K=8 | `build_sim5x.py` encodes the copy index base-4 ⟹ **every copy pair differs at exactly ONE column** | It describes the GENERATOR, not identifiability |
| H | The 7-family negative-control stratum establishes O1's specificity | **5 of 7 negatives contain REAL families** | Nothing about specificity is quotable from it |
| H | The committed 76.2% Soto member-detection baseline is O1's baseline | it is a `--cross-chrom` number, **~16 points** from the homology default (65.5% vs 81.8%) | Never re-established on the shipped recipe |
| H | Offline (non-pipeline) evaluation of a proposed O1 change is a valid test | it produced a **"40-family blob" that does not exist** in pipeline output (real worst fusion: 2) | Three confident wrong calls in one session; treat any proxy as a hypothesis generator (T8) |
| H | Score an edge rule by the purity of the best NPIP-containing component | a blind user must choose among **377–455 components** | Oracle-selected: unobtainable in the setting it claims to describe (T3) |
| H | LRRC37 demonstrates that junction co-threading works | after the strand fix **all 6 LRRC37 pairs co-thread ≥0.89** | A signal that says yes to everything cannot discriminate; a positives-only demo proves nothing |
| H | Edge-level / node-level pair scoring is valid evidence a boundary change will help | three predictions passed offline and failed e2e: 0.553→0.433; +0.023→−0.072; in-band 25%→51% → F1 0.704→0.401 | Anything that changes what a NODE IS cannot be judged on node-level metrics (T7) |
| H | Absolute size-accuracy figures ("44% correctly sized") are a property of our method | in-band 44% vs RefSeq but **37% vs CAT** (CAT has 8,312 genes on chr1+chr15 vs RefSeq's 5,511) | Reference-dependent; quote the DELTA between arms. Soto's gene set is CAT-derived |
| H | "Bases explained" measures how well a locus model fits its reads | baseline **84.8% vs co-threaded 68.3%** — the baseline scores higher because 76% of its reps are stubs | It rewards UNSPLICED models (T5); junction recall is the sound metric |
| H | Family-level P/R/F1 headlines can be quoted without stating the tie-break | **176 of 296** detected truth members are overlapped by more than one predicted copy; oracle tie-break gives F1 **0.580 vs 0.438** | Truth members are gene-sized SD intervals, our copies are transcript loci |
| H | The identity-floor change (0.70 → 0.93) improves the catalog | SOFT scorer **+0.050 F1**; STRICT scorer **−0.005** | Scorer-dependent: the verdict flips sign between two defensible rules |
| H | Family impurity measured by `union_ab_score.py` reflects real over-merging | **16 pairs** of Soto truth members from different families overlap ⟹ a fixed artifact floor of ~16 families in EVERY run | One correctly-placed copy is claimed by two families no matter what we do (sum@0.70: 54 credit-all vs 38 real) |
| H | E_r edges can be re-derived offline from a PAF to test a rule | the offline re-derivation gave **3 edges where the pipeline has 17** | "The fourth offline-proxy error of this line of work". Must use `-c -X --no-long-join` and the per-RECORD rule |
| H | The genomic-span union fires on 1 of 7 gorilla families = a 14% firing rate | **n = 7 is not a denominator** | The supportable claim is EXISTENTIAL (such a family exists) |
| H | Score tied-seed changes by Soto member recall | member recall alone is PARTITION-BLIND | Same trap as the seeding-relaxation over-merge; score RECOVERED and PHANTOM separately (T9) |
| H | Quote the NULL A / NULL B false-positive rates (0.27% per locus) as the detector's FPR | 780 loci scored twice, same code: whole-genome **1.01%** vs mini-reference **4.16%** | Substrate mismatch. Quote the substrate-MATCHED 3.29% (16/487) |
| H | Escaping junctions are canonical GT-AG, so they are real splice events | `ts:A` is '+' on **100%** of records (1,146,371/1,146,371); motif predicted by alignment strand at 99.73–99.99% | Both BAMs used `-ax splice:hq -uf`, which forces every N gap onto a canonical motif |
| H | Identity values from different minimap2 presets can be compared across families | on byte-identical sequence the preset alone moves identity by up to **0.17** (GWFAM227 0.9629 → 0.8766) | A cross-preset delta must never be quoted as a biological one (T13) |
| H | The pre-registered secondary rule gives a matched-vs-cross difference at p = 0.0312 | at the correct FAMILY unit **McNemar p = 0.125** (6 matched-only loci but only 4 families) | The locus unit is anticonservative. DO NOT QUOTE 0.0312 (T12) |
| H | Soto 2025 is a gorilla resource that can validate gorilla family calls | Soto is a HUMAN study (SD98 on CHM13, famCN from SGDP n=269), **zero gorilla-specific results** | Its exclusion rule drops any family also duplicated in great apes ⟹ orthogonal to gorilla BY CONSTRUCTION |
| H | Report the probe-calibration rate over an ORDERED denominator of candidate pairs (4,330) | corrected **2,165 UNORDERED** pairs — exactly 2× | minimap2 `-X` skips self and dual mappings, so a pair is reported at most once |
| H | Mask the copy with fewer reads in the excision positive control | (the pre-registered rule masks the copy with the LARGER start coordinate) | Masking by depth ties the choice to the outcome; the accepted rule is verified uncorrelated with depth/span/divergence |
| H | Masking a whole locus out of the reference is a positive control for `--absent-copies` | GWFAM47: reads DID collapse (1,002 @ copy2) but 1 transcript → 1 rep → **0 edges → 0 families** ⟹ discovery never runs | It tests the DIVERGENT route, which is unbuilt; only the COLLAPSED route is wired |
| H | Our RNA catalog achieves 81.8% recall of Soto family members | **46.1%** at ≥50% member coverage; family-level pure-and-whole **6/83 = 7.2%** | 81.8% is ANY-OVERLAP DETECTION, not family identification |
| H | famCN reproduction is validated by a median ratio of 0.997 | within 1 copy **55.0%**; rounded integer match **34.8%**; a CONSTANT predictor scores 1.0000 | Median-ratio is null-informative (T10) |
| H | Split our families by famCN and then score against Soto | Soto's own famCN passes the between-vs-within test 6/6 **because they grouped by famCN MAD<1** | famCN is CONSTITUTIVE of the ground truth; the 6/8 gain measures nothing |
| H | Our Soto family evaluation reports 100% precision | honest: sensitivity **76.2%** (276/362), precision **91.4%** (393/430) | It defined a "known copy" as one overlapping a Soto member and moved discoveries to a separate table |
| H | The Soto member recall number is 76.2% | 65.5% single-pass · 71.0% per-chrom · **76.2% = a union of 4 legs, 3 of them unrebuildable** | Three leg files are absent from the cache dir and the scorer prints them as MISSING without erroring |
| H | The gorilla known-gene-family benchmark is independent ground truth | truth is defined in HUMAN by gene symbol (CHM13 RefSeq) and mapped by orthology — 14 families / 94 loci | Human annotation projected, not gorilla truth |
| H | Evaluate an edge-ADDING change by whether member recall improves | PCDHB recall stayed **14/14** while pure families collapsed **12 → 9** | Recall is partition-blind (T9); fourth appearance of the over-merge-vs-reach tension |
| H | Simulate TSPY copies using each copy's own annotated exon boundaries | byte-identical copies get spliced lengths **1147/1140/1108** ⟹ the naive sim gave garbage n_decisive = 835 | Per-copy exons inject FAKE copy-specific junctions; apply ONE shared splice structure |
| H | The bipartite size-match shows DNA-side predicted sizes match truth (median 1.00) | median 1.00, IQR **1.00–1.00 by construction** | CIRCULAR — `--from-genome` seeds representative nodes FROM the Soto member windows. Never quote it |
| H | The aggregate status histogram is sufficient to validate a PSV-engine change | GSTM per-copy support shifted **1196/1255 → 242/2209** with copy calls unchanged | The histogram masked a large per-copy redistribution (T17) |
| H | 70.9% orthogonal confirmation is the catalog's precision | 100/141 mappable; the nucleotide re-check is **circular**, and the 41 unconfirmed have 0 genuine false merges | It is a FLOOR, not a precision |
| H | 89.2% DNA-confirmed (targeted genomic-span alignment) is orthogonal precision | 140/157 = 89.2%, robust 85.4–91.7% | Only PARTLY orthogonal — the genomic span contains the same exons. Report as a LOWER BOUND |
| H | An oracle-measured prize proves a downstream fix will realize that gain | oracle-first averted a **−114** catastrophe and corrected mechanism attribution 4× | An oracle prize is a CEILING; verify the killers exist in your own population (T18) |
| H | Gorilla ground truth can be transferred from the human annotation | **40.8% MAPQ 0** (113/277 human members); median target-span/query-length **8.01×** | A human resource supplies both the roster AND the denominator for a gorilla headline |
| M | The TBC1D3 subfamily agglomeration leaves the family graph unchanged | (no number — post-hoc agglomeration over frozen files cannot re-partition) | The check has NO POWER; it would pass for any method |
| M | The O1-vs-O2 roster agreement table shows the two objectives build the same object | bite is **1 of 25** regions; at defaults GSTM gives 4 catalog copies but `copy_assign` gives **0 families** on 6,031 reads | It mixed pre-M1 O1 output with post-M1 O2 |
| M | Quote 0/386 tied reads resolvable as the K=0 resolution rate | a single-panel decisive-test count on byte-identical references | Not a rate over any population |
| M | Use the 27-node control panel to test whether RNA edges are redundant with DNA | that panel is the **COMPLETE graph** ⟹ redundancy null power **ZERO** | (T11) |
| M | Measure the detector's false-positive rate against single-copy genes | (not applicable — a single-copy locus has no family, so the envelope is undefined) | The detector abstains on all of them; the rate is inapplicable, not low |
| M | On the 1,333-gene fixed set the MAD peak at 1.0 is real evidence of fidelity | on the full fixed set the curve is **flat 1.0–3.0 and peaks at 3.0** (ARI 0.657/0.750/0.722/0.712/0.687) | At best a coverage/quality trade-off |
| M | The gate-sensitivity sweep shows the gates are robust | recall band **81.2–82.3%** while copies move monotonically **376 → 513** | It measures RECALL robustness only; the gates are a copy-SEPARATION knob |
| M | Soto misses labelled "not-expressed" are genuinely silent copies | SPAG11A has **19 primary reads** yet is labelled not-expressed (mean depth 1.6× over 21 kb) | The label is mean-depth-based and conflates truly-silent with sparsely-expressed |
| M | The high-divergence per-copy read splits produced by the new DP default are validated | the real splits (242/2209) have **NO independent ground truth** | The ship rests only on provable optimality plus zero sim_eval regression |
| M | `validate_genomic_dna.py` can be pointed at the genome with a relative path | a relative genome path reports **0/N confirmed** rather than erroring | Silent failure |
| M | Any overlap > 0 is enough to call a bipartite match | **18/232** pairs had neither interval 50% covered; a 47 bp graze counted as a match AND a 17.8× over-merge | Fixed by requiring reciprocal containment |

## 7.2 REFUTED / RETRACTED

| R | Claim | Killing number | Why |
|---|---|---|---|
| H | 47.7% of reads carry ≥20 bp soft-clips, so flank sequence is abundant | same region `-F 2308`: **373 primaries, only 1.1%** have a ≥20 bp clip (~40× collapse) | The audit counted SECONDARIES, which are clipped by construction under `-Y` (T14) |
| H | 86% of true co-family pairs share no exon — a representative defect of our exon rule | whole-representative alignment gives the **same 14%**; the spanning-forest minimum is 10%, observed 18% | Do not compare a pairwise-completeness number to a connectivity requirement |
| H | `multifamily_p1p4.py`'s P1 "NO" verdicts show the definition is not seed-invariant | `symbol.startswith(family)` puts **TBC1D30/31/32** into "TBC1D3" | An annotation-string artifact |
| H | Offline regrouping of a fixed `copies.fa` is a LOWER bound on pipeline performance | offline predicted F1 **0.553**, pipeline gave **0.433**; recall INVERTED (0.302→0.399 predicted, 0.302→0.283 actual) | Exactly backwards — it is an UPPER bound: the floor changes which copies are emitted at all |
| H | `bench/soto/replay_er.py` gives trustworthy F1 predictions now the copies.fa flaw is fixed | replay predicted **0.745**, pipeline returned **0.595** (recall 0.686 vs 0.446) | Two independent defects: it uses connected components where the pipeline uses γ-quasi-cliques. SCREEN ONLY |
| H | Rank families for O3 recovery by the median over all copy PAIRS | GWFAM478's pair-median 0.9761 names an identity **NO COPY HAS** (0.9537/0.9978/0.9978/0.9542) | A median cannot see a maximum, while gate 5 is an ANY/MAX rule |
| H | We reproduce 53% of Soto's families EXACTLY, 220/415 | honest **104/491 = 21.2%**; gene-weighted 13.2%, pair-weighted 3.4% | The 53% restricted BOTH sides to their shared genes; 69/105 exact hits are trivial 2-gene pairs |
| H | The 376 copies flagged `off_benchmark_family` are false positives | of 125 off-benchmark families: **119 REAL non-Soto SD discoveries**, ~2 spurious, 4 ambiguous | Genuine FP is ~2 repeat-bridge families + ~22 giant-span readthroughs, not 376 |
| H | Score the ambiguous-read regime as reads with `mapq == 0` | collapsed-segdup reads carry **MAPQ 6–16** with secondaries | The correct split is "read multimaps (has secondary OR MAPQ0)" |
| H | Constrain found loci to be roughly the size of the expected ones | Soto truth median smallest/largest member **0.33** (73% below 0.5×) vs our predictions 0.49 | Truth families are MORE size-heterogeneous than our predictions; SRGAP2C 5.6 kb vs SRGAP2 208 kb is a TRUE pair |
| H | Padding loci ±10 kb shows homology continues past the 3′ boundary ⟹ TES truncation | **73% continue past 5′ vs 76% past 3′** — near-symmetric; 90th percentile runs to the 20,000 bp pad limit | Confounded: homology continuing past a boundary is EXPECTED inside an SD |
| H | Every empirical accuracy headline is validated (silver 99.9%, 70.9%, 89%, 75.1%) | all four are self-consistency: silver = minimap2's own placement, DNA-span ⊃ exons, 70.9% is an annotation FLOOR | The one orthogonal check (SEDEF/BISER) is the one not run |
| H | "Recovers truth families at AUC 0.982" / "MAPQ-0 is usually solved" / "Thm1–3 executable" | the AUC is circular against minimap2; the silver standard is **EMPTY in the hard regime**; the abstain certificate is Python-only | `--vg` emits ZERO copy-assignment by default (`psv_linkage.rs:292`) |
| M | Blind mode sits at chance (6 recovered vs a null of 6) | corrected locus-pair null median **10 [5,16]**, p(≥6) = 0.9644 ⟹ blind sits at the **3.6th percentile** | Observation in LOCUS-pair units, null in NODE-pair units, coincidentally both 6 (T12) |
| M | The blind-mode metric was pre-registered | no phase-1 scorer imports `blind_metric.py` and PREREG.md's mtime post-dates all blind output | Pre-registration applies to the chr1 arm ONLY |
| M | We use Soto's shared-exon criterion (`bedtools intersect -f 0.99`) | `RUSTLE_SHARED_EXON_MIN_BP` default is an **absolute 100 bp** | A 3 kb exon matching 300 bp of a shared domain passes |
| M | 12.4% genuine K=0 independently matches the 12% combinatorial census exactly | corrected to **9.7%** over a different family set | The 12.4% was the contaminated genomic-span figure; the exact agreement was a coincidence |
| M | os1 is a phantom tied-seed copy with 0 paralogs — tied-seed precision is 5/6 | os1↔os4 is a **99.3% identity inverted duplication** (NM 18/2,543, mapq60). Precision is **6/6** | A metric bug: the eval aligned to the ANNOTATION, but os4 is itself UNANNOTATED |
| M | Escaping junctions show a meaningful antisense excess | escapes agree with locus strand 37.3%/27.7% vs internal 97.65–99.38%, but **alignment strand predicts orientation at 99.73–99.99%** | Orientation is determined by `-uf` forcing, not biology |
| M | A `same_gene` landing means the locus is TRUNCATED (and there are zero new fusions) | **19,134/19,134** human NULL A extents are byte-identical to their own GFF gene feature | Truncation is impossible there by construction; every "truncation" is a landing in a different overlapping gene |
| M | After excising a copy, ask whether orphans went to the CATALOG'S sister copy (15.93%) | of 104 ABSORBED families only **48 (46.2%)** sent orphans to the catalog's own sister; **56 (53.8%)** to a different paralog | Wrong question — the receiving locus must be found EMPIRICALLY |
| M | Soto's S1A caption says ">1 exon", so use >1 exon as the SD98 gene rule | ">1 exon" gives **74.9%** recall and drops PPIAL4A-F; ">=1 exon" with `-f 1` gives **100.0%** | The METHODS text and the data both say ≥1 |
| M | Our 9.6% SD98 excess is due to a length filter, satellite masking or a stricter identity cut | ≥1 kb 107.2 · ≥5 kb 106.5 · satBases=0 72.4 · fracMatch ≥0.99 73.4 · **v1.0 = 97.8 exactly** | Only the v1.0/v2.0 assembly switch reproduces the number |
| M | `recompute_perchrom.sh` reproduces the genome-wide within-chrom detection exactly | it scored **26.2%** vs the headline; chr1 45→180, chr15 18→132 when switched to `--cross-chrom` | The headline used global Louvain; the script used per-chrom connected components |
| M | Build the gorilla known-family truth by merging PAF hits within 5 kb | spans up to **234 kb** and an uninterpretable 14% | Only the rebuilt one-locus-per-human-member truth is valid |
| M | Match a predicted copy to its planted copy by a positional window | copies closer than **W = 3000** are mis-assigned | Nearest-start matching is required |
| M | Within-family smallest/largest size concordance is non-circular truncation evidence | within-family **0.49** vs per-locus **0.54** — a coincidence between two different quantities | Self-retracted. The PAIRED per-locus 0.54 and the 32% boundary-limited figure still stand |
| M | The TSS extension makes size agreement worse (43% → 35%) | paired: 6/23 change, 4 worse / 2 better, **sign test p = 0.69**, Wilcoxon p = 0.16 | The matcher selects the SAME 23 members, so the test must be paired. `RUSTLE_TSS_SNAP` is opt-in for absence of benefit, not harm |
| M | The Rustle/StringTie FP gap comes from bad junctions, exons, colors or graph edges | introns 99.3/98.0; bundles 3,351/3,430 exact; nodes **96.3% byte-identical**; only-RU junctions = 0 | The foundation is near-identical; the gap is OVER-ENUMERATION of real introns |
| L | `family_def_dna_pr.json`'s recall figures (0.66 crossmap, 0.81 tied) | direct recompute gives **0.51**, and 0.81 is not reproducible at all | An unrecoverable subset definition. DISCARDED for slides |
| L | At identity 0.98 only ONE impure family remains | **24** impure families under credit-all (8 after best-overlap) | Artifact of `head -22` truncating the scorer's top-3 listing |
| L | 654 genes are lost before famCN is ever consulted | correct: **334/1,679 (19.9%)** fail at the graph + **160 (9.5%)** isolated BY the famCN split | The famCN MAD split is not innocent |
| L | The recompute lost 6 Soto members relative to the previous run | after the scorer fix recall **76.2% → 80.4%** per-chrom / 84.5% genome-wide, **0 regressions** | The expression-collapse leg's col-7 projection loci were not expanded |
| L | The sim5x PSV ladder allele design `BASES[(c+j)%4]` gives each copy private bases | copies 0 and 4 **always collide** | Must use base-4 digits `BASES[(c//4**j)%4]` |
| L | A boundary snap can only shorten a locus | measured **12 lengthen, 15 shorten** | `snap_boundary`'s near branch returns `fallback.min/max(peak)` |
| L | Fixing terminal edge wiring will close the parity gap | net **−4** | Rustle's extra terminal edges are read-backed REAL endpoints |

## 7.3 DEAD-END

| R | Claim / proposal | Killing number | Why |
|---|---|---|---|
| H | Run Soto's map-back on our merged 817-block SD98 BED with their verbatim `-N50 -p0.5` | **106 of 817** regions emitted ONLY their self-alignment = 31.4 Mb = 32.1% of all SD98 sequence; 69 entire families structurally invisible | Merging ~11k units into 817 blocks makes each block's self-alignment the primary, and `-p 0.5` deletes the paralog chain |
| H | "Soto used DNA segdups in human, so find his families in gorilla" | the paper is "213 **human-specific** gene families"; gorilla appears twice, as comparator | Human-specific BY CONSTRUCTION — FOXO3 was removed for being "also duplicated in other great apes" |
| H | Test our family detection on well-known human gene families using the human IsoSeq BAMs | **0 reads** at GSTM, HBB, PCDHB, KRT, TUBB, S100A in both human BAMs | Both are restricted to Soto regions; the benchmark had to move to gorilla |
| H | BISER (via mamba) will produce the orthogonal gorilla segdup map | its precompiled align step **SEGFAULTS (-11)** under WSL2; the unrefined search output gave a bogus 1/157 | Alu-repeat noise linking gene loci to repeats rather than paralogs |
| M | The pseudogene stratum is a valid negative control for O1 | a BED holding PPIA + all 51 of its hits gave **23 blocks, 0 edges, 0 families** | Incapable of failing: retrocopies are dispersed and are not independently transcribed |
| M | Correlate Soto's per-family gorilla copy number against our RNA-derived copy number | Soto's per-family gorilla CN is on their internal share, **not public** | What was run instead is a gene-set concordance, which is not a CN correlation |
| M | Parallelise the Soto recompute per FAMILY instead of per chromosome | per-family **1,042 copies / 82.6% inflated** vs per-chrom 225 copies | Per-family parallelism isolates co-located families and over-detects |
| M | The gorilla annotation can adjudicate family membership on the families that matter | `GGO_genomic.gff` names **1 of 19** NPIP genes and **0** TBC1D3; 46.6% of its genes are LOC ids | It names almost nothing on exactly the multi-copy families the thesis is about |
| M | Phylogenetic gene-tree/species-tree reconciliation supplies ground truth for recent families | reconciliation is hypersensitive to gene-tree errors and fails on near-identical recent paralogs | The field switched to SD detection + read-depth CN for that regime (valid for ANCIENT families) |
| M | A retained-intron filter (parity lever 1) will remove false positives | the +25 ceiling is unrealizable; the geometry predicate would kill **133 real isoforms** | Oracle averted −114 |
| L | The published family recalls can be reproduced from their reference catalog | `bench/family_rna_refine.tsv` was regenerated one day AFTER the showcase, shifting every family id (70→69, 549→547) | That reference catalog is GONE |
| L | A2's LOC-anonymisation arm shows robustness to label removal | **A2 == A0 in 81/81** families across 10 cells; `visible_symbol` has ZERO consumers | A2 is A0 on a random subsample and the anonymisation is entirely inert |
| L | Use plain bedtools OVERLAP instead of `-f 1` full containment for SD98 genes | mere overlap **88.5%** precision vs **90.1%** for containment | Containment is the correct rule per their methods |
| L | Color + mm_negative segmentation will close the parity gap | oracle-first bounding gave **prize 0** | Aborted before implementation |
| L | Donor-snapping (parity lever 2A) will remove false positives | **prize 0** — the donors are real | |

## 7.4 SUPERSEDED

| R | Claim | Replacement | Why |
|---|---|---|---|
| H | The m≥3 result clears a size-distribution-matched null, so the effect is hard to reach | a position-respecting cutter raises the null **124×** (0.003 → 0.372 of 22) | Genuinely size-matched but POSITION-BLIND, and chr1 truth is 80% clustered. Still clears (0/2000), margin overstated |
| H | Our replication of Soto's family partition agrees at ARI 0.720 | **0.514** (our famCN) / 0.544 (Soto's) on Soto's fixed n=1,793 universe | 0.720 was computed over a shrinking gene set |
| M | Family-definition precision is 0.64 | shipped exon-sum refine: **0.94** (1,709 TP / 103 FP) | 0.64 is the RAW conflict graph PRE-refine |
| M | `aln_frac` achieves AUC 0.915 at separating true families from false merges | **0.830** held out on the hard TP-vs-genuine-FP task | The truthbar cases are easy; the genuine FPs are the residual |
| M | "Pure-and-whole 6 of 83 families" is the catalog's family-level score | **12/83, 9/83 or 5/83** depending on the reading — not 6/83 | Does not reproduce and is definition-sensitive; every quote must state the definition |
| M | The co-located resolution census is 258 assignment-relevant pairs (77.5% resolvable) | exhaustive: **236 pairs / 97 families / 42,313 reads**, 87% resolvable / 13% K=0 | The old census was SAMPLED and used a 3-way bucket |
| M | The O3 detector achieves TPR 41.13% per arm | whole-genome held-out **TPR 0.2703** (10/37) [0.1540, 0.4298] | Every positive was scored on a family-only mini-reference ⟹ an UPPER BOUND |
| M | The detector's published hard-null FPR is 0.18% [0.07, 0.46] (4/2210) | clean-read value **0.00% (0/2210)** [0.00, 0.17] | 100% mis-chain ARTIFACT (each FP carries an N gap 10–30× its own locus). Restate 0.18% as an artifact CEILING |
| M | Our Soto replication reproduces the gene set at precision 90.1% | **96.19% (1793/1864)** | An `ID_\d+$` filter silently dropped all 114 `Unassigned_*` singleton rows |
| M | Family-level Soto recovery is X of 83 families | denominator **76** — 7 Soto families have a single member | Quote FOUND+isoform 48/76 = 63% |
| M | DNA mode wins — COMPLETE 87% / FOUND 100% vs RNA's 78% | with the purity guard applied like-for-like: **COMPLETE 46% / FOUND 57%** (RNA 67%) | The DNA numbers were measured without the purity guard |
| M | Same-individual gorilla data does not exist, so every gorilla number is cross-individual | the blind-spot audit measures on a **matched fibroblast BAM** (n = 1,537,238 reads) | The live limitation is one individual / one tissue, not the absence of matched data |
| L | Our v2.0 SD98 recovers 100% of Soto's autosomal paralogs (1831/1831) | **99.89% (1831/1833)** — S1E has 1,833, two of them liftover failures | |
| L | The docstring figure for region-level linking: a 178-gene family at ARI 0.20 | not reproducible; the correct control is **ARI 0.267 with a 1,102-gene component** | |
| L | Copy-level precision against Soto is 91.4% | **88.2%** at member level with candidates counted as calls | |

---
# §8 INFRASTRUCTURE & TOOLING

## 8.1 REFUTED / RETRACTED

| R | Claim | Killing number | Why |
|---|---|---|---|
| H | The test suite protects O1's defining parameters | **714/714 passing** under min→max, summed coverage and re-enabled asm20; only 9 of 246 O1-path tests run an aligner | 3 of 4 defining parameters change silently; the coverage denominator is untestable (all tests use equal-length sequences) |
| H | The minimap2 `-p 0.5` DNA bug has an RNA analogue, so `-N 500 -p 0.01` recovers edges | adding them to either tier gives a **BYTE-IDENTICAL PAF**; dropping `-X` collapses asm20 **2777 → 1053** edges | `-X` sets `-P`, so minimap2 skips primary/secondary selection entirely. `-X` is PROTECTIVE |
| H | `cargo test --lib` is sufficient to catch regressions | two broken tests hid a real inert-guard bug and a fixture deleted 16 days earlier | Use `--all-targets` (see also the compile break below) |
| H | Exact-A* in poasta gives a 2.4–2.8× speedup for `discover_psvs` | A/B: **0.0 / 0.5 / 12.7 s**, PCDHB marginally SLOWER, families/quant/assignments diff = 0 | The agent's benchmark was spurious; poasta's min-gap heuristic does not prune paralog alignments |
| H | The asymptotic audit's prediction that per-read assignment is the bottleneck | per-read assignment **~1.1 s** and `build_read_placements` 0.0 s vs POA diagnostic **~458 s** + poasta 69.7 s | Direct measurement with stage timers showed the audit was wrong |
| H | Auto-build the splice `.mmi` index so `RUSTLE_PROJECT_MMI` works transparently | a splice `.mmi` of a 3 Gb genome is **13.6 GB**, and the auto-build leaked it into `/tmp` | Made strictly OPT-IN against a persistent index (measured win when opt-in: 2.1×) |
| H | `~/winloci_scratch` is a physical Linux disk with plenty of space | writing a 90 GB BAM drove **C: to 100% with 2.8 GB free** | It sits on the WSL2 ext4 VHDX whose backing file lives on Windows C: |
| H | Deleting large files inside WSL frees the space back on C: | (structural) a WSL2 VHDX does not auto-shrink on delete | Only `wsl --shutdown` + `Optimize-VHD -Mode Full` reclaims it; `fstrim` needs sparseVhd=true |
| H | The genome-wide `rustle --vg` run crashed the terminal because of a rustle bug | one run peaks at **4.5 GB RSS**; the launcher ran 3 concurrent ⟹ ~13.5 GB + page cache > the 19 GB ceiling | WSL2 OOM: truncated GTFs with NO panic line (abrupt SIGKILL), memory fully reclaimed after |
| M | The `homology_catalog_never_touches_the_conflict_graph` guard proves O1 does not use E_c | it passed while `locus_unique_mapper_counts → reads_distinguish` flowed | It banned FOUR HAND-PICKED STRINGS; the leak uses different names and lives in a helper |
| M | `read_junction_support`'s junction counts validate the read-through filter | DAZ readthrough: **56 distinct junctions on primaries vs 154 with secondaries** | It skipped only `is_supplementary` while `bam_reads` holds secondaries (~63% under `-N 50`) (T14) |
| M | The tied-seed A/B result is a measured null for the O1 family catalog | `tied_seed` is read only inside `detect_and_assign`, and `DenovoConfig::from_env()` never sets it | A PLUMBING GAP, not a null — the lever cannot currently affect O1 at all |
| M | Widening the remap scope from a window to the genome can only lower best identity | **14 copies move UP, 9 flat, 6 DOWN** (to −0.1431) | Selection maximises BLOCK LENGTH, not identity: 1222 bp @0.8198 beats 1209 bp @0.9603 |
| M | `hunt_v4b.py`'s `best_target_identity` implements its docstring and mirrors gate 5 | it disagreed with its own docstring in **10 of 22** families; the driver passes a ~150 kb LOCAL window as `--fasta` | Differs from gate 5 in query object, target type, self-exclusion and fasta. 8/24 families tie bit-exactly |
| M | The shipped mis-chain guard (`is_giant_intron_mischain`) protects per-read statistics | the 827,011 bp mis-chain at GWFAM244:2 is carried by **97 primary reads — 32× the <3-read gate** | Blind BY DESIGN: it fires only on unsupported junctions. A POPULAR mis-chain is evidence of SYSTEMATIC mis-chaining |
| M | mGorGor1 is Kamilah's genome | BioSample **SAMN04003007 = GgoY_KB3781, donor "JIM", MALE**; the full 67 Mb chrY carries SRY/DDX3Y/UTY | Kamilah is a different, FEMALE gorilla (gorGor3/4/5/6). Fixed in 5 memory files |
| M | A v1 interval "cleaner" clustering unclipped read spans tightens the O3 windows | it **GREW 870/915 loci**, 39.80 → 73.44 Mb, >50 kb count 213 → 613 | Not cleaning — it silently changed what a node IS. Retained as `intervals_v1_UNCLIPPED_WRONG.*` |
| M | Gate-5 remap with minimap2 asm20 is adequate for checking a candidate exists | P2 simulation **4/4 pass** only after switching to `-cx splice` | asm20 fails on real multi-kb introns, so genuine copies could never chain |
| M | `family_detect::candidate_pairs` canonical k-mer LSH is a live decision step | deleting it leaves families/assignments **BYTE-IDENTICAL** (loses one log line and `.fallback.tsv`) | Diagnostic-only unless `--complete-core` is opted in |
| M | `FamilyGraph` is live via `consensus.rs`, so wiring `to_gfa()` gives real VG output | the struct is built only at `consensus.rs:284/310` — **tests** | Consensus-correction is not called in any live binary flow; `to_gfa()` would be dead code |
| M | The editing filter is a 4.5 s hotspot in per-region `copy_assign` | real cost **~0.1 s**; the true hotspot is `build_family_profiles`/`discover_psvs` at 4.4 s | A TIMER DOUBLE-COUNT — `t_psv` was reused cumulatively |
| M | minimap2 subprocess spawns are a pipeline bottleneck worth eliminating | spawns are **~5–15 ms each** | The bottleneck is PSV alignment, not process launch |
| M | Running BISER/SEDEF for the orthogonal segdup validation is blocked on WSL2 | SEDEF builds and runs end-to-end (fracMatch 0.985 on synthetic); the real limit is **~17–40 h at -j2**, and WSL2 rebooted at 832/1813 | BISER's −11 is a runtime crash in its Codon binary. Move to a cluster; partial output is resumable |
| M | FxHash is 26% faster than SipHash on HERC2 (13.423 s vs 9.867 s) | these regions are BIMODAL under a fixed binary (~10.0–10.5 s or ~13.0–14.1 s); interleaved n=5 gave FxHash **2/5, 3/5, 2/5** | A sampling artifact; SipHash was actually faster on HERC2 (0.958) |
| M | Correcting locus boundaries upward is cheap enough to run genome-wide | E_r input **5.6 Mb → 22.3 Mb (4.0×)**; the sensitive tier went ~18 min → **>70 min for ONE alignment** | Truncation was making the pipeline cheap; and the sensitive tier is the necessary one, so `SENSITIVE_ONLY` cannot help |
| L | Cloning `mapq_reads` so it outlives a drop is harmless | peak RSS **doubled to 15.5 GB** and thrashed for an hour on a 94%-secondary BAM | Compute from the borrow BEFORE the drop |
| L | `bench/FALSE_NEGATIVES.md` is a citation that never existed | committed **4586ba8**, deleted in 9b0814f; recovered in full, 100 lines, 13 measured false negatives | Fully recoverable from git |
| L | The 61-node rep-transcript row is UNTESTABLE because the s1a FASTA is gone | the file is at `o1_four/pooled.rep.fa`, **61 records**, same keys | Substituting the pooled substrate skipped the claim with teeth on a false premise |
| L | Add a `span > 1 kb` guard to `breakpoint_microhomology` | `is_rt_switch` reads only two short windows ENDING at each endpoint, never the intron between | Rejected as a misread and harmful; wide PSV brackets are already fine |
| L | `assignments.tsv` is self-describing — every decision is re-derivable from the file | **9,139 rows = 0.55%** are unrecoverable | `margin` is printed at 3 decimals while the rule tests `margin > 0` |
| L | Our step-5 famCN split numbers (`cov_split_results.txt`, `surr_mad*.tsv`, …) | **no generator script exists on disk** for any of the 4 artifacts | Not reproducible; renamed `.RETRACTED.txt` |
| L | The SRGAP2 variation-graph panel shows the full family span | the caption overstated the span **~8×**; abpoa asked 512 GiB on a 19 GB box | The panel was windowed (`WINDOW_BP=20000`); the caption now reports the actual window |

## 8.2 TRAP / TAUTOLOGICAL

| R | Claim | Killing number | Why |
|---|---|---|---|
| H | `HSA_genomic.gff` is a usable annotation for O1 truth | drops **17,024 loci = 29.1%** genome-wide, 32.5% of the panel, 92% of HERC2 | Gene-features-only. Use `Reference/chm13v2.0_RefSeq_full.gff.gz` |
| H | `cargo test --lib` is sufficient to validate a change to the O1 path | `--all-targets` **DID NOT COMPILE** for an extended period while `--lib` passed clean throughout | After the fix: 699 passed / 0 failed / 11 ignored over 15 binaries |
| H | Junction positions on a representative's exon-sum are comparable across strands | a '−' rep's exon-sum is stored reverse-complemented ⟹ every +/− pair scores **0.00** | Getting it wrong LOOKS like a real signal; it bit LRRC37 before the negative control existed |
| H | BAM CIGAR junction coordinates are comparable to the pipeline's `exons` column | it reported **0% junction recall and 100% non-canonical motifs** | SAM POS is 1-based, `exons` is 0-based half-open. Cost one full run |
| H | Use `copies.tsv` `n_reads` as the read-support source for the V4b hunt | undercounts on **105/105** copies, median **5.1×**, max 205×; disagrees at the threshold on 29/105 | Two families were wrongly excluded. Use the BAM |
| H | Count secondary alignments as witnesses when computing alt-read fraction | (review CRITICAL) secondaries double-count molecules | A `-F 2308` primary-only witness filter is required before `detect_hidden_copy` (T14) |
| H | The per-chrom recompute reported all units done successfully | a run where **ALL 21 units died** on an unmounted FASTA printed "all units done in 17s" | `rc=$?` was captured AFTER an echo/wc pipeline, and stale `*.copies.tsv` survived for the re-glob |
| M | Setting `RUSTLE_SHARED_EXON=1` alongside `--homology-genomic-span` tests the span substrate | **byte-identical md5** between arms | SHARED_EXON REPLACES the nucleotide runs, so the span sequences are built then discarded |
| M | Swapping the E_r coverage denominator while keeping the 0.50 floor is a valid test | **byte-identical output** — the tell that two edge-altering changes did nothing | The floor must travel WITH the denominator. Cost one 18-min run |
| M | The `additive_genomic_tier_fired` row reports whether the genomic-span leg fired | it is computed as `au_edges_genomic_new > 0` — "did the leg ADD an edge" | Mislabelled on its most quotable row: one file says `fired = false` beside a 230-record PAF on disk |
| M | Region-parallel sweeping (`--region-threads N`) gives ~N× speedup | 11-family contig **1.86×**; 8-contig set **1.19×** | The ceiling is the single heaviest family (~78–80 s), which is already internally parallel |
| M | `--dump-psv` can be used at genome scale | it accumulates all PSV matrices in memory across the parallel sweep | Use with small N |
| M | `<out>.allproj.tsv` can be diffed run-to-run to verify reproducibility | `gw_family_catalog.rs:414` iterates a HashMap with **no sort** before writing | Sorted content is identical, only ORDER differs (pre-existing, not fixed) |
| M | The parallel per-chrom conflict loop speeds up the gw catalog run you are doing | it is only reached by the DEFAULT path (or `--single-copy-baseline`) | `--cross-chrom` uses a different function and `--homology-primary` never reaches it |
| M | A non-empty `rustle.gtf` means the per-chrom job completed | OOM-killed jobs leave truncated GTFs that look complete (last line cut off, no completion line) | Use a `.done` sentinel — and test it with `-f`, not `-s`, because it is 0-byte |

## 8.3 DEAD-END

| R | Claim / proposal | Killing number | Why |
|---|---|---|---|
| H | Keep deleting modules to simplify the codebase | **all 9** remaining removal candidates verified NEEDS-REWIRE-FIRST (o2_materialize trio 2,549 L, asj_*, positional/hidden_copy) | The deletion strategy has hit its floor; remaining wins are CONSOLIDATION |
| M | Run `copy_assign` genome-wide with `collapsed_copies` via `split_locus_copies` | **8/56 regions in 30 min** vs 357 loci / 70 families in **19 s** for the direct detector | Too slow to be usable; replaced by a direct pileup detector on known copy loci |
| M | `--collapse-enumerate` can be used from the `detect_and_assign` / `copy_assign` path | it is inert-with-warning there | It only works on the `gw_family_catalog --homology-primary` path |
| M | `rna_only_edge_oracle.py` is runtime Python that must be ported to Rust | only one python spawn exists in the tree and it is `#[cfg(test)]` | A numpy/scipy feature-AUC harness consumed by NO rust code |
| M | `family_graph`'s `refine_by_minimizer_jaccard` is running in production | **0 production call sites** | Test-only, retired VG-assembler code. If the advisor asks to see it run, nothing is there |
| M | `copy_assign` can be run genome-wide in one process over the full BAM | the first full-BAM run **OOM'd (exit 144)**; the per-contig driver peaks at 2.84 GB on chr1 in 14 min | All contigs + an 11.7 GB BAM + POA over 6 threads exceed the 19 GB ceiling |
| M | Re-plumb the read loops to decode CRAM natively instead of transcoding | noodles' `reference_sequence_repository()` is `pub(crate)` and `records()` is self-referential | 6 read loops would need re-plumbing, and the pipeline reads the input in MULTIPLE passes |
| L | The `FamilyGraph` / `consensus.rs` exon layer can be serialized as the exon-presence graph | its `per_copy_cov` documentation references a function that **does not exist** | DEAD/test-only; v2 built the exon graph fresh |
| L | Wire an inline Rust DNA-support veto loader into the mosaic discriminator | the ref0 catalog covers **2.88%** of the genome and "absent" is unreliable | DNA is offline corroboration only; the prototype is the deliverable |
| L | The re-bundle clone (`RUSTLE_VG_UNION_BASELINE`) splits the over-collapsed mega-bundle | (structural) the clone reuses the already-merged span | A no-op; a real fix must strip secondaries and RE-DERIVE boundaries from primaries |
| L | The POA homology diagnostic is needed in cluster/sweep runs | skipping it: **537 s → 79 s (6.8×)**; it was ~85% of wall-clock | Purely diagnostic — families come from the de-tie conflict graph |
| L | The StringTie assembler stack and the network-flow machinery are part of the thesis | the 41-module assembler island is unreachable from the 5 thesis binaries, pinned by 12 tendril symbols | Explicitly retired as DEAD (`bench/RETIREMENT_MAP.tsv`, 89 modules mapped) |

## 8.4 SUPERSEDED

| R | Claim | Replacement | Why |
|---|---|---|---|
| H | The real Linux disk is `/dev/sde2` (1.8 TB, ~1 TB free) and the VHDX is `/dev/sdd` | reversed: the **VHDX is `/dev/sde`** (reports 757 G free — fiction) and the physical disk is **sdd2, 922 G real** | The VHDX's reported free space is bounded by C:'s free space |
| M | NPIP RNA edge count is 346 (from the panel's stored `summary.json`) | **339** under the current rule | 346 was computed at a coverage bar of 0.45 before `rustlib.py` changed on 08-11. Never mix pre-08-11 artifacts |
| M | The P2 densities printed by `seed_family_report.py` describe the shipped definition | it applied a PANEL-era rule (nmatch/blocklen ≥0.80, pre-M1 coverage) to a SHIPPED PAF | Every P2 density it printed is void; so is any "no alignment record" line from `graph_vs_graph_report.py` |
| M | The original `GGO.bam` is the right substrate for copy assignment | `GGO_mm.bam` — chr19 alone: 114,663 primary + **194,152 SECONDARY (63%)** | The original was the UNDERCOUNT (the N=5-cap problem) |
| M | minimap2's default secondary cap (`-N 5`) is adequate for family graphs | `-N 50 -p 0.1` heals arrays 5 → 11 copies with 0 added FP | The N=5 cap FRAGMENTS real >6-copy arrays. (The 5→11 heal is now non-reproducible — source BAM deleted) |
| L | Per-family `mmseqs` calls are the way to run the protein divergent tier | one batched `easy-search` is **158/158 BYTE-IDENTICAL** and essentially free | ~265 shell-outs made the tier ~3× slower than the whole default run |
| L | `WEB_MAX_DENSITY = 0.15` is the right de-novo web-density gate | raised to **0.30**; density is bimodal with nothing legitimate in [0.15, 0.30) | A real DROP: the 4 newly-dropped are large multi-chrom over-merges (DSFAM0 = 164 ZNF genes over 19 chromosomes) |
| L | The post-M1 control reference contains 24 copies / 6 families | recount: **27 copies / 7 families** | |

---
# §9 FRAMING & SCOPE

## 9.1 REFUTED / RETRACTED

| R | Claim | Killing number | Why |
|---|---|---|---|
| H | O1 is orthogonal to O2 — O1 uses reads only as SUPPORT, never as assignment | `locus_unique_mapper_counts` (read_conflict.rs:267) is reached at denovo_pipeline.rs:3116 **unconditionally**, not behind `--refine` | MAPQ>0 read counts decide how many nodes EXIST. The precise claim is "membership by SEQUENCE alone" |
| H | The read-conflict relation is a subset of segdups (E_c ⊆ E_a) | APOBEC3D/F is in E_c but SEDEF scores it **88.4% < the 90% Bailey cutoff** | E_a is INCOMPARABLE; the only unconditional containment is E_c ⊆ E_b (exonic homology) |
| H | The operational read-conflict relation is a subset of the protein relation (E_c ⊆ E_p) | witness: E_c \ E_p = non-coding/retrocopy read-conflicts | INCOMPARABLE in both directions. Only E_c ∩ {coding-both} ⊆ E_p holds |
| H | O1 uses no read evidence at all on the homology path | two read-derived inputs survive: locus discovery from primaries, and `distinguishing_uniq` | "No conflict graph" ≠ "no reads". `distinguishing_uniq` can never add a member or link two blocks |
| H | Predicted loci are 0.55× the true size — our loci are half the size they should be | spliced rep n=45 **median 0.77, in-band 69%** vs stub rep n=32 **median 0.25, in-band 28%** | Bimodal and entirely a single-exon-representative artifact. The size problem and the stub problem are ONE problem |
| H | "Transcripts that align poorly to the assembly" identifies "copies missing from the assembly" | the two loci with real signal are the bottom two of 913 (**0.9410 and 0.9511** vs median 0.9987) — and both are immunoglobulin | The identity signal and specificity were real; the CAUSE was immunology |
| H | Copy number = minimum path cover | `chi_h` is a greedy CHROMATIC COLOURING of the read-conflict graph | Path cover is the CMCPC junction-augmented generalisation — a different object. Reverted in glossary and figures |
| H | O1 and O2 are views of ONE graph object ("one framework" figure) | (structural) two graphs with different node types: homology over LOCI, variation over SEQUENCE | The PSV matrix, χ_H and copies-as-paths are views of the latter only. Reframed to "TWO graphs, one arrow" |
| H | The E_c read-conflict graph is "the PRINCIPLED family definition" | it contradicted `read_conflict.rs` and DEFINITIONS_FORMAL | The family is the E_r homology relation; conflict COUNTS copies. This docstring was the root cause of the advisor's "inconsistent approaches" complaint |
| M | The advisor's objection applies: T2T is an average, so any absent-copy signal is confounded | hom-alt drop **13.50×** with a het internal control at **1.096** | mGorGor1 is a haplotype-resolved assembly of ONE animal and the IsoSeq is the same CELL LINE. Live only for the CROSS arm and human work |
| M | The ID_63 flagship panel's red members illustrate the K=0 identifiability floor | those reds were "distinguishable-but-MERGED" — a genuine METHOD miss | Mislabelled as the floor; caught only by cross-referencing `soto_floor_decomposition.tsv`. Panel swapped to amylase |
| M | The cDNA all-vs-all + protein-homology pipeline (856 families) IS the family definition | the read-based definition is **196 families / 1,136 de-novo loci**; the 856-family pipeline is annotation-side | The definition uses ONLY IsoSeq reads. Once cDNA/genome homology became the backbone substrate, even orthogonality collapsed |
| M | The only hand-rolled k-mer steps are `family_detect::candidate_pairs` and `family_graph` | the real live sites are `family_rescue::canonical_kmer_set` (KMER=18) and `multi_repeat_bridge::locus_node_set` | Exactly ONE live k-mer site per default binary path. Wrong in an advisor-dangerous way |
| M | The flagship ASJ hits PSMD2 and DAXX are textbook splice-site variants | **0/475** anchors sit on a core GT-AG dinucleotide; PSMD2's donor GT and DAXX's minus-strand donor are INTACT | They are splice-REGION variants; the load-bearing result is the per-molecule allele→junction LINKAGE |
| M | PSMD2 and DAXX are the exemplars to lead the ASJ result with | PSMD2 **SOR = 10.45** (all 14 junction reads forward-strand), DAXX 7.08, both failing SOR<3.0 | A third independent reason to retire them. Lead with the genetic-core SET and DAZ1-vs-DAZL |

## 9.2 TRAP / TAUTOLOGICAL

| R | Claim | Killing number | Why |
|---|---|---|---|
| H | Use the shared-read multimap certificate (171/171 NPIP pairs) as the family EDGE rule | the census kept secondaries (`-N 50 -p 0.1 --secondary=yes`) and **97.4% of reads have a MAPQ=0 alignment** | E_c is the ambiguity oracle and belongs to O2; using it to build families destroys O1⊥O2. "Linked" means the aligner could not separate them |
| H | Randomized/sampled instance checking is adequate to validate an identifiability condition | sampling **missed the violating instances**; only exhaustive enumeration discovered Strong Separation | It gave false confidence |
| H | Use facility location (pay to open a site) to decide which loci exist | (architectural) facility location is O2 machinery | Using it to build loci imports assignment into O1 and breaks O1⊥O2. Explicitly marked "Do not." |
| H | Write up the DNA-side result as "0/915" — a null result about the method | Yoo/Rhie's 1–2 Mbp collapse per haplotype over the **1.1224%** these windows span predicts only **0.47–0.94** collapses | Underpowered BY CONSTRUCTION. Four bounds must travel with it (tandem-only 93.3%, ≥97.3% id over ≥5.8 kb, probe visibility ~0.17) |
| H | Promote a coherent pile of unmapped/unaligned reads to a reference-absent FAMILY | a genuinely absent copy yields a coherent unmapped pile in only **1/3 of cases** (54/162 ORPHANED) | Our nodes are genomic intervals: a read-only cluster has no interval, no λ, no seed-invariance. Both retracted candidates were exactly such clusters |
| H | We replicated Soto 2025 / this validates our RNA method | **3 of the 4 legs consume Soto's own files** (SD98 from their SEDEF track, genes from their CAT v4, famCN from identical WSSD windows) | Only the clustering computes anything, and no held-out set exists. Say "re-implementation concordance" |
| H | Present the variation graph as a DETECTOR of families | (structural) the VG only REPRESENTS what is given to it | Circular. The locked honest framing is "identifiability ceiling" |
| H | Pitch the method as "zero k-mers / all graph-to-graph" | the repeat-bridge gate needs a genome-wide repeat-FREQUENCY statistic ("recurs ≥8× genome-wide") | Alignment gives similarity, never frequency. Honest pitch: corner k-mers to the ONE irreducible role; "sequence-to-graph, NOT graph-to-graph" |
| H | Facility location / LP-rounding (the advisor's own 2016 formulation) is in the codebase | **0 occurrences in `src/`** | THEORY-ONLY; the bench check that exists uses a copy-INDEPENDENT weight. Do not sell "threading = graph optimization" |
| H | Frame the α-gate's abstention as DERIVED from an iterative joint estimator | the objective decomposes ⟹ per-read argmax IS the optimum; the EM changed **0/3,081** decisions | Explicitly forbidden framing. Say instead: an EM was run against the shipped rule and changed no evidenced decision |
| H | O1 and O2 can be presented as one decision rule | (scope) O1 = identification via E_r; O2 = assignment via E_c and the `de` criterion | E_c does NOT decide family membership. Flagged by the user as exactly what would irritate the advisor |
| H | Objective numbering can be quoted directly from any project document | "O3" is the EM in one document, allele-specific junctions in another, reference-absent copies now | Three incompatible schemes (2026-06-03 five objectives, 06-25 six, 08-07 three). Resolve the scheme before quoting |
| M | The IG/somatic-hypermutation confound is a novel trap our work discovered | Rodriguez 2021 (PLoS ONE 16(12):e0261374): V(D)J somatic deletions average **464 kb**; IGLoo removes somatic-haplotype reads explicitly | A known, named, guarded-against failure mode. Cite it, do not claim it |
| M | The 475 ASJ calls have collapsed-paralog masquerade RULED OUT | the control removed **0/475** | Non-binding: a control that removes nothing has excluded nothing (T20), and the separator was never run on the 44 LOC* calls |

## 9.3 DEAD-END

| R | Claim / proposal | Killing number | Why |
|---|---|---|---|
| H | Delete the unique-mapper-count branch in `distinct_locus_reps` to restore O1 ⊥ O2 | measured inert (0/15, 0/94 pairs, 0 firings over 451 chr1 loci) yet **not removable** | Its real domain is a co-located same-strand junction-less rep pair (true K=0), where nothing in the sequence separates them. Scope and monitor |
| H | Build O1 loci by max-weight bipartite assignment / facility location over all alignments | (architectural, agreed with the user) it makes O1 depend on per-read ASSIGNMENT | Support is a RELATION, assignment is a FUNCTION. The O1-safe answer is to abstain: keep a site if ≥3 reads share an intron chain by ANY alignment |
| H | Build one VG per family and recover copies+isoforms by constrained flow decomposition | Guitart/Eichler Fig S3: "the majority of TBC1D3 copies were **isolated nodes** with single-haplotype support" | Eichler's group built exactly this and ABANDONED it for phylogenetic groups; AS≥10 + phylogeny did the assignment |
| M | Discovering a reference-absent copy from transcriptome data is an established approach | **0 of ~30** primary papers; the closest analogue discovered copies by RE-ASSEMBLY | Accepted standards are S1 re-assemble / S2 depth+PSV / S3 unique peptides; transcript anomaly is not accepted evidence |
| M | The `--read-coherence` recall layer (the biggest measured win, +2,248 FSM) should headline | (scope) the advisor is explicitly not interested in a better assembler | Keep it only as the substrate producing per-molecule paths |

## 9.4 SUPERSEDED

| R | Claim | Replacement | Why |
|---|---|---|---|
| M | The thesis is a family-variation-graph co-assembly method | **three objectives, explicitly NOT an assembler**: O1 definition, O2 assignment, O3 detect+flag | Re-agreed with the advisor 2026-08-07. No flow decomposition or co-assembly remains in scope |
| M | The allele-specific-junction objective is a thesis deliverable and the headline | dropped by agreement — even though the 06-25 audit rated ASJ the **ONLY cleanly-attained** objective | De-scoped, not deleted. Do not present the reduced scope as though nothing was given up |
| M | 120 transversion-anchored ASJ calls are the unambiguously genetic headline set | genetic core **~77** (44 of the 120 are copy-confounded LOC* calls); airtight set **≈20** splice-proximal | "120" conflates not-edit with not-copy |

---

*End of register. Entries without a killing number were not admitted.*
