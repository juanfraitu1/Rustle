# Coupled Flow-Selection + Downstream Port — Design

**Status:** In execution (2026-06-11). Builds on the ST-faithful flip (6 gates shipped, chr19 291→264). Revert point: `RUSTLE_PRECISE=1` (byte-identical to 4705ab1). Goal: drive chr19 rustle-vs-ST divergence 264 → lower by porting flow-selection + downstream filtering *together* (the documented coupled wall).

## 1. The landscape (chr19, 264 divergent chains / 138 loci)

Built via exhaustive per-locus instrument-and-diff (23-agent workflow `wf_896f7150-f9e`, full classification in `/tmp/loci/classification.json`). Both tools' chr19 path/selection parity logs captured once (`/tmp/{st,ru}_chr19.jsonl`), sliced per-locus (`/tmp/loci/<idx>.txt`), classified by mechanism + portability.

| Category | Loci | Prize (chains) | Portability |
|---|---|---|---|
| GRAPH_DIFF | 51 | 111 | 30 loci/74 = WALL (junction-acceptance), 21 loci/37 coupled |
| CHECKTRF_ROUTING | 31 | 61 | ~all coupled |
| ALT_SPLICE_COMBO | 6 | 21 | coupled |
| RETAINED_INTRON | 15 | 18 | 11 coupled / 4 local |
| TERMINAL_OVEREXT | 11 | 17 | 10 coupled / 1 local |
| INCLUDED_DROP | 10 | 16 | 6 coupled / 4 local |
| FLOW_SWAP | 10 | 16 | coupled |
| OTHER/COV_SCALE | 4 | 4 | local |

**Wall (~82 prize / 32 loci, 31%):** GRAPH_DIFF rooted in junction acceptance — rustle accepts/rejects micro-exon splits, fuzzy 3-12bp donor/acceptor shifts, extra terminal/cassette nodes that ST doesn't. Matching = the documented precision catastrophe (rustle good-junctions not a clean subset/superset; st_shadow regresses; flow-port at ~12-chain ceiling). **NOT addressable at this level.**

**Addressable: ~182 prize / 106 loci** (165 coupled, 17 local).

## 2. The unifying root: cov-scale + checktrf routing

Two coupled roots drive most of the addressable divergence:

**(A) cov-scale mismatch (the dominant root).** ST's `pred->cov` is **tlen-proportional** (a long chain with few reads gets HIGH cov); rustle's flow-cov is **longcov/flux-proportional** (same chain gets LOW cov). Documented as "ST cov 60x higher" (`project_coverage_metrics_deviation`). This makes the SAME pairwise predicate diverge:
- **Verified at idx8** (graph-identical): ST `included_drop` fires because `n1.cov(6172) > n2.cov(7738)*DROP(0.5)=3869` ✓; rustle fails because flow-cov gives `n1.cov(1.02) > n2.cov(2.04)*0.5=1.02` ✗ (by 0.0006 — exactly 0.5× because lc1 vs lc2). The agent's "killer=None unconditional-containment branch" diagnosis was WRONG: `killer=None` is a payload artifact (ST's pairwise `included_drop` emit at rlink.cpp:19010 omits killer fields). Root is cov-scale, not a missing branch.
- Same root underlies **readthr** (rustle flow-cov 0.93-0.99 on chains ST scores >1.0 → readthr<1.0 floor misfires; idx77/80/91/134 pure cases) and **RI** thresholds.

**(B) checktrf routing + store-breadth.** ST routes alt-terminus/exon-skip chains via transfrag-bounded `get_trf_long` checktrf (terminating at read 3'); rustle's flow over-extracts them (locus-B mechanism) OR rustle's `checktrf_rescue` manufactures 8-65 surplus short/truncated fragments vs ST's 0-14. Coupled with downstream survival (the rescue is recall-load-bearing — blind drop over-kills 156 ST-shared isoforms per `project_post_flow_gate_pin`).

## 3. Ordered targets (synthesis `wf_896f7150-f9e`)

0. **checktrf_rescue store-breadth + survival** (prize 56, MED): clamp store-count to ST's `get_trf_long -L` breadth AND port RI/included_drop survival exemption (coupled — neither alone works).
1. **readthr flow-cov vs tlen-cov scale** (prize 27, MED): feed tlen-scaled cov to the readthr predicate (not lower the threshold).
2. **included_drop cov-scale** (prize 16, was mislabeled LOW→actually cov-scale MED): the `n1.cov > n2.cov*DROP` direction needs ST's tlen-cov.
3. **flow intron-combo selection order** (prize 25, HIGH): deepest flow-enumeration layer; prior re-ports regressed; defer.
4. **terminal extension termination + killer-coord cascade** (prize 17, MED).
5. **RI cross-locus killer + lowintron(bpcov)** (prize 11, MED): ST RI-kills using an external container whose intron is lowintron-flagged from **bpcov**; rustle's flow-cov never sets the flag / doesn't pool the cross-locus container. Distinct from cov-scale; well-localized (predcluster_st.rs).
6. **duplicate-component dedup** (prize 8, LOW): rustle processes some components twice (361 graphs emitted 2×, 56 with >3 nodes) — graphnode_list emit is in the main bundle path, so genuine double-process. Impact mechanism (halved flow → readthr misfire) needs bundling analysis.

## 3a. VIABILITY PROVEN (2026-06-11 experiment)

The coupled port is **viable** — demonstrated decisively. Implemented the readthr cov-scale exemption alone (lower `readthr_longcov_min` 2.0→1.0 under ST-faithful default, exempting longcov≥1 read-backed multi-exon flow chains from the readthr<1.0 kill). Result on chr19:
- **ST-only 76→65** (recovered 11 — the readthr exemption is correct ST-faithful behavior; ST keeps these read-backed chains).
- **RU-only 188→270** (+82 surplus) → total 264→335, net-negative ALONE. **Reverted.**
- **Decisive cross-check:** of the +82 surplus RU-only, **91% (75/82) are chains ST ALSO EXTRACTS but kills downstream** (matched against ST `path_extracted`); only 7 are chains ST never extracts. → Porting ST's downstream kills would remove ~75 and keep the 11 → **net convergence ~264→253 reachable** from this one coupled pair.

**Refined cov-scale understanding:** ST's `pred->cov` IS per-base normalized (`pred->cov /= abs(tlen)` at rlink.cpp:10858; the `path_extracted` cov ~2071 is *pre*-normalization). So the scale is closer than raw logs suggest; the readthr divergence is a boundary effect (rustle 0.95 vs ST ~1.0). **Filter order is the coupling:** rustle runs pairwise RI/included_drop (transcript_filter.rs:3210) BEFORE the readthr gate (8404), so the 75 surplus already survived rustle's pairwise filters — ST's pairwise RI/included_drop fired on them, rustle's didn't. The coupling lever is the **pairwise downstream comparison**, not readthr itself.

## 3b. COUPLED-PAIR ATTEMPT — inline shortcut FALSIFIED (2026-06-11)

Tried to implement the coupled pair as a single inline decision: readthr exemption + an inline containment guard (exempt a read-backed longcov≥1 chain from readthr ONLY IF it is not an intron-subset of a higher-cov sibling — RI & included both mean victim.introns ⊆ killer.introns). Result chr19: ST-only 76→66 (recovered 10), RU-only 188→**247** (+59 surplus still leaked) → total **313**, still net-negative. **Reverted.**

**Why it can't be shortcut:** the strict-subset guard caught only 23 of 82 surplus. The other 59 are NOT strict subsets, and cov-fraction does NOT separate them either (surplus median 0.049 vs recovered median 0.075, fully overlapping; recovered range 0.017–0.384). **The 11 ST-keeps and 75 ST-kills are feature-indistinguishable at the readthr-exemption point.** Critically, rustle's pairwise RI/included_drop runs BEFORE readthr and already returned "don't kill" for these 59 — so any inline re-check at readthr re-runs the SAME diverging decision. **The real gap is rustle's pairwise RI returning FALSE where ST's `retainedintron` returns TRUE** (confirmed example: 16-intron victim cov 0.98, killer cov 17.18 — ST RI-kills `0.98<0.1×17.18`, rustle keeps).

**Conclusion:** the coupled pair MUST be implemented as **pairwise RI/included_drop decision parity** (make rustle's `retainedintron_like` + killer-pairing + lowintron match ST's `retainedintron`), NOT an inline readthr guard. The divergence is in: lowintron flags (bpcov-undercount, transcript_filter.rs:1550), the killer pairing/predcluster grouping, or the cov comparison branch. Diagnose via the cross-tool RI instrumentation already in place: ST `pred_ri_eval` (every eval, ri_result + killer_nlow) vs rustle `RUSTLE_RI_TRACE` — diff ri_result per (victim,killer) pair on the 59 surplus. Then the readthr exemption becomes safe (the surplus dies in pairwise, the 11 survive). This is the multi-session pairwise port; the readthr exemption alone is the proven-correct flow-side half waiting for it.

## 3c. LOWINTRON ORACLE — pins the downstream half, bounds the prize (2026-06-11)

ST emits the exact per-intron lowintron mask via `pred_intron_low` (rlink.cpp:18536, `intron_low`+`chain`). Built an ST-exact lowintron oracle (4023 masks) and fed it to rustle via the wired `RUSTLE_LOWINTRON_ORACLE` hook (transcript_filter.rs:1551/1665; chain-key conventions match: rustle `{end+1}-{start}` == ST `{end+1}-{start-1}`). Three configs:
- **A. oracle alone:** chr19 264→**258** (RU-only 188→184, ST-only 76→74). ST-exact lowintron is a **clean net-positive −6** — the bpcov-undercount (rustle flags fewer low-introns than ST, e.g. example-2 killer rustle 3 vs ST `killer_nlow=7`) really does let surplus survive rustle's structural middle-exon RI (`ri=2`).
- **B. oracle + readthr exemption:** 327 (worse). Even with ST-exact lowintron the readthr surplus still floods +76 — lowintron removes only ~6 of the 82. So the readthr surplus is killed by ST via **isofrac/included_drop, NOT lowintron-RI**. The readthr exemption is unsalvageable as a standalone half.
- **C. cheap `RUSTLE_RI_ZERO_INTRON_LOW` approximation:** 267 (worse) — over-fires, loses 6 ST-only. The −6 needs *exact* lowintron, not the blunt zero-intron reinterpretation.

**Decisive conclusions:**
1. The **readthr+RI coupled pair does NOT yield a clean net-positive** — the surplus is feature-indistinguishable at the readthr point AND survives correct lowintron (killed by ST via isofrac/included_drop, a different cascade).
2. The one genuine net-positive lever found is **lowintron parity, bounded at −6** (264→258). It requires aligning rustle's **bpcov intron-coverage computation** with ST's (rustle computes exact 0.0 where ST computes tiny-nonzero → undercounts low-introns; transcript_filter.rs:1606-1618). The cheap approximations overshoot; the precise fix is the bpcov-computation alignment (the documented cov-computation lever), validatable against the oracle's 258.
3. The bigger addressable prizes (checktrf-routing 61) are untouched and remain the higher-leverage next targets — but they share the same cov-computation entanglement.

**Net for "what is possible at this level":** the coupled pair is more entangled than the 91% viability proof implied — the downstream half diverges from ST in the underlying **bpcov/cov computation**, not just in surgical filter rules. The clean, bounded win available is the lowintron-bpcov port (−6, oracle-validated). Larger gains need the coverage-computation alignment, which is the multi-session cov-scale port.

## 3d. LOWINTRON-BPCOV PORT — strand hypothesis FALSIFIED; root is the coverage substrate (2026-06-11)

Pursued the −6 lowintron win (oracle target 258). Found `build_lowintron_flags` (transcript_filter.rs:1575) receives a STRAND-SPECIFIC bpcov (`use_plus=strand!='-'`, pipeline.rs:13879) while ST reads `get_cov(1,...)` = the all-strand layer. Hypothesized strand was the root and fed the all-strand bpcov (`bpcov_stranded.all`) to lowintron (via thread-local) — **FALSIFIED: chr19 = 265, not 258** (recovered 1 ST-only but over-flagged +2 RU-only). All-strand for the whole pairwise = also 265.

**Measured the actual divergence** (ST `pred_intron_low` for the example-2 killer 23290512-23306475): ST flags intron 0 because its `introncov = 0.002315` — a **tiny non-zero from ~1 base of intronic coverage (a single read overhang)**; introns 1-3 are exactly 0 in ST too (and unflagged, same as rustle). So the divergence is NOT strand and NOT the exact-zero guard — it's that **ST's per-base coverage picks up sparse single-base intronic coverage (read overhangs / soft-clips near junctions) that rustle's bpcov accumulation drops** (rustle delta-encodes per-aligned-exon, contributing exactly 0 in clean introns).

**Definitive conclusion for the whole coupled port:** every addressable lever funnels through the **coverage substrate** — rustle's read→bpcov / flow-cov accumulation diverges from ST's at a fine grain (single-base intronic overhangs here; tlen-normalization and strand-layering elsewhere). Surgical filter/exemption ports cannot capture it (readthr coupled pair falsified §3b; lowintron strand fix falsified here). The real port is **aligning rustle's coverage accumulation with StringTie's**, which is a deep, broad-blast-radius change (bpcov feeds the whole pipeline) with modest per-lever payoff (−6 for lowintron). That is the genuine next-level work; it is NOT a single-session surgical win. The 6 shipped gates (chr19 264) remain the durable ST-faithful baseline; the oracle (`/tmp/st_lowintron.oracle`, target 258) is the validation harness for the eventual coverage-accumulation port.

## 3e. COVERAGE-ACCUMULATION PORT — fully diagnosed; two-part, broad-radius (2026-06-11)

Started the coverage port. Dumped rustle's COMPUTED lowintron masks + introncov genome-wide (`RUSTLE_LOWINTRON_DUMP`, single-threaded) and diffed against ST `pred_intron_low` (mask + `intron_covs`). Findings:

1. **When rustle has intronic coverage, it matches ST EXACTLY (ratio 1.000).** Not a scaling/tlen issue. The entire divergence is rustle computing coverage where ST has MORE / rustle exactly 0 where ST nonzero. Baseline: rustle UNDER-flags 1923 introns (96% are rustle exact-0), over-flags only 122. Masks differ on 1192/3334 shared chains.

2. **Part A — strand (PROVEN).** Root case at locus A (mini3): intron 17171361-17171490 in a `+` prediction has ST introncov=2.0, rustle 0. The 2 reads covering it are **flag 16 = MINUS strand** — antisense retained-intron reads. Rustle's bpcov here is strand-specific (`use_plus`, pipeline.rs:13879) so it drops them; ST's `get_cov(1,...)` = all-strand `bpcov[1]` counts them. Feeding the all-strand layer (`bpcov_stranded.all`) to lowintron improves mask agreement **1192→952** disagreements and raises introncov-match to 92%.

3. **Part B — base-level CIGAR accumulation (residual).** Even with all-strand, ST's introncov stays HIGHER than rustle's in **2840** intron cases (mostly small, e.g. ST 0.035 vs rustle 0; median ~0.02 = a few bases). ST's per-base coverage counts intronic bases rustle's per-aligned-exon delta-encoding (`Bpcov::from_reads`, bpcov.rs:222) drops — read overhangs / soft-clips / boundary bases near splice sites.

**Net:** all-strand alone = chr19 **265** (net −1 worse): it fixes ~240 chains' masks but the 6 oracle-relevant chains need BOTH parts, and all-strand introduces +2 over-flags. Reaching the oracle's 258 requires the strand change AND a base-level read→bpcov accumulation rewrite matching ST's CIGAR handling — **broad blast radius (bpcov feeds the whole pipeline) for a −6 ceiling**. Even the strand half alone is net-negative.

**VERDICT:** the coverage-substrate divergence is now precisely characterized and proven (strand + base-level CIGAR), with the diagnostic harness in place (`RUSTLE_LOWINTRON_DUMP` vs ST `pred_intron_low`; oracle target 258). But the full port is a major, broad-blast-radius coverage-accumulation rewrite whose payoff (−6) is small and whose partial forms are net-negative — it is not a clean incremental win and risks regressing the 1745 in-both chains. **The 6 shipped gates (chr19 264) remain the durable ST-faithful baseline.** Pursue the coverage port only if byte-exact StringTie parity is a hard requirement; then port `Bpcov::from_reads` to count intronic-overhang/soft-clip bases like ST's `add_read_to_cov`/`get_cov` AND use the all-strand layer for lowintron, validating each step against the dump + oracle.

## 3f. COVERAGE PORT EXECUTION — layer-by-layer; Part A solved, Part B = read-set completeness (2026-06-11)

Pursued the coverage-accumulation port (strict-parity requirement). Worked the layers with the `RUSTLE_LOWINTRON_DUMP` (computed mask + introncov, single-threaded) vs ST `pred_intron_low` (mask + `intron_covs`) harness:

- **Layer 1 — strand (SOLVED).** Root: rustle splits reads into strand-specific bundles, so a prediction's bpcov never sees opposite-strand (antisense) reads (proven: locus-A intron 17171361-17171490, the 2 covering reads are `ts:A:+` antisense in a `−` prediction). FIX: stash an all-strand bpcov built from `region_reads` (3-strand groupflow holds both strands) for `build_lowintron_flags` only (pipeline.rs:~15151 → thread-local → transcript_filter.rs). Introncov for antisense cases now matches ST (0→2.0/1.0). Mask disagreement 1192→995, introncov match ~86–91%.
- **Layer 2 — window/boundary (RULED OUT).** Tested introncov over left-shift / include-donor-base windows: the CURRENT window matches ST best (86% vs 9%/8%). The `cov_edge_add(...+1)` convention does NOT cause a net offset; the query region+divisor already match ST exactly.
- **Layer 3 — read weight / NH (RULED OUT).** Residual intron 32068767-32071937 (ru 2.607 vs ST 2.918): all 615 region alignments are PRIMARY (no secondary/supplementary), the M-covering reads are NH=1. Not multimapper-weight, not filtered-flag.
- **Layer 4 — read-set completeness (THE REMAINING LAYER).** The residual is ~986 base-coverage-units ST has that rustle lacks in that intron — i.e. **reads overlapping the bundle region that `region_reads` (keyed by exact `(chrom, b.start, b.end)`) does not include.** ST's bpcov is over ALL reads overlapping the genomic bundle; rustle's per-bundle/exact-span-keyed read set is a strict subset.

**Status:** Layer 1 implemented (correct, but net −1 at chr19 265 ALONE because Layers 2-4 unresolved → the 6 oracle chains need all layers). Layers 2-3 ruled out. **Layer 4 (read-set completeness) is the precise remaining work:** build the lowintron bpcov from ALL reads overlapping `[graph_bundle.start, graph_bundle.end]` regardless of bundle — e.g. a per-chromosome all-strand per-base coverage from the full read set (built once before bundling, queried in `build_lowintron_flags`), matching ST's genomic-bundle bpcov. Reverted to clean baseline 264. The harness (`RUSTLE_LOWINTRON_DUMP` vs `pred_intron_low`, oracle 258) is in place to drive Layer 4. This is a bounded but real architectural addition (a region/chromosome coverage structure); once it lands, Layers 1+4 together should reach the oracle's 258, then validate full suite + escape-hatch.

## 3g. LAYER 4 EXECUTED → ARCHITECTURAL CONCLUSION (2026-06-11)

Implemented Layer 4: a per-chromosome ALL-STRAND, ALL-READ per-base coverage built from every bundle's reads (`Bpcov::all_strand_from_read_refs`, prefix-sum pre-built, `Arc`-shared, stashed per-bundle for `build_lowintron_flags`; ~2.5 GB peak, fine). Result: chr19 **still 265**, NOT 258. The dump showed the chromosome-wide coverage **over-corrects**: introncov now ST<rustle in 2690 cases (was the reverse), i.e. it counts reads from neighbouring/cross-component bundles that StringTie's per-bundle `bpcov` excludes. And the original residual intron 32068767-32071937 is **unchanged at 2.607** (ST 2.918) even with ALL bundle reads — its missing coverage is **not in any rustle bundle** (pre-bundle-filtered), and rustle's read weight is `YC/YK ÷ NH` (bam.rs:624-629), another tag-interpretation layer.

**The full descent (each layer fixed/ruled-out revealed the next):**
1. strand (all-strand reads) — SOLVED, revealed →
2. window/boundary — RULED OUT, revealed →
3. read-weight at midpoint (NH=1) — RULED OUT, revealed →
4. bundle-completeness (chromosome-wide) — OVER- *and* UNDER-counts vs ST's per-bundle read set, revealed →
5. read-filtering + YC/YK/NH weighting — reads not in any bundle, differently weighted.

**CONCLUSION (architectural — systematic-debugging Phase 4.5, 3+ fixes each exposing a new layer):** StringTie's lowintron `get_cov(1,...)` reads its **per-bundle bpcov**, which is the product of ST's *entire* read pipeline — read loading, YC/YK/NH weighting, strand inference, and connected-component bundle definition. Rustle differs from ST at **every** one of these layers. Matching ST's exact lowintron coverage therefore requires **byte-exact parity of the whole read→bpcov pipeline** (loading + weighting + bundling), not a targeted coverage patch. That is a massive, broad-blast-radius undertaking for a **−6** chr19 payoff, and partial forms over/under-shoot (net-negative). **Verdict: the lowintron −6 win is not achievable without full read-pipeline parity; it is out of scope as a targeted fix.** All Layer-1..4 code reverted; baseline 264 clean; escape-hatch byte-identical. The oracle (258) and dump harness remain for any future full-pipeline-parity effort. **The 6 shipped gates (chr19 264) are the durable ST-faithful baseline and the recommended stopping point.**

## 4. Strategy

The cov-scale root (A) is central but **known-hard/reverted** (changing cov scale perturbs every threshold). The coupled approach: rather than change the EMITTED cov, port **ST's cov semantics into the filter-comparison basis only** (a tlen-scaled comparison quantity), gated behind `precise_mode()`. Test each target on the cheapest pure-graph-identical loci first, validate chr19 (`< 264`), revert on regression. Discipline unchanged: one change at a time, escape-hatch byte-identical, validate genome-wide before keeping.

## 5. Methodology + artifacts

- Work-list: `bench/divergent_loci.py` (cluster), `bench/slice_loci_digests.py` (per-locus digests).
- Validate: `target/release/rustle -L GGO_19.bam` vs `/tmp/st_all.gtf`, `bench/gtf_chain_diff.py` (baseline 264).
- Escape-hatch: `RUSTLE_PRECISE=1` byte-match 4705ab1.
- ⚠ Verify each agent diagnosis before implementing (idx8 lesson: over-optimistic LOCAL labels).
