# Copy assignment (O2): definition, decision rule, abstention, certificate

**O2's domain (user, 2026-10-19): AS-tied multimapping reads only**: a read with a rival alignment within 10 AS of its best (Eichler's within-10 rule; under splice:hq, 10 AS is two substitutions). A read with a clear best placement is already placed by minimap2; applying O2 to it is tautological. O2 goes beyond minimap2's coin toss among AS-tied alignments.

Compacted 2026-11-01. Full text of every section: `docs/archive/2026-11/copy_assignment_definition.md` (verbatim) and git history. Each entry ends with a pointer `(full: archive/2026-11/copy_assignment_definition.md:LINE)`.

## Current state (2026-11-03)

- **Scope (user, 2026-11-03):** reads minimap2 could not place = runner-up within 10 `ms` of the best, any AS, minus the K = 0 wall. Read a near tie on BOTH scores: ms and AS are one unit per record but a molecule's two margins differ by splice shape (r1508). The AS-gap class (spliced primary, unspliced higher-AS rival) is the aligner's (r1403/1404).
- **Two O2 code paths.** `copy_assign` (read-star PSV certificate, §1-§6) and the pile-based `copy_phasing run`. Default `copy_assign`: AS-tied gate at exact tie, `best_by_psv` then `best_by_duel`, Poisson-binomial certificate at alpha/(n-1), alpha 1e-3, genomic origin certificate, `--gtf-copy-set`.
- **Opt-in (measured 2026-11-03):** `--tie-ms W` adds the ms gate (`ms_best`, `ms_second`, `ms_margin` columns only under the flag; never a tie on ms 0; refuses an unindexed BAM). `--primary-on-family` (ON under `--tie-ms` since r1542, 2026-11-05; `--no-primary-on-family` is the escape): a tied molecule enters the gate only with a primary record among the family's windows; otherwise its best placement is outside the family and the PSV model, which knows only catalog copies, cannot test it (r1514: removes 3,465 of 4,529 held-out assignments at no cost to what was right). Older opt-ins: `--union-certificate` (r1103, r1481 FAIL), `--tie-structure`, `--certify-all-records`.
- **Truth and numbers:** KB3781 = Jim; `haplotype_truth.py` + `strict_truth.py` at slack 0. Panel (40 families): 284 right / 2 wrong vs primary 280 / 6 (r1509, post hoc truth). Held-out (454 families) with the guard: 129 / 22 vs 125 / 26; O2 right where the primary is wrong 5, the reverse 1 (r1514). Unverifiable rows (identical loci over the read) are their own class.
- **Failure classes seen:** (1) outside-catalog origin (uncatalogued parent gene; MCL148, MCL13) -> the guard; (2) shape mismatch (spliced read assigned to a single-exon copy; same families) -> a shape guard, not yet registered; (3) allele vs paralog (MCL183: the mat allele of copy 0 fits copy 1's reference at the PSVs; r1413 trap) -> needs haplotype-aware candidates.
- **Report, no default change:** r1484 PSV VCF (`bench/o3_graph/psv_vcf.py`; bcftools agrees 54.6% as run, 93.9% with MAPQ reset).
- **Prohibition (standing):** never claim O2 assigns better than minimap2 at scale; defend O2 on abstention and on not assigning reads the catalog cannot explain (§0, Reassignment §7).

## State as of 2026-10-07 (read this before the rest)

Body below last committed 2026-09-01, describes the binary of 2026-08-10/15; defaults changed 2026-09-09. Human (A119b, CHM13) and gorilla (KB3781/OR6737, mGorGor1) numbers are never pooled.

**(a) Shipped O2 at 2026-10-07.**
- Code: binary `src/copy_assign_bin.rs`; library `src/family.rs` `pub mod copy_assign` (`AssignParams`: `alpha: 1e-3`, `error_rate: 0.003`, `junction_weight: 5.0`, `origin_drop_indels: true`, `best_by_psv: true`, `best_by_duel: true`), `copy_assign_pipeline`, `copy_split::phred_err` (clamp `[1e-4, 0.25]`). Old paths `src/rustle/vg_family/...`, `src/bin/copy_assign.rs` no longer exist.
- Input = `families` roster `<id>.fam.copies.tsv` (`mcl_families --from-gtf --min-exonic-bp 1 --min-shared-exon-frac 0.60`); `gw_family_catalog` is LEGACY.
- Step 1 AS-tied gate (default ON 09-09): only molecules with >= 2 placements whose runner-up AS >= `--as-tie-ratio` x best (default `1.0` = exact tie); everything else dropped before the certificate. Ties outside every family are registered (`tie_outside_catalog`).
- Step 2 best copy = `best_by_psv`, then `best_by_duel` (escapes `--best-by-alignment`, `--no-best-by-duel`); r799: origin certificate is bk-only.
- Step 3 certificate: pairwise exact Poisson-binomial tail at `alpha/(n-1)` plus a GENOMIC origin certificate ignoring I/D (`--no-origin-drop-indels`); a read disagreeing with every candidate beyond error is `origin-rejected` ("matches no catalog copy"). `sole_candidate` ON.
- Step 4 abstention is a PSV-count step: abstain at m <= 1, assign at m >= 2 unless the read's bases conflict; at m = 1, p = e/3 = 0.001 = alpha (r1175, `docs/PREREG_o2_power_2026-09-29.md`).
- Step 5 outputs: `--gtf-copy-set` ON (undecided isoforms emitted once with `copies "A,B[,outside]"`); opt-in `--union-certificate` (r1103).
- Operating numbers (contested set only): human A119b MCL0 1,118 contested = 262 assigned / 512 tied / 344 ambiguous, excision sweep 247/262 = 94.3 % abstain; gorilla MCL1 33 = 4/28/1, MCL7 11 = 0/10/1, MCL58 40 = 2/38/0. Simulated truth human chr16: OWN 163/163, gorilla 0/0.

**(b) Superseded or refuted claims of the 09-01 body:** §0 22 % near-tie population (gate admits exact ties only); §1 argmax-over-logl with lowest-index tie-break (now PSV best + duel); §1 three statuses (`origin-rejected`, `tie_outside_catalog` added; r786-790); per-family decomposition refuted across families (no cross-family arbitration; chr16 sim PRIMARY 0.29 / ANY 0.09 vs aligner 0.49); alpha plateau measured pre-gate only; §2(a) every-record pool (now AS-tied only); §2(c) catalog input (now families table); §5 file:line refs dead; §7/§12 verification predates gate (not re-run); §9.1 `tie_invariant` nearly vacuous under the gate (r786, r1092, r1109); §9.3/§10 `recombinant_abstain` never connected; §11.6 primary-in-copy denominator superseded; §11.4 contradiction-rule mechanism superseded (dominant abstention now `origin-rejected` below 99.5 % identity, `docs/PSV_CEILING_2026-10-03.md`); §0 "no ground truth" partly superseded by simulated truth (HG002 real-read attempt covered 28 % of tied reads, C0 failed).

**(c) Open items (2026-10-07):** r1220 `--tie-outside-at-best` prereg not run (`as_tie_ratio` text 0.98 vs CLI 1.0); default-roster simulations UNDERPOWERED (r1252 chrY); cross-family arbitration (strict union certificate assigns 0 of 1,259 chr16 sim reads) vs expression prior, undecided (r1103); m = 1 boundary arbitrary; `assigned_outside_catalog` scoring; real-read truth needs the individual's own copies (HG002, not built); catalog hygiene (gorilla rRNA SM5/SM7/SM577, chrY satellite); "matches no catalog copy" is O3's regime. (full: archive/2026-11/copy_assignment_definition.md:5)

## INDEX

Sections, in order: §0 objective; NOTICE; 0 survival criteria; 1 Definition; 2 Input objects; 3 Feature space; 4 Certificate; 5 Constants; 6 Outputs; 7 Verification; 8 Reconciliation; 9 Known defects; 10 Orphaned modules; 11 Real data; 12 Summary; Reassignment (structural-anchor truth); Appendix A prereg. (full: archive/2026-11/copy_assignment_definition.md:58)

## §0 ⭐⭐ THE OBJECTIVE, RESTATED (2026-08-15) — read this before §1

> **O2 decides, for a read at a multi-copy locus, whether the evidence warrants assigning it to a copy at all, and abstains when it does not.** No 1/k.

- MAPQ-0 names the wrong population: inside 915 multi-copy loci MAPQ 0 is 0.0004 of primaries (genome-wide 0.0020). Ambiguity is real but hidden: 33.70% carry >= 1 secondary, AS margin <= 5% for 21.75% (superseded as the shipped population, see State (b)).
- Reassignment is not the contribution: a secondary fits better by `de` in 1.96% of reads-with-secondary, 12.16% of tight ties (<= 1%), net headroom ~0.1% of reads. **NEVER claim O2 assigns reads better than minimap2.**
- Defended by excision truth: divergence z-score AUC **0.7995** vs MAPQ **0.4944**; operating point TPR 0.5066 / FPR 0.0280; TPR collapses to 0.2404 where foreign reads > 50% of the pile (29 of 83 loci). minimap2 is near-optimal at WHICH copy, at chance about WHETHER; O2 supplies abstention.
- EM changed 0 of 3,081 evidenced decisions (objective decomposes). No non-circular truth for reassignment. (full: archive/2026-11/copy_assignment_definition.md:84)

## ⚠⚠ NOTICE — READ BEFORE QUOTING O2 ANYWHERE (2026-08-10)

`bench/DEFINITIONS_FORMAL.md` called O2 "max-weight facility location"; no shipped code solves it (0 grep hits; only `bench/mwca_approximation_check.py` on synthetic instances). The shipped rule equals facility location only in the degenerate all-open case. Do not quote "facility location" for any measured result. (full: archive/2026-11/copy_assignment_definition.md:162)

## 0. What the definition must survive

Three arbitrariness tests: rule-dependence (two gates in the binary, two in unwired modules), threshold-dependence (alpha has a measured plateau), circularity (the §9.1 headline metric fails it and is not quoted). (full: archive/2026-11/copy_assignment_definition.md:185)

## 1. Definition

Family `F` = ordered copies `c_1 ... c_n`; reads `R` overlapping them; each copy a feature vector, each read an observation vector (§3).

> **Objective.** `L(j) = Σ_{r ∈ R} log P(obs(r) | copy_{j(r)})`, emission = PSV term + junction term. No term couples two reads or two copies (no opening cost, no cardinality, no abundance prior in the hard path), so `L` decomposes.
>
> **Decision rule (per read, 08-10 binary).** `best(r) = argmax_i logl_i(r)`, ties to the lowest copy index (superseded 09-09 by `best_by_psv` + `best_by_duel`; argmax is the `--best-by-alignment` escape).
>
> **Abstention rule.** `thr = α / max(n−1, 1)` (Bonferroni over competitors).
>
> | outcome | condition |
> |---|---|
> | **`tied`** | `min_p(r) ≥ thr`: unresolvable in principle (K = 0 floor) |
> | **`assigned`** | `min_p(r) < thr` **and** `p(r) < thr` **and** `margin(r) > 0` |
> | **`ambiguous`** | `min_p(r) < thr` but the read's evidence does not clear it (or argmax not strict) |
>
> `margin = logl_best − logl_secondbest` (`+∞` for n = 1); `p(r)` = worst pairwise significance, `min_p(r)` = worst attainable; both exact Poisson-binomial upper tails (§4). `tied` is about the substrate, `ambiguous` about this molecule. **Never splits weight 1/k.**

**Emission ≠ decision:** `assigned_copy = best_copy` is printed on every row; filter on `status`.

| symbol | shipped | measured role |
|---|---|---|
| **α** | **1e-3** (CLI `--alpha`) | OPERATIVE, only constant in the gate. Plateau [1e-3, 1e-1]: assigned moves 0.384%; −1.81% at 1e-4, −5.30% at 1e-5 |
| `margin > 0` | strict MLE | sole decider on 7,709 / 1,648,460 rows = 0.468% |
| τ (`margin`) | 2.0 lib, 6.9 CLI | INERT (`use_margin_gate` false) |
| `error_rate` | 0.003 | only where the read has no per-base QV |
| `junction_weight` | 5.0 | ±5.0 nats ≈ 7 PSV matches; fixed, unswept |
| `boundary_tol` | 4 bp | OPEN, never swept |
| `junction_err` | 1e-4 | certificate ε for a junction; OPEN |
| `edit_rate` / `rna_editing_filter` | 0.2 / true | A↔G editing filter, certificate only |
| `PSV_MIN_ALLELE_READS` / `PSV_MIN_JUDGE_COV` | 2 / 4 | read-support PSV validation; < 4 reads kept unjudged |
| `iterative_prune` | false | refused under `--families` |
| `GATE_MIN_READS` / `QUANT_ERROR` | 3 / 0.01 | invariance certificate / soft EM only |

(full: archive/2026-11/copy_assignment_definition.md:199)

## 2. Input objects

(a) Reads: RECORDS in, MOLECULES out. 08-10 pool = every mapped record incl. secondary/supplementary (minimap2 `-Y`, full SEQ), observed in the frame of `best_overlap_copy`, MAPQ never gates (superseded: AS-tied only). One result per (molecule, family): >= 2 records `assigned` to distinct copies ⇒ `ambiguous`; else representative = max `n_decisive`, then margin, min p, lowest index. Arbitrating by margin would equal the primary flag on 323/323 SHARP and 1,629/2,912 = 55.9% panel-wide.
(b) Copy set historically re-derived inside `copy_assign` (GSTM: O1 4 copies vs 0 families).
(c) `--families` / `--copies-fa` consumes the O1 catalog and switches off every family-construction leg; `family_join.tsv` joins by catalog `tid`; RABL2 (5 contigs) refused. Bonferroni denominator = catalog `n`. (full: archive/2026-11/copy_assignment_definition.md:271)

## 3. The feature space

- **PSV columns** from copy-vs-copy alignment, validated by read pileup (dropped only if >= `PSV_MIN_JUDGE_COV` coverage and < 2 alleles reach `PSV_MIN_ALLELE_READS`).
- Since 2026-10-12 each copy is oriented against copy 0 (`is_reverse_homolog`, 15-mer sign test, alpha 1e-3); before, 10/56 KB3781 families got PSVs at most bases (723 vs 9.8 columns/kb); assignments unchanged.
- Engine default (2026-10-12, user decision) `semiglobal_msa_pair` (+2 / −3, gaps min(5 + 2k, 30), significant blocks Karlin-Altschul alpha 1e-3): held-out precision 0.993-0.998 vs 0.42, 5-20x faster; `RUSTLE_PSV_ENGINE=global` restores the old DP (capped 2^29 cells). Assignments and quant byte-identical either way.
- Per spanned decisive column each copy gets `log(1−e)` or `log(e/3)`, `e` = base QV or `error_rate`.
- **Junctions** mapped into the copy's spliced frame (`gen2off`); ±`junction_weight` within `boundary_tol`.
- `n_decisive` = decisive PSVs + decisive junctions; `combined` (PSV + junction) is the decision, `psv` reported for `junction_only`. (full: archive/2026-11/copy_assignment_definition.md:326)

## 4. The certificate

For (`best`, competitor `c`): distinguishing positions the read observes, each with `εⱼ` (`e/3` PSV, `junction_err` junction, `edit_rate` if editing-flagged); `k` support `best`.
- `p_{bc} = P(Σ Bernoulli(εⱼ) ≥ k)`, exact Poisson-binomial O(n²) DP.
- `attain_{bc} = Π εⱼ` (identifiability bound).
- `p(r) = max_c p_{bc}`, `min_p(r) = max_c attain_{bc}`; `min_p ≥ thr` ⟹ `tied`. `min_p` is relative to the same `best`. (full: archive/2026-11/copy_assignment_definition.md:368)

## 5. Every constant, and where it comes from

Table of 25 constants with 08-10 file:line (now dead, see State (a)): values as in §1 plus `PSV_BAND_MARGIN` 1024, `MOSAIC_EPS` 0.01, `MAX_MOSAIC_SITES` 250, `MosaicParams` (11 fields), `RefineParams` min_identity 0.80 / min_coverage 0.50 / sensitive_identity 0.60 (not exposed; moot under `--families`), `--min-copies` 2, `--win` 5,000,000, `--read-cap` 6000 (documented NO-OP). Env vars reaching O2 include `RUSTLE_PSV_*`, `RUSTLE_POSTERIOR_PRIOR`, `RUSTLE_VG_MOSAIC_*`; `RUSTLE_NO_RECOMBINANT_ABSTAIN` has no effect. (full: archive/2026-11/copy_assignment_definition.md:389)

## 6. Outputs — what every column means

- `assignments.tsv` (one row per read x family; `read_name` not unique across families): `family_id`, `assigned_copy` (= best, unconditional), `status`, `n_decisive`, `margin` (3 decimals), `p_value`, `min_p_value` (compare to `α/(n−1)`), `as_*` (reported, never decisive).
- `families.tsv`: `psv_cols`, `resolvable_psv` / `resolvable_j`, `junction_only`, `assigned_j`, `uniq` / `uniq_agree` (see §9.1).
- `quant.tsv`: `abundance` (soft EM, uniform prior at zero PSVs, never feeds hard assignment), `n_reads_hard`, `anchored_reads`, `tie_invariant` (`anchored_reads ≥ GATE_MIN_READS`), `junction_invariant` (>= 3 copy-specific junction reads); FALSE on both = copy leans on the primary flag.
- `famcn_readonly.tsv` (`chi_H`, `depth_cn`, `regime`); `family_join.tsv` under `--families` (aborts if a catalog copy is missing).
- Opt-in: `--posterior`, `--em`, `mosaic.tsv` (report-only), `--dump-psv`, `--phase` (a relabeling, §9.4). (full: archive/2026-11/copy_assignment_definition.md:426)

## 7. Verification — does the binary do what §1 says?

Rule re-implemented from the text (`verify_o2_rule.py`) on 171 runs, 1,648,460 rows: 9,139 (0.55%) undecidable (`margin == "0.000"`), **1,639,305 = 99.9990%** of decidable reproduced; 16 residual are printing ties. Alpha sweep: 1e-1 1,052,310; 1e-3 1,048,266 (−0.384%); 1e-5 −5.30%; 1e-50 −31.06%; uneven below 1e-6 (DAZ −43% by 1e-12, MAGEA −0.7%). Predates the 09-09 defaults, not re-run. (full: archive/2026-11/copy_assignment_definition.md:501)

## 8. Reconciliation with `bench/DEFINITIONS_FORMAL.md`

The real object is the per-read decomposed argmax with a Poisson-binomial abstention certificate, not facility location: the copy set S is fixed by O1, and facility location has no abstention (`tied`/`ambiguous` = 36.3% of reads). The three sites in `bench/DEFINITIONS_FORMAL.md` (lines 141, 225, 343) were corrected; `bench/THEORY.md` §5c left as theory. "assign-or-abstain, no 1/k" survives. (full: archive/2026-11/copy_assignment_definition.md:542)

## 9. Known defects (stated, not hidden)

- 9.1 `uniq_agree / uniq` is anti-correlated with correctness (denominator conditioned on the prediction; all-wrong primary flags move it 100.0% → 0.0% at unchanged accuracy; covers 13.49% of reads): do not quote. "63.7% assigned" WITHDRAWN (record-level); per molecule on the gorilla panel 2,097 / 7,121 = 29.4%.
- 9.2 `margin` printed at 3 decimals: 9,139 rows' status unrecoverable; fix not done.
- 9.3 Recombinant (mosaic) reads detected but never abstained; `recombinant_abstain.rs` unreachable.
- 9.4 `--phase` is a relabeling of `assignments.tsv`: symmetric difference 0 over 119,524 rows (10 run pairs).
- 9.5 Region-scoped (RABL2 refused); O1 definition parameters not exposed; X.4 (2026-08-11) shared the edge rule (`rule.tsv` diff EMPTY on 5 refine calls), clustering still not shared; `--read-cap` no-op. (full: archive/2026-11/copy_assignment_definition.md:574)

## 10. Orphaned modules — verdicts, with evidence

| module | verdict |
|---|---|
| `phasing.rs` (742 lines, 12 tests, 0 callers) | DELETED |
| `o2_materialize.rs`, `o2_columns.rs`, `o2_margin_gate.rs` | DELETE, not executed |
| `recombinant_abstain.rs` | KEEP AND CONNECT (not done; State (b): no module declared in `src/` now) |

25.2% of "O2 tests" guarded unreachable code. Tests 767 / 0 / 11 after deletion (−12); control panel 75 files no diff. (full: archive/2026-11/copy_assignment_definition.md:663)

## 11. O2 on REAL data — commitment rate + abstention certificate (2026-08-11, O2.3)

**Nothing here is an accuracy.** Gorilla OR6737 control panel, 5 multi-copy panels, N = 5,378 primaries (superseded denominator, State (b)).
- 11.1 assigned 0.3466 [0.3340, 0.3594] all; MAPQ 0 (n = 497, 496 MAGEA) 0.0302; MAPQ > 0 0.3788.
- 11.2 2,097 assigned rows, 0 certificate violations, worst `p/thr` 0.2666. `phred_err` clamps ε to [1e-4, 0.25] (Q40 floor), so single-feature quanta are 3.333e-5 (PSV) and 1e-4 (junction): a lone junction certifies only for n ≤ 10, a lone PSV for n ≤ 30. Independent recomputation from BAM QVs: 160/160 clear thr.
- 11.3 tied 45 / 4,885 = 0.0092 of evaluated; 13 / 488 = 0.0266 in MAPQ 0.
- 11.4 99.96% (4,615 rows) of abstentions were the molecule-level contradiction rule; dropping it would raise commitment 0.3466 → 0.8998 (abstention cost, not error).
- 11.5 493/493 unevaluated reads overlap no catalog copy. 11.6 one-sentence summary with guard rail. (full: archive/2026-11/copy_assignment_definition.md:736)

## 12. One-paragraph summary for a reader in a hurry

O2 assigns a long read to one paralog copy of a family fixed by O1, using PSV columns and copy-specific junctions; per-read argmax is the global optimum; reported only if an exact Poisson-binomial certificate clears `α/(n−1)` (α = 1e-3) with a strict argmax, else `tied` (K = 0 floor) or `ambiguous`; never 1/k; unit = molecule; rule reproduced on 99.9990% of decidable rows (08-10 binary). (full: archive/2026-11/copy_assignment_definition.md:909)

## Reassignment on real reads with structural-anchor truth (2026-08-15)

**Status: RUN, ATTACKED, UNDERPOWERED, route refuted** for the contested population. KB3781 matched individual, structural anchors, N = 53 reads (2 families, n_eff ≈ 3); prereg §8 abandons the head-to-head at n < 150. C1 constant baseline 53/53 = primary flag 53/53 (by label imbalance); O2 channel 4/53 = 0.0755, abstention 49/53. L1 removes 29/29 anchors in the contested families: anchoring and near-ties are mutually exclusive. Amendments: reassignment validated only in simulation; near-tie rate 23.83% in anchored two-copy families vs 0.17%; prohibition "O2 assigns better than minimap2" stands. Truth file sha256 `87d1119f06ef9ccac3ff7f750ceee4d558eb9e336eaa2c2d236723e35f556e80`. (full: archive/2026-11/copy_assignment_definition.md:926)

## Appendix A — the PRE-REGISTRATION, written before the run

Design + feasibility of the structural-anchor reassignment truth (2026-08-15, frozen prereg sha256 `0154e72a077dfa002eb72c418e5d566b0766f88c5ea76551e8a99c748de79791`): structural anchors not PSVs, feasibility on the 162 two-copy excision panel, protocol, scope, traps, artifacts. (full: archive/2026-11/copy_assignment_definition.md:1389)
