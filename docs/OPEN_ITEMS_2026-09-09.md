# Open items and loose ends — end of 2026-09-09

State at close: branch `dna-from-genome` @ `099766d`, **85 commits ahead of origin (unpushed)**, suite 871 / 0 / 11,
17 pre-registrations today (`docs/PREREG_*_2026-09-09.md`), register at row 804. Shipped defaults decided today:
AS-tied gate · `best_by_psv` · `origin_drop_indels` ON · `best_by_duel` ON · `gtf_copy_set` ON. Full pre-09-09
escape: `--no-as-tied-only --best-by-alignment --no-origin-drop-indels --no-best-by-duel --no-gtf-copy-set`.
Headline (human MCL0): 1,118 contested = 262 / 512 / 344; excision 247/262 = 94.3 %; GTF 488 placed / 14 dropped /
6 lifted / 40 undecided-with-set; per-tool phantoms 2 / 4 / 15 / 45.

Legend: **P1** = before the next advisor conversation · **P2** = next working week · **P3** = when convenient.
Effort in hours of work, compute excluded unless noted. "Evidence" = where the fact is recorded.

---

## A. O2 engine — the assignment itself

| # | item | evidence | why it matters | suggested action | prio / effort |
|---|---|---|---|---|---|
| A1 | **92 molecules lose their row when an unrelated copy is excised** (36 under copy 6, 26 under 24, 22 under 23, …; primaries at copies 2/22/12) | §6hg, PREREG 4c832450 outcome | Not a wrong call (absent = abstained) but an unexplained mechanism inside the gate/row emission; the excision sweep cannot be quoted as a mechanism table until it is traced | Trace one molecule (`5459_8000`, copy 2, n_cand 19) through the gate with `RUSTLE_STAR_DEBUG` under `excise_all/no6`; expect the tie-outside pass or the region binding to drop it; document or fix | P1 / 2 h |
| A2 | **The origin certificate is tested against `bk` only** (row 799): 27 human reads consistent with copy 12 but winning duels at copy 9 are labelled `origin_rejected` instead of `ambiguous` | §6hj, PREREG 819c1615 | Both abstain, but "origin-rejected" claims the read is foreign to the family; the contested denominator moved by 28 because of it | Pre-register: test the certificate against every candidate that passes; report `origin_consistent_set`; a read consistent with ≥ 1 candidate but not with `bk` becomes `ambiguous`, not rejected | P2 / 4 h |
| A3 | **Sibling blind spot** (row 796): the certificate cannot reject a copy ≲ 0.6 % diverged over a read; 5 of 262 assignments move confidently on excision; **19/262 rest on ≤ 2 columns** inside sibling groups (2/8 0.996, 10/11 0.992, 12/19 0.988) | §6hg | The only precision failure class; the number the advisor will probe | (a) add a per-row `sibling_identity` (max identity of the assigned copy to any other candidate over the read footprint) and a `columns_vs_sibling` count; (b) decide whether ≤ 2-column assignments between ≥ 0.985 siblings are reported as `assigned` or as a labelled weaker class. User's call | P1 / 4 h |
| A4 | **Tie width** (`--as-tie-ratio` < 1.0): "report both, pick later" — only the 0.98 admission count exists (42 % of the 2,278 blind-spot molecules) | §6hc addendum, PREREG 268fe69b | The advisor's "same AS" is exact; a near-tie band is the natural next widening and unmeasured | Pre-register 0.99 / 0.98 on human + MCL1; excision on the newly admitted; present both | P2 / 3 h + 1 h compute |
| A5 | `tie_invariant` is vacuous under the gate (row 786) | §6gy | Dead column / dead test | Remove the column and its test, or redefine as "every tie partner is a catalog copy" | P3 / 1 h |
| A6 | `assigned_outside_catalog`: the outside tie partner (EIF3C class, 8,944 human molecules) is never scored or named | §6gz, row 790 | 27 of 40 undecided isoforms say `copies "A,outside"` — "outside" is unnamed; naming the partner locus makes the set actionable | Record the outside placement's locus (chrom:start-end, and its RefSeq gene if any) in the gate registry; emit `copies "A,outside:chr16:28.7M(EIF3CL)"` | P2 / 3 h |
| A7 | `--admit-aligner-disagreement` and `--indel-psv` remain as OFF flags with measured refutations; `--indel-psv` has a known artifact class (terminal gaps, row 794) it does not exclude | §6hd, §6he | Code surface that will never ship; the register carries the result | Either delete both flags (keep the PREREGs) or leave with a doc comment pointing at the rows. User's call | P3 / 1 h |
| A8 | The contested denominator moves with the certificate (759 → 1,143 → 1,118 today) | §6hf, §6hj | A moving denominator invites the "conditioned on the prediction" trap | Standing rule (already in the docs): always quote the in-catalog AS-tied pool (1,643 human, 122 MCL1) beside the contested count | done as a rule; enforce in any new table |
| A9 | Under `best_by_duel` 4 non-contested MCL1 rows changed `bk`/margin (no status change) | §6hj P5 | Cosmetic; a byte-identity claim must exclude these columns | Note in the escape spec; nothing to do | P3 |

## B. The GTF — the deliverable

| # | item | evidence | why it matters | suggested action | prio / effort |
|---|---|---|---|---|---|
| B1 | **3 certified reads' isoforms fail to lift** to their copy (no transcript there, lift fails) — silently not emitted | §6hn/§6hp (`3 lifts failed`) | An assigned read with no transcript at its copy contradicts "the GTF O2 believes" | Fallback: place by the read's own read-star alignment to that copy (the columns are already computed) or emit at the source copy with `lift_failed "c"`; count them in the stderr line | P1 / 3 h |
| B2 | **Singleton chains never become transcripts** (`min_reads` 3): the 42 copy-22 O2-assigned reads — the SV-carrying isoform, realized inconsistently by the aligner — have no transcript in ANY tool's GTF; 220 of 290 contested molecules we do not carry are singletons | §6hh, row 798 | The hard-locus claim is undercut where it matters most (copy 22 is the SV copy) | Pre-register an "evidence-backed singleton" rule: a certificate-assigned read's chain is emitted even at support 1, tagged `support "1"`; measure the phantom table and the bakeoff again | P1 / 4 h |
| B3 | Old copy attributes (`copy_status`, `copy_votes`, `adjudicated_reads`, `matched_reads`) were computed over POST-gate reads only and still sit beside the new ones (`copies`, `placed_by`, `evidence_*`) | §6hp trap | Two vocabularies on one line; `copy_status "unadjudicated"` means "no AS-tied read", not "no read" | Rename the old ones `tied_*` or drop them under the default; keep them under `--no-gtf-copy-set` | P2 / 2 h |
| B4 | Two residual "phantoms" by the lift script's rule are LIFTED placements (their evidence is a certified read at another locus) | §6hn P1 | The measuring script and the emitter disagree on what evidence is | Teach `isoform_copy_lift.py` to count `placed_by "assigned_read"` as evidence; re-run the table (expect 0) | P2 / 1 h |
| B5 | One group differs between binary and post-processor ({15,16,17} vs {15},{16,17}) — fragment-choice sensitivity of the lift (`min_by_key(d)` over minimap2 fragments) | §6hp | Grouping is not fully deterministic across implementations | Use the best-scoring fragment per pair consistently (PAF column 10 matches) instead of nearest-distance; unit-test on that pair | P3 / 2 h |
| B6 | The lift-aware "carried" metric rewards duplicated placements (row 803); no single headline metric exists for the GTF | §6ho | The advisor will ask for ONE number | Define the composite: per isoform, (carried, address ⊆ evidence set, no phantom); report the three columns together, never "carried" alone; write it into `bench/isoform_bakeoff.py` as `--summary` | P1 / 2 h |
| B7 | Gorilla: partner copies (`member_status partner`) — are they placement targets? MCL1's copy-set run placed 272 with 1 phantom; partners were in `targets` | §6hp | If partners can receive placements, a "partner" address may be a non-member | Check `target_idx` vs `member_status`; decide whether partner rows are addresses or `outside` | P2 / 1 h |
| B8 | The GTF rule was developed on human MCL0 only; gorilla MCL1 ran under it but the held-back family (MCL58) has not | trap 15 | Same rule as for O2: validate where not developed | Run `--gtf` on MCL58, the phantom table (ours vs ours-raw), and the 2-assigned isoforms' placement | P1 / 1 h |

## C. Validation and benchmarks

| # | item | evidence | why it matters | suggested action | prio / effort |
|---|---|---|---|---|---|
| C1 | **Excision against competitors never run** (bakeoff PREREG P4): "only we abstain when the copy is gone" is untested for flair / StringTie / isoseq | O1_O2_OPEN §D | The strongest comparative claim, unmeasured | Reference-level excision: mask copy 2 in CHM13, realign the family reads, run each tool, count transcripts placed at the sibling; one copy is enough | P2 / 6 h + cluster compute |
| C2 | **Gorilla competitor arms are TESTIS** (`GGO_mm.bam`); only ours vs ours-raw is comparable there | §6hb | Cross-species replication of the bakeoff is missing | Re-run isoseq collapse / flair / StringTie on `GCA_029281585.2_flnc_mm.bam` (fibroblast) on the cluster (isoseq collapse: 3 min 38 s last time) | P2 / cluster |
| C3 | The excision sweep's pooled number (94.3 %) counts "row vanished" as abstain (A1) | §6hg | If A1 turns out to be a leak, the number changes | Blocked on A1 | — |
| C4 | Pooled isoform certificate refuted as a placer (71 %) but its per-isoform `pooled_copy`/`pooled_margin` attributes are only in the post-processor, not in the binary | §6hm | The attribute was promised as "report only" | Either add to the binary (cheap: the per-read pairwise stats exist) or drop the promise | P3 / 3 h |
| C5 | Support-matched comparisons (chain multiplicity ≥ 2) are a post-hoc control in two PREREGs | §6hh, §6ho | Honest but not pre-registered | Pre-register the support-matched form once, as the standing secondary table | P3 / 0.5 h |
| C6 | MCL2 (64 copies) is uninformative for O2 (1 AS-tied) and was the nominal hold-out for O1; MCL58 is O2's hold-out | §6hi | O1 still has no informative held-out family for the last rules (cross-family overlap, node floor) | Pick an O1 hold-out by an outcome-blind rule (e.g. most kept units among unused families) | P2 / 2 h |

## D. O1 / O3 carried over (from `O1_O2_OPEN_2026-09-09.md`, unchanged today)

| # | item | evidence | suggested action | prio |
|---|---|---|---|---|
| D1 | **G3**: the constants 0.70 / 0.30 / 300 bp / size 3 / 50 kb are unjustified; the inflation and Jaccard sweeps are the instrument | §6gu | Extend both sweeps to the identity and coverage floors; the code exists | P2 |
| D2 | Identity-weighted density (advisor criticism 5): shown to differ (0.818 vs 0.724 on TBC1D3), not validated, not swept, not priced on a held-out family | §6gw | Sweep + held-out family before it is claimed | P2 |
| D3 | Genome-wide alignment runtime: killed at 8 min for memory contention — the only regime where the Jaccard prefilter argument bites | §6gu | Re-run `minimap2 -x asm20 … seedgw/allgenes_gw.fa` under `/usr/bin/time` alone on the machine | P3 |
| D4 | **G4**: the O2 certificate is conditional on the column set (row 691) | §6eu | Now partly answered by the duel (§6hj) and the excision sweep; write the one-paragraph statement | P2 |
| D5 | **G5**: RNA admission of unannotated loci is not implemented (~4 % of Soto members lost at any-overlap); this is where O3 lives | §6ev | The real scientific gap in O1; a design decision, not a measurement | P2 (user) |
| D6 | O3 flag pass: pending user decisions — read-linked admission of partial paralogues, the chain rule for members whose reads fall outside the chain, the 44 nested unit pairs, the 102 annotated-no-unit loci | §6fm/§6fn | Decisions | P2 (user) |
| D7 | Two O1 definitions coexist in the docs (E_r/γ-quasi-clique vs MCL units) | FRAMING_AUDIT | One sentence in THESIS_OBJECTIVES choosing the shipped one | P1 (user) |

## E. Engineering hygiene

| # | item | suggested action | prio / effort |
|---|---|---|---|
| E1 | **85 commits unpushed** on `dna-from-genome` | `git push` (after `git rev-parse` sanity); `tools/stringtie` stays dirty by convention | P1 / 5 min |
| E2 | Five escape flags reproduce the pre-09-09 output | Add `--legacy-2026-09-08` that sets all five; keep a test that its outputs hash to `91081887` / `ff0b8f16` / `a4d0f5cd` on the fixtures | P2 / 2 h |
| E3 | The lift rule exists twice (Python `isoform_copy_lift.py` + Rust `LiftBlocks`) and the emission rule twice (`gtf_copy_set.py` + the binary) — drift risk | A `tests/` fixture where the binary's GTF must equal the post-processor's (small synthetic family, 3 copies, one phantom, one lift) | P2 / 3 h |
| E4 | `RUSTLE_COPYSET_DEBUG*` env vars | Fold into `--dump-star`-style flags or leave documented | P3 |
| E5 | The GTF emitter's `bam_reads` are post-gate; other attribute code paths may silently assume "every read" (the productivity ORF bar uses transcripts only — fine; `matched_reads` does not) | Grep every use of `read_chain`/`matched` in the emitter; document the invariant at the top of the block | P2 / 1 h |
| E6 | `MEMORY.md` is 206 lines (limit 200) | Compact the index; move detail to topic files | P3 / 0.5 h |
| E7 | Scratch drivers (`run_*.sh`, `count_astied.py`, `indel_dump_analysis.py`) live in the session scratchpad; only `bench/` is in the repo | Move the reusable ones (`count_astied.py` — the hold-out selection rule) into `bench/` | P2 / 0.5 h |
| E8 | `tool_bakeoff.py --restrict` was added but the comparison scripts use the unrestricted calls | Leave; note in the docstring | P3 |

## F. Meeting and documentation

| # | item | suggested action | prio |
|---|---|---|---|
| F1 | The five meeting artifacts (overlap, primaries-only, cores in a browser view, runtime vs Jaccard, TBC1D3 hierarchy, density vs similarity) predate today's numbers | Refresh the O2 numbers in the overlap/primaries artifacts (262/512/344; the copy-set GTF; the per-tool phantom table) | P1 / 2 h |
| F2 | `THESIS_OBJECTIVES.md` still says "NOT an assembler" | User's wording change (Part 0h) | P1 (user) |
| F3 | `docs/MEMORY_DIGEST.md` has not absorbed the 09-09 sections (§6gs–§6hp) | Append the digest lines | P3 / 1 h |
| F4 | The advisor dossier Part 0h now states the deliverable; the excision sweep, the sibling limit, and the pooled-certificate refutation deserve one paragraph each there | Write them | P1 / 1 h |
| F5 | A one-page "what the GTF says" reader guide (attributes: `copies`, `copies_undecided`, `placed_by`, `evidence_*`, `copy_status_final`, `lifted_from`) | Write `docs/GTF_COPY_SET_FORMAT.md` | P1 / 1 h |

## Suggested order

1. **E1 push**, then **A1** (trace the vanishing rows) and **B8** (GTF on the held-out family) — both cheap and both guard numbers already quoted.
2. **B1 + B2** (failed lifts; evidence-backed singletons) — the two places the deliverable still contradicts O2.
3. **A3** (sibling class labelling) and **B6** (the composite GTF metric) — the two questions the advisor will ask.
4. **F1, F4, F5** — refresh the meeting material with the final numbers and the GTF format guide.
5. **A2, A4, A6, B3** — the next pre-registered measurements.
6. **C1, C2** — the cluster jobs (competitor excision; fibroblast arms).
7. D-items as the user decides.
