# Design: StringTie-faithful retained-intron scope widening

Status: DESIGN (approved 2026-05-31). Next: implementation plan (writing-plans).
Grounding: needy-loci decision trace (`project_needy_loci_decision_trace`); over-enum floor
(`project_jfp_missr_characterization`); prior local-fix-vs-headline caveat (`project_future_parity_targets`).

## 1. Problem & insight

The needy-loci decision trace (2026-05-31) found that the over-enumeration FP gap is NOT mostly the
flow-depletion floor — it decomposes into 3 targetable hidden ST decisions. The highest-value/
lowest-risk one is **retained-intron**: StringTie fires `pred_kill reason=retained_intron stage=pairwise`
against a LOCAL high-coverage exonic killer (exon_cov ~190 vs intron_cov ~18), removing FP isoforms
whose intron is really a low-coverage dip inside a real exon.

Rustle ALREADY has the machinery — a pairwise killer/victim filter (`transcript_filter.rs:3063-3265`,
`n1`=higher-coverage killer, `n2`=overlapping victim) gated by `build_lowintron_flags`
(`transcript_filter.rs:1557`) + `retainedintron_like` (`:1666`). But it UNDER-FIRES: in RSTL.398 the
victim `.17` was only paired against a single full-length killer (30721740-30848012) and never against
the local one; in RSTL.184 the check fired 0×. So the gap is a SCOPE issue, one of:
- **(a) lowintron FLAGS:** RU doesn't mark the victim's intron as low-cov (the per-tx bpcov drop in
  `build_lowintron_flags` differs from ST's local-exonic comparison).
- **(b) killer-PAIRING scope:** the victim only pairs with killers in its `build_significant_overlap_adj`
  set that precede it in `ord` (`ord_pos[n2] <= oi` skip at `:3069`), excluding the right local killer.

**Insight:** RU owns the filter; this is a scope fix, not new code. RU also already has
`RUSTLE_LOWINTRON_ORACLE` (`transcript_filter.rs:1645-1659`) that forces ST's lowintron mask, and ST
already emits `pred_intron_low` — so the flags-vs-pairing diagnosis is nearly free. CAVEAT: prior
memory repeatedly shows local fixes don't move the genome-wide headline; the trace found ~3 FPs across
2 loci, so the genome-wide prize may be single-digit and/or offset by TP-cost (over-killing real
alt-isoforms). Hence oracle-first.

## 2. Goal & non-goals

- **Goal:** measure — oracle-first, default-OFF — whether widening RU's retained-intron check to ST's
  pairwise/local scope removes over-enum FPs genome-wide at acceptable Sn cost; diagnose which scope
  (flags vs pairing) under-fires; decide ship/shelve from the oracle.
- **Decision rule (user-chosen):** oracle-bound first, then decide; bound prize AND diagnose root cause.
- **Non-goals:** the donor-snap/transfrag-collapse lever (#2) and checktrf-drop (#3) from the trace
  (separate experiments if this clears); the flow-depletion floor; non-chr19 data.

## 3. Architecture — Phase 0 (two oracles) → abort gate → Phase 1 (targeted fix)

- **Phase 0a — total prize ceiling:** capture ST `pred_kill(reason=retained_intron)` genome-wide
  (chr19); attribute onto RU final chains; net = (RU-only-FP matched) − (RU-TP matched). The ceiling
  for any retained-intron-scope fix; the abort number.
- **Phase 0b — root-cause + flag-fixable prize (existing infra):** capture ST `pred_intron_low`
  genome-wide → build the `RUSTLE_LOWINTRON_ORACLE` map → run RU with it → measure FP-reduction/
  TP-cost/Sn. Diagnosis: 0b ≈ 0a → FLAGS; 0b ≪ 0a → PAIRING. Also an independent second prize estimate.
- **Abort gate:** net ≤ 0 or < ~5 → shelve Phase 1.
- **Phase 1 (gated) — targeted fix** behind `RUSTLE_RI_LOCAL` (default OFF): widen exactly the scope
  0b identified.

## 4. Components

### Phase 0a — `bench/retained_intron_prize.py` (new; + optional tiny ST emit enhancement)
- **Inputs:** ST `pred_kill` events genome-wide (`STRINGTIE_PARITY_LOG`, `STRINGTIE_PARITY_FILTER_CHROM=
  NC_073243.2`, all steps, NO range), filtered to `reason=retained_intron`; RU final GTF; ST final GTF;
  reference `/mnt/c/Users/jfris/Desktop/GGO_19.gtf`.
- **Matching:** each ST retained-intron kill → RU final chain by `(strand, span, nexons)`. If the
  `pred_kill` payload lacks the intron chain and span+nexons is ambiguous, enhance ST's `pred_kill`
  emit (`tools/stringtie/rlink.cpp`) to include the killed chain's introns (~5 lines, mirroring the
  `junction_raw` emit) for exact matching.
- **Output:** prize = matched RU-only-FP chains; TP-cost = matched RU-TP chains; net; the gate verdict.

### Phase 0b — `bench/build_lowintron_oracle.py` (new) + existing `RUSTLE_LOWINTRON_ORACLE`
- **Inputs:** ST `pred_intron_low` events genome-wide.
- **Build map:** convert ST's per-chain low-intron mask into the exact format `RUSTLE_LOWINTRON_ORACLE`
  consumes (chain-string key `"d1-a1,d2-a2,..."` using `exon[j-1].end+1 .. exon[j].start`, value
  `Vec<bool>`; see `transcript_filter.rs:1648-1657`).
- **Run:** `RUSTLE_LOWINTRON_ORACLE=<map> ./target/release/rustle -L GGO_19.bam -o /tmp/ru_riorc.gtf`
  → gffcompare vs reference + `bench/gtf_chain_diff.py` vs baseline.
- **Output:** FP-reduction, TP-cost, Sn; flags-vs-pairing diagnosis (vs 0a).

### Phase 1 — targeted fix (`RUSTLE_RI_LOCAL`, default OFF; only if gate clears)
- **If FLAGS:** widen `build_lowintron_flags` (`transcript_filter.rs:1557`) — flag a victim's intron as
  low when a LOCAL high-coverage prediction is exonic over it, not only via the per-tx bpcov drop.
- **If PAIRING:** widen the killer set at `transcript_filter.rs:3063-3068` so victims also pair against
  local high-coverage overlapping killers currently excluded by the `ord_pos[n2] <= oi` / significant-
  overlap restriction.
- Composes with existing `RUSTLE_RI_FILTER_OFF` / `RUSTLE_END_RI_*` toggles.

## 5. Data flow (fix ON)

flow → predictions → **[pairwise filter: widened retained-intron scope (`RUSTLE_RI_LOCAL`)]** → final GTF.

## 6. Validation (gates)

- **Phase 0a (abort gate):** net = RU-only-FP matched − RU-TP matched to ST retained-intron kills.
  Abort if ≤ 0 or < ~5.
- **Phase 0b:** `RUSTLE_LOWINTRON_ORACLE` FP-reduction/TP-cost/Sn; diagnosis + independent prize
  estimate. If 0a and 0b disagree wildly, investigate before trusting (segmentation-prize lesson).
- **Phase 1 (if reached):** genome-wide j-FP reduction + Sn vs 95.6/90.5. Ship only on a net F1 gain
  (Pr up, Sn not materially hurt); else opt-in + record the cost.
- **Default-unchanged guard** at every step (flags default OFF; baseline 95.6/90.5 verified).

## 7. Safety & exit

- 0a analysis-only; 0b uses existing default-OFF `RUSTLE_LOWINTRON_ORACLE`; any ST `pred_kill` emit
  enhancement is parity-logging-gated; Phase 1 behind `RUSTLE_RI_LOCAL` (default OFF). Default verified
  each step.
- **Dominant risk — TP-cost:** a widened retained-intron filter can kill genuine alt-isoforms with
  legitimately low-coverage introns. 0a measures it directly (RU-TP matched); Phase 1 re-checks Sn
  genome-wide. This is the main abort trigger.
- **Abort at 0a (plausible per the standing caveat):** record the retained-intron prize ceiling + the
  flags-vs-pairing diagnosis; shelve Phase 1. Settles whether the trace's local finding generalizes.
- **Proceed → Phase 1:** widen the identified scope; ship only on a confirmed F1 gain.
- **Honest expectation:** trace found ~3 FPs/2 loci; local fixes historically don't move the headline.
  0a may show a single-digit prize offset by TP-cost. Either outcome is a definitive answer.
