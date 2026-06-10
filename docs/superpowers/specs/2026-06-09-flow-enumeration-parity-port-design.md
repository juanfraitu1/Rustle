# Flow-enumeration parity port (parse_trflong_st convergence) — strategy design

**Date:** 2026-06-09
**Status:** Strategy design (approved approach; multi-session implementation deferred to fresh sessions)
**Builds on:** `docs/STRINGTIE_PARITY_FINDINGS.md` §6n/§6o, `project_flow_parity_scope`,
`project_post_flow_gate_pin`, `project_st_shadow_mode` memories; workflow wf_abb54bf1-34e.

## Problem & framing (read this first — the goal is subtle)

Default rustle emits 187 multi-intron chains StringTie does not. Two prior downstream fixes
(Approach A: bpcov in the retained-intron filter; Approach B: checktrf multinode store-gate) were
each built and **falsified by validation** — both over-kill recall because rustle's downstream
compensates for differences in its **flow**. The convergent conclusion: the only remaining lever is
flow-enumeration parity itself.

**Critical reframe (project_flow_parity_scope, MUST be understood):** at the distinct-extracted-chain
level rustle does NOT over-enumerate — rustle extracts 4029 chains, ST extracts 4023. The 2× seed
surplus (10526 vs 7982) collapses to duplicates. rustle and ST extract nearly the same COUNT but
DIFFERENT SETS: rustle-extra = 496 (21 real isoforms ST's flow MISSED + 475 FP), ST-extra = 490
(2 real). So **rustle's flow finds ~19 more real isoforms than ST's**; the 187 final rustle-only =
**32 TP (recall wins) + 155 FP**, inseparable.

**Therefore this port is explicitly ST-MIMICRY: converging rustle's flow to StringTie's will TRADE
~32 real-isoform recall for precision-vs-ST.** The user has chosen this with eyes open, to maximize
ST-parity precision. It is NOT a pure improvement vs the real annotation. This framing must survive
into the plan and the findings notes so the −32 TP is never mistaken for a regression-to-fix.

## Current state of the port (grounded 2026-06-09)

A substantial flow port already exists and is **default-off**:
- `src/rustle/parse_trflong_st.rs` (~60 KB): `long_max_flow_st` (:991), `parse_trflong_st` (:1077),
  gate `canonical_active()` (:1284 = `RUSTLE_PARSE_TRFLONG_ST_CANONICAL`).
- `src/rustle/path_extract.rs` canonical branches: back_to_source at :7351, fwd_to_sink at :7446,
  the `canonical` selector at :7163; plus `stringtie_exact()` sites (:3609, :4345).
- `stringtie_exact()` (`stringtie_parity.rs:173`) is **default-ON** — rustle's DEFAULT already runs
  many ST-faithful gated sites (junction merge filters, alt-junction demotion, fwd_to_sink past_seed
  gate, onpath_long reach, novel_splice_rescue).

**Measured (chain-diff vs /tmp/stP.gtf):**

| Config | Rustle-only |
|--------|-------------|
| default (`stringtie_exact` ON, canonical OFF) | 187 |
| `RUSTLE_PARSE_TRFLONG_ST_CANONICAL=1` | **223 (WORSE)** |

So the existing canonical port **regresses** — it is a non-converging partial attempt. This effort is
therefore **"debug a partial port that currently makes things worse,"** not "finish a near-complete
one." This matches the documented "shadow layers do not converge incrementally" finding.

## Goal

Under `canonical_active()` only (default byte-identical throughout), drive the **path_extracted chain
divergence** (rustle-extra 496 / ST-extra 490) toward 0 by converging `long_max_flow_st` and the
canonical back/fwd extraction to StringTie's `long_max_flow` (rlink.cpp:9856) + back_to_source /
fwd_to_sink, one first-divergence at a time. As the flow converges, the final-output gate
(`bench/gtf_chain_diff.py`) rustle-only should fall from 223 toward **~0** (ST-parity) — note this
removes ALL 187, BOTH the 155 FP AND the 32 recall-win TP, because an ST-faithful flow produces
neither. The 32-TP loss is the cost of ST-mimicry, incurred en route, not a floor. The true residual
floor is whatever remains from downstream-filtering differences once the flow matches (expected small).

## Method (the one that has actually worked here)

**Deterministic first-divergence tracing**, the method that cracked the original 187 analysis and
that the user's repeated guidance endorses ("StringTie is deterministic and available; divergences
MUST be replicable"). NOT rule-tuning, NOT threshold sweeps (every such attempt hit the wall).

Per iteration:
1. Run both tools with the parity log (`STRINGTIE_PARITY_LOG` / `RUSTLE_PARITY_LOG` +
   `RUSTLE_PARSE_TRFLONG_ST_CANONICAL=1`), capturing `path_extracted` + `parse_trflong_seed`.
2. Pick ONE divergent locus where canonical extracts a different path than ST (start with the loci
   driving the 223-vs-187 regression — chains canonical extracts that BOTH default-rustle and ST do
   not).
3. Trace the FIRST point inside `long_max_flow_st` / canonical back/fwd where the path or flux
   diverges from ST's `long_max_flow` at that locus (per-node flux, capacity, the back/fwd reach).
4. Fix that one divergence in `parse_trflong_st.rs` under the canonical gate. Default path untouched.
5. Re-measure at the chain-diff gate. Apply the convergence fork (below).

## Architecture & components (implementation is multi-session)

### Component 1 (FIRST sub-step) — Characterize why canonical regresses to 223
- Diff the canonical `path_extracted` chain set vs (a) ST's and (b) default-rustle's. Name the chains
  canonical adds (the +36) and removes. Are they ST-extra recoveries (good) or new canonical-only FPs
  (bad)? Output: `docs/STRINGTIE_PARITY_FINDINGS.md §6p` note + the worklist of first-divergence loci.
- **No porting before the regression is named.** (This is the only sub-step to execute first; it
  decides whether the existing port is close-but-broken or fundamentally off — feeding the
  debug-vs-rebuild decision the user deferred.)

### Component 2..N — Per-divergence convergence of long_max_flow_st
- Each is one first-divergence trace → one canonical-gated fix in `parse_trflong_st.rs` → one
  chain-diff measurement. Sized one locus/mechanism at a time; each its own session-scoped task.
- Where a divergence is a pure function (flux/capacity/reach given graph+path), add a Rust unit test
  pinning ST's behavior (pattern: the existing `long_max_flow_st` tests, if any, else new).

### Component final — Flip decision
- Once canonical rustle-only is near 0 (ST-parity) AND ST-extra-extracted ≈ 0, evaluate making
  canonical the default (it trades the 32 TP for precision-vs-ST). The flip is a SEPARATE, explicit
  decision with its own genome-wide validation — NOT automatic.

## Validation discipline (hard-won, non-negotiable)
- **Gate = `bench/gtf_chain_diff.py /tmp/ru_canon.gtf /tmp/stP.gtf`** (final-output rustle-only).
  Secondary: the `path_extracted` chain diff (rustle-extra / ST-extra) — the truer convergence signal.
- **Default-OFF regression guard every change:** default (canonical OFF) chain-diff must stay 187/104
  (±run-noise). Any default movement = the canonical gate leaked → fix immediately.
- **TP accounting via real annotation:** `bench/classify_checktrf_tpfp.py` (or equivalent) vs
  `../GGO_genomic.gff` NC_073243.2 — track the −32 TP explicitly so it is never mistaken for a bug.
- **F1 NOT the per-iteration signal** — the chain-divergence count is. Convergence is the goal.

## Pre-registered convergence fork (the risk management)
- **Converging** (each fix reduces rustle-extra-extracted toward ST, default stable) → continue the
  per-divergence port; it is the right lever.
- **Non-converging** (fixes are individually ST-faithful but the extracted-chain divergence does not
  fall, mirroring the shadow-layer history) → STOP. Decisive evidence the flow divergence is not
  closable by incremental faithful fixes, and the 187 is irreducible. Record and end the parity track.

Both outcomes are informative and worth the first few sessions. Defining the fork up front prevents an
open-ended sink — the same discipline that made Approaches A/B decisive instead of endless.

## Risk & abort
- **−32 TP is expected, not a bug** (the reframe). Abort is NOT triggered by losing those.
- **Non-convergence** is the real risk (three prior ST-faithful attempts did not converge
  incrementally). The pre-registered fork is the abort: if the extracted-chain divergence plateaus
  across several faithful fixes, stop.
- **Default leak** (a canonical change touching the default path) — caught by the default-OFF guard.

## Out of scope
- Default-mode behavior (all changes `canonical_active()`-gated, default byte-identical).
- Downstream gates (Approaches A/B — done, falsified, recorded).
- The genome-wide flip-to-default decision (separate, post-convergence).
- Non-`-L` / guided / VG modes.
