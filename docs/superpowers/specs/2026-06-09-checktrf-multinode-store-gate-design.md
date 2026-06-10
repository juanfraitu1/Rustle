# StringTie-faithful checktrf multi-node independent-store gate — design

**Date:** 2026-06-09
**Status:** IMPLEMENTED then FALSIFIED by validation → shipped opt-in (default-off). See
`docs/STRINGTIE_PARITY_FINDINGS.md` §6o. The default-on premise below was disproven: the gate
over-rejects 156 StringTie-shared real isoforms (ST-only 104→259) because rustle's checktrf rescue
is load-bearing for recall. Kept as the downstream complement to flow-enumeration parity.
**Builds on:** `docs/STRINGTIE_PARITY_FINDINGS.md` (§6m/§6n), `project_post_flow_gate_pin` memory,
the parse_trflong over-enumeration investigation (workflow wf_abb54bf1-34e).

## Problem & framing

Default rustle emits 187 multi-intron intron chains that StringTie does not (the "rustle-only"
precision gap on the chr19 GGO_19 parity dataset). A multi-agent investigation of the
`parse_trflong` over-enumeration split these 187 into two structurally distinct pools, verified
deterministically against rustle's `path_extracted` `source` field and the real RefSeq annotation
(`../GGO_genomic.gff`, NC_073243.2):

| Pool | Chains | TP / FP (vs real RefSeq, exact intron-chain match) | Nature |
|------|--------|------|--------|
| **flow-only** | 137 | 30 TP / 107 FP | path-selection in the flow loop — architectural, recall-risky |
| **checktrf-only** | 48 | **2 TP / 46 FP** | a clean structural gate divergence — ST-ground-truthed (gate-removable) |
| **checktrf+flow** | 2 | 0 TP / 2 FP | produced by BOTH sources — survive the gate via flow (not removable) |

The gate removes the **48 checktrf-only** chains. The 2 checktrf+flow chains are also produced by
a flow seed, so dropping their checktrf store leaves them emitted via flow — they are unaffected.

The four candidate flow-depletion/ordering/store-gate hypotheses were each **adversarially
refuted** (rustle's coverage ladders descend → depletion works; ST itself extracts near-duplicate
backbone paths). A prior downstream coverage fix (bpcov in the retained-intron gate) over-killed
119 true positives and was reverted (§6n). The flow pool is therefore **out of scope** for this
design and will be scoped separately.

The **checktrf pool** is a one-line, ST-faithful, F1-positive fix and is the subject of this spec.

## The divergence (verified 4 ways)

StringTie's `parse_trflong` checktrf handling (`tools/stringtie/rlink.cpp:10369`):

```cpp
if(!transfrag[t]->shortread && transfrag[t]->nodes.Count()>1) {
    // multi-node long-read: best_trf_match REDISTRIBUTION only
    ... if(tmatch.Count()) { donate coverage to matched preds; transfrag[t]->abundance=0; }
    // NO match → falls out of the block WITHOUT storing. Never becomes a prediction.
}
else {  // rlink.cpp:10413 — reachable ONLY for shortread OR single-node
    if(!eonly || transfrag[t]->guide) { ... store as independent prediction ... }
}
```

So ST's **independent-store branch is structurally unreachable for multi-node long-read
transfrags.** A multi-node long-read transfrag with no kept-path match is silently dropped.

Rustle (`src/rustle/path_extract.rs`): the multi-node long-read match/redistribute block exists at
9542-9765 (mirroring ST's `if(!shortread && nodes>1)`), but on the **no-match fall-through (9766)**
control drops into the independent-rescue block at **9779**, whose own comment admits the deviation:

```rust
// Independent rescue: store complete transfrag as its own prediction.
// Applies to: shortread/<=1-node, AND multi-node longread with no kept-path match.
```

That `AND multi-node longread with no kept-path match` clause is exactly what ST does not do.
These stores produce the 48 rustle-only chains (96% FP).

Verification performed:
1. **Source split** (rustle `path_extracted.source`): 137 flow-only + 48 checktrf-only +
   2 checktrf+flow of 187.
2. **TP/FP** (exact intron-chain vs `GGO_genomic.gff` NC_073243.2): the 48 checktrf-only =
   2 TP / 46 FP.
3. **Rustle mechanism**: the code comment at path_extract.rs:9780-9781.
4. **ST ground truth**: rlink.cpp:10369 if/else structure (read directly).

## Goal

Make rustle's checktrf control flow structurally match ST's: a multi-node long-read transfrag with
no kept-path match terminates (drop) instead of storing as an independent prediction. Ship
**default-on with an opt-out env**. Expected: rustle-only 187 → 139 (48 checktrf-only chains
removed), ST-only 104 unchanged, −46 FP / −2 TP vs real annotation (clean precision win).

## Architecture & components

A single `crate`-internal structural gate at the entry of the independent-rescue block in
`path_extract.rs`, plus opt-out env and trace instrumentation. No flow-loop changes.

### Component 1 — The gate (`path_extract.rs:9779`)

At the top of the independent-rescue block (after the no-match fall-through comment at 9778),
before the existing eonly skip at 9786, add:

```rust
// ST-faithful (rlink.cpp:10369/10413): a multi-node long-read transfrag with no kept-path
// match is redistribute-only in ST and NEVER stored as an independent prediction. Rustle
// historically fell through and stored it (the 48 checktrf rustle-only FP chains). Drop it
// here to match ST. Preserves: shortread/single-node rescue, guides, csr_triggered folds.
// Opt-out: RUSTLE_CHECKTRF_MULTINODE_RESCUE=1 restores the old store-it behavior.
if !is_shortread_tf
    && tf_nodes.len() > 1
    && !csr_triggered
    && !transfrags[t].guide
    && std::env::var_os("RUSTLE_CHECKTRF_MULTINODE_RESCUE").is_none()
{
    record_outcome!(t, SeedOutcome::ChecktrfMultinodeNoMatchDrop);
    emit_checktrf_result!(t, "multinode_no_match_drop", &tf_nodes);
    continue;
}
```

All guard inputs already exist in scope: `is_shortread_tf` (9540 = `!transfrags[t].longread`),
`tf_nodes` (9483, the transfrag's inner nodes), `csr_triggered` (9541), `transfrags[t].guide`.

### Component 2 — Abundance handling (ST-faithful: leave unchanged)

The gate does **not** set `transfrags[t].abundance = 0`. ST only zeroes abundance inside the
*matched* redistribution branch (rlink.cpp:~10410); on no-match multi-node it leaves abundance
untouched. Leaving it unchanged is the faithful behavior and avoids starving any later
rustle post-pass / nascent redistribution that consumes leftover abundance.

### Component 3 — Outcome + trace instrumentation

- Add a `SeedOutcome::ChecktrfMultinodeNoMatchDrop` variant (mirrors the existing
  `ChecktrfEonlySkip` / `ChecktrfRedistributed` variants) so the drop is recorded in the
  outcome accounting and visible in seed-decision dumps.
- Emit `checktrf_result` parity event `"multinode_no_match_drop"` (the `emit_checktrf_result!`
  macro is already used in this block) so the drop is traceable in the parity JSONL.

### Component 4 — Opt-out env

`RUSTLE_CHECKTRF_MULTINODE_RESCUE=1` restores the pre-fix behavior (store). Default-unset = new
ST-faithful behavior. Lets the before/after be toggled in one binary for the gate measurement and
provides a fast revert path if genome-wide validation surfaces a problem.

## Data flow

```
transfrag → checktrf deferred list → multi-node long-read match attempt (9542-9765)
  ├─ match found      → redistribute coverage, abundance=0, continue (UNCHANGED)
  └─ no match (9766)  → NEW GATE (9779):
       ├─ multi-node longread, non-guide, non-csr, env unset → DROP (continue)  ← the fix
       └─ otherwise (shortread / single-node / guide / csr / opt-out) → independent rescue store (UNCHANGED)
```

## Validation discipline

1. **Parity gate** (primary): `python3 bench/gtf_chain_diff.py /tmp/ruP.gtf /tmp/stP.gtf`
   → expect rustle-only 187 → 139, **ST-only 104 unchanged** (over-rejection alarm).
2. **TP/FP vs real annotation**: `/tmp/classify_tpfp.py` (exact intron-chain vs
   `../GGO_genomic.gff` NC_073243.2) → 48 checktrf-only chains removed; net −46 FP / −2 TP
   (the 2 checktrf+flow chains survive via flow and are unaffected).
3. **Default headline guard**: chr19 de-novo vs StringTie — precision up, sensitivity ~flat
   (−2 TP is minor). Confirm no unexpected movement (= leak).
4. **Opt-out byte-identical**: `RUSTLE_CHECKTRF_MULTINODE_RESCUE=1` reproduces current output
   exactly.
5. **2-TP safety check**: identify the 2 real-annotation chains lost; confirm they are NOT
   shipped strand-bundling (`RUSTLE_STRAND_PURE_MINORITY`) or read-chain recovery wins — if any
   overlap, the guard could double-suppress a shipped win (back off / exempt if so).
6. **Genome-wide** (gating, since default-on): run the gw harness (`gw_run.sh`, serial per
   OOM protocol) to confirm net-positive FSM across chromosomes before final sign-off.
7. **Unit test**: a focused Rust test pinning the gate predicate (multi-node longread non-guide
   non-csr → dropped; shortread / single-node / guide / csr → kept), following the existing
   predcluster/test patterns.

## Risk & abort

- **−2 real-annotation TP**: acceptable at 23:1 (46 FP : 2 TP) removed; the 2-TP safety check (step 5)
  guards against eroding a shipped recall win.
- **Post-pass interaction**: leaving abundance unchanged is ST-faithful, but verify no later
  rustle pass depends on the dropped transfrag becoming a stored pred.
- **Genome-wide regression**: default-on requires step 6 to pass. Concrete abort: if gw FSM is
  net-negative or ST-only rises materially, flip the default (make it opt-in) behind the same env
  and re-scope.

## Out of scope

- The **flow-139 pool** (107 FP / 30 TP): path-selection in the flow loop, architectural and
  recall-risky; the four refuted hypotheses establish it is not a simple gate. Scoped separately
  next (user decision: "checktrf now + scope flow next").
- Any flow-loop, depletion, ordering, or coverage-metric change.
- Downstream predcluster / retained-intron filter changes (the bpcov substitution already
  falsified, §6n).
