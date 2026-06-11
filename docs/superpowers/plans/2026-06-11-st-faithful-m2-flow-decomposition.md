# ST-Faithful Baseline — Milestone 2: Flow Isoform-Decomposition Co-Port

> **For agentic workers:** This continues `docs/superpowers/specs/2026-06-10-st-faithful-baseline-design.md` (Milestone 1 shipped: witness gate + SE-flux fix; mini3 3/16→1/7, chr19 291→269). Steps use checkbox (`- [ ]`).

**Goal:** Make rustle's long-read flow enumerate StringTie's full per-locus isoform set, converging the remaining 7 ST-only mini3 chains (and the genome-wide tail) that are *extracted-then-RI-killed* or *never-enumerated* because rustle's leaner flow produces fewer alt-terminus / alt-splice variants.

**Architecture:** Discovery-first. Read-level evidence (2026-06-11) proved the remaining ST-only are **assembly-invented sub-isoforms from the same reads** — StringTie's flow decomposes coverage into more isoforms (e.g. an alt-TES terminating at node 36062859-36063554) than rustle's. The existing partial `_st` flow port (`canonical_active()`) **regresses** (chr19 273→284, it *under*-extracts) — so a naive flip is wrong. The exact decomposition divergence must be pinned, then ported behind `precise_mode()`, accepting transient regression, converged via coupled filter stages, validated genome-wide at every step.

**Metric:** rustle-vs-StringTie TOTAL divergence (`gtf_chain_diff` vs `/tmp/st_all.gtf`), chr19 baseline now **269**. mini3 `bench/mini3/check.sh` (1/7). Escape-hatch `check_precise.sh` byte-identical to `4705ab1` — **non-negotiable** every task. Validate EVERY change on chr19 (7s) before keeping.

**Revert point:** `cab62d4` (current HEAD; Milestone-1 banked).

---

## Why this is hard (constraints learned)

- **7 targeted levers already regress** (ST-order, geometry-RI, endpoint-demotion, canonical, redist-scoping, RI-relax, hardend-truncate). The remainder is architectural (flow decomposition), not a gate.
- **The decomposition trades precision:** rustle's leaner flow = fewer FPs. Matching ST's richer enumeration risks rustle-only FPs; the win is only net-positive if the *coupled downstream filters* (RI/isofrac/readthr) also match ST so the extra isoforms that ST kills are killed by rustle too.
- **cov formula is NOT the issue** (shared-chain ratio 0.972, confirmed).

---

## Task 1 (DISCOVERY — gates everything): pin the flow-decomposition divergence

**Files:** read-only — `/tmp/m3_st3.jsonl` (ST), `/tmp/m3_r2.jsonl` (rustle), `tools/stringtie/rlink.cpp` (parse_trflong / long flow), `src/rustle/{max_flow.rs,path_extract.rs}`. Helpers `bench/mini3/{trace_query,list_divergent,nearmiss}.py`.

- [ ] **Step 1: Re-establish the worked example.** Locus B, the alt-TES pair: ST emits a path terminating at node `36062859-36063554` (STRG.2.1-class); rustle emits only the continuing path (`...-36063658`, RI-killed). Confirm via `path_extracted` (rustle) vs ST `path_extracted` that ST has the 36063554-terminus path and rustle does not.

- [ ] **Step 2: Instrument ST's path-generation fork.** In `tools/stringtie/rlink.cpp` parse_trflong / the long-read max-flow loop, add a `pd_emit("flow_path_gen", ...)` per extracted path capturing: seed transfrag, the chosen path's terminal node, the flux, and (critically) whether ST *also* spawned an alternate-terminus path at an internal branch. Rebuild ST (`cd tools/stringtie && make`; output must stay byte-identical — observational only). This reveals the exact decision where ST forks the alt-TES variant.

- [ ] **Step 3: Compare to rustle's path-generation.** Rustle already emits `path_extracted` + `node_flux`. Diff ST's `flow_path_gen` vs rustle's per-seed path set for locus B. Identify the precise rule: *when does ST emit an extra terminal variant that rustle does not?* (Hypotheses to test, not assume: ST re-seeds from internal high-flux nodes; ST's flow leaves residual that re-extracts as a shorter variant; ST splits at junction-acceptor boundaries.)

- [ ] **Step 4: Write the divergence finding** to `docs/STRINGTIE_PARITY_FINDINGS.md` (new §). The exact ST decision + source line + the rustle equivalent that must change. THIS determines Tasks 2-4 (do not pre-specify them — they are discovery-driven, like Milestone 1).

---

## Task 2 (PORT): port the decomposition decision behind `precise_mode()`

- [ ] Wrap the divergent path-generation decision: `if precise_mode() { today (lean) } else { ST-faithful (emit ST's variant set) }`. Use a new opt-out flag for the ST-faithful branch during bring-up (e.g. `RUSTLE_FLOW_DECOMP_ST`).
- [ ] Build; **escape-hatch `check_precise.sh` must stay green**; measure mini3 + chr19. Expect ST-only to DROP and rustle-only to RISE (the extra isoforms ST kills downstream). **Record both.** This regression is accepted *only if* Task 3 reclaims it.

## Task 3 (COUPLED CONVERGE): match ST's downstream filters on the new isoforms

- [ ] The extra isoforms must be killed by rustle's RI/isofrac/readthr exactly as ST kills them. Trace `pred_kill` on the new rustle-only chains vs ST `pred_kill`. Port the exact filter decision (the RI geometry / isofrac threshold) so the ST-killed extras are killed and the ST-kept extras survive.
- [ ] Build; escape-hatch green; mini3 + chr19. **Net divergence must improve below 269**, else revert the whole G-group (Tasks 2+3) — the co-port only ships if it generalizes.

## Task 4 (GATE): milestone validation

- [ ] mini3 toward 0/0 (or characterized residual); chr19 < 269; `check_precise.sh` green; suite 284/0 under `RUSTLE_PRECISE=1`; full-BAM smoke.
- [ ] Commit; update the design-doc status.

---

## Notes for the executor

- **Discovery (Task 1) is the gate.** Do NOT port before the exact divergence is pinned — every guess so far regressed.
- **Tasks 2+3 ship together or not at all** (the design's coupled-group rule). Flow-richer without matched filters = precision regression.
- ⚠ Never `pkill -f rustle`. ⚠ whole-genome OOMs — chr19 only. ST rebuild: `cd tools/stringtie && make` (DEBUG ok, output identical). *.bam/*.gtf gitignored → `git add -f`.
