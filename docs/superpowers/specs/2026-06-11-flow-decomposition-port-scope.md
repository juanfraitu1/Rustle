# Flow-Decomposition Port — Scope / Design

**Status:** SCOPE (design phase — no implementation until approved)
**Date:** 2026-06-11
**Branch:** vg/flow-capacity-apportionment · baseline HEAD b26eb16
**Author context:** distilled from the readthr/isofrac root-cause investigation (memory `project_coupled_port_viability`).

---

## 1. Problem statement (the unified root)

The entire chr19 rustle-vs-StringTie divergence — **both** the 186 Rustle-only and 80 ST-only multi-intron chains — reduces to **one** cause: rustle's max-flow **over-decomposes** each locus into more, thinner transcript paths than StringTie.

Worked example (locus STRG.225, NC_073243.2:32005068-32079598, strand −):

| quantity | rustle | StringTie | note |
|---|---|---|---|
| per-node `nodecov`, `noderate` | 560, 362.6 | 564, 362.6 | **MATCH** — substrate is sound |
| dominant-path flux | **119** | **143** | the seed of the cascade |
| flow paths through node 2 | 64 | 29 | over-decomposition |
| seeds processed | 152 | 73 | weak depletion (see §3) |
| transcripts emitted | 17 | 10 | |
| STRG.225.9 cov | **1.9** | **9.6** | minority starved |
| → fate | killed by faithful isofrac (0.01×usedcov=4.7) | kept | → ST-only divergence |

The downstream filters (`readthr`, `isofrac`, `retained_intron`) are **faithful ports** — verified line-for-line. They correctly kill chains that, *at rustle's under-apportioned coverage*, are genuinely underwater. The fix must therefore be in the flow decomposition, **not** the filters (patching filters = the knob-tuning `feedback_full_stringtie_mimicry` forbids).

### Causal chain (each link empirically established, see §7 evidence)
```
rustle dominant flux 119 (vs ST 143)        ← THE bottleneck; identical nodecov/noderate
  → weak depletion (drains less shared cov)
  → downstream seeds retain flow (152 vs 73; ST self-eliminates via zero_flux reject)
  → 17 thin paths vs 10 fat paths
  → minority isoforms get 1.9 cov vs 9.6
  → die to faithful isofrac/readthr  →  80 ST-only
  → the extra thin paths              →  186 Rustle-only
```

### What it is NOT (all falsified this session)
- **Not the bpcov/lowintron substrate** — `nodecov`/`noderate` match ST node-for-node. Categorically different from the lowintron dead-end.
- **Not the `!exact` guards** — default already runs guards-off (`stringtie_exact()` default-on); forcing them on (`RUSTLE_STRINGTIE_EXACT=0`) is *worse* (flux 110, div 227/123).
- **Not a seed-merge gap** — ST does not merge thin seeds (60% of ST seeds are single-read, == rustle's 57%).
- **Not `trthr`** — both eliminate the same 7 sub-1.0 seeds at the locus.

---

## 2. Goal & success criteria

**Goal:** make rustle's long-read flow decomposition concentrate like StringTie's — fatter dominant path → stronger depletion → fewer surviving seeds → minorities clear the faithful thresholds.

**Primary metric:** chr19 multi-intron divergence vs `/tmp/st_all.gtf` (`bench/gtf_chain_diff.py`), currently **186 RU-only / 80 ST-only**, target → 0.

**Local oracle:** rustle's dominant-path flux at STRG.225 → **143** (matches ST), path count → ~10, STRG.225.9 emitted at cov ~9.6.

**Hard invariant:** `RUSTLE_PRECISE=1` stays **byte-identical to commit 4705ab1**. The entire port lives behind `!precise_mode()`.

**Guardrail (tracked, not hard-gated):** the 1741 shared chains' `cov` is currently aligned with ST (Pearson 0.95). The port re-apportions flux, so this *will* move; per `feedback_full_stringtie_mimicry` we **accept transient regression per stage** but track shared-chain alignment as the health signal — a *faithful* port should keep/improve it, a *broken* one will wreck it.

---

## 3. The central unknown (Phase 0 — gating diagnostic)

Everything downstream depends on **why rustle's dominant-path flux is 119 where ST's is 143**, given identical `nodecov`/`noderate`. The flux is bounded by an internal junction **edge capacity** = Σ transfrag abundances on the bottleneck edge (all nodes have `nodecov ≫ 143`). Four candidate sources — Phase 0 must attribute the bottleneck to exactly one (or rank them):

1. **Graph node construction.** Confirmed difference: ST keeps node 21 = 32049072-32049150 whole; rustle splits it at junction donor 32049138 into 21+22, and has a node 43 ST lacks. Different nodes → different edges → different capacity matrix. (The flux plateau is constant *through* the split, hinting the bottleneck is source-side, but this must be confirmed, not assumed.)
2. **Transfrag keep-set.** Which transfrags enter the capacity (`cond_on_path`, `keeptr_gap_ok`, CHI_WIN gap test). A different kept subset on the bottleneck edge → different capacity.
3. **Transfrag abundances (read→transfrag granularity).** rustle may split the same reads into more/thinner transfrags, so the dominant transfrag on the bottleneck edge sums to 119 vs ST's 143. This is the documented "flow-input" layer (`project_flow_input_parity`) and the largest-blast-radius possibility.
4. **Ford-Fulkerson augmenting-path policy.** BFS order picks different augmenting paths, distributing the same total flow differently across paths.

**Phase 0 deliverable:** instrument both tools to dump, for the STRG.225 dominant path, the **bottleneck edge** (the min-residual edge that sets 119/143): its node pair, the list of contributing transfrags, and each transfrag's abundance. Compare rustle vs ST. Output: which of (1)–(4) is the root, with the specific divergent quantity.

**Why this gates everything:** the port's entry point and blast radius differ completely by case — (1) is graph-build parity (moderate, upstream), (2) is a keep-gate (small, local), (3) is the read→transfrag rewrite (largest), (4) is the FF solver (medium). We do **not** commit to an approach until Phase 0 answers this.

Tools already in place: `ST_TRACE_COV_NODES` (ST `[ST_COVNODE]`), `RUSTLE_TRACE_COV_NODES` (rustle `[COVNODE]`), `STRINGTIE_PARITY_LOG` / `RUSTLE_PARITY_LOG` (transfrag_define, parse_trflong_seed, seed_reject, graph_edge). Phase 0 likely needs a small **per-edge capacity trace** added to both `long_max_flow` implementations (rlink.cpp:8753 / max_flow.rs:1280).

---

## 4. Approach (contingent phases after Phase 0)

The port is **phased and contingent** — we execute the phase Phase 0 selects, validate against the oracle, then re-scope. Do not pre-commit to all phases.

- **Phase A — Graph node-construction parity** (if root = (1)). Port ST's node-boundary rule (when a junction donor/acceptor falls inside a node, ST's split policy). Re-measure dominant flux. Gate `!precise_mode()`.
- **Phase B — Transfrag keep-set parity** (if root = (2)). Align `keeptr`/CHI_WIN/`cond_on_path` so the bottleneck edge sums the same transfrags. Small, local.
- **Phase C — read→transfrag accumulation parity** (if root = (3)). The big one: port ST's `update_abundance` + `tr2no` collapse so reads merge into the same transfrag granularity. Broadest blast radius; touches all loci. Prior partial attempts ("canonical", `parse_trflong_st`) found this "close-but-broken" — Phase 0 must confirm it's worth it before entering.
- **Phase D — FF augmenting-path policy parity** (if root = (4)). Align rustle's BFS augmenting order/`max_fl` handling to ST's so flow concentrates on the dominant.

Each phase: failing-test-first against the node-flux oracle (rustle dominant flux → 143 at STRG.225), implement minimal change behind `!precise_mode()`, measure chr19 div + shared-chain cov alignment, commit if net-positive or transient-per-stage-acceptable.

---

## 5. Risks

- **Broad blast radius.** Any apportionment/graph change moves `cov` for all 1741 shared chains. Mitigation: the shared-chain cov-alignment metric (§2 guardrail) is computed every iteration; large degradation = the change is unfaithful, revert.
- **Prior "close-but-broken" history.** The `parse_trflong_st` / canonical flow port (`project_flow_parity_scope`) hit back/fwd-extension and free-sink-shortcut bugs and regressed vs annotation. Those were measured vs *annotation*; under the ST-faithful goal the metric flips to vs-ST, but the implementation hazards remain — reuse the existing `_st` scaffolding, don't rebuild.
- **Junction-acceptance entanglement.** If the root is the keep-set/graph (cases 1–2), it touches junction acceptance, where matching ST is a documented precision tradeoff (`project_junction_parity`). Under ST-faithful that's the point, but it will move precision-vs-annotation.
- **Whole-genome OOM** — per-chrom serial only; chr19 is the dev loop.

---

## 6. Validation harness (already built this session)

- **Oracle:** `ST_TRACE_COV_NODES=1 stringtie -L GGO_19.bam` → `/tmp/st_covnode.err` ([ST_COVNODE] per-node nodeflux/nodecov). rustle `RUSTLE_TRACE_COV_NODES=1`. Per-node flux diff at STRG.225 = ground truth.
- **Metric:** `rustle -L GGO_19.bam -o out.gtf; python3 bench/gtf_chain_diff.py out.gtf /tmp/st_all.gtf`.
- **Seed/transfrag parity:** `STRINGTIE_PARITY_LOG` vs `RUSTLE_PARITY_LOG` — seed count, transfrag_define, seed_reject reasons.
- **Shared-chain cov guardrail:** the inline `cov`-alignment probe (Pearson on 1741 shared chains; anchor regex to avoid the `longcov` match bug).
- **Escape-hatch:** `RUSTLE_PRECISE=1 bash bench/mini3/check_precise.sh` (must print "byte-matches 4705ab1") + `RUSTLE_PRECISE=1 cargo test --release --lib` (284/0).

---

## 7. Evidence index (where each claim is proven)
- Filters faithful + cov gap: memory `project_coupled_port_viability` "PRIMARY-LEVER ROOT CAUSE".
- Node-trace substrate match + over-decomposition: "NODE-LEVEL COV-TRACE DIFF".
- Guards falsified: "DOMINANT-FLUX GAP CHASED".
- Seed-merge probe (ST depletes, not merges): "ST-SEED-MERGE PROBE".
- Artifacts: `/tmp/{st,ru}_covnode.err`, `/tmp/st_parity.jsonl`, `/tmp/ru_parity.jsonl`, `/tmp/div_{ru,st}_only.txt`.

---

## 8. Open decisions for the user
1. **Entry commitment:** approve Phase 0 (bounded diagnostic, no source changes to ship) as the first and only currently-scoped step — then re-scope based on which of (1)–(4) it finds? (Recommended.)
2. **Risk posture:** accept transient shared-chain perturbation per stage (full ST-faithful, per `feedback_full_stringtie_mimicry`), with the cov-alignment metric as the revert signal? (Recommended.)
3. **Stop conditions:** if Phase 0 finds the root is (3) read→transfrag (broadest), do we proceed into it or stop at the characterization (it is the documented ceiling)?
