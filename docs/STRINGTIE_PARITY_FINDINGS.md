# StringTie Parity — Authoritative Findings & Ceilings

**Status:** current authority for the StringTie long-read parity effort (de novo, GGO_19 chr19).
**Last updated:** 2026-05-28. Supersedes the precision/parity-plan docs listed at the bottom.

This document is the single source of truth for: where Rustle stands vs StringTie, what the
remaining gaps are, what their *bounded ceilings* are, and which levers are exhausted vs open.
Detailed working trace: `docs/superpowers/specs/2026-05-28-pathpat-flow-parity-scope.md`.

## 1. Current operating point

De novo `rustle -L GGO_19.bam` vs `GGO_19.gtf` (StringTie's output), gffcompare intron-chain:

- **Sn 96.5 / Pr 90.7 / F1 93.5** (1974 transcripts).
- Rustle sits at a **higher-sensitivity / lower-precision** point than StringTie — and that is *partly a feature*: Rustle's flow recovers ~40 real reference isoforms StringTie does not extract.

**Shipped win:** `parity/isofrac-chain-dedup` branch (committed + pushed) — the isofrac multicov
**chain-dedup** fix: count each unique intron chain once in the isofrac denominator (StringTie
collapses dominant 3'-end variants before isofrac; Rustle now matches). +12 tx, Sn 96.1→96.5,
Pr 91.0→90.7, F1 +0.03; recovers STRG.267.3 / 343.2 / 15.1. Default ON, opt out via
`RUSTLE_ISOFRAC_CHAIN_DEDUP_OFF=1`.

## 2. The precision gap, fully characterized (~180 extra chains)

Every extra Rustle emits that StringTie does not was attributed (bench/pathpat_phase0*.py):

| bucket | count | mechanism |
|---|---:|---|
| filter divergence (co-extracted, ST kills) | ~97* | ST extracts the same chain then kills it; Rustle keeps it |
| alt-splice junction selection | ~48 | Rustle's flow generates a different/extra junction; near already-matched refs |
| subset/contained | ~32 | Rustle drops a junction (shorter chain) |
| genuine extra-junction over-generation | ~9 | true over-enumeration |

\* The ~97 "filter divergence" count is **inflated** by a 60bp coordinate-tolerance match in the
attribution; the genuinely filter-recoverable set is ~2 (see §3, bucket-A oracle).

## 3. Lever ceilings — what each fix could achieve (all oracle-bounded)

The reusable method: before building an expensive fix, feed StringTie's *exact* values in as an
**oracle** and measure the F1 ceiling. This killed two false starts cheaply.

| lever | ceiling (oracle) | verdict |
|---|---|---|
| **chain-dedup** | +0.03 (actual) | ✅ shipped — the one clean win |
| **bucket-A lowintron** (read-overhang bpcov) | **+0.1 F1** | ❌ falsified: oracle with ST's exact lowintron masks (93% match) gives only +0.1 Pr / 0 Sn / −2 j |
| **targeted FP-only suppression** | — | ❌ no discriminator: FP over-extractions (min_jct_mm med 4, longcov mostly 1) sit *inside* the real-isoform distribution; any gate that kills 31 FPs kills ~360 real isoforms |
| **flow path-enumeration rewrite** | **+1.1 F1** | ⚠️ open, but it's a **reshape**: Pr +4.2 / Sn −2.2, *loses ~40 real isoforms*; highest-risk, multi-session |

## 4. Abundance values (longcov / TPM / cov)

Measured on 1748 matched chains (bench, corrected GTF parsing — note ST writes `longcov` *before*
`cov`, so naive `cov "…"` regex mis-parses):

- **cov**: median ru/st 0.965, 95% within 2× — **well-aligned**. (The old "60× tlen-inflation" claim
  was largely a parsing artifact on final matched transcripts.)
- **TPM**: systematic **~1.26×**. Identical formula (`cov·1e6/Cov_Sum`); the gap is the denominator —
  Rustle sums final-transcript cov (clean); ST accumulates a per-bundle `sum_cov` over *pre-filter*
  predictions + gene-level sums (ST's own code: *"isn't this double counting?"*). Relative TPMs are
  preserved; ST's denominator is arguably *less* correct. Aligning = copying ST's messier accounting.
- **longcov**: 68% exact, 12% diverge (some extreme, e.g. ST STRG.26.1 longcov 1 vs Rustle 1211).
  **Tried** switching Rustle's longcov from `read_count` to `abundance` (ST's exact field): **no-op** —
  Rustle's transfrag *abundance itself* is 1211 vs ST's 1. The divergence is in read-to-transfrag
  assignment, upstream of longcov. Reverted.

## 5. The unifying conclusion

**Every remaining StringTie-parity gap — precision extras, longcov, TPM — traces to one foundation:
the read-to-isoform / read-to-transfrag assignment in the flow decomposition.**

- cov agrees on matched chains because the *path* is the same.
- longcov diverges where the flow attributes a different read mass to the chain's transfrag.
- the precision extras exist because the flow traverses alt-junction paths ST doesn't.
- the FP and real-minor isoforms are **indistinguishable** by any available feature (cov, longcov,
  min_jct_mm) — confirmed repeatedly.

There is **no cheap filter- or normalization-level lever left**. The only lever with a >1-F1 ceiling
is matching StringTie's read-assignment (the flow rewrite), and that is a sensitivity-for-precision
**reshape**, not a free win — a product decision about which operating point is desired, not a bug fix.

## 5b. Read-to-transfrag assignment — the fix foundation (infrastructure built 2026-05-28)

Since every gap traces to read-to-isoform assignment, the work starts one level down at
**read-to-transfrag** construction. Infrastructure: `bench/transfrag_parity_diff.py` joins both
tools' transfrag events (`transfrag_define` / `transfrag_pre_depl`, both carry abund + intron-chain)
by `(strand, chain)` and reports per-chain abundance divergence.

**What it establishes (full chr19):** 93% of common transfrag chains agree (79% exact). But **159
long multi-intron chains have Rustle over-attributing read mass** — e.g. abund **1270 vs ST 1**, and
`ru_read_count == ru_abund` (1270). ST creates a transfrag per read-pattern (few reads span a long
exact chain → abund 1); **Rustle aggregates partial/spanning reads onto the full chain.** Rustle also
has **6303 transfrag chains ST doesn't** (over-segmentation) vs 196 ST-only.

**Two distinct divergences (session 1 of the rewrite, 2026-05-28):**
1. **159 over-attributions** — Rustle abund 1270 vs ST 1 on long chains. `ru_read_count==abund`, and
   `RUSTLE_ST_TRIM=1` does NOT change it → **not the terminal trim, not loose matching** (Rustle
   matching is exact: `t.node_ids==key && t.pattern==pattern`, map_reads.rs:~1798). Drives the
   **longcov** divergence; mostly a *reporting* difference (both tools still make the chain dominant
   via flow). Lower F1 priority.
2. **6303 Rustle-only transfrag chains** — mostly **low-abundance (73% ≤1 read, 87% ≤2), median 7
   introns** chains ST never creates. Extra distinct chains → extra seeds → extra alt-junction paths
   → the **j-class precision extras**. **This is the F1-relevant divergence.**

**Over-segmentation magnitude (apples-to-apples, transfrag_pre_depl both sides):** Rustle **7335**
multi-intron transfrags vs ST **4383** → **3157 Rustle-only** chains, 205 ST-only (the
`transfrag_define` comparison inflated this to 6303 — use pre_depl both sides, now the diff default).

**Mechanism narrowed by elimination (session 2, 2026-05-28):** deep-dived a Rustle-only 4-intron
chain (`…15670251-15674428`, abund 1). ST's dominant there is longer (`…15670251-15681376,…`). The
weak junction `15670250-15674429` (2 reads) is **accepted by BOTH tools** — so NOT a junction-accept
divergence. Combined with trim ruled out (ST_TRIM no-op) and matching exact, the conclusion is:
**the same reads produce a different transfrag CHAIN in each tool** — the read→node-path construction
(`map_reads.rs::collect_read_nodes_exact` :250 / `split_read_segments` :1308) builds different node
paths per read than ST's `update_abundance` (rlink.cpp:4367). e.g. Rustle ends a chain at acceptor
15674428 where ST's read continues / forms a different chain.

**Session 3 (2026-05-28) — kill-split ruled out, mechanism still elusive.** The deep-dive read
(`m64076…/11600654`) splices through the weak junction AND continues to 15850439, yet Rustle emits a
truncated sub-chain ending at 15674428 → looked like read-splitting. But StringTie has NO
kill-and-split-reads-into-orphans logic (no `killed`/`orphan` in rlink.cpp read handling), and
`RUSTLE_NO_KILLED_SPLIT=1` (pass `None` for killed_junction_pairs at map_reads_to_graph,
pipeline.rs:13405) is a **no-op**: transfrag over-segmentation stays 7335/3157, F1 unchanged. So the
V99 kill-split is NOT the cause (or this isn't the call site building these transfrags — there are
others, e.g. sub-bundle pipeline.rs:15352).

**Four mechanisms now ruled out for the over-segmentation:** junction-acceptance (both accept),
terminal trim (ST_TRIM no-op), add_or_update matching (exact), killed-junction read-split (no-op).
The positive cause is still in read→node-path construction but unidentified.

**Session 4 (2026-05-28) — per-read trace DONE; root found but isolated alignment regresses.**
Added `RUSTLE_TRACE_READ=<substr>` (map_reads.rs, the DEFAULT `map_reads_to_graph_bundlenodes` path —
NOT the fallback `map_reads_to_graph`, which is why earlier gates were no-ops) and ST
`ST_TRACE_READ_START`. Traced read `/11600654`: Rustle splits it into seg0 (`…15674428-15674552`,
the truncated over-seg chain) + seg1 (orphan) at the `15674552→15681376` junction. That junction
(`15674526→15681377`, 1 read) is **rejected by Rustle (`mm_negative`) but accepted by ST**
(nreads=1). Globally Rustle rejects **10542** ST-accepted junctions with reason `mm_negative`
(graph_build.rs:834: `stat.mm < 0.0`; mm=−1 is a higherr/long-read DEMOTION MARKER, not absence of
support — junctions.rs:20 even warns this). When a read uses a rejected junction it splits → over-seg.

**BUT both isolated alignments REGRESS (env-gated, default off):**
- `RUSTLE_NO_KILLED_SPLIT=1` (don't orphan-split at killed junctions, bundlenodes path): Sn 96.5→94.3,
  Pr 90.7→83.5, j 120→270, over-seg 3157→4109. The orphan flag is PROTECTIVE.
- `RUSTLE_KEEP_MM_NEG=1` (accept mm<0 junctions with nreads_good>0, like ST): Sn→93.1, Pr→87.8, j→167,
  over-seg 3157→**3103 (only −54)**. So mm_negative explains only ~54 of the 3157, AND accepting them
  regresses (extra edges → spurious paths).

**Conclusion:** the per-read trace positively identified a root (mm_negative junction rejection → read
split), but (a) it explains a small fraction (~54/3157 — the over-seg is multi-causal) and (b) every
isolated ST-ward change regresses, because Rustle's strict junction filtering is **load-bearing /
protective** given its other divergences. Same interdependence as the whole effort: the
read-assignment/junction system can't be aligned piecemeal — only a coherent multi-stage change
(junction acceptance + flow + downstream filters together) would work, ceiling +1.1 F1 reshape (§3).
Gates `RUSTLE_NO_KILLED_SPLIT` / `RUSTLE_KEEP_MM_NEG` left as documented-regressive env opt-ins;
per-read trace `RUSTLE_TRACE_READ` / `ST_TRACE_READ_START` kept.

## 6. Instrumentation toolkit (kept, env-gated; default behavior unchanged)

- **Cross-tool parity diff**: `RUSTLE_PARITY_LOG` / `STRINGTIE_PARITY_LOG` (+ `_FILTER_CHROM/RANGE/STEPS`)
  emit JSONL decision events (path_extracted, pred_kill, checktrf_*, junction_accept, …) joined on
  `(step,chrom,start,end,strand)`.
- **Coverage/lowintron tracing**: `RUSTLE_TRACE_COV_NODES`, `RUSTLE_RI_TRACE`+`RUSTLE_TRACE_LOCUS`,
  `RUSTLE_BLF_TRACE`, `RUSTLE_JCTMM_DUMP`; StringTie side `ST_TRACE_COV_NODES`, `ST_ILMASK` (rlink.cpp).
- **Oracles**: `RUSTLE_LOWINTRON_ORACLE=<ST mask file>` (feed ST's lowintron in);
  `RUSTLE_LEFTOVER_REDIST`, `RUSTLE_RI_ZERO_INTRON_LOW`, `RUSTLE_LONGCOV_READCOUNT` (opt-in experiments).
- **Analysis scripts**: `bench/pathpat_phase0.py` / `phase0c.py` (extra attribution),
  `bench/pathpat_bucketB.py` (junction-shift characterization).
- **Build note (WSL):** StringTie submodule — `make clean release` (10s/run; DEBUG is ~9min). `/mnt/c`
  mtimes are flaky so `make` may skip recompiling edited `rlink.cpp`; `touch rlink.cpp` or `make clean`.
  Do NOT use `-o /dev/null` with StringTie (temp-dir error) — use a real path.

## Shadow Layer 1 — DONE (2026-05-29)

**Objective:** drive the `mm_negative` bucket (ST-accepted junctions Rustle rejected due to `mm < 0`)
to 0 under `RUSTLE_ST_SHADOW=1`.

**Gate mechanism:** `st_shadow()` predicate OR'd into `keep_mm_neg` at
`src/rustle/graph_build.rs` `filter_junctions_for_bundle` (commit f252e96, branch
`parity/isofrac-chain-dedup`). When shadow mode is on, any junction with `nreads_good > 0` that
would otherwise be demoted by the `mm < 0` marker is accepted, matching StringTie's behaviour
(which never uses the demotion marker as a hard reject).

**Results (chrom NC_073243.2, `-L` mode):**

| Condition | `mm_negative` | `absent_from_rustle` | `strand_mismatch` | Rustle accepted |
|-----------|--------------|---------------------|-------------------|-----------------|
| Shadow OFF (baseline) | 3355 | 6420 | 399 | 7303 |
| Shadow ON (`RUSTLE_ST_SHADOW=1`) | **0** | 6420 | 287 | 10770 |

Layer-1 target met: `LAYER-1 GATE (mm_negative bucket): 0 (target 0)`.

No residual `mm_negative` entries remain — the gate is clean (no nreads_good==0 cases needed
investigation; all previously-rejected mm<0 junctions are now accepted under shadow).

**Deferred to later layers:**
- `absent_from_rustle` (6420 unchanged) — extra junctions ST has that Rustle never sees; root is
  graph construction / bundle-boundary differences (Layer 2).
- `strand_mismatch` (399 → 287, partially improved as a side-effect of accepting more junctions)
  — strand assignment divergence, also Layer 2.

**F1 note:** F1 is NOT measured mid-stack per the implementation plan spec. Default mode
(shadow OFF) remains at Sn 96.5% / Pr 90.7%.

---

## 6b. Shadow Layer 2 — PARTIAL (2026-05-29, commit e6a2a05)

**Goal:** drive the `absent_from_rustle` bucket (junctions ST accepts that Rustle never even
*evaluates* — killed/demoted before the `junction_accept` emit at graph_build.rs:875) toward 0.

**Root cause (investigation):** the absent junctions are NOT a bundling/extraction/coordinate
problem — Rustle bundles them, sees the reads, extracts the junction, then kills/demotes them inside
`good_junc` (`src/rustle/killed_junctions.rs`). The dominant mechanism is a **Rustle-invented `mm=-1`
"all_bad" demotion** (every carrying read has a mismatch in its splice anchor) that StringTie's
`good_junc` (rlink.cpp:13700-13803) does NOT have — ST only sets a local `mismatch` flag and kills
via the long-intron gate. Secondary: Rustle's main witness check defaulted to a 10x long-read
multiplier vs ST's 100x (rlink.cpp:13744), i.e. Rustle was *more* aggressive at killing.

**Fix (shadow only):** under `st_shadow()`, (1) skip the two `higherr_low_support_bad` `mm=-1`
demotion sites (killed_junctions.rs ~716 and ~1239); (2) use ST's 100x witness multiplier.

**Effect (of ST-accepted junctions, junction_accept_diff.py):** `absent_from_rustle` **6420 → 4574**,
`strand_mismatch` **287 → 15**, Rustle-accepted **7303 → 12888**. Layer 1 (`mm_negative`) still 0.
Default (shadow OFF) unchanged at 96.5/90.7.

**Residual 4574 is NOT a good_junc-layer gap (floor reached for this layer):**
- ~735 witness (645 left + 90 right) + ~170 bad_long_intron kills are **ST-faithful** — ST's
  identical checks kill them too; their presence in ST's parity log is a pre-kill candidate /
  coordinate-strand artifact, irreducible at this layer.
- ~1222 `low_splice_frac` fire because **Rustle's `leftsupport`/`rightsupport`/`nreads` differ from
  ST's** (ST's values don't satisfy reason 6, so ST keeps). This is a read→transfrag support-accounting
  divergence → **Layer 3**, not good_junc. Gating the site-1371 defer was TESTED and REVERTED
  (net −10, non-monotonic: keeping them alive at good_junc doesn't make them ACCEPTED — they drop in
  canonicalization downstream).
- ~2361 die **pre-good_junc** (bundle formation / apply_higherr / isofrac mm=-1 sites at :1582/:1833)
  — an earlier stage.
- ~86 `min_support` — a small `nreads_good` accounting divergence (upstream).

**Conclusion:** Layer 2 (good_junc/junction acceptance) is exhausted at 6420→4574. The remaining
absences belong to Layer 3 (read→transfrag support accounting, the `low_splice_frac` 1222 + the
no-trace 2361) — consistent with the layered model. Coordinate-snap (29 HE_DEMOTE) deferred.

## 6c. Shadow Layer 3 — ENTRY MEASUREMENT (2026-05-29, not yet implemented)

Re-baselined `bench/transfrag_parity_diff.py` (transfrag_pre_depl both sides) on the current build:
- **Shadow OFF:** Rustle 7335 chains / ST 4383; Rustle-only **3157** (reproduces the prior baseline exactly).
- **Shadow ON (Layers 1+2):** Rustle 9774 / ST 4383; Rustle-only **5756** (WORSE), in-both 4178→4018.

This is the expected mid-stack regression and validates the shadow-mode premise: Layers 1+2 fixed the
junction SET (kept ~2600 more junctions), giving Rustle's *unchanged* read→transfrag CONSTRUCTION more
junctions to build distinct chains from → MORE over-segmentation. **Layer 3 cannot be judged by
transfrag count until the construction itself is made ST-faithful.** Default (shadow OFF) unchanged.

**Layer 3 scope:** port ST's `update_abundance` (rlink.cpp:4367 — one transfrag per read node-path)
into Rustle's read→node-path construction (`collect_read_nodes_exact` map_reads.rs:250,
`split_read_segments` map_reads.rs:1308, `add_or_update_transfrag` map_reads.rs:1635) under
`st_shadow()`. This is the architectural change the prior 4 read→transfrag sessions identified but
could not do in isolation (see [[abundance-value-alignment-longcov-tpm-cov-vs-stringtie]]). Gate:
transfrag_pre_depl Rustle-only → 0. It is a multi-session effort of its own — scope/plan before coding.

## 6d. Shadow Layer 2.5 — node-boundary divergence characterization (2026-05-29)

**Oracle (bench/node_parity_oracle.py):** within shared bundles under shadow L1+2, node-MATCH
bundles (1916/2171) carry only 40 Rustle-only transfrag chains vs 2383 in the 255 node-MISMATCH
bundles — node parity is the gating lever (60x concentration). PROCEED verdict. (Caveat: ~5353
additional Rustle-only chains live in non-shared *bundles* — a bundle-boundary divergence beneath
node parity; the real global gate stays transfrag_parity_diff.py Rustle-only.)

**Divergence class (bench/graphnode_diff.py, shadow ON, 255 mismatched shared bundles):**
- **SPLIT (222 bundles, 575 ST-nodes split) — dominant.** Rustle subdivides one ST node into 2+.
  Of 1514 internal split coords, **1466 (97%) sit at a coordinate both tools accept as a junction
  endpoint, but 1410 (93%) are at junctions ST DEMOTES** below its support thresholds (84% are
  single-read, nreads_good=1).
- MERGE 4, SHIFT 0, OTHER 33.

**Root cause — Layer 1 over-corrected.** ST's `junction_accept` parity emit (rlink.cpp:16811) runs
BEFORE a second `build_graphs` junction-cleanup pass (rlink.cpp:14065 `nreads_good < DROP/ERROR_PERC`
=5 → strand=0; 14232 `nreads_good<1.25*junctionthr` → mm=-1; 14245 coverage-consensus). Demoted
junctions are skipped in create_graph (rlink.cpp:3577) → no node boundary. Rustle's Layer-1 shadow
keeps these junctions (matching the *log*), so they create node boundaries ST never makes.

**Fix attempted (Task 4, commit eb9e72f):** in `compute_demoted_alt_coords`
(src/rustle/graph_build.rs:1257), under `st_shadow()`, suppress the node-boundary of junctions with
`nreads_good < 5.0` (DROP/ERROR_PERC), `!guide_match`, non-canonical splice (`consleft/consright != 1`),
and coord not used by a surviving junction (kept-coords guard). Default unchanged (shadow-gated; note
the default Intron-chain figure is run-to-run NONDETERMINISTIC ±0.1pp — 90.7↔90.8 — so the regression
guard must not be over-read at that precision).

**RESULT — consensus is a dead end for the genome-less benchmark; true lever is the mm_negative
graph-boundary, not consensus (commit ae17791, predicate corrected):**
- A buggy interim predicate (`consleft != 1`) demoted ALL weak junctions (consleft is always −1 in Rustle
  — genome consensus never wired into stats), giving graphnode 255→211 and transfrag Rustle-only
  5756→4486 but ST-only 365→**821** — i.e. it OVER-merged (traded over-seg for under-seg), NOT real parity.
- Corrected to `consleft == 0` to match ST exactly: ST's `leftcons` is `char` init **−1**, and ST demotes
  `if(!leftcons && nreads_good<5)` — in C `!(-1)` is FALSE, so ST demotes ONLY when consensus is computed
  AND non-canonical (consleft==0). **Our parity runs give NEITHER tool a genome** (`stringtie GGO_19.bam -L`,
  ST has no genome flag → `bdata->gseq` null → leftcons stays −1 → ST never demotes on consensus). So the
  faithful match is `==0`, which is INERT here (consleft always −1) → reverts to L1+2 (255 / 5756 / ST-only
  365). No regression; default 96.5/90.7 (ic nondeterministic ±0.1).
- **The 125→255 SPLIT inflation is NOT consensus.** It is Layer 1's `keep_mm_neg` (graph_build.rs:843)
  propagating into graph-node creation: Rustle keeps mm_negative junctions and creates a node boundary at
  each. StringTie marks them `mm=-1` then `strand=0` in build_graphs A5 (rlink.cpp:~14478) BEFORE
  create_graph, so ST keeps them in the junction_accept LOG (which Layer 1 matched) but EXCLUDES them from
  node boundaries. **TRUE LAYER-2.5 LEVER (next): under shadow, suppress the graph-node boundary of
  mm<0 junctions (e.g. add their coords to demoted_donors/acceptors in compute_demoted_alt_coords keyed on
  stats.mm<0), keeping the edge for read routing.** Open design question: does this conflict with Layer 1's
  read-split-fix intent? (Layer 1 kept mm_negative to stop reads splitting; boundary-suppression keeps the
  edge so routing is preserved — needs verification.) Consensus port (P1, wire consleft/consright from
  genome.rs:76) is only relevant if BOTH tools are run WITH `--genome-fasta` (../GGO.fasta available) — a
  separate parity-configuration decision.

**mm<0 boundary-suppression TESTED & REVERTED (2026-05-29):** added (under shadow) demotion of mm<0
junction coords in compute_demoted_alt_coords. Removed 327 extra node intervals (Rustle-extra 2020→1693)
but bundle gate barely moved (255→249) and transfrag Rustle-only got slightly WORSE (5756→5779, ST-only
unchanged). So this lever (like the nreads<5 and consensus ones) shifts the node set sideways without
converting mismatch→match or reducing transfrag over-seg. **CONCLUSION: three localized node-boundary
levers all fail to converge transfrag parity in isolation.** This is the architecture signal: graph-node
parity cannot be achieved by localized demotion rules measured against the transfrag gate alone — the
node and read→transfrag layers are interdependent (the project's core premise). PATH FORWARD options:
(a) implement Layer 3 (split/trim) AND node changes together, measuring only the combined transfrag gate
(true coherent-shadow approach); (b) genome-enabled parity (run both tools with --genome-fasta, wire P1
consensus). Layer 2.5 left at the safe inert consensus-corrected state (ae17791); no parity gain yet.

## 6e. Shadow Layer 3a — split-removal TESTED & REVERTED; truncations are UPSTREAM (2026-05-29)

Implemented ST's `get_fragment_pattern` split rule under shadow in `split_read_segments`
(map_reads.rs:1339): suppress Rustle's killed/BADJUNC split + single-node-fragment drop (sites at
map_reads.rs:652 and :1199), split only at non-contiguous nodes lacking a boundary-matching read junction.
Default byte-identical (shadow-gated). **Result: transfrag Rustle-only 5756 → 7206 (WORSE, +1450).**
Decisive breakdown: contiguous-sublist TRUNCATIONS of ST chains barely changed (2780 → 2722, −58) while
NON-truncation Rustle-only chains exploded (2976 → 4484, +1508); in-both flat (4018 → 4019). **So the
read SPLIT is NOT the cause of the truncations** — removing it doesn't collapse reads onto ST's full
chains, it lets them chain through more nodes into new divergent chains. **The truncations originate
upstream in read→node-path construction (`collect_read_nodes_exact` map_reads.rs:250 — which node set a
read maps to), not in the split (Layer 3a) nor the node boundaries (Layer 2.5).** Reverted; map_reads.rs
at committed state, no parity change.

## 6f. ARCHITECTURAL CONCLUSION (Layers 2.5+3) — 2026-05-29

The transfrag over-segmentation (5756 Rustle-only) is rooted in **read→node-path construction**
(`collect_read_nodes_exact`: the `rnode`/node-membership a read resolves to), which sits BELOW both the
graph node boundaries and the segment split. Across this effort we tried, each in isolation, and each
shifted the error mode WITHOUT converging the transfrag gate:
- Layer 2.5: weak-junction (`nreads<5`) demotion, splice-consensus demotion, `mm<0` boundary suppression.
- Layer 3a: ST-faithful split + single-node-drop removal.
Node-construction ↔ read-node-mapping ↔ transfrag chains are **circularly coupled** — the project's
"layers reinforce only together" premise, manifesting more strongly than hoped: the layers can't even be
*validated* incrementally because each isolated change moves the metric the wrong way. **The true bottom
of the stack is `collect_read_nodes_exact` (how a read's alignment becomes a node path) vs ST's
`get_read_pattern` (rlink.cpp:4041).** Closing transfrag parity requires matching that construction
(a per-read `rnode` vs `unique_nodes` diff is the entry point; `RUSTLE_TRACE_NODE_MAP` hook exists), and
likely doing so TOGETHER with the node-set and split as one coherent change measured only on the combined
transfrag gate. This is a major sub-effort. **Shipped & safe: Layers 1+2 (junction acceptance). Layer 2.5
inert/safe. No regression to the default (96.5/90.7, ±0.1 nondeterministic).**

> **CORRECTION (see §6g): the "root is collect_read_nodes_exact / circularly coupled" conclusion above
> is WRONG.** A per-read node-map trace proved Rustle's `collect_read_nodes_exact` is EQUIVALENT to ST's
> `get_read_pattern` (both build full paths). The real root is a single MISSING step — ST's transfrag
> containment-collapse — described next. It is NOT circular; it is cleanly localized.

## 6g. Shadow Layer 3 — ROOT FOUND: missing transfrag containment-collapse (2026-05-29)

Per-read node-map tracing (collect_read_nodes_exact vs get_read_pattern on truncation cases) DISPROVED the
node-mapping hypothesis: **both tools build identical, full-length node paths per read.** The truncated
Rustle-only chains are *real physically-short reads* (genuinely fewer junctions, verified by exact-CIGAR
match — e.g. 92 reads that truly stop after intron 13 of a 14-intron ST chain). 

**StringTie folds these contained reads into the longer compatible chain BEFORE the transfrag_pre_depl
snapshot, via the keeptrf containment-collapse loop (rlink.cpp:5588-5800, inside build_graphs):**
1. `transfrag.Sort(longtrCmp)` (5589) — most abundant/complete first, so full chains enter `keeptrf` first.
2. For each later transfrag, `compatible_long(t, …)` (5716; defined 5267) returns 1 (t1 has extra intron
   past t2), 2 (t2 extends past t1 → t1 is a prefix), or 3 (compatible ends).
3. Prefix `case 2` (5733): if non-guide and end-gaps `len[1]<ssdist && len[3]<ssdist`, set `included=true`,
   fold mass `keeptrf[t2].cov += abundance; keeptrf[t2].group.Add(t1)` — the prefix is NEVER added to
   keeptrf/trflong → absent from transfrag_pre_depl.
4. Guard (5780): a non-included transfrag lacking longstart/longend (not hardstart/hardend) → `weak=1`, dropped.

**Rustle has NO equivalent pre-seed collapse.** `transfrag_pre_depl` (path_extract.rs:6656-6690) emits every
`trflong_seed` directly; `parse_trflong` (path_extract.rs:5863) only orders/filters by weak/usepath, never
absorbs contained prefixes. So each distinct read path survives → over-segmentation.

**Truncation shape (confirms containment):** of 2780, 1018 are 5'-only, 1471 3'-only, 291 both-ends, 0
internal — all strict prefix/suffix/infix of an ST chain, exactly what compatible_long collapses.

**FIX (the real Layer 3): port the keeptrf containment-collapse loop under st_shadow()** into the
seed-selection stage (path_extract.rs around parse_trflong / the trflong_seed collection feeding
transfrag_pre_depl). Building blocks EXIST: `compatible_long(tf1,tf2,graph) -> (u8,[i64;4])`
(transfrag_process.rs:932, ST's exact return convention) + longstart/longend/weak/trflong_seed/usepath
fields on GraphTransfrag. New code ≈ abundance-sorted (longtrCmp) outer loop + case 1/2/3 absorb + weak
guard (~150-250 LOC mirroring rlink.cpp:5588-5800). Run BEFORE the transfrag_pre_depl emit and parse_trflong.
**Risk:** medium — exact longtrCmp tie-break + constants (ssdist/edgedist/DROP/longintronanchor) must match
ST; depends on longstart/longend/hardstart/hardend being set ST-faithfully (validate alongside).
**Impact ceiling:** all 2780 truncations (~49% of 5756 Rustle-only) are this ONE mechanism.
**ST trace note:** ST_TRACE_READ_START only fires in --merge mode (get_read_to_transfrag); the -L path
(get_read_pattern via get_fragment_pattern) has no trace hook — algorithm verified by reading + chain-count
deltas, not the hook.

---

## 7. Superseded documents

The following are superseded by this file for the precision/parity-gap analysis (kept for history):
`EXACT_STRINGTIE_PARITY_PLAN.md` (its "100% parity via guided mode" premise is circular — `-G GGO_19.gtf`
scored against `GGO_19.gtf`), `PRECISION_DIAGNOSTIC_FINDINGS.md`, `PRECISION_IMPROVEMENT_PLAN.md`,
`PRECISION_WORK_CONTINUATION.md`, `docs/STRINGTIE_PARITY_SYSTEMATIC.md`,
`docs/STRINGTIE_LONGREAD_PARITY_CONTINUATION_GUIDE.md`, `docs/SUBBUNDLE_PARITY_HANDOFF.md`.
