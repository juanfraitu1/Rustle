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

## 6h. DECISIVE: keeptrf collapse ALREADY EXISTS; shadow is F1-CATASTROPHIC on final output (2026-05-29)

Two corrections to §6g, from reading the Rustle code + measuring FINAL output (not the pre_depl proxy):

1. **Rustle ALREADY has the keeptrf containment-collapse** — `process_transfrags` (transfrag_process.rs
   ~2050-2280) is a close port of rlink.cpp:5588-5800 (guide handling, compatible_long case 1/2/3,
   leftdist/rightdist). The `transfrag_pre_depl` emit (path_extract.rs:6662) fires on EVERY `trflong_seed`
   WITHOUT the `weak==0 && usepath>=0` filter that parse_trflong applies — so it snapshots PRE-collapse
   seeds, while ST's transfrag_pre_depl is POST-collapse. **The 5756 "Rustle-only" was ~3.8x inflated by
   this staging mismatch** (collapse removes ~4250). So "port the keeptrf loop" is MOOT — it exists.

2. **FINAL-output measurement (the metric that matters), shadow ON vs ST vs truth:**

   | | transcripts | Rustle-only chains vs ST | Sn/Pr vs truth (GGO_19.gtf) |
   |---|---|---|---|
   | Default (shadow OFF) | 1973 | 219 | **96.5 / 90.8** |
   | Shadow ON (Layers 1+2) | 3049 | 1506 | **85.1 / 62.9** |

   **The shadow path is F1-CATASTROPHIC: precision 90.8 → 62.9.** Layers 1+2 ("ST-faithful junction
   acceptance") flood the graph with the weak/mm_negative junctions Rustle's default deliberately rejects;
   the existing keeptrf collapse + flow + filters cannot recover it → ~1100 spurious transcripts.

**CONCLUSION ON THE SHADOW STRATEGY.** The premise was "make every layer ST-faithful together; F1 only at
the end." We've now shown: (a) the upper layers (junction acceptance) make the FINAL output far WORSE
(90.8→62.9 Pr), because Rustle's downstream is tuned to Rustle's stricter junction set; (b) the layers do
NOT converge incrementally (each isolated change shifts the error mode); (c) the pieces we thought were
missing (keeptrf collapse) already exist. Reaching bit-exact parity would require making the ENTIRE
pipeline (flow allocation, all filters, abundance) ST-faithful simultaneously — a full reimplementation —
and the PAYOFF is only "Rustle can reproduce ST" (a proof), NOT better F1: Rustle's DEFAULT already matches
ST-level precision (90.8). The shadow mode is a research/proof instrument, not an F1 lever. **Shipped &
safe: Layers 1+2 (default-OFF, junction-acceptance parity proven). Default operating point untouched at
96.5/90.8. Tooling (junction_accept_diff, graphnode_diff, node_parity_oracle, transfrag_parity_diff) and
the full layer-by-layer root-cause map are the durable deliverables.**

## 6i. Layer 4 (FLOW) localized — the next port for full-pipeline parity (2026-05-29)

Decision: pursue full-pipeline reimplementation ("start as StringTie, then improve"). The HONEST gate is
final-output chain parity (`bench/gtf_chain_diff.py`): shadow ON = **1506 Rustle-only** final chains
(target 0). Localized via path_extracted / pred_kill diff:

| category | count | % |
|---|---|---|
| **FLOW divergence** (in Rustle path_extracted, ABSENT from ST path_extracted — ST never extracts) | **1387** | **92.1%** |
| FILTER divergence (in BOTH path_extracted; ST kills via pred_kill: retained_intron 88 / isofrac 42) | 113 | 7.5% |
| post-extraction | 6 | 0.4% |

**So 92% is FLOW, not filters.** Characterization of the 1506: median cov 1.08 (75% in [1,2) — coverage is
NOT the lever), 1–53 introns, 785 loci (240 over-enumeration hotspots), 67% alt-variants of an ST locus +
33% novel. Rustle flow-div = 1077 main-flow + 310 checktrf_rescue.

**Root cause:** ST `long_max_flow` (rlink.cpp:8471) tail (~8634-8665) SUBTRACTS allocated flow from every
participating transfrag's abundance, clamping to 0 (node-flow conservation). After a dominant path is
extracted, sibling seeds over shared nodes get flux=0 → fail the store gate (rlink.cpp:9807/9917) → fall to
checktrf, which in `!mixedMode` long-read mode rarely stores (rlink.cpp:9975). This enforces ~one path per
coverage unit. **Rustle's per-seed depletion (path_extract.rs main loop ~6507) does NOT reproduce ST's
flow-conservation strength → siblings keep nonzero flux → get stored; and Rustle's checktrf-rescue store
gates (path_extract.rs:10146/10162) are default-OFF → +310 ST drops.** (Confirmed not the edgecov toggle:
RUSTLE_NO_EDGECOV_DEPL=1 moved 1506→1505 — the divergence is the depletion MAGNITUDE / store gate inside
the flow.)

**NEXT = port Layer 4:** bring Rustle's per-seed flow depletion into agreement with ST's `long_max_flow`
abundance-subtraction-to-zero (rlink.cpp:8471, tail ~8634-8665) so post-dominant-path siblings return
flux=0 and demote to checktrf; pair with an ST-restrictive checktrf store gate. Targets the 1387 flow-div
(92%). The coverage/abundance divergence ("quantification extremely off": ST cov median ~5663 read×len
scale vs Rustle ~1.0 per-bp) is a SEPARATE Layer-6 concern, implicated only for the 113 filter-div.
Key files: ST rlink.cpp (parse_trflong 9693, long_max_flow 8471, store gate 9807/9917, checktrf 9975);
Rustle path_extract.rs (parse_trflong 5863, main loop ~6507, checktrf store 10248, gates 10146/10162),
global_flow.rs, max_flow.rs.

## 6j. Transfrag-PROVENANCE parity layer ADDED (2026-05-29, Rustle b0d13ed / ST fork 160e5e5)

Closed the provenance instrumentation gap. New default-OFF parity event `transfrag_collapse` emitted by
BOTH tools at the keeptrf containment-collapse: key (chrom, rep_start, rep_end, strand); payload
`rep_introns, group_cov, n_members, members="chain:abund;..."` — i.e. for each kept representative, which
transfrags were folded in and their abundances. Rustle: process_transfrags (transfrag_process.rs, after
keeptrf finalized; threaded bundle_chrom/strand from pipeline.rs:13747). ST: rlink.cpp:5803 keeptrf→trflong
loop. Diff tool `bench/transfrag_collapse_diff.py` joins by (strand, rep_introns). ALSO fixed the
`transfrag_pre_depl` STAGING mismatch (path_extract.rs:6662): now filters `trflong_seed && weak==0 &&
usepath>=0` so Rustle emits POST-collapse like ST (was pre-collapse — the source of the old 3.8x-inflated
5756 proxy). Default unchanged 96.5/90.7; additive/gated.

**BASELINE (shadow ON) — provenance divergence now directly visible:** Rustle **9776** keeptrf reps vs ST
**3792** (2.6x more — Rustle's consolidation is much WEAKER, keeps more distinct seeds); 2409 reps in both
(chains join byte-identically), 7367 Rustle-only reps, 1383 ST-only. On the 2409 COMMON reps: group_cov
matches within 5% only **50.1%** (1208), n_members matches **54.2%** — i.e. even where both tools keep the
same representative chain, they fold in DIFFERENT abundance/members half the time (e.g. ru_cov 585 vs
st_cov 3 on a 3-intron chain). This is the direct provenance view of the over-extraction: Rustle under-
consolidates (2.6x reps) AND mis-attributes group abundance. This is the gate for any future read→transfrag
/ keeptrf-collapse parity work (drive Rustle-only reps → 0 and group_cov/n_members match → 100%).

---

## Predcluster selection-parity (sub-project 1) — DONE, opt-in

Status: **DONE (2026-05-30), committed, default OFF.** The capability is the foundation for sub-project 2.

### What was built
`RUSTLE_PREDCLUSTER_ST=1` (default OFF) swaps Rustle's per-cluster winner selection for an
ST-faithful selection in `select_predictions_st` (`src/rustle/predcluster_st.rs`), dispatched
from `transcript_filter.rs::print_predcluster_with_summary_multi` (fully guarded; default-OFF path
byte-unaffected — the only non-test changes to existing files are the guarded dispatch + a few
`pub(crate)` visibility bumps). The flag predicate `stringtie_parity::st_predcluster()` is
unit-tested (`st_predcluster_default_off`), and each sub-stage has its own unit test (9 tests).

The selection runs ST's sub-stages in ST's order:
1. **Survival/coverage gate** — cov + read-through gate; drops low-cov / too-short chains, guides exempt.
2. **Pairwise containment** — ST `rlink.cpp:7363-7404` + `has_retained_intron`: contained-no-RI ⇒ kill
   unconditionally; contained-with-RI ⇒ kill only if `cov_i < cov_j`. Plus **`included_drop`**, which
   requires ST's **strict `included_pred`** (the coverage-free structural-inclusion variant), not Rustle's
   looser containment.
3. **Per-maxint isofrac** — ST `rlink.cpp:18734-18794` `longunder` loop, seeded per maximal-coverage
   interval; kills `cov < isofraclong·usedcov[s]`, guides exempt.
4. **Significant-overlap matrix gate** on the pairwise stage (`build_significant_overlap_matrix`,
   OvlTracker-like sweep) — recovers TPs that unconditional containment over-kills.

### Key findings
- **ST's predcluster retained-intron handling is lowintron-gated** at `rlink.cpp:18528 → 17117`,
  **NOT** at the 7363 site. RI demotion only fires for chains flagged low-intron; modeling it at the
  bare 7363 containment site over-kills.
- **`included_drop` needs ST's strict `included_pred`** (structural, coverage-free). Using Rustle's
  looser containment here mis-drops valid alternative chains.
- The **significant-overlap matrix gate** is load-bearing: without it, unconditional contained-no-RI
  kills remove true positives that ST keeps because their overlap is not "significant."

### Numbers (chain-level, vs `../GGO_19.gtf`; multi-intron chains keyed (strand, intron-tuple))
| Config | TP | FN | FP | Sn | Pr | F1 |
|---|---|---|---|---|---|---|
| **Baseline (flag OFF, shipped default)** | 1750 | 64 | 168 | 96.5 | 91.2 | **93.78** |
| **Flag-ON `RUSTLE_PREDCLUSTER_ST=1`** | 1734 | 80 | 155 | 95.7 | 92.1 | **93.65** |
| **Oracle ceiling (ST winners on matching clusters)** | 1762 | 52 | 123 | 97.0 | 93.5 | **95.27** |

Flag-ON is **deterministic** (3 identical runs, identical chain sets). Default (flag OFF) is unchanged.

### Selection-cluster convergence (candidate-matching clusters only, where Rustle's and ST's
`path_extracted` candidate sets are identical — selection isolated from extraction divergence;
495 matching clusters, 1328 reference chains in them):
| Config | TP | FN | FP | selection errors (FP+FN) |
|---|---|---|---|---|
| Baseline | 1306 | 22 | 63 | **85** |
| Flag-ON | 1295 | 33 | 43 | **76** |
| Oracle ceiling | 1318 | 10 | 18 | **28** |

Flag-ON drops FP **63 → 43 (−20)** but raises FN **22 → 33 (+11)**, a net **−9 selection errors**
(85 → 76), closing ~16% of the baseline→oracle gap (85→28). The FP reduction is the genuine selection
win; the FN rise is coverage-input divergence (below).

### Decision: **opt-in (default OFF).** Per the plan, flip to default-ON only if flag-ON F1 > baseline
93.78 with Pr up and Sn not hurt. Flag-ON F1 = **93.65 < 93.78** (Pr ↑ 91.2→92.1, FP ↓ 168→155, but Sn
↓ 96.5→95.7) ⇒ F1-neutral-to-slightly-negative standalone ⇒ **keep `RUSTLE_PREDCLUSTER_ST` opt-in.**
No default change. The capability is committed and is the foundation for sub-project 2.

### Coupling conclusion (why standalone is F1-neutral)
The ST selection rules reproduce ST's **OUTPUT** only when fed ST's **COVERAGE INPUTS**. Run on Rustle's
flow-based coverage, the same rules kill ~42 contested-minority chains because Rustle's flow cov is below
ST's on those chains (the per-maxint isofrac and cov gate see a different `usedcov`/`cov` ratio than ST
does). Those 42 residual selection-FN are **coverage divergence**, not selection bugs; another **38 FN are
never-extracted candidates** (ST extracts them, Rustle's flow never produces them). Both pools are
**sub-project 2 (candidate-extraction + coverage parity)**. The selection logic is correct (Pr exceeds
baseline, FP below baseline); realizing the **+1.49pp oracle ceiling (95.27)** requires sub-project 2 to
give selection ST-equivalent coverage and candidate sets.

---

## 6j. Mechanism audits + color/CGroup divergence (2026-05-30)

Three byte-level mechanism audits and one shipped fix on branch `parity/isofrac-chain-dedup`.

**SHIPPED (default, inert):** `get_min_start` tie-break `bundle_builder.rs:399` changed `<`→`<=`
so start-coordinate ties resolve to the higher `sno`, matching ST's pairwise cascade
(`rlink.cpp:627-668`). Provably inert on GGO_19 (272/272 tied cross-strand pairs overlap; 0
order-sensitive multi-group patterns). Bench unchanged 95.6/90.5 tx. Pure faithfulness win.

**Audits — all faithful (no reachable -L divergence):**
- *Three strands* (+/-/`.`): aligned across per-strand graphs, unstranded `neg_prop`
  dual-placement, junction strand, cross-strand demotion. 0 strand-flipped final chains. One
  cosmetic gap: Rustle reads only minimap2 `ts` (`bam.rs:333`); ST reads `XS` first
  (`GSam.cpp:338`) — invisible on minimap2 data, bites only XS-bearing BAMs.
- *compatible_long* `len[0..3]`: byte-faithful (`transfrag_process.rs:932` vs `rlink.cpp:5267`).
  Half-open↔inclusive convention (`Rustle.end = ST.end+1`) absorbs every `+1`; gates match
  (ssdist=25, edgedist=100, DROP=0.5). ret==2 work (febb383) did not desync it.
- *Color union-find / CGroup*: logic faithful (min-index rep, path compression,
  eqposcol/eqnegcol, set_strandcol). But output colors are NOT identical — see below.

**Key finding — we do NOT have ST's colors (interrupted colors):**
- NOT mm_negative (hypothesis refuted): the mm_negative reject (`graph_build.rs:845`) filters the
  junction from the graph EDGE list but never sets `JunctionStat.strand=0`. All 3467 RU
  mm_negative junctions keep `jstrand=±1`; zero RU junctions have `jstrand==0`.
- Real root: production call `pipeline.rs:14819` passes `junction_stats=None`, so the color-break
  test is `!has_good_left` (`bundle_builder.rs:137`, junction ∉ good_junctions) instead of ST's
  `strand==0` (`rlink.cpp:1069`). `has_good_left` is strictly narrower → Rustle over-breaks color.
  2510/2536 RU break sites (99%) land on a junction ST keeps with strand ±1.
- Impact tempered: union-find re-merges most breaks (only 11 of 125 boundary-mismatch bundles are
  true splits). Traced FPs: RSTL.103.1 (class p), RSTL.331.2 (class j); a third split benign.
- Entangled — literal ST-faithful regresses: `RUSTLE_SUBBUNDLE_USE_STATS=1` (ST's strand==0 break)
  worsens parity (bundles 3405→2569 under-segmenting; mismatches 125→155). The `None`/
  `!has_good_left` default is a tuned compensation: Rustle over-breaks color AND over-filters
  mm_negative edges; the two partially cancel. True faithfulness needs BOTH halves aligned, which
  lands back in the mm_negative floor (§6h, F1-catastrophic). Both-halves experiment scoped
  separately (design doc 2026-05-30).

### Phase 0 results (2026-05-30) — junction-set parity + segmentation prize-bound → ABORT

The scoped experiment's Phase 0 (analysis-only; new env-gated `junction_raw` + `good_junction`
emits in both tools, `bench/junction_set_diff.py`, `bench/segmentation_prize.py`) ran. Result:
**ABORT — Phase 1 (the two-halves color+mm_negative alignment) is shelved as F1-irrelevant.**

**Junction-set parity (the definitive "do we have ST's junctions?" answer):**
- Layer 1 RAW: shared=17464, only-RU=2, only-ST=1004 — the tools nearly agree on observed junctions
  (ST sees ~1004 RU never does).
- Layer 3 good_junction membership: shared=7307, **only-RU=0**, **only-ST=10170**. Rustle's kept set
  is a STRICT SUBSET of ST's; ST keeps 10170 junctions Rustle drops (and Rustle keeps none ST drops).
- Divergent (ST keeps, RU drops, both observed) = 10157: `not_in_ru_accept`=6403, `mm_negative`=3355
  (matches the known shadow Layer-1 figure), `strand_mismatch`=399. This is the support-accounting +
  mm_negative floor (§6h), not a new mechanism.

**Segmentation prize-bound (the abort gate):** net reachable FP+FN from perfect ST segmentation =
**0**. The tool first reported 6, but adversarial verification found (a) a load-bearing bug — the
split detector's containment filter dropped the real spanning RU bundle before the spanning guard
ran, miscounting single-exon side-bundles (which ST emits too) as fragments → 11 phantom splits; and
(b) causal prize 0 — every attributed FP is within-bundle over-enumeration (one, RSTL.284.2, is in
ST's own output) and every FN is alt-splice / splice-precision (e.g. a 4bp acceptor shift), none
caused by segmentation. After the fix (`bench/segmentation_prize.py`, commit 5801f51): **0 true
splits, 0 FP, 0 FN.**

**Conclusion:** Rustle's bundle segmentation already MATCHES ST's at every divergent locus — it has a
spanning bundle there; the "125 mismatched bundles" (graphnode_diff) are single-exon side-bundle
differences + spanning-overshoot, NOT fragmentation. So aligning colors/mm_negative to ST would fix
nothing in the final output. The color divergence (§6j) is confirmed **structurally negligible for
F1**. Residual at those loci is isoform over-enumeration + splice-site precision (orthogonal to
segmentation). The mm_negative junction divergence (10170 only-ST good junctions) is real but is the
known F1-catastrophic floor (§6h), not a segmentation lever. Phase 0 instrumentation
(`junction_raw`/`good_junction`) is committed and default-OFF for future use.

### Foundation parity audit (2026-05-30) — "same standpoint as ST?" = YES (4/5 layers measured)

Consolidated multi-agent audit (5 layers + adversarial verify) confirming Rustle starts from ST's
graph substrate before the known flow/transfrag divergence:
- **Bundle envelopes**: 3351/3430 exact (98.4% of RU); 25-bundle gap = benign within-locus
  merge/split reshaping; 10 wholesale-missing bundles all 0-tx / 0-ref. Substrate-equivalent.
- **Bundlenodes**: 3226/3351 byte-identical (96.3%); 125 mismatched = micro-node spacers + boundary
  jitter + junction-coupled splits. CAVEAT (adversarial refute): ≥1 GENUINE independent
  node-construction policy diff — ST materializes intron-interior nodes from opposing-strand
  junctions where RU leaves a gap (42 nodes; bundle 20117101-20120045 adds (20117599,20117788) under
  byte-identical good-junctions). DOWNSTREAM-INERT (identical final tx RSTL.32.1==STRG.76.1).
- **Graph edges**: now DIRECTLY MEASURED (2026-05-30) — `graph_edge` parity emit wired in both tools
  (RU pipeline.rs commit 769149f; ST rlink.cpp submodule b18c4b8; bench/edge_diff.py commit 30e3a77).
  Edges emitted as node-coordinate-pair tokens (FROM>TO, SRC/SNK for source/sink), comparable
  cross-tool. RESULT: on the 3226 identical-node shared bundles, edge sets identical on 3037/3226
  (94.2%). **Junction-spanning connectivity = same standpoint**: RU is missing exactly ONE junction
  edge vs ST on an identical-node bundle, and it is a CORRECT omission (alt-donor 16912473 absent from
  reference; ST builds+emits it as non-ref STRG.24.4, RU emits only the correct-donor chains — RU more
  precise). The other 60 RU-missing + 161 RU-extra junction edges are all on node-mismatched bundles
  (different decomposition, expected). The real edge-layer divergence is **source/sink (terminal)
  wiring**: RU has 2708 EXTRA source/sink edges (1308 on identical-node bundles, near-perfect superset,
  only 1 missing) — same coverage-drop threshold (ERROR_PERC*DROP=0.05) but RU uses a per-node
  coverage-ratio test + phantom-zero recursion (graph_build.rs:68-272) vs ST's structural-terminal +
  intra-exon cliff find_trims (rlink.cpp:1500-1601). This is the KNOWN coverage-trim / flow-depletion
  divergence (more permissive terminal wiring → more path endpoints → feeds over-enumeration), i.e. the
  edge layer hands off cleanly to exactly the downstream flow stage. NOT a connectivity bug.
- **Colors**: junction-clean spans identical 1609/1610 (99.94%); 124/125 partition mismatches span an
  only-ST good junction. No independent color-logic divergence (fully reducible to junction floor).
- **Junctions**: RAW = same standpoint — RU raw observes 6456/6456 ST-final + 6396/6396 reference
  introns (RU raw is a SUPERSET on the real-intron axis); only-ST raw=1003 are inert off-by-1
  strandless ghosts (0 in final, 0 in ref). good-set = strict subset (only-RU=0), the known floor.

**Propagation guarantee:** segmentation_prize = 0 true-split / 0 FP / 0 FN; every FN ref-chain has
all introns in RU good_junction (missing_jx=[]); every FP intron resolves to junctions good in BOTH
tools. NOT ONE final FP/FN is attributable to any substrate divergence — the problems are genuinely
downstream (flow / path-enum / transfrag selection) on a shared substrate. Verdict: NEAR_IDENTICAL
across all five layers (not literally byte-identical only due to the inert node-policy diff + raw
ghosts + strict-subset good-set). Lone instrumentation gap: graph edges (proved, not measured).

### Source/sink terminal-wiring experiment (2026-05-30) — ABORT, terminal edges are read-backed not FP-causing

The edge audit flagged terminal (source/sink) wiring as the first RU↔ST divergence (RU 2708 extra
SRC/SNK edges) and the on-ramp to flow. Oracle-first experiment (spec + plan 2026-05-30):

- **Phase 0a free probe:** `RUSTLE_COVLINK_RECURSE_ZERO_OFF=1` (disable RU's phantom-zero recursion)
  changed NOTHING — extra edges stayed 2708, F1 unchanged. The phantom recursion is inert; the extra
  edges come from RU's blanket per-node coverage-drop rule, not the recursion.
- **Phase 0b injection oracle** (`RUSTLE_TERMINAL_ORACLE`, default OFF; `terminal_oracle.rs` +
  `bench/terminal_oracle_report.py`): on identical-node bundles, override RU terminal edges with ST's
  captured set. Result: **NET −4 (FP_removed 1 − TP_lost 3 − FP_added 2), F1 93.89→93.78** (Tx
  95.6/90.5→95.4/90.4). Removing terminal edges to match ST hurt BOTH Sn and Pr → **ABORT**.
- **Verified robust** (adversarial): the oracle fires correctly (removes 1816 SS edges, 3289/3405
  bundles, 62 legit skips); only 518 net-dropped because two downstream stages RE-DERIVE terminal
  connectivity from read evidence after the hook — prune source/sink auto-attach
  (`graph_build.rs:3466-3483`, +121) and read-to-graph mapping (`map_reads.rs:1015`, +1177). The
  re-added edges are MORE read-supported (that's why they're re-established), so fuller convergence
  would lose MORE TPs, never reverse the sign.

**Conclusion:** RU's "extra" terminal edges are **read-backed real-transcript endpoints, not FP
sources**. Terminal wiring is NOT the over-enumeration lever and is NOT independently fixable —
terminal connectivity is re-derived from reads at the prune, read-mapping, AND flow
(`build_lr_edge_capacities`) stages. The FP floor remains the flow-depletion / coverage-divergence
issue (§6h, `project_coverage_metrics_deviation`), not graph-construction terminal wiring. The first
RU↔ST divergence point is correctness (read-backed endpoints), confirming the foundation substrate is
sound through to flow. `terminal_oracle.rs` + report tool are committed, default-OFF, reusable.

### Retained-intron scope oracle (2026-05-31) — PROCEED (+25 prize), fix = killer-PAIRING not flags

Needy-loci trace (§ memory `project_needy_loci_decision_trace`) found over-enum FPs are 3 targetable
hidden ST decisions; lever #1 (retained_intron) tested oracle-first (spec+plan 2026-05-31). Added ST
intron-chain emit to `pred_kill` + `pred_intron_low` (submodule 1bff461; byte-verified vs RU key).
- **Phase 0a prize:** of 1572 ST `retained_intron`-killed chains, 29 match RU final chains → 27 RU-only
  FP + 2 TP-cost → **NET +25**. First gate-clearing prize this session (addresses 27/158 RU-only FPs).
- **Phase 0b diagnosis:** forcing ST's lowintron FLAGS via the existing `RUSTLE_LOWINTRON_ORACLE`
  (oracle confirmed fired) → FP-reduction only **2**. ≪ 25 → the under-firing is the **killer-PAIRING
  scope** (`transcript_filter.rs:3063-3071`), NOT `build_lowintron_flags`. RU pairs victims only against
  a single full-length killer, not local high-coverage killers.
- **Phase 1 (pairing fix) FAILED — diagnosis corrected.** Precise locus trace (RSTL.398.17) showed the
  real mechanism is NEITHER flags NOR pairing-scope: the local high-cov killers ARE already paired (9 of
  them) and the lowintron flags ARE correct. The miss is the **`retainedintron_like` PREDICATE
  GEOMETRY** — victim's first exon lies fully INSIDE the killer's low-cov intron and ENDS BEFORE the
  killer's next exon (a spurious-donor-within-intron the `j==0` logic skips). The prescribed pairing
  widening realized only **−3 FP / −4 TP = net-negative** → REVERTED (commit bb871fa; `RUSTLE_RI_LOCAL`
  predicate kept, reusable). Default restored to 96.2/91.7, 95.6/90.5.
- **Real fix = a `retainedintron_like` predicate change** (add a `j==0` contained-first-exon kill) — a
  NEW mechanism beyond this plan, with a real OVER-KILL risk: it fires by geometry, not ST's exact
  killed list, and the −4 TP collateral in the pairing attempt warns that broadening the kill costs real
  isoforms. The +25 oracle ceiling (ST's EXACT kills: 27 FP + 2 TP) is the bound; a geometry predicate
  may kill a broader set → more TP-cost. Needs its own oracle-first sub-experiment (bound: does the
  contained-exon predicate realize ~+25 at ≤~2 TP, or over-kill?). Tools `bench/retained_intron_prize.py`
  (047e5de), `bench/build_lowintron_oracle.py` (d698d1a). Lesson reaffirmed: 0b ruled OUT flags but the
  "therefore pairing" inference was wrong — the per-locus trace found the true (predicate) cause.
- **FINAL: SHELVED — the +25 prize is NOT realizable by any RU-side predicate.** Two more oracles bounded
  the predicate space: (a) contained-exon GEOMETRY (`retained_intron_geometry_oracle.py`, ce8a0c0) →
  net **−114** (kills 133 real alt-TSS isoforms), captures 5/27 — oracle-first SAVED a catastrophic
  regression; (b) faithful CHAIN-SUBSET + dominance (`retained_intron_chainsubset_oracle.py`, ff87ebf) →
  best net **+0** (fires on nothing at dom≥20, or kills 9 TP/1 FP at dom≤5), captures 1-2/27. DECISIVE
  STRUCTURAL REASON: **18/27 of ST's FP victims have NO chain-superset killer in RU's FINAL transcript
  set** — ST kills them with predictions that exist only in its PRE-FILTER population and don't survive
  to RU's final output. The "114.8x dominance" was an ST pre-filter measurement that doesn't transfer.
  So retained-intron parity needs ST's upstream prediction POPULATION, not a downstream predicate →
  loops back to the flow/coverage-depletion floor (§6h, `project_coverage_metrics_deviation`). All
  retained-intron code reverted; flag + bench oracles kept. Lesson: an oracle PRIZE (matching ST's exact
  final-absent chains) is a CEILING, not proof a downstream fix can realize it — verify the killers
  exist in RU's population before building.

### Lever #2 — donor-snap (2A) + chimeric-suffix-rescue (2B) (2026-05-31)

Trace lever #2 ("transfrag-collapse snap/fold") mechanism-verified into TWO independent sub-levers,
NEITHER in transfrag-collapse (3rd mechanism mis-attribution caught by grounding): 2A = junction
higherr donor-snap scope gap (`junction_correction.rs:178`); 2B = chimeric_suffix_rescue extra RU-only
mechanism (`path_extract.rs:2992/9485`). Spec+plan 2026-05-31-prediction-population-snap-fold.
- **2B SHIPPED (opt-in, REALIZED):** `RUSTLE_CSR_FOLD` (predicate 49d889f, guard 299afa7) suppresses the
  rescue when the 5'-truncated suffix is contained in a kept full-length FLOW path (computable via
  `flow_kept_paths_len` + `out_idx`, no new state). Flag-ON: RU-only FP **158→155**, IC Pr 91.7→**91.8**,
  Tx Pr 90.5→**90.6**, **Sn flat 95.6, 0 TP cost**. Default OFF byte-faithful (96.2/91.7, 95.6/90.5).
  csr_classify: 3 CSR-FP / 0 CSR-TP (CSR recovers 0 reference TPs on GGO_19 — purely net-negative here).
  **FIRST realized over-enum reduction this session** (contrast: retained-intron was an unrealizable
  ceiling). Recommend default-flip pending approval. Tools `bench/csr_classify.py` (f7f386c).
- **2A (donor-snap) SHELVED — prize = 0** (`bench/donor_snap_prize.py`, 3eb6cb1). Of 158 RU-only FPs,
  NONE has a weak donor within 25bp of a stronger one using a NON-reference intron. The donors the trace
  flagged as snappable noise (e.g. 30800990 @ RSTL.398) are REAL reference donors; the FP isoforms
  (RSTL.398.5/.16) are wrong COMBINATIONS of real introns (chimeric chains), not spurious-donor
  artifacts. Snapping would destroy 29 real introns (25 isoform-TPs) for 0 prize. ABORT (cost≥prize, and
  net still negative under the narrowed canonical/guide gate). 4th oracle correction of a trace
  mechanism claim.

**Lever-#2 net: +3 FP realized (2B opt-in); 2A shelved.** Meta-insight reinforced across all over-enum
work: RU's FPs are wrong COMBINATIONS of REAL introns (chimeric path enumeration), not spurious building
blocks — every building-block-level fix (retained-intron predicate, donor-snap) finds 0 realizable prize
because the introns/donors are all real. The ONE realized win (2B) worked precisely because it suppressed
an EXTRA RU-only mechanism (the chimeric-suffix rescue) ST lacks — a clean behavioral difference, not a
population-matching problem. Realizable parity gains live where RU does something EXTRA vs ST (suppress
it), NOT where RU must reproduce ST's population-dependent path selection (= the flow floor).

---

## §6j — Layer 4 store-gate divergence (2026-06-09)

Grounds the next task. `RUSTLE_ST_SHADOW=1` shadow mode. Branch `vg/flow-capacity-apportionment`.
NO Rust changed in this task — characterization only.

### (a) Current gate baseline (HONEST gate, final emitted GTFs)
`python3 bench/gtf_chain_diff.py /tmp/ru_final.gtf /tmp/st_final.gtf`:
- Rustle: 2992 tx (2991 multi-intron chains, 1 single-exon)
- ST: 1845 tx (1821 multi-intron chains, 24 single-exon)
- multi-intron in both: 1497 — **Rustle-only: 1494** — **ST-only: 324**

(11-day-old memory baseline was ~1506/~322; the layer baseline shifted slightly to **1494 / 324**.
Use 1494 rustle-only as this layer's working baseline.)

### (b) Flow vs filter split
For each of the 1494 rustle-only final chains, checked whether its intron chain appears in ST's
`path_extracted` parity log (`/tmp/st.jsonl`). Coordinate convention: a GTF intron tuple
`(exon_end, next_exon_start)` maps to a `path_extracted` `introns` segment `(exon_end+1, next_exon_start-1)`
(offset `+1/-1`, verified empirically — 200/200 sample chains matched at exactly that offset).

- **FLOW divergence (chain ABSENT from ST `path_extracted`): 1379 (92.3%)** — ST's flow never even
  extracted a path with this chain. Matches the prior "~92% flow" estimate.
- **FILTER divergence (chain PRESENT in ST `path_extracted`): 115 (7.7%)** — of those, ST `pred_kill`'d
  49 (its chains appear in a `pred_kill` payload), 66 present-but-not-killed (likely collapsed into a
  different ST survivor or coordinate-edge mismatches; not this task's target).

This task targets the **1379 FLOW-divergence chains** — the surviving sibling seed-paths.

### Trace cases (flow-divergence)
All three confirmed absent from ST `path_extracted`. Signature: rustle stores **sibling seed-paths**
each carrying minimal **independent** mass (`flux≈1.0, longcov=1.0, entry_abund=1.0`):
- **Trace 1 (cleanest):** `NC_073243.2 (-) 29234259-29293356`, 1 intron `29236027-29292449`.
  Rustle `path_extracted`: `source=flow cov=1.0000 longcov=1.0 entry_abund=1.0 flux=1.0000 raw_flow=2.0 seed_tf=0 nexons=2`.
  ST: **zero `path_extracted` records anywhere in 29234000-29294000** — ST never extracts this path at all.
- **Trace 2:** `NC_073243.2 (-) 111697442-111734611`, 10-intron sibling (`seed_tf=29`), one of ~17 sibling
  seeds at this locus; the high-abundance seeds (`seed_tf=1, flux=36`) ARE shared with ST, but the thin
  siblings (`flux=1.0, raw_flow=16`) are rustle-only.
- **Trace 3:** `NC_073243.2 (+) 22704221-22707567`, 5-intron sibling (`seed_tf=6, flux=1.0, raw_flow=8.0`)
  living beside a `flux=3.0` primary (`seed_tf=1`) that ST keeps.

### (c) Rustle store site
`src/rustle/path_extract.rs` — inside `parse_trflong`'s seed loop. The store cascade:
1. `path_extract.rs:8879` — `if long_read_mode && !mixed_mode && !guide && flow_flux <= 0.0` → defer to
   checktrf (zero-flux), `continue`.
2. `path_extract.rs:8991` — `if coverage < min_cov_gate && !is_guide_pred` → defer to checktrf
   (`min_cov_gate` default **0.15**, env `RUSTLE_MIN_COV_GATE`; ST-parity-lite, NOT in ST), `continue`.
3. `path_extract.rs:9064` — `if config.eonly && !is_guide_pred` → skip.
4. **`path_extract.rs:9180` — `out.push(Transcript{…})` = the STORE.**

**Operative rustle store predicate:** store iff `flow_flux > 0.0` (8879) **AND** `coverage >= 0.15`
(8991) AND `!eonly` (9064). `coverage` and `flow_flux` are computed in the block at
`path_extract.rs:8377-8859`: `flow_flux` = the `flux` (`lf`) returned by `long_max_flow_st` /
`long_max_flow_seeded_with_used_pathpat` (long-read branch, `path_extract.rs:8539-8557`); `coverage` =
`Σ_j nodeflux[j]·noderate[path[j]]` summed over path nodes (`path_extract.rs:8741-8815`), mirroring ST's
`ecov`. The flow-divergence siblings have `flow_flux≈1.0, coverage≈1.0` and pass trivially.

### (d) StringTie store condition
`tools/stringtie/rlink.cpp`, `parse_trflong` (defined line 9758):
- `rlink.cpp:9856` — `flux = long_max_flow(...)` computes the path's max flow.
- `rlink.cpp:9872` — `if(flux || transfrag[t]->guide)` — **GATE A:** only non-zero-flux (or guide)
  paths proceed; else `tocheck` stays true → `checktrf.Add(t)` at `rlink.cpp:9044`. This is reached only
  after `back_to_source_fast_long` (9850) AND `fwd_to_sink_fast_long` (9854) both succeed.
- `rlink.cpp:9918-9950` — build exons; `cov += nodeflux[j]*noderate[path[j]]` (9923/9946).
- **`rlink.cpp:9982` — `if(g || (!eonly && len>=mintranscriptlen && cov>epsilon))` → `pred.Add(p)`
  (`rlink.cpp:10020`) = ST's STORE.** GATE B.

**`long_max_flow` capacity source:** `rlink.cpp:8614/8621` builds network capacity from the CURRENT
`transfrag[t]->abundance`, which is DEPLETED by every prior seed's extraction (`nodecov[path[j]] -=
nodeflux[j]` at `rlink.cpp:9926`, plus per-transfrag abundance consumption). So for a sibling seed
processed after its mass was consumed, `long_max_flow` returns `flux ≈ 0`.

### (e) The named divergence
**The store PREDICATE TEXT is essentially identical on both sides** (`flux>0 ∧ cov>ε ∧ len≥mintranscriptlen`,
rustle additionally requiring `coverage≥0.15`). **The divergence is NOT the predicate — it is the VALUE
of `flow_flux`/`coverage` that rustle's flow extraction returns for sibling seed-paths.**

- ST: for a thin sibling whose mass an earlier seed already depleted, `long_max_flow` (rlink.cpp:9856)
  returns **`flux = 0`** → ST fails GATE A (`rlink.cpp:9872`) → the path **never reaches the store at
  `rlink.cpp:9982`**; it is routed to `checktrf` and (in `!mixedMode` long mode) almost never stored
  there. Net: ST emits NO `path_extracted` for the chain → it is rustle-only.
- Rustle: for the same sibling, `long_max_flow_st` / `long_max_flow_seeded_with_used_pathpat`
  (`path_extract.rs:8539-8557`) returns **`flow_flux ≈ 1.0`** with **`coverage ≈ 1.0`** (INDEPENDENT,
  non-zero per-seed mass), so it passes 8879 + 8991 and is STORED at `path_extract.rs:9180`.

So ST's effective store gate carries an **implicit depletion-coupling requirement** that rustle's does
not realize: ST only stores a seed-path if `long_max_flow` over the **already-depleted** graph still
carries non-zero flow — i.e. the path must own residual, NOT-yet-claimed abundance mass. Rustle's
extraction preserves independent positive per-seed flux for the siblings (the prior depletion-port was
gate-neutral precisely because the depletion ran but the siblings' flux survived it). 1379 of 1494
rustle-only chains are these surviving independent-abundance siblings.

**Rustle-side inputs the next task will gate on** (the variables that differ from ST's `flux≈0`):
- `flow_flux` (`path_extract.rs:8859`, the `flux`/`lf` from the long-read max-flow call).
- `coverage` (`path_extract.rs:8378`, = `Σ nodeflux·noderate`).
- supporting: `raw_flow_sum_out` (`path_extract.rs:8377/9153`, logged as `raw_flow`), and the per-seed
  entry abundance `seed_entry_abund[idx]` (logged as `entry_abund`).
The next task must reproduce ST's depletion-coupled rejection: a sibling should be stored only if its
flux reflects residual UNCLAIMED mass after prior seeds, mirroring how ST's `long_max_flow` returns 0
once the supporting transfrag abundance has been consumed — NOT merely `flow_flux>0 ∧ coverage≥0.15`.

### Caveats / honesty
- The store SITES are pinned exactly on both sides (`path_extract.rs:9180` / `rlink.cpp:9982`, plus the
  upstream flux gate `rlink.cpp:9872`).
- I did NOT instrument ST's `long_max_flow` to print `flux=0` for the specific Trace-1 transfrag (ST's
  per-seed flux is not in the parity log); the "ST returns flux≈0 for depleted siblings" conclusion is a
  read of the capacity-from-`abundance` code path (rlink.cpp:8614/8621 + depletion at 9926) plus the
  EMPIRICAL fact that ST emits zero `path_extracted` for Trace 1's region. A confirmatory instrumented ST
  flux print for one trace transfrag is left to the next task if a tighter proof is wanted.

---

## §6j-diag — Layer 4 flux-survives-depletion diagnosis (2026-06-09, Task 2 Step 1)

**Verdict for Task 2: BLOCKED.** The premise that the rustle-only siblings survive because the
per-seed flux reads an **un-depleted capacity structure** is REFUTED by direct instrumentation. The
capacity IS already coupled to the live, shared, depleted abundance, and depletion propagates
correctly. The 92% FLOW-divergence is NOT a flux-coupling defect — it is an **upstream read→transfrag /
seed-construction divergence**: rustle forms long-read seeds that StringTie never forms at all.

### What carries the capacity (the structure the task asked to identify)
`max_flow.rs:1280` — `capacity[i][end_idx] += tf.abundance` (and the source/sink mirror at 1298-1299)
reads the **CURRENT value of the shared `transfrags[t_idx].abundance` slice**, gated by `tf.abundance > 0`
at `max_flow.rs:1189`. That slice is the SAME `&mut [GraphTransfrag]` threaded across every seed in the
loop (`path_extract.rs:6751 for (loop_pos, idx) in order` → `long_max_flow_st` at 8540, no per-seed
clone on the long-read non-guide path). The depletion tail (`max_flow.rs:1834 if !no_subtract`,
`update_transfrag_capacity` at 1930/1961) writes back into that same slice and is **structurally
identical to ST** (rlink.cpp:8692-8731: same `flow[n1][n2]>0 ∧ nodes[0]==path[i]` predicate, same
subtract-to-zero). The rustle-only protections (subseq-protect 1896, alt-splice-protect 1868,
seed-chord 1242) are all correctly **disabled under `st_shadow()`/`exact`** — verified
`RUSTLE_DEBUG_FLOW_SUBSEQ=1` → 0 hits. **So the flux is NOT independent of depletion by construction;
it reads the depleted pool already.**

### Empirical proof that depletion reaches and zeroes siblings (so coupling already works)
Trace-1 bundle `NC_073243.2:29200000-29400000` slice (49-tf bundle), `RUSTLE_ST_SHADOW=1
RUSTLE_DEPLETION_DIAG=1` + a temporary `SEED_CHAIN` print (since removed, `path_extract.rs` clean):
- Dominant seed idx=9 chain `[3,4,5,6,7,8,9,10,22,23,24]` flux=69.8 depletes the shared pool.
- Later overlapping seeds whose support it claimed: idx=1 `[5,6,7,8,9,10,22,23,24]`, idx=3, idx=8 all
  enter the loop at `ab_before=0.0 → flux=0.0`. **Depletion DID propagate through the shared slice.**
- The SURVIVING siblings idx=21 `[1,2,9,10,22,23,24]`, idx=22 `[3,4,9,10,22,23,24]`, idx=26
  `[3,4,9,22,23,24]`, idx=27 `[3,9,10,22,23,24]` each enter at `ab_before=1.0 → flux=1.0`. They are
  **distinct single-read transfrags** (each a different first-exon / skip-junction structure) that the
  dominant seed's flow never routed through, so `update_transfrag_capacity` never touched them. This is
  ST's algorithm exactly — ST seeds these same reads and also gives them independent flux.

### The REAL divergence (why ST emits zero for these chains): seeds ST never forms
Ran ST and rustle with the parity-decisions JSONL (`STRINGTIE_PARITY_LOG` / `RUSTLE_PARITY_LOG`) on the
same slice and diffed seeds vs stores:
- ST: 25 `parse_trflong_seed`, 25 `transfrag_pre_depl`, 25 `transfrag_collapse`, 17 `path_extracted`.
- rustle: 34 `parse_trflong_seed`, 34 `transfrag_pre_depl`, 22 `path_extracted`.
- Of the **9 rustle-only stored chains in the region, 8 have NO matching ST seed at all** — including
  **Trace 1 `29236027-29292449`**, which is rustle seed `f_idx=0 t_idx=0 entry_abund=1.0` but appears
  NOWHERE in ST's log (not seed, not pre_depl, not collapse). ST's read→transfrag stage
  (keeptrf containment-collapse, rlink.cpp:5805 `transfrag_collapse`) folds/drops these reads so they
  never become distinct `trflong` seeds; rustle keeps 34 vs ST's 25.
- The **1 of 9** that IS an ST seed (`29314688-29317778,29318188-29320347`, ST `entry_abund=1`, ST
  seeds it but does NOT store it) is the only genuine flux-coupling case in the region — ST depletes it
  via an earlier overlapping seed → flux=0 → drop; rustle stores it. This is the residual that a
  flux-coupling change could address, but it is the minority, not the 92%.

### Why this is the documented "everything-at-once" wall (DECISION = BLOCKED)
A change to the per-seed flux/store gate cannot reach the dominant 8/9 mechanism, because those chains
have **no ST analog to converge toward** — there is no ST seed/store whose `flux≈0` rustle should
reproduce; the entity simply does not exist on ST's side. Forcing flux=0 on rustle's extra seeds would
mean **suppressing the seeds themselves** (matching ST's keeptrf collapse magnitude / read→transfrag
folding), which is the upstream Layer-2/3 (read→transfrag, collapse) restructuring — not the Layer-4
flow store gate this task scoped. The siblings' flux is already depletion-coupled and faithful; the gap
is that rustle constructs ~36% more long-read seeds than ST. Escalate: the realizable lever is the
read→transfrag collapse layer, not depletion-coupled per-seed flux. (Tooling for the next pass:
`RUSTLE_DEPLETION_DIAG`, `RUSTLE_TF_INV_TSV`, and the `parity-decisions` JSONL seed/collapse diff above.)

---

## §6L — Default rustle-only (187) first-divergence stage distribution

**Date:** 2026-06-09. **Branch:** `vg/flow-capacity-apportionment`.
**Setup:** `rustle GGO_19.bam -L` vs `stringtie GGO_19.bam -L` (no shadow, no guided, no VG).
Both tools run with full parity logging (five stages: `junction_accept`, `graphnode_list`,
`transfrag_collapse`, `parse_trflong_seed`, `path_extracted`). The 187 rustle-only chains
were confirmed with `bench/gtf_chain_diff.py`. Analysis script: `bench/stage_divergence_187.py`.

### Stage distribution

Each of the 187 rustle-only multi-intron chains was classified by the EARLIEST pipeline stage
where rustle has it and StringTie does not.

| Stage        | Count |    % | Meaning |
|-------------|------:|------:|---------|
| `junction`  |     1 |  0.5% | ST didn't accept ≥1 junction; rustle did |
| `graph`     |     0 |  0.0% | — |
| `collapse`  |     0 |  0.0% | — |
| `seed`      |    61 | 32.6% | ST has no transfrag_collapse rep for this chain; rustle seeds it from an extra transfrag |
| `flow`      |    31 | 16.6% | ST seeds this chain but doesn't extract it via path_extracted; rustle does |
| `post_flow` |    94 | 50.3% | ST DOES emit this chain via path_extracted, then its pred_kill/filter stage removes it; rustle keeps it |
| `unknown`   |     0 |  0.0% | — |
| **TOTAL**   | **187** | | |

### Dominant stage: `post_flow` (94 chains, 50.3%)

ST actually extracts these chains in its flow stage (`path_extracted` event is present) but then
removes them before the final GTF. The pred_kill events were not captured in this run (not in the
filter-steps list), but the ST `path_extracted` payload reveals the mechanism: these chains have
**very low `longcov`**. Of the 89 post_flow chains that were also confirmed in `st_seeds` (a
slightly stricter sub-check):

- Median ST `longcov` = **1.0**
- 49/89 (55%) have `longcov ≤ 1`
- 71/89 (80%) have `longcov ≤ 2`
- 18/89 (20%) have `longcov ≥ 3`

ST's pred_kill / isofrac filter eliminates low-longcov transcripts; rustle's equivalent filter
passes them. This is consistent with the earlier (§5b) finding that rustle over-attributes read
mass to long chains (`ru_read_count == ru_abund` >> ST's 1).

**Three concrete examples (coords = path_extracted 1-based intron format, strand):**

1. `strand=-, 4 introns, first: 53304432-53306277,53306384-53307260,...`
   ST: cov=3956.2, longcov=**2**, entry_abund=2. Kept by rustle, killed by ST.

2. `strand=+, 14 introns, first: 110889996-110901426,110901623-110903773,...`
   ST: cov=3948.7, longcov=**1**, entry_abund=1. Kept by rustle, killed by ST.

3. `strand=-, 15 introns, first: 20534118-20534556,20534676-20534829,...`
   ST: cov=6665.3, longcov=**3**, entry_abund=3. Kept by rustle, killed by ST.

### Second stage: `seed` (61 chains, 32.6%)

These chains never appear in ST's `transfrag_collapse` set — ST does not build a transfrag
representative for them at all. Of the 71 seed-classified chains, **0/71 are in ST's collapse
set** (verified independently). This is the over-segmentation problem documented in §5b:
Rustle produces ~36% more long multi-intron transfrags than ST (7335 vs 4383 at
`transfrag_pre_depl`), creating seeds that ST never forms.

**Three concrete examples:**

1. `strand=-, 15 introns, first: 32006522-32009343,32009424-32011744,...`
   Not in ST collapse, not in ST seeds. Source: extra Rustle transfrag.

2. `strand=-, 25 introns, first: 22460560-22461440,22461607-22461696,...`
   Not in ST collapse, not in ST seeds.

3. `strand=-, 5 introns, first: 56745726-56746448,56746680-56747380,...`
   Not in ST collapse, not in ST seeds.

### Third stage: `flow` (31 chains, 16.6%)

ST seeds these chains (present in `parse_trflong_seed`) but flow depletion prevents extraction
(`path_extracted` absent for ST). This is the "flux=0 / depletion-coupled" mechanism described
at §6K — the chain exists as a seed but a competing path with higher flux depletes the shared
edges before the seed can be extracted.

### Single `junction`-stage case (1 chain)

`strand=-, 11 introns, first junction donor=21845195 acceptor=21852482`: ST's `junction_accept`
events do not include this junction (with either strand); rustle accepts it (mm=ok). This is an
isolated case; the junction may be near the mm_negative boundary.

### Accuracy caveats

1. **Graph-node coverage check**: The `graphnode_list` stage check uses a coordinate-overlap
   approximation (exon boundary in any node interval). Since 0 chains were classified as `graph`,
   this approximation had no effect on the final counts.
2. **ST transfrag_collapse truncation**: ST's parity log uses a ~480-char line buffer, truncating
   1081 of 7982 `transfrag_collapse` events mid-`rep_introns`. The analysis uses partial matches
   (dropping the truncated final segment) where possible. This may under-count ST's collapse set,
   but since 0 of 187 chains are classified as `collapse` and 0/71 seed chains matched any ST
   collapse entry, the truncation does not affect the main conclusions.
3. **15 bypass-collapse chains**: 15 rustle-only chains reach `path_extracted` without appearing
   in rustle's `transfrag_collapse` events (likely sub-bundle or checktrf paths). These are
   correctly classified by the seed/flow/post_flow checks at subsequent stages.
4. **`post_flow` pred_kill reason**: The specific ST filter reason (isofrac / predcluster /
   longcov threshold) was not captured in this run. The longcov distribution strongly suggests
   the `pred_intron_low`/isofrac check as the cause, consistent with §5b.

### Actionability

| Stage        | Realizable fix | Risk |
|-------------|---------------|------|
| `post_flow` (50.3%) | ~~longcov-floor~~ → §6m PINS this: match ST's coverage metric (lowintron+pred cov), NOT a filter | High — root is flow-coverage divergence; longcov floor over-kills (ST keeps 439 longcov==1) |
| `seed` (32.6%) | Reduce Rustle transfrag over-segmentation (Layer-2/3 read→transfrag) | High — architectural change; §5b documents prior attempts |
| `flow` (31, 16.6%) | Match ST's flow depletion / path selection | High — same root as over-segmentation |

The highest-leverage **filter-level** fix is at `post_flow` (longcov floor for long-chain
predictions), which requires no structural change to the flow. The `seed` + `flow` population
requires the multi-session read→transfrag rewrite documented in §5b.

---

## §6m — Pinned: which GATE drops the 94 `post_flow` chains (2026-06-09)

Direct trace of the 94 `post_flow` rustle-only chains against ST's instrumented log
(`/tmp/stP.jsonl`: 5199 `path_extracted` + 2100 `pred_kill`). Tooling: `/tmp/pin_gate.py`,
`/tmp/pin_gate2.py`, `/tmp/ri.py`, `/tmp/diffsets.py`, `/tmp/probe64.py`.

### Answer: NOT GATE A (flux) or GATE B (store cov). Downstream PAIRWISE selection.

All 94 chains appear in ST's `path_extracted` (so `long_max_flow` extracted them → **GATE A flux
passes**) and carry a high flux-`cov` (median 4941, min 729 → **GATE B `cov>epsilon` store passes**;
GATE B uses flux-cov, which is high). ST drops them **after** store, at the pairwise
selection / output-filter stage. The discriminating signal is **`longcov` (long-read count),
not flux-cov**: post_flow `longcov` median = 1 (50/94 exactly 1; entry_abund==longcov throughout)
vs ST-kept median = 3 — but **not a clean threshold** (ST keeps 439 chains at longcov==1, so an
absolute longcov floor over-kills).

### Mechanism split of the 94
- **30 chains** → explicit `pred_kill reason=retained_intron`, `stage=pairwise`, cov/killer_cov
  ratio **0.03–0.10** (28/30 ≤0.10). A retained-intron variant of an overlapping higher-cov
  spliced killer, retainer cov <~10% of killer. Named ST filter.
- **64 chains** → not in `pred_kill` (filtered at an uninstrumented stage). Low-cov near-equals:
  14 subset of an ST output chain, 17 one-junction-off, 5 superset, 23 other-overlap, 5 no overlap;
  **54/64 have longcov ≤2**. Same path-selection-on-shared-junctions pattern as §6 altsplice.

### Why rustle can't cleanly reproduce these (THE WALL, precision side)
Running rustle with `RUSTLE_PREDCLUSTER_ST=1` (the faithful ST selection port) on the same input:
rustle-only 187→175 (−12) but ST-only 104→**117 (+13 over-killed TPs)** — net wash. Critically it
kills only **2 of the 30** retained_intron-targeted chains while over-killing 15 shared TPs and
removing 33 *different* rustle-only FPs. **The port's RI/included_drop/isofrac decisions hinge on
ST's per-base `lowintron` mask + pred coverage, which rustle computes divergently** (the documented
flow-coverage divergence, `project_coverage_metrics_deviation`). So the selection *logic* is faithful;
the *coverage input* is not. No filter or threshold reproduces ST's 94 drops without the upstream
coverage metric.

### Conclusion (SUPERSEDED — see §6n for the precise root cause)
Closing the 94 `post_flow` requires **matching ST's flow enumeration**, not a coverage formula fix.
See §6n for the complete trace. The §6L "longcov floor" row is **downgraded**: a longcov floor
over-kills (ST keeps 439 longcov==1 chains).

---

## §6n — Coverage divergence root cause for the 30 RI kills (2026-06-09)

Traced the 30 retained_intron post_flow kills to their exact coverage divergence. Tooling:
`RUSTLE_LOWINTRON_TRACE` (transcript_filter.rs:1663), comparison scripts `/tmp/ru_cov.py`,
`/tmp/diffsets.py`, `/tmp/pin_gate2.py`.

### Finding: coverage divergence is flow-enumeration depletion, NOT a formula bug

For the specific killer case at locus `59398665-59575456` (26-exon, ST killer_cov=37.1):
- **Rustle has the exact same 26-exon chain** (same coords, same nexons) — RSTL.479.11, cov=**3.29**
- **But rustle also emits two extra chains** (27-exon cov=22.08 + another 26-exon cov=18.99)
- Those extra chains were extracted FIRST and depleted the shared node coverage budget
- When rustle extracts the ST-equivalent 26-exon killer, only leftover coverage remains → 3.29 vs 37.1 (11x lower)
- Victim/killer ratio in rustle = 2.98/3.29 = **0.904 >> 0.1** → RI kill doesn't fire

For another case (`31707587-31722117`, 8-exon, ST cov=37.2):
- Rustle has the 8-exon chain (RSTL.100.3 cov=2.98) but also a 9-exon chain (RSTL.100.1 cov=6.22)
  and RSTL.100.5 (8-exon, raw_flow_sum=9.85)
- `RUSTLE_COV_RAW_FLOW` raises raw killer cov to 9.85 vs ST's 37.2 — improvement is 1.6x,
  still 3.8x short of ST; ratio = 1.91/9.85 = 0.194 >> 0.1

**The 14 killer-present cases break down:**
- 2/14 kill correctly (ratios 0.034, 0.094 in rustle — both killer_cov large and close to ST)
- 12/14 fail: rustle's killer is 2–12x under-covered (depletion by extra paths)
- `RUSTLE_COV_RAW_FLOW` does NOT fix this — improvements 1.3–1.6x while divergences are 4–12x
- `RUSTLE_RI_USE_RAW_FLOW` is same (uses raw_flow_sum in the RI ratio check)

### Self-reinforcing loop: the precise mechanism

```
Rustle over-enumerates → extra paths extracted first
  → deplete shared nodecov
  → dominant "killer" chains get under-counted coverage (3.29 vs 37.1)
  → RI kills don't fire (victim/killer ratio too high)
  → victim chains survive → counted in the 94 post_flow FPs
```

The 94 post_flow FPs and the 30 RI coverage failures are the SAME problem: rustle's flow
over-enumeration. The coverage divergence IS the flow divergence. Fixing one requires fixing the other.

### Why predcluster_st doesn't help

`RUSTLE_PREDCLUSTER_ST=1` kills only 2/30 RI victims while over-killing 13 shared TPs — net wash.
It fails because:
1. 16/30 killers aren't in rustle at all (structural chain divergence — different path set)
2. 12/14 killer-present cases: killer too under-covered to trigger the 0.1 ratio (depletion)
3. The lowintron flags DO fire on the killers — `RUSTLE_LOWINTRON_TRACE` shows flagged introns —
   so the flag is not the problem

### Updated actionability for §6L post_flow row

| Stage | Old conclusion | Updated conclusion |
|-------|---------------|-------------------|
| `post_flow` (50.3%) | Match ST coverage metric | Fix flow over-enumeration (parse_trflong). Coverage is a symptom, not the root. |

### Approach A (bpcov killer coverage in RI gate) — TESTED and FALSIFIED (2026-06-09)

Hypothesis: substitute the killer's flow-depleted `coverage` with its undepleted per-base bpcov
average in the RI cov-ratio gate (`n2.coverage < frac * n1_cov`), so the kill fires despite
depletion. Implemented in `predcluster_st::retainedintron_st` + call site, tested under
`RUSTLE_PREDCLUSTER_ST=1` vs `../GGO_19.gtf`:

| Config | Rustle-only | ST-only | both |
|--------|-------------|---------|------|
| default (predcluster_st off) | 159 | 69 | 1745 |
| predcluster_st + bpcov-RI (both branches) | 113 | **188** | 1626 |
| predcluster_st + bpcov-RI (RI branch only) | 113 | **188** | 1626 |

**RESULT: over-killed 119 shared TPs** (ST-only 69→188). The RI branch alone produces the entire
regression. Root cause of the failure: raw bpcov averages ALL reads at the killer's exons,
including sibling-isoform reads at shared exons, so `n1_cov` is inflated far above ST's allocated
value. This makes `n2.coverage < 0.1 * n1_cov` true for legitimate alternative isoforms, not just
RI artifacts — the cov-ratio gate's discrimination is destroyed.

**Decisive conclusion:** NO per-transcript scalar in rustle recovers ST's allocated-undepleted
killer coverage: `tj.coverage` is depleted (3.29), bpcov-avg is sibling-inflated (~37 but
over-kills), `longcov` is a small read count (2-3, wrong units). ST's 37.1 only exists if the
extra sibling paths are never extracted. **The fix MUST be at flow enumeration (parse_trflong
over-enumeration), not the downstream RI filter.** Code reverted to clean (default byte-identical,
159/69 restored). Confirms §6n self-reinforcing-loop diagnosis from the fix side.

---

## §6o — Approach B (checktrf multinode store-gate) — TESTED and FALSIFIED (2026-06-09)

A multi-agent investigation (workflow wf_abb54bf1-34e) of the 187 default rustle-only chains split
them into 139 flow-only + 48 checktrf-only. The checktrf pool looked like a clean structural
divergence: StringTie's `parse_trflong` (rlink.cpp:10369/10413) handles a multi-node long-read
transfrag with no kept-path match by **redistribute-only** — its independent-store `else` branch is
structurally unreachable for such transfrags, so they are NEVER stored. Rustle stored them
(path_extract.rs:9779), producing the 48 checktrf-only chains (measured 2 TP / 46 FP vs RefSeq).

Built the ST-faithful gate (predicate `checktrf_multinode_no_match_drop` + `SeedOutcome` variant +
the gate; TDD + two-stage review). **Validation killed it.** Gate removed 201 chains:

| Config (vs /tmp/stP.gtf) | Rustle-only | ST-only |
|--------------------------|-------------|---------|
| default (pre-fix / opt-out) | 187 | 104 |
| gate ON (default-on, as built) | 157 | **259** |

The gate over-rejected **156 StringTie-shared real isoforms** (up to 37 introns) and cost **18
RefSeq TPs** — not the ~46 FPs it was scoped for. Root cause: rustle's checktrf rescue is
**load-bearing for recall** — it recovers 156 chains StringTie finds via FLOW but rustle's weaker
flow misses. Separability analysis (gate-time transfrag features, `/tmp` analysis):

| Feature | should-DROP (43 genuine FP) | should-KEEP (158 shared/TP) |
|---------|------|------|
| cov median | 1.01 (range 0.43–13.2) | 1.06 (range 0.37–619.8) |
| longcov median | 1.0 (30/43 ≤1) | 1.0 (85/158 ≤1) |
| entry_abund median | 1.0 | 1.0 |
| nexons median | 11 | 9 |

**No gate-time feature separates them** (cov low-ends identical; FPs are if anything MORE complex).
This is the SAME WALL as Approach A (§6n): a downstream gate "ST-faithful in isolation" over-kills
because rustle's UPSTREAM FLOW differs from ST's. My scoping error: measured only the rustle-ONLY
checktrf FPs, missed the 156 shared-chain collateral; the validation caught it.

**Disposition:** flipped to **opt-in, default-off** (`RUSTLE_CHECKTRF_MULTINODE_DROP=1`); default
byte-identical (187/104 restored). The gate is preserved as the **downstream complement** to
flow-enumeration parity: once rustle's flow produces those 156 chains directly, the checktrf
compensation is unnecessary and enabling this gate becomes correct + non-regressing.

**Convergent conclusion (Approaches A + B):** the 187 rustle-only chains are NOT closable by any
downstream gate — both the RI filter (bpcov) and the checktrf store-gate over-kill identically,
because rustle's checktrf/filtering compensates for flow gaps. **The only real lever is
flow-enumeration parity** (`parse_trflong`): the ~2× seed count / keeptrf-consolidation divergence
the workflow surfaced. That is the next (deferred, multi-session) target. Spec/plan:
`docs/superpowers/{specs,plans}/2026-06-09-checktrf-multinode-store-gate*`.

---

## §6p — Canonical flow-port regression characterized (Component 1, 2026-06-09)

Flow-parity Component 1 (spec `2026-06-09-flow-enumeration-parity-port-design.md`): characterize why
`RUSTLE_PARSE_TRFLONG_ST_CANONICAL=1` regresses (223 rustle-only / 144 ST-only, vs default 187/104).
Tooling: `bench/canonical_divergence.py`, path-level trace from the path_extracted logs.

### The regression is bidirectional
3-way chain-set diff (default / canonical / ST vs `../GGO_genomic.gff` NC_073243.2):
- canonical **ADDS 87** chains vs default: 16 ST-shared (converging) / **69 canonical-only-FP
  (regressing)** / 2 canonical-only-TP.
- canonical **REMOVES 91**: **56 ST-shared (recall lost → ST-only 104→144)** / 9 RefSeq-TP.
- canonical extracts **5883 path_extracted events** vs ST 5199 / default 5098 (same 10526 seeds) —
  it OVER-extracts paths.

### Root mechanism: long_max_flow_st under-depletes (ONE shared bug, not 49 causes)
The 71 canonical-only additions cluster across 49 loci (38 singletons, 11 multi). Path-level trace
at the top 3 loci (44156486-, 19115694-, 22459860-) shows the SAME mechanism at each: canonical
extracts MORE paths than BOTH ST and default, and the extras are **low-cov (<1.0) near-duplicate
sub-paths of the dominant backbone, clustered just at/below cov 1.0** (e.g. locus 22459860:
canonical 27 paths / 19 low-cov vs default 15/6 vs ST 11). ST and default both deplete these residual
sub-paths away; `long_max_flow_st` leaves ~1.0 residual flux and keeps extracting them (they pass the
min_cov_gate 0.15 store gate). So canonical is LESS parsimonious than both ST and default — a
depletion/stopping defect in `long_max_flow_st`, not a pervasive multi-cause divergence.

### VERDICT: close-but-broken → DEBUG the existing port
The 49-loci scatter is one shared depletion bug, so a single fix in `long_max_flow_st`
(`src/rustle/parse_trflong_st.rs:991`) — correct the abundance/nodecov depletion so dominant-backbone
sub-paths drop below the flux/store threshold like ST's `long_max_flow` (rlink.cpp:9856, 9926/9939
nodecov depletion) — should converge most of the 49 loci. Component 2 = trace + fix that depletion,
canonical-gated, measured at `bench/gtf_chain_diff.py /tmp/ru_canon.gtf /tmp/stP.gtf` (223 → toward
187 then below). Worklist (top regression loci, all strand `-` unless noted): 44156486-44307898 (6),
19115694-19174217 (4), 22459860-22469708 (4), 70716579-70761097 + (3), 27242275-27289018 (3),
97518383-97540599 (3), then the 11 multi-chain loci before the 38 singletons.

⚠ Node-flux instrumentation gap: `RUSTLE_COV_DEBUG` does not fire in canonical mode (canonical routes
through `long_max_flow_st`, which lacks that trace). Component 2 should add a COV_DEBUG-equivalent
inside `long_max_flow_st` to trace per-node depletion directly. Path-level tracing (above) sufficed
for the Component-1 verdict.

---

## 7. Superseded documents

The following are superseded by this file for the precision/parity-gap analysis (kept for history):
`EXACT_STRINGTIE_PARITY_PLAN.md` (its "100% parity via guided mode" premise is circular — `-G GGO_19.gtf`
scored against `GGO_19.gtf`), `PRECISION_DIAGNOSTIC_FINDINGS.md`, `PRECISION_IMPROVEMENT_PLAN.md`,
`PRECISION_WORK_CONTINUATION.md`, `docs/STRINGTIE_PARITY_SYSTEMATIC.md`,
`docs/STRINGTIE_LONGREAD_PARITY_CONTINUATION_GUIDE.md`, `docs/SUBBUNDLE_PARITY_HANDOFF.md`.
