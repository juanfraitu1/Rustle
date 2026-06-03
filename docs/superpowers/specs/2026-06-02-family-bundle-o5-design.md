# O5: family-guided "richer bundle" — sharing structural evidence across paralog copies

**Date:** 2026-06-02
**Status:** ⚠️ PREMISE CORRECTED 2026-06-02 — see banner below. Design retained for the *labelling/honesty* parts; the *structure-borrowing* parts are ALREADY IMPLEMENTED and default-on, so this design must NOT be implemented as written.
**Scope:** VG / `--vg`-gated only. Default de-novo (95.6/90.5) must stay byte-identical.
**Depends on:** [`2026-06-01-flow-capacity-apportionment-design.md`](2026-06-01-flow-capacity-apportionment-design.md) — O5 borrows STRUCTURE; that EM supplies ABUNDANCE. The two are complementary and O5 is unsafe without it (see *Over-enumeration*).

---

> ## ⚠️ CORRECTION (2026-06-02, after grounding against the code)
>
> This design assumed O5 structure-borrowing was unbuilt. **It is already implemented and default-on** in `--vg`:
> `vg::build_bundle_completion_nodes` (dark sibling exon spans → bundlenodes; gate `RUSTLE_VG_COMPLETION_OFF` unset = ON),
> `vg::build_bundle_borrow_junctions` (sibling junctions), `vg::build_bundle_borrow_coverage` (borrowed coverage),
> consumed at pipeline.rs:13089/13122/13194. The earlier "richer bundle = pool all alignments" prototype showed value only because it was a plain `rustle -L` run with **no `--vg` and no genome** — i.e. the existing machinery was never active in the comparison.
>
> **Empirical ground truth (re-run `--vg --genome-fasta <chrY>`):** RBMY c6 (2 genuine primary reads) still assembles to **0 transcripts** *with the completion machinery on*; c1–c5 give 16. When c6's structure was supplied as a `-G` guide, it emitted at **cov 0.000**. So the real O5 gap is on the **coverage/flow** side (structure is present; c6's reads don't produce flow through it), **not** the structure side this design targets.
>
> **Do not implement the structure-borrowing components below.** The still-valid, still-unbuilt parts are: the **honesty labelling** (`family_guided` attribute + low `capacity_confidence` on structurally-inferred transcripts) and the **no-fabrication eligibility gate**. The next action is the diagnostic in `## Re-grounded next step`, not this design's components.

## The O5 promise and the gap

`docs/THESIS_OBJECTIVES.md` O5 — *"share evidence across copies via the family graph"* — is the one objective still partial. O1 (define family), O2 (recover missed copies), O3 (EM assignment), O4 (per-copy isoforms) all land; O5's missing piece is concretely:

> Assembly is **per-bundle**. A paralog copy that is read-starved or whose reads are 5′-truncated assembles to **nothing or a fragment**, even though its sibling copies — which are 93–99.97% identical — show the full gene structure. The `FamilyGraph` *knows* the shared structure but never feeds it back into the per-copy assembly.

This is not hypothetical. **RBMY1 copy 6** (LOC129530242, `NC_073248.2:19,717,578–19,730,926`) is a real testis-expressed copy:

| | c6 default | c6 with family evidence (prototype) |
|---|---|---|
| reads | 2 (both 5′-truncated) | 2 own + family structure |
| transcripts | **0** (no source→sink path) | **full isoform**, spanning the gene |
| structure | — | matches siblings: 22/64 introns within 10 bp of a c5 intron; first introns `(579,1)`c5 ↔ `(579,2)`c6 — the expected per-copy slop, **not noise** |

So the evidence to assemble c6 *exists in the family*; it is simply not shared into c6's bundle. O5 = share it.

## Prototype (what proved the value, and its flaw)

The MVP — physically pool the family's 68 alignments at c6 (clear `SECONDARY`/`SUPPLEMENTARY` flags), assemble:

```
samtools view -H tspy.bam > c6pool.sam
samtools view tspy.bam NC_073248.2:19717578-19730926 \
  | awk 'BEGIN{OFS="\t"}{f=$2; if(and(f,256))f-=256; if(and(f,2048))f-=2048; $2=f; print}' >> c6pool.sam
samtools sort -o c6pool.bam c6pool.sam && samtools index c6pool.bam
rustle -L c6pool.bam -o c6pool.gtf
```

Result: c6 went 0 → a **full isoform set** (transcripts at 19,717,580–19,731,181). **Structure recovered, value proven.**

**The flaw — and it is the whole design constraint:** naive pooling *over-counts and over-enumerates*. c6 came out with 6 isoforms at coverage 11/10/8/4/3/2, but c6's own 2 reads support **one** isoform at coverage ~2. The pooled coverage is sibling reads echoing — the exact chimeric-isoform / double-count failure mode catalogued in `bench/tandem_duplicate_experiment.md` and the j-FP analyses. **Pooling reads for abundance is wrong.** Pooling structure while keeping abundance from the copy's own apportioned reads is right.

## Chosen design — structure-borrow, abundance-own

Two separable channels, mapped to two existing subsystems:

- **STRUCTURE (topology) ← `FamilyGraph`.** The family graph already holds the union of exon classes with **per-copy genome coordinates** (`ExonClass.per_copy_spans`, family_graph.rs:28) and junction edges with `family_support` (family_graph.rs). For a starved copy `c`, walk the family graph's `c`-path (the `path_for_copy` walk, family_graph.rs:80) to get **where c's homologous exons/junctions sit in c's own coordinates** — even exons c has zero reads on. This is a family-derived structural **guide**, analogous to a `-G` annotation but synthesized from the family rather than supplied externally.
- **ABUNDANCE (flow) ← the copy's own apportioned reads.** Coverage/capacity comes only from reads apportioned to `c` by the flow-capacity EM (`read.weight` after `run_fingerprint_em`, vg.rs:5060). A sibling read decisively assigned elsewhere contributes its **junction to the graph** (so the path exists) but **~0 weight to c's flow** (so c's abundance and isoform set stay honest).

The guide injects **edges (topology)**; the flow runs on **apportioned coverage**. That separation is what makes O5 safe: borrowed structure that c's own reads don't support carries ~0 flow and is filtered downstream — over-enumeration cannot survive the apportioned-coverage gate.

### Two regimes (the design only acts on the first)

1. **Starved / truncated copy** (c6: 2 reads, both 5′-short): the copy's reads cannot reach source→sink. The family guide **completes** the transfrag to full length. The 5′ extension is a **homology assertion** ("c6's transcript has the family's 5′ exons; we observed only its 3′ end"), *not* a coverage claim — and it is **labelled as such** (`family_guided "true"`, plus the `capacity_confidence` from the flow-capacity spec, which is low precisely because the 5′ exons have no anchored c6 coverage). Abundance = c6's 2 reads (~2), not 11.
2. **Already-resolved copy** (c1–c5): assembles fine without help. The guide is a no-op here — its edges duplicate edges the copy's own reads already provide, and apportioned flow is unchanged. **Untouched.**

### Honesty contract (non-negotiable, per project memory)

c6's **expression** (~2 reads) is *observed*; its **full-length structure** is *inferred from family homology*. The GTF must make that distinction legible:
- `family_guided "true"` on any transcript whose structure used borrowed exons/junctions the copy's own reads don't cover.
- `capacity_confidence` (flow-capacity spec) low on those transcripts (anchored coverage ≈ 0 on the borrowed exons).
- We do **not** report borrowed structure as independently observed, and we do **not** inflate abundance to the pooled coverage. This is the same discipline as the DAZ3 non-fabrication rule.

## Components

1. **Family structural guide.** New `vg.rs::family_structure_guide(fg: &FamilyGraph, cid: CopyId) -> Vec<(u64,u64)> /* exon spans in cid coords */` — wraps the existing `FamilyGraph::path_for_copy(cid)` walk (family_graph.rs:80) and reads `per_copy_spans` for each node on c's path. Junctions come from `JunctionEdge`s whose endpoints are both on c's path (`family_support`-ordered). No new graph build (the family graph already exists post-O1).
2. **Guide injection into the bundle.** In VG mode, before per-bundle assembly, for each family copy whose bundle is starved/truncated, inject the guide spans as **structural transfrags at zero own-coverage** into that bundle (mirrors how a `-G` reference guide is injected; reuse the guide-transfrag path rather than synthesising reads). Gate so the guide only *adds* edges — never removes or reweights the copy's own reads.
3. **Don't-prune-borrowed-junctions.** Guide-derived edges must survive graph pruning even at ~0 coverage (they carry topology, not abundance). Add a `from_family_guide: bool` marker on injected edges, exempting them from the low-coverage junction filter — and only them.
4. **Apportioned flow (reuse).** No new abundance code: the flow-capacity EM already sets `read.weight` per copy. Borrowed structure inherits ~0 flow automatically because no own-reads cover it. The starved copy's transcript abundance = Σ apportioned own-read weight (~2 for c6).
5. **Labelling.** `Transcript.family_guided: Option<bool>` (beside `capacity_confidence`); set true when the emitted path traverses ≥1 guide-only exon/junction. gtf.rs emit `family_guided "true";` when `Some(true)`.
6. **Starvation trigger.** A copy is "guide-eligible" iff (own primary reads ≥ 1) AND (default assembly yields 0 transcripts OR a transcript covering < X% of the family-consensus length). Threshold env `RUSTLE_VG_FAMILY_BUNDLE_MINFRAC` (default 0.6). A copy with **zero** own reads is **not** eligible (no observed expression → nothing to complete; that would be fabrication, the DAZ3 trap).

## Algorithm

**Phase A — eligibility (post-EM, pre-assembly).** For each family copy `c` with the family graph available: run the normal per-bundle path extraction; if c emits 0 transcripts or only a sub-`MINFRAC` fragment **and** c has ≥1 own primary read, mark c guide-eligible.

**Phase B — guide build.** `guide = family_structure_guide(fg, c)`: c's homologous exon spans (own coords) + family-supported junctions on c's path.

**Phase C — injection + re-assembly.** Inject `guide` as zero-own-coverage structural transfrags into c's bundle; mark injected edges `from_family_guide`. Re-run path extraction for c only. The copy's own apportioned reads now flow along a graph that reaches source→sink; truncated transfrags extend along guide edges.

**Phase D — abundance + labelling.** Flow runs on apportioned own-coverage (unchanged engine). Any emitted transcript traversing a guide-only exon/junction gets `family_guided = Some(true)` and (via flow-capacity spec) a low `capacity_confidence` on the borrowed segment. Abundance = apportioned own-read mass.

**Phase E — emit.** GTF carries `family_guided` + `capacity_confidence`. Abstain/low-confidence tagging from the flow-capacity spec applies unchanged.

## Why over-enumeration cannot survive (the key correctness argument)

The prototype's 6-isoform inflation came from **sibling reads contributing flow** to c6. In this design sibling reads contribute **topology only** (the guide is built from the family graph, not from injected sibling reads; and even if a sibling read overlaps c6, the EM apportions it ~0 weight at c6). Therefore:

- A family isoform c6 does **not** express → its exons get 0 apportioned c6 coverage → 0 flux through `max_flow` (BFS gates on residual, max_flow.rs:155; 0-capacity edge carries 0 flux) → filtered by the coverage gate. **It does not emit.**
- The one isoform c6 *does* express → its 3′ exons carry c6's 2-read apportioned flow; its 5′ exons are guide-completed at low/zero anchored coverage → emits **once**, at abundance ~2, flagged `family_guided`.

So the apportioned-flow channel is the over-enumeration guard. This is precisely why O5 is **gated together with** the flow-capacity EM and is unsafe without it.

## Gating

- New env `RUSTLE_VG_FAMILY_BUNDLE` (**default OFF** — opt-in prototype; promote to default-ON-in-VG only after the validation gates below pass and a benchmark shows net benefit, per the user's instruction).
- `RUSTLE_VG_FAMILY_BUNDLE_MINFRAC` (default 0.6), eligibility threshold.
- Entirely inside `config.vg_mode`; `Transcript.family_guided` written only in VG mode, emitted only when `Some`. **Default de-novo path never executes any of this** → 95.6/90.5 byte-identical (regression gate #1).
- Requires the flow-capacity EM (`RUSTLE_VG_ANCHOR_PRIOR`) active; if that flag is off, `FAMILY_BUNDLE` refuses to run (logs a warning) rather than risk pooled over-count.

## Test & validation plan

1. **Regression guard (first):** default de-novo on GGO_19, no `--vg` → 95.6/90.5, GTF byte-identical to HEAD. Assert no `family_guided` attr ever appears without `--vg`.
2. **RBMY c6 (primary success):** `--vg` + `RUSTLE_VG_FAMILY_BUNDLE=1` on the RBMY array → c6 emits **1** full-length isoform (not 6) at coverage ~2 (not 11), tagged `family_guided "true"` with low `capacity_confidence`; c1–c5 **byte-identical** to the no-flag VG run (guide is a no-op on resolved copies).
3. **Over-enumeration guard (the flaw test):** assert c6's emitted isoform **count == its own-read-supported count** (1), not the family-union count (6). This is the test the naive prototype fails and the design must pass.
4. **No-fabrication guard:** a synthetic family copy with **zero** reads stays empty (not guide-eligible) — guides complete *observed* copies, never manufacture silent ones. (DAZ3 stays silent under `--vg`.)
5. **Tandem synthetic (mamba `sim` env):** the `tandem_duplicate_experiment` synthetic — two near-identical copies, one starved — recovers the starved copy's structure at correct (low) abundance, with no chimeric junctions across copies.
6. **Genome-wide scan:** re-run `bench/paralog_secondary_scan`; the 89 `pure_spillover` candidates must **not** gain spurious guided transcripts (they have no decisive own reads → ineligible / abstained).

## Risks

1. **Guide injection bleeds into resolved copies** — a guide edge duplicating an own-read edge could perturb c1–c5 flux. Mitigate: guide is no-op when the copy already reaches source→sink (eligibility gate); test #2 asserts byte-identical c1–c5.
2. **Borrowed-junction pruning** — if the low-coverage junction filter drops guide edges, structure-completion fails (the original O5 family-guide-projection failure: cov=0). Mitigate: `from_family_guide` exemption (component 3); test #2 must show c6 full-length.
3. **Coordinate projection error** — `per_copy_spans` must give c's *own* coords, not the union span. Family graph already preserves this (family_graph.rs:26-28); assert each guide exon lies within c's genomic extent.
4. **Honesty regression** — a guided transcript emitted without the `family_guided`/low-`capacity_confidence` labels would misrepresent inferred structure as observed. Mitigate: emit is gated on the label being set; test asserts the attr present on any borrowed-exon transcript.
5. **Eligibility threshold tuning** — `MINFRAC` too low → never helps; too high → guides healthy copies. Default 0.6, env-exposed; scan (#6) is the false-positive guard.

## Promotion criteria (flag → permanent)

Make `RUSTLE_VG_FAMILY_BUNDLE` default-ON-in-VG only when **all** hold:
- Regression gate #1 byte-identical (no de-novo impact) — **hard requirement**.
- RBMY c6 emits the correct single full-length isoform at abundance ~2 with honesty labels (test #2/#3).
- No-fabrication + scan guards (#4/#6) show zero spurious guided transcripts on starved/silent copies.
- On a real multi-copy benchmark, family-bundle **adds** recovered real copies (O2/O4 sensitivity) without adding FPs (precision flat or up).

Until then it stays opt-in, documented in `docs/vg_genome_scoping.md` alongside the other VG envs.

## Relationship to the rest of VG

- **Flow-capacity EM** (sibling spec) is a hard dependency — it owns ABUNDANCE; O5 owns STRUCTURE. Neither is complete alone: EM keeps a starved copy's abundance honest but can't assemble it; family-bundle assembles it but would over-count without EM.
- **FamilyGraph** (O1) is the structure source — no new graph machinery.
- **`--read-chain`** is the orthogonal lever for surfacing a copy from raw junction-chains; family-bundle is the homology-completion lever. They can compose (read-chain to surface, family-bundle to complete), tested separately first.

## Re-grounded next step (2026-06-02)

The valuable, *unbuilt* question is narrow and empirical: **why does the existing default-on completion/borrow machinery fail to recover RBMY c6 (2 primary reads), and where in the path does c6 fall out?** Candidate failure points to instrument, in order:

1. **Family membership.** Does c6's 2-read bundle even enter a discovered family (`vg_families` after the quality-filter at pipeline.rs:10267)? If the quality filter drops the low-shared-read c6 copy, it never gets a `FamilyGraph` → no completion nodes built for it. *Check:* log whether c6's `bundle_idx` appears in `em_hmm_partitions` / `family_graphs`.
2. **Completion-node population.** Does `build_bundle_completion_nodes` produce spans for c6's `bundle_idx`? If c6 isn't a key in `bundle_completion_nodes`, the structure injection at pipeline.rs:13122 is a no-op for it. *Check:* log `bundle_completion_nodes.get(&c6_bundle_idx)`.
3. **Borrowed coverage vs gates.** If completion nodes ARE injected, does borrowed coverage (pipeline.rs:13194, with the `RUSTLE_VG_BORROW_FLOOR`/`_LEGACY` gates) put enough flow on c6's exons to seed a transcript, or is it below `readthr`/seed thresholds? The `o5.gtf` cov-0.000 guide result strongly implicates this: structure present, flow absent. *Check:* the per-exon borrowed coverage on c6 and the seed/readthr gate it must clear.

Only after locating the exact fall-out point should a fix be scoped. The likely shape (hypothesis, unconfirmed): c6 is excluded from the family at step 1 (its 2 reads share too little with siblings to pass the quality filter), so none of the structure machinery ever targets it — in which case the fix is in **family membership for low-shared-read copies**, not in new structure-borrowing. The `family_guided` labelling (this doc's one genuinely-unbuilt, still-valid component) attaches downstream of whichever fix lands.
