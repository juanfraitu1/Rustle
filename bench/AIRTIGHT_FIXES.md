# Adversarial review #2 — "make everything more airtight" (2026-06-28)

Second per-objective adversarial pass (after the defense-readiness audit). Each finding below was
**confirmed against the genome / by machine-check**, then fixed. SEDEF (O1 external precision) is excluded —
it runs on the cluster.

| # | Sev | Objective | Finding | Fix | Verified by |
|---|-----|-----------|---------|-----|-------------|
| **M1** | must | O3 ASJ | The flagship splice-mechanism claim — PSMD2/DAXX SNPs "on the canonical GT-AG dinucleotide", "creates/destroys the splice motif" — is **genome-FALSE by one base**. | Retracted; reframed as **splice-REGION (extended-consensus) variants**; anchors at donor−1 / exon boundary; **0/475 on a core dinucleotide**, GT-AG intact. | `bench/asj_motif_check.py` (re-derives from `GGO.fasta`, green) |
| **M2** | must | O2/O3 consistency | Stale strings: `copy_assign.py:487` "silver-standard accuracy"; `OBJECTIVES_STATUS.md` "120 genetic … thesis headline"; `DEFENSE_READINESS_AUDIT.md` "masquerade separator never run / airtight ≈ 20". | One-line corrections: silver = CIRCULAR; genetic core **~77**; O3 row marked **RESOLVED (P4)** + **CORRECTED (M1)**. | grep sweep clean |
| **H2** | high | THEORY | Theorem 4(ii) soundness silently assumed `origin(r) ∈ C`; false in the O4 reference-absent regime (a partial read can be confidently MISassigned to a wrong copy). | Made the **completeness precondition explicit** in statement + proof + scope note; added machine-checks **B6** (precondition necessity) and **B7** (recombinant-cover orthogonality). | `bench/bridge_theorem_check.py` (B1–B7 green); integrated into `copy_assignment_theory_checks.py` |
| **H3** | high | O2 | "Accuracy" is circular: the gate uses the same PSVs that defined the copies; silver = agreement with minimap2 primary. | Added **held-out-PSV cross-validation** (`copy_assign.py crossval`): assign on a disjoint TRAIN half of PSVs, confirm on the TEST half — no ground truth. | sim5x: **80% held-out confirmation = 1.6× / 3.2× / 6.4× over 1/K chance** (K=2/4/8) |
| **H4** | high | O4 | Genome-wide O4 re-indexes the 3 Gb genome on every candidate (minutes/call); no real bounded run had been done. | Added `RUSTLE_ABSENT_MMI` env-gate (pre-built splice `.mmi` target; byte-identical when unset) + bounded real `--absent-copies` run over the 82 principled-catalog regions. | compiles; bounded run → see below |

## M1 — the genome facts (0-based, as stored in `bench/asj_calls.tsv`)
- **PSMD2 (+)**: anchor `195406803` = **donor−1**; canonical donor `GT` at 0-based `[195406804,195406805]` is
  **intact** under both alleles (anchor is a different base). Linkage real: allele G → 14/14 (PSI 1.0), T → 0/18 (PSI 0).
- **DAXX (−)**: anchor `51042120` = the a-exclusive exon boundary; the minus-strand donor `GT` (forward `AC` at
  0-based `[51042118,51042119]`) is **intact**. Linkage: allele C → 9/9, A → 0/10.
- **Genome-wide: 0/475** called anchors fall on a core splice dinucleotide. The load-bearing O3 result is the
  per-molecule **allele→junction linkage**, NOT a dinucleotide-disruption mechanism.

## H2 — the precondition, made load-bearing
Theorem 4(ii) now reads: *under `origin(r) ∈ C`*, `δ(r) ≥ 1 ⟹` the unique consistent copy is the origin ⟹
correct assignment. **B6** drops the precondition and finds **2,616** instances where the gate confidently
misassigns a read whose origin is absent — and certifies the escape (a *full*-column read of an absent origin
is consistent with no copy in `C`, flagged novel; only *partial* reads are at risk). This is exactly the O4
two-stage-freeze rationale: abstain at stage 1, re-thread against `C ∪ {absent copies}`. **B7** exhibits a read
set with several distinct minimum covers of size ≥ 3 (the K≥3 recombination obstruction) on which the per-read
gate stays sound *within each fixed C* — Theorem 4 is per-read-given-C, not cover recovery.

## H3 — held-out-PSV cross-validation (`copy_assign.py crossval`)
```
 K reads_CV heldout_confirm%  chance(1/K)%  enrich  train_acc%(truth)
 2      200            80.0%         50.0%    1.6x             80.0%
 4      200            80.0%         25.0%    3.2x             80.0%
 8      200            80.0%         12.5%    6.4x             79.0%
```
The disjoint held-out PSVs — which played no role in the call — rank the TRAIN-assigned copy first at 80%,
up to **6.4× above chance**, with **no ground truth used**. This is the non-circular signal silver cannot give.
(Absolute rate ~80% not ~100% because each half carries only half the evidence; the enrichment is the point.)

## H4 — `.mmi` pre-index + bounded real run
- `RUSTLE_ABSENT_MMI` → a pre-built splice index (`minimap2 -x splice -d GGO.splice.mmi GGO.fasta`) is used as the
  remap target instead of re-reading the FASTA; **unset ⟹ byte-identical** (target = `fasta_path`).
- **Bounded real run (first real-data O4 pass with the corrected splice-preset gate).** One principled-catalog
  region `NC_073224.2:15670216-15791935` (~120 kb, 2,203 mapped reads), `--absent-copies --skip-poa-diagnostic`,
  `RUSTLE_ABSENT_MMI=GGO.splice.mmi`. Result: **0 reference-absent copies admitted** (no `AC_` in the
  assignments), **25 candidates all routed to DNA-needs** (fails-safe — none force-admitted as a copy). Rejection
  reasons: **15× ≥98% remap identity** (the MAPQ-0 paralog-leak / het the gate is built to catch), 6× host
  sequence unbuildable, 3× not min_p-distinct from host, 1× consensus unplaceable. Silver 188/188 (100%).
- **This empirically confirms the standing "mechanism-only / real-headroom ≈ 0 on GGO" finding on real data** with
  the *corrected* (`-x splice`) admission gate: the gate sees these reads-only clusters, remaps them, finds they
  match the reference at ≥98%, and correctly declines to invent a copy. 0 real reference-absent copies in GGO RNA;
  the architecture is sound and fails safe.
- **At-scale caveat:** ~10 min for this one region because minimap2 reloads the 13.6 GB `.mmi` on every candidate
  remap (25 here). The `.mmi` removes the *re-indexing* cost but not the per-call *load*; a full 82-region /
  genome-wide O4 pass should batch all candidate consensuses into ONE minimap2 invocation (next step) or run on
  the cluster. The env-gate is the enabling first step, not the complete at-scale solution.

## Reproduce
- `MINIFORGE python bench/asj_motif_check.py`            (M1 — genome-pinned, exits non-zero on drift)
- `python3 bench/bridge_theorem_check.py`                (H2 — B1–B7)
- `python3 bench/copy_assignment_theory_checks.py`       (H2 — full theory suite incl. B6/B7)
- `MINIFORGE python bench/copy_assign.py crossval`       (H3 — held-out-PSV CV)
- `RUSTLE_ABSENT_MMI=… target/release/copy_assign --absent-copies …`  (H4)
