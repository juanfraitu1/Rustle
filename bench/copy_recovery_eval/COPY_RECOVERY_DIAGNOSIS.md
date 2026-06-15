# Phase-0 diagnosis: why the PSV super-graph path doesn't recover starved copies

**Date:** 2026-06-14. Flagship: `LOC101129569` family on NC_086018.1 (~18.2–22.5 Mb, ~25 copies).
Missed copy `XM_055375231.2` @ 22178176–22192379 (headroom: mm=25, uniq=0).

## Setup
Slice `NC_086018.1:18000000-23000000` (10,833 reads, 1,152 secondary). `stringtie -L` guide (282 tx).
- baseline `rustle -G`: 284 tx
- PSV path `rustle --vg --vg-layer2 --vg-layer2-psv-linkage -G --genome-fasta`: 318 tx (+34)

## Key finding — copies are MIS-ASSEMBLED, not absent (rules out tree branches a/b/e)
At the flagship locus: **25 primary + 36 secondary** alignments; **6–8 transcripts already land
there** in both base and psv outputs. So:
- The locus is NOT empty and the family/graph DOES exist there → **NOT** (a) family-not-formed or
  (b) copy-not-a-path.
- `XM_055375231.2` is still class `MISS` in both base and psv → the *correct copy isoform is never
  extracted*; the locus is contaminated by sister-spillover reads (uniq=0 → all 25 primaries also
  map to sisters) and the wrong structure wins.
- PSV path emits **0** `vg_copy_recovery`/`PsvLinked`-tagged transcripts; its +34 are ordinary
  `flow`/`guide:` transcripts → the PSV recovery channel produced nothing for this copy.

## Narrowed to branch (c) or (d) — to pin next
- (c) genotyping/PSV columns: does `psv_columns_for_family` produce PSV columns separating this
  copy? does `genotype_family_reads` cast votes for it (using the 36 secondaries)? 
- (d) assignment/assembly/gate: does `assign_read_to_copy` pin reads to this copy? does
  `assemble_psv_isoforms` emit its chain? does `family_identifiability(fg,0.005,psv_family_min)`
  reject the family (E-gate)? is an assembled chain dropped by union-by-chain dedup?

## Reproduce
Slice + runs in `/tmp/diag/` (`s.bam`, `guide.gtf`, `base.gtf`, `psv.gtf`, `ref.gff3`,
`gc_base.*`, `gc_psv.*`). Trace target = `XM_055375231.2`.

## Implication for the design
Recovery is **"PSV-resolve the contaminated reads and extract the correct copy isoform from the
existing graph"** — squarely the super-graph/assignment approach, no clustering. The fix lives in
the genotype→assign→assemble→emit chain and/or the (E) gate, gated on `RUSTLE_VG_RECOVER_COPIES`.

## UPDATE — trace result (subagent) + STRAND-AWARE headroom correction

**Trace finding (on this flagship):** branch (d), failing at the **(E) gate**. `psv_columns_for_family`
returns **0 columns** because the `FamilyGraph` has **0 multi-copy (shared) exon nodes** — every node
is copy-private — so `family_identifiability` is false and the whole PSV pass is skipped (the
genotype/assign/assemble functions never run). Underneath: `XM_055375231.2` is **− strand**, but its
locus holds **25 `+` primary / 0 `−` primary** reads (the `+` sisters' cross-maps); the FamilyGraph
built there is `+` strand and the `−` target is **not even a member copy** (mixed-strand families
bail, family_graph.rs:498). DIAG instrumentation left in `layer2.rs` behind env `RUSTLE_DIAG_TARGET`
(default-off; lines ~550–753).

**=> This flagship is a STRAND-CONFOUNDED DUD** (a `−` copy with no own-strand reads — genuinely
unrecoverable). The strand-blind headroom (`mm=25`) counted wrong-strand sister cross-maps.

**STRAND-AWARE headroom (`/tmp/cre_guided/headroom_strand.py`, own-strand AS-decisive ≥2):**
- **27 GENUINE** own-strand-decisive recoverable copies (was 23 strand-blind), several strong:
  `XM_055374752.2` own=300, `XM_004063578.5` 67, `XM_055380753.2` 62, `GGTLC2` 42, `LOC129529430` 33,
  `LOC101138607` 25.
- 15 strand-confounded (wrong-strand only — the dud class, now excluded).
- 437 thin/tied (identifiability-capped), 18 unexpressed.

**Design implications:**
1. The recoverable target is the **27 genuine** copies; use the strand-aware list, not the 23.
2. PSV columns must come from a **reference-anchored alignment of the copies' exon sequences** (the
   `(E)` gate's "shared graph node" requirement is the blocker even for well-formed families).
3. Re-trace the gap on a GENUINE copy (e.g. `XM_055380753.2`, own=62) before fixing — the dud's
   strand-split may not be the genuine copies' failure mode.
