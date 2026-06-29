# P1–P4 (defense-readiness must-dos) — results

2026-06-28. Working through the DEFENSE_READINESS_AUDIT must-dos here (P0/SEDEF deferred to the cluster).

## P2 — O4 gate-5: `asm20` → `-x splice` (DONE, sim-verified)
The remap gate (`absent_copy.rs::remap_identity_minimap2`) shelled `minimap2 -cx asm20` with a SPLICED copy
consensus as query against the intron-bearing genome — `asm20` (non-spliced) cannot chain real multi-kb
introns, so a genuine multi-exon copy would fail to align and be wrongly routed to DnaNeeds. **Fixed to
`-cx splice`** (the cDNA-to-genome preset; `de:f:` then reflects exonic divergence, introns spliced out, not
counted as mismatch). 683 lib tests green; **`bench/absent_copy_sim.py` still 4/4 PASS** (AC_* copy admitted,
60 reads `status=assigned` to it, SIM_HET → DnaNeeds `<3 clusters`) — the fix does not regress the synthetic case.
Caveat carried: each remap call re-indexes the 3.5 GB genome (~1–2 min) → the ON path is perf-bound on real
data (an `.mmi` pre-index is the future optimization); the real-data admitted-copy attempt is below.

## P3 — non-circular O2 accuracy point + reconcile the 20%-vs-100% tables (DONE)
**Pinned `smoke_sim5x_ground_truth` in CI** (`copy_assign_pipeline.rs`): removed `#[ignore]`; it now runs in
the normal suite (early-returns harmlessly without `RUSTLE_SIM5X_DIR`) and ASSERTS the identifiability ladder
when the sim5x data is present: **K=0 → 100% tied / 0 assigned; K≥2 → acc|assigned ≥ 0.99**. Run with the data
(`RUSTLE_SIM5X_DIR=…/sim5x`) — **PASSES**. The ground-truth ladder:

| K | reads | PSV cols | resolvable% | acc\|assigned | acc\|argmax | tied% |
|---|---|---|---|---|---|---|
| 0 | 1000 | 0 | 0.0% | — | — | 100.0% |
| 1 | 1000 | 1 | 12.0% | 1.000 | 1.000 | 88.0% |
| 2 | 1000 | 2 | 20.0% | 1.000 | 1.000 | 80.0% |
| 4 | 1000 | 2 | 20.0% | 1.000 | 1.000 | 80.0% |

**Reconciliation: "20% vs 100% for K≥2" is NOT a data contradiction — it is a metric conflation.** At K≥2,
the **resolvable fraction is ~20%** (only reads that SPAN a distinguishing PSV can resolve; the synthetic reads
are short and PSVs are 2 sparse columns, so ~80% don't span one and are *correctly* Tied), while the
**accuracy on the assigned reads is 100%**. So the headline "K≥2 → 100%" is the *accuracy* (correct, defensible:
when the gate assigns, it is right), NOT the resolvable/assigned fraction. The honest statement is
**"K≥2 → 100% accuracy on the ~20% of reads that are resolvable; the rest abstain (Tied), not guessed."**
The ~20% is coverage-limited (read length × PSV density), not a method limit — on real GGO with longer reads /
denser PSVs the resolvable fraction is much higher (the definitive O2 = 75% assigned). Docs that state
"K≥2 → 100%" must read as accuracy, never as "all reads resolve."

This is a *non-circular* accuracy point: sim5x reads carry their TRUE copy in the read name (labeled ground
truth, not the circular silver standard), and the assertion (acc|assigned ≥ 0.99) is now enforced in CI.

## P4 — O3 masquerade separator on the LOC* calls (DONE)
The 120 transversion core (genetic ASJ) splits **76 at non-LOC genes + 44 at LOC\* paralog loci (18 distinct
genes)**. The audit's worry: at LOC\* loci the het-anchor "allele-specific" signal can be a paralog **copy**
masquerade (allele vs copy unresolvable from the het partition alone). I ran the clean separator —
`scan_gene_copy_specific_junctions` (`asj --mode copy`, partitions reads by PSV/COPY, not by het allele) — on
GGO_mm.bam over the 18 LOC\* gene windows:

- **17 of 18 LOC\* genes produce a copy-specific junction** (55 total, q<0.05 & |ΔPSI|≥0.3) → these are real
  multi-copy paralog loci where the allele-specific call is **copy-confounded (masquerade live, needs DNA)**.
- **1 of 18 (LOC101138206)** shows NO copy-specific junction → consistent with a **genuine within-gene het**.

**Honest two-count split (the deliverable):** of the 120 transversion "genetic core", **76 non-LOC calls are
the clean genetic core**; the **44 LOC\* calls (18 genes) are paralog-masquerade-suspect — 17/18 genes
copy-confounded, 1 clean.** So the defensible genetic-ASJ count is **~76 (non-LOC) + 1 clean LOC ≈ 77**, NOT
120, with the 44 LOC\* flagged copy-confounded. The **~20 splice-proximal dinucleotide calls** (PSMD2/DAXX-class)
remain mechanistically airtight regardless — their base-level motif disruption is independent of copy structure.
Data: `p4_loc.asj.tsv`, `p4_loc_regions.txt`. (The full 475/120 *recompute on GGO_mm.bam* — a genome-wide
single-mode rerun — is the remaining cluster-scale piece; the masquerade *separation* above is on GGO_mm.bam.)

## P1 — O2 on the principled conflict-graph catalog (DONE)
Ran the significance gate (`copy_assign --skip-poa-diagnostic`) on the **threshold-free de-tie conflict-graph
catalog** (`gw_conflict_catalog`, 82 same-chrom families → 106 detected within their spans, 206,186 reads;
102 min, all 82 regions):

| metric | principled conflict-graph catalog | (cf. annotation-refined co-located subset) |
|---|---|---|
| assigned | **63.9%** | 75.1% |
| ambiguous | 0.5% | 0.0% |
| certified-tied | **35.7%** | 24.8% |
| **of DECISIVE reads assigned** | **99.3%** | 99.9% |
| silver-standard | **99.8%** | 99.9% |

**The catalog story (answers the "killer question").** Report the **principled number as the genome-wide
headline: 63.9% assigned / 35.7% certified-tied** — the conflict-graph catalog is the elegant artifact (no
similarity threshold), so the headline and the principled method are now the **same object**, closing the
build-vs-run gap. The 75.1% is the *annotation-refined co-located subset* (cleaner because refinement drops the
Alu-bridge over-merges and harder families → fewer tied), and must be labeled as such, not as "the genome-wide
O2." **The DECISION RULE is identical on both** (99.3–99.9% of decisive reads assigned, silver ~99.8%, ~0%
ambiguous, no 1/k) — only the *tied* fraction moves (more genuinely unresolvable / K=0 families survive in the
unrefined catalog). So O2's defensible genome-wide claim: **"on the principled threshold-free catalog, 99.3% of
reads carrying any copy-distinguishing evidence are assigned with a calibrated certificate; 35.7% are
certified-tied (abstained, not guessed); silver 99.8%."** Data: `p1_conflict_o2.*`.
