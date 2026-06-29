# Adversarial review #3 — grounded in the committed literature (2026-06-28)

Re-run of the per-objective adversarial review, this time asking of each objective: *what do the prior-art
papers in memory do that we don't, what do they reveal we're missing, and what concretely can be done?* Each
finding is tagged with the paper, a verified code-state (**DONE / PARTIAL / GAP**, checked against the repo this
session), and a concrete action with rough **value × effort**.

Papers: Canzar 2016 (advisor), Soto 2025 Cell, Eichler/Guitart 2024 TBC1D3, Vollger 2019 SDA, Huang 2026
longcallR, Sahlin 2018 IsoCon, Zheng 2025 Clair3-RNA, family-detection prior-art (OrthoFinder/SEDEF).

---

## The high-leverage findings (what can actually be done)

### F1 — O4 positive control: mask a real copy, show re-admission  **[GAP → now feasible]  value HIGH · effort LOW**
- **Paper:** Soto 2025 — 668 paralogs (37%) sit in regions *missing or erroneous in GRCh38*; reference-absent
  copies are real at the DNA level, they're just a pre-T2T problem. Our bounded O4 run admitted **0** on the
  T2T-complete gorilla — which is either "machinery works, reference is complete" **or** "machinery is inert."
  We cannot currently tell which.
- **Code-state (verified):** GAP. Positive evidence is synthetic only (`absent_copy_sim.py` plants A1/A2
  haplotypes) + unit tests with mocked remap. `docs/VG_OBJECTIVES_AND_ROADMAP.md` Obj-2 lists the exact missing
  test: "mask one copy out of the reference, confirm the scan recovers its reads." Never run.
- **Action:** take one of the 4 protein-confirmed MHC copies (or a clean RABL2/DAZ copy), **delete it from the
  reference FASTA**, rebuild the `.mmi`, realign that locus's reads, run `--absent-copies`, and verify O4
  **re-admits it as an `AC_` copy** (not DnaNeeds). Decisive: turns "0 admitted" from ambiguous into "0 because
  the reference is complete, and here's proof the detector fires when a copy IS absent." The `.mmi` + bounded-run
  infrastructure from H4 makes this cheap now.

### F2 — O1 orthogonal, NON-circular family validation without SEDEF  **[PARTIAL/GAP]  value HIGH · effort MED**
- **Paper:** Soto 2025's family definition is DNA+CN: minimap2 map-back of SD regions → BEDTools intersect
  **-f 0.99** (keep mappings covering >99% of each exon = *shared exons*) → group genes whose family copy-number
  (famCN, read-depth) has mean-abs-deviation < 1. This is the field-standard recipe and needs **no SEDEF**.
- **Code-state (verified):** our O1 validation = Compara (only **12 mappable pairs**, 33% confirm), mmseqs
  protein (837 families), and minimap genomic-span homology (89.2% — but *partly circular*, the span contains the
  same exons). SEDEF is deferred to the cluster. There is **no famCN / shared-exon analog**.
- **Action:** implement the Soto recipe directly on the T2T gorilla genome we already have: map family-member
  exons back with `minimap2 -c --eqx -N50 -p0.5`, require ≥99% exon coverage for "shared exon," and add a
  read-depth famCN consistency check (MAD<1) per family. This is an **independent DNA/CN witness** for the RNA
  families that is *not* circular with the RNA conflict graph, and it's the exact method our own key reference
  (Soto) uses — strong to cite. Gets most of the SEDEF value without the cluster.

### F3 — O2 non-circular accuracy via the DNA-derived PSV catalog as a SUPERVISED prior  **[GAP]  value HIGH · effort MED-HIGH**
- **Paper:** Vollger SDA 2019 — DNA reads pre-phase PSVs→paralogs by correlation clustering; our identifiability
  theory proves the conditions under which that recovery is *exact*. The natural consequence (noted in memory):
  DNA **pre-phases PSVs→copies, turning the NP-hard RNA phasing into a supervised nearest-signature lookup.**
- **Code-state (verified):** GAP. The DNA-derived PSV catalog (24,256 pairs, 86% concordance) is **validation-
  only**; the per-read "signature decoder" (Phase 2 of the DNA-PSV-catalog spec) is explicitly *deferred*. RNA
  assignment uses only the read's own PSVs vs RNA-built profiles — never the DNA signatures as a prior.
- **Action:** wire Phase 2 — assign each RNA read to the DNA-defined copy whose PSV signature it matches
  (supervised), and report accuracy against the DNA labels. This is the **genuinely non-circular accuracy** the
  defense audit asked for (silver is circular; sim5x is synthetic; DNA labels are an independent oracle). Heavier
  than F1/F2, but it's the cleanest answer to "what's your real-data accuracy?"

### F4 — Theory: the Canzar-shaped capstone (LP-rounding approximation for the NP-hard cover)  **[GAP]  value HIGH (advisor) · effort HIGH**
- **Paper:** Canzar 2016 — his signature contribution is not just NP-hardness but an **LP relaxation + rounding
  with a provable approximation ratio** (0.19-approx, dependent rounding) for the multimapping-resolution
  problem cast as **maximum facility location**. We have Lemma 1 (MCC=χ(H)), Thm 1 (NP-hard), Thm 2–3 (exact
  under Strong Separation), Thm 4 (per-read bridge). We have **no approximation algorithm for the hard general
  case** — exactly the gap his own paper fills.
- **Action (research, not a quick fix):** cast copy-assignment as Canzar-style facility location (reads=clients,
  copies=facilities, conflict=incompatible PSV signatures) and give an **LP-rounding approximation with a ratio**
  for general MCC, or an approximation for the minimum copy-cover. This is the single most advisor-aligned move:
  it completes the "clean problem + hardness + *approximation guarantee*" arc in his exact taste. Flag as the
  thesis theory capstone, scoped separately.

---

## The cheaper robustness / citation wins

### F5 — O3 strand-bias (SOR) filter for ASJ  **[GAP]  value MED · effort LOW-MED**
- **Paper:** longcallR 2026 filters allele-specific junctions with **StrandOddsRatio ≥ 2** (plus Fisher + BH-FDR
  < 5%) to kill calls driven by reads from a single strand — a standard, expected robustness gate.
- **Code-state (verified):** GAP. Our ASJ gates are anchor-quality + balanced-het + min-span + |dPSI|≥0.3 +
  BH-FDR<0.05 + transversion/editing flags — **no strand-bias check**, and `AlignedRead` doesn't even carry
  `is_reverse` at the tally point (BAM flag read but never propagated). ~50–100 line change across 3–4 files.
- **Action:** thread read strand into the junction tally and add a SOR (or simpler per-strand balance) filter.
  Closes a robustness gap a reviewer who knows longcallR will immediately ask about; cheap and citable.

### F6 — O2 head-to-head: our significance gate vs the Eichler AS≥10 rule on the same reads  **[PARTIAL]  value MED · effort LOW**
- **Paper:** Eichler/Guitart 2024 — the assignment rule the advisor keeps citing: assign an Iso-Seq read to a
  paralog cluster iff its best minimap2 AS beats every other by **≥10**, else mark ambiguous and ignore
  (conservative, no 1/k).
- **Code-state (verified):** PARTIAL. We *calibrated* τ=6.9 to "≈ the Eichler AS≥10 operating point" but there
  is **no direct AS≥10-rule-vs-our-gate benchmark** on the same reads.
- **Action:** implement the literal AS≥10 rule and run it head-to-head against our significance gate on sim5x
  (labeled truth) + a GGO family. Expected story: same decisive calls, but our gate adds the per-read
  identifiability *certificate* (min_p) and abstains *with a reason* where AS≥10 silently discards. A clean
  one-figure exhibit that positions us as a principled generalization of his favourite rule.

### F7 — O2 IsoCon-style iterative *candidate* (copy) pruning  **[PARTIAL]  value MED · effort MED**
- **Paper:** IsoCon 2018 — after per-variant-position significance, it **iteratively removes non-significant
  candidates, reassigns their reads, and retests** until all surviving candidates are significant. We adopted the
  per-read least-significant-p test, but the *family layer* takes copies from the conflict graph and doesn't run
  this iterative copy-pruning (only the K=0/min_p=1 tie handles fully-identical copies).
- **Action:** add an IsoCon-style pass that merges/drops a copy not significantly distinct from its nearest
  neighbour and reassigns its reads. Connects directly to the K-frontier and to F3 (a copy with no DNA signature
  shouldn't be a separate copy). Tightens copy *count* (the precision axis the gate alone doesn't bound).

### F8 — Clair3-RNA DP/AD callable-region benchmarking (replace circular silver)  **[PARTIAL]  value MED · effort LOW-MED**
- **Paper:** Clair3-RNA 2025 — evaluate only in **callable ∩ adequate-coverage ∩ GIAB** regions, report by **DP
  (depth) + AD (allele depth)** at min 4×/10×; normalize for uneven RNA coverage. The editing filter we already
  shipped is from this paper.
- **Code-state:** PARTIAL — editing filter DONE; abundance/accuracy still leans on circular silver. The held-out
  CV (H3) is a good start.
- **Action:** report copy-assignment accuracy stratified by DP/AD in callable regions, as the honest non-circular
  evaluation frame (complements sim5x + the H3 held-out CV).

---

## What the literature CONFIRMS we already got right (defensive ammunition)

- **TBC1D3 (Eichler):** Minigraph/VG *failed* to group near-identical paralogs (isolated single-haplotype
  nodes) → they abandoned the graph for phylogenetic groups. This independently corroborates our VG-headroom
  refutations — the variation graph is **not** the copy-resolver. Cite it when defending the PSV/score approach.
- **SDA (Vollger):** hits the **same K=0 identifiability floor** we proved ("virtually identical duplications
  cannot be distinguished, need >100 kb reads"). Our identifiability theorem *explains* their empirical floor.
- **longcallR / Clair3-RNA:** both operate in the **uniquely-mappable** regime (longcallR drops HLA, MAPQ≥10);
  they assign to genes+haplotypes, never to family copies. The collapsed-paralog/multimapper regime we target is
  exactly the part they avoid — our niche is real and unoccupied.
- **IsoCon:** explicitly leaves "distinguish alleles from distinct copies" and "de-novo genome-wide discovery"
  as open — our ASJ (copy-vs-allele) + de-novo conflict-graph family definition are precisely those gaps.
- **SEDEF = Jaccard + chaining** = the advisor's own "maximize the Jaccardian + chaining" intuition is a
  published DNA SD-detection paradigm; our family work is the RNA-side analog.

---

## Prioritized "what can be done" (recommendation)

| # | Action | Objective | Value | Effort | Decisive? |
|---|--------|-----------|-------|--------|-----------|
| **F1** | Mask a real copy → show O4 re-admits it | O4 | HIGH | LOW | ✅ yes |
| **F2** | Soto famCN/shared-exon DNA validation (no SEDEF) | O1 | HIGH | MED | ✅ yes |
| **F6** | AS≥10 vs our gate head-to-head | O2 | MED | LOW | exhibit |
| **F5** | SOR strand-bias filter for ASJ | O3 | MED | LOW-MED | robustness |
| **F3** | DNA catalog as supervised prior (non-circular acc) | O2 | HIGH | MED-HIGH | ✅ yes |
| **F7** | IsoCon iterative copy-pruning | O2 | MED | MED | precision |
| **F8** | DP/AD callable-region benchmarking | O2 | MED | LOW-MED | honesty |
| **F4** | Canzar LP-rounding approximation (theory capstone) | Theory | HIGH | HIGH | advisor |

**Recommended next sitting:** F1 + F6 (both LOW effort, both decisive/exhibit, both directly answer a "is it
inert / is it just calibration?" doubt) — then F2 (the strongest non-circular O1 win available without the
cluster). F4 is the thesis theory capstone and should be scoped as its own effort.
