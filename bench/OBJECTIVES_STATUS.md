# Objectives status — what's attained, what's loose (grounded audit, 2026-06-25)

Honest assessment from a 6-objective audit against the actual artifacts (skeptical thesis-committee /
Canzar lens). **Verdict: not yet a clean "objectives attained" across the board — but the gaps are
NOT "the method doesn't work."** They are: over-claiming, built-but-unrun validations, an external
ground-truth gap, one fixable input flaw, and the (fundamental) DNA het-vs-copy wall.

## Per-objective

| # | objective | status | headline gap |
|---|---|---|---|
| 1 | RNA-level multi-copy family detection (~R ∩ ~B) | **PARTIAL→ conflict-graph catalog now RUN genome-wide** | ⭐ The principled de-tie READ-CONFLICT-GRAPH family definition (no similarity threshold) was RUN GENOME-WIDE for the first time (`gw_family_catalog` → `detect_conflict_catalog_genome_wide`): **82 clean families / 207 copies** (0 mixed-strand, 82/82 single-chrom, real 9/8/7/6/5-copy arrays), replacing the OLD `core_recip≥0.13` catalog (281, DNFAM0=728-member chr1→chrY over-merge). Closes L1. Still: 82 excludes cross-chrom families (colocated_families is same-chrom) + needs external (gorilla Compara) validation. bench/COPY_ASSIGN_RECOMPUTE.md |
| 2 | Copy assignment under ambiguity (PSV gate + AS-decisive) | **PARTIAL→ RECOMPUTED genome-wide on the COMPLETE BAM** | Canonical engine (L7, full PSV+junction, τ=ln((1−p)/p), assign-or-abstain) RAN GENOME-WIDE (via the `copy_assign` CLI; still `RUSTLE_VG_RECOVER_COPIES`-gated off default `--vg`) on `GGO_mm.bam`: **106 families / 206,186 reads on the principled threshold-free conflict-graph catalog (`gw_conflict_catalog`): 63.9% assigned / 0.5% ambiguous / 35.7% certified-tied / 99.3% of DECISIVE reads assigned / silver 99.8%** (silver = circular consistency check — agreement with minimap2's own primary placement, NOT accuracy; non-circular accuracy = sim5x labeled ladder K≥2 → acc\|assigned=1.000 on the ~20% of reads that are resolvable; K=0 → 100% Tied). Note: 75.1% cited elsewhere = annotation-refined co-located SUBSET, not genome-wide. Required 3 fixes that also corrected O1: minimap2 PSV-discovery (~100× faster), same-strand-only + disjoint-loci family gates (motif-validated). ⚠ OLD over-merged headlines (DSFAM817 90%, CAFAM0 213) RETIRED. bench/COPY_ASSIGN_RECOMPUTE.md |
| 3 | Allele-specific junctions | **ATTAINED** | the clean, committee-ready result — the natural headline |
| 4 | Reference-absent copies (this milestone) | **PARTIAL** | FP rate NOW QUANTIFIED (`o4_fp_bound.py`): raw hidden-copy flag fires on **7.39%** of definitionally-single-copy genes ≈ background 7.93% → raw flag is a non-specific SCREEN (het-dominated), not a copy detector; only the 4 protein-confirmed MHC candidates survive external check; copy-vs-allele needs DNA |
| 5 | Identifiability theorem (through-line) | **PARTIAL** | Theorems 1–4 proven + machine-checked B1–B7 (Thm 4 = bridge: production min_p gate is a sound per-read identifiability certifier for all K≥1, making theory load-bearing for the shipped method — under the explicit completeness precondition origin(r)∈C, machine-checked necessary [B6: dropping it → confident misassignment, the O4 hazard] and orthogonal to K≥3 cover non-uniqueness [B7]; RECOVER itself not run); Strong Separation is *sufficient not necessary* (true boundary = recombination-freeness, no closed form) |
| 6 | Cross-cutting: external validation + reproducibility | **PARTIAL** | validation leans on internal/circular validators; no single end-to-end pipeline |

## ⭐ Default-on / validated  vs  opt-in prototype  (the build-vs-run partition)

The single most important honesty disclosure (per the 2026-06-25 loose-ends audit, `LOOSE_ENDS_AUDIT.md`):
the principled artifacts the thesis is *about* are **implemented + unit-tested but OFF the default `--vg`
path and NOT run genome-wide**; the numbers shipped genome-wide came from the *older threshold methods*
they were meant to replace. State this up front; do not let a reviewer discover it.

| capability | default-on & validated at scale? | reality |
|---|---|---|
| StringTie-faithful assembly (baseline) | ✅ yes | the substrate; genome-wide, parity-tested |
| Allele-specific junctions (O3) | ✅ yes (mechanism) | the genuinely-attained result; ~77 genetic core (44 LOC* copy-confounded — masquerade separator RUN, 17/18 LOC* genes copy-specific; ~20 splice-proximal airtight) |
| de-tie conflict-graph **family definition** (O1) | ❌ panel-only (12 pairs) | shipped catalog uses `core_recip≥0.13` instead; conflict graph per-region, never genome-wide |
| **copy-assignment** canonical engine (O2) | ⚠ canonical declared + sim-validated; genome-wide blocked | `copy_assign::assign_read` (combined path) is canonical: full PSV+junction evidence, assign-or-abstain, principled τ=ln((1−p)/p) gate; sim5x labeled-truth acc\|assigned=1.000 @K≥1. Still default-off in `--vg` (`RUSTLE_VG_RECOVER_COPIES`) and genome-wide unrun (L4 deleted BAM) — a per-family/CLI + sim capability, not a default genome-wide output |
| Thm-3 disjoint-clique-union **abstain certificate** (O5) | ❌ Python-only | production `--vg` has no uniqueness certificate → silently assigns in the K≥3 recombinant regime the theorem exists to refuse |
| gene-conversion vs RT-switch **mosaic discriminator** | ❌ opt-in, never fired on real data | `RUSTLE_VG_MOSAIC_ON/_EMIT`; microhomology leg live but unobserved on GGO; DNA leg measured + deliberately not wired |
| read-coherence **terminal-exon trim** | ❌ opt-in, chr19-only | `RUSTLE_READCHAIN_TRIM_TERMINAL`; default `--read-coherence` GTF still carries terminal inflation |
| reference-absent **copy detector** (O4) | ⚠ runs, but unquantified | no FP rate; only the 4 protein-confirmed MHC candidates survive an external check |
| genome-wide PSV-linkage / injection-FP validations | ❌ built, never run | `gw_psvlink.sh` etc. exist with zero output artifacts |

## What is SOLID (the defensible core)

- **ASJ (O3) is genuinely attained end-to-end**: one phase-and-test engine → 475 single-anchor ASJ;
  the defensible genetic core is **~77** (full transversion set 120, minus 44 copy-confounded LOC\*
  calls — paralog masquerade, allele-vs-copy needs DNA), plus a multi-SNP superset and 146
  copy-specific differential junctions. Per-molecule allele→junction **linkage** is the load-bearing
  result (PSMD2 14/14 vs 0/18; DAZ1 vs DAZL dPSI 0.918, q=2.6e-151), and the genetic-vs-RNA-edit
  confound is controlled. **Mechanism caveat (genome-verified, `bench/asj_motif_check.py`):** the
  flagship anchors sit at donor−1 / the exon boundary — splice-REGION (extended-consensus) variants,
  NOT the invariant GT-AG dinucleotide (0/475 on a core dinucleotide; the dinucleotide is intact). The
  earlier "textbook splice-site / creates-destroys the motif" framing was genome-false and is retracted.
  **This is the thesis headline.**
- **The copy-assignment kill-test** is a principled, Canzar-aligned finding: votes ≡ LLR (16/16, monotone
  at flat error); the lever is the GATE (min_psv 3→1), not the scoring. Clean.
- **Family definition** is rigorous *where scoped*: error-model-derived constants (not fitted), perfect
  on the n=17 panel, three structural barriers that hold by construction; APOBEC3 correctly EXCLUDED
  despite being a Compara paralog (read-resolvable) = honesty discipline.
- **Identifiability theorem**: Theorems 1–4 formally proven + machine-checked (Thm 4 / bridge: production min_p gate is a sound per-read identifiability certifier for all K≥1, making the theory load-bearing for the shipped method; RECOVER itself not run); turns the
  negatives (winnowmap/short-read/aligner-invariance) into a contribution (the limit is identifiability).
- **Reference-absent milestone**: 4 endogenous MHC copies (protein-confirmed, contamination-ruled-out)
  + 15 multi-mapping-supported dispersed-paralog candidates, landing in MHC/PRDM9/ZNF as biology predicts.

## The loose ends, prioritized

**Tier 1 — blockers (undermine a core claim if unaddressed):**
1. **External ground-truth at scale (O1/O4/O6).** Every genome-wide validation leans on a circular
   universe (minimizer-Jaccard vs minimizer-Jaccard) or a self-built validator (the user's own mmseqs);
   Compara is a 12-pair human proxy. *Close:* frame protein-homology as corroboration (not truth);
   adjudicate the top 5–10 unvalidated high-read edges (OCLN~SEPTIN7, BCAS4~CCDC30) by split-read/synteny;
   ideally one external cross-check.
2. **The N=5 input flaw (O1).** minimap2's default secondary cap fragments REAL >6-copy families — exactly
   the DAZ/RBMY-class arrays the thesis cares about (verified: re-align `-N50 -p0.1` heals 5→11 copies,
   0 FP). *Close:* re-align GGO.bam and rebuild the conflict graph (~1–2 h), or emit >6-placement arrays
   as explicit incompleteness warnings.
3. **Genome-wide FP rate for reference-absent (O4).** The 73-candidate catalog has no measured specificity.
   *Close:* inject divergent synthetic copies at known loci, run the pipeline, measure sensitivity/precision.

**Tier 2 — built-but-never-run (cheap, high-value):**
4. **Run the genome-wide PSV-linkage validation (O2).** `gw_psvlink.sh` / `gw_psvlink_aggregate.py` exist
   but were never executed. *Close:* run per-chrom (watch the ~18 GB OOM); report PSV-net-new vs VG baseline.
5. **Validate per-copy abundance vs truth (O2).** `copy_abundance` (EM + CI) is emitted with zero accuracy
   check — a live fabrication risk (prior RBMY work showed confidence anti-correlated with identifiability).
   *Close:* sim5x abundance sweep, EM-estimate vs known per-copy fractions; else label "exploratory."

**Tier 3 — framing fixes (no computation, pre-empt the objection):**
6. **Re-scope O1**: call it "read-evidence-based recent-paralog / read-confusable copy detection" — a
   copy-assignment *substrate*, not an evolutionary taxonomy. State up front.
7. **State the DNA/parCN limit and the het-vs-copy wall up front** (O2/O4/O6) as the honest boundary, not
   defended under questioning. The 66 single-locus reference-absent candidates and the MHC copy-vs-allele
   are RNA-unresolvable by design — DNA parCN is the named resolver.
8. **Acknowledge Strong Separation is sufficient-not-necessary** (O5); scope the theorem claim to a
   conservative, provable sufficient condition with a polynomial certifier.

**Tier 4 — polish / coherence:**
9. **One reproducible end-to-end pipeline** (O6): `scripts/run_*.sh` BAM → families → assignments →
   reference-absent catalog, replacing scattered bench scripts on scratch artifacts.
10. Junction-decisive resolution wired into production (O2); DAZ MAPQ-0 ASJ re-check (O3); exon-union
    re-align all 145 families for the 90.3% bound (O5).

## The honest scope statement (what you can stand behind today)

*"At the RNA level, from great-ape long reads, we (a) detect read-confusable recent-paralog copy groups
— a copy-assignment substrate, not an evolutionary taxonomy (NB: the principled threshold-free
read-conflict-graph definition is validated on a panel and runs per-region; the genome-wide catalog
shipped to date still uses a similarity threshold and is being migrated to the conflict criterion); (b)
assign reads to individual copies up to a formally-characterised identifiability limit (Strong Separation,
a provable sufficient condition), with a calibrated decisive-margin gate; (c) model allele-specific and
copy-specific junctions per molecule without phasing — the fully-validated result; and (d) flag expressed
gene-family copies absent from the reference (confirmed in the MHC). The boundary throughout is
information-theoretic identifiability, which neither aligner choice nor read depth crosses; resolving
het-from-copy and absolute copy number requires DNA, which we do not claim from RNA."*

## Minimal closing sequence

1. Framing fixes (#6–8) — hours, no compute, pre-empt the biggest objections.
2. Run the built validations (#4, #5, #3) — 1–2 days, tooling exists, converts "shipped" → "shipped + attested."
3. Fix the N=5 input flaw (#2) — ~2 h, removes a silent undercount on the flagship arrays.
4. External-truth adjudication (#1) — days, caps how strongly genome-wide counts can be asserted.
5. One reproducible pipeline (#9) — days, makes it a thesis artifact not a script pile.
