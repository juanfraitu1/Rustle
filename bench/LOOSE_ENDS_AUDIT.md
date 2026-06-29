# What Is Still Loose — critical pipeline + objectives audit (2026-06-25)

*Adversarial multi-agent audit: 7 critical auditors (one per dimension) → every loose end independently
fact-checked for real / load-bearing → synthesis. 53 agents, **41 confirmed loose ends**. The two most
consequential structural claims (L1, default-off assignment) were additionally hand-verified in source.*

## ✅ Closed so far (2026-06-25, this session)

- **L16** — `promote_hidden_copies.py` verdict `"CONFIRMED-COPY" → "second-haplotype-candidate"`; MILESTONE
  headline + table relabeled with the copy-vs-allele caveat up front.
- **L11/L13** — `asj_findings.md` now LEADS with the **120 transversion** genetic core; the 475 is framed
  as the full candidate set (355 edit-confoundable transitions); the non-binding `frac_mq0` masquerade
  control (removed 0/475) is honestly disclosed with the within-gene-het-vs-paralog-locus TODO.
- **L1** — `denovo_families.py` header + `OBJECTIVES_STATUS.md` O1 row + scope statement now disclose the
  genome-wide catalog was built by the arbitrary `core_recip≥0.13` threshold, NOT the conflict graph.
- **L15 — O4 FP bound MEASURED** (`o4_fp_bound.py`): raw hidden-copy flag fires on **7.39%** of
  single-copy genes ≈ background → non-specific screen; only the 4 MHC survive external check.
- **L6** — `family_def_artifact_filter.py` docstring `(15) → (80)`.
- **Partition table** — `OBJECTIVES_STATUS.md` now has the default-on/validated vs opt-in/prototype table
  (the build-vs-run disclosure).

**Cheap doc-honesty batch (round 2):**
- **L2 — DENSITY GATE WIRED.** `family_split.rs::WEB_MAX_DENSITY` 0.15 → **0.30** (aligned with the
  validated DNA-manifest bar) + `denovo_family_split.py` matched; `Web` families are already excluded from
  copy-assignment (`denovo_pipeline.rs:129`), so this is a real DROP, not a tag. Measured: de-novo split
  695 family/3 web → **691/7**; the 4 newly-dropped are the large multi-chrom over-merges (DSFAM0 =
  164-member ZNF/19 chr, ZNF/APOL). `web_min_size` kept at 10 (vs manifest's 4) to protect small divergent
  families (MAGEB n=7). Tests updated (`classify_web_vs_family` + 629 suite green); shipped
  `denovo_families_split.tsv` re-labeled. **L5 (the cycle/self-recurrence leg) still de-novo-unwired.**
- **L18** — `copy_assignment_theory.md` §7.1: the 1026/1026 GGO figure relabeled as a CIRCULAR consistency
  check in the EASY (MAPQ>0) regime; the load-bearing evidence is the sim5x labeled-truth K-ladder.
- **L12** — `asj_findings.md`: the ~23 distal full-switch (|ΔPSI|=1.0) calls tagged "chimera-not-excluded"
  (no RT-switch guard in the ASJ path); local splice-proximal core unaffected.
- **L21** — `mosaic_discriminator/README.md`: disclosed that the GTF `gene_conversion "confirmed"`
  attribute + `[VG-MOSAIC]` report still label on recurrence alone (discriminator gates emission only).
- **L22** — `input_formats_and_ties.md` operational caveats (temp-BAM size, `$TMPDIR`, single-threaded,
  SIGKILL leak) + `main.rs` now prints the temp-BAM size/path.
- **L19** — already adequate (the exhaustive check is committed, reproducible, viol=0, L=3 scope stated).

**Copy-assignment batch (L7/L8/L9 — round 3):**
- **L7 — CANONICAL ENGINE declared.** `copy_assign::assign_read` via the `combined` pipeline path is THE
  engine: full long-read evidence (PSV columns + copy-specific junction chain — a test proves it
  out-resolves PSV-only), per-molecule/FLAIR-like, assign-or-abstain (no 1/k). The vote engine is its
  flat-error equivalent (kill-test 16/16); the CLI is the same scoring standalone. The margin is now a
  PRINCIPLED operating point **τ = ln((1−p)/p)** (`tau_from_p`/`p_from_tau`/`for_target_misassignment`,
  tested): τ=2.0≡p≈0.12 (recall) and τ≈6.9≡p=1e-3 (Eichler precision) are two choices of p, not arbitrary
  constants. Validated NON-CIRCULARLY on the sim5x labeled-truth ladder: **acc|assigned=1.000 at every
  K≥1, 100% abstain at K=0**. (631→637 tests green.)
- **L8 — abundance CI FIXED.** `copy_assign_pipeline.rs` now bases the CI on INFORMATIVE (decisive)
  coverage, not raw N: `n_eff = #reads with n_decisive≥1`; CI clamped to 0.5 and `=0.5` (full simplex)
  when `n_eff=0` (the K=0 unidentifiable regime). Test `abundance_ci_is_full_simplex_when_unidentifiable`.
- **L9 — silver reconciled.** FLAGSHIP DSFAM817 fixed: was "silver 3/3=100%", measured is **2/3=0.67**
  (circular + thin, 3 unique-mappers); the load-bearing evidence is now stated as the sim5x oracle, not
  the silver. DSFAM102 0/4 already raw.
- **L10 / genome-wide — BLOCKED on L4** (production GGO.bam missing/repointed). Documented as a
  per-family/CLI + sim-validated capability, not a default genome-wide output. `gw_psvlink.sh` built/unrun.

Remaining open (moderate, need steer): L3/L4 (external truth + reproducibility — RESTORE the deleted BAM,
which also unblocks L10/L20), L17 (port `is_disjoint_clique_union` abstain certificate to `--vg`), L20
(run mosaic on real data — needs the BAM). L2's deeper density-gate already wired; L5 (cycle leg) open.

## ⭐ O2 RECOMPUTE on the COMPLETE BAM (2026-06-26) — closes L10, advances L1/L4

Copy-assignment was RE-RUN genome-wide on `GGO_mm.bam` (the complete `-N50 -p0.1` multimapping; the old
`GGO.bam` was the undercount = L4). The recompute required + surfaced 3 fixes that also corrected the O1
family over-merge: (1) **minimap2 PSV-discovery** (~100× faster than poasta; the discrepancy exposed the
over-merges); (2) **same-strand-only** family gate (motif-validated: DSFAM817 = GT-AG/+ vs CT-AC/− =
antisense pair, NOT a family); (3) **disjoint-loci** de-dup (`prune_same_locus`: drops unspliced
read-throughs, tandem-safe). Genome-wide result: **58 real same-strand disjoint-loci families** (was 281
over-merged "families"), 144 copies, 978 collapsed recovered, **72.3% of 108k multimap reads assigned
(93.9% of decisive reads), 23% abstained at the identifiability floor (no 1/k), silver 99.8%**. ⚠ The OLD
flagship/headline numbers (DSFAM817 90%, CAFAM0 213) were on OVER-MERGED false families — RETIRED. See
`bench/COPY_ASSIGN_RECOMPUTE.md`. **L10 closed** (genome-wide run done on the right substrate). **L1
CLOSED**: the de-tie READ-CONFLICT-GRAPH family definition (no similarity threshold) was RUN GENOME-WIDE
(`gw_family_catalog` / `detect_conflict_catalog_genome_wide`, NEW): **82 clean families / 207 copies** (0
mixed-strand, 82/82 single-chrom, real 9/8/7/6/5-copy arrays) — replacing the OLD `core_recip≥0.13`
catalog (281 "families", DNFAM0 = 728-member chr1→chrY over-merge). 82 > the 58 lower-bound (confirming it
was a lower bound). Residual: 82 excludes CROSS-chrom families (colocated_families is same-chrom by design,
so RABL2-class cross-chrom pairs are not captured) + the principled catalog still wants external (gorilla
Compara/OrthoFinder) validation (L3).

## 1. Honest verdict (6 objectives + 2 interests)

**One objective is genuinely ATTAINED on mechanism: O3 (allele-specific junctions)** — the per-molecule
splice-dinucleotide cases (PSMD2 donor 14/14 vs 0/18, DAXX acceptor) are near-tautologically real and
need no external catalog. Everything else is **PARTIAL or LOOSE**, and the dominant pattern is a
**build-vs-run / principled-vs-shipped gap**: the elegant artifacts Canzar would value (the de-tie
conflict-graph family definition, the LLR copy-assignment engine, the Thm-3 disjoint-clique-union
abstain certificate, the mosaic discriminator) are implemented and unit-tested but are **not on the
default `--vg` path and were never run genome-wide**, while the numbers that *are* shipped genome-wide
were produced by the *older, threshold-based* methods those artifacts were meant to replace.

- **O1 (family definition, Interest I)** — the genome-wide deliverable is carried by an arbitrary
  similarity threshold (`T_CORE=0.13`) that produces a 728-member chr1→chrY over-merge — the exact
  failure the principled conflict-graph criterion exists to fix; the conflict graph never ran at scale.
- **O2 (copy-assignment, Interest II)** — headline numbers from a non-production engine; a post-hoc τ
  that governs no operative gate; a silver standard that is circular and empty in the hard MAPQ-0 regime.
- **O3 (ASJ)** — attained on mechanism; phrasing over-reaches (see L11/L13).
- **O4 (reference-absent)** — stands only on the 4 protein-BLAST-confirmed MHC copies; the larger pool
  has no measured FP rate and cannot separate an absent copy from a hyperdivergent allele.
- **O5 (theorem)** — mathematically sound and exhaustively checked in Python, but the proven
  RECOVER/abstain procedure is **not the shipped gate**.
- **O6 (cross-cutting)** — inherits the circular-validation and never-run-at-scale problems.

Net: the *science* is mostly right and the *math* is sound; the thesis is **"built but not yet
defensible at scale,"** and several headline phrasings over-reach their evidence.

## 2. Confirmed loose ends, by objective (ranked load-bearing × cheap-to-fix)

### O1 — Family definition (Interest I)
- **L1 (major).** Genome-wide catalog `denovo_families.tsv` (DNFAM0 = 728 members chr1→chrY) built by
  `core_recip ≥ T_CORE=0.13`; the de-tie conflict graph never ran genome-wide (`family_def_readconflict.tsv`
  = 12 pairs only; `fn recover`/`disjoint_clique` = 0 hits in src/). **Fix:** run `conflict_edges`
  genome-wide, or at minimum state the catalog's provenance. *(moderate)* **[hand-verified true]**
- **L2 (major).** "Family = clique" asserted in prose but the split path accepts density-0.167 Louvain
  blobs as families; the `web` tag (0.15) is decorative, nothing drops on it; artifact-filter's 0.30
  drop is on a *different* (DNA-manifest) pipeline. **Fix:** enforce/align a density gate on the de-novo
  split path, or drop the "clique" word and report the density distribution. *(cheap)*
- **L3 (major).** Genome-wide family validation is structurally circular (minimizer-Jaccard vs a
  minimap2-built universe; mmseqs is the user's own clustering); only external check = 12-pair human
  Compara at 33%. **Fix:** run gorilla Compara/OrthoFinder for one independent positive set; label the
  rest "self-built proxy" / "human spot-check." *(moderate)*
- **L4 (major).** The "N=5 cap heals 5→11 copies, 0 FP" result is non-reproducible today (GGO.bam
  symlink repointed to the uncapped BAM); the 416-family count's source BAM is deleted (stale `.bai`).
  **Fix:** re-run against an existing BAM; restore a genuinely `-N5` view. *(moderate)*
- L5 (minor). Topological artifact filter (15-mer self-recurrence + DROP) is cDNA-bench-only; the de-novo
  catalog gets only a `Web` flag. **Fix:** scope the "VG does triple duty" claim, or port the leg.
- L6 (cosmetic). `FAMILY_DEF_MAX_COPIES` docstring says 15, code default is 80. **Fix:** the one stale string.

### O2 — Copy-assignment under ambiguity (Interest II)
- **L7 (major).** Three uncoordinated engines; the flagship 94%/67% table comes from the Python
  prototype at margin 2.0, not the production vote engine (margin 1) nor the CLI LLR (6.9 nats); τ=6.9048
  matches neither operative gate. **Fix:** pick ONE canonical engine, validate genome-wide, present
  2.0/6.9 as chosen operating points, not a single "calibrated" τ. *(moderate)*
- **L8 (major).** Emitted per-copy abundance CI `1.96·√(θ(1−θ)/N)` shrinks with N exactly in the K=0
  regime the theory proves unidentifiable → false precision on a default user-facing output. **Fix:** make
  the CI track informative-PSV coverage, not raw N; report full-simplex uncertainty below threshold. *(moderate)*
- **L9 (major).** Copy-assignment silver standard is circular (truth = minimap2's own placement) and
  empty in the hard regime (denominators 3–15; DSFAM102 = 0/4); DSFAM817 listed as both 0.67 and 3/3.
  **Fix:** report raw fractions inline, reconcile the discrepancy, promote the sim5x read-name oracle as
  the load-bearing test. *(moderate)*
- **L10 (major).** Both engines are off the default `--vg` path (gated by `RUSTLE_VG_RECOVER_COPIES` /
  `RUSTLE_VG_LAYER2_PSV_LINKAGE` — **hand-verified: psv_linkage bails when unset**); the LLR engine is
  CLI/sim-only; `gw_psvlink.sh` has zero output artifacts. **Fix:** wire a genome-wide batch + report
  coverage, or state plainly that copy-assignment is a per-family/CLI + sim capability. *(moderate)*

### O3 — Allele-specific junctions (the one ATTAINED objective)
- **L11 (major).** The collapsed-paralog control is non-binding (frac_mq0 removes 0/475) and 36% of calls
  are LOC* paralog loci — frac_mq0 can't separate a het from a two-copy PSV in mappable flank. **Fix:**
  run `scan_gene_copy_specific_junctions` on LOC* loci; report within-gene-het vs LOC*-locus counts
  separately. *(cheap)*
- L12 (minor). No RT-switch/chimera guard in the ASJ path; the 23 distal full-switch calls are
  artifact-prone. **Fix:** annotate the distal/dPSI=1.0 class "chimera-not-excluded" (+ optional
  microhomology flag from copy_assign).
- L13 (minor). Only the 120 transversion calls are defensibly "genetic"; the headline 475 includes 355
  edit-confoundable transitions. **Fix:** lead with 120, present 475 as "all candidates."
- L14 (cosmetic). Per-call used/span counts dropped from committed TSVs. **Fix:** re-emit the 4 columns.

### O4 — Reference-absent copies
- **L15 (major, cheap) — ✅ CLOSED (`bench/o4_fp_bound.py`).** FP bound measured: the raw hidden-copy
  flag fires on **7.39% (828/11,206)** of definitionally-single-copy genes — ≈ the genome-wide
  background (7.93%). So the raw flag is a **non-specific screen** (het-dominated), not a copy detector;
  specificity comes from downstream filters and only the 4 protein-confirmed MHC candidates survive an
  external check. Documented in MILESTONE + OBJECTIVES_STATUS; the flagged pool is no longer asserted as
  a copy catalog.
- **L16 (major, cheap).** The detector can't separate an absent copy from a hyperdivergent allele
  (firewall calibrated on 1–2-SNP hets; all 4 MHC positives at 3.9–20.4% div = MHC allelic range), yet
  `promote_hidden_copies.py` prints "CONFIRMED-COPY" / "high" confidence. **Fix:** relabel
  "second-haplotype candidate"; test the firewall on a known multi-SNP het locus. *(cheap)*

### O5 — Identifiability theorem
- **L17 (major).** The proven RECOVER + disjoint-clique-union abstain certificate (refuses on K≥3
  recombination) exists only in Python; production's only guard is the column-count E-gate, so the K≥3
  recombinant regime gets silently assigned. **Fix:** port `is_disjoint_clique_union` as a per-family
  guard so the refusal fires in `--vg`, or scope the doc. *(moderate)*
- L18 (cheap). Thm-2 corroboration (1026/1026) is measured only in the MAPQ>0 easy regime with
  alignment-derived (circular) truth. **Fix:** relabel "consistent-with in the easy regime"; point the
  load-bearing claim at the sim5x labeled-truth ladder.
- L19 (cosmetic). The shipped exhaustive certificate is a bounded re-run (L=3); the full
  238,992-instance enumeration is cited but not committed. **Fix:** commit the regenerator, or down-state.

### Recent additions / cross-cutting
- **L20 (moderate).** The mosaic gene-conversion/RT-switch discriminator has never fired on real GGO
  data (validated by one synthetic-genome unit test; chr19 on/off byte-identical). **Fix:** run
  `MOSAIC_ON+EMIT` on a locus with surviving multi-copy families; report one real reclassification.
- L21 (minor). GTF `gene_conversion "confirmed"` + `[VG-MOSAIC]` stderr emitted on recurrence alone (the
  discriminator gates emission only). **Fix:** thread `classify_event` into the verdict, or correct the README.
- L22 (cheap). SAM/CRAM transcode is single-threaded, writes a multi-GB temp BAM to `/tmp`, leaks on
  SIGKILL; byte-identical verified chr19 only. **Fix:** document `$TMPDIR`/single-threaded, print temp size.

## 3. Top 5 next actions (cheap → "built" to "defensible")

1. **Stop over-labeling in two deliverable strings (cheap, kills two over-claims).** Relabel
   `promote_hidden_copies.py` "CONFIRMED-COPY" → "second-haplotype candidate" (L16), and lead O3 with the
   120 transversion calls / split the LOC*-locus ASJ count from within-gene-het (L11/L13). No recompute.
2. **Run the FP/specificity bound for O4 (cheap, unlocks the milestone).** Scan ~500 single-copy
   segdup-free genes; report the flag rate as an FP bound (L15). Until then assert only the 4 MHC copies.
3. **Disclose the family-catalog provenance + fix reproducibility (moderate, repairs the O1 headline).**
   State `denovo_families.tsv` came from T_CORE=0.13, report the de-tie size distribution alongside (L1),
   and re-derive the family count from a BAM that still exists (L4).
4. **Pick ONE copy-assignment engine and run it genome-wide (moderate, the O2 spine).** Resolve the
   three-engine/τ tangle (L7), make the abundance CI track informative-PSV coverage (L8), report silver as
   raw fractions with hard-regime denominators visible (L9).
5. **Close the theorem-vs-implementation seam (moderate, the O5 point Canzar will probe).** Port
   `is_disjoint_clique_union` as a per-family abstain guard so the refuse-when-non-unique guarantee fires
   in `--vg` (L17), or state plainly production runs a per-read calibrated-margin heuristic.

## 4. Claims that would NOT survive an external check (thesis-credibility risks)

- "We detect families genome-wide **without arbitrary similarity thresholds**" — the shipped catalog was
  produced *by* `T_CORE=0.13` and contains a 728-member chr1→chrY over-merge. (L1)
- "Recovers the truth families (AUC 0.982)" — circular; reproduces minimap2-defined families. (L3)
- "CONFIRMED [MHC] gene-conversion / CONFIRMED-COPY, high confidence" — can't distinguish an absent copy
  from a hyperdivergent MHC allele; DNA parCN never run. (L16)
- "475 ASJ, collapsed-paralog masquerade ruled out (0 collapsed)" — control removed 0/475; 36% are
  paralog loci; defensible headline = the 120 transversion / ~20 splice-proximal core. (L11/L13)
- "The difficult MAPQ-0 case is usually solved, validated" — circular silver standard, empty in the
  MAPQ-0 regime; 94%/67% from a Python prototype. (L7/L9)
- "Theorems 1–3 proven AND executable [in production]" — the abstain certificate is Python-only; the
  shipped `--vg` path silently assigns in the K≥3 regime the theorem exists to refuse. (L17)

## 5. "Shipped" features that are actually default-off and never validated genome-wide
(all confirmed gate-unset in source)
- `RUSTLE_VG_RECOVER_COPIES` — the PSV-linkage assignment path (**hand-verified: bails when unset**).
- `RUSTLE_VG_LAYER2_PSV_LINKAGE` (behind `RUSTLE_VG_LAYER2`) — `rustle --vg` emits zero copy-assignment by default.
- `RUSTLE_VG_MOSAIC_ON` / `_EMIT` — the gene-conversion/RT-switch discriminator; never observed firing on real data.
- `RUSTLE_READCHAIN_TRIM_TERMINAL` — terminal-exon inflation stays live by default in `--read-coherence`; trim measured chr19-only.
- The CLI LLR copy-assignment engine + `gw_psvlink.sh` / injection-FP experiments — built, zero output, never run.

**Cheapest honesty fix:** a one-table partition of "default-on / validated" vs "opt-in prototype" in
`OBJECTIVES_STATUS.md` (currently absent) — makes the build-vs-run gap a disclosed design choice rather
than a discovered weakness.
