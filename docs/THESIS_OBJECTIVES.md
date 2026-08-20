# Thesis objectives — multi-copy gene family assembly (VG / EM)

> ⚠⚠ **SUPERSEDED 2026-08-07. This file's five-objective VG/EM decomposition is NOT the current scope.**
> Agreed with the advisor, the thesis now has **THREE** objectives:
>
> | # | objective |
> |---|---|
> | **O1** | Define a multi-copy gene family topologically at the RNA level (quasi-clique in E_r; MCC = χ(H)) |
> | **O2** | Decide, for a read at a multi-copy locus, whether the evidence warrants assigning it to a COPY AT ALL — and abstain when it does not. Contested set = **alignment-score near-ties**, not MAPQ-0. No 1/k. |
> | **O3** | Detect + flag expressed transcript paths not explained by represented reference copies, **STRATIFIED by whether the orphaned reads have anywhere to go**. Detect-and-flag with a measured FPR per stratum; **completeness is never claimed**. (this was O4) |
>
> **The allele-specific-junction objective is DROPPED** — cut for time, and because it does not connect
> to the others. ⚠ It was the only objective the 2026-06-25 audit rated ATTAINED; the ASJ result itself
> is de-scoped, not retracted.
>
> ⚠ **"O3" denotes three different things across this repo's history** — the EM below, ASJ in the
> 06-25 audit, and reference-absent copies now. Resolve the scheme before quoting any objective number.
> The numbered sections below are kept for provenance only.
>
> ### ⭐ O3 RESTATED (2026-08-19) — stratify the target, bound each stratum
>
> O3 as originally posed — *"find copies missing from the genome"* — is **not achievable**, and the
> reason is measured, not conjectural. In the whole-genome excision control (one copy of 162 two-copy
> families deleted, matched IsoSeq) a deleted copy has **two fates**:
>
> | fate | rate | what the reads do |
> |---|---:|---|
> | **ORPHANED** | 33.3% | median **92.7% of reads unmapped** — detectable |
> | **ABSORBED** | 64.2% | reads land on the best paralogue at **1.75× depth**, concordance 0.967 — invisible to any unmapped-read method |
>
> **Expression is not the constraint** (99.34% of copies clear the floor); **where the reads go** is.
> ⚠ And O3's original target class — a collapsed paralogue — sits in the ABSORBED stratum.
>
> This is the same move that rescued O2: restate the target population and the claim, then bound each
> stratum honestly.
>
> **Target.** Copies absent from the assembly, **stratified by whether the orphaned reads have anywhere
> to go.** **Claim.** Detect-and-flag with a measured FPR, stated **per stratum**:
>
> | stratum × route | status | bound |
> |---|---|---|
> | unique sequence, **unmapped-read** route | **works** | **M ≤ 6.4** missing expressed copies |
> | paralogous sequence, **unmapped-read** route | ⚠⚠ **vacuous** | π = 1/35 = 0.0286, 0/26 at cov ≥ 0.8, **formally unbounded** |
> | paralogous sequence, **depth (S2)** route | **partial** | held-out **TPR 0.2703 / FPR 0.0200**; sensitivity set by **divergence, not abundance** (0.4500 above 0.01 divergence vs **0.0588** below) — ⚠ and **45.78% of positives lie below 0.01** |
>
> **Say this, not more:** *the instrument flags candidates with a measured false-positive rate in a
> named stratum, and has explicitly no power for unmapped-read detection in the collapsible stratum.*
> **Never claim completeness.** The signature is **UNMAPPED READS, not clipping** (34.53% pooled,
> MAPQ-60 before deletion) — no published collapse detector uses clipping.
>
> ⚠ **The one real candidate does not close:** STON1+GTF2A1L, ~116.7 kb absent from mGorGor1, 125
> near-full-length unmapped reads, gapless chromosome, 0 GFF lines, present in chimp and orangutan —
> but **single-copy, n = 1, p = 0.055, UNCONFIRMED**. It supports the instrument; it is not a result.
>
> ⭐ **The niche is empty**, which is the thesis value: nobody has found a reference-absent copy from
> transcriptome data. The field standard is S1 re-assemble / S2 depth+PSV / S3 peptides. A *bounded
> negative* is therefore itself publishable.
>
> ⚠ **The advisor's "the reference is an average" objection fails on PROVENANCE, not on measurement:**
> mGorGor1 is a haplotype-resolved assembly of **one animal** and the fibroblast IsoSeq is that
> animal's **own cell line**. Keep the two arguments separate — the objection dies on how the substrate
> was built, and the copy-number-polymorphism rate is **still unmeasured** (the 2026-08-19 mat-vs-pat
> run was uninformative: control floor 0.1512 vs signal 0.0278). ⚠⚠ And **do not quote the 8/9/9-vs-5/6/8
> haplotype deficit** — it does not reproduce; see [`o3_haplotype_cnv_result.md`](o3_haplotype_cnv_result.md).
>
> **Current O3 avenues and claim boundary (2026-08-16):** see
> [`docs/o3_missing_copy_evidence.md` §8](o3_missing_copy_evidence.md#8-possible-o3-avenues--decision-record-2026-08-16).
> Liftoff or a second genome is not required for the main experiment; natural RNA-only findings remain
> candidates unless independently validated with donor DNA.
>
> **Current O1 purity rules and expanded known-family graphs (2026-08-16):** see
> [`docs/o1_false_positive_rules.md`](o1_false_positive_rules.md) and the 19-graph
> [`expanded audit`](../bench/o1_expanded_family_audit/README.md). Soto SD membership is discovery
> evidence, not automatic gene-family membership; primary and audit graphs are emitted separately.


> Thesis-level narrative + status. The machine-generated synthetic scorecard lives in
> `bench/multi_copy_eval/OBJECTIVES_ASSESSMENT.md` (regenerated by `run_oracle.py`); this file is the
> human-maintained framing it rolls up into. Last revised 2026-06-01.

**Thesis claim:** assemble the isoforms of *individual paralogous copies* of multi-copy gene families from
long reads — recovering, quantifying, and resolving copies that primary-alignment assembly (StringTie) cannot.
The variation graph (`--vg`) shares evidence across copies; an EM assigns ambiguous reads to copies.

Status key: ✅ done · 🔄 partial / in progress · ⛔ blocked.

---

## O1 — Formally define a multi-copy gene family in the VG  ✅
**Foundational** — makes O2–O5 measurable (you cannot score copy recovery without defining a copy).

A multi-copy family is a **connected component of the read-multi-mapping graph** (a *splice-graph
generalization*: nodes = exon-classes carrying per-copy sequence, paths = copy × isoform) satisfying:
**M** ≥2 copies · **H** copies share a backbone (reads multi-map) · **X** ≥2 copies independently expressed
(edit-distance anchored) · **I** graded identifiability (full / partial / none), reported not assumed.

- Spec: `docs/superpowers/specs/2026-06-01-multicopy-family-definition-design.md` (jargon scrubbed for the advisor).
- Operational predicate `classify_family` shipped (`src/rustle/vg.rs`), emitted as `--vg` GTF attributes
  (`family_verdict`, `family_identifiability`, `family_n_copies/n_expressed`, `family_locus_rel`).
  Commits `d286400`..`6d91bf1`. Expression counted per **identifiability class** (`3609b60`), runs
  independent of EM eligibility (`6d91bf1`).
- Validated on **14 known families** (`bench/paralog_secondary_scan/validate_known_families.py`): DAZ / NBPF /
  RBMY / amylase → `family`; TSPY (97% ties) / MAGEA → `family_nonidentifiable`; SORD+LOC → spillover;
  β-defensin / protocadherin / PRAME correctly excluded (not expressed here).
- **Fresh O1 purity challenge (2026-08-16):** current Rustle was rerun on predeclared regional
  extracts from the original whole-genome-aligned GGO/HSA BAMs, without old node or family ids as
  inputs. It re-emits 124/133 audited loci and 72/75 independently named targets. Although 14/16
  unrelated conflicting-gene loci are real/reproducible emissions, only 1/16 rejoins the target;
  all nine re-emitted non-NBPF loci from the adversarial repeat bridge remain outside fresh NBPF.
  GOLGA2 is now separated as a documented broad-family/recent-subfamily outgroup rather than
  mis-scored as an unrelated false positive; an RNA identity-0.80 view removes it but damages MAGEA
  and NBPF, so this is a typed hierarchy rather than a new global threshold.
  The cost is explicit: 69/75 named targets land in the modal family (three not emitted, three
  split). The HSA run also discloses one node-pair decision delegated to O2 `reads_distinguish`,
  so this particular node set is not sequence-only. Direct Rustle tables, rule certificates,
  logs, and actual fresh `E_r` GFAs:
  `bench/o1_fresh_emission_validation/`; interpretation: `docs/o1_false_positive_rules.md`.
- **Deferred implementation:** emit the current RNA homology family as the broad family plus an
  opt-in, nested DNA-supported recent-copy subfamily (`RECENT_COPY`, `BROAD_ONLY`, or
  `DNA_UNRESOLVED`). The annotation-free algorithm, flag-off byte-identity requirement, output
  schema, GOLGA discriminator, and cross-family safety tests are specified in
  `docs/o1_hierarchical_family_followup.md`. Production Rustle does not yet emit these fields.
- **Provenance-model avenue:** represent loci as ordered paths through homologous duplication blocks,
  with separate RNA-homology, DNA-duplication, read-conflict, and optionally rooted ancestry edges.
  This can express the mosaic GOLGA2 + ITSN2-UTR origin of the chr15 GOLGA expansion without calling
  ITSN2 a GOLGA family member. A single genome yields an unrooted network; directional “ancestral”
  claims require outgroup sequence/synteny, without using the outgroup as the assembly reference.
  Formal model, five-family local pairwise-witness prototype, empirical concordance results, and
  proof-of-concept criteria: `docs/o1_duplication_provenance_model.md`. Durable typed tables and GFA
  projections: `bench/o1_provenance_witness_prototype/`. This is evidence that the representation
  separates coherent cores from repeat bridges; stable multi-locus block-class construction is
  still deferred and the current graphs remain `UNROOTED`. A deferred single-outgroup extension now
  specifies direct minimap2 block/flank alignment to both phased gorilla haplotypes, two-sided
  synteny rooting, explicit abstention states, and flag-off invariance. It uses no annotation
  projection and cannot change human family membership. Local ape assembly paths and WSL mount
  instructions: `docs/linuxdisk_data_access.md`. A GOLGA proof of concept now finds recurrent human
  intervals from GOLGA2 into 8 audited family loci and from ITSN2 into 6, while both proposed source
  loci have unique two-sided synteny in both gorilla haplotypes. These are retained as
  `ROOT_CANDIDATE_SINGLE_OUTGROUP` with `direction_status=UNROOTED`, because stable multi-locus
  block classes are not implemented. Evidence and rerun script:
  `bench/o1_outgroup_rooting_poc/`.
- **Open:** the §9 "report the range of consistent splits" quantification is deferred (shipped behavior:
  label non-identifiable + abstain from splitting); scope gate is `any_spliced` (lenient vs per-copy wording).

## O2 — Recover copies that primary-only assembly misses  🔄 capability shown · ⛔ no validated real-data recovery
Copies whose annotated transcripts are covered overwhelmingly by **secondary** alignments (their primary went
to a sibling paralog) — StringTie drops them; the VG *can instantiate a model* for them. No end-to-end run yet
demonstrates a dropped copy actually recovered on real data.

- Genome-wide scan (`bench/paralog_secondary_scan/`, GGO IsoSeq vs RefSeq): of 846 expressed multi-copy
  candidates, an edit-distance copy-anchor **classifies** ~93 expressed_real_copy (NBPF, LRPAP1-like, …) vs
  ~89 pure_spillover loci that naive secondary-use would *fabricate* (the honest counterexamples;
  SORD/LETM1/CES1 pseudogenes). This is *classification*, not assembly — no run shows `--vg` emitting a
  expressed_real_copy locus that primary-only output drops.
- **DAZ3 is a RETRACTED false positive (corrected 2026-06-02), not a recovery showcase.** Of its ~158–164
  +strand reads, the overwhelming majority are *secondary* alignments of DAZ1's reads (median NM ~6 at DAZ1
  vs ~88 at DAZ3); only ~2–6 reads genuinely prefer DAZ3. The earlier "5 isoforms at cov ~113" was phantom
  mass echoed from DAZ1. The committed anchored-prior + joint-strand EM collapses it to **cov 4.04 (2 tx,
  `low_confidence`, identifiability `none`)** — DAZ3 is real-but-very-lowly-expressed and the tool correctly
  **abstains**. (The cov-113 result reproduces only with the fix disabled:
  `RUSTLE_VG_JOINT_STRAND_EM=0 RUSTLE_VG_ANCHOR_PRIOR=0` → 5 tx, top cov 112.77.) DAZ3's value is now a
  **correctness** example — the EM refuses to fabricate a copy from a sibling's reads — guarded by
  `daz3_isoforms_max` in the oracle.
- **Open:** (1) the 93/89 rate is unreproducible from committed state (stalled verification, inputs
  uncommitted); (2) the decisive recovery test — mask a known copy out of the reference, confirm `--vg`
  recovers + assembles it — is not done; (3) no validated real-data O2 recovery exists.

## O3 — Assign ambiguous reads to the right copy (the EM)  ✅ synthetic
Fingerprint-EM: latent = copy of origin, parameters = per-copy abundance, per-read likelihood from
match/mismatch at distinguishing positions → the raw **ΔNM** anchor (T=2 decisive; ΔNM=0 non-identifiable).
Ties are **apportioned by the prior, never fabricated** — the answer to "a single SNP makes this trivial."

- Synthetic oracle: read-to-copy accuracy **100% / 100% decisive** (Obj 4).
- Derivation figures: `bench/paralog_secondary_scan/em_derivation*.png` (worked nucleotide example →
  E/M steps → convergence on DAZ-like numbers → the tie/decisive boundary).
- **Open:** real-data per-read validation on GGO is thin.

## O4 — Assemble each copy's distinct isoforms / structural variants  ✅ synthetic (capability) · 🔄 real-data
Per-copy isoform recovery — each copy gets its own source→sink paths through the family splice graph.

- Synthetic oracle: per-copy isoform **Sn/Pr 100/100** (Obj 3) — but this fixture is spatially
  separable, so the non-VG baseline ALSO scores 100/100; it certifies non-regression, not a
  StringTie-beating capability.
- **Discriminating benchmark built (2026-06-03):** `bench/tandem_attribution/o4o5_copy_benchmark.sh`
  (spec `docs/superpowers/specs/2026-06-03-o4o5-copy-resolution-benchmark.md`). Key finding:
  **gffcompare chain Sn/Pr cannot discriminate** — when copies have full-length reads the baseline
  always ties or beats VG. The honest discriminator is **copy ATTRIBUTION** (which copy each
  ambiguous read came from): on a merged-but-separable fixture VG attributes multimappers to their
  true source copy at **100% (default settings, id 0.97, 5/5 seeds)** while the baseline produces
  **no copy metric at all** (undefined, not worse). VG **abstains 100% at identical copies** (the
  DAZ limit — fabricates nothing). Measured calibration gap surfaced: at the **1-SNP boundary VG is
  overconfident** (dec_acc 0.44 at dec_frac 0.75) — fix-and-remeasure item.
- **Open:** real-data breadth (per-copy isoform recovery + copy attribution on curated GGO families);
  the synthetic truth is read-name prefixes.

## O5 — Share evidence across copies via the graph  🔄 partial · ⛔ one blocker
Cross-family borrowing (coverage + junction propagation) and the *intended* **structural-linkage channel** —
copy-specific *junctions* as distinguishing positions, adding to ΔNM to resolve reads SNPs alone cannot
(figure B's DAZ1−/DAZ3+ illustration). **Status caveat:** the junction channel (`lambda_j`) is currently
**inert/redundant** in the EM (ablation is a no-op — see `project_color_cgroup_parity`), and DAZ3 is
near-silent (~2–6 genuine reads, see O2), so this is the design concept, not a demonstrated resolution.

- Coverage/junction borrowing implemented in the VG pipeline.
- **The "intra-bundle splitter" blocker was a mirage (corrected 2026-06-01).** GOLGA6L7 is **antisense-silent**
  (0 own-strand reads; its 53+17 reads are +strand, belonging to an overlapping +strand lncRNA) — it correctly
  emits nothing; the "miss" was a gffcompare overlap artifact, not a tool bug. Checked on 3 genuine same-strand
  expressed paralog pairs (gap 15 bp–3 kb): they fail by **bundle rejection / partial assembly in separate
  bundles, never merge-collapse** — so an intra-bundle splitter has **no validated target**.
- **🔄 Real open gap:** paralog bundles (low-cov, family-secondary-stripped) get **rejected/partially
  assembled** rather than each copy emitted — a depth-aware-bundling / family-read-handling problem, not a
  splitter. This is the genuine O5 work item.

---

## Cross-cutting honest gaps
- Most scores are **synthetic**; real-data per-copy validation (O3/O4 on GGO) is the weakest link.
- O2's **93/89 rate** is classification-only and unreproducible from committed state; no validated real-data
  recovery exists, and the former DAZ3 showcase is a retracted false positive (see O2).
- O5's **GOLGA6L7 splitter** is unsolved.
- Default de-novo (non-`--vg`) headline held at **95.6 / 90.5** throughout — none of the VG work regressed it.
