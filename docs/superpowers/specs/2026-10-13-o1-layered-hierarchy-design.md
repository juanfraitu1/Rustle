# O1 as a layered hierarchy: sedef-rs SDs and duplicons, a Soto-style layer, our families — design (2026-10-13)

## Intent (agreed with the user)

**The user's words:**
- "do we now have the ability to run sedef locally with the sedef port? Could we include that in O1?"
- "the Soto replication pipeline should also be part of O1; ideally I want to include them as layers of it, but I think
  to do that I need to achieve more integration."

**Decisions** (each chosen by the user):
1. **Purpose: both.**
   - One pipeline computes every layer for any genome: SDs and duplicons from sedef-rs, Soto-style families, our MCL
     families and expansions.
   - Each layer is validated against its own reference.
   - O1's claims are judged against independent truth.
2. **Copy number:** assembly CN everywhere, so the method is the same for human and every ape. For human, read-depth CN
   is reported alongside to show what the switch costs.
3. **Soto layer:**
   - A Rust core, the SD98 exon-sharing edges plus copy-number clustering, lives inside the hierarchy module.
   - The existing Python replication is the test oracle.
4. **Architecture:** an `sd` stage, run once per genome and cached. The hierarchy flag computes every layer internally
   from that cache.

**Consequence, accepted.** Once Soto's method is a layer of ours, Soto's published families are no longer independent
truth for that layer. They give a replication check (S2), not a validation of O1.

**This supersedes** phases V and B of `docs/superpowers/specs/2026-10-13-close-o1-design.md` (amendment 9 there). That
spec's amendment-7 prerequisites are folded into section 3.

## What exists (verified 2026-10-13)

- **sedef-rs** (`/mnt/c/Users/jfris/Desktop/sd/sedef-rs`; README status 10-07):
  - ~2,700 lines of pure Rust in four stages: seed → chain → refine → decompose.
  - Outputs bp-exact SD pairs with strand and identity, elementary segments, and duplication families. 12 unit tests;
    deterministic.
  - Measured: chimp chr1 in 55 s at ≤ 3 GB. Whole-genome is extrapolated at ~12 min and was **never run**.
  - **Never compared with SEDEF.** Identities are k-mer lower bounds (no DP alignment); chaining is greedy;
    single-genome mode only.
  - Its release binary is not built here (the target dir is empty).
- **Soto replication** (`bench/soto/soto_replication.py`, ~2,700 lines of Python): `genesets → edges → cluster →
  score`.
  - With Soto's own inputs it reproduces the frozen `winloci_data/soto_replication/
    replicated_families_2334_median_finalhuman.tsv` (ARI 0.6959). Those inputs are SEDEF `final_human_clean.bed`,
    CAT v4 genes and S1C famCN, with the median MAD statistic.
  - Its header records that the famCN leg is circular against Soto.
  - Measured copy-number ladder (ARI / exact):

    | Copy number used | ARI | exact |
    |---|---|---|
    | sequence only | 0.7307 | 345 |
    | 268-sample SGDP read-depth (`famcn_ours_allwssd.tsv`) | 0.8855 | 375 |
    | S1C | 0.9698 | 479 |

    Copy number dominates the concordance.
- **Assembly copy number** (`bench/soto/parcn_assembly.py`; `docs/archive/2026-09/PREREG_soto_parcn_assembly_2026-09-29.md`,
  register 1171-1174):
  - It applies QuicK-mer2's own k-mer rule, as Soto ran it, counted exactly in complete assemblies: canonical k = 30,
    with the famCN analogue over the region's unmasked 30-mers.
- **Hierarchy flag** (`mcl_families --emit-hierarchy`): it reads a DupMasker colors BED (`--hierarchy-duplicons`) and
  a SEDEF pairs BED (`--hierarchy-sd`), and only *places* external families (`--hierarchy-external`).

## 1. Components and data flow

- **`sd` stage (new)**, once per genome, cached by the genome FASTA's content hash.
  - sedef-rs is vendored into Rustle as `src/sedef/`.
  - Run as `utilities sd-scan --genome FA --out PREFIX [-t N]`.
  - It writes `PREFIX.sd.pairs.tsv`, `PREFIX.sd.segments.tsv` and `PREFIX.sd.families.tsv`.
  - The Desktop repo stays the user's sandbox; changes are ported deliberately.
- **Hierarchy:** `mcl_families --emit-hierarchy --sd-cache PREFIX` computes the layers internally:
  - **SD blocks and complexes** from the sd pairs, replacing `--hierarchy-sd`.
  - **Duplicons** from the sd segments and families, replacing `--hierarchy-duplicons`. Expansions keep today's rule: a
    family's maximal duplicon signatures.
  - **Soto-style families** from the Rust core.
    - Edges join member loci whose exons share sequence through SD98 pairs (sd pairs at identity ≥ 0.98; the
      replication's exon rule).
    - They are clustered by assembly copy number with the replication's MAD statistic.
    - Assembly copy number is `parcn_assembly.py`'s famCN rule applied to each locus's exons and counted in the same
      genome. Its masking is the genome's own soft-masking; a FASTA without soft-masking is refused for CN, and that
      is reported.
  - **Our MCL families and superfamilies**, unchanged.
- The existing BED inputs and `--hierarchy-external` stay, for comparison runs and as fallbacks. The families stage
  does not change.
- **Outputs:** the hierarchy tables gain a `soto` level, with a `members.tsv` column and its own node rows.

## 2. Validation gates (each pre-registered before it runs)

- **V1 — SD calls, sedef-rs vs SEDEF.**
  - Whole-genome runs on CHM13 v2.0 and on gorilla (`GGO.fasta`).
  - References:
    - CHM13: `/mnt/linuxdisk/tmp/o1_close/lift/sedef_pairs.v2.bed` (`final_v1_clean.bed` lifted by the T2T chain;
      md5 7bd521edff9bd1e27b8141027da37a01);
    - gorilla: `winloci_data/GGO_sedef_final.bed`.
  - Compared base-pair-wise inside SEDEF's range (≥ 1 kb, ≥ 90% identity). Both sides of every pair are merged into
    SD regions per tool:
    - recall = SEDEF SD bases covered by sedef-rs SD regions / SEDEF SD bases;
    - precision = sedef-rs SD bases covered by SEDEF SD regions / sedef-rs SD bases.
  - Bar: recall and precision ≥ 0.80 on both genomes. Run time and peak memory are recorded.
  - **Fail:** sedef-rs is not adopted. The SD layer keeps external SEDEF where it exists, and the apes get none.
- **V1b — SD98.** The same comparison restricted to identity ≥ 0.98, since the Soto layer cuts there and sedef-rs
  identities are lower bounds.
  - Bar: ≥ 0.80.
  - **Fail:** the Soto layer takes SD98 pairs from SEDEF where it exists, and alignment-based identity refinement in
    sedef-rs becomes its own task.
- **V2 — duplicon layer (human only; DupMasker exists only for CHM13).**
  - Expansions built from sedef-rs duplicons agree with the DupMasker-built ones at ≥ 0.80 on a held-out family-hash
    half. Agreement is the bipartite measure of register 1263.
  - Soto families equal to an expansion strictly below a family (`expansion_not_family`) are no more than 10% below
    DupMasker's count.
  - **Fail:** sedef-rs gives the SD level only. Human expansions stay on DupMasker, and the apes have no expansion
    layer.
- **S1 — Soto core vs the Python oracle.** With Soto's own inputs, the Rust core reproduces
  `replicated_families_2334_median_finalhuman.tsv` as an identical partition (ARI 0.6959).
  - **Fail:** the Soto layer stops until the core matches.
- **CN1 — assembly copy number vs its oracle.** On the same regions and genome, the Rust famCN equals
  `parcn_assembly.py`'s values exactly. That includes resolved and unresolved status.
- **S2 — the Soto layer on our inputs (descriptive, no bar).**
  - Human A119b and testis: sedef-rs SD98, our loci and assembly CN.
  - Reported: concordance with Soto's published families (ARI, exact, matching over distinct member sets), next to
    S1's 0.6959.
  - Also reported: the same run with the 268-sample read-depth CN and with S1C, which shows what assembly CN costs.
  - No bar, because our loci are not Soto's gene set.

## 3. O1's claims against independent truth (each claim pre-registered; bars relative to nulls)

- **Truths:**
  - **famsim:** simulated families with known answers, re-run as a regression suite on every layer change.
  - **Liftoff copy pairs (every species):** gene copies Liftoff projects onto each assembly, independent of our edge
    rule; figure 8's framework.
  - **Ensembl Compara paralog groups (human):** reported, not barred, since they include far older duplications than
    our SD-era families.
  - **Yoo 2025 lineage-specific expansions (apes;** `winloci_data/yoo_2025`**):** scored at the expansion layer.
- **Scorer:** `bench/o1_hier/hier_match.py`, extended as close-O1 amendment 7 requires:
  - matching over distinct member sets;
  - the within-family null for any claim below the family level;
  - cross-catalog comparisons stratified by mapped set size.
- **Bars:**
  - **Family layer:** in each species, Liftoff copy pairs fall in the same family more often than under a
    size-matched null (p < 0.001), including on the held-out species.
  - **Expansion layer:** Yoo expansions recovered strictly below a family beat the within-family null (p < 0.001).
- **Held-out design.** Develop on human and chimp; hold gorilla out. Orangutan joins once its families stage has run.
  This doubles as the generalization pre-registration that was offered but never written.
- **Soto.** No longer truth for anything. Concordance with Soto's published families is reported for the Soto layer
  (S2) and for the MCL layer, labelled not independent.
- **Structure.** The layers must nest as stated (SD complex ⊇ superfamily ⊇ family ⊇ expansion). Every violation is
  counted and listed, with no tolerance.

## 4. Order, gates and scope

1. **Vendor sedef-rs** (`src/sedef/`, `utilities sd-scan`, its 12 tests ported). Run V1 and V1b on CHM13 and gorilla.
2. **The `sd` stage and the hierarchy reading its cache** for the SD and duplicon layers. Run V2 on human.
3. **The Soto core in Rust** (S1). Then assembly copy number (CN1), then S2 on human A119b and testis.
4. **The `soto` level in the hierarchy,** and the nesting check.
5. **O1 evaluation against the independent truths** (section 3).

Each phase gets its own implementation plan, written after the previous phase's results; the first plan covers steps
1-2. The gates above decide what each later phase builds on.

**Out of scope:**
- O2, O3, and the HG002 arm (waiting on the cluster BAM);
- duplicon IDs consistent across species (each genome gets its own sedef-rs families);
- alignment-based identity in sedef-rs (a separate task if V1b fails);
- population copy numbers for the apes.

## Deliverables

- `src/sedef/` and `utilities sd-scan`; the cached `sd` stage.
- SD, duplicon and `soto` layers computed inside `mcl_families --emit-hierarchy`.
- The Rust Soto core and assembly famCN, with oracle tests (S1, CN1).
- Pre-registrations and results for V1/V1b, V2, S1, CN1, S2 and the section-3 claims; register rows.
- `docs/THESIS_OBJECTIVES.md`'s O1 line updated to what the bars support.

## Results of phases 1-2 (2026-10-14)

V1 FAIL, V1b FAIL, V2 FAIL: `docs/PREREG_sedef_rs_layers_2026-10-13.md` (register 1278-1280). The
readings fixed there decide what phases 3-4 build on: the SD layer keeps SEDEF where it exists (V1); the Soto layer
takes SEDEF SD98 (V1b); human expansions stay on DupMasker and the apes have no expansion layer (V2). The readings
apply per layer (prereg, "Clarifications registered 2026-10-14"), so V2's fail does not decide the SD layer's source.

Consequences for section 3. With V2 failing, the apes have no expansion layer, so section 3's Yoo expansion bar (scored at
the expansion layer on apes) cannot be scored until an expansion source passes. With V1 and V1b failing, a genome without
SEDEF has no SD layer and no SD98 source. On disk only human (HSA) and gorilla (GGO) SEDEF tables exist, so chimp
(section 3's dev ape) has neither; the current scope is human and gorilla (user, 2026-10-13).
