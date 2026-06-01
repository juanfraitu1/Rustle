# Design: Haplotype-phasing read-to-copy assignment (replacing HMM-EM)

Status: DESIGN (approved 2026-05-31). Next: implementation plan (writing-plans).

Grounding (this session): the multi-copy oracle + DAZ verification; the genome-wide family
discovery (`bench/multi_copy_eval/full_em_families.tsv`); the phasing analysis on fam 175 / fam 214 /
DAZ; the HMM-footprint workflow (`hmm-removal-and-blockers`); memory `project_vg_wiring` (⚠ DAZ3
false-positive + NO-valid-read-reassignment-win sections, and the RESOLVED Obj-4 fingerprint-EM fix,
commit `9436e7e`). The phasing scorer (`score_read_exon_fingerprint`) and the evidence-adaptive gate
were already repaired in `9436e7e`; this design supplies the *sites* it was missing and makes phasing
the default.

## 1. Problem & the decisive evidence

The thesis claim — VG assigns multi-mapping long reads to their true paralog copy better than StringTie
(which assigns 1/n for short reads and discards multi-mappers entirely for long reads) — is correct in
principle but **the current implementation does the opposite**:

- On **fam 175** (a real co-expressed paralog pair, NC_073234.2:62.95M + 67.15M), the multi-mappers carry
  a strong haplotype signal: **~90 divergent sites per read** (median net `|NM_A − NM_B|` = 90). The
  truth (phased) per-copy read split is **copy B > copy A** (B≈43, A≈20). StringTie and rustle-baseline
  both get the ratio right (B>A). **rustle `--vg` INVERTS it** (A≈22 > B≈14).
- Root cause, verified: `build_exon_fingerprints` reports **`n_sites = 0`** for fam 175 and labels the
  copies *"identical"* — because the homologous exons of the two copies are grouped into **separate
  ExonClasses** (the k-mer-Jaccard merge drops below threshold at 5–7% divergence), so they are never
  compared. With **zero** sites the phasing EM has nothing to phase with → it falls back to a pileup
  prior / the HMM forward-DP path artifact → **mis-assigns**. This same mechanism fabricated DAZ3 (a
  copy with no genuinely-originating reads) from DAZ1's reads.
- The HMM-EM's discriminating signal is NOT base identity (DAZ reads fit DAZ1 at NM 6 vs DAZ3 at NM 88,
  yet the forward DP prefers DAZ3 by ~7900 log-units — a path-length/normalization artifact). So the HMM
  is not merely redundant; it is **actively wrong** on the co-expressed case (fam 175 inversion) and the
  starved case (DAZ3 fabrication).

**Insight:** the phasing signal is present and strong; the only thing missing is a correct **diagnostic-
site finder**. Once reads can be scored by their alleles at the true divergent positions, assignment is
decisive, explainable ("assign each read to the copy whose alleles it carries"), and correct — and the
HMM (which mis-assigns) should be removed.

## 2. Goal & non-goals

- **Goal:** enumerate the real divergent sites between paralog copies and assign each multi-mapping read
  to the copy whose haplotype it matches; make this the **default `--vg`** read-assignment method;
  **remove the HMM-EM**.
- **Decision rule:** the default de-novo operating point must never regress (96.2/91.7 intron-chain,
  95.6/90.5 tx). All changes are `--vg`-only and the default de-novo path is structurally untouched.
- **Non-goals:** the Obj-1 intra-bundle splitter; Obj-2 absent-copy validation (only the *re-routing* of
  its scorer off the HMM is in scope); per-copy isoform enumeration; expression-calibration end-to-end.

## 3. Approach (chosen)

Site-finding by **exon-pair synteny alignment** (Approach 1 of 3 considered — see §9). Match homologous
exons across each copy-pair by genomic synteny + minimizer anchoring (bypassing the broken ExonClass
merge), compare them base-by-base, emit divergent sites in genomic coordinates with both copies' alleles.
Feed those into the already-repaired phasing scorer. Make the phasing-EM the default `--vg` path and
delete the HMM.

Rejected alternatives: (2) fixing the k-mer merge — band-aid, the merge stays fragile and exon-relative
comparison still breaks on indels; (3) full transcript pairwise alignment — most general but most new
code and concatenation can fabricate junctions.

## 4. Architecture & data flow

```
--vg + --genome-fasta
  ├─ ingest (bundle.rs): thread genome ALWAYS in vg_mode w/ FASTA (today gated on --vg-snp)
  │     → read.mismatches populated on the default path
  ├─ family discovery (vg.rs, unchanged) → FamilyGraph (per-copy paths, per_copy_spans, per_copy_sequences)
  ├─ enumerate_diagnostic_sites(fg, n_copies)   ← NEW (replaces broken inner loop of build_exon_fingerprints)
  │     per copy-pair: synteny-match exons → align → divergent sites (genomic coords, both alleles)
  │     → returns ExonFingerprints { sites, per_copy_site_refs, n_copies, n_sites }
  ├─ phasing EM (run_fingerprint_em, repaired scorer): score each placement by read alleles at sites
  │     reads with >=1 covered divergent site → decisive haplotype assignment (adaptive gate)
  │     reads with 0 covered sites → apportion by per-copy expression prior (M-step)
  │     → updates BundleRead.weight
  └─ per-copy assembly with phased weights (unchanged downstream)
```

## 5. Components

### 5.1 `enumerate_diagnostic_sites` (NEW — src/rustle/vg.rs)
Replaces the inner per-ExonClass loop of `build_exon_fingerprints`. Signature unchanged in spirit:
`fn enumerate_diagnostic_sites(fg: &FamilyGraph, n_copies: usize) -> ExonFingerprints`.

Algorithm:
1. For each copy `c`, get its ordered exon list `(genomic_span, sequence)` from `recover_paralog_path(c)`
   + `per_copy_spans` + `per_copy_sequences`. Sequences are on the transcript strand (per the
   `ExonClass.per_copy_sequences` doc).
2. For each unordered copy-pair `(i, j)`:
   a. Determine relative orientation from the copies' strands (FamilyGraph nodes carry `strand`). If
      opposite (inverted paralogs, e.g. DAZ1−/DAZ3+), reverse one copy's exon ORDER and reverse-complement
      its sequences so the two are co-linear before matching.
   b. Match homologous exons: start from positional order; refine/validate each pairing by minimizer
      overlap (`minimizers()` in family_graph.rs) so an inserted/deleted whole exon in one copy doesn't
      shift all downstream pairings. Unmatched exons (no syntenic partner) contribute no sites.
   c. For each matched exon pair `(seqA @ spanA, seqB @ spanB)`:
      - Equal length → gapless compare: for each offset `o` with `seqA[o] != seqB[o]` and both ACGT, it
        is a divergent site at exon offsets `(o, o)`.
      - Different length → minimizer-anchor + banded fill to co-linearise, then compare within aligned
        blocks; emit sites only at aligned (non-indel) columns with offsets `(oA, oB)`.
   d. For each divergent site, map exon offsets → **genomic positions** accounting for strand:
        `posA = genomic_pos_of(spanA, oA, strandA)`, `posB = genomic_pos_of(spanB, oB, strandB)`,
      where `genomic_pos_of` accounts for whether the sequence is forward or reverse on the genome (the
      scorer matches against `read.exons` in genomic coords). Push:
        `copy_bases = [(i, baseA), (j, baseB)]` (alleles on the genome-forward strand, so they compare to
        `read.mismatches` which are genome-forward);
        `per_copy_site_refs[i].push((site_idx, posA, baseA))`,
        `per_copy_site_refs[j].push((site_idx, posB, baseB))`.
3. `n_copies > 2`: union sites across all pairs (a position divergent in any pair is diagnostic);
   per-copy site lists are then sorted + de-duplicated by genomic position (the existing collapse logic
   from `9436e7e` already de-dups per-copy fragment multiplicity; extend it to dedup identical
   `(ref_pos)` entries so a site found in multiple pairs is counted once per copy).

CRITICAL coordinate invariants (the fiddly part — the plan must unit-test these):
- `per_copy_site_refs[c]` positions must be in the SAME coordinate frame as `BundleRead.exons` and
  `BundleRead.mismatches` for a read aligned to copy `c` (genome-forward, 0-based, matching `ref_pos`
  used in `score_read_exon_fingerprint`).
- `copy_bases` alleles must be genome-forward (so `read_base == copy_base` is a valid comparison against
  `read.mismatches` values, which are query bases reported genome-forward by `extract_mismatches_vs_fasta`).
- For an inverted copy, the reverse-complement used for *matching* must NOT leak into the stored alleles
  — store the genome-forward base at the genome-forward position.

### 5.2 Phasing scorer (`score_read_exon_fingerprint`) — UNCHANGED
Already repaired in `9436e7e` (collapse-per-copy representative, evidence-adaptive gate
`eff_gap = min(gap_threshold, per_site_gap * n_sites_covered)`). It now receives real sites and needs no
change. The `build_exon_fingerprints` site-dedup collapse is retained.

### 5.3 Genome threading on the default path (src/rustle/pipeline.rs)
Today `vg_snp_genome` (threaded into `detect_bundles_from_bam_with_snp`) is gated on
`config.vg_mode && config.vg_snp && config.genome_fasta.is_some()` (pipeline.rs ~10039). Broaden to
`config.vg_mode && config.genome_fasta.is_some()` so `read.mismatches` are populated whenever VG runs
with a FASTA — the phasing EM needs them on the default path. (Non-VG runs still pass `None` and pay
nothing.) Keep the `--vg-snp` flag as a no-op alias / deprecate.

### 5.4 Default routing (src/rustle/pipeline.rs)
The dispatch (`if do_hmm { run_pre_assembly_em_hmm } else if build_graph { run_fingerprint_em } else {...}`)
becomes: in VG mode with a FASTA + family graph, always `run_fingerprint_em` (the phasing EM). Remove the
`do_hmm` branch and the `fit_profiles_in_place` call. Multimap-sequence collection
(`config.vg_multimap_sequences`, today gated on `vg_solver==On && !vg_no_hmm`) is no longer needed by the
phasing EM (it uses `read.mismatches` + the FamilyGraph sequences), so its HMM-only gating is removed.

### 5.5 Ambiguous reads (0 covered divergent sites)
Reads whose span covers no divergent site (genuinely conserved region) cannot be phased per-read. They are
apportioned by the per-copy expression prior the *phased* reads establish (the EM M-step
`log_priors[k] = ln(copy_total[k]/total_sum + 1e-3)`), i.e. proportional to each copy's phased abundance.
This is the honest "expression apportionment" (no invented per-read signal). The pileup-fallback wiring
added this session (run the EM instead of skipping when `n_sites == 0`) is retained for the all-conserved
family case.

### 5.6 HMM removal + module de-conflation
From the footprint analysis (`hmm-removal-and-blockers`):
- DELETE (~2,157 LOC): `vg_hmm/profile.rs`, `vg_hmm/scorer.rs`, `fit_profiles_in_place`,
  `run_pre_assembly_em_hmm` (vg.rs), `run_em_reweighting_hmm` (vg.rs), the HMM-rescue scoring in
  `vg_hmm/rescue.rs`, `ExonClass.profile` / `per_copy_profiles` fields, and the HMM env knobs
  (`RUSTLE_VG_HMM_*`).
- MOVE (shared, NOT HMM, ~2,330 LOC) out of `vg_hmm/` into a `vg_graph` module: `family_graph.rs`
  (FamilyGraph/ExonClass/NodeIdx/build_family_graph/recover_paralog_path/minimizers/refine_by_minimizer_jaccard),
  `positional.rs` (genome scan), the bundle-synthesis portion of `rescue.rs`, and the `RescueClass`
  framework in `diagnostic.rs`. Update `crate::vg_hmm::` references in vg.rs, pipeline.rs, path_extract.rs,
  transcript_filter.rs.
- DELETE the 4 diagnostic binaries (`dump_family_graph`, `roundtrip_family_graph`, `loo_assembly`,
  `leave_one_out_rescue`) and the HMM tests (`vg_hmm_profile.rs`, `vg_hmm_scorer.rs`, `vg_hmm_em.rs`, the
  HMM cases in `vg_hmm_rescue.rs`/`vg_hmm_family_graph.rs`).
- RE-ROUTE Obj-2 rescue (`vg_hmm/rescue.rs::run_rescue_in_memory` uses `forward_against_family`): switch
  candidate-read scoring to the k-mer prefilter + phasing identity already present (k-mer hit fraction);
  preserve the rescue scaffold's behavior (it is unvalidated either way — only the scoring primitive
  changes).

Sequencing: land the phasing site-finder + default-routing + validation FIRST (prove correctness), THEN
the HMM deletion + module move (mechanical, gated on the same green validation). This keeps a working
default at every commit.

## 6. Validation gate (every item must pass; encoded in the plan + the multi-copy oracle)

| Check | Target | How |
|---|---|---|
| fam 175 per-copy ratio | **B > A** (today VG inverts to A>B) | region run, compare cov(B) > cov(A) |
| DAZ | DAZ3 NOT fabricated (cov → ~0); reads stay at DAZ1 | DAZ subset run |
| synthetic fixture (Obj-4) | **100%** decisive accuracy | `run_oracle.py --fast --check` |
| fam 214, fam 136 | both copies, ratio matches phased truth | region runs |
| diagnostic-site count | fam 175 finds ~tens of sites (was 0) | `RUSTLE_VG_FP_SITE_TRACE` |
| default de-novo GGO_19 | **unchanged** 96.2/91.7, 95.6/90.5 | gffcompare |
| genome-wide `--vg` | no crash; not net-worse than baseline | full GGO run sanity |
| HMM removal | builds clean; all above still green | post-removal re-run |

A new oracle check `score_phasing_assignment` (real-family per-copy ratio vs phased truth, fam 175/214 as
fixtures-by-coordinate) is added to `bench/multi_copy_eval/` and to `expectations.json`.

## 7. Risks & isolation

- **Default de-novo isolation:** all changes are `--vg`-only; the de-novo path never executes the site-
  finder, the phasing EM, or any moved/deleted symbol on a non-VG run. Verified at each task.
- **Coordinate/strand bugs in the site-finder** (highest-risk): mitigated by dedicated unit tests on a
  synthetic 2-copy + inverted-copy fixture asserting exact genomic `ref_pos` and genome-forward alleles,
  and by the fam 175 / DAZ end-to-end gate.
- **Structurally-divergent copies** (whole-exon indels): minimizer anchoring matches surviving exons;
  unmatched exons contribute no sites; reads there fall to expression-apportionment — never mis-phased.
- **HMM removal breakage:** sequenced after the phasing gate is green; the footprint analysis confirms no
  non-VG caller depends on the deleted symbols; `RescueClass` (used by path_extract/transcript_filter) is
  shared and preserved.

## 8. What "done" looks like

`rustle --vg --genome-fasta` assigns multi-mappers by haplotype: fam 175 → B>A, DAZ → no DAZ3, synthetic
→ 100%, default de-novo unchanged; the HMM is gone; `vg_graph` holds the shared infrastructure; the
multi-copy oracle has a real-family phasing check that fails on regression.

## 9. Approaches considered (for the record)

1. **Exon-pair synteny alignment (CHOSEN):** match homologous exons by synteny+minimizers, compare base-
   by-base, genomic-coord sites. Fixes the real failure, self-contained, handles inversion.
2. **Fix the k-mer ExonClass merge:** least code, but the merge stays fragile (over-merge on repeats) and
   exon-relative comparison still breaks on indels — band-aid.
3. **Full transcript pairwise alignment:** most general (handles rearrangement) but most new code (full
   banded aligner + column→genomic mapping) and concatenation can fabricate junctions.
