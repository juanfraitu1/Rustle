# Haplotype-resolved isoform assembly (`--phase-isoforms`)

**Date:** 2026-06-07
**Status:** Design approved, pending implementation plan
**Axis:** A new dimension where rustle beats StringTie for long reads. StringTie is
*sequence-blind* — it uses only junctions and coverage, so it physically cannot tell
which haplotype/allele an isoform came from. Long reads span multiple heterozygous
SNVs; that phasing signal is exactly what ST discards. This feature outputs
haplotype-resolved isoform sets plus allele-specific isoform usage — a capability ST
structurally lacks.

## Goal

Given long reads (optionally aligned across a paralog family), produce:

1. Transcripts tagged with their haplotype (HP1/HP2) and phase set (PS).
2. An allele-specific isoform-usage table (`<output>.ase.tsv`) — per isoform, the
   read support split by haplotype, with flags for haplotype-specific isoforms and
   usage imbalance. **This is the headline scientific deliverable.**

All of this from rustle's own read sequence — no external phasing tool required —
while still honoring external HP tags when the BAM is pre-phased.

## Non-goals / boundaries

- **Default behavior is unchanged.** Feature is opt-in (`--phase-isoforms`). With the
  flag off, output is byte-identical to today (the invariant every prior feature in
  this codebase maintains).
- Not a general variant caller. Het detection exists only to drive phasing; we do not
  emit a VCF.
- Not polyploid-general in v1. We assume diploid (two haplotypes per copy). The MEC-DP
  is written for a bipartition; >2 haplotypes is a future extension.
- We do not attempt to phase across phase-set boundaries (regions never co-covered by a
  single read). Those become separate phase sets with locally-arbitrary HP labels.

## Key design insight: order the two axes

A read in a paralog-family bundle varies for two independent reasons:

- **Which copy** it came from — Paralog-Sequence-Variants (PSVs), i.e. positions fixed-
  different between copies.
- **Which haplotype** within that copy — germline heterozygous variants (hets).

These axes are orthogonal, so we resolve them in sequence rather than jointly:

1. **Copy assignment first** (reuse existing VG machinery: `compute_copy_ownership`,
   fingerprint-EM, decisive-evidence). This strips PSV variance.
2. **Per-copy MEC phasing second** (new). Within a single copy's read set the only
   remaining structured variance *is* het variance, so MEC-DP is clean and the
   het-vs-PSV discriminant problem disappears — it is solved by *ordering the steps*,
   not by a classifier.

A single-copy (non-family) locus is simply the trivial one-copy case, so single-copy and
multi-copy loci unify under one mechanism.

### Honest caveat — confidence inheritance

Per-copy phasing inherits the copy assignment's confidence. For decisively-resolved
copies it is sound; for thin / low-identifiability copies (the bounded regime documented
throughout this project's history) it would manufacture fake haplotypes from noise.
Therefore **per-copy phasing is gated on the copy's evidence certificate**: only copies
that clear a confidence threshold are phased; thin copies stay merged-unphased (HP=None)
rather than emitting spurious haplotypes. This is a deliberate detect-and-degrade design,
never fabricate.

## Architecture

Four units, sequenced. B0 is mostly wiring; B1 → B2 → A are the implementation order
with review checkpoints between.

### B0 — Copy partition (wiring)

- **Input:** a bundle with reads carrying `seq` / `mismatches` (requires the
  `--vg-snp`-style sequence+mismatch ingest path to be active) and, for family loci, the
  VG copy assignment already computed.
- **Output:** the bundle's reads partitioned into per-copy read sets, each annotated with
  the copy's evidence certificate (confidence). Non-family bundle → one trivial copy
  containing all reads.
- **Notes:** reuses existing `compute_copy_ownership` output. No new copy-assignment
  logic. Composes with `--vg`.

### B1 — Het detection (per copy read set)

- **Input:** one copy's read set; each read has `mismatches: Vec<(ref_pos, alt_base)>`
  (vs reference) and an aligned span.
- **Algorithm:** build a mismatch pileup. A read that matches reference at a site
  contributes no mismatch, so coverage must come from read spans:
  `MAF(pos) = (#reads with an alt allele at pos) / (#reads spanning pos)`.
  A candidate het is:
  - biallelic (one dominant alt allele),
  - `MAF ∈ [MIN_MAF, MAX_MAF]` (default ~0.25–0.75 — excludes errors near 0 and
    fixed/reference-artifact differences near 1.0),
  - ≥ `MIN_ALLELE_READS` per allele (default 3),
  - strand-balanced (both alleles seen on both strands, to drop strand-biased artifacts),
  - not inside a homopolymer / low-complexity window (the indel-artifact filter).
- **Output:** an ordered list of het sites for the copy, each with its two alleles.
- **Determinism:** sites sorted by reference position; allele tie-breaks by base.

### B2 — MEC-DP phasing (per copy read set)

- **Input:** the copy's reads + its het-site list from B1.
- **Algorithm:** whatshap-style Minimum-Error-Correction dynamic program for two
  haplotypes.
  - DP proceeds column-by-column over het sites; the DP state is the bipartition of the
    reads *active* (covering the current column) into the two haplotypes; cost is the
    number of allele observations that must be corrected.
  - **Coverage cap:** to bound the state space, active coverage per column is capped
    (default ~15) by deterministic read selection — longest / highest-quality reads
    first, ties broken by a stable key (read id / position). Dropped reads are phased by
    best-agreement assignment after the DP, not used to build the backbone.
  - Reads spanning 0 het sites → unphased (HP=None).
- **Phase sets:** het sites are linked only when co-covered by at least one read;
  connected components define phase sets. Each component gets a `ps_tag`.
- **Canonical labeling:** within each phase set, HP1 is the haplotype carrying the
  reference (or lexicographically lower) allele at the component's first het site, so
  labels are stable across runs.
- **Output:** `hp_tag` (1/2/None) and `ps_tag` written back onto each read.
- **External-tag precedence:** if ≥ `EXT_HP_FRAC` of the copy's reads already carry a
  BAM `hp_tag`, skip internal phasing for that copy and use the external tags as-is.
  One precedence rule, no blending.
- **Determinism:** fixed read-selection order, fixed column order, deterministic tie
  resolution in the DP traceback.

### A — Per-(copy, haplotype) assembly + outputs

- **Input:** the bundle with per-read `hp_tag`/`ps_tag` set.
- **Assembly:** for each (copy, haplotype) partition, run the existing assembly path via
  `split_bundle_by_phase` (extended so it splits per copy when a copy partition is
  present; unphased reads are duplicated into both haplotype sub-bundles, as the current
  prototype already does). Each resulting transcript is tagged with its HP and PS.
- **GTF output:** transcripts carry `haplotype "1"|"2"|"."` and `phase_set "<id>"`
  attributes.
- **ASE table (`<output>.ase.tsv`):** one row per emitted isoform —
  `chrom, locus, copy, transcript_id, phase_set, reads_hp1, reads_hp2,
  allele_specific (isoform observed on a single haplotype only), imbalance
  (|hp1-hp2| / total)`. Rows sorted by `(chrom, start, copy, phase_set,
  transcript_id)` for deterministic output.

## Data flow

```
BAM (seq + mismatches + optional HP/PS)
   │  --vg [family copy assignment]  --phase-isoforms
   ▼
B0 copy partition ──► per-copy read sets (+ confidence certificate)
   │  (gate: copy confidence ≥ threshold, else leave unphased)
   ▼
B1 het detection ──► het sites per copy
   ▼
B2 MEC-DP phasing ──► hp_tag / ps_tag per read   (or external HP tags if present)
   ▼
A  split_bundle_by_phase (per copy) ──► assemble each (copy, haplotype)
   ▼
HP/PS-tagged GTF   +   <output>.ase.tsv
```

## Error handling / edge cases

- **No hets in a copy:** copy emits one unphased transcript set (HP=None). Not an error.
- **Low-confidence copy:** skipped by the confidence gate; reads stay merged-unphased.
- **Single phase set vs multiple:** multiple phase sets get independent HP labelings;
  the ASE table is per-(transcript, phase_set) so cross-block labels are never conflated.
- **Reads dropped by the coverage cap:** still HP-assigned post-DP by best agreement;
  never silently lost.
- **`--phase-isoforms` without seq/mismatch ingest:** hard error with a clear message
  (the feature needs sequence data).

## Testing

TDD per unit (RED → GREEN), matching the repo convention.

- **B1 het detection:** unit tests — a planted het at MAF 0.5 is found; a 5% error site
  is rejected; a fixed (MAF ~1.0) difference is rejected; a homopolymer artifact is
  rejected; strand-biased site rejected.
- **B2 MEC-DP:** unit tests — two clean haplotypes recovered exactly; one corrected
  observation tolerated (MEC = 1); reads spanning no hets → unphased; coverage cap
  triggers deterministic selection; canonical labeling stable; two un-co-covered het
  groups → two phase sets.
- **A outputs:** unit tests — HP/PS attributes present in GTF; ASE table flags a
  haplotype-specific isoform; rows deterministically ordered.
- **Integration:** a synthetic phased bundle end-to-end → two haplotype isoform sets +
  ASE table.
- **Invariant:** with `--phase-isoforms` off, full output byte-identical to baseline;
  full lib suite stays green.

## Validation (beating StringTie, on data)

- **Synthetic:** badread reads (pacbio2021 HiFi profile, per prior sim-reads work, with
  the −strand re-orientation caveat) from a diploid locus with K planted SNVs and
  *differential isoform usage* (HP1 → isoform X, HP2 → isoform Y, plus a shared isoform
  Z). Success = recovery of both haplotypes and the haplotype-specific isoforms.
  External oracle = `whatshap haplotag` + per-HP assembly.
- **Real:** a het gene with known allele-specific expression. Report rustle-internal vs
  whatshap read-partition concordance (% of reads in the same bipartition), and contrast
  against StringTie's single blended isoform set (no HP, no ASE) — the qualitative
  "ST can't do this" comparison.

## Flags / config

- `--phase-isoforms` — master switch (default off).
- Tunables (env-overridable, sensible defaults): `MIN_MAF`, `MAX_MAF`,
  `MIN_ALLELE_READS`, coverage cap, copy-confidence gate threshold, `EXT_HP_FRAC`.
- Requires sequence + mismatch ingest (the `--vg-snp` path); composes with `--vg`.

## Open implementation questions (resolve during planning)

- Exact reuse surface of `compute_copy_ownership` for B0 (what it returns vs what the
  partitioner needs).
- Where `split_bundle_by_phase` hooks into the assembly loop, and how per-copy splitting
  layers on top of its current whole-bundle split.
- Default value of the copy-confidence gate (tie to the existing evidence-certificate
  thresholds rather than a new constant).
