# Assembly-based parCN (`parcn`) — design

## Goal

Close the **parCN** (paralog-specific copy number) gap — the thesis's known hard limit that RNA cannot
cross — for the copies that ARE resolvable, using **data already on disk** (the two phased mGorGor1
haplotype assemblies) and **no new sequencing**. Ship it as an **optional, cheap, standalone tool** that
consumes one catalog artifact and reports per-copy genomic dosage split by haplotype.

## Boundary — this is a demarcated assembly/DNA-side supplement

The core Rustle pipeline stays **RNA-exclusive**. `parcn` is a SEPARATE optional tool that exists only
because a T2T genome + phased maternal/paternal haplotypes happen to be available for this species. It:
- reads only `copies.fa` (a catalog artifact) + the two phased haplotype assemblies,
- writes only its own `<out>.parcn*.tsv`,
- **never** modifies, wires into, or is invoked by the RNA catalog / assignment path.
Nothing about the RNA-exclusive core changes; removing `parcn` leaves the pipeline byte-identical.

## Motivation

- famCN already exists (chi_H lower bound, `depth_cn`, genome projection). parCN does not: `depth_cn` is an
  RNA expression ratio, not per-copy genomic dosage.
- A phased diploid assembly IS the copy-number ground truth: each real paralog appears once per haplotype
  per genomic copy. Counting a paralog's assembly loci gives its parCN directly — no DNA depth, no
  read-depth model, no GC/mappability mask.
- The blocker was never compute — it was an assembly-level paralog identity to count against. The catalog's
  per-copy consensus sequences provide that identity (each copy's PRIVATE positions = its SUNs); the
  assembly provides the loci. Intersect them → parCN.
- ~82% of copies are SUN-identifiable (`bench/sun_identifiability.py`), so parCN is recoverable for the
  large majority today; the ~18% Tier-3 collapse is flagged for DNA (out of scope).

## Architecture

A new binary **`parcn`** (`src/bin/parcn.rs`) + a new library module **`src/rustle/vg_family/parcn.rs`**.
Single responsibility; consumes existing artifacts; touches no other output. Optional by construction.

### Inputs
- `--copies-fa <catalog>.copies.fa` — per-copy exon-sum (spliced) consensus sequences. Header format
  (as written by `gw_family_catalog`): `>{family_id}|{copy_idx}|{chrom}:{start}-{end}|{strand}|nexon={n}`.
- `--mat <mat.splice.mmi>` and `--pat <pat.splice.mmi>` — the two haplotype assemblies of **the same
  species/individual as the catalog** (gorilla catalog → gorilla mGorGor1; slots are generic so a human
  phased assembly substitutes when the catalog is human). **Splice-mode** minimap2 indexes, pre-built ONCE
  from the on-disk haplotype FASTAs (`minimap2 -x splice -d mat.splice.mmi mGorGor1.mat.cur.*.fasta.gz`); a
  raw FASTA is also accepted (minimap2 builds the index in-line). Splice mode is required because the query
  is an exon-sum consensus mapping across genomic introns.
- `--out <prefix>` — writes `<prefix>.parcn.tsv` + `<prefix>.parcn_families.tsv`.
- `--minimap2 <path>` (default `minimap2`), `--threads <n>` (default 4).

## Components (each independently testable)

### 1. Private positions (SUNs) from the copy consensuses — `sun_positions(copies) -> Vec<CopySun>`
Input: a family's copies (each `{copy_id, seq}`). For each copy B, pairwise-align B to every sibling copy
(the in-repo banded Gotoh aligner `banded_msa_pair`, exon-sum vs exon-sum, so no introns — a wide-enough
band handles the small copy-length differences), and collect B's **private positions**: offsets in B where
B's base is present and differs from the aligned base of EVERY sibling. Threshold-free. Also compute a
pairwise-identity max vs siblings to set the tier:
- **Tier-1 (SUN):** ≥1 private position.
- **Tier-2 (hap-unique, no SUN):** no private position but max sibling identity < 0.999 (a unique
  combination, no single pinning base).
- **Tier-3 (collapsed):** max sibling identity ≥ 0.999 (indistinguishable from another copy).
Returns per copy: `{copy_id, tier, private: Vec<(pos_in_copy, base)>}`.

### 2. Per-haplotype projection retaining the alignment — `project_with_cigar(...)`
For each haplotype (mat, then pat — one minimap2 index load each), project **all copies** (splice mode,
`min_identity 0.95`, `min_cov 0.90`, `known = {}` so NO locus is excluded — we WANT every genomic copy).
Unlike `project_families_batch` (which returns coordinate-only `CopyLocus`), this variant **retains the PAF
CIGAR** for each hit so the assembly base at any consensus position can be read. Group hits by family; union
per-copy hits and **dedup overlapping loci** (reciprocal-overlap ≥ 50%), keeping the highest-identity hit
(its `copy_id` = the locus's best copy, its CIGAR = the transfer alignment). Output per family per
haplotype: `Vec<Locus{ chrom, start, end, best_copy, identity, cigar }>`.

### 3. Read assembly bases at private positions — `assembly_bases_at(cigar, ref_start, positions) -> Vec<Option<u8>>`
Walk the projection CIGAR of a locus's best-copy hit (query = copy consensus, target = assembly). For each
requested consensus position, return the assembly base aligned to it (`None` if it falls in a
deletion/intron gap). Single hop: consensus-position → assembly base, directly from the alignment that
placed the locus.

### 4. Hybrid assignment — `assign_locus(locus, sun_positions, copies) -> Assignment`
- **Tier-1 (SUN, deterministic):** the locus's best copy is B; read the assembly bases at B's private
  positions (component 3). If the assembly carries B's private base at ≥1 private position → assign the
  locus to B, `method = SUN`.
- **Tier-2 (align fallback, flagged):** B has no private position (Tier-2/3 copy) OR the assembly did not
  carry B's private base → assign to the highest-identity copy at that locus, `method = align_fallback`,
  provided it beats the runner-up identity by ≥ 0.002; else →
- **Tier-3 (unresolved):** near-tie / no discriminator → `method = UNRESOLVED`, counted in the family's
  `n_unresolved_loci`, never assigned to a copy.

### 5. Tabulation + output — `write_parcn(...)`
`parCN[copy]` = # loci assigned to it; split `loci_mat` / `loci_pat`.

`<out>.parcn.tsv`
```
family_id  copy_id  sun_tier  loci_mat  loci_pat  parCN  assign_method
RBMY1A1    copy0    T1        1         1         2      SUN
RBMY1A1    copy3    T2        1         0         1      align_fallback
```
`<out>.parcn_families.tsv` (roll-up)
```
family_id  n_copies  famCN_diploid  n_unresolved_loci
RBMY1A1    6         11             1
```
`famCN_diploid = Σ parCN` over the family's copies.

## Data flow
```
copies.fa ─► parse (family, copy_id, seq)
              │
              ├─ sun_positions(copies) ─────────────────────────────┐  (per-copy private positions + tier)
              │                                                      │
              ├─ project_with_cigar(mat) ─► dedup loci(mat) ─► assign_locus ─┤
              └─ project_with_cigar(pat) ─► dedup loci(pat) ─► assign_locus ─┤
                                                                             ▼
                                                          write_parcn ─► parcn.tsv + parcn_families.tsv
```

## Error handling / graceful degradation
- minimap2 missing or non-zero exit on a haplotype → that haplotype contributes 0 loci (match
  `genome_projection`'s existing `Ok(None)` graceful-degradation contract), a WARNING is logged, the run
  continues (parCN then reflects only the haplotype that projected).
- A single-copy family (one record in `copies.fa`, no siblings) → no private positions to compute; its loci
  are reported with `sun_tier = NA`, `method = single_copy`, `parCN` = raw deduped locus count.
- A locus whose best-copy CIGAR cannot be parsed / is absent → dropped with a per-locus WARNING (counted in
  a `n_dropped` diagnostic printed to stderr, not in `n_unresolved`).

## Testing
1. **Unit — `sun_positions`:** planted copies (one with a private SNV, one identical to another) → assert
   the private position is found for the SNV copy (Tier-1), Tier-3 for the identical pair, correct tier for
   a unique-combination copy.
2. **Unit — `assembly_bases_at`:** a hand-built CIGAR with a deletion; assert the base at a pre-deletion
   consensus position is read correctly and a position inside the deletion returns `None`.
3. **Unit — `assign_locus`:** synthetic loci → Tier-1 SUN match (assembly carries the private base), Tier-2
   fallback with the ≥0.002 identity margin, and a near-tie → UNRESOLVED, each asserting the `method`.
4. **Unit — dedup:** overlapping per-copy loci (reciprocal-overlap ≥ 50%) collapse to one, keeping the
   highest-identity hit as `best_copy`.
5. **Integration — synthetic 2-haplotype:** two tiny FASTA "haplotypes", each carrying planted paralog
   copies with known private bases; run the binary end-to-end; assert `parcn.tsv` gives the planted per-copy
   `loci_mat`/`loci_pat` and the `famCN_diploid` roll-up.
6. **Real conservation + sanity:** on the **gorilla (GGO) catalog projected onto the gorilla mGorGor1
   haplotypes** (species-consistent) — assert `Σ parCN + n_unresolved_loci` = total distinct projected loci
   (conservation), and that a CN-stable gorilla family's `famCN_diploid ≈ 2×` its haploid catalog count
   (e.g. RBMY 6 haploid → ~12 diploid; DAZ, GSTM as secondary checks; report the actual number, don't
   hard-gate — real CN can deviate from exactly 2×).

## Non-goals (YAGNI)
- NOT resolving the Tier-3 collapsed residual — flagged as `n_unresolved_loci`, routed to DNA (out of scope).
- NOT changing famCN in the existing legs (`famcn.tsv` / `famcn_readonly.tsv` / `collapsed.tsv` unchanged).
- NOT wiring into or altering the RNA-exclusive core pipeline (see Boundary).
- NOT ingesting DNA reads or building a read-depth / GC / mappability model — that is the other route,
  deliberately separate.

## Success criteria
- A single `parcn` run on `copies.fa` + the two on-disk haplotypes emits per-copy `parCN` (mat/pat split)
  for the SUN/hap-resolvable copies, with the Tier-3 residual flagged, in minutes and within the box's RAM
  (two splice-index loads, one at a time).
- Deterministic SUN assignments carry a `SUN` method tag; fallback assignments are flagged `align_fallback`;
  nothing collapsed is silently assigned.
- Validation conservation holds and a CN-stable family's diploid famCN tracks ~2× its haploid annotation.
- The RNA-exclusive core is untouched (no new call sites, no changed output).
