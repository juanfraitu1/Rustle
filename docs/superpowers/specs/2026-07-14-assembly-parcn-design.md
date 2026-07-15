# Assembly-based parCN (`parcn`) — design

## Goal

Close the **parCN** (paralog-specific copy number) gap — the thesis's known hard limit that RNA cannot
cross — for the copies that ARE resolvable, using **data already on disk** (the two phased mGorGor1
haplotype assemblies) and **no new sequencing**. Ship it as an **optional, cheap, standalone tool** that
consumes the existing catalog and reports per-copy genomic dosage split by haplotype.

## Motivation

- famCN already exists (chi_H lower bound, `depth_cn`, genome projection). parCN does not: `depth_cn` is an
  RNA expression ratio, not per-copy genomic dosage.
- A phased diploid assembly IS the copy-number ground truth: each real paralog appears once per haplotype
  per genomic copy. Counting a paralog's assembly loci gives its parCN directly — no DNA depth, no
  read-depth model, no GC/mappability mask.
- The blocker was never compute (projection already runs on the box) — it was having an assembly-level
  paralog identity to count against. The RNA catalog provides that identity (each copy's SUN signature);
  the assembly provides the loci. Intersect them → parCN.
- ~82% of copies are SUN-identifiable (`bench/sun_identifiability.py`), so parCN is recoverable for the
  large majority today; the ~18% Tier-3 collapse is flagged for DNA (out of scope).

## Architecture

A new binary **`parcn`** + a new library module **`src/rustle/vg_family/parcn.rs`**. Single responsibility,
consumes existing artifacts, does not touch `gw_family_catalog` / `copy_assign` output. Optional by
construction: nothing runs it unless invoked.

### Inputs (all already produced by the pipeline, or on disk)
- `--copies-fa <catalog>.copies.fa` — per-copy exon-sum (spliced) consensus sequences (query for projection).
- `--psv-copies <copy_assign --dump-psv>.psv_copies.tsv` — `family_id, copy_index, copy_tid, alleles`
  (the per-copy PSV hap-vector; `alleles` is the copy's base at each PSV column, `.`=missing).
- `--psv-cols <…>.psv_cols.tsv` — `family_id, col_index, genome_pos` (PSV column → RNA-reference position;
  used to order columns and for provenance; assembly positions are derived by alignment transfer, below).
- `--mat <mGorGor1.mat.splice.mmi>` and `--pat <mGorGor1.pat.splice.mmi>` — the two haplotype assemblies of
  **the same species/individual as the catalog** (gorilla catalog → gorilla mGorGor1; the two slots are
  generic, so a human phased assembly could be substituted when the catalog is human). **Splice-mode**
  minimap2 indexes, pre-built ONCE from the on-disk haplotype FASTAs
  (`minimap2 -x splice -d mat.splice.mmi mGorGor1.mat.cur.*.fasta.gz`); the design also accepts the raw
  FASTA in these slots (minimap2 builds the index in-line). Splice mode is required because the query is an
  exon-sum (spliced) consensus mapping across genomic introns.
- `--out <prefix>` — writes `<prefix>.parcn.tsv` + `<prefix>.parcn_families.tsv`.

## Components (each independently testable)

### 1. SUN / private-column computation — `sun_columns(copy_alleles) -> HashMap<copy, Vec<col>>`
Pure. Input: the per-copy allele matrix (copies × PSV columns) from `psv_copies.tsv`. Column *c* is a **SUN
for copy X** iff X's allele at *c* is present and differs from every OTHER copy's allele at *c* (a private
singleton). Threshold-free. A copy with ≥1 SUN column is Tier-1; a copy whose full hap-row is unique but
has no SUN is Tier-2; a copy whose hap-row equals another's is Tier-3.

### 2. Per-haplotype projection + locus dedup — reuse `genome_projection::project_families_batch`
For each haplotype (mat, then pat — one minimap2 index load each), project **all copies** of all families
(`consensuses = [(copy_id, consensus)]`, `known = {}` so NO locus is excluded — we WANT every genomic
copy), `min_identity 0.95`, `min_cov 0.90`. Union the per-copy `CopyLocus` hits within a family and **dedup
overlapping loci** (reciprocal-overlap ≥ 50%) so each distinct genomic locus is counted once, regardless of
how many copies' consensuses hit it. Output per family per haplotype: a set of distinct genomic loci.

### 3. Column transfer — `locus_hapvector(assembly_fasta, locus, consensus, psv_cols) -> Vec<Option<u8>>`
For each distinct locus, fetch the assembly subsequence (`GenomeIndex::fetch_sequence` over the haplotype
FASTA), align **the consensus of the copy that projected best to that locus** to it (guarantees the columns
are spanned; the in-repo banded Gotoh PSV aligner, minimap2 CIGAR as fallback), and walk the alignment to
read the **assembly base at each PSV column position**
(column positions are in consensus coordinates; transferred to assembly coordinates through the CIGAR).
Result: the locus's own hap-vector in the same column order as the copies.

### 4. Hybrid assignment — `assign_locus(locus_hapvector, copy_alleles, sun_columns) -> Assignment`
- **Tier-1 (SUN, deterministic):** if the locus's alleles carry some copy X's private SUN allele(s) and no
  other copy's, assign the locus to X (`method = SUN`).
- **Tier-2 (align fallback, flagged):** else assign to the copy whose full hap-vector the locus matches best
  (max agreeing columns; must beat the runner-up by ≥1 column), `method = align_fallback`.
- **Tier-3 (unresolved):** no unique best (tie / hap-row shared) → `method = UNRESOLVED`, counted in a
  per-family `n_unresolved`, never assigned to a copy.

### 5. Tabulation + output
`parCN[copy]` = number of loci assigned to it; split `loci_mat` / `loci_pat`. Emit:

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
copies.fa ─┐
psv_copies ├─ sun_columns ──────────────┐
psv_cols  ─┘                            │
                                         ▼
copies.fa ──► project_families_batch(mat) ──► dedup loci(mat) ──► column transfer ──► assign ─┐
copies.fa ──► project_families_batch(pat) ──► dedup loci(pat) ──► column transfer ──► assign ─┤
                                                                                              ▼
                                                                          tabulate ─► parcn.tsv + roll-up
```

## Error handling / graceful degradation
- minimap2 missing or non-zero exit on a haplotype → that haplotype contributes 0 loci (match
  `project_family_copies`'s existing `Ok(None)` graceful-degradation contract), a WARNING is logged, and
  the run continues (parCN then reflects only the haplotype that projected).
- A copy absent from `psv_copies.tsv` (single-copy family, no PSV columns) → reported with `sun_tier=NA`,
  `parCN` = its raw deduped locus count, `method = single_copy` (no ambiguity to resolve).
- A locus whose consensus fails to align to the fetched subsequence → dropped with a per-locus WARNING
  (counted in a `n_dropped` diagnostic, not in `n_unresolved`).

## Testing
1. **Unit — `sun_columns`:** planted allele matrix; assert exactly the private columns are returned and the
   Tier ladder (T1/T2/T3) is correct including the shared-hap-row (T3) case.
2. **Unit — `assign_locus`:** synthetic locus hap-vectors → Tier-1 SUN match, Tier-2 fallback (with the
   ≥1-column margin), and Tier-3 tie → UNRESOLVED, each asserted with the expected `method`.
3. **Unit — column transfer:** a hand-built consensus + a mutated "assembly" subsequence with a known indel;
   assert the PSV column bases are read at the correct transferred positions (indel-shifted).
4. **Integration — synthetic 2-haplotype:** two tiny FASTA "haplotypes" each carrying planted paralog copies
   with known SUNs; run the binary end-to-end; assert `parcn.tsv` gives the planted per-copy `loci_mat`/
   `loci_pat` and `famCN_diploid`.
5. **Real conservation + sanity:** on the **gorilla (GGO) catalog projected onto the gorilla mGorGor1
   haplotypes** (species-consistent — that is what makes the count valid) — assert `Σ parCN +
   n_unresolved_loci` = total distinct projected loci (conservation), and that a CN-stable gorilla family's
   `famCN_diploid ≈ 2×` its haploid catalog count (e.g. RBMY 6 haploid → ~12 diploid; DAZ, GSTM as
   secondary checks; report the actual number, don't hard-gate — real CN can deviate from exactly 2×).

## Non-goals (YAGNI)
- NOT resolving the Tier-3 collapsed residual — flagged as `n_unresolved`, routed to DNA (out of scope).
- NOT changing famCN in the existing legs (`famcn.tsv` / `famcn_readonly.tsv` / `collapsed.tsv` unchanged).
- NOT calling SUNs from the assembly de-novo — paralog identity comes from the RNA catalog's PSV matrix.
- NOT ingesting DNA reads or building a read-depth / GC / mappability model — that is the other route,
  deliberately separate.

## Success criteria
- A single `parcn` run on the existing catalog + the two on-disk haplotypes emits per-copy `parCN`
  (mat/pat split) for the SUN/hap-resolvable copies, with the Tier-3 residual flagged, in minutes and well
  within the box's RAM (two splice-index loads, one at a time).
- Deterministic SUN assignments carry a `SUN` method tag; fallback assignments are flagged `align_fallback`;
  nothing collapsed is silently assigned.
- Validation conservation holds and a CN-stable family's diploid famCN tracks ~2× its haploid annotation.
