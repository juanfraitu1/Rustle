# Generalized projection (`--project-all-families`) — design

## Goal

Recover missed Soto family members that are near-identical to a **resolved sibling** but not to their
family's single best consensus — the "detected-family, no-projection" bucket (~53 members, ~10–20
realistically recoverable) — by projecting **every resolved copy's consensus** onto the genome, not just one
per family. Optional, precision-guarded, and honest about the copy-vs-allele line (these are **DNA-localized
parCN** localizations, not RNA-split members).

## Motivation (grounded, this session)

- Baseline de-novo recall 215/362 = 59.4 %; `--enumerate-copies` projection already localizes +36 missed
  members (→ 69.4 %), `--protein-tail` +6 RNA-split (→ 70.4 %), at 100 % precision (`bench/soto/
  SOTO_A119B_RECOVERY.md` Panel 4).
- `--enumerate-copies` projects **one consensus per family** — the *best-supported* copy's exon-sum
  (`gw_family_catalog.rs:316-330`, `project_families_batch(..., 0.80, 0.80, ...)`). A member near-identical
  to a *different* resolved copy, but below id≥0.98 to the family's best consensus, is never localized.
- Projecting **all** resolved copies' consensuses (all reps) casts a wider net at the same precision: a
  missed member's locus is hit by whichever resolved sibling is ≥98 % identical to it.
- The irreducible floor (5 silent + ~75–90 K=0 exon-identity) is untouched by this; the honest ceiling is
  ~76 %. This rule closes part of the recoverable gap, not the floor.

## Architecture

A new opt-in path in `gw_family_catalog`, reusing the existing batched projection. **No new module** —
a function added to `genome_projection.rs` (or `collapse_enumerate.rs`) plus wiring in the binary.

### Flag
`--project-all-families` (CLI) — **default OFF**. When OFF, output is **byte-identical** to today (no new
file, no changed columns). Requires `--homology-primary` (the path that produces `fams`), like
`--enumerate-copies`. Env escape `RUSTLE_PROJECT_ALL_FAMILIES=1` for parity with the other recall flags.

### Mechanism (when ON)
1. Build `consensuses = [("{family_id}|{copy_idx}", copy.exon_sum_seq)]` over **every** copy of every
   detected family (not one-per-family).
2. `project_families_batch(&consensuses, fasta, known, 0.98, 0.90, minimap2, threads)` — ONE minimap2 index
   load; id≥0.98, cov≥0.90 (the precise bucket, not the 0.80 totalCN bucket). `known` = each family's own
   resolved copy spans, so a copy projecting back onto its own already-catalogued locus is excluded.
3. Group hits by `family_id` (split the qname on the last `|`), union across the family's copies, and
   **dedup overlapping loci** (reciprocal-overlap ≥ 0.50 — the `parcn::dedup_loci` logic) so each distinct
   genomic locus counts once regardless of how many sibling consensuses hit it.
4. **Read-support gate:** for each surviving locus, count PRIMARY reads overlapping it in the BAM
   (`reads_in_region` over `chrom:start-end`, `-F 2308` semantics — primary only) and keep the locus only if
   `n_support ≥ 3`. This confirms a real *expressed* copy, not a projection artifact — the production
   precision guardrail (on the Soto benchmark id≥0.98 alone was already 100 %, but read-support generalizes
   the guarantee off-benchmark).

### Output — `<out>.allproj.tsv`
```
family_id  chrom  start  end  identity  n_support_reads  overlaps_existing_copy
GWFAM7     chr7   75976253  75991692  0.994  41  false
```
`overlaps_existing_copy` = whether the locus overlaps a copy already in `copies.tsv` (true = redundant with
an RNA-split copy; false = a genuinely new localized locus — the recall gain). No rows added to
`copies.tsv` / `famcn.tsv` (those stay byte-identical); this is a **separate, distinctly-labeled leg**.

## Components (each independently testable)
1. **`all_copy_consensuses(fams) -> Vec<(String, Vec<u8>)>`** — pure: flatten families to `("{fid}|{ci}",
   seq)` pairs over all copies. (Reuses `DenovoTranscript.seq`.)
2. **`known_from_fams(fams) -> HashMap<family_id, Vec<(chrom,start,end)>>`** — pure: each family's copy
   spans, for the `known` self-exclusion.
3. **`dedup + support-filter`** — reuse `parcn::dedup_loci` (or a local reciprocal-overlap merge) + a
   read-count over the BAM (`reads_in_region`), admit `n_support ≥ 3`.
4. **writer** — `format_allproj_row(family_id, locus, n_support, overlaps) -> String` + the guarded file
   writer in `gw_family_catalog` (only when `--project-all-families` AND non-empty).

## Error handling / graceful degradation
- minimap2 missing / non-zero exit → `project_families_batch` returns empty (existing `Ok(...)` contract) →
  no `allproj.tsv` written, WARN to stderr, run continues.
- A family with a single copy → its one consensus is projected (the `known` self-exclusion means it only
  yields loci OTHER than its own; often none — fine).
- BAM region query failure at a locus → treat as 0 support (locus dropped), WARN once.

## Testing
1. **Unit — `all_copy_consensuses` / `known_from_fams`:** hand-built families → correct `("fid|ci", seq)`
   pairs and per-family span map.
2. **Unit — dedup + support gate:** synthetic loci (some overlapping, some below 3 support) → correct
   survivors and `overlaps_existing_copy` flags.
3. **Integration — synthetic genome:** two copies of a family whose *best* consensus misses the 2nd locus at
   id≥0.98 but the 2nd copy's consensus hits it → `--project-all-families` localizes the 2nd locus,
   `--enumerate-copies` (best-only) does not. Asserts the recall mechanism end-to-end (minimap2-gated).
4. **Byte-identical OFF:** `gw_family_catalog` on GSTM/PCDHB/MAGEA/DAZ with the flag off → md5-identical
   `families.tsv`/`copies.tsv`/`famcn.tsv` to current.
5. **Soto A/B:** OFF vs ON on the Soto BAM. Report distinct NEW missed members covered (target +10–20 toward
   ~72–74 %), and precision = every admitted locus overlaps a real Soto member (0 spurious). Update Panel 4.

## Non-goals (YAGNI)
- NOT adding rows to `copies.tsv`/`famcn.tsv` or changing the RNA-split catalog (byte-identical OFF).
- NOT changing the family DEFINITION (this is a projection-recall leg, not a membership edge).
- NOT resolving K=0 exon-identity per-read (irreducible from RNA) or recovering silent members.
- NOT re-seeding reads (that is the separate Rule 4 / assignment-fed-seeding candidate).

## Success criteria
- OFF → byte-identical. ON → a distinctly-labeled `allproj.tsv` of id≥0.98 / ≥3-read genomic loci; on Soto,
  a measurable +10–20 newly-covered missed members at preserved 100 % precision, honestly framed as
  DNA-localized parCN. The RNA-exclusive catalog and family definition are untouched.
