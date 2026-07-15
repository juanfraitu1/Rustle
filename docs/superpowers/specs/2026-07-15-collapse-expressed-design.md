# Expressed-collapsed families (`--collapse-expressed`) — design

## Goal

Recover the ~29 "dead-family" Soto members — heavily EXPRESSED multi-copy families (ANKRD36B 4501 reads,
LIMS1 1652, TCAF2 901, …) that RNA drops because their copies are **exon-identical (K=0, zero PSVs)** and so
collapse to `<2` distinct loci. Report them as a distinct **expressed-collapsed family** class with
genome-projected **copy number** (per-read unresolvable → DNA/parcn), so nothing expressed is silently
dropped. A precision-guarded, opt-in sibling of `--collapse-enumerate`.

## Motivation (grounded, this session)

- Of the 101 members missed at 72.1%: 13 silent + 25 K=0-in-detected-families + **29 dead-family (whole
  family undetected despite ≥20 reads)** + 34 low-coverage. The 29 are the category-B K=0 families
  (ANKRD36B/LIMS1/TCAF2/NCF1/BAGE2/…): exon-identical, heavily expressed, **no PSVs**.
- `--collapse-enumerate` cannot catch them: its gate requires a `hidden_copy` PSV witness (a balanced 2nd
  haplotype), and **100%-identical copies have no PSV columns** → `detect_hidden_copy` never flags → it
  correctly abstains. These families need a **PSV-free** admission path.
- `depth_cn = E_fam/λ` is NOT a usable trigger — it is expression-confounded (a highly-expressed single-copy
  gene also has depth_cn ≫ 1). The real multi-copy evidence is the **genome projection**; the EEF1A1 guard
  is **per-locus read-support** (EEF1A1's dispersed processed pseudogenes are mostly silent + <98%).

## Design

### Trigger (PSV-free, EEF1A1-safe)
At `collapse-enumerate`'s existing drop point in `detect_homology_catalog_genome_wide` (a candidate that
collapsed to `<2` RNA-distinct loci and would be dropped), when `cfg.collapse_expressed` is on:
1. Project the candidate's consensus onto the genome (`project_family_copies`, id≥0.98, cov≥0.90, with the
   candidate's own span in `known` to exclude self).
2. For EACH projected locus, count PRIMARY reads over it (`reads_in_region`, `-F 2308`); keep loci with
   `n_support ≥ 3`.
3. **Admit iff ≥2 read-supported projection loci remain.** No `hidden_copy`/PSV requirement.

**Why this gate:** ANKRD36B's copies are real ≥98% genomic paralogs, all expressed → ≥2 read-supported loci
→ admit. EEF1A1's pseudogenes are silent (no reads) and often <98% → fail per-locus read-support → reject.
Highly-expressed single-copy genes project to 1 locus → fail ≥2 → reject. This is a STRONGER EEF1A1 guard
than `collapse-enumerate`'s PSV witness, and needs no PSVs (which these families lack).

### Relationship to `--collapse-enumerate`
Same drop point, same `project_family_copies`, same copy-NUMBER-only philosophy. Independent second path:
a dropped candidate is admitted by `--collapse-enumerate` (PSV witness + projection≥2) **and/or** by
`--collapse-expressed` (projection≥2 all read-supported). The two flags are independent; both default off.
`depth_cn` is reported as OPTIONAL metadata only when a `--lambda-file` is supplied — never a gate.

### Output — `<out>.expressed_collapsed.tsv`
Distinct file so the shipped `<out>.collapsed.tsv` stays byte-identical. Columns:
```
family_id  chrom  start  end  famCN  min_locus_reads  status  projection_loci
GWFAMe0    chr2   97950885  98048181  3  32  K0_COLLAPSED_EXPRESSED  chr2:...@0.998;...
```
`famCN` = # read-supported projection loci (+ the seed = seed-inclusive, matching `collapse-enumerate`'s
`famcn_from_projection`). `min_locus_reads` = the weakest admitted locus's support (transparency).
`projection_loci` = `chrom:start-end@identity` list. Never added to `copies.tsv`/`famcn.tsv`.

## Components (each independently testable)
1. **`admit_expressed_collapse(read_supported_loci: usize) -> bool`** — pure: `>= 2`.
2. **`readmit_locus_expressed(bam_path, chrom, lo, hi, consensus, genome, fasta_path, minimap2, threads, min_locus_reads) -> Option<ExpressedCollapsedFamily>`** — project → per-locus `reads_in_region` count → keep `n_support ≥ min_locus_reads` → admit if `≥2`. Reuses `project_family_copies` + `reads_in_region`.
3. **`ExpressedCollapsedFamily { chrom, start, end, famcn, min_locus_reads, projection: Vec<CopyLocus> }`** + **`format_expressed_collapsed_row`**.
4. **Wiring:** `detect_homology_catalog_genome_wide` returns a 3rd vec `expressed: Vec<ExpressedCollapsedFamily>` (extend the current `(out, collapsed)` tuple to `(out, collapsed, expressed)`); populated at the drop point only when `cfg.collapse_expressed`. `gw_family_catalog` gains `--collapse-expressed` (+ env `RUSTLE_COLLAPSE_EXPRESSED`) and writes `expressed_collapsed.tsv` (guarded by flag && non-empty).

## Error handling / graceful degradation
- minimap2 / BAM failure → the candidate stays dropped (`None`), matching `readmit_locus`.
- A projection locus whose `reads_in_region` errors → 0 support (dropped).
- `--collapse-expressed` requires `--homology-primary` (the path with the drop point), like `--collapse-enumerate`.

## Testing
1. **Unit — `admit_expressed_collapse`:** `0/1 → false`, `2/3 → true`.
2. **Unit — `format_expressed_collapsed_row`:** exact bytes incl `K0_COLLAPSED_EXPRESSED`, seed-inclusive famCN.
3. **Integration — synthetic 2-copy genome:** a family whose two exon-identical copies both carry reads →
   `readmit_locus_expressed` admits (≥2 read-supported loci); a variant where the 2nd locus is SILENT (no
   reads) → rejects (the EEF1A1 analogue). Minimap2-gated.
4. **Byte-identical OFF:** `gw_family_catalog` on GSTM/PCDHB/MAGEA/DAZ with the flag off → md5-identical
   `families.tsv`/`copies.tsv`/`famcn.tsv`/`collapsed.tsv`, no `expressed_collapsed.tsv`.
5. **Soto A/B:** OFF vs `--collapse-expressed` ON on the Soto BAM. Report dead-family members recovered as
   expressed-collapsed copy-number (target: a measurable share of the ~29 — ANKRD36B/LIMS1/TCAF2), and that
   every admitted family's projection loci overlap a real Soto family region (precision).
6. **EEF1A1 control:** `--collapse-expressed` on the EEF1A1-scoped BAM against the full FASTA → **NO**
   `expressed_collapsed.tsv` row (its pseudogenes are silent → fail per-locus read-support). If it admits,
   the read-support floor is raised and re-measured (the isolatable-flag rationale).

## Non-goals (YAGNI)
- NOT resolving these copies per-read (they are K=0 — physically impossible from RNA; that's the point).
- NOT changing `--collapse-enumerate`'s `collapsed.tsv` (byte-identical) or the family definition.
- NOT using `depth_cn` as a gate (expression-confounded); optional reported metadata only.
- NOT touching `copies.tsv`/`famcn.tsv`.

## Success criteria
- OFF → byte-identical. ON → `expressed_collapsed.tsv` recovering a measurable share of the ~29 dead-family
  members as flagged copy-number, EEF1A1 rejected, precision preserved (every admitted family overlaps a
  real Soto region). Moves the honest ceiling statement to: *every expressed Soto family is either resolved
  per-read or flagged as a K=0 copy-number family → DNA; nothing expressed is silently dropped.*
