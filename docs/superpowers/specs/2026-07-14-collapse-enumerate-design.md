# K=0-collapsed family enumeration (`--collapse-enumerate`) — design

## Goal

Recover near-identical multi-copy gene families that RNA collapses to fewer than 2 RNA-distinct loci —
which the pipeline currently **drops silently** at the `≥2-distinct-loci` gate — by **re-admitting** them
as *K=0-collapsed* families that report their **genome-projected copy number** plus a
*per-read-unresolvable → needs-DNA-parCN* flag. Behind a **default-off flag** so its effect is measurable in
isolation and can be deactivated without touching any other behaviour.

## Motivation (evidence, from the Soto/A119b benchmark)

- Precision is already 100% (245/245 detected copies real) and real-copy recall ≈ 90%.
- The residual misses (26 genuine, after excluding fragments/pseudogenes) are **uniformly the K=0
  near-identical collapse**: every one is 99–100% identical to a sibling (SPDYE8/9/13/15 = **100.00%**,
  ANKRD36 100%, TCAF 100%, LIMS 99.3%, NCF1 99.7%). No tunable gate helps (`--recover-copies`,
  `--protein-tail`, identity thresholds all confirmed inert — the collapse is intrinsic to RNA).
- **15 of the 26 (category B)** collapse *fully* to <2 RNA loci and are then removed by `min-copies=2`, so the
  whole family is lost from the catalog. RNA cannot resolve them *per read*, but their copy **number** IS
  recoverable by genome projection — the only thing missing is that we drop them before projecting.

## Design

### Flag
`RUSTLE_COLLAPSE_ENUMERATE=1` (env) and a `--collapse-enumerate` CLI flag on `gw_family_catalog` (and,
additively, `copy_assign`). **Default OFF.** When OFF: **byte-identical** to current output — no new families,
no changed columns, no new files.

### Hook point
In `denovo_pipeline`, the `<2-distinct-loci` / `copies.len() < min_copies` filter
(`colocated_families`, and the refine "homology component AND ≥2 distinct loci" gate) currently `continue`s
past such candidates. When the flag is ON, instead of discarding a single-locus (or sub-`min_copies`)
candidate whose reads are ambiguous, route it to the re-admission gate below.

### Re-admission gate — BOTH conditions required (this is the precision safeguard)
For a dropped candidate locus:
1. **Local collapse witness** — `hidden_copy::detect_hidden_copy` fires on that locus's reads: ≥12 balanced
   (0.20–0.60 alt-fraction) columns co-segregating on ≥5 reads = a coherent SECOND haplotype *at this locus*
   (not just MAPQ-0 ambiguity).
2. **Genome projection ≥2 loci** — `genome_projection::project_family_copies` of the candidate's consensus
   lands at **≥2 distinct genomic loci** at the famCN bucket (identity ≥0.98, cov ≥0.90).

Both fire → re-admit. Either fails → drop as today. The two-gate requirement is precisely what the earlier
**retired `collapse_gate`** lacked: that gate fired on EEF1A1 (χ(H)=7) because its MAPQ-0 reads map to
*dispersed processed pseudogenes* — ambiguity without a local second-haplotype witness *and* without a clean
≥2-locus genomic projection of a single collapsed consensus. Requiring the `hidden_copy` local witness AND the
projection is the discriminator that makes re-admission safe.

### Output
Re-admitted families are written distinctly, never as fabricated per-read copies:
- a `<out>.collapsed.tsv` (family_id, chrom, locus, `famCN` = projected copy number, `status = K0_COLLAPSED`,
  projection loci), and a `collapsed=true` marker in `<out>.families.tsv`.
- No rows added to `<out>.copies.tsv` / `.assignments.tsv` (there are no resolvable per-read copies).

## Testing

1. **Byte-identical OFF** — `gw_family_catalog` + `copy_assign` on GSTM/PCDHB/MAGEA/DAZ with the flag off →
   md5-identical `families.tsv` / `copies.tsv` / `assignments.tsv` to the current release.
2. **Unit tests** for the gate: (a) hidden-witness + ≥2 projection loci → admit; (b) hidden-witness only /
   projection only → reject; (c) an EEF1A1-style dispersed-pseudogene fixture → reject.
3. **Soto measurement** — `gw_family_catalog` on the scoped Soto BAM, flag ON vs OFF. Report: category-B
   K=0 families recovered as copy-number (target: a measurable share of the 15), precision (every re-admitted
   family's projection loci overlap a real Soto family region → 0 false re-admissions), and that all
   currently-resolved families and category-A families are UNCHANGED.

## Success criteria
- OFF → byte-identical (isolation guarantee).
- ON → recovers a measurable fraction of the 15 category-B collapsed Soto families as flagged copy-number,
  with **zero** false re-admissions (100% precision preserved), each carrying the K=0 / DNA-parCN flag.
- The flag toggles cleanly for A/B measurement and can be left off with no residual effect.

## Non-goals (YAGNI)
- NOT resolving K=0 copies per-read (physically impossible from RNA — that's the whole point).
- NOT category-A member-level recovery (within-already-detected-family projection) — a possible later
  follow-up, explicitly out of scope here.
- NOT changing any default behaviour or the meaning of existing columns.
