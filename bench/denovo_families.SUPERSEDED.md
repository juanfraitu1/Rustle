# ⚠ `denovo_families.tsv` is SUPERSEDED — do NOT cite it as the O1 family catalog

`bench/denovo_families.tsv` (1,130 "families") is the **OLD catalog built by the arbitrary similarity threshold
`core_recip ≥ 0.13`** — exactly the "with no arbitrary thresholds" claim the defense audit retired. It
**over-merges**: e.g. `DNFAM0` is a single 728-member family spanning chr1→chrY (unrelated genes bridged by
shared Alu/domain sequence).

**The principled O1 catalog is the threshold-free de-tie READ-CONFLICT-GRAPH catalog**
(`gw_family_catalog → detect_conflict_catalog_genome_wide`, persisted as
`/home/juanfra/winloci_scratch/gw_conflict_catalog.{families,copies}.tsv`): **82 same-chrom families / 207
copies**, no similarity threshold. This is the catalog O2 was run on (P1: 63.9% assigned / 35.7% certified-tied
genome-wide — `bench/P1_P4_RESULTS.md`).

`denovo_families.tsv` is kept only because a few legacy scripts (`copy_resolution_census.py`,
`denovo_families.py`, `dna_psv_catalog.py`, `gw_family_catalog.rs`) still read it for backward comparison. It is
**not** the family-definition result; cite `gw_conflict_catalog` (and, for external precision, the deferred
SEDEF segdup map — `bench/SEDEF_BUILD.md`).
