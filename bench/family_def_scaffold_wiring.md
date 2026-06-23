# Wiring the DNA-first scaffold into the --vg pipeline — what it revealed

Goal: promote the DNA-first scaffold (family_def_dna_scaffold.py) to the production pipeline
via `--family-manifest` (template-based family assembly) and re-measure completeness end-to-end.

## Mechanically wired
- Generated the manifest from cDNA all-vs-all paralog families: `make_dna_family_manifest.py`
  -> `/home/juanfra/winloci_scratch/dna_family_manifest.tsv` (1,460 families, 5,517 loci).
- Invocation (per-chrom, OOM-safe): `rustle --vg --family-manifest <tsv> --genome-fasta <fa>
  -L <bam> -o <gtf> --vg-report <tsv>`. The flag works; bundles are matched to manifest families.

## But end-to-end recovery does NOT work cleanly — two root causes

1. **The pipeline drops manifest families by default.** `--vg` does not ingest secondary
   alignments into bundles (deliberate gate, opt-in `RUSTLE_VG_INCLUDE_SECONDARY=1`), so the
   shared multimap reads that link copies are absent -> every family is `low_shared` -> 54->0
   on the test chrom NC_073240.2. The `--vg-family-min-shared` filter is a de-novo PRECISION
   gate, redundant+harmful against a DNA-validated manifest. Even with secondaries on + filters
   relaxed, only 5/54 formed -- and with 0 shared reads (no assignment possible).

2. **The DNA scaffold is ITSELF over-merged (the decisive finding).** The cDNA all-vs-all
   paralogy (id>=0.90, recip-cov>=0.30, connected components) produces mega-"families": max
   = 389 members; 19 families >20 members hold 23% of all 5,517 member genes. The manifest
   faithfully grabs every bundle for these -> the pipeline emits an over-merged 21-copy
   "family" spanning 88 Mb. So the "DNA backwards" ground truth has the SAME bridge /
   transitive-homology problem as RNA detection -- it is not a clean external truth.

## The deep conclusion
The bridge / over-merging problem is **fundamental to homology-based family definition** --
it appears in the DNA truth (389-member components), not just RNA. "DNA backwards" does not
escape it; it relocates it. So the completeness measurements (which used these DNA families)
are inflated by the mega-families, and the scaffold cannot be cleanly wired as-is.

## Path to realize the scaffold (future work)
1. **De-bridge the DNA families** before using as a scaffold: apply the same cleaning the RNA
   side uses (contiguous-core / reciprocal-coverage / clique structure) to the cDNA all-vs-all
   components, so the manifest carries clean families, not 389-member mega-components.
2. **Add a scaffold mode to --vg**: when a manifest is provided, bypass the de-novo
   `low_shared` precision filter (the manifest is the precision gate) and use the paused
   Phase-2 secondary side-index so shared reads link copies without contaminating bundles.
Until both are done, the DNA-first scaffold is a sound prototype-level idea whose production
realization is blocked by (a) the un-cleaned DNA truth and (b) de-novo-tuned pipeline mechanics.
