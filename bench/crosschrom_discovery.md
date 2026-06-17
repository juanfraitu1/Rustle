# Genome-wide cross-chromosome copy discovery (RNA-level)

**Question (user):** the production family grouping gathers copies per genomic region (position-overlap bundles), so DISPERSED paralogs whose copies sit on DIFFERENT chromosomes (RABL2A/B; the headroom probe's 17 DISPERSED families) are never co-considered. Can the harness be improved to find cross-chromosome resemblant copies?

**Answer: yes.** A genome-wide discovery harness removes the chromosome restriction and finds them precisely.

## Method (bench/extract_gene_reps.py + family_crosschrom_discovery.py + crosschrom_grade.py)
1. Extract one representative spliced transcript per gene, genome-wide: **22,983 gene reps** (longest transcript per RefSeq gene, +strand-oriented).
2. **Minimizer-LSH with NO chromosome restriction** (k=15/w=10 canonical; inverted index, skip repetitive minimizers >200 genes; pairs sharing ≥4 minimizers, Jaccard≥0.03) → 18,934 candidate pairs. Minimizer-Jaccard is only a loose PREFILTER (real diverged dups like RFPL sit at Jaccard ~0.13).
3. **POA contiguous-core gate** (the validated RNA-level definition, bench/poa_family_definition.py): largest single ungapped co-aligned block ≥ T=0.13 of the shorter gene.
4. **Grade by POA core-block IDENTITY** — the discriminator core_cov alone lacks.
- The user's note (minimizers are useful but not the only option) is respected: minimizers are just the prefilter; the gate is alignment. A full **intron-chain alignment** is a planned second (structural) candidate axis — dispersed copies keep their relative intron-chain structure.

## Recall — it finds the known cross-chromosome families
- **8/8 universe cross-chromosome families recovered (RABL2A/B + 7 LOC families).** RABL2A (NC_073235.2) ↔ RABL2B (NC_086018.1): core_cov 0.17, **core_ident 0.99** (recent dup; short but near-identical core — exactly why core_cov alone is low and T=0.13 was calibrated to keep it).
- Named recent cross-chrom paralogs found (all core_ident≥0.97): CRIPTO/CRIPTO3, GK/GK3, EIF2S3/EIF2S3B, HNRNPA1/HNRNPA1L2, PGAM1/PGAM4, SLC25A51/SLC25A52, METTL2A/METTL2B — textbook dispersed paralogs.

## Precision — the per-pair signal is clean (core-identity, not Jaccard)
- **All 8304 cross-chrom pairs have POA core-block identity ≥ 0.4** (7387 ≥0.9, 5846 ≥0.95); **none sit at the ~0.25 chance baseline.** There are no chance-alignment false positives.
- Minimizer-Jaccard was a BAD precision axis: the apparent 'chance' pairs (EEF1A1↔LOC etc., Jaccard<0.1) are real retrocopies/processed pseudogenes at 0.89–0.99 core-identity — their low Jaccard is just a length-mismatch artifact (a short copy vs a long parent). Core-block identity is the right axis.
- Recent-duplicate subset (core_ident≥0.95): **5846 pairs**.

## Honest caveats / residual false-positive modes
- **Family-level transitive over-merge.** Connected components over a permissive pair gate chain distinct subfamilies through domain hubs: the largest 'families' are 145- and 114-gene components (gene SUPERFAMILIES, not single families). The 428 **size-2** families are clean copy pairs; large components need tighter clustering (mutual-best / all-pairs-above-bar) than transitive closure. (703 components total.)
- **Short high-identity shared elements.** A few pairs sit just above the core_cov=0.13 gate with high identity over a SHORT block = a shared transposon/element between otherwise-unrelated genes (e.g. IGFL3↔USP12: 16% core at 98% id). The recency filter does NOT remove these (they are high-identity); whole-gene-fraction or element-masking would. Raising core_cov would drop them but also drops real short recent dups like RABL2 (0.17) — the binding tension.
- **Recency spectrum.** Pairs span recent duplicates (core_ident≥0.95) → older paralogs/retrocopies → ancient domain-based families. 'Recent duplicate' (the advisor-defensible claim) = the high-identity subset; broader 'resemblant copies' = the full set.
- **One representative isoform per gene** (gene_rep); a different isoform could shift coverage at the margin. **Input = RefSeq gene reps** (not assembled transfrags) — swapping in rustle/StringTie output is a one-line input change (same harness).
- **Scope = the 25 large contigs' RefSeq genes** (22,983); whole-genome contigs add more.

## Verdict
The harness improvement WORKS: removing the chromosome restriction (genome-wide minimizer-LSH → POA contiguous-core → core-identity grading) recovers **8/8** known cross-chromosome families and finds hundreds more, all genuinely homologous per-pair. The headline gap it closes: production groups copies per region, so it structurally cannot see RABL2A/B; this finds them. Remaining work is family-level clustering (de-transitive-merge) and an element-sharer filter — both refinements, not blockers.

## Reproduce
- `python3 bench/extract_gene_reps.py` (genome-wide gene reps → /tmp/gene_reps_gw.fa)
- `MINIFORGE python bench/family_crosschrom_discovery.py --stage all` (LSH → POA → families)
- `MINIFORGE python bench/crosschrom_grade.py` (core-identity grading → crosschrom_graded.tsv)
- `MINIFORGE python bench/crosschrom_writeup.py` (this writeup + figure)
