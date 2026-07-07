# Seed-and-extend family discovery with minimap2

## Goal
Recover multi-copy gene families that have **zero candidate edges** in the
`gw_xcbase` conflict-graph catalog, using **sequence alignment** (minimap2)
instead of k-mer clustering. The idea is simple and advisor-friendly: take each
existing family member as a seed, align it against genome-wide thin loci, and add
full-length homologous thin loci as new copies.

## Why not k-mer clustering?
The advisor prefers a mechanistic, alignment-based definition of homology.
minimap2 is an alignment tool; the k-mers it uses internally are not the
scientific criterion — the alignment statistics (identity, coverage) are.

## Method (`bench/seed_extend_minimap2.py`)

1. Load `gw_xcbase` family members as **seeds**.
2. Build a FASTA database of all genome-wide thin loci (support 1–2, spliced
   sequences) extracted by `bench/thin_locus_clustering.py`.
3. Run **one minimap2 -x asm20** job: all seeds vs. all thin loci.
4. Keep a hit only if:
   - minimap2 reports **`+` strand** (the transcription-orientation sequences
     align directly, not reverse-complement);
   - **identity ≥ 0.80**;
   - **coverage of shorter ≥ 0.40**;
   - **|n_exon_seed − n_exon_locus| ≤ 1** (similar exon count);
   - **spliced-length ratio ∈ [0.6, 1.7]** (similar transcript size).
5. Collapse rescued copies that share a junction (isoforms of the same locus).
6. Exclude a rescued copy if it is ≥ 50% contained in an existing seed span.
7. Merge rescued copies back into the `gw_xcbase` catalog.

## Results

| metric | `gw_xcbase` | `gw_seedext` | delta |
|---|---:|---:|---:|
| families | 265 | 265 | 0 |
| copies | 869 | 1,275 | **+406** |
| oracle multi-copy genes covered | 51 / 57 (0.895) | 54 / 57 (0.947) | **+3** |

Newly covered oracle genes:
- `LOC101141440`
- `LOC115934629`
- `LOC129534585`

Still missed:
- `LOC109029264` — cross-chrom candidate shares k-mers but minimap2 coverage is
  too low and POA does not confirm a full-length core.
- `UBE2Q2P16` — the UBE2Q2 family is absent from `gw_xcbase`, so there is no
  seed to extend from.
- `ZNF425` — no strong alignment to any existing `gw_xcbase` family member.

### Known-family windows (±100 kb)

| family | baseline | seed-ext | delta |
|---|---:|---:|---:|
| RABL2 | 3 | 6 | +3 |
| APOBEC3 | 0 | 0 | 0 |
| MAGEA | 8 | 15 | +7 |
| ANKRD18 | 10 | 22 | +12 |
| RGPD8 | 1 | 1 | 0 |
| ZNF92 | 0 | 0 | 0 |
| GSTM | 7 | 10 | +3 |
| HERC2 | 2 | 4 | +2 |

### Annotation purity of the 406 rescued copies

| label | count |
|---|---:|
| protein_coding | 313 |
| lncRNA | 38 |
| pseudogene | 42 |
| unannotated | 11 |
| other | 2 |

~85% of rescued copies overlap a gene annotation, and ~77% are protein-coding.

## Comparison with the POA rescue prototype

The k-mer+POA rescue (`bench/gw_rescue_prototype.py`) also recovers 3/6 missed
oracle genes with +479 copies. The two approaches agree on only 42 rescued loci;
they are largely complementary. The minimap2 seed-extend is more advisor-aligned
because the family-membership criterion is a visible alignment, not a k-mer
pre-filter.

## Files

- `bench/seed_extend_minimap2.py` — implementation
- `bench/eval_gw_rescue_full.py` — evaluation script
- `bench/seed_extend_minimap2.rescued.tsv` — 406 rescued-copy rows
- `/home/juanfra/winloci_scratch/gw_seedext.{families,copies}.tsv` — merged catalog

## Next steps

1. **Recover UBE2Q2 / ZNF425 families**: these need a seed. Options:
   - Use a curated annotation transcript of `UBE2Q2` as a one-time seed to find
     its pseudogene copies (annotation-assisted rescue).
   - Run a de novo pass that builds seeds from thin-locus pairs that align
     full-length to each other (still alignment-based, not k-mer clustering).
2. **Cross-chrom LOC109029264**: try lowering the coverage floor further or
   add POA confirmation after the minimap2 candidate step.
3. **Precision control**: require that a rescued copy has read support that is
   also multi-mapped to the seed family (read-confidence gate) before adding it.
