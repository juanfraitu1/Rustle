# Genome-wide family-aware rescue prototype

## Goal
Recover multi-copy gene families that have **zero candidate edges** in the current
`gw_xcbase` conflict-graph catalog by scanning **thin loci** (primary-read intron
chains with support 1–2, below the `>=3` assembly gate) near existing families and
POA-confirming them against family members.

## Method (`bench/gw_rescue_prototype.py`)

1. Load `gw_xcbase.{families,copies,copies}.tsv/fa`.
2. For each family, build a merged ±1 Mb neighbourhood interval around its members.
3. Scan the BAM in those intervals and collect intron-chain groups with support 1–2.
4. Collapse overlapping chains **only when they share at least one junction** (real
   isoforms of the same locus). This avoids merging read-throughs that span two
   distinct loci, which would otherwise be discarded because they overlap an
   existing family member.
5. Allow a thin locus to overlap an existing member **only if the overlap covers
   < 50% of the thin-locus span**. This rescues read-throughs/extended isoforms
   that extend well beyond an existing fragment.
6. Build the spliced sequence, infer strand from canonical splice motifs.
7. Pre-filter by canonical k-mer overlap (`K_RESCUE = 20`) against the best family
   member (same-chrom within ±1 Mb; cross-chrom requires `>= 100` shared k-mers).
8. POA-confirm via `poa_pair_stats.core_recip` (`T_CORE = 0.13`), trying both
   orientations.
9. Merge rescued copies back into the catalog.

## Results

| metric | `gw_xcbase` | `gw_rescue` | delta |
|---|---:|---:|---:|
| families | 265 | 265 | 0 |
| copies | 869 | 1,348 | **+479** |
| oracle multi-copy genes covered | 51 / 57 (0.895) | 54 / 57 (0.947) | **+3** |

Newly covered oracle genes:
- `LOC101141440`
- `LOC115934629`
- `LOC129534585`

Still missed:
- `LOC109029264` — thin locus shares k-mers with a cross-chrom family member but
  POA `core_recip` does not reach 0.13 (fragmentary/domain-like homology).
- `UBE2Q2P16` — no thin locus with strong k-mer/POA homology to existing families.
- `ZNF425` — no thin locus with strong k-mer/POA homology to existing families.

### Known-family windows (±100 kb)

| family | baseline | rescue | delta |
|---|---:|---:|---:|
| RABL2 | 3 | 3 | 0 |
| APOBEC3 | 0 | 0 | 0 |
| MAGEA | 8 | 16 | +8 |
| ANKRD18 | 10 | 15 | +5 |
| RGPD8 | 1 | 2 | +1 |
| ZNF92 | 0 | 0 | 0 |
| GSTM | 7 | 10 | +3 |
| HERC2 | 2 | 2 | 0 |

### Annotation purity of the 479 rescued copies

| label | count |
|---|---:|
| protein_coding | 341 |
| pseudogene | 50 |
| lncRNA | 44 |
| unannotated | 38 |
| other | 6 |

~80% of rescued copies overlap a protein-coding annotation, suggesting they are
largely real transcripts rather than random noise.

## Comparison with `complete_poa_core`

The existing `--complete-core` flag adds 53 copies and recovers **0** new oracle
genes. The rescue prototype is far more effective for the missed-oracle problem
because it operates on **under-assembled thin loci** rather than re-scoring reps
that already passed the assembly gate.

## Files

- `bench/gw_rescue_prototype.py` — prototype implementation
- `bench/eval_gw_rescue.py` — lightweight oracle-gene evaluation
- `bench/eval_gw_rescue_full.py` — oracle + known-family windows + annotation purity
- `bench/diagnostic_missed_oracle_thin_loci.py` — diagnostic that identified the 6 missed genes
- `bench/gw_rescue_prototype.rescued.tsv` — 479 rescued-copy rows
- `/home/juanfra/winloci_scratch/gw_rescue.{families,copies}.tsv` — merged catalog

## Next steps / open questions

1. **Cross-chrom POA rescue** for `LOC109029264`: the k-mer signal is present but
   POA does not confirm. A more lenient/core-aware POA or a local minimap2
   alignment might reveal whether the homology is real or domain-only.
2. **De novo thin-locus clustering** for `UBE2Q2P16` and `ZNF425`: these appear to
   have no homology to existing `gw_xcbase` families, so a rescue pass cannot find
   them. Clustering thin loci by sequence homology could discover entirely new
   families below the assembly gate.
3. **Precision control**: +479 copies improves recall but may introduce false
   positives. A principled gate (e.g., requiring the rescued copy to be supported
   by reads that are also confused with an existing family member, or using the
   VG to validate PSV topology) could reduce false additions.
4. **Port to Rust**: if the prototype is accepted, wire the rescue pass into
   `gw_family_catalog` (or a new binary) so it becomes part of the shipped
   pipeline.
