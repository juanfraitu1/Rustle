# Genome-grounded multi-copy gene family definition (read-independent, parallel)

`bench/genome_family_def.py` — defines multi-copy gene families directly on the haploid T2T genome
(GCF_029281585.2), in parallel across the genome, using inputs that **never see an RNA read**. This is
the independent oracle the adversarial reviews kept asking for: the RNA read-conflict family catalog can
be overlaid on it without circularity.

## Definition (the genome analog of the RNA read-conflict graph)

| | RNA read-conflict definition | **Genome definition (this)** |
|---|---|---|
| nodes | de-novo loci from the BAM | NCBI RefSeq gene loci (34,114) |
| edge | a read de-ties two loci | one SEDEF self-alignment segdup pair covers gene *u* on one side and gene *v* on the other |
| family | connected component, \|C\|≥2 | connected component, \|C\|≥2 |

The only membership rule is the standard "duplicated gene" criterion — a gene is *in* a segdup side iff
≥50 % of its span is covered (`COVER=0.50`). No tuned similarity threshold: SEDEF already enforces its own
identity floor on what counts as a segdup, and `|C|≥2` is the family bar. Edges and nodes are both
genome-only, so the catalog is non-circular w.r.t. the IsoSeq pipeline.

## Parallel structure

Edge emission is sharded by side-A contig over a multiprocessing pool (26 contigs), the genome analog of
`copy_assign --region-threads`; a single global union-find merges the per-contig edge lists, so
cross-chromosome paralog pairs are first-class. Whole genome (253 k segdup pairs × 34 k genes) runs in
**~5 s on 12 threads**.

## Two catalogs (built side by side)

- **A — annotated-gene families** (`genome_families_annotated.tsv`): nodes = genes; the family graph above.
- **B — de-novo duplication blocks** (`genome_families_blocks.tsv`, annotation-free): nodes = merged
  segdup footprints (units); edge = a segdup pair links two units; component ≥2 units, then *labeled* by the
  genes it covers. The graph is built without the annotation; genes only name the blocks.

Refined protein-coding view: `genome_families_protein_coding.tsv`.

## Segdup edge filter (default = raw/broad; Bailey SD(·) via `--bailey-sedef`)

This graph is an independent **oracle**, and its value is **recall** of real paralog families — so by
**default** `load_sedef_pairs` uses SEDEF's full ~50%-floor superset (all 253,029 pairs). The strict
Bailey-2002 predicate `SD(·)` (identity ≥ 90%, col 21; aln_len ≥ 1 kb, col 12; non-TE = uppercaseMatches/
aln_matches ≥ 0.50, col 28/29) is available via **`--bailey-sedef`** and keeps only **27,623 / 253,029
(10.9%)**. But `SD(·)` is a **narrower, incomparable oracle** (see `bench/SEGDUP_DEFINITION_FORMAL.md`):
it removes the repeat-*bridged* over-merge yet **drops divergent paralog families** that fall below the 90%
genomic-segdup cliff — verified losses: **CEACAM5/6/7, KRAB-ZNF (ZNF716…), IFITM1/2, ULBP1/3** all vanish
under `SD(·)` (the APOBEC3D/F = 88.4% case generalizes) — **and** it does *not* fix the transitive
single-linkage over-merge (max family still 317). That precision/recall trade is why the default is the
high-recall superset; the over-merge it admits is cosmetic noise filtered per-component when corroborating.

## Results

```
                     DEFAULT (raw superset)        --bailey-sedef (SD(·))
all-biotype          623 fam / 5162 genes / 261xc  446 / 3625 / 205
protein-coding       538 fam / 3340 genes / 223xc  390 / 2405 / 185
CATALOG B (de-novo)  285 blocks / 30,377 units      294 / 8,906
```

**Over-merge (default, honest).** The raw superset's mega-components are tRNA/snRNA/rRNA/repeat-lncRNA arrays
repeat-bridged together. `--bailey-sedef` removes the repeat-bridged edges but the transitive single-linkage
chaining **persists** (max 317): a mosaic SD sharing a valid ≥90% block with family X and another with Y
chains X–Y through two individually-valid edges. Per the formal note (§3.3), `SD(·)` is **necessary but not
sufficient**; the proper fix is the ≥2-distinct-loci / homology-component refinement the RNA side carries.

**Protein-coding size distribution** (default, 538 families): size 2 → 295, 3–4 → 124, 5–9 → 76, 10–49 → 35,
≥50 → 8 (median 2, p90 8). Mostly clean small paralog families with a short heavy tail (the repeat-bridged
mega-components, removable by `--bailey-sedef` or component refinement).

**Recovers textbook families from the genome alone** (named, small, clean): CEACAM5/6/7, **MAGEA1 /
MAGEA12 / CSAG1** (the same MAGE-A family validated on the RNA side), ZNF716/722/735/679 (KRAB-ZNF), PRSS1/2/3
(trypsinogen), IFITM1/2, ULBP1/3, LILRA2/5, FAM47A/B, RGPD8/PLGLB2. This is direct evidence the definition
identifies real multi-copy gene families without any expression data.

## Cross-checks (no extra edges)

- **Compara** (`compara_paralog_relation.json`): 3/112 multi-named families carry a Compara paralog
  relation — Compara coverage in that file is sparse (a curated subset), so this is a floor, not a rate.
- **A vs B granularity**: 580/623 (93 %) of A-families fall inside a single B-block — the gene graph and
  the block graph agree on grouping; B is coarser (merges adjacent families that share a duplication block).

## Outputs

`genome_families_annotated.tsv` · `genome_families_protein_coding.tsv` · `genome_families_blocks.tsv` ·
`genome_family_def.json` (stats).

## Cross-modal overlay vs the RNA catalog — done, with a circularity caveat

Overlaid catalog A on the RNA catalog (`denovo_families.tsv`): of within-RNA-family locus pairs, what
fraction are SEDEF-linked vs a distance/chromosome-matched null? Result (`bench/GENOME_RNA_OVERLAY.md`,
adversarially reviewed): precision 0.27 (excl. the DNFAM0 Alu over-merge), **204×** lift over a matched
null, conditioned recall 0.61. **But the adversarial panel caught a real circularity**: the RNA
"transcript" sequences the grouping scores are *reference-genome substrings* (`denovo_assemble_gate.py:58`
fetches `GGO.fasta` at read-defined exon coords), the same bytes SEDEF self-aligns — so this is a
**shared-reference-substrate** corroboration (read-derived = locus existence, expression, exon
architecture), **not** the independent oracle it first looked like. A genuinely cross-modal version
requires grouping loci by **read-consensus sequence** (PSV-carrying), or overlaying the shipped de-tie
`read_conflict.rs` graph (read-alignment-derived) instead of the POA-substring catalog.
