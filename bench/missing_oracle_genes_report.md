# Rescue investigation: three oracle genes missing from `gw_xcbase + seed-extend minimap2`

**Date:** 2026-07-07
**Working directory:** `/mnt/c/Users/jfris/Desktop/Rustle`
**Author:** Kimi Code subagent

## Executive summary

The current best catalog (`gw_xcbase` + `seed_extend_minimap2`) covers **54/57** diploid-CN oracle multi-copy genes. The three still-missing genes are:

- `LOC109029264` (NC_073234.2:24,978,855-24,995,090)
- `UBE2Q2P16`     (NC_073240.2:85,814,027-85,868,602)
- `ZNF425`        (NC_073230.2:158,066,726-158,091,246)

I exhausted the requested recovery passes (lenient seed-extend, annotation-free targeted rescue, annotation-seed rescue, structural checks). The bottom line:

- **`UBE2Q2P16` and `ZNF425`** appear to be **genuinely single-copy / non-segdup in the GGO assembly and read data**. No cross-chromosome or same-chromosome paralogs are detectable with any reasonable threshold, and they cannot be rescued without inventing synthetic evidence.
- **`LOC109029264`** has a **same-chromosome partial paralog ~20 kb upstream** (NC_073234.2:24,956-24,963 kb), but it is too short/fragmented to pass the full-length homology filter (`coverage-of-shorter >= 0.40`, `len-ratio >= 0.6`) that seed-extend uses. A catalog built without the cross-chromosome restriction (e.g. `gw_off`) does include it and covers `LOC109029264`, but that catalog loses ~37 other oracle genes, so it is not a net improvement.

**Recommendation:** None of the three can be recovered within the current `gw_xcbase + seed-extend` framework by reasonable parameter changes. If the user wants to push further on `LOC109029264`, the path is to relax the *same-chromosome* restriction (not the identity/coverage thresholds), but that is a catalog-design decision with a large precision/recall trade-off.

---

## 1. Reproduced the 54/57 metric

```bash
PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python -u bench/eval_gw_rescue.py \
  --baseline-copies /home/juanfra/winloci_scratch/gw_xcbase.copies.tsv \
  --baseline-families /home/juanfra/winloci_scratch/gw_xcbase.families.tsv \
  --test-copies /home/juanfra/winloci_scratch/gw_seedext.copies.tsv \
  --test-families /home/juanfra/winloci_scratch/gw_seedext.families.tsv
```

Result:

```
Oracle multi-copy genes: 57
Baseline covered: 51 (0.895)
Test covered:     54 (0.947)
Newly covered:    ['LOC101141440', 'LOC115934629', 'LOC129534585']
Lost:             []
Still missed:     ['LOC109029264', 'UBE2Q2P16', 'ZNF425']
Families:         265 -> 265 (+0)
Copies:           869 -> 1275 (+406)
```

The metric is defined in `bench/eval_gw_rescue.py`: a gene is "covered" iff some catalog copy lies on the same chromosome and overlaps the gene span. Cross-chromosome paralogs do **not** count.

---

## 2. Why the three genes are missing

### 2.1 Thin-locus evidence is present at all three loci

`bench/diagnostic_missed_oracle_thin_loci.py` reports strong primary-read intron-chain evidence within ±1 Mb of each gene:

| gene | chrom:span | sup1 | sup2 | sup3+ | chains |
|---|---|---|---|---|---|
| LOC109029264 | NC_073234.2:24,978,855-24,995,090 | 637 | 84 | 167 | 888 |
| UBE2Q2P16 | NC_073240.2:85,814,027-85,868,602 | 389 | 41 | 60 | 490 |
| ZNF425 | NC_073230.2:158,066,726-158,091,246 | 449 | 67 | 103 | 619 |

So the loci are **expressed and spliced** in the BAM; absence from the catalog is not due to missing read coverage.

### 2.2 None of the thin loci overlap any `gw_xcbase` seed

`bench/rescue_missing_oracle.py` aligned every thin locus overlapping the three genes against all 869 `gw_xcbase` seed copies. Result:

- For the majority of gene-overlapping loci, **minimap2 reported no alignment at all** to any seed.
- The few alignments that did appear were extremely weak (`identity < 0.72`, `coverage < 0.07`) and on unrelated families.

**Conclusion:** `gw_xcbase` simply contains no family homologous to these three genes. Because seed-extend only propagates from existing seeds, it cannot create a new family for them.

### 2.3 Structural reasons

| gene | annotated biotype | annotated paralogs in GGO? | SEDEF segdup partners? | same-chrom paralog in thin loci? |
|---|---|---|---|---|
| LOC109029264 | protein_coding (ovostatin homolog 2-like) | none with same Name | yes, short cross-chrom hits (~10 kb) to several chromosomes | yes, partial same-chrom locus at ~24.957 Mb |
| UBE2Q2P16 | transcribed_pseudogene | parent `UBE2Q2` is on same chromosome (NC_073240.2:76,716,711-76,776,552) | none | none |
| ZNF425 | protein_coding (zinc finger protein 425) | none with same Name | none | none |

Key observations:

- `UBE2Q2P16` is a **processed pseudogene** of `UBE2Q2`. Both are on the same chromosome. It has no detectable additional copies, so it is excluded by the `--cross-chrom` requirement of `gw_xcbase`.
- `ZNF425` is annotated once and has no SEDEF partners and no same-chrom paralog in the thin-locus data.
- `LOC109029264` is the only one with *any* detectable paralogy: SEDEF reports short (~10 kb) cross-chrom duplications to multiple chromosomes, plus a same-chrom partial locus ~20 kb upstream. None of these partners overlap the gene's own span, so they do not satisfy the coverage metric by themselves.

---

## 3. Recovery passes attempted

### 3.1 Lenient `seed_extend_minimap2` thresholds

Created `bench/seed_extend_minimap2_param.py` (parameterized copy of `bench/seed_extend_minimap2.py`) and swept:

```bash
# default-equivalent
PYTHONHASHSEED=0 python -u bench/seed_extend_minimap2_param.py \
  --out-prefix gw_seedext_default_check

# lenient
PYTHONHASHSEED=0 python -u bench/seed_extend_minimap2_param.py \
  --out-prefix gw_seedext_lenient_70_30 \
  --min-identity 0.70 --min-coverage 0.30 \
  --max-exon-diff 2 --min-len-ratio 0.5 --max-len-ratio 2.0

# very lenient
PYTHONHASHSEED=0 python -u bench/seed_extend_minimap2_param.py \
  --out-prefix gw_seedext_lenient_60_25 \
  --min-identity 0.60 --min-coverage 0.25 \
  --max-exon-diff 2 --min-len-ratio 0.5 --max-len-ratio 2.0

# ultra lenient
PYTHONHASHSEED=0 python -u bench/seed_extend_minimap2_param.py \
  --out-prefix gw_seedext_lenient_50_20 \
  --min-identity 0.50 --min-coverage 0.20 \
  --max-exon-diff 2 --min-len-ratio 0.5 --max-len-ratio 2.0
```

Results vs oracle genes:

| thresholds | rescued copies | oracle covered | newly covered | still missed | notes |
|---|---|---|---|---|---|
| 0.80 / 0.40 (default) | 406 | 54 | LOC101141440, LOC115934629, LOC129534585 | LOC109029264, UBE2Q2P16, ZNF425 | matches original |
| 0.70 / 0.30 | 545 | 53 | LOC101141440, LOC129534585 | + LOC115934629 lost | losing genes |
| 0.60 / 0.25 | 623 | 53 | LOC101141440, LOC129534585 | + LOC115934629 lost | losing genes |
| 0.50 / 0.20 | 687 | 53 | LOC101141440, LOC129534585 | + LOC115934629 lost | losing genes |

Lowering thresholds **did not recover any of the three target genes** and began to lose previously-covered genes (`LOC115934629`). The missing genes are not a threshold issue; they are a seed-availability issue.

### 3.2 Annotation-free targeted rescue

`bench/targeted_rescue_missing_oracle.py` used the thin loci overlapping each missed gene as seeds and aligned them against the full thin-locus database with permissive filters (`identity >= 0.50`, `coverage >= 0.30`, `max-exon-diff = 2`, `len-ratio 0.5-2.0`).

- `LOC109029264`: 43 candidates, **all same-chromosome isoforms** of the gene itself.
- `UBE2Q2P16`: 96 candidates, **all same-chromosome isoforms** of the gene itself.
- `ZNF425`: 7 candidates, **all same-chromosome isoforms** of the gene itself.

**No cross-chromosome paralogs were found.** The rescue therefore cannot create a multi-chromosome family for any of the three.

### 3.3 Annotation-seed rescue

`bench/annotation_seed_extend.py` (already present in the repo) uses the annotated gene transcript as a seed. It was run previously; its output `bench/annotation_seed_extend.rescued.tsv` contains only 11 rescued copies. For `LOC109029264` it rescued only the same-chromosome partial loci at ~24.957 Mb (which do **not** overlap the gene span and therefore do not satisfy the metric). It rescued nothing for `UBE2Q2P16` or `ZNF425`.

### 3.4 Alternative base catalogs already in scratch

Evaluated every available `gw_*` catalog in `/home/juanfra/winloci_scratch/`:

| catalog | oracle covered | LOC109029264 | UBE2Q2P16 | ZNF425 | comment |
|---|---|---|---|---|---|
| gw_xcbase | 51 | no | no | no | cross-chrom base |
| gw_seedext | 54 | no | no | no | best current |
| gw_rescue | 54 | no | no | no | rescue pass |
| gw_comp2 | 51 | no | no | no | complete-core |
| gw_refined | 47 | no | no | no | `--refine` only |
| gw_refined_seedext | 49 | no | no | no | `--refine` + seed-extend |
| gw_prottail | 48 | no | no | no | `--refine --protein_tail` |
| **gw_off** | **14** | **yes** | no | no | **same-chrom allowed, but loses 37 other genes** |
| gw_sig | 13 | yes | no | no | even stricter/off |

Only `gw_off` covers `LOC109029264`, by allowing same-chromosome multi-locus families. It loses so many other oracle genes that it is not a viable replacement.

### 3.5 What would it take to recover `LOC109029264`?

The same-chromosome partial paralog at NC_073234.2:24,957,014-24,963,117 is present in the thin-locus data. When used as a synthetic seed, minimap2 aligns it to the `LOC109029264` gene locus with:

- identity ~0.85-0.92
- coverage of the **shorter** (paralog) sequence ~0.45-0.71
- coverage of the **gene** locus ~0.15-0.25

Because the paralog is only ~800 bp while the gene locus is ~3,000 bp, the full-length homology filter rejects it. A rescue pass that added this partial paralog as a new seed *and* accepted partial overlaps could recover `LOC109029264`, but that is equivalent to switching from the cross-chromosome `gw_xcbase` to a same-chromosome-aware catalog — the trade-off is the 37 genes lost by `gw_off`.

---

## 4. Files created during this investigation

All new files are under `bench/` or in `/home/juanfra/winloci_scratch/` with descriptive prefixes. Main catalog files (`gw_xcbase.*`, `gw_seedext.*`) were **not** modified.

| file | purpose |
|---|---|
| `bench/rescue_missing_oracle.py` | Diagnostic: why the gene-overlapping thin loci fail to align `gw_xcbase` seeds |
| `bench/targeted_rescue_missing_oracle.py` | Annotation-free targeted rescue using gene-overlapping thin loci as seeds |
| `bench/seed_extend_minimap2_param.py` | Parameterized seed-extend for threshold sweeps |
| `bench/missing_oracle_genes_report.md` | This report |
| `/home/juanfra/winloci_scratch/gw_seedext_default_check.*` | Default-threshold seed-extend reproduction |
| `/home/juanfra/winloci_scratch/gw_seedext_lenient_70_30.*` | Lenient seed-extend output |
| `/home/juanfra/winloci_scratch/gw_seedext_lenient_60_25.*` | Very lenient seed-extend output |
| `/home/juanfra/winloci_scratch/gw_seedext_lenient_50_20.*` | Ultra lenient seed-extend output |
| `/home/juanfra/winloci_scratch/gw_off_seedext.*` | Seed-extend from the `gw_off` catalog (proof-of-concept) |

---

## 5. Recommendation

1. **Do not change `gw_xcbase` or `gw_seedext` thresholds** in hopes of recovering the three genes. The missing genes are not threshold false-negatives; they are absent from the seed set because the catalog contains no homologous family.

2. **`UBE2Q2P16` and `ZNF425` are effectively unrecoverable** from the current data under the current definition of a multi-copy family. They appear single-copy in the GGO assembly and have no SEDEF segdup partners or same-chrom paralogs. The oracle expectation that they are multi-copy may reflect a human/reference-centric annotation that does not hold for this GGO individual.

3. **`LOC109029264` is the only one with a plausible recovery path**, but it requires a catalog-design change, not a parameter tweak:
   - Allow same-chromosome multi-locus families (as in `gw_off`), which would recover `LOC109029264`.
   - Cost: many currently-covered cross-chromosome families are lost or split, dropping oracle coverage from 54 to 14 in the extreme `gw_off` setting.
   - A tuned hybrid (e.g., keep cross-chrom families from `gw_xcbase` and additionally emit high-confidence same-chrom partial families) could in principle recover `LOC109029264` while preserving most of the 54, but that is a non-trivial new catalog mode, not a seed-extend parameter change.

4. If the user wants to pursue `LOC109029264` further, the next experiment should be: run `gw_family_catalog` with `--cross-chrom` **and** `--complete-core` plus a relaxed POA-core threshold, or build a same-chromosome-aware supplemental catalog and merge it with `gw_seedext`. Both are larger changes than the requested recovery passes.

---

## Appendix: commands for reproducibility

```bash
# Reproduce 54/57 metric
PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python -u bench/eval_gw_rescue.py \
  --baseline-copies /home/juanfra/winloci_scratch/gw_xcbase.copies.tsv \
  --baseline-families /home/juanfra/winloci_scratch/gw_xcbase.families.tsv \
  --test-copies /home/juanfra/winloci_scratch/gw_seedext.copies.tsv \
  --test-families /home/juanfra/winloci_scratch/gw_seedext.families.tsv

# Diagnostic of missing-gene thin loci vs seeds
PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python -u bench/rescue_missing_oracle.py

# Targeted rescue
PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python -u bench/targeted_rescue_missing_oracle.py

# Lenient seed-extend sweep (example)
PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python -u bench/seed_extend_minimap2_param.py \
  --out-prefix gw_seedext_lenient_70_30 \
  --min-identity 0.70 --min-coverage 0.30 \
  --max-exon-diff 2 --min-len-ratio 0.5 --max-len-ratio 2.0
```
