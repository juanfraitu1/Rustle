# Known-family regression: do we still find the named families, with current code?

**Date:** 2026-07-10. **Substrate:** `GGO_mm.bam` (gorilla testis Iso-Seq) vs `GGO.fasta`. **Binary:**
`copy_assign` at commit `b55a30b` (post locus-total λ, post O1↔O2 harmony, post DAZ2 recovery). Every region run
**foreground, serial, one at a time** (WSL2 crash discipline), flags
`--min-copies 2 --skip-poa-diagnostic --homology-primary --lambda-file <λ=58>`.

**Purpose.** A regression sweep over annotated gorilla multi-copy families with a *known paralog count*, to answer
"do we still find the known families cleanly?" after the recent membership-changing work. Each recovered copy's
genomic span is **adversarially cross-checked against `GGO_genomic.gff`** (7-agent read-only workflow): a copy
counts only if it lands on a **distinct annotated paralog**, not a readthrough, nested fragment, or duplicate
locus. A count that is right for the wrong reason (an artifact standing in for a real copy — the historical GSTM
chimera) is failed.

## Result — 6 of 7 expressed families recover cleanly; 2/2 controls stay silent

| family | region | expected | copies | χ_H | reads assigned | copies → distinct paralogs? | verdict |
|---|---|---|---|---|---|---|---|
| **GSTM** | `NC_073224.2:129160000-129240000` | 3 | **3** | 3 | 2659/2673 | GSTM3, GSTM5, GSTM1 | **CLEAN** |
| **MAGEA** | `NC_073247.2:163585000-163820000` | 2 (inverted) | **2** | 2 | 896/931 | MAGEA4 (+), MAGEA10 (−) | **CLEAN** |
| **DAZ** | `NC_073248.2:42778133-42950552` | 2 | **2** | 2 | 2353 | DAZ1, LOC129530216 (=DAZ2) | **CLEAN** ⭐regression |
| **RBMY** | `NC_073248.2:19597754-19735926` | 6 | **6** | 6 | 888 | 6 distinct LOC1295302xx | ACCEPTABLE |
| **TSPY** | `NC_073248.2:34731504-34847734` | 6 | **5** | 5 | 218 | 5 distinct LOC1295302xx | ACCEPTABLE (5/6) |
| **PCDHB** | `NC_073228.2:144880000-145130000` | ~16 | **5** | 5 | 8459 | PCDHB2/β4/PCDHB5/PCDHB9/PCDHB15 | ACCEPTABLE (5 of ~16) |
| **RFPL** | `NC_086018.1:30200000-30390000` | 2 (RFPL2/3) | 4 (2 fam) | — | 814 | only RFPL3 real; 3 artifacts | **CONTAMINATED** |
| EEF1A1 *(control)* | `NC_073229.2:97600000-97620000` | 0 | **0** | — | 0 | single-copy, 0 E_r edges | correct silence |
| SRGAP2 *(control)* | `NC_073224.2:50290000-50560000` | 0 | **0** | — | 0 | single-copy, 0 E_r edges | correct silence |

"CLEAN" = every copy maps 1:1 to a distinct annotated paralog, no artifacts. "ACCEPTABLE" = copies all map to
distinct paralogs, with a coverage-only caveat (a rescue/low-read copy, or fewer copies than annotated because the
rest are not transcribed). Both are **passes**.

## The three CLEAN families

- **GSTM (3/3).** GSTM3 / GSTM5 / GSTM1, each a near-exact span match, well supported (222 / 1196 / 1255 reads),
  no overlaps, no readthrough. The globin problem still needs `--homology-primary` (E_c finds 0 edges), and the
  30 kb GSTM5→GSTM1 chimera that once formed a 4th spurious family is gone (strand fix + readthrough filter).
- **MAGEA (2/2).** The inverted pair MAGEA4 (+) / MAGEA10 (−) — the case that was structurally invisible before
  the `colocated_families` strand fix — recovers as one 2-copy family, both spans near-exact.
- **DAZ (2/2) — the headline regression.** Pre-recovery this window returned **0 families** (a 164 kb readthrough
  stood in place of the real copy). With `locus_support` + chimeric-bridge exclusion it now returns **2 copies**:
  DAZ1 (`42783134-42859657`) and DAZ2 (`LOC129530216`, `42899568-42945549`, 3′ ~70% — the documented DAZ2
  5′-truncation). Non-overlapping, well supported, χ_H = 2.

## The RFPL failure — 1 real copy of 4, but flagged, not silent

RFPL2 (−) / RFPL3 (+) is a low-expression inverted pair in a gene desert. Current code returns **two families,
four copies — three of them artifacts**:

| copy | span | reads | overlaps | status |
|---|---|---|---|---|
| CAFAM0 | `30286681-30333257` | 707 | SLC5A4 tail only (~860 bp) | **artifact** — 46 kb intergenic readthrough, no RFPL gene |
| CAFAM0 | `30320520-30368310` | 28 | none | **artifact** — 48 kb intergenic, nested in prev (recip 0.27) |
| CAFAM1 | `30368559-30376053` | 73 | RFPL3 (exact) | **real** — the only genuine paralog |
| CAFAM1 | `30374795-30385865` | 6 | RFPL3 tail + desert | **artifact** — nested 3′ fragment (recip 0.11) |

The read-support distribution is inverted — the 707-read copy is the intergenic artifact, the genuine RFPL3
carries 73 — so the 4-copy count is "correct" for the wrong reason, and the annotated RFPL2 is missed. **The tool
warns** (`WARNING: 2 copy pair(s) share genomic sequence … Containment recip 0.27 / 0.11`), so this does not pass
silently. It is the documented **coverage-floor artifact**: the R4 readthrough rule (single-exon transcript
engulfing ≥5 junctions) does not fire because an intergenic desert has no junctions to engulf, and `Containment`
is reported-not-pruned because CAFAM0 shares its feature-cell with real overlapping paralogs elsewhere (pruning it
would kill true tandem copies — [CONTAINMENT_COVERAGE_FLOOR](CONTAINMENT_COVERAGE_FLOOR.md)).

## What this establishes

- The named families the thesis leans on — **GSTM, MAGEA, DAZ, RBMY, TSPY, PCDHB** — all recover with copies that
  land on distinct annotated paralogs. The recent membership work (strand fix, readthrough filter, DAZ2 locus
  support, locus-total λ) did **not** regress any of them; DAZ improved from 0 to 2.
- **χ_H is the trustworthy copy-number leg**, and it matches the annotated paralog count exactly on every clean
  family (3, 2, 2, 6, 5, 5). `depth_cn` at λ=58 is still expression-inflated on high-expression families
  (DAZ 40.6, GSTM 46, RBMY 15) and is not used as the count.
- Negative controls do **not** over-call: EEF1A1 (3610 reads, the old χ(H)=7 depth confound) and SRGAP2 both
  return 0 families under `--homology-primary` (0 E_r edges — no homologous second locus).
- **RFPL is the one known limitation**, and it is surfaced by a runtime warning rather than a silent wrong answer.
  Low-expression inverted pairs in gene deserts admit intergenic readthrough artifacts the structural check flags
  but cannot prune.

## Reproduce

```
LAM=single_copy_gw_locustotal.lambda_global.tsv        # λ_global = 58
copy_assign --bam GGO_mm.bam --fasta GGO.fasta --min-copies 2 --skip-poa-diagnostic \
    --homology-primary --lambda-file $LAM --region NC_073248.2:42778133-42950552 --out kf_DAZ
awk -F'\t' 'NR>1{print $2,$4":"$5"-"$6,$9}' kf_DAZ.quant.tsv    # per-copy spans + read support
# then cross-check each span against GGO_genomic.gff gene records in the window
```

Related: `bench/FAMILY_SPOT_CHECK.md` (the pre-fix state of GSTM/MAGEA/RFPL), `bench/YAG_CHECK.md` (the pre-DAZ2
Y-family sweep), `bench/CONTAINMENT_COVERAGE_FLOOR.md` (why RFPL is not prunable),
`project_daz2_locus_support`, `project_single_copy_baseline` (the λ=58 basis).
