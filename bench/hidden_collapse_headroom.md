# Task (c): hidden collapsed-copy PSV signal at single-copy recall loci

Direct test of the headroom probe's undercount caveat: do read-coherence recall isoforms at
ANNOTATED single-copy loci actually sit on HIDDEN collapsed/cross-mapped paralog copies with
copy-discriminating PSV signal (annotation-free, called from the BAM pileup)?

Scanned single-copy recall loci with a linked PSV block (n_copies>=2 dumped). HIGH_COL=8.

## Verdict totals (loci with a linked PSV block)
| verdict | loci | recall isoforms | FSM |
|---|---|---|---|
| COLLAPSED_LIKE | 306 | 895 | 59 |
| AMBIGUOUS | 0 | 0 | 0 |
| HET_LIKE | 1397 | 4393 | 405 |

## n_coseg distribution (het-vs-collapse axis) — looking for a valley
```
 coseg   HET_LIKE  COLLAPSED
     0          0          0  
     1          0          0  
     2          0          0  
     3        707         46  ########################################
     4        313         29  ########################################
     5        173         16  ########################################
     6        113         10  ########################################
     7         91         11  ########################################
     8          0         52  ########################################
     9          0         32  ################################
    10          0         15  ###############
    11          0         21  #####################
    12          0         18  ##################
    13          0         10  ##########
    14          0         10  ##########
    15          0          9  #########
    16          0          4  ####
    17          0          5  #####
    18          0          2  ##
    19          0          1  #
   20+          0         15  ###############
```

## COLLAPSED_LIKE tiers (by evidence strength)
| tier | meaning | loci | recall isoforms | FSM |
|---|---|---|---|---|
| A_ge3groups | >=3 allele groups (not diploid het) | 22 | 158 | 6 |
| B_dense2copy | >=8 linked PSVs, 2 copies, low multimap | 167 | 437 | 40 |
| C_multimap | multimapping-driven (segdup/repeat possible) | 117 | 300 | 13 |

## Raw tier counts are NOT the headroom — confounds dominate (adversarial verification)
The adversarial workflow (5-mode methodology panel + hands-on BAM re-derivation, bench/wf_hidden_collapse_verify.js) confirmed every raw tier is dominated by a false-positive mode:
- **TIER B (≥8 linked cols, 2 copies) = diploid HETEROZYGOSITY.** A real 2nd genomic copy MUST
  multimap, but the TIER-B loci are UNIQUELY mapped (frac_mq0=0). The ≥8-column bar is mis-calibrated:
  scan windows average ~94 kb (whole genes, not transcripts), so a polymorphic gene trivially phases
  8–46 het SNPs (extreme cases: 46- and 30-coseg, both frac_mq0=0, MAPQ all 60, balanced minor frac).
- **TIER C (multimap) = mostly SEGDUP SPILLOVER.** frac_sec is structurally invalid: secondary
  records carry no SEQ → contribute zero alleles; the 'copies' are built from RESOLVED MAPQ-60
  primaries then OR-overridden by evidence-free spill-in reads.
- **RNA EDITING penetrates even TIER A.** Pure A>G/T>C transition spectra co-segregate as fake
  haplotypes and mint spurious ≥3-'copy' splits (36% of TIER A editing-suspect).

## Deterministic confound-controlled headroom
frac_mq0 bands over 306 COLLAPSED_LIKE: ≥0.3 (genuine local multimap)=2; [0.10,0.3)=4; <0.10 (uniquely mapped = het/editing/spillover)=300.
**Joint gate (frac_mq0≥0.3 AND n_coseg≥8) → 0 loci, 0 recall isoforms, 0 FSM.**
The 2 genuinely-multimapping loci all have n_coseg<8 (diploid-het floor, all 0 FSM) → not confidently collapse vs het-at-a-multimap-locus.

## VERDICT
- raw COLLAPSED_LIKE: 306 loci (895 recall isoforms, 59 FSM) — ALL confounds.
- naive 'TIER A+B' (BEFORE verification): 189 loci / 595 iso / 46 FSM — refuted (het + editing).
- **confound-controlled hidden-collapse headroom = 0 loci / 0 iso / 0 FSM.**
- **GO/NO-GO: NO-GO.** The direct-BAM scan, after confound control, finds 0 PSV-resolvable hidden
  collapse at single-copy recall loci — confirming the geometric probe's 0. The undercount caveat
  does NOT open real copy-resolution headroom.

Completeness gaps (honest, the detector cannot see): identical-sequence copies emit no PSVs (invisible); copies whose reads map to a separate locus/contig (RABL2/DAZ regime) never appear in this window; strand-blind (editing detected but not auto-gated); indel/STR copy differences discarded; no independent (segdup/Compara) paralog ground truth.
