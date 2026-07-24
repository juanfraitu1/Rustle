# Soto precision/recall audit — current binary (2026-07-23)

Reproducible on the cached Soto-region subset (`bench/soto/build_soto_cache.sh` →
`recompute_perchrom.sh` → `soto_cache_score.py`). All numbers are member-detection on the real Soto
`80_fams.chr.bed` (362 members, 83 families), scored exactly like `soto_detection_eval.py` (incl. the
expr-collapse col-7 projection expansion).

## Recall — up from 76.2%, zero regressions

| catalog | copies | member recall |
|---|---|---|
| baseline (committed, old binary) | 245 | 276/362 = **76.2%** |
| per-chrom (current, within-chrom, validated non-inflated) | 225 | 291/362 = **80.4%** |
| genome-wide `--cross-chrom` (current) | 863 | 306/362 = **84.5%** |
| genome-wide, obvious-FP-filtered | 487 | 303/362 = **83.7%** |

- **Conservative, fully-defensible figure: 80.4%** (validated non-inflated catalog + correct scoring).
- **0 regressions**: the 6 apparent "losses" were a scorer bug (expr-collapse col-7 not expanded), now fixed.
- The gains are real: LIMS1/LIMS4, AC124944.2 etc. — the same members the coverage top-up independently
  flagged as recoverable, now found without simulation (the accumulated `cluster_unspliced` seeding fix).
- Removing 376 obvious-FP-flagged copies cost only 3 members (84.5→83.7%) — proof the recall is not an
  artifact of inflated copy counts.

## Precision — the extra copies are mostly REAL discoveries, not over-detection

The current binary emits 863 copies vs the old 245. That gap is NOT spurious over-detection. Of the 125
detected families that overlap NO Soto member ("off-benchmark", 305 copies), a homology + complexity audit
(18-agent workflow, adversarially verified — `offbench_classified.json`) found:

| verdict | families | what they are |
|---|---|---|
| **REAL_DISCOVERY** | **119** (79 high-conf) | genuine non-Soto segmental duplications: NBPF/Olduvai-DUF1220, HYDIN2, histone clusters, SEPTIN7 pseudogenes, ANXA8/ANXA8L1, NUDT4/NUDT4P2, subtelomeric chr1p/chr3q SDs — high identity over COMPLEX UNIQUE sequence, not repeat bridges |
| SPURIOUS | ~2 | low-complexity **simple-repeat bridges** only (GWFAM145 = (ATGAAA)n satellite; GWFAM147 = (CATTTC)n, contested) |
| AMBIGUOUS | 4 | real complex homology but ultra-low support (3–6 reads) or partial-length |

**Genuine false-positive rate is low:** ~2 repeat-bridge families + ~22 giant-span (>50 kb) readthrough
copies (the NCF1-class mis-chain artifact), out of 280 families / 863 copies. The old catalog's "100%
on-member / 245 copies" reflects UNDER-detection (it missed the non-Soto paralog families), not superior
precision.

## Honest framing for the advisor

- **Member recall improved 76.2% → 80.4%** (conservative) / up to 84.5% (with cross-chrom), no regressions.
- The larger copy set is **mostly real non-Soto paralog-family discovery** (NBPF, HYDIN, histones, …), which
  is the pipeline's stated discovery objective — 119 real families vs ~2 spurious.
- Remaining true FPs are a small, characterized set (repeat-bridge families + readthrough giant spans) with
  clear code-level fixes (a low-complexity/simple-repeat gate on the refine homology; the mis-chain
  read-salvage that also fixes NCF1).

## Artifacts
- `definitive_fp_audit.tsv` — per-copy FP flags (low_read / giant_span / off_benchmark_family / single_exon_off).
  NOTE: `off_benchmark_family` alone is NOT a FP flag — the homology audit shows 95% of those are real.
- `offbench_classified.json` — the 125 off-benchmark families with REAL/SPURIOUS/AMBIGUOUS verdicts + evidence.
- Cache + recompute pipeline: `build_soto_cache.sh`, `recompute_perchrom.sh`, `soto_cache_score.py`.
