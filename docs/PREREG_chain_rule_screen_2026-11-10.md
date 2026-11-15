# PREREG r1619: an empirical rule on our own chains that beats the read-count curve (2026-11-10)

User question (2026-11-10): "can we try to find any other empirical rule that would give us an edge in precision and/or
recall over stringtie and flair in finding the reference chains?" Written after the tool baselines below, before any
feature was looked at.

## Baselines (multi-exon exact intron chains vs RefSeq; `r1619/curve.*.txt`, `bench/famcmp/curve.py`)
| | gorilla OR6737 chains / correct / precision | human A119b |
|---|---|---|
| ours (default) | 73,348 / 25,794 / 0.352 | 219,517 / 38,715 / 0.176 |
| StringTie | 67,491 / 23,248 / 0.344 | 211,865 / 31,781 / 0.150 |
| FLAIR | 138,125 / 24,293 / 0.176 | 626,983 / 39,230 / 0.063 |
| ours, reads >= 3 | 52,605 / 22,270 / 0.423 | 148,891 / 34,074 / 0.229 |
| StringTie, cov >= 2 | 46,307 / 18,842 / 0.407 | 154,976 / 25,520 / 0.165 |
Ours is already above StringTie's coverage curve on both samples at every cut looked at. FLAIR finds 515 more correct
chains on human at a third of the precision.

## Measure
A rule = a filter on the default chains computable from our own GTF and the genome (no other tool's output, no annotation).
Curve against curve: at the rule's correct count C, compare its precision with (a) our own reads-threshold curve and
(b) StringTie's cov curve, both linearly interpolated at C on the same sample. Dropped chains and dropped correct chains
are reported separately.

## Features screened (gorilla only; `bench/famcmp/chain_features.py`)
non-canonical splice motif; an intron within 10 bp of an intron with >= 5x its support (wobble); chain is a contiguous
run of a chain with more reads; an exon covers an intron of the locus with >= 5x the chain's reads (retention); minimum
intron support relative to the locus top; A-rich 20 bp downstream of the 3' end; reads; intron count. Read within read-count
strata, because RefSeq is incomplete and anything tied to expression looks good.
Near-dead already, not headlined: pair support (r1617), locus fraction alone (+0.01 over the reads curve).

## Bars (registered)
- Screen (gorilla): freeze ONE rule from the screen.
- Held out (A119b, `r1617/default/a.gtf`, no feature looked at there before the freeze): PASS if precision at the rule's
  correct count is >= 0.02 above our reads curve AND above StringTie's cov curve there, and the rule keeps >= 90% of the
  default correct chains.
- Second arm, only if the rule passes: rule + `RUSTLE_EXTEND_TRANSFRAGS=1` ends above default on both correct and precision.
Scope: multi-exon chains only.

## Screen on gorilla (2026-11-10, `r1619/feat.gorilla.tsv`, `r1619/rules.py`)
- Single features, share correct yes vs no: contiguous run of a stronger chain 48/1,639 (0.029) vs 0.359; retention 0.106
  (829 chains); wobble 0.192 vs 0.360 (about -0.05 within read strata); weakest intron < 5% of the locus top intron 0.160 vs
  0.428 (holds within strata: 0.120 vs 0.207 at 2 reads, 0.261 vs 0.665 at 11+); A-rich end -0.02 within strata; no
  non-canonical chains exist in the output.
- Against the reads curve at matched correct count every hand filter gains <= +0.020. A logistic model on all features
  (fitted on half the chains) gains +0.021 / +0.029 / +0.036 / +0.049 at 0.95 / 0.90 / 0.85 / 0.75 of the correct chains.
- A closed form matches the model (+0.022 / +0.027 / +0.035 / +0.045):
  **F = reads x sqrt(reads x introns / support of the locus's best-supported intron); a chain that is a contiguous run of a
  chain with more reads scores 0.** (Intron support = summed reads of our chains carrying it; locus = gene_id.)

## Frozen rule (before any A119b feature was computed)
Keep a chain when F >= 1.348 (the gorilla cut that keeps 0.92 of the correct chains). Held-out bars as registered above.

## Held out on A119b (2026-11-10, `r1619/feat.a119b.tsv`)
| | kept | correct kept | precision | reads curve at that correct count | dropped / correct among them |
|---|---|---|---|---|---|
| gorilla (screen) | 57,568 of 73,348 | 23,731 (0.920) | 0.412 | 0.387 (+0.025) | 15,780 / 2,063 |
| A119b (held out) | 153,130 of 219,517 | 34,090 (0.881) | 0.223 | 0.229 (-0.006) | 66,387 / 4,625 |
- **FAIL as registered** on both bars (gain < +0.02; correct kept 0.881 < 0.90). Still above StringTie's cov curve there
  (StringTie's whole output is 31,781 correct at 0.150), but so is the plain reads cut. Second arm not run.
- F curve vs reads curve on A119b at 0.95 / 0.90 / 0.85 / 0.75 / 0.50 of the correct chains: 0.202 vs 0.195, 0.217 vs 0.218,
  0.231 vs 0.239, 0.293 vs 0.283, 0.464 vs 0.396. The score helps only at the strict end (+0.07 at half the correct chains,
  +0.06 on gorilla); in the range that keeps most chains read count alone is as good.
- Post hoc on A119b, single features: contiguous run of a stronger chain 358 correct of 10,923 (0.033; gorilla 0.029);
  retention 307 of 6,814 (0.045); wobble 927 of 6,539 (0.142); weak intron 7,226 of 64,395 (0.112); all chains 0.176. The
  sub-chain class transfers and is what `PREFIX.final.gtf` already drops (r1498).
- Reading: where ours stands against StringTie and FLAIR is the baseline table (more correct chains than StringTie at
  higher precision on both samples; above its coverage curve at every cut). No rule from our own chain features moves the
  read-count curve by 0.02 on a held-out sample. Not adopted, nothing ported.
