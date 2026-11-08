# PREREG row 1580: damp non-SD edges before MCL (`--sd-damp`, exons + span link rule), held out on HG002 and gorilla

Written 2026-11-08, before any held-out run. The setting is chosen on A119b (dev) by the rule in "Setting" and appended
here before the held-out jobs start.

## Question

Does down-weighting homology edges that no segmental duplication (SD) supports make MCL's O1 families more precise at the
same recall, on samples the rule was not found on?

## Rule (the binary, `mcl_families --sd-damp-bed BED --sd-damp F`, uncommitted on 2026-11-08)

- Two graph nodes are **SD-linked** when ONE SEDEF pair holds, in one copy, >= half of node a's exon bases AND >= half of its
  span (first to last exon base), and the same for node b in the other copy (`SdPairs::links_exons`).
- Every edge between nodes that are not SD-linked has its weight multiplied by F; then MCL runs as usual.
- Everything else as the base families stage: `--min-exonic-bp 1 --min-shared-exon-frac 0.60 --emit-units`.

## Provenance (dev; all on A119b, so the rule is post hoc there)

- Step 0 (r1571/clu): pairs O1 gets right are 84% SD-linked, its false pairs 11%.
- The numpy prototype's big gain was a bug (`loci_in` saw only loci starting inside the SD copy). The exons-only rule gave
  +2-5 precision points at equal recall. The span half was added after diagnosing that bug: numpy, F = 0.25,
  I = 1.6 / 2.0 / 2.8 gives .867/.737, .891/.723, .912/.710 vs the plain numpy curve I = 3.2 .819/.734, 3.6 .847/.718,
  4.0 .850/.706: +5-6 points at equal recall.

## Setting (fixed on A119b binary runs r1578span before any held-out run)

Candidates: F = 0.25 at I = 1.6, 2.0, 2.8, and F = 0.5 at I = 2.0. Gate first: the binary's `sd_damp_edges_linked` must
be within 1% of the numpy span rule's 115,041 linked edges on the same graph (else the Rust span is not the prototype's and
the dev numbers do not transfer: stop). Then choose the candidate with the largest precision gain over the plain binary
curve (`pairscore.py` itself on PREFIX.fam.copies.tsv, the same scorer as the bar; linear interpolation on recall) among
candidates whose recall is >= the default (I = 2.8) recall - 0.04. The choice and its A119b numbers are appended below as "Setting chosen".

## Held-out samples

1. **HG002 (decisive).** `~/runs/hg002_kinnex/hg.families.gtf` on `hg002.pat.fa`. SDs: SEDEF 1.1 run here on
   `hg002.pat.fa` (`/anvme/.../sedef_build/hg002pat/final.bed`, default parameters, the same release as the CHM13 file). Known difference: both FASTAs are soft-masked
   but by different maskers (a 50 Mb slice: HG002 37% lowercase, CHM13 59%), so SEDEF may call a somewhat different SD set.
2. **Gorilla OR6737 (supporting).** `~/runs/gorilla_OR6737/or.families.gtf` on `GGO.fa` (mGorGor1 v2.0 RefSeq contigs).
   SDs: the user's SEDEF file for that assembly (`sedef_dl/gorilla/final.bed`, contigs checked: all present in GGO.fa).
   Truth maps human Compara by symbol onto GGO RefSeq genes, so it is small. If it has < 200 truth pairs, gorilla is
   reported descriptively and does not count toward the verdict.

## Arms per sample

- Plain curve: MCL at I = 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 5.0, no SD damping.
- SD arm: the chosen setting (the verdict).
- SD curve (reported, not the verdict): the chosen F at every plain-curve inflation, so the SD frontier can be compared
  with the plain frontier at every overlapping recall.
- Score: `pairscore.py --species hg002|gorilla`, level simiiformes, `PAIRSCORE_SAME_SYMBOL=ignore`, on PREFIX.fam.copies.tsv.

## Bar (per sample)

PASS when both hold:

1. Precision of the SD arm >= the plain curve's precision at the SD arm's recall (linear interpolation between the two
   plain points bracketing the SD arm's recall) + 0.03.
2. Recall of the SD arm >= recall of the plain default (I = 2.8) - 0.04. Why 0.04: on A119b the best candidates sit at
   -0.024 to -0.038 recall vs default (numpy), and a sample can shift the arm's recall by about a point while the gain at
   equal recall holds; clause 1 is the claim, clause 2 only stops the arm buying precision with a large recall loss.

If the SD arm's recall lies outside the plain curve's recall range, clause 1 cannot be judged: FAIL.

Verdict: PASS if HG002 passes. Gorilla is reported either way (counts as support only if it has >= 200 truth pairs).

## What follows

- PASS: propose `--sd-damp` as a documented opt-in flag (whether it becomes the default is the user's call). Next: SEDEF
  on new samples is a prerequisite, so cost is part of that decision (HG002 SEDEF wall is recorded).
- FAIL: register the row; the flag stays uncommitted or opt-in experimental at the user's choice.

## Setting chosen (2026-11-08, A119b binary, before any held-out run)

Gate: binary `sd_damp_edges_linked` 115,510 vs numpy 115,041 (+0.4%, within 1%): PASS.
`pairscore.py` (simiiformes, same-symbol ignored, 2,183 truth pairs) on A119b copies, plain binary curve:
I = 1.6 .689/.761, 1.8 .733/.753, 2.0 .750/.750, 2.4 .781/.743, 2.8 .802/.740 (default), 3.2 .837/.723, 4.0 .863/.704.

| candidate | precision / recall | plain at that recall | gain | recall vs default |
|---|---|---|---|---|
| F 0.25, I 1.6 | .883 / .736 | .810 | +0.073 | -0.004 |
| F 0.25, I 2.0 | .909 / .716 | .847 | +0.062 | -0.024 |
| **F 0.25, I 2.8** | **.925 / .718** | .844 | **+0.081** | -0.022 |
| F 0.5, I 2.0 | .860 / .737 | .808 | +0.052 | -0.003 |

Chosen: **`--sd-damp 0.25` at the default inflation 2.8** (largest gain, recall within 0.04). Held-out SD arm = `sd_I2.8`
in `r1580/heldout.sh SAMPLE 0.25`. Dev runs: `/anvme/workspace/v138da13-rustle/r1578span/`.
