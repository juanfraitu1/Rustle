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

## Amendment 1 (2026-11-08, before any HG002 run): Soto-style descriptive arms

Added at the advisor's request to compare with Soto et al. 2025 (families = SD98 genes sharing exons). NOT part of the
verdict of rows 1580 or 1581. In the same held-out jobs, at the default inflation: `soto_all_d01` (`--sd-damp 0.01`,
all SEDEF pairs: damping toward 0 keeps only SD-held edges, i.e. Soto's hard rule on de novo loci, the 0.60 shared-exon
gate already on), `soto98_d01` and `soto98_d25` (the same with SD98 pairs only, SEDEF fracMatch >= 0.98). A119b
reference: `/anvme/workspace/v138da13-rustle/r1580soto_dev/`. What Soto's method cannot be run as-is here: it clusters
annotated genes and uses population famCN (WSSD of ~270 human genomes), which de novo loci, gorilla and HG002 lack.

## Amendment 2 (2026-11-08, before any HG002 run): HG002 SEDEF, candidates over 5 Mb dropped

SEDEF on HG002 (paternal + maternal vs paternal, `sedef_build/hg002pat`) crashed every time on one alignment bucket
(`bucket_0338`, `vector::_M_default_append`) on a single 21.9 Mb seed candidate: chr9:47.8-69.7 Mb (pericentromeric
satellite) vs JAKCWS010000072.1. Rule, fixed before any SD file existed: in align buckets SEDEF did not finish, a
candidate with either side > 5 Mb is dropped (no real SD is that long; such hits are satellite arrays). Finished buckets are untouched.
Dropped rows are logged in `hg002pat/align.dropped_over5mb.log`. The dry run dropped 1 of 474 rows in bucket_0338. The
final count is reported with the verdict. The rule acts only on SEDEF input, before any family run, and is blind to truth.

## Amendment 3 (2026-11-08, before any HG002 family run): final drop count, and a 1 Mb cutoff for one stalled bucket

Under the 5 Mb rule three candidates were dropped, one each in buckets 0337, 0338 and 0364. All three are in the chr9
pericentromeric satellite (chr9:47-70 Mb) against chr9 itself or an unplaced contig. After that, bucket 0364 stalled
again for over 75 min at 100% CPU on its last candidate: chr9:59.6-61.3 Mb (1.69 Mb, the same satellite) vs
JAKCWS010000073.1. Every other candidate takes seconds. For the one bucket still unfinished (0364), the cutoff was lowered
to 1 Mb (`SEDEF_MAX_SIDE=1000000` in `sedef_build/job_resume.sh`). That drops only that row; the bucket's next-largest
candidate is 0.66 Mb. The other 999 buckets finished unchanged. In total 4 of ~474,000 candidates are dropped, all
chr9 satellite. Logs: `hg002pat/align.dropped_over5mb.run1.log` and `align.dropped_over5mb.log`. Blind to truth and to
any family result.

## Outcome (2026-11-08): PASS (HG002)

`r1580/hg002/`, pairscore on PREFIX.fam.copies.tsv. Gate: 145,409 edges SD-linked, 31,178 damped.
- Plain curve (I 1.6-5.0): .223/.396, .279/.379, .338/.370, .383/.343, .466/.345 (default), .491/.324, .559/.321,
  .570/.316, .640/.314.
- SD arm `sd_I2.8`: **.719 / .321**. Clause 1: the plain curve at recall .321 is .559 (I 3.6 lies exactly there); gain
  +0.160 >= 0.03. Clause 2: .321 >= .345 - 0.04 = .305. **PASS.**
- Damped curve (I 1.6-5.0): .590/.348, .622/.345, .656/.341, .693/.338, .719/.321, .747/.314, .765/.314, .795/.300,
  .820/.297.
- Soto-style arms (amendment 1, descriptive): soto_all_d01 .741/.304; soto98_d01 .826/.184; soto98_d25 .772/.229.
  SD98 identity read from column 22 of this build (see `docs/DATA.md`).
- Caveat: HG002's truth is small (143 hits at the default), and its plain precision is far below A119b's (.466 vs .802).
  The direction matches A119b and gorilla.
