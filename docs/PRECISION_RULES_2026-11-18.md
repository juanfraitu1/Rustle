# Precision rules on 12 samples, and the new chain score default (17-18 Nov 2026)

One page for what was tried after the ape samples arrived, what passed, what became default, and what was published.
Details and full tables are in the pre-registration documents named per row.

## Outcome
- **Default changed (user, 2026-11-18):** the chain score filter on `PREFIX.final.gtf` now uses `tools/data/chainscore_r1769.json`
  at cut 0.055. `RUSTLE_CHAINSCORE_MODEL=r1677` gives the earlier default (cut 0.061); `RUSTLE_CHAINSCORE=0` switches the filter off.
- **One rule passed its bar** (the r1769 model). Three feature families and a per-sample cut failed. Nothing else was ported.
- **Constraints set by the advisor during this work:** no samples beyond these 12 for rule finding; no annotation as guidance
  at run time.

## Samples
Human A119b (testis), HG002 (Kinnex, on CHM13); gorilla OR6737 (testis), KB3781 (fibroblast); chimpanzee, bonobo, Bornean and
Sumatran orangutan, each a testis sample and the reference individual's cell line. Truth = each genome's RefSeq, exact intron chain.

## What was tested
| row | question | test | result |
|---|---|---|---|
| 1769 | Does a chain score fitted on 11 samples of six species beat the two-sample r1677 score? | HG002 held out, same number of correct chains | PASS: +0.046 precision on HG002 at the fixed cut. Ported as `chainscore_r1769.json` |
| 1770 | Do BAM tags (second-best chain score, repeat length, minimiser density, clips, junction flanks) add precision? | leave one species out, 12 samples; bar: gain >= +0.005 on 10 of 12, none below -0.002 | FAIL: 7 of 12, median +0.006 |
| 1771 | Do read ends (cap G, clip content, shared end sites, end spread) add precision? | same | FAIL: 7 of 12, median +0.006. With row 1770: 8 of 12, median +0.010 |
| 1774 | Can the cut be set per sample without an annotation (decoy introns, score mixture)? | kept share within 0.015 of 0.95 | FAIL: fixed cut 5 of 12 within range, decoys unusable, mixture worse than the fixed cut |
| 1789 | Does splice-site sequence strength (PWM from the sample's own introns, pyrimidine tract, stronger site nearby) add precision? | same as 1770 | FAIL: 4 of 12 alone; with rows 1770 + 1771: 9 of 12, median +0.012 |
| 1792 | Diagnostic: are the wrong chains the score keeps in testis real transcripts RefSeq lacks? | replication in the other individual, relative to correct chains of equal support | Not supported: ratio >= 0.5 in 1 of 6 species; testis below cell line in 6 of 6 |

Documents: `PREREG_multisample_precision_2026-11-17.md`, `PREREG_bamtag_features_2026-11-17.md`,
`PREREG_readend_features_2026-11-17.md`, `PREREG_denovo_cut_2026-11-17.md`, `PREREG_splicesite_features_2026-11-18.md`,
`PREREG_testis_replication_2026-11-18.md`. Workspaces: `r1769_multisample_rules`, `r1770_bamtags`, `r1771_readends`,
`r1774_denovo_cut`, `r1789_splicesites`, `r1792_testis_replication` under `/anvme/workspace/v138da13-rustle`.

**Pattern in rows 1770, 1771, 1789.** The three families together gain +0.015 to +0.031 on HG002 and the six cell lines and
+0.001 to +0.004 on chimpanzee and both orangutan testis samples. Why is not known; row 1792 gives no support to "the testis
errors are real unannotated transcripts" (limit: a testis-specific transcript is not expected in a cell line).

**Why bambu is more precise.** It reports fewer transcripts. Cut to bambu's number of correct chains, Rustle is the more
precise on A119b (0.87 against 0.70 at bambu's default cut, 0.74 against 0.31 with every candidate) and OR6737 (0.78 against
0.68, 0.67 against 0.48); both samples were used to develop the score. bambu's own threshold is calibrated on an annotation,
which the advisor rules out.

## The new default
Share of correct default chains kept and precision of the kept chains.

| sample | no filter | r1677 at 0.061 (old default) | r1769 at 0.074 | r1769 at 0.055 (default) |
|---|---|---|---|---|
| A119b | 0.170 | 0.981, 0.199 | 0.932, 0.265 | 0.955, 0.242 |
| HG002 (not used to fit r1769) | 0.224 | 0.970, 0.259 | 0.966, 0.305 | 0.978, 0.280 |
| OR6737 | 0.347 | 0.954, 0.399 | 0.982, 0.385 | 0.989, 0.373 |
| KB3781 | 0.354 | 0.932, 0.388 | 0.972, 0.400 | 0.983, 0.385 |
| chimpanzee testis | 0.384 | 0.964, 0.447 | 0.978, 0.440 | 0.988, 0.426 |
| chimpanzee AG18354 | 0.337 | 0.945, 0.370 | 0.971, 0.392 | 0.982, 0.376 |
| bonobo testis | 0.315 | 0.949, 0.349 | 0.971, 0.353 | 0.983, 0.341 |
| bonobo PR00251 | 0.217 | 0.981, 0.261 | 0.936, 0.350 | 0.956, 0.319 |
| Bornean orangutan testis | 0.189 | 0.976, 0.226 | 0.938, 0.270 | 0.960, 0.251 |
| Bornean orangutan AG05252 | 0.302 | 0.965, 0.352 | 0.971, 0.370 | 0.983, 0.351 |
| Sumatran orangutan testis | 0.225 | 0.977, 0.276 | 0.953, 0.309 | 0.972, 0.289 |
| Sumatran orangutan AG06213 | 0.346 | 0.942, 0.387 | 0.978, 0.391 | 0.987, 0.377 |

- 0.055 is the loosest cut on a 0.005 grid at which every sample keeps at least 0.95 of its correct chains. The old default
  kept less than 0.95 on four samples.
- The two filters trade kept share against precision differently per sample, so this table does not rank them; the ranking
  comes from row 1769's matched-count test.
- Parity: on HG002 the driver with no environment variable drops 36,101 of 164,142 chains, identical to the analysis list;
  the old model by environment variable drops 26,316 as before (job 4480508). `cargo test --test pipeline_driver_dispatch chain` passes.
- Changed files: `tools/data/chainscore_r1769.json`, `tools/rustle_pipeline.sh` (`chainscore_model`), `tests/pipeline_driver_dispatch.rs`, `README.md`.

## HG002 in the tool comparison (scored 2026-11-18)
Primary alignments only, multi-exon chains against RefSeq on CHM13.

| tool | chains | match RefSeq | precision |
|---|---|---|---|
| Rustle | 160,977 | 36,631 | 0.228 |
| Rustle `final.gtf` (r1677 filter) | 129,404 | 35,313 | 0.273 |
| StringTie | 115,579 | 27,456 | 0.238 |
| FLAIR | 728,951 | 38,489 | 0.053 |
| IsoQuant | 67,703 | 16,087 | 0.238 |
| bambu, every candidate | 25,126 | 11,073 | 0.441 |
| bambu, default cut | 9,089 | 6,868 | 0.756 |

FLAIR's unguided collapse on the 158 GB read set needs about 440 thread-hours for its read-to-transcript alignment: with 18
threads it was at 60% after 14.4 h and was killed at the 24 h limit; the rerun with 72 threads took 7 h 43 min
(`artifact_refresh_2026-11-16/apes/hg002_fl_dn72.sbatch`; first-pass transcripts byte-identical between the two runs).

## New loci from tied secondary alignments, all 12 samples
Loci of the default assembly that the primary-only run, StringTie and FLAIR do not overlap (`bench/famcmp/seeded_loci.py`,
outputs in `artifact_refresh_2026-11-16/seeded/`).

| sample | new loci | on a RefSeq gene | spliced | exact RefSeq chain |
|---|---|---|---|---|
| A119b | 2,037 | 755 | 449 | 49 |
| HG002 | 954 | 406 | 280 | 19 |
| OR6737 | 446 | 299 | 232 | 70 |
| KB3781 | 305 | 250 | 163 | 28 |
| chimpanzee testis | 262 | 194 | 122 | 40 |
| chimpanzee AG18354 | 578 | 382 | 316 | 23 |
| bonobo testis | 429 | 309 | 286 | 94 |
| bonobo PR00251 | 776 | 343 | 215 | 16 |
| Bornean orangutan testis | 1,618 | 234 | 529 | 18 |
| Bornean orangutan AG05252 | 613 | 198 | 135 | 5 |
| Sumatran orangutan testis | 1,477 | 463 | 466 | 32 |
| Sumatran orangutan AG06213 | 331 | 207 | 176 | 11 |

## Published pages and the advisor's folder
- Assembler comparison, https://claude.ai/artifact/BKZMLTF4aqfVPD3yYjcYfd : HG002 added; the example loci tab removed; text
  cut to the interpretation; tab 3 shows the new-loci counts of every sample. Builder: `artifact_refresh_2026-11-16/build_all.py`.
- SQANTI3 categories, https://claude.ai/artifact/3neyKWVafDaKkDLqDYe5vh : HG002 added.
- Curated versus predicted RefSeq, https://claude.ai/artifact/1rHBdoUD7hbKtsCtbvJBQE : unchanged, A119b and OR6737 only.
- `advisor_share/artifacts_2026-11-17/` (synced to the lab share, 547 files, sizes verified): 12 samples with assemblies,
  annotations, score tables, SQANTI3 classifications and reads; `tables/new_loci/`; the pages; the scripts.

## Not done
- The pages and the advisor's folder still show `final.gtf` made with the r1677 filter, not the new default.
- The three failed feature families are not ported, not even as an opt-in.
- No sample outside these 12 has seen the r1769 model.
