# sedef-rs as the SD and duplicon layers: V1, V1b, V2 (2026-10-13)

**Pre-registered before any whole-genome sd-scan output is compared with SEDEF or DupMasker.** Spec:
`docs/superpowers/specs/2026-10-13-o1-layered-hierarchy-design.md`, section 2. Plan:
`docs/superpowers/plans/2026-10-13-o1-layers-phase-1-2.md`. Tool: `utilities sd-scan` (vendored sedef-rs, `src/sedef/`,
commit 96b75970), default parameters only; version key: `sedef-rs 0.1.0 vendored 2026-10-13 scan v1; seed k=16 window=2000 sketch=256 min_kmers=32 max_postings=100 min_shared=40; chain gap=2 min_score=20 merge_frac=0.5 min_cell_id=88 min_block_id=90; refine max_gap=400 diag_bin=50 min_span=500; decompose max_gap=200 diag_bin=50 min_seg=1000; max_pairs=1024; names=first-word`. V1 scorer: `bench/sedef/v1_compare.py`
(commit 9378b158).

**Not re-proposed.**
- Register 1263 (our own duplicon units, NOT ADOPTED at 0.638 held-out) is not re-run; V2 reuses its scorer and its
  split unchanged.
- SD or duplicon structure deciding family membership (register 676, 677/997, 817, 818) is not tried: sd-scan only
  annotates.

## Known before running (recorded now)

- Inverted SDs come out as window-sized pieces: on the toy genome, an 8 kb inverted copy gives four 2 kb pairs (`chain`
  links windows only along forward diagonals).
  - The base-pair union is unaffected except at a copy's ends.
  - The duplicons of an inverted SD are split at 2 kb.
- SD ends are clipped by about 1 kb (sandbox README caveat 6; the 2 kb window grid allows up to one window), and
  refined pairs stop one base short.
  Short SDs (1-5 kb) are the most exposed (recall).
- sd-scan identities are k-mer lower bounds; at V1b's 0.98 cut this costs recall.
- sd-scan writes est_id with one decimal, so at V1b's 0.98 cut a pair with est_id 97.95-97.99 rounds to 98.0 and is kept.
- sd-scan families are connected components of segment and >= 80%-overlap links. On chimp chr1 the largest holds 21.5% of
  all copies (3,750 of 17,426). Coarse duplicons in mosaic hubs are a risk for V2.
- The lifted CHM13 reference lost 925 of 151,758 SEDEF rows (rDNA rows and unliftable ends). Where their bases fall inside
  the universe, they count against sd-scan's precision.
- The vendored copy fixes two sandbox bugs that only multi-chromosome genomes show (`src/sedef/mod.rs`): chain merging
  across chromosome pairs, and family ids that followed HashMap order. Its runs are reproducible table for table
  (chimp chr1-3, two runs).

## V1 — SD calls, sd-scan vs SEDEF (both genomes)

**Inputs.**

| Genome | sd-scan run | SEDEF reference |
|---|---|---|
| CHM13 v2.0 | `chm13v2.0.fa` | `/mnt/linuxdisk/tmp/o1_close/lift/sedef_pairs.v2.bed` (v1.0 `final_v1_clean.bed`, lifted by the T2T chain; md5 7bd521edff9bd1e27b8141027da37a01) |
| Gorilla | `GGO.fasta` (mGorGor1 v2.0) | `winloci_data/GGO_sedef_final.bed` (md5 991ed34384832354c02df227e27f65b6) |

**Kept pairs.**
- SEDEF: column 21 (fracMatch) >= 0.90.
- sd-scan: status `ok` and est_id / 100 >= 0.90.
- Both: the shorter side >= 1,000 bp.

**Regions.** Both sides of every kept pair, merged per chromosome, then intersected with the universe:
- CHM13: the v2.0 bases inside the aligned blocks of `v1.0_to_v1.1_rdna_merged.chain`, i.e. bases with a v1.0 ancestor (chr1-22, X and M; not Y or the rDNA arrays);
- gorilla: every base of GGO.fasta.

**Measures.**
- recall = SEDEF region bases covered by sd-scan regions / SEDEF region bases;
- precision = sd-scan region bases covered by SEDEF regions / sd-scan region bases.

**Bar:** recall >= 0.80 AND precision >= 0.80, on CHM13 v2.0 AND on gorilla.

**Reading, fixed now.**
- **Pass:** sd-scan is the SD layer's source for every genome (the `sd` stage).
- **Fail:** sd-scan is not adopted. The SD layer keeps external SEDEF where it exists (human, gorilla); genomes without
  SEDEF get no SD layer.

## V1b — SD98

The same comparison, with 0.98 in place of 0.90 for both tools.

**Bar:** recall >= 0.80 AND precision >= 0.80, on both genomes.

**Reading.**
- **Pass:** the Soto layer (phase 3) takes its SD98 pairs from sd-scan.
- **Fail:** the Soto layer takes SEDEF SD98 where SEDEF exists, and alignment-based identity refinement in sedef-rs
  becomes its own task.

## Descriptive (V1 and V1b, no bar)

- Per chromosome: bases and coverage.
- Recall of the SEDEF pairs by identity bin and by length bin.
- Pair counts.
- Wall time and peak RSS of each sd-scan run (also the CHM13 v1.0 run, V2's input).

## V2 — duplicons, sd-scan vs DupMasker (human, F_O universe, CHM13 v1.0)

**Arms.** Both run with the binary named under **V2 binary** below; families are not recomputed from scratch, only
re-clustered from the same PAF.
- **DupMasker:**

      mcl_families --paf /mnt/linuxdisk/tmp/ourhom/loci.paf --gff /mnt/linuxdisk/tmp/ourhom/loci.gff3
        --min-exonic-bp 1 --min-shared-exon-frac 0.60 --dump-pairs --emit-hierarchy
        --hierarchy-duplicons winloci_data/duplicons/chm13.draft_v1.0_plus38Y_dupmasker_colors.bed
        --hierarchy-sd winloci_data/soto_replication/final_v1_clean.bed
        --hierarchy-external /mnt/linuxdisk/tmp/hier_flag/soto_external.tsv
        --out <arm dir>/ours

  This is the registered F_O command (10-09), re-run.
- **sd-scan:** the same command with `--hierarchy-duplicons <CHM13 v1.0 sd cache prefix>.sd.duplicons.tsv` (and its own
  `--out <arm dir>/ours`). The SD layer stays SEDEF, so only the duplicons differ.

**H0 (guard).** `ours.clusters.tsv`, `ours.loci.tsv` and `ours.pairs.tsv` are byte-identical between the arms.

**Scorer.** `bench/paf_duplicons/score_units.py`, unchanged (register 1263):
- the cover of a family (>= 2 members) = its expansions plus its `unique` class;
- family score = maximum-weight bipartite matching on Jaccard / max(#blocks);
- family-hash split: sha256 even = DEV, odd = HELD-OUT.
- Its `share_exonic_with_paf_unit` output is the sd-scan duplicons here.

**Bars.**
- **V2a:** mean family score over the HELD-OUT families >= 0.80.
- **V2b:** Soto families equal to an expansion strictly below a family, with sd-scan duplicons, >= ⌈0.9 × the DupMasker
  arm's count⌉ (10 on the 10-09 run, so >= 9).
  - Measured by `bench/o1_hier/hier_match.py MEMBERS soto_external.tsv OUT --fold LOCI --map key --null 0`,
    `exact_any.expansion_not_family`.
  - All Soto families; computed once per arm, at the held-out step.

**Procedure.**
- DEV is scored and written first.
- At most ONE amendment may follow, registered here before HELD-OUT is computed.
- HELD-OUT decides, and V2b is computed with it.

**Reading, fixed now.**
- **V2a and V2b hold:** sd-scan duplicons replace DupMasker as the expansion layer's source, for every genome.
  DupMasker stays available for comparison runs.
- **Otherwise:** sd-scan duplicons are not adopted; human expansions stay on DupMasker; apes have no expansion layer.

**Descriptive (no bar).**
- P2: sign test against the one-block null on structured families.
- P3: Soto EXPANSION placements per arm.
- P4: the two named fusions.
- Exact expansion sensitivity and precision; blocks per family.
- The NPIP family (MCL9) cover.
- Members' exonic bases carrying a duplicon, per arm.
- DEV and HELD-OUT side by side.
- Whether the re-run DupMasker arm equals the registered 10-09 files.

## Clarifications registered 2026-10-14, before any comparison

Added after this document's task review and before any sd-scan output was compared with SEDEF or DupMasker. No bar,
filter or measure changes.

- **Provenance.** The CHM13 v2.0 sd-scan run started at 10:36 on 2026-10-14, before this document's first commit
  (4e829791, 10:53), to overlap compute with other work; the gorilla run started at about 10:55 and the CHM13 v1.0 run
  (V2's input) at 11:49, both after it. Only the runs' own counts (pairs, families, run time, memory) were looked at.
  All three used the `utilities` binary built from branch `o1-layers` at ee6eb9c5, whose `src/sedef/` and
  `src/utilities.rs` equal commit 96b75970's, with the version key registered above. Starting the CHM13 v2.0 run before
  registration departs from spec section 2, which places the whole-genome runs inside the V1 gate and pre-registers each
  gate before it runs. The runs produced sd-scan's own outputs only; nothing was compared with SEDEF or DupMasker before
  the first registration commit (the plan's rule).
- **Readings apply per layer.** The SD layer follows V1; the Soto layer's SD98 source follows V1b; the expansion layer
  follows V2. V2's "Otherwise" reading concerns the duplicon and expansion layers only; the SD layer's source is whatever
  V1 decided. This narrows the spec's "sedef-rs gives the SD level only" (section 2, V2 fail), which presupposed a
  passing V1.
- **Genomes (md5 of the files as read).**
  - CHM13 v2.0: `/mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa`, 42b830df45df08ab63068e9128156e56.
  - Gorilla: `/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta`, 7a26fb029571b9e9450179769b2dc969.
  - CHM13 v1.0: `/mnt/linuxdisk/home/juanfraitu/winloci_data/soto_replication/t2t-chm13-v1.0.fa.gz`,
    67f9a15a0299e3d32df8e527a611db20 (24 records: chr1-22, X and M).
  - Relative `winloci_data/...` paths in this document mean `/mnt/linuxdisk/home/juanfraitu/winloci_data/...`.
- **V2's sd-scan arm** takes `--hierarchy-duplicons <prefix>.sd.duplicons.tsv`, where `<prefix>` is the line in
  `/mnt/linuxdisk/tmp/o1_layers/chm13v1.sd.prefix`, written by the driver's `sd` stage (`tools/rustle_pipeline.sh sd`)
  on the CHM13 v1.0 FASTA above.
- **V2 binary.** Both arms run with the one `mcl_families` binary built in Task 9 of
  `docs/superpowers/plans/2026-10-13-o1-layers-phase-1-2.md`. Its source differs from the registration commit's (the
  last commit that changes this document above `## Results`, named with the results) only by Task 9's hunks that add
  `--sd-cache` to `src/mcl_families.rs` and sd-scan rows to `parse_sd` in `src/family.rs`. Neither arm passes
  `--sd-cache` or feeds sd-scan rows to `parse_sd`; SEDEF rows parse as before. Its commit and the binary's sha256 are
  recorded with the results.
- **V2 comparator.** If the re-run DupMasker arm differs from the registered 10-09 files, the re-run arm is the
  comparator, and V2b's bar is ⌈0.9 × the re-run DupMasker arm's `expansion_not_family`⌉. If that count is 0, V2b
  cannot be decided and V2 fails.
- **H0 failure.** If `ours.clusters.tsv`, `ours.loci.tsv` or `ours.pairs.tsv` differ between the arms, V2 is not scored
  and the failure is reported. No V2 reading then applies: human expansions stay on DupMasker and the apes get no
  expansion layer, and the amendment may not be used to cure an H0 failure.
- **Task 9's A119b run is descriptive, not V2.** The plan's Task 9 prints, for human A119b, the hierarchy summary
  counts of a `--sd-cache` run beside the lifted-DupMasker run (members, families, expansions, duplicons owned, SD
  blocks and complexes). No family score, exact match or Soto placement is computed there.
- **The one amendment** may change only how the sd-scan arm's duplicon track is built from sd-scan's own outputs
  (segments, families). The bars, the scorer, the split, the DupMasker arm and V2b stay as registered.
- **V2a's number** is `P1_mean_family_score` from `score_units.py`'s HELD-OUT JSON, as printed (4 decimals); it
  passes when that value is at least 0.8000.
  - Family score = the sum of the matched Jaccards over the maximum-weight bipartite matching of the two covers' blocks
    (the DupMasker arm's and the sd-scan arm's), divided by the larger block count; a `unique` block can match only
    the other cover's `unique` block, and a non-`unique` block only a non-`unique` block
    (`score_units.py`, `family_score`).
  - Structured families, the one-block null (B1) and the `unique` class are as defined in
    `docs/PREREG_paf_duplicon_units_2026-10-10.md`; F_O is the registered DNA-level partition of the 2,334 loci of
    `/mnt/linuxdisk/tmp/ourhom/loci.paf` + `loci.gff3`, as defined in
    `docs/PREREG_hierarchy_duplicon_expansion_2026-10-07.md` ("Family layers"). Both arms re-cluster the same PAF with
    the V2 binary; the registered 10-09 partition is a descriptive reference, and H0 and the comparator rule govern any
    difference from it.
- **V2b's command** runs once per arm, with `<arm dir>` = `/mnt/linuxdisk/tmp/o1_layers/v2/dupmasker` and
  `/mnt/linuxdisk/tmp/o1_layers/v2/sdscan`: `MEMBERS` = `<arm dir>/ours.hierarchy.members.tsv`, `LOCI` =
  `<arm dir>/ours.loci.tsv`, `OUT` = `<arm dir>/key.json`, and `soto_external.tsv` =
  `/mnt/linuxdisk/tmp/hier_flag/soto_external.tsv`. The sd-scan arm's count is compared with the bar computed from
  the DupMasker arm's count. Each arm's `mcl_families` command writes `--out <arm dir>/ours`.
- **More known-before items.**
  - sd-scan uses soft-masked (lowercase) bases as k-mer bases; whether SEDEF's inputs were masked is not recorded here.
  - The one-decimal est_id rounding also applies at V1's 0.90 cut: 89.95-89.99 is kept.
  - CHM13 v2.0 has a chrY (HG002) and the v1.0 SEDEF reference has none. The X or autosome side of an sd-scan X-Y or
    autosome-Y pair therefore lies inside the universe with no possible SEDEF counterpart, a precision cost.
  - The DupMasker BED is `chm13.draft_v1.0_plus38Y` (6,306 GRCh38 chrY rows). F_O's loci lie on chr1-22 only, so those
    rows touch no F_O locus.
  - SD-end clipping: the sandbox README (caveat 6) says about 1 kb at the edges; the window grid allows up to one 2 kb
    window.

## Results

### V1 and V1b (2026-10-14)

Scored on 2026-10-14 with `bench/sedef/v1_compare.py` (commit 9378b158) and the four commands of the plan's Task 8, Step
3, after the provenance commit 092570c9 (the last commit that changes this document above `## Results`). The SEDEF
references are the registered files; their md5s were re-checked and equal the registered ones. The sd-scan inputs are
`<prefix>.sd.pairs.tsv` of the cache prefixes in `/mnt/linuxdisk/tmp/o1_layers/{chm13v2,ggo}.sd.prefix`. No bar, filter
or measure was changed after the numbers were seen. Outputs:
`/mnt/linuxdisk/tmp/o1_layers/v1_{chm13v2,ggo}_{0.90,0.98}.json`.

| genome | cut | SEDEF pairs | sd-scan pairs | SEDEF bp | sd-scan bp | recall | precision |
|---|---|---|---|---|---|---|---|
| CHM13 v2.0 | 0.90 | 44,734 | 59,692 | 237,098,879 | 102,510,802 | 0.3980 | 0.9204 |
| gorilla | 0.90 | 62,568 | 120,716 | 359,966,716 | 153,976,263 | 0.4067 | 0.9508 |
| CHM13 v2.0 | 0.98 | 3,323 | 27,332 | 131,531,378 | 63,761,476 | 0.3968 | 0.8185 |
| gorilla | 0.98 | 5,596 | 56,523 | 195,950,043 | 99,967,153 | 0.4236 | 0.8304 |

Pairs are the kept pairs (identity at or above the cut, shorter side >= 1,000 bp; sd-scan status `ok`). bp are the
merged region bases inside the universe (CHM13 v2.0: 3,044,903,725 bp; gorilla: 3,545,850,636 bp). Bases in both tools'
regions, the numerator of recall and of precision: CHM13 v2.0 94,354,434 (0.90) and 52,191,334 (0.98); gorilla
146,405,942 (0.90) and 83,009,680 (0.98).

**V1: FAIL. V1b: FAIL.** Both cuts miss the recall bar (0.80) on both genomes: recall 0.3980 (CHM13 v2.0) and 0.4067
(gorilla) at 0.90; 0.3968 and 0.4236 at 0.98. Precision clears 0.80 in all four rows (0.9204, 0.9508, 0.8185, 0.8304);
the bars were decided on the unrounded JSON values. Reading (fixed in advance):

- **V1 fail:** sd-scan is not adopted. The SD layer keeps external SEDEF where it exists (human, gorilla); genomes
  without SEDEF get no SD layer.
- **V1b fail:** the Soto layer takes SEDEF SD98 where SEDEF exists, and alignment-based identity refinement in sedef-rs
  becomes its own task.
- By the clarification "Readings apply per layer", this decides the SD layer and the Soto layer's SD98 source only; V2
  (the duplicon and expansion layers) is decided on its own bars.

Recall by SEDEF identity and length bin (descriptive, no bar; V1 = the 0.90 columns, V1b = the 0.98 columns). A cell is
the share of the SEDEF pairs' side bases inside the universe that sd-scan's regions at the same cut cover. Bins are
half-open [lo, hi); length is the shorter side of the pair. A base shared by overlapping SEDEF pairs counts once per
pair, so the bin bases (JSON `recall_bins`) add up to more than the merged SEDEF bp above: 1,967,262,036 (CHM13 v2.0)
and 2,507,676,446 (gorilla) at 0.90.

| SEDEF bin | CHM13 v2.0, 0.90 | gorilla, 0.90 | CHM13 v2.0, 0.98 | gorilla, 0.98 |
|---|---|---|---|---|
| identity 0.90-0.95 | 0.4350 | 0.3159 | - | - |
| identity 0.95-0.98 | 0.5266 | 0.5728 | - | - |
| identity 0.98-0.99 | 0.5323 | 0.7804 | 0.3527 | 0.6315 |
| identity >= 0.99 | 0.5312 | 0.7012 | 0.4601 | 0.6267 |
| length 1-5 kb | 0.3703 | 0.2180 | 0.2430 | 0.5144 |
| length 5-20 kb | 0.4582 | 0.4152 | 0.2553 | 0.6258 |
| length 20-100 kb | 0.5689 | 0.6383 | 0.4569 | 0.6390 |
| length >= 100 kb | 0.4521 | 0.5832 | 0.4037 | 0.6265 |

At 0.90 the pairs of >= 20 kb carry 77.5% (CHM13 v2.0) and 75.2% (gorilla) of the bin bases and the 1-5 kb pairs 5.9%
and 6.1%. The best single bin is 0.5689 (CHM13 v2.0) and 0.7804 (gorilla) at 0.90, 0.4601 and 0.6390 at 0.98: no bin
reaches 0.80 in any of the four columns.

Per chromosome (descriptive; full tables in the JSON `per_chrom`), over the sequences with at least 100,000 SEDEF bases
inside the universe:
- CHM13 v2.0, 0.90: 23 sequences; recall 0.1443 (chr18) to 0.5957 (chr16), no sequence reaches 0.80; precision 0.5274
  (chrX) to 0.9934 (chr8), below 0.80 on 6 of 23.
- CHM13 v2.0, 0.98: 23 sequences; recall 0.0555 (chr18) to 0.5749 (chr13), no sequence reaches 0.80; precision 0.3167
  (chr12) to 0.9752 (chr8), below 0.80 on 13 of 23.
- gorilla, 0.90: 25 sequences; recall 0.1736 (NC_073234.2) to 0.6392 (NC_073242.2), no sequence reaches 0.80; precision
  0.7835 (NC_073234.2) to 0.9966 (NC_073237.2), below 0.80 on 1 of 25.
- gorilla, 0.98: 25 sequences; recall 0.1883 (NC_073234.2) to 0.6802 (NC_073242.2), no sequence reaches 0.80; precision
  0.5737 (NC_073243.2) to 0.9528 (NC_073240.2), below 0.80 on 12 of 25.
- Left out of those ranges: CHM13 v2.0 chrM (7,615 SEDEF bases at 0.90, no sd-scan bases) and chrY (outside the
  universe); gorilla NC_011120.1, the mitochondrion (10,968 SEDEF bases at 0.90, no sd-scan bases). CHM13 v2.0 chrX has
  the lowest precision at 0.90 (0.5274); the chrY cost registered under "More known-before items" predicts a precision
  cost on X, but that split was not measured.

| run | wall s | peak RSS GiB | sequences | candidates | blocks | pairs (ok) | segments | families | duplicon rows |
|---|---|---|---|---|---|---|---|---|---|
| CHM13 v2.0 | 1111.86 | 10.31 | 25 | 600,273 | 189,574 | 189,574 (186,820) | 136,728 | 13,424 | 264,078 |
| gorilla | 3255.96 | 18.67 | 26 | 779,754 | 244,273 | 244,273 (240,384) | 209,924 | 15,404 | 401,916 |
| CHM13 v1.0 | 1080.44 | 9.44 | 24 | 544,839 | 178,980 | 178,980 (176,628) | 130,039 | 12,257 | 251,185 |

Wall s and peak RSS are `/usr/bin/time` around the whole `tools/rustle_pipeline.sh sd` call
(`bench/sedef/run_genome.sh`: 4 threads, 256 count shards, address space capped at 23 GiB by `ulimit -v 24117248`). Peak
RSS GiB = `peak_rss_kb` / 1,048,576, the unit of Task 4's projection: 10,812,312 (CHM13 v2.0), 19,578,244 (gorilla) and
9,900,128 (CHM13 v1.0) KiB. Stage seconds in `params.tsv` (seed / chain / refine / decompose): CHM13 v2.0 208.4 / 1.0 /
187.6 / 616.0; gorilla 229.1 / 1.6 / 336.5 / 2409.1; CHM13 v1.0 223.7 / 1.0 / 179.2 / 580.3. The driver logged the runs
at 10:36:09-10:54:34 (CHM13 v2.0), 10:55:02-11:49:10 (gorilla) and 11:49:36-12:07:29 (CHM13 v1.0) on 2026-10-14, as the
Provenance bullet states. The three cache prefixes are
`/mnt/linuxdisk/home/juanfraitu/sd_cache/sd/<md5>-2992c257/genome`; the md5s are the registered FASTA md5s, and the
`version_key` of each `params.tsv` equals the one registered above.

Observations (no bar):
- (a) CHM13 v1.0 has 24 sequences in `params.tsv` (chr1-22, X and M; chrM is included), not the 23 the plan expected.
  CHM13 v2.0 has 25 (the same plus Y) and gorilla 26.
- (b) The Task 4 memory projection (`bench/sedef/memory_probe.sh`, a power law through chimp chr1 and chr1-3)
  under-predicted both whole genomes. It gave 6.9 GiB (CHM13 v2.0) and 7.7 GiB (gorilla); the measured peaks are 10.31
  and 18.67 GiB, 1.5 and 2.4 times the projection (CHM13 v1.0 had no projection: 9.44 GiB). The gorilla peak RSS is 81%
  of the 23 GiB address-space cap (the cap limits virtual size, which is at least the RSS); all three runs finished with
  exit 0 and no allocation failure.
