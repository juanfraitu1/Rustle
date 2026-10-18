# Faster family tests, realistic reads, and what O2's evidence rests on (2026-10-16)

**User goal (2026-10-16):** "integrate O1, O2 and O3 more so it's faster and cheaper to test full families; also improve the
simulations to account for more realistic ones and determine if PSVs can also Include SUNs and probably some influence of
the cigar string, nm, de tags".

Everything below is development work on family sets 8 and 9 (rounds 7 and 8), both already seen. No claim is made here.
The pre-registered test on a fresh set is round 9 (`docs/PREREG_o3_graph_round9_2026-10-16.md`): FAIL under both
profiles. Section 9 holds the independent review's corrections; where they disagree with a section above, section 9 wins.

## 1. Faster, cheaper family tests

**Where the time went (one round-8 world, 8 families, 2,537 reads):**

| step | wall | peak memory |
|---|---|---|
| simulate | 0.6 s | 38 MB |
| map to the whole genome (index load) | 79 s | 14.5 GB |
| current method (copy_assign) | 16 s | 1.06 GB |
| O1 (o1seed) | 0.6 s | 31 MB |
| O2 (integ) | 5 s | 51 MB |
| O3 (detect) | 17 s | 98 MB |

REMOVED also needs a masked whole-genome index: 146 s and 19 GB.

**What changed:**
- **`localref.py`: a homology-closed local reference.** It holds every locus any isoform of the set aligns to (minimap2
  splice:hq, `-N 500 -p 0.01`, more permissive than the reads' `-N 50 -p 0.1`), plus each copy's span, padded 10 kb and merged.
  - Set 9's is 4.6 Mb in 75 regions, against the 3.1 Gb genome; its index is 22 MB.
  - Reads are mapped with the reads' exact minimap2 command and lifted back to genome coordinates (SA tags included), so
    every later step is unchanged.
  - It is built once per family set (38 s, the only whole-genome step). Its masked REMOVED variant takes 0.5 s.
- **Validation against whole-genome mapping (all six round-8 worlds):**
  - Unique placements (MAPQ >= 1) get 100% the same primary.
  - 91-93% of reads have the identical full alignment set.
  - Every MAPQ-0 primary that moved landed on an equal-score alignment that whole-genome mapping also reports (one exception
    per world). Those are coin flips, decided arbitrarily by either index.
- **`fastset.py`: one command per family set.**
  - It simulates every world.
  - It maps each distinct read sequence once (INTACT and OUTSIDE share their reads, CROSS and BRANCH contain INTACT's, SILENT
    is a subset): 3,284 distinct sequences for 12,712 world reads.
  - It runs the current method and O1 -> O2 -> O3 per world, two worlds at a time under the light lock, then scores the end
    to end and the bars, timing every step.
- **Round 8 rerun with `fastset.py`:**
  - 237 s at 1.08 GB peak, against about 655 s at 19.4 GB for round 8's chain, plus the one-time 38-s local reference.
  - Same verdict. Five worlds have identical bars; SILENT differs by 36 correct reads, a MAPQ-0 coin flip.

## 2. Realistic reads

**`calibrate.py`** measures real Iso-Seq alignments made with the same minimap2 command on the same CHM13 v2.0:
- **Errors and donor alleles:** primary MAPQ-60 reads inside single-isoform genes, so every difference is an error or the
  donor's allele.
  - A mismatch shared by >= 3 reads at >= 20% of coverage is the donor's allele (het below 80%, hom above).
  - Indels are classed by homopolymer context.
- **Ends:** 5' and 3' truncation in transcript bases.
- **Retained introns:** one aligned block spanning a whole annotated intron.
- **Expression:** the spread of per-copy primary counts within families whose copies can be told apart.

| per aligned base | A119b (lab testis; development) | public testis (ERR13885926; held back) | old simulator |
|---|---|---|---|
| substitution errors | 0.028% | 0.092% | 0.1% |
| insertion events | 0.087% | 0.046% | 0.03% |
| deletion events | 0.125% | 0.033% | 0.03% |
| indels in homopolymers | 64-65% | 60-64% | none by design |
| donor het / hom alleles | 0.027% / 0.033% | 0.014% / 0.036% | het 0.05%, hom none |
| reads reaching the 5' end | 41% | 56% | all (0-30 bp jitter) |
| retained intron (per intron) | 0.89% | 0.14% | none |
| per-copy expression log-sd | 0.86 | 1.16 | 0 (60 reads per copy) |

**`realworld.py --sim-profile`** draws:
- each read's 5'/3' truncation from the calibrated (5', 3') pairs;
- errors at the calibrated rates, with homopolymer-biased indels and calibrated lengths;
- the donor's hom alleles on both haplotypes and het alleles on one;
- retained introns;
- a lognormal per-copy expression.

Without a profile it is byte-identical to before (round 8's INTACT and CROSS reads reproduced exactly).

**Check against real reads.** Uniquely placed reads (MAPQ 60), `calibrate.py tags`:

| | de (median) | NM/kb | mismatches/kb | ins/kb | del/kb | indels in homopolymers |
|---|---|---|---|---|---|---|
| real A119b | 0.0024 | 2.65 | 1.42 | 0.89 | 1.22 | 64% |
| real public testis | 0.0017 | 1.78 | 1.39 | 0.47 | 0.38 | 61% |
| old simulator | 0.0019 | 1.93 | 1.67 | 0.30 | 0.34 | 18% |
| A119b profile | 0.0028 | 2.87 | 0.90 | 0.82 | 1.26 | 70% |

The two real libraries differ (indels 2-3x apart), so a profile is a library, not "real reads" in general. Results are
reported under both profiles.

## 3. What realistic reads change (integrated pipeline v1 vs the current method, `fastset.py`)

| set, profile | bars | where it fails |
|---|---|---|
| 9, old simulator | FAIL (round 8) | CROSS errors 1 vs 0 |
| 9, A119b | FAIL | B1: INTACT 926 vs 932, CROSS 1,011 vs 1,030; B2: REMOVED 8 vs 0 |
| 9, public | FAIL | B2: REMOVED 7 vs 0 |
| 8, old simulator | PASS | |
| 8, A119b | FAIL | B1: INTACT -93, REMOVED -70, SILENT -19; B2: CROSS 13 vs 8, REMOVED 17 vs 3, INTACT and SILENT 1 vs 0 |
| 8, public | FAIL | B2: CROSS 10 vs 8, REMOVED 8 vs 4 |

The integrated pipeline's advantage on catalogued copies was partly an artifact of idealized reads. It keeps its advantage
where the catalog is incomplete: OUTSIDE errors 24 vs 154 and BRANCH 1,194 vs 932 correct (set 9, A119b).

**Why:**
1. **O2 calls more reads `unexplained`:** 197 of 2,537 under A119b, against 8 with idealized reads. Of 139 examined:
   - 38 align to no path;
   - 51 have NM above delta x read length, because a homopolymer indel of 4+ bases or a retained intron counts per base;
   - 21 are clipped >= 50 bases;
   - 29 are short truncated reads where a few edits exceed delta x read length.
2. **REMOVED's catalogued errors are present-copy reads inside true absent-copy calls:** 7 of 7 on set 9 and 14 of 15 on
   set 8 under A119b. A119b's noisier reads widen the clusters.

## 4. Do PSVs include SUNs?

**Yes: a SUN is the PSV column where one copy's base differs and every other copy agrees.**
- Of set 9's 284 PSV columns, 196 are SUNs (69%), 77 shared and 11 two-copy pairs. The de novo families of 10-12 gave 78%.

From the current O2's read-star proofs (`copy_assign --dump-star`, which changes no call), each decisive column is classed
among the read's candidates:
- **About 70% of decisive columns are SUN-type:** the called copy's own SUN, or another copy's SUN the read rules out.
- **Complete catalogs (INTACT, CROSS): no call is wrong** under any profile, with or without an own-SUN column (0 of 827 to
  1,226 per world).
- **Missing copies are where every error is.**
  - OUTSIDE under idealized reads: 60 of 121 wrong calls rest on the called copy's own SUN. The uncatalogued copy carries that
    "private" base too.
  - **A SUN is private only relative to the copies considered.** Genome-wide SUNs (unique k-mers, Sudmant 2010) would guard
    against copies in the genome but missing from the catalog. Nothing in the reference guards against copies absent from it.

## 5. CIGAR, NM and de

- **de flags absent-copy reads.**
  - In REMOVED, the current O2's wrong calls sit in the reads' top de quartile: 85 of 254 against 26 of 354 in the bottom
    quartile (idealized), and 39/191 against 12/202 under A119b.
  - Within a complete catalog it predicts nothing, because there are no errors.
- **Soft clips (CIGAR S).**
  - O2's distance counts every clipped base as an edit. That made round 8's error: one sequencing error among 4 clipped
    bases scored 1 vs 5.
  - Aligning clipped ends against the path (semi-global) removes that error. But it also removes the inflated margins of a
    few correct single-difference calls: about -3 correct per world on set 9, +3 to +8 on set 8.
- **NM vs de: count indel events, not indel bases.**
  - Real errors are mostly homopolymer indels, and a retained intron is one event of many bases.
  - A gap-compressed distance (mismatches + indel events, like minimap2's `de`, with semi-global ends) recovers O2 on
    realistic reads:
    - set 9, A119b: INTACT 926 -> 953 correct, unexplained 197 -> 106;
    - set 8, A119b: INTACT 2,141 -> 2,215 correct, wrong 1 -> 0; CROSS 2,628 -> 2,719 correct, wrong 9 -> 1.
  - Under idealized reads it costs some calls on set 8 (INTACT 2,476 -> 2,407), with fewer errors (CROSS wrong 11 -> 3).
- **The read's own de cleans absent-copy calls.**
  - Present-copy reads inside true absent-copy calls sit within delta of their own locus: 51 of 55 at de <= delta. True
    absent-copy members do not: 45 of 1,154 at de <= delta.
  - A per-read genome check (`detect.py --member-genome-check`, arm I) drops them.

## 6. Integrated v2 (opt-in; defaults byte-identical)

**What it runs:** `o2graph.py run W integ 2 --distance gapcomp` + `detect.py W --member-genome-check`, or in one command
`fastset.py --o2-distance gapcomp --member-genome-check`.

**Development results** (bars per world, as round 8; the per-read column also classes a REMOVED copy's read as detectable by
its own `de`, `summary.py e2e --per-read-detectability`):

| set, profile | v1 | v2, copy-level | v2, per-read |
|---|---|---|---|
| 9, idealized | FAIL: CROSS errors 1 vs 0 | PASS | |
| 9, A119b | FAIL: INTACT, CROSS, REMOVED | FAIL: REMOVED errors 1 vs 0 | PASS |
| 9, public | FAIL: REMOVED errors 7 vs 0 | FAIL: REMOVED errors 1 vs 0 | PASS |
| 8, idealized | PASS | PASS | |
| 8, A119b | FAIL: 4 worlds | FAIL: INTACT 2,215 vs 2,234 correct; REMOVED errors 5 vs 3 | FAIL: INTACT 2,215 vs 2,234 correct |
| 8, public | FAIL: CROSS, REMOVED | PASS | PASS |

- **Set 9's one REMOVED error under both profiles** is a read of the removed USP17L.8 (copy median `de` 0.018) that is
  near-identical to USP17L.9 over its own truncated span. Per-read detectability classes it unknowable.
- **Set 8 under A119b reads, INTACT's shortfall (-19)** is O2's `unexplained` reads (SPRR2 37, AQP7 12), which the current
  method assigns.
  - A MAPQ-60 genome fallback for unexplained reads was measured and rejected.
  - Outside BRANCH it places 17 of 329 reads wrong, in worlds where the current method makes almost no errors.
  - In BRANCH it places 159 of 220 wrong: novel-exon reads land uniquely on their parent copy.
- **v2 against v1 under idealized reads:** fewer correct reads on set 8 (INTACT 2,407 vs 2,476), fewer errors (CROSS 3 vs 11,
  REMOVED 0 vs 13).

**Round 9** tests v2 on a fresh set: `docs/PREREG_o3_graph_round9_2026-10-16.md`.

## 7. Full families

Uncapped (`families.py build --cap 1000`; `/mnt/linuxdisk/tmp/o3_graph/families_full.json`):
- set 2: NPIPA 6, NPIPB 14, TBC1D3 10, GOLGA8 18 and four others; 69 copies;
- set 9: USP17L 59, TSPY 35 and six others; 116 copies.

Runs with `fastset.py` (integrated v2, A119b profile; post hoc, both sets already seen):

| set | copies | wall time | peak memory | heaviest steps (summed over worlds) | bars |
|---|---|---|---|---|---|
| 2 (NPIP, TBC1D3, GOLGA8, ...) | 69 | 1,508 s | 14.4 GB (local reference build only) | O3 1,215 s, current 457 s, O2 434 s | PASS |
| 9 (USP17L 59, TSPY 35, ...) | 116 | 1,203 s | 14.3 GB (local reference build only) | current 579 s, O2 543 s, O3 494 s | PASS |

- **Full families are feasible but no longer cheap.**
  - O3's `vg msga` family graphs grow steeply with copies and locus length.
  - copy_assign peaks at about 5 GB per world on set 9.
- **Since then:**
  - `fastset.py` shares family graphs between worlds with one catalog (`O3_GRAPH_CACHE`). INTACT, CROSS, BRANCH and SILENT
    share the catalog, so their graph is built once.
  - `--heavy-current` runs copy_assign under the heavy lock.
- These sets used the old family rule (section 8).

## 8. Domain sharers in the simulation's families (user, 2026-10-16)

The question: do the families hold genes that share only a domain with the others (over-representation)? Would requiring
links over >= 80-90% of the shorter transcript remove them?

**Families today** (`families.py`):
- genes are taken by name pattern;
- two are linked when their canonical transcripts align at >= 0.90 identity over >= half the shorter one;
- a family is the largest single-linkage component.

**Yardstick (`domains.py`).** Each gene is graded against its family by spliced alignment of its isoforms to the other genes'
genomic LOCI, in both directions:
- `whole` (each locus holds >= 80% of the other's transcript);
- `fragment` (one direction only);
- `partial` (50-80%);
- `segment` (< 50% either way: a domain sharer).

It uses the partner's genome, not its annotation. An earlier locus-to-locus version failed on old duplicates whose introns
diverged (SSX, AMY scored 0) and was replaced.

**Results over all ten sets' candidate families, uncapped (105 families with >= 2 genes):**

| link rule | genes | whole | fragment | partial | segment | genes dropped (of which whole + fragment) |
|---|---|---|---|---|---|---|
| canonical transcripts, >= 0.5 of the shorter (today) | 644 | 591 | 41 | 10 | 2 | |
| canonical transcripts, >= 0.8 | 588 | 545 | 33 | 4 | 0 | 56 (41) |
| canonical transcripts, >= 0.9 | 548 | 506 | 26 | 4 | 1 | 96 (73) |
| **locus: canonical transcript >= 0.8 inside the other gene's locus, one direction** | 670 | 617 | 53 | **0** | **0** | |

- **Domain sharers were rare.**
  - Only 2 of 644 genes (KIR2DL4, KRTAP9-7) were segments.
  - 10 were partial: NBPF7, HERC2P7, GUSBP3, SPANXN1/5, KRTAP5-4/-10, KRTAP4-2/-12/-16.
  - Over-representation was 1.9% of genes, not more.
- **Transcript links over >= 80% do remove them, at about five genuine copies per non-copy removed.**
  - The losses are copies whose canonical transcripts carry different exon sets; NPIPB's 14 whole copies fall to 7.
- **The locus rule removes all 12 and adds 38 copies the transcript rule missed:**
  - pseudogene copies such as LRRC37A and GSTA7P;
  - diverged-annotation copies such as NPIPB 14 -> 17, SSX 5 -> 8, PSG 7 -> 9.
- **Linking on any isoform instead of the canonical one inflates families** through short single-domain isoforms (LILR 2 -> 9,
  MT1 2 -> 9). The rule uses the canonical transcript for that reason.
- **The rebuilt simulation sets** (`/mnt/linuxdisk/tmp/o3_graph/families_locus.json`, each set's original cap and selection)
  change in every set.
  - Membership changes, for example GOLGA6L 4 -> 7, SPATA31 7 -> 5, MAGEA 5 -> 3.
  - Four first-8 selections shift: set 5 VCX -> CTAG, set 6 ANKRD20A -> CYP2D, set 8 FAM138 -> PRSS, set 10 KRTAP4 -> KRTAP9.

**Runs on the rebuilt sets** (integrated v2, v1 profiles; post hoc; per-read scoring):

| set | old families | rebuilt families |
|---|---|---|
| 10, A119b | FAIL: B1 INTACT -40, CROSS -37, REMOVED -25, SILENT -12 | FAIL: B1 INTACT -54, CROSS -51, SILENT -20; REMOVED +20 |
| 10, public | FAIL: B1 INTACT -18, CROSS -14, SILENT -10; B2 INTACT 5 vs 1, CROSS 5 vs 2, SILENT 2 vs 1 | FAIL: B1 INTACT -43, CROSS -38, SILENT -31; B2 INTACT 5 vs 1, CROSS 5 vs 2, SILENT 2 vs 1 (OUTSIDE errors 8 vs 547) |
| 9, A119b | PASS | FAIL: B2 CROSS 3 vs 0 |
| 8, A119b | FAIL: B1 INTACT -19 | FAIL: B1 INTACT -41, REMOVED -74, SILENT -50 |

- **Rebuilding the families changes no conclusion.** In every rebuilt set the integrated pipeline still has fewer errors in
  every world and more correct reads wherever the catalog is incomplete (OUTSIDE, BRANCH), and still loses correct reads on
  catalogued copies.
- **Single-read bars move.** For example, rebuilt set 9's 3 CROSS errors are reads of the newly added OR4F17 that carry OR4F5's
  annotated structure; they are near-identical to OR4F5.

## 9. Independent review of this branch and round 9 (2026-10-16)

No critical finding: the round-9 FAIL reproduces exactly. The corrections:

- **Public-profile errors (round 9):** they come from O1 dropping expressed copies with 4-5 reads, under its fixed 5-read
  support rule (PGA.2 had 2 uniquely placed reads). They were not an O2 effect. Recorded in the round-9 document and register
  row 1307.
- **Cost:** the per-world memory figures (about 1 GB) leave out the one-time local reference build: 38-75 s at about 14 GB, once
  per family set. A new family set still peaks near 14 GB (19.4 GB before); every later run of that set stays under 2 GB per
  world.
- **Read profiles, v2 (`profiles/*_v2.json`; v1 kept, and reproduces its reads byte for byte):**
  - **errors and donor alleles from every MAPQ-60 read** (`calibrate.py errors`). Single-isoform genes had under-counted the
    donor's het alleles:
    - A119b: substitution errors 0.41, het 0.58, hom 0.39 per kb, 1.38 mismatches per kb in all, against 0.89 in v1;
    - the simulator gives 1.33 per kb against 1.42 in real reads (v1 gave 0.90);
  - **expression spread** over copies with reads, less Poisson counting noise: A119b 0.86, public 1.03 (v1 1.16);
  - **truncation** redrawn while it would leave < 200 bases. v1 fell back to full length for 17-61% of short transcripts' reads.
- **Not simulated: the public library's soft clips (12-14% of bases)** are the Tn5 mosaic-end adapter (`AGATGTGTATAAGAGACAG`),
  left untrimmed at both ends of 99% of reads.
  - That is a preprocessing gap, not read noise. Such reads must be adapter-trimmed before any O2 that compares a whole read
    with paths.
  - A119b is trimmed: 0.3% clipped.
- **Local vs whole-genome mapping on realistic reads** (rebuilt set 10, A119b INTACT):
  - 99.87% of unique placements are identical; the 2 that differ are the same copy with a different first exon at the read's
    5' end;
  - 96.3% of reads have identical alignment sets;
  - every moved MAPQ-0 primary moved to an equal-score alignment.
- **Wording:**
  - the public library is called "held back" above, but it shaped v2's development;
  - its BAM was mapped without `-uf`;
  - the real-reads table mixes medians (de, NM) with pooled rates;
  - "complete catalogs: no wrong call" holds for sets 8-9, not round 9 (the current method misplaces FCGR1.0 reads);
  - a read's MAPQ-0 tie choice also depends on its name, so mapping by sequence moves tie choices (bars unchanged).
- **Numbers that came from scratch scripts:**
  - "51 of 55 / 45 of 1,154" (per-read genome check) was counted on the v1 REMOVED worlds: set 9 A119b, public and idealized,
    and set 8 A119b. The reviewer's 46/47 and 45/786 come from other worlds.
  - The MAPQ-60 fallback numbers (17/329, 159/220) were counted on the v2 development worlds; no output was saved.
- **Code fixes, each tested:**
  - `realworld.py` loads `o3_fewread/world.py` by path; a same-named `world.py` shadowed it in mixed test runs;
  - `fastset.py` refuses a local reference built for another family set;
  - `localref.py` lifts a same-contig mate.

**Post hoc runs added after the review** (integrated v2, per-read scoring):

| run | catalogued copies, correct vs current | errors | verdict |
|---|---|---|---|
| set 10, idealized reads (review item 6) | INTACT +8, CROSS +22, SILENT -1 | REMOVED 3 vs 1 | FAIL by single reads |
| set 10, A119b v1 (round 9) | INTACT -40, CROSS -37, SILENT -12 | fewer everywhere | FAIL |
| rebuilt set 10, A119b v2 | INTACT -63, CROSS -61, SILENT -22 | fewer everywhere | FAIL |
| rebuilt set 10, public v2 | INTACT -41, CROSS -41, SILENT -27 | INTACT 3 vs 2, SILENT 3 vs 2 | FAIL |

- **The effect of realistic reads is now measured on the held-out set itself.** Idealized reads give the integrated pipeline a
  small gain on catalogued copies; both calibrated libraries turn it into a 1-4% loss.
- **The corrected expression spread does what the review predicted.** Under public v2, O1 drops 0-2 expressed copies per world,
  against 2-4 under v1.
- **The pattern does not change.** In every realistic run, the integrated pipeline has far fewer errors wherever the catalog is
  incomplete (OUTSIDE: 3 vs 407-547) and loses correct reads on catalogued copies.

**Open:**
- O1's fixed 5-read rule for keeping a copy loses low-expressed copies under realistic expression. A rule on uniquely placed
  reads would keep them; it is the next O1 candidate.
- O2's divergence gate on realistic reads remains the catalogued-copy weakness.

## 10. Closing the remaining gaps (after the goal check, 2026-10-16)

**O1 keeps a copy with one uniquely placed read of its own** (`o1seed.py`). The 5-read primary count stays, so an identical
copy whose reads all tie is still kept. Over 138 worlds:

| profile | expressed copies dropped: 5-read rule / new rule | silent copies kept: 5-read rule / new rule |
|---|---|---|
| public | 68 / 14 | 35 / 35 |
| A119b | 11 / 11 | 82 / 84 |
| idealized | 2 / 2 | 22 / 22 |

The 14 expressed copies still dropped have only MAPQ-0 (tied) reads, so no read can single them out. End-to-end effect, on the same reads (rebuilt set 10, public v2 profile, integrated v2):
- **Expressed copies dropped:** 1-2 per world -> 0 in every world. The 8 silent copies are still dropped.
- **Integrated errors:** 3 -> 0 in each of INTACT, CROSS and SILENT, so B2 (errors <= current) now holds in all six worlds.
- **Correct reads:** +2 to +4 per world.
- **B1 still fails** on catalogued copies (INTACT -39, CROSS -39, SILENT -23): O2's divergence gate, not O1.

**Full families: faster and cheaper.** All of the following are exact: outputs are identical with and without them.
- **Family graphs shared** between worlds with one catalog (`O3_GRAPH_CACHE`).
- **O3 arm I only** (`detect.py --arms I`): the legacy arms' alignment and clustering are skipped; arm I's calls are
  identical.
- **A world's families in parallel** (`O3_PROCS`): identical with 1 or 2 processes, and to the previous code.
- **O2 scoring pruned** (gap-compressed distance): byte-identical.
- **4 threads for steps that run alone:** the shared mapping, REMOVED's mapping, and copy_assign under the heavy lock.
  copy_assign's output is identical at 2 and 4 threads.

| uncapped set | before (no cache) | cache + arm I only | + parallel families, 4-thread lone steps |
|---|---|---|---|
| 2 (69 copies) | 1,461 s, O3 1,215 s | 1,316 s, O3 770 s, peak 1.5 GB | **967 s**, O3 529 s, shared mapping 105 s (was 204 s), peak 2.2 GB, PASS |
| 9 (116 copies) | 1,145 s, O3 494 s | 1,027 s, O3 214 s, peak 5.0 GB (copy_assign alone) | |

Times exclude the one-time local reference build. All runs PASS their bars post hoc (A119b v1 profile).
- **Set 2: 1.5x faster than the first full-family run, O3 2.3x faster.**
- **What remains:**
  - copy_assign, the comparator rather than the integrated pipeline: about 60 s per world;
  - O2's path alignment: minimap2 of every hosted read against every path, about 80 s per world;
  - arm I's pooled per-host clustering: about 90 s per world.
- **Peak memory** stays under 2.2 GB per set once the local reference exists; set 9's largest copy_assign reaches 5 GB, alone
  under the heavy lock.

**Untrimmed-library adapters are now simulated and trimmed:**
- `calibrate.py adapters` measures them: the public library's 5' end is Tn5 mosaic end + G run + SMART primer (48 bp), its 3'
  end polyA + `GCATCG` + mosaic end (40 bp), on 98.8% of reads.
- Profile v3 (`public_testis_v3.json`) attaches them at the calibrated clip lengths, with sequencing errors.
- `trim.py` / `fastset.py --trim-adapters` removes them before mapping.

**The check (rebuilt set 10, public v3, integrated v2, per-read scoring):**

| reads | integrated, INTACT: correct / errors | current, INTACT: correct / errors | integrated O2 calls `unexplained` |
|---|---|---|---|
| adapters left on | 0 / 1,874 | 1,145 / 0 | 2,094 of 2,094 reads |
| adapters trimmed (`--trim-adapters`) | 1,334 / 0 | 1,368 / 2 | as without adapters |
| no adapters (v2, same O1 rule) | 1,327 / 0 | 1,366 / 2 | |

- **Left on, the integrated pipeline collapses.**
  - O2 calls every read unexplained.
  - O3 makes 8 false absent-copy calls in INTACT, a world with no absent copy, and sweeps the present copies' reads into them.
  - The current method aligns to the genome and clips the adapters, so it is barely affected.
- **Trimming restores it to the adapter-free level.** Every read is trimmed, removing 85.5 bases on average, and no adapter
  core is left at either end.
- **Rule:** an untrimmed library needs `--trim-adapters` (or lima/refine) before the integrated pipeline.

## 11. Integrated v3: the catalogued-read loss and the false calls (after the second goal check, 2026-10-16)

Development sets 8, 9 and 10 (all spent), simulated fresh with the v2 profiles and the new O1 rule. Every run is scored per
read (`--per-read-detectability`); "dev runs" are the 6 set x profile runs.
- Runs: `/mnt/linuxdisk/tmp/o3_graph/fast/fbdev/` and `fbdev_exact/` (`--union-exact`).
- Error lists: scratch scripts `errs.py`, `vs.py`, `why.py`, `crosserr.py` (session scratchpad), summarised here.

**O2: the genome fallback (`integl`).** It applies to the reads O2's paths leave `unexplained`.
- **Scoring:** every primary and secondary genome alignment, as mismatches + indel events.
- **What refuses:** an insertion or clip >= 30 bp (that is O3's evidence), and the read's best alignment lying off the
  family's seeded copies.
- **Placement:** within delta, at margin 2.
- **It removes the round-9 loss.** Correct reads, fallback minus no fallback: +29 to +90 per world (set 8 A119b INTACT
  2,132 -> 2,196), with no new error outside BRANCH. In BRANCH it adds +5 to +16 errors, reads of near-identical absent
  copies; current's BRANCH errors are 46-225.
- **Three other modes were measured on the same runs and rejected:**
  - `combined`: the genome on every read, with the guard. Errors drop (OUTSIDE 76 -> 0), but correct reads collapse
    where arrays have identical uncatalogued units (set 9 INTACT 935 -> 534).
  - `seeded`: the genome on explained reads, without the guard. It removes CROSS's structure-only errors, but INTACT loses
    13-47 reads.
  - `resolve`: the genome on O2's tied and thin reads. It adds 0-1 read per world.

**O3: three false-call mechanisms, found read by read.**
1. **The mosaic rule's nearest copy took the host's base wherever a single read covered a position** (a position covered by
   one read is no variable site). The cluster then disagreed with its own copy:
   - in set 9 LCE1, a cluster of LCE1.1's reads (one carrying a donor allele) was 0.025 from LCE1.1, so it became a false
     mosaic;
   - with the reads' own bases (`cluster_bases`, f080a30f), the distance is 0.0008.
2. **Mixed clusters:** reads of several present copies clustered together, a false mosaic (set 8 FAM138, once the
   consensus is honest).
3. **Present reads swept into a true call:** 2 USP17L reads in set 9 REMOVED.
- **Member attribution (`--member-attribution`, 12f8f2af) fixes all three.** An absent-copy PSV call keeps only the reads its
  consensus explains better than every single catalogued copy, by O2's margin of 2.
  - O3 errors fall to 0 in INTACT, CROSS and SILENT;
  - REMOVED's correct reads change by -2 to +13 per run.

**O2: the CROSS errors are failed projections.** All 17 CROSS errors of set 8 (both profiles) are reads of two structures
that O2's union never gave to the true copy:
- **WASH T0075974's last exon is 6 bp.** The spliced aligner cannot place it, so the projection had 4 exons, not 5, and was
  rejected.
- **DDX11L T0075982 projected,** but was merged with a structure within 10 bp of the copy's own.
- **`--union-exact` (36a430f8) fixes both.** A projection may lack a terminal exon < 30 bp, and structures merge only when
  their junctions match exactly. CROSS errors: 14 -> 1 and 10 -> 1.
- The simulator's projection (liftover) and its defaults are unchanged, so the test stays independent of O2's projection.

**v3 (fallback + attribution + `--union-exact`) against the current method, correct / errors:**

| world | set 8 A119b | set 8 public | set 9 A119b | set 9 public | set 10 A119b | set 10 public |
|---|---|---|---|---|---|---|
| INTACT | 2,191/0 vs 2,167/0 | 2,250/0 vs 2,153/1 | 935/0 vs 907/0 | 873/0 vs 828/0 | **1,136/0 vs 1,138/2** | 1,057/0 vs 1,051/2 |
| CROSS | 2,708/1 vs 2,551/8 | 2,788/1 vs 2,558/9 | 1,024/0 vs 1,001/0 | 963/0 vs 917/0 | 1,256/0 vs 1,248/3 | 1,159/0 vs 1,146/3 |
| SILENT | 1,917/0 vs 1,860/0 | 2,085/0 vs 1,909/1 | 678/0 vs 644/0 | 628/0 vs 587/0 | 702/0 vs 689/2 | 613/0 vs 602/2 |
| REMOVED | 2,137/1 vs 2,059/12 | 2,158/0 vs 2,072/7 | 798/0 vs 782/0 | 713/0 vs 698/0 | 1,049/2 vs 1,031/10 | 1,021/2 vs 963/4 |
| OUTSIDE | 2,232/76 vs 1,979/150 | 2,273/62 vs 2,000/120 | 760/23 vs 674/152 | 681/19 vs 605/123 | 1,095/3 vs 863/169 | 1,028/3 vs 789/167 |
| BRANCH | 2,505/104 vs 2,203/137 | 2,509/166 vs 2,176/191 | 1,196/36 vs 907/60 | 1,080/30 vs 828/46 | 1,320/66 vs 1,138/225 | 1,192/81 vs 1,051/200 |
| verdict | PASS | PASS | PASS | PASS | FAIL (B1 INTACT, -2) | PASS |

- **The remaining shortfall is structural.** Set 10 A119b INTACT loses 15 short reads (184-484 bp) whose 2-5 sequencing edits
  exceed delta's share of their length.
  - A looser gate would place reads that per-read scoring counts as detectable in REMOVED.
  - The current method makes 0 errors in REMOVED in set 9.
- **Faster and cheaper, every change exact (outputs byte-identical):**
  - **O2 on a full set 2 world: 83 -> 49 s.**
    - The clipped-end edit distance uses Myers/Hyyro bit vectors; it equals the DP on 3,000+ random cases.
    - CIGAR events are parsed by regex.
    - Families run in parallel (`O2_PROCS`, `--o2-procs 2`).
  - **O3 arm I on full set 2 worlds: 55-59 -> 32 s.**
    - Isoform maps are built lazily.
    - PSV clustering keeps its consensus current instead of recomputing it per comparison, and computes linkage on read bit
      sets; it equals the reference implementation on 60 random hosts each.
  - **O3's graph:** building is cached (section 10). The remaining O3 time is GraphAligner on the hosted reads (about a
    third) and the per-host pool alignments.
- **Round 10 pre-registers v3** on an eleventh set: `docs/PREREG_o3_graph_round10_2026-10-16.md`.

## 12. After round 10: v4 (2026-10-16, evening)

Round 10 (`docs/PREREG_o3_graph_round10_2026-10-16.md`, register 1308) failed. Its independent review corrected the first
failure's mechanism. Development since then uses round 10's worlds, now spent, plus sets 8-10, all with the v2 profiles; that
is nine set x profile runs including round 10's untrimmed-library run.

**1. Attribution no longer keeps a read that no catalogued copy competes for** (88c445d9). Round 10's false TAF11L call came
from the host window:
- TAF11L.44's +-2 kb window holds 225 bp of TAF11L.43, and TAF11L.33 reads aligned to that fragment;
- no copy's map covers it, so attribution kept them by default.

038890fd (`--member-attribution-delta`) was built on the wrong mechanism and is not used.

**2. Copies within delta of each other are undecidable** (`--twins`; c3113f61, 7b8092d1):
- O1 marks a new copy within delta of a kept curated copy as its twin (genomic spans, asm20 `de`). O2 then lets the twin
  take a read only by the margin.
- O1 keeps a curated copy that its support rule would drop when it lies within delta of a kept one. Round 10's PRR20.1 had
  all its reads tied, a MAPQ-0 tie-break dropped it, and its reads went to the twin as errors.

**3. Delta is error-free, reads are not** (`--error-floor`; f2789fdd, 2afdbc3d).
- Delta (0.00958) is the p99 exonic divergence between KB3781's haplotype assemblies (register 1199).
- Every read-level comparison with it compared an error-inclusive `de`. The world's sequencing-error floor is now subtracted
  first, as the median `de` of its MAPQ-60 primaries (round 10: A119b 0.0027-0.0031, public 0.0015-0.0018). It applies to:
  - O2's `unexplained` gate (paths and genome fallback): explained up to (delta + floor) x read length;
  - O3's genome check, call-level and per member;
  - the scorer's per-read detectability: a REMOVED read is detectable when its `de` minus the floor exceeds delta, for both
    systems.
- It removes both remaining failure mechanisms:
  - set 10's short reads whose sequencing edits exceeded delta's share of their length;
  - array clusters whose members' `de` (0.0103-0.0118) is sequencing noise above a near-identical copy.

**4. O1-hosted O3 (`--o1-hosts`, f45e0fe5) is measured and not adopted.** It removes the integrated pipeline's last use of
copy_assign's run, but loses REMOVED and BRANCH recall in set 9's arrays (REMOVED -5 and -20, BRANCH -55 and -53). It
passed 6 of the 9 runs.

**v4 = v3 + 1 + 2 + 3, against the current method.** Correct / errors, floor-corrected per-read scoring:

| world | r10 A119b | r10 public | set 8 A119b | set 8 public | set 9 A119b | set 9 public | set 10 A119b | set 10 public |
|---|---|---|---|---|---|---|---|---|
| INTACT | 737/0 vs 645/0 | 654/0 vs 576/0 | 2,221/0 vs 2,167/0 | 2,254/0 vs 2,153/1 | 944/0 vs 907/0 | 876/0 vs 828/0 | 1,146/0 vs 1,138/2 | 1,059/0 vs 1,051/2 |
| CROSS | 779/0 vs 689/0 | 694/0 vs 616/0 | 2,742/1 vs 2,551/8 | 2,792/1 vs 2,558/9 | 1,037/0 vs 1,001/0 | 966/0 vs 917/0 | 1,265/1 vs 1,248/3 | 1,161/0 vs 1,146/3 |
| OUTSIDE | 982/83 vs 828/142 | 919/72 vs 802/115 | 2,264/82 vs 1,979/150 | 2,276/77 vs 2,000/120 | 777/43 vs 674/152 | 690/34 vs 605/123 | 1,105/8 vs 863/169 | 1,030/5 vs 789/167 |
| REMOVED | 969/0 vs 879/2 | 919/0 vs 843/3 | 2,126/0 vs 2,059/7 | 2,159/0 vs 2,072/4 | 802/0 vs 782/0 | 710/0 vs 698/0 | 1,032/1 vs 1,031/4 | 1,018/0 vs 963/4 |
| BRANCH | 1,043/59 vs 682/154 | 946/49 vs 576/139 | 2,536/110 vs 2,203/137 | 2,513/168 vs 2,176/191 | 1,205/36 vs 907/60 | 1,083/31 vs 828/46 | 1,329/74 vs 1,138/225 | 1,194/84 vs 1,051/200 |
| SILENT | 476/0 vs 464/0 | 439/0 vs 424/0 | 1,929/0 vs 1,860/0 | 2,023/0 vs 1,909/1 | 676/0 vs 644/0 | 629/0 vs 587/0 | 700/0 vs 689/2 | 609/0 vs 602/2 |

- **Floor-corrected scoring: PASS in all 9 runs,** including the untrimmed-library one (not tabulated).
- **Round 10's scoring: PASS in 5.** The 4 failures are all REMOVED B2: the wider gate places reads whose `de` lies
  between delta and delta + floor, which that rule counts as errors.
  - set 8 public 9 vs 7;
  - set 9 A119b 1 vs 0;
  - set 10 A119b 33 vs 10;
  - set 10 public 13 vs 4.
- Runs: `/mnt/linuxdisk/tmp/o3_graph/fast/r11dev6/` (`bars_floor.txt`, `bars_r10rule.txt`), from r11dev3's O1 output.

## 13. Rounds 11-13: v4 passes; v5 and v5' fail (2026-10-16/17)

| round | pipeline | families | verdict | what decided it |
|---|---|---|---|---|
| 11 (r1309) | v4 | 8, 37 copies | **PASS, 12/12** | the INTACT/CROSS/SILENT margins are the current method's own O3 artefact (review) |
| 12 (r1310) | v5 = v4 + O1-hosted O3 (`--o1-hosts`, tiered) | 8, 104 copies (80-copy DUX4) | FAIL, both profiles | O2: 2 margin-2 coincidences in DUX4, 5 structural errors (retained introns absorbed by a paralog's path as one gap-compressed edit; a read DNA-identical at two copies). O3 cost: 791-1,355 s per set |
| 13 (r1311) | v5' = v5 + `--branch-genome-novel` | 8, 89 copies (43-copy OR7E, 3 retrocopy families) | FAIL (A119b); public 6/6 | one false O3 call on three haplotype-2 reads of MRGPRX.1 in four worlds: a mosaic call whose divergence is mostly one read's sequencing indels, kept by attribution's margin of 2 that two donor alleles meet (review). No O2 error on catalogued reads |

- **v5 makes the integrated pipeline independent of the current method.** O3's hosts come from O1's seeded copies, tiered:
  primary alignment on a seeded copy, then a secondary, then a +-1 kb window.
- **v5' makes it cheap on large arrays.** O3's branch arm clusters only reads whose own genome alignment leaves >= 50 bases
  unexplained. On round 12's worlds this cut O3's longest world from 699 s to <= 78 s. In round 13 the integrated steps took
  223 and 216 s per set of six worlds; the current method took 188 and 183 s.
- **v6 (`--distance gapstruct`) is not adopted.** Gap compression applied only below 30 bp still fails round 12, and it costs
  errors and correct reads in sets 8 and 10.
- **Open after round 13:**
  - round 12's structural O2 errors, which round 13's families did not provoke;
  - the allele-versus-copy confound in O3's genome check (section 14).

## 14. After round 13: v8, then v9 under a diploid donor (2026-10-17)

Round 13's review traced the false call to three things together:
- a mosaic distance that counts one long read's sequencing indels at positions only it covers;
- attribution's fixed margin of 2, which two donor alleles meet exactly;
- the per-read genome check, the last filter.

Development uses 16 spent runs: sets 8-10 and rounds 10-13, both profiles, plus the two untrimmed-library runs. v5' passes 13
of them. The rule, fixed before any result was read: a lever is adopted only if it removes round 13's false call or reduces
errors, and turns none of the 13 passing runs into a failure (`fast/r14dev/<lever>/`).

| lever | round 13 A119b | runs passing | turned a passing run into a failure | cost |
|---|---|---|---|---|
| O3 `--allele-discount` v7f: donor alleles leave the genome check | fixed | 13 | set 10 A119b (REMOVED 1,032 -> 1,024 vs 1,031) | low-divergence removed copies whose reads O2 partly assigns |
| O3 `--member-attribution-delta`: attribution margin max(2, delta x length) | fixed | 13 | set 10 A119b (1,027) | REMOVED -1 to -34 in 11 runs |
| **O3 `--mosaic-shared`: the mosaic distance only over positions >= 2 members cover** | **fixed** | **14** | **none** | **none: identical or better in every run** (set 10 A119b REMOVED +49, round 11 A119b +11) |
| **O2 `--o2-genome veto`: a path call overruled where the genome places the read confidently** | (no O2 error) | **15** with `--mosaic-shared` | **none** | CROSS -1 and -2 in set 8 |

**The allele discount took four versions,** each lost on a different removed copy:
- v7 learned alleles from every explained read and lost CTSLP.1. CTSLP.1's own reads, which O2 explains on its paralog, carry
  its PSVs.
- v7b counted candidates per seeded locus, and CTSLP.1's PSVs were diluted over its paralog's length.
- v7c grouped co-carried candidates and lost CXADRP.1 and AGAP.5.
- v7f learned only from path-assigned reads, taking a variant as an allele when its assigned carriers outnumber the calls'
  members. It still lost two small calls of within-delta copies that are detectable only through high-error reads.

It stays opt-in, not adopted. A population cannot always tell a present copy's allele from a low-divergence missing copy's
PSV.

**The genome veto took three versions:**
- v1 scored alignments with `locus_distance`, which leaves a clip under 30 bp free. A paralog's alignment that clipped a read's
  differing end then beat the read's own locus, and true KRTAP2.1 calls were moved.
- v2 charged clipped ends and long insertions in full (`veto_distance`). Now a short terminal exon the aligner missed at the
  true locus moved a correct DDX11L call (primary MAPQ 1).
- v3 adds a guard: the read's primary alignment must lie on the genome's winner at MAPQ >= 30.
  - Audit of v2 over 13 runs: every wrong veto had MAPQ 0-1, and every veto that fixed an error had MAPQ 53-60.
  - It fixes round 12's four CDRT15 retained-intron errors (round 12 public now passes), cuts OUTSIDE errors in 8 runs, and
    adds no error anywhere.

**v8 = v5' + `--mosaic-shared` + `--o2-genome veto`: 15 of 16 development runs pass.**
- The one failure is round 12 A119b, with 2 DUX4 margin-2 coincidences and 1 KRTAP2 read that is DNA-identical at two copies,
  in INTACT, CROSS and SILENT. The current method abstains on all of them, and no sequence rule separates them.

**The simulator** (`--sim-het-draw`, round 13 review I1): each copy's het and hom site counts are drawn (Binomial over its
exonic bases) instead of rounded. With rounding, the public donor was effectively haploid: 86 of 89 copies had no het site.

**Under drawn het sites, v8 meets the allele confound again.** v8 was run on round 13's families re-simulated with
`--sim-het-draw` (`fast/r14dev/hetdraw13/`, both profiles).
- **Public passes.**
- **A119b fails INTACT 5 vs 0 and CROSS 6 vs 2.**
  - Every error is an O3 call on FTH1P.6's haplotype-2 reads.
  - FTH1P.6 is a 553-bp retrocopy that drew 3 het sites within 131 bp. Its short haplotype-2 reads (203-258 bp, MAPQ 60)
    carry all three, and at that length 3 alleles plus their sequencing errors clear delta plus the floor.
  - On FTH1P.4's window, the 5 shortest form a PSV call with 3 linked unexplained alts.
  - A copy-level view puts the alleles at 3 / 553 = 0.0054, under delta: a haplotype, not a copy.
- **v8's rules do not reach it.**
  - The mosaic fix does not apply, since this is a linked-PSV call.
  - The genome check reads each read's local divergence.
  - The genome veto is O2.

**O3 `--haplotype-check`.** A PSV call is measured as the haplotype it is, not as its members alone. Development took three
versions:
- **v1, genome coordinates.** The substitutions >= 2 members share, then every primary carrying one, over the positions >= 2
  carriers cover.
  - It removes MRGPRX.1 and FTH1P.6, and passes 15 of 16 development runs.
  - It misses an array. On round 10's families re-simulated with drawn het sites, 3 haplotype-2 reads of TAF11L.10 (a 45-copy
    array) are called on TAF11L.44's window. Their MAPQ-0 primaries lie on three different units, so in genome coordinates they
    share nothing.
- **v2, the host window.** The family's explained reads are aligned to the call's window, those carrying the call's alts are
  added, and the extended cluster's nearest-copy distance is taken over positions >= 2 of them cover.
  - Selecting carriers by any alt pulled in 861 reads of 12 FTH1P copies, since a mosaic cluster's alts are other copies'
    bases. The mixture looked as distant as the call, and true REMOVED calls were lost.
- **v3 = v2 with carriers selected by the call's unexplained alts only** (the alts no copy's map carries).
  - FTH1P.6: 61 FTH1P.6 carriers, distance 0.0019.
  - TAF11L.10's call is removed too.
  - Every true call in those worlds' REMOVED is unchanged.

**v9 = v8 + `--haplotype-check` (v3).**
- **On the 16 spent runs: 15 pass**, as with v8. The failure is again round 12 A119b. v9's O3 rules leave 13 runs identical to
  v5' and improve 3: round 13 A119b's false call removed; REMOVED +2, +11 (round 11 A119b) and +45 (set 10 A119b).
- **On 8 family sets re-simulated with the diploid donor: 8 of 8 pass.** These are rounds 9, 10, 11 and 13's families, both
  profiles (`fast/r14dev/hetdraw*`, `hetfresh*`, `hetc_v9c`).
- **Cost:** the check aligns the family's explained reads to each host window with a candidate call. O3's longest world on
  round 13 A119b went from 26 to 58 s.
- Round 14 pre-registers v9 with the diploid donor (`docs/PREREG_o3_graph_round14_2026-10-17.md`).

## 15. After round 14: v10, the O2 proof, and the cost of full families (2026-10-17)

**v10 = v9 + O3 `--haplotype-check-genome`.** This brings back the genome-coordinate haplotype check beside the window one
(f867ec7c).
- **On the 16 spent runs: 15 pass,** as with v9. The failure is again round 12 A119b. Per world, the results are identical to
  v9 or have fewer OUTSIDE errors (round 12 public: 130 -> 106).
- **On round 14's 5 seeds and 8 diploid-donor sets: all 18 runs pass.** Under v9, A119b passed only 2 of the 4 extra seeds.
- **These sets are now spent.** Round 15 is the held-out test.

**NPIP as one family** (families.py group `mainstay`): NPIPA and NPIPB copies link into one component under the locus rule
(23 of 23), and split in two they placed NPIPB reads on NPIPA copies. On the merged mainstay set (seed 306), both profiles
pass every bar.

**The O2 proof** (`docs/PREREG_o2_proof_2026-10-17.md`, register 1314) passed in every cell.
- The integrated O2 made 0 wrong assignments in 22,631, on 8 held-out near-identical families.
- MAPQ >= 1 is the strong baseline: it is more powerful on weak reads and makes 55 errors.
- The same runs' end-to-end OUTSIDE worlds failed: O1 did not seed near-identical uncatalogued twins (C4, HSFY), and O2's twin
  rule gave twin reads to the curated copy (VCY).

**Where a full family's time goes** (mainstay CROSS world, O3_TIMING, one world):
- **vg graph construction was 121 of O3's 158 s with a cold cache.** The cache is keyed by the graph's inputs, and each world's
  O1 catalog differs, so most worlds miss.
- vg msga is not deterministic either: two builds of NPIP's graph differ in copy17's paths.
- **O2:** 32 s, of which minimap2 of 1,829 NPIP reads against 307 paths is 23 s.

**`--branch-linear` (O3, 44798db2): the branch step without the family graph.**
- Candidates are spliced-aligned (minimap2 splice:hq) to the family's loci as linear sequences.
- A read base is novel only where NO locus alignment explains it (`novel_union`). That is how the graph explains NPIP.14's
  125-253-bp insertions: other copies' alignments hold them.
- **Mainstay CROSS world:** calls identical to the graph version; O3 158 s -> 34 s.
- **On all 34 spent runs** (16 rounding runs, round 14's 10 seed runs, 8 diploid-donor sets):
  - every verdict is v10's;
  - 14 runs have identical bars;
  - 19 gain 1-9 correct BRANCH reads (two REMOVED worlds +1), with errors equal or one fewer.
- **Measured and not adopted:** `--branch-genome-segments`, which clusters the genome alignment's own unexplained pieces. It
  called 32 present NPIP reads on the mainstay CROSS world, whose pieces another NPIP copy explains.

**Exact caches for a set's worlds (`fastset --aln-cache`):**
- **O2's path alignments** are cached per read sequence (`O2_ALN_CACHE`, f2e59107). minimap2 aligns each read on its own;
  equal-score ties are broken by a hash of the read's name, which changes only the tp / s1 / s2 tags, never NM, the CIGAR or
  the coordinates O2's distances read.
  - **Mainstay, against the original tables:** INTACT cold, then CROSS warmed by INTACT, then CROSS warm. All three are
    byte-identical.
- **O3's catalog-isoform-to-host alignments** (`O3_ALN_CACHE`, 7839138c) are keyed by window, query names and sequences.
- **Not cached:** O3's read-to-host alignments. They keep only minimap2's primary, so a name-hash tie could change which one
  is kept.

## 16. After round 15: the O2 margin (2026-10-17, evening)

**Round 15's errors, read by read.** Every one is a FAM90A read, and the current method makes none of them.
1. **Path-set asymmetry (most of them).**
   - FAM90A.45's read-derived paths (an exon skip) fit reads of its siblings, and the siblings' own paths lack that
     structure.
   - Each such read lay 16 edits from copy 45's paths and 64-378 from every sibling's. On the genome, the siblings hold it at
     5 and copy 45 at 16.
   - Union projection adds a paralog's isoform to a copy only for an intron chain new to that copy, so read-derived
     structures and longer first exons can stay on one copy.
   - The fixed margin of 2 edits rewarded a copy for having the read's path, not for its sequence.
2. **A locus-fallback call on clipped bases** (FAM90A.51's read): copy 45's genome alignment clipped the 11 bases that tell
   the copies apart. locus_distance leaves a clip under 30 bases free, so copy 45 won by 2.
3. **A call no rival's paths could answer** (FAM90A.12's 201-bp read): only copy 11's paths aligned it, and eight loci hold
   it exactly.

**Round 12 A119b's standing failure** (the one development run v10 and v11 never passed) adds two more mechanisms:
4. **Indel margins against many rivals** (DUX4, 80 copies).
   - Two path calls won by 2 edits that were two 1-bp homopolymer indels, against 63 rivals.
   - copy_assign ties these reads: its PSVs leave indels out by default (`indel_psv: false`).
5. **A retained intron** (KRTAP2.0).
   - The simulator keeps an intron in 5.8% of reads (the A119b profile, from real reads).
   - Copy 1's long annotated first exon absorbs the intron as one gap-compressed insertion (2), and copy 0's paths clip it
     (85). Both loci hold the read exactly.

**The error model behind mechanism 4** (the v2 profiles, from real FLNC reads):
- substitutions 4.1e-4 per base;
- indel events 2.1e-3 per base, 63% of them inside homopolymer runs of 3 or more;
- donor alleles 1.0e-3 per base.

A homopolymer indel is about 20 times likelier than any one substitution, yet the gap-compressed distance counted both as one
edit.

**Four rules, each in O2's genome step** (`o2.integv`; O3's input, `o2.integ`, is unchanged). Each was measured on all 46
spent development runs (276 worlds: rounds 10-15, sets 8-10, round 14's seeds, the diploid-donor sets and the O2 proof's
cells), against v11's table, with the production code (`o2graph.py run W integ|integv`; the integ step with
`--o2-significance` reproduced every original table exactly). Adoption rule: no world that passed may fail.

| rule | mechanism | correct | errors | verdicts |
|---|---|---|---|---|
| `--veto-abstain` | 1 | 0 | -95 (27 worlds) | round 15 A119b 503 FAIL -> PASS |
| `--locus-pairwise` (with veto-abstain) | 2 | -5 | -6 | round 15 A119b 502 FAIL -> PASS |
| `--o2-significance` | 4 | -421 | -14 | none alone |
| `--structure-check` (with veto-abstain) | 3, 5 | -868 | -24 | fixes round 15 public 503 CROSS |
| **all four (v12)** | 1-5 | **-1,294 (-0.21%)** | **-139 (-2.1%)** | **36 -> 40 of 46 runs pass; no world newly fails** |

- **`--veto-abstain` (1bd6693a): an interval rule.**
  - A copy's genome cost lies between its locus_distance (a clip under 30 bp free) and its veto_distance (every clipped base
    paid).
  - A path call is abstained (`genome_contradicts`) when a rival's worst case beats the call's best case by the margin, and
    no confident single winner takes the read. The veto's MAPQ >= 30 rule still moves reads with one confident winner.
  - **The first form counted the clip-paid cost alone.** It abstained correct calls whose own alignment missed a short exon:
    DDX11L.6's 15-bp terminal exon, and DAZ.2's read end.
- **`--locus-clipped` (not adopted):** the same interval rule on locus-fallback calls.
  - It abstained correct calls wherever every copy's alignment clips the read's ends (HNRNPCL, TAF11L, FAM236).
  - It failed r14 A119b s205's B1.
- **`--locus-pairwise` (10208f8c).**
  - Against each rival, a copy pays its alignment's edits plus the read bases the rival aligns and it leaves out. Bases every
    copy clips are nobody's evidence.
  - The one read it loses (ANP.1's) is the ambiguity no alignment settles: the true copy's alignment clipped 14 bases that
    the rivals align with 6-9 edits.
- **`--o2-significance` (1bd6693a).**
  - O2's path distances, with each edit priced at -ln of its type's rate in the world's uniquely placed reads: a mismatch to
    one given base, a homopolymer indel per homopolymer run, another indel per base. The CIGAR's context is read against the
    path and genome sequences.
  - Round 12 A119b: 7.8 for a mismatch, 3.9 for a homopolymer indel, 7.3 for another indel.
  - A path call stands only when it beats every rival by ln(R / alpha), with R the rivals and alpha = 1e-3, as copy_assign's
    certificate divides alpha by its rivals.
  - **Calibration, round 12 A119b INTACT:** the gate took 39 calls, all DUX4 on two homopolymer indels (7.7 against a
    threshold of 11.0); 2 were wrong. That is 5% wrong, against the model's R e^-7.7 = 3%.
  - It adds 1-3 CPU seconds per family world (NPIP 1.4 s, GOLGA8 1.2 s).
- **`--structure-check` (10208f8c).**
  - A rival whose own locus (an alignment inside its padded span) explains the read better than its paths do, by the margin,
    lacks the read's structure among its paths.
  - The call is then compared with that rival on the genome, and abstained (`structure_only`) unless its best explanation
    wins there by the margin.
  - DAZ.2's rivals align 40-58 kb past their loci, so they are no evidence.
  - It also abstains correct calls that rest on annotation structure the genome supports by one edit or not at all. The
    simulator rewards those calls, since its reads come only from annotated isoforms. The rule keeps the pipeline's own
    margin of 2 rather than a tuned exception.

**The PSV question, as far as O2 goes.**
- A SUN-type difference (a substitution) now weighs 7.6-7.8 nats and a homopolymer indel 3.9-4.6. The CIGAR's edit types
  and the read's own error profile set the margin, where the fixed 2 edits treated every edit alike.
- `de` still flags absent-copy reads (section 5). Path alignments' mismatches come from NM less the indel bases (their
  CIGARs carry M); the genome's from =/X.

**v12 = v11 + all four.** Over the 276 worlds:
- v12 makes 602,563 correct calls and 6,478 errors, against v11's 603,857 and 6,617.
- The current method makes 556,431 and 24,141.
- **Four failing runs pass:** round 15 A119b 502 and 503, round 15 public 503, and round 12 A119b. Round 12 A119b is the one
  development run that v10 and v11 never passed.
- **The six still failing are the O2 proof's cells,** on their OUTSIDE twins (register 1315).
- **The thinnest B1 margins left are 3-6 reads** (integrated correct less current correct):
  - the O2 proof's public 402 INTACT and CROSS: 13 -> 3 and 15 -> 5;
  - het9 INTACT and SILENT: 6-8 -> 4-6;
  - set 10 public SILENT: 7 -> 5.
  B1 is where v12's cost could fail a held-out world.
- The structure check and the significance gate cost about 30-36 correct reads per error they remove. That is the price of
  holding errors at zero where the current method abstains.
- **Round 16 tests v12 on families it has not seen** (`docs/PREREG_o3_graph_round16_2026-10-17.md`).

**Round 16 (register 1318): FAIL on correctness (4 of 6 cells pass; A119b 601 and 602 fail); PASS on cost.**
- **B2 held in all 36 worlds.** The O2 margin removed every error where the current method makes none. Over the 36 worlds:
  1,833 errors against 8,796, and 147,583 correct calls against 137,426.
- **The three failing bars are B1, by 4-7 reads,** in A119b INTACT and SILENT of small families (3-7 copies).
- **The integrated pipeline's fixed thresholds are the cause, not v12's rules** (v11's genome step fails the cell as well):
  - `unexplained` short reads, whose edits exceed (delta + floor) x length;
  - `thin` margin-1 calls among few rivals, which copy_assign's certificate assigns.
- **The cost bar passed in CPU** (-23% and -30% on the mainstay set). On round 16's families the integrated stages are also
  faster on wall time (68-100 s against 118-138 s per cell).

**After round 16 (register 1319): v13 = v12 + `--struct-rescue --short-path-align`.**
- **What round 16's missing reads were:** reads the current method placed and v12 did not.
  - Short reads over the fixed divergence gate.
  - Thin margin-1 calls. Gap compression let the read's best path be an isoform lacking one of its exons, so the paralogous
    bases inside that exon went uncompared.
  - 150-200-bp reads with no path alignment at all: map-hifi's minimum DP score of 200 is out of their reach.
- **The adopted fixes, on 52 spent runs:**
  - Thin and tied reads are decided again on structure-aware distances (gapstruct: a long indel counts its length).
  - Reads with no path alignment are realigned with -s 40.
  - The result: 44 -> 46 runs pass, round 16's three failing worlds pass, no world newly fails; +4,289 correct, +95 errors in
    OUTSIDE and BRANCH.
- **Not adopted:**
  - A length-calibrated gate: it places reads of absent copies.
  - A significance-based margin: about one error per 1,000 rescued reads, where the current method makes none.
- **Round 17's families are name-free** (`seqfamilies.py`, user, 2026-10-17). They are every protein-coding gene, LOC genes
  included, linked by sequence identity (the locus rule) and by structure (the alignment carries one gene's junctions onto
  the other's introns).

**Round 17 (register 1320): FAIL on correctness (0 of 6 cells), PASS on cost.**
- **The setup:** name-free families (40, mostly pairs).
- **What held:** B1 in all 36 worlds; v13 made 124,264 correct calls against 82,008.
- **What failed:** OUTSIDE's B2 in every cell. When a pair's twin is uncatalogued and O1 does not seed it, O2 places the twin's
  reads on the one remaining copy (the O2 proof's open failure). The current method abstains on them.
- **Cost:** the integrated stages used 20-29% less CPU on the NPIP-sized set, and 43% less on round 17's cells.

**After round 17 (register 1321): v14 = v13 + O1 `--ref-copies`.** O1 searches the genome for copies the catalogue lacks (the
family rule on the world's homology-closed local reference); round 17's spent cells went 0 -> 5 passing.

**Round 18 (register 1322): FAIL on correctness (4 of 6 cells pass; both seed-803 cells fail), PASS on cost.**
- **What held:** B2 in all 36 worlds, OUTSIDE included (0-11 errors against 69-198). v14 made 130,722 correct calls against
  86,384, and 3,782 errors against 10,431.
- **What failed:** B1 in five worlds of seed 803, all from one family (OR4F21, OR4F16).
  - The genome search adds five OR4F unprocessed pseudogenes, 98.8-99.4% identical, which the simulator never expresses.
  - v14 abstains the reads tied with one, or one substitution from one. The current method, which compares catalogued copies
    only (`--no-as-tied-only`), places them on OR4F21, and in this simulation it is always right.
- **The same tie in another form failed `--twin-guard` (register 1321):** a read tied with an expressed uncatalogued twin and a
  read tied with a silent paralog are the same exact AS tie.
  - The current method wins the second kind because the simulation never expresses paralogs outside the catalogue.
  - Real A119b testis expresses these five pseudogenes: 41 primary reads, against 59 on the two coding copies.
- **Cost:** -3% and -13% CPU on the NPIP-sized set. O1's genome search now costs 141-147 CPU seconds there.

## 17. The goal's PSV question, answered by objective (2026-10-18)

The goal asks whether PSVs can include SUNs, and how the CIGAR, NM and de tags should count. The answers, from sections 4, 5
and 16 and rounds 16-18:

**O1, which copies exist.**
- **The CIGAR's introns are O1's structure evidence.** The genome search (`--ref-copies`, register 1321) takes a locus as a copy
  only when a catalogued transcript aligns to it at identity >= 0.90 over >= 0.8 with its exon chain carried: every N of the
  spliced alignment at one of the transcript's junctions.
- **A SUN is private only relative to the copies considered** (section 4). Completing the catalogue from the genome makes
  "private" mean private among the loci the family rule reaches. Round 18 showed the cost (register 1322): five OR4F
  pseudogenes made OR4F21's catalogue-private differences shared, and its reads tied.

**O2, which copy a read came from.**
- **SUNs are most of the evidence:** about 70% of the decisive PSV columns are SUN-type, the called copy's own SUN or a rival's
  SUN the read rules out (section 4). With a complete catalogue no call was wrong, with or without an own-SUN column.
- **The CIGAR's edit types set the margin** (`--o2-significance`, section 16). Each edit is priced at -ln of its rate in the
  world's uniquely placed reads:
  - a substitution 7.6-7.8 nats;
  - a homopolymer indel 3.9-4.6;
  - any other indel 7.3.
  - A call must beat every rival by ln(R / alpha). One SUN decides against one or two rivals; one homopolymer indel never.
- **Indels count as events, not bases** (section 5). NM counts a retained intron or a 4-bp homopolymer indel as many edits. A
  gap-compressed distance (mismatches + indel events, like de) recovered the calls NM lost.
  - A long indel (>= 30 bp) counts its length again (`--struct-rescue`), because it can hide the paralogous bases inside a
    skipped exon (SPAG.1, register 1319).
- **Soft clips bound the genome cost from both sides** (`--veto-abstain`): a clip under 30 bp is free at the lower bound and
  paid in full at the upper.

**O3, reads of copies the catalogue lacks.**
- **de flags them.** A read whose de, less the world's error floor, exceeds delta is detectable as an absent copy's read (the
  scorer's per-read detectability, section 5).
- **The read's own de cleans the calls:** present-copy reads inside a true absent-copy cluster sit within delta of their own
  locus (51 of 55; `--member-genome-check`).

**What it bought, measured on held-out families:** B2 (errors no more than the current method's) held in every world of rounds
16, 17 (except OUTSIDE, fixed by `--ref-copies`) and 18. Over round 18's 36 worlds: 3,782 errors against 10,431. The open bar
is B1 in worlds with silent paralogs (register 1322), which round 19 tests under the measured paralog expression.

