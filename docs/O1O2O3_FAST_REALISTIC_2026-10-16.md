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
