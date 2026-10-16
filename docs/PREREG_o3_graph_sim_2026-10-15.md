# O3 with a variation graph: three steps, tested on curated ideal scenarios (2026-10-15)

Pre-registration, written before any held-out world was built. Code: branch `o3-graph`, commit `883e54bc` (`bench/o3_graph/`).

**User goal (2026-10-15):** "proceed with your 3 step recommendation for giving the vg a bit more detection power; ... we
can use curated isoforms from refseq or cat to simulate ideal scenarios".

The recommendation followed the two-locus graph spike (`docs/PREREG_o3_family_best_2026-10-15.md`, Results, and memory).
It was: feed the graph into the current PSV-cluster method (`copy_assign --psv-clusters`, prototype
`bench/o3_unexplained/psvclust.py`) rather than build a separate detector.

## The three steps (as arms of `bench/o3_graph/detect.py`)

| arm | what it adds |
|---|---|
| CUR | the current method: `copy_assign --no-as-tied-only --psv-clusters`; a call = an `unexplained` cluster with est_div > delta |
| R0 | the prototype's per-host clustering against each copy's representative: CUR reproduced in Python |
| S1 | **step 1**: a cluster is explained against each copy's EVERY curated isoform, the isoform paths of the family's variation graph, not its one representative |
| S12 | S1 + **step 2**: branch calls. GraphAligner maps the family's hosted reads to its variation graph (`vg msga` of its copies' loci, every curated isoform added as a path with `vg rna`). A read with >= 50 bases on no graph path is a branch read: large insertions, dense error windows where the aligner forced novel sequence onto the graph, gaps, unaligned ends. Reads sharing a novel segment are linked; a group of >= 5 reads is a call. |
| S123 | S12 + **step 3**: a call is kept only when the genome does not explain its reads (see below) |
| S3 | R0 + step 3 alone |

Step 3 keeps a call only when the genome does not explain its reads:
- a PSV call: the median minimap2 `de` of its reads' primary alignments is > delta;
- a branch call: its reads' median unexplained bases on the genome (soft clips and dense error windows) are >= 50, or that
  median `de` is > delta.

delta = 0.00958, as in registers 1294-1296.

## Ideal scenarios (`bench/o3_graph/world.py`, famsim)

**Templates** (`templates.py`): curated CAT genes on CHM13 chr20 (`bakeoff/human_chr20/chr20_ref_cat.gtf`).
- Each has 4-14 exons, a span <= 40 kb, every intron GT..AG, and is not a transcript-assembly `MSTRG.*` gene.
- Each has >= 1 curated isoform that is an exon-skipping variant of the canonical and skips an internal cassette exon.
- 120 genes qualify. A seeded shuffle (seed 1) picks 24: the first 12 for development, the next 12 held out.
  - Development: ZNF217, SALL4, C20orf197, SLCO4A1-AS1, NSFL1C, GSS, AL031674.1, RBM38, FAM182A, AL160412.1, ARFGAP1, SNX5.
  - Held out: C20orf27, RBCK1, YTHDF1, POFUT1, PPDPF, LINC01370, SNRPB, AL139350.1, SAMD10, ZBP1, SPAG4, SNHG17.

**Copies and isoforms.**
- Every copy expresses the template's curated isoforms: the canonical, the MAIN isoform (skips the first cassette exon) and
  up to two more.
- The main isoform is the most expressed, so it is the catalog's representative, which lacks the cassette exon as real
  representatives do.
- 60 reads per expressed copy; famsim errors 0.1% substitution, 0.03% indel.

**Scenarios.** One world per scenario, one family per template; A = the template, B = A + 3% SNPs.

| scenario | copies | what should happen |
|---|---|---|
| NULL | A, B | no call |
| ALLELE | A, B, H = A + 0.4% (absent) | no call (an allele) |
| ABSENT@0.01 / @0.02 | A, B, E = A + 1% / 2% (absent) | E is called |
| BRANCH | A, B, E = A + a novel 150-bp exon + 0.2% SNPs (absent) | E is called: structure differs, sequence barely does |
| OUTSIDE@0.01 / @0.03 | A, B in the catalog, D = A + 1% / 3% in the genome only | no call: the genome explains D |
| TWIN_ALT | A, T = A + 15% SNPs in the cassette exon only, B | no call: T is a present copy |

**Scoring** (`score.py`, register 1294's rules):
- **Detected:** one call holds >= 0.8 of the absent copy's placed reads at purity >= 0.9.
- **False call:** a call whose reads are less than half the absent copy's. In a family with no absent copy, every call is
  false: NULL, ALLELE, OUTSIDE, TWIN_ALT.

## Known before running (recorded now)

**Development** (12 templates per scenario, seed 11), final code:

| scenario | CUR | R0 | S1 | S12 | S123 | S3 |
|---|---|---|---|---|---|---|
| ABSENT@0.01 detected | 8 | 8 | 8 | 8 | 7 | 7 |
| ABSENT@0.02 detected | 12 | 12 | 12 | 12 | 12 | 12 |
| BRANCH detected | 0 | 0 | 0 | 12 | 12 | 0 |
| OUTSIDE@0.01 / @0.03, families with a false call | 8 / 12 | 8 / 12 | 8 / 12 | 8 / 12 | 0 / 0 | 0 / 0 |
| TWIN_ALT, families with a false call | 4 | 4 | 0 | 0 | 0 | 0 |
| NULL, ALLELE, families with a false call | 0 | 0 | 0 | 0 | 0 | 0 |

The ABSENT@0.01 family S123 loses (SALL4) sits at delta: est_div 0.00979, the reads' median `de` 0.0095.

**Fixes made on development, before this prereg:**
- one world per scenario: in a mixed world, the same template's families are copies of each other;
- the novel exon goes into a long intron;
- unaligned ends, and dense error windows in GraphAligner's and minimap2's alignments, count as novel sequence. Both
  aligners force a random insert onto the reference as small mismatches and indels.
- novel segments are linked with short-sequence minimap2 parameters (the asm presets need an alignment score >= 200).

The current method and R0 make the same calls on development (the prototype = Rust, as on real data).

## Bars (held out: templates 13-24, seed 21; 12 families per scenario)

| id | bar |
|---|---|
| H1 | S123 detects >= 11 / 12 ABSENT@0.02 and >= 10 / 12 BRANCH, and ABSENT@0.01 >= CUR's count - 2 |
| H2 | S123 has a false call in <= 3 of the 60 families with no absent copy (NULL, ALLELE, OUTSIDE x2, TWIN_ALT) |
| H3 | S123 detects >= 8 more BRANCH families than CUR, and has >= 10 fewer false-call families than CUR over those 60 |
| H4 | each step is needed where it should be: (a) S1 has fewer TWIN_ALT false-call families than R0 (not measurable if R0 has none); (b) S12 detects >= 8 more BRANCH families than S1; (c) S123 has >= 10 fewer OUTSIDE false-call families (both rates pooled) than S12 |

## Readings (fixed now)

1. **H1-H4 pass.** Claimed, quoted with the counts: "In ideal scenarios built from curated isoforms, the graph-informed O3
   detector finds absent copies that differ only in structure, which the current method cannot. It also stops false calls
   from uncatalogued paralogs and from copies a representative makes look identical, and keeps the current detection of
   sequence-divergent absent copies; each step is needed for its case." Scope: simulation; real data not yet tested.
2. **H1 and H2 pass, H3 or H4 fails.** The detector works; the failing comparison is reported as such.
3. **H1 or H2 fails.** Not claimed.

Reported, not barred: every arm's counts per scenario; recall and purity of the calls; S3 (step 3 alone).

## Results

Added after the runs. Nothing above this heading changes once a statistic has been seen.

## Results (2026-10-15)

Run: code `23441fda` (as `883e54bc`; the prereg adds no code). Held-out worlds: templates 13-24, seed 21. Outputs in
`/mnt/linuxdisk/tmp/o3_graph/heldout/` (one directory per scenario, `score.json`). Each world takes about 30 s: build, map,
the current method, detector.

| held out (12 families each) | CUR | R0 | S1 | S12 | **S123** | S3 |
|---|---|---|---|---|---|---|
| ABSENT@0.01 detected | 8 | 8 | 8 | 8 | **8** | 8 |
| ABSENT@0.02 detected | 12 | 12 | 12 | 12 | **12** | 12 |
| BRANCH detected | 0 | 0 | 0 | 12 | **12** | 0 |
| OUTSIDE@0.01 / @0.03, families with a false call | 8 / 12 | 8 / 12 | 8 / 12 | 8 / 12 | **0 / 0** | 0 / 0 |
| TWIN_ALT, families with a false call | 9 | 9 | 0 | 0 | **0** | 0 |
| NULL, ALLELE, families with a false call | 0 | 0 | 0 | 0 | **0** | 0 |

- Families with a false call, of the 60 with no absent copy: CUR 29, S123 0.
- S123's detected calls: median recall 1.0, purity 1.0, minimum purity 1.0.
- CUR's BRANCH clusters (allele-level SNPs only, est_div < delta) are never `copy` calls.

| bar | value | |
|---|---|---|
| H1 | ABSENT@0.02 12 / 12, BRANCH 12 / 12, ABSENT@0.01 8 against CUR's 8 | pass |
| H2 | 0 of 60 | pass |
| H3 | BRANCH +12; false-call families -29 | pass |
| H4a | TWIN_ALT, R0 9 -> S1 0 | pass |
| H4b | BRANCH, S1 0 -> S12 12 | pass |
| H4c | OUTSIDE, S12 20 -> S123 0 | pass |

### Verdict (2026-10-15)

**Reading 1: claimed, in simulation.** In ideal scenarios built from curated CAT isoforms, the graph-informed O3 detector
(S123) does four things:
- it finds absent copies that differ only in structure: BRANCH 12 / 12, against 0 / 12 for the current method;
- it stops false calls from uncatalogued paralogs: OUTSIDE 20 -> 0;
- it stops false calls from copies a representative makes look identical: TWIN_ALT 9 -> 0;
- it keeps the current detection of sequence-divergent absent copies: 12 / 12 at 2%, 8 / 12 at 1% for every arm.

Each step is needed for its case: step 1 for TWIN_ALT, step 2 for BRANCH, step 3 for OUTSIDE.

**Scope and caveats.**
- These are famsim reads from curated isoforms in small worlds. Step 3's "genome" is the world's genome, not a whole
  genome.
- Real data's other failure modes are not modelled: reads from loci outside the window, other species, fragmented clusters.
  The 1% absent copies sit at delta (0.958%), and every arm misses the same 4 of 12:
  - three are pure clusters labelled allele-level (est_div 0.0069-0.0080);
  - the fourth (PPDPF) is a `copy` cluster holding 47 of the copy's 61 reads, under the 0.8 recall bar.
- Real data is next. On the two-locus spike the steps' real-data counterparts were visible:
  - MCL117's representatives hid every exonic PSV, as in TWIN_ALT;
  - most false `copy` clusters were other loci's reads, as in OUTSIDE.
