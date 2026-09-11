# PREREG: junction-aware greedy locus growth for de novo (reads-only) O1 boundary inference

**Written before implementation or measurement.** Copy committed to `docs/` with an md5 recorded below
before any run.

## Motivation

The de novo (annotation-free) O1 path (`gw_family_catalog`'s default, reads-only mode via
`denovo_pipeline.rs`/`family_detect.rs::collapse_loci_span_aware`) inherits a locus's boundary from ONE
picked representative transcript (`pick_locus_rep`). Diagnosed 08-03/04 (`project_locus_boundary_extent`):
single-exon-stub reps measure 0.29x true gene size vs 0.75x for spliced reps — "the size problem and the
single-exon-rep problem are one problem." Five independent fixes were tried and all failed end-to-end
because of a structural fact established 08-04 and never overturned: **boundaries and edges are coupled
through sequence sharing** — enlarging a locus toward its true size increases shared sequence with
neighbouring paralogs, creating more homology edges and merging components.

None of the five prior attempts used BOTH a gap signal and a relative-depth signal together as the
stopping rule. This session's `--min-boundary-fraction`/`--min-boundary-gap` (copy_assign.rs, §6hs)
validated exactly that compound signature for isoform-boundary outlier detection *within* an already-known
copy — depth alone failed there for the same reason `core-as-boundary` failed here (can't distinguish a
genuinely shallow real region from a readthrough tail); requiring gap AND low relative depth together
correctly separated ~30 smooth alternative-terminus points (kept) from isolated 5.5-139kb jumps (cut).

## Mechanism (NEW — not previously tried in this form)

Instead of summing/pooling every overlapping read into an aggregate span and then trimming (what
`locus_confident_extent`/`RUSTLE_LOCUS_DE_EXTENT` and `locus_core_bp` do), grow a locus's boundary
incrementally outward from its seed (the union of the junction-collapsed member transcripts' own exon
extents, via the existing `collapse_loci_span_aware_with_members`) one read-block at a time:

1. Seed = union of all member transcripts' own (start, end), i.e. the already-junction-validated
   isoform union — NOT the single picked rep's span.
2. Sort all primary-read genomic footprints on the locus's chromosome by position.
3. Walk outward (5' and 3' independently). At each step, before folding in the next read/block beyond the
   current boundary, compute: (a) the positional gap from the current boundary to that block's start,
   relative to the LOCAL median inter-read-start spacing already observed inside the current boundary
   (never a fixed global bp constant — the lesson from `--min-boundary-gap`'s first-pass failure); (b) the
   read count crossing that gap as a fraction of local depth just inside the boundary.
4. Stop growing in that direction the moment BOTH (a) the gap exceeds `k` times the local median spacing
   AND (b) the crossing read fraction is below `min_boundary_fraction` — mirroring the exact compound test
   already validated, reusing its already-calibrated constants (`min_boundary_fraction=0.10`) as the
   starting point rather than inventing new ones.
5. Otherwise fold the block in, update local statistics, continue.
6. No cross-locus homology reference is used anywhere in this rule — unlike the refuted `any-locus
   homology bounding` attempt, growth only ever consults this locus's own local read geometry.

Implemented as a new function `locus_growth_extent(bam_reads, transcripts, members_per_rep, k, min_frac)`
in `denovo_pipeline.rs`, parallel in shape to `locus_confident_extent`. New flag
`RUSTLE_LOCUS_GROWTH_EXTENT` (unset = off = byte-identical), consumed at the same site
`locus_confident_extent` is (mutually exclusive — if both are set, `RUSTLE_LOCUS_GROWTH_EXTENT` takes
precedence and a warning is printed). `k` and `min_frac` are read from
`RUSTLE_LOCUS_GROWTH_K` (default 3.0) / reuse `RUSTLE_MIN_BOUNDARY_FRACTION`-style default 0.10 if unset.

## Evaluation — TWO substrates, each scored on BOTH axes

**Human**: `winloci_data/A119b.t2t.bam` (confirmed full alignment, no region/downsample subset — direct
`minimap2 -ax splice:hq` to `chm13v2.0.fa`), chr1+chr15, truth = `bench/soto/80_fams.chr.bed` /
`80_fams.gene_preferred.bed` (RefSeq/CAT gene spans via `bench/soto/size_vs_refseq.py`), partition quality
via `bench/soto/partition_score.py` against Soto's 83 families.
⚠ Explicitly NOT `winloci_data/soto_reads.bam` — confirmed built with `samtools view -L 80_fams.chr.bed`,
a region-subset BAM with no flanking context, which would make any "does it correctly stop" measurement
vacuous (growth would look successful only because there is nothing outside the truth windows to grow
into). This is a substrate change from the original 08-03/04 experiments, so absolute numbers are not
directly comparable to that historical table; only the DELTA against a same-BAM shipped-baseline re-run
is comparable, and both will be reported.

**Gorilla**: `fibroblasts/GCA_029281585.2_flnc_mm.bam` (confirmed full alignment, no subset/downsample,
direct `minimap2 -ax splice:hq` to `GGO.fasta`), truth = `o1_oracle/npip31.regions` (31-locus NPIP
projection) plus the current default `mcl_families` catalog as the partition reference (today's fresh
274-cluster run, `scratchpad/rna_bp1_p9_repro.clusters.tsv`).

Metrics, both substrates: (1) size — median ratio to truth span, % in [0.5, 2.0]x (in-band fraction, never
the median alone — standing metric-trap rule), % <0.5x, % >2x; (2) partition/end-to-end — `partition_score.py`'s
member_recall / homogeneity / completeness / over_merge / split counts on human; NPIP recall (X/31) +
whether the fresh catalog's over-merge/fragmentation pattern (today's finding: 31/31 present but split
across 4 clusters) changes, on gorilla.

## Predictions, pre-committed

- **P1 (size)**: in-band fraction improves over the shipped representative-span baseline on the SAME BAM,
  without the over-extension blowup seen in `any-locus`/`confident-read-extent` (their >2x fraction hit
  26-37%). Pass: in-band fraction higher than shipped AND >2x fraction stays under 15%.
- **P2 (partition quality, the decisive test)**: does not reproduce the `RUSTLE_LOCUS_DE_EXTENT` collapse
  (F1 0.704→0.401, one family ballooning to 94 loci) or the `any-locus` regression (over_merge count roughly
  doubling). Pass, human: `partition_score.py`'s over_merge count does not increase by more than 20% over
  the shipped baseline, homogeneity does not fall by more than 5 points. Pass, gorilla: NPIP recall stays
  ≥31/31 and does not INCREASE cross-cluster fragmentation (today's 4-cluster split does not become worse).
- **Refutation criteria, stated in advance**: if in-band fraction improves but P2 fails on either substrate,
  this is refuted for the same structural reason as the prior five attempts (boundaries↔edges coupling),
  and the honest conclusion is that this coupling is not escapable by a smarter STOPPING rule alone — the
  next step would need to touch the edge-building rule itself, not the boundary rule.

## Explicit scope confirmation

No bipartite matching or facility-location step appears anywhere in the growth rule above (standing rule:
"NEVER build loci with facility location / bipartite matching" — construction only, `rep_quality.py`'s
own docstring: bipartite matching is a MEASUREMENT tool, "if it ever influences what a locus or family IS,
the method is broken"). Growth uses only local read-depth/gap geometry.

## Addendum (2026-09-10, later same session) — the promised md5, and status

This file's own text above (everything before this addendum) was never edited after being written and
committed (`d690715`); its md5 as pre-registered — and as of this addendum — is
**`8d5a63fb70875996d599c9f4218cd83c`**. The line-3 promise ("md5 recorded below") was left unfulfilled
until now, found by a 2026-09-10 audit; recorded here as an addendum rather than by inserting text into
the frozen body above, so the pre-registration itself stays untouched.

**Status: implemented, NOT measured.** `RUSTLE_LOCUS_GROWTH_EXTENT`/`locus_growth_extent()` shipped in
`denovo_pipeline.rs` (env-gated off, byte-identical when unset). A human chr1+chr15 baseline run was
started to gather the P1/P2 numbers above and killed mid-execution to prioritize the O2 loose-ends batch
(§6hw onward in `docs/o1_ledger.md`) — parked, not abandoned, not refuted. No P1/P2/refutation number in
this document has been measured yet; do not quote this mechanism as validated or as failed.
