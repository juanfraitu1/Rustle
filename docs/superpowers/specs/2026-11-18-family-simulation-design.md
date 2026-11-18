# Reads from the individual's own gene copies: a family simulation with exact provenance (design)

Date 2026-11-18. Status: draft 2 for approval (draft 1 was read by three independent reviewers; their findings are in
`W/famwise_design_2026-11-18/V_{builder,validity,codefit}.md` and are folded in). Follows
`2026-11-18-family-wise-assembly-design.md` (rows 1794-1796). W = `/anvme/workspace/v138da13-rustle`.
Rows claimed: 1797 (the simulation and whether it mimics), 1798 (what exact provenance is worth on it), 1799 (the same
on the simulated mouse and fly reads of row 1779).

## Intent

Rows 1794-1796 measured on real reads what true provenance is worth to assembly (x1.04 Rustle, x1.00 StringTie). Two
things limit that reading: the truth cannot judge a quarter of the near-tied reads, and where a read fits several copies
it cannot say which one emitted it. The user asked whether the way real families turn into reads can be simulated and
made to mimic real data. Here reads are drawn from HG002's own two haplotypes at real family genes, with real isoforms,
observed expression levels and a read profile measured on the real library, and aligned to CHM13 as real reads are.
Every read then has a known copy of origin.

## What the simulation can and cannot decide (read this first)

- **Free of assumption:** given the copy a read came from, where the aligner puts it. This depends on HG002's real
  sequence, the isoform, truncation and errors. It covers the reads the real truth cannot judge.
- **Assumed:** how many reads each copy emits. No measurement exists; every count on disk is an aligner placement.
  At identical copies the observed primaries are a fair toss (EIF3C / EIF3CL: 1,110 / 1,093), so a simulation that
  takes them as expression gives primary-only and true-origin inputs exchangeable reads, and R near 1 by construction.
- Therefore R is read as a **bracket over the unknown split** (three scenarios), together with a **ceiling** that no
  oracle can exceed. It is never presented as a measurement of the real split.

## Goals

1. A simulated HG002 family library whose every read has an exact origin (member gene, haplotype, isoform).
2. Row 1797: a registered check that it mimics the real library, including a replay of the real experiment of row 1794.
3. Row 1798: placement accuracy by read class with the exact truth; the ceiling and R in three expression scenarios.
4. By-product: the first measurement of how often the catalogue-free truth (rows 1775 / 1785) is right.
5. Row 1799: the same comparison on the simulated mouse and fly reads of row 1779 (the user's "then A").

## Non-goals

No assembler change, no default change. No instrument model (PBSIM3 + ccs). No reads from haplotype copies without a
reference counterpart. No template switching, internal priming, chimeras or allele-specific expression. The family-wise
layer's statuses are not scored on the simulation in this round.

## Anti-goals

- A simulation tuned until it gives a wanted answer: the read profile and the expression come from fixed measurements
  made before any simulated read is aligned, and rows 1797-1799 are registered before that.
- Fitted statistics passed off as validation: each statistic is declared calibration or validation.

## Approaches

- **A. The simulated mouse and fly reads of row 1779.** Real genomes and isoforms, source in the read name, aligned
  already. Reads come from the reference itself, errors uniform 1.6%, no truncation. Runs second (row 1799).
- **B. Reads from the individual's own haplotypes, aligned to the reference.** Chosen by the user; runs first.
- **C. An outside simulator for the sequencing step.** Only replaces B's read model; not installed; left out.

## What exists and is reused (imported, not edited)

- `bench/o3_graph/donorsim.py`: `read_paf`, `lift_exons`. `bench/o3_graph/realworld.py`: `splice`, `retain`,
  `draw_truncation`, `draw_edits`, `apply_edits`. `bench/o3_graph/calibrate.py`: `cigar_events`, `transcript_ends`,
  `retained_introns`, `classify_mismatches`. Their 52 tests pass on Fritz.
- `bench/o2_truth/neutral_truth.py` (`select`, `score`, `records`, `loci_of`), `tie_split.py`,
  `W/r1790_family_ties/family_ties.py`.
- `bench/famwise_eval/{origin_bam,region_chains,records}.py`; the job script of row 1794, copied into the new row
  directories with the directory as an argument (the original writes into its own row).
- Scout prototypes under `W/famwise_design_2026-11-18/sim_scratch/{mimic_spec,machinery}/` become tested scripts.

## Design

### 1. Copies and isoforms

- Members: the 419 truth families (1,459 member genes) of `bench/denovo_eval/famtruth.py`. Isoforms: every transcript of
  those genes in the CAT annotation (5,708); abutting exons merged.
- Counterpart on each haplotype (`hg002v1.1.mat`, `.pat`): the member's CHM13 span with 200 kb flanks, aligned with
  `minimap2 -x asm20 -c --secondary=no`. Among primary rows on the homologous chromosome, the row overlapping most of
  the member's own span (tie: the longer target span); MAPQ >= 20. Each isoform's exons go through that row
  (`lift_exons`, growth <= 50 bp); isoforms that do not lift are dropped on that haplotype and their weight goes to
  those that do; a copy needs one. A '-' row flips the transcript strand.
- Collapsed lifts: two members of one family whose lifts lie on the same haplotype strand and overlap by >= 50% of the
  shorter, and which do not overlap on CHM13, both have no counterpart on that haplotype (the repo's tested rule).
  Reviewer's measurement: 2.3% of real primaries, 39 of 60 such pairs sequence-equivalent, 13 not.
- Measured (recon, flank 200 kb, under a stricter all-isoforms rule): members holding 97.6% of the real family
  primaries lift onto at least one haplotype and 83.5% onto both; with a 20 kb flank 94.3%. HLA-DRB5 (13,077 primaries
  in the full library) is nearly the whole loss; HG002 appears to lack it.
- Deliverable: an accounting table per member and haplotype (row length, coverage of the member, `de`, MAPQ, isoforms
  lifted, and the reason when there is no counterpart).

### 2. Expression (calibration; assumed, see above)

Source: the seen HG002 quarter, all records (`W/r1796_famwise_origin/hg_seen/sel.bam`, 136,074 reads).

- A primary belongs to the member on its strand with the most overlapping exon bases (tie: lower gene id).
  Reads per member = its primaries. Depth = one quarter of the library.
- Isoform weights within a member: each spliced read whose junctions are exactly an isoform's chain or a contiguous
  part of it adds 1/k to each of the k isoforms it fits. Single-exon reads and reads fitting none carry no weight
  (27% of spliced primaries fit no annotated isoform: novel chains are not simulated). No weighted read: the canonical
  isoform (most exons, then longest span, then id).
- Haplotypes: even split over the haplotypes with a counterpart. Cell counts rounded (`round`).
- Scenarios, differing only in who emits the near-tied reads. Near-tied = two best places within 10 on `ms`.
  Members are grouped when a read has places on both within that band; a read with a non-member place in the band
  is left alone.
  - **observed**: every read at its primary's member.
  - **to-major**: a group's near-tied reads all go to its member with the most near-tied primaries (tie: gene id).
  - **to-minor**: all go to its member with the fewest (at least one).
  A member without a counterpart emits nothing; its reads are not moved to a sibling.
- For each scenario the simulated primaries per member are compared with the real ones. A scenario the real primaries
  contradict is said to be contradicted.

### 3. Read model (calibration, on genes outside the families)

Measured once on `W/reads/HSA/HG002_kinnex.t2t.bam`: protein-coding genes with one annotated transcript and two or more
exons, not family members, on the autosomes; primaries at MAPQ 60 on the gene's strand; a seeded sample of at most 50
per gene.

- Truncation: (5' fraction, 3' fraction) pairs by transcript length class (< 1, 1-2, 2-3, 3-5, >= 5 kb; an empty
  class takes the nearest). A read draws a pair from its isoform's class; redrawn up to 10 times if under 200 bases
  would remain. No extra end jitter.
- Errors: each simulated read takes the event rates of one real profiled read (substitutions, insertions, deletions
  per base), scaled to its length with stochastic rounding and no second Poisson draw; positions, homopolymer share
  and lengths as `draw_edits` places them. Substitutions that are HG002's own alleles are removed first on the full
  pileup of the gene (seen in >= 3 reads and >= 20% of the coverage). Indel alleles stay in the error rate.
- Retained introns: per-intron probability as measured. Strand: sense. FASTA output; names carry nothing and no `|`.
- One seed (2026) shared by all scenarios.
- **Control library "reference-drawn"**: the same model and the observed counts, but reads drawn from CHM13 itself.
  It is the negative control of row 1797.

### 4. Alignment and truths

- House command (`minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`) to CHM13; sorted, indexed.
- **Exact truth** (new script): a read's origin is its place on the contig and strand of, and sharing an exon base
  with, its source member; status `decidable`, or `no home record` when the read has no record there. Written for
  every mapped simulated read in the `reads.tsv` format, so `origin_bam.py` reads it unchanged.
- **Catalogue-free truth**: reads with a primary on a member (`neutral_truth.py select`, no sampling, then
  `samtools view -N`), then `score` with the repaired rule and its caps, family table, tie split: the real pipeline.
  Simulated reads whose primary lies on no member are counted per member.

### 5. Row 1797: does it mimic? (observed scenario; real = seen and held-out HG002 quarters)

Calibration statistics (reads per family, read length, error rates, share of reads with no near-perfect hit on their
own haplotypes: real 0.26-0.29) are reported without a verdict.

The bar, five items, all must hold:

- **F. Families.** For the ten families holding 79% of the seen quarter's near ties (named in the registration), the
  simulated near-tie share and the held-out quarter's have overlapping 95% intervals or differ by less than a factor
  1.5. At least 8 of 10.
- **S. Near-tie share without the four largest of those families** (fixed from the seen quarter): same test.
- **W. Primary off its origin among decidable near-tied reads** (catalogue-free truth; real 47 of 600 and 47 of 497):
  95% intervals overlap.
- **E. Completeness of simulated family reads** (out of sample: truncation is measured on other genes): both ends
  complete within 0.10 of the real 0.71.
- **A. Replay of row 1794** with the catalogue-free truth and the unchanged scripts: R within 0.10 of the real value
  for both tools (Rustle 1.04, StringTie 1.00).

The factor 1.5 is tied to the question: the dispute is over several-fold. If the bar fails, the library is called a
scenario, never a mimic of HG002, in every later sentence.

Controls on the bar itself: the reference-drawn library and KB3781's real statistics go through the same five items.
If the reference-drawn library also passes, the record says the bar cannot tell an individual's genome from the
reference, and passing is weak evidence.

Beside, no verdict: every other statistic of the scouting table (pooled near-tie share, lead 0 share, status mix,
ties that fit several copies, secondaries per read: all largely fixed by the fitted counts), the number of families
with disjoint intervals (real against real: 1 of 154), regions, and D over P.

Truth check (descriptive, with intervals; S2 complete, S1 the truth's own sample): of the catalogue-free truth's
`decidable` calls, the share at the exact origin; of its ties, the share containing it; of its `outside` calls, all are
false here by construction and are counted.

### 6. Row 1798: what is exact provenance worth here?

On every mapped simulated read of each scenario (no selection), from one alignment:
P (primaries only), D (all records, the default), O (each read only at its exact origin, moved if needed),
O-null (built like O with the origin set to the primary's own place). Rustle and StringTie on each.

- **Placement** (assumption-free part): share of reads with the primary at the origin, by lead class (clear, 1-10, 0)
  and by what the catalogue-free truth said (decidable, tie, could not judge).
- **Reference chains**: simulated isoform chains (multi-exon) carried in full by at least 2 drawn reads, at their
  CHM13 coordinates. **Affected chains**: those with at least one such read whose primary is not at its origin.
- **Ceiling** C = affected chains / affected chains P recovers. No input can give R above C.
- **R** = affected chains recovered, O over P; also by depth class (2-4, 5-19, 20+ reads), with and without the four
  largest near-tie families, and D over P.
- **False copies**: assembled multi-exon chains on the strand and exon span of a member that emitted no read, per input.
- O-null must give R between 0.98 and 1.02; otherwise the comparison carries a settings artefact and says so.

Readings (Rustle `PREFIX.gtf` and StringTie; affected chains; per scenario):
- R >= 2 for either tool: the advisor's expectation holds under that split.
- R <= 1.1 for both tools in all three scenarios: no split in the bracket gives several-fold.
- C < 2 in a scenario: two-fold is out of reach there for any method, at a quarter of the real depth.
- Anything else is reported as between. D above O is reported as what it is: the default borrowing a sibling's reads.

Prediction: C < 2 and R <= 1.1 in all three; under to-major and to-minor, false copies in P and D (O has none by
construction; the size is the finding). In the observed scenario R near 1 is by construction and is not evidence.

Beside: RefSeq-exact chains in the regions (`region_chains.py`, as row 1794); O built from the catalogue-free truth,
counted on the same affected chains.

### 7. Row 1799: mouse and fly

Exact truth from the read name and `truth.gtf`. Near-tied reads (two best places within 10 on `ms` or AS) are taken
with all their records; O = the primaries of all reads, with the near-tied reads placed at their origin; P = all
primaries; D = the full alignment. Reference chains = simulated transcripts with >= 2 reads; same ceiling, R and
readings; fly split by source type (transposon / gene). No mimic bar: these reads are drawn from the reference.

## Already seen (disclosed in the registration)

- Fly, a scout's format test: primary on the source for 20,761 of 40,885 near-tied reads; 825 of 1,248 clear-lead;
  20,474 reads moved, 759 regions; exact band equal to the truth for 12,446 of 41,741. No chain count was seen.
- Share of reads with the primary off its locus: fly 0.85%, mouse 0.48% (scout's quick count).
- Row 1346 (same machinery, HG00733): primaries right 0.619, oracle 0.970.
- Row 1795: false chains at a silent copy in P and D, never in the oracle.

## Limits

- No reads from haplotype copies without a reference counterpart and none from non-member loci: the `outside` class
  (about 4% of located near ties) is missing. One large family (G0006774) has near ties the simulation may not produce.
- Per-copy expression is assumed; the scenarios bracket it and do not test it.
- Novel isoforms are absent; every simulated chain is annotated. Depth classes show where this matters.
- Collapsed copies, members without a counterpart (2.4% of real primaries) and isoforms that do not lift are not
  simulated. Even haplotype split. One read model for all genes. The library holds family reads only.

## Left to the builder

Script layout under `bench/famwise_eval/` (driver for families / counts / counterparts / reads; read profile; exact
truth; mimic statistics; chain table), job scripts in the row directories, tests first, no fixed paths.

## Calls made by the author

- R reported as a bracket with a ceiling, not as a measurement (the reviewers' main point).
- Three expression scenarios and how each is built.
- Collapsed copies dropped rather than given two origins.
- Flank 200 kb; the row-covering-the-member rule; partial isoform lifts kept.
- Expression from the seen quarter; the unread half of the library stays unread.
- Read profile from protein-coding single-transcript genes; errors transplanted per read.
- The five bar items, the factor 1.5, and the reference-drawn control.
- Affected simulated chains as the count; O moves every read to its origin, not only near-tied ones.
- The layer's statuses are not scored here.
