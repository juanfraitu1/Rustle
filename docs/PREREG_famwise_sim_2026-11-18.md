# PREREG r1797 / r1798 / r1799: a family library with exact provenance, and what exact provenance is worth on it

Date 2026-11-18. Written before any simulated read was aligned. Design:
`docs/superpowers/specs/2026-11-18-family-simulation-design.md` (approved by the user 2026-11-18, including: R is reported
as a bracket with a ceiling; the row 1779 mouse and fly reads may be used for row 1799).
Workspace W = `/anvme/workspace/v138da13-rustle`; products in `W/r1797_famsim_mimic`, `W/r1798_famsim_oracle`,
`W/r1799_simreads_oracle`. Scripts: `bench/famwise_eval/{donor_reads,read_profile,exact_truth,sim_chains,mimic_stats,
read_stats,mimic_bar}.py`.

## What is assumed and what is not

Not assumed: given the copy a read came from, where the aligner puts it. Assumed: how many reads each copy emits (every
count on disk is an aligner placement). R on the simulation is therefore a bracket over three scenarios with a ceiling;
it is never a measurement of the real split. In the observed scenario R near 1 is by construction and is not evidence.

## The libraries (fixed before alignment)

- Members: the 419 truth families (1,459 genes) of `bench/denovo_eval/famtruth.py`; every CAT transcript (5,708).
- Counterparts on `hg002v1.1.mat` / `.pat`: member span with 200 kb flanks, `minimap2 -x asm20 -c --secondary=no`; the
  primary row on the homologous chromosome overlapping most of the member's span (tie: longer target span), MAPQ >= 20;
  exons lifted through it (growth <= 50 bp), isoforms that fail dropped, one needed. Two members of one family with
  lifts on the same haplotype strand overlapping >= 50% of the shorter, and not overlapping on CHM13: no counterpart
  for either on that haplotype.
- Expression from the seen HG002 quarter (`W/r1796_famwise_origin/hg_seen/sel.bam`): reads per member = primaries
  (member on the strand with most overlapping exon bases; tie lower gene id); isoform weights = 1/k per fitting
  isoform (junctions exactly the chain or a contiguous part), canonical isoform when no read fits; even haplotype split.
- Scenarios: **observed**; **to-major**; **to-minor** (near-tied = two best places within 10 on `ms`; members grouped
  by reads with places on both inside that band; a read with a non-member place in the band is left alone; the group's
  near-tied reads go to its member with the most / the fewest (>= 1) near-tied primaries, tie by gene id).
- Read model from the real library on protein-coding, single-transcript, multi-exon, non-member autosomal genes
  (MAPQ 60 primaries on the gene's strand, seeded sample of <= 50 per gene): truncation pairs by transcript length
  class (< 1, 1-2, 2-3, 3-5, >= 5 kb); each read takes one real read's event rates (substitution alleles of HG002
  removed on the full pileup), stochastic rounding; retained introns at the measured per-intron rate. Seed 2026.
- **Reference-drawn** control: the observed counts and the same model, reads drawn from CHM13 itself.
- Alignment: `minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes` to CHM13.
- Exact truth: a read's origin is its place on the contig and strand of, and sharing an exon base with, its source
  member; a read with no record there is `outside` (no home record).

## Row 1797: does the observed library mimic HG002?

Real values are the seen and held-out quarters (`r1785_truth_repair/hg`, `r1796_famwise_origin/hg_held`).
The bar: all five.

- **F.** The ten families with the most near ties in the seen quarter (79% of them): CHM13_G0021052, G0044653,
  G0006774, G0020650, G0049830, G0058367, G0020705, G0000527, G0055907, G0002961 (first member ids). For each, the
  simulated near-tie share and the held-out quarter's have overlapping Wilson 95% intervals or differ by less than a
  factor 1.5. At least 8 of 10.
- **S.** Pooled near-tie share without the first four of those families: same test.
- **W.** Primary off its origin among decidable near-tied reads under the catalogue-free truth (real 47 of 600 and
  47 of 497): the simulated Wilson interval overlaps the held-out quarter's.
- **E.** Share of selected reads with both ends complete against RefSeq (`read_stats.py`, strict rule; real 0.710 and
  0.709): within 0.10.
- **A.** Replay of row 1794 (catalogue-free truth with its caps, tie split, `origin_bam.py`, `region_chains.py` on
  RefSeq, the row's job steps unchanged): R within 0.10 of the held-out values, Rustle 1.04 and StringTie 1.00.

A failed bar: the library is called a scenario in every later sentence, never a mimic.
Controls: the reference-drawn library and KB3781's real statistics (items F, S, W where defined) go through the same
bar. If the reference-drawn library passes, the bar cannot tell the individual's genome from the reference and passing
is weak evidence.
Reported without verdict: calibration statistics (reads per family, read length, events per base, share with no
near-perfect hit on the own haplotypes: real 0.26-0.29); the other statistics of `mimic_stats.py`; families with
disjoint intervals (real against real 1 of 154); simulated against real primaries per member in each scenario.
Truth check (descriptive): of the catalogue-free `decidable` calls, the share at the exact origin; of its ties, the
share containing it; its `outside` calls (all false here by construction).

## Row 1798: what is exact provenance worth?

Per scenario, every mapped simulated read, one alignment: P (primaries), D (all records), O (each read only at its
exact origin; `origin_bam.py`), N (as O with the origin set to the primary's own place). Rustle (`PREFIX.gtf`;
P = `--seed-pool primary` on D, D = defaults, O and N = `--seed-pool all`) and StringTie (`-L`).

- Placement: share of reads with the primary at the origin by lead class (clear, 1-10, 0) and by the catalogue-free
  truth's status.
- Reference chains: multi-exon chains of simulated isoforms carried in full (every junction, >= 20 bases beyond the
  outer ones, no retained intron) by >= 2 drawn reads. Affected chains: those with >= 1 such read whose primary is not
  at its origin.
- Ceiling C = affected chains / affected chains recovered from P. R = affected chains recovered, O over P.
  Also D over P, N over P, by depth class (2-4, 5-19, 20+ full reads), and without the first four families of item F.
- False copies: assembled multi-exon chains on the strand and inside the exon span of a member that emitted no read
  and overlaps no emitting member on that strand.

Readings, per scenario, on affected chains:
- R >= 2 for either tool: the advisor's expectation holds under that split.
- R <= 1.1 for both tools in all three scenarios: no split in the bracket gives several-fold.
- C < 2: two-fold is out of reach in that scenario for any method at a quarter of the real depth.
- Otherwise: between. N over P outside 0.98-1.02: a settings artefact, said so, and R is read against N.
- D above O is reported as the default borrowing reads of a sibling copy.

Prediction: C < 2 and R <= 1.1 in all three scenarios; false copies in P and D under to-major and to-minor (O has none
by construction; the size is the finding).
Beside: RefSeq-exact chains in the regions as in row 1794; O from the catalogue-free truth on the same affected chains.

## Row 1799: mouse and fly (row 1779 reads)

Exact truth from the read name and `W/r1779_simulated_reads/{mmu,dme}/truth.gtf`. Reads whose two best places are
within 10 on `ms` or AS, with all records. P = all primaries; D = the full alignment; O = all primaries with the
near-tied reads placed at their origin. Reference chains = multi-exon simulated transcripts with >= 2 reads; affected,
C, R and readings as row 1798. Fly split by source type (transposon / gene). No mimic bar.

## Already seen

- Fly format test (scout): primary on the source for 20,761 of 40,885 near-tied reads, 825 of 1,248 clear-lead; 20,474
  reads moved, 759 regions; exact band equal to the truth for 12,446 of 41,741. No chain count.
- Primary off its locus: fly 0.85%, mouse 0.48% of reads (scout's quick count).
- Row 1346 (HG00733, same machinery): primaries right 0.619, oracle 0.970. Row 1795: false chains at a silent copy in
  P and D, never in the oracle.
- Recon for the design: counterpart accounting (97.6% of real primaries at members with a lift), isoform-fit and
  near-tie group counts on the seen quarter. No simulated read existed.

## Limits

No reads from haplotype copies without a reference counterpart, none from non-member loci (`outside` class missing);
per-copy expression assumed; novel isoforms absent; collapsed copies, members without a counterpart and isoforms that
do not lift are not simulated; even haplotype split; one read model; indel alleles counted as errors; family reads
only; the family-wise layer's statuses are not scored here.

## Amendments

Made while building, before any output was read (none changes a bar):
- Lead classes of the placement table are on AS between places (the truth scorer's strata); the scenarios' near-tie is on `ms`.
- The replay takes the selected reads' records from the one alignment by name instead of aligning them again.
- `sim_chains.py` tests a false copy on the member's span (first intron start to last intron end inside it).
- Read profile from the whole HG002 library (seed 2026), not the seen quarter.
- Row 1799: near-tied reads by the score table's record rule (two best records within 10 on `ms` or AS); the exact truth of
  every mapped read is computed in name-hash shards; every read counts as carrying its transcript in full (the simulator's
  truncation is not known per read); unmapped reads are left out.
- The truth check is a fourth table of `exact_truth.py --cf`.

Made after the result was read (post hoc, each beside the registered number, never instead of it):
- S: O without the records that share no exon base with the read's source member (`r1798_famsim_oracle/strict.sbatch`).
- The chain tables on `PREFIX.final.gtf` (after the chain score filter).
- C on attainable chains, and the corrected placement shares (final review's recount, `famwise_design_2026-11-18/F_result.md`).

## Result (2026-11-18)

Jobs 4482520-4482523 (libraries), 4482630 (post hoc sensitivity), 4482524 (fly), 4482552 (mouse; the first mouse job, 4482525, lost one truth shard to memory and was cancelled). Final review by two fresh
readers: every count below was rebuilt with independent code, 0 differences (`famwise_design_2026-11-18/F_code.md`,
`F_result.md`).

### The libraries

133,264 reads (observed), 133,226 (to-major), 133,272 (to-minor) from 431 / 424 / 425 members; 136,057 reference-drawn.
1,072 members have a counterpart on both haplotypes, 336 on one, 51 on none (2.0% of the seen quarter's reads, HLA-DRB5
2,425 of the 2,669). Not as registered in effect: 64 members (18,699 reads) emit from one haplotype only, 24 of them
autosomal (HLA-DRB1, EIF3C, EIF3CL, UGT2B17), so "even split" there is one haplotype at full dose; 56 members with one or
two real reads emit none (rounding). EIF3C and EIF3CL have 1,110 near-tied primaries each, so to-major and to-minor send
that group (39% of the near-tied reads) the same way.

### Row 1797: BAR FAILS (S, W, E). The library is a scenario, not a mimic.

| item | held-out HG002 | observed library | verdict | reference-drawn |
|---|---|---|---|---|
| F, ten families' near-tie share | | 8 of 10 | pass | 4 of 10, fail |
| S, pooled near-tie share without the first four | 0.0267 | 0.0522 | fail (factor 1.96) | 0.0418, fail |
| W, primary off origin, decidable near ties | 47 / 497 = 0.095 | 33 / 992 = 0.033 | fail | 0 / 732, fail |
| E, both ends complete | 0.709 | 0.448 | fail | 0.457, fail |
| A, replay of row 1794, R | 1.04 / 1.00 | 1.08 / 1.01 | pass | 1.09 / 1.00, pass |

The reference-drawn control fails F where the observed library passes, so F does tell HG002's genome from the reference.
KB3781: W 5 / 45, pass; the other items are not defined for it.
Without verdict: near-tied reads 8.6% of the selected reads (real 5.6%); read length p50 1,482 (real 1,715); events per
aligned base 0.0057 (real 0.0071); no near-perfect hit on the own haplotypes 0.21-0.26 (real 0.26-0.29); net 1.6% of the
reads have their primary at another member than the one that emitted them, 2.7% net difference from the real primaries.
What failed: the reads are cut shorter than the real ones (the end fractions were measured against the one CAT
transcript of single-transcript genes), and shorter reads tie more often.
Truth check: the catalogue-free truth's `decidable` call is the exact origin for 20,869 of 20,969 reads (99.5%); its tie
contains the origin for 6,772 of 6,800 (99.6%); all 1,367 `outside` calls are false, as constructed.

### Row 1798: BETWEEN. Not several-fold in any scenario; not nothing either. The prediction (R <= 1.1) was wrong.

Placement (observed; 133,262 mapped reads): primary at the origin for 99.05% of reads with a clear lead, 84.0% at a lead
of 1-10, 47.3% at a lead of 0; 5,002 reads (3.75%) off. 607 of those lie inside their own member (fragments of a
retained intron, which share no exon base): 4,395 reads (3.30%) are truly elsewhere, 3,655 at a sibling copy.

Affected chains recovered (Rustle `PREFIX.gtf`, before the chain score filter; StringTie `-L`):

| scenario | reference chains | affected | Rustle P / D / O / N | R | C | StringTie P / D / O / N | R | C |
|---|---|---|---|---|---|---|---|---|
| observed | 1,239 | 116 | 76 / 92 / 95 / 76 | 1.25 | 1.53 | 54 / 54 / 68 / 54 | 1.26 | 2.15 |
| to-major | 1,212 | 102 | 64 / 80 / 83 / 65 | 1.30 | 1.59 | 44 / 44 / 55 / 44 | 1.25 | 2.32 |
| to-minor | 1,225 | 104 | 67 / 80 / 83 / 67 | 1.24 | 1.55 | 49 / 49 / 59 / 49 | 1.20 | 2.12 |
| reference-drawn | 1,223 | 51 | 36 / 39 / 42 / 36 | 1.17 | 1.42 | 24 / 24 / 24 / 24 | 1.00 | 2.12 |

Readings as registered: R >= 2 in no scenario for either tool; R <= 1.1 for both tools in none; C < 2 for Rustle in all
three, not for StringTie. N over P is 1.00-1.02 on reference chains (genome-wide `--seed-pool all` adds 9% non-reference
chains in Rustle). D is below O (92 against 95): Rustle's defaults already get 16 of the 19 chains exact provenance adds.
- By depth (observed): 2-4 full reads, 31 affected chains: Rustle 6 -> 18, StringTie 3 -> 7; 5-19: 36 -> 43 and 27 -> 34;
  20+: 34 -> 34 and 24 -> 27. Several-fold exists only in the 2-4 class.
- Over all reference chains: Rustle 933 -> 955 (1.02), StringTie 737 -> 752 (1.02).
- Without the four largest families: 80 affected, Rustle 53 -> 66 (1.25), StringTie 45 -> 57 (1.27).
- Units: the 116 affected chains sit in 50 members of 26 families; Rustle's 19 gains in 13 members of 11 families.
- Post hoc, strict origin records (S): Rustle 95 / 82 / 82, StringTie 68 / 55 / 59: R unchanged within 0.02.
- Post hoc, after the chain score filter (`PREFIX.final.gtf`): P / D / O = 48 / 53 / 59 (R 1.23), 44 / 54 / 56 (1.27),
  46 / 52 / 57 (1.24); C 2.42 / 2.32 / 2.26; the 2-4 class is 2 -> 4.
- Post hoc, attainable chains (>= 2 records in O that carry the chain; 108 of 116): C 1.42 (Rustle), 2.00 (StringTie).
- Catalogue-free truth instead of the exact one (CF, other read set): Rustle 87 / 75 / 74 of the affected chains,
  StringTie 55 / 44 / 50.

False copies (chains inside a member that emitted no read; 1,016-1,022 such members):

| scenario | Rustle P / D / O / N | after the filter P / D / O | StringTie P / D / O / N |
|---|---|---|---|
| observed | 3 / 41 / 2 / 5 | 0 / 0 / 0 | 3 / 3 / 0 / 3 |
| to-major | 10 / 52 / 2 / 12 | 4 / 5 / 0 | 8 / 8 / 0 / 8 |
| to-minor | 16 / 65 / 2 / 18 | 7 / 10 / 0 | 19 / 19 / 0 / 19 |
| reference-drawn | 0 / 19 / 0 / 0 | 0 / 0 / 0 | 1 / 1 / 0 / 1 |

O's 2 are one place (HLA-G, from secondary records of HLA-A reads inside a merged origin locus); with strict origin
records O has 0 everywhere. Of D's 41, 26 rest on secondary records only; 34 have the exon pattern of a sibling's chain;
16 sit at members that are expressed in the real library and silent only here.
Beside: RefSeq-exact chains in the regions, O over P, 0.94 / 0.89 / 0.88 (Rustle) and 0.97 / 0.92 / 0.88 (StringTie):
below 1 because a chain assembled at a copy that emitted nothing still matches RefSeq. The measure of rows 1794-1796
credits false copies; on real reads that cannot be seen.

### Row 1799, fly: R 1.08 / 1.14 on 46 affected chains of 7,643

2,289,938 mapped reads. Primary at the origin: 99.98% (clear lead), 95.4% (lead 1-10), 31.0% (lead 0; 28,446 reads).
46 of 7,643 reference chains are affected: Rustle 38 / 41 / 41 (P / D / O; R 1.08, C 1.21), StringTie 36 / 36 / 41
(R 1.14, C 1.28). Transposon reads (239,017, where the lead-0 reads are) are single-exon and carry no chain.

### Row 1799, mouse: R 1.13 / 1.11 on 309 affected chains of 35,473

8,404,172 mapped reads. Primary at the origin: 99.88% (clear lead), 97.3% (lead 1-10), 57.5% (lead 0; 69,042 reads).
309 of 35,473 reference chains are affected: Rustle 217 / 244 / 246 (P / D / O; R 1.13, C 1.42), StringTie 201 / 201 / 224
(R 1.11, C 1.54). At 2-4 reads (55 affected chains): Rustle 11 -> 35, StringTie 11 -> 29; at 20+ reads (211): 182 -> 184
and 168 -> 173. Over all reference chains: 27,105 -> 27,135 and 23,583 -> 23,606 (1.00).

### What the three rows say together

Under every split of expression tried, and on two other species' simulated reads, knowing the true copy of every read
recovers 1.1 to 1.3 times the chains that misplaced reads touch, and 1.02 times the family chains overall; the touched
chains are 0.6% (fly), 0.9% (mouse) and 9% (HG002 scenario) of the chains. Several-fold appears only among chains with 2-4 reads. What
exact provenance removes is false copies: chains assembled at a copy that emitted no read.
The scenario has 1.5 times the real near-tie share but a lower share of wrong primaries among the near ties the
catalogue-free truth can decide (3.3% against 9.5%): which way the touched part is off is not known.
None of this measures the real split of expression between tied copies; that stays assumed.

## PREREG r1809: the same library with real read spans (written 2026-11-18, before any read of it was drawn or aligned)

The user's word: rerun the HG002 library with the truncation corrected, as a new row. Workspace `W/r1809_famsim_templates`.
This model was chosen AFTER rows 1797 and 1798 were read; everything in their Result sections counts as already seen.

What row 1797 got wrong: a real read was first fitted to an annotated isoform (as the whole chain or a part of it) and
the drawn read was then cut again with end fractions measured on other genes. Reads were cut twice, single-exon
fragments and unannotated chains were not drawn at all.

The one change. Every simulated read copies one real read of the seen quarter (`W/r1796_famwise_origin/hg_seen/sel.bam`):
the aligned blocks of its record at the member (ends and junctions as aligned; a block ends only at an intron).
- observed: the primary record of every read whose primary is on a member.
- to-major / to-minor: the members of row 1798's rule; a read that moves copies its record with the best `ms` on the
  member it moves to, and stays where it was when it has no record there (counted).
- The blocks are lifted onto a haplotype through the member's alignment row of row 1797 (each block's length within
  50 bp, order kept). The haplotype is drawn evenly (seed 2026, per read) among those where the member has a
  counterpart and the blocks lift; a read that lifts on none is not drawn (counted). Sequence = the lifted blocks,
  spliced, in the sense of the transcript; errors as in row 1797 (one real read's rates, stochastic rounding).
- Gone: isoform weights, the truncation pairs, the retained-intron rate. Unchanged: members, counterparts, the error
  model, the alignment, both truths, every script of the comparison.
- A read carries an isoform in full when its intron chain is the whole chain of an annotated isoform of its member and
  both end blocks are at least 20 bases. Reference chains, affected chains, C and R: as row 1798.
- Reference-drawn control: the same blocks spliced from CHM13 itself.

Bar: row 1797's five items, same thresholds, same held-out quarter. E is expected to pass nearly by construction (the
spans are those of the seen quarter; real against real passes every item), so E passing is not evidence; F, S, W and A
are. The library is called a mimic only if all five pass; otherwise a scenario, as before.

Comparison: row 1798's inputs P / D / O / N per scenario, plus S (O without the records that share no exon base with
the read's source member), registered this time; Rustle `PREFIX.gtf` and StringTie `-L`, `PREFIX.final.gtf` beside.
Readings on affected chains, per scenario, as row 1798: R >= 2 for either tool (S over P) = the advisor's expectation
under that split; R <= 1.1 for both tools in all three = no split gives several-fold; otherwise between; C < 2 = two-fold
out of reach. Depth classes, all reference chains, without the four largest families, false copies: reported as before.

Prediction: the bar passes E and A and fails W; R between 1.1 and 1.5 for both tools in all three scenarios; C < 2 for
Rustle; several-fold only at 2-4 full reads.

Limits beyond row 1798's: every read's structure is a real alignment to CHM13, so an alignment artefact of the real
read becomes true structure here; the origin of a read is still assumed to be where the real library's aligner put it.

### Amendments, row 1809

Post hoc (after the result was read; found by the reviewers):
- S keeps every primary record, also the 14-18 per library that share no exon base with the source member (the rule of
  row 1798's post hoc S), against the wording above. No read loses every record; S false copies are 0 either way.
- "Stays where it was (counted)": the count was not printed by the run. `donor_reads.py templates` prints it now
  (98 reads to-major, 344 to-minor; `templates.tsv` unchanged byte for byte).
- A second table judges "at origin" on the read's own template blocks instead of the best-AS place on the member
  (NBPF1, below). It is a sensitivity, not the registered rule.

### Result, row 1809 (2026-11-19; jobs 4482788-4482791, 31 min each)

A first submission (4482783-86) died in two seconds on a variable clash in the job script; no read was drawn.

Libraries. observed 132,321 reads from 484 members of 228 families (70,850 = 53.5% carry an annotated isoform in
full); to-major 132,304; to-minor 132,350; reference-drawn 136,074. 3,753 templates of 136,074 are not drawn: the
member has no counterpart on either haplotype (HLA-DRB5 2,425, FAM72C 229) or the blocks do not lift (HLA-DRB1 734).
11 members lose all their reads and count as silent. To-major moves 2,466 reads (98 stay), to-minor 3,246 (344 stay).
Median read length 1,700 (real 1,690-1,715; row 1797: 1,482).

**Bar: FAILS on W only; still a scenario, not a mimic.**

| item | held-out real | observed library | | reference-drawn |
|---|---|---|---|---|
| F near-tie share, ten families | | 9 of 10 | pass | 4 of 10, fail |
| S pooled near-tie share | 0.0267 | 0.0244 | pass | 0.0161, fail |
| W primary off origin, decidable near ties | 47 / 497 = 0.095 | 18 / 474 = 0.038 | **fail** | 1 / 210, fail |
| E both ends complete | 0.709 | 0.705 | pass (by construction) | 0.712 |
| A replay of row 1794, R | 1.04 / 1.00 | 1.06 / 1.01 | pass | 1.08 / 1.00 |

Placement (observed): the primary is at the origin for 99.68% of reads with a clear lead, 84.5% at a lead of 1-10,
50.5% at an exact tie (3,739 reads).

**Reading: BETWEEN as registered, and not distinguishable from 1. My prediction (R 1.1-1.5 everywhere) was wrong.**
Affected chains recovered, P / D / O / N / S:

| split | reference | affected | Rustle | S/P | C | StringTie | S/P | C |
|---|---|---|---|---|---|---|---|---|
| observed | 718 | 47 | 27 / 34 / 29 / 29 / 29 | 1.07 | 1.74 | 26 / 26 / 27 / 26 / 27 | 1.04 | 1.81 |
| to-major | 706 | 34 | 20 / 23 / 19 / 21 / 19 | 0.95 | 1.70 | 21 / 21 / 20 / 21 / 20 | 0.95 | 1.62 |
| to-minor | 711 | 41 | 25 / 28 / 24 / 26 / 24 | 0.96 | 1.64 | 19 / 19 / 23 / 19 / 23 | 1.21 | 2.16 |
| reference-drawn | 721 | 23 | 17 / 18 / 18 / 17 / 18 | 1.06 | 1.35 | 11 / 11 / 12 / 11 / 12 | 1.09 | 2.09 |

- R >= 2 nowhere. R <= 1.1 in five of the six cells; StringTie to-minor (19 -> 23) keeps the registered reading at
  "between". Every difference is one to four chains.
- Rustle's N over P is 1.07 / 1.05 / 1.04, outside the registered 0.98-1.02 band (a read's secondary inside its own
  locus counts under `--seed-pool all`), so R is read against N: 29 / 29, 19 / 21, 24 / 26 = 1.00, 0.90, 0.92.
- After the chain score filter (`PREFIX.final.gtf`), P -> S: 22 -> 26, 15 -> 15, 19 -> 19.
- 2-4 full reads (13 / 10 / 14 affected chains): Rustle 2 -> 5, 4 -> 3, 5 -> 6; StringTie 3 -> 3, 4 -> 4, 2 -> 3.
  Not several-fold, and too few chains to say.
- All reference chains: Rustle 519 -> 521, StringTie 433 -> 434 (observed).
- O equals S in every cell.

False copies (chains at members that emit no read), P / D / O / N / S:

| split | Rustle | Rustle after the filter (P / D / S) | StringTie |
|---|---|---|---|
| observed | 0 / 27 / 2 / 2 / 0 | 0 / 0 / 0 | 0 / 0 / 0 / 0 / 0 |
| to-major | 4 / 35 / 2 / 6 / 0 | 4 / 3 / 0 | 13 / 13 / 0 / 13 / 0 |
| to-minor | 11 / 47 / 2 / 13 / 0 | 8 / 13 / 0 | 14 / 14 / 0 / 14 / 0 |

9 / 8 / 9 of D's sit at FAM72C and HLA-DRB5, silent here only because their reads did not lift; O's and N's 2 are
HLA-DRB5. Without those two members D has 18 / 27 / 38. D's false copies are carried by secondary records only in
25 / 27 / 29 cases.

What the reviewers found (two recounts with their own code: 0 differences in 102,393 cells; every read of the
observed, to-major and reference-drawn libraries rebuilt byte for byte with an independent lift; notes in
`famwise_design_2026-11-18/sim_scratch/r1809_*`):

- **One affected chain is an artefact of the truth rule.** NBPF1 holds two near-identical units inside one member.
  Four full reads sit exactly on their template; the truth takes the other unit (AS + 1) as the origin and O moves
  them there. That one chain is lost by O and S for both tools in all three splits. With "at origin" judged on the
  template blocks (post hoc): observed Rustle 26 -> 29 (1.12), StringTie 25 -> 27 (1.08); to-major 19 -> 19, 20 -> 20;
  to-minor 24 -> 24, 18 -> 23 (1.28). No reading changes.
- **Why exact origin sometimes loses a chain** (14 lost cells = 6 chains): NBPF1 (6 cells); four chains that the
  assembler drops once the locus receives all its own reads (CBWD5 in Rustle, 136 -> 412 records on the locus; BOLA2B,
  FAM156A, ANKRD36B in StringTie); one chain that had lived on reads of the sibling copy (POLR2J). None is a missing,
  shorter or filtered origin record: all 1,245 off-origin full reads have a record at the origin and 1,232 carry the
  exact chain there.
- **Affected chains are hit hard, and P still finds them.** A typical affected chain has 48% of its full reads off
  origin (quartiles 0.33-0.63). In the observed split 24 of 47 are held up in P by reads of the other copy carrying
  the identical chain (EIF3C: 366 own + 352 other in P); in the other splits 21 and 23 chains keep five or more full
  reads at home.
- **C is not a reachable ceiling.** Of Rustle's 20 affected chains missing from P (observed), 13 are in no input at
  all; 9 of the 18 missing from O are a contiguous part of a deeper chain.
- **Why 718 reference chains and not 1,239.** Real spans carry fewer annotated isoforms in full: isoforms with >= 2
  full reads 1,268 -> 718. Of row 1798's 116 affected chains, 64 are not reference chains here (54 have no full read).
  Five families hold 31 of the 47 affected chains; EIF3C / EIF3CL holds 59% of the off-origin full reads.
- **Why D beats O in Rustle (34 vs 29, observed).** Of the 5 chains: 2 are the losses above, 2 need reads that are not
  at their origin (POLR2J: the registered "borrowing a sibling's reads"), 1 needs the member's own secondaries.

Prediction scored: bar passes E and A and fails W, right; C < 2 for Rustle, right (1.64-1.74); R 1.1-1.5 in all six
cells, wrong (one cell); several-fold at 2-4 full reads, not seen.

Limits of this row beyond those registered: 47 / 34 / 41 affected chains from 15 families; 16% of the templates have
under half their bases on the member's exons (real reads that are not transcripts of the annotated member); wrong
primaries among decidable near ties are still rarer than in the real library (3.8% vs 9.5%), and which way that moves
R is not known. Not checked: the to-minor library read by read; the assemblers' rules behind "dropped once the locus
gets its reads" (shown by perturbation runs on region BAMs, not by reading the code).

Rows 1798 and 1809 together: with reads drawn complete from every annotated isoform, exact origin recovers 1.2-1.3
times the affected chains; with real read spans it recovers 0.95-1.2 times, on fewer chains. Neither reaches 2.

