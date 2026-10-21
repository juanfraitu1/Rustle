# Pre-registration: O2 as a pool reducer in a 30-point intron-free window, on fresh substrates (2026-10-20)

**User (2026-10-20):**
- "the problem with nearly tied multimappers is that we cannot trust the primary, the primary and all the secondaries are a
  coin toss so that is why I am trying to either define a primary via PSVs and BAM tags or at least reduce the pool of
  possible confounded placements";
- "we could also go past the 10 AS difference from Eichler, if that is too close to the fully tied impossible cases";
- on the proposed test of "intron-free window 30 + the set rule" on fresh reads and simulations: "yes, please proceed".

## The rule under test (frozen at the commit that adds this file; `bench/o3_graph/o2window.py`)

1. **The window.** A read is in the domain when its alignments within 30 points of its best `ms` lie on two or more copies
   of one qualifying family, none off every copy, and minimap2's primary among them.
   - `ms` is minimap2's DP score without the intron charges. One substitution costs 5 points; "within 30" is `ms > best − 30`
     (as Eichler's 10 AS is a slack of 9).
   - Eichler's 10-AS domain is inside this one.
2. **Placements.** Each copy in the window is represented by its highest-`ms` alignment there. An alignment laying 50 or more
   bases on the exons of two copies (`stitched`) is no placement.
3. **The set** (`pareto_set`): the copies in the window that no other copy beats at every substitution column where the two
   differ, over the read bases every placement aligns. Copy r beats copy c when the read matches r and mismatches c at some
   column, and never the reverse. The set drops a copy only when the read contradicts it and never supports it against that
   rival.
4. **The single pick** (the tied copy with the fewest edits over the common window) is reported but carries no bar: on the
   gorilla development families it broke more than it fixed in every band.

## Development evidence (gorilla KB3781, development families; `docs/O2_WINDOW_SWEEP_2026-10-20.md`)

In the 30-point window, under the intron-free answer key: 9,482 scored reads; the set holds the truth on 9,476 (99.94%); it
is one copy for 8,340 (88.0%), right on 8,336 of them (99.95%); minimap2's primary is wrong on 529. The rule's constants
(30 points, 50 bp for a stitch, the Pareto comparison) were chosen there and are fixed here.

## Substrates: two, never seen by this rule

**S1, one genome (the user's regime for O2): simulations heldout24 + heldout25.** 36 worlds: samples A119b and public
× seeds 1401-1403 / 1501-1503 × INTACT, CROSS, SILENT; CHM13-only profiles (`*_v4_ref`), 21 + 14 families. Every copy is in
the reference. Truth = each read's simulated copy (`reads.truth.E0.tsv`). Copies: each world's `o2psv.copies.tsv`; genome:
its `genome.ref.fa`; reads: its `E0.bam`. These worlds were used before only for the corrector's S1 (register 1343) and the
regimen rounds, never for a window or a set rule.

**S2, real reads of another donor: HG002 Kinnex on CHM13** (`/mnt/linuxdisk/tmp/o3_graph/hg002proof/chm13.bam`, 633,939
candidate reads, `minimap2 -x splice:hq -uf --secondary=yes -N 100 -p 0.5`, no `--eqx`; match columns are read from the M
ops against the reference). Families: the 281 qualifying non-PAR families (`famident/qualify_nonpar.json`) of the pool
(`o2psv_dev/pool_copies.tsv`). Truth = the HG002 assembly's placement lifted to CHM13 (`hg002truth`: `truth.tsv` for the
10-AS-tied reads, `truth_untied.tsv` for the rest): class `unique` and `unique_elsewhere` → determined; `own_tie` →
undetermined; `no_counterpart`, `elsewhere` → off-copy; `poor_match` → left out.
- Caveat known in advance: HG002 is not the reference's donor, so every read may carry its own alleles (the O3 confound).
  That is the case where a set that does not guess should earn its keep.
- HG002 was development data for the corrector (register 1343) and the substrate of register 1329; neither used a window
  or a set.

Gorilla and human numbers are never pooled. The gorilla held-out half is not used: its SM392 chimeras motivated the
stitched rule, so it is not fresh for this test.

## Scoring (`o2window.py bars`, one window = 30)

A read is **scored** when its truth is determined and among the window's copies. Reported: scored reads, families, the
truth's statuses, reads whose truth lies outside the window, minimap2's primary right, the set holding the truth, singleton
sets and how many are right, minimap2's primary right on the singleton reads, the single pick's fixes and breaks, and the
same per band (reads first tied within 10, 20, 30 points).

## Bars, per substrate

- **P (power):** at least 300 scored reads from at least 11 families (S1: pooled over its 36 worlds, families distinct).
- **H (holds):** the set holds the truth on at least 99% of scored reads.
- **R (reduces):** the set is one copy for at least 50% of scored reads.
- **Q (quality):** of the singleton sets, at least 98% are right, and at least as many as minimap2's primary is right on
  those same reads.

**Verdict per substrate:** PASS = P, H, R and Q hold; FAIL = P holds and another bar fails; UNDERPOWERED = P fails. S1 is
the primary verdict (one genome); S2 is reported with its own verdict and never pooled with S1.

**A PASS on S1 would claim:** in one genome, O2 reduces the pool of a read tied within 30 intron-free points to the copies
its bases do not rule out; that pool keeps the true copy at least 99% of the time, is a single copy for at least half the
reads, and that single copy is right at least 98% of the time and at least as often as minimap2's primary. **A PASS on S2**
would extend that to a real donor whose copies differ from the reference's.

**It would not claim:** a single correct primary for every tied read (the set may keep several copies); copies absent from
the reference (off-copy: O3); exact ties (undetermined); anything about Eichler's 10-AS reads alone (they are a band here).

## Rules

- S1 and S2 are each prepared and scored once, with the code at this file's commit. Bugs found while running are logged in
  an amendment with the rerun; no constant changes.
- Heavy jobs in the foreground, one at a time, under 10 minutes per command.
- The metric-auditor checks both verdicts on a frozen copy before they are reported.
- Push only when asked.

## Amendments

### Amendment 1 (2026-10-20): two pooling bugs, fixed before any bar was read

Both surfaced on the first S1 scoring and touch no rule, constant or bar:
1. **Simulated read names repeat across worlds** (`s0000001` in every world), so the pooled truth table held 3,779 reads
   instead of 60,345 and pooled calls collided. Truth and calls tables of pooled worlds now key each read by its world
   (`keyed`, `world_tag`, `calls_tag`); a single calls table (S2) is unkeyed.
2. **The `bars` parser** read `--window 30` as two truth files when it followed `--truth`; it is now parsed first.

The first S1 bars output (3,779-read truth) was discarded unread beyond the traceback; S2's bars had not run. The
`prepare` step and its outputs are untouched. Both substrates were then scored once with the fixed code.

## Results (code b4594ea2; `personal/window/{s1,s2}/bars.json`; each substrate scored once)

### S1, one genome (36 simulated worlds, 29 families): **PASS**

| | scored | minimap2 right | set holds the truth | set of one copy | singletons right | minimap2 right on those reads |
|---|---|---|---|---|---|---|
| all, window 30 | 9,549 | 8,165 (85.5%) | 9,549 (100%) | 6,846 (71.7%) | 6,846 (100%) | 6,832 |
| first tied within 10 | 4,347 | 2,979 (68.5%) | 4,347 (100%) | 1,679 (38.6%) | 1,679 | 1,677 |
| first tied within 20 | 3,255 | 3,239 (99.5%) | 3,255 (100%) | 3,242 (99.6%) | 3,242 | 3,230 |
| first tied within 30 | 1,947 | 1,947 (100%) | 1,947 (100%) | 1,925 (98.9%) | 1,925 | 1,925 |

- P holds (9,549 reads, 29 families); H 100%; R 71.7%; Q 100% and 6,846 >= 6,832.
- The single pick: 12 fixes, 0 breaks.
- No read had its truth outside the window; no primary was stitched.
- In the band first tied within 10, the set is one copy for only 38.6% of reads: those are the near-identical copies, where
  the read often matches both copies over its whole length. The set keeps both rather than guessing, and minimap2 is wrong
  on 1,368 of those 4,347 reads.

### S2, HG002 Kinnex on CHM13 (61 families): **FAIL (H, Q)**

| | scored | minimap2 right | set holds the truth | set of one copy | singletons right | minimap2 right on those reads |
|---|---|---|---|---|---|---|
| all, window 30 | 17,881 | 15,798 (88.4%) | 17,353 (97.05%) | 11,355 (63.5%) | 10,852 (95.6%) | 10,836 |
| first tied within 10 | 5,272 | 3,341 (63.4%) | 4,771 (90.5%) | 1,379 (26.2%) | 900 (65.3%) | 896 |
| first tied within 20 | 2,865 | 2,770 (96.7%) | 2,843 (99.2%) | 1,971 (68.8%) | 1,952 (99.0%) | 1,952 |
| first tied within 30 | 9,744 | 9,687 (99.4%) | 9,739 (99.95%) | 8,005 (82.2%) | 8,000 (99.94%) | 7,988 |

- P holds (17,881 reads, 61 families); **H fails** at 97.05%; R holds at 63.5%; **Q fails** at 95.6% (the set's singletons
  are still right at least as often as minimap2 on the same reads, 10,852 >= 10,836).
- The single pick: 18 fixes, 59 breaks.
- Domain reads left unscored: undetermined 1,851 (HG002's own copies tie), poor_match 6,023, off-copy 376, missing 36;
  truth outside the window 93.
- **Both failures sit in the band first tied within 10.** There the set drops the truth on 501 reads (9.5%) and its
  singletons are right on 65%: a read carrying HG002's own allele of its copy, where that allele equals the other copy's
  base, contradicts its true copy at that column and sides with the rival. That is the donor-vs-reference confound (O3's
  case) the pre-registration named in advance.
- Bands 10-20 and 20-30 pass every bar on HG002 too.

## Audit (metric-auditor, 2026-10-20): S1 holds with caveats, S2 holds; the reading does not

- **Every number reproduces to the digit** (independent scorer, `/mnt/linuxdisk/tmp/audit_2026-10-20_o2window/`); the
  frozen code is byte-identical to b4594ea2; Amendment 1 changes only the pooling keys and the parser.
- **The bars do not test the set.** A PSV-free null — minimap2's primary when no rival is within 10 points, every window
  copy otherwise — passes both substrates (recomputed by me as well): S1 held 99.83%, singletons 54.5%, right 99.69%;
  S2 held 99.15%, singletons 70.5%, right 98.79%. The set beats that null on S1 (100% / 71.7% / 100%) and is **worse than
  it on S2 on every bar**. So S1's PASS certifies the 30-point window, not the Pareto set; a null arm was missing.
- **R is carried by the window.** In the Eichler band (first tied within 10) the set is one copy for 38.6% (S1) and 26.2%
  (S2) of reads; the 10-30 bands supply R, and there minimap2 is already right 99.5-100%.
- **The set's gain over minimap2 on S1 is 14 reads** (2 in the <10 band, 12 in 10-20; ZNF658 and GPR89A), 0 the reverse;
  on the other 1,370 reads minimap2 gets wrong the set keeps more than one copy.
- **S2's H failure is mostly one family, FCGR2A:** 383 of 528 dropped truths and 379 of 479 wrong singletons, all one
  pattern (truth .2, minimap2 and set .1; 381 are `unique_elsewhere` with donor_nm 0-1 and donor margin 274-321). Whether
  HG002's locus carries .1 alleles (the O3 confound) or the chain lift lands on the wrong paralog of the FCGR2A/2C/2B
  region is not established. Without FCGR2A, S2's "all" would pass; its <10 band would still fail (held 97.5%, singletons
  right 90.0%). SMN2 accounts for 9 drops.
- **"Holds the truth" is conditioned on the window:** 93 S2 reads whose truth lies outside the 30-point copies are
  excluded; 87 of them are singletons (the primary), all wrong. Counting them, S2 H = 96.55%.
- **Unscored singletons:** the set is one copy on 774 of the 1,851 reads HG002's own genome calls undecidable and on 4,288
  of the 6,023 poor-match reads. Q is on determinable reads only.
- S2's 10-20 band passes by 6 reads (22 failures against 28 allowed); minimap2 alone passes the 20-30 band.
- Stitched primaries: S1 0, S2 59 (11 among the dropped truths). S1 reads are not error-free (16% NM 0, mean de 0.0023),
  so H = 100% is favoured, not forced.
- S1 and S2 are fresh for this rule (constants from gorilla); but r1329 had shown that 73% of HG002's ties are donor ties
  before the bars were chosen.
- **Wording not to use:** "never drops the truth"; "the set does not guess"; "O2 reduces the pool" as a claim distinct
  from "abstain when a rival is within 10 points"; "the 10-30 bands pass on HG002 too" without the null; "the failures
  are the donor's own alleles" without naming FCGR2A and its unverified lift.

**What a next test needs:** a pre-registered null arm (the rule above), scored on the Eichler band alone, where the set
must beat the null on held truth and singletons; and the FCGR2A lift settled by aligning CHM13's FCGR2A.1 / .2 to HG002's
haplotypes.
