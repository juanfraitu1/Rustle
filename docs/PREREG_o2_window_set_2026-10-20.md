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

(none yet)
