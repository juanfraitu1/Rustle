# PREREG — O2 tie-outside registration at the molecule's best score (register row 1220), written before any run

## Finding being tested (row 1220, found by the o3_candidates final review)

`copy_assign`'s §6gz block (`src/bin/copy_assign.rs`, "which tied molecules have a tied placement OUTSIDE every supplied family
UNIT") marks an AS-tied molecule as `tie_outside_catalog` when one of its tied placements lies outside every target of the sweep,
and the mark later demotes an `Assigned` row of that molecule to `Tied` (`copy_assign_pipeline.rs`, the `is_tie_outside` test on the
assigned path). The block's notion of "tied" is the molecule's best AS **within the sweep** (`best_as` is folded over the sweep's
`bam_reads`), while the registry is process-wide: a family-less sweep of another contig registers a molecule from a pair that is tied
only among that contig's records. In the fig4/fig5 human O2 run (`/mnt/linuxdisk/tmp/rustle_figures/o2sim/human`, 290 families on
chr16, 25 contigs swept whole), 20 of the 163 marked molecules are marked ONLY by family-less sweeps of other contigs, and every one
of those 20 marks comes from a tied pair whose alignment score is BELOW the molecule's genome-wide best score. §6gz's intent is "the
molecule could belong to a locus we do not model"; a pair that is not at the molecule's best score does not say that.

## Rule under test

- **Current rule (R0):** a molecule is tie-outside when, in some sweep, a record of it with AS >= the molecule's best AS **in that
  sweep** lies outside every target of that sweep.
- **Proposed rule (R1):** the same, and the outside record's AS must be >= `as_tie_ratio` (0.98) x the molecule's **best AS over
  every record any sweep has seen** — the local best of the consuming sweep, the AS carried by each registered mark, and the
  genome-wide table when `RUSTLE_GTF_SECONDARY_AS_TABLE` is loaded (`global_best_as()`), whichever is largest.
- Mechanism (not a rule): `register_tie_outside_locus` carries the outside record's AS; the demotion test compares each mark's AS to
  the best known at assignment time. Visibility and order are unchanged from R0 (a mark registered by a later sweep demotes nothing in
  either rule); only the score condition changes.
- Switch: `--tie-outside-at-best` (default off for this test; byte-identical unset). Adopted = default on with the escape
  `--no-tie-outside-at-best`; the pre-09-09 escape set gains the flag.

## Substrate and arms

- Dev: the fig4/fig5 human O2 inputs (`figures/_o2.py`, the `o2` and `u2` invocations), run with the same binary under R0 and R1;
  all tables compared (`*.assignments.tsv`, `*.families.tsv`, `*.quant.tsv`, `*.famcn_readonly.tsv`).
- Held-out: the gorilla O2 figure inputs under `/mnt/linuxdisk/tmp/rustle_figures/o2sim/gorilla` (present).

## Pre-registered readouts and rules

- **T1 (scope of the change):** the set of molecules whose `tie_outside_catalog` flips R0 -> R1 is exactly a subset of the
  molecules all of whose outside pairs are below 0.98 x their best score (0 molecules gain the flag; none loses it while holding an
  outside pair at >= 0.98 x best).
- **T2 (effect on assignment):** the number of rows that change status (`Tied` -> `Assigned`) is reported; among them, the fraction
  whose assigned copy is the molecule's best-scoring placement >= 0.95 (R1 restores assignments the aligner already prefers).
- **T3 (held-out):** T1 holds on the gorilla inputs; T2's fraction is reported (no bar: the gorilla run has few contested rows).
- Adopt R1 (default) iff T1 and T2 hold on the dev inputs and T3 holds; otherwise keep R0 and record. The figure tables are
  regenerated and REPRODUCE.md gains a dated note either way. Register row(s) in the house format.

## Not changed

The 0.98 tie ratio, the assignment certificate, R3/R15 of the o3_candidates branch, the `~xchrom~` handling, sweep order.
