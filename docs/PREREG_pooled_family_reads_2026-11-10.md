# PREREG r1625: pooling a family's unassigned reads with copy_phasing (2026-11-10)

User (2026-11-10): "lets register the pooling of a family's ambiguous reads with copy phasing". Written before any lookup.
Scope note: this sits in the read-assignment layer the thesis advisor has not approved; it is an opt-in measurement, changes
no default and no product. AGENTS.md: the aligner's primary is reported beside every figure and is never called beaten at scale.

## Idea in data terms
`copy_assign` judges one molecule at a time against the catalog's copy-vs-copy columns and leaves most rows `ambiguous` or
`tied`. `copy_phasing run` piles every record on each copy, finds the variable columns of the pile (the sample's own variants
included, not only catalog PSVs), ATTRACTS reads that agree at every shared variable column, REPELS reads that disagree, and
assigns a whole cluster. A read with no deciding column of its own can then ride on cluster-mates that have one. What pooling
cannot add: for reads that share only catalog columns, agreeing with an assigned read says nothing new (stated now).

## Step 0 = lookup on existing files (no new run)
- Sample KB3781 (haplotype truth). Rows: `o2mstie/pof_panel.strict3_dump.tsv` + `pof_heldout.strict3_dump.tsv` with status
  `abstained` (r1514 arm, `--tie-ms 10 --primary-on-family`; these are the truth-decidable rows, truth = one catalog copy).
- Pooled call: `~/o2ties_real/kb.A.o2.tsv` (r1483 arm A, `copy_phasing run --semi --guard`, same catalog `kb.fam.copies.tsv`,
  run on 2026-11-01 and never scored against this truth). A call = `assigned` is a catalog copy's tid of the same family.
  `HOPELESS`, `ABSTAIN`, `new:N` and a missing row = no call.
- Aligner on the same rows: the dump's own verdict column.

## Bars
- S (stop): fewer than 100 abstained decidable rows have a `kb.A.o2.tsv` row for (read, family) -> STOP, underpowered.
- B1 reach: pooled calls on >= 10% of the abstained decidable rows.
- B2 precision: right / (right + wrong) >= 0.95 on those calls.
- B3 beside the aligner, on the rows with a pooled call: pooled right >= aligner right and pooled wrong <= aligner wrong.
- Reported without a bar: the same on panel and held-out families separately; calls on `new:N` clusters; how many pooled
  calls sit in clusters that contain >= 1 `copy_assign`-assigned read (the "rides on a mate" mechanism).
- PASS (B1-B3): register a fresh run with `--primary-on-family`-matched input and a second sample with truth (HG002).
  FAIL: closed; `copy_phasing` stays the O3 certificate it is.

## Outcome (2026-11-10, lookup only): FAIL on B3
| rows | abstained decidable | with a copy_phasing row | pooled calls | right / wrong | aligner on the same calls |
|---|---|---|---|---|---|
| all | 1,113 | 620 | 120 (10.8%) | 115 / 5 (0.958) | 118 / 2 |
| panel families | 788 | 366 | 113 | 108 / 5 | 111 / 2 |
| held-out families | 325 | 254 | 7 | 7 / 0 | 7 / 0 |
- S passes (620 rows). B1 PASS (10.8%, bar 10%). B2 PASS (0.958, bar 0.95). **B3 FAIL: 115 right / 5 wrong against the
  primary's 118 / 2 on the same rows.** Pooling calls reads the aligner already places, slightly worse.
- Of the 620 rows copy_phasing saw, 454 are its own ABSTAIN, 45 HOPELESS, 1 a `new:` cluster. Outside the panel it calls 7
  of 325 rows.
- "Rides on a mate": 0 of the 120 calls sit in a cluster of `kb.A.reads.tsv` that also holds a `copy_assign`-assigned read;
  that file lists only the candidate (new-haplotype) clusters, so this count does not measure the mechanism and is not used.
- Closed as registered: no fresh run. `copy_phasing` stays the O3 certificate.
